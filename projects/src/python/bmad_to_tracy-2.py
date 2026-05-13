#!/usr/bin/env python3
"""
bmad_to_tracy2.py

Translate the Bmad lattice subset used by racetrack_FOBOD_tuned.bmad into a
Tracy-2 style .lat file.

Supported Bmad constructs:
  - parameter[p0c] = ...
  - Drift, Quadrupole, SBend, Wiggler, RFCav element definitions
  - simple attribute lists, including Bmad continuation lines ending in '&'
  - Line = (...) definitions
  - line multipliers such as 32*FULL_CELL
  - reversed line references such as -HALF_CELL, emitted as inv(half_cell)
  - use, RING

Notes:
  - Bmad K1 maps to Tracy B_2 for quadrupoles and combined-function bends.
  - Bmad SBend angle maps to Tracy Phi = angle*180/pi.
  - Bmad wigglers are emitted as drifts by default, plus a commented CWIGGLER
    candidate when period/field data are present. This mirrors common Tracy-2
    workflows where the undulator is optically treated as a drift unless a
    local CWIGGLER model is explicitly enabled.
  - If an RFCav has `harmon` but no RF frequency, the script computes
    Frequency = harmon*c0/circumference from the selected use/ring line.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

C0 = 2.99792458e8


@dataclass
class Element:
    name: str
    kind: str
    attrs: Dict[str, str] = field(default_factory=dict)
    source: str = ""


@dataclass
class LineToken:
    name: str
    mult: int = 1
    inverse: bool = False
    raw: str = ""


@dataclass
class LineDef:
    name: str
    tokens: List[LineToken]
    source: str = ""


@dataclass
class Lattice:
    parameters: Dict[str, str] = field(default_factory=dict)
    elements: "OrderedDict[str, Element]" = field(default_factory=OrderedDict)
    lines: "OrderedDict[str, LineDef]" = field(default_factory=OrderedDict)
    use_line: Optional[str] = None
    warnings: List[str] = field(default_factory=list)


def strip_bmad_comment(line: str) -> str:
    """Remove Bmad ! comments. This subset does not support quoted ! strings."""
    return line.split("!", 1)[0]


def logical_bmad_lines(text: str) -> List[str]:
    """Join Bmad continuation lines ending in '&' after comment removal."""
    out: List[str] = []
    buf = ""
    for raw in text.splitlines():
        no_comment = strip_bmad_comment(raw).rstrip()
        if not no_comment.strip():
            if buf.strip():
                out.append(buf.strip())
                buf = ""
            continue
        cont = no_comment.endswith("&")
        if cont:
            no_comment = no_comment[:-1].rstrip()
            buf += no_comment + " "
        else:
            buf += no_comment
            if buf.strip():
                out.append(buf.strip())
            buf = ""
    if buf.strip():
        out.append(buf.strip())
    return out


def split_top_level_commas(s: str) -> List[str]:
    """Split on commas not inside (), [], or {}."""
    parts: List[str] = []
    start = 0
    depth = 0
    pairs = {"(": ")", "[": "]", "{": "}"}
    closers = set(pairs.values())
    for i, ch in enumerate(s):
        if ch in pairs:
            depth += 1
        elif ch in closers:
            depth = max(0, depth - 1)
        elif ch == "," and depth == 0:
            parts.append(s[start:i].strip())
            start = i + 1
    tail = s[start:].strip()
    if tail:
        parts.append(tail)
    return parts


def parse_attrs(parts: Iterable[str]) -> Dict[str, str]:
    attrs: Dict[str, str] = OrderedDict()
    for p in parts:
        if not p:
            continue
        if "=" not in p:
            # Keep bare flags rather than silently dropping them.
            attrs[p.strip().lower()] = "T"
            continue
        k, v = p.split("=", 1)
        attrs[k.strip().lower()] = v.strip()
    return attrs


def parse_line_expr(expr: str) -> List[LineToken]:
    """Parse the inside of Bmad Line=(...)."""
    tokens: List[LineToken] = []
    for raw in split_top_level_commas(expr):
        t = raw.strip()
        if not t:
            continue
        inverse = False
        if t.startswith("-"):
            inverse = True
            t = t[1:].strip()

        mult = 1
        m = re.match(r"^(\d+)\s*\*\s*(.+)$", t)
        if m:
            mult = int(m.group(1))
            t = m.group(2).strip()

        # Handle the rare form N*-NAME as well.
        if t.startswith("-"):
            inverse = not inverse
            t = t[1:].strip()

        tokens.append(LineToken(name=t, mult=mult, inverse=inverse, raw=raw))
    return tokens


def parse_bmad(path: Path) -> Lattice:
    lat = Lattice()
    text = path.read_text(encoding="utf-8-sig")
    for line in logical_bmad_lines(text):
        # parameter[p0c] = 1.25e9
        m = re.match(r"^parameter\s*\[\s*([^\]]+)\s*\]\s*=\s*(.+)$", line, re.I)
        if m:
            lat.parameters[m.group(1).strip().lower()] = m.group(2).strip()
            continue

        # bmad_com[...] = ... and other global knobs are intentionally ignored.
        if re.match(r"^\w+\s*\[", line):
            continue

        # use, RING
        m = re.match(r"^use\s*,\s*([A-Za-z_]\w*)\s*$", line, re.I)
        if m:
            lat.use_line = m.group(1).strip()
            continue

        if ":" not in line:
            lat.warnings.append(f"Ignored unrecognized line: {line}")
            continue

        name, rest = line.split(":", 1)
        name = name.strip()
        rest = rest.strip()
        line_match = re.match(r"^Line\s*=\s*\((.*)\)\s*$", rest, re.I)
        if line_match:
            lat.lines[name] = LineDef(name=name, tokens=parse_line_expr(line_match.group(1)), source=line)
            continue

        parts = split_top_level_commas(rest)
        if not parts:
            lat.warnings.append(f"Empty element definition for {name}")
            continue
        kind = parts[0].strip()
        attrs = parse_attrs(parts[1:])
        lat.elements[name] = Element(name=name, kind=kind, attrs=attrs, source=line)

    if lat.use_line is None and "RING" in lat.lines:
        lat.use_line = "RING"
    elif lat.use_line is None and "ring" in lat.lines:
        lat.use_line = "ring"
    return lat


def canon(name: str, preserve_case: bool = False) -> str:
    return name if preserve_case else name.lower()


def attr(attrs: Dict[str, str], *names: str, default: Optional[str] = None) -> Optional[str]:
    for n in names:
        if n.lower() in attrs:
            return attrs[n.lower()]
    return default


def as_float(s: Optional[str]) -> Optional[float]:
    if s is None:
        return None
    try:
        # Bmad files here use simple numeric expressions. Keep eval unavailable.
        return float(s)
    except ValueError:
        return None


def fmt_num(x: float) -> str:
    return f"{x:.12g}"


def fmt_attr_value(v: Optional[str], default: str = "0") -> str:
    if v is None or v == "":
        return default
    return v.strip()


def token_length(lat: Lattice, tok: LineToken, visiting: Optional[set] = None) -> float:
    return tok.mult * name_length(lat, tok.name, visiting=visiting)


def name_length(lat: Lattice, name: str, visiting: Optional[set] = None) -> float:
    if visiting is None:
        visiting = set()
    if name in visiting:
        lat.warnings.append(f"Recursive line-length dependency at {name}; length contribution set to 0")
        return 0.0
    if name in lat.elements:
        return element_length(lat.elements[name])
    if name in lat.lines:
        visiting.add(name)
        total = sum(token_length(lat, tok, visiting=visiting) for tok in lat.lines[name].tokens)
        visiting.remove(name)
        return total
    # Try case-insensitive fallback.
    emap = {k.lower(): k for k in lat.elements}
    lmap = {k.lower(): k for k in lat.lines}
    lk = name.lower()
    if lk in emap:
        return element_length(lat.elements[emap[lk]])
    if lk in lmap:
        return name_length(lat, lmap[lk], visiting=visiting)
    lat.warnings.append(f"Unknown line/element '{name}' while computing length; contribution set to 0")
    return 0.0


def element_length(el: Element) -> float:
    return as_float(attr(el.attrs, "l", default="0")) or 0.0


def p0c_to_energy_gev(p0c: Optional[str]) -> Optional[float]:
    x = as_float(p0c)
    if x is None:
        return None
    # Bmad p0c is commonly eV/c in files of this form.
    if abs(x) > 1e6:
        return x / 1e9
    return x


def tracy_token(tok: LineToken, preserve_case: bool = False) -> str:
    name = canon(tok.name, preserve_case)
    if tok.inverse:
        core = f"inv({name})"
    else:
        core = name
    if tok.mult != 1:
        return f"{tok.mult}*{core}"
    return core


def wrap_sequence(prefix: str, items: List[str], suffix: str = ";", width: int = 96) -> List[str]:
    """Wrap Tracy sequences after a prefix such as 'ring: '."""
    if not items:
        return [prefix.rstrip() + suffix]
    lines: List[str] = []
    current = prefix
    for i, item in enumerate(items):
        piece = item if i == 0 else ", " + item
        if len(current) + len(piece) + len(suffix) > width and current != prefix:
            lines.append(current.rstrip())
            current = "  " + item
        else:
            current += piece
    lines.append(current.rstrip() + suffix)
    return lines


def format_element(
    el: Element,
    lat: Lattice,
    circumference: Optional[float],
    args: argparse.Namespace,
) -> List[str]:
    n = canon(el.name, args.preserve_case)
    kind = el.kind.strip().lower()
    L = fmt_attr_value(attr(el.attrs, "l"), "0")

    if kind in {"drift"}:
        return [f"{n}: Drift, L = {L};"]

    if kind in {"quadrupole", "quad"}:
        b2 = fmt_attr_value(attr(el.attrs, "k1", "b_2"), "0")
        return [f"{n}: Quadrupole, L = {L}, B_2 = {b2}, N = Nquad;"]

    if kind in {"sbend", "rbend", "bend", "bending"}:
        phi = fmt_attr_value(attr(el.attrs, "angle", "phi"), "0")
        b2 = fmt_attr_value(attr(el.attrs, "k1", "b_2"), "0")
        text = f"{n}: Bending, L = {L}, Phi = ({phi})*180.0/pi, B_2 = {b2}, N = Nbend;"
        if len(text) <= 100:
            return [text]
        return [
            f"{n}: Bending, L = {L}, Phi = ({phi})*180.0/pi, B_2 = {b2},",
            "    N = Nbend;",
        ]

    if kind in {"sextupole", "sext"}:
        b3 = fmt_attr_value(attr(el.attrs, "k2", "b_3"), "0")
        return [f"{n}: Sextupole, L = {L}, B_3 = {b3}, N = Nsext;"]

    if kind in {"rfcav", "lcavity", "cavity"}:
        voltage = fmt_attr_value(attr(el.attrs, "voltage"), "0")
        harmonic = attr(el.attrs, "harmon", "har_num", "harnum", "harmonic")
        freq = attr(el.attrs, "frequency", "freq")
        if args.rf_frequency is not None:
            freq = fmt_num(args.rf_frequency)
        elif freq is None and not args.no_auto_rf_frequency:
            h = as_float(harmonic)
            if h is not None and circumference and circumference > 0:
                freq = fmt_num(h * C0 / circumference)
        har = fmt_attr_value(harmonic, "0")
        freq_part = f", Frequency = {freq}" if freq is not None else ""
        return [
            f"{n}: Cavity, L = {L}{freq_part}, Voltage = {voltage},",
            f"     HarNum = {har}, phase = 0, n = 1;",
        ]

    if kind in {"wiggler", "undulator"}:
        lines = [f"{n}: Drift, L = {L};"]
        l_period = as_float(attr(el.attrs, "l_period", "period"))
        bmax = attr(el.attrs, "b_max", "bmax")
        if l_period and l_period > 0:
            periods = int(round((as_float(attr(el.attrs, "l")) or 0.0) / l_period))
            cwig = f"{{ {n}: CWIGGLER, L={L}, B_MAX={fmt_attr_value(bmax, '0')}, PERIODS={periods}, SINUSOIDAL=1, HELICAL=1 }}"
        else:
            cwig = f"{{ {n}: source Bmad {el.kind}; attrs: " + ", ".join(f"{k}={v}" for k, v in el.attrs.items()) + " }}"
        lines.append(cwig)
        return lines

    if kind in {"marker"}:
        return [f"{n}: Marker;"]

    lat.warnings.append(f"Element '{el.name}' has unsupported Bmad type '{el.kind}'; emitted as Marker")
    return [f"{n}: Marker; {{ unsupported Bmad type: {el.kind} }}"]


def format_line(ld: LineDef, args: argparse.Namespace) -> List[str]:
    n = canon(ld.name, args.preserve_case)
    items = [tracy_token(tok, args.preserve_case) for tok in ld.tokens]
    return wrap_sequence(f"{n}: ", items)


def translate(lat: Lattice, args: argparse.Namespace) -> str:
    use_name = args.ring or lat.use_line
    circumference = name_length(lat, use_name) if use_name else None
    energy = args.energy_gev
    if energy is None:
        energy = p0c_to_energy_gev(lat.parameters.get("p0c"))

    out: List[str] = []
    out.append("define lattice;")
    out.append("ringtype = 1;")
    out.append("")
    if energy is not None:
        out.append(f"Energy = {fmt_num(energy)}; {{ Beam momentum [GeV]. }}")
    else:
        out.append("{ Energy was not found in Bmad parameter[p0c]; set Energy manually. }")
    out.append("")
    out.append("dP     = 1e-8;")
    out.append("CODeps = 1e-14;")
    out.append("")
    out.append(f"Nbend = {args.nbend};")
    out.append(f"Nquad = {args.nquad};")
    out.append(f"Nsext = {args.nsext};")
    out.append("")
    out.append("pi = 4.0*arctan(1.0);")
    out.append("c0 = 2.99792458e8;")
    if circumference and circumference > 0:
        out.append(f"{{ Translated from Bmad. Selected circumference = {fmt_num(circumference)} m. }}")
    out.append("")

    for el in lat.elements.values():
        out.extend(format_element(el, lat, circumference, args))
    out.append("")

    for ld in lat.lines.values():
        out.extend(format_line(ld, args))
    out.append("")

    if use_name:
        use_c = canon(use_name, args.preserve_case)
        if use_c != "ring" and "ring" not in {canon(k, args.preserve_case) for k in lat.lines}:
            out.append(f"ring: {use_c};")
        out.append(f"cell: {use_c if use_c == 'ring' else 'ring'}, symmetry = 1;")
    else:
        out.append("{ No Bmad 'use, ...' line found; set the final cell manually. }")
    out.append("")
    out.append("end;")
    out.append("")
    return "\n".join(out)


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Translate a supported subset of Bmad lattice syntax to Tracy-2 .lat syntax."
    )
    p.add_argument("input", type=Path, help="Input .bmad lattice file")
    p.add_argument("-o", "--output", type=Path, help="Output Tracy-2 .lat file; defaults to stdout")
    p.add_argument("--ring", help="Override the Bmad 'use' line/ring name used for circumference and final cell")
    p.add_argument("--energy-gev", type=float, help="Override Energy in GeV")
    p.add_argument("--rf-frequency", type=float, help="Override all RF cavity Frequency values in Hz")
    p.add_argument("--no-auto-rf-frequency", action="store_true", help="Do not compute RF frequency from harmonic*c0/circumference")
    p.add_argument("--nbend", type=int, default=10, help="Tracy integration slices for bends")
    p.add_argument("--nquad", type=int, default=10, help="Tracy integration slices for quadrupoles")
    p.add_argument("--nsext", type=int, default=2, help="Tracy integration slices for sextupoles")
    p.add_argument("--preserve-case", action="store_true", help="Keep Bmad element/line names instead of lowercasing them")
    p.add_argument("--quiet", action="store_true", help="Suppress translation report on stderr")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    lat = parse_bmad(args.input)
    text = translate(lat, args)
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)

    if not args.quiet:
        use_name = args.ring or lat.use_line
        circ = name_length(lat, use_name) if use_name else None
        print("bmad_to_tracy2 report", file=sys.stderr)
        print(f"  input:    {args.input}", file=sys.stderr)
        print(f"  output:   {args.output or '<stdout>'}", file=sys.stderr)
        print(f"  elements: {len(lat.elements)}", file=sys.stderr)
        print(f"  lines:    {len(lat.lines)}", file=sys.stderr)
        print(f"  use/cell: {use_name or '<none>'}", file=sys.stderr)
        if circ is not None:
            print(f"  circumference: {fmt_num(circ)} m", file=sys.stderr)
        if lat.warnings:
            print("  warnings:", file=sys.stderr)
            for w in dict.fromkeys(lat.warnings):
                print(f"    - {w}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
