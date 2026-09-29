# Tracking regression tests

`tracking_element_test` is a fast, deterministic test of Tracy's
single-particle tracking through individual elements and short fixed
lattices. It calls the tracking entry points directly; it reads no lattice
file and compares no output files.

## Running

```sh
make check                              # through Automake
./tracking_element_test                 # every named case, PASS or FAIL
TRACKING_TEST_MEASURE=1 ./tracking_element_test   # also print residuals
```

Automake keeps the program output in `tracking_element_test.log`.

## Layout

| File | Cases |
| --- | --- |
| `tracking_elements.cc` | single elements against analytic oracles, 2D frame transforms, fixed lattice, horizontal aperture loss |

`tracking_support.h` holds the assertions, bounds and fixtures the files
share. Each file ends in a `run_*_tests()` function that registers its cases; `tracking_element_test.cc` calls them in turn. Before every case
`run_test` resets the configuration flags to `Read_Lattice`'s defaults for a
ring, the mode lattices are normally tracked in; a case that needs another
mode sets it itself. Each case documents its oracle where it is defined. A new
group of cases goes in a new file.

## Conventions

Expected values are computed from the canonical equations, not captured from
Tracy's output. They are written against Tracy's own conventions:

- `(x, px, y, py, delta, ct)`, with `px`, `py` canonical momenta normalised to
  the reference momentum and `delta` the relative momentum deviation, so a
  drift advances `x` by `L px/(1+delta)`;
- the normal multipole expansion `B_y + i B_x = b_n (x + i y)^(n-1)`, with
  thin-lens kicks `Delta px = -B_y`, `Delta py = B_x`;
- `(delta, ct)` as a canonical pair with `delta` as the coordinate, since the
  drift gives `Delta ct = -dH/ddelta`; the symplecticity checks use this.

A convention error shared by an oracle and the implementation cancels and
passes. Validating the conventions themselves needs an independent code.

## What a passing case establishes

A pass establishes only what its case asserts:

- An analytic case compares tracking with equations evaluated in the test, so
  it checks strengths, signs and momentum dependence, within the conventions
  above.
- A consistency case (symplecticity, composition of transforms, convergence
  with slice count) checks structure. It still passes if both sides share a
  physics error.
- An order check fixes an exponent, not a coefficient. The sector bend's
  second-order coefficient, for example, is bounded but not computed.

## Bounds

Bounds come from an error model: truncation of an approximate oracle
(`C A^n`, with a companion case asserting the order `n`) or integration error
(derived from a slice-convergence case asserting the integrator's order).
Where no such model exists yet, the bound
is a fixed margin, usually ten, over the residual measured when it was set,
and its comment says "Provisional". Every bound's comment records that residual and its date.
Re-measure with `TRACKING_TEST_MEASURE=1` before changing one.

Every assertion rejects a nonfinite expected value, actual value, or bound,
and residual reductions keep a NaN rather than letting `std::max` drop it.

`expect_value`, `expect_state` and `expect_matrix` use a single
`tolerance = 1e-14` scaled by `max(1, |v|)` rather than a bound derived from
the case's amplitudes; the other assertions take an explicit bound.
