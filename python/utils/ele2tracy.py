import math
import numpy
import re
import StringIO
import sys

# Module to translate from ELEGANT to Tracy-2,3 lattice.


def get_index(tokens, token):
    return tokens.index(token) if token in tokens else None


def marker(line, tokens, decls):
    return  '%s: Marker;' % (tokens[0])

def marker_twiss(line, tokens, decls):
    return marker(line, tokens, decls) + ' { twiss }'

def marker_charge(line, tokens, decls):
    return marker(line, tokens, decls) + ' { %s }' % (line)

def marker_watch(line, tokens, decls):
    return marker(line, tokens, decls) + ' { watch }'

def marker_ematrix(line, tokens, decls):
    return marker(line, tokens, decls) + ' { ematrix }'

def marker_malign(line, tokens, decls):
    return marker(line, tokens, decls) + ' { malign }'

def bpm(line, tokens, decls):
    return  '%s: Beam Position Monitor;' % (tokens[0])

def drift(line, tokens, decls):
    loc_l = tokens.index('l')
    return '%s: Drift, L = %s;' % (tokens[0], get_arg(tokens[loc_l+1], decls))

def rcol(line, tokens, decls):
    return drift(line, tokens, decls) + ' { rcol }'

def bend(line, tokens, decls):
    # CSBEND is modeled by a symplectic integrator.
    # Angles are in [rad] for Elegant and degress for Tracy-2,3.
    loc_l = tokens.index('l')
    loc_phi = tokens.index('angle')
    loc_roll = get_index(tokens, 'tilt')
    loc_e1 = get_index(tokens, 'e1')
    loc_e2 = get_index(tokens, 'e2')
    loc_k = get_index(tokens, 'k1')
    loc_hgap = get_index(tokens, 'hgap')
    loc_fint = get_index(tokens, 'fint')
    loc_n = get_index(tokens, 'n_kicks')
    str = '%s: Bending, L = %s, T = (%s)*180.0/pi' % \
        (tokens[0], get_arg(tokens[loc_l+1], decls),
         get_arg(tokens[loc_phi+1], decls))
    if loc_roll: str += ', Roll = (%s)*180.0/pi' % \
            (get_arg(tokens[loc_roll+1], decls))
    if loc_e1: str += ', T1 = (%s)*180.0/pi' % \
            (get_arg(tokens[loc_e1+1], decls))
    if loc_e2: str += ', T2 = (%s)*180.0/pi' % \
                 (get_arg(tokens[loc_e2+1], decls))
    if loc_k:  str += ', K = %s' %  (get_arg(tokens[loc_k+1], decls))
    if loc_hgap and loc_fint:
        str += ', GAP = 4*%s*%s' % \
                (get_arg(tokens[loc_hgap+1], decls),
                 get_arg(tokens[loc_fint+1], decls))
    else:
        sys.stdout.write('bend: hgap or fint not defined\n')
    if False and loc_n != None:
        str += ', N = %s, Method = 4;' % (tokens[loc_n+1])
    else:
        str += ', N = Nbend, Method = 4;'
    return str

def quad(line, tokens, decls):
    # QUAD is modeled by a third order Taylor expansion.
    loc_l = tokens.index('l')
    loc_k = tokens.index('k1')
    str = '%s: Quadrupole, L = %s, K = %s, N = Nquad, Method = 4;' % \
        (tokens[0], get_arg(tokens[loc_l+1], decls),
         get_arg(tokens[loc_k+1], decls))
    return str

def sext(line, tokens, decls):
    loc_l = tokens.index('l')
    loc_k = tokens.index('k2')
    # The field expansion is a power series for Tracy-2,3 vs. a Taylor expansion
    # for Elegant; i.e., like in MAD-8.
    str = '%s: Sextupole, L = %s, K = %s/2.0, N = Nsext, Method = 4;' % \
        (tokens[0], get_arg(tokens[loc_l+1], decls),
         get_arg(tokens[loc_k+1], decls))
    return str

def oct1(line, tokens, decls):
    loc_l = tokens.index('l')
    loc_k = tokens.index('k3')
    # The field expansion is a power series for Tracy-2,3 vs. a Taylor expansion
    # for Elegant; i.e., like in MAD-8.
    str = '%s: Multipole, L = %s, HOM = (4, %s, 0.0), N = Nsext, Method = 4;' \
        % \
        (tokens[0], get_arg(tokens[loc_l+1], decls),
         get_arg(tokens[loc_k+1], decls))
    return str

def cavity(line, tokens, decls):
    loc_l = tokens.index('l')
    loc_f = tokens.index('freq')
    loc_v = tokens.index('volt')
    loc_phi = tokens.index('phase')
    # loc_entryf = tokens.index('end1_focus')
    # loc_exitf = tokens.index('end2_focus')
    str = '%s: Cavity, L = %s, Frequency = %s, Voltage = %s' % \
          (tokens[0], get_arg(tokens[loc_l+1], decls),
           get_arg(tokens[loc_f+1], decls),
           get_arg(tokens[loc_v+1], decls))
    if loc_phi: str += ', phi = 0*%s' % \
       (get_arg(tokens[loc_phi+1], decls))
    # if loc_entryf: str += ', rf_focus1 = %s' % \
    #    (get_arg(tokens[loc_entryf+1], decls))
    # if loc_exitf: str += ', rf_focus2 = %s' % \
    #    (get_arg(tokens[loc_exitf+1], decls))
    str +=', n = 1;'
    return str

def line(line, tokens, decls):
    n_elem = 10 # No of elements per line.
    n = len(tokens)
    tokens[2] = tokens[2].strip('(')
    tokens[n-1] = tokens[n-1].strip(')')
    str = '%s: ' % (tokens[0])
    if n >= n_elem: str += '\n  '
    for k in range(2, n):
        reverse = tokens[k].startswith('-')
        if reverse:
            tokens[k] = tokens[k].strip('-')
            str += 'inv('
        if (k-1) % (n_elem+1) == 0: str += '\n  '
        str += tokens[k]
        if reverse: str += ')'
        if k < n-1:
            str += ', '
        else:
            str += ';'
    return str


# Elegant -> Tracy-2,3 dictionary.
ele2tracy = {
    'charge'    : marker_charge,
    'mark'      : marker,
    'watch'     : marker_watch,
    'ematrix'   : marker_ematrix,
    'twiss'     : marker_twiss,
    'malign'    : marker_malign,
    'moni'      : bpm,
    'monitor'   : bpm,
    'drif'      : drift,
    'drift'     : drift,
    'edrift'    : drift,
    'csrdrif'   : drift,
    'rcol'      : drift,
    'kicker'    : drift,
    'bumper'    : drift,
    'csben'     : bend,
    'csbend'    : bend,
    'csrcsbend' : bend,
    'quad'      : quad,
    'kquad'     : quad,
    'sext'      : sext,
    'ksext'     : sext,
    'koct'      : oct1,
    'rfca'      : cavity,
    'line'      : line
    }


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass


def parse_rpnc(stack):
    # Reverse Polish Notation calculator.
    arg = stack.pop()
    if arg in '+-*/':
        op = arg
        arg = parse_rpnc(stack)
        if op == '*' and is_number(arg) and float(arg) < 0.0:
            return parse_rpnc(stack) + (' %s (%s)' % (op, arg))
        else:
            return parse_rpnc(stack) + (' %s %s' % (op, arg))
    else:
        return ('%s' % (arg))


def parse_decl(decl, decls):
    stack = (decl.strip('%')).split()
    lhs = stack.pop()
    op = stack.pop()
    decls.append(lhs)
    return ('%s = ' % (lhs)) + parse_rpnc(stack)


def get_arg(str, decls):
    if (str in decls) or (str.find('"') == -1):
        arg = str
    else:
        arg = parse_rpnc((str.strip('"')).split())
    return arg


def parse_definition(line, tokens, decls):
    n_elem = 10; # No of elements per line.

    for k in range(len(tokens)):
        # Remove white space; unless a string.
        if not tokens[k].startswith('"'):
            tokens[k] = re.sub('[\s]', '', tokens[k])
    try:
        str = ele2tracy[tokens[1]](line, tokens, decls)
    except KeyError:
        print '\n*** undefined token!'
        print line
        exit(1)
    return str


def parse_line(line, outf, decls):
    line_lc = line.lower()
    if not line_lc.rstrip():
        # Blank line.
        outf.write('\n')
    elif line_lc.startswith('!'):
        # Comment.
        outf.write('{ %s }\n' % (line.strip('!')))
    elif line_lc.startswith('%'):
        # Declaration.
        outf.write('%s;\n' % (parse_decl(line_lc.strip('%'), decls)))
    else:
        if line_lc.find(':') != -1:
            # Definition.
            tokens = re.split(r'[,:=]', line_lc)
            outf.write('%s\n' % (parse_definition(line_lc, tokens, decls)))
        else:
            print '\n*** undefined statement!'
            print line
            exit(1)

def prt_decl(outf):
    outf.write('define lattice; ringtype = 1;\n')
    outf.write('\nEnergy = 3.5; { Beam momentum [GeV]. }\n')
    outf.write('\ndP = 1e-8; CODeps = 1e-14;\n')
    outf.write('\nMeth = 4; Nbend = 10; Nquad = 10; Nsext = 2;\n')
    outf.write('\npi = 4.0*arctan(1.0); c0 = 2.99792458e8;\n\n')


def transl_file(file_name, decls):
    str = file_name.split('.')[0]+'.lat'
    inf = open(file_name, 'r')
    outf = open(str, 'w')
    prt_decl(outf)
    line = inf.readline()
    while line:
        line = line.strip('\r\n')
        # Remove trailing spaces.
        line = line.rstrip()
        if line.startswith('!'):
            # Comment.
            outf.write('{ %s }\n' % (line.strip('!')))
        else:
            # Skip characters after Line with Comment at end.
            [line, sep, tail] = line.partition('!')
            while line.endswith('&'):
                line = line.strip('&')
                line2 = (inf.readline()).strip('\r\n')
                # Remove trailing spaces.
                line2 = line2.rstrip()
                if line2.startswith('!'):
                    # Comment.
                    while line2.startswith('!'):
                        line2 = (inf.readline()).strip('\r\n')
                        # Remove trailing spaces.
                        line2 = line2.rstrip()
                line += line2
            parse_line(line, outf, decls)
        line = inf.readline()
    outf.write('\nline: ???;\n')
    outf.write('\ncell: line, symmetry = 1;\n')
    outf.write('\nend;\n')


home_dir = ''

decls = []

transl_file(home_dir+sys.argv[1], decls)
