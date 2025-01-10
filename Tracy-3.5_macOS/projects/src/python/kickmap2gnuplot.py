import math
import numpy
import sys

'''

Numpy:

  Arrays   a[i, j, k]

  Slices   a[i:j:k], i, i+k, i+2k,..., j

  Negative i, j means reverse: n+i, n+j; e.g. n = 10, a[-2:10], 8, 9.

  :        full range

'''

# Global constants.
X_ = 0; Y_ = 1; Z_ = 2


def printf(format, *args): sys.stdout.write(format % args)

def fprintf(outf, format, *args): outf.write(format % args)


def rd_n(inf):
    n = [0, 0]
    for k in range(0, 3): inf.readline()
    L = float(inf.readline().split()[0])
    inf.readline()
    n[X_] = int(inf.readline().split()[0])
    inf.readline()
    n[Y_] = int(inf.readline().split()[0])
    return [n, L]


def rd_theta(inf, n):
    theta = numpy.zeros((n[X_], n[Y_])); dy = []
    inf.readline(); inf.readline()
    dx = (numpy.array(inf.readline().split())).astype(numpy.float)
    for k in range(0, n[Y_]):
        tokens = inf.readline().split()
        dy.append(tokens[0]);
        for j in range(1, len(tokens)):
            theta[j-1, k] = float(tokens[j])
    dy = numpy.asanyarray(dy).astype(numpy.float)
    return [theta, dx, dy]


def rd_kick_map(file_name):
    inf  = open(file_name, 'r');

    [n, L] = rd_n(inf)
    theta = numpy.zeros((2, n[X_], n[Y_]))
    [theta[X_], dx, dy] = rd_theta(inf, n)
    [theta[Y_], dx, dy] = rd_theta(inf, n)

    inf.close();

    printf('\n  n     = [%d, %d]\n', n[X_], n[Y_])
    printf('  L [m] = %7.5f\n', L)

    return [theta, dx, dy, n]


def prt_gnuplot(file_name, theta, dx, dy, n):
    outf = open(file_name, 'w')

    for j in range(0, n[X_]):
        for k in range(0, n[Y_]):
            fprintf(outf, '  %12.5e %12.5e %12.5e %12.5e\n',
                    dx[j], dy[k], theta[X_, j, k], theta[Y_, j, k])
        fprintf(outf, '\n')

    outf.close()


home_dir  = '/home/ria34843/git_repos/tracy-3.5/projects/in/lattice/'
#home_dir  = '/home/johan/git_repos/tracy-3.5/projects/in/lattice/'
file_name = 'lattice_nsls-ii/w100_pole90mm_40div_7m.txt'

[theta, dx, dy, n] = rd_kick_map(home_dir+file_name)

prt_gnuplot('w100_km_gnuplot.out', theta, dx, dy, n)
