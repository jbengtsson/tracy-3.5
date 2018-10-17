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
    # Skip first three lines.
    inf.readline(); inf.readline(); inf.readline()
    L = float(inf.readline().split()[0])
    inf.readline()
    n[X_] = int(inf.readline().split()[0])
    inf.readline()
    n[Y_] = int(inf.readline().split()[0])
    return [n, L]

def rd_theta(inf, n_y):
    theta = []; dy = []
    inf.readline(); inf.readline()
    dx = (numpy.array(inf.readline().split())).astype(numpy.float)
    for k in range(0, n_y):
        tokens = inf.readline().split()
        dy.append(tokens[0]);
        for tx in tokens[1:]:
            theta.append(tx)
    dy = numpy.asanyarray(dy).astype(numpy.float)
    theta = numpy.asanyarray(theta).astype(numpy.float)
    return [theta, dx, dy]


def rd_kick_map(file_name):
    inf  = open(file_name, 'r');

    [n, L] = rd_n(inf)
    print L, n[X_], n[Y_]

    [theta_x, dx, dy] = rd_theta(inf, n[Y_])
    [theta_y, dx, dy] = rd_theta(inf, n[Y_])

    inf.close();

    print theta_x
    print theta_y

    printf('\n  n          = [%d, %d]\n', n[X_], n[Y_])
    print '  Shape{ theta_x } =', theta_x.shape
    print '  Shape{ theta_y } =', theta_y.shape

    # return [theta_x, x_min, dx, n]


def prt_gnuplot_3D_data(file_name, B, x_min, dx, n):
    outf = open(file_name, 'w')

    fprintf(outf,
            '#Bx [T], By [T], Bz [T] on 3D mesh: inmost loop vs X'
            ' (horizontal transverse position), outmost loop vs Z'
            ' (longitudinal position)\n')
    fprintf(outf, '#%12.5e #initial X position [m]\n', x_min[X_])
    fprintf(outf, '#%12.5e #step of X [m]\n', dx[X_])
    fprintf(outf, '#%1d #number of points vs X\n', n[X_])
    fprintf(outf, '#%12.5e #initial Y position [m]\n', x_min[Y_])
    fprintf(outf, '#%12.5e #step of Y [m]\n', dx[Y_])
    fprintf(outf, '#%1d #number of points vs Y\n', n[Y_])
    fprintf(outf, '#%12.5e #initial Z position [m]\n', x_min[Z_])
    fprintf(outf, '#%12.5e #step of Z [m]\n', dx[Z_])
    fprintf(outf, '#%1d #number of points vs Z\n', n[Z_])

    x = numpy.zeros(3)
    k = 0
    x[Z_] = x_min[Z_] - dx[Z_]
    while k < B[:, 0].size:
        x[Z_] += dx[Z_]
        x[Y_] = x_min[Y_] - dx[Y_]
        for i in range(0, n[Y_]):
            x[Y_] += dx[Y_]
            x[X_] = x_min[X_] - dx[X_]
            for j in range(0, 3):
                x[X_] += dx[X_]
                fprintf(outf, '   %13.5e %13.5e %13.5e\n',
                        B[k, 0], B[k, 1], B[k, 2])
            k += 1

    printf('  x_max      = [%12.5e, %12.5e, %12.5e]\n', x[X_], x[Y_], x[Z_])

    outf.close()


home_dir  = '/home/johan/git_repos/tracy-3.5/projects/in/lattice/'
file_name = 'w100_pole90mm_40div_7m.txt'

rd_kick_map(home_dir+file_name)

#prt_gnuplot_3D_data('w100_gnuplot_3D.out', B, x_min, dx, n)
