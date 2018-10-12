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


def rd_srw(file_name):
    inf  = open(file_name, 'r');

    x_min = numpy.zeros(3); dx = numpy.zeros(3); n = numpy.zeros(3)
    # Skipe 1st line.
    line = inf.readline()
    for k in range(0, 3):
        x_min[k] = inf.readline().split('#')[1]
        dx[k] = inf.readline().split('#')[1]
        n[k] = inf.readline().split('#')[1]

    B = []
    for line in inf:
        B.append([float(Bk) for Bk in line.split()])
    B = numpy.asanyarray(B)

    inf.close();

    return [B, x_min, dx, n]


def prt_srw_3D(file_name, B, x_min, dx, n):
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
    x[X_] = x_min[X_]
    for i in range(0, 2):
        x[Y_] = x_min[Y_]
        for j in range(0, 2):
            for k in range(0, B[:, 0].size):
                fprintf(outf, '   %13.5e %13.5e %13.5e\n',
                        B[k, 0], B[k, 1], B[k, 2])
            x[Y_] += dx[Y_]
        x[X_] += dx[X_]

    printf('  %7.5f %7.5f\n', x[X_], x[Y_])

    outf.close()


home_dir  = '/home/ria34843/git_repos/tracy-3.5/projects/in/lattice/'
file_name = '3pw_1p45dd_29_jb_2.dat'

[B, x_min, dx, n] = rd_srw(home_dir+file_name)

dx[X_]    = 5e-3;    dx[Y_]    = 2.5e-3
x_min[X_] = -dx[X_]; x_min[Y_] = -dx[Y_]
n[X_]     = 3;       n[Y_]     = 3;
prt_srw_3D('srw.out', B, x_min, dx, n)
