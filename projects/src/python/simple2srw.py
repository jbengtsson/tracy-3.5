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


def rd_simple(file_name):
    x_min = numpy.zeros(3); dx = numpy.zeros(3); n = [0, 0, 0]

    inf  = open(file_name, 'r');

    # Skip first two lines.
    inf.readline()
    inf.readline()

    B = [];  y0 = -1e30; z0 = -1e30
    for line in inf:
        tokens = line.split()
        [y, z] = [1e-3*float(tokens[1]), 1e-3*float(tokens[0])]
        B.append([float(tokens[3]), float(tokens[2]), float(tokens[4])])
        if y > y0:
            n[Y_] += 1
            if n[Y_] == 2:
                dx[Y_] = y - y0
            y0 = y
        if z > z0:
            n[Z_] += 1
            if n[Z_] == 1:
                x_min[Y_] = y;x_min[Z_] = z; 
            if n[Z_] == 2:
                dx[Z_] = z - z0
            z0 = z
        else:
            break

    for line in inf:
        tokens = line.split()
        y = 1e-3*float(tokens[1])
        B.append([float(tokens[3]), float(tokens[2]), float(tokens[4])])
        if y > y0:
            n[Y_] += 1
            y0 = y

    B = numpy.asanyarray(B)

    inf.close();

    printf('\n  n        = [%d, %d, %d]\n', n[X_], n[Y_], n[Z_])
    printf('  x_min    = [%12.5e, %12.5e, %12.5e]\n',
           x_min[X_], x_min[Y_], x_min[Z_])
    printf('  dx       = [%12.5e, %12.5e, %12.5e]\n', dx[X_], dx[Y_], dx[Z_])
    print '  Shape{B} =', B.shape

    return [B, x_min, dx, n]


def rd_field_map_csv(file_name):
    inf  = open(file_name, 'r');

    n = [int(k) for k in inf.readline().strip('\n\r').split(',')]
    inf.readline()

    B = numpy.zeros((3, n[X_], n[Y_], n[Z_]))
    for i in range(0, n[X_]):
        for k in range(0, n[Z_]):
            for j in range(0, n[Y_]):
                tokens = inf.readline().strip('\n\r').split(',')
                x = 1e-3*numpy.asanyarray(tokens[0:3]).astype(numpy.float)
                B[:, i, j, k] = \
                    numpy.asanyarray(tokens[3:6]).astype(numpy.float)
                if ((i == 0) and (j == 0) and (k == 0)):
                    x_min = x
                elif ((i == 1) and (j == 1) and (k == 1)):
                    dx = x - x_min

    inf.close();

    printf('\n  n        = [%d, %d, %d]\n', n[X_], n[Y_], n[Z_])
    printf('  x_min    = [%12.5e, %12.5e, %12.5e]\n',
           x_min[X_], x_min[Y_], x_min[Z_])
    printf('  dx       = [%12.5e, %12.5e, %12.5e]\n', dx[X_], dx[Y_], dx[Z_])
    print '  Shape{B} =', B.shape
    return [B, x_min, dx, n]


def prt_srw(file_name, B, x_min, dx, n):
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
    x[Z_] = x_min[Z_] - dx[Z_];
    for i in range(0, n[Z_]):
        x[Z_] += dx[Z_]
        x[Y_] = x_min[Y_] - dx[Y_]
        for j in range(0, n[Y_]):
            x[Y_] += dx[Y_]
            x[X_] = x_min[X_] - dx[X_]
            for k in range(0, n[X_]):
                x[X_] += dx[X_]
                fprintf(outf, '  %12.5e %12.5e %12.5e\n',
                        B[X_, k, j, i], B[Y_, k, j, i], B[Z_, k, j, i])

    outf.close()

    printf('  x_max    = [%12.5e, %12.5e, %12.5e]\n', x[X_], x[Y_], x[Z_])


home_dir  = '/home/ria34843/git_repos/tracy-3.5/projects/in/lattice/'
#home_dir  = '/home/johan/git_repos/tracy-3.5/projects/in/lattice/'

# file_name = '3pw_1p45dd_29_jb_2.dat'
# [B, x_min, dx, n] = rd_simple(home_dir+file_name)
# dx[X_] = 5e-3; x_min[X_] = -dx[X_]; n[X_] = 3;
# printf('\n  x_min      = [%12.5e, %12.5e, %12.5e]\n',
#        x_min[X_], x_min[Y_], x_min[Z_])
# prt_srw('w100.out', B, x_min, dx, n)

file_name = 'lattice_nsls-ii/w100v5_pole90mm_bxyz.csv'
[B, x_min, dx, n] = rd_field_map_csv(home_dir+file_name)
prt_srw('w100_srw.out', B, x_min, dx, n)
