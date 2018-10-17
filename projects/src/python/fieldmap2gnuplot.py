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


def rd_field_map(file_name):
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


def prt_gnuplot(file_name, B, x_min, dx, n):
    outf = open(file_name, 'w')

    x = x_min[X_];
    for j in range(0, n[X_]):
        z = x_min[Z_]
        for k in range(0, n[Z_]):
            fprintf(outf, '  %12.5e %12.5e %12.5e %12.5e %12.5e\n',
                    x, z, B[X_, j, n[Y_]/2, k], B[Y_, j, n[Y_]/2, k],
                    B[Z_, j, n[Y_]/2, k])
            z += dx[Z_]
        x += dx[X_]
        fprintf(outf, '\n')

    outf.close()


home_dir  = '/home/ria34843/git_repos/tracy-3.5/projects/in/lattice/'
file_name = 'lattice_nsls-ii/w100v5_pole90mm_bxyz.csv'

[B, x_min, dx, n] = rd_field_map(home_dir+file_name)

prt_gnuplot('w100_fm_gnuplot_3d.out', B, x_min, dx, n)
