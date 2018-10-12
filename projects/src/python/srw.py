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


def printf(format, *args): sys.stdout.write(format % args)

def fprintf(outf, format, *args): outf.write(format % args)


def rd_srw(file_name):
    inf  = open(file_name, 'r'); outf = open('srw.out', 'w')

    line = inf.readline()
    fprintf(outf, '%s', line)

    x_min = numpy.zeros(3); dx = numpy.zeros(3); n = numpy.zeros(3)
    for k in range(0, 3):
        x_min[k] = inf.readline().split('#')[1]
        dx[k] = inf.readline().split('#')[1]
        n[k] = inf.readline().split('#')[1]

    B = []
    for line in inf:
        B.append([float(Bk) for Bk in line.split()])
    B = numpy.asanyarray(B)

    inf.close(); outf.close()


def prt_srw_3D(file_name):
    outf = open(file_name, 'w')

    dx = 5e-3; dy = 2.5e-3

    fprintf(outf, '#%12.5e #initial X position [m]\n', x_min[0])
    fprintf(outf, '#%12.5e #step of X [m]\n', dx[0])
    fprintf(outf, '#%1d #number of points vs X\n', n[0])
    fprintf(outf, '#%12.5e #initial Y position [m]\n', x_min[1])
    fprintf(outf, '#%12.5e #step of Y [m]', dx[1])
    fprintf(outf, '#%1d #number of points vs Y\n', n[1])
    fprintf(outf, '#%12.5e #initial Z position [m]\n', x_min[2])
    fprintf(outf, '#%12.5e #step of Z [m]\n', dx[2])
    fprintf(outf, '#%1d #number of points vs Z\n', n[2])

    B = []
    for line in inf:
        B.append([float(Bk) for Bk in line.split()])
    B = numpy.asanyarray(B)

    x = -5e-3
    for i in range(0, 2):
        y = -2.5e-3
        for j in range(0, 2):
            for k in range(0, B[:, 0].size):
                fprintf(outf, '   %13.5e %13.5e %13.5e\n',
                        B[k, 0], B[k, 1], B[k, 2])
            y += dy
        x += dx

    printf('  %7.5f %7.5f\n', x, y)

    inf.close(); outf.close()


home_dir  = '/home/ria34843/git_repos/tracy-3.5/projects/in/lattice/'
file_name = '3pw_1p45dd_29_jb_2.dat'

rd_srw(home_dir+file_name)
