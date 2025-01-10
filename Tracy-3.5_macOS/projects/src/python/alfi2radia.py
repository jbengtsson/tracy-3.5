# Script to translate magnet kick map from Alfi's to RADIA format.

import math
import numpy as np
import re
import StringIO
import sys

def printf(format, *args): sys.stdout.write(format % args)

def fprintf(outf, format, *args): outf.write(format % args)

def rd_kick_map(file_name):
    theta = []
    f = open(file_name, 'r')
    # Skip 1st line.
    f.readline()
    with f:
        for l in f:
            if not l.isspace():  # is there anything on the line?
                vals = l.strip().split()
            assert len(vals) == 3
            for k in range(3):
                vals[k] = float(vals[k])
            # [cm] -> [m]. [(Gauss*cm)^2/cm] -> [(Tesla*m)^2/m].
            vals[0] *= 1e-2; vals[1] *= 1e-2; vals[2] *= 1e-12
            theta.append(vals)
            # printf('%6.3f %6.3f %16.8e\n', vals[0], vals[1], vals[2])
    return np.array(theta)

def prt_RADIA(outf, L, n_x, n_y, theta_x, theta_y):
    fprintf(outf, '# Author: A. Shahveh.\n');
    fprintf(outf, '# Main dipole.\n');
    fprintf(outf, '# Magnet Length [m]:\n');
    fprintf(outf, '%8.5f\n', L);
    fprintf(outf, '# Number of Horizontal Points:\n');
    fprintf(outf, '%2d\n', n_x);
    fprintf(outf, '# Number of Vertical Points\n');
    fprintf(outf, '%2d\n', 2*n_y);
    fprintf(outf, '# Horizontal 2nd Order Kick [T2m2]\n');
    fprintf(outf, 'START\n');
    # Print horizontal coordinates.
    for k in range(n_x):
        fprintf(outf, ' %12.5e', theta_x[k, 0])
    fprintf(outf, '\n')
    # Print vertical coordinate and kick vs. x.
    # Above mid-plane.
    for j in range(0, n_x*n_y-1, n_x):
        fprintf(outf, ' %12.5e', theta_x[n_x*n_y-1-j, 1])
        for k in range(n_x):
            fprintf(outf, ' %12.5e', -theta_x[n_x*n_y-1-j-k, 2])
        fprintf(outf, '\n')
    # Below mid-plane; assuming mid-plane symmetry.
    for j in range(0, n_x*n_y-1, n_x):
        fprintf(outf, ' %12.5e', -theta_x[j, 1])
        for k in range(n_x):
            fprintf(outf, ' %12.5e', theta_x[j+k, 2])
        fprintf(outf, '\n')

    fprintf(outf, '# Vertical 2nd Order Kick [T2m2]\n');
    fprintf(outf, 'START\n');
    # Print horizontal coordinates.
    for k in range(n_x):
        fprintf(outf, ' %12.5e', theta_y[k, 0])
    fprintf(outf, '\n')
    # Print vertical coordinate and kick vs. x.
    # Above mid-plane.
    for j in range(0, n_x*n_y-1, n_x):
        fprintf(outf, ' %12.5e', theta_y[n_x*n_y-1-j, 1])
        for k in range(n_x):
            fprintf(outf, ' %12.5e', theta_y[n_x*n_y-1-j-k, 2])
        fprintf(outf, '\n')
    # Below mid-plane; assuming mid-plane symmetry.
    for j in range(n_x, n_x*n_y-1, n_x):
        fprintf(outf, ' %12.5e', -theta_y[j, 1])
        for k in range(n_x):
            fprintf(outf, ' %12.5e', theta_y[j+k, 2])
        fprintf(outf, '\n')


home_dir = '/tmp/ria34843/U/DIAMOND/Model/'

file_name1 = 'dip_int_Thx_KickMap'
file_name2 = 'dip_int_Thy_KickMap'
file_name3 = 'dip_int_kickmap'

theta_x = rd_kick_map(home_dir+file_name1+'.txt')
theta_y = rd_kick_map(home_dir+file_name2+'.txt')
print theta_x[1]
print np.shape(theta_x)

outf = open(file_name3+'_radia.txt', 'w')
prt_RADIA(outf, 0.933, 37, 15, theta_x, theta_y)
