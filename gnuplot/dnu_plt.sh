#!/bin/sh

prm1=${1-1}
prm2=${2-0}

gnuplot << EOP

N = $prm1; ps = $prm2; pert = 1;

# MAX-VI: 1, SLS-2: 2.
case = 2;

f_s = 14; l_w = 2;
if (ps == 0) \
  set terminal x11; \
else if (ps == 1) \
  set terminal postscript enhanced color solid lw l_w "Times-Roman" f_s; \
  ext = "ps"; \
else if (ps == 2) \
  set terminal postscript eps enhanced color solid lw l_w "Times-Roman" f_s; \
  ext = "eps"; \
else if (ps == 3) \
  set terminal pdf enhanced color solid linewidth l_w font "Times-Roman f_s"; \
  ext = "pdf"; \
else if (ps == 4) \
  set term pngcairo enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "png";

if ((N == 1) && (case == 1)) \
  N_x = 102; N_y = 68; \
else if ((N == 1) && (case == 2)) \
  N_x = 39; N_y = 15; \
else if (N == 12) \
  N_x = 3; N_y = 1; \
else if (N == 20) \
  N_x = 5; N_y = 3;

if (case == 1) \
  x_min = 102.0; x_max = 102.5; y_min = 68.0; y_max = 68.5; \
else if (case == 2) \
  x_min = 39.0; x_max = 39.5; y_min = 15.0; y_max = 15.5;

# left adjusted labels
set key Left;

set grid;

set style line 1 lt 1 lw l_w lc rgb "blue";
set style line 2 lt 1 lw l_w lc rgb "dark-green";
set style line 3 lt 1 lw l_w lc rgb "red";
set style line 4 lt 1 lw l_w lc rgb "dark-orange";

set clabel "%5.2f"; set key left;

set palette rgbformulae 22, 13, -31 negative;

if (ps) set output "dnu.".(ext);

set multiplot;

set size 0.5, 0.5; set origin 0.0, 0.5;
set title "{/Symbol n}_x vs. A_{x,y}";
set xlabel "A_{x,y} [mm]"; set ylabel "{/Symbol n}_x";
if (!pert) \
  plot "dnu_dAx.out" using 1:(N*(N_x+\$5)) title "A_x" with lines ls 1, \
       "dnu_dAy.out" using 2:(N*(N_x+\$5)) title "A_y" with lines ls 3; \
else \
  plot "dnu_dAx.out" using 1:(N*(N_x+\$5)) title "A_x" with lines ls 1, \
       "dnu_dAy.out" using 2:(N*(N_x+\$5)) title "A_y" with lines ls 3, \
       "ptc/dnu_dAx_pert.out" using 1:(N*(N_x+\$3)) title "A_x (pert)" \
       with lines ls 2, \
       "ptc/dnu_dAy_pert.out" using 2:(N*(N_x+\$3)) title "A_y (pert)" \
       with lines ls 4;

set origin 0.0, 0.0;
set title "{/Symbol n}_y vs. A_{x,y}";
set xlabel "A_{x,y} [mm]"; set ylabel "{/Symbol n}_y"; \
if (!pert) \
  plot "dnu_dAx.out" using 1:(N*(N_y+\$6)) title "A_x" with lines ls 1, \
       "dnu_dAy.out" using 2:(N*(N_y+\$6)) title "A_y" with lines ls 3; \
else \
  plot "dnu_dAx.out" using 1:(N*(N_y+\$6)) title "A_x" with lines ls 1, \
       "dnu_dAy.out" using 2:(N*(N_y+\$6)) title "A_y" with lines ls 3, \
       "ptc/dnu_dAx_pert.out" using 1:(N*(N_y+\$4)) title "A_x (pert)" \
       with lines ls 2, \
       "ptc/dnu_dAy_pert.out" using 2:(N*(N_y+\$4)) title "A_y (pert)" \
       with lines ls 4;

set origin 0.5, 0.5;
set title "Chromaticity";
set xlabel "{/Symbol d} [%]"; set ylabel "{/Symbol n}_x";
set y2label "{/Symbol n}_y";
set ytics nomirror; set y2tics;
plot "chrom2.out" using 1:(N*\$2) title "{/Symbol n}_x" with lines ls 1, \
     "chrom2.out" using 1:(N*\$3) axis x1y2 title "{/Symbol n}_y" \
     with lines ls 3;

set origin 0.5, 0.0;

set style line 1 lw l_w lc rgb "red";
set style line 2 lw l_w lc rgb "dark-orange";
set style line 3 lw l_w lc rgb "blue";
set style line 4 lw l_w lc rgb "dark-green";
set style line 5 lw l_w lc rgb "purple";
set style line 6 lw l_w lc rgb "cyan";

set clabel "%4.1f"; set key left;
#unset clabel;
set noztics; unset colorbox;
set cbrange [0:1];

# x <-> horizontal, y <-> vertical, z <-> perpendicular to screen
set view 0, 0, 1, 1;

#set title "Norm of the Resonance Terms";
#set title "Norm of the Tune Shift Terms";
set title "Tune Footprint";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";

set parametric;

# Compute n for the different resonances:
#
#   n_x*nu_x + n_y*nu_y = n
#
#  nu_y = +/-(n-n_x*nu_y)/n_y
#
# The adjustment is the difference between x_min/x_max and y_min/y_max
# vs. the actual (nu_x,y) values

i10    = floor(x_min) + 1.0;
i01_1  = floor(y_min) + 1.0;
i01_2  = floor(y_min) + 2.0;

i20_1  = floor(x_min) + 1.5;
i20_2  = floor(x_min) + 0.5;
i02_1  = floor(y_min) + 0.5;
i02_2  = floor(y_min) + 1.5;
i02_3  = floor(y_min) + 2.5;

i11    = floor(x_min+y_max);
i1m1_1 = floor(x_min-y_min);
i1m1_2 = floor(x_min-y_min) - 1.0;

i30    = floor(x_min) + 4.0/3.0;
i12_1  = floor(x_min+2*y_max) - 2.0;
i12_2  = floor(x_min+2*y_max) - 1.0;
i12_3  = floor(x_min+2*y_max);
i1m2_1 = floor(x_min-2*y_min);
i1m2_2 = floor(x_min-2*y_min) - 1.0;
i1m2_3 = floor(x_min-2*y_min) - 2.0;
i1m2_4 = floor(x_min-2*y_min) + 1.0;

i40    = floor(x_min) + 1.25;
i04_1  = floor(y_min) + 0.25;
i04_2  = floor(y_min) + 1.25;
i04_3  = floor(y_min) + 1.75;
i04_4  = floor(y_min) + 2.25;
i22_1  = floor(2.0*x_min+2.0*y_max) - 2.0;
i22_2  = floor(2.0*x_min+2.0*y_max);
i2m2_1 = floor(2.0*x_min-2.0*y_min) - 1.0;
i2m2_2 = floor(2.0*x_min-2.0*y_min);
i2m2_3 = floor(2.0*x_min-2.0*y_min) + 1.0;

i14_1  = floor(x_min+4.0*y_max) - 5.0;
i14_2  = floor(x_min+4.0*y_max) - 4.0;
i14_3  = floor(x_min+4.0*y_max) - 3.0;
i14_4  = floor(x_min+4.0*y_max) - 2.0;
i14_5  = floor(x_min+4.0*y_max) - 1.0;
i14_6  = floor(x_min+4.0*y_max);
i32_1  = floor(3.0*x_min+2.0*y_max);
i32_2  = floor(3.0*x_min+2.0*y_max) + 1.0;
i32_3  = floor(3.0*x_min+2.0*y_max) + 2.0;
i1m4_1 = floor(x_min-4.0*y_min);
i1m4_2 = floor(x_min-4.0*y_min) - 1.0;
i1m4_3 = floor(x_min-4.0*y_min) - 2.0;
i1m4_4 = floor(x_min-4.0*y_min) - 3.0;
i1m4_5 = floor(x_min-4.0*y_min) - 4.0;
i1m4_6 = floor(x_min-4.0*y_min) - 5.0;
i3m2_1 = floor(3.0*x_min-2.0*y_min);
i3m2_2 = floor(3.0*x_min-2.0*y_min) + 1.0;
i3m2_3 = floor(3.0*x_min-2.0*y_min) + 2.0;

i21    = floor(2.0*x_min+y_min) + 2.0;
i2m1_1 = floor(2.0*x_min-y_min);
i2m1_2 = floor(2.0*x_min-y_min) + 1.0;

#i60    = floor(x_min);
i06    = floor(y_min) + 14.0/6.0;
#i51    = floor(5.0*x_min+1.0*y_max);
#i5m1   = floor(5.0*x_min-1.0*y_min);
i15    = floor(1.0*x_min+5.0*y_max);
i1m5_1 = floor(1.0*x_min-5.0*y_min) - 6;
i1m5_2 = floor(1.0*x_min-5.0*y_min) - 7;
#i24    = floor(2.0*x_min+4.0*y_max);
i2m4   = floor(2.0*x_min-4.0*y_min) - 5;
i42    = floor(4.0*x_min+2.0*y_max) + 2;
i4m2   = floor(4.0*x_min-2.0*y_min) - 1;
i33    = floor(3.0*x_min+3.0*y_max) + 1;
#i3m3   = floor(3.0*x_min-3.0*y_min) - 3;

#i6m4   = floor(6.0*x_min-4.0*y_min) - 2;

set xrange [x_min:x_max]; set yrange [y_min:y_max]; set zrange [0:1];
set urange [x_min:x_max]; set vrange [y_min:y_max];

splot "fmapdp_est.dat" using (N*(\$3+N_x)):(N*(\$4+N_y)):5 notitle \
      with points pt 2 lc rgb "blue" ps 0.5,\
      "fmap_est.dat" using (N*(\$3+N_x)):(N*(\$4+N_y)):5 notitle \
      with points pt 1 lc rgb "red" ps 0.5, \
      i10,   v,                  1.0 notitle with lines ls 1, \
      u,     i01_1,              1.0 notitle with lines ls 1, \
      u,     i01_2,              1.0 notitle with lines ls 1, \
                                                              \
      i20_1, v,                  1.0 notitle with lines ls 1, \
      i20_2, v,                  1.0 notitle with lines ls 1, \
      u,     i02_1,              1.0 notitle with lines ls 1, \
      u,     i02_2,              1.0 notitle with lines ls 1, \
      u,     i02_3,              1.0 notitle with lines ls 1, \
      u,     i11-u,              1.0 notitle with lines ls 2, \
      u,     u-i1m1_1,           1.0 notitle with lines ls 2, \
      u,     u-i1m1_2,           1.0 notitle with lines ls 2, \
                                                              \
      i30,   v,                  1.0 notitle with lines ls 1, \
      u,     (i12_1-u)/2.0,      1.0 notitle with lines ls 1, \
      u,     (i12_2-u)/2.0,      1.0 notitle with lines ls 1, \
      u,     (i12_3-u)/2.0,      1.0 notitle with lines ls 1, \
      u,     (u-i1m2_1)/2.0,     1.0 notitle with lines ls 1, \
      u,     (u-i1m2_2)/2.0,     1.0 notitle with lines ls 1, \
      u,     (u-i1m2_3)/2.0,     1.0 notitle with lines ls 1, \
      u,     (u-i1m2_4)/2.0,     1.0 notitle with lines ls 1, \
                                                              \
      i40, v,                    1.0 notitle with lines ls 3, \
      u,     i04_1,              1.0 notitle with lines ls 1, \
      u,     i04_2,              1.0 notitle with lines ls 1, \
      u,     i04_3,              1.0 notitle with lines ls 1, \
      u,     i04_4,              1.0 notitle with lines ls 1, \
      u,     (i22_1-2.0*u)/2.0,  1.0 notitle with lines ls 3, \
      u,     (i22_2-2.0*u)/2.0,  1.0 notitle with lines ls 3, \
      u,     (2.0*u-i2m2_1)/2.0, 1.0 notitle with lines ls 3, \
      u,     (2.0*u-i2m2_2)/2.0, 1.0 notitle with lines ls 3, \
      u,     (2.0*u-i2m2_3)/2.0, 1.0 notitle with lines ls 3, \
                                                              \
      u,     (i32_1-3.0*u)/2.0,  1.0 notitle with lines ls 4, \
      u,     (i32_2-3.0*u)/2.0,  1.0 notitle with lines ls 4, \
      u,     (i32_3-3.0*u)/2.0,  1.0 notitle with lines ls 4, \
      u,     (i14_1-u)/4.0,      1.0 notitle with lines ls 4, \
      u,     (i14_2-u)/4.0,      1.0 notitle with lines ls 4, \
      u,     (i14_3-u)/4.0,      1.0 notitle with lines ls 4, \
      u,     (i14_4-u)/4.0,      1.0 notitle with lines ls 4, \
      u,     (i14_5-u)/4.0,      1.0 notitle with lines ls 4, \
      u,     (i14_6-u)/4.0,      1.0 notitle with lines ls 4, \
      u,     (3.0*u-i3m2_1)/2.0, 1.0 notitle with lines ls 4, \
      u,     (3.0*u-i3m2_2)/2.0, 1.0 notitle with lines ls 4, \
      u,     (3.0*u-i3m2_3)/2.0, 1.0 notitle with lines ls 4, \
      u,     (u-i1m4_1)/4.0,     1.0 notitle with lines ls 4, \
      u,     (u-i1m4_2)/4.0,     1.0 notitle with lines ls 4, \
      u,     (u-i1m4_3)/4.0,     1.0 notitle with lines ls 4, \
      u,     (u-i1m4_4)/4.0,     1.0 notitle with lines ls 4, \
      u,     (u-i1m4_5)/4.0,     1.0 notitle with lines ls 4, \
      u,     (u-i1m4_6)/4.0,     1.0 notitle with lines ls 4
                                                              \
#      u,     i21-2.0*u,          1.0 notitle with lines ls 6, \
#      u,     2.0*u-i2m1_1,       1.0 notitle with lines ls 6, \
#      u,     2.0*u-i2m1_2,       1.0 notitle with lines ls 6; \

#                                                              \
#      u,     (4.0*u-i4m2)/2.0,   1.0 notitle with lines ls 5;
#      u,     i06,                1.0 notitle with lines ls 5, \
#      u,     (i15-u)/5.0,        1.0 notitle with lines ls 5, \
#      u,     (u-i1m5_1)/5.0,     1.0 notitle with lines ls 5, \
#      u,     (u-i1m5_2)/5.0,     1.0 notitle with lines ls 5, \
#      u,     (2.0*u-i2m4)/4.0,   1.0 notitle with lines ls 5, \
#      u,     (i42-4.0*u)/2.0,    1.0 notitle with lines ls 5, \
#      u,     (4.0*u-i4m2)/2.0,   1.0 notitle with lines ls 5, \
#      u,     (i33-3.0*u)/3.0,    1.0 notitle with lines ls 5;

#      u,     (6.0*u-i6m4)/4.0,   1.0 notitle with lines ls 5;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
