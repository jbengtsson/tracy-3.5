#!/bin/sh

prm1=${1-1}
prm2=${2-0}
prm3=${3-1}

gnuplot << EOP

ps = $prm2; scale = $prm3 == 1;

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
  set terminal pdf enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "pdf"; \
else if (ps == 4) \
  set term pngcairo enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "png";

if ($prm1 == 1) \
  N = 1; N_x = 101; N_y = 27; \
else if ($prm1 == 20) \
  N = 20; N_x = 5; N_y = 1;

#nu_x_min = 32.99; nu_x_max = 33.51; nu_y_min = 15.99; nu_y_max = 16.51;
#nu_x_min = 41.99; nu_x_max = 43.0; nu_y_min = 13.99; nu_y_max = 15.01;
#nu_x_min = 101.49; nu_x_max = 102.01; nu_y_min = 27.49; nu_y_max = 28.01;
nu_x_min = 100.9; nu_x_max = 101.6; nu_y_min = 26.9; nu_y_max = 27.6;

# x_min = -20.0; x_max = 20.0; y_min = -10.0; y_max = 10.0; \
# delta_min = -6.0; delta_max = 6.0;
x_min = -2.0; x_max = 2.0; y_min = -2.0; y_max = 2.0; \
delta_min = -5.1; delta_max = 5.1;

set grid;

set style line 1 lw 1 lc rgb "red";
set style line 2 lw 1 lc rgb "dark-orange";
set style line 3 lw 1 lc rgb "blue";
set style line 4 lw 1 lc rgb "dark-green";
set style line 5 lw 1 lc rgb "purple";
set style line 6 lw 1 lc rgb "cyan";

file_name = "`echo $TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name using (\$1/255):(\$2/255):(\$3/255);

#set cbrange [-10:-2];
set noztics; unset clabel;
set view map;
# To set y-axis to left side and avoid compression of color box.
unset pm3d;

set parametric;

# Compute n for the different resonances:
#
#   n_x*nu_x + n_y*nu_y = n
#
#  nu_y = +/-(n-n_x*nu_y)/n_y
#
# The adjustment is the difference between x_min/x_max and y_min/y_max
# vs. the actual (nu_x,y) values

i10    = floor(nu_x_min) + 1.0;
i01_1  = floor(nu_y_min) + 1.0;
i01_2  = floor(nu_y_min) + 2.0;

i20    = floor(nu_x_min) + 1.5;
i02_1  = floor(nu_y_min) + 0.5;
i02_2  = floor(nu_y_min) + 1.5;
i02_3  = floor(nu_y_min) + 2.5;

i11    = floor(nu_x_min+nu_y_max);
i1m1_1 = floor(nu_x_min-nu_y_min);
i1m1_2 = floor(nu_x_min-nu_y_min) - 1.0;

i30    = floor(nu_x_min) + 4.0/3.0;
i12_1  = floor(nu_x_min+2*nu_y_max) - 2.0;
i12_2  = floor(nu_x_min+2*nu_y_max) - 1.0;
i12_3  = floor(nu_x_min+2*nu_y_max);
i1m2_1 = floor(nu_x_min-2*nu_y_min);
i1m2_2 = floor(nu_x_min-2*nu_y_min) - 1.0;
i1m2_3 = floor(nu_x_min-2*nu_y_min) - 2.0;
i1m2_4 = floor(nu_x_min-2*nu_y_min) + 1.0;

i40    = floor(nu_x_min) + 1.25;
i04_1  = floor(nu_y_min) + 0.25;
i04_2  = floor(nu_y_min) + 1.25;
i04_3  = floor(nu_y_min) + 1.75;
i04_4  = floor(nu_y_min) + 2.25;
i22_1  = floor(2.0*nu_x_min+2.0*nu_y_max) - 2.0;
i22_2  = floor(2.0*nu_x_min+2.0*nu_y_max);
i2m2_1 = floor(2.0*nu_x_min-2.0*nu_y_min) - 1.0;
i2m2_2 = floor(2.0*nu_x_min-2.0*nu_y_min);
i2m2_3 = floor(2.0*nu_x_min-2.0*nu_y_min) + 1.0;

i14_1  = floor(nu_x_min+4.0*nu_y_max) - 5.0;
i14_2  = floor(nu_x_min+4.0*nu_y_max) - 4.0;
i14_3  = floor(nu_x_min+4.0*nu_y_max) - 3.0;
i14_4  = floor(nu_x_min+4.0*nu_y_max) - 2.0;
i14_5  = floor(nu_x_min+4.0*nu_y_max) - 1.0;
i14_6  = floor(nu_x_min+4.0*nu_y_max);
i32_1  = floor(3.0*nu_x_min+2.0*nu_y_max);
i32_2  = floor(3.0*nu_x_min+2.0*nu_y_max) + 1.0;
i32_3  = floor(3.0*nu_x_min+2.0*nu_y_max) + 2.0;
i1m4_1 = floor(nu_x_min-4.0*nu_y_min);
i1m4_2 = floor(nu_x_min-4.0*nu_y_min) - 1.0;
i1m4_3 = floor(nu_x_min-4.0*nu_y_min) - 2.0;
i1m4_4 = floor(nu_x_min-4.0*nu_y_min) - 3.0;
i1m4_5 = floor(nu_x_min-4.0*nu_y_min) - 4.0;
i1m4_6 = floor(nu_x_min-4.0*nu_y_min) - 5.0;
i3m2_1 = floor(3.0*nu_x_min-2.0*nu_y_min);
i3m2_2 = floor(3.0*nu_x_min-2.0*nu_y_min) + 1.0;
i3m2_3 = floor(3.0*nu_x_min-2.0*nu_y_min) + 2.0;

i21    = floor(2.0*nu_x_min+nu_y_min) + 2.0;
i2m1_1 = floor(2.0*nu_x_min-nu_y_min);
i2m1_2 = floor(2.0*nu_x_min-nu_y_min) + 1.0;

#i60    = floor(nu_x_min);
i06    = floor(nu_y_min) + 14.0/6.0;
#i51    = floor(5.0*nu_x_min+1.0*nu_y_max);
#i5m1   = floor(5.0*nu_x_min-1.0*nu_y_min);
i15    = floor(1.0*nu_x_min+5.0*nu_y_max);
i1m5_1 = floor(1.0*nu_x_min-5.0*nu_y_min) - 6;
i1m5_2 = floor(1.0*nu_x_min-5.0*nu_y_min) - 7;
#i24    = floor(2.0*nu_x_min+4.0*nu_y_max);
i2m4   = floor(2.0*nu_x_min-4.0*nu_y_min) - 5;
i42    = floor(4.0*nu_x_min+2.0*nu_y_max) + 2;
i4m2   = floor(4.0*nu_x_min-2.0*nu_y_min) - 1;
i33    = floor(3.0*nu_x_min+3.0*nu_y_max) + 1;
#i3m3   = floor(3.0*nu_x_min-3.0*nu_y_min) - 3;

#i6m4   = floor(6.0*nu_x_min-4.0*nu_y_min) - 2;

set urange [nu_x_min:nu_x_max]; set vrange [nu_y_min:nu_y_max];

if (ps) set output "fmap_1.".ext;

set multiplot;

#set pm3d at b map;
#set contour;
unset colorbox;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Tune Shift";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";
if (scale) set xrange [nu_x_min:nu_x_max]; set yrange [nu_y_min:nu_y_max]; 
splot "fmap.out" using \
      ((abs(\$3-int(\$3)) > 1e-6)? N*(N_x+\$3) : NaN): \
      ((abs(\$4-int(\$4)) > 1e-6)? N*(N_y+\$4) : NaN):7 \
      notitle w points pt 13 lt palette z, \
      i10,   v,                  1.0 notitle with lines ls 1, \
      u,     i01_1,              1.0 notitle with lines ls 1, \
      u,     i01_2,              1.0 notitle with lines ls 1, \
                                                              \
      i20,   v,                  1.0 notitle with lines ls 1, \
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
      u,     (u-i1m4_6)/4.0,     1.0 notitle with lines ls 4, \
                                                              \
      u,     i21-2.0*u,          1.0 notitle with lines ls 6, \
      u,     2.0*u-i2m1_1,       1.0 notitle with lines ls 6, \
      u,     2.0*u-i2m1_2,       1.0 notitle with lines ls 6;

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

set pm3d at b map;
#set contour;
unset colorbox;

set origin 0.0, 0.0;
set title "Diffusion Map";
set xlabel "A_x"; set ylabel "A_y";
if (scale) set xrange [x_min:x_max]; set yrange [y_min:y_max];
splot "fmap.out" using 1:2:((\$7 != -2.0)? \$7 : NaN) notitle lt palette z;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";

unset output;

# Workaround bug for multiplot.
reset;

set grid;

set style line 1 lw 1 lc rgb "red";
set style line 2 lw 1 lc rgb "dark-orange";
set style line 3 lw 1 lc rgb "blue";
set style line 4 lw 1 lc rgb "dark-green";
set style line 5 lw 1 lc rgb "purple";
set style line 6 lw 1 lc rgb "cyan";

file_name = "`echo $TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name using (\$1/255):(\$2/255):(\$3/255);

#set cbrange [-10:-2];
set noztics; unset clabel;
set view map;
# To set y-axis to left side and avoid compression of color box.
unset pm3d;

set parametric;

if (ps) set output "fmap_2.".ext;

set multiplot;

unset contour;
unset pm3d;
set colorbox;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Tune Shift";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";
if (scale) set xrange [nu_x_min:nu_x_max]; set yrange [nu_y_min:nu_y_max];
splot "fmapdp.out" using \
      ((abs(\$3-int(\$3)) > 1e-6)? N*(N_x+\$3) : NaN): \
      ((abs(\$4-int(\$4)) > 1e-6)? N*(N_y+\$4) : NaN):7 \
      notitle w points pt 13 lt palette z, \
      i10,   v,                  1.0 notitle with lines ls 1, \
      u,     i01_1,              1.0 notitle with lines ls 1, \
      u,     i01_2,              1.0 notitle with lines ls 1, \
                                                              \
      i20,   v,                  1.0 notitle with lines ls 1, \
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
      u,     (u-i1m4_6)/4.0,     1.0 notitle with lines ls 4, \
                                                              \
      u,     i21-2.0*u,          1.0 notitle with lines ls 6, \
      u,     2.0*u-i2m1_1,       1.0 notitle with lines ls 6, \
      u,     2.0*u-i2m1_2,       1.0 notitle with lines ls 6;

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

set pm3d at b map;
#set contour;
unset colorbox;

set origin 0.0, 0.0;
set title "Diffusion Map";
set xlabel "{/Symbol d} [%]"; set ylabel "A_x";
if (scale) set xrange [delta_min:delta_max]; set yrange [x_min:x_max];
splot "fmapdp.out" using 1:2:((\$7 != -2.0)? \$7 : NaN) notitle lt palette z;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
