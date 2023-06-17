#!/bin/sh

prm1=${1-""}
prm2=${2-0}

gnuplot << EOP

home_dir = "$prm1"; ps = $prm2;

c0      = 2.99792458e8
C       = 240.0
alpha_c = 7.038e-4

N   = 1
N_x = 17
N_y = 6

f_RF  = 499.6366302e6
f_rev = c0/C
h_RF  = f_RF/f_rev
print sprintf("\n  f_rev [MHz]     = %5.3f", 1e-6*f_rev)
print sprintf("  r_RF  [Mhz]     = %4.1f", 1e-6*f_RF)
print sprintf("  Harmonic number = %4.1f", h_RF)

f_beta(nu, N)  = (nu-N)*f_rev
df_RF(delta)   = -alpha_c*f_RF*delta

file1  = (home_dir)."dnu_dAx.out";
file12 = (home_dir)."dnu_dAx_pert.out"
file2  = (home_dir)."dnu_dAy.out";
file22 = (home_dir)."dnu_dAy_pert.out";
file3  = (home_dir)."chrom2.out";
file32 = (home_dir)."chrom2_pert.out";

f_s = 36; l_w = 2;
if (ps == 0) \
  set terminal qt 0 font "Sans, 9"; \
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
  ext = "png"; \
else if (ps == 5) \
  set term svg enhanced font "Times-Roman,f_s"; \
  ext = "svg";

# left adjusted labels
set key Left;

set grid;

set style line 1 lt 1 lw l_w lc rgb "blue";
set style line 2 lt 1 lw l_w lc rgb "dark-green";
set style line 3 lt 1 lw l_w lc rgb "red";
set style line 4 lt 1 lw l_w lc rgb "dark-orange";

set clabel "%5.2f"; set key left;

set palette rgbformulae 22, 13, -31 negative;

if (ps) set output (home_dir)."dnu_1.".(ext);

set title "{/Symbol n}_x vs. A_{x,y}";
set xlabel "A_{x,y} [mm]"; set ylabel "{/Symbol n}_x";
plot file1 using 1:(N*(N_x+\$5)) title "A_x" with lines ls 1, \
     file2 using 2:(N*(N_x+\$5)) title "A_y" with lines ls 3; \
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output (home_dir)."dnu_2.".(ext);

set title "{/Symbol n}_y vs. A_{x,y}";
set xlabel "A_{x,y} [mm]"; set ylabel "{/Symbol n}_y"; \
plot file1 using 1:(N*(N_y+\$6)) title "A_x" with lines ls 1, \
     file2 using 2:(N*(N_y+\$6)) title "A_y" with lines ls 3; \
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output (home_dir)."dnu_3.".(ext);

set title "Chromaticity";
set xlabel "{/Symbol D}f_{RF} [kHz]";
set ylabel  "{/Symbol D}f_{{/Symbol b},x} [kHz]";
set y2label "{/Symbol D}f_{{/Symbol b},y} [kHz]";
set ytics nomirror; set y2tics;
plot file3 using (1e-3*df_RF(1e-2*\$1)):(1e-3*f_beta(\$2, N_x)) \
     title "{/Symbol n}_x" with lines ls 1, \
     file3 using (1e-3*df_RF(1e-2*\$1)):(1e-3*f_beta(\$3, N_y)) \
     axis x1y2 title "{/Symbol n}_y" with lines ls 3; \
if (!ps) pause mouse "click on graph to cont.\n";
