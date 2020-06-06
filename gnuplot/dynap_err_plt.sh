#!/bin/sh

prm1=${1-""}
prm2=${2-0}

gnuplot << EOP

home_dir = "$prm1"; ps = $prm2; phys_app = 0;

file1 = (home_dir)."DA_bare_0.00.out";
file2 = (home_dir)."DA_real_0.00.out";
file3 = (home_dir)."DA_bare.out";
file4 = (home_dir)."DA_real.out";

f_s = 24; l_w = 2;
if (ps == 0) \
  set terminal qt; \
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

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue" ps 2 pt 1;
set style line 2 lt 1 lw 1 lc rgb "green" ps 2 pt 1;
set style line 3 lt 1 lw 1 lc rgb "red" ps 2 pt 1;

# draw projection of mechanical aperture
Ax = 17.5; Ay = 12.5;
beta_max_y = 25.5; beta_inj_y =  3.1;

if (phys_app) \
  x_hat = Ax; y_hat = Ay*sqrt(beta_inj_y/beta_max_y); \
  set arrow from -x_hat, 0.0 to -x_hat, y_hat nohead \
  lt 1 lw 1 lc rgb "black"; \
  set arrow from -x_hat, y_hat to  x_hat, y_hat nohead \
  lt 1 lw 1 lc rgb "black"; \
  set arrow from  x_hat, y_hat to  x_hat,   0.0 nohead \
  lt 1 lw 1 lc rgb "black";

if (ps) set output "dynap_err_1.".(ext);
set title "Dynamic Aperture\n";
set xlabel "x [mm]"; set ylabel "y [mm]";
#set xrange [-20:20];
#set yrange [0:4];
plot file1 using  1:2 title "bare" with linespoints ls 1, \
     file2 using  1:2 notitle with points ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

unset arrow;

if (ps) set output "dynap_err_2.".(ext);

#set multiplot;

#set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Horizontal Momentum Aperture\n";
set xlabel "{/Symbol d} [%]"; set ylabel "x^ [mm]";
set yrange [0:];
plot file3 using 1:5 title "bare" with linespoints ls 2, \
     file4 using 1:11:13 title "w errors" with errorbars ls 1, \
     file4 using 1:11 notitle with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "dynap_err_3.".(ext);

#set origin 0.0, 0.0;
set title "Vertical Momentum Aperture\n";
set xlabel "{/Symbol d} [%]"; set ylabel "y^ [mm]";
set yrange [0:];
plot file3 using 1:6 title "bare" with linespoints ls 2, \
     file4 using 1:14:16 title "w errors" with errorbars ls 3, \
     file4 using 1:14 notitle with lines ls 3;

#unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";

unset output;

# Workaround bug for multiplot.
reset;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue" ps 2 pt 1;
set style line 2 lt 1 lw 1 lc rgb "green" ps 2 pt 1;
set style line 3 lt 1 lw 1 lc rgb "red" ps 2 pt 1;

if (ps) set output "dynap_err_4.".(ext);

#set multiplot;

#set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Horizontal Momentum Acceptance\n";
set xlabel "{/Symbol d} [%]"; set ylabel "A_x [mm{/Symbol \327}mrad]";
set yrange [0:];
plot file3 using 1:3 title "bare" with linespoints ls 2, \
     file4 using 1:5:7 title "w errors" with errorbars ls 1, \
     file4 using 1:5 notitle with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "dynap_err_5.".(ext);

#set origin 0.0, 0.0;
set title "Vertical Momentum Acceptance\n";
set xlabel "{/Symbol d} [%]"; set ylabel "A_y [mm{/Symbol \327}mrad]";
set yrange [0:];
plot file3 using 1:4 title "bare" with linespoints ls 2, \
     file4 using 1:8:10 title "w errors" with errorbars ls 3, \
     file4 using 1:8 notitle with lines ls 3;

#unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
