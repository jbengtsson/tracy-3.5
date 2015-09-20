#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1; eps = 0; phys_app = 0;

ps = $prm1; eps = 0; norm = 0; mom_aper = 1; delta = "2.5"; phys_app = 0;

file1 = "dynap.out";
file2 = "dynap_dp".delta.".out";
file3 = "dynap_dp-".delta.".out";

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

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

if (ps) set output "dynap.ps";

if (!norm) \
  set title "Dynamic Aperture (bare lattice)"; \
  set xlabel "x [mm]"; set ylabel "y [mm]";
if (norm) \
  set title "Normalized Dynamic Aperture"; \
  set xlabel "x^"; set ylabel "y^"; \
  set xrange [-0.01:0.01]; set yrange[0.0:0.01];
#  set xrange [-0.012:0.012]; set yrange[0.0:0.012];

if (!mom_aper) \
  plot file1 using 1:2 title "{/Symbol d}=0" with linespoints ls 1;

if (mom_aper) \
  plot file1 using 1:2 title "{/Symbol d}=0.0%" with linespoints ls 1, \
       file2 using 1:2 title "{/Symbol d}=".delta."%" with linespoints ls 2, \
       file3 using 1:2 title "{/Symbol d}=-".delta."%" with linespoints  ls 3;

if (!ps) pause mouse "click on graph to cont.\n";
