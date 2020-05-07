#!/bin/sh

prm1=${1-""}
prm2=${2-0}
prm3=${3-3.0}

gnuplot << EOP

home_dir = "$prm1"; ps = $prm2; delta = sprintf("%3.1f", $prm3);
norm = 0; mom_aper = 1; phys_app = 0;

file1 = (home_dir)."dynap.out";
file2 = (home_dir)."dynap_dp".delta.".out";
file3 = (home_dir)."dynap_dp-".delta.".out";

f_s = 24; l_w = 2;
if (ps == 0) \
  set terminal qt 0 font "Sans, 9"; \
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
  ext = "png"; \
else if (ps == 5) \
  set term svg enhanced font "Times-Roman,f_s"; \
  ext = "svg";

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

if (ps) set output (home_dir)."dynap.".(ext);

if (!norm) \
  set title "Dynamic Aperture (bare lattice)"; \
  set xlabel "x [mm]"; set ylabel "y [mm]"; \
else \
  set title "Normalized Dynamic Aperture"; \
  set xlabel "x^"; set ylabel "y^"; \
  set xrange [-0.01:0.01]; set yrange[0.0:0.01];
#  set xrange [-0.012:0.012]; set yrange[0.0:0.012];

if (!mom_aper) \
  plot file1 using 1:2 title "{/Symbol d}=0" with linespoints ls 1; \
else \
  plot file1 using 1:2 title "{/Symbol d}=0.0%" with linespoints ls 1, \
       file2 using 1:2 title "{/Symbol d}=".delta."%" with linespoints ls 2, \
       file3 using 1:2 title "{/Symbol d}=-".delta."%" with linespoints  ls 3;

if (!ps) pause mouse "click on graph to cont.\n";

EOP
