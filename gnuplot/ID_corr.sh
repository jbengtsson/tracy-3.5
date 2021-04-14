#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

file_name_1 = "linlat1.out";
file_name_2 = "ID_corr_1.out";
file_name_3 = "ID_corr_res_1.out";

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
  ext = "png";

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "ID_corr_1.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol Db}_x/{/Symbol b}_x [%]";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot file_name_1 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_3 using 1:2 notitle with impulses ls 1;

set origin 0.0, 0.0;
set title "{/Symbol Db}_y/{/Symbol b}_y [%]";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot file_name_1 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_3 using 1:3 notitle with impulses ls 3;

unset multiplot;

if (!ps) pause -1;

if (ps) set output "ID_corr_2.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol Dn}_x";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot file_name_1 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_3 using 1:4 notitle with impulses ls 1;

set origin 0.0, 0.0;
set title "{/Symbol Dn}_y";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot file_name_1 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_3 using 1:5 notitle with impulses ls 3;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "ID_corr_3.ps"

set size 1.0, 1.0; set origin 0.0, 0.0;
set title "b_{2}L";
set xlabel "s [m]";
set ylabel "";
plot file_name_1 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_2 using 2:4 notitle with impulses ls 1;

unset multiplot;

if (!ps) pause mouse "click on graph to cont.\n";

EOP
