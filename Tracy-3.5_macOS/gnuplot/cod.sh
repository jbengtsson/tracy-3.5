#!/bin/sh

prm1=${1-0}
prm2=${2-"cod.out"}

gnuplot << EOP

ps        = $prm1;
file_name = "$prm2";


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

if (ps) set output "cod_1.".(ext);
set title "Horizontal Closed Orbit";
set xlabel "s [m]";
set ylabel "x [mm]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:9 notitle with impulses ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "cod_2.".(ext);
set title "Vertical Closed Orbit";
set xlabel "s [m]";
set ylabel "y [mm]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:10 notitle with impulses ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "cod_3.".(ext);
set title "Horizontal Corrector Strength";
set xlabel "s [m]";
set ylabel "{\Symbol t}_x [mrad]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:13 notitle with impulses ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "cod_4.".(ext);
set title "Vertical Corrector Strength";
set xlabel "s [m]";
set ylabel "{\Symbol t}_y [mrad]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:14 notitle with impulses ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
