#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

file_name_1 = "propagate_bc.out";
file_name_2 = "propagate_eps.out";

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

point_size = 0.1;
set style line 1 lt 1 lw 1 lc rgb "blue"  pointsize point_size;
set style line 2 lt 1 lw 1 lc rgb "green" pointsize point_size;
set style line 3 lt 1 lw 1 lc rgb "red"   pointsize point_size;

if (ps) set output "track_1.".(ext)

set key right top;
set title "Barycenter";
set xlabel "turn no";
 set ylabel "<x> [{/Symbol m}m]";
set y2label "<y> [{/Symbol m}m]";
set ytics nomirror;
set y2tics;
plot file_name_1 using 1:(1e6*\$2) title "<x>" with points ls 1, \
     file_name_1 using 1:(1e6*\$4) axis x1y2 title "<y>" with points ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "track_2.".(ext)

set multiplot;

set size 1.0, 0.5;
set origin 0.0, 0.5;
set key right top;
set title "Transverse Phase Space";
set xlabel "turn no";
set ylabel "{/Symbol e}_{x,y} [nm.rad]";
plot file_name_2 using 1:(1e9*\$2) title "{/Symbol e}_x" with points ls 1, \
     file_name_2 using 1:(1e9*\$4) title "{/Symbol e}_y" with points ls 3;

set origin 0.0, 0.0;
set key left top;
set title "Longitudinal Phase Space";
set xlabel "turn no";
set ylabel "{/Symbol s}_s [mm]";
set y2label "{/Symbol d} [1e-3]";
set ytics nomirror;
set y2tics;
plot file_name_2 using 1:(1e3*\$9) title "{/Symbol s}_s" with points ls 1, \
     file_name_2 using 1:(1e3*\$8) axis x1y2 title "{/Symbol d}" \
     with points ls 3;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";

EOP

