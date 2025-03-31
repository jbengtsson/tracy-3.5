#!/bin/sh

prm1=${1-0}
prm2=${2-"Deta_x"}

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
set style line 3 lt 1 lw 1 lc rgb "cyan";
set style line 4 lt 1 lw 1 lc rgb "red";

if (ps) set output file_name.".(ext);
set title "Linear Momentum Compaction - Driving Terms;
set xlabel "s [m]";
set ylabel "[m]";
set y2range [-1.5:20];
plot "cod.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name.".out" using 3:5 title "{/Symbol h}_x" with lines ls 1, \
     file_name.".out" using 3:6 title "{/Symbol h}'_x" with lines ls 2, \
     file_name.".out" using 3:7 title "d{/Symbol h}_x/d{/Symbol d}" with \
     lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
