#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

f_s = 24; l_w = 2;
# Enhanced is needed for Greek characters.
if (ps == 0) \
  set terminal qt 0 enhanced font "Sans, 9"; \
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

if (ps) set output "xi.".(ext)

set title "Chromaticity";
set xlabel "{/Symbol d} [%]"; set ylabel "{/Symbol n}_x";
set y2label "{/Symbol n}_y";
set ytics nomirror; set y2tics;
  plot "chrom2.out" using 1:2 title "{/Symbol n}_x" with lines ls 1, \
       "chrom2.out" using 1:3 axis x1y2 title "{/Symbol n}_y" \
       with lines ls 3;

if (!ps) pause mouse "click on graph to cont.\n";
