#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

f_s = 14; l_w = 2;
if (ps == 0) \
  set terminal x11; \
else if (ps == 1) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" f_s; \
  ext = "ps"; \
else if (ps == 2) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" f_s; \
  ext = "eps"; \
else if (ps == 3) \
  set terminal pdf enhanced color solid linewidth 2 font "Times-Roman f_s"; \
  ext = "pdf";

set grid;

if (ps) set output "Touschek.".ext
set title "Momentum Aperture";
set xlabel "s [m]";
set ylabel "{/Symbol d} [%]";
plot "mom_aper.out" using 2:3 notitle with fsteps ls 2, \
     "mom_aper.out" using 2:4 notitle with fsteps ls 2;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
