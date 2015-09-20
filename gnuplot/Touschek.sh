#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1; eps = 0;

font_size = 24; line_width = 2;
#font_size = 30; line_width = 2;
if (!ps) set terminal x11;
if (ps && eps) set terminal postscript eps enhanced color solid;
#  lw line_width "Times-Roman" font_size;
if (ps && !eps) set terminal postscript enhanced color solid;
#  lw line_width "Times-Roman" font_size;

set grid;

if (ps) set output "Touschek.ps"
set title "Momentum Aperture";
set xlabel "s [m]";
set ylabel "{/Symbol d} [%]";
plot "mom_aper.out" using 2:3 notitle with fsteps ls 2, \
     "mom_aper.out" using 2:4 notitle with fsteps ls 2;
if (!ps) pause mouse "click on graph to cont.\n";

EOP