#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

file_name_1 = "track.out"
file_name_2 = "fft.out"


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

set grid

point_size = 0.7
set style line 1 lt 1 lw 1 lc rgb "blue"  pointsize point_size
set style line 2 lt 1 lw 1 lc rgb "green" pointsize point_size
set style line 3 lt 1 lw 1 lc rgb "red"   pointsize point_size

if (ps) set output "fft_1.".(ext)

set multiplot

set size 0.5, 0.5
set origin 0.0, 0.5
set title "Horisontal Phase Space"
set xlabel "x [mm]"
set ylabel "p_x [mrad]"
plot "track.out" using (1e3*\$2):(1e3*\$3) notitle with points ls 1

set origin 0.5, 0.5
set title "Vertical Phase Space"
set xlabel "y [mm]"
set ylabel "p_y [mrad]"
plot "track.out" using (1e3*\$4):(1e3*\$5) notitle with points ls 3

set origin 0.5, 0.0
set title "Longitudinal Phase Space"
set xlabel "ct [mm]"
set ylabel "{/Symbol d [%]}"
plot "track.out" using (1e3*\$7):(1e2*\$6) notitle with points ls 2
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output "fft_2.".(ext)

set multiplot

set size 1.0, 0.5
set origin 0.0, 0.5
if (ps) set output "fft_2.".(ext)
set title "Time of Flight"
set xlabel "Turn Number"
set ylabel "ct [mm]"
plot "track.out" using 1:7 notitle with impulses ls 2

set size 1.0, 0.5
set origin 0.0, 0.0
set title "FFT of 2J_y"
set xlabel "{/Symbol n}_y"
set ylabel "A_y"
plot "fft.out" using 1:2 notitle with impulses ls 2
if (!ps) pause mouse "click on graph to cont.\n"

EOP
