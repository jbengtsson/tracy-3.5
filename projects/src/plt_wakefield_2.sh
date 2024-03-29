#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1

file_name = "moments.out"

f_s = 24
 l_w = 2
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

point_size = 0.1
set style line 1 lt 1 lw 1 lc rgb "blue"  pointsize point_size
set style line 2 lt 1 lw 1 lc rgb "green" pointsize point_size
set style line 3 lt 1 lw 1 lc rgb "red"   pointsize point_size

if (ps) set output "plt_wakefield_2_1.".(ext)

set multiplot

set size 1.0, 0.5
set origin 0.0, 0.5
set key right bottom
set title "<x> [{/Symbol m}m]"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:(1e6*\$2) notitle with linespoints ls 1

set origin 0.0, 0.0
set title "<y> [{/Symbol m}m]"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:(1e6*\$4) notitle with linespoints ls 3

unset multiplot
if (!ps) pause mouse "click on graph to cont.\n"


if (ps) set output "plt_wakefield_2_2.".(ext)

set multiplot

set origin 0.0, 0.5
set title "<s> [ps]"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:(1e12*\$7) notitle with linespoints ls 1

set origin 0.0, 0.0
set title "<{/Symbol d}>"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:6 notitle with linespoints ls 3

unset multiplot
if (!ps) pause mouse "click on graph to cont.\n"

exit

if (ps) set output "plt_wakefield_2_3.".(ext)

set multiplot

set size 1.0, 0.5
set origin 0.0, 0.5
set key right bottom
set title "{/Symbol s}_x [{/Symbol m}m]"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:(1e6*\$8) notitle with points ls 1

set origin 0.0, 0.0
set title "{/Symbol s}_y [{/Symbol m}m]"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:(1e6*\$10) notitle with points ls 3

unset multiplot
if (!ps) pause mouse "click on graph to cont.\n"


if (ps) set output "plt_wakefield_4.".(ext)

set multiplot

set origin 0.0, 0.5
set title "{/Symbol s}_s [ps]"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:(1e12*\$13) notitle with points ls 1

set origin 0.0, 0.0
set title "{/Symbol s_d}"
set xlabel "turn no"
set ylabel ""
plot file_name using 1:12 notitle with points ls 3

unset multiplot
if (!ps) pause mouse "click on graph to cont.\n"

EOP

