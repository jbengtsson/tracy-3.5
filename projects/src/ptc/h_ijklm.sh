#!/bin/sh

prm1=${1-0}
prm2=${2-"h_ijklm"}

gnuplot << EOP

ps        = $prm1;
file_name = "$prm2";
contour   = 0;

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

if (contour) \
  set nosurface; \
  set noztics;

#set cntrlabel format "%4.1f";
set key left;
#unset cntrlabel;
#set noztics;
unset colorbox;
#set cbrange [0:1];

#set cntrparam level 18;
set cntrparam level 25;

# x <-> horizontal, y <-> vertical, z <-> perpendicular to screen.
# rot_x, rot_z, scale, scale_z.
if (!contour) set view 65, 15, 1, 1;
if (contour) set view 0, 0, 1, 1;

set palette rgbformulae 22, 13, -31 negative;

if (ps) set output file_name."_1.".(ext);
if (ps) set output "tune_scan.ps"

set title "h_{ijklm}";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";

splot file_name.".dat" using 1:2:(log(abs(\$4))) notitle with lines lt \
      palette z;

if (!ps) pause mouse "click on graph to cont.\n";
