#!/bin/sh

prm1=${1-""}
prm2=${2-0}

echo $prm1
echo $prm2

gnuplot << EOP

ps = $prm2;
file_name = "$prm1";

f_s = 24; l_w = 2;
if (ps == 0) \
  set terminal x11; \
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

set view 65, 15, 1, 1;
set palette rgbformulae 22, 13, -31;
set clabel "%3.1f"; set key left;
unset colorbox;
set zlabel "";

if (ps) set output "kick_map.".(ext)
set title "Kick Map: {/Symbol q}_x [mrad]";
set xlabel "x [mm]"; set ylabel "y [mm]";
splot file_name using (1e3*\$1):(1e3*\$2):(1e3*\$3):2 notitle \
      w lines lt palette z;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
