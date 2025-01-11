#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

file_name = "FieldMap_pass.dat";

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

if (ps) set output "fm_track_1.".(ext)
set title "ID Trajectory";
set xlabel "s [m]"; set ylabel "x [mm]";
plot file_name using 2:(1e3*\$3) notitle w lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "fm_track_2.".(ext)
set title "ID Trajectory";
set xlabel "s [m]"; set ylabel "p_x [mrad]";
plot file_name using 2:(1e3*\$4) notitle w lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "fm_track_3.".(ext)
set title "ID Trajectory";
set xlabel "s [m]"; set ylabel "y [mm]";
plot file_name using 2:(1e3*\$5) notitle w lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "fm_track_4.".(ext)
set title "ID Trajectory";
set xlabel "s [m]"; set ylabel "p_y [mrad]";
plot file_name using 2:(1e3*\$6) notitle w lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "fm_track_5.".(ext)
set title "ID Trajectory";
set xlabel "s [m]"; set ylabel "delta";
plot file_name using 2:7 notitle w lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "fm_track_4.".(ext)
set title "ID Trajectory";
set xlabel "s [m]"; set ylabel "ct [mm]";
plot file_name using 2:(1e3*\$8) notitle w lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
