#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

file_name_1 = "beta_nl.out";
file_name_2 = "eta_nl.out";

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
set style line 4 lt 1 lw 1 lc rgb "cyan";

if (ps) set output "plt_nl_optics_1.".(ext);
set title "Linear Optics";
set xlabel "s [m]"; set ylabel "{/Symbol b} [m]";
set y2range [-2.0:20];
plot file_name_2 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_1 using 3:5 title "{/Symbol b}_x" with lines ls 1, \
     file_name_1 using 3:6 title "{/Symbol b}_y" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "plt_nl_optics_2.".(ext);
set title "Linear Dispersion";
set xlabel "s [m]"; set ylabel "{/Symbol h} [m]";
set y2range [-2.0:20];
plot file_name_2 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_2 using 3:5 title "{/Symbol h}_x" with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "plt_nl_optics_3.".(ext);
set title "Nonlinear Optics";
set xlabel "s [m]"; set ylabel "{/Symbol h} [m]";
set y2range [-2.0:20];
plot file_name_2 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_1 using 3:7 title "{/Symbol b}_x^{(1)}" with lines ls 1, \
     file_name_1 using 3:8 title "{/Symbol b}_y^{(1)}" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "plt_nl_optics_4.".(ext);
set title "Nonlinear Dispersion";
set xlabel "s [m]"; set ylabel "{/Symbol h} [m]";
set y2range [-2.0:20];
plot file_name_2 using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name_2 using 3:7 title "{/Symbol h}_x^{(2)}" with lines ls 1, \
     file_name_2 using 3:9 title "{/Symbol h}_x^{(3)}" with lines ls 2;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
