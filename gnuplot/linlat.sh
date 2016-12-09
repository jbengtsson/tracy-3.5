#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1; plt_nu = 1; plt_I5 = 1;

if (ps == 0) set terminal x11;
if (ps == 1) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 18;
if (ps == 2) \
  set terminal pdf enhanced color solid linewidth 2 font "Times-Roman, 18";

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps == 1) set output "linlat_1.eps"; \
else if (ps == 2) set output "linlat_1.pdf"
set title "Beta Functions";
set xlabel "s [m]"; set ylabel "{/Symbol b} [m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:6 title "{/Symbol b}_x" with lines ls 1, \
     "linlat.out" using 3:11 title "{/Symbol b}_y" with lines ls 3;
if (ps == 0) pause mouse "click on graph to cont.\n";

if (ps == 1) set output "linlat_2.eps"; \
else if (ps == 2) set output "linlat_2.pdf"
set title "Dispersion";
set xlabel "s [m]"; set ylabel "{/Symbol h} [m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:8 title "{/Symbol h}_x" with lines ls 1, \
     "linlat.out" using 3:13 title "{/Symbol h}_y" with lines ls 3;
if (ps == 0) pause mouse "click on graph to cont.\n";

if (ps == 1) set output "linlat_3.eps"; \
else if (ps == 2) set output "linlat_3.pdf"
if (plt_nu) \
  set title "Normalized Phase Advance"; \
  set xlabel "s [m]"; set ylabel "{/Symbol n}"; \
  set y2range [-1.5:20]; \
  plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
       lc rgb "black", \
       "linlat.out" using 3:7 title "{/Symbol n}_x" with lines ls 1, \
       "linlat.out" using 3:12 title "{/Symbol n}_y" with lines ls 3; \
  if (ps == 0) pause mouse "click on graph to cont.\n";

if (ps == 1) set output "linlat_4.eps"; \
else if (ps == 2) set output "linlat_4.pdf"
if (plt_nu) \
  set title "{/Symbol a}"; \
  set xlabel "s [m]"; set ylabel "{/Symbol a}"; \
  set y2range [-1.5:20]; \
  plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
       lc rgb "black", \
       "linlat.out" using 3:5 title "{/Symbol a}_x" with lines ls 1, \
       "linlat.out" using 3:10 title "{/Symbol a}_y" with lines ls 3; \
  if (ps == 0) pause mouse "click on graph to cont.\n";

if (ps == 1) set output "linlat_5.eps"; \
else if (ps == 2) set output "linlat_5.pdf"
set title "{/Symbol b}_{x,y}{/Symbol \264h}_x";
set xlabel "s [m]"; set ylabel "";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:(\$6*\$8) title "{/Symbol b}_x{/Symbol \264h}_x" \
     with lines ls 1, \
     "linlat.out" using 3:(\$11*\$8) title "{/Symbol b}_y{/Symbol \264h}_x" \
     with lines ls 3;
if (ps == 0) pause mouse "click on graph to cont.\n";

file_name = "`echo \$TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name using (\$1/255):(\$2/255):(\$3/255);
unset colorbox;
#set cbrange [-1.5:1.5];

if (ps == 1) set output "linlat_6.eps"; \
else if (ps == 2) set output "linlat_6.pdf"
set title "Floquet Space";
set xlabel "{/Symbol h}~"; set ylabel "{/Symbol h}'~";
plot "linlat.out" using 15:16:(abs(\$4)) notitle "{/Symbol n}_x" \
      with lines lt palette z;
if (ps == 0) pause mouse "click on graph to cont.\n";

EOP
