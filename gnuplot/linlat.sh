#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1; eps = 0; plt_nu = 1; plt_I5 = 1;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "linlat_1.ps"
set title "Beta Functions";
set xlabel "s [m]"; set ylabel "{/Symbol b} [m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:6 title "{/Symbol b}_x" with lines ls 1, \
     "linlat.out" using 3:11 title "{/Symbol b}_y" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_2.ps"
set title "Horizontal Dispersion";
set xlabel "s [m]"; set ylabel "{/Symbol h}_x [m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:8 notitle with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_3.ps"
set title "Vertical Dispersion";
set xlabel "s [m]"; set ylabel "{/Symbol h}_y [m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:13 notitle with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_4.ps";
if (plt_nu) \
  set title "{/Symbol m}"; \
  set xlabel "s [m]"; set ylabel ""; \
  set y2range [-1.5:20]; \
  plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
       lc rgb "black", \
       "linlat.out" using 3:7 title "{/Symbol n}_x" with lines ls 1, \
       "linlat.out" using 3:12 title "{/Symbol n}_y" with lines ls 3; \
  if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_5.ps";
if (plt_nu) \
  set title "{/Symbol a}"; \
  set xlabel "s [m]"; set ylabel ""; \
  set y2range [-1.5:20]; \
  plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
       lc rgb "black", \
       "linlat.out" using 3:5 title "{/Symbol a}_x" with lines ls 1, \
       "linlat.out" using 3:10 title "{/Symbol a}_y" with lines ls 3; \
  if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_6.ps"
set title "{/Symbol b}_{x,y}{/Symbol \264h}_x";
set xlabel "s [m]"; set ylabel "";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:(\$6*\$8) title "{/Symbol b}_x{/Symbol \264h}_x" \
     with lines ls 1, \
     "linlat.out" using 3:(\$11*\$8) title "{/Symbol b}_y{/Symbol \264h}_x" \
     with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

file_name = "`echo \$TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name using (\$1/255):(\$2/255):(\$3/255);
unset colorbox;
#set cbrange [-1.5:1.5];

if (ps) set output "linlat_7.ps";
set title "Floquet Space";
set xlabel "{/Symbol h}~"; set ylabel "{/Symbol h}'~";
plot "linlat.out" using 15:16:(abs(\$4)) notitle "{/Symbol n}_x" \
      with lines lt palette z;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
