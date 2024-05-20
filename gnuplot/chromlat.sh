#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

file_name = "linlat.out";

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

# Symbol characters broken for > \219.
# In particular, {/Symbol \264h} = eta.

if (ps) set output "chromlat_1.".(ext)
# set title "Linear Chromaticity: {/Symbol b}_{x,y}{/Symbol \264h}_x";
set title "Linear Chromaticity: {/Symbol b}_{x,y} x eta_x";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:(\$6*\$8) title "{/Symbol b}_x x eta_x" \
     with lines ls 1, \
     file_name using 3:(\$11*\$8) title "{/Symbol b}_y x eta_x" \
     with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

exit;

if (ps) set output "chromlat_2.".(ext)
set title "d{/Symbol b}_{x,y}/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:14 \
     title "d{/Symbol  b}_x/d{/Symbol d}" \
     with lines ls 1, \
     file_name using 3:15 \
     title "d{/Symbol  b}_y/d{/Symbol d}" \
     with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "chromlat_3.".(ext)
set title "Second Order Dispersion: d{/Symbol h}_x/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:16 notitle with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

exit;

if (ps) set output "chromlat_4.".(ext)
set title "Linear Coupling: sqrt({/Symbol b}_x{/Symbol \264b}_y)";
set xlabel "s [m]"; set ylabel "[m]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:6 notitle  with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "chromlat_5.".(ext)
set title "Second Order Chromaticity: d{/Symbol b}_{x,y}/d{/Symbol d\264h}_x";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:7 \
     title "d{/Symbol b}_x/d{/Symbol d\264h}_x" with lines ls 1, \
     file_name using 3:10 \
     title "d{/Symbol b}_y/d{/Symbol d\264h}_x" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "chromlat_6.".(ext)
set title "Second Order Chromaticity: {/Symbol b}_{x,y}{/Symbol \264}d{/Symbol h}_x/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:8 \
     title "{/Symbol b}_x{/Symbol \264}d{/Symbol h}_x/d{/Symbol d}" \
     with lines ls 1, \
     file_name using 3:11 \
     title "{/Symbol b}_y{/Symbol \264}d{/Symbol h}_x/d{/Symbol d}" \
     with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "chromlat_7.".(ext)
set title "Second Order Chromaticity: {/Symbol b}_{x,y}{/Symbol \264}d{/Symbol h}_x/d{/Symbol d} + {/Symbol \264}d{/Symbol h}_x/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:(\$7+\$8) title "Hor" with lines ls 1, \
     file_name using 3:(\$10+\$11) title "Ver" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "chromlat_8.".(ext)
set title "Chromaticity";
set xlabel "s [m]"; set ylabel "";
set y2range [-1.5:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:12 title "{/Symbol x}_x" with lines ls 1, \
     file_name using 3:13 title "{/Symbol x}_y" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
