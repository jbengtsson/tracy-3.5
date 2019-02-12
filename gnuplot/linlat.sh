#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1; plt_I5 = 1;

file_name = "linlat.out";

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

if (ps) set output "linlat_1.".(ext);
set title "Beta Functions";
set xlabel "s [m]"; set ylabel "{/Symbol b} [m]";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:6 title "{/Symbol b}_x" with lines ls 1, \
     file_name using 3:11 title "{/Symbol b}_y" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_2.".(ext);
set title "Dispersion";
set xlabel "s [m]"; set ylabel "{/Symbol h} [m]";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:8 title "{/Symbol h}_x" with lines ls 1, \
     file_name using 3:13 title "{/Symbol h}_y" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_3.".(ext);
set title "{/Symbol b}_{x,y}{/Symbol \264h}_x";
set xlabel "s [m]"; set ylabel "";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:(\$6*\$8) title "{/Symbol b}_x{/Symbol \264h}_x" \
     with lines ls 1, \
     file_name using 3:(\$11*\$8) title "{/Symbol b}_y{/Symbol \264h}_x" \
     with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

file_name_palette = "`echo \$TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name_palette \
  using (\$1/255):(\$2/255):(\$3/255);
unset colorbox;
#set cbrange [-2.0:1.5];

if (ps) set output "linlat_4.".(ext);
set title "{/ZapfChancery-MediumItalic H}_x({/Symbol h}_x\\\~, {/Symbol h}\'_x\\\~)";
set xlabel "{/Symbol h}_x\\\~ [10^{-3}]";
set ylabel "{/Symbol h}\'_x\\\~ [10^{-3}]";
set size square;
#set xrange [0:*]
plot file_name using (1e3*\$15):(1e3*\$16):(abs(\$4)) notitle "{/Symbol n}_x" \
     with lines lt palette z;
if (!ps) pause mouse "click on graph to cont.\n";

set size nosquare;

exit;

if (ps) set output "linlat_5.".(ext);
set title "{/Symbol h}_x [m]";
set xlabel "s [m]"; set ylabel "";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:8 notitle with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_6.".(ext);
set title "{/Symbol h}'_x";
set xlabel "s [m]"; set ylabel "";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:9 notitle with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_7.".(ext);
set title "eta{_x\\\~";
set xlabel "s [m]"; set ylabel "";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:15 title "eta_x\\\~" with lines ls 1, \
     file_name using 3:16 title "eta'_x\\\~" with lines ls 2;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_8.".(ext);
set title "{/ZapfChancery-MediumItalic H}_x(s)";
set xlabel "s [m]"; set ylabel "";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:(\$15**2+\$16**2) notitle with lines ls 1;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_9.".(ext);
set title "arg\\\{{/ZapfChancery-MediumItalic H}_x(s)\\\}";
set xlabel "s [m]"; set ylabel "";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:(-atan2(\$16, \$15)*180.0/pi) \
     title "arg\\\{curly\\\_H_x\\\}" with lines ls 1, \
     file_name using 3:(\$7*360.0) title "arg\\\{J_x\\\}" \
     with lines ls 2;
if (!ps) pause mouse "click on graph to cont.\n";

exit;

if (ps) set output "linlat_10.".(ext);
set title "Normalized Phase Advance";
set xlabel "s [m]"; set ylabel "{/Symbol n}";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:7 title "{/Symbol n}_x" with lines ls 1, \
     file_name using 3:12 title "{/Symbol n}_y" with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "linlat_11.".(ext);
set title "{/Symbol g}";
set xlabel "s [m]"; set ylabel "{/Symbol g}";
set y2range [-2.0:20];
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:((1.0+\$5**2)/\$6) title "{/Symbol g}_x" \
     with lines ls 1, \
     file_name using 3:((1.0+\$10**2)/\$11) title "{/Symbol g}_y" \
     with lines ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

exit;

if (ps) set output "linlat_12.".(ext);
 set title "{/Symbol a}"; \
set xlabel "s [m]"; set ylabel "{/Symbol a}"; \
set y2range [-2.0:20]; \
plot file_name using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 3:5 title "{/Symbol a}_x" with lines ls 1, \
     file_name using 3:10 title "{/Symbol a}_y" with lines ls 3; \
if (!ps) pause mouse "click on graph to cont.\n";

EOP
