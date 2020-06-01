#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

file_name = "drv_terms.out";

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

# left adjusted labels
#set key Left;
set clabel "%5.2f"; set key left;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";
set style line 4 lt 1 lw 1 lc rgb "cyan";
set style line 5 lt 1 lw 1 lc rgb "purple";
set style line 6 lt 1 lw 1 lc rgb "dark-orange";
set style line 7 lt 1 lw 1 lc rgb "dark-green";
set style line 8 lt 1 lw 1 lc rgb "yellow";

if (ps) set output "drv_terms_1.".(ext);
set title "2nd Order Chromatic  Terms";
set xlabel "s [m]"; set ylabel "|h_{ijklm}|";
set y2range [-2.0:20];
set ytics nomirror;
#set logscale y;
plot file_name using 2:3 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 2:4 title "h_{10002}" with steps ls 1, \
     file_name using 2:5 title "h_{20001}" with steps ls 2, \
     file_name using 2:6 title "h_{00201}" with steps ls 3;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "drv_terms_2.".(ext);
set title "2nd Order Geometric Terms";
set xlabel "s [m]"; set ylabel "|h_{ijklm}|";
plot file_name using 2:3 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 2:7 title "h_{30000}" with steps ls 1, \
     file_name using 2:8 title "h_{21000}" with steps ls 2, \
     file_name using 2:9 title "h_{10200}" with steps ls 3, \
     file_name using 2:10 title "h_{10020}" with steps ls 4, \
     file_name using 2:11 title "h_{10110}" with steps ls 5;
if (!ps) pause mouse "click on graph to cont.\n";

exit;

if (ps) set output "drv_terms_3.".(ext);
set title "3rd Order Geometric Terms";
set xlabel "s [m]"; set ylabel "|h_{ijklm}|";
set y2range [-2.0:20];
plot file_name using 2:3 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 2:12 title "h_{40000}" with steps ls 1, \
     file_name using 2:13 title "h_{31000}" with steps ls 2, \
     file_name using 2:14 title "h_{20200}" with steps ls 3, \
     file_name using 2:15 title "h_{11200}" with steps ls 4, \
     file_name using 2:16 title "h_{00400}" with steps ls 5, \
     file_name using 2:17 title "h_{20110}" with steps ls 6, \
     file_name using 2:18 title "h_{00310}" with steps ls 7, \
     file_name using 2:19 title "h_{20020}" with steps ls 8;
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "drv_terms_4.".(ext);
set title "4th Order Geometric Terms";
set xlabel "s [m]"; set ylabel "|h_{ijklm}|";
set y2range [-2.0:20];
plot file_name using 2:3 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     file_name using 2:20 title "h_{50000}" with steps ls 1, \
     file_name using 2:21 title "h_{41000}" with steps ls 2, \
     file_name using 2:22 title "h_{32000}" with steps ls 3, \
     file_name using 2:23 title "h_{30200}" with steps ls 4, \
     file_name using 2:24 title "h_{21200}" with steps ls 5, \
     file_name using 2:25 title "h_{10400}" with steps ls 6, \
     file_name using 2:26 title "h_{30110}" with steps ls 7, \
     file_name using 2:27 title "h_{21110}" with steps ls 8, \
     file_name using 2:28 title "h_{10310}" with steps ls 1, \
     file_name using 2:29 title "h_{30020}" with steps ls 2, \
     file_name using 2:30 title "h_{21020}" with steps ls 3, \
     file_name using 2:31 title "h_{10220}" with steps ls 4, \
     file_name using 2:32 title "h_{10130}" with steps ls 5, \
     file_name using 2:34 title "h_{10040}" with steps ls 6, \
     file_name using 2:35 title "h_{01040}" with steps ls 7;
if (!ps) pause mouse "click on graph to cont.\n";

EOP
