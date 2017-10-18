ps = 0;

file_name = "drv_terms.out";

f_s = 14; l_w = 2;
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
set style line 5 lt 1 lw 1 lc rgb "purple";
set style line 6 lt 1 lw 1 lc rgb "dark-orange";

if (ps) set output "drv_terms_1.".(ext);

set title "Chromatic Driving Terms";
set xlabel "s [m]"; set ylabel "Re{h_ijklm}";
plot file_name using 2:3  title "h_{10002}" with steps ls 1, \
     file_name using 2:4  title "h_{20001}" with steps ls 2, \
     file_name using 2:5  title "h_{00201}" with steps ls 3;

if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "drv_terms_2.".(ext);

set title "Driving Terms";
set xlabel "s [m]"; set ylabel "Re{h_ijklm}";
plot file_name using 2:6  title "h_{21000}" with steps ls 1, \
     file_name using 2:7  title "h_{10110}" with steps ls 2, \
     file_name using 2:8  title "h_{30000}" with steps ls 3, \
     file_name using 2:9  title "h_{10200}" with steps ls 4, \
     file_name using 2:10 title "h_{10020}" with steps ls 5;

if (!ps) pause mouse "click on graph to cont.\n";
