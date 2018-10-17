ps = 0; eps = 0; plt_nu = 1; plt_I5 = 1;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

if (ps) set output "linlat_1.ps"
set title "Beta Functions";
set xlabel "s [m]";
set ylabel "[m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt -1, \
     "linlat.out" using 3:6 title "{/Symbol b}_x" with lines 2, \
     "linlat.out" using 3:11 title "{/Symbol b}_y" with lines 1;
if (!ps) pause -1;

if (ps) set output "linlat_2.ps"
set title "Horizontal Dispersion";
set xlabel "s [m]";
set ylabel "{/Symbol h}_x [m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt -1, \
     "linlat.out" using 3:8 notitle with lines 2;
if (!ps) pause -1;

if (ps) set output "linlat_3.ps"
set title "Vertical Dispersion";
set xlabel "s [m]";
set ylabel "{/Symbol h}_y [m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt -1, \
     "linlat.out" using 3:13 notitle with lines 2;
if (!ps) pause -1;

if (ps) set output "linlat_4.ps";
if (plt_nu) \
  set title "Phase Advance"; \
  set xlabel "s [m]"; \
  set ylabel ""; \
  set y2range [-1.5:20]; \
  plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt -1, \
       "linlat.out" using 3:7 title "{/Symbol n}_x" with lines 2, \
       "linlat.out" using 3:12 title "{/Symbol n}_y" with lines 1; \
  if (!ps) pause -1;

if (ps) set output "linlat_5.ps";
if (plt_I5) \
  set title "I5"; \
  set xlabel "s [m]"; \
  set ylabel ""; \
  set y2range [-1.5:20]; \
  plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt -1, \
       "linlat.out" using 3:15 title "{/Symbol n}_x" with lines 2; \
  if (!ps) pause -1;
