ps = 0; eps = 0; plt_nu = 1; plt_I5 = 1;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "linlat_err_1.ps"
set title "Beta Functions (bare)";
set xlabel "s [m]";
set ylabel "[m]";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat.out" using 3:6 title "{/Symbol b}_x" with lines ls 1, \
     "linlat.out" using 3:11 title "{/Symbol b}_y" with lines ls 3;
if (!ps) pause -1;

if (ps) set output "linlat_err_2.ps"
set title "Beta Functions (w IDs)";
set xlabel "s [m]";
set ylabel "[m]";
set y2range [-1.5:20];
#plot "linlat_IDs.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
#     lc rgb "black", \
#     "linlat_IDs.out" using 3:6 title "{/Symbol b}_x" with lines ls 1, \
#     "linlat_IDs.out" using 3:11 title "{/Symbol b}_y" with lines ls 3;
#if (!ps) pause -1;

if (ps) set output "linlat_err_3.ps"
set title "Beta Functions (corrected)";
set xlabel "s [m]";
set ylabel "[m]";
set y2range [-1.5:20];
#plot "linlat_IDs_corr.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
#     lc rgb "black", \
#     "linlat_IDs_corr.out" using 3:6 title "{/Symbol b}_x" with lines ls 1, \
#     "linlat_IDs_corr.out" using 3:11 title "{/Symbol b}_y" with lines ls 3;
#if (!ps) pause -1;

if (ps) set output "linlat_err_4.ps"
set title "Horizontal Dispersion (corrected)";
set xlabel "s [m]";
set ylabel "{/Symbol h}_y [m]";
set y2range [-1.5:20];
#plot "linlat_IDs_corr.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
#     lc rgb "black", \
#     "linlat_IDs_corr.out" using 3:8 notitle with lines ls 1;
#if (!ps) pause -1;

if (ps) set output "linlat_err_5.ps"
set title "Beta Functions (/w errors)";
set xlabel "s [m]";
set ylabel "[m]";
set y2range [-1.5:20];
plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat_err.out" using 3:6 title "{/Symbol b}_x" with lines ls 1, \
     "linlat_err.out" using 3:11 title "{/Symbol b}_y" with lines ls 3;
if (!ps) pause -1;

if (ps) set output "linlat_err_6.ps"
set title "Vertical Dispersion (/w errors)";
set xlabel "s [m]";
set ylabel "{/Symbol h}_x [m]";
set y2range [-1.5:20];
plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "linlat_err.out" using 3:13 notitle with lines ls 3;
if (!ps) pause -1;

