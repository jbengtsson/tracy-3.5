ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "ID_corr_1.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol Db}_x/{/Symbol b}_x [%]";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "ID_corr_res.out" using 1:2 notitle with impulses ls 1;

set origin 0.0, 0.0;
set title "{/Symbol Db}_y/{/Symbol b}_y [%]";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "ID_corr_res.out" using 1:3 notitle with impulses ls 3;

unset multiplot;

if (!ps) pause -1;

if (ps) set output "ID_corr_2.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol Dn}_x";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "ID_corr_res.out" using 1:4 notitle with impulses ls 1;

set origin 0.0, 0.0;
set title "{/Symbol Dn}_y";
set xlabel "s [m]";
set ylabel "";
set y2range [-1.5:20];
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "ID_corr_res.out" using 1:5 notitle with impulses ls 3;

unset multiplot;

if (!ps) pause -1;

if (ps) set output "ID_corr_3.ps"

set size 1.0, 1.0; set origin 0.0, 0.0;
set title "b_{2}L";
set xlabel "s [m]";
set ylabel "";
plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "ID_corr.out" using 2:4 notitle with impulses ls 1;

unset multiplot;

if (!ps) pause -1;
