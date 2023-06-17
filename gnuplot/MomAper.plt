ps = 0; eps = 0; stat = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "MomAper.ps"
set title "Momentum Aperture";
#set title "Momentum Aperture ({/Symbol t} = 2.6 hrs ())";
set xlabel "s [m]";
set ylabel "{/Symbol d} [%]";
set y2range [-20:20];

plot "linlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "mom_aper.out" using 2:3 notitle with lines ls 1, \
     "mom_aper.out" using 2:4 notitle with lines ls 1;

if (!ps) pause -1;
