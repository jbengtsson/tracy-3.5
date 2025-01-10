ps = 0; eps = 1;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "sigma_1.ps"
set title "RMS Beam Size";
set xlabel "s [m]"; set ylabel "[mm]";
set y2range [-1.5:20];
plot "sigma.out" using 3:4 axis x1y2 notitle with fsteps 3, \
     "sigma.out" using 3:5 title "x" with lines 2, \
     "sigma.out" using 3:7 title "y" with lines 1;
if (!ps) pause -1;

if (ps) set output "sigma_2.ps"
set title "RMS Beam Angular Divergence";
set xlabel "s [m]"; set ylabel "[mrad]";
set y2range [-1.5:20];
plot "sigma.out" using 3:4 axis x1y2 notitle with fsteps 3, \
     "sigma.out" using 3:6 title "x" with lines 2, \
     "sigma.out" using 3:8 title "y" with lines 1;
if (!ps) pause -1;
