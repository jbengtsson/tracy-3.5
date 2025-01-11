ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

set yrange[-50:50];

if (ps) set output "aperture_1.ps"
set title "Horizontal Mechanical Aperture";
set xlabel "s [m]";
set ylabel "x [mm]";
plot "aperture.out" using 1:2 notitle with fsteps 2, \
     "aperture.out" using 1:3 notitle with fsteps 2;
if (!ps) pause -1;

if (ps) set output "aperture_2.ps"
set title "Vertical Mechanical Aperture";
set xlabel "s [m]";
set ylabel "y [mm]";
plot "aperture.out" using 1:4 notitle with fsteps 2, \
     "aperture.out" using 1:5 notitle with fsteps 2;
if (!ps) pause -1;
