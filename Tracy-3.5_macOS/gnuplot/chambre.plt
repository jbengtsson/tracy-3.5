ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "chambre_1.ps"
set title "Horizontal Mechanical Aperture";
set xlabel "s [m]";
set ylabel "x [mm]";
plot "chambre.out" using 3:4 notitle with fsteps 2, \
     "chambre.out" using 3:5 notitle with fsteps 2;
if (!ps) pause -1;

if (ps) set output "chambre_2.ps"
set title "Vertical Mechanical Aperture";
set xlabel "s [m]";
set ylabel "y [mm]";
set yrange [0:];
plot "chambre.out" using 3:6 notitle with fsteps 2;
if (!ps) pause -1;
