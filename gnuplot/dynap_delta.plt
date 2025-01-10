ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "dynap_delta.ps"
set title "Dynamical Aperture";
set xlabel "x [mm]";
set ylabel "y [mm]";
#set xrange[-40:40];
#set yrange[-20:20];
plot "dynap_0.dat" using 1:2 title "No Cavity/Rad" with linespoints 1, \
     "dynap_delta_0.dat" using 1:2 title "{/Symbol d}=0" with linespoints 2, \
     "dynap_delta_1.dat" using 1:2 title "{/Symbol d}=1%" with linespoints 3, \
     "dynap_delta_2.dat" using 1:2 title "{/Symbol d}=2%" with linespoints 4;
if (!ps) pause -1;
