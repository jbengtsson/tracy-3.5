ps = 1; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "track_1.ps"
set title "Horizontal Phase Space";
set xlabel "x [mm]";
set ylabel "px [mrad]";
#plot "track.dat" using 2:3 notitle with points 1;
#if (!ps) pause -1;

if (ps) set output "track_2.ps"
set title "Vertical Phase Space";
set xlabel "y [mm]";
set ylabel "py [mrad]";
#plot "track.dat" using 4:5 notitle with points 1;
#if (!ps) pause -1;

if (ps) set output "track_3.ps"
set title "Longitudinal Phase Space";
set xlabel "[fsec]";
set ylabel "{/Symbol d} [%]";
plot "iso_0.dat" using 3:2 notitle with points 1, \
     "iso_1.dat" using 3:2 notitle with points 2;
if (!ps) pause -1;
