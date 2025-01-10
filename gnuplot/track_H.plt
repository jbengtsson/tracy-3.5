ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "track_H.ps"
set title "H";
set xlabel "Turn No";
set ylabel "[%]";
#set yrange[0:*];
plot "track_H.dat" using 1:2 notitle with impulses 3;
if (!ps) pause -1;

if (ps) set output "track_H.ps"
set title "H";
set xlabel "Turn No";
set ylabel "[%]";
#set yrange[0:*];
plot "track_H.dat" using 1:3 notitle with impulses 2;
if (!ps) pause -1;
