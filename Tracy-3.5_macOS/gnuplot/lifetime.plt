ps = 1; eps = 1;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid lw 2 "Times-Roman" 25;
if (ps && eps) set terminal postscript enhanced color solid lw 2 "Times-Roman" 25;
#if (ps && !eps) set terminal postscript enhanced color solid;
#if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "lifetimes.ps"
set title "Touschek lifetime vs. CPMU vertical half gap";
set xlabel "half gap size [mm]";
set ylabel "{/Symbol t} (hrs)";
#set yrange [0:];

plot "lifetimes.dat" using 1:2 notitle with linespoints lw 2 pointsize 3 lc rgb "red";


if (!ps) pause -1;
