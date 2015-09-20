ps = 0; eps = 1;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "eta_y.ps"
set title "Vertical Dispersion";
set xlabel "s [m]";
set ylabel "{/Symbol d} [mm]";

#set yrange [0:];
#plot "momentum_acceptance.CM.dat" using 1:(100*$2) notitle with lines 2, \
#     "Touschek.out" using 2:3 notitle with lines 3;
#plot "Touschek.out" using 1:2 notitle with lines 3, \
#     "Touschek.out" using 1:3 notitle with lines 3, \
#     "Touschek.out" using 1:($2+$4) notitle 5, \
#     "Touschek.out" using 1:($2-$4) notitle 5, \
#     "Touschek.out" using 1:($3+$5) notitle 5, \
#     "Touschek.out" using 1:($3-$5) notitle 5;

plot "eta_y.out" using 2:5 notitle with lines ls 3;

if (!ps) pause -1;
