ps = 0; eps = 0; plt_nu = 1; plt_I5 = 1;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "dlinlat_1.ps"
set title "Beta Beat";
set xlabel "s [m]";
set ylabel "[%]";
set y2range [-1.5:20];
plot "dlinlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "dlinlat.out" using 3:5 title "{/Symbol Db}_x/{/Symbol b}_x" \
     with lines ls 1, \
     "dlinlat.out" using 3:7 title "{/Symbol Db}_y/{/Symbol b}_y" \
     with lines ls 3;
if (!ps) pause -1;

if (ps) set output "dlinlat_2.ps";
if (plt_nu) \
  set title "Normalized Phase Beat"; \
  set xlabel "s [m]"; \
  set ylabel ""; \
  set y2range [-1.5:20]; \
  plot "dlinlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
       lc rgb "black", \
       "dlinlat.out" using 3:6 title "{/Symbol n}_x" with lines ls 1, \
       "dlinlat.out" using 3:8 title "{/Symbol n}_y" with lines ls 3; \
  if (!ps) pause -1;
