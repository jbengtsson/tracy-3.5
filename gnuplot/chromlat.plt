ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "chromlat_1.ps"
set title "Linear Chromaticity: {/Symbol b}_{x,y}{/Symbol \264h}_x";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:5 title "{/Symbol b}_x{/Symbol \264h}_x" \
     with lines ls 1, \
     "chromlat.out" using 3:9 title "{/Symbol b}_y{/Symbol \264h}_x" \
     with lines ls 3;
if (!ps) pause -1;

if (ps) set output "chromlat_2.ps"
set title "Linear Coupling: sqrt({/Symbol b}_x{/Symbol \264b}_y)";
set xlabel "s [m]"; set ylabel "[m]";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:6 notitle  with lines ls 1;
if (!ps) pause -1;

if (ps) set output "chromlat_3.ps"
set title "Second Order Chromaticity: d{/Symbol b}_{x,y}/d{/Symbol d\264h}_x";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:7 \
     title "d{/Symbol b}_x/d{/Symbol d\264h}_x" with lines ls 1, \
     "chromlat.out" using 3:10 \
     title "d{/Symbol b}_y/d{/Symbol d\264h}_x" with lines ls 3;
if (!ps) pause -1;

if (ps) set output "chromlat_4.ps"
set title "Second Order Chromaticity: {/Symbol b}_{x,y}{/Symbol \264}d{/Symbol h}_x/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:8 \
     title "{/Symbol b}_x{/Symbol \264}d{/Symbol h}_x/d{/Symbol d}" \
     with lines ls 1, \
     "chromlat.out" using 3:11 \
     title "{/Symbol b}_y{/Symbol \264}d{/Symbol h}_x/d{/Symbol d}" \
     with lines ls 3;
if (!ps) pause -1;

if (ps) set output "chromlat_5.ps"
set title "Second Order Chromaticity: {/Symbol b}_{x,y}{/Symbol \264}d{/Symbol h}_x/d{/Symbol d} + {/Symbol \264}d{/Symbol h}_x/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m^2]";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:($7+$8) title "Hor" with lines ls 1, \
     "chromlat.out" using 3:($10+$11) title "Ver" with lines ls 3;
if (!ps) pause -1;

if (ps) set output "chromlat_6.ps"
set title "Chromaticity";
set xlabel "s [m]"; set ylabel "";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:12 title "{/Symbol x}_x" with lines ls 1, \
     "chromlat.out" using 3:13 title "{/Symbol x}_y" with lines ls 3;
if (!ps) pause -1;

if (ps) set output "chromlat_7.ps"
set title "d{/Symbol b}_{x,y}/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m]";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:14 \
     title "d{/Symbol  b}_x/d{/Symbol d}" \
     with lines ls 1, \
     "chromlat.out" using 3:15 \
     title "d{/Symbol  b}_y/d{/Symbol d}" \
     with lines ls 3;
if (!ps) pause -1;

if (ps) set output "chromlat_8.ps"
set title "Second Order Dispersion: d{/Symbol h}_x/d{/Symbol d}";
set xlabel "s [m]"; set ylabel "[m]";
set y2range [-1.5:20];
plot "chromlat.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "chromlat.out" using 3:16 notitle with lines ls 1;
if (!ps) pause -1;
