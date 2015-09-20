ps = 0; eps = 0;

font_size = 30; line_width = 2;
if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid \
  lw line_width "Times-Roman" font_size;
if (ps && eps) \
  set terminal postscript eps enhanced color solid \
  lw line_width "Times-Roman" font_size;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "cod_1.ps"
set title "Horizontal Closed Orbit";
set xlabel "s [m]";
set ylabel "x [mm]";
set y2range [-1.5:20];
plot "cod.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "cod.out" using 3:9 notitle with impulses ls 1;
if (!ps) pause -1;

if (ps) set output "cod_2.ps"
set title "Vertical Closed Orbit";
set xlabel "s [m]";
set ylabel "y [mm]";
set y2range [-1.5:20];
plot "cod.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "cod.out" using 3:10 notitle with impulses ls 3;
if (!ps) pause -1;

if (ps) set output "cod_3.ps"
set title "Horizontal Corrector Strength";
set xlabel "s [m]";
set ylabel "{\Symbol t}_x [mrad]";
set y2range [-1.5:20];
plot "cod.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "cod.out" using 3:13 notitle with impulses ls 1;
if (!ps) pause -1;

if (ps) set output "cod_4.ps"
set title "Vertical Corrector Strength";
set xlabel "s [m]";
set ylabel "{\Symbol t}_y [mrad]";
set y2range [-1.5:20];
plot "cod.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "cod.out" using 3:14 notitle with impulses ls 3;
if (!ps) pause -1;
