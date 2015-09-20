ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "cod_rms_1.ps"
set title "Horizontal RMS Orbit \n \
({/Symbol D}x_{rms}=100 {/Symbol m}m, {/Symbol D}y_{rms}=100 {/Symbol m}m)";
set xlabel "s [m]"; set ylabel "[mm]";
set yrange [-0.01:*];
set y2range [-1.5:20];
plot "cod_rms.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "cod_rms.out" using 3:5:5:7 notitle with errorbars ls 1;
if (!ps) pause -1;

if (ps) set output "cod_rms_2.ps"
set title "Vertical RMS Orbit \n \
({/Symbol D}x_{rms}=100 {/Symbol m}m, {/Symbol D}y_{rms}=100 {/Symbol m}m)";
set xlabel "s [m]"; set ylabel "[mm]";
set yrange [-0.015:*];
set y2range [-1.5:20];
plot "cod_rms.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "cod_rms.out" using 3:8:8:10 notitle with errorbars ls 3;
if (!ps) pause -1;
