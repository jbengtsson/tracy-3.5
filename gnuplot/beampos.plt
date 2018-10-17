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

if (ps) set output "beampos_1.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Horizontal Trajectory";
set xlabel "s [m]";
set ylabel "x [mm]";
set y2range [-1.5:20];
plot "beampos.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beampos.out" using 3:(1e3*$5) notitle with lines ls 1;

set origin 0.0, 0.0;
set title "Vertical Trajectory";
set xlabel "s [m]";
set ylabel "y [mm]";
set y2range [-1.5:20];
plot "beampos.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beampos.out" using 3:(1e3*$6) notitle with lines ls 3;

unset multiplot;
if (!ps) pause(-1);

if (ps) set output "beampos_2.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Horizontal Trajectory";
set xlabel "s [m]";
set ylabel "x [mm]";
set y2range [-1.5:20];
plot "beampos_corr.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
    "beampos_corr.out" using 3:(1e3*$5) notitle with lines ls 1;

set origin 0.0, 0.0;
set title "Vertical Trajectory";
set xlabel "s [m]";
set ylabel "y [mm]";
set y2range [-1.5:20];
plot "beampos_corr.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beampos_corr.out" using 3:(1e3*$6) notitle with lines ls 3;

unset multiplot;
if (!ps) pause(-1);
