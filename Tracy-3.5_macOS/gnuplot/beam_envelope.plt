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

if (ps) set output "beam_envelope_1.ps"
set title "{/Symbol e}^*_x";
set xlabel "s [m]"; set ylabel "[nm{/Symbol \327}rad]";
set yrange [0:];
set y2range [-1.5:20];

plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beam_envelope.out" using 3:(1e9*sqrt($4*$5-$6*$6)) \
     notitle with lines ls 1;

if (!ps) pause -1;

if (ps) set output "beam_envelope_2.ps"
set title "{/Symbol e}^*_y";
set xlabel "s [m]"; set ylabel "[nm{/Symbol \327}rad]";
set yrange [0:];
set y2range [-1.5:20];

plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beam_envelope.out" using 3:(1e9*sqrt($7*$8-$9*$9)) notitle \
     with lines ls 1;

if (!ps) pause -1;

if (ps) set output "beam_envelope_3.ps"
set title "{/Symbol s}_x";
set xlabel "s [m]"; set ylabel "[{/Symbol m}m]";
set yrange [0:];
set y2range [-1.5:20];

plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beam_envelope.out" using 3:(1e6*sqrt($4)) notitle with lines ls 1;

if (!ps) pause -1;

if (ps) set output "beam_envelope_4.ps"
set title "{/Symbol s}_{p_x}";
set xlabel "s [m]"; set ylabel "[{/Symbol m}rad]";
set yrange [0:];
set y2range [-1.5:20];

plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beam_envelope.out" using 3:(1e6*sqrt($5)) notitle with lines ls 1;

if (!ps) pause -1;

if (ps) set output "beam_envelope_5.ps"
set title "{/Symbol s}_y";
set xlabel "s [m]"; set ylabel "[{/Symbol m}m]";
set yrange [0:];
set y2range [-1.5:20];

plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beam_envelope.out" using 3:(1e6*sqrt($7)) notitle with lines ls 3;

if (!ps) pause -1;

if (ps) set output "beam_envelope_6.ps"
set title "{/Symbol s}_{p_y}";
set xlabel "s [m]"; set ylabel "[{/Symbol m}rad]";
set yrange [0:];
set y2range [-1.5:20];

plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beam_envelope.out" using 3:(1e6*sqrt($8)) notitle with lines ls 3;

if (!ps) pause -1;

if (ps) set output "beam_envelope_7.ps"
set title "Transverse Coupling Angle";
set xlabel "s [m]"; set ylabel "[{/Symbol \260}]";
set autoscale y;
set y2range [-1.5:20];

plot "linlat_err.out" using 3:4 axis x1y2 notitle with fsteps lt 1 lw 1 \
     lc rgb "black", \
     "beam_envelope.out" using 3:10 notitle with lines ls 1;

if (!ps) pause -1;
