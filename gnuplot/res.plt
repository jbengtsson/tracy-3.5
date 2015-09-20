ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

set parametric;

if (ps) set output "res.ps"
set title "Resonances";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";

set xrange [-0.1:1.1]; set yrange [-0.1:1.1];
#set xrange [0.0:0.5]; set yrange [0.5:1.0];
#set xrange [0.35:0.43]; set yrange [0.58:0.65];
#set xrange [0.80:0.86]; set yrange [0.38:0.45];

plot 0, t title "1" with lines 1, \
     1, t notitle with lines 1, \
     t, 0 notitle with lines 1, \
     t, 1 notitle with lines 1, \
\
     0.0/2.0, t title "2" with lines 2, \
     1.0/2.0, t notitle with lines 2, \
     2.0/2.0, t notitle with lines 2, \
     t, 0.0/2.0 notitle with lines 2, \
     t, 1.0/2.0 notitle with lines 2, \
     t, 2.0/2.0 notitle with lines 2, \
     t, t notitle with lines 2, \
     t, 1-t notitle with lines 2, \
\
     1.0/3.0, t title "3" with lines 3, \
     2.0/3.0, t notitle with lines 3, \
     3.0/3.0, t notitle with lines 3, \
     t, 2*t notitle with lines 3, \
     t, -1+2*t notitle with lines 3, \
     t, 1-2*t notitle with lines 3, \
     t, 2-2*t notitle with lines 3, \
\
     0.0/4.0, t title "4" with lines 4, \
     1.0/4.0, t notitle with lines 4, \
     2.0/4.0, t notitle with lines 4, \
     3.0/4.0, t notitle with lines 4, \
     4.0/4.0, t notitle with lines 4, \
     t, 0.0/4.0 notitle with lines 4, \
     t, 1.0/4.0 notitle with lines 4, \
     t, 2.0/4.0 notitle with lines 4, \
     t, 3.0/4.0 notitle with lines 4, \
     t, 4.0/4.0 notitle with lines 4, \
     2*t, 2*t notitle with lines 4, \
     2*t, 1-2*t notitle with lines 4, \
\
     0.0/5.0, t title "5" with lines 5, \
     1.0/5.0, t notitle with lines 5, \
     2.0/5.0, t notitle with lines 5, \
     3.0/5.0, t notitle with lines 5, \
     4.0/5.0, t notitle with lines 5, \
     5.0/5.0, t notitle with lines 5, \
     t, 0.0/5.0 notitle with lines 5, \
     t, 1.0/5.0 notitle with lines 5, \
     t, 2.0/5.0 notitle with lines 5, \
     t, 3.0/5.0 notitle with lines 5, \
     t, 4.0/5.0 notitle with lines 5, \
     t, 5.0/5.0 notitle with lines 5, \
     3*t, 2*t notitle with lines 5, \
     3*t, 1-2*t notitle with lines 5, \
     t, 4*t notitle with lines 5, \
     t, -1+4*t notitle with lines 5, \
     t, -2+4*t notitle with lines 5, \
     t, -3+4*t notitle with lines 5, \
     t, 1-4*t notitle with lines 5, \
     t, 2-4*t notitle with lines 5, \
     t, 3-4*t notitle with lines 5, \
     t, 4-4*t notitle with lines 5, \
\
     0.4, 0.6 notitle with points 3, \
     0.377, 0.611 notitle with points 3, \
\
     "dnu_dAx.out" using 5:6 title "{/Symbol nu}(A_x)" with linespoints 2, \
     "dnu_dAy.out" using 5:6 title "{/Symbol nu}(A_y)" with linespoints 1;

if (!ps) pause -1;
