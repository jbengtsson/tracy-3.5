ps = 0; eps = 1; action_angle = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

if (ps) set output "dnu_dA_1.ps"

set title "{/Symbol n}_x vs. A_{x,y}";
if (!action_angle) \
  set xlabel "A_{x,y} [mm]"; set ylabel ""; \
  plot "dnu_dAx.out" using 1:5 title "A_x" with linespoints 2, \
       "dnu_dAy.out" using 2:5 title "A_y" with linespoints 1;
#  plot "dnu_dAx.out" using 1:5 title "A_x (num.)" with linespoints ls 3, \
#       "/home/bengtsson/Thor/src/dnu_dAx.dat" using 1:3 \
#       title "A_x (pert. theory)" with linespoints ls 3, \
#       "dnu_dAy.out" using 2:5 title "A_y (num.l)" with linespoints ls 1, \
#       "/home/bengtsson/Thor/src/dnu_dAy.dat" using 1:3 \
#       title "A_y (pert. theory)" with linespoints ls 4;
if (action_angle) \
  set xlabel "A_{x,y} [mm]"; set ylabel "{/Symbol n}_x"; \
  set ytics nomirror; set y2tics; \
  plot "dnu_dAx.out" using 3:5 title "J_x" with linespoints 2;
#  plot "dnu_dAx.out" using 3:5 title "J_x" with linespoints ls 1, \
#       "dnu_dAy.out" using 4:5 title "J_y" with linespoints ls 3;

if (!ps) pause -1;

if (ps) set output "dnu_dA_2.ps"

  set title "Amplitude Dependent Tune Shift";
if (!action_angle) \
  set xlabel "A_{x,y} [mm]"; set ylabel "{/Symbol n}_y"; \
  plot "dnu_dAx.out" using 1:6 title "A_x" with linespoints 2, \
       "dnu_dAy.out" using 2:6 title "A_y" with linespoints 1;
#  plot "dnu_dAx.out" using 1:6 title "A_x (num.)" with linespoints  ls1, \
#       "/home/bengtsson/Thor/src/dnu_dAx.dat" using 1:4 \
#       title "A_x (pert. theory)" with linespoints 1, \
#       "dnu_dAy.out" using 2:6 title "A_y (num.)" with linespoints ls 3, \
#       "/home/bengtsson/Thor/src/dnu_dAy.dat" using 1:4 \
#       title "A_y (pert. theory)" with linespoints ls 3;
if (action_angle) \
  set xlabel "A_{x,y} [mm]"; set ylabel "{/Symbol n}_y"; \
  plot "dnu_dAx.out" using 3:6 title "J_x" with linespoints 2, \
       "dnu_dAy.out" using 4:6 title "J_y" with linespoints 1;

if (!ps) pause -1;
