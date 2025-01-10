ps = 0; J_phi = 0; f_rf  = 0;

f_s = 24; l_w = 2;
if (ps == 0) \
  set terminal x11; \
else if (ps == 1) \
  set terminal postscript enhanced color solid lw l_w "Times-Roman" f_s; \
  ext = "ps"; \
else if (ps == 2) \
  set terminal postscript eps enhanced color solid lw l_w "Times-Roman" f_s; \
  ext = "eps"; \
else if (ps == 3) \
  set terminal pdf enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "pdf"; \
else if (ps == 4) \
  set term pngcairo enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "png";

set grid;

point_size = 0.7;
set style line 1 lt 1 lw 1 lc rgb "blue"  pointsize point_size;
set style line 2 lt 1 lw 1 lc rgb "green" pointsize point_size;
set style line 3 lt 1 lw 1 lc rgb "red"   pointsize point_size;

if (ps) set output "track.".(ext)

set multiplot;

set size 0.5, 0.5; set origin 0.0, 0.5;
if (!J_phi) \
  set title "Horizontal Phase Space"; \
  set xlabel "x [mm]"; set ylabel "p_x [mrad]"; \
  plot "track.out" using 2:3 notitle with points ls 1; \
else \
  set title "Horizontal Action-Angle Variables"; \
  set xlabel "{/Symbol f}_x [rad]"; set ylabel "J_x [mm.mrad]"; \
  set yrange [0:]; \
  plot "dJ.out" using 3:2 notitle with points ls 1;

set origin 0.5, 0.5;
if (!J_phi) \
  set title "Vertical Phase Space"; \
  set xlabel "y [mm]"; set ylabel "p_y [mrad]"; \
  plot "track.out" using 4:5 notitle with points ls 3; \
else \
  set title "Vertical Action-Angle Variables"; \
  set xlabel "{/Symbol f}_y [rad]"; set ylabel "J_y [mm.mrad]"; \
  set yrange [0:]; \
  plot "dJ.out" using 5:4 notitle with points ls 3;

set origin 0.5, 0.0;
set title "Longitudinal Phase Space";
set ylabel "{/Symbol d} [%]";
set yrange [*:*];
if (f_rf) \
  set xlabel "{/Symbol f} [{/Symbol \260}]"; \
else \
  set xlabel "ct [mm]";
plot "track.out" using 7:6 notitle with points ls 2;

unset multiplot;
if (!ps) pause -1;
