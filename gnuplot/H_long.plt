ps = 0;

f_s = 14; l_w = 2;
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

set nosurface; set contour base; set noztics; set key off; unset colorbox;
# x <-> horizontal, y <-> vertical, z <-> perpendicular to screen
# rot_x, rot_z, scale, scale_z
set view 0, 0, 1, 1;
set palette rgbformulae 22, 13, -31;

if (ps) set output "H_long.".(ext);

set cntrparam level 75;
set title "Longitudinal Phase Space to O({/Symbol a}_4)"
set xlabel "{/Symbol f} [{/Symbol \260}]"; set ylabel "{/Symbol d} [%]";
splot "H_long.dat" using 1:2:3 notitle with lines lt palette z;
if (!ps) pause mouse "click on graph to cont.\n";
