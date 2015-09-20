ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set nosurface; set contour base; set noztics; set key off; unset colorbox;
# x <-> horizontal, y <-> vertical, z <-> perpendicular to screen
# rot_x, rot_z, scale, scale_z
set view 0, 0, 1, 1;
set palette rgbformulae 22, 13, -31;

if (ps) set output "H_long.ps"

set cntrparam level 20;
set title "Longitudinal Phase Space to O({/Symbol a}_4)"
set xlabel "{/Symbol f} [{/Symbol \260}]"; set ylabel "{/Symbol d} [%]";
splot "H_long.dat" using 1:2:3 notitle with lines lt palette z;
if (!ps) pause(-1);
