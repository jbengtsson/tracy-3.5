ps = 0; contour = 0;

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
  set terminal pdf enhanced color solid linewidth l_w font "Times-Roman f_s"; \
  ext = "pdf"; \
else if (ps == 4) \
  set term pngcairo enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "png";

set grid;

set style line 1 lt 1 lw l_w lc rgb "blue";
set style line 2 lt 1 lw l_w lc rgb "dark-green";
set style line 3 lt 1 lw l_w lc rgb "red";
set style line 4 lt 1 lw l_w lc rgb "dark-orange";

file_name = "`echo $TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name using ($1/255):($2/255):($3/255);

#set noztics; unset clabel;
#set view map;
# To set y-axis to left side and avoid compression of color box.
#unset pm3d;

#set cntrparam level 18;
set cntrparam level 25;

# x <-> horizontal, y <-> vertical, z <-> perpendicular to screen
# rot_x, rot_z, scale, scale_z
if (!contour) \
  set view 65, 15, 1, 1; \
else \
  set view 0, 0, 1, 1;

#set surface
#set hidden3d
set pm3d
#set contour;

#set parametric;

#set pm3d at b map;
#unset colorbox;

if (ps) set output "quad_scan_1.".ext
set view 66, 22, 1, 1;
set title "Horizontal Tune"
set xlabel "k_{QF}"; set ylabel "k_{BH}"; set zlabel "{/Symbol n}_x"
splot "quad_scan.out" using 1:2:3 notitle w lines lt palette z;
#if (!ps) pause(-1);
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "quad_scan_2.".ext
set view 36, 33, 1, 1;
set title "Vertical Tune"
set xlabel "k_{QF}"; set ylabel "k_{BH}"; set zlabel "{/Symbol n}_y"
splot "quad_scan.out" using 1:2:4 notitle w lines lt palette z;
#if (!ps) pause(-1);
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "quad_scan_3.".ext
set view 43, 36, 1, 1;
set title "Horizontal Linear Chromaticity"
set xlabel "k_{QF}"; set ylabel "k_{BH}"; set zlabel "{/Symbol x}_x"
splot "quad_scan.out" using 1:2:5 notitle w lines lt palette z;
#if (!ps) pause(-1);
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "quad_scan_4.".ext
set view 57, 48, 1, 1;
set title "Vertical Linear Chromaticity"
set xlabel "k_{QF}"; set ylabel "k_{BH}"; set zlabel "{/Symbol x}_y"
splot "quad_scan.out" using 1:2:6 notitle w lines lt palette z;
#if (!ps) pause(-1);
if (!ps) pause mouse "click on graph to cont.\n";

if (ps) set output "quad_scan_5.".ext
set view 64, 33, 1, 1;
set title "Horizontal Emittance"
set xlabel "k_{QF}"; set ylabel "k_{BH}";
set zlabel "{/Symbol e}_x\n[pm{\264}rad]"
splot "quad_scan.out" using 1:2:7 notitle w lines lt palette z;
#if (!ps) pause(-1);
if (!ps) pause mouse "click on graph to cont.\n";
