ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "track_map_1.ps"
set title "Normalized Horizontal Phase Space";
set xlabel "x";
set ylabel "p_x";
plot "track_map.dat" using 1:2:7:8 notitle with vector linetype 2;
if (!ps) pause -1;

if (ps) set output "track_map_2.ps"
set title "Normalized Vertical Phase Space";
set xlabel "y";
set ylabel "p_y";
plot "track_map.dat" using 3:4:9:10 notitle with vector linetype 2;
if (!ps) pause -1;

if (ps) set output "track_map_3.ps"
set title "Longitudinal Phase Space";
set xlabel "{/Symbol f} [{/Symbol \260}]";
set ylabel "{/Symbol d} [%]";
plot "track_map.dat" using 5:6:11:12 notitle with vector linetype 2;
if (!ps) pause -1;
