ps = 0; eps = 1; stat = 1;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "../out/MomAper.ps"
set title "Momentum Aperture";
set xlabel "s [m]";
set ylabel "{/Symbol d}";
#set yrange [0:];

if (stat) plot "MomAper.out" using 2:3:5 notitle with errorlines 3, \
          "MomAper.out" using 2:4:6 notitle with errorlines 3;

if(!stat) plot "MomAper.out" using 2:3 notitle with lines 3, \
          "MomAper.out" using 2:4 notitle with lines 3;

if (!ps) pause -1;
