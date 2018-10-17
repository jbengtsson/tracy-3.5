ps = 1; eps = 1; stat = 1; gap = 0; dw = 1;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
#if (ps && !eps) set terminal postscript enhanced color solid;
#if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "MomAper_mult.ps"
set title "Momentum Aperture";
set xlabel "s [m]";
set ylabel "{/Symbol d}";
#set yrange [-.035:.035];
if(gap) set yrange [-.05:.05];
if(dw) set yrange [-.07:.07];

if(gap)	plot "linlat.out" using 3:(0.003*$4) notitle with fsteps lt 1 lw 1 \
           lc rgb "black", \
          "0.1mm/MomAper.out" using 2:3 title "0.1 mm" with lines lc rgb "red", \
          "0.1mm/MomAper.out" using 2:4 notitle with lines lc rgb "red", \
	   "1mm/MomAper.out" using 2:3 title "1mm" with lines lc rgb "orange", \
           "1mm/MomAper.out" using 2:4 notitle with lines lc rgb "orange", \
 	   "2mm/MomAper.out" using 2:3 title "2mm" with lines lc rgb "yellow", \
           "2mm/MomAper.out" using 2:4 notitle with lines lc rgb "yellow", \
           "2.5mm/MomAper.out" using 2:3 title "2.5mm" with lines lc rgb "green", \
           "2.5mm/MomAper.out" using 2:4 notitle with lines lc rgb "green", \
           "5mm/MomAper.out" using 2:3 title "5mm" with lines lc rgb "blue", \
           "5mm/MomAper.out" using 2:4 notitle with lines lc rgb "blue", \
           "10mm/MomAper.out" using 2:3 title "10mm" with lines lc rgb "purple", \
           "10mm/MomAper.out" using 2:4 notitle with lines lc rgb "purple";


if(dw) plot "linlat.out" using 3:(0.003*$4) notitle with fsteps lt 1 lw 1 \
           lc rgb "black", \
          "dw1/MomAper.out" using 2:3 title "1 dw" with lines lc rgb "red", \
          "dw1/MomAper.out" using 2:4 notitle with lines lc rgb "red", \
          "dw2/MomAper.out" using 2:3 title "2 dw" with lines lc rgb "blue", \
          "dw2/MomAper.out" using 2:4 notitle with lines lc rgb "blue", \
	  "dw5/MomAper.out" using 2:3 title "5 dw" with lines lc rgb "green", \
          "dw5/MomAper.out" using 2:4 notitle with lines lc rgb "green";

if (!ps) pause -1;
