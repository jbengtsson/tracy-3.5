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

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";
set style line 4 lt 1 lw 1 lc rgb "dark-orange";

if (ps) set output "ibs.".ext;

set title "Impact of IBS: {/Symbol e}_{x,y} = 16 pm.rad" \
          . ", {/Symbol s_d} = 1{/Symbol \264}10^{-3}, I_b = 5 nC";
set xlabel "Bunch Length {/Symbol s}_b [cm]";
set ylabel "{/Symbol e}_{x,y} [pm.rad]";
set y2label "{/Symbol s_d} [10^{-3}]";
set ytics nomirror; set y2tics;
plot "ibs.out" using 1:2 title "{/Symbol e}_{x,y}" with lines ls 1, \
     "ibs.out" using 1:(1e3*$7) axis x1y2 title "{/Symbol s_d}" \
     with lines ls 2;

if (!ps) pause mouse "click on graph to cont.\n";
