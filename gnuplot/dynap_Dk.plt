ps = 0; eps = 1;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set palette rgbformulae 22, 13, -31;
unset colorbox;

if (ps) set output "dynap_Dk.ps"

if (1) \
  set title "Dynamic Aperture vs. {/Symbol D}b_2/b_2"; \
  set xlabel "x [mm]"; set ylabel "y [mm]"; \
  plot "dynap.out" using 1:2 title "Bare Lattice" \
       with linespoints lt palette frac 0.0, \
       "1e-4/dynap.out" using 1:2 title "{/Symbol D}b_2/b_2=1e-4" \
       with linespoints lt palette frac 0.2, \
       "2e-4/dynap.out" using 1:2 title "{/Symbol D}b_2/b_2=2e-4" \
       with linespoints lt palette frac 0.4, \
       "5e-4/dynap.out" using 1:2 title "{/Symbol D}b_2/b_2=5e-4" \
       with linespoints lt palette frac 0.7, \
       "1e-3/dynap.out" using 1:2 title "{/Symbol D}b_2/b_2=1e-3" \
       with linespoints lt palette frac 0.8, \
       "2e-3/dynap.out" using 1:2 title "{/Symbol D}b_2/b_2=2e-3" \
       with linespoints lt palette frac 1.0;

if (0) \
  set title "Dynamic Aperture vs. {/Symbol D}b_n/b_{2,3}"; \
  set xlabel "x [mm]"; set ylabel "y [mm]"; \
  plot "dynap_bare.out" using 1:2 title "Bare Lattice" with linespoints 3, \
       "dynap_mag.out" using 1:2 title "w. Multipole Errors" \
       with linespoints 1;


if (!ps) pause -1;
