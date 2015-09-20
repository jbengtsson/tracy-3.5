ps = 0; eps = 1;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set palette rgbformulae 22, 13, -31;
unset colorbox;

if (ps) set output "dynap_cod.ps"
set title "Dynamic Aperture vs. Orbit in the Sextupoles";
set xlabel "x [mm]"; set ylabel "y [mm]";

plot "dynap.out" using 1:2 title "Bare Lattice" \
     with linespoints lt palette frac 0.0, \
     "10e-6/dynap.out" using 1:2 title "{/Symbol s}_{x,y}=10e-6" \
     with linespoints lt palette frac 0.2, \
     "20e-6/dynap.out" using 1:2 title "{/Symbol s}_{x,y}=20e-6" \
     with linespoints lt palette frac 0.4, \
     "50e-6/dynap.out" using 1:2 title "{/Symbol s}_{x,y}=50e-6" \
     with linespoints lt palette frac 0.8, \
     "100e-6/dynap.out" using 1:2 title "{/Symbol s}_{x,y}=100e-6" \
     with linespoints  lt palette frac 1.0

if (!ps) pause -1;
