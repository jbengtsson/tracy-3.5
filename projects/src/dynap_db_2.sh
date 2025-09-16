#!/bin/sh

prm1=${1:-""}
prm2=${2:-0}

gnuplot << EOP

home_dir = "$prm§"
ps       = $prm2

file_name = home_dir."linlat"

f_s = 14
l_w = 2

# Enhanced is needed for Greek characters.
if (ps == 0) {
  set terminal qt enhanced font "DejaVu Sans,12"
} else if (ps == 1) {
  set terminal postscript enhanced color solid lw l_w font "Times-Roman,".f_s
  ext = "ps"
} else if (ps == 2) {
  set terminal postscript eps enhanced color solid lw l_w font \
 "Times-Roman,".f_s
  ext = "eps"
} else if (ps == 3) {
  set terminal pdfcairo enhanced color solid lw l_w font "Times-Roman,".f_s
  ext = "pdf"
} else if (ps == 4) {
  set terminal pngcairo enhanced color solid lw l_w font "Times-Roman,".f_s
  ext = "png"
}

set grid

set style line 1 lt 1 lw 1 lc rgb "green"
set style line 2 lt 1 lw 1 lc rgb "blue"
set style line 3 lt 1 lw 1 lc rgb "cyan"
set style line 4 lt 1 lw 1 lc rgb "purple"
set style line 5 lt 1 lw 1 lc rgb "red"

if (ps) set output "dynap_db_2.ps"

set title "Dynamic Aperture vs. {/Symbol D}b_2/b_2"
set xlabel "x [mm]"
set ylabel "y [mm]"
plot "dynap_0.00e+00.out" using 1:2 title "Bare Lattice" \
     with linespoints ls 1, \
     "dynap_2.50e-04.out" using 1:2 title "{/Symbol D}b_2/b_2=2.5e-4" \
     with linespoints ls 2, \
     "dynap_5.00e-04.out" using 1:2 title "{/Symbol D}b_2/b_2=5.0e-4" \
     with linespoints ls 3, \
     "dynap_1.00e-03.out" using 1:2 title "{/Symbol D}b_2/b_2=1.0e-3" \
     with linespoints ls 4, \
     "dynap_2.50e-03.out" using 1:2 title "{/Symbol D}b_2/b_2=2.5e-3" \
     with linespoints ls 5

if (!ps) pause mouse "click on graph to cont.\n"

EOP
