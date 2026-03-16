#!/bin/sh

prm1=${1:-""}
prm2=${2:-0}

gnuplot << EOP

home_dir = "$prm1"
ps       = $prm2

file_name_1 = home_dir."tbt.out"
file_name_2 = home_dir."sls.out"
file_name_3 = home_dir."sls_fft.out"
file_name_4 = home_dir."tbt_phase_space.out"

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

set grid;

set style line 1 lt 1 lw 1 pt 1 lc rgb "blue";
set style line 2 lt 1 lw 1 pt 1 lc rgb "dark-green";
set style line 3 lt 1 lw 1 pt 1 lc rgb "red";
set style line 4 lt 1 lw 1 pt 1 lc rgb "dark-orange";


if (ps) set output "fft_sls_1.".ext;

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol b}_x";
set xlabel "s [m]"; set ylabel "[m]";
plot file_name_1 using 2:3:5 notitle with errorbars ls 1, \
     file_name_1 using 2:15 notitle with lines ls 1;

set origin 0.0, 0.0;
set title "{/Symbol b}_y";
set xlabel "s [m]"; set ylabel "[m]";
plot file_name_1 using 2:6:8 notitle with errorbars ls 3, \
     file_name_1 using 2:16 notitle with lines ls 3;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output "fft_sls_2.".ext;

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol Dm}_x/2{/Symbol p}";
set xlabel "s [m]"; set ylabel "";
plot file_name_1 using 2:9:11 notitle with errorbars ls 1;

set origin 0.0, 0.0;
set title "{/Symbol Dm}_y/2{/Symbol p}";
set xlabel "s [m]"; set ylabel "";
plot file_name_1 using 2:12:14 notitle with errorbars ls 3;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output "fft_sls_3.".ext;

set multiplot;

set size 0.5, 0.5; set origin 0.0, 0.5;
set title "Horizontal Position";
set xlabel "n"; set ylabel "[mm]";
plot file_name_2 using 1:(1e3*\$2) notitle with impulses ls 1;

set origin 0.5, 0.5;
set title "Vertical Position";
set xlabel "n"; set ylabel "[mm]";
plot file_name_2 using 1:(1e3*\$3) notitle with impulses ls 3;

set origin 0.0, 0.0;
set title "FFT of Horizontal Position";
set xlabel "{/Symbol n}"; set ylabel "";
plot file_name_3 using 2:(1e3*\$3) notitle with impulses ls 1;

set origin 0.5, 0.0;
set title "FFT of Vertical Position";
set xlabel "{/Symbol n}"; set ylabel "";
plot file_name_3 using 2:(1e3*\$4) notitle with impulses ls 3;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output "fft_sls_4.".ext;

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Hor Floquet Space";
set xlabel "x"; set ylabel "p_x";
plot file_name_4 using (1e3*\$2):(1e3*\$3) notitle w points lt 3 \
     pointsize 0.5;

set origin 0.0, 0.0;
set title "Ver Floquet Space";
set xlabel "y"; set ylabel "p_y";
plot file_name_4 using (1e3*\$4):(1e3*\$5) notitle w points lt 1 \
     pointsize 0.5;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output "fft_sls_5.".ext;

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "2J_x";
set xlabel "{Symbol f}_x"; set ylabel "";
set logscale y
set format y "10^{%L}"
plot file_name_4 using 1:6 notitle w points lt 3;

set origin 0.0, 0.0;
set title "2J_y";
set xlabel "{Symbol f}_y"; set ylabel "";
plot file_name_4 using 1:7 notitle w points lt 1;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n"

EOP
