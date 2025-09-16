#!/bin/sh

prm1=${1:-""}
prm2=${2:-0}

gnuplot << EOP

home_dir = "$prm1"
ps       = $prm2

# file_name = home_dir."mom_aper.out"
file_name = home_dir."touschek.out"

f_s = 14
l_w = 2

# Enhanced is needed for Greek characters.
if (ps == 0) {
  set terminal qt enhanced font "DejaVu Sans,12""
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

if (ps) set output "touschek.".(ext)
set title "Momentum Aperture"
set xlabel "s [m]"
set ylabel "{/Symbol d} [%]"
#set yrange [-5.2 : 5.2]
plot file_name using 2:3 notitle with fsteps ls 2, \
     file_name using 2:4 notitle with fsteps ls 2
if (!ps) pause mouse "click on graph to cont.\n"

EOP
