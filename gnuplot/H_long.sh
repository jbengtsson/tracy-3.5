#!/bin/sh

prm1=${1-0}

gnuplot << EOP

ps = $prm1;

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

set nosurface; set contour base; set noztics; set key off; unset colorbox
# x <-> horizontal, y <-> vertical, z <-> perpendicular to screen
# rot_x, rot_z, scale, scale_z

set view map
set palette rgbformulae 22, 13, -31
set cntrparam level 50

if (ps) set output "H_long.".(ext)

set title "Longitudinal Phase Space to O({/Symbol a}_4)"
# Greek letters doesn't work for terminal.
set xlabel "phase [deg]"
# Degree symbol does not work for qt terminal.
# set xlabel "{/Symbol f} [{/Symbol \260}]"
set ylabel "{/Symbol d} [%]"
splot "H_long.dat" using 1:2:3 notitle with lines lt palette z
if (!ps) pause mouse "click on graph to cont.\n"

EOP
