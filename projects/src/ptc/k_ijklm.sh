#!/bin/sh

prm1=${1-0}
prm2=${2-"k_ijklm"}

gnuplot << EOP

ps        = $prm1
file_name = "$prm2"
contour   = 1

f_s = 24
l_w = 2
# Enhanced is needed for Greek characters.
if (ps == 0) \
  set terminal qt 0 enhanced font "Sans, 9" \
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

set grid

set style line 1 lt 1 lw 1 lc rgb "blue"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "red"

# file_name_jet = "`echo $TRACY_LIB`/gnuplot/jet"
# Load 64-color palette for Jet
# set palette model RGB file file_name_jet using (\$1/255):(\$2/255):(\$3/255)

set cntrparam level 25
set view map
set pm3d

set macros
JET="define (0 0 0 0.5, 1./8 0 0 1, 3./8 0 1 1, 5./8 1 1 0, 7./8 1 0 0, 1 0.5 0 0)"

MATLAB = "defined (0  0.0 0.0 0.5, \
                   1  0.0 0.0 1.0, \
                   2  0.0 0.5 1.0, \
                   3  0.0 1.0 1.0, \
                   4  0.5 1.0 0.5, \
                   5  1.0 1.0 0.0, \
                   6  1.0 0.5 0.0, \
                   7  1.0 0.0 0.0, \
                   8  0.5 0.0 0.0 )"


set palette @JET

if (ps) set output file_name."_1.".(ext)
set title "k_{22000}"
set xlabel "{/Symbol n}_x"
 set ylabel "{/Symbol n}_y"
splot file_name.".dat" using 1:2:(log(\$5)) notitle with pm3d
if (!ps) pause mouse "click on graph to cont.\n"

# set cbrange [-20:-10]

if (ps) set output file_name."_2.".(ext)
set title "k_{11110}"
set xlabel "{/Symbol n}_x"
 set ylabel "{/Symbol n}_y"
splot file_name.".dat" using 1:2:(log(\$6)) notitle with pm3d
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output file_name."_3.".(ext)
set title "k_{00220}"
set xlabel "{/Symbol n}_x"
set ylabel "{/Symbol n}_y"
splot file_name.".dat" using 1:2:(log(\$7)) notitle with pm3d
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output file_name."_4.".(ext)
set title "k_{ijklm} RMS"
set xlabel "{/Symbol n}_x"
 set ylabel "{/Symbol n}_y"
splot file_name.".dat" using 1:2:(log(\$8)) notitle with pm3d
if (!ps) pause mouse "click on graph to cont.\n"

unset cbrange

if (ps) set output file_name."_5.".(ext)
set title "k_{11002}"
set xlabel "{/Symbol n}_x"
 set ylabel "{/Symbol n}_y"
splot file_name.".dat" using 1:2:(log(\$9)) notitle with pm3d
if (!ps) pause mouse "click on graph to cont.\n"

if (ps) set output file_name."_6.".(ext)
set title "k_{00112}"
set xlabel "{/Symbol n}_x"
 set ylabel "{/Symbol n}_y"
splot file_name.".dat" using 1:2:(log(\$10)) notitle with pm3d
if (!ps) pause mouse "click on graph to cont.\n"
