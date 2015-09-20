ps = 0; eps = 0; scn = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;
set angles degrees;
set polar;

if (ps) set output "h_20001.ps"
set title "h_{20001}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_20001.out" using 8:7 notitle with impulses 2, \
     "h_20001_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_00201.ps"
set title "h_{00201}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_00201.out" using 8:7 notitle with impulses 2, \
     "h_00201_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_10002.ps"
set title "h_{10002}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10002.out" using 8:7 notitle with impulses 2, \
     "h_10002_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_21000.ps"
set title "h_{21000}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_21000.out" using 8:7 notitle with impulses 2, \
     "h_21000_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_30000.ps"
set title "h_{30000}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_30000.out" using 8:7 notitle with impulses 2, \
     "h_30000_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_10110.ps"
set title "h_{10110}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10110.out" using 8:7 notitle with impulses 2, \
     "h_10110_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_10020.ps"
set title "h_{10020}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10020.out" using 8:7 notitle with impulses 2, \
     "h_10020_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_10200.ps"
set title "h_{10200}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10200.out" using 8:7 notitle with impulses 2, \
     "h_10200_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;


if (ps) set output "h_12003000.ps"
set title "h_{12003000}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_12003000.out" using 13:14 notitle with impulses 2, \
     "h_12003000_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_21003000.ps"
set title "h_{21003000}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_21003000.out" using 13:14 notitle with impulses 2, \
     "h_21003000_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_30000111.ps"
set title "h_{30000111}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_30000111.out" using 13:14 notitle with impulses 2, \
     "h_30000111_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{10112100}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10112100.out" using 13:14 notitle with impulses 2, \
     "h_10112100_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{10021020}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10021020.out" using 13:14 notitle with impulses 2, \
     "h_10021020_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{10201200}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10201200.out" using 13:14 notitle with impulses 2, \
     "h_10201200_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{21000120}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_21000120.out" using 13:14 notitle with impulses 2, \
     "h_21000120_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{10110120}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10110120.out" using 13:14 notitle with impulses 2, \
     "h_10110120_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_10022100{}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10022100.out" using 13:14 notitle with impulses 2, \
     "h_10022100_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{30000102}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_30000102.out" using 13:14 notitle with impulses 2, \
     "h_30000102_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{10021011}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10021011.out" using 13:14 notitle with impulses 2, \
     "h_10021011_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{30000120}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_30000120.out" using 13:14 notitle with impulses 2, \
     "h_30000120_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{10202100}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10202100.out" using 13:14 notitle with impulses 2, \
     "h_10202100_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{10111020}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_10111020.out" using 13:14 notitle with impulses 2, \
     "h_10111020_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{01201011}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_01201011.out" using 13:14 notitle with impulses 2, \
     "h_01201011_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{01111020}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_01111020.out" using 13:14 notitle with impulses 2, \
     "h_01111020_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;

if (ps) set output "h_.ps"
set title "h_{01201020}"; set xlabel "Re(h)"; set ylabel "Im(h)";
plot "h_01201020.out" using 13:14 notitle with impulses 2, \
     "h_01201020_sum.out" using 2:1 notitle with linespoints 1;
if (!ps) pause -1;
