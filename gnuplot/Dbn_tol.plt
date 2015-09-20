# DA vs. b_n.

ps = 0; eps = 0;

home_dir = "/home/bengtsson/Results/Sep_2010/Thor/2.0_2.0/ID_tol/3/";

B_rho = 10.0;

#font_size = 30; line_width = 2;
font_size = 24; line_width = 2;
if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript portrait enhanced color solid \
  lw line_width "Times-Roman" font_size;
if (ps && eps) \
  set terminal postscript eps enhanced color solid \
  lw line_width "Times-Roman" font_size;

set grid;

set style line 1 lt 1 lw line_width lc rgb "blue";
set style line 2 lt 1 lw line_width lc rgb "dark-green";
set style line 3 lt 1 lw line_width lc rgb "red";
set style line 4 lt 1 lw line_width lc rgb "dark-orange";

set key right bottom;

if (ps) set output "Dbn_tol_1.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "DA vs. RMS Error in Long Straights (13)";
set xlabel "({/Symbol D}B_3L) [T/m]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_rms_tol_+3.dw_90" \
     using (B_rho*$2):4 title "({/Symbol D}B_3L)" with points ls 1, \
     home_dir."Dbn_rms_tol_-3.dw_90" \
     using (B_rho*$3):4 title "({/Symbol D}A_3L)" with points ls 2;

set origin 0.0, 0.0; 
set title "DA vs. RMS Error in Long Straights (13)";
set xlabel "({/Symbol D}B_4L) [T/m^{-2}]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_rms_tol_+4.dw_90" \
     using (B_rho*$2):4 title "({/Symbol D}B_4L)" with points ls 1, \
     home_dir."Dbn_rms_tol_-4.dw_90" \
     using (B_rho*$3):4 title "({/Symbol D}A_4L)" with points ls 2;

unset multiplot;
if (!ps) pause -1;

if (ps) set output "Dbn_tol_2.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "DA vs. RMS Error in Short Straights (15)";
set xlabel "({/Symbol D}B_3L) [T/m]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_rms_tol_+3.ivu_20" \
     using (B_rho*$2):4 title "({/Symbol D}B_3L)" with points ls 1, \
     home_dir."Dbn_rms_tol_-3.ivu_20" \
     using (B_rho*$3):4 title "({/Symbol D}A_3L)" with points ls 2;

set origin 0.0, 0.0; 
set title "DA vs. RMS Error in Short Straights (15)";
set xlabel "({/Symbol D}B_4L) [T/m^{-2}]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_rms_tol_+4.ivu_20" \
     using (B_rho*$2):4 title "({/Symbol D}B_4L)" with points ls 1, \
     home_dir."Dbn_rms_tol_-4.ivu_20" \
     using (B_rho*$3):4 title "({/Symbol D}A_4L)" with points ls 2;

unset multiplot;
if (!ps) pause -1;

if (ps) set output "Dbn_tol_3.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "DA vs. Sys Error in Long Straights (13)";
set xlabel "({/Symbol D}B_3L) [T/m]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_sys_tol_+3.dw_90" \
     using (B_rho*$2):4 title "({/Symbol D}B_3L)" with linespoints ls 1, \
     home_dir."Dbn_sys_tol_-3.dw_90" \
     using (B_rho*$3):4 title "({/Symbol D}A_3L)" with linespoints ls 2;

set origin 0.0, 0.0; 
set title "DA vs. Sys Error in Long Straights (13)";
set xlabel "({/Symbol D}B_4L) [T/m^{-2}]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_sys_tol_+4.dw_90" \
     using (B_rho*$2):4 title "({/Symbol D}B_4L)" with linespoints ls 1, \
     home_dir."Dbn_sys_tol_-4.dw_90" \
     using (B_rho*$3):4 title "({/Symbol D}A_4L)" with linespoints ls 2;

unset multiplot;
if (!ps) pause -1;

if (ps) set output "Dbn_tol_4.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "DA vs. Sys Error in Short Straights (15)";
set xlabel "({/Symbol D}B_3L) [T/m]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_sys_tol_+3.ivu_20" \
     using (B_rho*$2):4 title "({/Symbol D}B_3L)" with linespoints ls 1, \
     home_dir."Dbn_sys_tol_-3.ivu_20" \
     using (B_rho*$3):4 title "({/Symbol D}A_3L)" with linespoints ls 2;

set origin 0.0, 0.0; 
set title "DA vs. Sys Error in Short Straights (15)";
set xlabel "({/Symbol D}B_4L) [T/m^{-2}]"; set ylabel "DA [mm^2]";
set yrange [0:];
set logscale x;
plot home_dir."Dbn_sys_tol_+4.ivu_20" \
     using (B_rho*$2):4 title "({/Symbol D}B_4L)" with linespoints ls 1, \
     home_dir."Dbn_sys_tol_-4.ivu_20" \
     using (B_rho*$3):4 title "({/Symbol D}A_4L)" with linespoints ls 2;

unset multiplot;
if (!ps) pause -1;
