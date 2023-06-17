ps = 1; eps = 1; rms=1; sxt=1;

 R=0.03;
# R = 0.04;
S=1.0;
if(sxt) S=R;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid lw 2;
if (ps && eps) set terminal postscript eps enhanced color solid lw 2;

set grid;

if (ps && rms) set output "all_rms.ps"
if (ps && !rms) set output "all_sys.ps"
if (rms && sxt) set title "Dynamic Aperture vs. RMS Error in Sextupoles";
if (rms && !sxt) set title "Dynamic Aperture vs. RMS Error in Quadrupoles";
if (!rms && sxt) set title "Dynamic Aperture vs.Systematic Error in Sextupoles";
if (!rms && !sxt) set title "Dynamic Aperture vs.Systematic Error in Quadrupoles";
set label "R=%2.2g m",R at graph .9,.6 center;
set logscale x;
if (sxt) set xtics (1e-7, 1e-6,1e-5,1e-4,1e-3,2e-3,1e-2,1e-1)
if (!sxt) set xtics (1e-6, 1e-5,1e-4,4e-4,1e-3,1e-2,1e-1)
if (sxt) set xlabel "{/Symbol D}B^{(3)}_n";
if (!sxt) set xlabel "{/Symbol D}B^{(2)}_n";
set ylabel "DA [{mm}^2]";

if (!sxt && !rms) \
plot    "../out/N2Db3_sys" using ($4*R/S):7 title "B_3" with linespoints 3, \
	"../out/N2Da3_sys" using ($4*R/S):7 title "A_3" with linespoints 4, \
	"../out/N2Db4_sys" using ($4*R*R/S):7 title "B_4" with linespoints 5, \
	"../out/N2Da4_sys" using ($4*R*R/S):7 title "{A_4" with linespoints 6, \
	"../out/N2Db5_sys" using ($4*R*R*R/S):7 title "B_5" with linespoints 7, \
	"../out/N2Da5_sys" using ($4*R*R*R/S):7 title "A_5" with linespoints 8, \
	"../out/N2Db6_sys" using ($4*R*R*R*R/S):7 title "B_6" with linespoints 9, \
	"../out/N2Da6_sys" using ($4*R*R*R*R/S):7 title "A_6" with linespoints 10;
;
if (sxt && !rms) \
plot    "../out/N3Db3_sys" using ($4*R/S):7 title "B_3" with linespoints 3, \
	"../out/N3Da3_sys" using ($4*R/S):7 title "A_3" with linespoints 4, \
	"../out/N3Db4_sys" using ($4*R*R/S):7 title "B_4" with linespoints 5, \
	"../out/N3Da4_sys" using ($4*R*R/S):7 title "{A_4" with linespoints 6, \
	"../out/N3Db5_sys" using ($4*R*R*R/S):7 title "B_5" with linespoints 7, \
	"../out/N3Da5_sys" using ($4*R*R*R/S):7 title "A_5" with linespoints 8, \
	"../out/N3Db6_sys" using ($4*R*R*R*R/S):7 title "B_6" with linespoints 9, \
	"../out/N3Da6_sys" using ($4*R*R*R*R/S):7 title "A_6" with linespoints 10;

if (!sxt && rms) \
plot    "../out/N2Db3_rms" using ($4*R/S):7 title "B_3" with linespoints 3, \
	"../out/N2Da3_rms" using ($4*R/S):7 title "A_3" with linespoints 4, \
	"../out/N2Db4_rms" using ($4*R*R/S):7 title "B_4" with linespoints 5, \
	"../out/N2Da4_rms" using ($4*R*R/S):7 title "{A_4" with linespoints 6, \
	"../out/N2Db5_rms" using ($4*R*R*R/S):7 title "B_5" with linespoints 7, \
	"../out/N2Da5_rms" using ($4*R*R*R/S):7 title "A_5" with linespoints 8, \
	"../out/N2Db6_rms" using ($4*R*R*R*R/S):7 title "B_6" with linespoints 9, \
	"../out/N2Da6_rms" using ($4*R*R*R*R/S):7 title "A_6" with linespoints 10;
;
if (sxt && rms) \
plot    "../out/N3Db3_rms" using ($4*R/S):7 title "B_3" with linespoints 3, \
	"../out/N3Da3_rms" using ($4*R/S):7 title "A_3" with linespoints 4, \
	"../out/N3Db4_rms" using ($4*R*R/S):7 title "B_4" with linespoints 5, \
	"../out/N3Da4_rms" using ($4*R*R/S):7 title "{A_4" with linespoints 6, \
	"../out/N3Db5_rms" using ($4*R*R*R/S):7 title "B_5" with linespoints 7, \
	"../out/N3Da5_rms" using ($4*R*R*R/S):7 title "A_5" with linespoints 8, \
	"../out/N3Db6_rms" using ($4*R*R*R*R/S):7 title "B_6" with linespoints 9, \
	"../out/N3Da6_rms" using ($4*R*R*R*R/S):7 title "A_6" with linespoints 10;
