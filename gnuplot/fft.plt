ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "fft_1.ps"
set title "FFT of 2J_x";
set xlabel "{/Symbol n}_x";
set ylabel "A_x";
plot "track_fft.fft" using 1:2 notitle with impulses ls 2;
if (!ps) pause -1;

if (ps) set output "fft_2.ps"
set title "FFT of 2J_y";
set xlabel "{/Symbol n}_y";
set ylabel "A_y";
plot "track_fft.fft" using 1:4 notitle with impulses ls 3;
if (!ps) pause -1;

exit;

if (ps) set output "fft_3.ps"
set title "FFT of  {/Symbol p}_y";
set xlabel "{/Symbol n}_x";
set ylabel "A_x";
plot "track.fft" using 1:3 notitle with impulses ls 2;
if (!ps) pause -1;

if (ps) set output "fft_3.ps"
set title "FFT of {/Symbol p}_y";
set xlabel "{/Symbol n}_y";
set ylabel "A_y";
plot "track.fft" using 1:5 notitle with impulses ls 3;
if (!ps) pause -1;
