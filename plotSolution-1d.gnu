#!/bin/gnuplot

set key autotitle columnhead
set terminal svg size 3840,2160 enhanced background rgb "#D0E0D0" font "Cantarell,42"
set size 1.02,1.02

set colorsequence default

#set key outside
set key bottom right
set key font "Cantarell,36"
set key autotitle columnhead

set grid xtics mxtics ytics mytics

set xtics
set mxtics 5
set ytics
set mytics 5
set ytics nomirror
set y2tics
set my2tics 5

set xlabel "X Coordinate"

set style line  1 lt 1        lc rgb "#000000" lw 15 pt 1 ps 1.0 # Continuous thick black line for initial solution at SPs

set style line 10 lt 1        lc rgb "#9D0208" lw  7 pt 1 ps 1.5 # Continuous medium orange line for current iteration solution at SPs
set style line  6 lt 1        lc rgb "#8f1402" lw  7 pt 1 ps 1.5 # Dark red continuous line for solution at SPs in next iteration

set style line  3 lt 1        lc rgb "#069af3" lw  7 pt 1 ps 1.0 # Light Blue points for interpolated solution at face FP 1
set style line  4 lt 1        lc rgb "#1d5dec" lw  7 pt 1 ps 1.0 # Bark Blue points for interpolated solution at face FP 2

set style line 30 lt 1        lc rgb "#DC2F02" lw  7 pt 1 ps 1.0 # Total Residue
set style line 31 lt 1 dt "_" lc rgb "#FFBA08" lw  4 pt 1 ps 0.8 # Residue Component 1 - Source Term
set style line 32 lt 1 dt "-" lc rgb "#FAA307" lw  4 pt 1 ps 0.8 # Residue Component 2 - Discontinuus Flux
set style line 33 lt 1 dt "." lc rgb "#F48C06" lw  4 pt 1 ps 0.8 # Residue Component 3 - Continuous Flux

unset y2label
set ytics mirror

set ylabel "Density"
set y2label "Density Residual"
set title "Density plot - Initial vs Current" font "Cantarell,60"
set output "gnuplot-density-debug.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:3 with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:3 with linespoints ls 10 axis x1y1 title "Solution at SP", \
     "res_sp_gnuplt_sorted.dat"  using 2:3 with linespoints ls 30 axis x1y2 title "Residue", \
     "res_src_gnuplt_sorted.dat" using 2:3 with linespoints ls 31 axis x1y2 title "Residue - Src Term", \
     "res_dsc_gnuplt_sorted.dat" using 2:3 with linespoints ls 32 axis x1y2 title "Residue - Dsc Flux", \
     "res_cnt_gnuplt_sorted.dat" using 2:3 with linespoints ls 33 axis x1y2 title "Residue - Cnt Flux", \
     "sol_fp1_gnuplt_sorted.dat" using 2:3 with      points ls  3 axis x1y1 title "Solution at FP1", \
     "sol_fp2_gnuplt_sorted.dat" using 2:3 with      points ls  4 axis x1y1 title "Solution at FP2"

set ylabel "Momentum"
set y2label "Momentum Residual"
set title "Momentum plot - Initial vs Current"
set output "gnuplot-momentum-debug.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:4 with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:4 with linespoints ls 10 axis x1y1 title "Solution at SP", \
     "res_sp_gnuplt_sorted.dat"  using 2:4 with linespoints ls 30 axis x1y2 title "Residue", \
     "res_src_gnuplt_sorted.dat" using 2:4 with linespoints ls 31 axis x1y2 title "Residue - Src Term", \
     "res_dsc_gnuplt_sorted.dat" using 2:4 with linespoints ls 32 axis x1y2 title "Residue - Dsc Flux", \
     "res_cnt_gnuplt_sorted.dat" using 2:4 with linespoints ls 33 axis x1y2 title "Residue - Cnt Flux", \
     "sol_fp1_gnuplt_sorted.dat" using 2:4 with      points ls  3 axis x1y1 title "Solution at FP1", \
     "sol_fp2_gnuplt_sorted.dat" using 2:4 with      points ls  4 axis x1y1 title "Solution at FP2"

set ylabel "Energy"
set y2label "Energy Residual"
set title "Energy plot - Initial vs Current"
set output "gnuplot-energy-debug.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:5 with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:5 with linespoints ls 10 axis x1y1 title "Solution at SP", \
     "res_sp_gnuplt_sorted.dat"  using 2:5 with linespoints ls 30 axis x1y2 title "Residue", \
     "res_src_gnuplt_sorted.dat" using 2:5 with linespoints ls 31 axis x1y2 title "Residue - Src Term", \
     "res_dsc_gnuplt_sorted.dat" using 2:5 with linespoints ls 32 axis x1y2 title "Residue - Dsc Flux", \
     "res_cnt_gnuplt_sorted.dat" using 2:5 with linespoints ls 33 axis x1y2 title "Residue - Cnt Flux", \
     "sol_fp1_gnuplt_sorted.dat" using 2:5 with      points ls  3 axis x1y1 title "Solution at FP1", \
     "sol_fp2_gnuplt_sorted.dat" using 2:5 with      points ls  4 axis x1y1 title "Solution at FP2"

unset y2label
set ytics mirror

set ylabel "Density"
set title "Density plot - Analytical vs Numeric" font "Cantarell,60"
set output "gnuplot-density.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:3             with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:3             with lines       ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Momentum"
set title "Momentum plot - Analytical vs Numeric"
set output "gnuplot-momentum.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:4             with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:4             with lines       ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Energy"
set title "Energy plot - Analytical vs Numeric"
set output "gnuplot-energy.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:5             with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:5             with lines       ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Velocity"
set title "Velocity plot - Initial vs Current"
set output "gnuplot-velocity.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:($4/$3)       with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:($4/$3)       with lines       ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Pressure"
set title "Pressure plot - Initial vs Current"
set output "gnuplot-pressure.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:6             with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:6             with linespoints ls 10 axis x1y1 title "Solution at SP"

set ylabel "Temperature"
set title "Temperature plot - Analytical vs Numeric"
set output "gnuplot-temperature.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:($6/($3*287)) with lines ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:($6/($3*287)) with lines ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Mach"
set title "Mach plot - Initial vs Current"
set output "gnuplot-mach.svg"
plot "sol_sp_gnuplt_initial.dat" using 2:7             with lines       ls  1 axis x1y1 title "Initial", \
     "sol_sp_gnuplt_sorted.dat"  using 2:7             with linespoints ls 10 axis x1y1 title "Solution at SP", \
     "sol_fp1_gnuplt_sorted.dat" using 2:7             with      points ls  3 axis x1y1 title "Solution at FP1", \
     "sol_fp2_gnuplt_sorted.dat" using 2:7             with      points ls  4 axis x1y1 title "Solution at FP2"

set ylabel "Density Flux"
set title "Density Flux" font "Cantarell,60"
set output "gnuplot-density-flux.svg"
plot "flx_fp1_gnuplt_sorted.dat" using 2:3             with linespoints ls  3 axis x1y1 title "Flux at FP1", \
     "flx_fp2_gnuplt_sorted.dat" using 2:3             with linespoints ls  4 axis x1y1 title "Flux at FP2"
