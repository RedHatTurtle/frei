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

set style line  1 lt 1        lc rgb "#000000" lw 15 pt 1 ps 1.0 # Continuous thick black line for reference solution at SPs

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
set yrange [0:1.05]
set title "Density plot - Analytical vs Numeric" font "Cantarell,60"
set output "gnuplot-density.svg"
plot "analytical/ana010.txt" using ($1+0.5):($2/10) with lines ls  1 axis x1y1 title "Analytical", \
     "sol_sp_gnuplt_sorted.dat"        using 2:3              with lines ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Momentum"
set yrange [0:1.05]
set title "Momentum plot - Analytical vs Numeric"
set output "gnuplot-momentum.svg"
plot "analytical/ana010.txt" using ($1+0.5):($2*$3/10) with lines ls  1 axis x1y1 title "Analytical", \
     "sol_sp_gnuplt_sorted.dat"        using 2:4                 with lines ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Energy"
set yrange [0:2.625]
set title "Energy plot - Analytical vs Numeric"
set output "gnuplot-energy.svg"
plot "analytical/ana010.txt" using ($1+0.5):($4/0.4+0.5*$2*$3*$3)/10 with lines ls  1 axis x1y1 title "Analytical", \
     "sol_sp_gnuplt_sorted.dat"        using 2:5                               with lines ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Velocity"
set yrange [0:1.05]
set title "Velocity plot - Analytical vs Numeric"
set output "gnuplot-velocity.svg"
plot "analytical/ana010.txt" using ($1+0.5):($3) with lines ls  1 axis x1y1 title "Analytical", \
     "sol_sp_gnuplt_sorted.dat"        using 2:($4/$3)     with lines ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Pressure"
set yrange [0:1.05]
set title "Pressure plot - Analytical vs Numeric"
set output "gnuplot-pressure.svg"
plot "analytical/ana010.txt" using ($1+0.5):($4/10) with lines ls  1 title "Analytical", \
     "sol_sp_gnuplt_sorted.dat"        using 2:6              with lines ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Temperature"
set yrange [0:300]
set title "Temperature plot - Analytical vs Numeric"
set output "gnuplot-temperature.svg"
plot "analytical/ana010.txt" using ($1+0.5):($4/(10*$2*287)) with lines ls  1 title "Analytical", \
     "sol_sp_gnuplt_sorted.dat"        using 2:($6/($3*287))          with lines ls 10 axis x1y1 title "Solution at SPs"

set ylabel "Mach"
set yrange [0:1.05]
set title "Mach plot - Analytical vs Numeric"
set output "gnuplot-mach.svg"
plot "analytical/ana010.txt" using ($1+0.5):($5) with lines ls  1 axis x1y1 title "Analytical", \
     "sol_sp_gnuplt_sorted.dat"        using 2:7           with lines ls 10 axis x1y1 title "Solution at SPs"
