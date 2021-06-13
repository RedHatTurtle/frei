#!/usr/bin/gnuplot

set key autotitle columnhead
set terminal svg size 1280,1280 enhanced background rgb "#95a3a6" font "Cantarell,20"
set grid

set xtics
set ytics
set y2tics
set ytics nomirror

set xlabel "X Coordinate"

set style line 1 lt 1 lc rgb "#000000" lw 6 pt 1 ps 1.0
set style line 2 lt 1 lc rgb "#fe4b03" lw 3 pt 1 ps 1.5
set style line 3 lt 1 lc rgb "#1d5dec" lw 3 pt 1 ps 1.0
set style line 4 lt 1 lc rgb "#069af3" lw 3 pt 1 ps 1.0
set style line 5 lt 1 lc rgb "#ffff81" lw 2 pt 1 ps 0.8

#set xrange [-5:5]
#set yrange [0:1.1]
#set y2range [0:1000]

set ylabel "Normalized Density"
set y2label "Density Residual"
set title "Density plot - Analytical vs Numeric" font "Cantarell,32"
set output "gnuplot-density.svg"
plot "sol_sp_gnuplt_0_sorted.dat"  using 2:3 with lines       ls 1 axis x1y1 title "Analytical Density", \
     "res_sp_gnuplt_1_sorted.dat"  using 2:3 with linespoints ls 5 axis x1y2 title "Sol Pts Density Residual", \
     "sol_sp_gnuplt_1_sorted.dat"  using 2:3 with linespoints ls 2 axis x1y1 title "Sol Pts Density", \
     "sol_fp1_gnuplt_1_sorted.dat" using 2:3 with      points ls 3 axis x1y1 title "Flx Pts 1 Density", \
     "sol_fp2_gnuplt_1_sorted.dat" using 2:3 with      points ls 4 axis x1y1 title "Flx Pts 2 Density"

set ylabel "Normalized Velocity"
set y2label "Velocity Residual"
set title "Velocity plot - Analytical vs Numeric"
set output "gnuplot-velocity.svg"
plot "sol_sp_gnuplt_0_sorted.dat"  using 2:4 with lines       ls 1 axis x1y1 title "Analytical Velocity", \
     "res_sp_gnuplt_1_sorted.dat"  using 2:4 with linespoints ls 5 axis x1y2 title "Sol Pts Velocity Residual", \
     "sol_sp_gnuplt_1_sorted.dat"  using 2:4 with linespoints ls 2 axis x1y1 title "Calculated Velocity", \
     "sol_fp1_gnuplt_1_sorted.dat" using 2:4 with      points ls 3 axis x1y1 title "Flx Pts 1 Velocity", \
     "sol_fp2_gnuplt_1_sorted.dat" using 2:4 with      points ls 4 axis x1y1 title "Flx Pts 2 Velocity"

set ylabel "Normalized Energy"
set y2label "Energy Residual"
set title "Energy plot - Analytical vs Numeric"
set output "gnuplot-energy.svg"
plot "sol_sp_gnuplt_0_sorted.dat"  using 2:5 with lines       ls 1 axis x1y1 title "Analytical Energy", \
     "res_sp_gnuplt_1_sorted.dat"  using 2:5 with linespoints ls 5 axis x1y2 title "Sol Pts Energy Residual", \
     "sol_sp_gnuplt_1_sorted.dat"  using 2:5 with linespoints ls 2 axis x1y1 title "Sol Pts Energy", \
     "sol_fp1_gnuplt_1_sorted.dat" using 2:5 with      points ls 3 axis x1y1 title "Flx Pts 1 Energy", \
     "sol_fp2_gnuplt_1_sorted.dat" using 2:5 with      points ls 4 axis x1y1 title "Flx Pts 2 Energy"

set ylabel "Normalized Pressure"
set title "Pressure plot - Analytical vs Numeric"
set output "gnuplot-pressure.svg"
plot "sol_sp_gnuplt_0_sorted.dat"  using 2:6 with lines       ls 1 axis x1y1 title "Analytical Pressure", \
     "sol_sp_gnuplt_1_sorted.dat"  using 2:6 with linespoints ls 2 axis x1y1 title "Sol Pts Pressure", \
     "sol_fp1_gnuplt_1_sorted.dat" using 2:6 with      points ls 3 axis x1y1 title "Flx Pts 1 Pressure", \
     "sol_fp2_gnuplt_1_sorted.dat" using 2:6 with      points ls 4 axis x1y1 title "Flx Pts 2 Pressure"

set ylabel "Mach"
set title "Mach plot - Analytical vs Numeric"
set output "gnuplot-mach.svg"
plot "sol_sp_gnuplt_0_sorted.dat"  using 2:7 with lines       ls 1 axis x1y1 title "Analytical Mach", \
     "sol_sp_gnuplt_1_sorted.dat"  using 2:7 with linespoints ls 2 axis x1y1 title "Sol Pts Mach", \
     "sol_fp1_gnuplt_1_sorted.dat" using 2:7 with      points ls 3 axis x1y1 title "Flx Pts 1 Mach", \
     "sol_fp2_gnuplt_1_sorted.dat" using 2:7 with      points ls 4 axis x1y1 title "Flx Pts 2 Mach"
