#!/bin/sh

sort -gk2 sol_sp_gnuplt-0.dat                > sol_sp_gnuplt_reference.dat
sort -gk2 sol_sp_gnuplt-$1.dat               > sol_sp_gnuplt_sorted.dat

gnuplot -c plotSolution-1d.gnu
