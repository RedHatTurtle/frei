#!/bin/sh

sort -gk2 sol_sp_gnuplt-0.dat                > sol_sp_gnuplt_init.dat
sort -gk2 sol_sp_gnuplt-$1.dat               > sol_sp_gnuplt_sorted.dat

gnuplot -c shocktube-1d.gnu
