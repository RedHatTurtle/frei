#!/bin/sh

sort -gk2 sol_sp_gnuplt-0.dat                > sol_sp_gnuplt_initial.dat
sort -gk2 sol_sp_gnuplt-$1.dat               > sol_sp_gnuplt_sorted.dat
sort -gk2 sol_fp1_gnuplt-$1.dat              > sol_fp1_gnuplt_sorted.dat
sort -gk2 sol_fp2_gnuplt-$1.dat              > sol_fp2_gnuplt_sorted.dat
sort -gk2 flx_fp1_gnuplt-iter_$1-stage_1.dat > flx_fp1_gnuplt_sorted.dat
sort -gk2 flx_fp2_gnuplt-iter_$1-stage_1.dat > flx_fp2_gnuplt_sorted.dat
sort -gk2 res_src_gnuplt-iter_$1-stage_1.dat > res_src_gnuplt_sorted.dat
sort -gk2 res_dsc_gnuplt-iter_$1-stage_1.dat > res_dsc_gnuplt_sorted.dat
sort -gk2 res_cnt_gnuplt-iter_$1-stage_1.dat > res_cnt_gnuplt_sorted.dat
sort -gk2 res_sp_gnuplt-$1.dat               > res_sp_gnuplt_sorted.dat
sort -gk2 res_tot_gnuplt-iter_$1-stage_1.dat > res_sp_gnuplt_sorted.dat

gnuplot -c shocktube-1d.gnu
