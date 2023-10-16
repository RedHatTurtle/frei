#! /bin/sh

mkdir -p ./bin

FLAG_DBG="--warnings --detailed-errors"
FLAG_OPT="--fast"
FLAG_NOBLAS="--set blasImpl=off --set lapackImpl=off"
FLAG_SYSBLAS="-lcblas -llapacke"

export CHPL_RT_NUM_THREADS_PER_LOCALE=1

clear;
echo
echo "Compiling benchmarks"
module purge;
module load chapel-1.32.0-sl;
echo "  - Opt, single-locale, no  BLAS"
chpl -o bin/inverse-bench-sl-opt-noblas.chapel $FLAG_OPT $FLAG_NOBLAS  inverse-bench.chpl inverse.chpl
echo
echo "  - Opt, single-locale, Sys BLAS"
chpl -o bin/inverse-bench-sl-opt-blas.chapel   $FLAG_OPT $FLAG_SYSBLAS inverse-bench.chpl inverse.chpl
echo
module purge;
module load chapel-1.32.0-ml;
echo "  - Opt,  multi-locale, no  BLAS"
chpl -o bin/inverse-bench-ml-opt-noblas.chapel $FLAG_OPT $FLAG_NOBLAS  inverse-bench.chpl inverse.chpl
echo
echo "  - Opt,  multi-locale, Sys BLAS"
chpl -o bin/inverse-bench-ml-opt-blas.chapel   $FLAG_OPT $FLAG_SYSBLAS inverse-bench.chpl inverse.chpl
echo
echo
echo "############################################################"
echo
echo "Running Benchmarks"
echo "############################################################"
echo "  - Opt, single-locale, no BLAS"
./bin/inverse-bench-sl-opt-noblas.chapel      -nl 1 --n=64 --testTime=10
echo
echo "############################################################"
echo "  - Opt, single-locale, Sys BLAS"
./bin/inverse-bench-sl-opt-blas.chapel -nl 1 --n=64 --testTime=10
echo
echo "############################################################"
echo "  - Opt, multi-locale,  no BLAS"
./bin/inverse-bench-ml-opt-noblas.chapel      -nl 1 --n=64 --testTime=10
echo
echo "############################################################"
echo "  - Opt, multi-locale,  Sys BLAS"
./bin/inverse-bench-ml-opt-blas.chapel -nl 1 --n=64 --testTime=10
echo
echo "############################################################"
