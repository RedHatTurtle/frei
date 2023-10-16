#! /bin/sh

mkdir -p ./bin

FLAG_DBG="--warnings --detailed-errors"
FLAG_OPT="--fast"
FLAG_NOBLAS="--set blasImpl=off --set lapackImpl=off"
FLAG_SYSBLAS="-lcblas -llapacke"

export CHPL_RT_NUM_THREADS_PER_LOCALE=2

clear;
module purge;
module load chapel-1.32.0-sl;
echo
echo "Compiling tests"
echo "  - Dbg, single-locale, no  BLAS"
chpl -o bin/inverse-test-sl-dbg-noblas.chapel $FLAG_DBG $FLAG_NOBLAS  inverse-test.chpl  inverse.chpl
echo
echo "  - Dbg, single-locale, Sys BLAS"
chpl -o bin/inverse-test-sl-dbg-blas.chapel   $FLAG_DBG $FLAG_SYSBLAS inverse-test.chpl  inverse.chpl
echo
module purge;
module load chapel-1.32.0-ml;
echo "  - Dbg,  multi-locale, no  BLAS"
chpl -o bin/inverse-test-ml-dbg-noblas.chapel $FLAG_DBG $FLAG_NOBLAS  inverse-test.chpl  inverse.chpl
echo
echo "  - Dbg,  multi-locale, Sys BLAS"
chpl -o bin/inverse-test-ml-dbg-blas.chapel  $FLAG_DBG $FLAG_SYSBLAS inverse-test.chpl  inverse.chpl
echo
echo
echo "############################################################"
echo
echo "Running Tests"
echo "############################################################"
echo "  - Dbg, single-locale, no BLAS";
./bin/inverse-test-sl-dbg-noblas.chapel -nl 1 --n=3
echo
echo "############################################################"
echo "  - Dbg, single-locale, Sys BLAS";
./bin/inverse-test-sl-dbg-blas.chapel   -nl 1 --n=3
echo
echo "############################################################"
echo "  - Dbg, multi-locale,  no BLAS";
./bin/inverse-test-ml-dbg-noblas.chapel -nl 1 --n=3
echo
echo "############################################################"
echo "  - Dbg, multi-locale,  Sys BLAS";
./bin/inverse-test-ml-dbg-blas.chapel   -nl 1 --n=3
echo
echo "############################################################"
