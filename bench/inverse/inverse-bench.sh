#! /bin/sh

export CHPL_RT_NUM_THREADS_PER_LOCALE=1

echo "Compiling benchmarks"
chpl -o inverse-bench-dbg               --set blasImpl=off --set lapackImpl=off  inverse-bench.chpl
chpl -o inverse-bench-dbg-lapack        -lcblas            -llapacke             inverse-bench.chpl
chpl -o inverse-bench-opt        --fast --set blasImpl=off --set lapackImpl=off  inverse-bench.chpl
chpl -o inverse-bench-opt-lapack --fast -lcblas            -llapacke             inverse-bench.chpl
echo
echo "Running Benchmarks"
./inverse-bench-dbg        --n=10 --testTime=10
./inverse-bench-dbg-lapack --n=10 --testTime=10
./inverse-bench-opt        --n=10 --testTime=10
./inverse-bench-opt-lapack --n=10 --testTime=10

