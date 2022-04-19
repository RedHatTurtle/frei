#!/bin/bash

# System BLAS/LAPACK
PATH_TO_CBLAS_DIR="/usr/include"
PATH_TO_BLAS_LIBS="/usr/lib64"
PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
PATH_TO_LIBGFORTRAN="/usr/lib64"
PATH_TO_LAPACK_BINARIES="/usr/lib64"

# AMD Libraries
#PATH_TO_CBLAS_DIR="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/include"
#PATH_TO_BLAS_LIBS="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/lib64"
#PATH_TO_LAPACKE_INCLUDE_DIR="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/include"
#PATH_TO_LIBGFORTRAN="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/lib"
#PATH_TO_LAPACK_BINARIES="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/lib"

# Intel Libraries
#PATH_TO_CBLAS_DIR="/opt/intel/oneapi/mkl/latest/include"
#PATH_TO_BLAS_LIBS="/opt/intel/oneapi/mkl/latest/lib/intel64"
#PATH_TO_LAPACKE_INCLUDE_DIR="/opt/intel/oneapi/mkl/latest/include"
#PATH_TO_LIBGFORTRAN="/opt/intel/oneapi/mkl/latest/lib/intel64"
#PATH_TO_LAPACK_BINARIES="/opt/intel/oneapi/mkl/latest/lib/intel64"

# Build all tests and save log
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Interpolation Bench:"
chpl -o bench/interpBench.chapel                                  \
     --fast                                                       \
     -L/usr/lib64/gcc/x86_64-suse-linux/11                        \
     -I$PATH_TO_CBLAS_DIR                                         \
     -L$PATH_TO_BLAS_LIBS -lcblas                                 \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR                               \
     -L$PATH_TO_LIBGFORTRAN -lgfortran                            \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas        \
     --main-module InterpBench bench/interpBench.chpl
    2>&1 | tee bench/interpBench-build.log
if [ $? -eq 0 ]; then
  echo -e "\nSuccess"
else
  echo -e "\nFailed"
fi
echo
echo
echo "------------------------------------------------------------"
echo
echo -e "Running Tests..."
echo
# Run tests and output to file
./bench/interpBench &> bench/interpBench.log
echo -e "Done"
echo
echo "------------------------------------------------------------"
echo
