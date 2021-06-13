#!/bin/bash

##################################################
###   Debug Build                              ###
##################################################

PATH_TO_CBLAS_DIR="/usr/include"
PATH_TO_BLAS_LIBS="/usr/lib64"
PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
PATH_TO_LIBGFORTRAN="/usr/lib64"
PATH_TO_LAPACK_BINARIES="/usr/lib64"

if test -f "frei_dbg-build.log" ; then
    mv frei_dbg-build.log frei_dbg-build.log.old
fi

#    --set lapackImpl=off                        \
#    --set blasImpl=off                          \

chpl -o frei_dbg                                 \
     --warnings                                  \
     --warn-unstable                             \
     -I/usr/include                              \
     -L/usr/lib64 -lcblas                        \
     -I/usr/include                              \
     -L/usr/lib64 -lgfortran                     \
     -L/usr/lib64 -llapacke -llapack -lcblas     \
     --main-module "FREI" src/FREI.chpl          \
                          src/fr.chpl            \
                          src/temporal.chpl      \
                          src/output.chpl        \
                          src/boundary.chpl      \
                          src/frmesh.chpl        \
                          src/mesh.chpl          \
                          src/gmesh.chpl         \
                          src/riemann.chpl       \
                          src/flux.chpl          \
                          src/correction.chpl    \
                          src/interpolation.chpl \
                          src/polynomials.chpl   \
                          src/sourceterm.chpl    \
                          src/init.chpl          \
                          src/config.chpl        \
                          src/input.chpl         \
                          src/parameters.chpl    \
                          src/testing.chpl       |
    tee frei_dbg-build.log

##################################################
###   Optimized / Production Build             ###
##################################################

PATH_TO_CBLAS_DIR="/usr/include"
PATH_TO_BLAS_LIBS="/usr/lib64"
PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
PATH_TO_LIBGFORTRAN="/usr/lib64"
PATH_TO_LAPACK_BINARIES="/usr/lib64"

if test -f "frei_opt-build.log" ; then
    mv frei_opt-build.log frei_opt-build.log.old
fi
chpl -o frei_opt                                 \
     --fast                                      \
     -I$PATH_TO_CBLAS_DIR                        \
     -L$PATH_TO_BLAS_LIBS -lcblas                \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR              \
     -L$PATH_TO_LIBGFORTRAN -lgfortran           \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas     \
     --main-module "FREI" src/FREI.chpl          \
                          src/fr.chpl            \
                          src/temporal.chpl      \
                          src/output.chpl        \
                          src/boundary.chpl      \
                          src/frmesh.chpl        \
                          src/mesh.chpl          \
                          src/gmesh.chpl         \
                          src/riemann.chpl       \
                          src/flux.chpl          \
                          src/correction.chpl    \
                          src/interpolation.chpl \
                          src/polynomials.chpl   \
                          src/sourceterm.chpl    \
                          src/init.chpl          \
                          src/config.chpl        \
                          src/input.chpl         \
                          src/parameters.chpl    \
                          src/testing.chpl       |
    tee frei_opt-build.log

##################################################
###   Optimized / Production Build             ###
##################################################

PATH_TO_CBLAS_DIR="/opt/intel/oneapi/mkl/latest/include"
PATH_TO_BLAS_LIBS="/opt/intel/oneapi/mkl/latest/lib/intel64"
PATH_TO_LAPACKE_INCLUDE_DIR="/opt/intel/oneapi/mkl/latest/include"
PATH_TO_LIBGFORTRAN="/opt/intel/oneapi/mkl/latest/lib/intel64"
PATH_TO_LAPACK_BINARIES="/opt/intel/oneapi/mkl/latest/lib/intel64"

#if test -f "frei_opt_mkl-build.log" ; then
#    mv frei_opt_mkl-build.log frei_opt_mkl-build.log.old
#fi
#chpl -o frei_opt_mkl                             \
#     --fast                                      \
#     -I$PATH_TO_CBLAS_DIR                        \
#     -L$PATH_TO_BLAS_LIBS -lblas                 \
#     -I$PATH_TO_LAPACKE_INCLUDE_DIR              \
#     -L$PATH_TO_LIBGFORTRAN -lgfortran           \
#     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lrefblas \
#     --set blasImpl=mkl                          \
#     --set lapackImpl=mkl                        \
#     --main-module "FREI" src/FREI.chpl          \
#                          src/fr.chpl            \
#                          src/temporal.chpl      \
#                          src/output.chpl        \
#                          src/boundary.chpl      \
#                          src/frmesh.chpl        \
#                          src/mesh.chpl          \
#                          src/gmesh.chpl         \
#                          src/riemann.chpl       \
#                          src/flux.chpl          \
#                          src/correction.chpl    \
#                          src/interpolation.chpl \
#                          src/polynomials.chpl   \
#                          src/sourceterm.chpl    \
#                          src/init.chpl          \
#                          src/config.chpl        \
#                          src/input.chpl         \
#                          src/parameters.chpl    \
#                          src/testing.chpl       |
#    tee frei_opt-build.log
