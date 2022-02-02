#!/bin/bash

FREI_DIR=$(pwd)
GIT_HASH=$(git rev-parse HEAD)
GIT_HASH_SHORT=$(git rev-parse --short HEAD)
COMPILE_TIME=$(date -u +'%Y-%m-%d %H:%M:%S UTC')
GIT_BRANCH=$(git branch | grep "^\*" | sed 's/^..//')
SOURCES_HASH=$(sha1sum ./src/*.chpl | sha1sum | head -c 40)
BUILD_SCRIPT_HASH=$(sha1sum ./build_frei.sh | head -c 40)
CHANGES=``

git diff --quiet HEAD ./src; SOURCE_CHANGES=$?
git diff --quiet HEAD ./build_frei.sh; BUILD_SCRIPT_CHANGES=$?
if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
    echo "Compiling FREI"
    echo "   - Commit:" $GIT_HASH
    echo "   - Branch:" $GIT_BRANCH
    echo "   - Time  :" $COMPILE_TIME
elif [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -ne 0 ] ; then
    CHANGES=$CHANGES".bld-${BUILD_SCRIPT_HASH::6}"
    echo "Compiling FREI with modified build script"
    echo "   - Commit:" $GIT_HASH
    echo "   - Branch:" $GIT_BRANCH
    echo "   - Time  :" $COMPILE_TIME
    echo "Build script hash:" $BUILD_SCRIPT_HASH
else
    CHANGES=$CHANGES".src-${SOURCES_HASH::6}"
    if [ $BUILD_SCRIPT_CHANGES -ne 0 ] ; then
        CHANGES=$CHANGES".bld-${BUILD_SCRIPT_HASH::6}"
    fi
    echo "Compiling FREI with modified sources based on"
    echo "   - Commit:" $GIT_HASH
    echo "   - Branch:" $GIT_BRANCH
    echo "   - Time  :" $COMPILE_TIME
    echo "Build script hash:" $BUILD_SCRIPT_HASH
    echo "Source files hash:" $SOURCES_HASH
fi

EXTENSION="$GIT_HASH_SHORT$CHANGES"

if [ ! -d "build" ]; then
    echo "Creating build directory"
    mkdir build
fi
echo "------------------------------------------------------------"

##################################################
###   Debug Build                              ###
##################################################
PATH_TO_CBLAS_DIR="/usr/include"
PATH_TO_BLAS_LIBS="/usr/lib64"
PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
PATH_TO_LIBGFORTRAN="/usr/lib64"
PATH_TO_LAPACK_BINARIES="/usr/lib64"

echo "(1/3) Building Debug version of Frei..."
echo
chpl -o build/frei_dbg.$EXTENSION                          \
     --warnings                                            \
     --warn-unstable                                       \
     -I$PATH_TO_CBLAS_DIR                                  \
     -L$PATH_TO_BLAS_LIBS -lcblas                          \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
     -L$PATH_TO_LIBGFORTRAN -lgfortran                     \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
     --main-module "FREI" src/FREI.chpl                    \
                          src/fr.chpl                      \
                          src/temporal.chpl                \
                          src/output.chpl                  \
                          src/boundary.chpl                \
                          src/frmesh.chpl                  \
                          src/mapping.chpl                 \
                          src/mesh.chpl                    \
                          src/gmesh.chpl                   \
                          src/riemann.chpl                 \
                          src/flux.chpl                    \
                          src/correction.chpl              \
                          src/limiter.chpl                 \
                          src/projection.chpl              \
                          src/quadrature.chpl              \
                          src/interpolation.chpl           \
                          src/polynomials.chpl             \
                          src/sourceterm.chpl              \
                          src/init.chpl                    \
                          src/ringleb.chpl                 \
                          src/config.chpl                  \
                          src/input.chpl                   \
                          src/parameters.chpl              \
                          src/testing.chpl                 \
    2>&1 | tee build/frei_dbg.$EXTENSION.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "OK"
    ln -srf build/frei_dbg.$EXTENSION ./build/frei_dbg
    if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
        ln -srf build/frei_dbg.$EXTENSION ./build/frei_dbg.$GIT_BRANCH
    else
        ln -srf build/frei_dbg.$EXTENSION ./build/frei_dbg.$GIT_BRANCH-alpha
    fi
else
    echo "FAILED"
fi
echo "------------------------------------------------------------"

##################################################
###   Optimized / Production Build             ###
##################################################
echo "(2/3) Building Optimized version of Frei..."
echo
chpl -o build/frei_opt.$EXTENSION                          \
     --fast                                                \
     -I$PATH_TO_CBLAS_DIR                                  \
     -L$PATH_TO_BLAS_LIBS -lcblas                          \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
     -L$PATH_TO_LIBGFORTRAN -lgfortran                     \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
     --main-module "FREI" src/FREI.chpl                    \
                          src/fr.chpl                      \
                          src/temporal.chpl                \
                          src/output.chpl                  \
                          src/boundary.chpl                \
                          src/frmesh.chpl                  \
                          src/mapping.chpl                 \
                          src/mesh.chpl                    \
                          src/gmesh.chpl                   \
                          src/riemann.chpl                 \
                          src/flux.chpl                    \
                          src/correction.chpl              \
                          src/limiter.chpl                 \
                          src/projection.chpl              \
                          src/quadrature.chpl              \
                          src/interpolation.chpl           \
                          src/polynomials.chpl             \
                          src/sourceterm.chpl              \
                          src/init.chpl                    \
                          src/ringleb.chpl                 \
                          src/config.chpl                  \
                          src/input.chpl                   \
                          src/parameters.chpl              \
                          src/testing.chpl                 \
    2>&1 | tee build/frei_opt.$EXTENSION.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "OK"
    ln -srf build/frei_opt.$EXTENSION ./build/frei_opt
    if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
        ln -srf build/frei_opt.$EXTENSION ./build/frei_opt.$GIT_BRANCH
    else
        ln -srf build/frei_opt.$EXTENSION ./build/frei_opt.$GIT_BRANCH-alpha
    fi
else
    echo "FAILED"
fi
echo "------------------------------------------------------------"

##################################################
###   Optimizes Intel MKL Build                ###
##################################################
PATH_TO_CBLAS_DIR="/opt/intel/oneapi/mkl/latest/include"
PATH_TO_BLAS_LIBS="/opt/intel/oneapi/mkl/latest/lib/intel64"
PATH_TO_LAPACKE_INCLUDE_DIR="/opt/intel/oneapi/mkl/latest/include"
PATH_TO_LIBGFORTRAN="/opt/intel/oneapi/mkl/latest/lib/intel64"
PATH_TO_LAPACK_BINARIES="/opt/intel/oneapi/mkl/latest/lib/intel64"

echo "(3/3) Building Optimized Intel MKL version of Frei..."
echo
chpl -o build/frei_mkl.$EXTENSION                          \
     --fast                                                \
     --set blasImpl=mkl                                    \
     --set lapackImpl=mkl                                  \
     -I$PATH_TO_CBLAS_DIR                                  \
     -L$PATH_TO_BLAS_LIBS -lblas                           \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
     -L$PATH_TO_LIBGFORTRAN -lgfortran                     \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
     --main-module "FREI" src/FREI.chpl                    \
                          src/fr.chpl                      \
                          src/temporal.chpl                \
                          src/output.chpl                  \
                          src/boundary.chpl                \
                          src/frmesh.chpl                  \
                          src/mapping.chpl                 \
                          src/mesh.chpl                    \
                          src/gmesh.chpl                   \
                          src/riemann.chpl                 \
                          src/flux.chpl                    \
                          src/correction.chpl              \
                          src/limiter.chpl                 \
                          src/projection.chpl              \
                          src/quadrature.chpl              \
                          src/interpolation.chpl           \
                          src/polynomials.chpl             \
                          src/sourceterm.chpl              \
                          src/init.chpl                    \
                          src/ringleb.chpl                 \
                          src/config.chpl                  \
                          src/input.chpl                   \
                          src/parameters.chpl              \
                          src/testing.chpl                 \
    2>&1 | tee build/frei_mkl.$EXTENSION.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "OK"
    ln -srf build/frei_mkl.$EXTENSION ./build/frei_mkl
    if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
        ln -srf build/frei_mkl.$EXTENSION ./build/frei_mkl.$GIT_BRANCH
    else
        ln -srf build/frei_mkl.$EXTENSION ./build/frei_mkl.$GIT_BRANCH-alpha
    fi
else
    echo "FAILED"
fi
echo "------------------------------------------------------------"
