#!/bin/bash

if [ "${-#*i}" == "$-" ]; then
    Color_Off='\033[0m'       # Text Reset
    BRed='\033[1;31m'         # Bold Red
    BGreen='\033[1;32m'       # Bold Green
fi

# Initialise option flags
BUILD_GENERIC="true"
BUILD_INTEL="false"
BUILD_AMD="false"
BUILD_OPENBLAS="false"
BUILD_DBG="true"
BUILD_OPT="true"

if [[ $@ =~ --amd      ]]; then BUILD_AMD="true"; BUILD_GENERIC="false"; fi
if [[ $@ =~ --intel    ]]; then BUILD_INTEL="true"; BUILD_GENERIC="false"; fi
if [[ $@ =~ --openblas ]]; then BUILD_OPENBLAS="true"; BUILD_GENERIC="false"; fi
if [[ $@ =~ --generic  ]]; then BUILD_GENERIC="true"; fi
if [[ $@ =~ --all      ]]; then BUILD_GENERIC="true"; BUILD_OPENBLAS="true"; BUILD_AMD="true"; BUILD_INTEL="true"; fi
if [[ $@ =~ --no-dbg   ]]; then BUILD_DBG="false"; fi
if [[ $@ =~ --no-opt   ]]; then BUILD_OPT="false"; fi

FREI_DIR=$(pwd)
GIT_HASH=$(git rev-parse HEAD)
GIT_HASH_SHORT=$(git rev-parse --short HEAD)
COMPILE_TIME=$(date -u +"%Y-%m-%d %H:%M:%S UTC")
GIT_BRANCH=$(git branch | grep "^\*" | sed "s/^..//")
SOURCES_HASH=$(sha1sum ./src/*.chpl | sha1sum | head -c 40)
BUILD_SCRIPT_HASH=$(sha1sum ./build_frei.sh | head -c 40)
CHANGES=``

git diff --quiet HEAD ./src; SOURCE_CHANGES=$?
git diff --quiet HEAD ./build_frei.sh; BUILD_SCRIPT_CHANGES=$?
if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
    echo -e "Compiling FREI unaltered from"
    echo -e "   - Commit:" $GIT_HASH
    echo -e "   - Branch:" $GIT_BRANCH
    echo -e "   - Time  :" $COMPILE_TIME
elif [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -ne 0 ] ; then
    CHANGES=$CHANGES".bld-${BUILD_SCRIPT_HASH::6}"
    echo -e "Compiling FREI with modified build script"
    echo -e "   - Commit:" $GIT_HASH
    echo -e "   - Branch:" $GIT_BRANCH
    echo -e "   - Time  :" $COMPILE_TIME
    echo -e "Build script hash:" $BUILD_SCRIPT_HASH
else
    CHANGES=$CHANGES".src-${SOURCES_HASH::6}"
    if [ $BUILD_SCRIPT_CHANGES -ne 0 ] ; then
        CHANGES=$CHANGES".bld-${BUILD_SCRIPT_HASH::6}"
    fi
    echo -e "Compiling FREI with modified sources based on"
    echo -e "   - Commit:" $GIT_HASH
    echo -e "   - Branch:" $GIT_BRANCH
    echo -e "   - Time  :" $COMPILE_TIME
    echo -e "Build script hash:" $BUILD_SCRIPT_HASH
    echo -e "Source files hash:" $SOURCES_HASH
fi

EXTENSION="$GIT_HASH_SHORT$CHANGES"
GIT_BRANCH=${GIT_BRANCH//\//-}

if [ ! -d "build" ]; then
    echo -e "Creating build directory"
    mkdir build
fi
echo -e "------------------------------------------------------------"

##################################################
###   System native BLAS/Lapack based build    ###
##################################################

if [[ $BUILD_GENERIC == "true" ]]; then
    # System BLAS/LAPACK
    #PATH_TO_LIBGFORTRAN="/usr/lib64"
    PATH_TO_CBLAS_DIR="/usr/include"
    PATH_TO_BLAS_LIBS="/usr/lib64/blas"
    PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
    PATH_TO_LAPACK_BINARIES="/usr/lib64/lapack"

    if [[ $BUILD_DBG == "true" ]]; then
        ##################################################
        ###   Debug Build                              ###
        ##################################################
        echo -e "(1/2) Building Debug generic version of Frei..."
        echo
        chpl -o build/frei_dbg.$EXTENSION                          \
             --warnings                                            \
             --warn-unstable                                       \
             --detailed-errors                                     \
             -I$PATH_TO_CBLAS_DIR                                  \
             -L$PATH_TO_BLAS_LIBS -lcblas                          \
             -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
             -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_dbg.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_dbg.$EXTENSION ./build/frei_dbg
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_dbg.$EXTENSION ./build/frei_dbg.$GIT_BRANCH
            else
                ln -srf build/frei_dbg.$EXTENSION ./build/frei_dbg.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
    fi

    if [[ $BUILD_OPT == "true" ]]; then
        ##################################################
        ###   Optimized / Production Build             ###
        ##################################################
        echo -e "(2/2) Building Optimized generic version of Frei..."
        echo
        chpl -o build/frei_opt.$EXTENSION                          \
             --fast                                                \
             -I$PATH_TO_CBLAS_DIR                                  \
             -L$PATH_TO_BLAS_LIBS -lcblas                          \
             -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
             -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_opt.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_opt.$EXTENSION ./build/frei_opt
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_opt.$EXTENSION ./build/frei_opt.$GIT_BRANCH
            else
                ln -srf build/frei_opt.$EXTENSION ./build/frei_opt.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
    fi
fi


##################################################
###   OpenBLAS based build                     ###
##################################################

if [[ $BUILD_OPENBLAS == "true" ]]; then
    # OpenBLAS libraries
    #PATH_TO_LIBGFORTRAN="/usr/lib64"
    PATH_TO_CBLAS_DIR="/usr/include"
    PATH_TO_OPENBLAS_LIBS="/usr/lib64/openblas-serial/"
    PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
    PATH_TO_LAPACK_BINARIES="/usr/lib64/lapack"

    if [[ $BUILD_DBG == "true" ]]; then
        ##################################################
        ###   Debug Build                              ###
        ##################################################
        echo -e "(1/2) Building Debug OpenBLAS version of Frei..."
        echo
        chpl -o build/frei_dbg_openblas.$EXTENSION                 \
             --warnings                                            \
             --warn-unstable                                       \
             --detailed-errors                                     \
             -L$PATH_TO_OPENBLAS_LIBS -lopenblas                   \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_dbg_openblas.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_dbg_openblas.$EXTENSION ./build/frei_dbg_openblas
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_dbg_openblas.$EXTENSION ./build/frei_dbg_openblas.$GIT_BRANCH
            else
                ln -srf build/frei_dbg_openblas.$EXTENSION ./build/frei_dbg_openblas.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
    fi

    if [[ $BUILD_OPT == "true" ]]; then
        ##################################################
        ###   Optimized / Production Build             ###
        ##################################################
        echo -e "(2/2) Building Optimized OpenBLAS version of Frei..."
        echo
        chpl -o build/frei_opt_openblas.$EXTENSION                 \
             --fast                                                \
             -L$PATH_TO_OPENBLAS_LIBS -lopenblas                   \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_opt_openblas.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_opt_openblas.$EXTENSION ./build/frei_opt_openblas
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_opt_openblas.$EXTENSION ./build/frei_opt_openblas.$GIT_BRANCH
            else
                ln -srf build/frei_opt_openblas.$EXTENSION ./build/frei_opt_openblas.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
    fi
fi

##################################################
###   Intel MKL based Build                    ###
##################################################

if [[ $BUILD_INTEL == "true" ]]; then
    # Intel Libraries
    #PATH_TO_INTEL_LIBGFORTRAN="/opt/intel/oneapi/mkl/latest/lib/intel64"
    PATH_TO_MKL_INCLUDES="/opt/intel/oneapi/mkl/latest/include"
    PATH_TO_MKL_LIBS="/opt/intel/oneapi/mkl/latest/lib/intel64"

    if [[ $BUILD_DBG == "true" ]]; then
        echo -e "(1/2) Building Debug Intel MKL version of Frei..."
        echo
        chpl -o build/frei_dbg_mkl.$EXTENSION                      \
             --warnings                                            \
             --warn-unstable                                       \
             --detailed-errors                                     \
             --set blasImpl=mkl                                    \
             --set lapackImpl=mkl                                  \
             -I$PATH_TO_MKL_INCLUDES                               \
             -L$PATH_TO_MKL_LIBS -lblas -lcblas -llapack -llapacke \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_dbg_mkl.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_dbg_mkl.$EXTENSION ./build/frei_dbg_mkl
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_dbg_mkl.$EXTENSION ./build/frei_dbg_mkl.$GIT_BRANCH
            else
                ln -srf build/frei_dbg_mkl.$EXTENSION ./build/frei_dbg_mkl.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
        fi

    if [[ $BUILD_OPT == "true" ]]; then
        echo -e "(2/2) Building Optimized Intel MKL version of Frei..."
        echo
        chpl -o build/frei_opt_mkl.$EXTENSION                      \
             --fast                                                \
             --set blasImpl=mkl                                    \
             --set lapackImpl=mkl                                  \
             -I$PATH_TO_MKL_INCLUDES                               \
             -L$PATH_TO_MKL_LIBS -lblas -lcblas -llapack -llapacke \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_opt_mkl.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_opt_mkl.$EXTENSION ./build/frei_opt_mkl
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_opt_mkl.$EXTENSION ./build/frei_opt_mkl.$GIT_BRANCH
            else
                ln -srf build/frei_opt_mkl.$EXTENSION ./build/frei_opt_mkl.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
    fi
fi

##################################################
###   AMD AOCL LibM based Build                ###
##################################################

if [[ $BUILD_AMD == "true" ]]; then
    # AMD Libraries
    #PATH_TO_AMD_LIBGFORTRAN="/opt/AMD/aocl/aocl-linux-aocc-4.0.0/lib"
    PATH_TO_AMD_INCLUDES="/opt/AMD/aocl/aocl-linux-aocc-4.0.0/include"
    PATH_TO_AMD_LIBS="/opt/AMD/aocl/aocl-linux-aocc-4.0.0/lib"

    if [[ $BUILD_DBG == "true" ]]; then
        echo -e "(1/2) Building Debug AMD AOCL LibM version of Frei..."
        echo
        chpl -o build/frei_dbg_amd.$EXTENSION                      \
             --warnings                                            \
             --warn-unstable                                       \
             --detailed-errors                                     \
             -I$PATH_TO_AMD_INCLUDES                               \
             -L$PATH_TO_AMD_LIBS -lblas -lcblas -llapack -llapacke \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_dbg_amd.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_dbg_amd.$EXTENSION ./build/frei_dbg_amd
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_dbg_amd.$EXTENSION ./build/frei_dbg_amd.$GIT_BRANCH
            else
                ln -srf build/frei_dbg_amd.$EXTENSION ./build/frei_dbg_amd.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
    fi

    if [[ $BUILD_OPT == "true" ]]; then
        echo -e "(2/2) Building Optimized AMD AOCL LibM version of Frei..."
        echo
        chpl -o build/frei_opt_amd.$EXTENSION                      \
             --fast                                                \
             -I$PATH_TO_AMD_INCLUDES                               \
             -L$PATH_TO_AMD_LIBS -lblas -lcblas -llapack -llapacke \
             --main-module "FREI" src/FREI.chpl                    \
                                  src/fr.chpl                      \
                                  src/temporal.chpl                \
                                  src/output.chpl                  \
                                  src/error.chpl                   \
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
                                  src/functions/determinant.chpl   \
                                  src/functions/sorttuple.chpl     \
                                  src/parameters.chpl              \
                                  src/testing.chpl                 \
            2>&1 | tee build/frei_opt_amd.$EXTENSION.log
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo -e "${BGreen}Done${Color_Off}"
            ln -srf build/frei_opt_amd.$EXTENSION ./build/frei_opt_amd
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                ln -srf build/frei_opt_amd.$EXTENSION ./build/frei_opt_amd.$GIT_BRANCH
            else
                ln -srf build/frei_opt_amd.$EXTENSION ./build/frei_opt_amd.$GIT_BRANCH-alpha
            fi
        else
            echo -e "${BRed}Fail${Color_Off}"
        fi
        echo -e "------------------------------------------------------------"
    fi
fi
