#!/bin/bash

# Define color variables for colored output if in interactive shell
if [ "${-#*i}" == "$-" ]; then
    Color_Off='\033[0m'       # Text Reset
    BRed='\033[1;31m'         # Bold Red
    BGreen='\033[1;32m'       # Bold Green
    UBlue='\033[4;34m'        # Underline Blue
fi

# Default sucessful and failed build messages
MSG_DONE="${BGreen}Done${Color_Off}"
MSG_FAIL="${BRed}Fail${Color_Off}"

BUILD_DIR="./build"
BUILD_LOG_DIR="./build/buildLogs"
if [ ! -d "build" ]; then
    echo
    echo -e "Creating build directory:      $BUILD_DIR"
    mkdir -p $BUILD_DIR
fi
if [ ! -d "build/buildLogs" ]; then
    echo -e "Creating build logs directory: $BUILD_LOG_DIR"
    mkdir -p $BUILD_LOG_DIR
fi
echo -e "------------------------------------------------------------"

FREI_DIR=$(pwd)
GIT_HASH=$(git rev-parse HEAD)
GIT_HASH_SHORT=$(git rev-parse --short HEAD)
COMPILE_TIME=$(date -u +"%Y-%m-%d %H:%M:%S UTC")
GIT_BRANCH=$(git branch | grep "^\*" | sed "s/^..//")
SOURCES_HASH=$(sha1sum ./src/*.chpl | sha1sum | head -c 40)
BUILD_SCRIPT_HASH=$(sha1sum ./build_frei.sh | head -c 40)
CHANGES=""

# Check if:
# - We're compiling an unmodified commit
# - The source code changed from the base commit
# - The build script changed from the base commit
git diff --quiet HEAD ./src; SOURCE_CHANGES=$?
git diff --quiet HEAD ./build_frei.sh; BUILD_SCRIPT_CHANGES=$?
if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
    echo
    echo -e "Compiling unmodified FREI from"
    echo -e "   - Commit:" $GIT_HASH
    echo -e "   - Branch:" $GIT_BRANCH
    echo -e "   - Time  :" $COMPILE_TIME
elif [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -ne 0 ] ; then
    CHANGES=$CHANGES".bld-${BUILD_SCRIPT_HASH::6}"
    echo
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
    echo
    echo -e "Compiling FREI with modified sources based on"
    echo -e "   - Commit:" $GIT_HASH
    echo -e "   - Branch:" $GIT_BRANCH
    echo -e "   - Time  :" $COMPILE_TIME
    echo -e "Build script hash:" $BUILD_SCRIPT_HASH
    echo -e "Source files hash:" $SOURCES_HASH
fi
echo -e "------------------------------------------------------------"

# Replace "/" with "_" on branch names so they can be user on file names
GIT_BRANCH=${GIT_BRANCH//\//_}
VERS_HASH="$GIT_HASH_SHORT$CHANGES"

##################################################
###   List of FREI Source Files                ###
##################################################

SRC_FILES="                          \
    src/FREI.chpl                    \
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
    src/dimensional.chpl             \
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
    src/testing.chpl                 "

##################################################
###   Configure Build Flags                    ###
##################################################

# Configure locale mode
FLAG_SL="--local"
FLAG_ML="--no-local"

if [[ $@ =~ --sl       ]]; then
    FLAG_LOCALE=$FLAG_SL
    BLD_LOCALE="sl";
    VER_LOCALE="${UBlue}Single Locale${Color_Off}"
elif [[ $@ =~ --ml       ]]; then
    FLAG_LOCALE=$FLAG_ML
    BLD_LOCALE="ml";
    VER_LOCALE="${UBlue}Multi Locale${Color_Off}"
else
    FLAG_LOCALE=$FLAG_SL
    BLD_LOCALE="sl";
    VER_LOCALE="${UBlue}Single Locale${Color_Off}"
fi

# Set build flag for each build type
FLAG_DBG="--warnings --warn-unstable --detailed-errors"
FLAG_DBG="--warnings --detailed-errors"
FLAG_OPT="--fast"

if [[ $@ =~ --dbg      ]]; then
    BLD_FLAG=$FLAG_DBG
    BLD_MODE="dbg"
    VER_MODE="${UBlue}Debug${Color_Off}"
elif [[ $@ =~ --opt      ]]; then
    BLD_FLAG=$FLAG_OPT
    BLD_MODE="opt"
    VER_MODE="${UBlue}Optimized${Color_Off}"
else
    BLD_FLAG=$FLAG_OPT
    BLD_MODE="opt"
    VER_MODE="${UBlue}Optimized${Color_Off}"
fi

# Configure Lapack distribution
BUILD_GENERIC="true"
BUILD_AMD="false"
BUILD_INTEL="false"
BUILD_OPENBLAS="false"

if [[ $@ =~ --amd      ]]; then BUILD_AMD="true";      BUILD_GENERIC="false"; fi
if [[ $@ =~ --intel    ]]; then BUILD_INTEL="true";    BUILD_GENERIC="false"; fi
if [[ $@ =~ --openblas ]]; then BUILD_OPENBLAS="true"; BUILD_GENERIC="false"; fi
if [[ $@ =~ --blis     ]]; then BUILD_BLIS="true";     BUILD_GENERIC="false"; fi
if [[ $@ =~ --generic  ]]; then BUILD_GENERIC="true";  fi

##################################################
###   Generate binary and links names          ###
##################################################

BUILD_NAME="frei_${BLD_LOCALE}_${BLD_MODE}.${VERS_HASH}"
LINK_NAME="frei_${BLD_LOCALE}_${BLD_MODE}"
if [[ "$GIT_BRANCH" != *"HEAD detached"* ]]; then
    LINK_NAME="${LINK_NAME}.${GIT_BRANCH}"
fi

echo
echo BUILD DIR  = $BUILD_DIR
echo BUILD NAME = $BUILD_NAME
echo LINK  NAME = $LINK_NAME

##################################################
###   System native BLAS/Lapack based build    ###
##################################################

if [[ $BUILD_GENERIC == "true" ]]; then
    # System BLAS/LAPACK
    VER_LAPACK="${UBlue}System LAPACK${Color_Off}"
    #PATH_TO_LIBGFORTRAN="/usr/lib64"
    #PATH_TO_BLAS_DIR="/usr/include"
    #PATH_TO_BLAS_LIBS="/usr/lib64/blas"
    #PATH_TO_LAPACK_INCLUDE_DIR="/usr/include"
    #PATH_TO_LAPACK_BINARIES="/usr/lib64/lapack"

    echo
    echo -e "Building $VER_LOCALE $VER_MODE $VER_LAPACK version of Frei..."
    echo
    chpl $BLD_FLAG                                          \
        -I$PATH_TO_BLAS_DIR                                 \
        -L$PATH_TO_BLAS_LIBS -lcblas                        \
        -I$PATH_TO_LAPACK_INCLUDE_DIR                       \
        -L$PATH_TO_LAPACK_BINARIES -llapacke                \
        --main-module "FREI" $SRC_FILES                     \
        -o $BUILD_DIR/$BUILD_NAME                           \
        2>&1 | tee $BUILD_LOG_DIR/$BUILD_NAME.log
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        # Short link with build version specs only
        ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/${LINK_NAME%%\.*}

        # Longer link with branch name if not a detached HEAD
        if [[ "$GIT_BRANCH" != *"HEAD detached"* ]]; then
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                # Link with branch name if code and build scripts unchanged
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME
            else
                # Link with branch name and "-alpha" for edited code or build script
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME-alpha
            fi
        fi

        echo
        echo -e $MSG_DONE
    else
        echo
        echo -e $MSG_FAIL
    fi
    echo -e "------------------------------------------------------------"
fi


##################################################
###   OpenBLAS based build                     ###
##################################################

if [[ $BUILD_OPENBLAS == "true" ]]; then
    # OpenBLAS libraries
    VER_LAPACK="${UBlue}OpenBLAS${Color_Off}"
    #PATH_TO_LIBGFORTRAN="/usr/lib64"
    #PATH_TO_OPENBLAS_DIR="/usr/include/openblas"
    #PATH_TO_OPENBLAS_LIBS="/usr/lib64/openblas-serial"

    echo
    echo -e "Building $VER_LOCALE $VER_MODE $VER_LAPACK version of Frei..."
    echo
    chpl $FLAG_LOCALE $BLD_FLAG                             \
        -I$PATH_TO_OPENBLAS_DIR                             \
        -L$PATH_TO_OPENBLAS_LIBS -lopenblas                 \
        --main-module "FREI" $SRC_FILES                     \
        -o $BUILD_DIR/$BUILD_NAME                           \
        2>&1 | tee $BUILD_LOG_DIR/$BUILD_NAME.log
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        # Short link with build version specs only
        ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/${LINK_NAME%%\.*}

        # Longer link with branch name if not a detached HEAD
        if [[ "$GIT_BRANCH" != *"HEAD detached"* ]]; then
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                # Link with branch name if code and build scripts unchanged
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME
            else
                # Link with branch name and "-alpha" for edited code or build script
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME-alpha
            fi
        fi

        echo
        echo -e $MSG_DONE
    else
        echo
        echo -e $MSG_FAIL
    fi
    echo -e "------------------------------------------------------------"
fi

##################################################
###   BLIS based build                         ###
##################################################

if [[ $BUILD_BLIS == "true" ]]; then
    # BLIS libraries
    VER_LAPACK="${UBlue}BLIS${Color_Off}"
    # echo "Bulding Blis"; ./configure --enable-cblas auto; make showconfig; make -j; make check -j; sudo make install
    #PATH_TO_LIBGFORTRAN="/usr/lib64"
    #PATH_TO_BLIS_DIR="/usr/local/include/blis"
    #PATH_TO_BLIS_LIBS="/usr/local/lib"
    #PATH_TO_LAPACK_INCLUDE_DIR="/usr/include"
    #PATH_TO_LAPACK_BINARIES="/usr/lib64/lapack"

    echo
    echo -e "Building $VER_LOCALE $VER_MODE $VER_LAPACK version of Frei..."
    echo
    chpl $FLAG_LOCALE $BLD_FLAG                             \
        -I$PATH_TO_BLIS_DIR                                 \
        -L$PATH_TO_BLIS_LIBS -lblis                         \
        -I$PATH_TO_LAPACK_INCLUDE_DIR                       \
        -L$PATH_TO_LAPACK_BINARIES -llapacke                \
        --main-module "FREI" $SRC_FILES                     \
        -o $BUILD_DIR/$BUILD_NAME                           \
        2>&1 | tee $BUILD_LOG_DIR/$BUILD_NAME.log
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        # Short link with build version specs only
        ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/${LINK_NAME%%\.*}

        # Longer link with branch name if not a detached HEAD
        if [[ "$GIT_BRANCH" != *"HEAD detached"* ]]; then
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                # Link with branch name if code and build scripts unchanged
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME
            else
                # Link with branch name and "-alpha" for edited code or build script
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME-alpha
            fi
        fi

        echo
        echo -e $MSG_DONE
    else
        echo
        echo -e $MSG_FAIL
    fi
    echo -e "------------------------------------------------------------"
fi

##################################################
###   Intel MKL based Build                    ###
##################################################

if [[ $BUILD_INTEL == "true" ]]; then
    # Intel Libraries
    VER_LAPACK="${UBlue}Intel MKL${Color_Off}"
    #PATH_TO_INTEL_LIBGFORTRAN="/opt/intel/oneapi/mkl/latest/lib/intel64"
    PATH_TO_MKL_INCLUDES="/opt/intel/oneapi/mkl/latest/include"
    PATH_TO_MKL_LIBS="/opt/intel/oneapi/mkl/latest/lib/intel64"

    echo
    echo -e "Building $VER_LOCALE $VER_MODE $VER_LAPACK version of Frei..."
    echo
    chpl $FLAG_LOCALE $BLD_FLAG                             \
        --set blasImpl=mkl                                  \
        --set lapackImpl=mkl                                \
        -I$PATH_TO_MKL_INCLUDES                             \
        -L$PATH_TO_MKL_LIBS                                 \
        -lmkl_intel_lp64 -lmkl_sequential -lmkl_core        \
        -lpthread -lm -ld                                   \
        --main-module "FREI" $SRC_FILES                     \
        -o $BUILD_DIR/$BUILD_NAME                           \
        2>&1 | tee $BUILD_LOG_DIR/$BUILD_NAME.log
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        # Short link with build version specs only
        ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/${LINK_NAME%%\.*}

        # Longer link with branch name if not a detached HEAD
        if [[ "$GIT_BRANCH" != *"HEAD detached"* ]]; then
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                # Link with branch name if code and build scripts unchanged
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME
            else
                # Link with branch name and "-alpha" for edited code or build script
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME-alpha
            fi
        fi

        echo
        echo -e $MSG_DONE
    else
        echo
        echo -e $MSG_FAIL
    fi
    echo -e "------------------------------------------------------------"
fi

##################################################
###   AMD AOCL LibM based Build                ###
##################################################

if [[ $BUILD_AMD == "true" ]]; then
    # AMD Libraries
    VER_LAPACK="${UBlue}AMD AOCL LibM${Color_Off}"
    #PATH_TO_AMD_LIBGFORTRAN="/opt/AMD/aocl/aocl-linux-aocc-4.0/lib"
    #PATH_TO_AMD_INCLUDES="/opt/AMD/aocl/aocl-linux-aocc-4.0/include"
    #PATH_TO_AMD_LIBS="/opt/AMD/aocl/aocl-linux-aocc-4.0/lib"
    PATH_TO_AMD_INCLUDES="/opt/AMD/aocl/aocl-linux-gcc-4.0/include"
    PATH_TO_AMD_LIBS="/opt/AMD/aocl/aocl-linux-gcc-4.0/lib"
    #PATH_TO_LAPACK_INCLUDE_DIR="/usr/include"
    #PATH_TO_LAPACK_BINARIES="/usr/lib64/lapack"

    echo
    echo -e "Building $VER_LOCALE $VER_MODE $VER_LAPACK version of Frei..."
    echo
    chpl $FLAG_LOCALE $BLD_FLAG                             \
        -I$PATH_TO_AMD_INCLUDES                             \
        -L$PATH_TO_AMD_LIBS -lblis -lamdlibm -lamdlibm      \
        -I$PATH_TO_LAPACK_INCLUDE_DIR                       \
        -L$PATH_TO_LAPACK_BINARIES -llapacke                \
        --main-module "FREI" $SRC_FILES                     \
        -o $BUILD_DIR/$BUILD_NAME                           \
        2>&1 | tee $BUILD_LOG_DIR/$BUILD_NAME.log
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        # Short link with build version specs only
        ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/${LINK_NAME%%\.*}

        # Longer link with branch name if not a detached HEAD
        if [[ "$GIT_BRANCH" != *"HEAD detached"* ]]; then
            if [ $SOURCE_CHANGES -eq 0 ] && [ $BUILD_SCRIPT_CHANGES -eq 0 ] ; then
                # Link with branch name if code and build scripts unchanged
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME
            else
                # Link with branch name and "-alpha" for edited code or build script
                ln -srf $BUILD_DIR/$BUILD_NAME $BUILD_DIR/$LINK_NAME-alpha
            fi
        fi

        echo
        echo -e $MSG_DONE
    else
        echo
        echo -e $MSG_FAIL
    fi
    echo -e "------------------------------------------------------------"
fi
