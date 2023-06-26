#!/bin/bash

if [ "${-#*i}" == "$-" ]; then
    Color_Off='\033[0m'       # Text Reset
    BRed='\033[1;31m'         # Bold Red
    BGreen='\033[1;32m'       # Bold Green
fi

CBLAS_DIR="/usr/include"
BLAS_LIBS="/usr/lib64/lapack"
LAPACKE_INCLUDE_DIR="/usr/include"
LIBGFORTRAN="/usr/lib64/lapack"
LAPACK_BINARIES="/usr/lib64/lapack"

# AMD Libraries
#CBLAS_DIR="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/include"
#BLAS_LIBS="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/lib64"
#LAPACKE_INCLUDE_DIR="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/include"
#LIBGFORTRAN="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/lib"
#LAPACK_BINARIES="/opt/AMD/aocl/aocl-linux-aocc-3.1.0/lib"

# Intel Libraries
#CBLAS_DIR="/opt/intel/oneapi/mkl/latest/include"
#BLAS_LIBS="/opt/intel/oneapi/mkl/latest/lib/intel64"
#LAPACKE_INCLUDE_DIR="/opt/intel/oneapi/mkl/latest/include"
#LIBGFORTRAN="/opt/intel/oneapi/mkl/latest/lib/intel64"
#LAPACK_BINARIES="/opt/intel/oneapi/mkl/latest/lib/intel64"

# Build all tests and save log
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Polynomials Tests:"
chpl -o tests/polynomials_tests                                   \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Polynomials src/polynomials.chpl               \
                               src/parameters.chpl                \
                               src/testing.chpl                   \
    2>&1 | tee tests/polynomials_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Quadrature Tests:"
chpl -o tests/quadrature_tests                                    \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Quadrature src/quadrature.chpl                 \
                              src/polynomials.chpl                \
                              src/parameters.chpl                 \
                              src/testing.chpl                    \
    2>&1 | tee tests/quadrature_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Interpolation Tests:"
chpl -o tests/interpolation_tests                                 \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Interpolation src/interpolation.chpl           \
                                 src/polynomials.chpl             \
                                 src/mapping.chpl                 \
                                 src/frmesh.chpl                  \
                                 src/temporal.chpl                \
                                 src/flux.chpl                    \
                                 src/mesh.chpl                    \
                                 src/gmesh.chpl                   \
                                 src/input.chpl                   \
                                 src/functions/determinant.chpl   \
                                 src/functions/sorttuple.chpl     \
                                 src/parameters.chpl              \
                                 src/testing.chpl                 \
    2>&1 | tee tests/interpolation_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Projection Tests:"
chpl -o tests/projection_tests                                    \
     --warnings                                                   \
     --warn-unstable                                              \
     -I $CBLAS_DIR                                                \
     -L $BLAS_LIBS -lcblas                                        \
     -I $LAPACKE_INCLUDE_DIR                                      \
     -L $LIBGFORTRAN -lgfortran                                   \
     -L $LAPACK_BINARIES -llapacke -llapack -lcblas               \
     --main-module Projection src/projection.chpl                 \
                              src/quadrature.chpl                 \
                              src/polynomials.chpl                \
                              src/parameters.chpl                 \
                              src/testing.chpl                    \
    2>&1 | tee tests/projection_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Correction Tests:"
chpl -o tests/correction_tests                                    \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Correction src/correction.chpl                 \
                              src/interpolation.chpl              \
                              src/polynomials.chpl                \
                              src/input.chpl                      \
                              src/mapping.chpl                    \
                              src/frmesh.chpl                     \
                              src/temporal.chpl                   \
                              src/flux.chpl                       \
                              src/mesh.chpl                       \
                              src/gmesh.chpl                      \
                              src/config.chpl                     \
                              src/functions/determinant.chpl      \
                              src/functions/sorttuple.chpl        \
                              src/parameters.chpl                 \
                              src/testing.chpl                    \
    2>&1 | tee tests/correction_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Dimensional Tests:"
chpl -o tests/dimensional_tests                                   \
     --warnings                                                   \
     --warn-unstable                                              \
     --detailed-errors                                            \
     --main-module Dimensional     src/dimensional.chpl           \
                                   src/flux.chpl                  \
                                   src/mesh.chpl                  \
                                   src/gmesh.chpl                 \
                                   src/input.chpl                 \
                                   src/functions/determinant.chpl \
                                   src/functions/sorttuple.chpl   \
                                   src/parameters.chpl            \
                                   src/testing.chpl               \
    2>&1 | tee tests/dimensional_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Gmesh Tests:"
chpl -o tests/gmesh_tests                                         \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Gmesh src/gmesh.chpl                           \
                         src/mesh.chpl                            \
                         src/frmesh.chpl                          \
                         src/temporal.chpl                        \
                         src/flux.chpl                            \
                         src/mapping.chpl                         \
                         src/interpolation.chpl                   \
                         src/polynomials.chpl                     \
                         src/input.chpl                           \
                         src/functions/determinant.chpl           \
                         src/functions/sorttuple.chpl             \
                         src/parameters.chpl                      \
                         src/testing.chpl                         \
    2>&1 | tee tests/gmesh_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Mesh Tests:"
chpl -o tests/mesh_tests                                          \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Mesh src/mesh.chpl                             \
                        src/gmesh.chpl                            \
                        src/polynomials.chpl                      \
                        src/interpolation.chpl                    \
                        src/mapping.chpl                          \
                        src/frmesh.chpl                           \
                        src/temporal.chpl                         \
                        src/flux.chpl                             \
                        src/input.chpl                            \
                        src/functions/determinant.chpl            \
                        src/functions/sorttuple.chpl              \
                        src/parameters.chpl                       \
                        src/testing.chpl                          \
    2>&1 | tee tests/mesh_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building FRMesh Tests:"
chpl -o tests/frmesh_tests                                        \
     --warnings                                                   \
     --warn-unstable                                              \
     -I $CBLAS_DIR                                                \
     -L $BLAS_LIBS -lcblas                                        \
     -I $LAPACKE_INCLUDE_DIR                                      \
     -L $LIBGFORTRAN -lgfortran                                   \
     -L $LAPACK_BINARIES -llapacke -llapack -lcblas               \
     --main-module FRMesh src/frmesh.chpl                         \
                          src/temporal.chpl                       \
                          src/flux.chpl                           \
                          src/interpolation.chpl                  \
                          src/polynomials.chpl                    \
                          src/mesh.chpl                           \
                          src/gmesh.chpl                          \
                          src/mapping.chpl                        \
                          src/config.chpl                         \
                          src/input.chpl                          \
                          src/functions/determinant.chpl          \
                          src/functions/sorttuple.chpl            \
                          src/parameters.chpl                     \
                          src/testing.chpl                        \
    2>&1 | tee tests/frmesh_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Mapping Tests:"
chpl -o tests/mapping_tests                                       \
     --warnings                                                   \
     --warn-unstable                                              \
     -I $CBLAS_DIR                                                \
     -L $BLAS_LIBS -lcblas                                        \
     -I $LAPACKE_INCLUDE_DIR                                      \
     -L $LIBGFORTRAN -lgfortran                                   \
     -L $LAPACK_BINARIES -llapacke -llapack -lcblas               \
     --main-module Mapping src/mapping.chpl                       \
                           src/interpolation.chpl                 \
                           src/polynomials.chpl                   \
                           src/frmesh.chpl                        \
                           src/temporal.chpl                      \
                           src/flux.chpl                          \
                           src/mesh.chpl                          \
                           src/gmesh.chpl                         \
                           src/input.chpl                         \
                           src/functions/determinant.chpl         \
                           src/functions/sorttuple.chpl           \
                           src/parameters.chpl                    \
                           src/testing.chpl                       \
     2>&1 | tee tests/mapping_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Flux Tests:"
chpl -o tests/flux_tests                                          \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Flux src/flux.chpl                             \
                        src/parameters.chpl                       \
                        src/testing.chpl                          \
    2>&1 | tee tests/flux_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Riemann Tests:"
chpl -o tests/riemann_tests                                       \
     --warnings                                                   \
     --warn-unstable                                              \
     -I $CBLAS_DIR                                                \
     -L $BLAS_LIBS -lcblas                                        \
     -I $LAPACKE_INCLUDE_DIR                                      \
     -L $LIBGFORTRAN -lgfortran                                   \
     -L $LAPACK_BINARIES -llapacke -llapack -lcblas               \
     --main-module Riemann src/riemann.chpl                       \
                           src/interpolation.chpl                 \
                           src/polynomials.chpl                   \
                           src/flux.chpl                          \
                           src/input.chpl                         \
                           src/frmesh.chpl                        \
                           src/temporal.chpl                      \
                           src/flux.chpl                          \
                           src/mesh.chpl                          \
                           src/gmesh.chpl                         \
                           src/mapping.chpl                       \
                           src/functions/determinant.chpl         \
                           src/functions/sorttuple.chpl           \
                           src/parameters.chpl                    \
                           src/testing.chpl                       \
    2>&1 | tee tests/riemann_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo "Building Ringleb Tests:"
chpl -o tests/ringleb_tests                                       \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Ringleb src/ringleb.chpl                       \
                           src/interpolation.chpl                 \
                           src/polynomials.chpl                   \
                           src/input.chpl                         \
                           src/frmesh.chpl                        \
                           src/temporal.chpl                      \
                           src/flux.chpl                          \
                           src/mesh.chpl                          \
                           src/gmesh.chpl                         \
                           src/mapping.chpl                       \
                           src/functions/determinant.chpl         \
                           src/functions/sorttuple.chpl           \
                           src/parameters.chpl                    \
                           src/testing.chpl                       \
   2>&1 | tee tests/ringleb_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Init Tests:"
chpl -o tests/init_tests                                          \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Init src/init.chpl                             \
                        src/interpolation.chpl                    \
                        src/polynomials.chpl                      \
                        src/ringleb.chpl                          \
                        src/flux.chpl                             \
                        src/config.chpl                           \
                        src/input.chpl                            \
                        src/frmesh.chpl                           \
                        src/temporal.chpl                         \
                        src/mesh.chpl                             \
                        src/gmesh.chpl                            \
                        src/mapping.chpl                          \
                        src/functions/determinant.chpl            \
                        src/functions/sorttuple.chpl              \
                        src/parameters.chpl                       \
                        src/testing.chpl                          \
    2>&1 | tee tests/init_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Boundary Tests:"
chpl -o tests/boundary_tests                                      \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module Boundary src/boundary.chpl                     \
                            src/interpolation.chpl                \
                            src/polynomials.chpl                  \
                            src/init.chpl                         \
                            src/ringleb.chpl                      \
                            src/flux.chpl                         \
                            src/frmesh.chpl                       \
                            src/temporal.chpl                     \
                            src/mesh.chpl                         \
                            src/gmesh.chpl                        \
                            src/mapping.chpl                      \
                            src/config.chpl                       \
                            src/input.chpl                        \
                            src/functions/determinant.chpl        \
                            src/functions/sorttuple.chpl          \
                            src/parameters.chpl                   \
                            src/testing.chpl                      \
    2>&1 | tee tests/boundary_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building Limiter Tests:"
chpl -o tests/limiter_tests                                       \
     --warnings                                                   \
     --warn-unstable                                              \
     -I $CBLAS_DIR                                                \
     -L $BLAS_LIBS -lcblas                                        \
     -I $LAPACKE_INCLUDE_DIR                                      \
     -L $LIBGFORTRAN -lgfortran                                   \
     -L $LAPACK_BINARIES -llapacke -llapack -lcblas               \
     --main-module Limiter src/limiter.chpl                       \
                           src/projection.chpl                    \
                           src/quadrature.chpl                    \
                           src/polynomials.chpl                   \
                           src/functions/sorttuple.chpl           \
                           src/parameters.chpl                    \
                           src/testing.chpl                       \
    2>&1 | tee tests/limiter_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Building FR Tests:"
chpl -o tests/fr_tests                                            \
     --warnings                                                   \
     --warn-unstable                                              \
     --main-module FR src/fr.chpl                                 \
                      src/testing.chpl                            \
    2>&1 | tee tests/fr_build.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo -e "${BGreen}Done${Color_Off}"
else
  echo -e "${BRed}Fail${Color_Off}"
fi
echo
echo "------------------------------------------------------------"
echo
echo -e "Running Tests..."
echo
# Run tests and output to file
./tests/polynomials_tests   &> tests/polynomials_tests.log
./tests/quadrature_tests    &> tests/quadrature_tests.log
./tests/interpolation_tests &> tests/interpolation_tests.log
./tests/projection_tests    &> tests/projection_tests.log

./tests/correction_tests    &> tests/correction_tests.log
./tests/dimensional_tests   &> tests/dimensional_tests.log

./tests/gmesh_tests         &> tests/gmesh_tests.log
./tests/mesh_tests          &> tests/mesh_tests.log
./tests/frmesh_tests        &> tests/frmesh_tests.log
./tests/mapping_tests       &> tests/mapping_tests.log

./tests/flux_tests          &> tests/flux_tests.log
./tests/riemann_tests       &> tests/riemann_tests.log

./tests/ringleb_tests       &> tests/ringleb_tests.log
./tests/init_tests          &> tests/init_tests.log
./tests/boundary_tests      &> tests/boundary_tests.log
./tests/limiter_tests       &> tests/limiter_tests.log

./tests/fr_tests            &> tests/fr_tests.log
echo -e "Done"
echo
echo "------------------------------------------------------------"
echo
