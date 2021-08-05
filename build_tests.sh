#!/bin/bash

PATH_TO_CBLAS_DIR="/usr/include"
PATH_TO_BLAS_LIBS="/usr/lib64"
PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
PATH_TO_LIBGFORTRAN="/usr/lib64"
PATH_TO_LAPACK_BINARIES="/usr/lib64"

# Build all tests and save log
echo "Building Polynomials Tests:"
chpl -o polynomials_tests                          \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Polynomials src/polynomials.chpl src/testing.chpl src/parameters.chpl |
    tee polynomials_build.log
echo "done"
echo
echo "Building Quadrature Tests:"
chpl -o quadrature_tests                           \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Quadrature src/quadrature.chpl src/polynomials.chpl src/testing.chpl src/parameters.chpl |
    tee quadrature_build.log
echo "done"
echo
echo "Building Interpolation Tests:"
chpl -o interpolation_tests                        \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Interpolation src/interpolation.chpl src/polynomials.chpl src/parameters.chpl src/testing.chpl |
    tee interpolation_build.log
echo "done"
echo
echo "Building Projection Tests:"
chpl -o projection_tests                           \
     --warnings                                    \
     --warn-unstable                               \
     -I/usr/include                                \
     -L/usr/lib64 -lcblas                          \
     -I/usr/include                                \
     -L/usr/lib64 -lgfortran                       \
     -L/usr/lib64 -llapacke -llapack -lcblas       \
     --main-module Projection src/projection.chpl  \
                              src/quadrature.chpl  \
                              src/polynomials.chpl \
                              src/testing.chpl     \
                              src/parameters.chpl  |
    tee projection_build.log
echo "done"
echo
echo "Building Correction Tests:"
chpl -o correction_tests                           \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Correction src/correction.chpl src/polynomials.chpl src/testing.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/config.chpl src/parameters.chpl |
    tee correction_build.log
echo "done"
echo
echo "Building Gmesh Tests:"
chpl -o gmesh_tests                                \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Gmesh src/gmesh.chpl src/parameters.chpl |
    tee gmesh_build.log
echo "done"
echo
echo "Building Mesh Tests:"
chpl -o mesh_tests                                 \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Mesh src/mesh.chpl src/gmesh.chpl src/testing.chpl src/parameters.chpl |
    tee mesh_build.log
echo "done"
echo
echo "Building FRMesh Tests:"
chpl -o frmesh_tests                               \
     --warnings                                    \
     --warn-unstable                               \
     --main-module FRMesh src/frmesh.chpl src/polynomials.chpl src/mesh.chpl src/gmesh.chpl src/testing.chpl src/config.chpl src/input.chpl src/parameters.chpl |
    tee frmesh_build.log
echo "done"
echo
echo "Building Flux Tests:"
chpl -o flux_tests                                 \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Flux src/flux.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/parameters.chpl |
    tee flux_build.log
echo "done"
echo
echo "Building Riemann Tests:"
chpl -o riemann_tests                              \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Riemann src/riemann.chpl src/flux.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/parameters.chpl |
    tee riemann_build.log
echo "done"
echo
echo "Building Init Tests:"
chpl -o init_tests                                 \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Init src/init.chpl src/flux.chpl src/config.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/parameters.chpl |
    tee init_build.log
echo "done"
echo
echo "Building Boundary Tests:"
chpl -o boundary_tests                             \
     --warnings                                    \
     --warn-unstable                               \
     --main-module Boundary src/boundary.chpl src/init.chpl src/flux.chpl src/mesh.chpl src/gmesh.chpl src/config.chpl src/input.chpl src/parameters.chpl |
    tee boundary_build.log
echo "done"
echo
echo "Building Limiter Tests:"
chpl -o limiter_tests                              \
     --warnings                                    \
     --warn-unstable                               \
     -I/usr/include                                \
     -L/usr/lib64 -lcblas                          \
     -I/usr/include                                \
     -L/usr/lib64 -lgfortran                       \
     -L/usr/lib64 -llapacke -llapack -lcblas       \
     --main-module Limiter src/limiter.chpl src/projection.chpl src/quadrature.chpl src/polynomials.chpl src/parameters.chpl src/testing.chpl |
    tee limiter_build.log
echo "done"
echo
echo "Building FR Tests:"
chpl -o fr_tests                                   \
     --warnings                                    \
     --warn-unstable                               \
     --main-module FR src/fr.chpl |
    tee fr_build.log
echo "done"
echo

echo "Running Tests"
# Run tests and output to file
./polynomials_tests   &> polynomials_tests.log
./quadrature_tests    &> quadrature_tests.log
./interpolation_tests &> interpolation_tests.log
./projection_tests    &> projection_tests.log

./correction_tests    &> correction_tests.log

./gmesh_tests         &> gmesh_tests.log
./mesh_tests          &> mesh_tests.log
./frmesh_tests        &> frmesh_tests.log

./flux_tests          &> flux_tests.log
./riemann_tests       &> riemann_tests.log

./init_tests          &> init_tests.log
./boundary_tests      &> boundary_tests.log
./limiter_tests       &> limiter_tests.log

./fr_tests            &> fr_tests.log
echo "Done"
