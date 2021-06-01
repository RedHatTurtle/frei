#!/bin/bash

# Build all tests and save log
echo "Building Polynomials Tests:"
chpl -o polynomials_tests --main-module Polynomials src/polynomials.chpl src/testing.chpl src/parameters.chpl |
    tee polynomials_build.log
echo "done"
echo
echo "Building Interpolation Tests:"
chpl -o interpolation_tests --main-module Interpolation src/interpolation.chpl src/polynomials.chpl src/parameters.chpl src/testing.chpl |
    tee interpolation_build.log
echo "done"
echo
echo "Building Correction Tests:"
chpl -o correction_tests --main-module Correction src/correction.chpl src/polynomials.chpl src/testing.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/config.chpl src/parameters.chpl |
    tee correction_build.log
echo "done"
echo
echo "Building Gmesh Tests:"
chpl -o gmesh_tests --main-module Gmesh src/gmesh.chpl src/parameters.chpl |
    tee gmesh_build.log
echo "done"
echo
echo "Building Mesh Tests:"
chpl -o mesh_tests --main-module Mesh src/mesh.chpl src/gmesh.chpl src/testing.chpl src/parameters.chpl |
    tee mesh_build.log
echo "done"
echo
echo "Building FRMesh Tests:"
chpl -o frmesh_tests --main-module FRMesh src/frmesh.chpl src/polynomials.chpl src/mesh.chpl src/gmesh.chpl src/testing.chpl src/config.chpl src/input.chpl src/parameters.chpl |
    tee frmesh_build.log
echo "done"
echo
echo "Building Flux Tests:"
chpl -o flux_tests --main-module Flux src/flux.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/parameters.chpl |
    tee flux_build.log
echo "done"
echo
echo "Building Riemann Tests:"
chpl -o riemann_tests --main-module Riemann src/riemann.chpl src/flux.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/parameters.chpl |
    tee riemann_build.log
echo "done"
echo
echo "Building Init Tests:"
chpl -o init_tests --main-module Init src/init.chpl src/flux.chpl src/config.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/parameters.chpl |
    tee init_build.log
echo "done"
echo
echo "Building Boundary Tests:"
chpl -o fr_tests --main-module Boundary src/boundary.chpl src/init.chpl src/flux.chpl src/mesh.chpl src/gmesh.chpl src/config.chpl src/input.chpl src/parameters.chpl |
    tee boundary_build.log
echo "done"
echo
echo "Building FR Tests:"
chpl -o fr_tests --main-module FR src/fr.chpl |
    tee fr_build.log
echo "done"
echo

echo "Running Tests"
# Run tests and output to file
./polynomials_tests   &> polynomials_tests.log
./interpolation_tests &> interpolation_tests.log
./correction_tests    &> correction_tests.log

./gmesh_tests         &> gmesh_tests.log
./mesh_tests          &> mesh_tests.log
./frmesh_tests        &> frmesh_tests.log

./flux_tests          &> flux_tests.log
./riemann_tests       &> riemann_tests.log

./init_tests          &> init_tests.log
./boundary_tests      &> boundary_tests.log

./fr_tests            &> fr_tests.log
echo "Done"
