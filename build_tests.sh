#!/bin/bash

# Build all tests and save log
echo "Building Polynomials Tests:"
chpl -o polynomials_tests --main-module Polynomials src/polynomials.chpl src/testing.chpl src/parameters.chpl |
    tee polynomials_build_test.log
echo "done"
echo
echo "Building Correction Tests:"
chpl -o correction_tests --main-module Correction src/correction.chpl src/polynomials.chpl src/testing.chpl src/input.chpl src/mesh.chpl src/gmesh.chpl src/config.chpl src/parameters.chpl |
    tee correction_build_test.log
echo "done"
echo
echo "Building Interpolation Tests:"
chpl -o interpolation_tests --main-module Interpolation src/interpolation.chpl src/polynomials.chpl src/parameters.chpl src/testing.chpl |
    tee interpolation_build_test.log
echo "done"
echo
echo "Building Gmesh Tests:"
chpl -o gmesh_tests --main-module Gmesh src/gmesh.chpl src/parameters.chpl |
    tee gmesh_build_test.log
echo "done"
echo
echo "Building Mesh Tests:"
chpl -o mesh_tests --main-module Mesh src/mesh.chpl src/gmesh.chpl src/testing.chpl src/parameters.chpl |
    tee mesh_build_test.log
echo "done"
echo
echo "Building FlowVars Tests:"
chpl -o flowvars_tests --main-module FlowVars src/flowVars.chpl src/mesh.chpl src/gmesh.chpl src/testing.chpl src/config.chpl src/input.chpl src/parameters.chpl |
    tee flowvars_build_test.log
echo "done"
echo

echo "Running Tests"
# Run tests and output to file
./polynomials_tests   &> polynomials_test.log
./correction_tests    &> correction_test.log
./interpolation_tests &> interpolation_test.log

./gmesh_tests         &> gmesh_test.log
./mesh_tests          &> mesh_test.log
./flowvars_tests      &> flowvars_test.log
echo "Done"
