#!/bin/bash

# Build all tests and save log
chpl -o polynomials_tests --main-module Polynomials src/polynomials.chpl src/testing.chpl src/parameters.chpl |
    tee polynomials_build_test.log
echo
chpl -o correction_tests --main-module Correction src/correction.chpl src/polynomials.chpl src/testing.chpl src/parameters.chpl |
    tee correction_build_test.log
echo
chpl -o interpolation_tests --main-module Interpolation src/interpolation.chpl src/polynomials.chpl src/parameters.chpl src/testing.chpl |
    tee interpolation_build_test.log
echo
chpl -o gmesh_tests --main-module Gmesh src/gmesh.chpl src/parameters.chpl |
    tee gmesh_build_test.log
echo
chpl -o mesh_tests --main-module Mesh src/mesh.chpl src/gmesh.chpl src/testing.chpl src/parameters.chpl |
    tee mesh_build_test.log

# Run tests and output to file
./polynomials_tests   &> polynomials_test.log
./correction_tests    &> correction_test.log
./interpolation_tests &> interpolation_test.log

./gmesh_tests         &> gmesh_test.log
./mesh_tests          &> mesh_test.log

