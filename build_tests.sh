#!/bin/bash

# Build all tests and save log
chpl -o polynomials_tests --main-module Polynomials polynomials.chpl testing.chpl parameters.chpl |
    tee polynomials_build_test.log
chpl -o correction_tests --main-module Correction correction.chpl polynomials.chpl testing.chpl parameters.chpl |
    tee correction_build_test.log
chpl -o interpolation_tests --main-module Interpolation interpolation.chpl polynomials.chpl parameters.chpl testing.chpl |
    tee interpolation_build_test.log

chpl -o gmesh_tests --main-module Gmesh gmesh.chpl parameters.chpl |
    tee gmesh_build_test.log
chpl -o mesh_tests --main-module Mesh mesh.chpl testing.chpl |
    tee mesh_build_test.log

# Run tests and output to file
./polynomials_tests   &> polynomials_test.log
./correction_tests    &> correction_test.log
./interpolation_tests &> interpolation_test.log

./gmesh_tests         &> gmesh_test.log
./mesh_tests          &> mesh_test.log

