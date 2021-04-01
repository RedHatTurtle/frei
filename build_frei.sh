#!/bin/bash

#if test -f "frei_opt-build.log" ; then
#    mv frei_opt-build.log frei_opt-build.log.old
#fi
#chpl -o frei_opt --main-module "FREI" src/FREI.chpl          \
#                                      src/input.chpl         \
#                                      src/config.chpl        \
#                                      src/correction.chpl    \
#                                      src/polynomials.chpl   \
#                                      src/interpolation.chpl \
#                                      src/mesh.chpl          \
#                                      src/gmesh.chpl         \
#                                      src/testing.chpl       \
#                                      src/parameters.chpl    |
#    tee frei_opt-build.log


if test -f "frei_dbg-build.log" ; then
    mv frei_dbg-build.log frei_dbg-build.log.old
fi

chpl -o frei_dbg --main-module "FREI" src/FREI.chpl          \
                                      src/input.chpl         \
                                      src/config.chpl        \
                                      src/polynomials.chpl   \
                                      src/correction.chpl    \
                                      src/interpolation.chpl \
                                      src/mesh.chpl          \
                                      src/gmesh.chpl         \
                                      src/testing.chpl       \
                                      src/parameters.chpl    |
    tee frei_dbg-build.log

