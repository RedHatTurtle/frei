#!/bin/bash

#if test -f "frei_opt-build.log" ; then
#    mv frei_opt-build.log frei_opt-build.log.old
#fi
#chpl -o frei_opt --main-module "FREI" src/FREI.chpl          \
#                                      src/fr.chpl            \
#                                      src/frmesh.chpl        \
#                                      src/mesh.chpl          \
#                                      src/gmesh.chpl         \
#                                      src/riemann.chpl       \
#                                      src/flux.chpl          \
#                                      src/correction.chpl    \
#                                      src/interpolation.chpl \
#                                      src/polynomials.chpl   \
#                                      src/config.chpl        \
#                                      src/input.chpl         \
#                                      src/parameters.chpl    \
#                                      src/testing.chpl       |
#    tee frei_opt-build.log


#if test -f "frei_dbg-build.log" ; then
#    mv frei_dbg-build.log frei_dbg-build.log.old
#fi
#
#chpl -o frei_dbg --main-module "FREI" src/FREI.chpl          \
#                                      src/fr.chpl            \
#                                      src/frmesh.chpl        \
#                                      src/mesh.chpl          \
#                                      src/gmesh.chpl         \
#                                      src/riemann.chpl       \
#                                      src/flux.chpl          \
#                                      src/correction.chpl    \
#                                      src/interpolation.chpl \
#                                      src/polynomials.chpl   \
#                                      src/config.chpl        \
#                                      src/input.chpl         \
#                                      src/parameters.chpl    \
#                                      src/testing.chpl       |
#    tee frei_dbg-build.log

if test -f "frei-build.log" ; then
    mv frei-build.log frei-build.log.old
fi

chpl --set lapackImpl=off                        \
     --set blasImpl=off                          \
     -o frei                                     \
     --main-module "FREI" src/FREI.chpl          \
                          src/fr.chpl            \
                          src/output.chpl        \
                          src/frmesh.chpl        \
                          src/mesh.chpl          \
                          src/gmesh.chpl         \
                          src/riemann.chpl       \
                          src/flux.chpl          \
                          src/correction.chpl    \
                          src/interpolation.chpl \
                          src/polynomials.chpl   \
                          src/config.chpl        \
                          src/input.chpl         \
                          src/parameters.chpl    \
                          src/testing.chpl       |
    tee frei-build.log

