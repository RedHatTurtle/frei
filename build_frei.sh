#!/bin/bash

##################################################
###   Debug Build                              ###
##################################################

if test -f "frei_dbg-build.log" ; then
    mv frei_dbg-build.log frei_dbg-build.log.old
fi

chpl --set lapackImpl=off                        \
     --set blasImpl=off                          \
     --warnings                                  \
     --warn-unstable                             \
     -o frei_dbg                                 \
     --main-module "FREI" src/FREI.chpl          \
                          src/fr.chpl            \
                          src/output.chpl        \
                          src/boundary.chpl      \
                          src/frmesh.chpl        \
                          src/mesh.chpl          \
                          src/gmesh.chpl         \
                          src/riemann.chpl       \
                          src/flux.chpl          \
                          src/correction.chpl    \
                          src/interpolation.chpl \
                          src/polynomials.chpl   \
                          src/sourceterm.chpl    \
                          src/init.chpl          \
                          src/config.chpl        \
                          src/input.chpl         \
                          src/parameters.chpl    \
                          src/testing.chpl       |
    tee frei_dbg-build.log

##################################################
###   Optimized / Prduction Build              ###
##################################################

if test -f "frei_opt-build.log" ; then
    mv frei_opt-build.log frei_opt-build.log.old
fi
chpl --set lapackImpl=off                        \
     --set blasImpl=off                          \
     --fast                                      \
     -o frei_opt                                 \
     --main-module "FREI" src/FREI.chpl          \
                          src/fr.chpl            \
                          src/output.chpl        \
                          src/boundary.chpl      \
                          src/frmesh.chpl        \
                          src/mesh.chpl          \
                          src/gmesh.chpl         \
                          src/riemann.chpl       \
                          src/flux.chpl          \
                          src/correction.chpl    \
                          src/interpolation.chpl \
                          src/polynomials.chpl   \
                          src/sourceterm.chpl    \
                          src/init.chpl          \
                          src/config.chpl        \
                          src/input.chpl         \
                          src/parameters.chpl    \
                          src/testing.chpl       |
    tee frei_opt-build.log
