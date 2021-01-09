chpl -o frei --main-module "Control" frei.chpl          \
                                     correction.chpl    \
                                     polynomials.chpl   \
                                     interpolation.chpl \
                                     testing.chpl       \
                                     input.chpl         \
                                     mesh.chpl          \
                                     gmesh.chpl         \
                                     parameters.chpl    |
    tee frei_build_dbg.log



chpl -o frei --main-module "Control" frei.chpl          \
                                     correction.chpl    \
                                     polynomials.chpl   \
                                     interpolation.chpl \
                                     testing.chpl       \
                                     input.chpl         \
                                     mesh.chpl          \
                                     gmesh.chpl         \
                                     parameters.chpl    |
    tee frei_build_opt.log

