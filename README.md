# FREE
*Flux Reconstruction Experimental Edition*
*by RedHatTurtle (Fábio Mallaco Moreira)*

## Intro
This is my first attempt at writing a functional Chapel code. The code is based on a computational project from ITA's
graduate class CC-299 on High-Resolution methods for CFD. The code attempts to simulate a 1D shock-tube experiemnt with
several 2nd or higher order spatial discretization methods, including the Flux Reconstructiom method as described by
H.T. Huyhn from NASA Glenn.

## Compiling and running instructions

```
chpl -o frade --main-module "Control" frade.chpl
./frade --inputFile=inputExample.toml
```
