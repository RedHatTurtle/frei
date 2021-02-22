prototype module Solution {
  use Input;

  class flux_reconstruction_c
  {
    var nDims : int;
    var nVars : int;

    // Domains
    var  xyzSP_d : domain(rank=2, idxType=int);     // {1..nSPs, 1..nDims}
    var  xyzFP_d : domain(rank=2, idxType=int);     // {1..nSPs, 1..nDims}

    var  solSP_d : domain(rank=2, idxType=int);     // {1..nSPs, 1..nVars}
    var  solFP_d : domain(rank=3, idxType=int);     // {1..nFPs, 1..nVars}
    var  flxFP_d : domain(rank=2, idxType=int);     // {1..nFPs, 1..nVars}
    var dSolSP_d : domain(rank=2, idxType=int);     // {1..nSPs, 1..nVars, 1..nVars}
    var dSolFP_d : domain(rank=4, idxType=int);     // {1..2, 1..nFPs, nVars, nVars}

    var cellSPidx_d : domain(rank=2, idxType=int);  // {1..nCells, 1..3}
    var faceFPidx_d : domain(rank=2, idxType=int);  // {1..nFaces, 1..3}

    // FR solver variables
    var xyzSP : [xyzSP_d] real;
    var xyzFP : [xyzFP_d] real;

    var oldSolSP : [ solSP_d] real;     // Backup of the solution at the beginning of residue calculation
    var    solSP : [ solSP_d] real;     // Conserved variables at SPs
    var    solFP : [ solFP_d] real;     // Conserved variables at FPs (0-Left / 1-right)
    var    flxFP : [ flxFP_d] real;     // Unique convective flux at FPs
    var   dSolSP : [dSolSP_d] real;     // Gradient, at the SPs, of the discontinuous solution interpolation
    var   dSolFP : [dSolFP_d] real;     // Gradient, at the FPs, of the discontinuous flux reconstruction
    var    resSP : [ solSP_d] real;     // conserved variables residual

    var cellSPidx : [cellSPidx_d] int;  // Index of the first SP and number of SPs of a cell
    var faceFPidx : [faceFPidx_d] int;  // Index of the first FP and number of FPs of a face
  }
}
