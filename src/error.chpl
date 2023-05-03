module ErrorCalc
{
  import FRMesh.fr_mesh_c;

  var meshVol : real = 0.0;

  proc error_calc(problemType: int, frMesh : fr_mesh_c) : [] (string, real)
  {
    var errors_d : domain(1);
    var errors : [errors_d] (string, real);

    use Parameters.ParamInput;
    select problemType
    {
      when IC_RINGLEB do
        return error_ringleb(frMesh);
    }

    return errors;
  }

  proc error_ringleb(frMesh : fr_mesh_c) : [] (string, real)
  {
    use Parameters.ParamMesh;
    use Quadrature;
    use Ringleb;
    use Flux;

    var error_d : domain(1);
    var errors  : [error_d] (string,real);

    var l1ErrorAbs, l1ErrorRel : [1..frMesh.nVars] real = 0;
    var l2ErrorAbs, l2ErrorRel : [1..frMesh.nVars] real = 0;
    var lfErrorAbs, lfErrorRel : [1..frMesh.nVars] real = 0;
    var l1EntrErrorAbs, l2EntrErrorAbs, lfEntrErrorAbs : real;
    var l1EntrErrorRel, l2EntrErrorRel, lfEntrErrorRel : real;

    // Calculate mesh length/area/volume
    if (meshVol <= 0.0) then
      for cellIdx in frMesh.cellList.domain {
        var unit : [1..frMesh.cellSPidx[cellIdx, 2]] real = 1.0;
        var jacobian = frMesh.jacSP[frMesh.cellSPidx[cellIdx, 1].. #frMesh.cellSPidx[cellIdx, 2]];
        meshVol += integrate( unit, jacobian, TOPO_QUAD, frMesh.solOrder);
      }

    // Loop through cells
    for cellIdx in frMesh.cellList.domain
    {
      // Get the reference values for the points in this cell
      var solRef : [1..frMesh.nVars, 1..frMesh.cellSPidx[cellIdx, 2]] real;
      for spIdx in 1..frMesh.cellSPidx[cellIdx, 2] do
        solRef[.., spIdx] = ringleb_sol(frMesh.xyzSP[frMesh.cellSPidx[cellIdx, 1]+spIdx-1, 1..2]);

      var entrRef : real = entropy_cv(solRef[.., 1]);

      // Calculate the errors at each point
      var jacobian   = frMesh.jacSP[frMesh.cellSPidx[cellIdx, 1].. #frMesh.cellSPidx[cellIdx, 2]];
      var ptErrorAbs = abs(frMesh.solSP[1..4, frMesh.cellSPidx[cellIdx, 1].. #frMesh.cellSPidx[cellIdx, 2]] - solRef[1..4, 1..frMesh.cellSPidx[cellIdx, 2]]);
      var ptErrorRel = ptErrorAbs/abs(solRef[1..4, 1..frMesh.cellSPidx[cellIdx, 2]]);

      var ptEntrErrorAbs : [1..frMesh.cellSPidx[cellIdx, 2]] real;
      for spIdx in 1..frMesh.cellSPidx[cellIdx, 2] do
        ptEntrErrorAbs[spIdx] = abs(entropy_cv(frMesh.solSP[1..4, frMesh.cellSPidx[cellIdx, 1]+spIdx-1]) - entrRef);

      var ptEntrErrorRel = ptEntrErrorAbs/abs(entrRef);

      // Absolute errors
      {
        // L2 Errors
        {
          //Conserved Variables
          for varIdx in 1..frMesh.nVars
          {
            var ptErrorAbs2 = ptErrorAbs[varIdx, ..]**2;
            l2ErrorAbs[varIdx] += integrate( ptErrorAbs2, jacobian, TOPO_QUAD, frMesh.solOrder );
          }

          //Entropy
          var ptEntrErrorAbs2 = ptEntrErrorAbs**2;
          l2EntrErrorAbs += integrate( ptEntrErrorAbs2, jacobian, TOPO_QUAD, frMesh.solOrder );
        }
        // Lf Errors
        {
          //Conserved Variables
          for varIdx in 1..frMesh.nVars
          {
            var ptErrorAbsM : real = max reduce (ptErrorAbs[varIdx, ..]);
            lfErrorAbs = max(lfErrorAbs, ptErrorAbsM);
          }

          //Entropy
          var ptEntrErrorAbsM : real = max reduce (ptEntrErrorAbs);
          lfEntrErrorAbs = max(lfEntrErrorAbs, ptEntrErrorAbsM);
        }
        // L1 Errors
        {
          //Conserved Variables
          for varIdx in 1..frMesh.nVars
          {
            var ptErrorAbs1 = ptErrorAbs[varIdx, ..];
            l1ErrorAbs[varIdx] += integrate( ptErrorAbs1, jacobian, TOPO_QUAD, frMesh.solOrder );
          }

          //Entropy
          var ptEntrErrorAbs1 = ptEntrErrorAbs;
          l1EntrErrorAbs += integrate( ptEntrErrorAbs1, jacobian, TOPO_QUAD, frMesh.solOrder );
        }
      }

      // Relative errors
      {
        // L2 Errors
        {
          //Conserved Variables
          for varIdx in 1..frMesh.nVars
          {
            var ptErrorRel2 = ptErrorRel[varIdx, ..]**2;
            l2ErrorRel[varIdx] += integrate( ptErrorRel2, jacobian, TOPO_QUAD, frMesh.solOrder );
          }

          //Entropy
          var ptEntrErrorRel2 = ptEntrErrorRel**2;
          l2EntrErrorRel += integrate( ptEntrErrorRel2, jacobian, TOPO_QUAD, frMesh.solOrder );
        }
        // Lf Errors
        {
          //Conserved Variables
          for varIdx in 1..frMesh.nVars
          {
            var ptErrorRelM : real = max reduce (ptErrorRel[varIdx, ..]);
            lfErrorRel = max(lfErrorRel, ptErrorRelM);
          }

          //Entropy
          var ptEntrErrorRelM : real = max reduce (ptEntrErrorRel);
          lfEntrErrorRel = max(lfEntrErrorRel, ptEntrErrorRelM);
        }
        // L1 Errors
        {
          //Conserved Variables
          for varIdx in 1..frMesh.nVars
          {
            var ptErrorRel1 = ptErrorRel[varIdx, ..];
            l1ErrorRel[varIdx] += integrate( ptErrorRel1, jacobian, TOPO_QUAD, frMesh.solOrder );
          }

          //Entropy
          var ptEntrErrorRel1 = ptEntrErrorRel;
          l1EntrErrorRel += integrate( ptEntrErrorRel1, jacobian, TOPO_QUAD, frMesh.solOrder );
        }
      }
    }


    // Absolute errors
    {
      //Conserved Variables
      for varIdx in 1..frMesh.nVars
      {
        error_d = {1..errors.domain.high+1};
        try! errors[errors.domain.high] = ("L2(Error(Sol[%1i]))".format(varIdx), l2ErrorAbs[varIdx]/meshVol);
      }
      // Entropy
      error_d = {1..errors.domain.high+1};
      try! errors[errors.domain.high] = ("L2(Error(%6s))".format("Entr"), l2EntrErrorAbs/meshVol);

      //Conserved Variables
      for varIdx in 1..frMesh.nVars
      {
        error_d = {1..errors.domain.high+1};
        try! errors[errors.domain.high] = ("Lf(Error(Sol[%1i]))".format(varIdx), lfErrorAbs[varIdx]/meshVol);
      }
      // Entropy
      error_d = {1..errors.domain.high+1};
      try! errors[errors.domain.high] = ("Lf(Error(%6s))".format("Entr"), lfEntrErrorAbs/meshVol);

      //Conserved Variables
      for varIdx in 1..frMesh.nVars
      {
        error_d = {1..errors.domain.high+1};
        try! errors[errors.domain.high] = ("L1(Error(Sol[%1i]))".format(varIdx), l1ErrorAbs[varIdx]/meshVol);
      }
      // Entropy
      error_d = {1..errors.domain.high+1};
      try! errors[errors.domain.high] = ("L1(Error(%6s))".format("Entr"), l1EntrErrorAbs/meshVol);
    }

    // Relative errors
    {
      //Conserved Variables
      for varIdx in 1..frMesh.nVars
      {
        error_d = {1..errors.domain.high+1};
        try! errors[errors.domain.high] = ("L2(%%Error(Sol[%1i]))".format(varIdx), l2ErrorRel[varIdx]/meshVol);
      }
      // Entropy
      error_d = {1..errors.domain.high+1};
      try! errors[errors.domain.high] = ("L2(%%Error(%6s))".format("Entr"), l2EntrErrorRel/meshVol);

      //Conserved Variables
      for varIdx in 1..frMesh.nVars
      {
        error_d = {1..errors.domain.high+1};
        try! errors[errors.domain.high] = ("Lf(%%Error(Sol[%1i]))".format(varIdx), lfErrorRel[varIdx]/meshVol);
      }
      // Entropy
      error_d = {1..errors.domain.high+1};
      try! errors[errors.domain.high] = ("Lf(%%Error(%6s))".format("Entr"), lfEntrErrorRel/meshVol);

      //Conserved Variables
      for varIdx in 1..frMesh.nVars
      {
        error_d = {1..errors.domain.high+1};
        try! errors[errors.domain.high] = ("L1(%%Error(Sol[%1i]))".format(varIdx), l1ErrorRel[varIdx]/meshVol);
      }
      // Entropy
      error_d = {1..errors.domain.high+1};
      try! errors[errors.domain.high] = ("L1(%%Error(%6s))".format("Entr"), l1EntrErrorRel/meshVol);
    }

    // Return some stats for the command line logs
    return errors;
  }
}
