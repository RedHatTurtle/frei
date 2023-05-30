module Temporal
{
  /*
    This module contains the time iterating functions.
    At the moment no methods that require storing intermediate stages are implemented and there is no defined adequate
    way to store them.
  */

  proc time_advance(oldSol : [] real, sol : [oldSol.domain] real, res : [oldSol.domain] real, dt : real,
      thisStage : int, timeScheme : int) : [sol.domain] real
  {
    // Advance the solution in time
    // Select and call the selected time marching scheme

    use Parameters.ParamInput;

    var newSol : [sol.domain] real;

    select timeScheme
    {
      when TIME_EULER
      {}
      when TIME_RK_CLASSIC
      {}
      when TIME_TVDRK_O2S2 do
        newSol = sspRK_o2sN(oldSol, sol, res, dt, thisStage, 2);
      when TIME_TVDRK_O2S3 do
        newSol = sspRK_o2sN(oldSol, sol, res, dt, thisStage, 3);
      when TIME_TVDRK_O2S4 do
        newSol = sspRK_o2sN(oldSol, sol, res, dt, thisStage, 4);
      when TIME_TVDRK_O2SN
      {
        // Placeholder hard coded stage number
        var nStages = 5;

        newSol = sspRK_o2sN(oldSol, sol, res, dt, thisStage, nStages);
      }
      when TIME_TVDRK_O3S3 do
        newSol = sspRK_o3s3(oldSol, sol, res, dt, thisStage);
      when TIME_TVDRK_O3S4 {}
      when TIME_TVDRK_O3S5 {}
      when TIME_TVDRK_O4S5 {}
    }

    return newSol;
  }

  proc max_wave_speed(consVars : [] real) : real
  {
    use LinearAlgebra;
    import Flux.sound_speed_cv;
    import Flux.velocity_magnitude_cv;

    // Calculate the maximum wave speed at all points in this cell
    const velocity   : real = velocity_magnitude_cv(consVars);
    const soundSpeed : real = sound_speed_cv(consVars);

    const maxWaveSpeed : real = soundSpeed + velocity;

    return maxWaveSpeed;
  }

  proc max_wave_speed_array(consVarsArray : [] real) : real
  {
    use LinearAlgebra;
    import Flux.sound_speed_cv;
    import Flux.velocity_magnitude_cv;

    var maxWaveSpeed : real = 0.0;

    // Calculate the maximum wave speed at all points in this cell
    for spIdx in consVarsArray.domain.dim(1)
    {
      // Pick maximum wave speed at any SP
      maxWaveSpeed = max( maxWaveSpeed, max_wave_speed(consVarsArray[.., spIdx]));
    }

    return maxWaveSpeed;
  }

  proc adjust_cfl()
  {
    // Dynamic CFL setting mainly for implicit time stepping schemes
    // Proposed features:
    //   0. Constant CFL
    //   1. Linear CFL ramp
    //   2. Exponential CFL ramp
    //   3. Backtrack ramping if residue or solution diverges
  }

  proc sspRK_o2sN(oldSol : [] real, sol : [oldSol.domain] real, res : [oldSol.domain] real, dt : real, thisStage : int,
      nStages : int) : [sol.domain] real
  {
    var newSol : [sol.domain] real;

    if thisStage < nStages then
      newSol = sol - 1.0/(nStages-1.0)*res*dt;
    else if thisStage == nStages then
      newSol = (oldSol + sol*(nStages-1.0) - res*dt)/nStages;

    return newSol;
  }

  proc sspRK_o3s3(oldSol : [] real, sol : [oldSol.domain] real, res : [oldSol.domain] real, dt : real, thisStage : int)
      : [sol.domain] real
  {
    var newSol : [sol.domain] real;

    if thisStage == 1 then
      newSol = sol - res*dt;
    else if thisStage == 2 then
      newSol = (3.0/4.0)*oldSol + (1.0/4.0)*sol - (1.0/4.0)*res*dt;
    else if thisStage == 3 then
      newSol = (1.0/3.0)*oldSol + (2.0/3.0)*sol - (2.0/3.0)*res*dt;

    return newSol;
  }

  proc main()
  {
    use IO.FormattedIO;
    use Testing;

    // Solver a simple random ODE and validate the solution

  }
}
