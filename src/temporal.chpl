prototype module Temporal_Methods
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

  proc time_step(cfl : real, consVars : [] real) : real
  {
    import Flux.pressure_cv;
    import Flux.temperature_cv;
    import Flux.sound_speed_cv;

    // Calculate time step size for variable time stepping / constant CFL time marching

    var idxRho : int   = consVars.domain.dim(0).low;           // First element is density
    var idxMom : range = consVars.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = consVars.domain.dim(0).high;          // Last element is energy

    var velocity   : real = norm(consVars[idxMom]) / consVars[idxRho];
    var soundSpeed : real = sound_speed_cv(consVars);

    var maxEigenvalue = soundSpeed + velocity;

    // Placeholder values
    var clength = 1;
    var timeStep : real;

    //if ( isViscous == 1 ) then
    //    getTimeStep = min( clength**2 * reynolds * consVars[inxRho] * dt / viscosity_cv( consVars ), clength * h /
    //        maxEigenvalue);
    //else
        timeStep = cfl * clength / maxEigenvalue;

    return timeStep;
  }

  proc adjust_cfl()
  {
    // Dynamic CFL setting mainly for implicit time stepping schemes
    // Proposed features:
    //   1. Linear CFL ramp
    //   2. Exponential CFL rampo
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
