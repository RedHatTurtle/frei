prototype module SourceTerm
{
  use Random;
  use UnitTest;
  use Parameters;

  proc source_term(xyz : [] real, consVars : [] real, eqSet : int) : [consVars.domain] real
  {
    use Parameters.ParamInput;
    import Init.nozzle_area;

    var sourceTerm : [consVars.domain] real;

    select eqSet
    {
      when EQ_QUASI_1D_EULER
      {
        for point in consVars.domain.dim(0)
        {
          // There could be a switch for different area functions if there are any. For this a new input variable might
          // be needed.
          var area  : real = nozzle_area(xyz[point, 1]);
          var dArea : real = nozzle_area_deriv(xyz[point, 1]);

          sourceTerm[point, 1..3] = quasi_1d_euler(consVars[point, 1..3], area, dArea);
        }
      }
    }

    return sourceTerm;
  }

  proc quasi_1d_euler(consVars : [1..3] real, area : real, dArea : real) : [1..3] real
  {
    import Flux.pressure_cv;

    var sourceTerm : [1..3] real;
    var p : real = pressure_cv(consVars);

    sourceTerm[1] = consVars[2];
    sourceTerm[2] = consVars[2]*consVars[2]/consVars[1];
    sourceTerm[3] = consVars[2]/consVars[1]*(consVars[3] + p);

    sourceTerm = - (dArea/area) * sourceTerm;

    return sourceTerm;
  }
}
