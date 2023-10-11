module Ringleb
{
    {
      // Lazy parameters
      //qTop = 0.2
      //kMin = 0.2
      //kMax = 0.4

      // Subsonic parameters
      //qTop = 0.35
      //kMin = 0.40
      //kMax = 0.60

      // Barth parameters
      //qTop = 0.43
      //kMin = 0.60
      //kMax = 0.98

      // Wolf parameters
      //qTop = 0.30
      //kMin = 0.40
      //kMax = 0.80

      // ZJ Wang parameters
      //qTop = 0.50
      //kMin = 0.70
      //kMax = 1.20

      // HO-Workshop parameters
      //qTop = 0.50
      //kMin = 0.70
      //kMax = 1.50
    }

  proc ringleb_sol(xy : [1..2] real) : [] real
  {
    import Input.fGamma;

    var vel   : real = ringleb_xy_to_vel(xy);
    var dens  : real = ringleb_dens(ringleb_c(vel));
    var k     : real = 1.0/sqrt(1.0/(2.0*vel*vel) + dens*(ringleb_j(ringleb_c(vel))/2.0 + xy[1]));
    //var theta : real = arcsin(vel/k);

    var consVars : [1..4] real;

    consVars[1] = dens;
    consVars[2] = dens * vel*sqrt(max(1.0 - (vel/k)**2, 0.0)) * sgn(xy[2]);
    consVars[3] =-dens * vel*min(1.0, vel/k);
    consVars[4] = dens * vel*vel/2.0 + ringleb_pres(ringleb_c(vel))/(fGamma-1);

    return consVars;
  }

  proc ringleb_kq_to_xy(k : real, q : real) : [1..2] real
  {
    import Input.fGamma;

    var gamm1 : real = fGamma-1.0;

    var c   : real = ringleb_c(q);
    var rho : real = ringleb_dens(c);
    var j   : real = ringleb_j(c);

    var xy : [1..2] real;

    xy[1] = (1.0/q**2 - 2.0/k**2)/(2.0*rho) + j/2.0;
    xy[2] = sqrt(1.0-(q/k)**2)/(k*rho*q);
    if isNan(xy[2]) then xy[2] = 0.0;

    return xy;
  }

  proc ringleb_xy_to_vel(xy : [1..2] real) : real
  {
    use Parameters.ParamConstants;

    const maxIter : int = 128;

    // Bisection initial values
    var qLo : real = EPS4;
    var qHi : real = ringleb_k_max() - EPS12;

    var cLo : real = ringleb_c(qLo);
    var cHi : real = ringleb_c(qHi);

    var lhsLo : real = (xy[1] + ringleb_j(cLo)/2.0)**2 + xy[2]**2 - (1/(2*ringleb_dens(cLo)*qLo**2))**2;
    var lhsHi : real = (xy[1] + ringleb_j(cHi)/2.0)**2 + xy[2]**2 - (1/(2*ringleb_dens(cHi)*qHi**2))**2;

    // Iteration variables
    var qNew : real = (qLo+qHi)/2.0;
    var cNew, lhsNew : real;

    for k in 1..maxIter
    {
      // Evaluate LHS
      cNew = ringleb_c(qNew);
      lhsNew = (xy[1] + ringleb_j(cNew)/2.0)**2 + xy[2]**2 - (1/(2*ringleb_dens(cNew)*qNew**2))**2;

      // Update solution bounds
      {
        if sgn(lhsNew) * sgn(lhsLo) == 1
        {
          qLo = qNew;
          lhsLo = lhsNew;
        }
        if sgn(lhsNew) * sgn(lhsHi) == 1
        {
          qHi = qNew;
          lhsHi = lhsNew;
        }
      }

      // Check for convergente
      if qNew == (qLo+qHi)/2.0 then
        return qNew;
      else
        qNew = (qLo+qHi)/2.0;
    }

    writeln();
    writeln("RINGLEB SOLUTION NOT CONVERGED");
    writeln("    xy = ", xy);
    writef("    qLo = %22.15er | lhsLo = %22.15er\n", qLo, lhsLo);
    writef("    qHi = %22.15er | lhsHi = %22.15er\n", qHi, lhsHi);
    writef("    qNew= %22.15er | lhsnew= %22.15er\n", qNew, lhsNew);
    writef("    qNex= %22.15er\n", (qLo+qHi)/2.0);

    return -1.0;
  }

  proc ringleb_c(q : real) : real
  {
    import Input.fGamma;

    var c : real = sqrt(1.0-(fGamma-1.0)*(q**2)/2.0);
    if isNan(c) then c = 0.0;

    return c;
  }

  proc ringleb_dens(c : real) : real
  {
    import Input.fGamma;

    return c**(2/(fGamma-1));
  }

  proc ringleb_pres(c : real) : real
  {
    import Input.fGamma;

    return c**(2*fGamma/(fGamma-1))/fGamma;
  }

  proc ringleb_j(c : real) : real
  {
    import Math.log;
    import Input.fGamma;

    return 1/c + 1/(3*c**3) + 1/(5*c**5) - log((1.0+c)/(1.0-c))/2.0;
  }

  proc ringleb_k_max() : real
  {
    import Input.fGamma;

    return 4.0/(fGamma+1.0);
  }

  proc main()
  {
    var xy : [1..2] real;

    xy = [ 0.2, 0.0];
    writeln("xy = ", xy, ", sol = ", ringleb_sol(xy));
    xy = [ 0.0, 0.0];
    writeln("xy = ", xy, ", sol = ", ringleb_sol(xy));
    xy = [-0.2, 0.0];
    writeln("xy = ", xy, ", sol = ", ringleb_sol(xy));
    xy = [-0.2, 0.5];
    writeln("xy = ", xy, ", sol = ", ringleb_sol(xy));
    xy = [-0.2, 1.0];
    writeln("xy = ", xy, ", sol = ", ringleb_sol(xy));
    xy = [-0.2,-1.0];
    writeln("xy = ", xy, ", sol = ", ringleb_sol(xy));
  }
}
