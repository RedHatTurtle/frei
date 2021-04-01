prototype module FR
{
  use Parameters.Input;

  const fGamma : real = 1.4;

  proc initial_conditions(IC : int, xyz : real, sol : real)
  {
    select IC
    {
      when IC_SHOCKTUBE
      {
        for i in xyz.dom(0)
        {
          if xyz[i,1] < 0.4 then
          {
            sol[i,1] = 1.0;
            sol[i,2] = 1.0;
            sol[i,3] = 1.0;
          }
          else
          {
            sol[i,1] = 1.0;
            sol[i,2] = 1.0;
            sol[i,3] = 1.0;
          }
        }
      }
      when IC_1D_NOZZLE
      {
        const xT : real = 0.4;
        const aE : real = 0.4;
        const aT : real = 0.2;
        var areaCrit : real = nozzle_area(xT); // Smooth transonic flow = Throat Area

        for i in xyz.dom(0)
        {
          const rho0 = 1;
          const p0 = 1;
          const t0 = 1;

          var area = nozzle_ratio(xyz[i,1]);
          var mach = nozzle_mach(area);

          var aux = 1 + (fGamma-1.0)*mach**2/2;

          var rho = rho0 * aux**(-1/(fGamma-1));
          var p = p0 * aux**(-fGamma/(fGamma-1));
          var t = t0 * aux**(-1);

          var a = sqrt(fGamma*p/rho);

          sol[i,1] = rho;  // Density
          sol[i,2] = rho*mach*a;  // Momentum
          sol[i,3] = p*(1/(fGamma-1) + fGamma*mach**2);  // Energy
        }
      }
      otherwise {}
    }
  }

  proc nozzle_ratio(in x : real) : real
  {
    // Nozzle shape (from "I do Like CFD", vol 1 by Katate Masatsuka (Hiroaki Nishikawa)_)
    // A(x) = 25/9(Ae-At)(x-0.4)^2 + At
    // aE = Exit area
    // aT = Throat area
    const xT : real = 0.4;
    const aE : real = 0.4;
    const aT : real = 0.2;

    return (25.0/9.0*(aE-aT)*(x-xT)**2 + aT)/aT;
  }

  proc nozzle_mach_sub(in areaRatio : real) : real
  {
    use IO;

    var machLo : real = 0.0;
    var machHi : real = 1.0;
    var mach : real = 0.5*(machLo + machHi);

    while (mach != machLo) && (mach != machHi)
    {
      var res = (1.0/mach)*((2.0+(fGamma-1.0)*mach**2)/(fGamma+1.0))**((fGamma+1.0)/(fGamma-1.0));

      if res > areaRatio**2 then
        machLo = mach;
      else
        machHi = mach;

      mach = 0.5*(machLo + machHi);
    }

    return mach;
  }

  proc nozzle_mach_sup(in areaRatio : real) : real
  {
    var machLo : real =   1.0;
    var machHi : real = 100.0;
    var mach : real = 0.5*(machLo + machHi);

    while (mach != machLo) && (mach != machHi)
    {
      var res = (1.0/mach)*((2.0+(fGamma-1.0)*mach**2)/(fGamma+1.0))**((fGamma+1.0)/(fGamma-1.0));

      if res < areaRatio**2 then
        machLo = mach;
      else
        machHi = mach;

      mach = 0.5*(machLo + machHi);
    }

    return mach;
  }

  proc main()
  {
    use IO;

    const iMax = 20;
    var x : real;

    writeln("1D Nozzle Flow");
    writeln("Point #, X-Coord, Area Ratio,    Mach");
    for i in 0..iMax
    {
      x = i/iMax:real;

      if (x < 0.4) then
        writeln("%7i, %7.4dr, %10.4dr, %7.4dr".format( i, x, nozzle_ratio(x), nozzle_mach_sub(nozzle_ratio(x)) ));
      else
        writeln("%7i, %7.4dr, %10.4dr, %7.4dr".format( i, x, nozzle_ratio(x), nozzle_mach_sup(nozzle_ratio(x)) ));
    }
  }
}
