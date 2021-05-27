prototype module Init
{
  use Parameters.ParamInput;
  use Input;
  use Config;

  // Nozzle shape from "I do Like CFD", vol 1 by Katate Masatsuka (Hiroaki Nishikawa), page 248
  // A(x) = 25/9(Ae-At)(x-xT)^2 + At
  // aE = Exit area
  // aT = Throat area
  // xT = Throat location
  param xThroat    : real = 0.4;
  param exitArea   : real = 0.4;
  param throatArea : real = 0.2;

  proc initial_condition(IC : int, const ref xyz : [] real) : [] real
  {
    var sol : [xyz.domain.dim(0), 1..3] real;
    select IC
    {
      when IC_SHOCKTUBE
      {
        var xMin : real = 0.0;
        var xMax : real = 1.0;

        var rhoHi : real = 2.0;
        var eHi   : real = 2.0;
        var rhoLo : real = 1.0;
        var eLo   : real = 1.0;

        for i in xyz.domain.dim(0)
        {
          if xyz[i,1] < 0.5*(xMin+xMax) then
          {
            sol[i,1] = rhoHi;
            sol[i,2] = 0.0;
            sol[i,3] = eHi;
          }
          else
          {
            sol[i,1] = rhoLo;
            sol[i,2] = 0.0;
            sol[i,3] = eLo;
          }
        }
      }
      when IC_1D_NOZZLE_SUBSONIC
      {
        const dens0 : real = 1.0;
        const pres0 : real = 1.0;
        const temp0 : real = 1.0;

        // Set the exit pressure to uniquely define the solution
        const presExit : real = 0.97;
        const machExit : real = sqrt(2/(fGamma-1)*((presExit/pres0)**(-(fGamma-1)/fGamma)-1));

        const aux : real = ( (2+(fGamma-1)*machExit**2)/(fGamma+1) )**((fGamma+1)/(fGamma-1)) / machExit**2;

        var aCrit : real = exitArea / sqrt(aux);

        for i in xyz.domain.dim(0)
        {
          var area : real = nozzle_area(xyz[i,1]);

          var mach : real = nozzle_mach_sub(area/aCrit);

          var aux : real = 1.0 + 0.5*(fGamma-1.0)*mach**2;

          var dens : real = dens0 * aux**(-1/(fGamma-1));
          var pres : real = pres0 * aux**(-fGamma/(fGamma-1));
          var temp : real = temp0 * aux**(-1);

          var a : real = sqrt(fGamma*pres/dens);
          var vel : real = a*mach;

          sol[i,1] = dens;
          sol[i,2] = dens*vel;
          sol[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
        }
      }
      when IC_1D_NOZZLE_SMOOTH_TRANSONIC
      {
        const dens0 : real = 1.0;
        const pres0 : real = 1.0;
        const temp0 : real = 1.0;

        var aCrit : real = nozzle_area(xThroat);

        for i in xyz.domain.dim(0)
        {
          var mach : real;
          var area = nozzle_area(xyz[i,1]);

          if xyz[i,1] < xThroat then
            // Solve for subsonic conditions
            mach = nozzle_mach_sub(area/throatArea);
          else
            // Solve for supersonic conditions
            mach = nozzle_mach_sup(area/throatArea);

          var aux : real = 1.0 + 0.5*(fGamma-1.0)*mach**2;

          var dens : real = dens0 * aux**(-1/(fGamma-1));
          var pres : real = pres0 * aux**(-fGamma/(fGamma-1));
          var temp : real = temp0 * aux**(-1);

          var a : real = sqrt(fGamma*pres/dens);
          var vel : real = a*mach;

          sol[i,1] = dens;
          sol[i,2] = dens*vel;
          sol[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
        }
      }
      when IC_1D_NOZZLE_SHOCKED_TRANSONIC
      {
        const dens0 : real = 1.0;
        const pres0 : real = 1.0;
        const temp0 : real = 1.0;

        // Start with calculating the ideal solution
        for i in xyz.domain.dim(0)
        {
          var mach : real;
          var area = nozzle_area(xyz[i,1]);

          if xyz[i,1] < xThroat then
            // Solve for subsonic conditions
            mach = nozzle_mach_sub(area/throatArea);
          else
            // Solve for supersonic conditions
            mach = nozzle_mach_sup(area/throatArea);

          var aux : real = 1.0 + 0.5*(fGamma-1.0)*mach**2;

          var dens : real = dens0 * aux**(-1/(fGamma-1));
          var pres : real = pres0 * aux**(-fGamma/(fGamma-1));
          var temp : real = temp0 * aux**(-1);

          var a : real = sqrt(fGamma*pres/dens);
          var vel : real = a*mach;

          sol[i,1] = dens;
          sol[i,2] = dens*vel;
          sol[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
        }

        // Now calculate the subsonic solution based on the exit pressure
        {
          param presExit : real = 0.9;
          const machExit : real = sqrt(2/(fGamma-1)*((presExit/pres0)**(-(fGamma-1)/fGamma)-1));
          const aux : real = ( (2+(fGamma-1)*machExit**2)/(fGamma+1) )**((fGamma+1)/(fGamma-1)) / machExit**2;
          const areaExitCrit : real = exitArea / sqrt(aux);

          for i in xyz.domain.dim(0) by -1
          {
            var area : real = nozzle_area(xyz[i,1]);

            // Solve for supersonic conditions
            var machIdeal : real = nozzle_mach_sup(area/throatArea);
            // Calculate mach after a hipothetical shock
            var machPostShock : real = sqrt( (1+0.5*(fGamma-1)*(machIdeal)**2)/(fGamma*machIdeal**2-0.5*(fGamma-1)) );
            // Calculate thesubsonic mach from the exit pressure
            var machSub : real = nozzle_mach_sub(area/areaExitCrit);

            // Check if a shock links the supersonic solution to the subsonic one based on the exit pressure
            if machPostShock > machSub
            {
              // If it doesn´t then propagate the subsonic solution backwards
              var aux : real = 1.0 + 0.5*(fGamma-1.0)*machSub**2;

              var dens : real = dens0 * aux**(-1/(fGamma-1));
              var pres : real = pres0 * aux**(-fGamma/(fGamma-1));
              var temp : real = temp0 * aux**(-1);

              var a : real = sqrt(fGamma*pres/dens);
              var vel : real = a*machSub;

              sol[i,1] = dens;
              sol[i,2] = dens*vel;
              sol[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
            }
            else
              // If the shock fits then keep the ideal solution on the rest of the nozzle
              break;
          }
        }
      }
    }

    return sol;
  }

  proc nozzle_area(in x : real) : real
  {
    // Nozzle shape from "I do Like CFD", vol 1 by Katate Masatsuka (Hiroaki Nishikawa), page 248
    // A(x) = 25/9(Ae-At)(x-xT)^2 + At
    // aE = Exit area
    // aT = Throat area
    // xT = Throat location

    return 25.0/9.0*(exitArea-throatArea)*(x-xThroat)**2 + throatArea;
  }

  proc nozzle_mach_sub(in areaRatio : real) : real
  {
    var machA : real = 1.0;
    var machB : real = 0.0;
    var mach  : real = 0.5*(machA + machB);

    while (mach != machA) && (mach != machB)
    {
      var res = (1.0/mach**2)*((2.0+(fGamma-1.0)*mach**2)/(fGamma+1.0))**((fGamma+1.0)/(fGamma-1.0));

      if res > areaRatio**2 then
        machB = mach;
      else
        machA = mach;

      mach = 0.5*(machA + machB);
    }

    return mach;
  }

  proc nozzle_mach_sup(in areaRatio : real) : real
  {
    var machA : real =   1.0;
    var machB : real = 100.0;
    var mach  : real = 0.5*(machA + machB);

    while (mach != machA) && (mach != machB)
    {
      var res = (1.0/mach**2)*((2.0+(fGamma-1.0)*mach**2)/(fGamma+1.0))**((fGamma+1.0)/(fGamma-1.0));

      if res > areaRatio**2 then
        machB = mach;
      else
        machA = mach;

      mach = 0.5*(machA + machB);
    }

    return mach;
  }

  proc entropy_wave_1d(x : real) : real
  {}

  proc main()
  {
    use IO;
    use Parameters.ParamInput;
    use Flux;

    const nNodes = 21;
    var xyz : [1..nNodes, 1..3] real;
    var sol : [1..nNodes, 1..3] real;

    writeln("1D Meshes");

    for i in xyz.domain.dim(0)
    {
      xyz[i, 1] = (i-1.0)/(nNodes-1.0);
      xyz[i, 2] = 0.0;
      xyz[i, 3] = 0.0;
    }

    writeln();
    writeln("1D Shock Tube");
    sol = initial_condition(IC_SHOCKTUBE, xyz);
    writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3]");
    for i in xyz.domain.dim(0) do
      writeln("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er".format(i, xyz[i, 1], xyz[i, 2], xyz[i, 3],
                                                                                    sol[i, 1], sol[i, 2], sol[i, 3]));

    writeln();
    writeln("1D Nozzle Flow - Subsonic");
    sol = initial_condition(IC_1D_NOZZLE_SUBSONIC, xyz);
    writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3],   Pressure,       Mach");
    for i in xyz.domain.dim(0) do
      writeln("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er".format(i, xyz[i, 1],
            xyz[i, 2], xyz[i, 3], sol[i, 1], sol[i, 2], sol[i, 3], pressure_cv(sol[i,..]), mach_cv(sol[i,..])));

    writeln();
    writeln("1D Nozzle Flow - Smooth Transonic");
    sol = initial_condition(IC_1D_NOZZLE_SMOOTH_TRANSONIC, xyz);
    writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3],   Pressure,       Mach");
    for i in xyz.domain.dim(0) do
      writeln("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er".format(i, xyz[i, 1],
            xyz[i, 2], xyz[i, 3], sol[i, 1], sol[i, 2], sol[i, 3], pressure_cv(sol[i,..]), mach_cv(sol[i,..])));

    writeln();
    writeln("1D Nozzle Flow - Shocked Transonic");
    sol = initial_condition(IC_1D_NOZZLE_SHOCKED_TRANSONIC, xyz);
    writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3],   Pressure,       Mach");
    for i in xyz.domain.dim(0) do
      writeln("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er".format(i, xyz[i, 1],
            xyz[i, 2], xyz[i, 3], sol[i, 1], sol[i, 2], sol[i, 3], pressure_cv(sol[i,..]), mach_cv(sol[i,..])));
  }
}
