module Init
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

  proc flow_condition(familySubType : int, familyParameters : [] real, xyz : [] real) : [] real
  {
    import Dimensional.scales;
    use Ringleb;

    var flowVars : [xyz.domain.dim(0), 1..Input.nEqs] real;

    select familySubType
    {
      when IC_SINUSOIDAL
      {
        for i in xyz.domain.dim(0) do
          flowVars[i,1] = sin(xyz[i,1]*pi);
      }
      when IC_GAUSSPULSE
      {
        for i in xyz.domain.dim(0)
        {
          var a : real = 1/4;
          var b : real = 3/4;

          if (xyz[i,1] > a && xyz[i,1] < b) then
            flowVars[i,1] = exp(1/((xyz[i,1]-a)*(xyz[i,1]-b))+(abs(a-b)/2)**(-2.0));
          else
            flowVars[i,1] = 0.0;;
        }
      }
      when IC_SHOCKTUBE
      {
        const xMin : real = 0.0;
        const xMax : real = 1.0;

        ref  rhoHi = familyParameters[1];
        ref  pHi   = familyParameters[2];
        ref  rhoLo = familyParameters[3];
        ref  pLo   = familyParameters[4];

        for i in xyz.domain.dim(0)
        {
          if xyz[i,1] < 0.5*(xMin+xMax) then
          {
            flowVars[i,1] = rhoHi;
            flowVars[i,2] = 0.0;
            flowVars[i,3] = pHi/(fGamma-1.0);
          }
          else
          {
            flowVars[i,1] = rhoLo;
            flowVars[i,2] = 0.0;
            flowVars[i,3] = pLo/(fGamma-1.0);
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

          flowVars[i,1] = dens;
          flowVars[i,2] = dens*vel;
          flowVars[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
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

          flowVars[i,1] = dens;
          flowVars[i,2] = dens*vel;
          flowVars[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
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

          flowVars[i,1] = dens;
          flowVars[i,2] = dens*vel;
          flowVars[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
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

              flowVars[i,1] = dens;
              flowVars[i,2] = dens*vel;
              flowVars[i,3] = pres/(fGamma-1) + 0.5*dens*vel**2;
            }
            else
              // If the shock fits then keep the ideal solution on the rest of the nozzle
              break;
          }
        }
      }
      when IC_RINGLEB
      {
        for i in xyz.domain.dim(0) do
          flowVars[i,1..4] = scales!.dim2non_cv(ringleb_sol(xyz[i,1..2]));
      }
      when IC_GENERIC_MEANFLOW
      {
        for ptIdx in xyz.domain.dim(0) do
          flowVars[ptIdx, 1..Input.nEqs] = scales!.dim2non_cv(familyParameters[1..Input.nEqs]);
      }
    }

    return flowVars;
  }

  proc nozzle_area(x : real) : real
  {
    // Nozzle shape from "i do like cfd", vol 1 by Katate Masatsuka (Hiroaki Nishikawa), page 248
    // a(x) = 25/9(ae-at)(x-xt)^2 + at
    // ae = exit area
    // at = throat area
    // xt = throat location

    return 25.0/9.0*(exitArea-throatArea)*(x-xThroat)**2 + throatArea;
  }

  proc nozzle_area_deriv(x : real) : real
  {
    // Nozzle shape from "i do like cfd", vol 1 by Katate Masatsuka (Hiroaki Nishikawa), page 248
    // a(x) = 25/9(ae-at)(x-xt)^2 + at
    // a'(x) = 50/9(ae-at)(x-xt)
    // ae = exit area
    // at = throat area
    // xt = throat location

    return 50.0/9.0*(exitArea-throatArea)*(x-xThroat);
  }

  proc nozzle_mach_sub(areaRatio : real) : real
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

  proc nozzle_mach_sup(areaRatio : real) : real
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

  proc freestream(familyParameters : [] real) : real
  {
    // Calculate the conserved variables given:
    // - Mach
    // - Alpha (Angle of Attack)
    // - Beta (Sideslip)
    // - Pressure
    // - Temperature
  }

  proc main()
  {
    use IO;
    use Parameters.ParamInput;
    use Flux;

    const nNodes = 21;
    param fGamma = 1.4;
    var xyz : [1..nNodes, 1..3] real;
    var sol : [1..nNodes, 1..3] real;

    writeln();
    writeln("1D Meshes");
    {
      // Generate a 1D mesh
      for i in xyz.domain.dim(0)
      {
        xyz[i, 1] = (i-1.0)/(nNodes-1.0);
        xyz[i, 2] = 0.0;
        xyz[i, 3] = 0.0;
      }

      writeln();
      writeln("1D Shock Tube");
      Input.nEqs = 3;
      var famlParameters : [1..4] real = [1.0, 1.0, 0.1, 0.1];
      sol = flow_condition(IC_SHOCKTUBE, famlParameters, xyz);
      writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3]");
      for i in xyz.domain.dim(0) do
        writef("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er\n", i, xyz[i, 1], xyz[i, 2], xyz[i, 3],
            sol[i, 1], sol[i, 2], sol[i, 3]);

      writeln();
      writeln("1D Nozzle Flow - Subsonic");
      Input.nEqs = 3;
      sol = flow_condition(IC_1D_NOZZLE_SUBSONIC, [0.0], xyz);
      writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3],   Pressure,       Mach");
      for i in xyz.domain.dim(0) do
        writef("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er\n", i, xyz[i, 1], xyz[i, 2],
            xyz[i, 3], sol[i, 1], sol[i, 2], sol[i, 3], pressure_cv(sol[i,..], fGamma), mach_cv(sol[i,..], fGamma));

      writeln();
      writeln("1D Nozzle Flow - Smooth Transonic");
      Input.nEqs = 3;
      sol = flow_condition(IC_1D_NOZZLE_SMOOTH_TRANSONIC, [0.0], xyz);
      writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3],   Pressure,       Mach");
      for i in xyz.domain.dim(0) do
        writef("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er\n", i, xyz[i, 1], xyz[i, 2],
            xyz[i, 3], sol[i, 1], sol[i, 2], sol[i, 3], pressure_cv(sol[i,..], fGamma), mach_cv(sol[i,..], fGamma));

      writeln();
      writeln("1D Nozzle Flow - Shocked Transonic");
      Input.nEqs = 3;
      sol = flow_condition(IC_1D_NOZZLE_SHOCKED_TRANSONIC, [0.0], xyz);
      writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3],   Pressure,       Mach");
      for i in xyz.domain.dim(0) do
        writef("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er\n", i, xyz[i, 1], xyz[i, 2],
            xyz[i, 3], sol[i, 1], sol[i, 2], sol[i, 3], pressure_cv(sol[i,..], fGamma), mach_cv(sol[i,..], fGamma));
    }

    writeln();
    writeln("2D Meshes");
    {
      // Generate Ringleb Mesh
    }
  }
}
