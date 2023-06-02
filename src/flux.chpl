module Flux
{
  proc pressure(dens : real, temp : real) : real
  {
    import Input.fR;

    var pressure : real = fR*dens*temp;

    return pressure;
  }

  proc pressure_cv(cons : [] real) : real
  {
    import Input.fGamma;
    import LinearAlgebra.dot;

    var idxDens : int   = cons.domain.dim(0).low;           // First element is density
    var idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    var pressure : real = (fGamma-1.0)*(cons[idxEner] - 0.5*dot(cons[idxMomV],cons[idxMomV])/cons[idxDens]);

    return pressure;
  }

  proc temperature(dens : real, pres : real) : real
  {
    import Input.fR;

    var temperature : real = pres/dens/fR;

    return temperature;
  }

  proc temperature_cv(cons : [] real) : real
  {
    import Input.fR;

    var idxDens : int   = cons.domain.dim(0).low;           // First element is density

    var pres : real = pressure_cv(cons);

    var temperature : real = pres/(cons[idxDens]*fR);

    return temperature;
  }

  proc velocity_magnitude_cv(cons : [] real) : real
  {
    import LinearAlgebra.norm;
    import LinearAlgebra.normType;

    var idxDens : int   = cons.domain.dim(0).low;           // First element is density
    var idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    var velM : real = norm(cons[idxMomV], normType.norm2)/cons[idxDens];

    return velM;
  }

  proc entropy_cv(cons : [] real) : real
  {
    import Input.fGamma;

    var idxDens : int   = cons.domain.dim(0).low;          // First element is density

    var pressure : real = pressure_cv(cons);

    var entropy : real = pressure/(cons[idxDens]**fGamma);

    return entropy;
  }

  proc internal_energy_cv(cons : [] real) : real
  {
    import Input.fCv;
    import Input.fR;

    var p : real = pressure_cv(cons);

    var internalEnergy : real = fCv*p/fR;

    return internalEnergy;
  }

  proc enthalpy_cv(cons : [] real) : real
  {
    var enerInt : real = internal_energy_cv(cons);
    var pres    : real = pressure_cv(cons);

    var enthalpy : real = enerInt + pres;

    return enthalpy;
  }

  proc enthalpy_stagnation_cv(cons : [] real) : real
  {
    var idxEner : int   = cons.domain.dim(0).high;         // Last element is energy

    var pres    : real = pressure_cv(cons);

    var enthalpyStagnation : real = cons[idxEner] + pres;

    return enthalpyStagnation;
  }

  proc sound_speed(dens : real, pres : real) : real
  {
    import Input.fGamma;

    var a : real = sqrt(fGamma*pres/dens);

    return a;
  }

  proc sound_speed(temp : real) : real
  {
    import Input.fGamma;
    import Input.fR;

    var a : real = sqrt(fGamma*fR*temp);

    return a;
  }

  proc sound_speed_cv(cons : [] real) : real
  {
    import Input.fGamma;

    var idxDens : int   = cons.domain.dim(0).low;           // First element is density
    var idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    var p : real = pressure_cv(cons);

    return sqrt(fGamma*p/cons[idxDens]);
  }

  proc sound_speed_pv(prim : [] real) : real
  {
    import Input.fGamma;

    var idxDens : int   = prim.domain.dim(0).low;           // First element is density
    var idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxPres : int   = prim.domain.dim(0).high;          // Last element is the pressure

    return sqrt(fGamma*prim[idxPres]/prim[idxDens]);
  }

  proc mach_cv(cons : [] real) : real
  {
    import LinearAlgebra.norm;

    var idxDens : int   = cons.domain.dim(0).low;           // First element is density
    var idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    var velV : [idxMomV] real = cons[idxMomV]/cons[idxDens];

    var mach : real = norm(velV) / sound_speed_cv(cons);

    return mach;
  }

  proc density(temp : real, pres : real) : real
  {
    import Input.fR;

    var dens : real = pres/(temp*fR);

    return dens;
  }

  proc energy_pv(prim : [] real) : real
  {
    import Input.fGamma;
    import LinearAlgebra.dot;

    var idxDens : int   = prim.domain.dim(0).low;           // First element is density
    var idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxPres : int   = prim.domain.dim(0).high;          // Last element is pressure

    var ener : real = prim[idxPres]/(fGamma-1) + 0.5*prim[idxDens]*dot(prim[idxVelV], prim[idxVelV]);

    return ener;
  }

  proc cons2prim(cons : [] real) : [] real
  {
    var idxDens : int   = cons.domain.dim(0).low;           // First element is density
    var idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the momentum vector
    var idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    var idxVelV : range = idxMomV;
    var idxPres : int   = idxEner;

    var prim : [cons.domain] real;

    prim[idxDens] = cons[idxDens];
    prim[idxVelV] = cons[idxMomV]/cons[idxDens];
    prim[idxPres] = pressure_cv(cons);

    return prim;
  }

  proc prim2cons(prim : [] real) : [] real
  {
    var idxDens : int   = prim.domain.dim(0).low;           // First element is density
    var idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocity vector
    var idxPres : int   = prim.domain.dim(0).high;          // Last element is pressure

    var idxMomV : range = idxVelV;
    var idxEner : int   = idxPres;

    var cons : [prim.domain] real;

    cons[idxDens] = prim[idxDens];
    cons[idxMomV] = prim[idxDens]*prim[idxVelV];
    cons[idxEner] = energy_pv(prim);

    return cons;
  }

  proc convection_flux_cv_1d(cons : real) : real
  {
    import Input.convectionSpeed;

    return cons*convectionSpeed[1];
  }

  proc convection_flux_cv(cons : real) : [] real
  {
    import Input.convectionSpeed;
    import Input.nEqs;

    return cons*convectionSpeed[1..nEqs];
  }

  proc burgers_flux_cv_1d(cons : real) : real
  {
    return 0.5*cons**2.0;
  }

  proc euler_flux_cv_1d(cons : [1..3] real) : [1..3] real
  {
    var euler_flux_cv : [cons.domain] real;
    var p : real = pressure_cv(cons);

    euler_flux_cv[1] = cons[2];
    euler_flux_cv[2] = cons[2]*cons[2]/cons[1] + p;
    euler_flux_cv[3] = cons[2]/cons[1]*(cons[3] + p);

    return euler_flux_cv;
  }

  proc euler_flux_cv(cons : [] real) : [] real
  {
    var idxDens : int   = cons.domain.dim(0).low;           // First element is density
    var idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    var euler_flux_cv : [idxMomV-1, cons.domain.dim(0)] real;

    var velV : [idxMomV] real = cons[idxMomV] / cons[idxDens];
    var pres : real = pressure_cv(cons);

    for idxDir in idxMomV-1
    {
      euler_flux_cv[idxDir, idxDens] = cons[idxDir+1];
      euler_flux_cv[idxDir, idxMomV] = cons[idxDir+1]*velV[idxMomV];
      euler_flux_cv[idxDir, idxEner] = velV[idxDir+1]*(cons[idxEner] + pres);

      euler_flux_cv[idxDir, idxDir+1] += pres;
    }

    return euler_flux_cv;
  }

  proc euler_flux_pv(prim : [] real) : [] real
  {
    use Input;
    import LinearAlgebra.dot;

    var idxDens : int   = prim.domain.dim(0).low;           // First element is density
    var idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocity components
    var idxPres : int   = prim.domain.dim(0).high;          // Last element is pressure

    var euler_flux_pv : [idxVelV-1, prim.domain.dim(0)] real;

    for idxDir in idxVelV-1
    {
      euler_flux_pv[idxDir, idxDens] = prim[idxDens+idxDir]*prim[idxDens];
      euler_flux_pv[idxDir, idxVelV] = prim[idxDens+idxDir]*prim[idxDens]*prim[idxVelV];
      euler_flux_pv[idxDir, idxPres] = prim[idxDens+idxDir]*(prim[idxPres]*fGamma/(fGamma-1.0)+0.5*prim[idxDens]*dot(prim[idxVelV],prim[idxVelV]));

      euler_flux_pv[idxDir, idxDir+1] += prim[idxPres];
    }

    return euler_flux_pv;
  }

  proc visc_flux_cv(cons : [] real, consGrad : [] real) : [] real
  {
  }

  proc visc_flux_pv(prim : [] real, primGrad : [] real) : [] real
  {
  }

  proc main()
  {
    // Flow defining variables
    param temp : real = 288.15;
    param pres : real = 101_325.0;
    param mach : real = 0.6;
    param alph : real = 0.03490658503988659153847381536977225426885743777083; // 2.0 degrees in Radians
    param beta : real = 0.00872664625997164788461845384244306356721435944271; // 0.5 degree  in Radians

    // Calculated with 50 digits of precision
    param dens : real =           1.22504581154325263871323932326963696202399523031669; // Pres/(fR*Temp)
    param aSpd : real =         340.28762770443822699513075710728002306623909230587339; // sqrt(fGamma*fR*Temp)
    param velM : real =         204.17257662266293619707845426436801383974345538352403; // Mach*aSpd
    param ener : real =     278_846.40000000000000000000000000000000000000000000000000; // Pres/(fGamma-1) + 0.5*Dens*velM^2
    param entr : real =      76_261.15011769449782489957624473444208242962290561088006; // Pres/(Dens**fGamma)
    param enthStag : real = 380_171.40000000000000000000000000000000000000000000000000; // Ener + Pres
    param enerInt : real =  253_312.50000000000000000000000000000000000000000000000000; // fCv*Pres/fR
    param enth : real =     354_637.50000000000000000000000000000000000000000000000000; // EnerInt + Pres

    const velV1d : [1..1] real = [velM                                                              ];
    const velV2d : [1..2] real = [velM*cos(alph)          , velM*sin(alph)                          ];
    const velV3d : [1..3] real = [velM*cos(alph)*cos(beta), velM*sin(alph)*cos(beta), velM*sin(beta)];

    const cons1d : [1..3] real = [dens, dens*velM                                                                        , ener];
    const cons2d : [1..4] real = [dens, dens*velM*cos(alph)          , dens*velM*sin(alph)                               , ener];
    const cons3d : [1..5] real = [dens, dens*velM*cos(alph)*cos(beta), dens*velM*sin(alph)*cos(beta), dens*velM*sin(beta), ener];

    const prim1d : [1..3] real = [dens, velM                                                              , pres];
    const prim2d : [1..4] real = [dens, velM*cos(alph)          , velM*sin(alph)                          , pres];
    const prim3d : [1..5] real = [dens, velM*cos(alph)*cos(beta), velM*sin(alph)*cos(beta), velM*sin(beta), pres];

    writef("Conserved variables, Density (kg/m^3), XYZ Momentum components (kg*m/s*m^3), Energy (J/m^3):\n");
    writef("  1D - Test Reference: %.8ht\n", cons1d);
    writef("            prim2cons: %.8ht\n", prim2cons(prim1d));
    writef("  2D - Test Reference: %.8ht\n", cons2d);
    writef("            prim2cons: %.8ht\n", prim2cons(prim2d));
    writef("  3D - Test Reference: %.8ht\n", cons3d);
    writef("            prim2cons: %.8ht\n", prim2cons(prim3d));
    writef("\n");
    writef("Primitive Variables, Density (kg/m^3), XYZ Velocity components (m/s), Pressure (Pa):\n");
    writef("  1D - Test Reference: %.8ht\n", prim1d);
    writef("            cons2prim: %.8ht\n", cons2prim(cons1d));
    writef("  2D - Test Reference: %.8ht\n", prim2d);
    writef("            cons2prim: %.8ht\n", cons2prim(cons2d));
    writef("  3D - Test Reference: %.8ht\n", prim3d);
    writef("            cons2prim: %.8ht\n", cons2prim(cons3d));
    writef("\n");
    writef("Pressure (Pa) functions:\n");
    writef("      Test Reference: %19.12er\n", pres);
    writef("  From Dens and Temp: %19.12er\n", pressure(dens, temp));
    writef("  1D - From ConsVars: %19.12er\n", pressure_cv(cons1d));
    writef("  2D - From ConsVars: %19.12er\n", pressure_cv(cons1d));
    writef("  3D - From ConsVars: %19.12er\n", pressure_cv(cons1d));
    writef("\n");
    writef("Temperature (K) functions:\n");
    writef("      Test Reference: %19.12er\n", temp);
    writef("  From Dens and Pres: %19.12er\n", temperature(dens, pres));
    writef("  1D - from ConsVars: %19.12er\n", temperature_cv(cons1d));
    writef("  2D - from ConsVars: %19.12er\n", temperature_cv(cons2d));
    writef("  3D - from ConsVars: %19.12er\n", temperature_cv(cons3d));
    writef("\n");
    writef("Density (kg/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", dens);
    writef("  From Temp and Pres: %19.12er\n", density(temp, pres));
    writef("\n");
    writef("Energy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", ener);
    writef("  1D - From PrimVars: %19.12er\n", energy_pv(prim1d));
    writef("  2D - From PrimVars: %19.12er\n", energy_pv(prim2d));
    writef("  3D - From PrimVars: %19.12er\n", energy_pv(prim3d));
    writef("\n");
    writef("Internal Energy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", enerInt);
    writef("  1D - From ConsVars: %19.12er\n", internal_energy_cv(cons1d));
    writef("  2D - From ConsVars: %19.12er\n", internal_energy_cv(cons2d));
    writef("  3D - From ConsVars: %19.12er\n", internal_energy_cv(cons3d));
    writef("\n");
    writef("Entropy () functions:\n");
    writef("      Test Reference: %19.12er\n", entr);
    writef("  1D - From ConsVars: %19.12er\n", entropy_cv(cons1d));
    writef("  2D - From ConsVars: %19.12er\n", entropy_cv(cons2d));
    writef("  3D - From ConsVars: %19.12er\n", entropy_cv(cons3d));
    writef("\n");
    writef("Enthalpy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", enth);
    writef("  1D - From ConsVars: %19.12er\n", enthalpy_cv(cons1d));
    writef("  2D - From ConsVars: %19.12er\n", enthalpy_cv(cons2d));
    writef("  3D - From ConsVars: %19.12er\n", enthalpy_cv(cons3d));
    writef("\n");
    writef("Stagnation Enthalpy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", enthStag);
    writef("  1D - From ConsVars: %19.12er\n", enthalpy_stagnation_cv(cons1d));
    writef("  2D - From ConsVars: %19.12er\n", enthalpy_stagnation_cv(cons2d));
    writef("  3D - From ConsVars: %19.12er\n", enthalpy_stagnation_cv(cons3d));
    writef("\n");
    writef("Speed of Sound (m/s) and Mach (non-dimensional) functions:\n");
    writef("  Speed of Sound from Dens and Pres: %19.12er,     from Temp: %19.12er\n",
        sound_speed(dens, pres), sound_speed(temp));
    writef("  1D - Speed of Sound from ConsVars: %19.12er, from PrimVars: %19.12er - Mach: %19.12er\n",
        sound_speed_cv(cons1d), sound_speed_pv(prim1d), mach_cv(cons1d));
    writef("  2D - Speed of Sound from ConsVars: %19.12er, from PrimVars: %19.12er - Mach: %19.12er\n",
        sound_speed_cv(cons2d), sound_speed_pv(prim2d), mach_cv(cons2d));
    writef("  3D - Speed of Sound from ConsVars: %19.12er, from PrimVars: %19.12er - Mach: %19.12er\n",
        sound_speed_cv(cons3d), sound_speed_pv(prim3d), mach_cv(cons3d));
    writef("\n");
    writef("Euler flux functions:\n");
    writef("  1D Fluxes:\n");
    writef("    1D  Euler Flux (ConsVars): %.8t\n", euler_flux_cv_1d(cons1d));
    writef("    Gen Euler Flux (ConsVars): %.8t\n", euler_flux_cv(cons1d));
    writef("    Gen Euler Flux (PrimVars): %.8t\n", euler_flux_pv(prim1d));
    writef("\n");
    writef("  2D Fluxes:\n");
    writef("    Gen Euler Flux (ConsVars): %.8t\n", euler_flux_cv(cons2d)[1,..]);
    writef("                               %.8t\n", euler_flux_cv(cons2d)[2,..]);
    writef("    Gen Euler Flux (PrimVars): %.8t\n", euler_flux_pv(prim2d)[1,..]);
    writef("                               %.8t\n", euler_flux_pv(prim2d)[2,..]);
    writef("\n");
    writef("  3D Fluxes:\n");
    writef("    Gen Euler Flux (ConsVars): %.8t\n", euler_flux_cv(cons3d)[1,..]);
    writef("                               %.8t\n", euler_flux_cv(cons3d)[2,..]);
    writef("                               %.8t\n", euler_flux_cv(cons3d)[3,..]);
    writef("    Gen Euler Flux (PrimVars): %.8t\n", euler_flux_pv(prim3d)[1,..]);
    writef("                               %.8t\n", euler_flux_pv(prim3d)[2,..]);
    writef("                               %.8t\n", euler_flux_pv(prim3d)[3,..]);
    //writef("\n");
    //writef("Generic Viscous flux function (conserver vars):\n");
    //writef("  Viscous 1D Flux:\n", visc_flux_cv(cons1d));
    //writef("  Viscous 2D Flux:\n", visc_flux_cv(cons2d));
    //writef("  Viscous 3D Flux:\n", visc_flux_cv(cons3d));
    //writef("\n");
    //writef("Generic Viscous flux function (primitive vars):\n");
    //writef("  Viscous 1D Flux:\n", visc_flux_pv(prim1d));
    //writef("  Viscous 2D Flux:\n", visc_flux_pv(prim2d));
    //writef("  Viscous 3D Flux:\n", visc_flux_pv(prim3d));
  }
}
