module Flux
{
  proc pressure(dens : real, temp : real, fR : real) : real
  {
    const pres : real = fR*dens*temp;

    return pres;
  }

  proc pressure_cv(const cons : [] real, fGamma : real) : real
  {
    import LinearAlgebra.dot;

    const idxDens : int   = cons.domain.dim(0).low;           // First element is density
    const idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    const pres : real = (fGamma-1.0)*(cons[idxEner] - 0.5*dot(cons[idxMomV],cons[idxMomV])/cons[idxDens]);

    return pres;
  }

  proc temperature(dens : real, pres : real, fR : real) : real
  {
    const temp : real = pres/dens/fR;

    return temp;
  }

  proc temperature_cv(const cons : [] real, fGamma : real, fR : real) : real
  {
    const idxDens : int   = cons.domain.dim(0).low;           // First element is density

    const pres : real = pressure_cv(cons, fGamma);
    const temp : real = pres/cons[idxDens]/fR;

    return temp;
  }

  proc velocity_magnitude_cv(const cons : [] real) : real
  {
    import LinearAlgebra.norm;
    import LinearAlgebra.normType;

    const idxDens : int   = cons.domain.dim(0).low;           // First element is density
    const idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    const velM : real = norm(cons[idxMomV], normType.norm2)/cons[idxDens];

    return velM;
  }

  proc entropy_cv(const cons : [] real, fGamma : real) : real
  {
    const idxDens : int   = cons.domain.dim(0).low;          // First element is density

    const pres : real = pressure_cv(cons, fGamma);
    const entr : real = pres/(cons[idxDens]**fGamma);

    return entr;
  }

  proc internal_energy_cv(const cons : [] real, fGamma : real, fR : real, fCv : real) : real
  {
    const pres    : real = pressure_cv(cons, fGamma);
    const enerInt : real = fCv*pres/fR;

    return enerInt;
  }

  proc enthalpy_cv(const cons : [] real, fGamma : real, fR : real, fCv : real) : real
  {
    const enerInt : real = internal_energy_cv(cons, fGamma, fR, fCv);
    const pres    : real = pressure_cv(cons, fGamma);
    const enth    : real = enerInt + pres;

    return enth;
  }

  proc enthalpy_stagnation_cv(const cons : [] real, fGamma : real) : real
  {
    const idxEner : int   = cons.domain.dim(0).high;         // Last element is energy

    const pres     : real = pressure_cv(cons, fGamma);
    const enthStag : real = cons[idxEner] + pres;

    return enthStag;
  }

  proc sound_speed(dens : real, pres : real, fGamma : real) : real
  {
    const aSpd : real = sqrt(fGamma*pres/dens);

    return aSpd;
  }

  proc sound_speed_temp(temp : real, fGamma : real, fR : real) : real
  {
    const aSpd : real = sqrt(fGamma*fR*temp);

    return aSpd;
  }

  proc sound_speed_cv(const cons : [] real, fGamma : real) : real
  {
    const idxDens : int   = cons.domain.dim(0).low;           // First element is density
    const idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    const pres : real = pressure_cv(cons, fGamma);

    return sqrt(fGamma*pres/cons[idxDens]);
  }

  proc sound_speed_pv(const prim : [] real, fGamma : real) : real
  {
    const idxDens : int   = prim.domain.dim(0).low;           // First element is density
    const idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxPres : int   = prim.domain.dim(0).high;          // Last element is the pressure

    return sqrt(fGamma*prim[idxPres]/prim[idxDens]);
  }

  proc mach_cv(const cons : [] real, fGamma : real) : real
  {
    import LinearAlgebra.norm;

    const idxDens : int   = cons.domain.dim(0).low;           // First element is density
    const idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    const velV : [idxMomV] real = cons[idxMomV]/cons[idxDens];
    const mach : real = norm(velV) / sound_speed_cv(cons, fGamma);

    return mach;
  }

  proc density(temp : real, pres : real, fR : real) : real
  {
    const dens : real = pres/temp/fR;

    return dens;
  }

  proc energy_pv(const prim : [] real, fGamma : real) : real
  {
    import LinearAlgebra.dot;

    const idxDens : int   = prim.domain.dim(0).low;           // First element is density
    const idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxPres : int   = prim.domain.dim(0).high;          // Last element is pressure

    const ener : real = prim[idxPres]/(fGamma-1) + 0.5*prim[idxDens]*dot(prim[idxVelV], prim[idxVelV]);

    return ener;
  }

  proc cons2prim(const cons : [] real, fGamma : real) : [] real
  {
    const idxDens : int   = cons.domain.dim(0).low;           // First element is density
    const idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the momentum vector
    const idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    const idxVelV : range = idxMomV;
    const idxPres : int   = idxEner;

    var prim : [cons.domain] real;

    prim[idxDens] = cons[idxDens];
    prim[idxVelV] = cons[idxMomV]/cons[idxDens];
    prim[idxPres] = pressure_cv(cons, fGamma);

    return prim;
  }

  proc prim2cons(const prim : [] real, fGamma : real) : [] real
  {
    const idxDens : int   = prim.domain.dim(0).low;           // First element is density
    const idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocity vector
    const idxPres : int   = prim.domain.dim(0).high;          // Last element is pressure

    const idxMomV : range = idxVelV;
    const idxEner : int   = idxPres;

    var cons : [prim.domain] real;

    cons[idxDens] = prim[idxDens];
    cons[idxMomV] = prim[idxDens]*prim[idxVelV];
    cons[idxEner] = energy_pv(prim, fGamma);

    return cons;
  }

  proc convection_flux_cv_1d(cons : real, convVelV : real) : real
  {
    return cons*convVelV;
  }

  proc convection_flux_cv(cons : real, convVelV : real) : [] real
  {
    return cons*convVelV;
  }

  proc burgers_flux_cv_1d(cons : real) : real
  {
    return 0.5*cons**2.0;
  }

  proc euler_flux_cv_1d(const cons : [1..3] real, fGamma : real) : [1..3] real
  {
    const pres : real = pressure_cv(cons, fGamma);

    var euler_flux_cv : [cons.domain] real;

    euler_flux_cv[1] = cons[2];
    euler_flux_cv[2] = cons[2]*cons[2]/cons[1] + pres;
    euler_flux_cv[3] = cons[2]/cons[1]*(cons[3] + pres);

    return euler_flux_cv;
  }

  proc euler_flux_cv(const cons : [] real, fGamma : real) : [] real
  {
    const idxDens : int   = cons.domain.dim(0).low;           // First element is density
    const idxMomV : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxEner : int   = cons.domain.dim(0).high;          // Last element is energy

    var euler_flux_cv : [idxMomV-1, cons.domain.dim(0)] real;

    const velV : [idxMomV] real = cons[idxMomV] / cons[idxDens];
    const pres : real = pressure_cv(cons, fGamma);

    for idxDir in idxMomV-1
    {
      euler_flux_cv[idxDir, idxDens] = cons[idxDir+1];
      euler_flux_cv[idxDir, idxMomV] = cons[idxDir+1]*velV[idxMomV];
      euler_flux_cv[idxDir, idxEner] = velV[idxDir+1]*(cons[idxEner] + pres);

      euler_flux_cv[idxDir, idxDir+1] += pres;
    }

    return euler_flux_cv;
  }

  proc euler_flux_pv(const prim : [] real, fGamma : real) : [] real
  {
    import LinearAlgebra.dot;

    const idxDens : int   = prim.domain.dim(0).low;           // First element is density
    const idxVelV : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocity components
    const idxPres : int   = prim.domain.dim(0).high;          // Last element is pressure

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

  proc visc_flux_cv(const cons : [] real, const consGrad : [] real) : [] real
  {
  }

  proc visc_flux_pv(const prim : [] real, const primGrad : [] real) : [] real
  {
  }

  proc main()
  {
    // Fluid properties
    param gasR     : real = 8.314462618;            // Ideal Gas Constant (2018 CODATA), J/(mol*K)
    param fMolMass : real = 2.8966e-2;              // Molar Mass of the fluid mixture, kg/mol
    param fGamma   : real = 1.4;                    // Ideal diatomic gas specific heat ratio
    param fVisc    : real = 1.78792e-5;             // Dynamic viscosity coefficient, kg/(m*s)
    param fKappa   : real = 2.52992e-2;             // Thermal Conductivity, W/(m*K)
    param fR       : real = gasR/fMolMass;          // Calculate the specific gas constant for this fluid, J/(kg*K)
    param fCp      : real = fR/(1.0-1.0/fGamma);    // Calculate the heat coefficients at constant pressure Cp, J/(kg*K)
    param fCv      : real = fR/(fGamma-1.0);        // Calculate the heat coefficients at constant volume Cv, J/(kg*K)
    param fPrandtl : real = fVisc*fCp/fKappa;       // Calculate the Prandtl number of the fluid

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
    param enth : real =     354_637.50000000000000000000000000000000000000000000000000; // EnerInt + Pres
    param enthStag : real = 380_171.40000000000000000000000000000000000000000000000000; // Ener + Pres
    param enerInt : real =  253_312.50000000000000000000000000000000000000000000000000; // fCv*Pres/fR

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
    writef("            prim2cons: %.8ht\n", prim2cons(prim1d, fGamma));
    writef("  2D - Test Reference: %.8ht\n", cons2d);
    writef("            prim2cons: %.8ht\n", prim2cons(prim2d, fGamma));
    writef("  3D - Test Reference: %.8ht\n", cons3d);
    writef("            prim2cons: %.8ht\n", prim2cons(prim3d, fGamma));
    writef("\n");
    writef("Primitive Variables, Density (kg/m^3), XYZ Velocity components (m/s), Pressure (Pa):\n");
    writef("  1D - Test Reference: %.8ht\n", prim1d);
    writef("            cons2prim: %.8ht\n", cons2prim(cons1d, fGamma));
    writef("  2D - Test Reference: %.8ht\n", prim2d);
    writef("            cons2prim: %.8ht\n", cons2prim(cons2d, fGamma));
    writef("  3D - Test Reference: %.8ht\n", prim3d);
    writef("            cons2prim: %.8ht\n", cons2prim(cons3d, fGamma));
    writef("\n");
    writef("Pressure (Pa) functions:\n");
    writef("      Test Reference: %19.12er\n", pres);
    writef("  From Dens and Temp: %19.12er\n", pressure(dens, temp, fR));
    writef("  1D - From ConsVars: %19.12er\n", pressure_cv(cons1d, fGamma));
    writef("  2D - From ConsVars: %19.12er\n", pressure_cv(cons1d, fGamma));
    writef("  3D - From ConsVars: %19.12er\n", pressure_cv(cons1d, fGamma));
    writef("\n");
    writef("Temperature (K) functions:\n");
    writef("      Test Reference: %19.12er\n", temp);
    writef("  From Dens and Pres: %19.12er\n", temperature(dens, pres, fR));
    writef("  1D - from ConsVars: %19.12er\n", temperature_cv(cons1d, fGamma, fR));
    writef("  2D - from ConsVars: %19.12er\n", temperature_cv(cons2d, fGamma, fR));
    writef("  3D - from ConsVars: %19.12er\n", temperature_cv(cons3d, fGamma, fR));
    writef("\n");
    writef("Density (kg/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", dens);
    writef("  From Temp and Pres: %19.12er\n", density(temp, pres, fR));
    writef("\n");
    writef("Energy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", ener);
    writef("  1D - From PrimVars: %19.12er\n", energy_pv(prim1d, fGamma));
    writef("  2D - From PrimVars: %19.12er\n", energy_pv(prim2d, fGamma));
    writef("  3D - From PrimVars: %19.12er\n", energy_pv(prim3d, fGamma));
    writef("\n");
    writef("Internal Energy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", enerInt);
    writef("  1D - From ConsVars: %19.12er\n", internal_energy_cv(cons1d, fGamma, fR, fCv));
    writef("  2D - From ConsVars: %19.12er\n", internal_energy_cv(cons2d, fGamma, fR, fCv));
    writef("  3D - From ConsVars: %19.12er\n", internal_energy_cv(cons3d, fGamma, fR, fCv));
    writef("\n");
    writef("Entropy () functions:\n");
    writef("      Test Reference: %19.12er\n", entr);
    writef("  1D - From ConsVars: %19.12er\n", entropy_cv(cons1d, fGamma));
    writef("  2D - From ConsVars: %19.12er\n", entropy_cv(cons2d, fGamma));
    writef("  3D - From ConsVars: %19.12er\n", entropy_cv(cons3d, fGamma));
    writef("\n");
    writef("Enthalpy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", enth);
    writef("  1D - From ConsVars: %19.12er\n", enthalpy_cv(cons1d, fGamma, fR, fCv));
    writef("  2D - From ConsVars: %19.12er\n", enthalpy_cv(cons2d, fGamma, fR, fCv));
    writef("  3D - From ConsVars: %19.12er\n", enthalpy_cv(cons3d, fGamma, fR, fCv));
    writef("\n");
    writef("Stagnation Enthalpy (J/m^3) functions:\n");
    writef("      Test Reference: %19.12er\n", enthStag);
    writef("  1D - From ConsVars: %19.12er\n", enthalpy_stagnation_cv(cons1d, fGamma));
    writef("  2D - From ConsVars: %19.12er\n", enthalpy_stagnation_cv(cons2d, fGamma));
    writef("  3D - From ConsVars: %19.12er\n", enthalpy_stagnation_cv(cons3d, fGamma));
    writef("\n");
    writef("Speed of Sound (m/s) and Mach (non-dimensional) functions:\n");
    writef("  Speed of Sound from Dens and Pres: %19.12er,     from Temp: %19.12er\n",
        sound_speed(dens, pres, fGamma), sound_speed_temp(temp, fGamma, fR));
    writef("  1D - Speed of Sound from ConsVars: %19.12er, from PrimVars: %19.12er - Mach: %19.12er\n",
        sound_speed_cv(cons1d, fGamma), sound_speed_pv(prim1d, fGamma), mach_cv(cons1d, fGamma));
    writef("  2D - Speed of Sound from ConsVars: %19.12er, from PrimVars: %19.12er - Mach: %19.12er\n",
        sound_speed_cv(cons2d, fGamma), sound_speed_pv(prim2d, fGamma), mach_cv(cons2d, fGamma));
    writef("  3D - Speed of Sound from ConsVars: %19.12er, from PrimVars: %19.12er - Mach: %19.12er\n",
        sound_speed_cv(cons3d, fGamma), sound_speed_pv(prim3d, fGamma), mach_cv(cons3d, fGamma));
    writef("\n");
    writef("Euler flux functions:\n");
    writef("  1D Fluxes:\n");
    writef("    1D  Euler Flux (ConsVars): %.8t\n", euler_flux_cv_1d(cons1d, fGamma));
    writef("    Gen Euler Flux (ConsVars): %.8t\n", euler_flux_cv(cons1d, fGamma));
    writef("    Gen Euler Flux (PrimVars): %.8t\n", euler_flux_pv(prim1d, fGamma));
    writef("\n");
    writef("  2D Fluxes:\n");
    writef("    Gen Euler Flux (ConsVars): %.8t\n", euler_flux_cv(cons2d, fGamma)[1,..]);
    writef("                               %.8t\n", euler_flux_cv(cons2d, fGamma)[2,..]);
    writef("    Gen Euler Flux (PrimVars): %.8t\n", euler_flux_pv(prim2d, fGamma)[1,..]);
    writef("                               %.8t\n", euler_flux_pv(prim2d, fGamma)[2,..]);
    writef("\n");
    writef("  3D Fluxes:\n");
    writef("    Gen Euler Flux (ConsVars): %.8t\n", euler_flux_cv(cons3d, fGamma)[1,..]);
    writef("                               %.8t\n", euler_flux_cv(cons3d, fGamma)[2,..]);
    writef("                               %.8t\n", euler_flux_cv(cons3d, fGamma)[3,..]);
    writef("    Gen Euler Flux (PrimVars): %.8t\n", euler_flux_pv(prim3d, fGamma)[1,..]);
    writef("                               %.8t\n", euler_flux_pv(prim3d, fGamma)[2,..]);
    writef("                               %.8t\n", euler_flux_pv(prim3d, fGamma)[3,..]);
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
