module Dimensional
{
  type scales_t = unmanaged scales_c?;
  var scales : scales_t;

  proc init_scales(lengRef, velMRef, tempRef, presRef) do
    Dimensional.scales = new unmanaged scales_c(lengRef, velMRef, tempRef, presRef);

  class scales_c
  {
    // Directly defined scales
    const lengScale : real = 1.0; // Metre
    const velMScale : real = 1.0; // Metres per second
    const tempScale : real = 1.0; // Kelvin
    const presScale : real = 1.0; // Pascal = kg/(m*s^2)

    // Other fundamental units' scales
    const timeScale : real = 1.0; // Second

    // Heat capacity scale
    const heatScale : real = 1.0; // J/(kg*K) = m^2/(s^2*K)

    // Mass Scale
    const massScale : real = 1.0; // Kilogram

    // Conserved variables' scales
    const densScale : real = 1.0; // kg/m^3
    const momMScale : real = 1.0; // kg/(m^2*s)
    const enerScale : real = 1.0; // J/m^3 = kg/(m*s^2)

    // Viscosity related scales (not used yet)
    //const viscScale : real = 1.0;
    //const reynScale : real = 1.0; // Not really a scale since Reynolds is non-dimensional
    //const pranScale : real = 1.0; // Also not really a scale since Prandtl is non-dimensional

    proc init(lengRef, velMRef, tempRef, presRef)
    {
      // The free parameters of the scaling are directly defined. These variables were selected because they are the
      // ones most commonly used to define an open flow / atmospheric condition.
      this.lengScale = lengRef;
      this.velMScale = velMRef;
      this.tempScale = tempRef;
      this.presScale = presRef;

      // To preserve the validity of regular physical equations without any added constants all scales must be
      // dimensionally compatible. For example, specifying a velocity and a length scale implicitly defines a unique
      // time scale that preserves the validity of the equation $Velocity = C * Distance / Time$ with the constant C=1.
      // Therefore to prevent the need for extra parameters in the functions for these equations we enforce the
      // dimensional compatibility of the scales. Naturally non-dimensional parameters, like the heat capacity ratio,
      // don't need to be scaled.

      // Time scale defined by the length and velocity scales as explained in example above
      this.timeScale = this.lengScale/this.velMScale;

      // Heat capacity scales according to length, time and temperature
      // i.e. Heat Capacity = leng^2 / (time^2 * temp)
      this.heatScale = this.velMScale * this.velMScale / this.tempScale;

      // Mass scale must be compatible with the prescribed pressure scale
      // I.e. Pres = force / leng^2 = mass / (leng * time^2)
      this.massScale = this.presScale * this.lengScale * this.timeScale * this.timeScale;

      // Calculate conserved variables' scales. Momentum and Energy scales are different from their
      // typical units (kg*m/s and J) because volume specific quantities are used.
      this.densScale = this.presScale / this.velMScale / this.velMScale;
      this.momMScale = this.densScale * this.velMScale;
      this.enerScale = this.presScale;

      // Viscosity coefficient
      //this.viscScale = this.presScale * this.timeScale;

      // Momentum and heat diffusivity parameters. These are not dimensional numbers but Reynolds
      // appears in the non-dimensional form of the Navier-Stokes equations. And Prandtl might show
      // up somewhere in the future too.

      // Reynolds = Dens * Velocity * Length / Dynamic Viscosity
      //this.reynScale = this.densScale * this.lengScale * this.velMScale / this.viscScale;

      // Prandtl = Cp * Dynamic Viscosity / Thermal Conductivity
      //this.pranScale = this.heatScale * this.viscScale / this.condScale;
    }

    // Dim -> Non
    proc dim2non_time(time :    real) :    real do return time/timeScale;
    proc dim2non_leng(leng :    real) :    real do return leng/lengScale;

    proc dim2non_dens(dens :    real) :    real do return dens/densScale;
    proc dim2non_momv(momV : [] real) : [] real do return momV/momMScale;
    proc dim2non_momm(momM :    real) :    real do return momM/momMScale;
    proc dim2non_ener(ener :    real) :    real do return ener/enerScale;

    proc dim2non_velv(velV : [] real) : [] real do return velV/velMScale;
    proc dim2non_velm(velM :    real) :    real do return velM/velMScale;
    proc dim2non_temp(temp :    real) :    real do return temp/tempScale;
    proc dim2non_pres(pres :    real) :    real do return pres/presScale;

    proc dim2non_heat(heat :    real) :    real do return heat/heatScale;

    // Non -> Dim
    proc non2dim_time(time :    real) :    real do return time*timeScale;
    proc non2dim_leng(leng :    real) :    real do return leng*lengScale;

    proc non2dim_dens(dens :    real) :    real do return dens*densScale;
    proc non2dim_momv(momV : [] real) : [] real do return momV*momMScale;
    proc non2dim_momm(momM :    real) :    real do return momM*momMScale;
    proc non2dim_ener(ener :    real) :    real do return ener*enerScale;

    proc non2dim_velv(velV : [] real) : [] real do return velV*velMScale;
    proc non2dim_velm(velM :    real) :    real do return velM*velMScale;
    proc non2dim_temp(temp :    real) :    real do return temp*tempScale;
    proc non2dim_pres(pres :    real) :    real do return pres*presScale;

    proc non2dim_heat(heat :    real) :    real do return heat*heatScale;

    proc dim2non_cv(consDim : [] real) : [] real
    {
      var consNon : [consDim.domain] real;

      var idxDens : int   = consDim.domain.dim(0).low;        // First element is density
      var idxMomV : range = consDim.domain.dim(0).expand(-1); // Intermediary elements are the momentum vector components
      var idxEner : int   = consDim.domain.dim(0).high;       // Last element is energy

      consNon[idxDens] = consDim[idxDens] / this.densScale;
      consNon[idxMomV] = consDim[idxMomV] / this.momMScale;
      consNon[idxEner] = consDim[idxEner] / this.enerScale;

      return consNon;
    }

    proc non2dim_cv(consNon : [] real) : [] real
    {
      var consDim : [consNon.domain] real;

      var idxDens : int   = consNon.domain.dim(0).low;        // First element is density
      var idxMomV : range = consNon.domain.dim(0).expand(-1); // Intermediary elements are the momentum vector components
      var idxEner : int   = consNon.domain.dim(0).high;       // Last element is energy

      consDim[idxDens] = consNon[idxDens] * this.densScale;
      consDim[idxMomV] = consNon[idxMomV] * this.momMScale;
      consDim[idxEner] = consNon[idxEner] * this.enerScale;

      return consDim;
    }

    proc dim2non_pv(primDim : [] real) : [] real
    {
      var primNon : [primDim.domain] real;

      var idxDens : int   = primDim.domain.dim(0).low;        // First element is density
      var idxVelV : range = primDim.domain.dim(0).expand(-1); // Intermediary elements are the velocity vector components
      var idxPres : int   = primDim.domain.dim(0).high;       // Last element is energy

      primNon[idxDens] = primDim[idxDens] / this.densScale;
      primNon[idxVelV] = primDim[idxVelV] / this.velMScale;
      primNon[idxPres] = primDim[idxPres] / this.enerScale;

      return primNon;
    }

    proc non2dim_pv(primNon : [] real) : [] real
    {
      var primDim : [primNon.domain] real;

      var idxDens : int   = primNon.domain.dim(0).low;        // First element is density
      var idxVelV : range = primNon.domain.dim(0).expand(-1); // Intermediary elements are the velocity vector components
      var idxPres : int   = primNon.domain.dim(0).high;       // Last element is energy

      primDim[idxDens] = primNon[idxDens] * this.densScale;
      primDim[idxVelV] = primNon[idxVelV] * this.velMScale;
      primDim[idxPres] = primNon[idxPres] * this.enerScale;

      return primDim;
    }
  }

  proc main()
  {
    use Testing;

    // Fluid properties
    param gasR     : real = 8.314462618;            // Ideal Gas Constant (2018 CODATA), J/(mol*K)
    param fMolMass : real = 2.8966e-2;              // Molar Mass of the fluid mixture, kg/mol
    param fGamma   : real = 1.4;                    // Ideal diatomic gas specific heat ratio
    param fR       : real = gasR/fMolMass;          // Calculate the specific gas constant for this fluid, J/(kg*K)
    param fCp      : real = fR/(1.0-1.0/fGamma);    // Calculate the heat coefficients at constant pressure Cp, J/(kg*K)
    param fCv      : real = fR/(fGamma-1.0);        // Calculate the heat coefficients at constant volume Cv, J/(kg*K)
  //param fVisc    : real = 1.78792e-5;             // Dynamic viscosity coefficient, kg/(m*s)
  //param fKappa   : real = 2.52992e-2;             // Thermal Conductivity, W/(m*K)
  //param fPrandtl : real = fVisc*fCp/fKappa;       // Calculate the Prandtl number of the fluid

    // Define scale calculation parameters
    param lengRef : real = 2.0;
    param velMRef : real = 241.92454098943761000787319395449783147814277224947322;
    param tempRef : real = 216.6;
    param presRef : real = 17_933.0;

    // Define test flow conditions (~ Standard atmosphere @ 12.5km)
    // https://www.pdas.com/atmosTable2SI.html
    param time : real = 1.0e-4;
    param leng : real =      1.0;
    param temp : real =    216.6;
    param pres : real = 17_933.0;
    param mach : real = 0.82;
    param heat : real = fR;
    param alph : real = 0.03490658503988659153847381536977225426885743777083; // 2.0 degrees in Radians
    param beta : real = 0.00872664625997164788461845384244306356721435944271; // 0.5 degree  in Radians
    //param visc : real = 1.422e-5 // kg/(m*s)

    // Calculate remaining required properties with a 50 digit calculator
    param dens : real =       0.28843558377998645594115862622493950794900043071807; // Pres/(fR*Temp)
    param aSpd : real =     295.02992803589952439984535848109491643675947835301612; // sqrt(fGamma*fR*Temp)
    param velM : real =     241.92454098943761000787319395449783147814277224947322; // Mach*aSpd
    param momM : real =      69.77964621099369921629054180976936089541572335796040; // Dens*Mach*aSpd
    param ener : real =  53_273.20444000000000000000000000000000000000000000000000; // Pres/(fGamma-1) + 0.5*Dens*velM^2
    const velV1d : [1..1] real = velM*[1                                                  ];
    const velV2d : [1..2] real = velM*[cos(alph)          , sin(alph)                     ];
    const velV3d : [1..3] real = velM*[cos(alph)*cos(beta), sin(alph)*cos(beta), sin(beta)];
    const momV1d : [1..1] real = momM*[1                                                  ];
    const momV2d : [1..2] real = momM*[cos(alph)          , sin(alph)                     ];
    const momV3d : [1..3] real = momM*[cos(alph)*cos(beta), sin(alph)*cos(beta), sin(beta)];

    // Build conserved and primitive variables arrays
    const consVars1dDim : [1..3] real = [dens, dens*velV1d[1], ener];
    const consVars2dDim : [1..4] real = [dens, dens*velV2d[1], dens*velV2d[2], ener];
    const consVars3dDim : [1..5] real = [dens, dens*velV3d[1], dens*velV3d[2], dens*velV3d[3], ener];
    const primVars1dDim : [1..3] real = [dens, velV1d[1], pres];
    const primVars2dDim : [1..4] real = [dens, velV2d[1], velV2d[2], pres];
    const primVars3dDim : [1..5] real = [dens, velV3d[1], velV3d[2], velV3d[3], pres];

    // Reference correct ratios
    param timeRatio : real = 0.01209622704947188050039365969772489157390713861247;
    param lengRatio : real = 0.5;
    param densRatio : real = 0.94136000000000000000000000000000000000000000000001;
    param tempRatio : real = 1.0;
    param presRatio : real = 1.0;
    param velMRatio : real = 1.0;
    param momMRatio : real = densRatio*velMRatio;
    param heatRatio : real = 1.06229285289368573128239993201325741480411319792640;
    param enerRatio : real = 2.97068000000000000000000000000000000000000000000000;
    const velV1dNon : [1..1] real = [1                                                  ];
    const velV2dNon : [1..2] real = [cos(alph)          , sin(alph)                     ];
    const velV3dNon : [1..3] real = [cos(alph)*cos(beta), sin(alph)*cos(beta), sin(beta)];
    const momV1dNon : [1..1] real = densRatio*velMRatio*[1                                                  ];
    const momV2dNon : [1..2] real = densRatio*velMRatio*[cos(alph)          , sin(alph)                     ];
    const momV3dNon : [1..3] real = densRatio*velMRatio*[cos(alph)*cos(beta), sin(alph)*cos(beta), sin(beta)];

    // Build conserved and primitive variables arrays
    const consVars1dNon : [1..3] real = [densRatio, densRatio*velV1dNon[1], enerRatio];
    const consVars2dNon : [1..4] real = [densRatio, densRatio*velV2dNon[1], densRatio*velV2dNon[2], enerRatio];
    const consVars3dNon : [1..5] real = [densRatio, densRatio*velV3dNon[1], densRatio*velV3dNon[2], densRatio*velV3dNon[3], enerRatio];
    const primVars1dNon : [1..3] real = [densRatio, velV1dNon[1], presRatio];
    const primVars2dNon : [1..4] real = [densRatio, velV2dNon[1], velV2dNon[2], presRatio];
    const primVars3dNon : [1..5] real = [densRatio, velV3dNon[1], velV3dNon[2], velV3dNon[3], presRatio];

    // Initialize the scales object with the most common reference values
    init_scales(lengRef, velMRef, tempRef, presRef);

    // Testing single variable methods
    writef("Function Name | Reference           | Result              | Abs Error | Rel Error\n");
    writef("--------------|---------------------|---------------------|-----------|----------\n");
    writef(" Dim2Non Time | %19.12er | %19.12er | %9.2er | %9.2er\n",   timeRatio, scales!.dim2non_time(time     ),      error(timeRatio, scales!.dim2non_time(time     )), relative_error(timeRatio, scales!.dim2non_time(time     )));
    writef("         Leng | %19.12er | %19.12er | %9.2er | %9.2er\n",   lengRatio, scales!.dim2non_leng(leng     ),      error(lengRatio, scales!.dim2non_leng(leng     )), relative_error(lengRatio, scales!.dim2non_leng(leng     )));
    writef("         Dens | %19.12er | %19.12er | %9.2er | %9.2er\n",   densRatio, scales!.dim2non_dens(dens     ),      error(densRatio, scales!.dim2non_dens(dens     )), relative_error(densRatio, scales!.dim2non_dens(dens     )));
    writef("         Momm | %19.12er | %19.12er | %9.2er | %9.2er\n",   momMRatio, scales!.dim2non_momm(momM     ),      error(momMRatio, scales!.dim2non_momm(momM     )), relative_error(momMRatio, scales!.dim2non_momm(momM     )));
    writef("         Ener | %19.12er | %19.12er | %9.2er | %9.2er\n",   enerRatio, scales!.dim2non_ener(ener     ),      error(enerRatio, scales!.dim2non_ener(ener     )), relative_error(enerRatio, scales!.dim2non_ener(ener     )));
    writef("         Velm | %19.12er | %19.12er | %9.2er | %9.2er\n",   velMRatio, scales!.dim2non_velm(velM     ),      error(velMRatio, scales!.dim2non_velm(velM     )), relative_error(velMRatio, scales!.dim2non_velm(velM     )));
    writef("         Temp | %19.12er | %19.12er | %9.2er | %9.2er\n",   tempRatio, scales!.dim2non_temp(temp     ),      error(tempRatio, scales!.dim2non_temp(temp     )), relative_error(tempRatio, scales!.dim2non_temp(temp     )));
    writef("         Pres | %19.12er | %19.12er | %9.2er | %9.2er\n",   presRatio, scales!.dim2non_pres(pres     ),      error(presRatio, scales!.dim2non_pres(pres     )), relative_error(presRatio, scales!.dim2non_pres(pres     )));
    writef("         Heat | %19.12er | %19.12er | %9.2er | %9.2er\n",   heatRatio, scales!.dim2non_heat(heat     ),      error(heatRatio, scales!.dim2non_heat(heat     )), relative_error(heatRatio, scales!.dim2non_heat(heat     )));
    writef("--------------|---------------------|---------------------|-----------|----------\n");
    writef(" Non2Dim Time | %19.12er | %19.12er | %9.2er | %9.2er\n",   time     , scales!.non2dim_time(timeRatio),      error(time     , scales!.non2dim_time(timeRatio)), relative_error(time     , scales!.non2dim_time(timeRatio)));
    writef("         Leng | %19.12er | %19.12er | %9.2er | %9.2er\n",   leng     , scales!.non2dim_leng(lengRatio),      error(leng     , scales!.non2dim_leng(lengRatio)), relative_error(leng     , scales!.non2dim_leng(lengRatio)));
    writef("         Dens | %19.12er | %19.12er | %9.2er | %9.2er\n",   dens     , scales!.non2dim_dens(densRatio),      error(dens     , scales!.non2dim_dens(densRatio)), relative_error(dens     , scales!.non2dim_dens(densRatio)));
    writef("         Momm | %19.12er | %19.12er | %9.2er | %9.2er\n",   momM     , scales!.non2dim_momm(momMRatio),      error(momM     , scales!.non2dim_momm(momMRatio)), relative_error(momM     , scales!.non2dim_momm(momMRatio)));
    writef("         Ener | %19.12er | %19.12er | %9.2er | %9.2er\n",   ener     , scales!.non2dim_ener(enerRatio),      error(ener     , scales!.non2dim_ener(enerRatio)), relative_error(ener     , scales!.non2dim_ener(enerRatio)));
    writef("         Velm | %19.12er | %19.12er | %9.2er | %9.2er\n",   velM     , scales!.non2dim_velm(velMRatio),      error(velM     , scales!.non2dim_velm(velMRatio)), relative_error(velM     , scales!.non2dim_velm(velMRatio)));
    writef("         Temp | %19.12er | %19.12er | %9.2er | %9.2er\n",   temp     , scales!.non2dim_temp(tempRatio),      error(temp     , scales!.non2dim_temp(tempRatio)), relative_error(temp     , scales!.non2dim_temp(tempRatio)));
    writef("         Pres | %19.12er | %19.12er | %9.2er | %9.2er\n",   pres     , scales!.non2dim_pres(presRatio),      error(pres     , scales!.non2dim_pres(presRatio)), relative_error(pres     , scales!.non2dim_pres(presRatio)));
    writef("         Heat | %19.12er | %19.12er | %9.2er | %9.2er\n",   heat     , scales!.non2dim_heat(heatRatio),      error(heat     , scales!.non2dim_heat(heatRatio)), relative_error(heat     , scales!.non2dim_heat(heatRatio)));

    // Test array methods
    writef("\n");
    writef("Velocity Vector, Dim -> Non :\n");
    writef("  1D - Test Reference: %.8ht\n", velV1dNon);
    writef("         dim2non_velv: %.8ht\n", scales!.dim2non_velv(velV1d));
    writef("  2D - Test Reference: %.8ht\n", velV2dNon);
    writef("         dim2non_velv: %.8ht\n", scales!.dim2non_velv(velV2d));
    writef("  3D - Test Reference: %.8ht\n", velV3dNon);
    writef("         dim2non_velv: %.8ht\n", scales!.dim2non_velv(velV3d));
    writef("\n");
    writef("Momentum Vectors, Dim -> Non :\n");
    writef("  1D - Test Reference: %.8ht\n", momV1dNon);
    writef("         dim2non_momv: %.8ht\n", scales!.dim2non_momv(momV1d));
    writef("  2D - Test Reference: %.8ht\n", momV2dNon);
    writef("         dim2non_momv: %.8ht\n", scales!.dim2non_momv(momV2d));
    writef("  3D - Test Reference: %.8ht\n", momV3dNon);
    writef("         dim2non_momv: %.8ht\n", scales!.dim2non_momv(momV3d));
    writef("\n");
    writef("Velocity Vectors, Non -> Dim :\n");
    writef(" 1D -  Test Reference: %.8ht\n", velV1d);
    writef("         non2dim_velv: %.8ht\n", scales!.non2dim_velv(velV1dNon));
    writef(" 2D -  Test Reference: %.8ht\n", velV2d);
    writef("         non2dim_velv: %.8ht\n", scales!.non2dim_velv(velV2dNon));
    writef(" 3D -  Test Reference: %.8ht\n", velV3d);
    writef("         non2dim_velv: %.8ht\n", scales!.non2dim_velv(velV3dNon));
    writef("\n");
    writef("Momentum Vectors, Non -> Dim :\n");
    writef(" 1D -  Test Reference: %.8ht\n", momV1d);
    writef("         non2dim_momv: %.8ht\n", scales!.non2dim_momv(momV1dNon));
    writef(" 2D -  Test Reference: %.8ht\n", momV2d);
    writef("         non2dim_momv: %.8ht\n", scales!.non2dim_momv(momV2dNon));
    writef(" 3D -  Test Reference: %.8ht\n", momV3d);
    writef("         non2dim_momv: %.8ht\n", scales!.non2dim_momv(momV3dNon));

    // Test array methods
    writef("\n");
    writef("Conserved Variables, Dim -> Non :\n");
    writef("  1D - Test Reference: %.8ht\n", consVars1dNon);
    writef("           dim2non_cv: %.8ht\n", scales!.dim2non_cv(consVars1dDim));
    writef("  2D - Test Reference: %.8ht\n", consVars2dNon);
    writef("           dim2non_cv: %.8ht\n", scales!.dim2non_cv(consVars2dDim));
    writef("  3D - Test Reference: %.8ht\n", consVars3dNon);
    writef("           dim2non_cv: %.8ht\n", scales!.dim2non_cv(consVars3dDim));
    writef("\n");
    writef("Primitive Variables, Dim -> Non :\n");
    writef("  1D - Test Reference: %.8ht\n", primVars1dNon);
    writef("           dim2non_pv: %.8ht\n", scales!.dim2non_pv(primVars1dDim));
    writef("  2D - Test Reference: %.8ht\n", primVars2dNon);
    writef("           dim2non_pv: %.8ht\n", scales!.dim2non_pv(primVars2dDim));
    writef("  3D - Test Reference: %.8ht\n", primVars3dNon);
    writef("           dim2non_pv: %.8ht\n", scales!.dim2non_pv(primVars3dDim));
    writef("\n");
    writef("Conserved Variables, Non -> Dim :\n");
    writef(" 1D -  Test Reference: %.8ht\n", consVars1dDim);
    writef("           non2dim_cv: %.8ht\n", scales!.non2dim_cv(consVars1dNon));
    writef(" 2D -  Test Reference: %.8ht\n", consVars2dDim);
    writef("           non2dim_cv: %.8ht\n", scales!.non2dim_cv(consVars2dNon));
    writef(" 3D -  Test Reference: %.8ht\n", consVars3dDim);
    writef("           non2dim_cv: %.8ht\n", scales!.non2dim_cv(consVars3dNon));
    writef("\n");
    writef("Primitive Variables, Non -> Dim :\n");
    writef(" 1D -  Test Reference: %.8ht\n", primVars1dDim);
    writef("           non2dim_pv: %.8ht\n", scales!.non2dim_pv(primVars1dNon));
    writef(" 2D -  Test Reference: %.8ht\n", primVars2dDim);
    writef("           non2dim_pv: %.8ht\n", scales!.non2dim_pv(primVars2dNon));
    writef(" 3D -  Test Reference: %.8ht\n", primVars3dDim);
    writef("           non2dim_pv: %.8ht\n", scales!.non2dim_pv(primVars3dNon));

    // Dump the entire object
    writef("\n");
    writef("Scales object dump:\n");
    writef("  %.8ht\n", scales);
  }
}
