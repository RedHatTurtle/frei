module Boundary
{
  use Random;
  use UnitTest;
  import LinearAlgebra.norm;
  import Mesh.faml_r;
  import Input.fGamma;
  import Input.fR;

  proc boundary(hostConsVars : [] real, faml : faml_r, xyz : [] real, nrm : [] real) : [hostConsVars.domain] real
  {
    use Parameters.ParamInput;
    import Ringleb.ringleb_sol;

    var ghstConsVars : [hostConsVars.domain] real = 0;

    select faml.bocoType
    {
      when BC_TYPE_FLOW {}
      when BC_TYPE_INFLOW do
        select faml.bocoSubType
        {
          when BC_SUBTYPE_GENERIC_INFLOW {}
          when BC_SUBTYPE_SUB_INFLOW do
            ghstConsVars = sub_inflow(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_SUP_INFLOW {}
          when BC_SUBTYPE_MDOT_INFLOW {}
        }
      when BC_TYPE_OUTFLOW do
        select faml.bocoSubType
        {
          when BC_SUBTYPE_GENERIC_OUTFLOW {}
          when BC_SUBTYPE_SUB_OUTFLOW do
            ghstConsVars = sub_outflow(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_SUP_OUTFLOW do
            ghstConsVars = sup_outflow(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_MDOT_OUTFLOW {}
        }
      when BC_TYPE_OPENING do
        select faml.bocoSubType
        {
          when BC_SUBTYPE_RIEMANN do
            ghstConsVars = riemann(hostConsVars, faml.bocoProperties, nrm);
          when BC_SUBTYPE_GENERIC_FREEFLOW {}
          when BC_SUBTYPE_FREESTREAM {}
          when BC_SUBTYPE_FIXED {}
        }
      when BC_TYPE_WALL do
        select faml.bocoSubType
        {
          when BC_SUBTYPE_SLIP_WALL do
            ghstConsVars = symmetry(hostConsVars, nrm);
          when BC_SUBTYPE_EULER_WALL do
            ghstConsVars = symmetry(hostConsVars, nrm);
          when BC_SUBTYPE_ADIABATIC_WALL do
            ghstConsVars = wall_adiabatic(hostConsVars);
          when BC_SUBTYPE_ISOTHERMAL_WALL {}
          when BC_SUBTYPE_DEFAULT_WALL do
            ghstConsVars = symmetry(hostConsVars, nrm);
        }
      when BC_TYPE_SPECIAL do
        select faml.bocoSubType
        {
          when BC_SUBTYPE_SYMMETRY {}
          when BC_SUBTYPE_PERIODIC do
            ghstConsVars = periodic(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_DIRICHLET do
            ghstConsVars = dirichlet(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_1D_NOZZLE_CRIT_INFLOW do
            ghstConsVars = nozzle_ideal_inflow(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_1D_NOZZLE_SUBSONIC_INFLOW do
            ghstConsVars = nozzle_subsonic_inflow(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_1D_NOZZLE_SHOCKED_INFLOW do
            ghstConsVars = nozzle_shocked_inflow(hostConsVars, faml.bocoProperties);
          when BC_SUBTYPE_RINGLEB do
            ghstConsVars = 2*ringleb_sol(xyz[1..2]) - hostConsVars;
          when BC_SUBTYPE_MMS_DIRICHLET {}
        }
    }

    return ghstConsVars;
  }

  proc symmetry(hostConsVars : [] real, nrm : [] real) : [hostConsVars.domain] real
  {
    var idxDens : int   = hostConsVars.domain.dim(0).low;        // First element is density
    var idxMom  : range = hostConsVars.domain.dim(0).expand(-1); // Intermediary elements are the momentum vector
    var idxEner : int   = hostConsVars.domain.dim(0).high;       // Last element is energy

    var uniNrm : [nrm.domain] real = nrm/norm(nrm);

    var ghstConsVars : [hostConsVars.domain] real;

    // Preserve density and energy
    ghstConsVars[idxDens] = hostConsVars[idxDens];
    ghstConsVars[idxEner] = hostConsVars[idxEner];

    // Reflect momentum relative to the face normal
    ghstConsVars[idxMom]  = hostConsVars[idxMom] - 2.0*dot(hostConsVars[idxMom], uniNrm)*uniNrm;

    return ghstConsVars;
  }

  proc wall_adiabatic(hostConsVars : [] real) : [hostConsVars.domain] real
  {
    var idxDens : int   = hostConsVars.domain.dim(0).low;        // First element is density
    var idxMom  : range = hostConsVars.domain.dim(0).expand(-1); // Intermediary elements are the velocities
    var idxEner : int   = hostConsVars.domain.dim(0).high;       // Last element is energy

    var ghstConsVars : [hostConsVars.domain] real;

    ghstConsVars[idxDens] =  hostConsVars[idxDens];
    ghstConsVars[idxMom]  = -hostConsVars[idxMom] ;
    ghstConsVars[idxEner] =  hostConsVars[idxEner];

    return ghstConsVars;
  }

  proc sub_inflow(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Subsonic inflow boundary condition:
    //   - Fix density to free stream
    //   - Fix velocities to free stream
    //   - Extrapolate pressure
    //
    // Boco properties:
    //   - Free stream density
    //   - Free stream velocity components

    import LinearAlgebra.dot;
    import Input.fGamma;
    import Flux.pressure_cv;
    import Flux.temperature_cv;

    var ghstConsVars : [hostConsVars.domain] real;

    const idxRho : int   = hostConsVars.domain.dim(0).low;           // First element is density
    const idxMom : range = hostConsVars.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxEne : int   = hostConsVars.domain.dim(0).high;          // Last element is energy

    // Pre-compute the free stream specific kinetic energy
    const eneKin = 0.5 * bocoProperties[idxRho] * dot(bocoProperties[idxMom], bocoProperties[idxMom]);

    // Set all variables but pressure to free stream
    ghstConsVars = bocoProperties[hostConsVars.domain.dim(0)];

    // Get the interior pressure
    var presHost = pressure_cv( hostConsVars );

    // Compute the total energy using interior pressure and free stream specific kinetic energy
    ghstConsVars[idxEne] = presHost/(fGamma-1) + eneKin;

    return ghstConsVars;
  }

  proc sup_inflow(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Supersonic inflow boundary condition:
    //   - Fix all variables to free stream
    //
    // Boco properties:
    //   - Free stream density
    //   - Free stream velocity components
    //   - Free stream static pressure

    var ghstConsVars : [hostConsVars.domain] real;

    ghstConsVars = bocoProperties[hostConsVars.domain.dim(0)];

    return ghstConsVars;
  }

  proc sub_outflow(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Subsonic outflow boundary condition:
    //   - Extrapolate density
    //   - Extrapolate velocities
    //   - Fix pressure to freestream
    //
    // Boco properties:
    //   - Free stream static pressure

    import LinearAlgebra.norm;
    import Input.fGamma;

    var ghstConsVars : [hostConsVars.domain] real;

    const idxRho : int   = hostConsVars.domain.dim(0).low;           // First element is density
    const idxMom : range = hostConsVars.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    const idxEne : int   = hostConsVars.domain.dim(0).high;          // Last element is energy

    // Unit normal and tangent at this flux point
    // faceNormUnit = faceNorm / norm(faceNorm);

    // Calculate the normal velocity in the boundary
    // var velNorm = dot(faceNormalUnit, hostConsVars(idxMom)/hostConsVars(idxRho) );
    // var soundSpeed = sqrt( max( fGamma * hostConsVars[idxEne]/hvfp[idxRho] , rndoff ) );

    // Extrapolate Density and Momentum
    ghstConsVars[idxRho] = hostConsVars[idxRho];
    ghstConsVars[idxMom] = hostConsVars[idxMom];

    // If the normal velocity from the interior says this is an inflow boundary, set the ghost velocity to the
    // magnitude of the interior velocity multiplied by the unit normal
    // if (veln < 0) then
    //   ghostConsVars[idxMom] = norm( hostConsVars(idxMom) ) * faceNormUnit[..];

    // Check if the normal velocity is supersonic
    // if (abs(vn_int)/aspd_int >= one) then
         // If it's supersonic set the ghost pressure to the total pressure
         // ghstPres = total_pressure_cv(bocoConsVars);
    // else
         // If it's subsonic set the ghost pressure to the BoCo pressure
         var ghstPres = bocoProperties[1];
    // end if

    // Calculate Energy based on the extrapolated kinetic energy and the BoCo pressure
    var hostEneKin = 0.5*dot(hostConsVars[idxMom], hostConsVars[idxMom])/hostConsVars[idxRho];
    ghstConsVars[idxEne] = ghstPres/(fGamma-1.0) + hostEneKin;

    return ghstConsVars;
  }

  proc sup_outflow(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Supersonic outflow boundary condition:
    //   - Extrapolate all variables
    //
    // Boco properties:
    //   - None

    var ghstConsVars : [hostConsVars.domain] real;

    ghstConsVars = hostConsVars;

    return ghstConsVars;
  }

  proc riemann(hostConsVars : [] real, bocoProperties : [] real, nrm : [] real) : [hostConsVars.domain] real
  {
    var idxDens : int   = hostConsVars.domain.dim(0).low;        // First element is density
    var idxMom  : range = hostConsVars.domain.dim(0).expand(-1); // Intermediary elements are the velocities
    var idxEner : int   = hostConsVars.domain.dim(0).high;       // Last element is energy

    var uniNrm : [nrm.domain] real = nrm/norm(nrm);

    var ghstConsVars : [hostConsVars.domain] real;

    // Internal mesh condition
    var densInt : real = hostConsVars[idxDens];
    var  velInt : [idxMom] real = hostConsVars[idxMom]/densInt;
    var enerInt : real = hostConsVars[idxEner];

    var presInt : real = pressure_cv(hostConsVars);
    var    aInt : real = sqrt( fGamma * presInt / densInt );

    // Compute the internal invariants
    var velNrmInt : real = dot(velInt, uniNrm);
    var riemannP  : real = velNrmInt + 2.0*aInt/(fGamma-1);

    // External boundary condition
    var presExt : real = bocoProperties[idxDens];
    var  velExt : [idxMom] real = bocoProperties[idxMom];
    var tempExt : real = bocoProperties[idxEner];

    var densExt : real = presExt / ( fR * tempExt );
    var    aExt : real = sqrt( fGamma * fR * tempExt );

    // Compute the external invariants
    var velNrmExt : real = dot(velExt, uniNrm);
    var riemannM  : real = velNrmExt - 2.0*aExt/(fGamma-1);

    // If the flow is supersonic then use only internal or external variables to calculate the invariants
    if ( abs(velNrmInt/aInt) > 1.0 ) then
      if ( velNrmInt < 0.0 ) then
        riemannP = velNrmExt + 2.0*aExt/(fGamma-1);
      else
        riemannM = velNrmInt - 2.0*aInt/(fGamma-1);

    // Compute normal velocity component and speed of sound at the interface using the invariants
    var velNrm : real = (riemannP + riemannM)/2.0;
    var a      : real = (riemannP - riemannM)*(fGamma-1)/4.0;
    var vel    : [idxMom] real;
    var entr   : real;

    // Calculate the boundary properties based on the normal flow direction
    if ( velNrm < 0.0 )
    { // This is an inflow case
      vel = velExt + (velNrm-velNrmExt)*uniNrm;
      entr = (aExt**2)/(fGamma*densExt**(fGamma-1));
    }
    else
    { // This is an outflow case
      vel = velInt + (velNrm-velNrmInt)*uniNrm;
      entr = (aInt**2)/(fGamma*densInt**(fGamma-1));
    }

    // Compute the ghost variables
    var dens : real = ((a**2)/(fGamma*entr))**(1/(fGamma-1));
    var pres : real = dens*a**2/fGamma;
    var ener : real = pres/(fGamma-1) + dens*norm(vel)/2.0;

    ghstConsVars[idxDens] = 2*dens     - densInt;
    ghstConsVars[idxMom ] = 2*dens*vel - densInt*velInt;
    ghstConsVars[idxEner] = 2*ener     - enerInt;

    return ghstConsVars;
  }

  proc dirichlet(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Dirichlet boundary condition:
    //   - Fix all variables to defined parameters
    //
    // Boco properties:
    //   - Free stream density
    //   - Free stream momentum components
    //   - Free stream energy

    var ghstConsVars : [hostConsVars.domain] real;

    ghstConsVars = bocoProperties[hostConsVars.domain.dim(0)];

    return ghstConsVars;
  }

  proc nozzle_ideal_inflow(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Use cases:
    //   - Calculate inlet mach and density for ideal solution. Requires the area ratio.
    //   - Calculate inlet mach and density for subsonic solution. Requires an area ratio with the critical area defined
    //     by the exit pressure.

    // Get data from the initial condition
    use Parameters.ParamInput;
    import Init.flow_condition;

    var xyz : [1..1, 1..1] real = 0;

    // Local area + Critical area -> Local Mach -> Primitive Variables
    // Local Pressure -> Local Mach -> Primitive Variables

    var ghstConsVars : [hostConsVars.domain] real = reshape(flow_condition(IC_1D_NOZZLE_SMOOTH_TRANSONIC, [0.0], xyz),
        hostConsVars.domain);

    return ghstConsVars;
  }

  proc nozzle_subsonic_inflow(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Use cases:
    //   - Calculate inlet mach and density for ideal solution. Requires the area ratio.
    //   - Calculate inlet mach and density for subsonic solution. Requires an area ratio with the critical area defined
    //     by the exit pressure.

    // Get data from the initial condition
    use Parameters.ParamInput;
    import Init.flow_condition;

    var xyz : [1..1, 1..1] real = 0;

    // Local area + Critical area -> Local Mach -> Primitive Variables
    // Local Pressure -> Local Mach -> Primitive Variables

    var ghstConsVars : [hostConsVars.domain] real = reshape(flow_condition(IC_1D_NOZZLE_SUBSONIC, [0.0], xyz),
        hostConsVars.domain);

    return ghstConsVars;
  }

  proc nozzle_shocked_inflow(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Use cases:
    //   - Calculate inlet mach and density for ideal solution. Requires the area ratio.
    //   - Calculate inlet mach and density for subsonic solution. Requires an area ratio with the critical area defined
    //     by the exit pressure.

    // Get data from the initial condition
    use Parameters.ParamInput;
    import Init.flow_condition;

    var xyz : [1..1, 1..1] real = 0;

    // Local area + Critical area -> Local Mach -> Primitive Variables
    // Local Pressure -> Local Mach -> Primitive Variables

    var ghstConsVars : [hostConsVars.domain] real = reshape(flow_condition(IC_1D_NOZZLE_SHOCKED_TRANSONIC, [0.0], xyz),
        hostConsVars.domain);

    return ghstConsVars;
  }

  proc periodic(hostConsVars : [] real, bocoProperties : [] real) : [hostConsVars.domain] real
  {
    // Get the family name of the corresponding BoCo
    // What criteria can be used reliably to match FPs?
    // It seems safe to assume that the boundaries have the same orientation and lenght
    //   1. Require a specified "translation vector"
    //   2. Add to that a rotatation axis and angle

    // Pre-Processing:
    //  - Build a list of matching FPs

    // 1. Get index of corresponding FP
    // 2. Copy conserved variables - frMesh.solFP[meshFP, side, var]
    // 3. Copy interpolated discontinuous flux - frMesh.flxFP[meshFP, side, var]

    return 0.0;
  }

  proc main()
  {}
}
