prototype module Riemann
{
  proc upwind_1d(consL : [1..3] real, consR : [1..3] real, nrm : [1..1] real) : [1..3] real
  {
    use Flux;

    var upwind : [consL.domain] real;

    if nrm[1] < 0 then
      upwind = convection_flux_cv_1d(consR);
    else
      upwind = convection_flux_cv_1d(consL);

    return upwind;
  }

  proc rusanov_1d(consL : [1..3] real, consR : [1..3] real, nrm : [1..1] real) : [1..3] real
  {
    use Input;
    import Flux.pressure_cv;
    import Flux.euler_flux_cv;

    var idxDens : int   = consL.domain.dim(0).low;           // First element is density
    var idxMom  : range = consL.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEner : int   = consL.domain.dim(0).high;          // Last element is energy

    // Primitive variables
    //  Left state
    var densL : real = consL[idxDens];
    var velL  : real = consL[idxMom]/consL[idxDens];
    var presL : real = pressure_cv(consL);
    var enthL : real = ( consL[idxEner] + presL ) / densL;
    //  Right state
    var densR : real = consR[idxDens];
    var velR  : real = consR[idxMom]/consR[idxDens];
    var presR : real = pressure_cv(consR);
    var enthR : real = ( consR[idxEner] + presR ) / densR;

    // Compute the Roe Averages
    var rt  : real = sqrt(densR/densL);
    var vel : real = (velL+rt*velR)/(1.0+rt);
    var enth: real = (enthL+rt*enthR)/(1.0+rt);
    var a   : real = sqrt( (fGamma-1.0)*(enth-0.5*vel*vel) );

    var smax : real = abs(vel) + a;

    // Compute the average flux.
    var rusanov : [1..3] real = 0.5*( euler_flux_cv(consL)[1,1..3] + euler_flux_cv(consR)[1,1..3] - nrm[1]*smax*(consR-consL));

    return rusanov;
  }

  proc roe_1d(consL : [1..3] real, consR : [1..3] real, nrm : [1..1] real) : [1..3] real
  {
    use Input;
    use LinearAlgebra;
    import Flux.pressure_cv;
    import Flux.euler_flux_cv;

    var idxDens : int   = consL.domain.dim(0).low;        // First element is density
    var idxMom  : range = consL.domain.dim(0).expand(-1); // Intermediary elements are the velocities
    var idxEner : int   = consL.domain.dim(0).high;       // Last element is energy

    var uniNrm : [nrm.domain] real = nrm/norm(nrm, normType.norm2);

    //Primitive and other variables.
    //  Left state
    var densL   : real = consL[1];
    var velL    : real = consL[2]/consL[1];
    var velNrmL : real = velL * uniNrm[1];
    var presL   : real = pressure_cv(consL);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( consL[3] + presL ) / densL;
    //  Right state
    var densR   : real = consR[1];
    var velR    : real = consR[2]/consR[1];
    var velNrmR : real = velR * uniNrm[1];
    var presR   : real = pressure_cv(consR);
    var aR      : real = sqrt(fGamma*presR/densR);
    var enthR   : real = ( consR[3] + presR ) / densR;

    //First compute the Roe Averages
    var RT     : real = sqrt(densR/densL);
    var dens   : real = RT*densL;
    var vel    : real = (velL+RT*velR)/(1.0+RT);
    var velNrm : real = (velNrmL+RT*velNrmR)/(1.0+RT);
    var enth   : real = (enthL+RT*enthR)/(1.0+RT);
    var a      : real = sqrt( (fGamma-1.0)*(enth-0.5*vel*vel) );

    //Differences in primitive variables.
    var dDens   : real =   densR - densL;
    var dVel    : real =    velR - velL;
    var dVelNrm : real = velNrmR - velNrmL;
    var dPres   : real =   presR - presL;

    //Wave strength (Characteristic Variables).
    var dV : [1..3] real;
    dV[1] = 0.5*(dPres-dens*a*dVelNrm)/(a*a);
    dV[2] = -( dPres/(a*a) - dDens );
    dV[3] = 0.5*(dPres+dens*a*dVelNrm)/(a*a);

    //Absolute values of the wave speeds (Eigenvalues)
    var ws : [1..3] real;
    ws[1] = abs(velNrm - a);
    ws[2] = abs(velNrm    );
    ws[3] = abs(velNrm + a);

    // Harten's Entropy Fix JCP(1983), 49, pp357-393. Only for the nonlinear fields.
    // It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
    if ( ws(1) < 0.2 ) then ws(1) = (ws(1)**2)/0.4+0.1;
    if ( ws(3) < 0.2 ) then ws(3) = (ws(3)**2)/0.4+0.1;

    //Right eigenvectors
    var R : [1..3,1..3] real;
    R[1,1] = 1.0;
    R[2,1] = vel  - a;
    R[3,1] = enth - a*velNrm;

    R[1,2] = 1.0;
    R[2,2] = vel;
    R[3,2] = 0.5*vel*vel;

    R[1,3] = 1.0;
    R[2,3] = vel  + a;
    R[3,3] = enth + a*velNrm;

    // Compute the let and right fluxes
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2] + ws[3]*dV[3]*R[.., 3];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc roe_2d(consL : [1..4] real, consR : [1..4] real, nrm : [1..2] real) : [1..4] real
  {
    use Input;
    use LinearAlgebra;
    import Flux.pressure_cv;
    import Flux.euler_flux_cv;

    var idxDens : int   = consL.domain.dim(0).low;        // First element is density
    var idxMom  : range = consL.domain.dim(0).expand(-1); // Intermediary elements are the velocities
    var idxEner : int   = consL.domain.dim(0).high;       // Last element is energy

    var uniNrm : [nrm.domain] real = nrm/norm(nrm, normType.norm2);

    //Primitive and other variables.
    //  Left state
    ref densL   : real = consL[idxDens];
    ref enerL   : real = consL[idxEner];
    var velL    : [idxMom-1] real = consL[idxMom]/densL;
    var velNrmL : real = dot[velL, uniNrm];
    var presL   : real = pressure_cv(consL);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL;

    //  Right state
    ref densR   : real = consR[idxDens];
    ref enerR   : real = consR[idxEner];
    var velR    : [idxMom-1] real = consR[idxMom]/densR;
    var velNrmR : real = dot[velR, uniNrm];
    var presR   : real = pressure_cv(consR);
    var aR      : real = sqrt(fGamma*presR/densR);
    var enthR   : real = ( enerR + presR ) / densR;

    // Compute the Roe Averages
    var RT     : real = sqrt( densR/densL );
    var dens   : real = RT*densL;
    var vel    : [idxMom-1] real = (velL + RT*velR)/(1.0 + RT);
    var velNrm : real = dot(vel, uniNrm);
    var enth   : real = (enthL + RT*enthR)/(1.0 + RT);
    var a      : real = sqrt( (fGamma-1.0)*(enth-0.5*dot(vel, vel)) );

    //Differences in primitive variables
    var dDens   = densR   - densL;
    var dVelNrm = velNrmR - velNrmL;
    var dPres   = presR   - presL;
    var dVel    = velR    - velL;

    // Wave strengths (Characteristic Variables).
    var dV : [1..4] real;
    dV[1] =  (dPres - dens*a*dVelNrm)/(2.0*a**2);
    dV[2] =   dDens - dPres/(a**2);
    dV[3] =  (dPres + dens*a*dVelNrm)/(2.0*a**2);
    dV[4] =   dens;
    //dV[4] =   dens*dVelTng/a;

    // Absolute values of the wave speeds (Eigenvalues)
    var ws : [1..4] real;
    ws[1] = abs(velNrm - a);
    ws[2] = abs(velNrm    );
    ws[3] = abs(velNrm + a);
    ws[4] = abs(velNrm    );

    // Harten's Entropy Fix JCP(1983), 49, pp357-393. Only for the nonlinear fields.
    // It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
    if ( ws[1] < 0.2 ) then ws[1] = (ws[1]**2)/0.4 + 0.1;
    if ( ws[3] < 0.2 ) then ws[3] = (ws[3]**2)/0.4 + 0.1;

    // Right eigenvectors
    var R : [1..4, 1..4] real;

    // Negative acoustic wave
    R[idxDens, 1] = 1.0;
    R[idxMom , 1] = vel  - a*uniNrm;
    R[idxEner, 1] = enth - a*velNrm;

    // Entropy wave
    R[idxDens, 2] = 1.0;
    R[idxMom , 2] = vel;
    R[idxEner, 2] = dot(vel, vel)/2.0;

    // Positive acoustic wave
    R[idxDens, 3] = 1.0;
    R[idxMom , 3] = vel  + a*uniNrm;
    R[idxEner, 3] = enth + a*velNrm;

    // Shear wave
    R[idxDens, 4] = 0.0;
    R[idxMom , 4] = dVel - dVelNrm*uniNrm;
    R[idxEner, 4] = dot(vel, dVel) - velNrm*dVelNrm;

    // Compute the let and right fluxes
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2]
                                   + ws[3]*dV[3]*R[.., 3] + ws[4]*dV[4]*R[.., 4];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc roe_3d(consL : [1..5] real, consR : [1..5] real, nrm : [1..3] real) : [1..5] real
  {
    use Input;
    use LinearAlgebra;
    import Flux.pressure_cv;
    import Flux.euler_flux_cv;

    var idxDens : int   = consL.domain.dim(0).low;        // First element is density
    var idxMom  : range = consL.domain.dim(0).expand(-1); // Intermediary elements are the velocities
    var idxEner : int   = consL.domain.dim(0).high;       // Last element is energy

    var uniNrm : [nrm.domain] real = nrm/norm(nrm, normType.norm2);

    //Primitive and other variables.
    //  Left state
    ref densL   : real = consL[idxDens];
    ref enerL   : real = consL[idxEner];
    var velL    : [idxMom-1] real = consL[idxMom]/densL;
    var velNrmL : real = dot[velL, uniNrm];
    var presL   : real = pressure_cv(consL);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL;

    //  Right state
    ref densR   : real = consR[idxDens];
    ref enerR   : real = consR[idxEner];
    var velR    : [idxMom-1] real = consR[idxMom]/densR;
    var velNrmR : real = dot[velR, uniNrm];
    var presR   : real = pressure_cv(consR);
    var aR      : real = sqrt(fGamma*presR/densR);
    var enthR   : real = ( enerR + presR ) / densR;

    // Compute the Roe Averages
    var RT     : real = sqrt( densR/densL );
    var dens   : real = RT*densL;
    var vel    : [idxMom-1] real = (velL + RT*velR)/(1.0 + RT);
    var velNrm : real = dot(vel, uniNrm);
    var enth   : real = (enthL + RT*enthR)/(1.0 + RT);
    var a      : real = sqrt( (fGamma-1.0)*(enth-0.5*dot(vel, vel)) );

    //Differences in primitive variables
    var dDens   = densR   - densL;
    var dVelNrm = velNrmR - velNrmL;
    var dPres   = presR   - presL;
    var dVel    = velR    - velL;

    // Wave strengths (Characteristic Variables).
    var dV : [1..4] real;
    dV[1] =  (dPres - dens*a*dVelNrm)/(2.0*a**2);
    dV[2] =   dDens - dPres/(a**2);
    dV[3] =  (dPres + dens*a*dVelNrm)/(2.0*a**2);
    dV[4] =   dens;

    // Absolute values of the wave speeds (Eigenvalues)
    var ws : [1..4] real;
    ws[1] = abs(velNrm - a);
    ws[2] = abs(velNrm    );
    ws[3] = abs(velNrm + a);
    ws[4] = abs(velNrm    );

    // Harten's Entropy Fix JCP(1983), 49, pp357-393. Only for the nonlinear fields.
    // It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
    if ( ws[1] < 0.2 ) then ws[1] = (ws[1]**2)/0.4 + 0.1;
    if ( ws[3] < 0.2 ) then ws[3] = (ws[3]**2)/0.4 + 0.1;

    // Right eigenvectors
    var R : [1..5, 1..4] real;

    // Negative acoustic wave
    R[idxDens, 1] = 1.0;
    R[idxMom , 1] = vel  - a*uniNrm;
    R[idxEner, 1] = enth - a*velNrm;

    // Entropy wave
    R[idxDens, 2] = 1.0;
    R[idxMom , 2] = vel;
    R[idxEner, 2] = dot(vel, vel)/2.0;

    // Positive acoustic wave
    R[idxDens, 3] = 1.0;
    R[idxMom , 3] = vel  + a*uniNrm;
    R[idxEner, 3] = enth + a*velNrm;

    // Combined shear waves
    R[idxDens, 4] = 0.0;
    R[idxMom , 4] = dVel - dVelNrm*uniNrm;
    R[idxEner, 4] = dot(vel, dVel) - velNrm*dVelNrm;

    // Compute the let and right fluxes
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2]
                                   + ws[3]*dV[3]*R[.., 3] + ws[4]*dV[4]*R[.., 4];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc rusanov(consL : [] real, consR : [] real, nrm : [] real) : [] real
  {
    // Generic 1D/2D/3D Rusanov solver
  }

  proc roe(consL : [] real, consR : [] real, nrm : [] real) : [] real
  {
    use Input;
    use LinearAlgebra;
    import Flux.pressure_cv;
    import Flux.euler_flux_cv;

    var idxDens : int   = consL.domain.dim(0).low;        // First element is density
    var idxMom  : range = consL.domain.dim(0).expand(-1); // Intermediary elements are the velocities
    var idxEner : int   = consL.domain.dim(0).high;       // Last element is energy

    var uniNrm : [nrm.domain] real = nrm/norm(nrm, normType.norm2);

    // Primitive and other variables.
    //   Left state
    ref densL   : real = consL[idxDens];
    var velL    : [idxMom-1] real = consL[idxMom]/densL;
    ref enerL   : real = consL[idxEner];

    var velNrmL : real = dot[velL, uniNrm];
    var presL   : real = pressure_cv(consL);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL;

    //   Right state
    ref densR   : real = consR[idxDens];
    var velR    : [idxMom-1] real = consR[idxMom]/densR;
    ref enerR   : real = consR[idxEner];

    var velNrmR : real = dot[velR, uniNrm];
    var presR   : real = pressure_cv(consR);
    var aR      : real = sqrt(fGamma*presR/densR);
    var enthR   : real = ( enerR + presR ) / densR;

    // Compute the Roe Averages
    var RT     : real = sqrt( densR/densL );
    var dens   : real = RT*densL;
    var vel    : [idxMom-1] real = (velL + RT*velR)/(1.0 + RT);
    var velNrm : real = dot(vel, uniNrm);
    var enth   : real = (enthL + RT*enthR)/(1.0 + RT);
    var a      : real = sqrt( (fGamma-1.0)*(enth-0.5*dot(vel, vel)) );

    // Differences in primitive variables
    var dDens   = densR   - densL;
    var dVelNrm = velNrmR - velNrmL;
    var dPres   = presR   - presL;
    var dVel    = velR    - velL;

    // Wave strengths (Characteristic Variables).
    var dV : [1..4] real;
    dV[1] =  (dPres - dens*a*dVelNrm)/(2.0*a**2);
    dV[2] =   dDens - dPres/(a**2);
    dV[3] =  (dPres + dens*a*dVelNrm)/(2.0*a**2);
    dV[4] =   dens;

    // Absolute values of the wave speeds (Eigenvalues)
    var ws : [1..4] real;
    ws[1] = abs(velNrm - a);
    ws[2] = abs(velNrm    );
    ws[3] = abs(velNrm + a);
    ws[4] = abs(velNrm    );

    // Harten's Entropy Fix JCP(1983), 49, pp357-393. Only for the nonlinear fields.
    // It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
    if ( ws[1] < 0.2 ) then ws[1] = (ws[1]**2)/0.4 + 0.1;
    if ( ws[3] < 0.2 ) then ws[3] = (ws[3]**2)/0.4 + 0.1;

    // Right eigenvectors
    var R : [consL.domain.dim(0), 1..4] real;

    // Negative acoustic wave
    R[idxDens, 1] = 1.0;
    R[idxMom , 1] = vel  - a*uniNrm;
    R[idxEner, 1] = enth - a*velNrm;

    // Entropy wave
    R[idxDens, 2] = 1.0;
    R[idxMom , 2] = vel;
    R[idxEner, 2] = dot(vel, vel)/2.0;

    // Positive acoustic wave
    R[idxDens, 3] = 1.0;
    R[idxMom , 3] = vel  + a*uniNrm;
    R[idxEner, 3] = enth + a*velNrm;

    // Combined shear waves
    R[idxDens, 4] = 0.0;
    R[idxMom , 4] = dVel - dVelNrm*uniNrm;
    R[idxEner, 4] = dot(vel, dVel) - velNrm*dVelNrm;

    // Compute the let and right fluxes
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2] + ws[3]*dV[3]*R[.., 3];
    if idxEner > 3 then diss += ws[4]*dV[4]*R[.., 4];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc hll(consL : [] real, consR : [] real, nrm : [] real) : [] real
  {
    // The HLL (Harten-Lax-van Leer)
  }

  proc hllc(consL : [] real, consR : [] real, nrm : [] real) : [] real
  {
    // The HLLC (Harten-Lax-van Leer-Contact) solver was introduced by Toro. It restores the missing Rarefaction wave by
    // some estimates, like linearisations, these can be simple but also more advanced exists like using the Roe average
    // velocity for the middle wave speed. They are quite robust and efficient but somewhat more diffusive.
    //
    // Toro, E. F.; Spruce, M.; Speares, W. (1994), "Restoration of the contact surface in the HLL-Riemann solver"
    // Shock Waves, 4 (1): 25–34, doi:10.1007/BF01414629
    //
    // Quirk, J. J. (1994), "A contribution to the great Riemann solver debate"
    // Int. J. Numer. Methods Fluids, 18 (6): 555–574, doi:10.1002/fld.1650180603
  }

  proc rotated_rhll(consL : [] real, consR : [] real, nrm : [] real) : [] real
  {
    // These solvers were introduced by Nishikawa and Kitamura, in order to overcome the carbuncle problems of the Roe
    // solver and the excessive diffusion of the HLLE solver at the same time. They developed robust and accurate Riemann
    // solvers by combining the Roe solver and the HLLE/Rusanov solvers: they show that being applied in two orthogonal
    // directions the two Riemann solvers can be combined into a single Roe type solver (the Roe solver with modified
    // wave speeds). In particular, the one derived from the Roe and HLLE solvers, called Rotated-RHLL solver, is
    // extremely robust (carbuncle free for all possible test cases on both structured and unstructured grids) and
    // accurate (as accurate as the Roe solver for the boundary layer calculation).
    //
    // Nishikawa, H.; Kitamura, K. (2008), "Very simple, carbuncle-free, boundary-layer-resolving, rotated-hybrid Riemann solvers"
    // J. Comput. Phys., 227 (4): 2560–2581, doi:10.1016/j.jcp.2007.11.003
  }

  proc godunov(consL : [] real, consR : [] real, nrm : [] real) : [] real
  {
  }

  proc main()
  {
    use LinearAlgebra;
    import Flux.euler_flux_cv;

    var nrm1 : [1..3] real = [+0.612372436, +0.353553391, +0.707106781];
    var nrm2 : [1..3] real = [-0.612372436, -0.353553391, -0.707106781];


    writeln();
    writeln("1D Tests");
    {
      var cons1dL : [1..3] real = [1.225, 250.1160830494510, 278846.40];
      var cons1dR : [1..3] real = [1.125, 265.2881676121270, 259260.30];

      var uniNrm1 : [1..1] real = nrm1[1..1]/norm(nrm1[1..1]);
      var uniNrm2 : [1..1] real = nrm2[1..1]/norm(nrm2[1..1]);

      writef("  Conserved variables:\n");
      writef("    Left  %t: %.4ht\n", cons1dL.domain, cons1dL);
      writef("    Right %t: %.4ht\n", cons1dR.domain, cons1dR);

      writef("\n");
      writef("  Face normal unit vectors:\n");
      writef("    Normal 1 %t: %.4ht\n", nrm1[1..1].domain, nrm1[1..1]);
      writef("    Normal 2 %t: %.4ht\n", nrm2[1..1].domain, nrm2[1..1]);

      writef("\n");
      writef("  Euler flux:\n");
      writef("    Left  %t: %.4ht\n", euler_flux_cv(cons1dL).domain,
                                      euler_flux_cv(cons1dL)[1,..] );
      writef("    Right %t: %.4ht\n", euler_flux_cv(cons1dR).domain,
                                      euler_flux_cv(cons1dR)[1,..] );

      writef("\n");
      writef("  Normal 1:\n");
      writef("    Left  Euler Flux %t: %.4ht\n", dot( nrm1[1..1], euler_flux_cv(cons1dL)).domain,
                                                 dot( nrm1[1..1], euler_flux_cv(cons1dL))       );
      writef("    Right Euler Flux %t: %.4ht\n", dot( nrm1[1..1], euler_flux_cv(cons1dR)).domain,
                                                 dot( nrm1[1..1], euler_flux_cv(cons1dR))       );
      writef("    1D  Roe     Flux %t: %.4ht\n", roe_1d(cons1dR, cons1dL, nrm1[1..1]    ).domain,
                                                 roe_1d(cons1dR, cons1dL, nrm1[1..1]    )       );
      writef("    Gen Roe 1D  Flux %t: %.4ht\n", roe(   cons1dR, cons1dL, nrm1[1..1]    ).domain,
                                                 roe(   cons1dR, cons1dL, nrm1[1..1]    )       );

      writef("  Normal 2:\n");
      writef("    Left  Euler Flux %t: %.4ht\n", dot( nrm2[1..1], euler_flux_cv(cons1dL)).domain,
                                                 dot( nrm2[1..1], euler_flux_cv(cons1dL))       );
      writef("    Right Euler Flux %t: %.4ht\n", dot( nrm2[1..1], euler_flux_cv(cons1dR)).domain,
                                                 dot( nrm2[1..1], euler_flux_cv(cons1dR))       );
      writef("    1D  Roe     Flux %t: %.4ht\n", roe_1d(cons1dR, cons1dL, nrm2[1..1]    ).domain,
                                                 roe_1d(cons1dR, cons1dL, nrm2[1..1]    )       );
      writef("    Gen Roe 1D  Flux %t: %.4ht\n", roe(   cons1dR, cons1dL, nrm2[1..1]    ).domain,
                                                 roe(   cons1dR, cons1dL, nrm2[1..1]    )       );

      writef("\n");
      writef("  Normal 1:\n");
      writef("    Left  Euler Flux %t: %.4ht\n", outer(nrm1[1..1], dot( nrm1[1..1], euler_flux_cv(cons1dL))).domain,
                                                 outer(nrm1[1..1], dot( nrm1[1..1], euler_flux_cv(cons1dL)))[1, ..]);
      writef("    Right Euler Flux %t: %.4ht\n", outer(nrm1[1..1], dot( nrm1[1..1], euler_flux_cv(cons1dR))).domain,
                                                 outer(nrm1[1..1], dot( nrm1[1..1], euler_flux_cv(cons1dR)))[1, ..]);
      writef("    1D  Roe     Flux %t: %.4ht\n", outer(nrm1[1..1], roe_1d(cons1dR, cons1dL, nrm1[1..1]    )).domain,
                                                 outer(nrm1[1..1], roe_1d(cons1dR, cons1dL, nrm1[1..1]    ))[1, ..]);
      writef("    Gen Roe 1D  Flux %t: %.4ht\n", outer(nrm1[1..1], roe(   cons1dR, cons1dL, nrm1[1..1]    )).domain,
                                                 outer(nrm1[1..1], roe(   cons1dR, cons1dL, nrm1[1..1]    ))[1, ..]);

      writef("  Normal 2:\n");
      writef("    Left  Euler Flux %t: %.4ht\n", outer(nrm2[1..1], dot( nrm2[1..1], euler_flux_cv(cons1dL))).domain,
                                                 outer(nrm2[1..1], dot( nrm2[1..1], euler_flux_cv(cons1dL)))[1, ..]);
      writef("    Right Euler Flux %t: %.4ht\n", outer(nrm2[1..1], dot( nrm2[1..1], euler_flux_cv(cons1dR))).domain,
                                                 outer(nrm2[1..1], dot( nrm2[1..1], euler_flux_cv(cons1dR)))[1, ..]);
      writef("    1D  Roe     Flux %t: %.4ht\n", outer(nrm2[1..1], roe_1d(cons1dL, cons1dR, nrm2[1..1]    )).domain,
                                                 outer(nrm2[1..1], roe_1d(cons1dL, cons1dR, nrm2[1..1]    ))[1, ..]);
      writef("    Gen Roe 1D  Flux %t: %.4ht\n", outer(nrm2[1..1], roe(   cons1dL, cons1dR, nrm2[1..1]    )).domain,
                                                 outer(nrm2[1..1], roe(   cons1dL, cons1dR, nrm2[1..1]    ))[1, ..]);
    }

    writeln();
    writeln("2D Tests");
    {
      var cons2dL : [1..4] real = [1.225, 249.1643158413380, 21.79905299130990, 278846.40];
      var cons2dR : [1..4] real = [1.125, 264.2786660416750, 23.12138729040020, 259260.30];

      var uniNrm1 : [1..2] real = nrm1[1..2]/norm(nrm1[1..2]);
      var uniNrm2 : [1..2] real = nrm2[1..2]/norm(nrm2[1..2]);

      writef("  Conserved variables:\n");
      writef("    Left  %t: %.4ht\n", cons2dL.domain, cons2dL);
      writef("    Right %t: %.4ht\n", cons2dR.domain, cons2dR);

      writef("\n");
      writef("  Face normal unit vectors:\n");
      writef("    Normal 1 %t: %.4ht\n", nrm1[1..2].domain, nrm1[1..2]);
      writef("    Normal 2 %t: %.4ht\n", nrm2[1..2].domain, nrm2[1..2]);

      writef("\n");
      writef("  Euler flux:\n");
      writef("    Left  %t: %.4ht, %.4ht\n", euler_flux_cv(cons2dL).domain,
                                             euler_flux_cv(cons2dL)[1,..] ,
                                             euler_flux_cv(cons2dL)[2,..] );
      writef("    Right %t: %.4ht, %.4ht\n", euler_flux_cv(cons2dR).domain,
                                             euler_flux_cv(cons2dR)[1,..] ,
                                             euler_flux_cv(cons2dR)[2,..] );

      writef("\n");
      writef("  Normal 1:\n");
      writef("    Left  Euler Flux %t: %.4ht\n", dot( nrm1[1..2], euler_flux_cv(cons2dL)).domain,
                                                 dot( nrm1[1..2], euler_flux_cv(cons2dL))       );
      writef("    Right Euler Flux %t: %.4ht\n", dot( nrm1[1..2], euler_flux_cv(cons2dR)).domain,
                                                 dot( nrm1[1..2], euler_flux_cv(cons2dR))       );
      writef("    2D  Roe     Flux %t: %.4ht\n", roe_2d(cons2dR, cons2dL, nrm1[1..2]    ).domain,
                                                 roe_2d(cons2dR, cons2dL, nrm1[1..2]    )       );
      writef("    Gen Roe 2D  Flux %t: %.4ht\n", roe(   cons2dR, cons2dL, nrm1[1..2]    ).domain,
                                                 roe(   cons2dR, cons2dL, nrm1[1..2]    )       );

      writef("\n");
      writef("  Normal 2:\n");
      writef("    Left  Euler Flux %t: %.4ht\n", dot( nrm2[1..2], euler_flux_cv(cons2dL)).domain,
                                                 dot( nrm2[1..2], euler_flux_cv(cons2dL))       );
      writef("    Right Euler Flux %t: %.4ht\n", dot( nrm2[1..2], euler_flux_cv(cons2dR)).domain,
                                                 dot( nrm2[1..2], euler_flux_cv(cons2dR))       );
      writef("    2D  Roe     Flux %t: %.4ht\n", roe_2d(cons2dL, cons2dR, nrm2[1..2]    ).domain,
                                                 roe_2d(cons2dL, cons2dR, nrm2[1..2]    )       );
      writef("    Gen Roe 2D  Flux %t: %.4ht\n", roe(   cons2dL, cons2dR, nrm2[1..2]    ).domain,
                                                 roe(   cons2dL, cons2dR, nrm2[1..2]    )       );

      writef("\n");
      writef("  Normal 1:\n");
      writef("    Left  Euler Flux %t: %.4ht, %.4ht\n", outer(nrm1[1..2], dot( nrm1[1..2], euler_flux_cv(cons2dL))).domain,
                                                        outer(nrm1[1..2], dot( nrm1[1..2], euler_flux_cv(cons2dL)))[1, ..],
                                                        outer(nrm1[1..2], dot( nrm1[1..2], euler_flux_cv(cons2dL)))[2, ..]);
      writef("    Right Euler Flux %t: %.4ht, %.4ht\n", outer(nrm1[1..2], dot( nrm1[1..2], euler_flux_cv(cons2dR))).domain,
                                                        outer(nrm1[1..2], dot( nrm1[1..2], euler_flux_cv(cons2dR)))[1, ..],
                                                        outer(nrm1[1..2], dot( nrm1[1..2], euler_flux_cv(cons2dR)))[2, ..]);
      writef("    2D  Roe     Flux %t: %.4ht, %.4ht\n", outer(nrm1[1..2], roe_2d(cons2dR, cons2dL, nrm1[1..2]    )).domain,
                                                        outer(nrm1[1..2], roe_2d(cons2dR, cons2dL, nrm1[1..2]    ))[1, ..],
                                                        outer(nrm1[1..2], roe_2d(cons2dR, cons2dL, nrm1[1..2]    ))[2, ..]);
      writef("    Gen Roe 2D  Flux %t: %.4ht, %.4ht\n", outer(nrm1[1..2], roe(   cons2dR, cons2dL, nrm1[1..2]    )).domain,
                                                        outer(nrm1[1..2], roe(   cons2dR, cons2dL, nrm1[1..2]    ))[1, ..],
                                                        outer(nrm1[1..2], roe(   cons2dR, cons2dL, nrm1[1..2]    ))[2, ..]);
      writef("  Normal 2:\n");
      writef("    Left  Euler Flux %t: %.4ht, %.4ht\n", outer(nrm2[1..2], dot( nrm2[1..2], euler_flux_cv(cons2dL))).domain,
                                                        outer(nrm2[1..2], dot( nrm2[1..2], euler_flux_cv(cons2dL)))[1, ..],
                                                        outer(nrm2[1..2], dot( nrm2[1..2], euler_flux_cv(cons2dL)))[2, ..]);
      writef("    Right Euler Flux %t: %.4ht, %.4ht\n", outer(nrm2[1..2], dot( nrm2[1..2], euler_flux_cv(cons2dR))).domain,
                                                        outer(nrm2[1..2], dot( nrm2[1..2], euler_flux_cv(cons2dR)))[1, ..],
                                                        outer(nrm2[1..2], dot( nrm2[1..2], euler_flux_cv(cons2dR)))[2, ..]);
      writef("    2D  Roe     Flux %t: %.4ht, %.4ht\n", outer(nrm2[1..2], roe_2d(cons2dL, cons2dR, nrm2[1..2]    )).domain,
                                                        outer(nrm2[1..2], roe_2d(cons2dL, cons2dR, nrm2[1..2]    ))[1, ..],
                                                        outer(nrm2[1..2], roe_2d(cons2dL, cons2dR, nrm2[1..2]    ))[2, ..]);
      writef("    Gen Roe 2D  Flux %t: %.4ht, %.4ht\n", outer(nrm2[1..2], roe(   cons2dL, cons2dR, nrm2[1..2]    )).domain,
                                                        outer(nrm2[1..2], roe(   cons2dL, cons2dR, nrm2[1..2]    ))[1, ..],
                                                        outer(nrm2[1..2], roe(   cons2dL, cons2dR, nrm2[1..2]    ))[2, ..]);
    }

    writeln();
    writeln("3D Tests");
    {
      var cons3dL : [1..5] real = [1.225, 249.0125316723220, 17.79514529902850, 7.098538138029130, 278846.40];
      var cons3dR : [1..5] real = [1.125, 264.1176746188930, 23.12138729040020, 9.223192434062800, 259260.30];

      var uniNrm1 : [1..3] real = nrm1[1..3]/norm(nrm1[1..3]);
      var uniNrm2 : [1..3] real = nrm2[1..3]/norm(nrm2[1..3]);

      writef("  Conserved variables:\n");
      writef("    Left  %t: %.4ht\n", cons3dL.domain, cons3dL);
      writef("    Right %t: %.4ht\n", cons3dR.domain, cons3dR);

      writef("\n");
      writef("  Face normal unit vectors:\n");
      writef("    Normal 1 %t: %.4ht\n", nrm1[1..3].domain, nrm1[1..3]);
      writef("    Normal 2 %t: %.4ht\n", nrm2[1..3].domain, nrm2[1..3]);

      writef("\n");
      writef("  Euler flux:\n");
      writef("    Left  %t: %.4ht, %.4ht, %.4ht\n", euler_flux_cv(cons3dL).domain,
                                                    euler_flux_cv(cons3dL)[1,..] ,
                                                    euler_flux_cv(cons3dL)[2,..] ,
                                                    euler_flux_cv(cons3dL)[3,..] );
      writef("    Right %t: %.4ht, %.4ht, %.4ht\n", euler_flux_cv(cons3dR).domain,
                                                    euler_flux_cv(cons3dR)[1,..] ,
                                                    euler_flux_cv(cons3dR)[2,..] ,
                                                    euler_flux_cv(cons3dR)[3,..] );

      writef("\n");
      writeln("  Normal 1:");
      writef("    Left  Euler Flux %t: %.4ht\n", dot( nrm1[1..3], euler_flux_cv(cons3dL)).domain,
                                                 dot( nrm1[1..3], euler_flux_cv(cons3dL))       );
      writef("    Right Euler Flux %t: %.4ht\n", dot( nrm1[1..3], euler_flux_cv(cons3dR)).domain,
                                                 dot( nrm1[1..3], euler_flux_cv(cons3dR))       );
      writef("    2D  Roe     Flux %t: %.4ht\n", roe_3d(cons3dL, cons3dR, nrm1[1..3]    ).domain,
                                                 roe_3d(cons3dL, cons3dR, nrm1[1..3]    )       );
      writef("    Gen Roe 2D  Flux %t: %.4ht\n", roe(   cons3dL, cons3dR, nrm1[1..3]    ).domain,
                                                 roe(   cons3dL, cons3dR, nrm1[1..3]    )       );

      writeln("  Normal 2:");
      writef("    Left  Euler Flux %t: %.4ht\n", dot( nrm2[1..3], euler_flux_cv(cons3dL)).domain,
                                                 dot( nrm2[1..3], euler_flux_cv(cons3dL))       );
      writef("    Right Euler Flux %t: %.4ht\n", dot( nrm2[1..3], euler_flux_cv(cons3dR)).domain,
                                                 dot( nrm2[1..3], euler_flux_cv(cons3dR))       );
      writef("    2D  Roe     Flux %t: %.4ht\n", roe_3d(cons3dL, cons3dR, nrm2[1..3]    ).domain,
                                                 roe_3d(cons3dL, cons3dR, nrm2[1..3]    )       );
      writef("    Gen Roe 2D  Flux %t: %.4ht\n", roe(   cons3dL, cons3dR, nrm2[1..3]    ).domain,
                                                 roe(   cons3dL, cons3dR, nrm2[1..3]    )       );

      writef("\n");
      writeln("  Normal 1:");
      writef("    Left  Euler Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL))).domain,
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL)))[1, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL)))[2, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL)))[3, ..]);
      writef("    Right Euler Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR))).domain,
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR)))[1, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR)))[2, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR)))[3, ..]);
      writef("    3D  Roe     Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    )).domain,
                                                               outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    ))[1, ..],
                                                               outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    ))[2, ..],
                                                               outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    ))[3, ..]);
      writef("    Gen Roe 3D  Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    )).domain,
                                                               outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    ))[1, ..],
                                                               outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    ))[2, ..],
                                                               outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    ))[3, ..]);

      writeln("  Normal 2:");
      writef("    Left  Euler Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL))).domain,
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL)))[1, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL)))[2, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dL)))[3, ..]);
      writef("    Right Euler Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR))).domain,
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR)))[1, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR)))[2, ..],
                                                               outer(nrm1[1..3], dot( nrm1[1..3], euler_flux_cv(cons3dR)))[3, ..]);
      writef("    3D  Roe     Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    )).domain,
                                                               outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    ))[1, ..],
                                                               outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    ))[2, ..],
                                                               outer(nrm1[1..3], roe_3d(cons3dR, cons3dL, nrm1[1..3]    ))[3, ..]);
      writef("    Gen Roe 3D  Flux %t: %.4ht, %.4ht, %.4ht\n", outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    )).domain,
                                                               outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    ))[1, ..],
                                                               outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    ))[2, ..],
                                                               outer(nrm1[1..3], roe(   cons3dR, cons3dL, nrm1[1..3]    ))[3, ..]);
    }
  }
}
