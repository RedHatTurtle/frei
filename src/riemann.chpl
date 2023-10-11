module Riemann
{
  proc upwind_1d(const ref consL : [1..3] real, const ref consR : [1..3] real, const ref nrm : [1..1] real) : [1..3] real
  {
    use Flux;
    import Input.convectionSpeed;

    var upwind : [consL.domain] real;

    if nrm[1] < 0 then
      upwind = convection_flux_cv_1d(consR, convectionSpeed);
    else
      upwind = convection_flux_cv_1d(consL, convectionSpeed);

    return upwind;
  }

  proc rusanov_1d(const ref consL : [1..3] real, const ref consR : [1..3] real, const ref nrm : [1..1] real) : [1..3] real
  {
    import Input.fGamma;
    import Flux.pressure_cv;
    import Flux.euler_flux_cv;

    var idxDens : int   = consL.domain.dim(0).low;           // First element is density
    var idxMom  : range = consL.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEner : int   = consL.domain.dim(0).high;          // Last element is energy

    // Primitive variables
    //  Left state
    var densL : real = consL[idxDens];
    var velL  : real = consL[idxMom]/consL[idxDens];
    var presL : real = pressure_cv(consL, fGamma);
    var enthL : real = ( consL[idxEner] + presL ) / densL;
    //  Right state
    var densR : real = consR[idxDens];
    var velR  : real = consR[idxMom]/consR[idxDens];
    var presR : real = pressure_cv(consR, fGamma);
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

  proc roe_1d(const ref consL : [1..3] real, const ref consR : [1..3] real, const ref nrm : [1..1] real) : [1..3] real
  {
    use LinearAlgebra;
    import Input.fGamma;
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
    var presL   : real = pressure_cv(consL, fGamma);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( consL[3] + presL ) / densL;
    //  Right state
    var densR   : real = consR[1];
    var velR    : real = consR[2]/consR[1];
    var velNrmR : real = velR * uniNrm[1];
    var presR   : real = pressure_cv(consR, fGamma);
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
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL, fGamma));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR, fGamma));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2] + ws[3]*dV[3]*R[.., 3];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc roe_2d(const ref consL : [1..4] real, const ref consR : [1..4] real, const ref nrm : [1..2] real) : [1..4] real
  {
    use LinearAlgebra;
    import Input.fGamma;
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
    var presL   : real = pressure_cv(consL, fGamma);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL;

    //  Right state
    ref densR   : real = consR[idxDens];
    ref enerR   : real = consR[idxEner];
    var velR    : [idxMom-1] real = consR[idxMom]/densR;
    var velNrmR : real = dot[velR, uniNrm];
    var presR   : real = pressure_cv(consR, fGamma);
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
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL, fGamma));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR, fGamma));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2]
                                   + ws[3]*dV[3]*R[.., 3] + ws[4]*dV[4]*R[.., 4];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc roe_3d(const ref consL : [1..5] real, const ref consR : [1..5] real, const ref nrm : [1..3] real) : [1..5] real
  {
    use LinearAlgebra;
    import Input.fGamma;
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
    var presL   : real = pressure_cv(consL, fGamma);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL;

    //  Right state
    ref densR   : real = consR[idxDens];
    ref enerR   : real = consR[idxEner];
    var velR    : [idxMom-1] real = consR[idxMom]/densR;
    var velNrmR : real = dot[velR, uniNrm];
    var presR   : real = pressure_cv(consR, fGamma);
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
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL, fGamma));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR, fGamma));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2]
                                   + ws[3]*dV[3]*R[.., 3] + ws[4]*dV[4]*R[.., 4];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc rotated_rhll_2d(const ref consL : [1..4] real, const ref consR : [1..4] real, const ref nrm : [1..2] real) : [1..4] real
  {
    use LinearAlgebra;
    import Input.fGamma;
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
    var velVL   : [idxMom-1] real = consL[idxMom]/densL;
    var velNrmL : real = dot[velVL, uniNrm];
    var presL   : real = pressure_cv(consL, fGamma);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL;

    //  Right state
    ref densR   : real = consR[idxDens];
    ref enerR   : real = consR[idxEner];
    var velVR   : [idxMom-1] real = consR[idxMom]/densR;
    var velNrmR : real = dot[velVR, uniNrm];
    var presR   : real = pressure_cv(consR, fGamma);
    var aR      : real = sqrt(fGamma*presR/densR);
    var enthR   : real = ( enerR + presR ) / densR;

    // Compute the Roe Averages
    var RT     : real = sqrt( densR/densL );
    var dens   : real = RT*densL;
    var velV   : [idxMom-1] real = (velVL + RT*velVR)/(1.0 + RT);
    var enth   : real = (enthL + RT*enthR)/(1.0 + RT);
    var a      : real = sqrt( (fGamma-1.0)*(enth-0.5*dot(velV, velV)) );

    //--------------------------------------------------------------------------------
    // Define nrmA and nrmB, and compute alpha1 and alpha2: (4.2) in the original paper.
    //
    // Note: nrma and nrmB may need to be frozen at some point during a steady calculation to fully make it converge.
    //       For time-accurate calculation, it is not necessary.
    // Note: For a boundary face, you may want to set nrmB = nrm and nrmA = Tangent to force use the Roe flux.

    // Calculate the magnitude of the velocity difference vector
    var dVelV    = velVR-velVL;
    var absDVelV = norm(dVelV, normType.norm2);

    // If it's large enough it can be used as the first normal
    var uniNrmA : [nrm.domain] real;
    if absDVelV > 1.0e-12 {
        // Normal to shock or tangent to shear
        uniNrmA = dVelV/absDVelV;
    } else {
        // A face tangent vector is used rotating the face normal
        uniNrmA[1] = -uniNrm[2];
        uniNrmA[2] =  uniNrm[1];
    }

    // Calculate the angle between the face normal and nrmA
    var alphaA = dot(uniNrm, uniNrmA);

    // Make alphaA always positive.
    if sgn(alphaA) < 0 {
        uniNrmA = sgn(alphaA) * uniNrmA;
        alphaA  = sgn(alphaA) * alphaA;
    }

    //----------------------------------------------------
    // Compute the wave speed estimates for the HLL part, following Einfeldt:
    // B. Einfeldt, On Godunov-type methods for gas dynamics, SIAM Journal on Numerical Analysis 25 (2) (1988) 294–318.
    //
    // Note: HLL is actually applied to n1, but this is all we need to incorporate HLL. See JCP2008 paper.

    var velNrmA  = dot(velV, uniNrmA);
    var SRp = max( 0.0, velNrmA + a, dot(velVR, uniNrmA) + aR ); // Maximum wave speed estimate
    var SLm = min( 0.0, velNrmA - a, dot(velVL, uniNrmA) - aL ); // Minimum wave speed estimate

    // This is the only place where uniNrmA=(nx1,ny1,nz1), except from calculating uniNrmB right bellow.
    //----------------------------------------------------

    // nrmB = direction perpendicular to nrmA.
    //     Note: There are infinitely many choices for this vector. The best choice may be discovered in future.
    // The paper authors suggest employing the formula (4.4) in the paper:
    //     NrmB = (nrmAxn)xnrmA / |(n1xn)xn1|    ('x' is the vector product.)

    // Author's implementation/suggestion
    var uniNrmB : [nrm.domain] real;
    uniNrmB[1] = -uniNrmA[2];
    uniNrmB[2] =  uniNrmA[1];
    uniNrmB = uniNrmB/norm(uniNrmB, normType.norm2);

    // Calculate the angle between the face normal and nrmB
    var alphaB = dot(uniNrm, uniNrmB);

    // Make alphaB always positive.
    if sgn(alphaB) < 0 {
        uniNrmB = sgn(alphaB) * uniNrmB;
        alphaB  = sgn(alphaB) * alphaB;
    }

    // Now we are going to compute the Roe flux with n2 as the normal with modified wave speeds (5.12).
    // NOTE: The Roe flux here is computed without tangent vectors. See "I do like CFD, VOL.1" for details: page 57,
    //       Equation (3.6.31).

    // Recalculate normal velocities using uniNrmB
    var velNrmBL = dot(velVL, uniNrmB);
    var velNrmBR = dot(velVR, uniNrmB);
    var velNrmB  = dot(velV , uniNrmB);

    // Differences in primitive variables
    var dDens    = densR    - densL;
    var dPres    = presR    - presL;
    var dVelNrmB = velNrmBR - velNrmBL;

    // Wave strengths (Characteristic Variables).
    var dV : [1..4] real;
    dV[1] =  (dPres - dens*a*dVelNrmB)/(2.0*a**2);
    dV[2] =   dDens - dPres/(a**2);
    dV[3] =  (dPres + dens*a*dVelNrmB)/(2.0*a**2);
    dV[4] =   dens;

    // Absolute values of the wave speeds (Eigenvalues)
    var eig : [1..4] real;
    eig[1] = velNrmB - a; // Left moving acoustic wave
    eig[2] = velNrmB;     // Entropy wave
    eig[3] = velNrmB + a; // Right moving acoustic wave
    eig[4] = velNrmB;     // Shear wave

    var ws : [1..4] real = abs(eig);

    // Harten's Entropy Fix JCP(1983), 49, pp357-393. Only for the nonlinear fields.
    // It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
    if ( ws[1] < 0.2 ) then ws[1] = (ws[1]**2)/0.4 + 0.1;
    if ( ws[3] < 0.2 ) then ws[3] = (ws[3]**2)/0.4 + 0.1;

    // Combine the wave speeds for Rotated-RHLL: Eq.(5.12) in the original JCP2008 paper.
    ws = ws*alphaB - (2.0*alphaA*SRp*SLm + alphaB*(SRp+SLm)*eig)/(SRp-SLm);

    // Right eigenvectors
    var R : [1..5, 1..4] real;

    // Negative acoustic wave
    R[idxDens, 1] = 1.0;
    R[idxMom , 1] = velV - a*uniNrmB;
    R[idxEner, 1] = enth - a*velNrmB;

    // Entropy wave
    R[idxDens, 2] = 1.0;
    R[idxMom , 2] = velV;
    R[idxEner, 2] = dot(velV, velV)/2.0;

    // Positive acoustic wave
    R[idxDens, 3] = 1.0;
    R[idxMom , 3] = velV + a*uniNrmB;
    R[idxEner, 3] = enth + a*velNrmB;

    // Combined shear waves
    R[idxDens, 4] = 0.0;
    R[idxMom , 4] = dVelV - dVelNrmB*uniNrmB;
    R[idxEner, 4] = dot(velV, dVelV) - velNrmB*dVelNrmB;

    // Compute the let and right fluxes
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL, fGamma));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR, fGamma));

    // Calculate the dissipation term
    var roeDiss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2]
                                      + ws[3]*dV[3]*R[.., 3] + ws[4]*dV[4]*R[.., 4];

    // Compute the final 2D Rotated-RHLL flux
    var rotRHLL : [consL.domain] real = (fluxL*SRp - fluxR*SLm)/(SRp-SLm) - roeDiss/2.0;

    return rotRHLL;
  }

  proc rotated_rhll_3d(const ref consL : [1..5] real, const ref consR : [1..5] real, const ref nrm : [1..3] real) : [1..5] real
  {
    use LinearAlgebra;
    import Input.fGamma;
    import Flux.pressure_cv;
    import Flux.euler_flux_cv;

    var idxDens : int   = consL.domain.dim(0).low;        // First element is density
    var idxMom  : range = consL.domain.dim(0).expand(-1); // Intermediary elements are the velocities
    var idxEner : int   = consL.domain.dim(0).high;       // Last element is energy

    var uniNrm : [1..3] real = nrm/norm(nrm, normType.norm2);

    //Primitive and other variables.
    //  Left state
    ref densL   : real = consL[idxDens];
    ref enerL   : real = consL[idxEner];
    var velVL   : [idxMom-1] real = consL[idxMom]/densL;
    var velNrmL : real = dot[velVL, uniNrm];
    var presL   : real = pressure_cv(consL, fGamma);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL;

    //  Right state
    ref densR   : real = consR[idxDens];
    ref enerR   : real = consR[idxEner];
    var velVR   : [idxMom-1] real = consR[idxMom]/densR;
    var velNrmR : real = dot[velVR, uniNrm];
    var presR   : real = pressure_cv(consR, fGamma);
    var aR      : real = sqrt(fGamma*presR/densR);
    var enthR   : real = ( enerR + presR ) / densR;

    // Compute the Roe Averages
    var RT     : real = sqrt( densR/densL );
    var dens   : real = RT*densL;
    var velV   : [idxMom-1] real = (velVL + RT*velVR)/(1.0 + RT);
    var enth   : real = (enthL + RT*enthR)/(1.0 + RT);
    var a      : real = sqrt( (fGamma-1.0)*(enth-0.5*dot(velV, velV)) );

    //--------------------------------------------------------------------------------
    // Define nrmA and nrmB, and compute alpha1 and alpha2: (4.2) in the original paper.
    //
    // Note: nrma and nrmB may need to be frozen at some point during a steady calculation to fully make it converge.
    //       For time-accurate calculation, it is not necessary.
    // Note: For a boundary face, you may want to set nrmB = nrm and nrmA = Tangent to force use the Roe flux.

    // Calculate the magnitude of the velocity difference vector
    var dVelV    = velVR-velVL;
    var absDVelV = norm(dVelV, normType.norm2);

    // If it's large enough it can be used as the first normal
    var uniNrmA : [1..3] real;
    if absDVelV > 1.0e-12 then
        // Normal to shock or tangent to shear
        uniNrmA = dVelV/absDVelV;
    else {
        // A face tangent vector is used rotating the two largest components of uniNrm
        var nx2 = uniNrm[1]*uniNrm[1];
        var ny2 = uniNrm[2]*uniNrm[2];
        var nz2 = uniNrm[3]*uniNrm[3];
        if (nx2 <= ny2 && nx2 <= nz2) {
            uniNrmA[1] =  uniNrm[1];
            uniNrmA[1] =  0.0;          // Authors' implementation
            uniNrmA[2] = -uniNrm[3];
            uniNrmA[3] =  uniNrm[2];
        } else if (ny2 <= nx2 && ny2 <= nz2) {
            uniNrmA[1] = -uniNrm[3];
            uniNrmA[2] =  uniNrm[2];
            uniNrmA[2] =  0.0;          // Authors' implementation
            uniNrmA[3] =  uniNrm[1];
        } else if (nz2 <= nx2 && nz2 <= ny2) {
            uniNrmA[1] = -uniNrm[2];
            uniNrmA[2] =  uniNrm[1];
            uniNrmA[3] =  uniNrm[3];
            uniNrmA[3] =  0.0;          // Authors' implementation
        } else {
            writeln("Error defining Normal A in the Rotated-RHLL Riemann Flux");
        }
    }

    // Calculate the angle between the face normal and nrmA
    var alphaA = dot(uniNrm, uniNrmA);

    // Make alphaA always positive.
    if sgn(alphaA) < 0 {
        uniNrmA = sgn(alphaA) * uniNrmA;
        alphaA  = sgn(alphaA) * alphaA;
    }

    //----------------------------------------------------
    // Compute the wave speed estimates for the HLL part, following Einfeldt:
    // B. Einfeldt, On Godunov-type methods for gas dynamics, SIAM Journal on Numerical Analysis 25 (2) (1988) 294–318.
    //
    // Note: HLL is actually applied to n1, but this is all we need to incorporate HLL. See JCP2008 paper.

    var velNrmA  = dot(velV, uniNrmA);
    var SLm = min( 0.0, velNrmA - a, dot(velVL, uniNrmA) - aL ); // Minimum wave speed estimate
    var SRp = max( 0.0, velNrmA + a, dot(velVR, uniNrmA) + aR ); // Maximum wave speed estimate

    // This is the only place where uniNrmA=(nx1,ny1,nz1), except from calculating uniNrmB right bellow.
    //----------------------------------------------------

    // nrmB = direction perpendicular to nrmA.
    //     Note: There are infinitely many choices for this vector. The best choice may be discovered in future.
    // The paper authors suggest employing the formula (4.4) in the paper:
    //     NrmB = (nrmAxn)xnrmA / |(n1xn)xn1|    ('x' is the vector product.)

    //var uniNrmB : [nrm.domain] real = cross(uniNrmA, uniNrm);
    var uniNrmB : [nrm.domain] real = cross(cross(uniNrmA, uniNrm), uniNrmA); // Author's implementation/suggestion
    uniNrmB = uniNrmB/norm(uniNrmB, normType.norm2);

    // Calculate the angle between the face normal and nrmB
    var alphaB = dot(uniNrm, uniNrmB);

    // Make alphaB always positive.
    if sgn(alphaB) < 0 {
        uniNrmB = sgn(alphaB) * uniNrmB;
        alphaB  = sgn(alphaB) * alphaB;
    }

    // Now we are going to compute the Roe flux with n2 as the normal with modified wave speeds (5.12).
    // NOTE: The Roe flux here is computed without tangent vectors. See "I do like CFD, VOL.1" for details: page 57,
    //       Equation (3.6.31).

    // Recalculate normal velocities using uniNrmB
    var velNrmBL = dot(velVL, uniNrmB);
    var velNrmBR = dot(velVR, uniNrmB);
    var velNrmB  = dot(velV , uniNrmB);

    // Differences in primitive variables
    var dDens    = densR    - densL;
    var dPres    = presR    - presL;
    var dVelNrmB = velNrmBR - velNrmBL;

    // Wave strengths (Characteristic Variables).
    var dV : [1..4] real;
    dV[1] =  (dPres - dens*a*dVelNrmB)/(2.0*a**2);
    dV[2] =   dDens - dPres/(a**2);
    dV[3] =  (dPres + dens*a*dVelNrmB)/(2.0*a**2);
    dV[4] =   dens;

    // Absolute values of the wave speeds (Eigenvalues)
    var eig : [1..4] real;
    eig[1] = velNrmB - a; // Left moving acoustic wave
    eig[2] = velNrmB;     // Entropy wave
    eig[3] = velNrmB + a; // Right moving acoustic wave
    eig[4] = velNrmB;     // Shear wave

    var ws : [1..4] real = abs(eig);

    // Harten's Entropy Fix JCP(1983), 49, pp357-393. Only for the nonlinear fields.
    // It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
    if ( ws[1] < 0.2 ) then ws[1] = (ws[1]**2)/0.4 + 0.1;
    if ( ws[3] < 0.2 ) then ws[3] = (ws[3]**2)/0.4 + 0.1;

    // Combine the wave speeds for Rotated-RHLL: Eq.(5.12) in the original JCP2008 paper.
    ws = ws*alphaB - (2.0*alphaA*SRp*SLm + alphaB*(SRp+SLm)*eig)/(SRp-SLm);

    // Right eigenvectors
    var R : [1..5, 1..4] real;

    // Negative acoustic wave
    R[idxDens, 1] = 1.0;
    R[idxMom , 1] = velV - a*uniNrmB;
    R[idxEner, 1] = enth - a*velNrmB;

    // Entropy wave
    R[idxDens, 2] = 1.0;
    R[idxMom , 2] = velV;
    R[idxEner, 2] = dot(velV, velV)/2.0;

    // Positive acoustic wave
    R[idxDens, 3] = 1.0;
    R[idxMom , 3] = velV + a*uniNrmB;
    R[idxEner, 3] = enth + a*velNrmB;

    // Combined shear waves
    R[idxDens, 4] = 0.0;
    R[idxMom , 4] = dVelV - dVelNrmB*uniNrmB;
    R[idxEner, 4] = dot(velV, dVelV) - velNrmB*dVelNrmB;

    // Compute the let and right fluxes
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL, fGamma));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR, fGamma));

    // Calculate the dissipation term
    var roeDiss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2]
                                      + ws[3]*dV[3]*R[.., 3] + ws[4]*dV[4]*R[.., 4];

    // Compute the final 3D Rotated-RHLL flux
    var rotRHLL : [consL.domain] real = (fluxL*SRp - fluxR*SLm)/(SRp-SLm) - roeDiss/2.0;

    return rotRHLL;
  }

  proc rusanov(const ref consL : [] real, const ref consR : [] real, const ref nrm : [] real) : [] real
  {
    // Generic 1D/2D/3D Rusanov solver
  }

  proc roe(const ref consL : [] real, const ref consR : [] real, const ref nrm : [] real) : [] real
  {
    use LinearAlgebra;
    import Input.fGamma;
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
    var presL   : real = pressure_cv(consL, fGamma);
    var aL      : real = sqrt(fGamma*presL/densL);
    var enthL   : real = ( enerL + presL ) / densL; // Mass specific Stagnation/Total Enthalpy

    //   Right state
    ref densR   : real = consR[idxDens];
    var velR    : [idxMom-1] real = consR[idxMom]/densR;
    ref enerR   : real = consR[idxEner];

    var velNrmR : real = dot[velR, uniNrm];
    var presR   : real = pressure_cv(consR, fGamma);
    var aR      : real = sqrt(fGamma*presR/densR);
    var enthR   : real = ( enerR + presR ) / densR; // Mass specific Stagnation/Total Enthalpy

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
    var fluxL : [consL.domain] real = dot(uniNrm, euler_flux_cv(consL, fGamma));
    var fluxR : [consR.domain] real = dot(uniNrm, euler_flux_cv(consR, fGamma));

    // Calculate the dissipation term
    var diss : [consL.domain] real = ws[1]*dV[1]*R[.., 1] + ws[2]*dV[2]*R[.., 2] + ws[3]*dV[3]*R[.., 3];
    if idxEner > 3 then diss += ws[4]*dV[4]*R[.., 4];

    // Compute the Roe flux
    var roe : [consL.domain] real = (fluxL + fluxR - diss)/2.0;

    return roe;
  }

  proc hll(const ref consL : [] real, const ref consR : [] real, const ref nrm : [] real) : [] real
  {
    // The HLL (Harten-Lax-van Leer)
  }

  proc hllc(const ref consL : [] real, const ref consR : [] real, const ref nrm : [] real) : [] real
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

  proc rotated_rhll(const ref consL : [] real, const ref consR : [] real, const ref nrm : [] real) : [] real
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

  proc godunov(const ref consL : [] real, const ref consR : [] real, const ref nrm : [] real) : [] real
  {
  }

  proc main()
  {
    use LinearAlgebra;
    import Flux.euler_flux_cv;

    param fGamma : real = 1.4;

    var nrm1 : [1..3] real = [+0.612372436, +0.353553391, +0.707106781];
    var nrm2 : [1..3] real = [-0.612372436, -0.353553391, -0.707106781];
    var nrm3 : [1..3] real = [+1.13651    , +0.699182   , +0.1];

    writeln();
    writeln("1D Tests:");
    {
      var consL1 : [1..3] real = [1.225, 250.1160830494510, 278846.40];
      var consR1 : [1..3] real = [1.125, 265.2881676121270, 259260.30];
      var consL2 : [1..3] real = [1.225, 250.1160830494510, 278846.40];
      var consR2 : [1..3] real = [1.125, 265.2881676121270, 259260.30];

      var uniNrm1 : [1..1] real = nrm1[1..1]/norm(nrm1[1..1]);
      var uniNrm2 : [1..1] real = nrm2[1..1]/norm(nrm2[1..1]);

      writef("\n");
      writef("  Normal 1:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 1 %t: %+.8ht\n",    nrm1[1..1].domain,    nrm1[1..1]);
      writef("      UniNrm 1 %t: %+.8ht\n", uniNrm1[1..1].domain, uniNrm1[1..1]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  %t: %+.8ht\n", consL1.domain, consL1);
      writef("      Right %t: %+.8ht\n", consR1.domain, consR1);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm1[1..1], euler_flux_cv(consL1, fGamma)).domain,
                                                  dot( uniNrm1[1..1], euler_flux_cv(consL1, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm1[1..1], euler_flux_cv(consR1, fGamma)).domain,
                                                  dot( uniNrm1[1..1], euler_flux_cv(consR1, fGamma))       );
      writef("    1D  Roe     Flux %t: %+.8ht\n", roe_1d(consR1, consL1, nrm1[1..1]).domain,
                                                  roe_1d(consR1, consL1, nrm1[1..1])       );
      writef("    Gen Roe 1D  Flux %t: %+.8ht\n", roe(   consR1, consL1, nrm1[1..1]).domain,
                                                  roe(   consR1, consL1, nrm1[1..1])       );

      writef("\n");
      writef("  Normal 2:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 2 %t: %+.8ht\n",    nrm2[1..1].domain,    nrm2[1..1]);
      writef("      UniNrm 2 %t: %+.8ht\n", uniNrm2[1..1].domain, uniNrm2[1..1]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  %t: %+.8ht\n", consL2.domain, consL2);
      writef("      Right %t: %+.8ht\n", consR2.domain, consR2);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm2[1..1], euler_flux_cv(consL2, fGamma)).domain,
                                                  dot( uniNrm2[1..1], euler_flux_cv(consL2, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm2[1..1], euler_flux_cv(consR2, fGamma)).domain,
                                                  dot( uniNrm2[1..1], euler_flux_cv(consR2, fGamma))       );
      writef("    1D  Roe     Flux %t: %+.8ht\n", roe_1d(consR2, consL2, nrm2[1..1]).domain,
                                                  roe_1d(consR2, consL2, nrm2[1..1])       );
      writef("    Gen Roe 1D  Flux %t: %+.8ht\n", roe(   consR2, consL2, nrm2[1..1]).domain,
                                                  roe(   consR2, consL2, nrm2[1..1])       );
    }

    writeln();
    writeln("------------------------------------------------------------------------------------------------------------------------");
    writeln();
    writeln("2D Tests:");
    {
      var consL1 : [1..4] real = [1.225, 249.1643158413380, 21.79905299130990, 278846.40];
      var consR1 : [1..4] real = [1.125, 264.2786660416750, 23.12138729040020, 259260.30];
      var consL2 : [1..4] real = [1.225, 249.1643158413380, 21.79905299130990, 278846.40];
      var consR2 : [1..4] real = [1.125, 264.2786660416750, 23.12138729040020, 259260.30];
      var consL3 : [1..4] real = [0.714259, 0.356919, 0.0124639, 1.3648];
      var consR3 : [1..4] real = [0.714259, 0.356919, 0.0124639, 1.3648];

      var uniNrm1 : [1..2] real = nrm1[1..2]/norm(nrm1[1..2]);
      var uniNrm2 : [1..2] real = nrm2[1..2]/norm(nrm2[1..2]);
      var uniNrm3 : [1..2] real = nrm3[1..2]/norm(nrm3[1..2]);

      writef("\n");
      writef("  Normal 1:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 1 %t: %+.8ht\n",    nrm1[1..2].domain,    nrm1[1..2]);
      writef("      UniNrm 1 %t: %+.8ht\n", uniNrm1[1..2].domain, uniNrm1[1..2]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  1 %t: %+.8ht\n", consL1.domain, consL1);
      writef("      Right 1 %t: %+.8ht\n", consR1.domain, consR1);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm1[1..2], euler_flux_cv(consL1, fGamma)).domain,
                                                  dot( uniNrm1[1..2], euler_flux_cv(consL1, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm1[1..2], euler_flux_cv(consR1, fGamma)).domain,
                                                  dot( uniNrm1[1..2], euler_flux_cv(consR1, fGamma))       );
      writef("    2D  Roe     Flux %t: %+.8ht\n", roe_2d(         consL1, consR1, nrm1[1..2]).domain,
                                                  roe_2d(         consL1, consR1, nrm1[1..2])       );
      writef("    Gen Roe 2D  Flux %t: %+.8ht\n", roe(            consL1, consR1, nrm1[1..2]).domain,
                                                  roe(            consL1, consR1, nrm1[1..2])       );
      writef("    2D  RotRHLL Flux %t: %+.8ht\n", rotated_rhll_2d(consL1, consR1, nrm1[1..2]).domain,
                                                  rotated_rhll_2d(consL1, consR1, nrm1[1..2])       );

      writef("\n");
      writef("  Normal 2:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 2 %t: %+.8ht\n",    nrm2[1..2].domain,    nrm2[1..2]);
      writef("      UniNrm 2 %t: %+.8ht\n", uniNrm2[1..2].domain, uniNrm2[1..2]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  2 %t: %+.8ht\n", consL2.domain, consL2);
      writef("      Right 2 %t: %+.8ht\n", consR2.domain, consR2);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm2[1..2], euler_flux_cv(consL2, fGamma)).domain,
                                                  dot( uniNrm2[1..2], euler_flux_cv(consL2, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm2[1..2], euler_flux_cv(consR2, fGamma)).domain,
                                                  dot( uniNrm2[1..2], euler_flux_cv(consR2, fGamma))       );
      writef("    2D  Roe     Flux %t: %+.8ht\n", roe_2d(         consL2, consR2, nrm2[1..2]).domain,
                                                  roe_2d(         consL2, consR2, nrm2[1..2])       );
      writef("    Gen Roe 2D  Flux %t: %+.8ht\n", roe(            consL2, consR2, nrm2[1..2]).domain,
                                                  roe(            consL2, consR2, nrm2[1..2])       );
      writef("    2D  RotRHLL Flux %t: %+.8ht\n", rotated_rhll_2d(consL2, consR2, nrm2[1..2]).domain,
                                                  rotated_rhll_2d(consL2, consR2, nrm2[1..2])       );

      writef("\n");
      writef("  Test 3:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 3 %t: %+.8ht\n",    nrm3[1..2].domain,    nrm3[1..2]);
      writef("      UniNrm 3 %t: %+.8ht\n", uniNrm3[1..2].domain, uniNrm3[1..2]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  3 %t: %+.8ht\n", consL3.domain, consL3);
      writef("      Right 3 %t: %+.8ht\n", consR3.domain, consR3);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm3[1..2], euler_flux_cv(consL3, fGamma)).domain,
                                                  dot( uniNrm3[1..2], euler_flux_cv(consL3, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm3[1..2], euler_flux_cv(consR3, fGamma)).domain,
                                                  dot( uniNrm3[1..2], euler_flux_cv(consR3, fGamma))       );
      writef("    2D  Roe     Flux %t: %+.8ht\n", roe_2d(         consL3, consR3, nrm3[1..2]).domain,
                                                  roe_2d(         consL3, consR3, nrm3[1..2])       );
      writef("    Gen Roe 2D  Flux %t: %+.8ht\n", roe(            consL3, consR3, nrm3[1..2]).domain,
                                                  roe(            consL3, consR3, nrm3[1..2])       );
      writef("    2D  RotRHLL Flux %t: %+.8ht\n", rotated_rhll_2d(consL3, consR3, nrm3[1..2]).domain,
                                                  rotated_rhll_2d(consL3, consR3, nrm3[1..2]));
    }

    writeln();
    writeln("------------------------------------------------------------------------------------------------------------------------");
    writeln();
    writeln("3D Tests:");
    {
      var consL1 : [1..5] real = [1.225, 249.0125316723220, 17.79514529902850, 7.098538138029130, 278846.40];
      var consR1 : [1..5] real = [1.125, 264.1176746188930, 23.12138729040020, 9.223192434062800, 259260.30];
      var consL2 : [1..5] real = [1.225, 249.0125316723220, 17.79514529902850, 7.098538138029130, 278846.40];
      var consR2 : [1..5] real = [1.125, 264.1176746188930, 23.12138729040020, 9.223192434062800, 259260.30];
      var consL3 : [1..5] real = [0.714259, 0.356919, 0.0124639, 0.005, 1.3648];
      var consR3 : [1..5] real = [0.714259, 0.356919, 0.0124639, 0.005, 1.3648];

      var uniNrm1 : [1..3] real = nrm1[1..3]/norm(nrm1[1..3]);
      var uniNrm2 : [1..3] real = nrm2[1..3]/norm(nrm2[1..3]);
      var uniNrm3 : [1..3] real = nrm3[1..3]/norm(nrm3[1..3]);

      writef("\n");
      writef("  Test 1:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 1 %t: %+.8ht\n",    nrm1[1..3].domain,    nrm1[1..3]);
      writef("      UniNrm 1 %t: %+.8ht\n", uniNrm1[1..3].domain, uniNrm1[1..3]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  1 %t: %+.8ht\n", consL1.domain, consL1);
      writef("      Right 1 %t: %+.8ht\n", consR1.domain, consR1);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm1[1..3], euler_flux_cv(consL1, fGamma)).domain,
                                                  dot( uniNrm1[1..3], euler_flux_cv(consL1, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm1[1..3], euler_flux_cv(consR1, fGamma)).domain,
                                                  dot( uniNrm1[1..3], euler_flux_cv(consR1, fGamma))       );
      writef("    3D  Roe     Flux %t: %+.8ht\n", roe_3d(         consL1, consR1, nrm1[1..3]).domain,
                                                  roe_3d(         consL1, consR1, nrm1[1..3])       );
      writef("    Gen Roe 3D  Flux %t: %+.8ht\n", roe(            consL1, consR1, nrm1[1..3]).domain,
                                                  roe(            consL1, consR1, nrm1[1..3])       );
      writef("    3D  RotRHLL Flux %t: %+.8ht\n", rotated_rhll_3d(consL1, consR1, nrm1[1..3]).domain,
                                                  rotated_rhll_3d(consL1, consR1, nrm1[1..3])       );

      writef("\n");
      writef("  Test 2:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 2 %t: %+.8ht\n",    nrm2[1..3].domain,    nrm2[1..3]);
      writef("      UniNrm 2 %t: %+.8ht\n", uniNrm2[1..3].domain, uniNrm2[1..3]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  2 %t: %+.8ht\n", consL2.domain, consL1);
      writef("      Right 2 %t: %+.8ht\n", consR2.domain, consR1);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm2[1..3], euler_flux_cv(consL1, fGamma)).domain,
                                                  dot( uniNrm2[1..3], euler_flux_cv(consL1, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm2[1..3], euler_flux_cv(consR1, fGamma)).domain,
                                                  dot( uniNrm2[1..3], euler_flux_cv(consR1, fGamma))       );
      writef("    3D  Roe     Flux %t: %+.8ht\n", roe_3d(         consL2, consR2, nrm2[1..3]).domain,
                                                  roe_3d(         consL2, consR2, nrm2[1..3])       );
      writef("    Gen Roe 3D  Flux %t: %+.8ht\n", roe(            consL2, consR2, nrm2[1..3]).domain,
                                                  roe(            consL2, consR2, nrm2[1..3])       );
      writef("    3D  RotRHLL Flux %t: %+.8ht\n", rotated_rhll_3d(consL2, consR2, nrm2[1..3]).domain,
                                                  rotated_rhll_3d(consL2, consR2, nrm2[1..3])       );

      writef("\n");
      writef("  Test 3:\n");
      writef("    Face normal unit vectors:\n");
      writef("      Normal 3 %t: %+.8ht\n",    nrm2[1..3].domain,    nrm3[1..3]);
      writef("      UniNrm 3 %t: %+.8ht\n", uniNrm2[1..3].domain, uniNrm3[1..3]);
      writef("\n");
      writef("    Conserved variables:\n");
      writef("      Left  3 %t: %+.8ht\n", consL3.domain, consL1);
      writef("      Right 3 %t: %+.8ht\n", consR3.domain, consR1);
      writef("\n");
      writef("    Left  Euler Flux %t: %+.8ht\n", dot( uniNrm3[1..3], euler_flux_cv(consL1, fGamma)).domain,
                                                  dot( uniNrm3[1..3], euler_flux_cv(consL1, fGamma))       );
      writef("    Right Euler Flux %t: %+.8ht\n", dot( uniNrm3[1..3], euler_flux_cv(consR1, fGamma)).domain,
                                                  dot( uniNrm3[1..3], euler_flux_cv(consR1, fGamma))       );
      writef("    3D  Roe     Flux %t: %+.8ht\n", roe_3d(         consL3, consR3, nrm3[1..3]).domain,
                                                  roe_3d(         consL3, consR3, nrm3[1..3])       );
      writef("    Gen Roe 3D  Flux %t: %+.8ht\n", roe(            consL3, consR3, nrm3[1..3]).domain,
                                                  roe(            consL3, consR3, nrm3[1..3])       );
      writef("    3D  RotRHLL Flux %t: %+.8ht\n", rotated_rhll_3d(consL3, consR3, nrm3[1..3]).domain,
                                                  rotated_rhll_3d(consL3, consR3, nrm3[1..3])       );

    }
  }
}
