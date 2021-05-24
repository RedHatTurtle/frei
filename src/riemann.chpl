prototype module Riemann
{
  proc rusanov_1d(uL : [1..3] real, uR : [1..3] real) : [1..3] real
  {
    use Input;
    import Flux.pressure_cv;
    import Flux.invs_flux_cv;

    var idxRho : int   = uL.domain.dim(0).low;           // First element is density
    var idxMom : range = uL.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = uL.domain.dim(0).high;          // Last element is energy

    var gm1 : real = fGamma - 1.0;
    var gp1 : real = fGamma + 1.0;

    // Primitive variables
    //  Left state
    var rhoL : real = uL[1];
    var vL   : real = uL[2]/uL[1];
    var pL   : real = pressure_cv(uL);
    var hL   : real = ( uL[3] + pL ) / rhoL;
    //  Right state
    var rhoR : real = uR[1];
    var vR   : real = uR[2]/uR[1];
    var pR   : real = pressure_cv(uR);
    var hR   : real = ( uR[3] + pR ) / rhoR;

    // Compute the Roe Averages
    var rt  : real = sqrt(rhoR/rhoL);
    var rho : real = rt*rhoL;
    var v   : real = (vL+rt*vR)/(1.0+rt);
    var h   : real = (hL+rt*hR)/(1.0+rt);
    var a   : real = sqrt( (fGamma-1.0)*(h-0.5*v*v) );

    var smax : real = abs(v) + a;

    // Compute the average flux.
    var rusanov : [1..3] real = 0.5*( invs_flux_cv(uL)[1,1..3] + invs_flux_cv(uR)[1,1..3] - smax*(uR-uL));

    return rusanov;
  }

  proc roe_1d(uL : [1..3] real, uR : [1..3] real) : [1..3] real
  {
    use Input;
    import Flux.pressure_cv;
    import Flux.invs_flux_cv;

    var gm1 : real = fGamma - 1.0;
    var gp1 : real = fGamma + 1.0;

    //Primitive and other variables.
    //  Left state
    var rhoL : real = uL[1];
    var vL   : real = uL[2]/uL[1];
    var pL   : real = pressure_cv(uL);
    var aL   : real = sqrt(fGamma*pL/rhoL);
    var HL   : real = ( uL[3] + pL ) / rhoL;
    //  Right state
    var rhoR : real = uR[1];
    var vR   : real = uR[2]/uR[1];
    var pR   : real = pressure_cv(uR);
    var aR   : real = sqrt(fGamma*pR/rhoR);
    var HR   : real = ( uR[3] + pR ) / rhoR;

    //First compute the Roe Averages
    var RT  : real = sqrt(rhoR/rhoL);
    var rho : real = RT*rhoL;
    var v   : real = (vL+RT*vR)/(1.0+RT);
    var H   : real = (HL+RT*HR)/(1.0+RT);
    var a   : real = sqrt( (fGamma-1.0)*(H-0.5*v*v) );

    //Differences in primitive variables.
    var drho = rhoR - rhoL;
    var du   =   vR - vL;
    var dP   =   pR - pL;

    //Wave strength (Characteristic Variables).
    var dV : [1..3] real;
    dV[1] =  0.5*(dP-rho*a*du)/(a*a);
    dV[2] = -( dP/(a*a) - drho );
    dV[3] =  0.5*(dP+rho*a*du)/(a*a);

    //Absolute values of the wave speeds (Eigenvalues)
    var ws : [1..3] real;
    ws[1] = abs(v-a);
    ws[2] = abs(v  );
    ws[3] = abs(v+a);

    //Modified wave speeds for nonlinear fields (to remove expansion shocks).
    //There are various ways to implement an entropy fix. This is just one
    //example.
    var Da : real = max(0.0, 4.0*((vR-aR)-(vL-aL)) );
    if (ws(1) < 0.5*Da) then ws(1) = ws(1)*ws(1)/Da + 0.25*Da;
    Da = max(0.0, 4.0*((vR+aR)-(vL+aL)) );
    if (ws(3) < 0.5*Da) then ws(3) = ws(3)*ws(3)/Da + 0.25*Da;

    //Right eigenvectors
    var R : [1..3,1..3] real;
    R[1,1] = 1.0;
    R[2,1] = v - a;
    R[3,1] = H - v*a;

    R[1,2] = 1.0;
    R[2,2] = v;
    R[3,2] = 0.5*v*v;

    R[1,3] = 1.0;
    R[2,3] = v + a;
    R[3,3] = H + v*a;

    //Compute the average flux.
    var roe : [1..3] real = 0.5*( invs_flux_cv(uL)[1,1..3] + invs_flux_cv(uR)[1,1..3] );

    //!Add the matrix dissipation term to complete the Roe flux.
    for j in 1..3 do {
      for k in 1..3 do {
        roe[j] = roe[j] - 0.5*ws[k]*dV[k]*R[j,k];
      }
    }

    return roe;
  }

  proc main()
  {
    import Flux.invs_flux_cv;

    var cons1dL : [1..3] real = [1.225, 250.1160830494510, 278846.40];
    var cons2dL : [1..4] real = [1.225, 249.9637190895680, 8.728925415626780, 278846.40];
    var cons3dL : [1..5] real = [1.225, 249.9637190895680, 8.728925415626780, 0.0, 278846.40];

    var cons1dR : [1..3] real = [1.5, 300.1160830494510, 298846.40];
    var cons2dR : [1..4] real = [1.5, 300.9637190895680, 9.728925415626780, 298846.40];
    var cons3dR : [1..5] real = [1.5, 300.9637190895680, 9.728925415626780, 0.0, 298846.40];

    writeln("Conserverd variables:");
    writeln("  Left:");
    writeln("1D: ", cons1dL);
    writeln("2D: ", cons2dL);
    writeln("3D: ", cons3dL);
    writeln("  Right:");
    writeln("1D: ", cons1dR);
    writeln("2D: ", cons2dR);
    writeln("3D: ", cons3dR);
    writeln();
    writeln("Left  Flux:   ", invs_flux_cv(cons1dL));
    writeln("Right Flux:   ", invs_flux_cv(cons1dR));
    writeln("Runanov Flux: ", rusanov_1d(cons1dL, cons1dR));
  }
}
