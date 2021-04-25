prototype module Flux
{
  proc pressure(cons : [] real ) : real
  {
    use Input;
    import LinearAlgebra.dot;

    var idxRho : int   = cons.domain.dim(0).low;           // First element is density
    var idxMom : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = cons.domain.dim(0).high;          // Last element is energy

    var pressure : real = (fGamma-1.0)*(cons[idxEne] - 0.5*dot(cons[idxMom],cons[idxMom])/cons[1]);

    return pressure;
  }

  proc invs_flux_cv_1d(u : [1..3] real ) : [1..3] real
  {
    import LinearAlgebra.dot;

    var idxRho : int   = u.domain.dim(0).low;           // First element is density
    var idxMom : range = u.domain.dim(0).expand(-1);    // Intermediary elements are the momentum components
    var idxEne : int   = u.domain.dim(0).high;          // Last element is energy

    var invs_flux_cv : [u.domain] real;
    var p : real = pressure(u);

    invs_flux_cv[1] = u[2];
    invs_flux_cv[2] = u[2]*u[2]/u[1] + p;
    invs_flux_cv[3] = u[2]/u[1]*(u[3] + p);

    return invs_flux_cv;
  }

  proc invs_flux_cv(cons : [] real ) : [] real
  {
    import LinearAlgebra.dot;

    var idxRho : int   = cons.domain.dim(0).low;           // First element is density
    var idxMom : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = cons.domain.dim(0).high;          // Last element is energy

    var invs_flux_cv : [idxMom-1, cons.domain.dim(0)] real;

    var vel : [idxMom] real = cons[idxMom] / cons[idxRho];
    var p   : real = pressure(cons);

    for i in idxMom-1
    {
      invs_flux_cv[i, idxRho] = cons[i+1];
      invs_flux_cv[i, idxMom] = cons[i+1]*vel[idxMom];
      invs_flux_cv[i, idxEne] = vel[i+1]*(cons[idxEne] + p);

      invs_flux_cv[i, i+1] += p;
    }

    return invs_flux_cv;
  }

  proc invs_flux_pv(prim : [] real ) : [] real
  {
    use Input;
    import LinearAlgebra.dot;

    var idxRho : int   = prim.domain.dim(0).low;           // First element is density
    var idxVel : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocity components
    var idxPre : int   = prim.domain.dim(0).high;          // Last element is pressure

    var invs_flux_pv : [idxVel-1, prim.domain.dim(0)] real;

    for i in idxVel-1
    {
      invs_flux_pv[i, idxRho] = prim[idxRho+i]*prim[idxRho];
      invs_flux_pv[i, idxVel] = prim[idxRho+i]*prim[idxRho]*prim[idxVel];
      invs_flux_pv[i, idxPre] = prim[idxRho+i]*(prim[idxPre]*fGamma/(fGamma-1.0)+0.5*prim[idxRho]*dot(prim[idxVel],prim[idxVel]));

      invs_flux_pv[i, i+1] += prim[idxPre];
    }

    return invs_flux_pv;
  }

  proc main()
  {
    var cons1d : [1..3] real = [1.225, 250.1160830494510, 278846.40];
    var cons2d : [1..4] real = [1.225, 249.9637190895680, 8.728925415626780, 278846.40];
    var cons3d : [1..5] real = [1.225, 249.9637190895680, 8.728925415626780, 0.0, 278846.40];

    var prim1d : [1..3] real = [1.225, 204.1763943260830, 101325.0];
    var prim2d : [1..4] real = [1.225, 204.0520155833210, 7.125653400511660, 101325.0];
    var prim3d : [1..5] real = [1.225, 204.0520155833210, 7.125653400511660, 0.0, 101325.0];

    writeln("Conserverd variables:");
    writeln("1D: ", cons1d);
    writeln("2D: ", cons2d);
    writeln("3D: ", cons3d);
    writeln();
    writeln("Primitive Variables:");
    writeln("1D: ", prim1d);
    writeln("2D: ", prim2d);
    writeln("3D: ", prim3d);
    writeln();
    writeln("Inviscid 1D Flux:\n", invs_flux_cv_1d(cons1d));
    writeln();
    writeln("Generic inviscid flux function (conserver vars):");
    writeln("  Inviscid 1D Flux:\n", invs_flux_cv(cons1d));
    writeln("  Inviscid 2D Flux:\n", invs_flux_cv(cons2d));
    writeln("  Inviscid 3D Flux:\n", invs_flux_cv(cons3d));
    writeln();
    writeln("Generic inviscid flux function (primitive vars):");
    writeln("  Inviscid 1D Flux:\n", invs_flux_pv(prim1d));
    writeln("  Inviscid 2D Flux:\n", invs_flux_pv(prim2d));
    writeln("  Inviscid 3D Flux:\n", invs_flux_pv(prim3d));
  }
}
