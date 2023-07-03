module Limiter
{
  use UnitTest;
  use Set;

  proc scaled_marker_projection(solPoly : [] real, jacobian : [] real, cellTopo : int, solDegree : int, projDegree : int) : real
  {
    use Parameters.ParamMesh;
    use Projection;
    use Quadrature;
    import Math.half_pi;

    /*
      Resolution Indicator
    */

    // Project the density poly to a lower degree $P_m$
    var projPoly : [solPoly.domain] real = project_poly(nodalPoly = solPoly,
                                                        cellTopo = cellTopo,
                                                        polyDegree = solDegree,
                                                        projDegree = projDegree);

    var sol2Poly : [solPoly.domain] real = solPoly**2.0;
    var dif2Poly : [solPoly.domain] real = (solPoly - projPoly)**2.0;

    // Integrate $(P_n)^2$ and $(P_n-P_m)^2$ in the cell
    //l2_p = sum( rho_p(1:np) * rho_p(1:np) * metrics_dt(nc)%jac(1:np) * &
    //            std_elem(this_geom,this_order)%wts(1:np) )
    var l2_p    : real = integrate(nodalPoly = sol2Poly,
                                   jacobian = jacobian,
                                   elemTopo = cellTopo,
                                   nodeCnt = solDegree);

    //l2_pm1 = sum( rho_pm1(1:np) * rho_pm1(1:np) * metrics_dt(nc)%jac(1:np) * &
    //              std_elem(this_geom,this_order)%wts(1:np) )
    var l2_pm1  : real = integrate(nodalPoly = dif2Poly,
                                   jacobian = jacobian,
                                   elemTopo = cellTopo,
                                   nodeCnt = solDegree);

    // Calculate the Log10 of the ratio between the quantities above
    //marker(nc) = log10( l2_pm1 / l2_p )
    var marker : real = log10(l2_pm1 / l2_p);

    // Apply a sigmoid transformation to this value and calculate a "scaled" marker
    //psi0 = -one*(c1 + c2*log10(real(cutoff_order,kind=wp)))
    param dpsi : real = 0.5;
    param c1 : real = 4.0;
    param c2 : real = 4.25;
    var psi0 : real = -1.0*(c1 + c2*log10((solDegree-1):real));
    var scaledMarker : real;

    // If the scaled marker is over a threshold then mark the cell as problematic
    if (marker <= psi0-dpsi) then
      scaledMarker = 0.0;
    else if (marker >= psi0+dpsi) then
      scaledMarker = 1.0;
    else
      scaledMarker = 0.5*(1.0 + sin(half_pi*(marker-psi0)/dpsi));

    return scaledMarker;
  }

  proc troubled_cell_marker(solPoly : [] real, jacobian : [] real, cellTopo : int, solDegree : int) : int
  {
    use Parameters.ParamMesh;
    use Projection;
    use Quadrature;
    import Math.half_pi;
    import Math.log10;
    import Math.sin;

    /*
      Resolution Indicator
    */

    var refPoly : [solPoly.domain] real = solPoly;
    var stableDegree : int;

    for projDegree  in 0..solDegree-1 by -1
    {
      // Project the density poly to a lower degree $P_m$
      var projPoly : [solPoly.domain] real = project_poly(nodalPoly = solPoly,
                                                          cellTopo = cellTopo,
                                                          polyDegree = solDegree,
                                                          projDegree = projDegree);

      var ref2Poly : [solPoly.domain] real = refPoly**2.0;
      var dif2Poly : [solPoly.domain] real = (refPoly - projPoly)**2.0;

      // Integrate $(P_n)^2$ and $(P_n-P_m)^2$ in the cell
      //l2_p = sum( rho_p(1:np) * rho_p(1:np) * metrics_dt(nc)%jac(1:np) * &
      //            std_elem(this_geom,this_order)%wts(1:np) )
      var l2_ref  : real = integrate(nodalPoly = ref2Poly,
                                     jacobian = jacobian,
                                     elemTopo = cellTopo,
                                     nodeCnt = solDegree);

      //l2_pm1 = sum( rho_pm1(1:np) * rho_pm1(1:np) * metrics_dt(nc)%jac(1:np) * &
      //              std_elem(this_geom,this_order)%wts(1:np) )
      var l2_diff : real = integrate(nodalPoly = dif2Poly,
                                     jacobian = jacobian,
                                     elemTopo = cellTopo,
                                     nodeCnt = solDegree);

      // Calculate the Log10 of the ratio between the quantities above
      //marker(nc) = log10( l2_pm1 / l2_p )
      var marker : real = log10(l2_diff / l2_ref);

      // Apply a sigmoid transformation to this value and calculate a "scaled" marker
      //psi0 = -one*(c1 + c2*log10(real(cutoff_order,kind=wp)))
      param dpsi : real = 0.5;
      param c1 : real = 4.0;
      param c2 : real = 4.25;
      var psi0 : real = -1.0*(c1 + c2*log10((projDegree):real));
      var scaledMarker : real;

      // If the scaled marker is over a threshold then mark the cell as problematic
      if (marker <= psi0-dpsi) then
        scaledMarker = 0.0;
      else if (marker >= psi0+dpsi) then
        scaledMarker = 1.0;
      else
        scaledMarker = 0.5*(1.0 + sin(half_pi*(marker-psi0)/dpsi));

      // Check if the reference solution is stable and either stop and output the stable degree or keep projecting down
      if scaledMarker <= 1.0e-8
      {
        stableDegree = projDegree+1;
        break;
      }
      else
      {
        refPoly = projPoly;
      }
    }

    return stableDegree;
  }

  proc projection_limiter(solPoly : [] real, cellTopo : int, solDegree : int, projDegree : int)
  {
    use Projection;
    /*
       Project problematic solution interpolation to to a lower order
    */

    return project_poly(nodalPoly = solPoly, cellTopo = cellTopo, polyDegree = solDegree, projDegree = projDegree);
  }

  proc main()
  {
    // Generate a smooth and a discontinuous test functions to validate the marker and limiter

    use Random;
    use Parameters.ParamTest;
    use Parameters.ParamMesh;
    use Polynomials;
    use Projection;
    use Quadrature;
    import Math.half_pi;
    import Math.exp;
    import Math.log;
    import Math.sin;

    // Create a set with the cell topologies contained in the hypothetical test mesh
    var cellTopos : set(int);
    cellTopos.add(TOPO_LINE);
    cellTopos.add(TOPO_QUAD);

    var minDegree : int = 0;
    var maxDegree : int = 5;

    init_polyProj(minDegree, maxDegree, cellTopos);
    init_quadratureWeights(minDegree, maxDegree, cellTopos);

    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    // Line Mode
    writeln();
    writeln("------------------------------------------------------------");
    writeln();
    writeln("Line Cell");
    for polyDegree in minDegree..maxDegree
    {
      writeln();
      writef("%i-Degree Interpolation:\n", polyDegree);

      var nodeDistLine : [1..polyDegree+1] real = nodes_legendre_gauss(polyDegree+1);

      var nodeCnt   : int = polyDegree+1;
      var jacobian  : [1..nodeCnt] real = 0.5;
      var nodalSin  : [1..nodeCnt] real;
      var nodalExp  : [1..nodeCnt] real;
      var nodalLog  : [1..nodeCnt] real;
      var nodalStep : [1..nodeCnt] real;

      for nodeIdx in 1..nodeCnt
      {
        var xi  : real = nodeDistLine[nodeIdx];

        nodalSin[nodeIdx] = sin(xi*half_pi);
        nodalExp[nodeIdx] = exp(xi);
        nodalLog[nodeIdx] = log(1.0+xi);

        if xi > 0.0 then
          nodalStep[nodeIdx] = 1.0;
        else
          nodalStep[nodeIdx] = 0.0;
      }

      writeln();
      writef("  Sin Stable degree: %2i\n",    troubled_cell_marker( nodalSin, jacobian , TOPO_LINE, polyDegree));
      writef("  Original Solution: %{ ht}\n", nodalSin);
      writef("  Limited Solution:  %{ ht}\n", projection_limiter(   nodalSin, TOPO_LINE, polyDegree ,
                                              troubled_cell_marker( nodalSin, jacobian , TOPO_LINE, polyDegree)));
      writeln();
      writef("  Exp Stable degree: %2i\n",    troubled_cell_marker( nodalExp, jacobian , TOPO_LINE, polyDegree));
      writef("  Original Solution: %{ ht}\n", nodalExp);
      writef("  Limited Solution:  %{ ht}\n", projection_limiter(   nodalExp, TOPO_LINE, polyDegree ,
                                              troubled_cell_marker( nodalExp, jacobian , TOPO_LINE, polyDegree)));
      writeln();
      writef("  Log Stable degree: %2i\n",    troubled_cell_marker( nodalLog, jacobian , TOPO_LINE, polyDegree));
      writef("  Original Solution: %{+ht}\n", nodalLog);
      writef("  Limited Solution:  %{ ht}\n", projection_limiter(   nodalLog, TOPO_LINE, polyDegree ,
                                              troubled_cell_marker( nodalLog, jacobian , TOPO_LINE, polyDegree)));
      writeln();
      writef("  Step Stable degree: %2i\n",   troubled_cell_marker(nodalStep, jacobian, TOPO_LINE, polyDegree));
      writef("  Original Solution: %{ ht}\n", nodalStep);
      writef("  Limited Solution:  %{ ht}\n", projection_limiter(  nodalStep, TOPO_LINE, polyDegree ,
                                              troubled_cell_marker(nodalStep, jacobian, TOPO_LINE, polyDegree)));
    }

    // Quad Mode
    writeln();
    writeln("------------------------------------------------------------");
    writeln();
    writeln("Quad Cell");
    for polyDegree in minDegree..maxDegree
    {
      writeln();
      writef("%i-Degree Interpolation:\n", polyDegree);

      var nodeDistLine : [1..polyDegree+1] real = nodes_legendre_gauss(polyDegree+1);

      var nodeCnt   : int = (polyDegree+1)**2;
      var jacobian  : [1..nodeCnt] real = 0.25;
      var nodalSin  : [1..nodeCnt] real;
      var nodalExp  : [1..nodeCnt] real;
      var nodalLog  : [1..nodeCnt] real;
      var nodalStep : [1..nodeCnt] real;

      for nodeIdx in 1..nodeCnt
      {
        var xi  : real = nodeDistLine[(nodeIdx-1)/(polyDegree+1)+1];
        var eta : real = nodeDistLine[(nodeIdx-1)%(polyDegree+1)+1];

        nodalSin[nodeIdx] = sin(xi*half_pi) * sin(eta*half_pi);
        nodalExp[nodeIdx] = exp(xi) * exp(eta);
        nodalLog[nodeIdx] = log(1.0+xi) * log(1.0+eta);

        if xi+2*eta > 0.0 then
          nodalStep[nodeIdx] = 1.0;
        else
          nodalStep[nodeIdx] = 0.0;
      }

      writeln();
      writef("  Sin Stable degree: %2i\n",    troubled_cell_marker( nodalSin , jacobian, TOPO_QUAD, polyDegree));
      writef("  Original Solution: %{ ht}\n", nodalSin);
      writef("  Limited Solution:  %{ ht}\n", projection_limiter(   nodalSin , TOPO_QUAD, polyDegree ,
                                              troubled_cell_marker( nodalSin , jacobian, TOPO_QUAD, polyDegree)));
      writeln();
      writef("  Exp Stable degree: %2i\n",    troubled_cell_marker( nodalExp , jacobian, TOPO_QUAD, polyDegree));
      writef("  Original Solution: %{ ht}\n", nodalExp);
      writef("  Limited Solution:  %{ ht}\n", projection_limiter(   nodalExp , TOPO_QUAD, polyDegree ,
                                              troubled_cell_marker( nodalExp , jacobian, TOPO_QUAD, polyDegree)));
      writeln();
      writef("  Log Stable degree: %2i\n",    troubled_cell_marker( nodalLog , jacobian, TOPO_QUAD, polyDegree));
      writef("  Original Solution: %{ ht}\n", nodalLog);
      writef("  Limited Solution:  %{ ht}\n", projection_limiter(   nodalLog , TOPO_QUAD, polyDegree ,
                                              troubled_cell_marker( nodalLog , jacobian, TOPO_QUAD, polyDegree)));
      writeln();
      writef("  Step Stable degree: %2i\n",    troubled_cell_marker(nodalStep, jacobian, TOPO_QUAD, polyDegree));
      writef("  Original Solution:  %{ ht}\n", nodalStep);
      writef("  Limited Solution:   %{ ht}\n", projection_limiter(  nodalStep, TOPO_QUAD, polyDegree ,
                                               troubled_cell_marker(nodalStep, jacobian, TOPO_QUAD, polyDegree)));
    }
  }
}
