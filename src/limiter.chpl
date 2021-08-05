prototype module Limiter
{
  use UnitTest;
  use Set;

  proc scaled_marker_projection(solPoly : [] real, cellTopo : int, solDegree : int, projDegree : int) : real
  {
    use Parameters.ParamMesh;
    use Projection;
    use Quadrature;

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
                                   elemTopo = cellTopo,
                                   nodeCnt = solDegree);

    //l2_pm1 = sum( rho_pm1(1:np) * rho_pm1(1:np) * metrics_dt(nc)%jac(1:np) * &
    //              std_elem(this_geom,this_order)%wts(1:np) )
    var l2_pm1  : real = integrate(nodalPoly = dif2Poly,
                                   elemTopo = cellTopo,
                                   nodeCnt = solDegree);

    // Calculate the Log10 of the ratio between the quatities above
    //marker(nc) = log10( l2_pm1 / l2_p )
    var marker : real = log10(l2_pm1 / l2_p);

    // Apply a sigmoid tranformation to this value and calculate a "scaled" marker
    //psi0 = -one*(c1 + c2*log10(real(cutoff_order,kind=wp)))
    param dpsi : real = 0.5;
    param c1 : real = 4.0;
    param c2 : real = 4.25;
    var psi0 : real = -1.0*(c1 + c2*log10((solDegree-1):real));
    var scaledMarker : real;

    // If the scaled marker is over a threshhold then mark the cell as problematic
    if (marker <= psi0-dpsi) then
      scaledMarker = 0.0;
    else if (marker >= psi0+dpsi) then
      scaledMarker = 1.0;
    else
      scaledMarker = 0.5*(1.0 + sin(half_pi*(marker-psi0)/dpsi));

    return scaledMarker;
  }

  proc troubled_cell_marker(solPoly : [] real, cellTopo : int, solDegree : int) : int
  {
    use Parameters.ParamMesh;
    use Projection;
    use Quadrature;

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
                                     elemTopo = cellTopo,
                                     nodeCnt = solDegree);

      //l2_pm1 = sum( rho_pm1(1:np) * rho_pm1(1:np) * metrics_dt(nc)%jac(1:np) * &
      //              std_elem(this_geom,this_order)%wts(1:np) )
      var l2_diff : real = integrate(nodalPoly = dif2Poly,
                                     elemTopo = cellTopo,
                                     nodeCnt = solDegree);

      // Calculate the Log10 of the ratio between the quatities above
      //marker(nc) = log10( l2_pm1 / l2_p )
      var marker : real = log10(l2_diff / l2_ref);

      // Apply a sigmoid tranformation to this value and calculate a "scaled" marker
      //psi0 = -one*(c1 + c2*log10(real(cutoff_order,kind=wp)))
      param dpsi : real = 0.5;
      param c1 : real = 4.0;
      param c2 : real = 4.25;
      var psi0 : real = -1.0*(c1 + c2*log10((projDegree):real));
      var scaledMarker : real;

      // If the scaled marker is over a threshhold then mark the cell as problematic
      if (marker <= psi0-dpsi) then
        scaledMarker = 0.0;
      else if (marker >= psi0+dpsi) then
        scaledMarker = 1.0;
      else
        scaledMarker = 0.5*(1.0 + sin(half_pi*(marker-psi0)/dpsi));

      // Check if the reference solution is stable and either stop and output the stable degree or keep projectiong down
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
    /*
       Project problematic soulution interpolation to to a lower order
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

    // Create a set with the cell topologies contained in the hypothetics test mesh
    var cellTopos : set(int);
    cellTopos.add(TOPO_LINE);  // Add Line element to the set

    init_polyProj(1, 9, cellTopos);
    init_quadratureWeights(1, 9, cellTopos);

    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    for nodeCnt in 2..4
    {
      writef("%i Node Interpolation\n", nodeCnt);

      var nodes     : [1..nodeCnt] real = nodes_legendre_gauss(nodeCnt);
      var nodalSin  : [1..nodeCnt] real;
      var nodalExp  : [1..nodeCnt] real;
      var nodalLog  : [1..nodeCnt] real;
      var nodalStep : [1..nodeCnt] real;

      var cellSize : real = 0.1;

      for i in nodes.domain
      {
        nodalSin[i] = sin(((randStreamSeeded.getNth(0)-0.5) + nodes[i]*cellSize/2.0)*half_pi);
        nodalExp[i] = exp(((randStreamSeeded.getNth(1)-0.5) + nodes[i]*cellSize/2.0));
        nodalLog[i] = log((randStreamSeeded.getNth(2) + nodes[i]*cellSize/2.0));
        if nodes[i] < 0.0 then nodalStep[i] = 0.0; else nodalStep[i] = 1.0;
      }

      writef("   Limiter: %{8.4dr}  | Sin:  %{ ht}\n", troubled_cell_marker(nodalSin, TOPO_LINE, nodeCnt-1), nodalSin);
      writef("   Limited Polynomial         %{ ht}\n", projection_limiter(nodalSin, TOPO_LINE, nodeCnt-1, nodeCnt-2));
      writef("   Limiter: %{8.4dr}  | Exp:  %{ ht}\n", troubled_cell_marker(nodalExp, TOPO_LINE, nodeCnt-1), nodalExp);
      writef("   Limited Polynomial:        %{ ht}\n", projection_limiter(nodalExp, TOPO_LINE, nodeCnt-1, nodeCnt-2));
      writef("   Limiter: %{8.4dr}  | Log:  %{ ht}\n", troubled_cell_marker(nodalLog, TOPO_LINE, nodeCnt-1), nodalLog);
      writef("   Limited Polynomial         %{ ht}\n", projection_limiter(nodalLog, TOPO_LINE, nodeCnt-1, nodeCnt-2));
      writef("   Limiter: %{8.4dr}  | Step: %{ ht}\n", troubled_cell_marker(nodalStep, TOPO_LINE, nodeCnt-1), nodalStep);
      writef("   Limited Polynomial         %{ ht}\n", projection_limiter(nodalStep, TOPO_LINE, nodeCnt-1, troubled_cell_marker(nodalStep, TOPO_LINE, nodeCnt-1)));
      writef("\n");
    }
  }
}
