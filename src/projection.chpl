/*
  This module initializes and manages the tools required to project orthogonal polynomials to lower orders
*/
module Projection
{
  use UnitTest;
  use Set;

  class projection_coefficients_c
  {
    var coefs_d : domain(3); // {to_order, proj_matrix_row, proj_matrix_col}
    var coefs: [coefs_d] real;
  }

  // Define type for the interpolation structure. Add "?" to allow default initialization to nil.
  type projection_coefficients_t = unmanaged projection_coefficients_c?;

  // Domains
  var polyProj_d : domain(2*int); // {elem_topo, point_dist}

  // Coefficient structures
  var polyProj : [polyProj_d] projection_coefficients_t;

  //////////////////////////////////
  //   Initialization Procedure   //
  //////////////////////////////////

  proc init_polyProj(minPolyDegree : int, maxPolyDegree : int, cellTopos : set(int))
  {
    // We can build a matrix M[i,j] = P_j(x_i) that multiplied by a vector alpha[j] containing a linear combination of the orthogonal
    // polynomial basis results in a array y[x_i] with the value of the interpolated polynomial at the point x_i.
    //
    //  alpha = {  alpha_0 ,  alpha_1 ,  alpha_2 , ... ,  alpha_n }
    //      y = {   y[x_1] ,   y[x_2] ,   y[x_3] , ... ,  y[x_n]  }
    //
    // In 1D M degenerates into a Vandermonde matrix
    //    M_x = { P_0(x_1) , P_1(x_1) , P_2(x_1) , ... , P_n(x_1)
    //            P_0(x_2) , P_1(x_2) , P_2(x_2) , ... , P_n(x_2)
    //            P_0(x_3) , P_1(x_3) , P_2(x_3) , ... , P_n(x_3)
    //              ...    ,   ...    ,   ...    , ... ,   ...
    //            P_0(x_n) , P_1(x_n) , P_2(x_n) , ... , P_n(x_n) }
    //
    // matmul( M[i,j] , alpha[i] ) = y[j]
    //
    // Store inverted M so that alpha = matmul( M^(-1) , j[j] ) is easy to calculate

    use Time;
    use Parameters.ParamMesh;
    use Polynomials;
    use Quadrature;

    writeln();
    writeln("Initializing Polynomial Projection matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minPolyDegree);
    writeln("    Maximum Polynomial Degree: ", maxPolyDegree);
    var stopwatch : Timer;
    stopwatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for polyDegree in minPolyDegree..maxPolyDegree do
        polyProj_d.add((cellTopo, polyDegree));

    // Check if the quadrature weights structure exists for this element and initialize it if not
    for cellTopo in cellTopos do
      for polyDegree in minPolyDegree..maxPolyDegree do
        if quadratureWeights.domain.contains((cellTopo, polyDegree))
        {
          if quadratureWeights[(cellTopo, polyDegree)] == nil then
            init_quadratureWeights(0, maxPolyDegree, cellTopos);
        }
        else
        {
          init_quadratureWeights(0, maxPolyDegree, cellTopos);
        }

    for (cellTopo, polyDegree) in polyProj.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = polyDegree+1;

          polyProj[(cellTopo, polyDegree)] = new projection_coefficients_t({0..polyDegree, 1..spCnt, 1..spCnt})!;

          // Evaluate each component "n" of the polynomial basis at the interpolation points "i"
          var nodalBasis : [1..spCnt, 1..spCnt] real = nodal_basis_line(polyDegree, polyDegree);

          for projDegree in 0..polyDegree by -1 do
            polyProj[(cellTopo, polyDegree)]!.coefs[projDegree, .., ..] = proj_matrix(nodalBasis[.., 1..projDegree+1],
                quadratureWeights[(cellTopo, polyDegree)]!.weights);
        }
        when TOPO_TRIA
        {
          var spCnt : int = ((polyDegree+1)*(polyDegree+2))/2;

          var nodalBasis : [1..spCnt, 0..(polyDegree+1)**2-1] real;
        }
        when TOPO_QUAD
        {
          var spCnt : int = (polyDegree+1)**2;

          polyProj[(cellTopo, polyDegree)] = new projection_coefficients_t({0..polyDegree, 1..spCnt, 1..spCnt})!;

          // Evaluate each component "n" of the polynomial basis at the interpolation points (i,j)
          var nodalBasis : [1..spCnt, 0..(polyDegree+1)**2-1] real = nodal_basis_quad(polyDegree, polyDegree);

          for projDegree in 0..polyDegree by -1 do
            polyProj[(cellTopo, polyDegree)]!.coefs[projDegree, .., ..] = proj_matrix(nodalBasis[.., 0..(projDegree+1)**2-1],
                quadratureWeights[(cellTopo, polyDegree)]!.weights);
        }
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at projection initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", stopwatch.elapsed(TimeUnits.milliseconds));
  }

  //////////////////////////
  //   Projection Basis   //
  //////////////////////////

  proc nodal_basis_line(polyDegree : int, interpDegree : int)
  {
    // Calculate nodal values in several points for each component of the basis
    use Polynomials;

    var nodeCnt : int = interpDegree+1;
    var modeCnt : int = polyDegree+1;

    var nodeDistLine : [1..polyDegree+1] real = nodes_legendre_gauss(polyDegree+1);

    var nodalBasis : [1..nodeCnt, 0..modeCnt-1] real;

    // Calculate nodal values in several points for each component of the basis
    for nodeIdx in 1..nodeCnt do // Loop through nodes
      for modeIdx in 0..polyDegree do // Loop through modes
        nodalBasis[nodeIdx, modeIdx] = eval_legendre_poly(modeIdx, nodeDistLine[nodeIdx]);

    return nodalBasis;
  }

  proc nodal_basis_tria(polyDegree : int, interpDegree : int)
  {
    // Calculate nodal values in several points for each component of the basis
    use Polynomials;

    var nodeCnt : int = ((interpDegree+1)*(interpDegree+2))/2;
    var modeCnt : int = ((polyDegree+1)*(polyDegree+2))/2;

    var nodeDistLine : [1..polyDegree+1] real = nodes_legendre_gauss(polyDegree+1);

    var nodalBasis   : [1..nodeCnt, 0..modeCnt-1] real;

    // Loop through nodes
    for nodeIdx in 1..nodeCnt
    {
      // Loop through modes
      for diag in 0..modeCnt-1 do
        for n in 0..diag
        {
          var modeIdx : int = (diag*(diag+1))/2 + n;
          //nodalBasis[nodeIdx, modeIdx] = eval_legendre_poly(diag-n, nodeDistLine[])
          //                              *eval_legendre_poly(     n, nodeDistLine[]);
        }
    }

    return nodalBasis;
  }

  proc nodal_basis_quad(polyDegree : int, interpDegree : int)
  {
    // Calculate nodal values in several points for each component of the basis
    use Polynomials;

    var nodeCnt : int = (interpDegree+1)**2;
    var modeCnt : int = (polyDegree+1)**2;

    var nodeDistLine : [1..polyDegree+1] real = nodes_legendre_gauss(polyDegree+1);

    var nodalBasis   : [1..nodeCnt, 0..modeCnt-1] real;

    // Loop through nodes
    for nodeIdx in 1..nodeCnt
    {
      // Loop through modes
      for degree in 0..polyDegree
      {
        var modeIdx : int = (degree+1)**2-1;
        nodalBasis[nodeIdx, modeIdx] = eval_legendre_poly(degree, nodeDistLine[(nodeIdx-1)/(polyDegree+1)+1])
                                      *eval_legendre_poly(degree, nodeDistLine[(nodeIdx-1)%(polyDegree+1)+1]);

        for aux in 0..degree-1
        {
          modeIdx = degree**2 + 2*aux;
          nodalBasis[nodeIdx, modeIdx] = eval_legendre_poly(degree, nodeDistLine[(nodeIdx-1)/(polyDegree+1)+1])
                                        *eval_legendre_poly(   aux, nodeDistLine[(nodeIdx-1)%(polyDegree+1)+1]);

          modeIdx = degree**2 + 2*aux + 1;
          nodalBasis[nodeIdx, modeIdx] = eval_legendre_poly(   aux, nodeDistLine[(nodeIdx-1)/(polyDegree+1)+1])
                                        *eval_legendre_poly(degree, nodeDistLine[(nodeIdx-1)%(polyDegree+1)+1]);
        }
      }
    }

    return nodalBasis;
  }

  proc nodal_basis_tetr() {}
  proc nodal_basis_pyra() {}
  proc nodal_basis_pris() {}
  proc nodal_basis_hexa() {}

  ///////////////////////////
  //   Projection Matrix   //
  ///////////////////////////

  proc proj_matrix(basis : [] real, weights : [] real, normalizedBasis : bool = false)
  {
    use LinearAlgebra;

    // basis : [1..spCnt, 1..modeCnt] real
    //      Matrix of the basis vectors of the projection space in nodal form
    //      Values of each polynomial mode at each SP
    //
    // weights : [1..sps] real
    //      Gaussian quadrature weights for the abscissa used to represent the basis vectors

    // ProjectionMatrix(B) = B * (B^T * W * B)^{-1} * B^T * W
    //      Where W is a diagonal matrix with the quadrature weights


    // B^T * W
    var BtW : [basis.T.domain] real;
    for i in weights.domain do
      BtW[..,i] = basis[i,..] * weights[i];

    // If the basis is not normalized then normalize projection
    if !normalizedBasis
    {
      // Inverse of the internal product of the basis vectors with themselves. This should be a diagonal matrix with the
      // inverse of the norms of the basis vectors unless the basis is not orthogonal.
      var invNorm : [basis.domain.dim(1), basis.domain.dim(1)] real = inv(dot(BtW, basis));

      // Apply the normalization
      BtW = dot(invNorm, BtW);
    }

    // Pre multiply the projection matrix by the nodal form of the basis so that the projection result shall be in the
    // same form.
    var projMatrix : [weights.domain.dim(0), weights.domain.dim(0)] real = dot(basis, BtW);

    return projMatrix;
  }

  /////////////////////////////
  //   Projection Function   //
  /////////////////////////////

  proc project_poly(nodalPoly : [] real, cellTopo : int, polyDegree : int, projDegree : int)
  {
    use Parameters.ParamMesh;
    use LinearAlgebra;

    var projection : [nodalPoly.domain] real;

    select cellTopo
    {
      when TOPO_LINE
      {
        projection = dot(polyProj[(cellTopo, polyDegree)]!.coefs[projDegree, .., ..], nodalPoly.T);
      }
      when TOPO_TRIA {}
      when TOPO_QUAD
      {
        projection = dot(polyProj[(cellTopo, polyDegree)]!.coefs[projDegree, .., ..], nodalPoly.T);
      }
      when TOPO_TETR {}
      when TOPO_PYRA {}
      when TOPO_PRIS {}
      when TOPO_HEXA {}
      otherwise do writeln("Unsupported mesh element found at polynomial projection.");
    }

    return projection;
  }

  ///////////////////////////////
  //   Module Test Procedure   //
  ///////////////////////////////

  proc main()
  {
    /*
       Generate a random n-th degree polynomial by way of a random linear combination of Legendre Polynomials
    */

    use Random;
    use Set;
    use IO.FormattedIO;
    use Testing;
    use Parameters.ParamMesh;
    use Parameters.ParamTest;
    use Quadrature;
    use Polynomials;
    use LinearAlgebra;

    var randStream = new RandomStream(real);
    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    var minPolyDegree : int = 0;
    var maxPolyDegree : int = 9;

    // Create a set with the cell topologies contained in the hypothetical test mesh
    var cellTopos : set(int);
    cellTopos.add(TOPO_LINE);
    cellTopos.add(TOPO_QUAD);

    // Initialize projection matrix structure
    writeln();
    writeln("Initializing projection module for FR (polyProj):");
    init_polyProj(minPolyDegree, maxPolyDegree, cellTopos);

    // Print initialized projection structure
    writeln();
    writeln("Interpolation initialized structure for FR (polyProj):");
    writeln(polyProj);
    writeln();

    writeln();
    writeln("Test projecting randomly generated polynomials");

    // Define test polynomial parameters
    var polyDegree : int = 9;

    writeln();
    writeln("------------------------------------------------------------");
    writeln();
    writeln("Topo Line:");
    {
      var nodeCnt : int = polyDegree+1;
      var nodes   : [1..nodeCnt] real = nodes_legendre_gauss(nodeCnt);

      // Generate random polynomial in modal form
      var modalPoly : [0..polyDegree] real;
      randStreamSeeded.fillRandom(modalPoly);

      // Calculate nodal values in several points for each component of the basis
      var nodalBasis : [1..nodeCnt, 0..polyDegree] real;
      nodalBasis = nodal_basis_line(polyDegree, polyDegree);

      // Calculate the nodal representation of the random test polynomial
      var nodalPoly : [1..nodeCnt] real = dot(nodalBasis, modalPoly);

      // Write out test parameters
      {
        writeln();
        writef("Nodes: ");
        for i in 1..nodeCnt do
          writef(" %{13.6er}", nodes[i]);
        writef("\n");

        writef("Modes: ");
        for i in 0..polyDegree do
          writef(" %{13.6er}", modalPoly[i]);
        writef("\n");

        writeln();
        writef("Nodal basis:\n");
        for j in 0..polyDegree
        {
          writef("   Legendre(%2i): ", j);
          for i in 1..nodeCnt do
            writef(" %{13.6er}", nodalBasis[i,j]);
          writef("\n");
        }
      }

      // Test projections to every degree down to 0
      for projDegree in (0..polyDegree by -1)
      {
        // Print nodal values of projected polynomial
        writeln();
        writef("Projection of a degree %i polynomial to %i degree:\n", polyDegree, projDegree);

        {
          var nodalSum : [1..polyDegree+1] real = dot(nodalBasis[1..polyDegree+1, 0..projDegree], modalPoly[0..projDegree]);

          var nodalProj : [1..polyDegree+1] real = project_poly(nodalPoly = nodalPoly,
                                                         cellTopo  = TOPO_LINE,
                                                         polyDegree = polyDegree,
                                                         projDegree = projDegree);

          writef("   Node | Sum of modes     | Nodal Projection | Absolute Error | Relative Error\n");

          for i in 1..polyDegree+1 do
            writef("   %4i | %{16.8er} | %{16.8er} | %{14.6er} | %{14.6er}\n", i, nodalSum[i], nodalProj[i],
                error(nodalSum[i], nodalProj[i]), relative_error(nodalSum[i], nodalProj[i]));
        }
      }
    }

    writeln();
    writeln("------------------------------------------------------------");
    writeln();
    writeln("Topo Quad:");
    {
      var nodeCnt : int = (polyDegree+1)**2;
      var modeCnt : int = (polyDegree+1)**2;

      var nodeDistLine : [1..polyDegree+1] real = nodes_legendre_gauss(polyDegree+1);

      // Generate random polynomial in modal form
      var modalPoly : [0..(polyDegree+1)**2-1] real;
      randStreamSeeded.fillRandom(modalPoly);

      // Calculate nodal values in several points for each component of the basis
      var nodalBasis : [1..nodeCnt, 0..modeCnt-1] real;
      nodalBasis = nodal_basis_quad(polyDegree, polyDegree);

      // Calculate the nodal representation of the random test polynomial
      var nodalPoly : [1..nodeCnt] real = dot(nodalBasis, modalPoly);

      // Write out test parameters
      {
        writeln();
        writef("Nodes:\n");
        for i in 1..polyDegree+1 by -1
        {
          for j in 1..polyDegree+1 do
            writef("   (%{7.4dr}, %{7.4dr})", nodeDistLine[j], nodeDistLine[i]);
          writef("\n");
        }

        writeln();
        writef("Mode weights:\n");
        for n in 0..(polyDegree+1)**2-1 do
          writef(" (%2i,%2i) | %{13.6er}\n", n/(polyDegree+1), n%(polyDegree+1), modalPoly[n]);

        writeln();
        writef("Nodal basis:\n");
        for majorDegree in 0..polyDegree
        {
          for n in 0..majorDegree-1
          {
            writef("   Legendre(%2i,%2i):", majorDegree, n);
            for i in 1..polyDegree+1 by -1
            {
              for j in 1..polyDegree+1 do
                writef(" %{13.6er}", nodalBasis[(i-1)*(polyDegree+1)+j, majorDegree**2 + 2*n]);
              writef("\n                   ");
            }
            writeln();

            writef("   Legendre(%2i,%2i):", n, majorDegree);
            for i in 1..polyDegree+1 by -1
            {
              for j in 1..polyDegree+1 do
                writef(" %{13.6er}", nodalBasis[(i-1)*(polyDegree+1)+j, majorDegree**2 + 2*n + 1]);
              writef("\n                   ");
            }
            writeln();
          }
          writef("   Legendre(%2i,%2i):", majorDegree, majorDegree);
          for i in 1..polyDegree+1 by -1
          {
            for j in 1..polyDegree+1 do
              writef(" %{13.6er}", nodalBasis[(i-1)*(polyDegree+1)+j, (majorDegree+1)**2-1]);
            writef("\n                   ");
          }
          writeln();
        }
      }

      // Test projections to every degree down to 0
      for projDegree in (0..polyDegree by -1)
      {
        // Print nodal values of projected polynomial
        writeln();
        writef("Projection of a degree %i polynomial to %i degree:\n", polyDegree, projDegree);

        {
          var nodalSum : [1..nodeCnt] real = dot(nodalBasis[1..nodeCnt, 0..(projDegree+1)**2-1],
                                                 modalPoly[0..(projDegree+1)**2-1]  );

          var nodalProj : [1..nodeCnt] real = project_poly(nodalPoly  = nodalPoly,
                                                           cellTopo   = TOPO_QUAD,
                                                           polyDegree = polyDegree,
                                                           projDegree = projDegree);

          writef("   Node    | Sum of modes     | Nodal Projection | Absolute Error | Relative Error\n");
          for ij in 1..nodeCnt do
            writef("   (%2i,%2i) | %{16.8er} | %{16.8er} | %{14.6er} | %{14.6er}\n",
                (ij-1)/(polyDegree+1)+1,
                (ij-1)%(polyDegree+1)+1,
                nodalSum[ij],
                nodalProj[ij],
                error(nodalSum[ij], nodalProj[ij]),
                relative_error(nodalSum[ij], nodalProj[ij])
            );
          writef("\n");
        }
      }
    }
  }
}
