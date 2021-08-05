/*
  This module initializes and manages the tools required to project orthogonal polynomials to lower orders
*/
prototype module Projection
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

    use Parameters.ParamMesh;
    use Polynomials;
    use Quadrature;

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for polyDegree in minPolyDegree..maxPolyDegree do
        polyProj_d.add((cellTopo, polyDegree));

    for (cellTopo, polyDegree) in polyProj.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = polyDegree+1;
          var nodes   : [1..spCnt] real = nodes_legendre_gauss(spCnt);
          var weights : [1..spCnt] real = weights_legendre_gauss(spCnt);

          polyProj[(cellTopo, polyDegree)] = new projection_coefficients_t({0..polyDegree, 1..spCnt, 1..spCnt})!;

          // calculate nodal values in several points for each component of the basis
          var nodalbasis : [1..spCnt, 1..spCnt] real;
          for i in 1..spCnt do
            for j in 1..spCnt do
              nodalbasis[i,j] = eval_legendre_poly(j-1, nodes[i]);

          for projDegree in 0..polyDegree by -1 do
            polyProj[(cellTopo, polyDegree)]!.coefs[projDegree, .., ..] = proj_matrix(nodalbasis[.., 1..projDegree+1], weights);
        }
        when TOPO_TRIA {}
        when TOPO_QUAD {}
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at projection initialization.");
      }
    }
  }

  proc proj_matrix(basis : [] real, weights : [] real, normalizedBasis : bool = false)
  {
    use LinearAlgebra;

    // basis : [1..spCnt, 1..modeCnt] real
    //      Matrix of the basis vectors of the projection space in nodal form
    //
    // weights : [1..sps] real
    //      Gaussian quadrature weights for the abscissa used to represent the basis vectors

    // ProjectionMatrix(B) = B * (B^T * W * B)^{-1} * B^T * W


    // B^T * W
    var BtW : [basis.T.domain] real;
    for i in weights.domain do
      BtW[..,i] = basis[i,..] * weights[i];

    // If the basis is not normalized then normalize projection
    if !normalizedBasis
    {
      // Inverse of the tnternal product of the basis vectors with themselfes. This shoud be a diagonal matrix with the
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

  proc project_poly(nodalPoly : [1..polyDegree+1] real, cellTopo : int, polyDegree : int, projDegree : int)
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
      when TOPO_QUAD {}
      when TOPO_TETR {}
      when TOPO_PYRA {}
      when TOPO_PRIS {}
      when TOPO_HEXA {}
      otherwise do writeln("Unsupported mesh element found at polynomial projection.");
    }

    return projection;
  }

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

    // Create a set with the cell topologies contained in the hypothetics test mesh
    var cellTopos : set(int);
    cellTopos.add(2);  // Add Line element to the set

    // Intitialize projection matrix structure
    writeln();
    writeln("Interpolation initialized structure for FR (polyProj):");
    writeln();
    init_polyProj(minPolyDegree, maxPolyDegree, cellTopos);

    // Define test polynomial parameters
    var polyDegree : int = 9;
    var spCnt : int = polyDegree+1;

    var nodes   : [1..spCnt] real = nodes_legendre_gauss(spCnt);
    var weights : [1..spCnt] real = weights_legendre_gauss(spCnt);

    // Generate random polynomial in modal form
    var modalPoly : [0..polyDegree] real;
    randStreamSeeded.fillRandom(modalPoly);

    // Calculate nodal values in several points for each component of the basis
    var nodalBasis : [1..spCnt, 0..polyDegree] real;
    for i in 1..spCnt do
      for j in 0..polyDegree do
        nodalBasis[i,j] = eval_legendre_poly(j, nodes[i]);

    // Calculate the nodal representation of the random test polynomial
    var nodalPoly : [1..spCnt] real = dot(nodalBasis, modalPoly);

    // Write out test parameters
    writef("Nodes: ");
    for i in 1..spCnt do
      writef(" %{13.6er}", nodes[i]);
    writef("\n");

    writef("Modes: ");
    for i in 0..polyDegree do
      writef(" %{13.6er}", modalPoly[i]);
    writef("\n\n");

    writef("Nodal basis:\n");
    for j in 0..polyDegree
    {
      writef("   Legendre(%2i): ", j);
      for i in 1..spCnt do
        writef(" %{13.6er}", nodalBasis[i,j]);
      writef("\n");
    }
    writef("\n");

    // Test projections to every degree down to 0
    for projDegree in (0..polyDegree by -1)
    {
      var nodalSum : [1..spCnt] real = dot(nodalBasis[1..spCnt, 0..projDegree], modalPoly[0..projDegree]);

      var nodalProj : [1..spCnt] real = project_poly(nodalPoly = nodalPoly,
                                                     cellTopo = TOPO_LINE,
                                                     polyDegree = polyDegree,
                                                     projDegree = projDegree);

      // Print nodal values of projected polynomial
      writef("Projection of a degree %i polynomial to %i degree\n", polyDegree, projDegree);
      writef("Node |   Sum of modes   | Nodal Projection | Absolute Error | Relative Error\n");
      for i in 1..spCnt do
        writef("%4i | %{16.8er} | %{16.8er} | %{14.6er} | %{14.6er}\n", i, nodalSum[i], nodalProj[i],
            error(nodalSum[i], nodalProj[i]), relative_error(nodalSum[i], nodalProj[i]));
      writef("\n");
    }
  }
}
