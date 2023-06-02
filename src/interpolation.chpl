module Interpolation
{
  use Random;
  use UnitTest;
  use Set;

  class interpolation_coefficients_c
  {
    var coefs_d : domain(2); // {nFPs, nSPs}
    var coefs: [coefs_d] real;
  }

  class derivation_coefficients_c
  {
    var coefs_d : domain(3); // {nFPs, nDims, nSPs}
    var coefs: [coefs_d] real;
  }

  // Define type for the interpolation structure. Add "?" to allow default initialization to nil.
  type interpolation_coefficients_t = unmanaged interpolation_coefficients_c?;
  type derivation_coefficients_t    = unmanaged derivation_coefficients_c?;

  // Domains
  var sp2fpInterp_d   : domain(2*int);
  var sp2spDeriv_d    : domain(2*int);
  var sp2nodeInterp_d : domain(2*int);

  // Coefficient structures
  var sp2fpInterp   : [sp2fpInterp_d] interpolation_coefficients_t; // Cell to face interpolation
  var sp2spDeriv    : [sp2spDeriv_d]  derivation_coefficients_t;    // Cell to face derivative
  var sp2nodeInterp : [sp2spDeriv_d]  interpolation_coefficients_t; // Cell to node interpolation

  // Alternative coefficient structures.
  //    1- It's assumed interpolation is always from the cellTopo SPs to somewhere else.
  //    2- Interpolating to a lower order is possible but obviously does not preserve order of accuracy.
  //    3- fromPtIdx is typically the cell SP index
  //    4-   toPtIdx can be either the cell vertices, the face FPs, or an arbitrary distribution of points in the edge.
  //
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toCell[          solOrder ].coefs[toPtIdx, fromPtIdx]
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toFace[(faceIdx, solOrder)].coefs[toPtIdx, fromPtIdx]
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toEdge[(edgeIdx, solOrder)].coefs[toPtIdx, fromPtidx]
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toNode[(nodeIdx, solOrder)].coefs[toPtIdx, fromPtidx]
  //
  // Ex: interpolate[(TOPO_PRIS, 4)].toFace[(1, 4)].

  //////////////////////////////////
  //   Initialization Procedure   //
  //////////////////////////////////

  proc init_sp2fpInterp(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Polynomials;

    writeln();
    writeln("Initializing SP -> FP Interpolation matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var initWatch : stopwatch;
    initWatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2fpInterp_d.add((cellTopo, interpOrder));

    // Calculate all relevant coefficients
    for (cellTopo, interpOrder) in sp2fpInterp.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = interpOrder+1;
          var fpCnt : int = 2;

          sp2fpInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..fpCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spDistLine : [1..spCnt] real = nodes_legendre_gauss(spCnt);
          var fpDistLine : [1..fpCnt] real = [-1.0, 1.0];

          for fp in 1..fpCnt do
            sp2fpInterp[(cellTopo, interpOrder)]!.coefs[{fp.. #1, 1..spCnt}] = reshape(
                eval_LagrangePoly1D_array(fpDistLine[fp], spDistLine), {fp.. #1, 1..spCnt});
        }
        when TOPO_TRIA {}
        when TOPO_QUAD
        {
          var spCnt : int = (interpOrder+1)**2;
          var fpCnt : int = (interpOrder+1)*4;

          sp2fpInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..fpCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
          var fpDistLine : [1..2] real = [-1.0, 1.0];

          for cellFace in 1..4 do
            for faceFP in 1..interpOrder+1
            {
              var  xiFP : real;
              var etaFP : real;

              // Get xi and eta depending on the cell face
              select cellFace
              {
                when 1
                {
                  xiFP  = spDistLine[faceFP];
                  etaFP = fpDistLine[1];
                }
                when 2
                {
                  xiFP  = fpDistLine[2];
                  etaFP = spDistLine[faceFP];
                }
                when 3
                {
                  xiFP  = spDistLine[interpOrder+2-faceFP];
                  etaFP = fpDistLine[2];
                }
                when 4
                {
                  xiFP  = fpDistLine[1];
                  etaFP = spDistLine[interpOrder+2-faceFP];
                }
              }

              for sp in 1..spCnt
              {
                var i : int = (sp-1)%(interpOrder+1)+1;
                var j : int = (sp-1)/(interpOrder+1)+1;
                var fp : int = faceFP+(cellFace-1)*(interpOrder+1);

                sp2fpInterp[(cellTopo, interpOrder)]!.coefs[fp, sp] = eval_LagrangePoly1D( xiFP, i, spDistLine)
                                                                     *eval_LagrangePoly1D(etaFP, j, spDistLine);
              }
            }
        }
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", initWatch.elapsed()*1000);
  }

  proc init_sp2spDeriv(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Polynomials;

    writeln();
    writeln("Initializing SP -> SP' Differentiation matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var initWatch : stopwatch;
    initWatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2spDeriv_d.add((cellTopo, interpOrder));

    // Calculate all relevant coefficients
    for (cellTopo, interpOrder) in sp2spDeriv.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = interpOrder+1;

          sp2spDeriv[(cellTopo, interpOrder)] = new derivation_coefficients_t({1..spCnt, 1..1, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spLoc : [1..spCnt] real = nodes_legendre_gauss(spCnt);

          for sp in 1..spCnt do
            sp2spDeriv[(cellTopo, interpOrder)]!.coefs[{sp.. #1, 1.. #1, 1..spCnt}] =
                  reshape(eval_DLagrangeDx_array(spLoc[sp], spLoc), {sp.. #1, 1.. #1, 1..spCnt});
        }
        when TOPO_TRIA {}
        when TOPO_QUAD
        {
          var spCnt : int = (interpOrder+1)**2;

          sp2spDeriv[(cellTopo, interpOrder)] = new derivation_coefficients_t({1..spCnt, 1..2, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);

          for toSP in 1..spCnt
            {
              var  xi : real = spDistLine[(toSP-1)%(interpOrder+1)+1];
              var eta : real = spDistLine[(toSP-1)/(interpOrder+1)+1];

              for fromSP in 1..spCnt
              {
                var i : int = (fromSP-1)%(interpOrder+1)+1;
                var j : int = (fromSP-1)/(interpOrder+1)+1;

                sp2spDeriv[(cellTopo, interpOrder)]!.coefs[toSP, 1, fromSP] = eval_DLagrangeDx(xi, i, spDistLine)
                                                                             *eval_LagrangePoly1D(eta, j, spDistLine);

                sp2spDeriv[(cellTopo, interpOrder)]!.coefs[toSP, 2, fromSP] = eval_LagrangePoly1D(xi, i, spDistLine)
                                                                             *eval_DLagrangeDx(eta, j, spDistLine);
              }
            }
        }
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", initWatch.elapsed()*1000);
  }

  proc init_sp2nodeInterp(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Polynomials;

    writeln();
    writeln("Initializing SP -> Node Differentiation matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var initWatch : stopwatch;
    initWatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2nodeInterp_d.add((cellTopo, interpOrder));

    // Calculate all relevant coefficients
    for (cellTopo, interpOrder) in sp2nodeInterp.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = interpOrder+1;
          var nodeCnt : int = 2;//elem_vertices(cellTopo);

          sp2nodeInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..nodeCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spDistLine : [1..spCnt] real = nodes_legendre_gauss(spCnt);
          var nodeDistLine : [1..nodeCnt] real = [-1.0, 1.0];

          for node in 1..nodeCnt do
            sp2nodeInterp[(cellTopo, interpOrder)]!.coefs[{node.. #1, 1..spCnt}] = reshape(
                eval_LagrangePoly1D_array(nodeDistLine[node], spDistLine), {node.. #1, 1..spCnt});
        }
        when TOPO_TRIA {}
        when TOPO_QUAD
        {
          var spCnt : int = (interpOrder+1)**2;
          var nodeCnt : int = 4;//elem_vertices(cellTopo);

          sp2nodeInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..nodeCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
          var nodeDistLine : [1..2] real = [-1.0, 1.0];

          for cellNode in 1..nodeCnt
          {
            var  xiFP : real;
            var etaFP : real;

            select cellNode
            {
              when 1
              {
                xiFP  = nodeDistLine[1];
                etaFP = nodeDistLine[1];
              }
              when 2
              {
                xiFP  = nodeDistLine[2];
                etaFP = nodeDistLine[1];
              }
              when 3
              {
                xiFP  = nodeDistLine[2];
                etaFP = nodeDistLine[2];
              }
              when 4
              {
                xiFP  = nodeDistLine[1];
                etaFP = nodeDistLine[2];
              }
            }

            for sp in 1..spCnt
            {
              var i : int = (sp-1)%(interpOrder+1)+1;
              var j : int = (sp-1)/(interpOrder+1)+1;

              sp2nodeInterp[(cellTopo, interpOrder)]!.coefs[cellNode, sp] = eval_LagrangePoly1D( xiFP, i, spDistLine)
                                                                           *eval_LagrangePoly1D(etaFP, j, spDistLine);
            }
          }
        }
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", initWatch.elapsed()*1000);
  }

  //////////////////////////////////////////
  //   Lagrange Interpolation functions   //
  //////////////////////////////////////////

  proc eval_LagrangePoly1D(x : real, k : int, xi : [] real) : real
  {
    // Evaluate the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[] at the point x.
    //
    //               nfp
    //              _____
    //               | |     x  - xi(i)
    //    L_k (x) =  | |  -------------
    //               i=1  xi(k) - xi(i)
    //              i/=k
    //
    //  x  : Coordinate at which to evaluate the k-th Lagrange basis polynomial
    //  k  : Index of the basis vector being evaluated
    //  xi : Array of all the interpolation nodes coordinates
    //
    // The k-th basis vector of a Lagrange basis is a polynomial function that evaluates to 0 at all nodes xi[] except
    // at xi[k] where it evaluates to 1.

    var return_value : real = 1.0;

    for i in xi.domain do
      if i != k then
        return_value *= (x-xi[i]) / (xi[k]-xi[i]);

    return return_value;
  }

  proc eval_LagrangePoly1D_array(x : real, xi : [] real) : [] real
  {
    // Get the values of the 1D Lagrange polynomials at point x
    //
    //  x  : Coordinate at which to evaluate the k-th Lagrange basis polynomial
    //  xi : Array of all the interpolation nodes coordinates
    //
    // This evaluates the Lagrange polynomial at the location x. The Lagrange polynomial corresponds to a function that
    // is 0.0 at all iterpolation points xi(i) for i=1,nfp except at xi(k) where the function is equal to 1.

    var return_value : [xi.domain] real;

    for k in xi.domain do
      return_value[k] = eval_LagrangePoly1D(x, k, xi);

    return return_value;
  }

  proc eval_DLagrangeDx(x : real, k : int, xi : [] real) : real
  {
    // Evaluate the derivative of the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[] at the
    // point x.
    //
    //                    nfp                  nfp
    //                    ___                 _____
    //     d              \          1         | |     x  - xi(j)
    //    -- [L_k (x)] =  /__  -------------   | |  -------------
    //    dx                   xi(k) - xi(i)   j=1  xi(k) - xi(j)
    //                    i=1                 j/=i
    //                   i/=k                 j/=k
    //
    //  x  : Coordinate at which to evaluate the derivative of the k-th Lagrange basis polynomial
    //  k  : Index of the basis vector being evaluated
    //  xi : Array of all the interpolation nodes coordinates

    var evalDLagrangeDx : real = 0.0;
    var xiSliced : [xi.domain.low.. #(xi.size-1)] real;

    // Populate xiSliced with all nodes except the k-th
    for i in xiSliced.domain do
      if i<k then
        xiSliced[i] = xi[i];
      else
        xiSliced[i] = xi[i+1];

    for i in xiSliced.domain
    {
      var aux = xiSliced[i];
      xiSliced[i] = xi[k];
      evalDLagrangeDx += eval_LagrangePoly1D(x, i, xiSliced) / (xiSliced(i)-aux);
      xiSliced[i] = aux;
    }

    return evalDLagrangeDx;
  }

  proc eval_DLagrangeDx_array(x : real, xi : [] real) : [] real
  {
    // Evaluate the derivative of the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[] at the
    // point x.
    //
    //  x  : Coordinate at which to evaluate the derivative of the k-th Lagrange basis polynomial
    //  xi : Array of all the interpolation nodes coordinates

    var evalDLagrangeDxArray : [xi.domain] real;

    for k in xi.domain do
      evalDLagrangeDxArray[k] = eval_DLagrangeDx(x, k, xi);

    return evalDLagrangeDxArray;
  }

  proc eval_D2LagrangeDx2(x : real, k : int, xi : real) : real
  {
    // Evaluate the second derivative of the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[]
    // at the point x.
    //
    //                      nfp                nfp                  nfp
    //                      ___                ___                 _____
    //     d²               \         1        \          1         | |     x  - xi(n)
    //    --- [L_k (x)]  =  /__ -------------  /__  -------------   | |  -------------
    //    dx²                   xi(k) - xi(l)       xi(k) - xi(m)   n=1  xi(k) - xi(n)
    //                      l=1                m=1                 n/=l
    //                     l/=k               m/=l                 n/=m
    //                                        m/=k                 n/=k
    //
    //  x  : Coordinate at which to evaluate the second derivative of the k-th Lagrange basis polynomial
    //  k  : Index of the basis vector being evaluated
    //  xi : Array of all the interpolation nodes coordinates

    var evalD2LagrangeDx2 : real = 0.0;
    var xiSliced : [xi.domain.low.. #(xi.size-1)] real;

    // Populate xiSliced with all nodes except the k-th
    for i in xiSliced.domain do
      if i<k then
        xiSliced[i] = xi[i];
      else
        xiSliced[i] = xi[i+1];

    for l in xi.domain do
      if l != k then
        evalD2LagrangeDx2 += eval_DLagrangeDx(x, k, xiSliced) / (xi(k)-xi(l));

    return evalD2LagrangeDx2;
  }

  proc eval_D2LagrangeDx2_array(x : real, k : int, xi : real) : [xi.domain] real
  {
    // Get value of second derivative of k-th Lagrange polynomial at point x.
    //
    //  x  : Coordinate at which to evaluate the second derivative of the k-th Lagrange basis polynomial
    //  xi : Array of all the interpolation nodes coordinates

    var evalD2LagrangeDx2Array : [xi.domain] real = 0.0;

    for k in xi.domain do
      evalD2LagrangeDx2Array[k] = eval_D2LagrangeDx2(x, k, xi[xi.domain]);

    return evalD2LagrangeDx2Array;
  }

  proc eval_LagrangePoly2D(x : [1..2] real, k : [1,.2] int, xi : [] real, eta : [] real) : real
  {
    return  eval_LagrangePoly1D(x[1], k[1],  xi[])
           *eval_LagrangePoly1D(x[2], k[2], eta[]);
  }

  proc eval_LagrangePoly3D(x : [1..3] real, k : [1..2] int, xi : [] real, eta : [] real, zeta : [] real) : real
  {
    return  eval_LagrangePoly1D(x[1], k[1],   xi[])
           *eval_LagrangePoly1D(x[2], k[2],  eta[])
           *eval_LagrangePoly1D(x[3], k[3], zeta[]);
  }

  ///////////////////////////////
  //   Module Test Procedure   //
  ///////////////////////////////

  proc main()
  {
    use IO;
    use Testing;
    use Parameters.ParamTest;
    use Parameters.ParamMesh;
    use Polynomials;

    var randStream = new RandomStream(real);
    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    var interpolation : real;
    var interpDeriv   : real;
    var node          : [0..9] real;
    var x             : [0..9] real;
    var y_node        : [0..9] real = 0;
    var y_x           : [0..9] real = 0;
    var coef          : [0..9] real;
    var basis         : [0..9, 0..9] real;
    var basisDeriv    : [0..9, 0..9] real;

    writeln();

    // Get Chebyshev roots to use as interpolation nodes
    node = nodes_legendre_gauss(10);
    writeln("Interpolation nodes (Legendre-Gauss):         ", node);
    for i in 0..9 do
      node[i] = -cos( half_pi * (2*i+1)/10 );
    writeln("Normalized interpolation nodes [-1.0 to 1.0]: ", node);
    writeln();

    // Get random interpolation targets [-1.0, 1.0] and random polynomials
    x[0] = 0;
    x[9] = 1;
    randStreamSeeded.fillRandom(x[1..8]);
    randStreamSeeded.fillRandom(coef);
    x = 2.0*x-1.0;
    coef = 2.0*coef -1.0;
    writeln("Interpolation Targets   = ", x);
    writeln("Polynomial Coefficients = ", coef);
    writeln();

    // Create a set with the cell topologies contained in the hypothetical test mesh
    var cellTopos : set(int);
    cellTopos.add(TOPO_LINE);
    cellTopos.add(TOPO_QUAD);

    // Cycle through supported mesh element topologies
    for elemTopo in cellTopos
    {
      // Generate a random n-th degree test polynomial
      // Evaluate the test polynomial at the SPs
      // Interpolate the solution at the FPs
      // Evaluate the test polynomial at the FPs and compare the values
    }

    // Calculate the Lagrange Basis
    writeln("Interpolation Basis");
    for i in 0..9 {
      basis[i,..] = eval_LagrangePoly1D_array(x[i], node);
      writef("For x_%i = %7.4dr : %.4t\n", i, x[i], basis[i,..]);
    }
    writeln();

    for i in 0..9 {
      for j in 0..9 do
        y_node[j] += (coef[i]) * (node[j] ** i);
      for j in 0..9 {
        y_x[j] += (coef[i]) * (x[j] ** i);
        interpolation = 0.0;
        for k in 0..9 do
          interpolation += basis[j,k]*y_node[k];
        writef( "y_%i (%7.4dr): %7.4dr   Interpolation: %7.4dr   Error: %10.3er   Relative: %10.3er\n", i, x[j], y_x[j],
            interpolation, error(y_x(j), interpolation), relative_error(y_x(j), interpolation) );
      }
      writeln();
    }

    // Calculate the Lagrange Derivative Basis
    writeln("Interpolation Derivative Basis");
    for i in 0..9 {
      basisDeriv[i,..] = eval_DLagrangeDx_array(x[i], node);
      writef("For x_%i = %7.4dr : %.4t", i, x[i], basisDeriv[i,..]);
    }
    writeln();

    y_node = 0.0;
    y_x = 0.0;
    for i in 0..9 {
      for j in 0..9 do
        y_node[j] += (coef[i]) * (node[j] ** i);
      for j in 0..9 {
        y_x[j] += (i:real * coef[i]) * (x[j] ** (i-1));
        interpDeriv = 0.0;
        for k in 0..9 do
          interpDeriv += basisDeriv[j,k]*y_node[k];
        writef( "y_%i' (%7.4dr): %7.4dr   Interpolation Deriv: %7.4dr   Error: %10.3er   Relative: %10.3er\n", i, x[j], y_x[j],
            interpDeriv, error(y_x(j), interpDeriv), relative_error(y_x(j), interpDeriv) );
      }
      writeln();
    }

    var minOrder : int = 0;
    var maxOrder : int = 9;

    // Calculate the FR structures
    writeln();
    writeln("Interpolation initialized structure for FR (sp2fpInterp):");
    writeln();

    init_sp2fpInterp(minOrder, maxOrder, cellTopos);
    writeln(sp2fpInterp);
    writeln();

    writeln();
    writeln("Interpolation derivative initialized structure for FR (sp2spDeriv):");
    writeln();

    init_sp2spDeriv(minOrder, maxOrder, cellTopos);
    writeln(sp2spDeriv);
    writeln();

    writeln();
    writeln("Node interpolation initialized structure for FR (sp2nodeInterp):");
    writeln();

    init_sp2nodeInterp(minOrder, maxOrder, cellTopos);
    writeln(sp2nodeInterp);
    writeln();

    // Test the FR structures
    writeln();
    writeln("Testing the interpolation structure for FR(sp2fpInterp):");
    writeln();
    for interpOrder in minOrder..maxOrder
    {
      var fpCnt : int = 2;
      var spCnt : int = interpOrder+1;

      // Interpolation abscissa
      var spLoc : [1..spCnt] real = nodes_legendre_gauss(interpOrder+1);
      // Interpolation targets
      var fpLoc : [1..fpCnt] real = [-1.0, 1.0];

      // Generate random polynomial of degree interpOrder
      var coef : [0..interpOrder] real;
      randStreamSeeded.fillRandom(coef);
      coef = 2.0*coef -1.0;

      // Function evaluated at interpolation abscissa
      var yDirectSP : [1..spCnt] real = 0.0;
      for spIdx in spLoc.domain do
        for k in coef.domain do
          yDirectSP[spIdx] += coef[k] * spLoc[spIdx]**k;

      for fpIdx in fpLoc.domain
      {
        // Calculate values of this polynomial at the interpolation abscissa
        var yDirectFP : real = 0.0;
        var yInterpFP : real = 0.0;

        for k in coef.domain
        {
          yDirectFP += coef[k] * fpLoc[fpIdx]**k;
          yInterpFP += sp2fpInterp[(2, interpOrder)]!.coefs[fpIdx, k+1] * yDirectSP[k+1];
        }

        writef("Degree %2i interpolation | y(%7.4dr) = %7.4dr | interp(%7.4dr) = %7.4dr | Abs Erro = %10.3er | Rel Error = %10.3er\n",
            interpOrder, fpLoc[fpIdx], yDirectFP, fpLoc[fpIdx], yInterpFP,
            error(yDirectFP, yInterpFP), relative_error(yDirectFP, yInterpFP));
      }
      writeln();
    }

    writeln();
    writeln("Testing the interpolation derivative structure for FR(sp2spDeriv):");
    writeln();
    for interpOrder in minOrder..maxOrder
    {
      var spCnt : int = interpOrder+1;

      // Interpolation abscissa
      var spLoc : [1..spCnt] real = nodes_legendre_gauss(interpOrder+1);

      // Generate random polynomial of degree interpOrder
      var coef : [0..interpOrder] real;
      randStreamSeeded.fillRandom(coef);
      coef = 2.0*coef -1.0;

      // Function evaluated at interpolation abscissa
      var yDirectSP : [1..spCnt] real = 0.0;
      for spIdx in spLoc.domain do
        for k in coef.domain do
          yDirectSP[spIdx] += coef[k] * spLoc[spIdx]**k;

      for spIdx in spLoc.domain
      {
        // Calculate values of this polynomial at the interpolation abscissa
        var dyDirectSP : real = 0.0;
        var dyInterpSP : real = 0.0;

        for k in coef.domain
        {
          if k != 0 then
            dyDirectSP += k * coef[k] * spLoc[spIdx]**(k-1);
          dyInterpSP += sp2spDeriv[(2, interpOrder)]!.coefs[spIdx, 1, k+1] * yDirectSP[k+1];
        }

        writef("Degree %2i interpolation | y(%7.4dr) = %6.4dr | interp(%7.4dr) = %7.4dr | Abs Erro = %10.3er | Rel Error = %10.3er\n",
            interpOrder, spLoc[spIdx], dyDirectSP, spLoc[spIdx], dyInterpSP,
            error(dyDirectSP, dyInterpSP), relative_error(dyDirectSP, dyInterpSP));
      }
      writeln();
    }
  }
}
