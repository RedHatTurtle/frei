prototype module Interpolation
{
  use Random;
  use UnitTest;
  use Set;

  class interpolation_coefficients_c
  {
    var coefs_d : domain(2); // {nFPs, nSPs}
    var coefs: [coefs_d] real;
  }

  // Define type for the interpolation structure. Add "?" to allow default initialization to nil.
  type interpolation_coefficients_t = unmanaged interpolation_coefficients_c?;

  // Domains
  var sp2fpInterp_d : domain(2*int);
  var sp2spDeriv_d  : domain(2*int);

  // Coefficient structures
  var sp2fpInterp : [sp2fpInterp_d] interpolation_coefficients_t;
  var sp2spDeriv  : [sp2spDeriv_d]  interpolation_coefficients_t;

  proc init_sp2fpInterp(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Parameters.ParamMesh;
    use Polynomials;

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2fpInterp_d.add((cellTopo, interpOrder));

    for (cellTopo, interpOrder) in sp2fpInterp.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : range = 1..interpOrder;
          var fpCnt : range = 1..2;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spLoc : [spCnt] real = nodes_legendre_gauss(interpOrder);
          var fpLoc : [fpCnt] real = [-0.5, 0.5];

          sp2fpInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({fpCnt, spCnt})!;

          for fp in fpCnt do
            sp2fpInterp[(cellTopo, interpOrder)]!.coefs[{fp..fp,spCnt}] =
                  reshape(eval_LagrangePoly1D_array(fpLoc[fp], spLoc), {fp..fp, spCnt});
        }
        when TOPO_TRIA {}
        when TOPO_QUAD {}
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    // Calculate all relevant coefficients
  }

  proc init_sp2spDeriv(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Parameters.ParamMesh;
    use Polynomials;

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2spDeriv_d.add((cellTopo, interpOrder));

    for (cellTopo, interpOrder) in sp2spDeriv.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : range = 1..interpOrder;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spLoc : [spCnt] real = nodes_legendre_gauss(interpOrder);

          sp2spDeriv[(cellTopo, interpOrder)] = new interpolation_coefficients_t({spCnt, spCnt})!;

          for sp in spCnt do
            sp2spDeriv[(cellTopo, interpOrder)]!.coefs[{sp..sp,spCnt}] =
                  reshape(eval_DLagrangeDx_array(spLoc[sp], spLoc), {sp..sp, spCnt});
        }
        when TOPO_TRIA {}
        when TOPO_QUAD {}
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    // Calculate all relevant coefficients
  }

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

  proc main()
  {
    use IO;
    use Testing;
    use Parameters.ParamTest;
    use Polynomials;

    var randStream = new RandomStream(real);
    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    var interpolation      : real;
    var interpDeriv : real;
    var node          : [0..9] real;
    var x             : [0..9] real;
    var y_node        : [0..9] real = 0;
    var y_x           : [0..9] real = 0;
    var coef          : [0..9] real;
    var basis         : [0..9, 0..9] real;
    var basisDeriv    : [0..9, 0..9] real;
    var cellTopos     : set(int);

    writeln();

    // Get Chebyshev roots to use as interpolation nodes
    node = nodes_legendre_gauss(10);
    writeln("Interpolation nodes (Legendre-Gauss):         ", node);
    for i in 0..9 do
      node[i] = -cos( half_pi * (2*i+1)/10 )/2;
    writeln("Normalized interpolation nodes [-0.5 to 0.5]: ", node);
    writeln();

    // Get random interpolation targets [-0.5,0.5] and random polynomials
    x[0] = 0;
    x[9] = 1;
    randStreamSeeded.fillRandom(x[1..8]);
    randStreamSeeded.fillRandom(coef);
    x = x-0.5;
    coef = 2*coef -1;
    writeln("Interpolation Targets   = ", x);
    writeln("Polynomial Coefficients = ", coef);
    writeln();

    // Calculate the Lagrange Basis
    writeln("Interpolation Basis");
    for i in 0..9 {
      basis[i,..] = eval_LagrangePoly1D_array(x[i], node);
      writeln("For x_%i = %8.5dr : ".format(i, x[i]), basis[i,..]);
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
        writeln( "y_%i (%6.3dr): %7.4dr   Interpolation: %7.4dr   Error: %11.3er   Relative: %11.3er".format(i, x[j],
              y_x[j], interpolation, error(y_x(j), interpolation), relative_error(y_x(j), interpolation)) );
      }
      writeln();
    }

    // Calculate the Lagrange Derivative Basis
    writeln("Interpolation Derivative Basis");
    for i in 0..9 {
      basisDeriv[i,..] = eval_DLagrangeDx_array(x[i], node);
      writeln("For x_%i = %8.5dr : ".format(i, x[i]), basisDeriv[i,..]);
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
        writeln( "y_%i' (%6.3dr): %7.4dr   Interpolation Deriv: %7.4dr   Error: %11.3er   Relative: %11.3er".format(i,
              x[j], y_x[j], interpDeriv, error(y_x(j), interpDeriv), relative_error(y_x(j), interpDeriv)) );
      }
      writeln();
    }

    // Calculate the FR structures
    writeln();
    writeln("Interpolation initialized structure for FR:");
    writeln();

    cellTopos.add(2);

    init_sp2fpInterp(1, 9, cellTopos);
    writeln(sp2fpInterp);
    writeln();

    writeln();
    writeln("Interpolation derivative initialized structure for FR:");
    writeln();

    init_sp2spDeriv(1, 9, cellTopos);
    writeln(sp2spDeriv);
    writeln();
  }
}
