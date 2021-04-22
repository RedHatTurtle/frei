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

  // Perhaps it might be useful in the future to have a sparse domain with the cell topologies present in the mesh
  var sp2fpInterp_d : domain(2*int);
  //var sp2fpInterp_d : domain(2); // {elemType, interpOrder}

  var sp2fpInterp : [sp2fpInterp_d] interpolation_coefficients_t;

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

  proc eval_LagrangePoly1D( in k : int, in x : real, in xi : [] real ) : real
  {
    // Get the value of the 1D Lagrange polynomial at point x
    //
    //                  nfp
    //   d             _____
    //  -- [L_k (x)] =  | |   x   - xi(i)
    //  dx              | |  -------------
    //                  i=1  xi(k) - xi(i)
    //                 i/=k
    //
    // K  : The value of the Lagrange polynomial is 1 at xi(k)
    // X  : Location at which to evaluate the k-th Lagrange polynomial
    // XI : Array of all the interpolation points
    //
    // This evaluates the Lagrange polynomial at the location x.
    // The Lagrange polynomial corresponds to a function that
    // is 0.0 at all interpolation points xi(i) for i=1,nfp
    // except at xi(k) where the function is equal to 1.

    var return_value : real = 1;
    var factors : domain(int) = xi.domain;

    for i in (factors - k) do
      return_value *= (x-xi[i]) / (xi[k]-xi[i]);

    return return_value;

  } // eval_LagrangePoly1D

  proc eval_LagrangePoly1D_array( in x : real, in xi : [] real) : [] real
  {
    // Get the values of the 1D Lagrange polynomials at point x
    //
    // X  : Location at which to evaluate the k-th Lagrange polynomial
    // XI : Array of all the interpolation points
    //
    // This evaluates the Lagrange polynomial at the location x.
    // The Lagrange polynomial corresponds to a function that
    // is 0.0 at all iterpolation points xi(i) for i=1,nfp
    // except at xi(k) where the function is equal to 1.

    var return_value : [xi.domain] real;

    for k in xi.domain do
      return_value[k] = eval_LagrangePoly1D(k, x, xi);

    return return_value;

  } // eval_LagrangePoly1D_array

  proc eval_DLagrangeDx( in k : int, in x : real, in xi : [] real) : real
  {
    // Get the value of the derivative of the
    // k-th Lagrange polynomial at the point x.
    //
    //                  nfp                nfp
    //   d              ___               _____
    //  -- [L_k (x)] =  \         1        | |   x   - x(j)
    //  dx              /__  -----------   | |  -----------
    //                       x(k) - x(i)   j=1  x(k) - x(j)
    //                  i=1               j/=i
    //                 i/=k               j/=k
    //
    // k  : index where the value at xi(k) is 1
    // X  : x location at which to evaluate the derivative
    //      of the k-th Lagrange polynomial
    // XI : array of all the interpolation points

    var return_value : real = 0.0;

    if (xi.size == 1) then
      return 0;

    for i in xi.domain - k do
      return_value += eval_LagrangePoly1D(i,x,xi[xi.domain-i]) / (xi(k)-xi(i));

    return return_value;

  } // eval_DLagrangeDx

  proc eval_D2LagrangeDx2(m : int, x : real, xi : real) : real
  {
    // Get value of second derivative of
    // m-th Lagrange polynomial at point x.
    //
    //                   nfp                nfp                nfp
    //  d2               ___                ___               _____
    // --- [L_m (x)]  =  \         1        \         1        | |   x   - x(j)
    // dx2               /__  -----------   /__  -----------   | |  -----------
    //                        x(m) - x(l)        x(m) - x(n)   j=1  x(m) - x(j)
    //                   l=1                n=1               j/=l
    //                  l/=m               n/=l               j/=n
    //                                     n/=m               j/=m
    //
    // m  : index where the value at xi(m) is 1
    // X  : x location at which to evaluate the second
    //      derivative of the m-th Lagrange polynomial
    // XI : array of all the interpolation points

//    //.. Local Scalars ..
//    var l, n, nfp : int;
//
//    //.. Local Arrays ..
//    int, dimension(1:size(xi)) :: imask
//    logical(lk), dimension(1:size(xi)) :: lmask
//
//    nfp = size(xi)
//
//    imask(:) = (/ (l,l=1,nfp) /)
//
//    return_value = 0.0
//
//    if (nfp > 2) {
//
//      for l in [1..nfp] {
//        if (l == m) cycle
//        for n in [1..nfp] {
//          if ( (n == m) .or. (n == l) ) cycle
//          lmask = (imask/=m .and. imask/=l .and. imask/=n)
//          return_value = return_value + product( x    -xi(1:nfp) , lmask ) / &
//                                        product( xi(m)-xi(1:nfp) , lmask ) / &
//                                        (xi(m) - xi(n))
//        }
//        return_value = return_value / (xi(m) - xi(l))
//      }
//    }
  }

  proc main()
  {
    use IO;
    use Testing;
    use Parameters.ParamTest;
    use Polynomials;

    var randStream = new RandomStream(real);
    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    var node : [0..9] real;
    var x : [0..9] real;
    var y_node : [0..9] real = 0;
    var y_x : [0..9] real = 0;
    var coef : [0..9] real;
    var basis : [0..9, 0..9] real;
    var interpolation : real;
    var cellTopos : set(int);

    writeln();

    // Get Chebyshev roots to use as interpolation nodes
    node = nodes_legendre_gauss(10);
    writeln("Interpolation nodes: ", node);
    for i in 0..9 do
      node[i] = -cos( half_pi * (2*i+1)/10 )/2;
    writeln("Normalized interpolation nodes: ", node);
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
        interpolation = 0;
        for k in 0..9 do
          interpolation += basis[j,k]*y_node[k];
        writeln( "y_%i (%6.3dr): %7.4dr   Interpolation: %7.4dr   Error: %11.3er   Relative: %11.3er".format(i, x[j], y_x[j], interpolation, error(y_x(j), interpolation), relative_error(y_x(j), interpolation)) );
      }
      writeln();
    }

    writeln();
    writeln("Interpolation initialized structure for FR:");
    writeln();

    cellTopos.add(2);

    init_sp2fpInterp(1,9,cellTopos);
    writeln(sp2fpInterp);
    writeln();

    //writeln( "y_%i (%6.3dr): %7.4dr   Interpolation: %7.4dr   Error: %11.3er   Relative: %11.3er".format(i, x[j], y_x[j], interpolation, error(y_x(j), interpolation), relative_error(y_x(j), interpolation)) );
    writeln();
  }
}
