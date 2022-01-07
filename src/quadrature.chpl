prototype module Quadrature
{
  use UnitTest;
  use Set;

  class quadrature_weights_c
  {
    var weights_d : domain(1); // {quadratureAbscissa}
    var weights: [weights_d] real;
  }

  // Define type for the interpolation structure. Add "?" to allow default initialization to nil.
  type quadrature_weights_t = unmanaged quadrature_weights_c?;

  // Domains
  var quadratureWeights_d : domain(2*int); // {elemTopo, quadratureDegree}

  // Coefficient structures
  var quadratureWeights : [quadratureWeights_d] quadrature_weights_t;

  //////////////////////////////////
  //   Initialization Procedure   //
  //////////////////////////////////

  proc init_quadratureWeights(minPolyDegree : int, maxPolyDegree : int, elemTopos : set(int))
  {
    use Set;
    use Parameters.ParamMesh;
    use Polynomials;

    // Add all combination of cell topology and interpolation order to the domain
    for elemTopo in elemTopos do
      for quadratureDegree in minPolyDegree..maxPolyDegree do
        quadratureWeights_d.add((elemTopo, quadratureDegree));

    for (elemTopo, quadratureDegree) in quadratureWeights.domain
    {
      select elemTopo
      {
        when TOPO_LINE
        {
          var nodeCnt : int = quadratureDegree+1;

          quadratureWeights[(elemTopo, quadratureDegree)] = new quadrature_weights_t({1..nodeCnt})!;

          quadratureWeights[(elemTopo, quadratureDegree)]!.weights = weights_legendre_gauss(nodeCnt);
        }
        when TOPO_TRIA {}
        when TOPO_QUAD {}
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at quadrature initialization.");
      }
    }
  }

  ///////////////////////////////
  //   Integration Procedure   //
  ///////////////////////////////

  proc integrate(nodalPoly : [] real, elemTopo : int, nodeCnt : int) : real
  {
    use LinearAlgebra;

    return dot(quadratureWeights[(elemTopo, nodeCnt)]!.weights, nodalPoly);
  }

  /////////////////////////////////////////
  //   Weight set finding functions      //
  /////////////////////////////////////////

  proc weights_legendre_gauss(n : int) : [1..n] real
  {
    /*
      COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
      N  = ORDER OF THE FORMULA
      CS = 0.0ES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
      DZ = VECTOR OF THE DERIVATIVES AT THE 0.0ES, DZ(I), I=1,N
      WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
    */

    use Polynomials;

    var xi : [1..n] real = nodes_legendre_gauss(n);
    var dy : [1..n] real;
    var weights : [1..n] real;

    for i in 1..n do
      dy[i] = eval_legendre_poly_dy(n, xi[i]);

    for i in 1..n do
      weights(i) = 2.0 / ((1.0 - xi(i)**2) * dy(i)**2);

    return weights;
  }

  proc weights_legendre_gauss_lobatto(n : int, vn : real, out weights : real)
  {
//    /*
//      COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
//      N  = ORDER OF THE FORMULA
//      ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
//      VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
//      WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
//    */
//
//    //.. Formal Arguments ..
//    int     ,                   intent(in)  :: n;
//    real, dimension(1:n+1), intent(in)  :: vn;
//    real, dimension(1:n+1), intent(out) :: weights;
//
//    //.. Local Scalars ..
//    int      :: i;
//    real :: c;
//
//    c = 2.0;
//    if (n > 0) then
//      c = c / (n*n+1) : real;
//
//    for i in 1..n+1 do
//      weights(i) = c / (vn(i)*vn(i));
  } // weights_legendre_gauss_lobatto

  proc weights_dgbook_jacobi_gauss_lobatto(n : int, alpha : real, beta : real, xi : real, out weights : real)
  {
//    /*
//      COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
//      N  = ORDER OF THE FORMULA
//      alpha  = PARAMETER > -1
//      beta  = PARAMETER > -1
//      ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
//      WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
//    */
//
//    //.. Formal Arguments ..
//    int,                        intent(in)  :: n;
//    real,                   intent(in)  :: alpha;
//    real,                   intent(in)  :: beta;
//    real, dimension(1:n+1), intent(in)  :: xi;
//    real, dimension(1:n+1), intent(out) :: weights;
//
//    //.. Local Scalars ..
//    int  :: i;
//    real :: a1, a2, b1, b2, ab, ab1, ab2;
//    real :: c, c1, c2, c3, ri, rn, su;
//    real :: y, z, dy, d2y, scl;
//
//    weights(:) = 0.0;
//
//    if (n > 1) {
//      a1  = alpha + 1.0;
//      a2  = alpha + 2.0;
//      b1  = beta  + 1.0;
//      b2  = beta  + 2.0;
//      ab  = alpha + beta;
//      ab1 = alpha + beta + 1.0;
//      ab2 = alpha + beta + 2.0;
//
//      c = (2.0**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2);
//
//      rn = real(n,kind=lp);
//      c  = c * (2.0*rn + ab) / (rn + ab1);
//      c1 = c * a1 / (b2 * ab2);
//      c2 = c * b1 / (a2 * ab2);
//      c3 = 0.5 * c * a1 * b1;
//
//      for i in 2..n-1 {
//        ri = i-1 : real;;
//        c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2));
//        c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2));
//        c3 = c3*(ri+a1)*(ri+b1) / ((ri+2.0)*(ri+ab1));
//      }
//
//      su = 0.0;
//      for i in 2..n do
//        su = su + xi(i);
//
//      scl = ( (1.0-xi(1))**a1 ) * ( (1.0+xi(1))**b1 );
//      weights(1)   = scl*c1*(rn-1.0-su);
//
//      scl = ( (1.0-xi(n+1))**a1 ) * ( (1.0+xi(n+1))**b1 );
//      weights(n+1) = scl*c2*(rn-1.0+su);
//
//      for i in 2..n {
//        eval_jacobi_poly(n  ,alpha,beta,xi(i),y,dy,d2y);
//        eval_jacobi_poly(n-1,alpha,beta,xi(i),z,dy,d2y);
//        scl = ( (1.0-xi(i))**a1 ) * ( (1.0+xi(i))**b1 );
//        weights(i) = -scl*c3/y/dy;
//      }
//    }
  }

  ///////////////////////////////
  //   Module Test Procedure   //
  ///////////////////////////////

  proc main()
  {
    /*
       Generate a random n-th degree polynomial and compare numerical quadrature to algebraic integration
    */

    use Random;
    use Set;
    use IO.FormattedIO;
    use Testing;
    use Parameters.ParamTest;
    use Parameters.ParamMesh;
    use Polynomials;
    use LinearAlgebra;

    var randStream = new RandomStream(real);
    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    var minQuadratureDegree : int = 0;
    var maxQuadratureDegree : int = 9;

    // Create a set with the cell topologies contained in the hypothetics test mesh
    var cellTopos : set(int);
    cellTopos.add(TOPO_LINE);  // Add Line element to the set

    // Intitialize projection matrix structure
    writeln();
    writeln("Quadrature initialized structure for FR (quadratureWeights)");
    writeln();
    init_quadratureWeights(minQuadratureDegree, maxQuadratureDegree, cellTopos);
    writeln(quadratureWeights);
    writef("\n");

    writef("Poly Degree | Algebraic Integral |  Simple Quadrature, Abs Error, Rel Error | Gauss Quadrature  , Abs Error, Rel Error\n");
    for testPolyDegree in minQuadratureDegree..maxQuadratureDegree
    {
      // Generate random polynomial coefficients
      var testPoly : [0..testPolyDegree] real;
      randStreamSeeded.fillRandom(testPoly);
      writef("%11i", testPolyDegree);

      // Algebraic Integration
      var algebraicIntegral : real = 0.0;
      for i in 0..testPolyDegree do
        algebraicIntegral += (testPoly[i]/(i+1)*1.0**(i+1)) - (testPoly[i]/(i+1)*(-1.0)**(i+1));

      writef(" | %{18.11er}", algebraicIntegral);

      // Simple Quadrature
      {
        var nodeCnt : int = testPolyDegree+1;
        var nodes         : [1..nodeCnt] real = nodes_legendre_gauss(nodeCnt);
        var nodalTestPoly : [1..nodeCnt] real = 0.0;
        for i in 1..nodeCnt do
          for j in 0..testPolyDegree do
            nodalTestPoly[i] += testPoly[j]*nodes[i]**j;
        writef(" | %{18.11er}, %{9.2er}, %{9.2er}", integrate(nodalTestPoly, TOPO_LINE, nodeCnt-1),
            error(algebraicIntegral, integrate(nodalTestPoly, TOPO_LINE, nodeCnt-1)),
            relative_error(algebraicIntegral, integrate(nodalTestPoly, TOPO_LINE, nodeCnt-1)));
      }

      // Gauss quadrature
      {
        var nodeCnt : int = testPolyDegree/2+1;
        var nodes         : [1..nodeCnt] real = nodes_legendre_gauss(nodeCnt);
        var nodalTestPoly : [1..nodeCnt] real = 0.0;
        for i in 1..nodeCnt do
          for j in 0..testPolyDegree do
            nodalTestPoly[i] += testPoly[j]*nodes[i]**j;
        writef(" | %{18.11er}, %{9.2er}, %{9.2er}", integrate(nodalTestPoly, TOPO_LINE, nodeCnt-1),
            error(algebraicIntegral, integrate(nodalTestPoly, TOPO_LINE, nodeCnt-1)),
            relative_error(algebraicIntegral, integrate(nodalTestPoly, TOPO_LINE, nodeCnt-1)));
      }
      writef("\n");
    }
  }
}
