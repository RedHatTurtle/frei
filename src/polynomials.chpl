prototype module Polynomials
{
  use Random;
  use UnitTest;

  proc eval_jacobi_poly (in  n: int , in  alpha: real, in beta: real, in x: real,
                         out y: real, out dy: real, out d2y: real)
  {

    var yp, dyp, d2yp : real;
    var ym, dym, d2ym : real;
    var ab, di : real;
    var c0, c1, c2, c3, c4 : real;

    // Initialize 0-th order Jacobi polynomial y=1
    y   = 1.0;
    dy  = 0.0;
    d2y = 0.0;

    // If the polynomial order is greater than 0 initialize the 1st degree Jacobi polynomial
    if (n > 0) {
      ab  = alpha + beta;

      y   = 0.5*(ab + 2.0)*x + 0.5*(alpha - beta);
      dy  = 0.5*(ab + 2.0);
      d2y = 0.0;
    }

    // Use the recurrence relation to get the higher order Jacobi polynomials
    if (n > 1) {

      // Set previous as 0-th order Jacobi polynomial
      yp   = 1.0;
      dyp  = 0.0;
      d2yp = 0.0;

      // Iterate through the recurrence relation
      for i in 2..n do {
        di = i : real;
        c0 = 2.0*di + ab;
        c1 = 2.0 * di * (di+ab) * (c0-2.0);
        c2 = (c0-1.0) * (c0-2.0) * c0;
        c3 = (c0-1.0) * (alpha-beta) * ab;
        c4 = 2.0 * c0 * (di+alpha-1.0) * (di+beta-1.0);

        ym = y;
        y  = ((c2*x + c3)*y - c4*yp)/c1;
        yp = ym;

        dym = dy;
        dy  = ((c2*x + c3)*dy - c4*dyp + c2*yp)/c1;
        dyp = dym;

        d2ym = d2y;
        d2y  = ((c2*x + c3)*d2y - c4*d2yp + 2.0*c2*dyp)/c1;
        d2yp = d2ym;
      }
    }
  }

  proc normalized_jacobi_poly(in n : int, in alpha : real, in beta : real, in x : real) : real
  {
//    //.. Local Scalars ..
//    var i : int;
//    var ab, ab1, ab2, ab3, amb, a1, b1 : real;
//    var gamma0, gamma1 : real;
//    var aold, anew, bnew, h1, ri : real;
//    var sqrt_numer, sqrt_denom : real;
//
//    //.. Local Arrays ..
//    var pl : [1..n+1] real;
//
//    ab = alpha + beta;
//    amb = alpha - beta;
//    a1 = alpha + 1.0;
//    b1 = beta + 1.0;
//    ab1 = ab + 1.0;
//    ab2 = ab + 2.0;
//    ab3 = ab + 3.0;
//
//    gamma0 = 2.0**ab1 / ab1 * gammafun(a1) * gammafun(b1) / gammafun(ab1);
//
//    pl(1) = 1.0/sqrt(gamma0);
//
//    if (n == 0) {
//      return pl(1);
//    } else {
//
//      gamma1 = a1 * b1 / ab3 * gamma0;
//
//      pl(2) = 0.5*( ab2*x + amb ) / sqrt(gamma1);
//
//      aold = 2.0 / ab2 * sqrt( a1*b1/ab3 );
//
//      for i in 1..n-1 {
//        ri = real(i,kind=lp);
//        h1 = 2.0*ri + ab;
//        sqrt_numer = (ri+1.0)*(ri+ab1)*(ri+a1)*(ri+b1);
//        sqrt_denom = (h1+1.0)*(h1+3.0);
//        anew = sqrt( sqrt_numer / sqrt_denom ) * 2.0 / (h1+2.0);
//        bnew = (beta*beta - alpha*alpha) / h1 / (h1+2.0);
//        pl(i+2) = 1.0/anew*( -aold*pl(i) + (x-bnew)*pl(i+1) );
//        aold = anew;
//      }
//
//      return pl(n+1);
//    }
  } // normalized_jacobi_poly

  proc eval_legendre_poly (in n:int, in x:real,
                           out y:real, out dy:real, out d2y:real)
  {
    eval_jacobi_poly(n, 0, 0, x, y, dy, d2y);
  }

  proc eval_radau_poly (in n:int, in x:real,
                           out y:real, out dy:real, out d2y:real)
  {
    var y_n, dy_n, d2y_n : real;
    var y_p, dy_p, d2y_p : real;

    eval_legendre_poly(  n, x, y_n, dy_n, d2y_n);
    eval_legendre_poly(n-1, x, y_p, dy_p, d2y_p);

    y   = 0.5*(   y_n -  y_p ) * (-1)**n;
    dy  = 0.5*(  dy_n - dy_p ) * (-1)**n;
    d2y = 0.5*( d2y_n -d2y_p ) * (-1)**n;
  }

  proc grad_normalized_jacobi_poly(n : int, alpha : real, beta : real, x : real) : real
  {
//    //.. Local Scalars ..
//    var rn, a1, b1, ab1 : real;
//
//    if (n == 0) {
//      return 0.0;
//    } else {
//      rn  = n : real;
//      a1  = alpha + 1.0;
//      b1  = beta  + 1.0;
//      ab1 = alpha + b1;
//
//      return sqrt(rn*(rn+ab1)) * normalized_jacobi_poly(n-1,a1,b1,x);
//    }
  }

  proc nodes_uniform(in n : int) : [1..n] real
  {
    var xi : [1..n] real;

    for i in 1..n do
      xi[i] = 2*(i/(n+1):real)-1;

    return xi;
  }

  proc nodes_uniform(in   n : int,
                     out xi : [1..n] real)
  {
     xi = nodes_uniform(n);
  }

  proc nodes_uniform_lobatto(in n : int) : [1..n] real
  {
    var xi : [1..n] real;

    xi[1] = -1;
    xi[n] =  1;
    xi[2..n-1] = nodes_uniform(n-2);

    return xi;
  }

  proc nodes_uniform_lobatto(in   n : int,
                             out xi : [1..n] real)
  {
    xi = nodes_uniform_lobatto(n);
  }

  proc nodes_jacobi_gauss(in n : int, in alpha : real, in beta : real,
                          out xi : [1..n] real, out dxi : [1..n] real)
  {
    use Parameters.ParamConstants;

    //.. Local Scalars ..
    var i, it : int;
    var konst, x, y, dy, d2y : real;

    if (n > 0)
    {
      xi(1)  = (beta - alpha) / (alpha + beta + 2.0);
      dxi(1) = 0.5*(alpha + beta) + 1.0;
    }

    if (n > 1)
    {
      konst = 0.5*PI / (2.0*n + alpha + beta + 1.0);

      for i in 1..n
      {
        x  = -cos( konst * (4.0*i + alpha + beta - 1.0) );

        for it in 1..16
        {
          eval_jacobi_poly(n, alpha, beta, x, y, dy, d2y);

          if (abs(y) <= EPS16) then
            break;

          x = x - y/dy;
        }

        if (abs(x) <= EPS17) then
          x = 0.0;

        xi(i)  = x;
        dxi(i) = dy;
      }

      for i in 1..n/2 {
        xi (n-i+1) = -xi (i);
       //dxi(n-i+1) = -dxi(i);
      }
    }
  }

  proc nodes_legendre_gauss(in  n  : int,
                            out xi : [1..n] real, out dxi : [1..n] real)
  {
    nodes_jacobi_gauss(n, 0, 0, xi, dxi);
  }

  proc nodes_legendre_gauss(in  n  : int) : [1..n] real
  {
    var xi, dxi : [1..n] real;
    nodes_jacobi_gauss(n, 0, 0, xi, dxi);
    return xi;
  }

  proc nodes_legendre_gauss_lobatto(n,xi,vn)
  {
//    //.. Formal Arguments ..
//    int     ,                   intent(in)  :: n;
//    real, dimension(1:n+1), intent(out) :: xi;
//    real, dimension(1:n+1), intent(out) :: vn;
//
//    //.. Local Scalars ..
//    int      :: i, it, n2;
//    real :: sn, x, y, dy, d2y;
//
//    //.. Local Constants ..
//    real, parameter :: pi   = 3.14159265358979323846;
//
//    if (n == 0) {
//      xi(1) = 0.0;
//      vn(1) = 1.0;
//    } else {
//      n2 = (n-1)/2;
//      sn = real(2*n - 4*n2 - 3,kind=lp);
//
//      xi(1) = -1.0;
//      vn(1) =  sn;
//
//      xi(n+1) = 1.0;
//      vn(n+1) = 1.0;
//    }
//
//    if (n > 1) {
//      x = 0.0;
//      eval_legendre_poly(n,x,y,dy,d2y);
//      xi(n2+2) = x;
//      vn(n2+2) = y;
//    }
//
//    if (n > 2) {
//      for i in 1..n2 {
//        x = cos( pi * real(i,kind=lp) / real(n,kind=lp) )
//
//        for it in 1..16 {
//          eval_legendre_poly(n,x,y,dy,d2y)
//          x = x - dy/d2y
//        }
//
//        xi(i+1)   = -x
//        xi(n-i+1) =  x
//
//        vn(i+1)   = y*sn
//        vn(n-i+1) = y
//      }
//    }
  }

  proc weights_legendre_gauss(n,xi,dy,weights)
  {
//    //****************************************************************
//    //   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
//    //   N  = ORDER OF THE FORMULA
//    //   CS = 0.0ES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
//    //   DZ = VECTOR OF THE DERIVATIVES AT THE 0.0ES, DZ(I), I=1,N
//    //   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
//    //****************************************************************
//
//    //.. Formal Arguments ..
//    int     ,                 intent(in)  :: n;
//    real, dimension(1:n), intent(in)  :: xi;
//    real, dimension(1:n), intent(in)  :: dy;
//    real, dimension(1:n), intent(out) :: weights;
//
//    //.. Local Scalars ..
//    int :: i;
//
//    weights(:) = 2.0;
//
//    for i in 1..n do
//      weights(i) = 2.0 / ((1.0 - xi(i)*xi(i)) * dy(i) * dy(i));
  }

  proc weights_legendre_gauss_lobatto(in n : int, in vn : real,
                                      out weights : real)
  {
//    //**********************************************************************
//    //   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
//    //   N  = ORDER OF THE FORMULA
//    //   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
//    //   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
//    //   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
//    //**********************************************************************
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

  proc weights_dgbook_jacobi_gauss_lobatto(in n : int, in alpha : real, in beta : real, in xi : real,
                                           out weights : real)
  {
//    //********************************************************************
//    //   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
//    //   N  = ORDER OF THE FORMULA
//    //   alpha  = PARAMETER > -1
//    //   beta  = PARAMETER > -1
//    //   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
//    //   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
//    //********************************************************************
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

  proc derivative_matrix_legendre_gauss(n,a)
  {
//    //***********************************************************************
//    //  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
//    //  LEGENDRE GAUSS NODES
//    //***********************************************************************
//
//    //.. Formal Arguments ..
//    int     , intent(in) :: n; // ORDER OF LEGENDRE GAUSS POLYNOMIAL
//    real, intent(in) :: a; // REAL VALUE TO DISTINGUISH
//
//    //.. Function Result ..
//    real, dimension(1:n+1,1:n+1) :: return_value; // DERIVATIVE MATRIX
//
//    //.. Local Scalars ..
//    int  :: i,j,np;
//    real :: vi,zi,vj,zj;
//
//    //.. Local Arrays ..
//    real, dimension(1:n+1) :: cs; // VECTOR OF THE 0.0ES
//    real, dimension(1:n+1) :: dz; // LEGENDRE DERIVATIVES AT THE 0.0ES
//
//    np = n + 1;
//
//    return_value(1,1) = 0.0;
//
//    if (np <= 1) then
//      return;
//
//    // Compute the 0.0s of the Legendre Gauss polynomial of
//    // order 'n' and the value of the derivative at each 0.0
//
//    nodes_legendre_gauss(np,cs,dz);
//
//    for i in 1..np {
//      vi = dz(i);
//      zi = cs(i);
//      for j in 1..np {
//        if (i /= j) {
//          vj = dz(j);
//          zj = cs(j);
//          return_value(i,j) = vi/(vj*(zi-zj));
//        } else {
//          return_value(i,i) = zi/(1.0-zi*zi);
//        }
//      }
//    }
  }

  proc derivative_matrix_legendre_gauss_lobatto(n,a)
  {
//    //***********************************************************************
//    //  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
//    //  LEGENDRE GAUSS-LOBATTO NODES
//    //***********************************************************************
//
//    //.. Formal Arguments ..
//    int,  intent(in) :: n; // ORDER OF LEGENDRE GAUSS-LOBATTO POLYNOMIAL
//    real, intent(in) :: a; // REAL VALUE TO DISTINGUISH BETWEEN DP AND QP VERSIONS
//
//    //.. Function Result ..
//    real, dimension(1:n+1,1:n+1) :: return_value; // DERIVATIVE MATRIX
//
//    //.. Local Scalars ..
//    int  :: i,j,np;
//    real :: ei,ej,vi,vj,dn,c;
//
//    //.. Local Arrays ..
//    real, dimension(1:n+1) :: et; // VECTOR OF THE NODES
//    real, dimension(1:n+1) :: vn; // VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES
//
//    np = n + 1;
//
//    // Compute the 0.0s of the Legendre Gauss Lobatto polynomial of
//    // order 'n' and the value of the Legendre Gauss polynomial of the
//    // same order at each 0.0
//
//    nodes_legendre_gauss_lobatto(n,et,vn);
//
//    for i in 1..np {
//      vi = vn(i);
//      ei = et(i);
//      for j in 1..np {
//        if (i /= j) {
//          vj = vn(j);
//          ej = et(j);
//          return_value(i,j) = vi/(vj*(ei-ej));
//        } else {
//          return_value(i,i) = 0.0;
//        }
//      }
//    }
//
//    if (np > 1) {
//      dn = real(np-1,kind=lp);
//      c  = 0.25 * dn * (dn+1.0);
//      return_value(1,1) = -c;
//      return_value(np,np) = c;
//    }
  }

  proc main()
  {
    use IO.FormattedIO;
    use Testing;

    var a, b : real;
    var x, y, dy, ddy : real;

    a = -0.5;
    b = -0.5;

    for order in 1..9 by 1 {
      writeln( "Jacobi: n=%3i a=%5.2dr b=%5.2dr".format(order, a, b));
      for i in 0..order+2 by 1 {
        x = (i/(order+2):real)*2-1;
        eval_jacobi_poly(order, a, b, x, y, dy, ddy);
        writeln( "x=%7.4dr   y=%24.20dr   dy=%25.20dr   ddy=%26.20dr".format(x, y, dy, ddy) );
      }
      writeln();

      writeln( "Legendre: n=%3i".format(order));
      for i in 0..order+2 by 1 {
        x = (i/(order+2):real)*2-1;
        eval_legendre_poly(order, x, y, dy, ddy);
        writeln( "x=%7.4dr   y=%24.20dr   dy=%25.20dr   ddy=%26.20dr".format(x, y, dy, ddy) );
      }
      writeln();

      writeln( "Radau: n=%3i".format(order));
      for i in 0..order+2 by 1 {
        x = (i/(order+2):real)*2-1;
        eval_radau_poly(order, x, y, dy, ddy);
        writeln( "x=%7.4dr   y=%24.20dr   dy=%25.20dr   ddy=%26.20dr".format(x, y, dy, ddy) );
      }
      writeln();

      var xi, dxi : [1..order] real;
      nodes_uniform(order, xi);
      writeln("Uniformly spaced Nodes: xi=[", xi, "]");
      writeln();
      nodes_uniform_lobatto(order, xi);
      writeln("Uniformly spaced Nodes with edges: xi=[", xi, "]");
      writeln();
      nodes_jacobi_gauss(order, a, b, xi, dxi);
      writeln("Jacobi-Gauss Nodes: xi=[", xi, "]");
      writeln("    dxi=[", dxi, "]" );
      writeln();
      nodes_legendre_gauss(order, xi, dxi);
      writeln("Legendre-Gauss Nodes: xi=[", xi, "]");
      writeln("    dxi=[", dxi, "]" );
      writeln();
      //nodes_legendre_gauss_lobatto(order, a, b, xi, dxi);
      //writeln("Legendre-Gauss-Lobatto Nodes: xi=[", xi, "]");
      //writeln("    dxi=[", dxi, "]" );
      //writeln();

      writeln();
      writeln();
      writeln();
    }
  }

}
