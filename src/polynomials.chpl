module Polynomials
{
  use Random;
  use UnitTest;

  /////////////////////////////////////////
  //   Polynomial evaluation functions   //
  /////////////////////////////////////////

  proc eval_jacobi_poly (n: int , alpha: real, beta: real, x: real, out y: real, out dy: real, out d2y: real)
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

  proc normalized_jacobi_poly(n : int, alpha : real, beta : real, x : real) : real
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

  proc eval_legendre_poly (n:int, x:real, out y:real, out dy:real, out d2y:real)
  {
    eval_jacobi_poly(n, 0, 0, x, y, dy, d2y);
  }

  proc eval_legendre_poly (n:int, x:real) : real
  {
    var y, dy, d2y : real;

    eval_jacobi_poly(n, 0, 0, x, y, dy, d2y);

    return y;
  }

  proc eval_legendre_poly_dy (n:int, x:real) : real
  {
    var y, dy, d2y : real;

    eval_jacobi_poly(n, 0, 0, x, y, dy, d2y);

    return dy;
  }

  proc eval_radau_poly (n:int, x:real, out y:real, out dy:real, out d2y:real)
  {
    var y_n, dy_n, d2y_n : real;
    var y_p, dy_p, d2y_p : real;

    eval_legendre_poly(  n, x, y_n, dy_n, d2y_n);
    eval_legendre_poly(n-1, x, y_p, dy_p, d2y_p);

    y   = 0.5*(   y_n -  y_p ) * (-1)**n;
    dy  = 0.5*(  dy_n - dy_p ) * (-1)**n;
    d2y = 0.5*( d2y_n -d2y_p ) * (-1)**n;
  }

  /////////////////////////////////////////
  //   Node/Root set finding functions   //
  /////////////////////////////////////////

  proc nodes_uniform(n : int) : [1..n] real
  {
    var xi : [1..n] real;

    for i in 1..n do
      xi[i] = 2*(i/(n+1):real)-1;

    return xi;
  }

  proc nodes_uniform(n : int, out xi : [1..n] real)
  {
     xi = nodes_uniform(n);
  }

  proc nodes_uniform_lobatto(n : int) : [1..n] real
  {
    var xi : [1..n] real;

    xi[1] = -1;
    xi[n] =  1;
    xi[2..n-1] = nodes_uniform(n-2);

    return xi;
  }

  proc nodes_uniform_lobatto(n : int, out xi : [1..n] real)
  {
    xi = nodes_uniform_lobatto(n);
  }

  proc nodes_jacobi_gauss(n : int, alpha : real, beta : real, out xi : [1..n] real, out dxi : [1..n] real)
  {
    use Parameters.ParamConstants;
    import Math.cos;

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

  proc nodes_legendre_gauss(n : int, out xi : [1..n] real, out dxi : [1..n] real)
  {
    nodes_jacobi_gauss(n, 0, 0, xi, dxi);
  }

  proc nodes_legendre_gauss(n  : int) : [1..n] real
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

  /////////////////////////////////////////
  //   Derivative evaluation functions   //
  /////////////////////////////////////////

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

  /////////////////////////////////////////
  //   Test main function                //
  /////////////////////////////////////////

  proc main()
  {
    use IO.FormattedIO;
    use Testing;

    var a, b : real;
    var x, y, dy, ddy : real;

    a = -0.5;
    b = -0.5;

    for order in 1..9 by 1 {
      writef( "Jacobi: n=%2i a=%5.2dr b=%5.2dr\n", order, a, b);
      for i in 0..order+2 by 1 {
        x = (i/(order+2):real)*2-1;
        eval_jacobi_poly(order, a, b, x, y, dy, ddy);
        writef("x=%7.4dr   y=%24.20dr   dy=%25.20dr   ddy=%26.20dr\n", x, y, dy, ddy);
      }
      writeln();

      writef( "Legendre: n=%2i\n", order);
      for i in 0..order+2 by 1 {
        x = (i/(order+2):real)*2-1;
        eval_legendre_poly(order, x, y, dy, ddy);
        writef("x=%7.4dr   y=%24.20dr   dy=%25.20dr   ddy=%26.20dr\n", x, y, dy, ddy);
      }
      writeln();

      writef( "Radau: n=%2i\n", order);
      for i in 0..order+2 by 1 {
        x = (i/(order+2):real)*2-1;
        eval_radau_poly(order, x, y, dy, ddy);
        writef("x=%7.4dr   y=%24.20dr   dy=%25.20dr   ddy=%26.20dr\n", x, y, dy, ddy);
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
