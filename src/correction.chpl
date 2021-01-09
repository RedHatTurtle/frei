prototype module Correction
{
  use Random;
  use UnitTest;

  // In GFR the CorrectionMatrix_Edge procedure returns the value of the selected

  proc correction_dg (in n : int, in x : real) : real
  {
    // g_{DG}(n,x) = R_{R,n}(x) = \frac{(-1)^n}{2}(P_n-P_{n-1})

    use Polynomials;

    var y, dy, ddy : real;

    eval_radau_poly(n, x, y, dy, ddy);

    return y;
  }

  proc correction_ga (in n : int, in x : real) : real
  {
    // g_ga(n,x) defined by roots and normalized with g_ga(-1)=1

    use Polynomials;

    var xi, tmp : [1..n] real;
    var y_norm, y_x : real;
    var y, dy, ddy : real;

    eval_legendre_poly(n-1, x, y, dy, ddy);
    y = ((-1)**n)*y*(x-1)/2;

    return y;
  }

  proc correction_g2 (in n : int, in x : real) : real
  {
    // g_2(n,x) = \frac{n-1}{2n-1}*R_{R,n} + \frac{n}{2n-1}*R_{R,n-1}

    use Polynomials;

    var y_n, dy_n, ddy_n : real;
    var y_p, dy_p, ddy_p : real;

    eval_radau_poly(  n, x, y_n, dy_n, ddy_n);
    eval_radau_poly(n-1, x, y_p, dy_p, ddy_p);

    return y_n*( (n-1):real/(2*n-1):real ) + y_p*( n:real/(2*n-1):real );
  }

  proc main()
  {
    use Testing;
    use Polynomials;

    for i in 1..9{
      var x : [0..i+1] real;
      x[0] = -1;
      x[i+1] = 1;
      x[1..i] = nodes_legendre_gauss(i);
      for j in 0..i+1 {
        writeln( "n: %2i,   x = %7.4dr,   g_dg = %8.5dr,   g_ga = %8.5dr,   g_2 = %8.5dr"
            .format(i, x[j], correction_dg(i, x[j]), correction_ga(i, x[j]), correction_g2(i, x[j])) );
      }
      writeln();
    }
  }
}
