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

  proc correction_ga (in n : int, in x : real) : real
  {
    // 
  }

  proc main()
  {
    use Testing;

    var x : real;

    for i in 0..9 {
      for j in 0..i {
        x = j:real/i:real - 0.5;
        writeln( "n: %2i,   x = %7.4dr,   g_dg = %8.5dr,   g_2 = %8.5dr".format(i, x, correction_dg(i, x), correction_g2(i, x)) );
      }
      writeln();
    }
  }
}
