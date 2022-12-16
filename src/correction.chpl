module Correction
{
  // In GFR there is correction structure initialization procedure that then calls individual procedures for each cell
  //    geometry. Ex: CorrectionMatrix_Edge, CorrectionMatrix_Quad, CorrectionMatrix_Tri.
  //    The edge procedure selects the correction function depending on the input file and the other procedures are
  //    based on the values returned by this one. Ex: Using tensor product of the edge correction function to build the
  //    quadrilateral and hexahedral correction functions.
  //
  // The GFR structure for storing these correction function value has the following format:
  //    correct(cell_geom,cell_order)%mat(sp_id,fp_id)

  use Random;
  use UnitTest;
  use Set;

  class flux_correction_c
  {
    var correction_d : domain(2); // {nFPs, nSPs}
    var correction : [correction_d] real;
  }

  // Define type for the flux correction structure. Add "?" to allow default initialization to nil.
  type flux_correction_t = unmanaged flux_correction_c?;

  // Perhaps it might be useful in the future to have a sparse domain with the cell topologies present in the mesh
  //var flux_correction_d : sparse domain(2);

  var flux_correction_d : domain(2*int); // {cellType, interpOrder}

  var flux_correction : [flux_correction_d] flux_correction_t;

  //////////////////////////////////
  //   Initialization Procedure   //
  //////////////////////////////////

  proc init_correction(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Time;
    import Input.frScheme;
    use Parameters.ParamInput;
    use Parameters.ParamMesh;
    use Polynomials;

    writeln();
    writeln("Initializing Correction matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var stopwatch : Timer;
    stopwatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        flux_correction_d.add((cellTopo, interpOrder));

    for (cellTopo, interpOrder) in flux_correction.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = interpOrder;
          var fpCnt : int = 2;

          flux_correction[(cellTopo, interpOrder)] = new flux_correction_t({1..fpCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh uses has SPs on Legendre roots. xD
          var spDistLine : [1..interpOrder] real = nodes_legendre_gauss(interpOrder);

          var correctionCoefsLine : [1..interpOrder] real;
          select frScheme
          {
            when FR_DG do
              correctionCoefsLine = correction_dg_deriv(interpOrder, spDistLine);
            when FR_GA do
              correctionCoefsLine = correction_ga_deriv(interpOrder, spDistLine);
            when FR_G2 do
              correctionCoefsLine = correction_g2_deriv(interpOrder, spDistLine);
          }

          for (fpIdx, spIdx) in {1..1, 1..spCnt} do
            flux_correction[(cellTopo, interpOrder)]!.correction[fpIdx, spIdx] = correctionCoefsLine[spIdx];

          // Invert the correction derivative for right side FP
          for (fpIdx, spIdx) in {2..2, 1..spCnt} do
            flux_correction[(cellTopo, interpOrder)]!.correction[fpIdx, spIdx] = -1.0*correctionCoefsLine[spCnt - (spIdx-1)];

        }
        when TOPO_TRIA {}
        when TOPO_QUAD
        {
          var spCnt : int = (interpOrder)**2;
          var fpCnt : int = (interpOrder)*4;

          flux_correction[(cellTopo, interpOrder)] = new flux_correction_t({1..fpCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh uses has SPs on Legendre roots.
          var spDistLine : [1..interpOrder] real = nodes_legendre_gauss(interpOrder);

          var correctionCoefsLine : [1..interpOrder] real;
          select frScheme
          {
            when FR_DG do
              correctionCoefsLine = correction_dg_deriv(interpOrder, spDistLine);
            when FR_GA do
              correctionCoefsLine = correction_ga_deriv(interpOrder, spDistLine);
            when FR_G2 do
              correctionCoefsLine = correction_g2_deriv(interpOrder, spDistLine);
          }

          for (cellFace, faceFP, lineSP) in {1..4, 1..interpOrder, 1..interpOrder}
          {
            var cellFP : int = faceFP+(cellFace-1)*(interpOrder);

            select cellFace
            {
              when 1
              {
                var cellSP : int = (lineSP-1)*(interpOrder) + faceFP;
                flux_correction[(cellTopo, interpOrder)]!.correction[cellFP,cellSP] = +correctionCoefsLine[lineSP];
              }
              when 2
              {
                var cellSP : int = faceFP*(interpOrder) - (lineSP-1);
                flux_correction[(cellTopo, interpOrder)]!.correction[cellFP,cellSP] = -correctionCoefsLine[lineSP];
              }
              when 3
              {
                var cellSP : int = (spCnt+1)-faceFP - (lineSP-1)*(interpOrder);
                flux_correction[(cellTopo, interpOrder)]!.correction[cellFP,cellSP] = -correctionCoefsLine[lineSP];
              }
              when 4
              {
                var cellSP : int = (spCnt+1)-faceFP*(interpOrder) + (lineSP-1);
                flux_correction[(cellTopo, interpOrder)]!.correction[cellFP,cellSP] = +correctionCoefsLine[lineSP];
              }
            }
          }
        }
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at flux correction initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", stopwatch.elapsed(TimeUnits.milliseconds));
  }

  //////////////////////////////
  //   Correction Functions   //
  //////////////////////////////

  proc correction_dg (in n : int, in x : real) : real
  {
    // g_{DG}(n,x) = R_{R,n}(x) = \frac{(-1)^n}{2}(P_n-P_{n-1})

    use Polynomials;

    var y, dy, ddy : real;

    eval_radau_poly(n, x, y, dy, ddy);

    return y;
  }

  proc correction_dg_deriv (in n : int, in x : real) : real
  {
    // g_{DG}(n,x) = R_{R,n}(x) = \frac{(-1)^n}{2}(P_n-P_{n-1})

    use Polynomials;

    var y, dy, ddy : real;

    eval_radau_poly(n, x, y, dy, ddy);

    return dy;
  }

  proc correction_ga (in n : int, in x : real) : real
  {
    // g_ga(n,x) defined by roots and normalized with g_ga(-1)=1

    use Polynomials;

    var y, dy, ddy : real;

    eval_legendre_poly(n-1, x, y, dy, ddy);
    y = ((-1)**n)*y*(x-1)/2;

    return y;
  }

  proc correction_ga_deriv (in n : int, in x : real) : real
  {
    // g_ga(n,x) defined by roots and normalized with g_ga(-1)=1

    use Polynomials;

    var y, dy, ddy : real;

    eval_legendre_poly(n-1, x, y, dy, ddy);

    // This comes from applyng the product rule to the expression on correction_ga
    dy = ((-1)**n)*(dy*(x-1.0)/2.0 + y/2.0);

    return dy;
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

  proc correction_g2_deriv (in n : int, in x : real) : real
  {
    // g_2(n,x) = \frac{n-1}{2n-1}*R_{R,n} + \frac{n}{2n-1}*R_{R,n-1}

    use Polynomials;

    var y_n, dy_n, ddy_n : real;
    var y_p, dy_p, ddy_p : real;

    eval_radau_poly(  n, x, y_n, dy_n, ddy_n);
    eval_radau_poly(n-1, x, y_p, dy_p, ddy_p);

    return dy_n*( (n-1):real/(2*n-1):real ) + dy_p*( n:real/(2*n-1):real );
  }

  ///////////////////////////////
  //   Module Test Procedure   //
  ///////////////////////////////

  proc main()
  {
    use IO;
    use Parameters.ParamMesh;
    use Testing;
    use Polynomials;

    for i in 1..9{
      var x : [0..i+2] real;
      x[0] = -1;
      x[i+2] = 1;
      x[1..i+1] = nodes_uniform(i+1);
      for j in 0..i+2 {
        writef( "n: %2i, x = %7.4dr, g_dg = %9.5dr, g'_dg = %9.5dr, g_ga = %9.5dr, g'_ga = %9.5dr, g_2 = %9.5dr, g'_2 = %9.5dr\n",
            i, x[j], correction_dg(i, x[j]), correction_dg_deriv(i, x[j]), correction_ga(i, x[j]), correction_ga_deriv(i, x[j]),
            correction_g2(i, x[j]), correction_g2_deriv(i, x[j]) );
      }
      writeln();
    }

    for i in 1..9{
      var x : [0..i+1] real;
      x[0] = -1;
      x[i+1] = 1;
      x[1..i] = nodes_legendre_gauss(i);
      for j in 0..i+1 {
        writef( "n: %2i, x = %7.4dr, g_dg = %9.5dr, g'_dg = %9.5dr, g_ga = %9.5dr, g'_ga = %9.5dr, g_2 = %9.5dr, g'_2 = %9.5dr\n",
            i, x[j], correction_dg(i, x[j]), correction_dg_deriv(i, x[j]), correction_ga(i, x[j]), correction_ga_deriv(i, x[j]),
            correction_g2(i, x[j]), correction_g2_deriv(i, x[j]) );
      }
      writeln();
    }

    // Create a set with the cell topologies contained in the hypothetics test mesh
    var cellTopos : set(int);
    cellTopos.add(TOPO_LINE);
    cellTopos.add(TOPO_QUAD);

    var minOrder : int = 0;
    var maxOrder : int = 9;

    // Calculate the FR structures
    writeln();
    writeln("Flux correction initialized structure for FR:");
    writeln();

    init_correction(minOrder, maxOrder, cellTopos);
    writeln();
    writeln(flux_correction);
    writeln();
  }
}
