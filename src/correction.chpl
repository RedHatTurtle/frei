prototype module Correction
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

  class flux_correction_c
  {
    var correction_d : domain(2); // {nFPs, nSPs}
    var correction : [correction_d] real;
  }

  // Define type for the flux correction structure. Add "?" to allow default initialization to nil.
  type flux_correction_t = unmanaged flux_correction_c?;

  // Perhaps it might be useful in the future to have a sparse domain with the cell topologies present in the mesh
  //var flux_correction_d : sparse domain(2);

  var flux_correction_d : domain(2); // {elemType, interpOrder}

  var flux_correction : [flux_correction_d] flux_correction_t;

  proc init_correction()
  {
    use Parameters.Input;
    use Parameters.Mesh;
    use Polynomials;

    // Allocate flux correction structure
    var elemTopos : range = 2..2;
    var solOrders : range = 1..9;
    flux_correction_d = {elemTopos, solOrders};

    for (elemTopo, solOrder) in flux_correction.domain
    {
      select elemTopo
      {
        when TOPO_LINE
        {
          var spCnt : range = 1..solOrder;
          var fpCnt : range = 1..2;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh uses has SPs on Legendre roots. xD
          var spLoc : [spCnt] real = nodes_legendre_gauss(solOrder);
          var fpLoc : [fpCnt] real = [-0.5, 0.5];

          flux_correction[elemTopo, solOrder] = new flux_correction_t({fpCnt, spCnt})!;

          for (fp,sp) in {fpCnt,spCnt} do
            select frScheme
            {
                when FR_DG do
                  flux_correction[elemTopo, solOrder]!.correction[fp,sp] = correction_dg(solOrder, spLoc[sp]);
                when FR_GA do
                  flux_correction[elemTopo, solOrder]!.correction[fp,sp] = correction_ga(solOrder, spLoc[sp]);
                when FR_G2 do
                  flux_correction[elemTopo, solOrder]!.correction[fp,sp] = correction_g2(solOrder, spLoc[sp]);
            }
        }
        when TOPO_TRIA {}
        when TOPO_QUAD {}
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise do writeln("Unsupported mesh element found at flux correction initialization.");
      }
    }
  }

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
    use IO;
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

    writeln();
    writeln("Flux correction initialized structure for FR:");
    writeln();

    init_correction();
    writeln(flux_correction);
    writeln();
  }
}
