import BLAS;
private param usingBLAS = BLAS.header != '';

config const n : int = 3;

proc main()
{
    use Time;
    use Random;
    use LinearAlgebra;
    use Inverse;

    // This breaks the linear algebra module
    //var testDoms : [1..2] domain(2) = [{0..<n, 0..<n},
    //                                   {0..<n, 1..n},
    //                                   {1..n , 0..<n},
    //                                   {1..n , 1..n}];

    var testDoms : [1..2] domain(2) = [{0..<n, 0..<n},
                                       {1..n , 1..n}];

    for matrix_d in testDoms
    {
      var matrix : [matrix_d] real;

      // Generate random test matrix
      fillRandom(matrix, 17);

      // Make the matrix slightly ill conditioned
      //matrix /= 1000.0;
      //matrix += 1.0;

      writeln("Test domain: ", matrix.domain);

      writeln("  - Matrix times Inverse:");
      {
        const inverse = Inverse.inv(matrix);
        writeln("    Matrix  domain: ", matrix.domain);
        writeln("    Inverse domain: ", inverse.domain);
        writeln("    Product domain: ", dot(matrix, inverse).domain);
        for lin in matrix.domain.dim(0)
        {
          for col in matrix.domain.dim(1) do
            writef("%10.4dr", matrix[lin, col]);
          writef("  |  ");
          for col in matrix.domain.dim(1) do
            writef("%24.16er", inverse[lin, col]);
          writef("  |  ");
          for col in matrix.domain.dim(1) do
            writef("%16.8er", dot(matrix, inverse)[lin, col]);

          writeln();
        }
      }
      writeln();
      writeln("  - Matrix times Inverse_dot:");
      {
        const inverse = Inverse.inv_dot(matrix);
        writeln("    Matrix  domain: ", matrix.domain);
        writeln("    Inverse domain: ", inverse.domain);
        writeln("    Product domain: ", dot(matrix, inverse).domain);
        for lin in matrix.domain.dim(0)
        {
          for col in matrix.domain.dim(1) do
            writef("%10.4dr", matrix[lin, col]);
          writef("  |  ");
          for col in matrix.domain.dim(1) do
            writef("%24.16er", inverse[lin, col]);
          writef("  |  ");
          for col in matrix.domain.dim(1) do
            writef("%16.8er", dot(matrix, inverse)[lin, col]);

          writeln();
        }
      }
      if usingBLAS then
      {
        writeln();
        writeln("  - Matrix times Lapack Inverse:");
        {
          const inverse = LinearAlgebra.inv(matrix);
          writeln("    Matrix  domain: ", matrix.domain);
          writeln("    Inverse domain: ", inverse.domain);
          writeln("    Product domain: ", dot(matrix, inverse).domain);
          for lin in matrix.domain.dim(0)
          {
            for col in matrix.domain.dim(1) do
              writef("%10.4dr", matrix[lin, col]);
            writef("  |  ");
            for col in matrix.domain.dim(1) do
              writef("%24.16er", inverse[lin, col]);
            writef("  |  ");
            for col in matrix.domain.dim(1) do
              writef("%16.8er", dot(matrix, inverse)[lin, col]);

            writeln();
          }
        }
      }
      writeln();
    }
}
