module Determinant
{
  proc det(const ref matrix : [] real) : real
  {
    import LinearAlgebra.det;

    var jacobian : real;

    if matrix.size == 1 then
      jacobian = matrix[1,1];
    else if matrix.size == 4 then
      jacobian = matrix[1,1]*matrix[2,2] - matrix[1,2]*matrix[2,1];
    else if matrix.size == 9 then
      jacobian = matrix[1,1]*(matrix[2,2]*matrix[3,3] - matrix[2,3]*matrix[3,2])
                +matrix[1,2]*(matrix[2,3]*matrix[3,1] - matrix[2,1]*matrix[3,3])
                +matrix[1,3]*(matrix[2,1]*matrix[3,2] - matrix[2,2]*matrix[3,1]);
    else
      jacobian = det(matrix);

    return jacobian;
  }
}
