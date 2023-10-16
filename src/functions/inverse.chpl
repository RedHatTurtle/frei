module Inverse
{
  import BLAS;
  private param usingBLAS = BLAS.header != '';

  proc inv_triu(const ref U : [?U_d] real, unit_diag : bool = false) : [U_d] real
  {
    const n : int = U_d.dim(0).size;
    var U_inv : [U_d] real;

    for j in 0..<n do
    {
      if unit_diag then
        U_inv(j, j) = 1.0;
      else
        U_inv(j, j) = 1.0 / U(j, j);

      if j > 0 then
        forall i in 0..<j
        {
          for k in 0..<j do
            U_inv[  i, j] += U_inv[i, k]*U[k, j];
          U_inv[    i, j] *= -U_inv( j, j);
        }
    }

    return U_inv;
  }

  proc inv_tril(const ref L : [?L_d] real, unit_diag : bool = false) : [L_d] real
  {
    const n : int = L_d.dim(0).size;
    var L_inv : [L_d] real;

    for j in 0..<n by -1 do
    {
      if unit_diag then
        L_inv[j, j] = 1.0;
      else
        L_inv[j, j] = 1.0 / L[j, j];

      if j+1 < n then
        forall i in j+1..<n
        {
          for k in j+1..<n do
            L_inv[i, j] += L_inv[i, k]*L[k, j];
          L_inv[i, j] *= -L_inv(j, j);
        }
    }

    return L_inv;
  }

  proc inv_triu_dot(const ref U : [?U_d] real, unit_diag : bool = false) : [U_d] real
  {
    const n : int = U_d.dim(0).size;
    var U_inv : [U_d] real;

    for j in 0..<n do
    {
      if unit_diag then
        U_inv(j, j) = 1.0;
      else
        U_inv(j, j) = 1.0 / U(j, j);

      if j > 0 then
        U_inv[0..<j, j] = -U_inv(j, j) * dot(U_inv[0..<j, 0..<j], U[0..<j, j]);
    }



    return U_inv;
  }

  proc inv_tril_dot(const ref L : [?L_d] real, unit_diag : bool = false) : [L_d] real
  {
    const n : int = L_d.dim(0).size;
    var L_inv : [L_d] real;

    for j in 0..<n by -1 do
    {
      if unit_diag then
        L_inv[j, j] = 1.0;
      else
        L_inv[j, j] = 1.0 / L[j, j];

      if j+1 < n then
        L_inv[j+1..<n, j] = -L_inv(j, j) * dot(L_inv[j+1..<n, j+1..<n], L[j+1..<n, j]);
    }



    return L_inv;
  }

  proc inv(const ref M : [?M_d] real) : [M_d] real
  {
    const n : int = M_d.dim(0).size;

    // Allocate return matrix
    var mInv : [M_d.dim(1), M_d.dim(0)] real;

    // Calculate LU factorization
    const (luM, iPiv) = lu(M.reindex(0..<M_d.shape(0), 0..<M_d.shape(1)));

    if (!usingBLAS || M_d.size <= 36**2) then
    {
      // Invert each LU factor
      const uInv = inv_triu(U=luM);
      const lInv = inv_tril(L=luM, unit_diag = true);

      // Multiply inverted factors and pivot
      forall i in 0..<n do
        for j in 0..<n do
          for k in 0..<n do
            mInv[M_d.dim(0).orderToIndex(i),M_d.dim(1).orderToIndex(iPiv[j])] += uInv[i,k] * lInv[k,j];
    }
    else
    {
      // Invert each LU factor
      var uInv : [M_d] real = inv_triu_dot(U=luM);
      var lInv : [M_d] real = inv_tril_dot(L=luM, unit_diag = true);

      // Multiply inverted factors and pivot
      for j in 0..<n do
        mInv[.., M_d.dim(1).orderToIndex(iPiv[j])] = dot(uInv, lInv[.., j]);
    }

    return mInv;
  }

  proc inv_dot(const ref M : [?M_d]) : [M_d] real
  {
    const n : int = M_d.dim(0).size;

    // Allocate return matrix
    var mInv : [M_d.dim(1), M_d.dim(0)] real;

    // Calculate LU factorization
    const (luM, iPiv) = lu(M.reindex(0..<M_d.shape(0), 0..<M_d.shape(1)));

    // Invert each LU factor
    var uInv = inv_triu_dot(U=luM);
    var lInv = inv_tril_dot(L=luM, unit_diag = true);

    // Multiply inverted factors and pivot
    for j in 0..<n do
      mInv[.., M_d.dim(1).orderToIndex(iPiv[j])] = dot(uInv, lInv[.., j]);

    return mInv;
  }
}
