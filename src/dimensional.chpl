prototype module Dimensional
{
  proc dim2nondim
  {
    use Input;
    use Mesh;

    rhoLow  = rhoLow  / rhoRef;
    pLow    = pLow    / pRef;
    rhoHigh = rhoHigh / rhoRef;
    pHigh   = pHigh   / pRef;
  }

  proc nondim2dim
  {
    use Input;

    rhoLow  = rhoLow  * rhoRef;
    pLow    = pLow    * pRef;
    rhoHigh = rhoHigh * rhoRef;
    pHigh   = pHigh   * pRef;
  }
}
