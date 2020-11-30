prototype module Solution {
  use Input;

  var nSolPts : int(64);
  var nFlxPts : int(64);

  var uSolPt : [1..nSolPts, 1..nEqs] real(64);
  var rSolPt : [1..nSolPts, 1..nEqs] real(64);
  var uFlxPt : [1..nFlxPts, 1..nEqs] real(64);
  var fFlxPt : [1..nFlxPts, 1..nEqs] real(64);
}
