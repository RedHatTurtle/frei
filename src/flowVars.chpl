prototype module Solution {
  use Input;

  var nSolPts : int;
  var nFlxPts : int;

  var uSolPt : [1..nSolPts, 1..nEqs] real;
  var rSolPt : [1..nSolPts, 1..nEqs] real;
  var uFlxPt : [1..nFlxPts, 1..nEqs] real;
  var fFlxPt : [1..nFlxPts, 1..nEqs] real;
}
