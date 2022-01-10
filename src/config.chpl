prototype module Config
{
  use Input;
  use Parameters;

  // Derived data
  //var nEqs  : int = 1;

  //var fpLoc : int = 1;
  //var spLoc : int = 1;

  enum MESH_TYPE   {GENERATE, GMESH2, GMESH4, CGNS};
  enum MESH_SCHEME {UNIFORM, RANDOM, RINGLEB};

  proc configure_solver() {}
}
