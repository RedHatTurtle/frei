/* Documentation for FREI */
module FREI
{
  //Runtime constants
  config const inputFile : string = "input.toml";

  proc main() {
    use Parameters.ParamInput;
    use Config;
    use Input;
    use Interpolation;
    use Output;
    use Gmesh;
    use Mesh;
    use FRMesh;
    use FR;

    var iteration : int = 0;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();

    // 3. Read / define mesh
    var gmesh2 = new unmanaged gmesh2_c(nNodes=7, nElements=8, nFamilies=3);
    gmesh2.random1D(-1,1);

    // 4. Initialize the solver, pre calculate stuff
    init_interpolation(minOrder,maxOrder);

    // 5. Convert input mesh to solver mesh
    var fr_mesh = new unmanaged fr_mesh_c(nDims=1, nVars=3);
    fr_mesh.import_gmesh2(gmesh2);
    fr_mesh.allocate_vars();
    fr_mesh.set_points_locations();

    // Save mesh file in internal format

    // Initialize solution
    fr_mesh.solSP = initial_condition(IC_SHOCKTUBE, fr_mesh.xyzSP);

    // Save restart file

    // Output initial state
    iterOutput(iteration, fr_mesh);

    // Solve flow
    for iteration in 0..maxIter
    {
      // Iterate Solver (single or multiple time steps)
        // Interpolate solution to FPs
        // Calculate numerical flux
        // Calculate internal correction
        // Calculate interface correction

      // Print solver status / log

      // Save restart file

      // Output solution state
      iterOutput(0, fr_mesh);

      // Check input changes
    }

    // Output the final solution
    iterOutput(iteration, fr_mesh);
  }
}
