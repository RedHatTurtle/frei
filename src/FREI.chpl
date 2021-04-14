/* Documentation for FREI */
module FREI
{
  //Runtime constants
  config const inputFile : string = "input.toml";

  proc main() {
    use Parameters.ParamInput;
    use Config;
    use Input;
    use Gmesh;
    use Mesh;
    use FRMesh;
    use FR;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();

    // 3. Read / define mesh
    var gmesh2 = new unmanaged gmesh2_c(nNodes=7, nElements=8, nFamilies=3);
    gmesh2.random1D(-1,1);

    // 4. Initialize the solver, pre calculate stuff

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
    //fr_mesh.output_gnuplot();

    // Solve flow
      // Iterate Solver (single or multiple time steps
      // Print solver status
      // Output solution state
      // Save restart file
      // Check input changes

    // Output the final solution
  }
}
