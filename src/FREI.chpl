/* Documentation for FREI */
module FREI
{
  //Runtime constants
  config const inputFile : string = "input.toml";

  proc main() {
    use Input;
    use Dimensional;
    use Spatial_Methods;

    // Read inputa data
    indat(inputFile);

    // Process input data
    dim2nondim;

    // Initialize the problem
    //allocate_variables;
    //init_mesh;
    //init_solution;

    // Output the initial state
    //output_init;

    // Start the solver itrations
    //solver;

    // Output the final state
    //output_final;
  }
}
