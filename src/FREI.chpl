/* Documentation for FREI */
module FREI
{
  //Runtime constants
  config const inputFile : string = "input.toml";

  proc main() {
    use Input;
    use Config;

    // Read input data
    indat(inputFile);

    // Process input data into config
    configure_solver();

    // Initialize the problem
    //allocate_variables;
    //init_mesh;
    //dim2nondim;
    //init_solution;

    // Output the initial state
    //output_init;

    // Start the solver iterations
    //solver;

    // Output the final state
    //output_final;
  }
}
