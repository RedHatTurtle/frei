module Control {
  //Runtime constants
  config const inputFile : string = "input.toml";

  proc main() {
    use Input;

    // Read inputa data
    indat(inputFile);

    // Process input data
    //nonDimensionalization;

    // Initialize the problem
    //allocateVariable;
    //initMesh;
    //initSollution;

    // Output the initial state
    //outputInit;

    // Start the solver itrations
    //solver;

    // Output the final state
    //outputFinal;
  }
}

module Input {
  proc indat(fileName : string) {
    use TOML;

    var tomlFile : file;
    try {
      tomlFile = open(fileName, iomode.r);
    } catch e : FileNotFoundError {
      writeln("Critical Error: Input file not found.");
      writeln("Stopping Execution immediately.");
    } catch {
      writeln("Unknown Error opening input file.");
      writeln("Stopping Execution immediately.");
    }

    var tomlData = parseToml(tomlFile);

    writeln();
    writeln("Printing input file used:");
    writeln(tomlData);
    writeln("End of input file.");
    writeln();
  }
}
