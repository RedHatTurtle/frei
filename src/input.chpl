prototype module Input
{
  use Random;
  use UnitTest;
  use Parameters.ParamInput;

  //parPhysics
  var nDims              : int = 1;
  var eqSet              : int = EQ_EULER;
  var initialConditions  : int = IC_SHOCKTUBE;
  var boundaryConditions : int = BC_DIRICHLET;

  //parFluid
  var fGamma : real = 1.4;           // Set the ratio of heat coefficients cp/cv
  var fR     : real = 287.0;         // Set the gas constant

  //parMesh
  var meshGen       : bool = false;
  var xMin          : real = -1.0;
  var xMax          : real =  1.0;
  var nCells        : int = 1000;
  var meshingScheme : int = MESH_UNIFORM;        // How to divide the domain into cells

  //parSpatial
  var spatialScheme     : int = SPATIAL_BEAMWARMING;
  var dissipationScheme : int = DISS_NONE;    // Type of numerical dissipation added

  //parFR
  var iOrder   : int = 3;                // Solution polynomial interpolation order
  var minOrder : int = 3;                // Minimum interpolation order to calculate coefficients
  var maxOrder : int = 3;                // Maximum interpolation order to calculate coefficients
  var distSP   : int = PTS_LEGENDRE;     // Distribution of SPs for SD method
  var frScheme : int = FR_DG;

  //parTime
  var timeScheme : int = TIME_RK_CLASSIC;
  var timeStep   : real = 1e-4;

  //parOutput
  var ioIter  : int = 0;             // Number of iterations between output dumps
  var maxIter : int = 0;             // Maximum number of iterations
  var ioTime  : real =-0.01;         // Time interval between output dumps
  var maxTime : real = 1.00;         // Maximum time to simulate

  //parRef
  var rhoRef  : real = 1.0;          // Reference density for non-dimensionalization
  var pRef    : real = 1.0;          // Reference pressure for non-dimensionalization

  //parInit
  var rhoLow  : real = 1.0;          // Density on the LOW pressure side of the shock tube
  var pLow    : real = 1.0;          // Pressure on the LOW pressure side of the shock tube
  var rhoHigh : real = 5.0;          // Density on the HIGH pressure side of the shock tube
  var pHigh   : real = 5.0;          // Pressure on the HIGH pressure side of the shock tube

  //////////////////////////////////////////////////////////////////////////////

  // Derived data
  var nEqs    : int = 1;
  var nDOF    : int = 1;
  var nGhosts : int = 1;
  var nPoints : int = nCells+1;

  //////////////////////////////////////////////////////////////////////////////

  proc indat(fileName : string)
  {
    use IO;
    use SysError;
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
    writeln("################################################################################");
    writeln("###   Printing input file used                                               ###");
    writeln("################################################################################");
    writeln();

    writeln(tomlData);

    writeln();
    writeln("################################################################################");
    writeln("###   End of input file                                                      ###");
    writeln("################################################################################");
    writeln();

    try {
      eqSet = tomlData["parPhysics"]!["eqSet"]!.i : int;
      initialConditions  = tomlData["parPhysics"]!["initialConditions"]!.i : int;
      boundaryConditions = tomlData["parPhysics"]!["boundaryConditions"]!.i : int;

      fGamma = tomlData["parFluid"]!["fGamma"]!.re : real;
      fR     = tomlData["parFluid"]!["fR"]!.re : real;

      xMin          = tomlData["parMesh"]!["xMin"]!.re : real;
      xMax          = tomlData["parMesh"]!["xMax"]!.re : real;
      nCells        = tomlData["parMesh"]!["nCells"]!.i : int;
      meshingScheme = tomlData["parMesh"]!["meshingScheme"]!.i : int;

      spatialScheme     = tomlData["parSpatial"]!["spatialScheme"]!.i : int;
      dissipationScheme = tomlData["parSpatial"]!["dissipationScheme"]!.i : int;

      iOrder = tomlData["parFR"]!["iOrder"]!.i : int;
      minOrder = tomlData["parFR"]!["minOrder"]!.i : int;
      maxOrder = tomlData["parFR"]!["maxOrder"]!.i : int;
      distSP = tomlData["parFR"]!["distSP"]!.i : int;
      distSP = tomlData["parFR"]!["frScheme"]!.i : int;

      timeStep = tomlData["parTime"]!["timeScheme"]!.i : int;
      timeStep = tomlData["parTime"]!["timeStep"]!.re  : real;

      ioIter  = tomlData["parOutput"]!["ioIter"]!.i : int;
      maxIter = tomlData["parOutput"]!["maxIter"]!.i : int;
      ioTime  = tomlData["parOutput"]!["ioTime"]!.re : real;
      maxTime = tomlData["parOutput"]!["maxTime"]!.re : real;

      rhoRef  = tomlData["parRef"]!["rhoRef"]!.re : real;
      pRef    = tomlData["parRef"]!["pRef"]!.re : real;

      rhoLow  = tomlData["parInit"]!["rhoLow"]!.re : real;
      pLow    = tomlData["parInit"]!["pLow"]!.re : real;
      rhoHigh = tomlData["parInit"]!["rhoHigh"]!.re : real;
      pHigh   = tomlData["parInit"]!["pHigh"]!.re : real;
    } catch {
      write("Error reading input file.");
    }

    writeln();
    writeln("################################################################################");
    writeln("###   Finished reading input file                                            ###");
    writeln("################################################################################");
    writeln();

    nPoints = nCells + 1;

    select eqSet {
      when EQ_CONVECTION   do nEqs=1;
      when EQ_INVBURGERS   do nEqs=1;
      when EQ_DIFFUSION    do nEqs=1;
      when EQ_LINBURGERS   do nEqs=1;
      when EQ_VISBURGERS   do nEqs=1;
      when EQ_EULER        do nEqs=3;
      when EQ_NAVIERSTOKES do nEqs=3;
    }
  }
}
