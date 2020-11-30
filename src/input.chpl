prototype module Input
{
  use Parameters;

  //parPhysics
  var eqset              : int(8) = eqEuler;
  var initialConditions  : int(8) = icShockTube;
  var boundaryConditions : int(8) = bcDirichlet;

  //parFluid
  var fGamma : real(64) = 1.4;           // Set the ratio of heat coefficients cp/cv
  var fR     : real(64) = 287.0;         // Set the gas constant

  //parMesh
  var xMin          : real(64) = -1.0;
  var xMax          : real(64) =  1.0;
  var nCells        : int(64) = 1000;
  var meshingScheme : int(8) = meshUniform;        // How to divide the domain into cells

  //parSpatial
  var spatialScheme     : int(8) = spatialBeamWarming;
  var dissipationScheme : int(8) = dissNone;    // Type of numerical dissipation added

  //parFR
  var iOrder : int(8) = 3;               // Solution polynomial interpolation order
  var distSP : int(8) = ptsLegendre;     // Distribution of SPs for SD method
  var frScheme : int(8) = fr_gDG;

  //parTime
  var timeScheme : int(8) = timeRK_Classic;
  var timeStep   : real(64) = 1e-4;

  //parOutput
  var ioIter  : int(64) = 0;             // Number of iterations between output dumps
  var maxIter : int(64) = 0;             // Maximum number of iterations
  var ioTime  : real(64) =-0.01;         // Time interval between output dumps
  var maxTime : real(64) = 1.00;         // Maximum time to simulate

  //parRef
  var rhoRef  : real(64) = 1.0;          // Reference density for non-dimensionalization
  var pRef    : real(64) = 1.0;          // Reference pressure for non-dimensionalization

  //parInit
  var rhoLow  : real(64) = 1.0;          // Density on the LOW pressure side of the shock tube
  var pLow    : real(64) = 1.0;          // Pressure on the LOW pressure side of the shock tube
  var rhoHigh : real(64) = 5.0;          // Density on the HIGH pressure side of the shock tube
  var pHigh   : real(64) = 5.0;          // Pressure on the HIGH pressure side of the shock tube

  //////////////////////////////////////////////////////////////////////////////

  // Derived data
  var nEqs    : int(8) = 1;
  var nDOF    : int(8) = 1;
  var nGhosts : int(8) = 1;
  var nPoints : int(64) = nCells+1;

  //////////////////////////////////////////////////////////////////////////////

  proc indat(fileName : string)
  {
    use TOML;
    use Mesh;

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
    writeln(tomlData);
    writeln("################################################################################");
    writeln("###   End of input file                                                      ###");
    writeln("################################################################################");
    writeln();

    try {
      eqset = tomlData["parPhysics"]["eqset"].i : int(8);
      initialConditions  = tomlData["parPhysics"]["initialConditions"].i : int(8);
      boundaryConditions = tomlData["parPhysics"]["boundaryConditions"].i : int(8);

      fGamma = tomlData["parFluid"]["fGamma"].re : real(64);
      fR     = tomlData["parFluid"]["fR"].re : real(64);

      xMin          = tomlData["parMesh"]["xMin"].re : real(64);
      xMax          = tomlData["parMesh"]["xMax"].re : real(64);
      nCells        = tomlData["parMesh"]["nCells"].i : int(64);
      meshingScheme = tomlData["parMesh"]["meshingScheme"].i : int(8);

      spatialScheme     = tomlData["parSpatial"]["spatialScheme"].i : int(8);
      dissipationScheme = tomlData["parSpatial"]["dissipationScheme"].i : int(8);

      iOrder = tomlData["parFR"]["iOrder"].i : int(8);
      distSP = tomlData["parFR"]["distSP"].i : int(8);
      distSP = tomlData["parFR"]["frScheme"].i : int(8);

      timeStep = tomlData["parTime"]["timeScheme"].i : int(8);
      timeStep = tomlData["parTime"]["timeStep"].re  : real(64);

      ioIter  = tomlData["parOutput"]["ioIter"].i : int(64);
      maxIter = tomlData["parOutput"]["maxIter"].i : int(64);
      ioTime  = tomlData["parOutput"]["ioTime"].re : real(64);
      maxTime = tomlData["parOutput"]["maxTime"].re : real(64);

      rhoRef  = tomlData["parRef"]["rhoRef"].re : real(64);
      pRef    = tomlData["parRef"]["pRef"].re : real(64);

      rhoLow  = tomlData["parInit"]["rhoLow"].re : real(64);
      pLow    = tomlData["parInit"]["pLow"].re : real(64);
      rhoHigh = tomlData["parInit"]["rhoHigh"].re : real(64);
      pHigh   = tomlData["parInit"]["pHigh"].re : real(64);
    } catch {
      write("Error reading input file.");
    }
    writeln();
    writeln("################################################################################");
    writeln("###   Finished reading input file                                            ###");
    writeln("################################################################################");
    writeln();

    nPoints = nCells + 1;

    select eqset {
      when eqConvection   do nEqs=1;
      when eqInvBurgers   do nEqs=1;
      when eqDiffusion    do nEqs=1;
      when eqLinBurgers   do nEqs=1;
      when eqVisBurgers   do nEqs=1;
      when eqEuler        do nEqs=3;
      when eqNavierStokes do nEqs=3;
    }
  }
}
