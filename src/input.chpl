module Input
{
  use Random;
  use UnitTest;
  use Parameters.ParamInput;
  import Mesh.faml_r;

  //parPhysics
  var nDims            : int = 3;
  var eqSet            : int = EQ_EULER;

  //parConvection
  var convectionSpeed : [1..3] real = [1.0, 0.0, 0.0];

  //parDiffusion
  var diffusionCoef : real = 1.0;

  // parFluid
  // SI defined constant as the product of the Avogadro number and the Boltzman constant. Calculated with values defined
  // in the 2018 CODATA report
  param gasR   : real = 8.314462618; // Ideal Gas Constant (2018 CODATA), J/(mol*K)
  // Gatley et al., A Twenty-First Century Molar Mass for Dry Air - DOI:10.1080/10789669.2008.10391032
  var fMolMass : real = 2.8966e-2;   // Molar Mass of the fluid mixture, kg/mol
  // Assuming an Ideal Gas
  var fGamma   : real = 1.4;         // Set the ratio of heat coefficients, Cp/Cv
  // Standard values for dry air from https://www.engineeringtoolbox.com/dry-air-properties-d_973.html
  var fVisc    : real = 1.78792e-5;  // Dynamic viscosity coefficient, kg/(m*s)
  var fKappa   : real = 2.52992e-2;  // Thermal Conductivity, W/(m*K)

  //parMesh
  var meshFormat    : int = MESH_GENERATE;
  var meshFileName  : string = "mesh.msh";
  var xMin          : real = -1.0;
  var xMax          : real =  1.0;
  var nCells        : int = 1000;
  var meshingScheme : int = MESH_GEN_UNIFORM;        // How to divide the domain into cells

  //parSpatial
  var spatialScheme     : int = SPATIAL_FR;
  var dissipationScheme : int = DISS_NONE;    // Type of numerical dissipation added
  var limiterScheme     : int = LIMITER_NONE; // Select limiter implementation

  //parFR
  var iOrder   : int = 3;                // Solution polynomial interpolation order
  var minOrder : int = 3;                // Minimum interpolation order to calculate coefficients
  var maxOrder : int = 3;                // Maximum interpolation order to calculate coefficients
  var distSP   : int = PTS_LEGENDRE;     // Distribution of SPs for SD method
  var frScheme : int = FR_DG;

  //parTime
  var timeScheme : int = TIME_TVDRK_O2S3;
  var timeStepMethod : int = DT_LOCAL_CFL;
  var timeStep   : real = 1e-6;
  var cfl : real = 1.0;

  //parOutput
  var ioIter  : int =   100;         // Number of iterations between output dumps
  var maxIter : int = 10000;         // Maximum number of iterations
  var ioTime  : real =-0.01;         // Time interval between output dumps
  var maxTime : real = 1.00;         // Maximum time to simulate
  var outError: int =     0;         // Calculate and output solution error

  //parRef
  var rhoRef  : real = 1.0;          // Reference density for non-dimensionalization
  var pRef    : real = 1.0;          // Reference pressure for non-dimensionalization

  //parFamilies
  var nFaml : int = 0;
  var famlList_d : domain(1);
  var famlList   : [famlList_d] faml_r;

  //////////////////////////////////////////////////////////////////////////////

  // Derived data
  var nEqs    : int = 5;
  var nDOF    : int = 1;
  var nGhosts : int = 1;
  var nPoints : int = nCells+1;

  var fR        : real = gasR/fMolMass;    // Calculate the specific gas constant for this fluid, J/(kg*K)
  var fCp       : real = fR/(1-1/fGamma);  // Calculate the heat coefficients at constant pressure Cp, J/(kg*K)
  var fCv       : real = fR/(fGamma-1);    // Calculate the heat coefficients at constant volume Cv, J/(kg*K)
  var fPrandtl  : real = fVisc*fCp/fKappa; // Calculate the Prandtl number of the fluid

  //////////////////////////////////////////////////////////////////////////////

  proc indat(fileName : string)
  {
    use IO;
    use OS;
    use TOML;

    var tomlFile : file;
    var tomlData : shared Toml?;

    // Open input config file
    try {
      writeln("Opening input file");
      tomlFile = open(fileName, ioMode.r);
    } catch e : FileNotFoundError {
      writeln("Critical Error: Input file not found.");
      writeln("Stopping Execution immediately.");
    } catch {
      writeln("Unknown Error opening input file.");
      writeln("Stopping Execution immediately.");
    }

    // Parse input config file contents into TOML object
    try {
      writeln("Parsing input file");
      tomlData = parseToml(tomlFile);
    } catch {
      writeln("Critical Error: Invalid TOML data");
    }

    // Dump parsed input config file
    {
      writeln();
      writeln("################################################################################");
      writeln("###   Input file dump                                                        ###");
      writeln("################################################################################");
      writeln(tomlData);
      writeln("################################################################################");
      writeln("###   End of input file                                                      ###");
      writeln("################################################################################");
      writeln();
    }

    // Copy configuration from TOML object to module variables
    try {
      writeln("Configuring solver");

      // parPhysics
      nDims             = tomlData!["parPhysics"]!["nDims"]!.i : int;
      eqSet             = tomlData!["parPhysics"]!["eqSet"]!.i : int;

      if eqSet == EQ_CONVECTION then
        convectionSpeed[1..3] = tomlData!["parConvection"]!["convectionSpeed"]!.arr[0..2]!.re : real;

      if eqSet == EQ_DIFFUSION then
        diffusionCoef = tomlData!["parDiffusion"]!["diffusionCoef"]!.re : real;

      // parFluid
      if (eqSet == EQ_EULER || eqSet == EQ_NAVIERSTOKES || eqSet == EQ_QUASI_1D_EULER)
      {
        fMolMass = tomlData!["parFluid"]!["fMolMass"]!.re : real;
        fGamma   = tomlData!["parFluid"]!["fGamma"]!.re : real;
        fVisc    = tomlData!["parFluid"]!["fVisc"]!.re : real;
        fKappa   = tomlData!["parFluid"]!["fKappa"]!.re : real;
      }

      // parMesh
      meshFormat    = tomlData!["parMesh"]!["meshFormat"]!.i : int;
      if meshFormat == MESH_GENERATE
      {
        xMin          = tomlData!["parMesh"]!["xMin"]!.re : real;
        xMax          = tomlData!["parMesh"]!["xMax"]!.re : real;
        nCells        = tomlData!["parMesh"]!["nCells"]!.i : int;
        meshingScheme = tomlData!["parMesh"]!["meshingScheme"]!.i : int;
      }
      else
      {
        meshFileName = tomlData!["parMesh"]!["meshFile"]!.s : string;
      }

      // parSpatial
      spatialScheme     = tomlData!["parSpatial"]!["spatialScheme"]!.i : int;
      dissipationScheme = tomlData!["parSpatial"]!["dissipationScheme"]!.i : int;
      limiterScheme     = tomlData!["parSpatial"]!["limiterScheme"]!.i : int;

      // parFR
      iOrder   = tomlData!["parFR"]!["iOrder"]!.i : int;
      minOrder = tomlData!["parFR"]!["minOrder"]!.i : int;
      maxOrder = tomlData!["parFR"]!["maxOrder"]!.i : int;
      distSP   = tomlData!["parFR"]!["distSP"]!.i : int;
      frScheme = tomlData!["parFR"]!["frScheme"]!.i : int;

      // parTime
      timeScheme = tomlData!["parTime"]!["timeScheme"]!.i : int;
      timeStepMethod = tomlData!["parTime"]!["timeStepMethod"]!.i : int;
      timeStep   = tomlData!["parTime"]!["timeStep"]!.re  : real;
      cfl        = tomlData!["parTime"]!["cfl"]!.re  : real;

      // parOutput
      ioIter   = tomlData!["parOutput"]!["ioIter"]!.i : int;
      maxIter  = tomlData!["parOutput"]!["maxIter"]!.i : int;
      ioTime   = tomlData!["parOutput"]!["ioTime"]!.re : real;
      maxTime  = tomlData!["parOutput"]!["maxTime"]!.re : real;
      outError = tomlData!["parOutput"]!["outError"]!.i : int;

      // parRef
      rhoRef   = tomlData!["parRef"]!["rhoRef"]!.re : real;
      pRef     = tomlData!["parRef"]!["pRef"]!.re : real;

      // parFamilies
      nFaml    = tomlData!["parFamilies"]!["nFamilies"]!.i : int;
      famlList_d = 1..nFaml;

      for famlIdx in famlList.domain
      {
        famlList[famlIdx].name           = tomlData!["parFamilies"]![famlIdx:string]!["familyName"]!.s : string;
        famlList[famlIdx].nDim           = tomlData!["parFamilies"]![famlIdx:string]!["familyDimension"]!.i : int;
        famlList[famlIdx].bocoType       = tomlData!["parFamilies"]![famlIdx:string]!["familyType"]!.i : int;
        famlList[famlIdx].bocoSubType    = tomlData!["parFamilies"]![famlIdx:string]!["familySubType"]!.i : int;
        for propIdx in tomlData!["parFamilies"]![famlIdx:string]!["familyParameters"]!.arr.domain do
          famlList[famlIdx].bocoProperties[propIdx+1] =
            tomlData!["parFamilies"]![famlIdx:string]!["familyParameters"]!.arr[propIdx]!.re : real;
      }
    } catch {
      write("Critical Error: Invalid solver configuration");
    }

    // Calculate some basic derived configurations
    select eqSet {
      when EQ_CONVECTION     do nEqs=nDims;
      when EQ_INVBURGERS     do nEqs=1;
      when EQ_DIFFUSION      do nEqs=nDims;
      when EQ_LINBURGERS     do nEqs=nDims;
      when EQ_VISBURGERS     do nEqs=1;
      when EQ_EULER          do nEqs=nDims+2;
      when EQ_NAVIERSTOKES   do nEqs=nDims+2;
      when EQ_QUASI_1D_EULER do nEqs=nDims+2;
      otherwise
      {
        writeln("Undefined number or equations \"nEqs\"");
      }
    }

    // Calculate derived fluid properties from input
    fR       = gasR/fMolMass;
    fCp      = fR/(1-1/fGamma);
    fCv      = fR/(fGamma-1);
    fPrandtl = fVisc*fCp/fKappa;

    // Mesh Generation parameters
    nPoints = nCells + 1;
  }
}
