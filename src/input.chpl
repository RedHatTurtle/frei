prototype module Input
{
  use Random;
  use UnitTest;
  use Parameters.ParamInput;
  import Mesh.faml_r;

  //parPhysics
  var nDims            : int = 1;
  var eqSet            : int = EQ_EULER;

  //parConvection
  var convectionSpeed : [1..3] real = [1.0, 0.0, 0.0];

  //parDiffusion
  var diffusionCoef : real = 1.0;

  //parFluid
  var fGamma : real =    1.4;        // Set the ratio of heat coefficients Cp/Cv
  var fCp    : real = 1004.5;        // Set the heat coefficients at constant pressure Cp
  var fCv    : real =  717.5;        // Set the heat coefficients at constant volume Cv
  var fR     : real =  287.0;        // Set the specific gas constant J/(kg*K)

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

  //parFR
  var iOrder   : int = 3;                // Solution polynomial interpolation order
  var minOrder : int = 3;                // Minimum interpolation order to calculate coefficients
  var maxOrder : int = 3;                // Maximum interpolation order to calculate coefficients
  var distSP   : int = PTS_LEGENDRE;     // Distribution of SPs for SD method
  var frScheme : int = FR_DG;

  //parTime
  var timeScheme : int = TIME_TVDRK_O2S3;
  var timeStep   : real = 1e-6;

  //parOutput
  var ioIter  : int =   100;         // Number of iterations between output dumps
  var maxIter : int = 10000;         // Maximum number of iterations
  var ioTime  : real =-0.01;         // Time interval between output dumps
  var maxTime : real = 1.00;         // Maximum time to simulate

  //parRef
  var rhoRef  : real = 1.0;          // Reference density for non-dimensionalization
  var pRef    : real = 1.0;          // Reference pressure for non-dimensionalization

  //parFamilies
  var nFaml : int = 1;
  var famlList_d : domain(1);
  var famlList   : [famlList_d] faml_r;

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
    var tomlData : unmanaged Toml?;

    try {
      writeln("Opening input file");
      tomlFile = open(fileName, iomode.r);
    } catch e : FileNotFoundError {
      writeln("Critical Error: Input file not found.");
      writeln("Stopping Execution immediately.");
    } catch {
      writeln("Unknown Error opening input file.");
      writeln("Stopping Execution immediately.");
    }

    try {
      writeln("Parsing input file");
      tomlData = parseToml(tomlFile);
    } catch {
      writeln("Critical Error: Invalid TOML data");
    }

    writeln();
    writeln("################################################################################");
    writeln("###   Input file dump                                                        ###");
    writeln("################################################################################");
    writeln(tomlData);
    writeln("################################################################################");
    writeln("###   End of input file                                                      ###");
    writeln("################################################################################");
    writeln();

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
        fGamma = tomlData!["parFluid"]!["fGamma"]!.re : real;
        fCp    = tomlData!["parFluid"]!["fCp"]!.re : real;
        fCv    = tomlData!["parFluid"]!["fCv"]!.re : real;
        fR     = tomlData!["parFluid"]!["fR"]!.re : real;
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

      // parSpacial
      spatialScheme     = tomlData!["parSpatial"]!["spatialScheme"]!.i : int;
      dissipationScheme = tomlData!["parSpatial"]!["dissipationScheme"]!.i : int;

      // parFR
      iOrder   = tomlData!["parFR"]!["iOrder"]!.i : int;
      minOrder = tomlData!["parFR"]!["minOrder"]!.i : int;
      maxOrder = tomlData!["parFR"]!["maxOrder"]!.i : int;
      distSP   = tomlData!["parFR"]!["distSP"]!.i : int;
      frScheme = tomlData!["parFR"]!["frScheme"]!.i : int;

      // parTime
      timeScheme = tomlData!["parTime"]!["timeScheme"]!.i : int;
      timeStep   = tomlData!["parTime"]!["timeStep"]!.re  : real;

      // parOutput
      ioIter   = tomlData!["parOutput"]!["ioIter"]!.i : int;
      maxIter  = tomlData!["parOutput"]!["maxIter"]!.i : int;
      ioTime   = tomlData!["parOutput"]!["ioTime"]!.re : real;
      maxTime  = tomlData!["parOutput"]!["maxTime"]!.re : real;

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

    nPoints = nCells + 1;

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
  }
}
