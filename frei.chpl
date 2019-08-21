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

module Parameters {
  // Equation Sets
  param eqConvection   : int(8) = 1;    // Convection Eq            du/dt + c*(du/dx) = 0
  param eqInvBurgers   : int(8) = 2;    // Inviscid Burgers Eq      du/dt + u*(du/dx) = 0
  param eqDiffusion    : int(8) = 3;    // Diffusion Eq             du/dt             + k*(ddu/dxx) = 0
  param eqLinBurgers   : int(8) = 4;    // Linear Burgers Eq        du/dt + c*(du/dx) + k*(ddu/dxx) = 0
  param eqVisBurgers   : int(8) = 5;    // Viscous Burgers Eq       du/dt + u*(du/dx) + k*(ddu/dxx) = 0
  param eqEuler        : int(8) = 6;    // Euler Eq
  param eqNavierStokes : int(8) = 7;    // Navier-Stokes Eq

  // Initial Conditions
  param icSinusoidal      : int(8) = 11;    // Sinusoidal Wave
  param icGaussPulse      : int(8) = 12;    // Gaussian Pulse Wave
  param icEllipticalPulse : int(8) = 13;    // Elliptic Pulse Wave
  param icSquare          : int(8) = 14;    // Square Wave
  param icMixedWave       : int(8) = 15;    // Mixed Wave
  param icShockTube       : int(8) = 61;    // Shock Tube

  // Boundary Conditions
  param bcPeriodic  : int(8) = 1;    // Periodic
  param bcDirichlet : int(8) = 2;    // Dirichlet

  // Meshing Scheme
  param meshUniform    : int(8) = 1;    // 
  param meshSinusoidal : int(8) = 2;    // 
  param meshRandom     : int(8) = 3;    // 

  // Spatial Schemes
  param spatialBeamWarming      : int(8) =  1;    // 
  param spatialLaxWendroff      : int(8) =  2;    // 
  param spatialMacCormak        : int(8) =  3;    // 

  param spatialStegerWarming_o1 : int(8) =  4;    // 
  param spatialStegerWarming_o2 : int(8) =  5;    // 

  param spatialStegerWarmingO1  : int(8) =  4;    // 
  param spatialStegerWarmingO2  : int(8) =  5;    // 

  param spatial_Steger_Warming_o1 : int(8) =  4;    // 
  param spatial_Steger_Warming_o2 : int(8) =  5;    // 

  param spatialVanLeer_o1       : int(8) =  6;    // 
  param spatialVanLeer_o2       : int(8) =  7;    // 
  param spatialAUSM_o1          : int(8) =  8;    // 
  param spatialAUSMplus_o1      : int(8) =  9;    // 
  param spatialRoe              : int(8) = 10;    // 
  param spatialFR               : int(8) = 11;    // 

  // Time Schemes
  param timeEuler      : int(8) = 0;
  param timeRK_Classic : int(8) = 1;
  param timeTVDRK_o2s2 : int(8) = 2;    // Flux Reconstruction
  param timeTVDRK_o2s3 : int(8) = 3;    // Roe
  param timeTVDRK_o2s4 : int(8) = 4;    // 
  param timeTVDRK_o2sN : int(8) = 5;    // 
  param timeTVDRK_o3s3 : int(8) = 6;    // 
  param timeTVDRK_o3s4 : int(8) = 7;    // 
  param timeTVDRK_o3s5 : int(8) = 8;    // 
  param timeTVDRK_o4s5 : int(8) = 9;    // 

  // Inviscid Numerical Flux Schemes
  param fluxRusanov : int(8) = 1;    // 
  param fluxRoe     : int(8) = 2;    // 
  param fluxHLL     : int(8) = 3;    // 
  param fluxHLLC    : int(8) = 3;    // 
  param fluxRHLL    : int(8) = 3;    // 

  // Viscou Numerical Flux Scheme
  param viscBR1 : int(8) = 1;    // 
  param viscBR2 : int(8) = 2;    // 
  param viscLDG : int(8) = 3;    // 

  // Dissipation Scheme
  param dissNone    : int(8) = 0;    // 
  param dissSecond  : int(8) = 1;    // 
  param dissFourth  : int(8) = 2;    // 
  param dissJameson : int(8) = 3;    // 

  // Point distributions
  param ptsUniform          : int(8) = 1;    // Uniform
  param ptsLegendre         : int(8) = 2;    // Gauss-Legendre
  param ptsLegendreLobatto  : int(8) = 3;    // Gauss-Legendre-Lobatto
  param ptsChebyshev        : int(8) = 4;    // Gauss-Chebyshev
  param ptsChebyshevLobatto : int(8) = 5;    // Gauss-Chebyshev-Lobatto

  // FR correction functions
  param fr_gDG    : int(8) = 1;
  param fr_g2     : int(8) = 2;
  param fr_gGauss : int(8) = 3;
  param fr_gSD    : int(8) = 4;
  param fr_g3     : int(8) = 5;

  param frGDG    : int(8) = 1;
  param frG2     : int(8) = 2;
  param frGGauss : int(8) = 3;
  param frGSD    : int(8) = 4;
  param frG3     : int(8) = 5;

  param frGdg    : int(8) = 1;
//param frG2     : int(8) = 2;
  param frGgauss : int(8) = 3;
  param frGsd    : int(8) = 4;
//param frG3     : int(8) = 5;
}

module Input {
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

module Mesh {
  var nBoundaries : int(8);
  var nFaces : int(64);
  var nCells : int(64);

  var allFaces : range = 1..nFaces;
  var allCells : range = 1..nCells;

  //////////////////////////////////////////////////////////////////////////////

  record bFace_ {
    var cell : int(64);
    var bcType : int(8);
    var cell_face : int(8); // 1 for left, 2 for right.
    var fp : int(64);
  }

  var bFace : [1..nBoundaries] bFace_;

  //////////////////////////////////////////////////////////////////////////////

  record face_cell_ {
    var cell : int(64);
    var cell_face : int(8);
    var fp : int(64);
  }

  record face_ {
    var x : real(64);
    var left  : face_cell_;
    var right : face_cell_;
  }

  var face : [1..nFaces] face_;

  //////////////////////////////////////////////////////////////////////////////

  record cell_face_ {
    var idx  : int(64);
    var side : int(8);
    var fp   : int(8);
  }

  record cell_ {
    var begSP : int(64);
    var endSP : int(64);
    var face  : [1..2] cell_face_;
  }

  var cell : [1..nCells] cell_;

  //////////////////////////////////////////////////////////////////////////////
}

module Solution {
  var nSolPts : int(64);
  var nFlxPts : int(64);

  var uSolPt : [1..nSolPts, 1..nEqs] real(64);
  var rSolPt : [1..nSolPts, 1..nEqs] real(64);
  var uFlxPt : [1..nFlxPts, 1..nEqs] real(64);
  var fFlxPt : [1..nFlxPts, 1..nEqs] real(64);
}
