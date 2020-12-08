prototype module Parameters
{
  // Equation Sets
  param eqConvection   : int = 1;    // Convection Eq            du/dt + c*(du/dx) = 0
  param eqInvBurgers   : int = 2;    // Inviscid Burgers Eq      du/dt + u*(du/dx) = 0
  param eqDiffusion    : int = 3;    // Diffusion Eq             du/dt             + k*(ddu/dxx) = 0
  param eqLinBurgers   : int = 4;    // Linear Burgers Eq        du/dt + c*(du/dx) + k*(ddu/dxx) = 0
  param eqVisBurgers   : int = 5;    // Viscous Burgers Eq       du/dt + u*(du/dx) + k*(ddu/dxx) = 0
  param eqEuler        : int = 6;    // Euler Eq
  param eqNavierStokes : int = 7;    // Navier-Stokes Eq

  // Initial Conditions
  param icSinusoidal      : int = 11;    // Sinusoidal Wave
  param icGaussPulse      : int = 12;    // Gaussian Pulse Wave
  param icEllipticalPulse : int = 13;    // Elliptic Pulse Wave
  param icSquare          : int = 14;    // Square Wave
  param icMixedWave       : int = 15;    // Mixed Wave
  param icShockTube       : int = 61;    // Shock Tube

  // Boundary Conditions
  param bcPeriodic  : int = 1;    // Periodic
  param bcDirichlet : int = 2;    // Dirichlet

  // Meshing Scheme
  param meshUniform    : int = 1;    // 
  param meshSinusoidal : int = 2;    // 
  param meshRandom     : int = 3;    // 

  // Spatial Schemes
  param spatialBeamWarming      : int =  1;    // 
  param spatialLaxWendroff      : int =  2;    // 
  param spatialMacCormak        : int =  3;    // 

  param spatialStegerWarming_o1 : int =  4;    // 
  param spatialStegerWarming_o2 : int =  5;    // 

  param spatialStegerWarmingO1  : int =  4;    // 
  param spatialStegerWarmingO2  : int =  5;    // 

  param spatial_Steger_Warming_o1 : int =  4;    // 
  param spatial_Steger_Warming_o2 : int =  5;    // 

  param spatialVanLeer_o1       : int =  6;    // 
  param spatialVanLeer_o2       : int =  7;    // 
  param spatialAUSM_o1          : int =  8;    // 
  param spatialAUSMplus_o1      : int =  9;    // 
  param spatialRoe              : int = 10;    // 
  param spatialFR               : int = 11;    // 

  // Time Schemes
  param timeEuler      : int = 0;
  param timeRK_Classic : int = 1;
  param timeTVDRK_o2s2 : int = 2;    // Flux Reconstruction
  param timeTVDRK_o2s3 : int = 3;    // Roe
  param timeTVDRK_o2s4 : int = 4;    // 
  param timeTVDRK_o2sN : int = 5;    // 
  param timeTVDRK_o3s3 : int = 6;    // 
  param timeTVDRK_o3s4 : int = 7;    // 
  param timeTVDRK_o3s5 : int = 8;    // 
  param timeTVDRK_o4s5 : int = 9;    // 

  // Inviscid Numerical Flux Schemes
  param fluxRusanov : int = 1;    // 
  param fluxRoe     : int = 2;    // 
  param fluxHLL     : int = 3;    // 
  param fluxHLLC    : int = 3;    // 
  param fluxRHLL    : int = 3;    // 

  // Viscou Numerical Flux Scheme
  param viscBR1 : int = 1;    // 
  param viscBR2 : int = 2;    // 
  param viscLDG : int = 3;    // 

  // Dissipation Scheme
  param dissNone    : int = 0;    // 
  param dissSecond  : int = 1;    // 
  param dissFourth  : int = 2;    // 
  param dissJameson : int = 3;    // 

  // Point distributions
  param ptsUniform          : int = 1;    // Uniform
  param ptsLegendre         : int = 2;    // Gauss-Legendre
  param ptsLegendreLobatto  : int = 3;    // Gauss-Legendre-Lobatto
  param ptsChebyshev        : int = 4;    // Gauss-Chebyshev
  param ptsChebyshevLobatto : int = 5;    // Gauss-Chebyshev-Lobatto

  // FR correction functions
  param fr_gDG    : int = 1;
  param fr_g2     : int = 2;
  param fr_gGauss : int = 3;
  param fr_gSD    : int = 4;
  param fr_g3     : int = 5;

  param frGDG    : int = 1;
  param frG2     : int = 2;
  param frGGauss : int = 3;
  param frGSD    : int = 4;
  param frG3     : int = 5;

  param frGdg    : int = 1;
  param frG2     : int = 2;
  param frGgauss : int = 3;
  param frGsd    : int = 4;
  param frG3     : int = 5;
}
