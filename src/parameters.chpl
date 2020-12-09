prototype module Parameters
{
  prototype module Input
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

    // Parametric internal 1D meshing methods
    param meshUniform    : int = 1;    // 
    param meshRandom     : int = 2;    // 

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
  }

  prototype module Gmesh
  {

    // Gmesh element topologies
    //param GMESH_TOPO_PNT     : int = 1;
    //param GMESH_TOPO_LIN     : int = 2;
    //param GMESH_TOPO_TRI     : int = 3;
    //param GMESH_TOPO_QUA     : int = 4;
    //param GMESH_TOPO_TET     : int = 5;
    //param GMESH_TOPO_PYR     : int = 6;
    //param GMESH_TOPO_PRI     : int = 7;
    //param GMESH_TOPO_HEX     : int = 8;
    //param GMESH_TOPO_POLYG   : int = 9;
    //param GMESH_TOPO_POLYH   : int = 10;
    //param GMESH_TOPO_XFEM    : int = 11;
    //param GMESH_TOPO_MINI    : int = 12;
    //param GMESH_TOPO_TRIH    : int = 13;

    // Gmesh Element Types
    //
    //Dimension | Shape         | Linear    | Quadratic           | Cubic                        | Quartic
    //---------------------------------------------------------------------------------------------------------------------------
    //  1-D     |   Line        |   LIN_2   |   LIN_3             |   LIN_4                      |   LIN_5
    //  2-D     |   Triangle    |   TRI_3   |   TRI_6             |   TRI_9    TRI_10            |   TRI_12     TRI_15
    //          |   Quadrangle  |   QUAD_4  |   QUAD_8   QUAD_9   |   QUAD_12  QUAD_16           |   QUAD_P4_16 QUAD_25
    //  3-D     |   Tetrahedron |   TETRA_4 |   TETRA_10          |   TETRA_16 TETRA_20          |   TETRA_22   TETRA_34 TETRA_35
    //          |   Pyramid     |   PYRA_5  |   PYRA_13  PYRA_14  |   PYRA_21  PYRA_29  PYRA_30  |   PYRA_P4_29 PYRA_50  PYRA_55
    //          |   Pentahedron |   PENTA_6 |   PENTA_15 PENTA_18 |   PENTA_24 PENTA_38 PENTA_40 |   PENTA_33   PENTA_66 PENTA_75
    //          |   Hexahedron  |   HEXA_8  |   HEXA_20  HEXA_27  |   HEXA_32  HEXA_56  HEXA_64  |   HEXA_44    HEXA_98  HEXA_125

    param GMESH_PNT      : int = 15;

    //     Line:                 Line3:          Line4:
    //
    //       v
    //       ^
    //       |
    //       |
    // 0-----+-----1 --> u   0----2----1     0---2---3---1

    param GMESH_LIN_2    : int =   1;
    param GMESH_LIN_3    : int =   8;
    param GMESH_LIN_4    : int =  26;
    param GMESH_LIN_5    : int =  27;
    param GMESH_LIN_6    : int =  28;
    param GMESH_LIN_7    : int =  62;
    param GMESH_LIN_8    : int =  63;
    param GMESH_LIN_9    : int =  64;
    param GMESH_LIN_10   : int =  65;
    param GMESH_LIN_11   : int =  66;
    param GMESH_LIN_B    : int =  67;
    param GMESH_LIN_C    : int =  70;
    param GMESH_LIN_1    : int =  84;
    param GMESH_LIN_SUB  : int = 134;

    // Triangle:               Triangle6:          Triangle9/10:          Triangle12/15:
    //
    // v
    // ^                                                                   2
    // |                                                                   | \
    // 2                       2                    2                      9   8
    // |`\                     |`\                  | \                    |     \
    // |  `\                   |  `\                7   6                 10 (14)  7
    // |    `\                 5    `4              |     \                |         \
    // |      `\               |      `\            8  (9)  5             11 (12) (13) 6
    // |        `\             |        `\          |         \            |             \
    // 0----------1 --> u      0-----3----1         0---3---4---1          0---3---4---5---1

    param GMESH_TRI_3    : int =   2;
    param GMESH_TRI_6    : int =   9;
    param GMESH_TRI_9    : int =  20;
    param GMESH_TRI_10   : int =  21;
    param GMESH_TRI_12   : int =  22;
    param GMESH_TRI_15   : int =  23;
    param GMESH_TRI_15I  : int =  24;
    param GMESH_TRI_21   : int =  25;
    param GMESH_TRI_28   : int =  42;
    param GMESH_TRI_36   : int =  43;
    param GMESH_TRI_45   : int =  44;
    param GMESH_TRI_55   : int =  45;
    param GMESH_TRI_66   : int =  46;
    param GMESH_TRI_18   : int =  52;
    param GMESH_TRI_21I  : int =  53;
    param GMESH_TRI_24   : int =  54;
    param GMESH_TRI_27   : int =  55;
    param GMESH_TRI_30   : int =  56;
    param GMESH_TRI_B    : int =  68;
    param GMESH_TRI_1    : int =  85;
    param GMESH_TRI_SUB  : int = 135;
    param GMESH_TRI_MINI : int = 138;
    param GMESH_TRIH_4   : int = 140;

    // Quadrangle:            Quadrangle8:            Quadrangle9:
    //
    //       v
    //       ^
    //       |
    // 3-----------2          3-----6-----2           3-----6-----2
    // |     |     |          |           |           |           |
    // |     |     |          |           |           |           |
    // |     +---- | --> u    7           5           7     8     5
    // |           |          |           |           |           |
    // |           |          |           |           |           |
    // 0-----------1          0-----4-----1           0-----4-----1

    param GMESH_QUA_4    : int =  3;
    param GMESH_QUA_9    : int = 10;
    param GMESH_QUA_8    : int = 16;
    param GMESH_QUA_16   : int = 36;
    param GMESH_QUA_25   : int = 37;
    param GMESH_QUA_36   : int = 38;
    param GMESH_QUA_12   : int = 39;
    param GMESH_QUA_16I  : int = 40;
    param GMESH_QUA_20   : int = 41;
    param GMESH_QUA_49   : int = 47;
    param GMESH_QUA_64   : int = 48;
    param GMESH_QUA_81   : int = 49;
    param GMESH_QUA_100  : int = 50;
    param GMESH_QUA_121  : int = 51;
    param GMESH_QUA_24   : int = 57;
    param GMESH_QUA_28   : int = 58;
    param GMESH_QUA_32   : int = 59;
    param GMESH_QUA_36I  : int = 60;
    param GMESH_QUA_40   : int = 61;
    param GMESH_QUA_1    : int = 86;

    // Tetrahedron:                          Tetrahedron10:
    //
    //                    v
    //                  .
    //                ,/
    //               /
    //            2                                     2
    //          ,/|`\                                 ,/|`\
    //        ,/  |  `\                             ,/  |  `\
    //      ,/    '.   `\                         ,6    '.   `5
    //    ,/       |     `\                     ,/       8     `\
    //  ,/         |       `\                 ,/         |       `\
    // 0-----------'.--------1 --> u         0--------4--'.--------1
    //  `\.         |      ,/                 `\.         |      ,/
    //     `\.      |    ,/                      `\.      |    ,9
    //        `\.   '. ,/                           `7.   '. ,/
    //           `\. |/                                `\. |/
    //              `3                                    `3
    //                 `\.
    //                    ` w

    param GMESH_TET_4    : int =   4;
    param GMESH_TET_10   : int =  11;
    param GMESH_TET_20   : int =  29;
    param GMESH_TET_56   : int =  31;
    param GMESH_TET_22   : int =  32;
    param GMESH_TET_28   : int =  33;
    param GMESH_TET_84   : int =  71;
    param GMESH_TET_120  : int =  72;
    param GMESH_TET_165  : int =  73;
    param GMESH_TET_220  : int =  74;
    param GMESH_TET_286  : int =  75;
    param GMESH_TET_34   : int =  79;
    param GMESH_TET_35   : int =  30;
    param GMESH_TET_40   : int =  80;
    param GMESH_TET_46   : int =  81;
    param GMESH_TET_52   : int =  82;
    param GMESH_TET_58   : int =  83;
    param GMESH_TET_1    : int =  87;
    param GMESH_TET_SUB  : int = 136;
    param GMESH_TET_16   : int = 137;
    param GMESH_TET_MINI : int = 139;

    // Pyramid:                     Pyramid13:                   Pyramid14:
    //
    //                4                            4                            4
    //              ,/|\                         ,/|\                         ,/|\
    //            ,/ .'|\                      ,/ .'|\                      ,/ .'|\
    //          ,/   | | \                   ,/   | | \                   ,/   | | \
    //        ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
    //      ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
    //    ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
    //  ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
    // 0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
    //  `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
    //    `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
    //      `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \
    //        `\.'      `\     `\`          `\.'             `\`         `\.'             `\`
    //           1----------------2            1--------8-------2           1--------8-------2
    //                     `\
    //                       u

    param GMESH_PYR_5    : int =   7;
    param GMESH_PYR_14   : int =  14;
    param GMESH_PYR_13   : int =  19;
    param GMESH_PYR_30   : int = 118;
    param GMESH_PYR_55   : int = 119;
    param GMESH_PYR_91   : int = 120;
    param GMESH_PYR_140  : int = 121;
    param GMESH_PYR_204  : int = 122;
    param GMESH_PYR_285  : int = 123;
    param GMESH_PYR_385  : int = 124;
    param GMESH_PYR_29   : int = 126;
    param GMESH_PYR_37   : int = 127;
    param GMESH_PYR_21   : int = 125;
    param GMESH_PYR_45   : int = 128;
    param GMESH_PYR_53   : int = 129;
    param GMESH_PYR_61   : int = 130;
    param GMESH_PYR_69   : int = 131;
    param GMESH_PYR_1    : int = 132;

    // Prism:                      Prism15:               Prism18:
    //
    //            w
    //            ^
    //            |
    //            3                       3                      3
    //          ,/|`\                   ,/|`\                  ,/|`\
    //        ,/  |  `\               12  |  13              12  |  13
    //      ,/    |    `\           ,/    |    `\          ,/    |    `\
    //     4------+------5         4------14-----5        4------14-----5
    //     |      |      |         |      8      |        |      8      |
    //     |    ,/|`\    |         |      |      |        |    ,/|`\    |
    //     |  ,/  |  `\  |         |      |      |        |  15  |  16  |
    //     |,/    |    `\|         |      |      |        |,/    |    `\|
    //    ,|      |      |\        10     |      11       10-----17-----11
    //  ,/ |      0      | `\      |      0      |        |      0      |
    // u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
    //     |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
    //     |,/         `\|         |,/         `\|        |,/         `\|
    //     1-------------2         1------9------2        1------9------2

    param GMESH_PRI_6    : int =   6;
    param GMESH_PRI_18   : int =  13;
    param GMESH_PRI_15   : int =  18;
    param GMESH_PRI_1    : int =  89;
    param GMESH_PRI_40   : int =  90;
    param GMESH_PRI_75   : int =  91;
    param GMESH_PRI_126  : int = 106;
    param GMESH_PRI_196  : int = 107;
    param GMESH_PRI_288  : int = 108;
    param GMESH_PRI_405  : int = 109;
    param GMESH_PRI_550  : int = 110;
    param GMESH_PRI_24   : int = 111;
    param GMESH_PRI_33   : int = 112;
    param GMESH_PRI_42   : int = 113;
    param GMESH_PRI_51   : int = 114;
    param GMESH_PRI_60   : int = 115;
    param GMESH_PRI_69   : int = 116;
    param GMESH_PRI_78   : int = 117;

    // Hexahedron:             Hexahedron20:          Hexahedron27:
    //
    //        v
    // 3----------2            3----13----2           3----13----2
    // |\     ^   |\           |\         |\          |\         |\
    // | \    |   | \          | 15       | 14        |15    24  | 14
    // |  \   |   |  \         9  \       11 \        9  \ 20    11 \
    // |   7------+---6        |   7----19+---6       |   7----19+---6
    // |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
    // 0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
    //  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
    //   \ |     \  \ |         10 |        12|        10 |  21    12|
    //    \|      w  \|           \|         \|          \|         \|
    //     4----------5            4----16----5           4----16----5

    param GMESH_HEX_8    : int =   5;
    param GMESH_HEX_27   : int =  12;
    param GMESH_HEX_20   : int =  17;
    param GMESH_HEX_1    : int =  88;
    param GMESH_HEX_64   : int =  92;
    param GMESH_HEX_125  : int =  93;
    param GMESH_HEX_216  : int =  94;
    param GMESH_HEX_343  : int =  95;
    param GMESH_HEX_512  : int =  96;
    param GMESH_HEX_729  : int =  97;
    param GMESH_HEX_1000 : int =  98;
    param GMESH_HEX_32   : int =  99;
    param GMESH_HEX_44   : int = 100;
    param GMESH_HEX_56   : int = 101;
    param GMESH_HEX_68   : int = 102;
    param GMESH_HEX_80   : int = 103;
    param GMESH_HEX_92   : int = 104;
    param GMESH_HEX_104  : int = 105;
  }

  prototype module CGNS
  {
    // CGNS Element Types   -   https://cgns.github.io/CGNS_docs_current/sids/conv.html
    //
    //Dimension | Shape         | Linear    | Quadratic           | Cubic                        | Quartic
    //---------------------------------------------------------------------------------------------------------------------------
    //  0-D     |   Point       |   NODE    |   NODE              |   NODE                       |   NODE
    //  1-D     |   Line        |   BAR_2   |   BAR_3             |   BAR_4                      |   BAR_5
    //  2-D     |   Triangle    |   TRI_3   |   TRI_6             |   TRI_9    TRI_10            |   TRI_12     TRI_15
    //          |   Quadrangle  |   QUAD_4  |   QUAD_8   QUAD_9   |   QUAD_12  QUAD_16           |   QUAD_P4_16 QUAD_25
    //  3-D     |   Tetrahedron |   TETRA_4 |   TETRA_10          |   TETRA_16 TETRA_20          |   TETRA_22   TETRA_34 TETRA_35
    //          |   Pyramid     |   PYRA_5  |   PYRA_13  PYRA_14  |   PYRA_21  PYRA_29  PYRA_30  |   PYRA_P4_29 PYRA_50  PYRA_55
    //          |   Pentahedron |   PENTA_6 |   PENTA_15 PENTA_18 |   PENTA_24 PENTA_38 PENTA_40 |   PENTA_33   PENTA_66 PENTA_75
    //          |   Hexahedron  |   HEXA_8  |   HEXA_20  HEXA_27  |   HEXA_32  HEXA_56  HEXA_64  |   HEXA_44    HEXA_98  HEXA_125
    param CGNS_NODE       : int;

    param CGNS_BAR_2      : int;
    param CGNS_BAR_3      : int;
    param CGNS_BAR_4      : int;
    param CGNS_BAR_5      : int;

    param CGNS_TRI_3      : int;
    param CGNS_TRI_6      : int;
    param CGNS_TRI_9      : int;
    param CGNS_TRI_10     : int;
    param CGNS_TRI_12     : int;
    param CGNS_TRI_15     : int;

    param CGNS_QUAD_4     : int;
    param CGNS_QUAD_8     : int;
    param CGNS_QUAD_9     : int;
    param CGNS_QUAD_12    : int;
    param CGNS_QUAD_16    : int;
    param CGNS_QUAD_P4_16 : int;
    param CGNS_QUAD_25    : int;

    param CGNS_TETRA_4    : int;
    param CGNS_TETRA_10   : int;
    param CGNS_TETRA_16   : int;
    param CGNS_TETRA_20   : int;
    param CGNS_TETRA_22   : int;
    param CGNS_TETRA_34   : int;
    param CGNS_TETRA_35   : int;

    param CGNS_PYRA_5     : int;
    param CGNS_PYRA_13    : int;
    param CGNS_PYRA_14    : int;
    param CGNS_PYRA_21    : int;
    param CGNS_PYRA_29    : int;
    param CGNS_PYRA_30    : int;
    param CGNS_PYRA_P4_29 : int;
    param CGNS_PYRA_50    : int;
    param CGNS_PYRA_55    : int;

    param CGNS_PENTA_6    : int;
    param CGNS_PENTA_15   : int;
    param CGNS_PENTA_18   : int;
    param CGNS_PENTA_24   : int;
    param CGNS_PENTA_38   : int;
    param CGNS_PENTA_40   : int;
    param CGNS_PENTA_33   : int;
    param CGNS_PENTA_66   : int;
    param CGNS_PENTA_75   : int;

    param CGNS_HEXA_8     : int;
    param CGNS_HEXA_20    : int;
    param CGNS_HEXA_27    : int;
    param CGNS_HEXA_32    : int;
    param CGNS_HEXA_56    : int;
    param CGNS_HEXA_64    : int;
    param CGNS_HEXA_44    : int;
    param CGNS_HEXA_98    : int;
    param CGNS_HEXA_125   : int;
  }
}
