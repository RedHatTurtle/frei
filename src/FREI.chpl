/* Documentation for FREI */
module FREI
{
  // Run-time constants
  config const inputFile : string = "input.toml";
  config const inputMesh : string = "mesh.mesh";

  proc main() {
    use IO;
    use Time;
    use Parameters.ParamInput;
    use Parameters.ParamMesh;
    use Config;
    use Input;
    use Dimensional;
    use Flux;
    use Riemann;
    use Interpolation;
    use Output;
    use ErrorCalc;
    use Gmesh;
    use Mesh;
    use FRMesh;
    use Boundary;
    use Correction;
    use Init;
    use Quadrature;
    use Projection;
    use Limiter;
    use FR;
    use LinearAlgebra;
    use SourceTerm;
    use Temporal;

    // Declare timing variables
    var initTime      : real = 0.0;
    var iterTime      : real = 0.0;
    var timeStepTime  : real = 0.0;
    var stabilizeTime : real = 0.0;
    var iterWatch     : stopwatch;
    var totalWatch    : stopwatch;
    totalWatch.start();

    var residueWatch : stopwatch;
    var residueTime  : real = 0.0;
    residueWatch.start();

    var mainWatch  : stopwatch;
    var srcTermTime : real = 0.0;
    var dscFluxTime : real = 0.0;
    var cntFluxTime : real = 0.0;
    mainWatch.start();

    var dscFluxWatch : stopwatch;
    var dscFluxTime1 : real = 0.0;
    var dscFluxTime2 : real = 0.0;
    var dscFluxTime3 : real = 0.0;
    var dscFluxTime4 : real = 0.0;
    dscFluxWatch.start();

    var cntFluxWatch : stopwatch;
    var cntFluxTime1 : real = 0.0;
    var cntFluxTime2 : real = 0.0;
    var cntFluxTime3 : real = 0.0;
    var cntFluxTime4 : real = 0.0;
    cntFluxWatch.start();

    var jumpCorrectionWatch : stopwatch;
    var riemTime : real = 0.0;
    var jumpTime : real = 0.0;
    var corrTime : real = 0.0;
    jumpCorrectionWatch.start();

    var iteration : int = 0;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();
    init_scales(lengRef=Input.lengRef, velMRef=Input.velMRef, tempRef=Input.tempRef, presRef=Input.presRef);
    var timeStepStages : int;
    select timeScheme
    {
      when TIME_TVDRK_O2S2 do
       timeStepStages = 2;
      when TIME_TVDRK_O2S3 do
       timeStepStages = 3;
      when TIME_TVDRK_O2S4 do
       timeStepStages = 4;
      when TIME_TVDRK_O2SN do
       timeStepStages = 5;
      when TIME_TVDRK_O3S3 do
       timeStepStages = 3;
      when TIME_TVDRK_O3S4 do
       timeStepStages = 4;
      when TIME_TVDRK_O3S5 do
       timeStepStages = 5;
      when TIME_TVDRK_O4S5 do
       timeStepStages = 5;
      otherwise do
       timeStepStages = 1;
    }

    // 3. Read / define mesh
    var gmesh2 = new unmanaged gmesh2_c();
    select Input.meshFormat
    {
      when MESH_GENERATE
      {
        select Input.meshingScheme
        {
          when MESH_GEN_UNIFORM do
            gmesh2.uniform1D(Input.nCells, Input.xMin, Input.xMax);
          when MESH_GEN_RANDOM do
            gmesh2.random1D(Input.nCells, Input.xMin, Input.xMax);
        }
      }
      when MESH_GMESH do
        gmesh2.read_gmesh_file(Input.meshFileName);
      when MESH_CGNS {}
    }

    // 4. Convert input mesh to solver mesh
    var frMesh = new unmanaged fr_mesh_c(nDims=gmesh2.mesh_dimension(), nVars=Input.nEqs, solOrder=Input.iOrder);
    frMesh.import_gmesh2(gmesh2);   // Convert mesh to native format
    frMesh.set_families(famlList);  // Get families data from input file and write to mesh

    // 5. Initialize FR mesh
    frMesh.allocate_fr_vars();      // Allocate SP and FP solution/flux/residue arrays
    frMesh.set_points_locations();  // Calculate coordinate transformations and point coordinates
    frMesh.build_cell_char_leng();  // Calculate the characteristic length of the mesh cells

    // 6. Save mesh file in internal format

    // 7. Initialize the FR solver, pre calculate coefficients and stuff
    init_sp2fpInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2spDeriv(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2nodeInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_quadratureWeights(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_polyProj(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_correction(Input.minOrder+1, Input.maxOrder+1, frMesh.cellTopos);

    // 8. Initialize solution
    forall cellIdx in frMesh.cellList.domain
    {
      ref familyIdx = frMesh.cellList[cellIdx].family;

      ref familyType = frMesh.famlList[familyIdx].bocoType;
      ref familySubType = frMesh.famlList[familyIdx].bocoSubType;
      ref familyParameters = frMesh.famlList[familyIdx].bocoProperties;

      ref cellSPini = frMesh.cellSPidx[cellIdx, 1];
      ref cellSPcnt = frMesh.cellSPidx[cellIdx, 2];
      ref thisCell = frMesh.cellList[cellIdx];

      frMesh.solSP[.., cellSPini.. #cellSPcnt] = flow_condition(familySubType                           ,
                                                                familyParameters                        ,
                                                                frMesh.xyzSP[cellSPini.. #cellSPcnt, ..]).T;

      // Interpolate solution to FPs
      for cellFace in thisCell.faces.domain
      {
        ref faceIdx  : int = thisCell.faces[cellFace];
        ref faceSide : int = thisCell.sides[cellFace];
        ref thisFace = frMesh.faceList[faceIdx];

        for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
        {
          var cellFP : int;
          if faceSide == 1 then
            cellFP = (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
          else
            cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

          frMesh.solFP[meshFP, faceSide, ..] = dot(frMesh.solSP[.., cellSPini.. #cellSPcnt]                              ,
                                                   sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, ..]);
        }

        // Check if the face's right neighbor is a Boundary Condition
        if frMesh.faceList[faceIdx].cells[2] < 0
        {
          // Yep, it is, lets get some local iteration variables
          ref faceFPini : int = frMesh.faceFPidx[faceIdx, 1];
          ref faceFPcnt : int = frMesh.faceFPidx[faceIdx, 2];

          ref thisBoco = frMesh.bocoList[-frMesh.faceList[faceIdx].cells[2]];
          ref thisFaml = frMesh.famlList[thisBoco.family];

          // Iterate through the FPs on this face
          for meshFP in faceFPini.. #faceFPcnt
          {
            // Calculate the boundary condition using the solution at the left neighbor´s corresponding FP
            frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml             ,
                                                            frMesh.xyzFP[meshFP, ..], frMesh.nrmFP[meshFP, ..]);
          }
        }
      }
    }

    // 9. Output initial solution
    iterOutput(iteration, frMesh, flagNormals = true);

    // 10. Stabilize initial solution
    if Input.limiterScheme != LIMITER_NONE
    {
      // Loop through cells
      forall cellIdx in frMesh.cellList.domain
      {
        ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
        ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

        for varIdx in 1..frMesh.nVars
        {
          var stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                        jacobian = frMesh.jacSP[cellSPini.. #cellSPcnt]       ,
                                                        cellTopo = frMesh.cellList[cellIdx].elemTopo()        ,
                                                        solDegree = frMesh.solOrder                           );

          if stableDegree < frMesh.solOrder then
            frMesh.solSP[varIdx, cellSPini.. #cellSPcnt] = projection_limiter(solPoly = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                                              cellTopo = frMesh.cellList[cellIdx].elemTopo()        ,
                                                                              solDegree = frMesh.solOrder                           ,
                                                                              projDegree = stableDegree                             );
        }
      }
    }

    // 11. Save first restart file

    // 12. Initialize convergence monitoring variables
    var l2SolDeltaAbsIni : [1..frMesh.nVars] real;
    var convergenceLogFile : file;
    var       errorLogFile : file;
    try! {
      convergenceLogFile = open("convergence.dat", ioMode.cw);
            errorLogFile = open(      "error.dat", ioMode.cw);
    } catch {
      try! stdout.writeln("Unknown Error opening convergence log file.");
      try! stderr.writeln("Unknown Error opening convergence log file.");
    }
    var convergenceLogWriter = try! convergenceLogFile.writer();
    var       errorLogWriter = try!       errorLogFile.writer();

    writeln();
    initTime = mainWatch.elapsed();
    writef("Initialization Time: %10.2dr ms\n", initTime*1000);
    writeln();
    writef("Start Iterating\n");

    // Main: Solve flow
    iterWatch.start();
    for iteration in 1..Input.maxIter
    {
      iterWatch.clear();

      // Save initial solution
      frMesh.oldSolSP = frMesh.solSP;

      // Iterate RK stages
      for stage in 1..timeStepStages
      {
        // Calculate residue for this iteration
        {
          // The residue has 3 components:
          //   1. Continuous Flux
          //   2. Discontinuous Flux
          //   3. Source terms
          //
          // The residual array is reset in the time stepping procedure
          residueWatch.clear();

          // Component 1: Source Term
          {
            mainWatch.clear();

            if Input.eqSet == EQ_QUASI_1D_EULER then
              forall spIdx in frMesh.resSP.domain.dim(1) do
                frMesh.resSP[.., spIdx.. #1] += -source_term(frMesh.xyzSP[spIdx.. #1, ..],
                                                             frMesh.solSP[.., spIdx.. #1],
                                                             Input.eqSet                 )
                                               * frMesh.jacSP[spIdx];

            srcTermTime += mainWatch.elapsed();
          }

          // Component 2: Discontinuous Flux
          {
            mainWatch.clear();

            // Calculate flux at SPs and it's divergence
            forall cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              ref thisCell = frMesh.cellList[cellIdx];
              var cellTopo : int = thisCell.elemTopo();

              // Allocate temporary flux array
              var flxSP : [1..frMesh.nDims, 1..frMesh.nVars, 1..cellSPcnt] real;

              // Step 1: Calculate fluxes at SPs
              //dscFluxWatch.clear();
              for meshSP in cellSPini.. #cellSPcnt do
                select Input.eqSet
                {
                  when EQ_CONVECTION do
                    flxSP[ 1, .., meshSP-cellSPini+1] = convection_flux_cv_1d(frMesh.solSP[.., meshSP], Input.convectionSpeed);
                  when EQ_INVBURGERS do
                    flxSP[ 1, .., meshSP-cellSPini+1] = burgers_flux_cv_1d(frMesh.solSP[.., meshSP]);
                  when EQ_QUASI_1D_EULER do
                    flxSP[ 1, .., meshSP-cellSPini+1] = euler_flux_cv_1d(frMesh.solSP[.., meshSP], Input.fGamma);
                  when EQ_EULER do
                    flxSP[.., .., meshSP-cellSPini+1] = euler_flux_cv(frMesh.solSP[.., meshSP], Input.fGamma);
                }
              //dscFluxTime1 += dscFluxWatch.elapsed();

              // Step 2: Interpolate fluxes to FPs and save the FP normal flux
              //dscFluxWatch.clear();
              for cellFace in thisCell.faces.domain
              {
                // Get loop variables
                ref faceIdx  : int = thisCell.faces[cellFace];
                ref faceSide : int = thisCell.sides[cellFace];
                ref thisFace = frMesh.faceList[faceIdx];

                // Iterate though all FPs on this face
                for meshFP in frMesh.faceFPidx[faceIdx, 1].. #frMesh.faceFPidx[faceIdx, 2]
                {
                  // Allocate temporary flux array
                  var flx : [1..2] real;

                  var faceFP = meshFP+1-frMesh.faceFPidx[faceIdx, 1];
                  var uniNrm = frMesh.nrmFP[meshFP, ..]/norm(frMesh.nrmFP[meshFP, ..]);

                  var cellFP : int;
                  if faceSide == 1 then
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                  else
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

                  for varIdx in 1..frMesh.nVars
                  {
                    flx[1] = dot( flxSP[ 1, varIdx, 1..cellSPcnt]                            ,
                                  sp2fpInterp[(cellTopo, frMesh.solOrder)]!.coefs[cellFP, ..]);
                    flx[2] = dot( flxSP[ 2, varIdx, 1..cellSPcnt]                            ,
                                  sp2fpInterp[(cellTopo, frMesh.solOrder)]!.coefs[cellFP, ..]);
                    // Doesn't work even though it seems like it should
                    //flx = dot( flxSP[.., varIdx, 1..cellSPcnt]                            ,
                    //           sp2fpInterp[(cellTopo, frMesh.solOrder)]!.coefs[cellFP, ..]);

                    frMesh.flxFP[meshFP, faceSide, varIdx] = dot(flx, uniNrm);
                  }
                }
              }
              //dscFluxTime2 += dscFluxWatch.elapsed();

              // Step 3: Convert fluxes from physical to computational domain
              //dscFluxWatch.clear();
              for meshSP in cellSPini.. #cellSPcnt
              {
                // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determinant
                var flxsp = flxSP[.., .., meshSP-cellSPini+1];
                flxSP[.., .., meshSP-cellSPini+1] = dot( frMesh.metSP[meshSP, .., ..], flxsp)*frMesh.jacSP[meshSP];
              }
              //dscFluxTime3 += dscFluxWatch.elapsed();

              // Step 4: Calculate flux divergence
              //dscFluxWatch.clear();
              for cellSP in 1..cellSPcnt
              {
                var meshSP = cellSPini + cellSP - 1;

                for dimIdx in 1..frMesh.nDims
                {
                  var flxsp = flxSP[dimIdx, .., 1..cellSPcnt];
                  frMesh.resSP[.., meshSP] += dot(flxsp, sp2spDeriv[(cellTopo, frMesh.solOrder)]!.coefs[cellSP, dimIdx, ..]);
                }
              }
              //dscFluxTime4 += dscFluxWatch.elapsed();
            }
            dscFluxTime += mainWatch.elapsed();
          }

          // Component 3: Continuous Flux
          {
            mainWatch.clear();

            // Step 1: Interpolate solution to FPs
            cntFluxWatch.clear();
            forall cellIdx in frMesh.cellList.domain
            {
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              ref thisCell = frMesh.cellList[cellIdx];

              forall cellFace in thisCell.faces.domain
              {
                ref faceIdx  : int = thisCell.faces[cellFace];
                ref faceSide : int = thisCell.sides[cellFace];
                ref thisFace = frMesh.faceList[faceIdx];

                forall meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var cellFP : int;

                  if faceSide == 1 then
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                  else
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

                  frMesh.solFP[meshFP, faceSide, ..] = dot(frMesh.solSP[.., cellSPini.. #cellSPcnt]                              ,
                                                           sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, ..]);
                }
              }
            }
            cntFluxTime1 += cntFluxWatch.elapsed();

            // Step 2: Apply boundary conditions
            cntFluxWatch.clear();
            forall faceIdx in frMesh.faceList.domain
            {
              // Get loop variables
              ref faceFPini : int = frMesh.faceFPidx[faceIdx, 1];
              ref faceFPcnt : int = frMesh.faceFPidx[faceIdx, 2];
              ref thisFace = frMesh.faceList[faceIdx];

              // Check if the face´s right neighbor is a Boundary Condition
              if frMesh.faceList[faceIdx].cells[2] < 0
              {
                // Yep, it is, lets get some local iteration variables
                ref thisBoco = frMesh.bocoList[-frMesh.faceList[faceIdx].cells[2]];
                ref thisFaml = frMesh.famlList[thisBoco.family];

                // Iterate through the FPs on this face
                forall meshFP in faceFPini.. #faceFPcnt
                {
                  // Calculate the boundary condition using the solution at the left neighbor´s corresponding FP
                  frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml             ,
                                                                  frMesh.xyzFP[meshFP, ..], frMesh.nrmFP[meshFP, ..]);
                }
              }
            }
            cntFluxTime2 += cntFluxWatch.elapsed();

            // Step 3: Calculate interface correction
            cntFluxWatch.clear();
            forall cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              ref thisCell = frMesh.cellList[cellIdx];
              var cellTopo : int = thisCell.elemTopo();

              for cellFace in thisCell.faces.domain
              {
                ref faceIdx  : int = thisCell.faces[cellFace];
                ref faceSide : int = thisCell.sides[cellFace];
                ref thisFace = frMesh.faceList[faceIdx];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var jump : [1..frMesh.nVars] real;

                  // Operation 1: Calculate Riemann flux at the FP
                  //jumpCorrectionWatch.clear();
                  select Input.eqSet
                  {
                    when EQ_CONVECTION do
                      jump = upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_INVBURGERS do
                      jump = upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_QUASI_1D_EULER do
                      jump = roe_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_EULER do
                      jump = roe(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                  }
                  //riemTime += jumpCorrectionWatch.elapsed();

                  // Operation 2: Calculate jump at a FP and convert it to the physical domain
                  //jumpCorrectionWatch.clear();
                  {
                    // Calculate the flux jump = -1*(local_flux) + numerical_flux
                    jump -= frMesh.flxFP[meshFP, faceSide, ..];

                    // Convert fluxes from physical to computational domain.
                    // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determinant
                    jump[..] = jump[..] * norm(frMesh.nrmFP[meshFP, ..], normType.norm2);

                    select cellTopo
                    {
                      when TOPO_LINE
                      {
                        if      faceSide == 1 && cellFace == 1 then
                          jump *= -1;
                        else if faceSide == 2 && cellFace == 2 then
                          jump *= -1;
                      }
                      when TOPO_QUAD
                      {
                        if      faceSide == 1 && (cellFace == 1 || cellFace == 4) then
                          jump *= -1;
                        else if faceSide == 2 && (cellFace == 2 || cellFace == 3) then
                          jump *= -1;
                      }
                    }
                  }
                  //jumpTime += jumpCorrectionWatch.elapsed();

                  // Operation 3: Apply the correction to the residue matrix
                  //jumpCorrectionWatch.clear();
                  {
                    // For 1D each face has 1 FP therefore the FP and the Face have the same index Relative to it's
                    // position in the cell
                    var cellFP : int;
                    var faceFP : int = meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                    if faceSide == 1 then
                      cellFP = (cellFace-1)*(frMesh.solOrder+1) +  faceFP;
                    else
                      cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] + 1 - faceFP);

                    // The correction function was calculated in the computational domain already, therefore no
                    // transformation is required.
                    frMesh.resSP[.., cellSPini.. #cellSPcnt] += outer(jump[..]                          ,
                        flux_correction[(cellTopo, frMesh.solOrder+1)]!.correction[cellFP, 1..cellSPcnt]);
                  }
                  //corrTime += jumpCorrectionWatch.elapsed();
                }
              }
            }
            cntFluxTime3 += cntFluxWatch.elapsed();

            cntFluxTime += mainWatch.elapsed();
          }

          residueTime += residueWatch.elapsed();
        }

        // Advance RK Stage
        {
          mainWatch.clear();

          frMesh.calc_time_step();

          // Loop through cells
          forall cellIdx in frMesh.cellList.domain
          {
            ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
            ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

            // Convert the residual from the computational to the physical domain
            forall meshSP in cellSPini.. #cellSPcnt do
              frMesh.resSP[.., meshSP] /= frMesh.jacSP[meshSP];

            // Update solution
            frMesh.solSP[.., cellSPini.. #cellSPcnt] = time_advance(frMesh.oldSolSP[.., cellSPini.. #cellSPcnt],
                                                                    frMesh.solSP[   .., cellSPini.. #cellSPcnt],
                                                                    frMesh.resSP[   .., cellSPini.. #cellSPcnt],
                                                                    frMesh.cellTimeStep[cellIdx],
                                                                    stage, Input.timeScheme                     );
          }

          // Zero out residue
          frMesh.resSP = 0.0;

          timeStepTime += mainWatch.elapsed();
        }

        // Stabilize Solution
        {
          mainWatch.clear();

          if Input.limiterScheme != LIMITER_NONE
          {
            // Loop through cells
            forall cellIdx in frMesh.cellList.domain
            {
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

              for varIdx in 1..frMesh.nVars
              {
                var stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                              jacobian = frMesh.jacSP[cellSPini.. #cellSPcnt]       ,
                                                              cellTopo = frMesh.cellList[cellIdx].elemTopo()        ,
                                                              solDegree = frMesh.solOrder                           );

                if stableDegree < frMesh.solOrder then
                  frMesh.solSP[varIdx, cellSPini.. #cellSPcnt] = projection_limiter(solPoly    = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                                                    cellTopo   = frMesh.cellList[cellIdx].elemTopo()         ,
                                                                                    solDegree  = frMesh.solOrder                             ,
                                                                                    projDegree = stableDegree                                );
              }
            }
          }

          stabilizeTime += mainWatch.elapsed();
        }
      }

      // IO
      {
        // Save restart file

        // Check if we should write the solution this iteration
        if iteration % ioIter == 0 then
          iterOutput(iteration, frMesh);

        // Calculate solution error and write to error log
        if Input.outError > 0
        {
          var errors = error_calc(Input.outError, frMesh);
          print_log(errorLogWriter, iteration, errors);
        }

        // Calculate and print convergence metrics
        {
          // Calculate solution delta from previous iteration
          const l1SolDeltaAbs : [1..frMesh.nVars] real = [varIdx in 1..frMesh.nVars]
                + reduce abs(frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]);
          const l2SolDeltaAbs : [1..frMesh.nVars] real = [varIdx in 1..frMesh.nVars] sqrt(
                + reduce    (frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..])**2);
          const lfSolDeltaAbs : [1..frMesh.nVars] real = [varIdx in 1..frMesh.nVars]
              max reduce abs(frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]);
          const l1SolDeltaRel : [1..frMesh.nVars] real = [varIdx in 1..frMesh.nVars]
                + reduce abs((frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]) / frMesh.oldSolSP[varIdx, ..]);
          const l2SolDeltaRel : [1..frMesh.nVars] real = [varIdx in 1..frMesh.nVars] sqrt(
                + reduce    ((frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]) / frMesh.oldSolSP[varIdx, ..] )**2);
          const lfSolDeltaRel : [1..frMesh.nVars] real = [varIdx in 1..frMesh.nVars]
              max reduce abs((frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]) / frMesh.oldSolSP[varIdx, ..]);

          // Save deltas from first iterations as reference
          if iteration == 1 then
            l2SolDeltaAbsIni = l2SolDeltaAbs;

          // Output full state to log file
          log_convergence(convergenceLogWriter, iteration, l1SolDeltaAbs, l2SolDeltaAbs, lfSolDeltaAbs,
                                                           l1SolDeltaRel, l2SolDeltaRel, lfSolDeltaRel);

          // Output summarized convergence metrics to stdOut
          writef("Iteration %9i | Time %{ 10.2dr}ms | Log10(L2(ΔSol)/L2(ΔSol0)) = %{ 7.4dr}",
              iteration, iterWatch.elapsed()*1000, log10(norm(l2SolDeltaAbs)/norm(l2SolDeltaAbsIni)));

          if iteration % ioIter == 0 then
            writef(" | Solution file saved\n");
          else
            writef("\n");
        }

        // Check if input file changed
        {}

      }
    }

    /////////////////////////////
    // Output and Finalization //
    /////////////////////////////

    // Output the final solution
    iterOutput(iteration, frMesh);

    var totalTime : real = totalWatch.elapsed();
    writeln();
    writef("Time splits:\n");
    writef("- Init      : %11.2dr ms - %4.1dr%% of Run-Time\n",      initTime*1000,      initTime/totalTime   *100);
    writef("- Residue   : %11.2dr ms - %4.1dr%% of Run-Time\n",   residueTime*1000,   residueTime/totalTime   *100);
    writef("  - Src Term: %11.2dr ms - %4.1dr%% of Run-Time\n",   srcTermTime*1000,   srcTermTime/totalTime   *100);
    writef("  - Dsc Flux: %11.2dr ms - %4.1dr%% of Run-Time\n",   dscFluxTime*1000,   dscFluxTime/totalTime   *100);
    //writef("    - Step 1: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime1*1000,  dscFluxTime1/dscFluxTime *100);
    //writef("    - Step 2: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime2*1000,  dscFluxTime2/dscFluxTime *100);
    //writef("    - Step 3: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime3*1000,  dscFluxTime3/dscFluxTime *100);
    //writef("    - Step 4: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime4*1000,  dscFluxTime4/dscFluxTime *100);
    writef("  - Cnt Flux: %11.2dr ms - %4.1dr%% of Run-Time\n",   cntFluxTime*1000,   cntFluxTime/totalTime   *100);
    writef("    - Step 1: %11.2dr ms - %4.1dr%% of Cnt Flux\n",  cntFluxTime1*1000,  cntFluxTime1/cntFluxTime *100);
    writef("    - Step 2: %11.2dr ms - %4.1dr%% of Cnt Flux\n",  cntFluxTime2*1000,  cntFluxTime2/cntFluxTime *100);
    writef("    - Step 3: %11.2dr ms - %4.1dr%% of Cnt Flux\n",  cntFluxTime3*1000,  cntFluxTime3/cntFluxTime *100);
    //writef("      - Op 1: %11.2dr ms - %4.1dr%% of St3 Flux\n",      riemTime*1000,      riemTime/cntFluxTime3*100);
    //writef("      - Op 2: %11.2dr ms - %4.1dr%% of St3 Flux\n",      jumpTime*1000,      jumpTime/cntFluxTime3*100);
    //writef("      - Op 3: %11.2dr ms - %4.1dr%% of St3 Flux\n",      corrTime*1000,      corrTime/cntFluxTime3*100);
    writef("- Stabilize : %11.2dr ms - %4.1dr%% of Run-Time\n", stabilizeTime*1000, stabilizeTime/totalTime   *100);
    writef("- Time-Step : %11.2dr ms - %4.1dr%% of Run-Time\n",  timeStepTime*1000,  timeStepTime/totalTime   *100);
    writef("---------------------------------------------------------\n");
    var sumTime : real = initTime + srcTermTime + dscFluxTime + cntFluxTime + stabilizeTime + timeStepTime;
    writef("  Sum       : %11.2dr ms - %4.1dr%% of Run-Time\n",  sumTime*1000,  sumTime/totalTime*100);
    writef("  Run-time  : %11.2dr ms\n", totalTime*1000);

    writeln();
    writeln("Fin");
  }
}
