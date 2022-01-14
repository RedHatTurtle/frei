/* Documentation for FREI */
prototype module FREI
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
    use Flux;
    use Riemann;
    use Interpolation;
    use Output;
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
    use Temporal_Methods;

    // Declare timing variables
    const timeUnit = TimeUnits.milliseconds;
    var initTime : real = 0.0;
    var iterTime : real = 0.0;
    var residueTime : real = 0.0;
    var timeStepTime : real = 0.0;
    var stabilizeTime : real = 0.0;
    var iterTimer  : Timer;
    var totalTimer : Timer;
    totalTimer.start();

    var stopwatch  : Timer;
    var srcTermTime : real = 0.0;
    var dscFluxTime : real = 0.0;
    var cntFluxTime : real = 0.0;
    stopwatch.start();

    var dscFluxWatch : Timer;
    var dscFluxTime1 : real = 0;
    var dscFluxTime2 : real = 0;
    var dscFluxTime3 : real = 0;
    var dscFluxTime4 : real = 0;
    dscFluxWatch.start();

    var cntFluxWatch : Timer;
    var cntFluxTime1 : real = 0;
    var cntFluxTime2 : real = 0;
    var cntFluxTime3 : real = 0;
    var cntFluxTime4 : real = 0;
    cntFluxWatch.start();

    var iteration : int = 0;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();
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

    // 6. Save mesh file in internal format

    // 7. Initialize the FR solver, pre calculate coefficients and stuff
    init_sp2fpInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2spDeriv(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2nodeInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_quadratureWeights(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_polyProj(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_correction(Input.minOrder+1, Input.maxOrder+1, frMesh.cellTopos);

    // 8. Initialize solution
    for cellIdx in frMesh.cellList.domain
    {
      ref familyIdx = frMesh.cellList[cellIdx].family;

      ref familyType = frMesh.famlList[familyIdx].bocoType;
      ref familySubType = frMesh.famlList[familyIdx].bocoSubType;
      ref familyParameters = frMesh.famlList[familyIdx].bocoProperties;

      ref cellSPini = frMesh.cellSPidx[cellIdx, 1];
      ref cellSPcnt = frMesh.cellSPidx[cellIdx, 2];
      ref thisCell = frMesh.cellList[cellIdx];

      frMesh.solSP[cellSPini.. #cellSPcnt, ..] = flow_condition(familySubType,
                                                                familyParameters,
                                                                frMesh.xyzSP[cellSPini.. #cellSPcnt, ..]);

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

          frMesh.solFP[meshFP, faceSide, ..] = dot(sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, ..],
                                                   frMesh.solSP[cellSPini..#cellSPcnt,..]                       );
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
            frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml,
                                                            frMesh.xyzFP[meshFP, ..], frMesh.nrmFP[meshFP, ..]);
          }
        }
      }
    }

    // 9. Output initial solution
    iterOutput(iteration, frMesh);

    // 10. Stabilize initial solution
    {
      // Loop through cells
      for cellIdx in frMesh.cellList.domain
      {
        ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
        ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

        for varIdx in 1..frMesh.nVars
        {
          var stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                        cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                        solDegree = Input.iOrder);

          if stableDegree < Input.iOrder then
            frMesh.solSP[cellSPini.. #cellSPcnt, varIdx] = projection_limiter(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                                              cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                                              solDegree = Input.iOrder,
                                                                              projDegree = stableDegree);
        }
      }
    }

    // 11. Save first restart file

    // 12. Initialize convergence monitoring variables
    var l2DeltaIni           : [1..frMesh.nVars] real;
    var l2RelativeDeltaIni   : [1..frMesh.nVars] real;
    var convergenceLog : file;
    try {
      convergenceLog = open("convengence.dat" , iomode.cw);
    } catch {
      stdout.writeln("Unknown Error opening convergence log file.");
      stderr.writeln("Unknown Error opening convergence log file.");
    }
    var convergenceLogChan = convergenceLog.writer();

    writeln();
    initTime = stopwatch.elapsed(timeUnit);
    writef("Stopwatch - Init    : %10.2dr ms\n", initTime);
    writef("Start Iterating\n");

    // Main: Solve flow
    iterTimer.start();
    for iteration in 1..Input.maxIter
    {
      iterTimer.clear();

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

          // Component 1: Source Term
          {
            stopwatch.clear();

            if Input.eqSet == EQ_QUASI_1D_EULER then
              for spIdx in frMesh.resSP.domain.dim(0) do
                frMesh.resSP[spIdx..#1, ..] += -source_term(frMesh.xyzSP[spIdx..#1, ..],
                                                           frMesh.solSP[spIdx..#1, ..],
                                                           Input.eqSet                )
                                              * frMesh.jacSP[spIdx];

            srcTermTime += stopwatch.elapsed(timeUnit);
          }

          // Component 2: Discontinuous Flux
          {
            stopwatch.clear();

            // Calculate flux at SPs and it's divergence
            for cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              ref thisCell = frMesh.cellList[cellIdx];

              // Allocate temporary flux array
              var flxSP : [cellSPini.. #cellSPcnt, 1..frMesh.nDims, 1..frMesh.nVars] real;

              // Step 1: Calculate fluxes
              dscFluxWatch.clear();
              for meshSP in cellSPini.. #cellSPcnt do
                select Input.eqSet
                {
                  when EQ_CONVECTION do
                    flxSP[meshSP, 1, ..] = convection_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_INVBURGERS do
                    flxSP[meshSP, 1, ..] = burgers_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_QUASI_1D_EULER do
                    flxSP[meshSP, 1, ..] = euler_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_EULER do
                    flxSP[meshSP, .., ..] = euler_flux_cv(frMesh.solSP[meshSP, ..]);
                }
              dscFluxTime1 += dscFluxWatch.elapsed(timeUnit);

              // Step 2: Interpolate fluxes to FPs and save the FP normal flux
              dscFluxWatch.clear();
              for cellFace in thisCell.faces.domain
              {
                // Get loop variables
                ref faceIdx  : int = thisCell.faces[cellFace];
                ref faceSide : int = thisCell.sides[cellFace];
                ref thisFace = frMesh.faceList[faceIdx];

                // Allocate temporary flux array
                var flx : [1..frMesh.faceFPidx[faceIdx, 2], 1..2, 1..frMesh.nVars] real;

                // Iterate though all FPs on this face
                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var cellFP : int;
                  if faceSide == 1 then
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                  else
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

                  var uniNrm : [frMesh.nrmFP[meshFP, ..].domain] real = frMesh.nrmFP[meshFP, ..]/norm(frMesh.nrmFP[meshFP, ..]);

                  for varIdx in 1..frMesh.nVars
                  {
                    flx[meshFP+1-frMesh.faceFPidx[faceIdx, 1], 1, varIdx] =
                        dot( sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, ..],
                             flxSP[cellSPini..#cellSPcnt, 1, varIdx]                               );
                    flx[meshFP+1-frMesh.faceFPidx[faceIdx, 1], 2, varIdx] =
                        dot( sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, ..],
                             flxSP[cellSPini..#cellSPcnt, 2, varIdx]                   );

                    frMesh.flxFP[meshFP, faceSide, varIdx] = dot(uniNrm, flx[meshFP+1-frMesh.faceFPidx[faceIdx, 1], .., varIdx]);
                  }
                }
              }
              dscFluxTime2 += dscFluxWatch.elapsed(timeUnit);

              // Step 3: Convert fluxes from physical to computational domain
              dscFluxWatch.clear();
              for meshSP in cellSPini.. #cellSPcnt
              {
                // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determinant
                var jInv : [frMesh.metSP[meshSP, .., ..].domain] real = frMesh.metSP[meshSP, .., ..]**(-1);

                if frMesh.nDims == 2
                {
                  jInv[1,1] =  frMesh.metSP[meshSP, 2, 2];
                  jInv[1,2] = -frMesh.metSP[meshSP, 1, 2];
                  jInv[2,1] = -frMesh.metSP[meshSP, 2, 1];
                  jInv[2,2] =  frMesh.metSP[meshSP, 1, 1];
                }

                flxSP[meshSP, .., ..] = dot(jInv, reshape(flxSP[meshSP, .., ..], flxSP[meshSP, .., ..].domain));
              }
              dscFluxTime3 += dscFluxWatch.elapsed(timeUnit);

              // Step 4: Calculate flux divergence
              dscFluxWatch.clear();
              for cellSP in 1..cellSPcnt
              {
                var meshSP = cellSPini + cellSP - 1;

                for dimIdx in 1..frMesh.nDims
                {
                  var coefs : [sp2spDeriv[(thisCell.elemTopo(), Input.iOrder)]!.coefs[cellSP, dimIdx, ..].domain] real
                             = sp2spDeriv[(thisCell.elemTopo(), Input.iOrder)]!.coefs[cellSP, dimIdx, ..];
                  var flxsp : [flxSP[cellSPini..#cellSPcnt, dimIdx, ..].domain] real
                             = flxSP[cellSPini..#cellSPcnt, dimIdx, ..];

                  frMesh.resSP[meshSP, ..] += dot(coefs, flxsp);
                }
              }
              dscFluxTime4 += dscFluxWatch.elapsed(timeUnit);
            }
            dscFluxTime += stopwatch.elapsed(timeUnit);
          }

          // Component 3: Continuous Flux
          {
            stopwatch.clear();

            // Interpolate solution to FPs
            cntFluxWatch.clear();
            for cellIdx in frMesh.cellList.domain
            {
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              ref thisCell = frMesh.cellList[cellIdx];

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

                  frMesh.solFP[meshFP, faceSide, ..] = dot(sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, ..],
                                                           frMesh.solSP[cellSPini..#cellSPcnt,..]                       );
                }
              }
            }
            cntFluxTime1 += cntFluxWatch.elapsed(timeUnit);

            // Apply boundary conditions
            cntFluxWatch.clear();
            for faceIdx in frMesh.faceList.domain
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
                for meshFP in faceFPini.. #faceFPcnt
                {
                  // Calculate the boundary condition using the solution at the left neighbor´s corresponding FP
                  frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml,
                                                                  frMesh.xyzFP[meshFP, ..], frMesh.nrmFP[meshFP, ..]);
                }
              }
            }
            cntFluxTime2 += cntFluxWatch.elapsed(timeUnit);

            // Calculate interface correction
            cntFluxWatch.clear();
            for cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              ref thisCell = frMesh.cellList[cellIdx];

              for cellFace in thisCell.faces.domain
              {
                ref faceIdx  : int = thisCell.faces[cellFace];
                ref faceSide : int = thisCell.sides[cellFace];
                ref thisFace = frMesh.faceList[faceIdx];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  // For 1D each face has 1 FP therefore the FP and the Face have the same index Relative to it's
                  // position in the cell
                  var cellFP : int;
                  var faceFP : int = meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                  if faceSide == 1 then
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                  else
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

                  // Calculate the flux jump = -1*(local_flux) + numerical_flux
                  var jump : [1..frMesh.nVars] real = -frMesh.flxFP[meshFP, faceSide, ..];
                  select Input.eqSet
                  {
                    when EQ_CONVECTION do
                      jump += upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_INVBURGERS do
                      jump += upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_QUASI_1D_EULER do
                      jump += roe_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_EULER do
                      jump += roe(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                  }

                  // Convert fluxes from physical to computational domain.
                  // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determinant
                  jump[..] = jump[..] * norm(frMesh.nrmFP[meshFP, ..], normType.norm2);

                  select thisCell.elemTopo()
                  {
                    when TOPO_LINE
                    {
                      if faceSide == 1 && cellFace == 1 then
                        jump = -jump;
                      if faceSide == 2 && cellFace == 2 then
                        jump = -jump;
                    }
                    when TOPO_QUAD
                    {
                      if faceSide == 1 && (cellFace == 1 || cellFace == 4) then
                        jump = -jump;
                      if faceSide == 2 && (cellFace == 2 || cellFace == 3) then
                        jump = -jump;
                    }
                  }

                  // The correction function was calculated in the computational domain already, therefore no
                  // transformation is required.
                  frMesh.resSP[cellSPini.. #cellSPcnt, ..] += outer(
                      flux_correction[(thisCell.elemTopo(), Input.iOrder+1)]!.correction[cellFP, 1..cellSPcnt],
                      jump[..]);
                }
              }
            }
            cntFluxTime3 += cntFluxWatch.elapsed(timeUnit);

            cntFluxTime += stopwatch.elapsed(timeUnit);
          }

          residueTime += iterTimer.elapsed(timeUnit);
        }

        // Advance RK Stage
        {
          stopwatch.clear();

          // Loop through cells
          for cellIdx in frMesh.cellList.domain
          {
            ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
            ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

            // Calculate dt for this cell
            var dt : real = Input.timeStep;
            //if variableTimeStep then
            //  dt = time_step();

            // Convert the residual from the computational to the physical domain
            for meshSP in cellSPini.. #cellSPcnt do
              frMesh.resSP[meshSP, ..] /= frMesh.jacSP[meshSP];


            // Update solution
            frMesh.solSP[cellSPini.. #cellSPcnt, ..] = time_advance(frMesh.oldSolSP[cellSPini.. #cellSPcnt, ..],
                                                                    frMesh.solSP[cellSPini.. #cellSPcnt, ..],
                                                                    frMesh.resSP[cellSPini.. #cellSPcnt, ..],
                                                                    dt, stage, Input.timeScheme);

            // Zero out residue
            frMesh.resSP = 0.0;
          }

          timeStepTime += stopwatch.elapsed(timeUnit);
        }

        // Stabilize Solution
        {
          stopwatch.clear();

          if Input.limiterScheme != LIMITER_NONE
          {
            // Loop through cells
            for cellIdx in frMesh.cellList.domain
            {
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

              for varIdx in 1..frMesh.nVars
              {
                var stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                              cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                              solDegree = Input.iOrder);

                if stableDegree < Input.iOrder then
                  frMesh.solSP[cellSPini.. #cellSPcnt, varIdx] = projection_limiter(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                                                    cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                                                    solDegree = Input.iOrder,
                                                                                    projDegree = stableDegree);
              }
            }
          }

          stabilizeTime += stopwatch.elapsed(timeUnit);
        }
      }

      // Print solver status / log

      // Save restart file

      // Calculate and print convergence metrics
      {
        var l1Delta              : [1..frMesh.nVars] real;
        var l2Delta              : [1..frMesh.nVars] real;
        var lInfDelta            : [1..frMesh.nVars] real;
        var l1RelativeDelta      : [1..frMesh.nVars] real;
        var l2RelativeDelta      : [1..frMesh.nVars] real;
        var lInfRelativeDelta    : [1..frMesh.nVars] real;

        // Calculate solution delta from previous iteration
        for varIdx in 1..frMesh.nVars
        {
          l1Delta[varIdx]           =      + reduce     (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx]);
          l2Delta[varIdx]           = sqrt(+ reduce     (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx])**2);
          lInfDelta[varIdx]         = max    reduce  abs(frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx]);
          l1RelativeDelta[varIdx]   =      + reduce     (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx]);
          l2RelativeDelta[varIdx]   = sqrt(+ reduce     (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx])**2);
          lInfRelativeDelta[varIdx] = max    reduce abs((frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx])
                                                         /frMesh.oldSolSP[.., varIdx]);
        }

        // Save values from first iterations as reference
        if iteration == 1
        {
          l2DeltaIni         = l2Delta;
          l2RelativeDeltaIni = l2RelativeDelta;
        }

        // Output summarized convergence metrics to stdOut
        writef("Iteration %9i | Time %{ 10.2dr}ms | Log10(L2(ΔSol)/L2(ΔSol0)) = %{ 7.4dr}", iteration,
            iterTimer.elapsed(timeUnit), log10(norm(l2Delta)/norm(l2DeltaIni)));

        // Output full state to log file
        log_convergence(convergenceLogChan, iteration, l1Delta, l2Delta, lInfDelta, l1RelativeDelta, l2RelativeDelta, lInfRelativeDelta);

        if iteration % ioIter == 0 then
          writef(" | Saving solution file\n");
        else
          writef("\n");
      }

      // Check if we should write the solution this iteration
      if iteration % ioIter == 0 then
        iterOutput(iteration, frMesh);

      // Check if input file changed
    }

    // Output the final solution
    //iterOutput(iteration, frMesh);

    var totalTime : real = totalTimer.elapsed(timeUnit);
    writeln();
    writef("Time splits:\n");
    writef("- Init      : %11.2dr ms - %4.1dr%% of Run-Time\n",      initTime,      initTime/totalTime*100);
    writef("- Residue   : %11.2dr ms - %4.1dr%% of Run-Time\n",   residueTime,   residueTime/totalTime*100);
    writef("  - Src Term: %11.2dr ms - %4.1dr%% of Run-Time\n",   srcTermTime,   srcTermTime/totalTime*100);
    writef("  - Dsc Flux: %11.2dr ms - %4.1dr%% of Run-Time\n",   dscFluxTime,   dscFluxTime/totalTime*100);
    writef("    - Step 1: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime1,  dscFluxTime1/dscFluxTime*100);
    writef("    - Step 2: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime2,  dscFluxTime2/dscFluxTime*100);
    writef("    - Step 3: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime3,  dscFluxTime3/dscFluxTime*100);
    writef("    - Step 4: %11.2dr ms - %4.1dr%% of Dsc Flux\n",  dscFluxTime4,  dscFluxTime4/dscFluxTime*100);
    writef("  - Cnt Flux: %11.2dr ms - %4.1dr%% of Run-Time\n",   cntFluxTime,   cntFluxTime/totalTime*100);
    writef("    - Step 1: %11.2dr ms - %4.1dr%% of Cnt Flux\n",  cntFluxTime1,  cntFluxTime1/cntFluxTime*100);
    writef("    - Step 2: %11.2dr ms - %4.1dr%% of Cnt Flux\n",  cntFluxTime2,  cntFluxTime2/cntFluxTime*100);
    writef("    - Step 3: %11.2dr ms - %4.1dr%% of Cnt Flux\n",  cntFluxTime3,  cntFluxTime3/cntFluxTime*100);
    writef("- Stabilize : %11.2dr ms - %4.1dr%% of Run-Time\n", stabilizeTime, stabilizeTime/totalTime*100);
    writef("- Time-Step : %11.2dr ms - %4.1dr%% of Run-Time\n",  timeStepTime,  timeStepTime/totalTime*100);
    writef("---------------------------------------------------------\n");
    var sumTime : real = initTime + srcTermTime + dscFluxTime + cntFluxTime + stabilizeTime + timeStepTime;
    writef("  Sum       : %11.2dr ms - %4.1dr%% of Run-Time\n",  sumTime,  sumTime/totalTime*100);
    writef("  Run-time  : %11.2dr ms\n", totalTime);

    writeln();
    writeln("Fin");
  }
}
