/* Documentation for FREI */
prototype module FREI
{
  //Runtime constants
  config const inputFile : string = "input.toml";

  proc main() {
    use Parameters.ParamInput;
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
    use FR;
    use LinearAlgebra;
    use SourceTerm;
    use Temporal_Methods;

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
    var gmesh2 = new unmanaged gmesh2_c(nNodes=Input.nCells+1, nElements=Input.nCells+2, nFamilies=Input.nFaml);

    select Input.meshingScheme
    {
      when MESH_UNIFORM do
        gmesh2.uniform1D(Input.xMin, Input.xMax);
      when MESH_RANDOM do
        gmesh2.random1D(Input.xMin, Input.xMax);
    }

    // 5. Convert input mesh to solver mesh
    var frMesh = new unmanaged fr_mesh_c(nDims=1, nVars=Input.nEqs, solOrder=Input.iOrder);
    frMesh.import_gmesh2(gmesh2);   // Convert mesh to native format
    frMesh.set_families(famlList);  // Get families data from input file and write to mesh
    frMesh.allocate_fr_vars();      // Allocate SP and FP solution/flux/residue arrays
    frMesh.set_points_locations();  // Calculate coordinate trasnformations and point coordinates

    // 4. Initialize the solver, pre calculate stuff
    init_sp2fpInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2spDeriv(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_correction(Input.minOrder+1, Input.maxOrder+1, frMesh.cellTopos);

    // Save mesh file in internal format

    // Initialize solution
    frMesh.solSP = initial_condition(Input.initialCondition, frMesh.xyzSP);

    // Save restart file

    // Output initial state
    iterOutput(iteration, frMesh);

    writeln("Start Iterating");

    // Solve flow
    for iteration in 1..Input.maxIter
    {
      // Zero out residue
      frMesh.resSP = 0.0;

      // Save initial solution
      frMesh.oldSolSP = frMesh.solSP;

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
            for spIdx in frMesh.resSP.domain.dim(0) do
              frMesh.resSP[spIdx..#1, ..] = -source_term(frMesh.xyzSP[spIdx..#1, ..],
                                                         frMesh.solSP[spIdx..#1, ..],
                                                         Input.eqSet                )
                                            * frMesh.jacSP[spIdx];
          }

          // Component 2: Discontinuous Flux
          {
            // Calculate flux at SPs and it´s divergence
            for cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              var thisCell = frMesh.cellList[cellIdx];
              var cellSPini = frMesh.cellSPidx[cellIdx, 1];
              var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

              // Allocate temporary flux array
              var flxSP : [cellSPini.. #cellSPcnt, 1..frMesh.nVars] real;

              // Calculate fluxes
              for meshSP in cellSPini.. #cellSPcnt do
                select Input.eqSet
                {
                  when EQ_CONVECTION do
                    flxSP[meshSP, ..] = convection_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_INVBURGERS do
                    flxSP[meshSP, ..] = burgers_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_EULER do
                    flxSP[meshSP, ..] = euler_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_QUASI_1D_EULER do
                    flxSP[meshSP, ..] = euler_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                }

              // Interpolate fluxes to FPs
              for cellFace in thisCell.faces.domain
              {
                var faceIdx  = thisCell.faces[cellFace];
                var thisFace = frMesh.faceList[faceIdx];
                var faceSide = thisCell.sides[cellFace];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var cellFP = cellFace;
                  frMesh.flxFP[meshFP, faceSide, ..] = dot(sp2fpInterp[(thisCell.elemTopo(), iOrder)]!.coefs(cellFP, ..),
                                                           flxSP[cellSPini..#cellSPcnt,..]                              );
                }
              }

              // Convert fluxes from physical to computational domain.
              // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determiant
              for meshSP in cellSPini.. #cellSPcnt do
                flxSP[meshSP, ..] = dot(flxSP[meshSP, ..], frMesh.metSP[meshSP, 1, 1]**(-1))*frMesh.jacSP[meshSP];

              // Calculate flux divergence
              for cellSP in 1..cellSPcnt
              {
                var meshSP = cellSPini + cellSP - 1;
                frMesh.resSP[meshSP,..] += dot(sp2spDeriv[(thisCell.elemTopo(), iOrder)]!.coefs(cellSP, ..),
                                               flxSP[cellSPini..#cellSPcnt,..]                             );
              }
            }
          }

          // Component 3: Continuous Flux
          {
            // Interpolate solution to FPs
            for cellIdx in frMesh.cellList.domain
            {
              var thisCell = frMesh.cellList[cellIdx];
              var cellSPini = frMesh.cellSPidx[cellIdx, 1];
              var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

              for cellFace in thisCell.faces.domain
              {
                var faceIdx  = thisCell.faces[cellFace];
                var thisFace = frMesh.faceList[faceIdx];
                var faceSide = thisCell.sides[cellFace];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var cellFP = cellFace;
                  frMesh.solFP[meshFP, faceSide, ..] = dot(sp2fpInterp[(thisCell.elemTopo(), iOrder)]!.coefs(cellFP, ..),
                                                           frMesh.solSP[cellSPini..#cellSPcnt,..]                       );
                }
              }
            }

            // Apply boundary conditions
            for faceIdx in frMesh.faceList.domain
            {
              // Get loop variables
              var thisFace = frMesh.faceList[faceIdx];
              var faceFPini = frMesh.faceFPidx[faceIdx, 1];
              var faceFPcnt = frMesh.faceFPidx[faceIdx, 2];

              // Check if the face´s right neighbor is a Boundary Condition
              if frMesh.faceList[faceIdx].cells[2] < 0
              {
                // Yep, it is, lets get some local iteration variables
                var bocoIdx = -frMesh.faceList[faceIdx].cells[2];
                var thisBoco = frMesh.bocoList[bocoIdx];

                var famlIdx = thisBoco.family;
                var thisFaml = frMesh.famlList[famlIdx];

                // Iterate through the FPs on this face
                for meshFP in faceFPini.. #faceFPcnt
                {
                  // Calculate the Boundary Condition using the solution at the left neighbor´s corresponding FP
                  frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml);
                }
              }
            }

            // Calculate interface correction
            for cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              var thisCell = frMesh.cellList[cellIdx];
              var cellSPini = frMesh.cellSPidx[cellIdx, 1];
              var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

              for cellFace in thisCell.faces.domain
              {
                var faceIdx  = thisCell.faces[cellFace];
                var thisFace = frMesh.faceList[faceIdx];
                var faceSide = thisCell.sides[cellFace];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  // For 1D each face has 1 FP therefore the FP and the Face have the same index relative to it's
                  // position in the cell
                  var cellFP = cellFace;

                  // Calculate the flux jump = -1*(local_flux) + numerical_flux
                  var jump : [1..frMesh.nVars] real = -frMesh.flxFP[meshFP, faceSide, ..];
                  select Input.eqSet
                  {
                    when EQ_CONVECTION do
                      jump += upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, 1]);
                    when EQ_INVBURGERS do
                      jump += upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, 1]);
                    when EQ_QUASI_1D_EULER do
                      jump += roe_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..]);
                    when EQ_EULER do
                      jump += roe_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..]);
                  }

                  // Convert fluxes from physical to computational domain.
                  // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determiant
                  jump[..] = dot(jump[..], frMesh.metFP[meshFP, faceSide, 1, 1]**(-1))*frMesh.jacFP[meshFP, faceSide];

                  // The correction function was calculated in the computational domain already, therefore no
                  // transformation is required.
                  frMesh.resSP[cellSPini.. #cellSPcnt, ..] += outer(
                      flux_correction[(thisCell.elemTopo(), iOrder+1)]!.correction[cellFP, 1..cellSPcnt],
                      jump[..]);
                }
              }
            }
          }
        }

        // Advance RK Stage
        {
          // Loop through cells
          for cellIdx in frMesh.cellList.domain
          {
            var cellSPini = frMesh.cellSPidx[cellIdx, 1];
            var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

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
          }
        }
      }

      // Print solver status / log

      // Save restart file

      // Output solution state
      if iteration % ioIter == 0 then
        iterOutput(iteration, frMesh);

      // Check input changes
    }

    // Output the final solution
    //iterOutput(iteration, frMesh);

    writeln();
    writeln("Fin.");
  }
}
