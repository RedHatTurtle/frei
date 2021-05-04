/* Documentation for FREI */
module FREI
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
    use FR;
    use LinearAlgebra;

    var iteration : int = 0;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();

    // 3. Read / define mesh
    var gmesh2 = new unmanaged gmesh2_c(nNodes=7, nElements=8, nFamilies=3);
    gmesh2.random1D(-1,1);

    // 5. Convert input mesh to solver mesh
    var frMesh = new unmanaged fr_mesh_c(nDims=1, nVars=3, solOrder=iOrder-1);
    frMesh.import_gmesh2(gmesh2);
    frMesh.allocate_vars();
    frMesh.set_points_locations();

    // 4. Initialize the solver, pre calculate stuff
    init_sp2fpInterp(minOrder, maxOrder, frMesh.cellTopos);

    // Save mesh file in internal format

    // Initialize solution
    frMesh.solSP = initial_condition(IC_SHOCKTUBE, frMesh.xyzSP);

    // Save restart file

    // Output initial state
    iterOutput(iteration, frMesh);

    // Solve flow
    for iteration in 0..maxIter
    {
      // Iterate Solver (single or multiple time steps)
      {
        // Interpolate solution to FPs
        for cellIdx in frMesh.cellList.domain
        {
          var thisCell = frMesh.cellList[cellIdx];
          for cellFace in thisCell.faces.domain
          {
            var faceIdx  = thisCell.faces[cellFace];
            var thisFace = frMesh.faceList[faceIdx];
            var faceSide = thisCell.sides[cellFace];
            for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
            {
              var cellFP = cellFace;
              var cellSPini = frMesh.cellSPidx[cellIdx, 1];
              var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];
              frMesh.solFP[meshFP, faceSide, 1] = dot(sp2fpInterp[(thisCell.elemTopo(), iOrder)]!.coefs(cellFP, ..),
                                                      frMesh.solSP[cellSPini..#cellSPcnt,1]                                           );
            }
          }
        }

        // Calculate numerical flux
        for face in frMesh.faceList.domain
        {
          for faceFP in frMesh.faceFPidx[face,1].. #frMesh.faceFPidx[face,2]
          {
            // Riemann flux
          }
        }

        // Calculate internal correction
        for cell in frMesh.cellList
        {}

        // Calculate interface correction
        for cell in frMesh.cellList
        {}
      }
      // Print solver status / log

      // Save restart file

      // Output solution state
      iterOutput(iteration, frMesh);

      // Check input changes
    }

    // Output the final solution
    iterOutput(iteration, frMesh);
  }
}
