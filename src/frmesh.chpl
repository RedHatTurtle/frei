prototype module FRMesh {
  use Input;
  use Mesh;

  class fr_mesh_c : mesh_c
  {
    const nVars  : int;
    var solOrder : int;

    // Domains
    var  xyzSP_d : domain(rank=2, idxType=int);    // {1..nSPs, 1..nDims}
    var  xyzFP_d : domain(rank=2, idxType=int);    // {1..nFPs, 1..nDims}

    var  metSP_d : domain(rank=3, idxType=int);    // {1..nSPs, 1..nDims, 1..nDims}       nDims+1 if moving mesh
    var  jacSP_d : domain(rank=1, idxType=int);    // {1..nSPs}                           nDims+1 if moving mesh

    var  solSP_d : domain(rank=2, idxType=int);    // {1..nSPs, 1..nVars}
    var  solFP_d : domain(rank=3, idxType=int);    // {1..nFPs, 1..2, 1..nVars}
    var  flxFP_d : domain(rank=3, idxType=int);    // {1..nFPs, 1..2, 1..nVars}
    var dSolSP_d : domain(rank=3, idxType=int);    // {1..nSPs, 1..nVars, 1..nVars}
    var dSolFP_d : domain(rank=4, idxType=int);    // {1..nFPs, 1..2, 1..nVars, 1..nVars}

    var cellSPidx_d : domain(rank=2, idxType=int); // {1..nCells, 1..2}
    var faceFPidx_d : domain(rank=2, idxType=int); // {1..nFaces, 1..2}

    // FR solver variables
    var xyzSP : [xyzSP_d] real;
    var xyzFP : [xyzFP_d] real;

    var metSP : [metSP_d] real;   // First degree metric terms, Jacobian matrix
    var jacSP : [jacSP_d] real;   // Jacobian
    var nrmFP : [xyzFP_d] real;   // Face normal vector pointing from the left (side 1) to the right (side 2) cell

    var oldSolSP : [ solSP_d] real;     // Backup of the solution at the beginning of residue calculation
    var    solSP : [ solSP_d] real;     // Conserved variables at SPs
    var    solFP : [ solFP_d] real;     // Conserved variables at FPs (1-Left / 2-right)
    var    flxFP : [ flxFP_d] real;     // Discontinuous flux at FPs (1-Left / 2-right)
    var   dSolSP : [dSolSP_d] real;     // Gradient, at the SPs, of the discontinuous solution interpolation
    var   dSolFP : [dSolFP_d] real;     // Gradient, at the FPs, of the discontinuous flux reconstruction
    var    resSP : [ solSP_d] real;     // conserved variables residual

    var cellSPidx : [cellSPidx_d] int;  // Index of the first SP and number of SPs of a cell
    var faceFPidx : [faceFPidx_d] int;  // Index of the first FP and number of FPs of a face

    proc allocate_fr_vars()
    {
      var nSPs : int;
      var nFPs : int;

      cellSPidx_d = {this.cellList_d.dim(0), 1..2};
      faceFPidx_d = {this.faceList_d.dim(0), 1..2};

      for cell in this.cellList_d
      {
        var cellSPcnt = n_cell_sps(this.cellList(cell).elemTopo(), this.solOrder);

        // Build the Cell -> SP map
        cellSPidx[cell, 1] = nSPs+1;    // First SP of this cell
        cellSPidx[cell, 2] = cellSPcnt; // Number of SPs on this cell (redundant info but convenient for programming)

        // Update the mesh total SP count
        nSPs += cellSPcnt;

        // Resize the SP data arrays
        xyzSP_d  = {1..nSPs, 1..this.nDims};
        metSP_d  = {1..nSPs, 1..this.nDims, 1..this.nDims};
        jacSP_d  = {1..nSPs};
        solSP_d  = {1..nSPs, 1..this.nVars};
        dSolSP_d = {1..nSPs, 1..this.nVars, 1..this.nDims};
      }

      for face in this.faceList_d
      {
        var faceFPcnt = n_face_fps(this.faceList(face).elemTopo(), this.solOrder);

        // Build the Face -> FP map
        faceFPidx[face, 1] = nFPs+1;    // First FP of this face
        faceFPidx[face, 2] = faceFPcnt; // Number of FPs on this face (redundant info but convenient for programming)

        nFPs += faceFPcnt;

        // Resize the FP data arrays
        xyzFP_d  = {1..nFPs, 1..this.nDims};
        solFP_d  = {1..nFPs, 1..2, 1..this.nVars};
        flxFP_d  = {1..nFPs, 1..2, 1..this.nVars};
        dSolFP_d = {1..nFPs, 1..2, 1..this.nVars, 1..this.nDims};
      }
    }

    proc set_points_locations()
    {
      use Parameters.ParamMesh;
      use Polynomials;
      use LinearAlgebra;
      use Mapping;
      use Set;

      init_mapping(minOrder=this.solOrder, maxOrder=this.solOrder, this.cellTypes|this.faceTypes);
      init_mapping_metrics(minOrder=this.solOrder, maxOrder=this.solOrder, this.cellTypes|this.faceTypes);

      for cellIdx in this.cellList_d
      {
        // Get the list of nodes that define this cell
        var elemNodes : [this.cellList[cellIdx].nodes_d] int = this.cellList[cellIdx].nodes;

        // Get the coordinates of these nodes
        var xyzMshNodes : [1..this.nDims, 1..elem_nodes(this.cellList[cellIdx].elemType)] real;
        for nodeIdx in elemNodes.domain do
          xyzMshNodes[.., nodeIdx] = this.nodeList[elemNodes[nodeIdx]].xyz[1..this.nDims];

        var cellType : 2*int = (this.cellList[cellIdx].elemType, this.solOrder);

        //////////////////////////
        ///   SP coordinates   ///
        //////////////////////////

        for spIdx in 1..#cellSPidx[cellIdx, 2] do
          this.xyzSP[cellSPidx[cellIdx, 1]+spIdx-1, ..] = dot( mapping[cellType]!.coefs[spIdx, ..],
                                                               xyzMshNodes[..,..].T               );

        /////////////////////////////////////
        ///   Calculate mapping metrics   ///
        /////////////////////////////////////

        for rstDim in 1..this.nDims do
          for spIdx in 1..#cellSPidx[cellIdx, 2] do
            this.metSP[cellSPidx[cellIdx, 1]+spIdx-1, .., rstDim] = dot( mappingMetrics[cellType]!.coefs[rstDim, spIdx, ..],
                                                                         xyzMshNodes[..,..].T                              );

        //////////////////////////////
        ///   Calculate Jacobian   ///
        //////////////////////////////

        // Calculate the Jacobian at SPs
        for spIdx in 1..#cellSPidx[cellIdx, 2] do
          this.jacSP[cellSPidx[cellIdx, 1]+spIdx-1] = determinant(this.metSP[cellSPidx[cellIdx, 1]+spIdx-1, .., ..]);
      }

      for faceIdx in this.faceList_d
      {
        // Get the list of nodes that define this face
        var elemNodes : [this.faceList[faceIdx].nodes_d] int = this.faceList[faceIdx].nodes;

        // Get the coordinates of these nodes
        var xyzMshNodes : [1..this.nDims, 1..elem_nodes(this.faceList[faceIdx].elemType)] real;
        for nodeIdx in elemNodes.domain do
          xyzMshNodes[.., nodeIdx] = this.nodeList[elemNodes[nodeIdx]].xyz[1..this.nDims];

        var faceType : 2*int = (this.faceList[faceIdx].elemType, this.solOrder);

        //////////////////////////
        ///   FP coordinates   ///
        //////////////////////////

        for faceFPIdx in 1..faceFPidx[faceIdx, 2] do
          this.xyzFP[faceFPidx[faceIdx, 1]+faceFPIdx-1, ..] = dot( mapping[faceType]!.coefs[faceFPIdx, ..],
                                                                   xyzMshNodes[..,..].T               );

        //////////////////////
        ///   FP normals   ///
        //////////////////////

        // For node faces in 1D meshes check if normal needs to be reversed
        var leftCell : int = this.faceList[faceIdx].cells[1];

        // Check in which position relative to the left cell this face is on
        var reverse  : real = 1.0;
        var leftFace : int = this.cellList[leftCell].faces[1];
        if (this.nDims == 1) && (faceIdx == leftFace) then
          reverse = -1.0;

        for faceFPIdx in 1..faceFPidx[faceIdx, 2]
        {
          // Initialize metrics so that for 1D and 2D meshes the cross product results in a correct normal
          var metricsFP : [1..2, 1..3] real = reshape([0, 1, 0,
                                                       0, 0, 1], {1..2, 1..3});

          for rstDim in 1..elem_dimension_type(this.faceList[faceIdx].elemType) do
            metricsFP[rstDim, 1..this.nDims] = dot( mappingMetrics[faceType]!.coefs[rstDim, faceFPIdx, ..],
                                                    xyzMshNodes[..,..].T                               );

          // Get face normal through the cross product and slice it to the appropriate dimension vector
          this.nrmFP[faceFPidx[faceIdx, 1]+faceFPIdx-1, ..] = reverse * cross( reshape(metricsFP[1, 1..3], {1..3}),
                                                                               reshape(metricsFP[2, 1..3], {1..3}) )[1..this.nDims];
        }
      }
    }
  }

  proc determinant(metrics : [] real) : real
  {
    use LinearAlgebra;

    var jacobian : real;

    if metrics.size == 1 then
      jacobian = metrics[1,1];
    else if metrics.size == 4 then
      jacobian = metrics[1,1]*metrics[2,2] - metrics[1,2]*metrics[2,1];
    else if metrics.size == 9 then
      jacobian = metrics[1,1]*(metrics[2,2]*metrics[3,3] - metrics[2,3]*metrics[3,2])
                +metrics[1,2]*(metrics[2,3]*metrics[3,1] - metrics[2,1]*metrics[3,3])
                +metrics[1,3]*(metrics[2,1]*metrics[3,2] - metrics[2,2]*metrics[3,1]);
    else
      jacobian = det(metrics);

    return jacobian;
  }

  proc n_cell_sps(in elemTopo : int, in solOrder) : int
  {
    use Parameters.ParamMesh;

    select elemTopo
    {
      when TOPO_NODE do return 1;
      when TOPO_LINE do return (solOrder+1);
      when TOPO_TRIA do return (solOrder+1)*(solOrder+2)/2;
      when TOPO_QUAD do return (solOrder+1)**2;
      when TOPO_TETR do return (solOrder+1)*(solOrder+2)*(solOrder+3)/6;
      when TOPO_PYRA do return (solOrder+1)*(solOrder+2)*(2*solOrder+3)/6;
      when TOPO_PRIS do return (solOrder+1)*(solOrder+1)*(solOrder+2)/2;
      when TOPO_HEXA do return (solOrder+1)**3;
      otherwise return -1;
    }
  }

  proc n_cell_fps(in elemTopo : int, in solOrder) : int
  {
    use Parameters.ParamMesh;

    select elemTopo
    {
      when TOPO_LINE do return 2;
      when TOPO_TRIA do return 3*(solOrder+1);
      when TOPO_QUAD do return 4*(solOrder+1);
      when TOPO_TETR do return 4*(solOrder+1)*(solOrder+2)/2;
      when TOPO_PYRA do return 4*(solOrder+1)*(solOrder+2)/2 +1*(solOrder+1)**2;
      when TOPO_PRIS do return 2*(solOrder+1)*(solOrder+2)/2 +3*(solOrder+1)**2;
      when TOPO_HEXA do return 6*(solOrder+1)**2;
      otherwise return -1;
    }
  }

  proc n_face_fps(in elemTopo : int, in solOrder) : int
  {
    use Parameters.ParamMesh;

    select elemTopo
    {
      when TOPO_NODE do return 1;
      when TOPO_LINE do return (solOrder+1);
      when TOPO_TRIA do return (solOrder+1)*(solOrder+2)/2;
      when TOPO_QUAD do return (solOrder+1)**2;
      otherwise return -1;
    }
  }

  proc main()
  {
    use Gmesh;

    var test_gmesh2 = new unmanaged gmesh2_c();
    test_gmesh2.random1D(nCells=6, xMin=-1.0, xMax=1.0);

    // Get number of physical dimensions from mesh or input
    var test_mesh = new unmanaged mesh_c(nDims=1);
    test_mesh.import_gmesh2(test_gmesh2);

    writeln("Test 1: Random 1D mesh - Gmesh => Native:");
    writeln(test_mesh);
    writeln();

    // Get number of physical dimensions from mesh or input
    var test_frmesh = new unmanaged fr_mesh_c(nDims=1, nVars=3, solOrder=2);
    test_frmesh.import_gmesh2(test_gmesh2);
    test_frmesh.allocate_fr_vars();

    writeln("Test 2: Allocate FR variables");
    writeln("nVarsSP  = ", test_frmesh.nVars);
    writeln("solOrder = ", test_frmesh.solOrder);
    writeln();
    writeln("SP Coordinates     ", test_frmesh.xyzSP.domain);
    writeln("FP Coordinates     ", test_frmesh.xyzFP.domain);
    writeln("SP Metric Terms    ", test_frmesh.metSP.domain);
    writeln("SP Jacobian        ", test_frmesh.jacSP.domain);
    writeln("FP Normal          ", test_frmesh.nrmFP.domain);
    writeln("SP Solution (old)  ", test_frmesh.oldSolSP.domain);
    writeln("SP Solution        ", test_frmesh.solSP.domain);
    writeln("FP Solution        ", test_frmesh.solFP.domain);
    writeln("FP Flux            ", test_frmesh.flxFP.domain);
    writeln("SP d/Solution      ", test_frmesh.dSolSP.domain);
    writeln("FP d/Solution      ", test_frmesh.dSolFP.domain);
    writeln("SP Residue         ", test_frmesh.resSP.domain);
    writeln("SP Sparse Idx      ", test_frmesh.cellSPidx.domain);
    writeln("FP Sparse Idx      ", test_frmesh.faceFPidx.domain);
    writeln();

    // Calculate cell transformations. SPs and FPs locations, metric terms and Jacobians.
    test_frmesh.set_points_locations();

    writeln("Test 3: Calculate cell coordinate transformations and metric terms");
    writeln(test_frmesh);
    writeln();
  }
}

