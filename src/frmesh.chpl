module FRMesh {
  use Gmesh;
  use Input;
  use Mesh;
  use BlockDist;
  //use Block1DDist;
  use ReplicatedDist;

  class fr_mesh_c : mesh_c
  {
    const nVars  : int;
    const solOrder : int;

    const nSPs : int;
    const nFPs : int;

    // Array sizing domains
    const cellSPidxSpace : domain(rank=2, idxType=int); // {1..nCells, 1..2}
    const faceFPidxSpace : domain(rank=2, idxType=int); // {1..nFaces, 1..2}
    const xyzSPSpace     : domain(rank=2, idxType=int); // {1..nSPs, 1..nDims}
    const xyzFPSpace     : domain(rank=2, idxType=int); // {1..nFPs, 1..nDims}
    const metSPSpace     : domain(rank=3, idxType=int); // {1..nSPs, 1..nDims, 1..nDims}       nDims+1 if moving mesh
    const jacSPSpace     : domain(rank=1, idxType=int); // {1..nSPs}                           nDims+1 if moving mesh
    const solSPSpace     : domain(rank=2, idxType=int); // {1..nVars, 1..nSPs}
    const solFPSpace     : domain(rank=3, idxType=int); // {1..nFPs, 1..2, 1..nVars}
    const flxFPSpace     : domain(rank=3, idxType=int); // {1..nFPs, 1..2, 1..nVars}

    // Sincgle Locale Domains
    //var cellSPidx_d : cellSPidxSpace ;
    //var faceFPidx_d : faceFPidxSpace ;
    //var xyzSP_d     : xyzSPSpace     ;
    //var xyzFP_d     : xyzFPSpace     ;
    //var metSP_d     : metSPSpace     ;
    //var jacSP_d     : jacSPSpace     ;
    //var solSP_d     : solSPSpace     ;
    //var solFP_d     : solFPSpace     ;
    //var flxFP_d     : flxFPSpace     ;

    // BlockDist Domains
    var cellSPidx_d = cellSPidxSpace dmapped Block(boundingBox=cellSPidxSpace);
    var faceFPidx_d = faceFPidxSpace dmapped Block(boundingBox=faceFPidxSpace);
    var xyzSP_d     = xyzSPSpace     dmapped Block(boundingBox=xyzSPSpace    );
    var xyzFP_d     = xyzFPSpace     dmapped Block(boundingBox=xyzFPSpace    );
    var metSP_d     = metSPSpace     dmapped Block(boundingBox=metSPSpace    );
    var jacSP_d     = jacSPSpace     dmapped Block(boundingBox=jacSPSpace    );
    var solSP_d     = solSPSpace     dmapped Block(boundingBox=solSPSpace    );
    var solFP_d     = solFPSpace     dmapped Block(boundingBox=solFPSpace    );
    var flxFP_d     = flxFPSpace     dmapped Block(boundingBox=flxFPSpace    );

    // BlockDist Domains
    //var cellSPidx_d = cellSPidxSpace dmapped Block1D(boundingBox=cellSPidxSpace);
    //var faceFPidx_d = faceFPidxSpace dmapped Block1D(boundingBox=faceFPidxSpace);
    //var xyzSP_d     = xyzSPSpace     dmapped Block1D(boundingBox=xyzSPSpace    );
    //var xyzFP_d     = xyzFPSpace     dmapped Block1D(boundingBox=xyzFPSpace    );
    //var metSP_d     = metSPSpace     dmapped Block1D(boundingBox=metSPSpace    );
    //var jacSP_d     = jacSPSpace     dmapped Block1D(boundingBox=jacSPSpace    );
    //var solSP_d     = solSPSpace     dmapped Block1D(boundingBox=solSPSpace    );
    //var solFP_d     = solFPSpace     dmapped Block1D(boundingBox=solFPSpace    );
    //var flxFP_d     = flxFPSpace     dmapped Block1D(boundingBox=flxFPSpace    );

    // FR solver variables

    // Point indexing arrays
    var cellSPidx : [cellSPidx_d] int;  // Index of the first SP and number of SPs of a cell
    var faceFPidx : [faceFPidx_d] int;  // Index of the first FP and number of FPs of a face

    // Physical coordinates
    var xyzSP : [xyzSP_d] real;
    var xyzFP : [xyzFP_d] real;

    // Transformation metrics
    var metSP : [metSP_d] real;   // First degree metric terms, Jacobian matrix
    var jacSP : [jacSP_d] real;   // Jacobian
    var nrmFP : [xyzFP_d] real;   // Face normal vector pointing from the left (side 1) to the right (side 2) cell

    // Flow Variables
    var oldSolSP : [ solSP_d] real;     // Backup of the solution at the beginning of residue calculation
    var    solSP : [ solSP_d] real;     // Conserved variables at SPs
    var    solFP : [ solFP_d] real;     // Conserved variables at FPs (1-Left / 2-right)
    var    flxFP : [ flxFP_d] real;     // Discontinuous flux at FPs (1-Left / 2-right)
    // For future viscous flow implementation
    //var   dSolSP : [dSolSP_d] real;     // Gradient, at the SPs, of the discontinuous solution interpolation
    //var   dSolFP : [dSolFP_d] real;     // Gradient, at the FPs, of the discontinuous flux reconstruction
    var    resSP : [ solSP_d] real;     // Conserved variables residual

    proc init(mesh : gmesh2_c, nVars : int, solOrder : int)
    {
      super.init(mesh);

      this.nVars = nVars;
      this.solOrder = solOrder;

      var spCnt : int;
      const cellTopoCnt = mesh.cell_topo_cnt();
      for cellTopo in cellTopoCnt.keys() do
        spCnt = try! cellTopoCnt[cellTopo] * n_cell_sps(cellTopo, this.solOrder);
      this.nSPs = spCnt;
      //writeln("Mesh allocated with ", this.nSPs, " SPs");

      var fpCnt : int;
      const faceTopoCnt = mesh.face_topo_cnt();
      for faceTopo in faceTopoCnt.keys() do
        fpCnt = try! faceTopoCnt[faceTopo] * n_face_fps(faceTopo, this.solOrder);
      this.nFPs = fpCnt;
      //writeln("Mesh allocated with ", this.nFPs, " FPs");

      this.cellSPidxSpace = {1..this.nCells, 1..2};
      this.faceFPidxSpace = {1..this.nFaces, 1..2};

      this.xyzSPSpace     = {1..this.nSPs, 1..nDims};
      this.xyzFPSpace     = {1..this.nFPs, 1..nDims};
      this.metSPSpace     = {1..this.nSPs, 1..nDims, 1..nDims};
      this.jacSPSpace     = {1..this.nSPs};
      this.solSPSpace     = {1..this.nVars, 1..nSPs};
      this.solFPSpace     = {1..this.nFPs, 1..2, 1..nVars};
      this.flxFPSpace     = {1..this.nFPs, 1..2, 1..nVars};

      this.cellSPidx_d = this.cellSPidxSpace dmapped Block(boundingBox=this.cellSPidxSpace);
      this.faceFPidx_d = this.faceFPidxSpace dmapped Block(boundingBox=this.faceFPidxSpace);

      this.xyzSP_d     = this.xyzSPSpace     dmapped Block(boundingBox=this.xyzSPSpace    );
      this.xyzFP_d     = this.xyzFPSpace     dmapped Block(boundingBox=this.xyzFPSpace    );
      this.metSP_d     = this.metSPSpace     dmapped Block(boundingBox=this.metSPSpace    );
      this.jacSP_d     = this.jacSPSpace     dmapped Block(boundingBox=this.jacSPSpace    );
      this.solSP_d     = this.solSPSpace     dmapped Block(boundingBox=this.solSPSpace    );
      this.solFP_d     = this.solFPSpace     dmapped Block(boundingBox=this.solFPSpace    );
      this.flxFP_d     = this.flxFPSpace     dmapped Block(boundingBox=this.flxFPSpace    );
    }

    override proc import_gmesh2(gmesh : gmesh2_c)
    {
      super.import_gmesh2(gmesh);
    }

    proc allocate_fr_vars()
    {
      var spCnt : int = 0;
      var fpCnt : int = 0;

      cellSPidx_d = {this.cellList_d.dim(0), 1..2};
      faceFPidx_d = {this.faceList_d.dim(0), 1..2};

      for cellIdx in this.cellList.domain
      {
        var cellSPcnt = n_cell_sps(this.cellList(cellIdx).elemTopo(), this.solOrder);

        // Build the Cell -> SP map
        cellSPidx[cellIdx, 1] = spCnt+1;    // First SP of this cell
        cellSPidx[cellIdx, 2] = cellSPcnt; // Number of SPs on this cell (redundant info but convenient for programming)

        // Update the mesh total SP count
        spCnt += cellSPcnt;

        // Resize the SP data arrays
        xyzSP_d  = {1..spCnt, 1..this.nDims};
        metSP_d  = {1..spCnt, 1..this.nDims, 1..this.nDims};
        jacSP_d  = {1..spCnt};
        solSP_d  = {1..this.nVars, 1..spCnt};

        // For future viscous flow implementation
        //dSolSP_d = {1..spCnt, 1..this.nVars, 1..this.nDims};
      }

      for faceIdx in this.faceList.domain
      {
        var faceFPcnt = n_face_fps(this.faceList(faceIdx).elemTopo(), this.solOrder);

        // Build the Face -> FP map
        faceFPidx[faceIdx, 1] = fpCnt+1;    // First FP of this face
        faceFPidx[faceIdx, 2] = faceFPcnt; // Number of FPs on this face (redundant info but convenient for programming)

        fpCnt += faceFPcnt;

        // Resize the FP data arrays
        xyzFP_d  = {1..fpCnt, 1..this.nDims};
        solFP_d  = {1..fpCnt, 1..2, 1..this.nVars};
        flxFP_d  = {1..fpCnt, 1..2, 1..this.nVars};

        // For future viscous flow implementation
        //dSolFP_d = {1..fpCnt, 1..2, 1..this.nVars, 1..this.nDims};
      }
    }

    proc set_points_locations()
    {
      use Parameters.ParamMesh;
      use Polynomials;
      use LinearAlgebra;
      use Mapping;
      use Set;
      import Inverse.inv;
      import Determinant.det;

      init_mapping(minOrder=this.solOrder, maxOrder=this.solOrder, this.cellTypes|this.faceTypes);
      init_mapping_metrics(minOrder=this.solOrder, maxOrder=this.solOrder, this.cellTypes|this.faceTypes);

      this.xyzSP = 0.0;
      this.metSP = 0.0;
      this.jacSP = 0.0;
      this.xyzFP = 0.0;
      this.nrmFP = 0.0;

      forall cellIdx in this.cellList.domain
      {
        // Get loop variables
        const cellSPini = this.cellSPidx[cellIdx, 1];
        const cellSPcnt = this.cellSPidx[cellIdx, 2];

        // Get the list of nodes that define this cell
        const cellNodes : [this.cellList[cellIdx].nodes_d] int = this.cellList[cellIdx].nodes;

        // Get the coordinates of these nodes
        var xyzMshNodes : [1..this.nDims, 1..elem_nodes(this.cellList[cellIdx].elemType)] real;
        for nodeIdx in cellNodes.domain do
          xyzMshNodes[.., nodeIdx] = this.nodeList[cellNodes[nodeIdx]].xyz[1..this.nDims];

        const cellType : 2*int = (this.cellList[cellIdx].elemType, this.solOrder);

        for cellSPidx in 1..cellSPcnt
        {
          const meshSPidx = cellSPini + cellSPidx - 1;

          // Calculate the physical coordinates of this cell's SPs
          //this.xyzSP[meshSPidx, ..] = dot( mapping[cellType]!.coefs[cellSPidx, ..],
          //                                 xyzMshNodes[.., ..].T                  );
          for dimIdx in 1..this.nDims do
            for cellNodeIdx in xyzMshNodes.domain.dim(1) do
              this.xyzSP[meshSPidx, dimIdx] += mapping[cellType]!.coefs[cellSPidx, cellNodeIdx]
                                              *xyzMshNodes[dimIdx, cellNodeIdx];

          // Calculate the derivatives of the physical coordinates by the computational coordinates
          // [[x_xi, x_eta, x_zeta],
          //  [y_xi, y_eta, y_zeta],
          //  [z_xi, z_eta, z_zeta]]
          //for rstDim in 1..this.nDims do
          //  this.metSP[meshSPidx, .., rstDim] = dot( mappingMetrics[cellType]!.coefs[rstDim, cellSPidx, ..],
          //                                           xyzMshNodes[..,..].T                                  );
          for physDimIdx in 1..this.nDims do
            for compDimIdx in 1..this.nDims do
              for cellNodeIdx in xyzMshNodes.domain.dim(1) do
                this.metSP[meshSPidx, physDimIdx, compDimIdx] += mappingMetrics[cellType]!.coefs[compDimIdx, cellSPidx, cellNodeIdx]
                                                                *xyzMshNodes[physDimIdx, cellNodeIdx];

          // Calculate the Jacobian at SPs
          this.jacSP[meshSPidx] = det(this.metSP[meshSPidx, .., ..]);

          // Invert the metric matrix since we only use that
          // [[  xi_x,   xi_y,   xi_z],
          //  [ eta_x,  eta_y,  eta_z],
          //  [zeta_x, zeta_y, zeta_z]]
          //this.metSP[meshSPidx, .., ..] = inv(this.metSP[meshSPidx, .., ..]);
          this.metSP[meshSPidx, .., ..] = inv(this.metSP[meshSPidx, .., ..]);
        }
      }

      forall faceIdx in this.faceList.domain
      {
        // Get loop variables
        const faceFPini = this.faceFPidx[faceIdx, 1];
        const faceFPcnt = this.faceFPidx[faceIdx, 2];

        // Get the list of nodes that define this face
        const faceNodes : [this.faceList[faceIdx].nodes_d] int = this.faceList[faceIdx].nodes;

        // Get the coordinates of these nodes
        var xyzMshNodes : [1..this.nDims, 1..elem_nodes(this.faceList[faceIdx].elemType)] real;
        for nodeIdx in faceNodes.domain do
          xyzMshNodes[.., nodeIdx] = this.nodeList[faceNodes[nodeIdx]].xyz[1..this.nDims];

        const faceType : 2*int = (this.faceList[faceIdx].elemType, this.solOrder);

        // The face normal is calculated based on the left side cell (which is, by convention, never a boundary)
        const leftCell : int = this.faceList[faceIdx].cells[1];

        // Depending on the position of this face in relation to the left cell the normal mingh need to be reversed
        const leftFace : int = this.cellList[leftCell].faces[1];
        const reverse  : real = if (this.nDims == 1) && (faceIdx == leftFace)
          then -1.0
          else +1.0;

        for faceFPidx in 1..faceFPcnt
        {
          const meshFPidx = faceFPini + faceFPidx - 1;

          // Calculate the physical coordinates of this face's FPs
          //this.xyzFP[meshFPidx, ..] = dot( mapping[faceType]!.coefs[faceFPidx, ..],
          //                                 xyzMshNodes[..,..].T                   );
          for dimIdx in 1..this.nDims do
            for faceNodeIdx in xyzMshNodes.domain.dim(1) do
              this.xyzFP[meshFPidx, dimIdx] += mapping[faceType]!.coefs[faceFPidx, faceNodeIdx]
                                              *xyzMshNodes[dimIdx, faceNodeIdx];

          //////////////////////
          ///   FP normals   ///
          //////////////////////

          var metricsFP : [1..2, 1..3] real;

          // Initialize metrics so that for 1D and 2D meshes the cross product results in a correct normal
          if this.nDims == 1 then
            metricsFP = reshape([0, 1, 0,
                                 0, 0, 1], {1..2, 1..3});
          else if this.nDims == 2 then
            metricsFP = reshape([0, 0, 0,
                                 0, 0, 1], {1..2, 1..3});

          //for compDim in 1..elem_dimension_type(this.faceList[faceIdx].elemType) do
          //  metricsFP[compDim, 1..this.nDims] = dot( mappingMetrics[faceType]!.coefs[compDim, faceFPidx, ..],
          //                                           xyzMshNodes[..,..].T                                   );
          for physDimIdx in 1..this.nDims do
            for compDimIdx in 1..this.nDims-1 do
              for cellNodeIdx in xyzMshNodes.domain.dim(1) do
                metricsFP[compDimIdx, physDimIdx] += mappingMetrics[faceType]!.coefs[compDimIdx, faceFPidx, cellNodeIdx]
                                                    *xyzMshNodes[physDimIdx, cellNodeIdx];

          // Get face normal through the cross product and slice it to the appropriate dimension vector
          this.nrmFP[meshFPidx, ..] = reverse * cross( reshape(metricsFP[1, 1..3], {1..3}),
                                                       reshape(metricsFP[2, 1..3], {1..3}) )[1..this.nDims];
        }
      }
    }

    proc calc_time_step()
    {
      use Parameters.ParamInput;
      import Temporal.max_wave_speed;
      import Input.timeStepMethod;
      import Input.timeStep;

      // Option 1: More legible code
      forall cellIdx in this.cellList.domain do
        select timeStepMethod
        {
          when DT_GLOBAL_CONST do cellTimeStep[cellIdx] = Input.timeStep;
          when DT_GLOBAL_CFL   do cellTimeStep[cellIdx] = time_step_cell(cellIdx);
          when DT_LOCAL_CFL    do cellTimeStep[cellIdx] = time_step_cell(cellIdx);
        }

      this.minTimeStep = min reduce(cellTimeStep);

      if timeStepMethod == DT_GLOBAL_CFL then
        this.cellTimeStep = this.minTimeStep;
    }

    proc time_step_cell(cellIdx : int)
    {
      import Temporal.max_wave_speed_array;
      import Input.cfl;

      ref cellSPini : int = this.cellSPidx[cellIdx, 1];
      ref cellSPcnt : int = this.cellSPidx[cellIdx, 2];

      const spWaveSpeed : real = max_wave_speed_array(this.solSP[.., cellSPini.. #cellSPcnt]);

      const timeStep : real = cfl * this.cellCharLeng[cellIdx] / spWaveSpeed;

      return timeStep;
    }
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
    test_mesh.build_cell_char_leng();
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
    // For future viscous flow implementation
    //writeln("SP d/Solution      ", test_frmesh.dSolSP.domain);
    //writeln("FP d/Solution      ", test_frmesh.dSolFP.domain);
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
