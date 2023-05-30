module FRMesh {
  use Gmesh;
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

    var  solSP_d : domain(rank=2, idxType=int);    // {1..nVars, 1..nSPs}
    var  solFP_d : domain(rank=3, idxType=int);    // {1..nFPs, 1..2, 1..nVars}
    var  flxFP_d : domain(rank=3, idxType=int);    // {1..nFPs, 1..2, 1..nVars}
    // For future viscous flow implementation
    //var dSolSP_d : domain(rank=3, idxType=int);    // {1..nSPs, 1..nVars, 1..nVars}
    //var dSolFP_d : domain(rank=4, idxType=int);    // {1..nFPs, 1..2, 1..nVars, 1..nVars}

    var cellSPidx_d : domain(rank=2, idxType=int); // {1..nCells, 1..2}
    var faceFPidx_d : domain(rank=2, idxType=int); // {1..nFaces, 1..2}

    // FR solver variables

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

    var cellSPidx : [cellSPidx_d] int;  // Index of the first SP and number of SPs of a cell
    var faceFPidx : [faceFPidx_d] int;  // Index of the first FP and number of FPs of a face

    override proc import_gmesh2(gmesh : gmesh2_c)
    {
      use Gmesh;

      // Copy nodes
      this.nNodes = gmesh.nodes.domain.dim(0).high;
      this.nodeList_d = {1..this.nNodes};
      for node in this.nodeList_d do
        this.nodeList[node].xyz = gmesh.nodes[node,1..3];

      // Get family names and dimensions. These are used to sort the elements.
      this.nFamls = gmesh.families.domain.dim(0).high;
      this.famlList_d = {1..this.nFamls};
      for faml in this.famlList_d
      {
        this.famlList[faml].name = gmesh.families[faml].name;
        this.famlList[faml].nDim = gmesh.families[faml].nDim;
      }
      var cellDim : int = max reduce this.famlList[..].nDim;
      var bocoDim : int = cellDim - 1;

      // Sort elements into cells and boundaries and copy them
      for element in gmesh.elements.domain
      {
        var elemDim     : int;
        var elemFamlIdx : int;

        // Search though families for this element's family.
        for famlIdx in gmesh.families.domain
        {
          // Gmesh family tags are not globally unique, they are unique within families with the same number of spatial
          // dimensions therefore we search for a match of both criteria.
          if gmesh.elements[element].tags[1]   == gmesh.families[famlIdx].tag
          && gmesh.elements[element].elemDim() == gmesh.families[famlIdx].nDim
          {
            elemDim = gmesh.families[famlIdx].nDim;
            elemFamlIdx = famlIdx;
          }
        }

        // Add the mesh element to either the cell or boco list depending on the element/family dimension
        select elemDim
        {
          when cellDim
          {
            // Increment cell count with new element
            this.nCells += 1;
            // Resize domain to expand the array
            this.cellList_d = {1..this.nCells};
            // Fill up the cell properties. Maybe this should be an initializer?
            this.cellList[this.nCells].nodes_d  = gmesh.elements[element].nodes_d;
            this.cellList[this.nCells].nodes    = gmesh.elements[element].nodes;
            this.cellList[this.nCells].elemType = elem_type_gmsh2mesh(gmesh.elements[element].elemType);
            this.cellList[this.nCells].family   = elemFamlIdx;
          }
          when bocoDim
          {
            // Increment boco count with new element
            this.nBocos += 1;
            // Resize domain to expand the array
            this.bocoList_d = {1..this.nBocos};
            // Fill up the boco properties. Maybe this should be an initializer?
            this.bocoList[this.nBocos].nodes_d  = gmesh.elements[element].nodes_d;
            this.bocoList[this.nBocos].nodes    = gmesh.elements[element].nodes;
            this.bocoList[this.nBocos].elemType = elem_type_gmsh2mesh(gmesh.elements[element].elemType);
            this.bocoList[this.nBocos].family   = elemFamlIdx;
          }
          otherwise do
          {
            writeln("Found element of unexpected dimension in mesh.");
            writeln("   Maximum dimension element found: ", cellDim);
            writeln("   Assumed cell dimension: ", cellDim);
            writeln("   Assumed boco dimension: ", bocoDim);
            writeln("Problematic element dimension: ", elemDim, " ID: ", element);
            writeln();
          }
        }
      }

      // Build the face list for Riemann solver iteration
      this.build_face_list();

      // Build the sets of cell and face element types and topologies present in this mesh
      this.build_elem_sets();
    }

    override proc build_face_list()
    {
      use Parameters.ParamMesh;
      use Parameters.ParamGmesh;
      import SortTuple.sort_tuple;

      // Build face list
      var faceVerts : [1..6] 4*int;
      var faceMap_d : domain(4*int);
      var faceMap : [faceMap_d] int;

      // Add all boundaries to the face map
      for boco in this.bocoList_d
      {
        // Each boundary defines a face. Build the face nodes tuple based on the element topology.
        select elem_topology(this.bocoList[boco].elemType)
        {
          when TOPO_NODE do faceVerts[1] = (this.bocoList[boco].nodes[1], 0, 0, 0);
          when TOPO_LINE do faceVerts[1] = (this.bocoList[boco].nodes[1], this.bocoList[boco].nodes[2], 0, 0);
          when TOPO_TRIA do faceVerts[1] = (this.bocoList[boco].nodes[1], this.bocoList[boco].nodes[2],
                                            this.bocoList[boco].nodes[3], 0);
          when TOPO_QUAD do faceVerts[1] = (this.bocoList[boco].nodes[1], this.bocoList[boco].nodes[2],
                                            this.bocoList[boco].nodes[3], this.bocoList[boco].nodes[4]);
          otherwise {}
        }

        // Create the face element on the face list
        // Increment the face counter with the new face and expand the face list domain
        this.nFaces += 1;
        this.faceList_d = {1..this.nFaces};

        // Add the face vertices to the face matching list and store the face index given to this face.
        faceMap_d.add(sort_tuple(faceVerts[1]));
        faceMap[sort_tuple(faceVerts[1])] = this.nFaces;

        // Fill up the face properties. First cell to contain a face is put of the right side of the face.
        this.faceList[this.nFaces].cells[2] = -boco;
        this.faceList[this.nFaces].elemType = this.bocoList[boco].elemType;

        // Fill the right side neighbor ID, boundaries are always on the right side of a face.
        // Boundaries have negative indexes so they can be easily distinguished from mesh cells.
        this.bocoList[boco].face = this.nFaces;
      }

      // Add faces from cells to the face map and perform the matching
      for cellIdx in this.cellList_d
      {
        // Allocate this cell's local face list
        this.cellList[cellIdx].faces_d = {1..elem_faces(elem_topology(this.cellList[cellIdx].elemType))};

        // Build the face nodes tuple
        select elem_topology(this.cellList[cellIdx].elemType)
        {
          when TOPO_LINE
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], 0, 0, 0);
            faceVerts[2] = (this.cellList[cellIdx].nodes[2], 0, 0, 0);
          }
          when TOPO_TRIA
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2], 0, 0);
            faceVerts[2] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3], 0, 0);
            faceVerts[3] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[1], 0, 0);
          }
          when TOPO_QUAD
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2], 0, 0);
            faceVerts[2] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3], 0, 0);
            faceVerts[3] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[4], 0, 0);
            faceVerts[4] = (this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[1], 0, 0);
          }
          when TOPO_TETR
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[2],
                            0);
            faceVerts[2] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[4],
                            0);
            faceVerts[3] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[4],
                            0);
            faceVerts[4] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[4],
                            0);
          }
          when TOPO_PYRA
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[2]);
            faceVerts[2] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[5],
                            0);
            faceVerts[3] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[5],
                            0);
            faceVerts[4] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[5],
                            0);
            faceVerts[5] = (this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[5],
                            0);
          }
          when TOPO_PRIS
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[5],
                            this.cellList[cellIdx].nodes[4]);
            faceVerts[2] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[6],
                            this.cellList[cellIdx].nodes[5]);
            faceVerts[3] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[4],
                            this.cellList[cellIdx].nodes[6]);
            faceVerts[4] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[2],
                            0);
            faceVerts[5] = (this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[5], this.cellList[cellIdx].nodes[6],
                            0);
          }
          when TOPO_HEXA
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[2]);
            faceVerts[2] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[6],
                            this.cellList[cellIdx].nodes[5]);
            faceVerts[3] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[7],
                            this.cellList[cellIdx].nodes[6]);
            faceVerts[4] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[8],
                            this.cellList[cellIdx].nodes[7]);
            faceVerts[5] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[5], this.cellList[cellIdx].nodes[8],
                            this.cellList[cellIdx].nodes[4]);
            faceVerts[6] = (this.cellList[cellIdx].nodes[5], this.cellList[cellIdx].nodes[6], this.cellList[cellIdx].nodes[7],
                            this.cellList[cellIdx].nodes[8]);
          }
          otherwise {}
        }

        // Iterate through this cell's faces
        for cellFaceIdx in this.cellList[cellIdx].faces.domain
        {
          // Check if this face is already in the map
          if faceMap_d.contains(sort_tuple(faceVerts[cellFaceIdx]))
          { // This a mapped face! :D
            var faceIdx : int = faceMap[sort_tuple(faceVerts[cellFaceIdx])];

            // Save the cell ID as the left neighbor of the face
            this.faceList[faceIdx].cells[1] = cellIdx;

            // Fill face nodes fro mthe left neighbors so that face are consistently oriented
            this.faceList[faceIdx].nodes_d = {1..elem_nodes(this.faceList[this.nFaces].elemType)};
            this.faceList[faceIdx].nodes   = this.cellList[cellIdx].nodes[elem_face_nodes(this.cellList[cellIdx].elemType,
                                                                                              cellFaceIdx                    )];

            // Save the face index and the side of the face this cell is on to the cell record
            this.cellList[cellIdx].faces[cellFaceIdx] = faceIdx;
            this.cellList[cellIdx].sides[cellFaceIdx] = 1;
          }
          else
          { // This isn't a mapped face. But don't panic! D:

            // Increment the face counter with the new face and expand the face list domain
            this.nFaces += 1;
            this.faceList_d = {1..this.nFaces};

            // Add the tuple if nodes that define this face to the map and save this face's face index
            faceMap_d.add(sort_tuple(faceVerts[cellFaceIdx]));
            faceMap[sort_tuple(faceVerts[cellFaceIdx])] = this.nFaces;

            // Fill up the face properties. First cell to contain a face is put of the right side of the face.
            this.faceList[this.nFaces].cells[2] = cellIdx;
            this.faceList[this.nFaces].elemType = elem_face_type(this.cellList[cellIdx].elemType, cellFaceIdx);

            // Save the face index and the side of the face this cell is on to the cell record
            this.cellList[cellIdx].faces[cellFaceIdx] = this.nFaces;
            this.cellList[cellIdx].sides[cellFaceIdx] = 2;
          }
        }
      }
      // Check if all faces have been correctly identified?
    }

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
        solSP_d  = {1..this.nVars, 1..nSPs};

        // For future viscous flow implementation
        //dSolSP_d = {1..nSPs, 1..this.nVars, 1..this.nDims};
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

        // For future viscous flow implementation
        //dSolFP_d = {1..nFPs, 1..2, 1..this.nVars, 1..this.nDims};
      }
    }

    proc set_points_locations()
    {
      use Parameters.ParamMesh;
      use Polynomials;
      use LinearAlgebra;
      use Mapping;
      use Set;
      import Determinant.determinant;

      init_mapping(minOrder=this.solOrder, maxOrder=this.solOrder, this.cellTypes|this.faceTypes);
      init_mapping_metrics(minOrder=this.solOrder, maxOrder=this.solOrder, this.cellTypes|this.faceTypes);

      for cellIdx in this.cellList_d
      {
        // Get loop variables
        ref cellSPini = this.cellSPidx[cellIdx, 1];
        ref cellSPcnt = this.cellSPidx[cellIdx, 2];

        // Get the list of nodes that define this cell
        var elemNodes : [this.cellList[cellIdx].nodes_d] int = this.cellList[cellIdx].nodes;

        // Get the coordinates of these nodes
        var xyzMshNodes : [1..this.nDims, 1..elem_nodes(this.cellList[cellIdx].elemType)] real;
        for nodeIdx in elemNodes.domain do
          xyzMshNodes[.., nodeIdx] = this.nodeList[elemNodes[nodeIdx]].xyz[1..this.nDims];

        var cellType : 2*int = (this.cellList[cellIdx].elemType, this.solOrder);

        for cellSP in 1..cellSPcnt
        {
          var meshSP = cellSPini + cellSP - 1;

          // Calculate the physical coordinates of this cell's SPs
          this.xyzSP[meshSP, ..] = dot( mapping[cellType]!.coefs[cellSP, ..],
                                        xyzMshNodes[..,..].T               );

          // Calculate the derivatives of the physical coordinates by the computational coordinates
          // {{x_xi, x_eta, x_zeta},
          //  {y_xi, y_eta, y_zeta},
          //  {z_xi, z_eta, z_zeta}}
          for rstDim in 1..this.nDims do
            this.metSP[meshSP, .., rstDim] = dot( mappingMetrics[cellType]!.coefs[rstDim, cellSP, ..],
                                                  xyzMshNodes[..,..].T                              );

          // Calculate the Jacobian at SPs
          this.jacSP[meshSP] = determinant(this.metSP[meshSP, .., ..]);

          // Invert the metric matrix since we only use that
          // {{  xi_x,   xi_y,   xi_z},
          //  { eta_x,  eta_y,  eta_z},
          //  {zeta_x, zeta_y, zeta_z}}
          this.metSP[meshSP, .., ..] = inv(this.metSP[meshSP, .., ..]);
        }
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
