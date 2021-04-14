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
    var  metFP_d : domain(rank=4, idxType=int);    // {1..nFPs, 1..2, 1..nDims, 1..nDims} nDims+1 if moving mesh
    var  jacSP_d : domain(rank=1, idxType=int);    // {1..nSPs}                           nDims+1 if moving mesh
    var  jacFP_d : domain(rank=2, idxType=int);    // {1..nFPs, 1..2}                     nDims+1 if moving mesh

    var  solSP_d : domain(rank=2, idxType=int);    // {1..nSPs, 1..nVars}
    var  solFP_d : domain(rank=3, idxType=int);    // {1..nFPs, 1..2, 1..nVars}
    var  flxFP_d : domain(rank=2, idxType=int);    // {1..nFPs, 1..nVars}
    var dSolSP_d : domain(rank=3, idxType=int);    // {1..nSPs, 1..nVars, 1..nVars}
    var dSolFP_d : domain(rank=4, idxType=int);    // {1..nFPs, 1..2, 1..nVars, 1..nVars}

    var cellSPidx_d : domain(rank=2, idxType=int); // {1..nCells, 1..2}
    var faceFPidx_d : domain(rank=2, idxType=int); // {1..nFaces, 1..2}

    // FR solver variables
    var xyzSP : [xyzSP_d] real;
    var xyzFP : [xyzFP_d] real;

    var metSP : [metSP_d] real;   // First degree metric terms, Jacobian matrix
    var metFP : [metFP_d] real;   // First degree metric terms, Jacobian matrix
    var jacSP : [jacSP_d] real;   // Jacobian
    var jacFP : [jacFP_d] real;   // Jacobian

    var oldSolSP : [ solSP_d] real;     // Backup of the solution at the beginning of residue calculation
    var    solSP : [ solSP_d] real;     // Conserved variables at SPs
    var    solFP : [ solFP_d] real;     // Conserved variables at FPs (0-Left / 1-right)
    var    flxFP : [ flxFP_d] real;     // Unique convective flux at FPs
    var   dSolSP : [dSolSP_d] real;     // Gradient, at the SPs, of the discontinuous solution interpolation
    var   dSolFP : [dSolFP_d] real;     // Gradient, at the FPs, of the discontinuous flux reconstruction
    var    resSP : [ solSP_d] real;     // conserved variables residual

    var cellSPidx : [cellSPidx_d] int;  // Index of the first SP and number of SPs of a cell
    var faceFPidx : [faceFPidx_d] int;  // Index of the first FP and number of FPs of a face

    proc allocate_vars()
    {
      var nSPs : int;
      var nFPs : int;

      cellSPidx_d = {this.cellList_d.dim(0), 1..2};
      faceFPidx_d = {1..this.nFaces, 1..2};

      for cell in this.cellList_d
      {
        cellSPidx[cell, 1] = nSPs+1;
        cellSPidx[cell, 2] = n_cell_sps(this.cellList(cell).elemTopo(), this.solOrder);

        nSPs += cellSPidx[cell, 2];

        xyzSP_d  = {1..nSPs, 1..this.nDims};

        metSP_d    = {1..nSPs, 1..this.nDims, 1..this.nDims};
        jacSP_d    = {1..nSPs};

        solSP_d  = {1..nSPs, 1..this.nVars};
        dSolSP_d = {1..nSPs, 1..this.nVars, 1..this.nDims};
      }

      for face in this.faceList_d
      {
        faceFPidx[face, 1] = nFPs+1;
        faceFPidx[face, 2] = n_face_fps(this.faceList(face).elemTopo(), this.solOrder);

        nFPs += faceFPidx[face, 2];

        xyzFP_d  = {1..nFPs, 1..this.nDims};

        metFP_d    = {1..nFPs, 1..2, 1..this.nDims, 1..this.nDims};
        jacFP_d    = {1..nFPs, 1..2};

        solFP_d  = {1..nFPs, 1..2, 1..this.nVars};
        flxFP_d  = {1..nFPs, 1..this.nVars};
        dSolFP_d = {1..nFPs, 1..2, 1..this.nVars, 1..this.nDims};
      }
    }

    proc set_points_locations()
    {
      use Parameters.ParamMesh;
      use Polynomials;
      use LinearAlgebra;

      for cell in this.cellList_d
      {
        var nodes_d : domain(rank=1, idxType=int) = this.cellList[cell].nodes_d;
        var nodes : [nodes_d] int = this.cellList[cell].nodes;          // List of nodes that define this cell
        var xyzStdSPs   : [1..this.nDims, 1..cellSPidx[cell, 2]] real;  // List of SPs that belong to this cell

        var xyzMshNodes : [1..this.nDims, nodes_d.dim(0)] real; // Physical coordinates of the nodes
        var xyzStdNodes : [1..this.nDims, nodes_d.dim(0)] real; // Computational coordinates of the nodes

        var coefficients : [1..this.nDims, nodes_d.dim(0)] real; // Coefficients of the transformation polynomial
        var exponents    : [1..this.nDims, nodes_d.dim(0)] int;  // Exponents of the transformation polynomial
        var mat : [nodes_d.dim(0), nodes_d.dim(0)] real = 1; // Coordinate transformation eq system matrix

        for i in nodes_d do
            xyzMshNodes[1..this.nDims,i] = this.nodeList[nodes[i]].xyz[1..this.nDims];

        // The B vector holds one of the physical coordinates in the physical domain of each of the mesh reference points.
        // Ex: B_x = { x[p_1], x[p_2], x[p_3], ... , x[p_n] }

        // The matrix M is composed of lines with the powers of the computational (reference) domain coordinates of a
        //    given reference point.
        // Each line of the matrix has the powers of the computational coordinates of a reference point

        select this.cellList[cell].elemTopo()
        {
          when TOPO_LINE
          {
            //    A_x = { coef_x[xi^0], coef_x[xi^1], coef_x[xi^2], ... , coef_x[xi^n] }
            //    B_x = {      x[p_1] ,      x[p_2] ,      x[p_3] , ... ,      x[p_n]  }
            //
            //    M_x = {   xi[p_1]^0 ,   xi[p_1]^1 ,   xi[p_1]^2 , ... ,    xi[p_1]^n
            //              xi[p_2]^0 ,   xi[p_2]^1 ,   xi[p_2]^2 , ... ,    xi[p_2]^n
            //              xi[p_3]^0 ,   xi[p_3]^1 ,   xi[p_3]^2 , ... ,    xi[p_3]^n
            //                 ...    ,      ...    ,      ...    , ... ,       ...
            //              xi[p_n]^0 ,   xi[p_n]^1 ,   xi[p_n]^2 , ... ,    xi[p_n]^n }

            // Get the coordinates of the standard/reference/computational domain nodes
            // Vertices
            xyzStdNodes[1, 1] = -1;
            xyzStdNodes[1, 2] = +1;
            // Edge nodes
            for node in nodes_d do
              xyzStdNodes[1, 2..nodes_d.high-2] = nodes_uniform(nodes_d.high-2);

            // Get the exponents
            for j in nodes_d do exponents[1, j] = j-1;

            // SP locations
            xyzStdSPs[1, 1..this.solOrder] = nodes_legendre_gauss(this.solOrder);
          }
          when TOPO_TRIA {}
          when TOPO_QUAD
          {
            // Get the coordinates of the standard/reference/computational domain nodes
            // Vertices
            xyzStdNodes[1, 1] = -1;
            xyzStdNodes[2, 1] = +1;
            xyzStdNodes[1, 2] = -1;
            xyzStdNodes[2, 2] = +1;
            xyzStdNodes[1, 3] = -1;
            xyzStdNodes[2, 3] = +1;
            xyzStdNodes[1, 4] = -1;
            xyzStdNodes[2, 4] = +1;
            // Edge nodes
            for node in nodes_d do
              xyzStdNodes[1, 2..nodes_d.high-2] = nodes_uniform(nodes_d.high-2);

            // Get the exponents
            for node in nodes_d
            {
              exponents[1, node] = (node-1) / nodes_d.high**0.5:int; // Xi  exponent
              exponents[2, node] = (node-1) % nodes_d.high**0.5:int; // Eta exponent
            }
          }
          when TOPO_TETR {}
          when TOPO_PYRA {}
          when TOPO_PRIS {}
          when TOPO_HEXA {}
          otherwise {}
        }

        // Build the matrix
        for (i, j) in mat.domain do
          for k in 1..this.nDims do
            mat[i, j] *= xyzStdNodes[k, i]**exponents[k, j];

        //////////////////////
        // Solve the system //
        //////////////////////

        // Using Chapel Linear Algebra Library
        for dim in 1..this.nDims do
          coefficients[dim, ..] = solve(reshape(mat,{nodes_d.dim(0)-1, nodes_d.dim(0)-1}),
                                        reshape(xyzMshNodes[1,..], {nodes_d.dim(0)-1})  );

        // Future 1: Using Lapack QR factorization
        //geqrf(matrix_order: lapack_memory_order = 101, a: [] real(32), tau: [] real(32));


        // Future 2: Using interpolation and Lagrange derivative

        /////////////////////////////////////////////////////////
        ///   Calculate transformation Metrics and Jacobian   ///
        /////////////////////////////////////////////////////////

        // Calculate metric terms for SPs
        for (sp, mshDim, stdDim) in {1..#cellSPidx[cell, 2], 1..this.nDims, 1..this.nDims} do
            this.metSP[cellSPidx[cell, 1]+sp-1, mshDim, stdDim] = coefficients[mshDim, 2];

        // Get faces from this cell
        for face in this.cellList[cell].faces do
          // Get FPs from this face
          for fp in faceFPidx[face,1]..#faceFPidx[face,2] do
            // Find out if we write to the Left or Right FPs of the face
            if this.faceList[face].cells[1] == cell then
              this.metFP[fp, 1, .., ..] = coefficients[1, 2];
            else if this.faceList[face].cells[2] == cell then
              this.metFP[fp, 2, .., ..] = coefficients[1, 2];
            else
              writeln("Error: Inconsistent mesh data found");

        // Calculate the Jacobian at SPs
        for sp in 1..#cellSPidx[cell, 2] do
          this.jacSP[cellSPidx[cell, 1]+sp-1] = this.metSP[cellSPidx[cell, 1]+sp-1, 1, 1];

        // Calculate the Jacobian at FPs
        for face in this.cellList[cell].faces do
          for fp in faceFPidx[face,1]..#faceFPidx[face,2] do
              if this.faceList[face].cells[1] == cell then
                this.jacFP[fp, 1] = this.metFP[fp, 1, 1, 1];
              else if this.faceList[face].cells[2] == cell then
                this.jacFP[fp, 2] = this.metFP[fp, 2, 1, 1];

        //////////////////////////////
        ///   Points coordinates   ///
        //////////////////////////////

        // Calculate the Jacobian at SPs
        for sp in 1..#cellSPidx[cell, 2] do
          this.xyzSP[cellSPidx[cell, 1]+sp-1, 1] = coefficients[1,1]+xyzStdSPs[1,sp]*coefficients[1,2];

        // Calculate the Jacobian at FPs
        for face in this.cellList[cell].faces do
          for fp in faceFPidx[face,1]..#faceFPidx[face,2] do
              if this.cellList[cell].faces[1] == face then
                this.xyzFP[fp, 1] = coefficients[1,1]-1*coefficients[1,2];
              else if this.cellList[cell].faces[2] == face then
                this.xyzFP[fp, 1] = coefficients[1,1]+1*coefficients[1,2];
      }
    }
  }

  proc n_cell_sps(in elemTopo : int, in solOrder) : int
  {
    use Parameters.ParamMesh;

    select elemTopo
    {
      when TOPO_LINE do return (solOrder+1);
      when TOPO_TRIA do return (solOrder+1)*(solOrder+2)/2;
      when TOPO_QUAD do return (solOrder+1)**2;
      when TOPO_TETR do return (solOrder+1)*(solOrder+2)*(solOrder+3)/6;
      when TOPO_PYRA do return (solOrder+1)*(solOrder+2)*(2*solOrder+1)/6;
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

    var test_gmesh2 = new unmanaged gmesh2_c(nNodes=7, nElements=8, nFamilies=3);
    test_gmesh2.random1D(-1,1);

    // Get number of physical dimensions from mesh or input
    var test_mesh = new unmanaged mesh_c(nDims=1);
    test_mesh.import_gmesh2(test_gmesh2);

    writeln("Test 1: Random 1D mesh - Gmesh => Native:");
    writeln(test_mesh);
    writeln();

    // Get number of physical dimensions from mesh or input
    var test_frmesh = new unmanaged fr_mesh_c(nDims=1, nVars=3, solOrder=2);
    test_frmesh.import_gmesh2(test_gmesh2);
    test_frmesh.allocate_vars();

    writeln("Test 2: Allocate FR variables");
    writeln("nVarsSP  = ", test_frmesh.nVars);
    writeln("solOrder = ", test_frmesh.solOrder);
    writeln();
    writeln("SP Coordinates     ", test_frmesh.xyzSP.domain);
    writeln("FP Coordinates     ", test_frmesh.xyzFP.domain);
    writeln("SP Metric Terms    ", test_frmesh.metSP.domain);
    writeln("FP Metric Terms    ", test_frmesh.metFP.domain);
    writeln("SP Jacobian        ", test_frmesh.jacSP.domain);
    writeln("FP Jacobian        ", test_frmesh.jacFP.domain);
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

