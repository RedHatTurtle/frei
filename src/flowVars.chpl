prototype module FlowVars {
  use Input;
  use Mesh;

  class fr_mesh_c : mesh_c
  {
    const nVars : int;
    var solOrder : int;

    // Domains
    var  xyzSP_d : domain(rank=2, idxType=int);     // {1..nSPs, 1..nDims}
    var  xyzFP_d : domain(rank=2, idxType=int);     // {1..nSPs, 1..nDims}

    var  solSP_d : domain(rank=2, idxType=int);     // {1..nSPs, 1..nVars}
    var  solFP_d : domain(rank=3, idxType=int);     // {1..nFPs, 1..2, 1..nVars}
    var  flxFP_d : domain(rank=2, idxType=int);     // {1..nFPs, 1..nVars}
    var dSolSP_d : domain(rank=3, idxType=int);     // {1..nSPs, 1..nVars, 1..nVars}
    var dSolFP_d : domain(rank=4, idxType=int);     // {1..nFPs, 1..2, 1..nVars, 1..nVars}

    var cellSPidx_d : domain(rank=2, idxType=int);  // {1..nCells, 1..2}
    var faceFPidx_d : domain(rank=2, idxType=int);  // {1..nFaces, 1..2}

    // FR solver variables
    var xyzSP : [xyzSP_d] real;
    var xyzFP : [xyzFP_d] real;

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

      cellSPidx_d = {1..this.nCells, 1..2};
      faceFPidx_d = {1..this.nFaces, 1..2};

      for cell in this.cellList_d
      {
        cellSPidx[cell, 1] = nSPs+1;
        cellSPidx[cell, 2] = n_cell_sps(this.cellList(cell).elemTopo(), this.solOrder);

        nSPs += cellSPidx[cell, 2];

        xyzSP_d  = {1..nSPs, 1..this.nDims};
        solSP_d  = {1..nSPs, 1..this.nVars};
        dSolSP_d = {1..nSPs, 1..this.nVars, 1..this.nDims};

      }

      for face in this.faceList_d
      {
        faceFPidx[face, 1] = nFPs+1;
        faceFPidx[face, 2] = n_face_fps(this.faceList(face).elemTopo(), this.solOrder);

        nFPs += faceFPidx[face, 2];

        xyzFP_d  = {1..nFPs, 1..this.nDims};
        solFP_d  = {1..nFPs, 1..2, 1..this.nVars};
        flxFP_d  = {1..nFPs, 1..this.nVars};
        dSolFP_d = {1..nFPs, 1..2, 1..this.nVars, 1..this.nDims};
      }
    }
  }

  proc n_cell_sps(in elemTopo : int, in solOrder) : int
  {
    use Parameters.Mesh;

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
    use Parameters.Mesh;

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
    use Parameters.Mesh;

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

    writeln("Test 1: Random 1D mesh - Gmesh:");
    writeln(test_gmesh2);
    writeln();

    // Get number of physical dimensions from mesh or input
    var test_mesh = new unmanaged mesh_c(nDims=1);
    test_mesh.import_gmesh2(test_gmesh2);

    writeln("Test 2: Random 1D mesh - Gmesh => Native:");
    writeln(test_mesh);
    writeln();

    writeln("Test 3: Cell and Face counts by topology:");
    writeln("Cell counts:");
    writeln(test_mesh.cell_count());
    writeln("Face counts:");
    writeln(test_mesh.face_count());
    writeln();

    // Get number of physical dimensions from mesh or input
    var test_frmesh = new unmanaged fr_mesh_c(nDims=1, nVars=3, solOrder=2);
    test_frmesh.import_gmesh2(test_gmesh2);
    test_frmesh.allocate_vars();

    writeln("Test 4: Initialize FR variables");
    writeln(test_frmesh);
    writeln();
  }
}
