prototype module Gmesh
{
  use Random;
  use UnitTest;

  // Planned supported mesh formats
  //
  // - 1D Internal generated mesh  (Frei 1.0)
  // - 1D/2D Gmesh version 2/4     (Frei 2.0)
  // - 2D Internal generated Gmesh (Frei 2.x)
  // - 3D Gmesh version 2/4        (Frei 3.0)
  // - 3D CGNS                     (Frei 4.0)
  // - 2D Internal generated CGNS  (Frei 4.x)
  // - 3D Internal generated CGNS  (Frei 4.x)
  // - 3D Internal generated Gmesh (Frei 4.x)
  //
  // It would probably be better to have a separate module for each IO standard (CGNS, Gmesh, Plot3D, etc). At the
  //    moment shit's just pilled up here for compactness since most of it is just general concepts and sketches.
  //
  // The idea behind having a separate class for each IO format it to isolate the conversion from the solver mesh format
  //    to the IO format from the reading and writing of the IO files. It might waste to much memory for just code
  //    organization's sake and should be reviewed in the future.

  class gmesh2_c
  {
    var nodes_d    : domain(rank=2, idxType=int);
    var elements_d : domain(rank=1, idxType=int);
    var families_d : domain(rank=1, idxType=int);

    var nodes    : [nodes_d] real;
    var elements : [elements_d] gmesh_element_r;
    var families : [families_d] gmesh_family_r;

    proc random1D(nCells : int, xMin : real = -1.0, xMax : real= 1.0)
    {
      use Random;
      use Parameters.ParamTest;
      use Parameters.ParamGmesh;
      import Sort.sort;

      var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

      // Allocate mesh elements
      this.nodes_d = {1..nCells+1, 1..3};
      this.elements_d = {1..nCells+2};
      this.families_d = {1..3};

      var nodePermutation : [this.nodes.domain.dim(0)] int;
      var elemPermutation : [this.elements.domain] int;
      var x = this.nodes[..,1];
      var cells : [1..nCells, 1..2] int;

      // Get random values from xMin to xMax
      x[1] = 0;
      x[2] = 1;
      randStreamSeeded.fillRandom(x[3..]);
      x = x*(xMax-xMin) + xMin;

      // Sort nodes to build elements and generate a random permutation
      sort(x);
      permutation(nodePermutation, RANDOM_SEED);
      permutation(elemPermutation, RANDOM_SEED);

      // Fill element list with non overlapping elements oriented from left to right
      for i in 1..nCells
      {
        cells[i,1] = i;
        cells[i,2] = i+1;
      }

      // Commit values to object

      // Set the boundary and internal families
      this.families[1].tag  = 1;
      this.families[1].nDim = 1;
      this.families[1].name = "flow";
      this.families[2].tag  = 2;
      this.families[2].nDim = 0;
      this.families[2].name = "left";
      this.families[3].tag  = 3;
      this.families[3].nDim = 0;
      this.families[3].name = "right";

      // Fill node list in randomized order
      for i in this.nodes.domain.dim(0) do
        this.nodes[nodePermutation[i],1] = x[i];

      // Fill element list with the internal elements in random order
      for i in this.elements.domain.dim(0).expand(-1)
      {
        this.elements[elemPermutation[i]].elemType = GMESH_LIN_2;
        this.elements[elemPermutation[i]].setNodes();
        this.elements[elemPermutation[i]].tags_d = {1..1};
        this.elements[elemPermutation[i]].tags[1] = 1; // Family
        this.elements[elemPermutation[i]].nodes[1] = nodePermutation[cells[i-1,1]];
        this.elements[elemPermutation[i]].nodes[2] = nodePermutation[cells[i-1,2]];
      }

      // Add left boundary point to elements list
      this.elements[elemPermutation[1]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[1]].setNodes();
      this.elements[elemPermutation[1]].tags_d = {1..1};
      this.elements[elemPermutation[1]].tags[1] = 2; // Family
      this.elements[elemPermutation[1]].nodes[1] = nodePermutation[1];

      // Add right boundary point to elements list
      this.elements[elemPermutation[nCells+2]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[nCells+2]].setNodes();
      this.elements[elemPermutation[nCells+2]].tags_d = {1..1};
      this.elements[elemPermutation[nCells+2]].tags[1] = 3; // Family
      this.elements[elemPermutation[nCells+2]].nodes[1] = nodePermutation[this.nodes.domain.dim(0).high];
    }

    proc uniform1D(nCells : int, xMin: real = -1.0, xMax: real = 1.0)
    {
      use Parameters.ParamGmesh;
      import Sort.sort;

      // Allocate mesh elements
      this.nodes_d = {1..nCells+1, 1..3};
      this.elements_d = {1..nCells+2};
      this.families_d = {1..3};

      var x = this.nodes[..,1];
      var cells : [1..nCells, 1..2] int;

      // Get uniform nodes from xMin to xMax
      var step : real = (xMax - xMin)/nCells;
      for node in this.nodes.domain.dim(0) do
        x[node] = xMin + (node-1)*step;

      // Fill element list with non overlapping elements oriented from left to right
      for i in 1..nCells
      {
        cells[i,1] = i;
        cells[i,2] = i+1;
      }

      // Commit values to object

      // Set the boundary and internal families
      this.families[1].tag  = 1;
      this.families[1].nDim = 1;
      this.families[1].name = "flow";
      this.families[2].tag  = 2;
      this.families[2].nDim = 0;
      this.families[2].name = "left";
      this.families[3].tag  = 3;
      this.families[3].nDim = 0;
      this.families[3].name = "right";

      // Fill node list
      for i in this.nodes.domain.dim(0) do
        this.nodes[i,1] = x[i];

      // Fill element list with the internal elements
      for i in this.elements.domain.dim(0).expand(-1)
      {
        this.elements[i].elemType = GMESH_LIN_2;
        this.elements[i].setNodes();
        this.elements[i].tags_d = {1..1};
        this.elements[i].tags[1] = 1; // Family
        this.elements[i].nodes[1] = cells[i-1,1];
        this.elements[i].nodes[2] = cells[i-1,2];
      }

      // Add left boundary point to elements list
      this.elements[1].elemType = GMESH_PNT_1;
      this.elements[1].setNodes();
      this.elements[1].tags_d = {1..1};
      this.elements[1].tags[1] = 2; // Family
      this.elements[1].nodes[1] = 1;

      // Add right boundary point to elements list
      this.elements[nCells+2].elemType = GMESH_PNT_1;
      this.elements[nCells+2].setNodes();
      this.elements[nCells+2].tags_d = {1..1};
      this.elements[nCells+2].tags[1] = 3; // Family
      this.elements[nCells+2].nodes[1] = this.nodes.domain.dim(0).high;
    }

    proc read_gmesh_file() {}
    proc write_gmesh_file() {}
  }

  class gmesh4_c
  {
    var nNodes : int;
    var nElements : int;
    var nFamilies : int;

    var nodes    : [1..nNodes, 1..3] real;
    var elements : [1..nElements] gmesh_element_r;
    var families : [1..nFamilies] gmesh_family_r;

    proc random1D(xMin: real, xMax: real)
    {
      use Random;
      use Parameters.ParamTest;
      use Parameters.ParamGmesh;
      import Sort.sort;

      var randStreamSeeded = new RandomStream(real, RANDOM_SEED);
      var x = this.nodes[..,1];
      var cells : [1..this.nElements-2, 1..2] int;
      var nodePermutation : [1..this.nNodes] int;
      var elemPermutation : [1..this.nElements] int;

      // Get random values from xMin to xMax
      x[1] = 0;
      x[2] = 1;
      randStreamSeeded.fillRandom(x[3..]);
      x = x*(xMax-xMin) + xMin;

      // Sort nodes to build elements and generate a random permutation
      sort(x);
      permutation(nodePermutation, RANDOM_SEED);
      permutation(elemPermutation, RANDOM_SEED);

      // Fill element list with non overlapping elements oriented from left to right
      for i in 1..this.nElements-2 {
        cells[i,1] = i;
        cells[i,2] = i+1;
      }

      // Commit values to object

      // Set the boundary and internal families
      this.families[1].nDim = 1;
      this.families[1].name = "flow";
      this.families[2].nDim = 0;
      this.families[2].name = "left";
      this.families[3].nDim = 0;
      this.families[3].name = "right";

      // Fill node list in randomized order
      for i in 1..this.nNodes do
        this.nodes[nodePermutation[i],1] = x[i];

      // Fill element list with the internal elements in random order
      for i in 2..this.nElements-1 {
        this.elements[elemPermutation[i]].elemType = GMESH_LIN_2;
        this.elements[elemPermutation[i]].setNodes();
        this.elements[elemPermutation[i]].tags[1] = 1; // Family
        this.elements[elemPermutation[i]].nodes[1] = nodePermutation[cells[i-1,1]];
        this.elements[elemPermutation[i]].nodes[2] = nodePermutation[cells[i-1,2]];
      }

      // Add left boundary point to elements list
      this.elements[elemPermutation[1]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[1]].setNodes();
      this.elements[elemPermutation[1]].tags[1] = 2; // Family
      this.elements[elemPermutation[1]].nodes[1] = nodePermutation[1];

      // Add right boundary point to elements list
      this.elements[elemPermutation[this.nElements]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[this.nElements]].setNodes();
      this.elements[elemPermutation[this.nElements]].tags[1] = 3; // Family
      this.elements[elemPermutation[this.nElements]].nodes[1] = nodePermutation[this.nNodes];
    }

    proc read_gmesh_file() {}
    proc write_gmesh_file() {}
  }

  record gmesh_element_r
  {
    var elemType : int;
    var tags_d   : domain(rank=1, idxType=int);
    var tags     : [tags_d] int;
    var nodes_d  : domain(rank=1, idxType=int);
    var nodes    : [nodes_d] int;

    proc setNodes()
    {
      use Parameters.ParamGmesh;

      select this.elemType {
        // Point
        when GMESH_PNT_1    do this.nodes_d = {1..1   };
        // Line
        when GMESH_LIN_2    do this.nodes_d = {1..2   };
        when GMESH_LIN_3    do this.nodes_d = {1..3   };
        when GMESH_LIN_4    do this.nodes_d = {1..4   };
        when GMESH_LIN_5    do this.nodes_d = {1..5   };
        when GMESH_LIN_6    do this.nodes_d = {1..6   };
        when GMESH_LIN_7    do this.nodes_d = {1..7   };
        when GMESH_LIN_8    do this.nodes_d = {1..8   };
        when GMESH_LIN_9    do this.nodes_d = {1..9   };
        when GMESH_LIN_10   do this.nodes_d = {1..10  };
        // Lagrange Triangles
        when GMESH_TRI_3    do this.nodes_d = {1..3   };
        when GMESH_TRI_6    do this.nodes_d = {1..6   };
        when GMESH_TRI_10   do this.nodes_d = {1..10  };
        when GMESH_TRI_15   do this.nodes_d = {1..15  };
        when GMESH_TRI_21   do this.nodes_d = {1..21  };
        when GMESH_TRI_28   do this.nodes_d = {1..28  };
        when GMESH_TRI_36   do this.nodes_d = {1..36  };
        when GMESH_TRI_45   do this.nodes_d = {1..45  };
        when GMESH_TRI_55   do this.nodes_d = {1..55  };
        // Lagrange Quadrilaterals
        when GMESH_QUA_4    do this.nodes_d = {1..4   };
        when GMESH_QUA_9    do this.nodes_d = {1..9   };
        when GMESH_QUA_16   do this.nodes_d = {1..16  };
        when GMESH_QUA_25   do this.nodes_d = {1..25  };
        when GMESH_QUA_36   do this.nodes_d = {1..36  };
        when GMESH_QUA_49   do this.nodes_d = {1..49  };
        when GMESH_QUA_64   do this.nodes_d = {1..64  };
        when GMESH_QUA_81   do this.nodes_d = {1..81  };
        when GMESH_QUA_100  do this.nodes_d = {1..100 };
        // Lagrange Tetrahedra
        when GMESH_TET_4    do this.nodes_d = {1..4   };
        when GMESH_TET_10   do this.nodes_d = {1..10  };
        when GMESH_TET_20   do this.nodes_d = {1..20  };
        when GMESH_TET_35   do this.nodes_d = {1..35  };
        when GMESH_TET_56   do this.nodes_d = {1..56  };
        when GMESH_TET_84   do this.nodes_d = {1..84  };
        when GMESH_TET_120  do this.nodes_d = {1..120 };
        when GMESH_TET_165  do this.nodes_d = {1..165 };
        when GMESH_TET_220  do this.nodes_d = {1..220 };
        // Lagrange Pyramids
        when GMESH_PYR_5    do this.nodes_d = {1..5   };
        when GMESH_PYR_14   do this.nodes_d = {1..14  };
        when GMESH_PYR_30   do this.nodes_d = {1..30  };
        when GMESH_PYR_55   do this.nodes_d = {1..55  };
        when GMESH_PYR_91   do this.nodes_d = {1..91  };
        when GMESH_PYR_140  do this.nodes_d = {1..140 };
        when GMESH_PYR_204  do this.nodes_d = {1..204 };
        when GMESH_PYR_285  do this.nodes_d = {1..285 };
        when GMESH_PYR_385  do this.nodes_d = {1..385 };
        // Lagrange Prisms
        when GMESH_PRI_6    do this.nodes_d = {1..6   };
        when GMESH_PRI_18   do this.nodes_d = {1..18  };
        when GMESH_PRI_40   do this.nodes_d = {1..40  };
        when GMESH_PRI_75   do this.nodes_d = {1..75  };
        when GMESH_PRI_126  do this.nodes_d = {1..126 };
        when GMESH_PRI_196  do this.nodes_d = {1..196 };
        when GMESH_PRI_288  do this.nodes_d = {1..288 };
        when GMESH_PRI_405  do this.nodes_d = {1..405 };
        when GMESH_PRI_550  do this.nodes_d = {1..550 };
        // Lagrange Hexahedra
        when GMESH_HEX_8    do this.nodes_d = {1..8   };
        when GMESH_HEX_27   do this.nodes_d = {1..27  };
        when GMESH_HEX_64   do this.nodes_d = {1..64  };
        when GMESH_HEX_125  do this.nodes_d = {1..125 };
        when GMESH_HEX_216  do this.nodes_d = {1..216 };
        when GMESH_HEX_343  do this.nodes_d = {1..343 };
        when GMESH_HEX_512  do this.nodes_d = {1..512 };
        when GMESH_HEX_729  do this.nodes_d = {1..729 };
        when GMESH_HEX_1000 do this.nodes_d = {1..1000};
      }
    }
  }

  record gmesh_family_r
  {
    var tag  : int;
    var nDim : int;
    var name : string;
  }

  proc main()
  {
    {
      writeln("Test 1: Random 1D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.random1D(nCells=6, xMin=-1, xMax=2);
      writeln(test_gmesh2);
      writeln();
    }

    {
      writeln("Test 2: Uniform 1D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.uniform1D(nCells=6, xMin=-1, xMax=2);
      writeln(test_gmesh2);
      writeln();
    }
  }
}
