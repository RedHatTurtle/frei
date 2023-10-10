module Gmesh
{
  use Map;
  use Random;
  use UnitTest;

  config const testMesh : string = "test-cases/gmesh-1d-test.msh";

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
      for node in this.nodes.domain.dim(0) do
        x[node] = ((nCells+1-node)*xMin + (node-1)*xMax)/nCells;

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
        this.elements[i].tags[1] = 1; // Family Idx
        this.elements[i].nodes[1] = cells[i-1,1];
        this.elements[i].nodes[2] = cells[i-1,2];
      }

      // Add left boundary point to elements list
      this.elements[1].elemType = GMESH_PNT_1;
      this.elements[1].setNodes();
      this.elements[1].tags_d = {1..1};
      this.elements[1].tags[1] = 2; // Family Idx
      this.elements[1].nodes[1] = 1;

      // Add right boundary point to elements list
      this.elements[nCells+2].elemType = GMESH_PNT_1;
      this.elements[nCells+2].setNodes();
      this.elements[nCells+2].tags_d = {1..1};
      this.elements[nCells+2].tags[1] = 3; // Family Idx
      this.elements[nCells+2].nodes[1] = this.nodes.domain.dim(0).high;
    }

    proc read_gmesh_file(meshFileName : string)
    {
      use IO;
      use OS;

      enum section {Main, MeshFormat, PhysicalNames, Nodes, Elements, Periodic, NodeData, ElementData,
                                ElementNodeData, InterpolationScheme};

      var meshFile : file;

      // Open the mesh file for reading only
      try {
        writeln();
        writeln("Opening Gmesh2 formated mesh file");
        meshFile = open(meshFileName, ioMode.r);
      } catch e : FileNotFoundError {
        writeln("Critical Error: Mesh file not found.");
      } catch {
        writeln("Unknown Error opening mesh file.");
      }

      // Open file reader to gmesh file
      var meshFileReader = try! meshFile.reader();

      // Set initial state
      var state = section.Main;
      var nLines : int = 0;
      var lineIdx : int = 1;

      // Start reading and parsing file
      for line in meshFileReader.lines()
      {
        select state
        {
          when section.Main
          {
            // Update state with the section starting at this line
            select line
            {
              when "$MeshFormat\n" do
              {
                write("Verifying Mesh Format  ...");
                state = section.MeshFormat;
              }
              when "$PhysicalNames\n"
              {
                write("Reading Physical Names ...");
                state = section.PhysicalNames;
              }
              when "$Nodes\n"
              {
                write("Reading Mesh Nodes     ...");
                state = section.Nodes;
              }
              when "$Elements\n"
              {
                write("Reading Mesh Elements  ...");
                state = section.Elements;
              }
              when "$Periodic\n" do
                state = section.Periodic;
              when "$NodeData\n" do
                state = section.NodeData;
              when "$ElementData\n" do
                state = section.ElementData;
              when "$ElementNodeData\n" do
                state = section.ElementNodeData;
              when "$InterpolationScheme\n" do
                state = section.InterpolationScheme;
              otherwise do
                writeln("Unexpected line on main section");
            }
          }
          when section.MeshFormat
          {
            if line == "$EndMeshFormat\n"
            {
                // Reset state to main section
                writeln(" done");
                state = section.Main;
                nLines = 0;
            }
            else if line != "2.2 0 8\n"
            {
              // Validate if gmesh format is compatible with this class
              writeln("Unsuported gmesh format version");
              writeln("Expected '2.2 0 8'");
            }
          }
          when section.PhysicalNames
         {
            if line == "$EndPhysicalNames\n"
            {
                // Reset state to main section
                writeln(" done");
                state = section.Main;
                nLines = 0;
                lineIdx = 1;
            }
            else if nLines == 0
            {
              // If it's the first line of the section get the number of Physical Namas and allocate families
              try! {
                nLines = line : int;
              }
              catch {
                writeln("Failed to cast `Physical Names Count` to int");
              }
              this.families_d = {1..nLines};
            }
            else
            {
              var physicalDim  : int;
              var physicalTag  : int;
              var physicalName : string;
              var valueIdx : int = 1;

              for value in line.split()
              {
                 if valueIdx == 1 then
                    try! this.families[lineIdx].nDim = value:int;
                 if valueIdx == 2 then
                    try! this.families[lineIdx].tag = value:int;
                 if valueIdx == 3 then
                    try! this.families[lineIdx].name = value.strip("\"");

                 valueIdx += 1;
              }

              lineIdx += 1;
            }
          }
          when section.Nodes
          {
            if line == "$EndNodes\n"
            {
                // Reset state to main section
                writeln(" done");
                state = section.Main;
                nLines = 0;
                lineIdx = 1;
            }
            else if nLines == 0
            {
              // If it's the first line of the section get the number of Physical Namas and allocate families
              try! {
                nLines = line : int;
              }
              catch {
                writeln("Failed to cast `Node Count` to int");
              }
              this.nodes_d = {1..nLines, 1..3};
            }
            else
            {
              var nodeIdx   : int;
              var nodeCoord : [1..3] real;

              var valueIdx  : int = 1;
              for value in line.split()
              {
                 if valueIdx == 1 then
                    try! nodeIdx = value : int;
                 if valueIdx == 2 then
                    try! nodeCoord[1] = value : real;
                 if valueIdx == 3 then
                    try! nodeCoord[2] = value : real;
                 if valueIdx == 4 then
                    try! nodeCoord[3] = value : real;

                 valueIdx += 1;
              }

              this.nodes[nodeIdx, 1..3] = nodeCoord[1..3];

              lineIdx += 1;
            }
          }
          when section.Elements
          {
            if line == "$EndElements\n"
            {
                // Reset state to main section
                writeln(" done");
                state = section.Main;
                nLines = 0;
                lineIdx = 1;
            }
            else if nLines == 0
            {
              // If it's the first line of the section get the number of Physical Namas and allocate families
              try! {
                nLines = line : int;
              }
              catch {
                writeln("Failed to cast `Number of Elements` to int");
              }
              this.elements_d = {1..nLines};
            }
            else
            {
              var elemIdx   : int;

              var valueIdx  : int = 1;
              for value in line.split()
              {
                 // Element Number
                 if valueIdx == 1 then
                    try! elemIdx = value : int;
                 // Element Type
                 if valueIdx == 2
                 {
                    try! this.elements[elemIdx].elemType = value : int;
                    this.elements[elemIdx].setNodes();
                 }
                 // Number of Tags
                 if valueIdx == 3 then
                    try! this.elements[elemIdx].tags_d = {1..value:int};
                 // Tags
                 if ((valueIdx > 3) && (valueIdx < 4+this.elements[elemIdx].tags.domain.dim(0).high)) then
                    try! this.elements[elemIdx].tags[valueIdx-3] = value : int;
                 // Node List
                 if valueIdx >= 4+this.elements[elemIdx].tags.domain.dim(0).high then
                   try! this.elements[elemIdx].nodes[valueIdx-(3+this.elements[elemIdx].tags_d.dim(0).high)] = value : int;

                 valueIdx += 1;
              }

              lineIdx += 1;
            }
          }
          when section.Periodic
          {
            if line == "$EndPeriodic\n" then
              state = section.Main;
          }
          when section.NodeData
          {
            if line == "$EndNodeData\n" then
              state = section.Main;
          }
          when section.ElementData
          {
            if line == "$EndElementData\n" then
              state = section.Main;
          }
          when section.ElementNodeData
          {
            if line == "$EndElementNodeData\n" then
              state = section.Main;
          }
          when section.InterpolationScheme
          {
            if line == "$EndInterpolationScheme\n" then
              state = section.Main;
          }
        }
      }

      writeln();
      writeln("Families found:");
      for family in this.families do
        writef("   %1i-D Elements, Tag: %2i, Name: %s\n", family.nDim, family.tag, family.name);
      writeln();
    }

    proc write_gmesh_file() {}

    proc mesh_dimension() do return max reduce this.families[..].nDim;

    proc node_cnt() : int do return this.nodes.domain.dim(0).size;

    proc face_cnt() : int
    {
      var nFaces : int = 0;

      for elem in this.elements
      {
        if elem.elemDim() == this.mesh_dimension() then
          nFaces += elem.elemFaces();
        else
          nFaces += 1;
      }

      return nFaces/2;
    }

    proc cell_cnt() : int
    {
      var cellCnt : int = 0;

      // Get cell element dimension
      const cellDim : int = this.mesh_dimension();

      // Iterate through elements searching for cells
      for elem in this.elements do
        if elem.elemDim() == cellDim then
          cellCnt += 1;

      return cellCnt;
    }

    proc boco_cnt() : int
    {
      var bocoCnt : int = 0;

      // Get boco element dimension
      const bocoDim : int = this.mesh_dimension()-1;

      // Iterate through elements searching for boundary faces
      for elem in this.elements do
        if elem.elemDim() == bocoDim then
          bocoCnt += 1;

      return bocoCnt;
    }

    proc faml_cnt() : int do return this.families.domain.dim(0).size;

    proc face_topo_cnt() : map
    {
      use Map;

      var faceCnt : map(keyType = int, valType=int);

      for elem in this.elements
      {
        if elem.elemDim() == this.mesh_dimension()
        {
          // Loop though and count each face of this cell
          for faceTopo in elem.elemFaceTopo()
          {
            if !faceCnt.contains(faceTopo) then
              faceCnt.add(faceTopo, 0);

            try! faceCnt[faceTopo] += 1 ;
          }
        }
        else if elem.elemDim() == this.mesh_dimension()-1
        {
          // This is a boundary face, add it to the count
          if !faceCnt.contains(elem.elemTopo()) then
            faceCnt.add(elem.elemTopo(), 0);

          try! faceCnt[elem.elemTopo()] += 1 ;
        }
      }

      for faceTopo in faceCnt.keys() do
        try! faceCnt[faceTopo] /= 2 ;

      return faceCnt;
    }

    proc cell_topo_cnt() : map
    {
      use Map;

      var cellCnt : map(keyType = int, valType=int);

      // Get cell element dimension
      const cellDim : int = this.mesh_dimension();

      // Iterate through elements searching for cells
      for elem in this.elements do
        if elem.elemDim() == cellDim then
        {
          // Add topology to map if not already there
          if !cellCnt.contains(elem.elemTopo()) then
            cellCnt.add(elem.elemTopo(), 0);

          // Count the cell
          try! cellCnt[elem.elemTopo()] += 1 ;
        }

      return cellCnt;
    }

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

    proc elemDim()   : int do return elem_dim(this.elemType);
    proc elemTopo()  : int do return elem_topo(this.elemType);
    proc elemFaces() : int do return elem_faces(this.elemTopo());
    proc elemFaceTopo() : [] int do return elem_face_topo(this.elemTopo());

    proc ref setNodes()
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

  proc elem_dim(const ref elemType : int) : int
  {
    use Parameters.ParamGmesh;

    select elemType {
      when GMESH_PNT_1    do return 0;
      when GMESH_LIN_2    do return 1;
      when GMESH_LIN_3    do return 1;
      when GMESH_LIN_4    do return 1;
      when GMESH_LIN_5    do return 1;
      when GMESH_LIN_6    do return 1;
      when GMESH_LIN_7    do return 1;
      when GMESH_LIN_8    do return 1;
      when GMESH_LIN_9    do return 1;
      when GMESH_LIN_10   do return 1;
      when GMESH_LIN_11   do return 1;
      when GMESH_TRI_3    do return 2;
      when GMESH_TRI_6    do return 2;
      when GMESH_TRI_10   do return 2;
      when GMESH_TRI_15   do return 2;
      when GMESH_TRI_21   do return 2;
      when GMESH_TRI_28   do return 2;
      when GMESH_TRI_36   do return 2;
      when GMESH_TRI_45   do return 2;
      when GMESH_TRI_55   do return 2;
      when GMESH_TRI_66   do return 2;
      when GMESH_QUA_4    do return 2;
      when GMESH_QUA_9    do return 2;
      when GMESH_QUA_16   do return 2;
      when GMESH_QUA_25   do return 2;
      when GMESH_QUA_36   do return 2;
      when GMESH_QUA_49   do return 2;
      when GMESH_QUA_64   do return 2;
      when GMESH_QUA_81   do return 2;
      when GMESH_QUA_100  do return 2;
      when GMESH_QUA_121  do return 2;
      when GMESH_PYR_5    do return 3;
      when GMESH_PYR_14   do return 3;
      when GMESH_PYR_30   do return 3;
      when GMESH_PYR_55   do return 3;
      when GMESH_PYR_91   do return 3;
      when GMESH_PYR_140  do return 3;
      when GMESH_PYR_204  do return 3;
      when GMESH_PYR_285  do return 3;
      when GMESH_PYR_385  do return 3;
      when GMESH_PRI_6    do return 3;
      when GMESH_PRI_18   do return 3;
      when GMESH_PRI_40   do return 3;
      when GMESH_PRI_75   do return 3;
      when GMESH_PRI_126  do return 3;
      when GMESH_PRI_196  do return 3;
      when GMESH_PRI_288  do return 3;
      when GMESH_PRI_405  do return 3;
      when GMESH_PRI_550  do return 3;
      when GMESH_HEX_8    do return 3;
      when GMESH_HEX_27   do return 3;
      when GMESH_HEX_64   do return 3;
      when GMESH_HEX_125  do return 3;
      when GMESH_HEX_216  do return 3;
      when GMESH_HEX_343  do return 3;
      when GMESH_HEX_512  do return 3;
      when GMESH_HEX_729  do return 3;
      when GMESH_HEX_1000 do return 3;
      otherwise return -1;
    }
  }

  proc elem_topo(const ref elemType : int) : int
  {
    use Parameters.ParamGmesh;

    select elemType {
      when GMESH_PNT_1   do return GMESH_PNT;
      when GMESH_LIN_2   do return GMESH_LIN;
      when GMESH_LIN_3   do return GMESH_LIN;
      when GMESH_LIN_4   do return GMESH_LIN;
      when GMESH_LIN_5   do return GMESH_LIN;
      when GMESH_TRI_3   do return GMESH_TRI;
      when GMESH_TRI_6   do return GMESH_TRI;
      when GMESH_TRI_10  do return GMESH_TRI;
      when GMESH_TRI_15  do return GMESH_TRI;
      when GMESH_QUA_4   do return GMESH_QUA;
      when GMESH_QUA_9   do return GMESH_QUA;
      when GMESH_QUA_16  do return GMESH_QUA;
      when GMESH_QUA_25  do return GMESH_QUA;
      when GMESH_TET_4   do return GMESH_TET;
      when GMESH_TET_10  do return GMESH_TET;
      when GMESH_TET_20  do return GMESH_TET;
      when GMESH_TET_35  do return GMESH_TET;
      when GMESH_PYR_5   do return GMESH_PYR;
      when GMESH_PYR_14  do return GMESH_PYR;
      when GMESH_PYR_30  do return GMESH_PYR;
      when GMESH_PYR_55  do return GMESH_PYR;
      when GMESH_PRI_6   do return GMESH_PRI;
      when GMESH_PRI_18  do return GMESH_PRI;
      when GMESH_PRI_40  do return GMESH_PRI;
      when GMESH_PRI_75  do return GMESH_PRI;
      when GMESH_HEX_8   do return GMESH_HEX;
      when GMESH_HEX_27  do return GMESH_HEX;
      when GMESH_HEX_64  do return GMESH_HEX;
      when GMESH_HEX_125 do return GMESH_HEX;
      otherwise return -1;
    }
  }

  proc elem_faces(const ref elemTopo : int) : int
  { // Assuming this mesh element is a cell how many faces does it have
    use Parameters.ParamGmesh;

    select elemTopo {
      when GMESH_PNT do return 0; // Vertex
      when GMESH_LIN do return 2; // Edge
      when GMESH_TRI do return 3; // Triangle
      when GMESH_QUA do return 4; // Quadrilateral
      when GMESH_TET do return 4; // Tetrahedron
      when GMESH_PYR do return 5; // Pyramid
      when GMESH_PRI do return 5; // Prism
      when GMESH_HEX do return 6; // Hexahedron
      otherwise return -1;
    }
  }

  proc elem_face_topo(const ref elemTopo : int) : [] int
  { // Assuming this mesh element is a cell how many faces does it have
    use Parameters.ParamGmesh;

    select elemTopo {
      when GMESH_LIN do return [GMESH_PNT, GMESH_PNT]; // Edge
      when GMESH_TRI do return [GMESH_LIN, GMESH_LIN, GMESH_LIN]; // Triangle
      when GMESH_QUA do return [GMESH_LIN, GMESH_LIN, GMESH_LIN, GMESH_LIN]; // Quadrilateral
      when GMESH_TET do return [GMESH_TRI, GMESH_TRI, GMESH_TRI, GMESH_TRI, ]; // Tetrahedron
      when GMESH_PYR do return [GMESH_TET, GMESH_TRI, GMESH_TRI, GMESH_TRI, GMESH_TRI, ]; // Pyramid
      when GMESH_PRI do return [GMESH_TRI, GMESH_TRI, GMESH_QUA, GMESH_QUA, GMESH_QUA, ]; // Prism
      when GMESH_HEX do return [GMESH_QUA, GMESH_QUA, GMESH_QUA, GMESH_QUA, GMESH_QUA, GMESH_QUA, ]; // Hexahedron
      otherwise return [-1];
    }
  }

  proc main()
  {
    {
      writeln("Test 1: Random 1D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.random1D(nCells=6, xMin=-1, xMax=2);
      writeln("Test Gmesh query functions:");
      writeln("  - Mesh Dim: ", test_gmesh2.mesh_dimension());
      writeln("  - Node Cnt: ", test_gmesh2.node_cnt());
      writeln("  - Face Cnt: ", test_gmesh2.face_cnt());
      writeln("  - Cell Cnt: ", test_gmesh2.cell_cnt());
      writeln("  - Boco Cnt: ", test_gmesh2.boco_cnt());
      writeln("  - Faml Cnt: ", test_gmesh2.faml_cnt());
      writeln("  - Face Topo Cnt: ", test_gmesh2.face_topo_cnt());
      writeln("  - Cell Topo Cnt: ", test_gmesh2.cell_topo_cnt());
      writeln();
      writeln("Dump Gmesh instance:");
      writeln(test_gmesh2);
      writeln();
    }

    {
      writeln("Test 2: Uniform 1D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.uniform1D(nCells=6, xMin=-1, xMax=2);
      writeln("Test Gmesh query functions:");
      writeln("  - Mesh Dim: ", test_gmesh2.mesh_dimension());
      writeln("  - Node Cnt: ", test_gmesh2.node_cnt());
      writeln("  - Face Cnt: ", test_gmesh2.face_cnt());
      writeln("  - Cell Cnt: ", test_gmesh2.cell_cnt());
      writeln("  - Boco Cnt: ", test_gmesh2.boco_cnt());
      writeln("  - Faml Cnt: ", test_gmesh2.faml_cnt());
      writeln("  - Face Topo Cnt: ", test_gmesh2.face_topo_cnt());
      writeln("  - Cell Topo Cnt: ", test_gmesh2.cell_topo_cnt());
      writeln();
      writeln("Dump Gmesh instance:");
      writeln(test_gmesh2);
      writeln();
    }

    {
      writeln("Test 3: Read 2D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.read_gmesh_file(testMesh);
      writeln("Test Gmesh query functions:");
      writeln("  - Mesh Dim: ", test_gmesh2.mesh_dimension());
      writeln("  - Node Cnt: ", test_gmesh2.node_cnt());
      writeln("  - Face Cnt: ", test_gmesh2.face_cnt());
      writeln("  - Cell Cnt: ", test_gmesh2.cell_cnt());
      writeln("  - Boco Cnt: ", test_gmesh2.boco_cnt());
      writeln("  - Faml Cnt: ", test_gmesh2.faml_cnt());
      writeln("  - Face Topo Cnt: ", test_gmesh2.face_topo_cnt());
      writeln("  - Cell Topo Cnt: ", test_gmesh2.cell_topo_cnt());
      writeln();
      writeln("Dump Gmesh instance:");
      writeln(test_gmesh2);
      writeln();
    }
  }
}
