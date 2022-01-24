module Mesh
{
  use Gmesh;
  use Random;
  use UnitTest;
  use Set;

  // Internal mesh structure

  record node_r
  {
    // Arrays
    var xyz : [1..3] real;
  }

  record edge_r
  {
    // Geometric properties
    var elemType : int;

    // Arrays
    var nodes : [1..2] int;

    proc elemTopo() : int
    {
      return elem_topology(this.elemType);
    }
  }

  record face_r
  {
    // Geometric properties
    var elemType : int;

    // Array sizing domains
    var nodes_d : domain(1);
    //var edges_d : domain(1);

    // Arrays
    var nodes : [nodes_d] int;
    //var edges : [edges_d] int;
    var cells : [1..2] int;

    proc elemTopo() : int
    {
      return elem_topology(this.elemType);
    }
  }

  record cell_r
  {
    // Geometric properties
    var elemType : int;

    // Array sizing domains
    var nodes_d : domain(1);
    //var edges_d : domain(1);
    var faces_d : domain(1);

    // Arrays
    var nodes : [nodes_d] int;
    //var edges : [edges_d] int;
    var faces : [faces_d] int;
    var sides : [faces_d] int; // Side of the face this cell is on (1-Left / 2-Right)

    var family : int;

    proc elemTopo() : int
    {
      return elem_topology(this.elemType);
    }
  }

  record boco_r
  {
    // Geometric properties
    var elemType : int;

    // Array sizing domains
    var nodes_d : domain(1);
    //var edges_d : domain(1);

    // Arrays
    var nodes : [nodes_d] int;
    //var edges : [edges_d] int;
    var face : int;

    var family : int;

    proc elemTopo() : int
    {
      return elem_topology(this.elemType);
    }
  }

  record faml_r
  {
    var name : string;
    var nDim : int;

    // Boco definition
    var bocoType, bocoSubType : int;
    var bocoProperties : [1..9] real;
  }

  class mesh_c
  {
    const nDims : int;
    var nNodes : int;
    var nEdges : int;
    var nFaces : int;
    var nCells : int;
    var nBocos : int;
    var nFamls : int;

    // Lists of types of elements in this mesh
    var faceTopos : set(int); // Gmesh element topology/shape. Ex: point, line, triangle, hexahedron
    var faceTypes : set(int); // CGNS element type, element topology + high-order spec. Ex: tetrahedron with 4 or 10 nodes
    var cellTopos : set(int); // Gmesh element topology/shape. Ex: point, line, triangle, hexahedron
    var cellTypes : set(int); // CGNS element type, element topology + high-order spec. Ex: tetrahedron with 4 or 10 nodes

    // Array sizing domains
    var nodeList_d : domain(rank=1, idxType=int);
    var edgeList_d : domain(rank=1, idxType=int);
    var faceList_d : domain(rank=1, idxType=int);
    var cellList_d : domain(rank=1, idxType=int);
    var bocoList_d : domain(rank=1, idxType=int);
    var famlList_d : domain(rank=1, idxType=int);

    // Arrays
    var nodeList : [nodeList_d] node_r;
    var edgeList : [edgeList_d] edge_r;
    var faceList : [faceList_d] face_r;
    var cellList : [cellList_d] cell_r;
    var bocoList : [bocoList_d] boco_r;
    var famlList : [famlList_d] faml_r;

    proc import_gmesh2(gmesh : gmesh2_c)
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

        // Search families for tag and get the elements dimension from family.
        // Gmesh family tags are not globally unique, they are unique within families with the same number of spatial dimensions.
        for famlIdx in gmesh.families.domain do
          if gmesh.elements[element].tags[1] == gmesh.families[famlIdx].tag
          && gmesh.elements[element].elemDim() == gmesh.families[famlIdx].nDim
          {
            elemDim = gmesh.families[famlIdx].nDim;
            elemFamlIdx = famlIdx;
          }

        // Add to cell or boco list depending on the element dimension
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
      this.face_list_builder();

      // Build the sets of cell and face element types and topologies present in this mesh
      this.elem_set_builder();
    }

    proc face_list_builder()
    {
      use Parameters.ParamMesh;
      use Parameters.ParamGmesh;

      // Build face list
      var faceNodes : [1..6] 4*int;
      var faceMap_d : domain(4*int);
      var faceMap : [faceMap_d] int;

      // Add all boundaries to the face map
      for boco in this.bocoList_d
      {
        // Each boundary defines a face. Build the face nodes tuple based on the element type.
        select elem_topology(this.bocoList[boco].elemType)
        {
          when TOPO_NODE do faceNodes[1] = (this.bocoList[boco].nodes[1], 0, 0, 0);
          when TOPO_LINE do faceNodes[1] = (this.bocoList[boco].nodes[1], this.bocoList[boco].nodes[2], 0, 0);
          when TOPO_TRIA do faceNodes[1] = (this.bocoList[boco].nodes[1], this.bocoList[boco].nodes[2],
                                            this.bocoList[boco].nodes[3], 0);
          when TOPO_QUAD do faceNodes[1] = (this.bocoList[boco].nodes[1], this.bocoList[boco].nodes[2],
                                            this.bocoList[boco].nodes[3], this.bocoList[boco].nodes[4]);
          otherwise {}
        }

        // Create the face element on the face list
        // Increment cell count with new element
        this.nFaces += 1;
        // Resize domain to expand the array
        this.faceList_d = {1..this.nFaces};
        // Fill up the face properties. Maybe this should be an initializer?
        select elem_topology(this.bocoList[boco].elemType)
        {
          when TOPO_NODE
          {
            this.faceList[this.nFaces].elemType = this.bocoList[boco].elemType;
            this.faceList[this.nFaces].nodes_d = {1..1};
            this.faceList[this.nFaces].nodes[1] = faceNodes[1][0];
          }
          when TOPO_LINE
          {
            this.faceList[this.nFaces].elemType = this.bocoList[boco].elemType;
            this.faceList[this.nFaces].nodes_d = {1..2};
            this.faceList[this.nFaces].nodes[1] = faceNodes[1][0];
            this.faceList[this.nFaces].nodes[2] = faceNodes[1][1];
          }
          when TOPO_TRIA
          {
            this.faceList[this.nFaces].elemType = this.bocoList[boco].elemType;
            this.faceList[this.nFaces].nodes_d = {1..3};
            this.faceList[this.nFaces].nodes[1] = faceNodes[1][0];
            this.faceList[this.nFaces].nodes[2] = faceNodes[1][1];
            this.faceList[this.nFaces].nodes[3] = faceNodes[1][2];
          }
          when TOPO_QUAD
          {
            this.faceList[this.nFaces].elemType = this.bocoList[boco].elemType;
            this.faceList[this.nFaces].nodes_d = {1..4};
            this.faceList[this.nFaces].nodes[1] = faceNodes[1][0];
            this.faceList[this.nFaces].nodes[2] = faceNodes[1][1];
            this.faceList[this.nFaces].nodes[3] = faceNodes[1][2];
            this.faceList[this.nFaces].nodes[4] = faceNodes[1][3];
          }
          otherwise {}
        }

        // Fill the right side neighbor ID. Bocos have negative indexes so they can be easily distinguished
        this.faceList[this.nFaces].cells[2] = -boco;
        this.bocoList[boco].face = this.nFaces;

        // Store the face ID with for when we find the left side neighbor
        faceMap_d.add(sort_tuple(faceNodes[1]));
        faceMap[sort_tuple(faceNodes[1])] = this.nFaces;
      }

      // Add faces from cells to the face map and perform the matching
      for cell in this.cellList_d
      {
        // Allocate the local face array
        this.cellList[cell].faces_d = {1..elem_faces(elem_topology(this.cellList[cell].elemType))};

        // Build the face nodes tuple
        select elem_topology(this.cellList[cell].elemType)
        {
          when TOPO_LINE
          {
            faceNodes[1] = (this.cellList[cell].nodes[1], 0, 0, 0);
            faceNodes[2] = (this.cellList[cell].nodes[2], 0, 0, 0);
          }
          when TOPO_TRIA
          {
            faceNodes[1] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[2], 0, 0);
            faceNodes[2] = (this.cellList[cell].nodes[2], this.cellList[cell].nodes[3], 0, 0);
            faceNodes[3] = (this.cellList[cell].nodes[3], this.cellList[cell].nodes[1], 0, 0);
          }
          when TOPO_QUAD
          {
            faceNodes[1] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[2], 0, 0);
            faceNodes[2] = (this.cellList[cell].nodes[2], this.cellList[cell].nodes[3], 0, 0);
            faceNodes[3] = (this.cellList[cell].nodes[3], this.cellList[cell].nodes[4], 0, 0);
            faceNodes[4] = (this.cellList[cell].nodes[4], this.cellList[cell].nodes[1], 0, 0);
          }
          when TOPO_TETR
          {
            faceNodes[1] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[3], this.cellList[cell].nodes[2],
                            0);
            faceNodes[2] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[2], this.cellList[cell].nodes[4],
                            0);
            faceNodes[3] = (this.cellList[cell].nodes[2], this.cellList[cell].nodes[3], this.cellList[cell].nodes[4],
                            0);
            faceNodes[4] = (this.cellList[cell].nodes[3], this.cellList[cell].nodes[1], this.cellList[cell].nodes[4],
                            0);
          }
          when TOPO_PYRA
          {
            faceNodes[1] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[4], this.cellList[cell].nodes[3],
                            this.cellList[cell].nodes[2]);
            faceNodes[2] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[2], this.cellList[cell].nodes[5],
                            0);
            faceNodes[3] = (this.cellList[cell].nodes[2], this.cellList[cell].nodes[3], this.cellList[cell].nodes[5],
                            0);
            faceNodes[4] = (this.cellList[cell].nodes[3], this.cellList[cell].nodes[4], this.cellList[cell].nodes[5],
                            0);
            faceNodes[5] = (this.cellList[cell].nodes[4], this.cellList[cell].nodes[1], this.cellList[cell].nodes[5],
                            0);
          }
          when TOPO_PRIS
          {
            faceNodes[1] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[2], this.cellList[cell].nodes[5],
                            this.cellList[cell].nodes[4]);
            faceNodes[2] = (this.cellList[cell].nodes[2], this.cellList[cell].nodes[3], this.cellList[cell].nodes[6],
                            this.cellList[cell].nodes[5]);
            faceNodes[3] = (this.cellList[cell].nodes[3], this.cellList[cell].nodes[1], this.cellList[cell].nodes[4],
                            this.cellList[cell].nodes[6]);
            faceNodes[4] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[3], this.cellList[cell].nodes[2],
                            0);
            faceNodes[5] = (this.cellList[cell].nodes[4], this.cellList[cell].nodes[5], this.cellList[cell].nodes[6],
                            0);
          }
          when TOPO_HEXA
          {
            faceNodes[1] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[4], this.cellList[cell].nodes[3],
                            this.cellList[cell].nodes[2]);
            faceNodes[2] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[2], this.cellList[cell].nodes[6],
                            this.cellList[cell].nodes[5]);
            faceNodes[3] = (this.cellList[cell].nodes[2], this.cellList[cell].nodes[3], this.cellList[cell].nodes[7],
                            this.cellList[cell].nodes[6]);
            faceNodes[4] = (this.cellList[cell].nodes[3], this.cellList[cell].nodes[4], this.cellList[cell].nodes[8],
                            this.cellList[cell].nodes[7]);
            faceNodes[5] = (this.cellList[cell].nodes[1], this.cellList[cell].nodes[5], this.cellList[cell].nodes[8],
                            this.cellList[cell].nodes[4]);
            faceNodes[6] = (this.cellList[cell].nodes[5], this.cellList[cell].nodes[6], this.cellList[cell].nodes[7],
                            this.cellList[cell].nodes[8]);
          }
          otherwise {}
        }

        for cellFaceIdx in this.cellList[cell].faces.domain
        {
          // Check if this face is already in the map
          if faceMap_d.contains(sort_tuple(faceNodes[cellFaceIdx]))
          {
            // This a mapped face

            // Save the face ID to the cell element
            this.cellList[cell].faces[cellFaceIdx] = faceMap[sort_tuple(faceNodes[cellFaceIdx])];
            // Save the side of the face this cell is on
            this.cellList[cell].sides[cellFaceIdx] = 1;
            // Save the cell ID as the left neighbor of the face
            this.faceList[faceMap[sort_tuple(faceNodes[cellFaceIdx])]].cells[1] = cell;
          }
          else
          {
            //If it's not then create the new face

            // Increment cell count with new element
            this.nFaces += 1;
            // Resize domain to expand the array
            this.faceList_d = {1..this.nFaces};

            // Save the face ID in the cell
            this.cellList[cell].faces[cellFaceIdx] = this.nFaces;
            // Save the side of the face this cell is on
            this.cellList[cell].sides[cellFaceIdx] = 2;

            // Fill up the face properties. Maybe this should be an initializer?
            select elem_topology(this.cellList[cell].elemType)
            {
              when TOPO_LINE
              {
                this.faceList[this.nFaces].elemType = TYPE_NODE;
                this.faceList[this.nFaces].nodes_d = {1..1};
                this.faceList[this.nFaces].nodes[1] = faceNodes[cellFaceIdx][0];
              }
              when TOPO_TRIA
              {
                this.faceList[this.nFaces].elemType = TYPE_LINE_2;

                // Reverse the order of the face nodes since this node set is from the right side cell
                this.faceList[this.nFaces].nodes_d = {1..2};
                this.faceList[this.nFaces].nodes[1] = faceNodes[cellFaceIdx][1];
                this.faceList[this.nFaces].nodes[2] = faceNodes[cellFaceIdx][0];
              }
              when TOPO_QUAD
              {
                this.faceList[this.nFaces].elemType = TYPE_LINE_2;

                // Reverse the order of the face nodes since this node set is from the right side cell
                this.faceList[this.nFaces].nodes_d = {1..2};
                this.faceList[this.nFaces].nodes[1] = faceNodes[cellFaceIdx][1];
                this.faceList[this.nFaces].nodes[2] = faceNodes[cellFaceIdx][0];
              }
              when TOPO_TETR {}
              otherwise {}
            }
            this.faceList[this.nFaces].cells[2] = cell;

            // Save nodes that define it in the map
            faceMap_d.add(sort_tuple(faceNodes[cellFaceIdx]));
            faceMap[sort_tuple(faceNodes[cellFaceIdx])] = this.nFaces;
          }
        }
      }
    }

    proc set_families(inputFamlList)
    {
      // Loop through the mesh families
      for meshFaml in this.famlList
      {
        // Check the input file for the corresponding family
        for inputFaml in inputFamlList
        {
          if meshFaml.name == inputFaml.name && meshFaml.nDim == inputFaml.nDim
          {
            meshFaml.name           = inputFaml.name;
            meshFaml.nDim           = inputFaml.nDim;
            meshFaml.bocoType       = inputFaml.bocoType;
            meshFaml.bocoSubType    = inputFaml.bocoSubType;
            meshFaml.bocoProperties = inputFaml.bocoProperties;
            continue;
          }
        }
      }
    }

    proc elem_set_builder()
    {
      for cell in this.cellList
      {
        cellTypes.add(cell.elemType);
        cellTopos.add(cell.elemTopo());
      }

      for face in this.faceList
      {
        faceTypes.add(face.elemType);
        faceTopos.add(face.elemTopo());
      }
    }

    proc cell_count()
    {
      use Parameters.ParamMesh;
      var cell_count_d : domain(rank=2, idxType=int) = {1..1,1..2};
      var cell_count : [cell_count_d] int = 0;
      var topo : int;

      // Initialize list with the first element
      cell_count[1,1] = elem_topology(this.cellList[1].elemType);

      // Go over cell list and count the number of elements of each topology
      for cell in this.cellList
      {
        topo = elem_topology(cell.elemType);

        // Search the count array for this topology
        for i in cell_count_d.dim(0)
        {
          if cell_count[i,1] == topo then
            // If there is a line for this topology then just add 1 do the count
            cell_count[i,2] += 1;
          else if i == cell_count_d.dim(0).high
          {
            // If not and this is the end of the array then add a new line and start counting
            cell_count_d = {1..cell_count_d.dim(0).high+1, cell_count_d.dim(1)};
            cell_count[i+1,1] = topo;
            cell_count[i+1,2] = 1;
          }
        }
      }

      return cell_count;
    }

    proc face_count()
    {
      use Parameters.ParamMesh;
      var face_count_d : domain(rank=2, idxType=int) = {1..1,1..2};
      var face_count : [face_count_d] int = 0;
      var topo : int;

      // Initialize list with the first element
      face_count[1,1] = elem_topology(this.faceList[1].elemType);

      // Go over face list and count the number of elements of each topology
      for face in this.faceList
      {
        topo = elem_topology(face.elemType);

        // Search the count array for this topology
        for i in face_count_d.dim(0)
        {
          if face_count[i,1] == topo then
            // If there is a line for this topology then just add 1 do the count
            face_count[i,2] += 1;
          else if i == face_count_d.dim(0).high
          {
            // If not and this is the end of the array then add a new line and start counting
            face_count_d = {1..face_count_d.dim(0).high+1, face_count_d.dim(1)};
            face_count[i+1,1] = topo;
            face_count[i+1,2] = 1;
          }
        }
      }

      return face_count;
    }
  }

  proc elem_type_gmsh2mesh(in elemTypeGmsh : int) : int
  {
    use Parameters.ParamGmesh;
    use Parameters.ParamMesh;

    select elemTypeGmsh {
      when GMESH_PNT_1   do return TYPE_NODE;
      when GMESH_LIN_2   do return TYPE_LINE_2;
      when GMESH_LIN_3   do return TYPE_LINE_2;
      when GMESH_LIN_4   do return TYPE_LINE_2;
      when GMESH_LIN_5   do return TYPE_LINE_2;
      when GMESH_TRI_3   do return TYPE_TRIA_3;
      when GMESH_TRI_6   do return TYPE_TRIA_3;
      when GMESH_TRI_10  do return TYPE_TRIA_3;
      when GMESH_TRI_15  do return TYPE_TRIA_3;
      when GMESH_QUA_4   do return TYPE_QUAD_4;
      when GMESH_QUA_9   do return TYPE_QUAD_4;
      when GMESH_QUA_16  do return TYPE_QUAD_4;
      when GMESH_QUA_25  do return TYPE_QUAD_4;
      when GMESH_TET_4   do return TYPE_TETR_4;
      when GMESH_TET_10  do return TYPE_TETR_4;
      when GMESH_TET_20  do return TYPE_TETR_4;
      when GMESH_TET_35  do return TYPE_TETR_4;
      when GMESH_PYR_5   do return TYPE_PYRA_5;
      when GMESH_PYR_14  do return TYPE_PYRA_5;
      when GMESH_PYR_30  do return TYPE_PYRA_5;
      when GMESH_PYR_55  do return TYPE_PYRA_5;
      when GMESH_PRI_6   do return TYPE_PRIS_6;
      when GMESH_PRI_18  do return TYPE_PRIS_6;
      when GMESH_PRI_40  do return TYPE_PRIS_6;
      when GMESH_PRI_75  do return TYPE_PRIS_6;
      when GMESH_HEX_8   do return TYPE_HEXA_8;
      when GMESH_HEX_27  do return TYPE_HEXA_8;
      when GMESH_HEX_64  do return TYPE_HEXA_8;
      when GMESH_HEX_125 do return TYPE_HEXA_8;
      otherwise
      {
        writeln("Error converting element types from Gmesh to Native");
        writeln("   Gmesh Type: ", elemTypeGmsh);
        return -1;
      }
    }
  }

  proc elem_dimension(in elemTopo : int) : int
  {
    use Parameters.ParamMesh;

    select elemTopo {
      when TOPO_NODE do return 0;
      when TOPO_LINE do return 1;
      when TOPO_TRIA do return 2;
      when TOPO_QUAD do return 2;
      when TOPO_TETR do return 3;
      when TOPO_PYRA do return 3;
      when TOPO_PRIS do return 3;
      when TOPO_HEXA do return 3;
      otherwise return -1;
    }
  }

  proc elem_dimension_type(in elemType : int) : int
  {
    use Parameters.ParamMesh;

    select elemType {
      when TYPE_NODE     do return 0;
      when TYPE_LINE_2   do return 1;
      when TYPE_LINE_3   do return 1;
      when TYPE_LINE_4   do return 1;
      when TYPE_LINE_5   do return 1;
      when TYPE_TRIA_3   do return 2;
      when TYPE_TRIA_6   do return 2;
      when TYPE_TRIA_10  do return 2;
      when TYPE_TRIA_15  do return 2;
      when TYPE_QUAD_4   do return 2;
      when TYPE_QUAD_9   do return 2;
      when TYPE_QUAD_16  do return 2;
      when TYPE_QUAD_25  do return 2;
      when TYPE_TETR_4   do return 3;
      when TYPE_TETR_10  do return 3;
      when TYPE_TETR_20  do return 3;
      when TYPE_TETR_35  do return 3;
      when TYPE_PYRA_5   do return 3;
      when TYPE_PYRA_14  do return 3;
      when TYPE_PYRA_30  do return 3;
      when TYPE_PYRA_55  do return 3;
      when TYPE_PRIS_6   do return 3;
      when TYPE_PRIS_18  do return 3;
      when TYPE_PRIS_40  do return 3;
      when TYPE_PRIS_75  do return 3;
      when TYPE_HEXA_8   do return 3;
      when TYPE_HEXA_27  do return 3;
      when TYPE_HEXA_64  do return 3;
      when TYPE_HEXA_125 do return 3;
      otherwise return -1;
    }
  }

  proc elem_topology(in elemType : int) : int
  {
    use Parameters.ParamMesh;

    select elemType {
      when TYPE_NODE     do return TOPO_NODE;
      when TYPE_LINE_2   do return TOPO_LINE;
      when TYPE_LINE_3   do return TOPO_LINE;
      when TYPE_LINE_4   do return TOPO_LINE;
      when TYPE_LINE_5   do return TOPO_LINE;
      when TYPE_TRIA_3   do return TOPO_TRIA;
      when TYPE_TRIA_6   do return TOPO_TRIA;
      when TYPE_TRIA_10  do return TOPO_TRIA;
      when TYPE_TRIA_15  do return TOPO_TRIA;
      when TYPE_QUAD_4   do return TOPO_QUAD;
      when TYPE_QUAD_9   do return TOPO_QUAD;
      when TYPE_QUAD_16  do return TOPO_QUAD;
      when TYPE_QUAD_25  do return TOPO_QUAD;
      when TYPE_TETR_4   do return TOPO_TETR;
      when TYPE_TETR_10  do return TOPO_TETR;
      when TYPE_TETR_20  do return TOPO_TETR;
      when TYPE_TETR_35  do return TOPO_TETR;
      when TYPE_PYRA_5   do return TOPO_PYRA;
      when TYPE_PYRA_14  do return TOPO_PYRA;
      when TYPE_PYRA_30  do return TOPO_PYRA;
      when TYPE_PYRA_55  do return TOPO_PYRA;
      when TYPE_PRIS_6   do return TOPO_PRIS;
      when TYPE_PRIS_18  do return TOPO_PRIS;
      when TYPE_PRIS_40  do return TOPO_PRIS;
      when TYPE_PRIS_75  do return TOPO_PRIS;
      when TYPE_HEXA_8   do return TOPO_HEXA;
      when TYPE_HEXA_27  do return TOPO_HEXA;
      when TYPE_HEXA_64  do return TOPO_HEXA;
      when TYPE_HEXA_125 do return TOPO_HEXA;
      otherwise return -1;
    }
  }

  proc elem_order(in elemType : int) : int
  {
    use Parameters.ParamMesh;

    select elemType {
      when TYPE_NODE     do return  1; // Vertex
      when TYPE_LINE_2   do return  1; // 1st order Edge
      when TYPE_LINE_3   do return  2; // 2nd order Edge
      when TYPE_LINE_4   do return  3; // 3rd order Edge
      when TYPE_LINE_5   do return  4; // 4th order Edge
      when TYPE_TRIA_3   do return  1; // 1st order Triangle
      when TYPE_TRIA_6   do return  2; // 2nd order Triangle
      when TYPE_TRIA_10  do return  3; // 3rd order Triangle
      when TYPE_TRIA_15  do return  4; // 4th order Triangle
      when TYPE_QUAD_4   do return  1; // 1st order Quadrilateral
      when TYPE_QUAD_9   do return  2; // 2nd order Quadrilateral
      when TYPE_QUAD_16  do return  3; // 3rd order Quadrilateral
      when TYPE_QUAD_25  do return  4; // 4th order Quadrilateral
      when TYPE_TETR_4   do return  1; // 1st order Tetrahedron
      when TYPE_TETR_10  do return  2; // 2nd order Tetrahedron
      when TYPE_TETR_20  do return  3; // 3rd order Tetrahedron
      when TYPE_TETR_35  do return  4; // 4th order Tetrahedron
      when TYPE_PYRA_5   do return  1; // 1st order Pyramid
      when TYPE_PYRA_14  do return  2; // 2nd order Pyramid
      when TYPE_PYRA_30  do return  3; // 3rd order Pyramid
      when TYPE_PYRA_55  do return  4; // 4th order Pyramid
      when TYPE_PRIS_6   do return  1; // 1st order Prism
      when TYPE_PRIS_18  do return  2; // 2nd order Prism
      when TYPE_PRIS_40  do return  3; // 3rd order Prism
      when TYPE_PRIS_75  do return  4; // 4th order Prism
      when TYPE_HEXA_8   do return  1; // 1st order Hexahedron
      when TYPE_HEXA_27  do return  2; // 2nd order Hexahedron
      when TYPE_HEXA_64  do return  3; // 3rd order Hexahedron
      when TYPE_HEXA_125 do return  4; // 4th order Hexahedron
      otherwise return -1;
    }
  }

  proc elem_vertices(in elemTopo : int) : int
  {
    use Parameters.ParamMesh;

    select elemTopo {
      when TOPO_NODE do return 1; // Vertex
      when TOPO_LINE do return 2; // Edge
      when TOPO_TRIA do return 3; // Triangle
      when TOPO_QUAD do return 4; // Quadrilateral
      when TOPO_TETR do return 4; // Tetrahedron
      when TOPO_PYRA do return 5; // Pyramid
      when TOPO_PRIS do return 6; // Prism
      when TOPO_HEXA do return 8; // Hexahedron
      otherwise return -1;
    }
  }

  proc elem_nodes(in elemType : int) : int
  {
    use Parameters.ParamMesh;

    select elemType {
      when TYPE_NODE     do return   1; // Vertex
      when TYPE_LINE_2   do return   2; // 1st order Edge
      when TYPE_LINE_3   do return   3; // 2nd order Edge
      when TYPE_LINE_4   do return   4; // 3rd order Edge
      when TYPE_LINE_5   do return   5; // 4th order Edge
      when TYPE_TRIA_3   do return   3; // 1st order Triangle
      when TYPE_TRIA_6   do return   6; // 2nd order Triangle
      when TYPE_TRIA_10  do return  10; // 3rd order Triangle
      when TYPE_TRIA_15  do return  15; // 4th order Triangle
      when TYPE_QUAD_4   do return   4; // 1st order Quadrilateral
      when TYPE_QUAD_9   do return   9; // 2nd order Quadrilateral
      when TYPE_QUAD_16  do return  16; // 3rd order Quadrilateral
      when TYPE_QUAD_25  do return  25; // 4th order Quadrilateral
      when TYPE_TETR_4   do return   4; // 1st order Tetrahedron
      when TYPE_TETR_10  do return  10; // 2nd order Tetrahedron
      when TYPE_TETR_20  do return  20; // 3rd order Tetrahedron
      when TYPE_TETR_35  do return  35; // 4th order Tetrahedron
      when TYPE_PYRA_5   do return   5; // 1st order Pyramid
      when TYPE_PYRA_14  do return  14; // 2nd order Pyramid
      when TYPE_PYRA_30  do return  30; // 3rd order Pyramid
      when TYPE_PYRA_55  do return  55; // 4th order Pyramid
      when TYPE_PRIS_6   do return   6; // 1st order Prism
      when TYPE_PRIS_18  do return  18; // 2nd order Prism
      when TYPE_PRIS_40  do return  40; // 3rd order Prism
      when TYPE_PRIS_75  do return  75; // 4th order Prism
      when TYPE_HEXA_8   do return   8; // 1st order Hexahedron
      when TYPE_HEXA_27  do return  27; // 2nd order Hexahedron
      when TYPE_HEXA_64  do return  64; // 3rd order Hexahedron
      when TYPE_HEXA_125 do return 125; // 4th order Hexahedron
      otherwise return -1;
    }
  }

  proc elem_edges(in elemTopo : int) : int
  {
    use Parameters.ParamMesh;

    select elemTopo {
      when TOPO_NODE do return  0; // Vertex
      when TOPO_LINE do return  1; // Edge
      when TOPO_TRIA do return  3; // Triangle
      when TOPO_QUAD do return  4; // Quadrilateral
      when TOPO_TETR do return  6; // Tetrahedron
      when TOPO_PYRA do return  8; // Pyramid
      when TOPO_PRIS do return  9; // Prism
      when TOPO_HEXA do return 12; // Hexahedron
      otherwise return -1;
    }
  }

  proc elem_faces(in elemTopo : int) : int
  {
    use Parameters.ParamMesh;

    select elemTopo {
      when TOPO_NODE do return 0; // Vertex
      when TOPO_LINE do return 2; // Edge
      when TOPO_TRIA do return 3; // Triangle
      when TOPO_QUAD do return 4; // Quadrilateral
      when TOPO_TETR do return 4; // Tetrahedron
      when TOPO_PYRA do return 5; // Pyramid
      when TOPO_PRIS do return 5; // Prism
      when TOPO_HEXA do return 6; // Hexahedron
      otherwise return -1;
    }
  }

  proc sort_tuple(in tuple : 4*int) : 4*int
  {
    if tuple[0] < tuple[2] then tuple[0] <=> tuple[2];
    if tuple[1] < tuple[3] then tuple[1] <=> tuple[3];
    if tuple[0] < tuple[1] then tuple[0] <=> tuple[1];
    if tuple[2] < tuple[3] then tuple[2] <=> tuple[3];
    if tuple[1] < tuple[2] then tuple[1] <=> tuple[2];
    return tuple;
  }

  proc main()
  {
    use Gmesh;

    var test_gmesh2 = new unmanaged gmesh2_c();
    test_gmesh2.random1D(nCells=6, xMin=-1.0, xMax=1.0);

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
  }
}
