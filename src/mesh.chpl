module Mesh
{
  use Gmesh;
  use Random;
  use UnitTest;
  use Set;
  use BlockDist;
  //use Block1DDist;
  use ReplicatedDist;

  //////////////////
  //  Mesh Class  //
  //////////////////

  class mesh_c
  {
    // Mesh basic parameters
    const nDims : int;

    // Mesh mutable parameters
    const nNodes : int = 0;
    const nCells : int = 0;
    const nFaces : int = 0;
    const nBocos : int = 0;
    const nEdges : int = 0;
    const nFamls : int = 0;

    // Array sizing domains
    const nodeSpace : domain(rank=1, idxType=int) = {1..0};
    const cellSpace : domain(rank=1, idxType=int) = {1..0};
    const faceSpace : domain(rank=1, idxType=int) = {1..0};
    const bocoSpace : domain(rank=1, idxType=int) = {1..0};
    const edgeSpace : domain(rank=1, idxType=int) = {1..0};
    const famlSpace : domain(rank=1, idxType=int) = {1..0};

    // Single Locale Domains
    //var nodeList_d : domain(rank=1, idxType=int);// dmapped nodeList_dist;
    //var cellList_d : domain(rank=1, idxType=int);// dmapped cellList_dist;
    //var faceList_d : domain(rank=1, idxType=int);// dmapped faceList_dist;
    //var bocoList_d : domain(rank=1, idxType=int);// dmapped bocoList_dist;
    //var edgeList_d : domain(rank=1, idxType=int);// dmapped edgeList_dist;

    // Block Dist Domains
    var nodeList_d = nodeSpace dmapped Block(boundingBox=nodeSpace);
    var cellList_d = cellSpace dmapped Block(boundingBox=cellSpace);
    var faceList_d = faceSpace dmapped Block(boundingBox=faceSpace);
    var bocoList_d = bocoSpace dmapped Block(boundingBox=bocoSpace);
    var edgeList_d = edgeSpace dmapped Block(boundingBox=edgeSpace);

    // Block1D Dist Domains
    //var nodeList_d = nodeSpace dmapped Block1D(boundingBox=nodeSpace);
    //var cellList_d = cellSpace dmapped Block1D(boundingBox=cellSpace);
    //var faceList_d = faceSpace dmapped Block1D(boundingBox=faceSpace);
    //var bocoList_d = bocoSpace dmapped Block1D(boundingBox=bocoSpace);
    //var edgeList_d = edgeSpace dmapped Block1D(boundingBox=edgeSpace);

    var famlList_d : domain(rank=1, idxType=int) dmapped Replicated();

    // Arrays
    var nodeList : [nodeList_d] node_r;
    var cellList : [cellList_d] cell_r;
    var faceList : [faceList_d] face_r;
    var bocoList : [bocoList_d] boco_r;
    var edgeList : [edgeList_d] edge_r;
    var famlList : [famlList_d] faml_r;

    // Lists of types of elements in this mesh
    var faceTopos : set(int); // Gmesh element topology/shape. Ex: point, line, triangle, hexahedron
    var faceTypes : set(int); // CGNS element type, element topology + high-order spec. Ex: tetrahedron with 4 or 10 nodes
    var cellTopos : set(int); // Gmesh element topology/shape. Ex: point, line, triangle, hexahedron
    var cellTypes : set(int); // CGNS element type, element topology + high-order spec. Ex: tetrahedron with 4 or 10 nodes

    // Variable time step variables
    var cellTimeStep : [cellList_d] real; // Current calculated time-step
    var cellCharLeng : [cellList_d] real; // Cell characteristic length
    var minTimeStep : real;
    var minCharLeng : real;

    proc init(nDims : int, nNodes : int, nEdges : int, nFaces : int, nCells : int, nBocos : int, nFamls : int)
    {
      this.nDims = nDims;

      this.nNodes = nNodes;
      this.nCells = nCells;
      this.nFaces = nFaces;
      this.nBocos = nBocos;
      this.nEdges = nEdges;
      this.nFamls = nFamls;

      this.nodeSpace = {1..#this.nNodes};
      this.cellSpace = {1..#this.nCells};
      this.faceSpace = {1..#this.nFaces};
      this.bocoSpace = {1..#this.nBocos};
      this.edgeSpace = {1..#this.nEdges};
      this.famlSpace = {1..#this.nFamls};

      this.nodeList_d = this.nodeSpace dmapped Block(boundingBox=this.nodeSpace);
      this.cellList_d = this.cellSpace dmapped Block(boundingBox=this.cellSpace);
      this.faceList_d = this.faceSpace dmapped Block(boundingBox=this.faceSpace);
      this.bocoList_d = this.bocoSpace dmapped Block(boundingBox=this.bocoSpace);
      //this.nodeList_d = this.nodeSpace dmapped Block1D(boundingBox=this.nodeSpace);
      //this.cellList_d = this.cellSpace dmapped Block1D(boundingBox=this.cellSpace);
      //this.faceList_d = this.faceSpace dmapped Block1D(boundingBox=this.faceSpace);
      //this.bocoList_d = this.bocoSpace dmapped Block1D(boundingBox=this.bocoSpace);

      this.edgeList_d = if (nEdges == 0) then {1..#this.nEdges} dmapped Block(boundingBox={1..1})
                                         else {1..#this.nEdges} dmapped Block(boundingBox={1..#this.nEdges});
      //this.edgeList_d = if (nEdges == 0) then {1..#this.nEdges} dmapped Block1D(boundingBox={1..1})
      //                                   else {1..#this.nEdges} dmapped Block1D(boundingBox={1..#this.nEdges});

      this.famlList_d = {1..#this.nFamls} dmapped Replicated();
    }

    proc init(mesh : gmesh2_c)
    {
      use Gmesh;

      // Initialize basic constants
      this.nDims = mesh.mesh_dimension();

      // Initialize mesh element counts
      this.nNodes = mesh.nodes.domain.dim(0).size;
      this.nCells = mesh.cell_cnt();
      this.nFaces = mesh.face_cnt();
      this.nBocos = mesh.boco_cnt();
      this.nEdges = 0;
      this.nFamls = mesh.faml_cnt();

      // Initialize non-distributed domains
      this.nodeSpace = {1..#this.nNodes};
      this.cellSpace = {1..#this.nCells};
      this.faceSpace = {1..#this.nFaces};
      this.bocoSpace = {1..#this.nBocos};
      this.edgeSpace = {1..#this.nEdges};
      this.famlSpace = {1..#this.nFamls};

      // Initialize distributed domains
      this.nodeList_d = this.nodeSpace dmapped Block(boundingBox=this.nodeSpace);
      this.cellList_d = this.cellSpace dmapped Block(boundingBox=this.cellSpace);
      this.faceList_d = this.faceSpace dmapped Block(boundingBox=this.faceSpace);
      this.bocoList_d = this.bocoSpace dmapped Block(boundingBox=this.bocoSpace);
      this.edgeList_d = if (nEdges == 0) then this.edgeSpace dmapped Block(boundingBox={1..1})
                                         else this.edgeSpace dmapped Block(boundingBox=this.edgeSpace);
      this.famlList_d = {1..#this.nFamls} dmapped Replicated();
    }

    proc import_gmesh2(gmesh : gmesh2_c)
    {
      use Gmesh;

      // Copy nodes
      this.nodeList_d = {1..this.nNodes};
      for node in this.nodeList_d do
        this.nodeList[node].xyz = gmesh.nodes[node,1..3];

      // Get family names and dimensions. These are used to sort the elements.
      this.famlList_d = {1..this.nFamls};
      forall loc in Locales do
        on loc do
          for faml in this.famlList_d
          {
            this.famlList[faml].name = gmesh.families[faml].name;
            this.famlList[faml].nDim = gmesh.families[faml].nDim;
          }

      var cellDim : int = max reduce this.famlList[..].nDim;
      var bocoDim : int = cellDim - 1;

      var cellIdx : int = 0;
      var bocoIdx : int = 0;

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
            cellIdx += 1;
            // Fill up the cell properties. Maybe this should be an initializer?
            this.cellList[cellIdx].nodes_d  = gmesh.elements[element].nodes_d;
            this.cellList[cellIdx].nodes    = gmesh.elements[element].nodes;
            this.cellList[cellIdx].elemType = elem_type_gmsh2mesh(gmesh.elements[element].elemType);
            this.cellList[cellIdx].family   = elemFamlIdx;
          }
          when bocoDim
          {
            // Increment boco count with new element
            bocoIdx += 1;
            // Fill up the boco properties. Maybe this should be an initializer?
            this.bocoList[bocoIdx].nodes_d  = gmesh.elements[element].nodes_d;
            this.bocoList[bocoIdx].nodes    = gmesh.elements[element].nodes;
            this.bocoList[bocoIdx].elemType = elem_type_gmsh2mesh(gmesh.elements[element].elemType);
            this.bocoList[bocoIdx].family   = elemFamlIdx;
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

    proc build_face_list()
    {
      use Parameters.ParamMesh;
      use Parameters.ParamGmesh;
      import SortTuple.sort_tuple;

      // Build face list
      var faceVerts : [1..6] 4*int;
      var faceMap_d : domain(4*int);
      var faceMap : [faceMap_d] int;

      var faceIdx : int = 0;

      // Add all boundaries to the face map
      for bocoIdx in this.bocoList.domain
      {
        // Each boundary defines a face. Build the face nodes tuple based on the element topology.
        select elem_topology(this.bocoList[bocoIdx].elemType)
        {
          when TOPO_NODE do faceVerts[1] = (this.bocoList[bocoIdx].nodes[1], 0, 0, 0);
          when TOPO_LINE do faceVerts[1] = (this.bocoList[bocoIdx].nodes[1], this.bocoList[bocoIdx].nodes[2], 0, 0);
          when TOPO_TRIA do faceVerts[1] = (this.bocoList[bocoIdx].nodes[1], this.bocoList[bocoIdx].nodes[2],
                                            this.bocoList[bocoIdx].nodes[3],                              0 );
          when TOPO_QUAD do faceVerts[1] = (this.bocoList[bocoIdx].nodes[1], this.bocoList[bocoIdx].nodes[2],
                                            this.bocoList[bocoIdx].nodes[3], this.bocoList[bocoIdx].nodes[4]);
          otherwise {}
        }

        // Increment the face index for the current new face
        faceIdx += 1;

        // Fill the face properties
        this.faceList[faceIdx].elemType = this.bocoList[bocoIdx].elemType;

        // Fill the right side neighbor ID, boundaries are always on the right side of a face.
        // Boundaries have negative indexes so they can be easily distinguished from mesh cells.
        this.faceList[faceIdx].cells[2] = -bocoIdx;
        this.bocoList[bocoIdx].face     = faceIdx;

        // Add the face vertices to the face matching list and store the face index given to this face.
        faceMap_d.add(sort_tuple(faceVerts[1]));
        faceMap[sort_tuple(faceVerts[1])] = faceIdx;
      }

      // Add faces from cells to the face map and perform the matching
      for cellIdx in this.cellList.domain
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
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[2],                              0 );
            faceVerts[2] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2],
                            this.cellList[cellIdx].nodes[4],                              0 );
            faceVerts[3] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[4],                              0 );
            faceVerts[4] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[1],
                            this.cellList[cellIdx].nodes[4],                              0 );
          }
          when TOPO_PYRA
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[4],
                            this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[2]);
            faceVerts[2] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2],
                            this.cellList[cellIdx].nodes[5],                              0 );
            faceVerts[3] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[5],                              0 );
            faceVerts[4] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[4],
                            this.cellList[cellIdx].nodes[5],                              0 );
            faceVerts[5] = (this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[1],
                            this.cellList[cellIdx].nodes[5],                              0 );
          }
          when TOPO_PRIS
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2],
                            this.cellList[cellIdx].nodes[5], this.cellList[cellIdx].nodes[4]);
            faceVerts[2] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[6], this.cellList[cellIdx].nodes[5]);
            faceVerts[3] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[1],
                            this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[6]);
            faceVerts[4] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[2],                              0 );
            faceVerts[5] = (this.cellList[cellIdx].nodes[4], this.cellList[cellIdx].nodes[5],
                            this.cellList[cellIdx].nodes[6],                              0 );
          }
          when TOPO_HEXA
          {
            faceVerts[1] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[4],
                            this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[2]);
            faceVerts[2] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[2],
                            this.cellList[cellIdx].nodes[6], this.cellList[cellIdx].nodes[5]);
            faceVerts[3] = (this.cellList[cellIdx].nodes[2], this.cellList[cellIdx].nodes[3],
                            this.cellList[cellIdx].nodes[7], this.cellList[cellIdx].nodes[6]);
            faceVerts[4] = (this.cellList[cellIdx].nodes[3], this.cellList[cellIdx].nodes[4],
                            this.cellList[cellIdx].nodes[8], this.cellList[cellIdx].nodes[7]);
            faceVerts[5] = (this.cellList[cellIdx].nodes[1], this.cellList[cellIdx].nodes[5],
                            this.cellList[cellIdx].nodes[8], this.cellList[cellIdx].nodes[4]);
            faceVerts[6] = (this.cellList[cellIdx].nodes[5], this.cellList[cellIdx].nodes[6],
                            this.cellList[cellIdx].nodes[7], this.cellList[cellIdx].nodes[8]);
          }
          otherwise {}
        }

        // Iterate through this cell's faces
        for cellFaceIdx in this.cellList[cellIdx].faces.domain
        {
          // Check if this face is already in the map
          if faceMap_d.contains(sort_tuple(faceVerts[cellFaceIdx]))
          { // This a mapped face! :D

            // Get the face ID
            const faceID : int = faceMap[sort_tuple(faceVerts[cellFaceIdx])];

            // Save the cell ID as the left neighbor of the face
            this.faceList[faceID].cells[1] = cellIdx;

            // Fill face nodes from the left neighbors so that face are consistently oriented
            this.faceList[faceID].nodes_d = {1..elem_nodes(this.faceList[faceID].elemType)};
            this.faceList[faceID].nodes   = this.cellList[cellIdx].nodes[elem_face_nodes(this.cellList[cellIdx].elemType,
                                                                                         cellFaceIdx                    )];

            // Save the face ID and the side of the face this cell is on to the cell record
            // The second element to map to a face is always defined to be on the left side of the face
            this.cellList[cellIdx].faces[cellFaceIdx] = faceID;
            this.cellList[cellIdx].sides[cellFaceIdx] = 1;
          }
          else
          { // This isn't a mapped face. But don't panic! D:

            // Increment the face index for the current new face
            faceIdx += 1;

            // Fill up the face properties. First cell to contain a face is put of the right side of the face.
            this.faceList[faceIdx].cells[2] = cellIdx;
            this.faceList[faceIdx].elemType = elem_face_type(this.cellList[cellIdx].elemType, cellFaceIdx);

            // Save the face index and the side of the face this cell is on to the cell record
            // Boundaries are assumed to be always on the right side of a face and they are added first so we follow the same
            // convention and put the first cell that contains a face in it'sright side.
            this.cellList[cellIdx].faces[cellFaceIdx] = faceIdx;
            this.cellList[cellIdx].sides[cellFaceIdx] = 2;

            // Add the tuple of nodes that define this face to the map and save this face's face index
            faceMap_d.add(sort_tuple(faceVerts[cellFaceIdx]));
            faceMap[sort_tuple(faceVerts[cellFaceIdx])] = faceIdx;
          }
        }
      }
      // Check if all faces have been correctly identified?
    }

    proc build_elem_sets()
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

    proc build_cell_char_leng()
    {
      forall cellIdx in this.cellList.domain do
        cellCharLeng[cellIdx] = elem_char_leng(this.cellList[cellIdx].elemTopo(), this.cellList[cellIdx].nodes);

      this.minCharLeng = min reduce (cellCharLeng);
    }

    proc set_families(inputFamlList)
    {
      // Since the families list is replicated on all locales we need to do this on all of them
      forall loc in Locales do
        on loc do
          // Loop through the mesh's families
          label MESH_LIST for meshFaml in this.famlList
          {
            // Look for this family in the input
            label INPUT_LIST for inputFaml in inputFamlList
            {
              if meshFaml.name == inputFaml.name && meshFaml.nDim == inputFaml.nDim
              {
                // When we find a match copy the configuration
                meshFaml.name           = inputFaml.name;
                meshFaml.nDim           = inputFaml.nDim;
                meshFaml.bocoType       = inputFaml.bocoType;
                meshFaml.bocoSubType    = inputFaml.bocoSubType;
                meshFaml.bocoProperties = inputFaml.bocoProperties;
                continue MESH_LIST;
              }
            }
            writeln("ERROR: Couldn't find family ", meshFaml.name, "in input config");
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

    proc elem_char_leng(elemTopo : int, elemNodes : [] int) : real
    {
      use Parameters.ParamMesh;

      select elemTopo
      {
        when TOPO_NODE do return 0.0;

        when TOPO_LINE do return line_leng( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz );

        when TOPO_TRIA do return tria_min_height( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                                  this.nodeList[elemNodes[3]].xyz );

        when TOPO_QUAD do return quad_min_height( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                                  this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz );

        when TOPO_TETR do return  tetr_min_height( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                                   this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz );

        when TOPO_PYRA do return  pyra_min_height( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                                   this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz,
                                                   this.nodeList[elemNodes[5]].xyz );

        when TOPO_PRIS do return  pris_min_height( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                                   this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz,
                                                   this.nodeList[elemNodes[5]].xyz, this.nodeList[elemNodes[6]].xyz );

        when TOPO_HEXA do return  hexa_min_height( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                                   this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz,
                                                   this.nodeList[elemNodes[5]].xyz, this.nodeList[elemNodes[6]].xyz,
                                                   this.nodeList[elemNodes[7]].xyz, this.nodeList[elemNodes[8]].xyz );
        otherwise {
          writeln("Error calculating mesh element characteristic length");
          writeln("Invalid element topology: ", elemTopo);
          return -1.0;
        }
      }
    }

    proc elem_size (elemTopo : int, elemNodes : [] int) : real
    {
      use Parameters.ParamMesh;

      select elemTopo
      {
        when TOPO_NODE do return 0.0;

        when TOPO_LINE do return line_leng( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz );

        when TOPO_TRIA do return tria_area( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                            this.nodeList[elemNodes[3]].xyz );

        when TOPO_QUAD do return quad_area( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                            this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz );

        when TOPO_TETR do return  tetr_vol( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                            this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz );

        when TOPO_PYRA do return  pyra_vol( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                            this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz,
                                            this.nodeList[elemNodes[5]].xyz );

        when TOPO_PRIS do return  pris_vol( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                            this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz,
                                            this.nodeList[elemNodes[5]].xyz, this.nodeList[elemNodes[6]].xyz );

        when TOPO_HEXA do return  hexa_vol( this.nodeList[elemNodes[1]].xyz, this.nodeList[elemNodes[2]].xyz,
                                            this.nodeList[elemNodes[3]].xyz, this.nodeList[elemNodes[4]].xyz,
                                            this.nodeList[elemNodes[5]].xyz, this.nodeList[elemNodes[6]].xyz,
                                            this.nodeList[elemNodes[7]].xyz, this.nodeList[elemNodes[8]].xyz );
        otherwise {
          writeln("Error calculating mesh element size");
          writeln("Invalid element topology: ", elemTopo);
          return -1.0;
        }
      }
    }
  }

  ///////////////////////
  //  Mesh Components  //
  ///////////////////////

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

  /////////////////////////////////
  //  Element Probing Functions  //
  /////////////////////////////////

  proc elem_type_gmsh2mesh(in elemTypeGmsh : int) : int
  {
    use Parameters.ParamGmesh;
    use Parameters.ParamMesh;

    select elemTypeGmsh {
      when GMESH_PNT_1   do return TYPE_NODE;
      when GMESH_LIN_2   do return TYPE_LINE_2;
      when GMESH_LIN_3   do return TYPE_LINE_3;
      when GMESH_LIN_4   do return TYPE_LINE_4;
      when GMESH_LIN_5   do return TYPE_LINE_5;
      when GMESH_TRI_3   do return TYPE_TRIA_3;
      when GMESH_TRI_6   do return TYPE_TRIA_6;
      when GMESH_TRI_10  do return TYPE_TRIA_10;
      when GMESH_TRI_15  do return TYPE_TRIA_15;
      when GMESH_QUA_4   do return TYPE_QUAD_4;
      when GMESH_QUA_9   do return TYPE_QUAD_9;
      when GMESH_QUA_16  do return TYPE_QUAD_16;
      when GMESH_QUA_25  do return TYPE_QUAD_25;
      when GMESH_TET_4   do return TYPE_TETR_4;
      when GMESH_TET_10  do return TYPE_TETR_10;
      when GMESH_TET_20  do return TYPE_TETR_20;
      when GMESH_TET_35  do return TYPE_TETR_35;
      when GMESH_PYR_5   do return TYPE_PYRA_5;
      when GMESH_PYR_14  do return TYPE_PYRA_14;
      when GMESH_PYR_30  do return TYPE_PYRA_30;
      when GMESH_PYR_55  do return TYPE_PYRA_55;
      when GMESH_PRI_6   do return TYPE_PRIS_6;
      when GMESH_PRI_18  do return TYPE_PRIS_18;
      when GMESH_PRI_40  do return TYPE_PRIS_40;
      when GMESH_PRI_75  do return TYPE_PRIS_75;
      when GMESH_HEX_8   do return TYPE_HEXA_8;
      when GMESH_HEX_27  do return TYPE_HEXA_27;
      when GMESH_HEX_64  do return TYPE_HEXA_64;
      when GMESH_HEX_125 do return TYPE_HEXA_125;
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

  proc elem_degree(in elemType : int) : int
  {
    use Parameters.ParamMesh;

    select elemType {
      when TYPE_NODE     do return  0; // Nodes have no degree
      when TYPE_LINE_2   do return  1; // 1st degree Edge
      when TYPE_LINE_3   do return  2; // 2nd degree Edge
      when TYPE_LINE_4   do return  3; // 3rd degree Edge
      when TYPE_LINE_5   do return  4; // 4th degree Edge
      when TYPE_TRIA_3   do return  1; // 1st degree Triangle
      when TYPE_TRIA_6   do return  2; // 2nd degree Triangle
      when TYPE_TRIA_10  do return  3; // 3rd degree Triangle
      when TYPE_TRIA_15  do return  4; // 4th degree Triangle
      when TYPE_QUAD_4   do return  1; // 1st degree Quadrilateral
      when TYPE_QUAD_9   do return  2; // 2nd degree Quadrilateral
      when TYPE_QUAD_16  do return  3; // 3rd degree Quadrilateral
      when TYPE_QUAD_25  do return  4; // 4th degree Quadrilateral
      when TYPE_TETR_4   do return  1; // 1st degree Tetrahedron
      when TYPE_TETR_10  do return  2; // 2nd degree Tetrahedron
      when TYPE_TETR_20  do return  3; // 3rd degree Tetrahedron
      when TYPE_TETR_35  do return  4; // 4th degree Tetrahedron
      when TYPE_PYRA_5   do return  1; // 1st degree Pyramid
      when TYPE_PYRA_14  do return  2; // 2nd degree Pyramid
      when TYPE_PYRA_30  do return  3; // 3rd degree Pyramid
      when TYPE_PYRA_55  do return  4; // 4th degree Pyramid
      when TYPE_PRIS_6   do return  1; // 1st degree Prism
      when TYPE_PRIS_18  do return  2; // 2nd degree Prism
      when TYPE_PRIS_40  do return  3; // 3rd degree Prism
      when TYPE_PRIS_75  do return  4; // 4th degree Prism
      when TYPE_HEXA_8   do return  1; // 1st degree Hexahedron
      when TYPE_HEXA_27  do return  2; // 2nd degree Hexahedron
      when TYPE_HEXA_64  do return  3; // 3rd degree Hexahedron
      when TYPE_HEXA_125 do return  4; // 4th degree Hexahedron
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
    // This will probably only be used in 3D meshes but is available for all element topologies for consistency
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
  { // Assuming this mesh element is a cell how many faces does it have
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

  proc elem_face_type(elemType : int, elemFaceIdx : int = 0) : int
  {
    use Parameters.ParamMesh;

    select elemType
    {
      when TYPE_LINE_2   do return TYPE_NODE;
      when TYPE_LINE_3   do return TYPE_NODE;
      when TYPE_LINE_4   do return TYPE_NODE;
      when TYPE_LINE_5   do return TYPE_NODE;
      when TYPE_TRIA_3   do return TYPE_LINE_2;
      when TYPE_TRIA_6   do return TYPE_LINE_3;
      when TYPE_TRIA_10  do return TYPE_LINE_4;
      when TYPE_TRIA_15  do return TYPE_LINE_5;
      when TYPE_QUAD_4   do return TYPE_LINE_2;
      when TYPE_QUAD_9   do return TYPE_LINE_3;
      when TYPE_QUAD_16  do return TYPE_LINE_4;
      when TYPE_QUAD_25  do return TYPE_LINE_5;
      when TYPE_TETR_4   do return TYPE_TRIA_3;
      when TYPE_TETR_10  do return TYPE_TRIA_6;
      when TYPE_TETR_20  do return TYPE_TRIA_10;
      when TYPE_TETR_35  do return TYPE_TRIA_15;
      when TYPE_PYRA_5
      {
        select elemFaceIdx
        {
          when 1    do return TYPE_QUAD_4;
          otherwise do return TYPE_TRIA_3;
        }
      }
      when TYPE_PYRA_14
      {
        select elemFaceIdx
        {
          when 1    do return TYPE_QUAD_9;
          otherwise do return TYPE_TRIA_6;
        }
      }
      when TYPE_PYRA_30
      {
        select elemFaceIdx
        {
          when 1    do return TYPE_QUAD_16;
          otherwise do return TYPE_TRIA_10;
        }
      }
      when TYPE_PYRA_55
      {
        select elemFaceIdx
        {
          when 1    do return TYPE_QUAD_25;
          otherwise do return TYPE_TRIA_15;
        }
      }
      when TYPE_PRIS_6
      {
        select elemFaceIdx
        {
          when 4    do return TYPE_TRIA_3;
          when 5    do return TYPE_TRIA_3;
          otherwise do return TYPE_QUAD_4;
        }
      }
      when TYPE_PRIS_18
      {
        select elemFaceIdx
        {
          when 4    do return TYPE_TRIA_6;
          when 5    do return TYPE_TRIA_6;
          otherwise do return TYPE_QUAD_9;
        }
      }
      when TYPE_PRIS_40
      {
        select elemFaceIdx
        {
          when 4    do return TYPE_TRIA_10;
          when 5    do return TYPE_TRIA_10;
          otherwise do return TYPE_QUAD_16;
        }
      }
      when TYPE_PRIS_75
      {
        select elemFaceIdx
        {
          when 4    do return TYPE_TRIA_15;
          when 5    do return TYPE_TRIA_15;
          otherwise do return TYPE_QUAD_25;
        }
      }
      when TYPE_HEXA_8   do return TYPE_QUAD_4;
      when TYPE_HEXA_27  do return TYPE_QUAD_9;
      when TYPE_HEXA_64  do return TYPE_QUAD_16;
      when TYPE_HEXA_125 do return TYPE_QUAD_25;
      otherwise return -1;
    }
  }

  proc elem_face_nodes(elemType : int, cellFaceIdx : int) : [] int
  {
    use Parameters.ParamMesh;

    var faceType    : int = elem_face_type(elemType, cellFaceIdx);
    var faceNodes   : [1..elem_nodes(faceType)] int;

    select elemType {
      when TYPE_LINE_2 do
        faceNodes = cellFaceIdx;
      when TYPE_TRIA_3 do
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%3+1];
      when TYPE_TRIA_6
      {
        // Corner nodes
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%3+1];
        // Mid-edge nodes
        for faceNodeIdx in 3..elem_nodes(faceType) do
          faceNodes[faceNodeIdx] = 3 + (cellFaceIdx-1)*(elem_degree(elemType)-1) + (faceNodeIdx-2);
      }
      when TYPE_TRIA_10
      {
        // Corner nodes
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%3+1];
        // Mid-edge nodes
        for faceNodeIdx in 3..elem_nodes(faceType) do
          faceNodes[faceNodeIdx] = 3 + (cellFaceIdx-1)*(elem_degree(elemType)-1) + (faceNodeIdx-2);
      }
      when TYPE_TRIA_15
      {
        // Corner nodes
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%3+1];
        // Mid-edge nodes
        for faceNodeIdx in 3..elem_nodes(faceType) do
          faceNodes[faceNodeIdx] = 3 + (cellFaceIdx-1)*(elem_degree(elemType)-1) + (faceNodeIdx-2);
      }
      when TYPE_QUAD_4 do
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%4+1];
      when TYPE_QUAD_9
      {
        // Corner nodes
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%4+1];
        // Mid-edge nodes
        for faceNodeIdx in 3..elem_nodes(faceType) do
          faceNodes[faceNodeIdx] = 4 + (cellFaceIdx-1)*(elem_degree(elemType)-1) + (faceNodeIdx-2);
      }
      when TYPE_QUAD_16
      {
        // Corner nodes
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%4+1];
        // Mid-edge nodes
        for faceNodeIdx in 3..elem_nodes(faceType) do
          faceNodes[faceNodeIdx] = 4 + (cellFaceIdx-1)*(elem_degree(elemType)-1) + (faceNodeIdx-2);
      }
      when TYPE_QUAD_25
      {
        // Corner nodes
        faceNodes[1..2] = [cellFaceIdx, cellFaceIdx%4+1];
        // Mid-edge nodes
        for faceNodeIdx in 3..elem_nodes(faceType) do
          faceNodes[faceNodeIdx] = 4 + (cellFaceIdx-1)*(elem_degree(elemType)-1) + (faceNodeIdx-2);
      }
      when TYPE_TETR_4
      {
        select faceType
        {
          when 1 do faceNodes = [1, 3, 2]; // Zeta min
          when 2 do faceNodes = [1, 2, 4]; // Eta min
          when 3 do faceNodes = [2, 3, 4]; // Xi + Eta + Zeta = Const
          when 4 do faceNodes = [3, 1, 4]; // Xi min
        }
      }
      when TYPE_TETR_10
      {
        select faceType
        {
          //                     Corners | Mid-Edge
          when 1 do faceNodes = [ 1, 3, 2, 7,  6,  5]; // Zeta min
          when 2 do faceNodes = [ 1, 2, 4, 5,  9,  8]; // Eta min
          when 3 do faceNodes = [ 2, 3, 4, 6, 10,  9]; // Xi + Eta + Zeta = Const
          when 4 do faceNodes = [ 3, 1, 4, 7,  8, 10]; // Xi min
        }
      }
      when TYPE_TETR_20
      {
        select faceType
        {
          //                     Corners | Mid-Edge              | Mid-Face
          when 1 do faceNodes = [ 1, 3, 2, 10,  9,  8,  7,  6,  5, 17]; // Zeta min
          when 2 do faceNodes = [ 1, 2, 4,  5,  6, 13, 14, 12, 11, 18]; // Eta min
          when 3 do faceNodes = [ 2, 3, 4,  7,  8, 15, 16, 14, 13, 19]; // Xi + Eta + Zeta = Const
          when 4 do faceNodes = [ 3, 1, 4,  9, 10, 11, 12, 16, 15, 20]; // Xi min
        }
      }
      when TYPE_TETR_35
      {
        select faceType
        {
          //                     Corners | Mid-Edge                          | Mid-Face
          when 1 do faceNodes = [ 1, 3, 2, 13, 12, 11, 10,  9,  8,  7,  6,  5, 23, 24, 25]; // Zeta min
          when 2 do faceNodes = [ 1, 2, 4,  5,  6,  7, 17, 18, 19, 16, 15, 14, 26, 27, 28]; // Eta min
          when 3 do faceNodes = [ 2, 3, 4,  8,  9, 10, 20, 21, 22, 19, 18, 17, 29, 30, 31]; // Xi + Eta + Zeta = Const
          when 4 do faceNodes = [ 3, 1, 4, 11, 12, 13, 14, 15, 16, 22, 21, 20, 32, 33, 34]; // Xi min
        }
      }
      when TYPE_PYRA_5
      {
        select faceType
        {
          when 1 do faceNodes = [1, 4, 3, 2]; // Zeta min
          when 2 do faceNodes = [   1, 2, 5]; // Eta min
          when 3 do faceNodes = [   2, 3, 5]; //
          when 4 do faceNodes = [   3, 4, 5]; //
          when 5 do faceNodes = [   4, 1, 5]; // Xi min
        }
      }
      when TYPE_PYRA_14
      {
        select faceType
        {
          //                       Corners | Mid-Edge      | Mid-Face
          when 1 do faceNodes = [1, 4, 3, 2,  9,  8,  7,  6, 14]; // Zeta min
          when 2 do faceNodes = [   1, 2, 5,      6, 11, 10    ]; // Eta min
          when 3 do faceNodes = [   2, 3, 5,      7, 12, 11    ]; //
          when 4 do faceNodes = [   3, 4, 5,      8, 13, 12    ]; //
          when 5 do faceNodes = [   4, 1, 5,      9, 10, 13    ]; // Xi min
        }
      }
      when TYPE_PYRA_30
      {
        select faceType
        {
          //                       Corners | Mid-Edge                      | Mid-Face
          when 1 do faceNodes = [1, 4, 3, 2, 12, 13, 11, 10,  9,  8,  7,  6, 22, 25, 24, 23]; // Zeta min
          when 2 do faceNodes = [   1, 2, 5,          6,  7, 16, 16, 15, 14,             26]; // Eta min
          when 3 do faceNodes = [   2, 3, 5,          8,  9, 18, 19, 17, 16,             27]; //
          when 4 do faceNodes = [   3, 4, 5,         10, 11, 20, 21, 19, 18,             28]; //
          when 5 do faceNodes = [   4, 1, 5,         12, 13, 14, 15, 21, 22,             29]; // Xi min
        }
      }
      when TYPE_PYRA_55
      {
        select faceType
        {
          //                       Corners | Mid-Edge                                      | Mid-Face
          when 1 do faceNodes = [1, 4, 3, 2, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6, 30, 31, 32, 33, 34, 35, 36, 37, 38]; // Zeta min
          when 2 do faceNodes = [   1, 2, 5,              6,  7,  8, 21, 22, 23, 20, 19, 18,                         39, 40, 41]; // Eta min
          when 3 do faceNodes = [   2, 3, 5,              9, 10, 11, 24, 25, 26, 23, 22, 21,                         42, 43, 44]; //
          when 4 do faceNodes = [   3, 4, 5,             12, 13, 14, 27, 28, 29, 26, 25, 24,                         45, 46, 47]; //
          when 5 do faceNodes = [   4, 1, 5,             15, 16, 17, 18, 19, 20, 29, 28, 27,                         48, 49, 50]; // Xi min
        }
      }
      when TYPE_PRIS_6
      {
        select faceType
        {
          when 1 do faceNodes = [1, 2, 5, 4]; // Eta min
          when 2 do faceNodes = [2, 3, 6, 5]; // Xi + Eta = Const
          when 3 do faceNodes = [3, 1, 4, 6]; // Xi min
          when 4 do faceNodes = [   1, 3, 2]; // Zeta min
          when 5 do faceNodes = [   4, 5, 6]; // Zeta max
        }
      }
      when TYPE_PRIS_18
      {
        select faceType
        {
          //                       Corners | Mid-Edge      | Mid-Face
          when 1 do faceNodes = [1, 2, 5, 4,  7, 11, 13, 10, 16]; // Eta min
          when 2 do faceNodes = [2, 3, 6, 5,  8, 12, 14, 11, 17]; // Xi + Eta = Const
          when 3 do faceNodes = [3, 1, 4, 6,  9, 10, 15, 12, 18]; // Xi min
          when 4 do faceNodes = [   1, 3, 2,      9,  8,  7    ]; // Zeta min
          when 5 do faceNodes = [   4, 5, 6,     13, 14, 15    ]; // Zeta max
        }
      }
      when TYPE_PRIS_40
      {
        select faceType
        {
          //                       Corners | Mid-Edge                      | Mid-Face
          when 1 do faceNodes = [1, 2, 5, 4,  7,  8, 15, 16, 20, 19, 14, 13, 26, 27, 28, 29]; // Eta min
          when 2 do faceNodes = [2, 3, 6, 5,  9, 10, 17, 18, 22, 21, 16, 15, 30, 31, 32, 33]; // Xi + Eta = Const
          when 3 do faceNodes = [3, 1, 4, 6, 11, 12, 13, 14, 24, 23, 18, 17, 34, 35, 36, 37]; // Xi min
          when 4 do faceNodes = [   1, 3, 2,         12, 11, 10,  9,  8,  7,             25]; // Zeta min
          when 5 do faceNodes = [   4, 5, 6,         19, 20, 21, 22, 23, 24,             38]; // Zeta max
        }
      }
      when TYPE_PRIS_75
      {
        select faceType
        {
          //                       Corners | Mid-Edge                                      | Mid-Face
          when 1 do faceNodes = [1, 2, 5, 4,  7,  8,  9, 19, 20, 21, 27, 26, 25, 18, 17, 16, 37, 38, 39, 40, 41, 42, 43, 44, 45]; // Eta min
          when 2 do faceNodes = [2, 3, 6, 5, 10, 11, 12, 22, 23, 24, 30, 29, 28, 21, 20, 19, 46, 47, 48, 49, 50, 51, 52, 53, 54]; // Xi + Eta = Const
          when 3 do faceNodes = [3, 1, 4, 6, 13, 14, 15, 16, 17, 18, 33, 32, 31, 24, 23, 22, 55, 56, 57, 58, 59, 60, 61, 62, 63]; // Xi min
          when 4 do faceNodes = [   1, 3, 2,             15, 14, 13, 12, 11, 10,  9,  8,  7,                         34, 35, 36]; // Zeta min
          when 5 do faceNodes = [   4, 5, 6,             25, 26, 27, 28, 29, 30, 31, 32, 33,                         64, 65, 66]; // Zeta max
        }
      }
      when TYPE_HEXA_8
      {
        select faceType
        {
          when 1 do faceNodes = [1, 4, 3, 2]; // Zeta min
          when 2 do faceNodes = [1, 2, 6, 5]; // Eta min
          when 3 do faceNodes = [2, 3, 7, 6]; // Xi min
          when 4 do faceNodes = [3, 4, 8, 7]; // Zeta max
          when 5 do faceNodes = [1, 5, 8, 4]; // Eta max
          when 6 do faceNodes = [5, 6, 7, 8]; // Xi max
        }
      }
      when TYPE_HEXA_27
      {
        select faceType
        {
          //                       Corners | Mid-Edge      | Mid-Face
          when 1 do faceNodes = [1, 4, 3, 2, 12, 11, 10,  9, 21]; // Zeta min
          when 2 do faceNodes = [1, 2, 6, 5,  9, 14, 17, 13, 22]; // Eta min
          when 3 do faceNodes = [2, 3, 7, 6, 10, 15, 18, 14, 23]; // Xi min
          when 4 do faceNodes = [3, 4, 8, 7, 11, 16, 19, 15, 24]; // Zeta max
          when 5 do faceNodes = [1, 5, 8, 4, 13, 20, 16, 12, 25]; // Eta max
          when 6 do faceNodes = [5, 6, 7, 8, 17, 18, 19, 20, 26]; // Xi max
        }
      }
      when TYPE_HEXA_64
      {
        select faceType
        {
          //                       Corners | Mid-Edge                      | Mid-Face
          when 1 do faceNodes = [1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10,  9, 33, 36, 35, 34]; // Zeta min
          when 2 do faceNodes = [1, 2, 6, 5,  9, 10, 19, 20, 26, 25, 18, 17, 37, 38, 39, 40]; // Eta min
          when 3 do faceNodes = [2, 3, 7, 6, 11, 12, 21, 22, 28, 27, 20, 19, 41, 42, 43, 44]; // Xi min
          when 4 do faceNodes = [3, 4, 8, 7, 13, 14, 23, 24, 30, 29, 22, 21, 45, 46, 47, 48]; // Zeta max
          when 5 do faceNodes = [1, 5, 8, 4, 17, 18, 32, 31, 24, 23, 15, 16, 49, 50, 51, 52]; // Eta max
          when 6 do faceNodes = [5, 6, 7, 8, 25, 26, 27, 28, 29, 30, 31, 32, 53, 54, 55, 56]; // Xi max
        }
      }
      when TYPE_HEXA_125
      {
        select faceType
        {
          //                       Corners | Mid-Edge                                      | Mid-Face
          when 1 do faceNodes = [1, 4, 3, 2, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9, 45, 46, 47, 48, 49, 50, 51, 52, 53]; // Zeta min
          when 2 do faceNodes = [1, 2, 6, 5,  9, 10, 11, 24, 25, 26, 35, 34, 33, 23, 22, 21, 54, 55, 56, 57, 58, 59, 60, 61, 62]; // Eta min
          when 3 do faceNodes = [2, 3, 7, 6, 12, 13, 14, 27, 28, 29, 38, 37, 36, 26, 25, 24, 63, 64, 65, 66, 67, 68, 69, 70, 71]; // Xi min
          when 4 do faceNodes = [3, 4, 8, 7, 15, 16, 17, 30, 31, 32, 41, 40, 39, 29, 28, 27, 72, 73, 74, 75, 76, 77, 78, 79, 80]; // Zeta max
          when 5 do faceNodes = [1, 5, 8, 4, 21, 22, 23, 44, 43, 42, 32, 31, 30, 18, 19, 20, 81, 82, 83, 84, 85, 86, 87, 88, 89]; // Eta max
          when 6 do faceNodes = [5, 6, 7, 8, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 90, 91, 92, 93, 94, 95, 96, 97, 98]; // Xi max
        }
      }
      otherwise do return [-1];
    }

    return faceNodes;
  }

  ///////////////////////////////////
  //  Mesh Element Size Functions  //
  ///////////////////////////////////

  proc elem_size(elemTopo : int, xyz : [] real) : real
  {
    use Parameters.ParamMesh;

    select elemTopo
    {
      when TOPO_NODE do return 0.0;
      when TOPO_LINE do return line_leng( xyz[1, 1..3], xyz[2, 1..3] );
      when TOPO_TRIA do return tria_area( xyz[1, 1..3], xyz[2, 1..3], xyz[3, 1..3] );
      when TOPO_QUAD do return quad_area( xyz[1, 1..3], xyz[2, 1..3], xyz[3, 1..3], xyz[4, 1..3] );
      when TOPO_TETR do return  tetr_vol( xyz[1, 1..3], xyz[2, 1..3], xyz[3, 1..3], xyz[4, 1..3] );
      when TOPO_PYRA do return  pyra_vol( xyz[1, 1..3], xyz[2, 1..3], xyz[3, 1..3], xyz[4, 1..3], xyz[5, 1..3] );
      when TOPO_PRIS do return  pris_vol( xyz[1, 1..3], xyz[2, 1..3], xyz[3, 1..3],
                                          xyz[4, 1..3], xyz[5, 1..3], xyz[6, 1..3] );
      when TOPO_HEXA do return  hexa_vol( xyz[1, 1..3], xyz[2, 1..3], xyz[3, 1..3], xyz[4, 1..3],
                                          xyz[5, 1..3], xyz[6, 1..3], xyz[7, 1..3], xyz[8, 1..3] );
      otherwise {
        writeln("Error calculating mesh element size");
        writeln("Invalid element topology: ", elemTopo);
        return -1.0;
      }
    }
  }

  proc line_leng( vert1 : [1..3] real, vert2 : [1..3] real ) : real
  {
    // Use the matrix formula to calculate the area of a triangle

    use LinearAlgebra;

    const edge12 : [1..3] real = vert2 - vert1;

    return norm(edge12, normType.norm2);
  }

  proc tria_area( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real ) : real
  {
    // Use the matrix formula to calculate the area of a triangle

    use LinearAlgebra;

    const edge12 : [1..3] real = vert2 - vert1;
    const edge13 : [1..3] real = vert3 - vert1;

    return 0.5*norm(cross(edge12, edge13), normType.norm2);
  }

  proc quad_area( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                  vert4 : [1..3] real                                           ) : real
  {
    // Split the quadrilateral into 2 triangles and use the matrix formula to calculate their areas

    use LinearAlgebra;

    const edge12 : [1..3] real = vert2 - vert1;
    const diag13 : [1..3] real = vert3 - vert1;
    const edge14 : [1..3] real = vert4 - vert1;

    const area1 : [1..3] real = cross(edge12, diag13);
    const area2 : [1..3] real = cross(diag13, edge14);

    return 0.5*(norm( area1, normType.norm2 ) + norm( area2, normType.norm2 ));
  }

  proc tetr_vol( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                 vert4 : [1..3] real                                           ) : real
  {
    import Determinant.determinant;

    var matrix : [1..3, 1..3] real;
    var volume : real = 0.0;

    matrix[1, 1..3] = vert2[1..3] - vert1[1..3];
    matrix[2, 1..3] = vert3[1..3] - vert1[1..3];
    matrix[3, 1..3] = vert4[1..3] - vert1[1..3];

    volume = determinant(matrix)/6.0;

    return volume;
  }

  proc pyra_vol( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                 vert4 : [1..3] real, vert5 : [1..3]                           ) : real
  {
    import Determinant.determinant;

    var matrix : [1..3, 1..3] real;
    var volume : real = 0.0;

    //writeln();
    //writeln("Calculating Pyramid Cell Volume");

    matrix[1, 1..3] = vert2[1..3] - vert1[1..3];
    matrix[2, 1..3] = vert3[1..3] - vert1[1..3];
    matrix[3, 1..3] = vert5[1..3] - vert1[1..3];
    //writeln("  Tetra1: ", determinant(matrix)/6.0);
    volume = determinant(matrix);

    matrix[1, 1..3] = vert3[1..3] - vert1[1..3];
    matrix[2, 1..3] = vert4[1..3] - vert1[1..3];
    matrix[3, 1..3] = vert5[1..3] - vert1[1..3];
    //writeln("  Tetra2: ", determinant(matrix)/6.0);
    volume += determinant(matrix);

    return volume/6.0;
  }

  proc pris_vol( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                 vert4 : [1..3] real, vert5 : [1..3] real, vert6 : [1..3] real ) : real
  {
    import Determinant.determinant;

    var matrix : [1..3, 1..3] real;
    var volume : real = 0.0;

    //writeln();
    //writeln("Calculating Prism Cell Volume");

    matrix[1, 1..3] = vert2[1..3] - vert1[1..3];
    matrix[2, 1..3] = vert3[1..3] - vert1[1..3];
    matrix[3, 1..3] = vert4[1..3] - vert1[1..3];
    //writeln("  Tetra1: ", determinant(matrix)/6.0);
    volume = determinant(matrix);

    matrix[1, 1..3] = vert4[1..3] - vert6[1..3];
    matrix[2, 1..3] = vert3[1..3] - vert6[1..3];
    matrix[3, 1..3] = vert5[1..3] - vert6[1..3];
    //writeln("  Tetra2: ", determinant(matrix)/6.0);
    volume += determinant(matrix);

    matrix[1, 1..3] = vert4[1..3] - vert5[1..3];
    matrix[2, 1..3] = vert3[1..3] - vert5[1..3];
    matrix[3, 1..3] = vert2[1..3] - vert5[1..3];
    //writeln("  Tetra3: ", determinant(matrix)/6.0);
    volume += determinant(matrix);

    return volume/6.0;
  }

  proc hexa_vol( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                 vert4 : [1..3] real, vert5 : [1..3] real, vert6 : [1..3] real,
                 vert7 : [1..3] real, vert8 : [1..3] real                      ) : real
  {
    import Determinant.determinant;

    var matrix : [1..3, 1..3] real;
    var volume : real = 0.0;

    //writeln();
    //writeln("Calculating Hexahedron Cell Volume");

    // First z=0 base tetra
    matrix[1, 1..3] = vert2[1..3] - vert1[1..3];
    matrix[2, 1..3] = vert4[1..3] - vert1[1..3];
    matrix[3, 1..3] = vert5[1..3] - vert1[1..3];
    //writeln("  Tetra1: ", determinant(matrix)/6.0);
    volume = determinant(matrix);

    // Second z=0 base tetra
    matrix[1, 1..3] = vert2[1..3] - vert3[1..3];
    matrix[2, 1..3] = vert7[1..3] - vert3[1..3];
    matrix[3, 1..3] = vert4[1..3] - vert3[1..3];
    //writeln("  Tetra2: ", determinant(matrix)/6.0);
    volume += determinant(matrix);

    // First z=1 base tetra
    matrix[1, 1..3] = vert2[1..3] - vert6[1..3];
    matrix[2, 1..3] = vert5[1..3] - vert6[1..3];
    matrix[3, 1..3] = vert7[1..3] - vert6[1..3];
    //writeln("  Tetra3: ", determinant(matrix)/6.0);
    volume += determinant(matrix);

    // Second z=1 base tetra
    matrix[1, 1..3] = vert5[1..3] - vert8[1..3];
    matrix[2, 1..3] = vert4[1..3] - vert8[1..3];
    matrix[3, 1..3] = vert7[1..3] - vert8[1..3];
    //writeln("  Tetra4: ", determinant(matrix)/6.0);
    volume += determinant(matrix);

    // Internal volume left
    matrix[1, 1..3] = vert2[1..3] - vert5[1..3];
    matrix[2, 1..3] = vert4[1..3] - vert5[1..3];
    matrix[3, 1..3] = vert7[1..3] - vert5[1..3];
    //writeln("  Tetra5: ", determinant(matrix)/6.0);
    volume += determinant(matrix);

    return volume/6.0;
  }

  ////////////////////////////////////////////
  //  Cell Characteristic Length Functions  //
  ////////////////////////////////////////////

  proc tria_min_height( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real ) : real
  {
    // Calculate the shortest height of a triangle from the coordinates of it's vertices

    use LinearAlgebra;

    // Get the double of the triangle area using the shoelace formula
    var area : real = tria_area(vert1, vert2, vert3);

    const edge12 : [1..3] real = vert2 - vert1;
    const edge23 : [1..3] real = vert3 - vert2;
    const edge31 : [1..3] real = vert1 - vert3;

    // Calculate the length of all sides and pick the largest
    var baseMax : real = max( norm(edge12, normType.norm2) ,
                              norm(edge23, normType.norm2) ,
                              norm(edge31, normType.norm2) );

    // The minimum height is the area divided by the largest base/side
    return 2.0*area/baseMax;
  }

  proc quad_min_height( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                        vert4 : [1..3] real                                           ) : real
  {
    // Calculate the minimum height of the quadrilateral

    use LinearAlgebra;

    // Calculate the length of each side
    const edge12 : [1..3] real = vert2-vert1;
    const edge23 : [1..3] real = vert3-vert2;
    const edge34 : [1..3] real = vert4-vert3;
    const edge41 : [1..3] real = vert1-vert4;

    var len12 : real = norm(edge12, normType.norm2);
    var len23 : real = norm(edge23, normType.norm2);
    var len34 : real = norm(edge34, normType.norm2);
    var len41 : real = norm(edge41, normType.norm2);

    // Get the area of each subtriangle
    var area234 : real = tria_area(vert2, vert3, vert4);
    var area341 : real = tria_area(vert3, vert4, vert1);
    var area412 : real = tria_area(vert4, vert1, vert2);
    var area123 : real = tria_area(vert1, vert2, vert3);

    // Calculate the distance of the 2 opposing vertices to each side and pick the maximum
    // Then select the shortest if these 4 heights
    var heightMin : real = 2.0*min( max(area412, area123)/len12 ,
                                    max(area234, area123)/len23 ,
                                    max(area234, area341)/len34 ,
                                    max(area341, area412)/len41 );

    return heightMin;
  }

  proc tetr_min_height( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                        vert4 : [1..3] real                                           ) : real
  {
    // Volume = (1/3) * BaseArea * Height
    // Height = 3 * Volume / BaseArea

    var base : real = max(tria_area(vert2, vert3, vert4),
                          tria_area(vert1, vert4, vert3),
                          tria_area(vert1, vert2, vert4),
                          tria_area(vert1, vert3, vert2) );

    return 3.0*tetr_vol(vert1, vert2, vert3, vert4)/base;
  }

  proc pyra_min_height( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                        vert4 : [1..3] real, vert5 : [1..3] real                      ) : real
  {
    return 3.0*min( // Face 1 heights
                    tetr_vol(vert1, vert2, vert4, vert5)/tria_area(vert1, vert4, vert2),

                    // Face 2 heights
                    max( tetr_vol(vert1, vert5, vert2, vert3),
                         tetr_vol(vert1, vert5, vert2, vert4) )/tria_area(vert1, vert2, vert5),

                    // Face 3 heights
                    max( tetr_vol(vert2, vert5, vert3, vert1),
                         tetr_vol(vert2, vert5, vert3, vert4) )/tria_area(vert2, vert3, vert5),

                    // Face 4 heights
                    max( tetr_vol(vert3, vert5, vert4, vert1),
                         tetr_vol(vert3, vert5, vert4, vert2) )/tria_area(vert3, vert4, vert5),

                    // Face 5 heights
                    max( tetr_vol(vert4, vert5, vert1, vert2),
                         tetr_vol(vert4, vert5, vert1, vert3) )/tria_area(vert4, vert1, vert5)
                  );
  }

  proc pris_min_height( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                        vert4 : [1..3] real, vert5 : [1..3] real, vert6 : [1..3] real ) : real
  {
    return 3.0*min( // Face 1 heights
                    max( tetr_vol(vert1, vert4, vert2, vert3),
                         tetr_vol(vert1, vert4, vert2, vert6) )/tria_area(vert1, vert2, vert4),

                    // Face 2 heights
                    max( tetr_vol(vert2, vert5, vert3, vert1),
                         tetr_vol(vert2, vert5, vert3, vert4) )/tria_area(vert2, vert3, vert5),

                    // Face 3 heights
                    max( tetr_vol(vert1, vert3, vert4, vert2),
                         tetr_vol(vert1, vert3, vert4, vert5) )/tria_area(vert1, vert4, vert3),

                    // Face 4 heights
                    max( tetr_vol(vert1, vert2, vert3, vert4),
                         tetr_vol(vert1, vert2, vert3, vert5),
                         tetr_vol(vert1, vert2, vert3, vert6) )/tria_area(vert1, vert3, vert2),

                    // Face 5 heights
                    max( tetr_vol(vert4, vert6, vert5, vert1),
                         tetr_vol(vert4, vert6, vert5, vert2),
                         tetr_vol(vert4, vert6, vert5, vert3) )/tria_area(vert4, vert5, vert6)
                  );
  }

  proc hexa_min_height( vert1 : [1..3] real, vert2 : [1..3] real, vert3 : [1..3] real,
                        vert4 : [1..3] real, vert5 : [1..3] real, vert6 : [1..3] real,
                        vert7 : [1..3] real, vert8 : [1..3] real                      ) : real
  {
    return 3.0*min( // Face 1 heights
                    max( tetr_vol(vert1, vert2, vert4, vert5),
                         tetr_vol(vert1, vert2, vert4, vert6),
                         tetr_vol(vert1, vert2, vert4, vert7),
                         tetr_vol(vert1, vert2, vert4, vert8) )/tria_area(vert1, vert2, vert4),

                    // Face 2 heights
                    max( tetr_vol(vert1, vert5, vert2, vert4),
                         tetr_vol(vert1, vert5, vert2, vert8),
                         tetr_vol(vert1, vert5, vert2, vert7),
                         tetr_vol(vert1, vert5, vert2, vert3) )/tria_area(vert1, vert2, vert5),

                    // Face 3 heights
                    max( tetr_vol(vert1, vert4, vert5, vert2),
                         tetr_vol(vert1, vert4, vert5, vert3),
                         tetr_vol(vert1, vert4, vert5, vert7),
                         tetr_vol(vert1, vert4, vert5, vert6) )/tria_area(vert1, vert4, vert5),

                    // Face 4 heights
                    max( tetr_vol(vert5, vert8, vert6, vert1),
                         tetr_vol(vert5, vert8, vert6, vert4),
                         tetr_vol(vert5, vert8, vert6, vert3),
                         tetr_vol(vert5, vert8, vert6, vert2) )/tria_area(vert5, vert6, vert8),

                    // Face 5 heights
                    max( tetr_vol(vert4, vert3, vert8, vert1),
                         tetr_vol(vert4, vert3, vert8, vert2),
                         tetr_vol(vert4, vert3, vert8, vert6),
                         tetr_vol(vert4, vert3, vert8, vert5) )/tria_area(vert4, vert3, vert8),

                    // Face 6 heights
                    max( tetr_vol(vert2, vert6, vert3, vert1),
                         tetr_vol(vert2, vert6, vert3, vert5),
                         tetr_vol(vert2, vert6, vert3, vert8),
                         tetr_vol(vert2, vert6, vert3, vert4) )/tria_area(vert2, vert3, vert6)
                  );
  }

  /////////////////////////
  //  Testing Procedure  //
  /////////////////////////

  proc main()
  {
    use Gmesh;
    use Testing;
    use Parameters.ParamMesh;

    writeln();
    writeln("--------------------------------------------------------------------------------");
    writeln();
    writeln("Test 1: Generate a random 1D mesh - Gmesh:");

    var test_gmesh2 = new unmanaged gmesh2_c();
    test_gmesh2.random1D(nCells=6, xMin=-1.0, xMax=1.0);
    writeln();
    writeln(test_gmesh2);

    writeln();
    writeln("--------------------------------------------------------------------------------");
    writeln();
    writeln("Test 2: Convert random mesh from Gmesh => Native:");

    // Get number of physical dimensions from mesh or input
    var test_mesh = new unmanaged mesh_c(nDims=1);
    test_mesh.import_gmesh2(test_gmesh2);
    test_mesh.build_cell_char_leng();
    writeln();
    writeln(test_mesh);

    writeln();
    writeln("--------------------------------------------------------------------------------");
    writeln();
    writeln("Test 3: Cell and Face counts by topology:");
    writeln();
    writeln("Cell counts:");
    writeln(test_mesh.cell_count());
    writeln("Face counts:");
    writeln(test_mesh.face_count());

    writeln();
    writeln("--------------------------------------------------------------------------------");
    writeln();
    writeln("Test 4: Calculate mesh element length/area/volume:");
    {
      // Vertices and reference values calculated with Mathematica script "MeshElementSize.m"
      // All vertices on quadrilateral faces must be coplanar
      const vert1 : [1..3] real = [  3.7974683544303797468354430379746835443037974683544e-2,
                                     1.1004126547455295735900962861072902338376891334250e-2,
                                     2.1691973969631236442516268980477223427331887201735e-3 ];
      const vert2 : [1..3] real = [  1.0472440944881889763779527559055118110236220472441e+0,
                                     1.5174506828528072837632776934749620637329286798179e-3,
                                    -1.9830028328611898016997167138810198300283286118980e-2 ];
      const vert3 : [1..3] real = [  1.0052413242807748377966353727258180499742236362837e+0,
                                     1.9781486945054348981743425265972475239538496550704e+0,
                                    -3.7933372255520781960666895242258112418885586133838e-2 ];
      const vert4 : [1..3] real = [  2.3738872403560830860534124629080118694362017804154e-2,
                                     1.9934138309549945115257958287596048298572996706915e+0,
                                    -1.6597510373443983402489626556016597510373443983402e-2 ];
      const vert5 : [1..3] real = [ -1.0118043844856661045531197301854974704890387858347e-2,
                                    -4.1916167664670658682634730538922155688622754491018e-2,
                                     3.0031201248049921996879875195007800312012480499220e+0 ];
      const vert6 : [1..3] real = [  9.7091855204967826520077260335877485377610528085053e-1,
                                    -5.1536506888283435473522093887582675766656712870591e-2,
                                     3.0041651072071832178859815959600809484731556420752e+0 ];
      const vert7 : [1..3] real = [  8.9474675213241119537897024534668256086146853085362e-1,
                                     1.8199116052025052422032508877600218964698897674390e+0,
                                     4.4081237748435055697809484989569817863820657771572e+0 ];
      const vert8 : [1..3] real = [ -4.6639390399799550175081311359315688855140302186572e-2,
                                     1.8351629033273894051855318360721711249111735321387e+0,
                                     4.4116391275403188403028577535296135156280577507533e+0 ];

      // Reference values
      const LineLeng : real = 1.0095537166582577686828933158331283089967683434117;
      const TriaArea : real = 1.0006111071724500240791894523022772118776921443961;
      const QuadArea : real = 1.9706018705641267330507965937128075578107481314163;
      const TetrVol  : real = 1.0001214782656461611348893435338889368568025784744;
      const PyraVol  : real = 1.9696375961994772192617584288865276601214484332378;
      const PrisVol  : real = 3.3755561111825761236808689900659946291972459913852;
      const HexaVol  : real = 7.0829400120134946303608705716380832379163203260756;

      writeln();
      writef("Using the topology specific funtions:\n");
      writef("Elem Topo  | Reference           | Calculated          | Abs Error | Rel Error\n");
      writef("  LineLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", LineLeng, line_leng(vert1, vert2),
                                                               error(LineLeng, line_leng(vert1, vert2)),
                                                      relative_error(LineLeng, line_leng(vert1, vert2)));
      writef("  TriaArea | %19.12er | %19.12er | %9.2er | %9.2er\n", TriaArea, tria_area(vert1, vert2, vert4),
                                                               error(TriaArea, tria_area(vert1, vert2, vert4)),
                                                      relative_error(TriaArea, tria_area(vert1, vert2, vert4)));
      writef("  QuadArea | %19.12er | %19.12er | %9.2er | %9.2er\n", QuadArea, quad_area(vert1, vert2, vert3, vert4),
                                                               error(QuadArea, quad_area(vert1, vert2, vert3, vert4)),
                                                      relative_error(QuadArea, quad_area(vert1, vert2, vert3, vert4)));
      writef("  TetrVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", TetrVol , tetr_vol( vert1, vert2, vert4, vert5),
                                                               error(TetrVol , tetr_vol( vert1, vert2, vert4, vert5)),
                                                      relative_error(TetrVol , tetr_vol( vert1, vert2, vert4, vert5)));
      writef("  PyraVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", PyraVol , pyra_vol( vert1, vert2, vert3, vert4, vert5),
                                                               error(PyraVol , pyra_vol( vert1, vert2, vert3, vert4, vert5)),
                                                      relative_error(PyraVol , pyra_vol( vert1, vert2, vert3, vert4, vert5)));
      writef("  PrisVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", PrisVol ,
                                 pris_vol( vert1, vert2, vert4, vert5, vert6, vert8),
                 error(PrisVol , pris_vol( vert1, vert2, vert4, vert5, vert6, vert8)),
        relative_error(PrisVol , pris_vol( vert1, vert2, vert4, vert5, vert6, vert8)));
      writef("  HexaVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", HexaVol ,
                                 hexa_vol( vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8),
                 error(HexaVol , hexa_vol( vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8)),
        relative_error(HexaVol , hexa_vol( vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8)));

      var xyzLine : [1..2, 1..3] real;
      xyzLine[1, 1..3] = vert1;
      xyzLine[2, 1..3] = vert2;

      var xyzTria : [1..3, 1..3] real;
      xyzTria[1, 1..3] = vert1;
      xyzTria[2, 1..3] = vert2;
      xyzTria[3, 1..3] = vert4;

      var xyzQuad : [1..4, 1..3] real;
      xyzQuad[1, 1..3] = vert1;
      xyzQuad[2, 1..3] = vert2;
      xyzQuad[3, 1..3] = vert3;
      xyzQuad[4, 1..3] = vert4;

      var xyzTetr : [1..4, 1..3] real;
      xyzTetr[1, 1..3] = vert1;
      xyzTetr[2, 1..3] = vert2;
      xyzTetr[3, 1..3] = vert4;
      xyzTetr[4, 1..3] = vert5;

      var xyzPyra : [1..6, 1..3] real;
      xyzPyra[1, 1..3] = vert1;
      xyzPyra[2, 1..3] = vert2;
      xyzPyra[3, 1..3] = vert3;
      xyzPyra[4, 1..3] = vert4;
      xyzPyra[5, 1..3] = vert5;

      var xyzPris : [1..6, 1..3] real;
      xyzPris[1, 1..3] = vert1;
      xyzPris[2, 1..3] = vert2;
      xyzPris[3, 1..3] = vert4;
      xyzPris[4, 1..3] = vert5;
      xyzPris[5, 1..3] = vert6;
      xyzPris[6, 1..3] = vert8;

      var xyzHexa : [1..8, 1..3] real;
      xyzHexa[1, 1..3] = vert1;
      xyzHexa[2, 1..3] = vert2;
      xyzHexa[3, 1..3] = vert3;
      xyzHexa[4, 1..3] = vert4;
      xyzHexa[5, 1..3] = vert5;
      xyzHexa[6, 1..3] = vert6;
      xyzHexa[7, 1..3] = vert7;
      xyzHexa[8, 1..3] = vert8;

      writeln();
      writef("Using the generic elem_size function:\n");
      writef("Elem Topo  | Reference           | Calculated          | Abs Error | Rel Error\n");
      writef("  LineLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", LineLeng, elem_size(TOPO_LINE, xyzLine),
                                                               error(LineLeng, elem_size(TOPO_LINE, xyzLine)),
                                                      relative_error(LineLeng, elem_size(TOPO_LINE, xyzLine)));
      writef("  TriaArea | %19.12er | %19.12er | %9.2er | %9.2er\n", TriaArea, elem_size(TOPO_TRIA, xyzTria),
                                                               error(TriaArea, elem_size(TOPO_TRIA, xyzTria)),
                                                      relative_error(TriaArea, elem_size(TOPO_TRIA, xyzTria)));
      writef("  QuadArea | %19.12er | %19.12er | %9.2er | %9.2er\n", QuadArea, elem_size(TOPO_QUAD, xyzQuad),
                                                               error(QuadArea, elem_size(TOPO_QUAD, xyzQuad)),
                                                      relative_error(QuadArea, elem_size(TOPO_QUAD, xyzQuad)));
      writef("  TetrVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", TetrVol , elem_size(TOPO_TETR, xyzTetr),
                                                               error(TetrVol , elem_size(TOPO_TETR, xyzTetr)),
                                                      relative_error(TetrVol , elem_size(TOPO_TETR, xyzTetr)));
      writef("  PyraVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", PyraVol , elem_size(TOPO_PYRA, xyzPyra),
                                                               error(PyraVol , elem_size(TOPO_PYRA, xyzPyra)),
                                                      relative_error(PyraVol , elem_size(TOPO_PYRA, xyzPyra)));
      writef("  PrisVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", PrisVol , elem_size(TOPO_PRIS, xyzPris),
                                                               error(PrisVol , elem_size(TOPO_PRIS, xyzPris)),
                                                      relative_error(PrisVol , elem_size(TOPO_PRIS, xyzPris)));
      writef("  HexaVol  | %19.12er | %19.12er | %9.2er | %9.2er\n", HexaVol , elem_size(TOPO_HEXA, xyzHexa),
                                                               error(HexaVol , elem_size(TOPO_HEXA, xyzHexa)),
                                                      relative_error(HexaVol , elem_size(TOPO_HEXA, xyzHexa)));

      writeln();
      writef("Using the mesh class methods:\n");
      writef("Test mesh cell | elem_size()\n");
      for cellIdx in test_mesh.cellList.domain do
        writef("  Cell %2i size | %19.12er\n", cellIdx,
          test_mesh.elem_size(test_mesh.cellList[cellIdx].elemTopo(), test_mesh.cellList[cellIdx].nodes)
        );
    }

    writeln();
    writeln("--------------------------------------------------------------------------------");
    writeln();
    writeln("Test 5: Calculate mesh element shortest passthrough distance:");
    {
      // Vertices and reference values calculated with Mathematica script "MeshElementCharLeng.m"
      // All vertices on quadrilateral faces must be coplanar
      const vert1 : [1..3] real = [  3.7974683544303797468354430379746835443037974683544e-2,
                                     1.1004126547455295735900962861072902338376891334250e-2,
                                     2.1691973969631236442516268980477223427331887201735e-3 ];
      const vert2 : [1..3] real = [  1.0472440944881889763779527559055118110236220472441e+0,
                                     1.5174506828528072837632776934749620637329286798179e-3,
                                    -1.9830028328611898016997167138810198300283286118980e-2 ];
      const vert3 : [1..3] real = [  1.0052413242807748377966353727258180499742236362837e+0,
                                     1.9781486945054348981743425265972475239538496550704e+0,
                                    -3.7933372255520781960666895242258112418885586133838e-2 ];
      const vert4 : [1..3] real = [  2.3738872403560830860534124629080118694362017804154e-2,
                                     1.9934138309549945115257958287596048298572996706915e+0,
                                    -1.6597510373443983402489626556016597510373443983402e-2 ];
      const vert5 : [1..3] real = [ -1.0118043844856661045531197301854974704890387858347e-2,
                                    -4.1916167664670658682634730538922155688622754491018e-2,
                                     3.0031201248049921996879875195007800312012480499220e+0 ];
      const vert6 : [1..3] real = [  9.7091855204967826520077260335877485377610528085053e-1,
                                    -5.1536506888283435473522093887582675766656712870591e-2,
                                     3.0041651072071832178859815959600809484731556420752e+0 ];
      const vert7 : [1..3] real = [  8.9474675213241119537897024534668256086146853085362e-1,
                                     1.8199116052025052422032508877600218964698897674390e+0,
                                     4.4081237748435055697809484989569817863820657771572e+0 ];
      const vert8 : [1..3] real = [ -4.6639390399799550175081311359315688855140302186572e-2,
                                     1.8351629033273894051855318360721711249111735321387e+0,
                                     4.4116391275403188403028577535296135156280577507533e+0 ];

      // Reference values
      const LineCharLeng : real = 1.0095537166582577686828933158331283089967683434117e+0;
      const TriaCharLeng : real = 8.9361432579282310476118100755966934326463688804176e-1;
      const QuadCharLeng : real = 1.0090861045298942923515333317876425756813824966488e+0;
      const TetrCharLeng : real = 8.4495650133317678516378493575905710319629504752832e-1;
      const PyraCharLeng : real = 9.4482836877434877107756908549626460820538193643161e-1;
      const PrisCharLeng : real = 8.9229626157848543100704494613623558312367332736129e-1;
      const HexaCharLeng : real = 1.0079388647889108935072665768907527132658799181828e+0;

      writeln();
      writef("Using the topology specific funtions:\n");
      writef("Elem Topo      | Reference           | Calculated          | Abs Error | Rel Error\n");
      writef("  LineCharLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", LineCharLeng, line_leng(vert1, vert2),
                                                                   error(LineCharLeng, line_leng(vert1, vert2)),
                                                          relative_error(LineCharLeng, line_leng(vert1, vert2)));
      writef("  TriaCharLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", TriaCharLeng, tria_min_height(vert1, vert2, vert4),
                                                                   error(TriaCharLeng, tria_min_height(vert1, vert2, vert4)),
                                                          relative_error(TriaCharLeng, tria_min_height(vert1, vert2, vert4)));
      writef("  QuadCharLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", QuadCharLeng, quad_min_height(vert1, vert2, vert3, vert4),
                                                                   error(QuadCharLeng, quad_min_height(vert1, vert2, vert3, vert4)),
                                                          relative_error(QuadCharLeng, quad_min_height(vert1, vert2, vert3, vert4)));
      writef("  TetrCharLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", TetrCharLeng, tetr_min_height(vert1, vert2, vert4, vert5),
                                                                   error(TetrCharLeng, tetr_min_height(vert1, vert2, vert4, vert5)),
                                                          relative_error(TetrCharLeng, tetr_min_height(vert1, vert2, vert4, vert5)));
      writef("  PyraCharLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", PyraCharLeng, pyra_min_height(vert1, vert2, vert3, vert4, vert5),
                                                                   error(PyraCharLeng, pyra_min_height(vert1, vert2, vert3, vert4, vert5)),
                                                          relative_error(PyraCharLeng, pyra_min_height(vert1, vert2, vert3, vert4, vert5)));
      writef("  PrisCharLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", PrisCharLeng,
                                     pris_min_height(vert1, vert2, vert4, vert5, vert6, vert8),
                 error(PrisCharLeng, pris_min_height(vert1, vert2, vert4, vert5, vert6, vert8)),
        relative_error(PrisCharLeng, pris_min_height(vert1, vert2, vert4, vert5, vert6, vert8)));
      writef("  HexaCharLeng | %19.12er | %19.12er | %9.2er | %9.2er\n", HexaCharLeng,
                                     hexa_min_height(vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8),
                 error(HexaCharLeng, hexa_min_height(vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8)),
        relative_error(HexaCharLeng, hexa_min_height(vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8)));

      writeln();
      writef("Using the mesh class methods:\n");
      writef("Test mesh cell | elem_char_leng()\n");
      for cellIdx in test_mesh.cellList.domain do
        writef("  Cell %2i size | %19.12er\n", cellIdx,
          test_mesh.elem_char_leng(test_mesh.cellList[cellIdx].elemTopo(), test_mesh.cellList[cellIdx].nodes)
        );
    }
  }
}
