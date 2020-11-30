prototype module Mesh
{
  use Random;
  use UnitTest;

  // Internal mesh structure

  class node_c
  {
    // Geometric properties
    var xyz : real(32);
  }

  class edge_c
  {
    // Geometric properties
    var geom_order : uint(8);
    var nodes : [2] uint(32);
  }

  class face_c
  {
    // Geometric properties
    var geom_type : uint(8);
    var geom_order : uint(8);
    var nodes : [elem_nodes(geom_type)] uint(32);
    var edges : [elem_edges(geom_type)] uint(32);

    // Numerical properties
    var sol_order : uint(8);
    var first_fp : uint(32);
    var num_fp : uint(16);
  }

  class cell_c
  {
    // Geometric properties
    var geom_type : uint(8);
    var geom_order : uint(8);
    var nodes : [elem_nodes(geom_type)] uint(32);
    var edges : [elem_edges(geom_type)] uint(32);
    var faces : [elem_faces(geom_type)] uint(32);

    // Numerical properties
    var sol_order : uint(8);
    var first_sp : uint(32);
    var num_sp : uint(16);
  }

  class mesh_1d_c
  {
    var nDims : uint(8) = 1;
    var nNodes : uint(32);
    var nFaces : uint(32);
    var nCells : uint(32);

    var nodeList : [nNodes] node_c;
    var faceList : [nFaces] face_c;
    var cellList : [nCells] cell_c;
  }

  class mesh_2d_c
  {
    var nDims : uint(8) = 2;
    var nNodes : uint(32);
    var nFaces : uint(32);
    var nCells : uint(32);

    var nodeList : [nNodes] node_c;
    var faceList : [nFaces] face_c;
    var cellList : [nCells] cell_c;
  }

  class mesh_3d_c
  {
    var nDims : uint(8) = 3;
    var nNodes : uint(32);
    var nEdges : uint(32);
    var nFaces : uint(32);
    var nCells : uint(32);

    var nodeList : [nNodes] node_c;
    var edgeList : [nEdges] edge_c;
    var faceList : [nFaces] face_c;
    var cellList : [nCells] cell_c;
  }

  private proc elem_nodes(in elem_type : uint[8]) : uint[8]
  {
    select elem_type {
      when 0 do return 1; // Vertex
      when 1 do return 2; // Edge
      when 2 do return 3; // Triangle
      when 3 do return 4; // Quadrilateral
      when 4 do return 4; // Tetrahedron
      when 5 do return 5; // Pyramid
      when 6 do return 6; // Prism
      when 7 do return 8; // Hexahedron
    }
  }

  private proc elem_edges(in elem_type : uint[8]) : uint[8]
  {
    select elem_type {
      when 0 do return  0; // Vertex
      when 1 do return  1; // Edge
      when 2 do return  3; // Triangle
      when 3 do return  4; // Quadrilateral
      when 4 do return  6; // Tetrahedron
      when 5 do return  8; // Pyramid
      when 6 do return  9; // Prism
      when 7 do return 12; // Hexahedron
    }
  }

  private proc elem_faces(in elem_type : uint[8]) : uint[8]
  {
    select elem_type {
      when 0 do return 0; // Vertex
      when 1 do return 0; // Edge
      when 2 do return 1; // Triangle
      when 3 do return 1; // Quadrilateral
      when 4 do return 4; // Tetrahedron
      when 5 do return 5; // Pyramid
      when 6 do return 5; // Prism
      when 7 do return 6; // Hexahedron
    }
  }

}
