prototype module Mesh
{
  use Random;
  use UnitTest;

  // Internal mesh structure

  class node_c
  {
    // Geometric properties
    var xyz : real;
  }

  class edge_c
  {
    // Geometric properties
    var geom_order : int;
    var nodes : [2] int;
  }

  class face_c
  {
    // Geometric properties
    var geom_type : int;
    var geom_order : int;
    var nodes : [elem_nodes(geom_type)] int;
    var edges : [elem_edges(geom_type)] int;

    // Numerical properties
    var sol_order : int;
    var first_fp : int;
    var num_fp : int;
  }

  class cell_c
  {
    // Geometric properties
    var geom_type : int;
    var geom_order : int;
    var nodes : [elem_nodes(geom_type)] int;
    var edges : [elem_edges(geom_type)] int;
    var faces : [elem_faces(geom_type)] int;

    // Numerical properties
    var sol_order : int;
    var first_sp : int;
    var num_sp : int;
  }

  class mesh_1d_c
  {
    var nDims : int = 1;
    var nNodes : int;
    var nFaces : int;
    var nCells : int;

    var nodeList : [nNodes] node_c;
    var faceList : [nFaces] face_c;
    var cellList : [nCells] cell_c;
  }

  class mesh_2d_c
  {
    var nDims : int = 2;
    var nNodes : int;
    var nFaces : int;
    var nCells : int;

    var nodeList : [nNodes] node_c;
    var faceList : [nFaces] face_c;
    var cellList : [nCells] cell_c;
  }

  class mesh_3d_c
  {
    var nDims : int = 3;
    var nNodes : int;
    var nEdges : int;
    var nFaces : int;
    var nCells : int;

    var nodeList : [nNodes] node_c;
    var edgeList : [nEdges] edge_c;
    var faceList : [nFaces] face_c;
    var cellList : [nCells] cell_c;
  }

  private proc elem_nodes(in elem_type : int) : int
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

  private proc elem_edges(in elem_type : int) : int
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

  private proc elem_faces(in elem_type : int) : int
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
