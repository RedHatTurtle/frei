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

  class gmesh_c
  {
    var nNodes : int;
    var nElements : int;
    var nFamilies : int;

    var nodes    : [1..nNodes, 1..3] real;
    var elements : [1..nElements] gmesh_element_r;
    var families : [0..nFamilies-1] gmesh_family_r;

    proc random1D(xMin: real, xMax: real)
    {
      use Random;
      use Parameters.Tests;
      use Parameters.Gmesh;
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
      this.families[0].nDim = 1;
      this.families[0].name = "flow";
      this.families[1].nDim = 0;
      this.families[1].name = "left";
      this.families[2].nDim = 0;
      this.families[2].name = "right";

      // Fill node list in randomized order
      for i in 1..this.nNodes do
        this.nodes[nodePermutation[i],1] = x[i];

      // Fill element list with the internal elements in random order
      for i in 2..this.nElements-1 {
        this.elements[elemPermutation[i]].elemType = GMESH_LIN_2;
        this.elements[elemPermutation[i]].setNodes;
        this.elements[elemPermutation[i]].tags[1] = 0;
        this.elements[elemPermutation[i]].nodes[1] = nodePermutation[cells[i-1,1]];
        this.elements[elemPermutation[i]].nodes[2] = nodePermutation[cells[i-1,2]];
      }

      // Add left boundary point to elements list
      this.elements[elemPermutation[1]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[1]].setNodes;
      this.elements[elemPermutation[1]].tags[1] = 1;
      this.elements[elemPermutation[1]].nodes[1] = nodePermutation[1];

      // Add right boundary point to elements list
      this.elements[elemPermutation[this.nElements]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[this.nElements]].setNodes;
      this.elements[elemPermutation[this.nElements]].tags[1] = 2;
      this.elements[elemPermutation[this.nElements]].nodes[1] = nodePermutation[this.nNodes];
    }
  }

  record gmesh_element_r
  {
    var elemType : int;
    var nTags : int = 1;
    var tags : [1..1] int;

    var nodes_d : domain(rank=1, idxType=int);
    var nodes : [nodes_d] int;

    proc setNodes
    {
      use Parameters.Gmesh;

      select this.elemType {
        when GMESH_PNT_1  do this.nodes_d = {1..1};
        when GMESH_LIN_2  do this.nodes_d = {1..2};
        when GMESH_LIN_3  do this.nodes_d = {1..3};
        when GMESH_LIN_4  do this.nodes_d = {1..4};
        when GMESH_TRI_3  do this.nodes_d = {1..3};
        when GMESH_TRI_6  do this.nodes_d = {1..6};
        when GMESH_TRI_9  do this.nodes_d = {1..9};
        when GMESH_TRI_10 do this.nodes_d = {1..10};
        when GMESH_QUA_4  do this.nodes_d = {1..4};
        when GMESH_QUA_9  do this.nodes_d = {1..9};
        when GMESH_QUA_8  do this.nodes_d = {1..8};
        when GMESH_QUA_16 do this.nodes_d = {1..16};
      }
    }
  }

  record gmesh_family_r
  {
    var nDim : int;
    var name : string;
  }

  proc main()
  {
    use Random;

    var mesh1 = new gmesh_c(nNodes=7, nElements=8, nFamilies=3);
    mesh1.random1D(-1,2);

    writeln("Test 1: Random 1D mesh:");
    writeln(mesh1);
  }
}
