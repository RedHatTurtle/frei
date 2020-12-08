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
  // It would probably be better to have a separate module for each IO standard (CGNS, Gmesh, Plot3D, etc). At the moment shit's
  // just pilled up here for compactness since most of it is just general concepts and sketches.
  //
  // The idea behind having a separate class for each IO format it to isolate the conversion from the solver mesh format to the IO
  //    format from the reading and writing of the IO files. It might waste to much memory for just code organization's sake and
  //    should be reviewed in the future. 

  class gmesh_c
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
      use Parameters.Gmesh;
      import Sort.sort;

      var seed : int = 47;
      var randStreamSeeded = new RandomStream(real, seed);
      var x = this.nodes[..,1];
      var elem : [1..this.nElements, 1..2] int;
      var nodePermutation : [1..this.nNodes] int;
      var elemPermutation : [1..this.nElements] int;

      // Get random values from xMin to xMax
      x[1] = 0;
      x[2] = 1;
      randStreamSeeded.fillRandom(x[3..]);
      x = x*(xMax-xMin) + xMin;

      // Sort nodes to build elements and generate a random permutation
      sort(x);
      permutation(nodePermutation, seed);
      permutation(elemPermutation, seed);

      // Fill element list with no overlapping elements oriented from left to right
      for i in 1..this.nElements {
        elem[i,1] = i;
        elem[i,2] = i+1;
      }

      // Set the boundary and internal families
      this.families[1].nDim = 0;
      this.families[1].name = "left";
      this.families[2].nDim = 0;
      this.families[2].name = "right";
      this.families[3].nDim = 1;
      this.families[3].name = "flow";

      // Commit values to object
      for i in 1..this.nNodes do
        this.nodes[nodePermutation[i],1] = x[i];

      for i in 1..this.nElements {
        this.elements[elemPermutation[i]].elemType = GMESH_LIN_2;
        this.elements[elemPermutation[i]].setNodes;
        this.elements[elemPermutation[i]].nodes[1] = nodePermutation[elem[i,1]];
        this.elements[elemPermutation[i]].nodes[2] = nodePermutation[elem[i,2]];
      }
    }
  }

  record gmesh_element_r
  {
    var elemType : int;
    var nTags : int;
    var tags : [1..nTags] int;

    var nodeDomain : domain(rank=1, idxType=int);
    var nodes : [nodeDomain] int;

    proc setNodes
    {
      use Parameters.Gmesh;

      select this.elemType {
        when GMESH_LIN_2  do this.nodeDomain = {1..2};
        when GMESH_LIN_3  do this.nodeDomain = {1..3};
        when GMESH_LIN_4  do this.nodeDomain = {1..4};
        when GMESH_TRI_3  do this.nodeDomain = {1..3};
        when GMESH_TRI_6  do this.nodeDomain = {1..6};
        when GMESH_TRI_9  do this.nodeDomain = {1..9};
        when GMESH_TRI_10 do this.nodeDomain = {1..10};
        when GMESH_QUA_4  do this.nodeDomain = {1..4};
        when GMESH_QUA_9  do this.nodeDomain = {1..9};
        when GMESH_QUA_8  do this.nodeDomain = {1..8};
        when GMESH_QUA_16 do this.nodeDomain = {1..16};
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

    var mesh1 = new gmesh_c(nNodes=7, nElements=6, nFamilies=3);
    mesh1.random1D(-1,2);

    writeln("Test 1: Random 1D mesh:");
    writeln(mesh1);
  }
}
