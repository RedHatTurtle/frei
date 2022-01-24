module CGNS
{
  use Random;
  use UnitTest;

  // Planned supported mesh formats
  //
  // - 1D Internal generated mesh  (Frei 1.0)
  // - 1D/2D cgns version 2/4     (Frei 2.0)
  // - 2D Internal generated cgns (Frei 2.x)
  // - 3D cgns version 2/4        (Frei 3.0)
  // - 3D CGNS                     (Frei 4.0)
  // - 2D Internal generated CGNS  (Frei 4.x)
  // - 3D Internal generated CGNS  (Frei 4.x)
  // - 3D Internal generated cgns (Frei 4.x)
  //
  // It would probably be better to have a separate module for each IO standard (CGNS, cgns, Plot3D, etc). At the moment shit's
  // just pilled up here for compactness since most of it is just general concepts and sketches.
  //
  // The idea behind having a separate class for each IO format it to isolate the conversion from the solver mesh format to the IO
  //    format from the reading and writing of the IO files. It might waste to much memory for just code organization's sake and
  //    should be reviewed in the future. 

  class cgns_mesh_c
  {
    var cgnsVersion   : string;
    var cellDimension : int;
    var physDimension : int;
    var cgnsZones     : [1..1] owned cgns_zone_c;
  }

  class cgns_zone_c
  {
    var vertexSize         : int;
    var cellSize           : int;
    var vertexSizeBoundary : int;

    var gridCoordinates    : [vertexSize, physDimension] real;

    // Options
    //   - Irregular Array with variable size in second dimension
    //   - Separate arrays for each element type
    //   - Array of class with individual allocation for node list
    //   - Array of arrays with individual allocation for node list
    var elementList : int;
  }

  record cgns_element_r
  {
    var elemType : int;
    var nTags : int;
    var tags : [1..nTags] int;

    var nodeDomain : domain(rank=1, idxType=int);
    var nodes : [nodeDomain] int;

    proc setNodes
    {
      use Parameters.cgns;

      select this.elemType {
        when CGNS_BAR_2  do this.nodeDomain = {1..2};
        when CGNS_BAR_3  do this.nodeDomain = {1..3};
        when CGNS_BAR_4  do this.nodeDomain = {1..4};
        when CGNS_TRI_3  do this.nodeDomain = {1..3};
        when CGNS_TRI_6  do this.nodeDomain = {1..6};
        when CGNS_TRI_9  do this.nodeDomain = {1..9};
        when CGNS_TRI_10 do this.nodeDomain = {1..10};
        when CGNS_QUA_4  do this.nodeDomain = {1..4};
        when CGNS_QUA_9  do this.nodeDomain = {1..9};
        when CGNS_QUA_8  do this.nodeDomain = {1..8};
        when CGNS_QUA_16 do this.nodeDomain = {1..16};
      }
    }
  }

  record cgns_family_r
  {
    var nDim : int;
    var name : string;
  }

  proc main()
  {
    use Random;

    var mesh1 = new cgns_c(nNodes=7, nElements=6, nFamilies=3);
    mesh1.random1D(-1,2);

    writeln("Test 1: Random 1D mesh:");
    writeln(mesh1);
  }
}
