prototype module InputMesh
{
  use Random;
  use UnitTest;

  // Planned supported mesh formats
  //
  // - 1D Internal generated mesh (Frei v1.0)
  // - 1D/2D Gmesh version 2/4    (Frei v2.0)
  // - 2D Internal generated mesh (Frei > v2.0)
  // - 3D Gmesh version 2/4       (Frei v3.0)
  // - 3D CGNS                    (Frei v4.0)
  // - 3D Internal generated mesh (Frei > v4.0)

  class input_gmesh_c
  {
    var nNodes : int;
    var nElements : int;
    var nFamilies : int;

    var nodes    : [1..nNodes] real;
    var elements : [1..nElements] int;
    var families : [1..nFamilies] string;
  }

  proc main()
  {
    use Testing;

    var mesh0 = new input_mesh_c();
    var mesh1 = new input_mesh_c(nNodes=3, nElements=2, nFamilies=1);
    var mesh2 = new input_mesh_c(nNodes=7, nElements=5, nFamilies=4);

    writeln(mesh0);
    writeln(mesh1);
    writeln(mesh2);
  }
}
