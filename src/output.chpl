prototype module Output
{
  import FRMesh.fr_mesh_c;

  proc iterOutput(nIter : int, fr_mesh : fr_mesh_c)
  {
    use Path;

    var stringIter : string = nIter:string;
    var outputDir  : string = curDir;

    // Create directory for this iteration output

    // Write all selected output files
    output_gnuplot(outputDir, stringIter, fr_mesh.xyzSP, fr_mesh.solSP);

    // Delete old output directories if necessary
  }

//  proc timeOutput(in curTime : real)
//  {
//    var stringTime : string;
//    var dirName : string;
//
//    // Define the string format to write the current time on the file/dir name string
//    // Build string with the file/directory name
//    write(stringTime,'(f0.4)') curTime;
//    dirName = 'time'+trim(stringTime);
//
//    call execute_command_line('mkdir -p '+trim(dirName) );
//    call stdOut(dirName);
//  }

  proc output_gnuplot(outputDir : string, fileSulfix : string, xyz : [] real, vars : [] real)
  {
    use FileSystem;
    use IO;
    use Path;
    use Flux;

    const GNUPlotIOStyle = new iostyle(min_width_columns=15, showpointzero=0, showplus=1, precision=6, realfmt=2);
    param fileRoot : string = "sol_gnuplt";
    param fileExt  : string = ".dat";
    param nameSep  : string = "_";

    var fileName   : string = fileRoot + nameSep + fileSulfix + fileExt;
    var outputFile : file;

    try {
      outputFile = open(outputDir + pathSep + fileName , iomode.cw, style=GNUPlotIOStyle);
    } catch {
      stdout.writeln("Unknown Error opening output file.");
      stderr.writeln("Unknown Error opening output file.");
    }

    var outputChan = outputFile.writer();

    for dof in xyz.domain.dim(0) do
      outputChan.writeln(dof, xyz[dof,..], vars[dof,..], pressure_cv(vars[dof,..]), mach_cv(vars[dof,..]));

    outputChan.close();
    outputFile.close();
  }
}

