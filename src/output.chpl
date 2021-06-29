prototype module Output
{
  import FRMesh.fr_mesh_c;

  proc iterOutput(nIter : int, fr_mesh : fr_mesh_c)
  {
    use IO;
    use Path;

    var stringIter : string = nIter:string;
    var outputDir  : string = curDir;

    // Create directory for this iteration output

    // Write all selected output files
    writeln("Writing output of iteration %7i".format(nIter));
    output_gnuplot(outputDir, "sol_sp_gnuplt", stringIter, fr_mesh.xyzSP, fr_mesh.solSP, true, true);
    output_gnuplot(outputDir, "res_sp_gnuplt", stringIter, fr_mesh.xyzSP, fr_mesh.resSP);
    output_gnuplot(outputDir, "sol_fp1_gnuplt", stringIter, fr_mesh.xyzFP, fr_mesh.solFP[..,1,..]);
    output_gnuplot(outputDir, "sol_fp2_gnuplt", stringIter, fr_mesh.xyzFP, fr_mesh.solFP[..,2,..]);

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

  proc output_gnuplot(outputDir : string, fileRoot : string, fileSulfix : string, xyz : [] real, vars : [] real,
      flagPressure : bool = false, flagMach : bool = false)
  {
    use FileSystem;
    use IO;
    use Path;
    use Flux;

    const GNUPlotIOStyle = new iostyle(min_width_columns=15, showpointzero=0, showplus=1, precision=6, realfmt=2);
    //param fileRoot : string = "sol_gnuplt";
    param fileExt  : string = ".dat";
    param nameSep  : string = "-";

    var fileName   : string = fileRoot + nameSep + fileSulfix + fileExt;
    var outputFile : file;

    try {
      outputFile = open(outputDir + pathSep + fileName , iomode.cw, style=GNUPlotIOStyle);
    } catch {
      stdout.writeln("Unknown Error opening output file.");
      stderr.writeln("Unknown Error opening output file.");
    }

    var outputChan = outputFile.writer();

    // Channel style:
    //   binary = 0,
    //   byteorder = 1,
    //   str_style = -65280,
    //   min_width_columns = 15,
    //   max_width_columns = 4294967295,
    //   max_width_characters = 4294967295,
    //   max_width_bytes = 4294967295,
    //   string_start = 34,
    //   string_end = 34,
    //   string_format = 0,
    //   bytes_prefix = 98,
    //   base = 0,
    //   point_char = 46,
    //   exponent_char = 101,
    //   other_exponent_char = 112,
    //   positive_char = 43,
    //   negative_char = 45,
    //   i_char = 105,
    //   prefix_base = 1,
    //   pad_char = 32,
    //   showplus = 1,
    //   uppercase = 0,
    //   leftjustify = 0,
    //   showpoint = 0,
    //   showpointzero = 0,
    //   precision = 6,
    //   realfmt = 2,
    //   complex_style = 0,
    //   array_style = 0,
    //   aggregate_style = 0,
    //   tuple_style = 0

    // Write header
    outputChan.writef("      #DOF");
    for dimIdx in xyz.domain.dim(1) do
      outputChan.writef("    Coordinate-%i", dimIdx);
    for varIdx in vars.domain.dim(1) do
      outputChan.writef("      Variable-%i", varIdx);
    if flagPressure then
      outputChan.writef("        Pressure");
    if flagMach then
      outputChan.writef("            Mach");
    outputChan.writeln();

    // Write values
    for dof in xyz.domain.dim(0)
    {
      outputChan.writef("%{10i}", dof);
      for dimIdx in xyz.domain.dim(1) do
        outputChan.writef("%{+16.6er}", xyz[dof,dimIdx]);
      for varIdx in vars.domain.dim(1) do
        outputChan.writef("%{+16.6er}", vars[dof,varIdx]);
      if flagPressure then
        outputChan.writef("%{+16.6er}", pressure_cv(vars[dof,..]));
      if flagMach then
        outputChan.writef("%{+16.6er}", mach_cv(vars[dof,..]));
      outputChan.writeln();
    }

    outputChan.close();
    outputFile.close();
  }
}

