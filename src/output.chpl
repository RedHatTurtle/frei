prototype module Output
{
  use IO;
  import FRMesh.fr_mesh_c;

  proc log_convergence(convergenceLogChan : channel, iteration : int, lInfDelta : [] real, l2Delta : [] real, l1Delta : [] real,
                       lInfRelativeDelta : [] real, l2RelativeDelta : [] real, l1RelativeDelta : [] real)
  {
    use IO;

    const header : [1..6] string = ["Linf(ΔSol[%1i])",   "L2(ΔSol[%1i])",   "L1(ΔSol[%1i])",
                                    "Linf(%%ΔSol[%1i])", "L2(%%ΔSol[%1i])", "L1(%%ΔSol[%1i])"];

    // Write header on first iteration
    if iteration == 1
    {
      convergenceLogChan.writef("%10s", "Iteration");
      for metric in header do
        for varIdx in l2Delta.domain do
          convergenceLogChan.writef("%16s", metric.format(varIdx));
      convergenceLogChan.writef("\n");
    }

    convergenceLogChan.writef("%10i", iteration);
    for varIdx in l2Delta.domain do
      convergenceLogChan.writef("%{ 16.8er}", lInfDelta[varIdx]);
    for varIdx in l2Delta.domain do
      convergenceLogChan.writef("%{ 16.8er}", l2Delta[varIdx]);
    for varIdx in l2Delta.domain do
      convergenceLogChan.writef("%{ 16.8er}", l1Delta[varIdx]);
    for varIdx in l2Delta.domain do
      convergenceLogChan.writef("%{ 16.8er}", lInfRelativeDelta[varIdx]);
    for varIdx in l2Delta.domain do
      convergenceLogChan.writef("%{ 16.8er}", l2RelativeDelta[varIdx]);
    for varIdx in l2Delta.domain do
      convergenceLogChan.writef("%{ 16.8er}", l1RelativeDelta[varIdx]);
    convergenceLogChan.writef("]\n");

    convergenceLogChan.flush();
  }

  // Select the adequate output routines that need to be run
  proc iterOutput(nIter : int, frMesh : fr_mesh_c)
  {
    use IO;
    use Path;

    var stringIter : string = nIter:string;
    var outputDir  : string = curDir;

    // Create directory for this iteration output

    // Write all selected output files
    if frMesh.nDims == 1
    {
      output_gnuplot(outputDir, "sol_sp_gnuplt", stringIter, frMesh.xyzSP, frMesh.solSP, true, true);
      output_gnuplot(outputDir, "res_sp_gnuplt", stringIter, frMesh.xyzSP, frMesh.resSP);
    }
    if frMesh.nDims == 2 then
      output_tecplot_dat(outputDir, "sol_ho_tecplot", stringIter, frMesh);

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
    outputChan.writef("       DOF");
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
        outputChan.writef("  %{ 14.7er}", xyz[dof,dimIdx]);
      for varIdx in vars.domain.dim(1) do
        outputChan.writef("  %{ 14.7er}", vars[dof,varIdx]);
      if flagPressure then
        outputChan.writef("  %{ 14.7er}", pressure_cv(vars[dof,..]));
      if flagMach then
        outputChan.writef("  %{ 14.7er}", mach_cv(vars[dof,..]));
      outputChan.writeln();
    }

    outputChan.close();
    outputFile.close();
  }

  proc output_tecplot_dat(outputDir : string, fileRoot : string, fileSulfix : string,
      frMesh : fr_mesh_c,
      flagPressure : bool = true, flagTemperature : bool = false, flagMach : bool = true, flagEntropy : bool = false,
      flagDouble : bool = false)
  {
    use FileSystem;
    use IO;
    use Path;
    use Flux;
    use Interpolation;
    use Parameters.ParamMesh;

    //const TecplotIOStyle = new iostyle(min_width_columns=15, showpointzero=0, showplus=1, precision=6, realfmt=2);
    //param fileRoot : string = "sol_gnuplt";
    param fileExt  : string = ".dat";
    param nameSep  : string = "-";

    var fileName   : string = fileRoot + nameSep + fileSulfix + fileExt;
    var outputFile : file;

    try {
      outputFile = open(outputDir + pathSep + fileName , iomode.cw);
    } catch {
      stdout.writeln("Unknown Error opening output file.");
      stderr.writeln("Unknown Error opening output file.");
    }

    var outputChan = outputFile.writer();

    // File Header
    {
      var tecplotTitle : string = "FREI OUTPUT";
      var tecplotFileType : string = "FULL";
      var tecplotVariables : string = "\"PointIdx\"";

      var dimName : [1..3] string = ["\"X\"", "\"Y\"", "\"Z\""];
      //if input.eqSet == EULER
      //{
      //  if frMesh.nDims == 3 then
      //    var varName : [1..frMesh.nVars] string = ["\"Density\"", "\"Momentum-X\"", "\"Energy\""];
      //  if frMesh.nDims == 4 then
          var varName : [1..frMesh.nVars] string = ["\"Density\"", "\"Momentum-X\"", "\"Momentum-Y\"", "\"Energy\""];
      //  if frMesh.nDims == 5 then
      //    var varName : [1..frMesh.nVars] string = ["\"Density\"", "\"Momentum-X\"", "\"Momentum-Y\"", "\"Momentum-Z\"", "\"Energy\""];
      //}

      for dimIdx in 1..frMesh.nDims do
        tecplotVariables += ", "+dimName[dimIdx];

      for varIdx in 1..frMesh.nVars do
        tecplotVariables += ", "+varName[varIdx];

      if flagPressure    then tecplotVariables += ", \"Pressure\"";
      if flagTemperature then tecplotVariables += ", \"Temperature\"";
      if flagMach        then tecplotVariables += ", \"Mach\"";
      if flagEntropy     then tecplotVariables += ", \"Entropy\"";

      outputChan.writef("TITLE = %\"S\n", tecplotTitle);
      outputChan.writef("FILETYPE = %\"S\n", tecplotFileType);
      outputChan.writef("VARIABLES = %s\n", tecplotVariables);
    }

    // Initially every family is combined into one single zone, but in the future more zones should be added
    // correnpoding to each mesh family.
    //
    // Ex 1, Ringleb Flow: Zone 1 = Inlet
    //                     Zone 2 = Outler
    //                     Zone 3 = Inner Wall
    //                     Zone 4 = Outter Wall
    //                     Zone 5 = Mesh Interior / Flow
    //
    // Ex 2, Airfoil Flow: Zone 1 = Farfield
    //                     Zone 2 = Airfoil
    //                     Zone 3 = Mesh Interior / Flow

    // Finite Element Zone Record (Mesh Interior / Flow)
    {
      // Control Line
      outputChan.writef("ZONE\n");

      // Zone data
      var spCnt : int = frMesh.xyzSP.domain.dim(0).high;
      var fpCnt : int = frMesh.xyzFP.domain.dim(0).high;
      var nodeCnt : int = frMesh.nodeList.domain.dim(0).high;
      var pointCnt : int = spCnt + fpCnt + nodeCnt;
      var elemCnt : int = frMesh.cellList.domain.dim(0).high*(frMesh.solOrder+2)**2;
      var realFormat : string = "  %{ 14.7er}";

      // Zone Header
      {
        var tecplotZoneTitle : string = "Flow";

        outputChan.writef("T = %\"S", tecplotZoneTitle);
        outputChan.writef(", ZONETYPE = FEQUADRILATERAL");
        outputChan.writef(", DATAPACKING = POINT"); // Default is BLOCK
        outputChan.writef(", NODES = %i, ELEMENTS = %i", pointCnt, elemCnt);
        outputChan.writef(", NV = 1"); // Position of the variable with the node index, assumed to be ordered
        outputChan.writef("\n");

        if flagDouble // Defualt is SINGLE
        {
          outputChan.writef("DT = (DOUBLE, DOUBLE, DOUBLE ...)");
          realFormat = "  %{ 22.15er}";
          outputChan.writef("\n");
        }
        //outputChan.writef(" VARLOCATION = (NODAL, NODAL, NODAL ...)"); // Default is NODAL
        //outputChan.writef("\n");
      }

      // Zone Data
      {
        var pointIdxFormat : string = "  %" + ceil(log10(pointCnt+1)):int:string + "i";

        // Loop through SPs
        for spIdx in frMesh.xyzSP.domain.dim(0)
        {
          outputChan.writef(pointIdxFormat, spIdx);

          // Loop through spacial coordinates
          for dimIdx in 1..frMesh.nDims do
            outputChan.writef(realFormat, frMesh.xyzSP[spIdx, dimIdx]);

          // Loop through conserved variables
          for varIdx in 1..frMesh.nVars do
            outputChan.writef(realFormat, frMesh.solSP[spIdx, varIdx]);

          if flagPressure    then outputChan.writef(realFormat, pressure_cv(frMesh.solSP[spIdx, ..]));
          if flagTemperature then outputChan.writef(realFormat, temperature_cv(frMesh.solSP[spIdx, ..]));
          if flagMach        then outputChan.writef(realFormat, mach_cv(frMesh.solSP[spIdx, ..]));
          if flagEntropy     then outputChan.writef(realFormat, entropy_cv(frMesh.solSP[spIdx, ..]));

          outputChan.writef("\n");
        }

        // Loop through FPs
        for fpIdx in frMesh.xyzFP.domain.dim(0)
        {
          outputChan.writef(pointIdxFormat, spCnt + fpIdx);

          // Loop through spacial coordinates
          for dimIdx in 1..frMesh.nDims do
            outputChan.writef(realFormat, frMesh.xyzFP[fpIdx, dimIdx]);

          // Loop through conserved variables
          for varIdx in 1..frMesh.nVars do
            outputChan.writef(realFormat, (frMesh.solFP[fpIdx, 1, varIdx]+frMesh.solFP[fpIdx, 2, varIdx])/2.0);

          if flagPressure    then outputChan.writef(realFormat, pressure_cv(frMesh.solFP[fpIdx, 1, ..]));
          if flagTemperature then outputChan.writef(realFormat, temperature_cv(frMesh.solFP[fpIdx, 1, ..]));
          if flagMach        then outputChan.writef(realFormat, mach_cv(frMesh.solFP[fpIdx, 1, ..]));
          if flagEntropy     then outputChan.writef(realFormat, entropy_cv(frMesh.solFP[fpIdx, 1, ..]));

          outputChan.writef("\n");
        }

        // Calculate the average solution on the nodes
        var cellCnt : [frMesh.nodeList.domain] int = 0;
        var solNode : [frMesh.nodeList.domain.dim(0), 1..frMesh.nVars] real = 0.0;
        for cellIdx in frMesh.cellList.domain
        {
          ref thisCell = frMesh.cellList[cellIdx];

          for cellNodeIdx in thisCell.nodes.domain
          {
            var meshNodeIdx : int = thisCell.nodes[cellNodeIdx];
            cellCnt[meshNodeIdx] += 1;
            solNode[meshNodeIdx, ..] += dot(sp2nodeInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellNodeIdx, ..],
                                        frMesh.solSP[frMesh.cellSPidx[cellIdx, 1]..#frMesh.cellSPidx[cellIdx, 2], ..]);
          }
        }

        // Loop through Nodes
        for nodeIdx in frMesh.nodeList.domain
        {
          outputChan.writef(pointIdxFormat, spCnt + fpCnt + nodeIdx);

          solNode[nodeIdx, ..] = solNode[nodeIdx, ..]/cellCnt[nodeIdx];

          // Loop through spacial coordinates
          for dimIdx in 1..frMesh.nDims do
            outputChan.writef(realFormat, frMesh.nodeList[nodeIdx].xyz[dimIdx]);

          // Loop through conserved variables
          for varIdx in 1..frMesh.nVars do
            outputChan.writef(realFormat, solNode[nodeIdx, varIdx]);

          if flagPressure    then outputChan.writef(realFormat, pressure_cv(solNode[nodeIdx, ..]));
          if flagTemperature then outputChan.writef(realFormat, temperature_cv(solNode[nodeIdx, ..]));
          if flagMach        then outputChan.writef(realFormat, mach_cv(solNode[nodeIdx, ..]));
          if flagEntropy     then outputChan.writef(realFormat, entropy_cv(solNode[nodeIdx, ..]));

          outputChan.writef("\n");
        }
      }

      // Connectivity Data
      {
        for cellIdx in frMesh.cellList.domain do
          select frMesh.cellList[cellIdx].elemTopo()
          {
            when TOPO_QUAD
            {
              var connectivityFormat : string = "  %" + ceil(log10(pointCnt+1)):int:string + "i"
                                               +"  %" + ceil(log10(pointCnt+1)):int:string + "i"
                                               +"  %" + ceil(log10(pointCnt+1)):int:string + "i"
                                               +"  %" + ceil(log10(pointCnt+1)):int:string + "i"
                                               +"\n";

              // Internal sub-cells
              for subCell in 1..frMesh.solOrder**2
              {
                var lin : int = (subCell-1)/(frMesh.solOrder);
                var col : int = (subCell-1)%(frMesh.solOrder);

                var sp1Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin  )*(frMesh.solOrder+1) + col;
                var sp2Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin  )*(frMesh.solOrder+1) + col+1;
                var sp3Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin+1)*(frMesh.solOrder+1) + col+1;
                var sp4Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin+1)*(frMesh.solOrder+1) + col;

                outputChan.writef(connectivityFormat, sp1Idx, sp2Idx, sp3Idx, sp4Idx);
              }

              // Sub-cells on edges, loop through faces
              for cellFaceIdx in 1..4
              {
                // Get the global mesh index of this face
                var meshFaceIdx = frMesh.cellList[cellIdx].faces[cellFaceIdx];

                for subCell in 1..frMesh.solOrder
                {
                  // Get the global index of two consecutive FPs on this face
                  var meshFP1Idx : int;
                  var meshFP2Idx : int;
                  if frMesh.cellList[cellIdx].sides[cellFaceIdx] == 1
                  {
                    meshFP1Idx = frMesh.faceFPidx[meshFaceIdx, 1] + (subCell-1);
                    meshFP2Idx = frMesh.faceFPidx[meshFaceIdx, 1] +  subCell;
                  }
                  else
                  {
                    meshFP1Idx = frMesh.faceFPidx[meshFaceIdx, 1] + frMesh.faceFPidx[meshFaceIdx, 2] -  subCell;
                    meshFP2Idx = frMesh.faceFPidx[meshFaceIdx, 1] + frMesh.faceFPidx[meshFaceIdx, 2] - (subCell+1);
                  }

                  // Now get the SPs
                  var sp1Idx : int;
                  var sp2Idx : int;
                  select cellFaceIdx
                  {
                    when 1
                    {
                      sp1Idx = frMesh.cellSPidx[cellIdx, 1] +  subCell;
                      sp2Idx = frMesh.cellSPidx[cellIdx, 1] + (subCell-1);
                    }
                    when 2
                    {
                      sp1Idx = frMesh.cellSPidx[cellIdx, 1] + (subCell  )*(frMesh.solOrder+1) + frMesh.solOrder;
                      sp2Idx = frMesh.cellSPidx[cellIdx, 1] + (subCell-1)*(frMesh.solOrder+1) + frMesh.solOrder;
                    }
                    when 3
                    {
                      sp1Idx = frMesh.cellSPidx[cellIdx, 1] + (frMesh.solOrder+1)*(frMesh.solOrder+1) - (subCell+1);
                      sp2Idx = frMesh.cellSPidx[cellIdx, 1] + (frMesh.solOrder+1)*(frMesh.solOrder+1) -  subCell;
                    }
                    when 4
                    {
                      sp1Idx = frMesh.cellSPidx[cellIdx, 1] + (frMesh.solOrder - subCell  )*(frMesh.solOrder+1);
                      sp2Idx = frMesh.cellSPidx[cellIdx, 1] + (frMesh.solOrder - subCell+1)*(frMesh.solOrder+1);
                    }
                  }

                  var fp1PntIdx : int = spCnt + meshFP1Idx;
                  var fp2PntIdx : int = spCnt + meshFP2Idx;
                  var sp1PntIdx : int = sp1Idx;
                  var sp2PntIdx : int = sp2Idx;

                  outputChan.writef(connectivityFormat, fp1PntIdx, fp2PntIdx, sp1PntIdx, sp2PntIdx);
                }
              }

              // Sub-cells on corners
              for cellNodeIdx in 1..4
              {
                // Get the global index of the corner node
                var meshNodeIdx : int = frMesh.cellList[cellIdx].nodes[cellNodeIdx];

                // Get the global index of the neighboring FP counterclockwise to the node
                var meshFace1Idx = frMesh.cellList[cellIdx].faces[cellNodeIdx];
                // Get first or last FP of the face depending on orientation
                var meshFP1Idx : int;
                if frMesh.cellList[cellIdx].sides[cellNodeIdx] == 1 then
                  meshFP1Idx = frMesh.faceFPidx[meshFace1Idx, 1]; // First FP from face
                else
                  meshFP1Idx = frMesh.faceFPidx[meshFace1Idx, 1] + frMesh.faceFPidx[meshFace1Idx, 2] - 1; // Last FP from face

                // Get the global index of the neighboring FP clockwise to the node
                var meshFace2Idx : int;
                var side : int;
                var meshFP2Idx : int;
                if cellNodeIdx == 1
                {
                  meshFace2Idx = frMesh.cellList[cellIdx].faces[4];
                  side = frMesh.cellList[cellIdx].sides[4];
                }
                else
                {
                  meshFace2Idx = frMesh.cellList[cellIdx].faces[cellNodeIdx-1];
                  side = frMesh.cellList[cellIdx].sides[cellNodeIdx-1];
                }
                // Get first or last FP of the face depending on orientation
                if side == 1 then
                  meshFP2Idx = frMesh.faceFPidx[meshFace2Idx, 1] + frMesh.faceFPidx[meshFace2Idx, 2] - 1; // Last FP from face
                else
                  meshFP2Idx = frMesh.faceFPidx[meshFace2Idx, 1]; // First FP from face

                var spIdx : int;
                select cellNodeIdx
                {
                  when 1 do spIdx = frMesh.cellSPidx[cellIdx, 1];
                  when 2 do spIdx = frMesh.cellSPidx[cellIdx, 1] + frMesh.solOrder;
                  when 3 do spIdx = frMesh.cellSPidx[cellIdx, 1] + frMesh.cellSPidx[cellIdx, 2] - 1;
                  when 4 do spIdx = frMesh.cellSPidx[cellIdx, 1] + (frMesh.solOrder+1)*(frMesh.solOrder);
                }

                var nodePntIdx : int = spCnt + fpCnt + meshNodeIdx;
                var fp1PntIdx  : int = spCnt + meshFP1Idx;
                var spPntIdx   : int = spIdx;
                var fp2PntIdx  : int = spCnt + meshFP2Idx;

                outputChan.writef(connectivityFormat, nodePntIdx, fp1PntIdx, spPntIdx, fp2PntIdx);
              }
            }
          }
      }

      // Zone Footer
      {}
    }

    // Text Record
    {}

    // Geometry Record
    {}

    // Custom Labels Record
    {}

    // Data Set Auxiliary Data Record
    {}

    // Variable Auxiliary Data Record
    {}
  }

  proc output_tecplot_plt(outputDir : string, fileRoot : string, fileSulfix : string, xyz : [] real, vars : [] real,
      flagPressure : bool = false, flagMach : bool = false)
  {}

  proc output_tecplot_szplt(outputDir : string, fileRoot : string, fileSulfix : string, xyz : [] real, vars : [] real,
      flagPressure : bool = false, flagMach : bool = false)
  {}
}

