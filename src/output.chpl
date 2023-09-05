module Output
{
  use IO;
  import FRMesh.fr_mesh_c;

  proc log_convergence(convergenceLogWriter : fileWriter, resToggle : bool = false,  iteration : int,
      l1Abs : [] real, l2Abs : [] real, lfAbs : [] real, l1Rel : [] real, l2Rel : [] real, lfRel : [] real)
  {
    use IO;

    const header : [1..6] string = if (resToggle) then
      ["L2(Res[%1i])",  "Lf(Res[%1i])",  "L1(Res[%1i])",  "L2(%%Res[%1i])",  "Lf(%%Res[%1i])",  "L1(%%Res[%1i])"]
    else
      ["L2(ΔSol[%1i])", "Lf(ΔSol[%1i])", "L1(ΔSol[%1i])", "L2(%%ΔSol[%1i])", "Lf(%%ΔSol[%1i])", "L1(%%ΔSol[%1i])"];

    // Open file writer and write to log
    try {
      // Write header on first iteration
      if iteration == 1
      {
        convergenceLogWriter.writef("%10s", "Iteration");
        for metric in header do
          for varIdx in l2Abs.domain do
            convergenceLogWriter.writef("%16s", metric.format(varIdx));
        convergenceLogWriter.writef("\n");
      }

      convergenceLogWriter.writef("%10i", iteration);
      for varIdx in l2Abs.domain do
        convergenceLogWriter.writef("%{ 16.8er}", l2Abs[varIdx]);
      for varIdx in lfAbs.domain do
        convergenceLogWriter.writef("%{ 16.8er}", lfAbs[varIdx]);
      for varIdx in l1Abs.domain do
        convergenceLogWriter.writef("%{ 16.8er}", l1Abs[varIdx]);
      for varIdx in l2Rel.domain do
        convergenceLogWriter.writef("%{ 16.8er}", l2Rel[varIdx]*100);
      for varIdx in lfRel.domain do
        convergenceLogWriter.writef("%{ 16.8er}", lfRel[varIdx]*100);
      for varIdx in l1Rel.domain do
        convergenceLogWriter.writef("%{ 16.8er}", l1Rel[varIdx]*100);
      convergenceLogWriter.writef("\n");

      convergenceLogWriter.flush();
    }
    catch {
      writeln("Failed to write to convergence log");
    }
  }

  proc print_log(logWriter : fileWriter, iteration : int, logVar : [] (string, real))
  {
    use IO;

    // Open file writer and write to log
    try {
      // Write header on first iteration
      if iteration == 1
      {
        logWriter.writef("%10s", "Iteration");
        for varIdx in logVar.domain do
          logWriter.writef(" %19s", logVar[varIdx](0));
        logWriter.writef("\n");
      }

      logWriter.writef("%10i", iteration);
      for varIdx in logVar.domain do
        logWriter.writef(" %{ 19.12er}", logVar[varIdx](1));
      logWriter.writef("\n");

      logWriter.flush();
    }
    catch {
      writeln("Failed to write log file log");
    }
  }

  // Select the adequate output routines that need to be run
  proc iterOutput(nIter : int, frMesh : fr_mesh_c, flagNormals : bool = false)
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
    if frMesh.nDims == 2 {
      output_fr_tecplot_dat(outputDir, "sol_ho_tecplot", stringIter, frMesh, flagNormals = flagNormals);
      output_fr_avg_tecplot_dat(outputDir, "sol_lo_tecplot", stringIter, frMesh);
    }

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
    import Input.fGamma;

    //param fileRoot : string = "sol_gnuplt";
    param fileExt  : string = ".dat";
    param nameSep  : string = "-";

    var fileName   : string = fileRoot + nameSep + fileSulfix + fileExt;
    var outputFile : file;

    try {
      outputFile = open(outputDir + pathSep + fileName , ioMode.cw);
    } catch {
      try! stdout.writeln("Unknown Error creating/opening output file.");
      try! stderr.writeln("Unknown Error creating/opening output file.");
    }

    // Open file writer and write gnuplot output
    try {
      var outputWriter = outputFile.writer();

      // Write header
      outputWriter.writef("       DOF");
      for dimIdx in xyz.domain.dim(1) do
        outputWriter.writef("    Coordinate-%i", dimIdx);
      for varIdx in vars.domain.dim(1) do
        outputWriter.writef("      Variable-%i", varIdx);
      if flagPressure then
        outputWriter.writef("        Pressure");
      if flagMach then
        outputWriter.writef("            Mach");
      outputWriter.writeln();

      // Write values
      for dof in xyz.domain.dim(0)
      {
        outputWriter.writef("%{10i}", dof);
        for dimIdx in xyz.domain.dim(1) do
          outputWriter.writef("  %{ 14.7er}", xyz[dof, dimIdx]);
        for varIdx in vars.domain.dim(0) do
          outputWriter.writef("  %{ 14.7er}", vars[varIdx, dof]);
        if flagPressure then
          outputWriter.writef("  %{ 14.7er}", pressure_cv(vars[.., dof], fGamma));
        if flagMach then
          outputWriter.writef("  %{ 14.7er}", mach_cv(vars[.., dof], fGamma));
        outputWriter.writeln();
      }

      outputWriter.close();
    }
    catch {
      writeln("Failed to write solution to GNUPlot file");
    }

    try! outputFile.close();
  }

  proc output_fr_tecplot_dat(outputDir : string, fileRoot : string, fileSulfix : string,
      frMesh : fr_mesh_c,
      flagPressure : bool = true, flagTemperature : bool = false, flagMach : bool = true, flagEntropy : bool = true,
      flagNormals : bool = false, flagDouble : bool = false)
  {
    use FileSystem;
    use IO;
    use Path;
    use Flux;
    use Interpolation;
    use Parameters.ParamMesh;
    import Math.log10;
    import Dimensional.scales;
    import Mesh.elem_vertices;
    import Input.fGamma;
    import Input.fR;

    //param fileRoot : string = "sol_gnuplt";
    param fileExt  : string = ".dat";
    param nameSep  : string = "-";

    var fileName   : string = fileRoot + nameSep + fileSulfix + fileExt;
    var outputFile : file;

    try {
      outputFile = open(outputDir + pathSep + fileName , ioMode.cw);
    } catch {
      try! stdout.writeln("Unknown Error creating/opening output file.");
      try! stderr.writeln("Unknown Error creating/opening output file.");
    }

    try {
      var outputWriter = outputFile.writer();

      // File Header
      {
        var tecplotTitle : string = "FREI OUTPUT";
        var tecplotFileType : string = "FULL";
        var tecplotVariables : string = "\"PointIdx\"";

        var dimName : [1..3] string = ["\"X\"", "\"Y\"", "\"Z\""];
        var nrmName : [1..3] string = ["\"Normal-X\"", "\"Normal-Y\"", "\"Normal-Z\""];
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

        if flagNormals then for dimIdx in 1..frMesh.nDims do
          tecplotVariables += ", "+nrmName[dimIdx];

        for varIdx in 1..frMesh.nVars do
          tecplotVariables += ", "+varName[varIdx];

        if flagPressure    then tecplotVariables += ", \"Pressure\"";
        if flagTemperature then tecplotVariables += ", \"Temperature\"";
        if flagMach        then tecplotVariables += ", \"Mach\"";
        if flagEntropy     then tecplotVariables += ", \"Entropy\"";

        outputWriter.writef("TITLE = %\"S\n", tecplotTitle);
        outputWriter.writef("FILETYPE = %\"S\n", tecplotFileType);
        outputWriter.writef("VARIABLES = %s\n", tecplotVariables);
      }

      // Initially every family is combined into one single zone, but in the future more zones should be added
      // corresponding to each mesh family.
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
        outputWriter.writef("ZONE\n");

        // Zone data
        var spCnt : int = frMesh.xyzSP.domain.dim(0).high;
        var fpCnt : int = frMesh.xyzFP.domain.dim(0).high;
        var elemCnt : int = frMesh.cellList.domain.dim(0).high*(frMesh.solOrder+2)**2;
        var realFormat : string = "  %{ 14.7er}";

        // Calculate the average solution on the vertices and index them
        var vertCnt : int = 0;
        var nodeCellCnt : [frMesh.nodeList.domain] int = 0;
        var nodeVertMap : [frMesh.nodeList.domain] int = 0;
        var solNode : [frMesh.nodeList.domain.dim(0), 1..frMesh.nVars] real = 0.0;
        for cellIdx in frMesh.cellList.domain
        {
          ref thisCell = frMesh.cellList[cellIdx];

          for cellNodeIdx in 1..elem_vertices(thisCell.elemTopo())
          {
            var meshNodeIdx : int = thisCell.nodes[cellNodeIdx];

            nodeCellCnt[meshNodeIdx] += 1;

            // If this is the first time we loop though this vertex then index it
            if nodeCellCnt[meshNodeIdx] == 1
            {
              vertCnt += 1;
              nodeVertMap[meshNodeIdx] = vertCnt;
            }

            solNode[meshNodeIdx, ..] += dot(frMesh.solSP[.., frMesh.cellSPidx[cellIdx, 1].. #frMesh.cellSPidx[cellIdx, 2]],
                                            sp2nodeInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellNodeIdx, ..]);
          }
        }
        var pointCnt : int = spCnt + fpCnt + vertCnt;

        // Zone Header
        {
          var tecplotZoneTitle : string = "Flow";

          outputWriter.writef("T = %\"S", tecplotZoneTitle);
          outputWriter.writef(", ZONETYPE = FEQUADRILATERAL");
          outputWriter.writef(", DATAPACKING = POINT"); // Default is BLOCK
          outputWriter.writef(", NODES = %i, ELEMENTS = %i", pointCnt, elemCnt);
          outputWriter.writef(", NV = 1"); // Position of the variable with the node index, assumed to be ordered
          outputWriter.writef("\n");

          if flagDouble // Default is SINGLE
          {
            outputWriter.writef("DT = (DOUBLE, DOUBLE, DOUBLE ...)");
            realFormat = "  %{ 22.15er}";
            outputWriter.writef("\n");
          }
          //outputWriter.writef(" VARLOCATION = (NODAL, NODAL, NODAL ...)"); // Default is NODAL
          //outputWriter.writef("\n");
        }

        // Zone Data
        {
          var pointIdxFormat : string = "  %" + ceil(log10(pointCnt+1)):int:string + "i";

          // Loop through SPs
          for spIdx in frMesh.xyzSP.domain.dim(0)
          {
            outputWriter.writef(pointIdxFormat, spIdx);

            // Loop through spatial coordinates
            for dimIdx in 1..frMesh.nDims do
              outputWriter.writef(realFormat, frMesh.xyzSP[spIdx, dimIdx]);

            if flagNormals then for dimIdx in 1..frMesh.nDims do
              outputWriter.writef(realFormat, 0.0);

            // Dimensionalize the conserved variables
            const dimConsVars : [1..frMesh.nVars] real = scales!.non2dim_cv(frMesh.solSP[.., spIdx]);

            // Loop through the dimensionalized conserved variables
            for varIdx in 1..frMesh.nVars do
              outputWriter.writef(realFormat, dimConsVars[varIdx]);

            if flagPressure    then outputWriter.writef(realFormat,    pressure_cv(dimConsVars, fGamma));
            if flagTemperature then outputWriter.writef(realFormat, temperature_cv(dimConsVars, fGamma, fR));
            if flagMach        then outputWriter.writef(realFormat,        mach_cv(dimConsVars, fGamma));
            if flagEntropy     then outputWriter.writef(realFormat,     entropy_cv(dimConsVars, fGamma));

            outputWriter.writef("\n");
          }

          // Loop through FPs
          for fpIdx in frMesh.xyzFP.domain.dim(0)
          {
            outputWriter.writef(pointIdxFormat, spCnt + fpIdx);

            // Average and dimensionalize variables
            const avgSol : [1..frMesh.nVars] real = (frMesh.solFP[fpIdx, 1, ..] + frMesh.solFP[fpIdx, 2, ..])/2.0;
            const dimAvgConsVars : [1..frMesh.nVars] real = scales!.non2dim_cv(avgSol);

            // Loop through spatial coordinates
            for dimIdx in 1..frMesh.nDims do
              outputWriter.writef(realFormat, frMesh.xyzFP[fpIdx, dimIdx]);

            if flagNormals then for dimIdx in 1..frMesh.nDims do
              outputWriter.writef(realFormat, frMesh.nrmFP[fpIdx, dimIdx]);

            // Loop through conserved variables
            for varIdx in 1..frMesh.nVars do
              outputWriter.writef(realFormat, dimAvgConsVars[varIdx]);

            if flagPressure    then outputWriter.writef(realFormat,    pressure_cv(dimAvgConsVars, fGamma));
            if flagTemperature then outputWriter.writef(realFormat, temperature_cv(dimAvgConsVars, fGamma, fR));
            if flagMach        then outputWriter.writef(realFormat,        mach_cv(dimAvgConsVars, fGamma));
            if flagEntropy     then outputWriter.writef(realFormat,     entropy_cv(dimAvgConsVars, fGamma));

            outputWriter.writef("\n");
          }

          // Loop through Nodes
          for nodeIdx in frMesh.nodeList.domain
          {
            // Skip nodes that aren't vertices
            if nodeCellCnt[nodeIdx] == 0 then continue;

            outputWriter.writef(pointIdxFormat, spCnt + fpCnt + nodeVertMap[nodeIdx]);

            // Average and dimensionalize variables
            const avgSol : [1..frMesh.nVars] real = solNode[nodeIdx, ..]/nodeCellCnt[nodeIdx];
            const dimAvgConsVars : [1..frMesh.nVars] real = scales!.non2dim_cv(avgSol);

            // Loop through spatial coordinates
            for dimIdx in 1..frMesh.nDims do
              outputWriter.writef(realFormat, frMesh.nodeList[nodeIdx].xyz[dimIdx]);

            if flagNormals then for dimIdx in 1..frMesh.nDims do
              outputWriter.writef(realFormat, 0.0);

            // Loop through conserved variables
            for varIdx in 1..frMesh.nVars do
              outputWriter.writef(realFormat, dimAvgConsVars[varIdx]);

            if flagPressure    then outputWriter.writef(realFormat,    pressure_cv(dimAvgConsVars, fGamma));
            if flagTemperature then outputWriter.writef(realFormat, temperature_cv(dimAvgConsVars, fGamma, fR));
            if flagMach        then outputWriter.writef(realFormat,        mach_cv(dimAvgConsVars, fGamma));
            if flagEntropy     then outputWriter.writef(realFormat,     entropy_cv(dimAvgConsVars, fGamma));

            outputWriter.writef("\n");
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
                for subCell in 1..(frMesh.solOrder)**2
                {
                  var lin : int = (subCell-1)/(frMesh.solOrder);
                  var col : int = (subCell-1)%(frMesh.solOrder);

                  var sp1Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin  )*(frMesh.solOrder+1) + col;
                  var sp2Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin  )*(frMesh.solOrder+1) + col+1;
                  var sp3Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin+1)*(frMesh.solOrder+1) + col+1;
                  var sp4Idx : int = frMesh.cellSPidx[cellIdx, 1] + (lin+1)*(frMesh.solOrder+1) + col;

                  outputWriter.writef(connectivityFormat, sp1Idx, sp2Idx, sp3Idx, sp4Idx);
                }

                // Sub-cells on edges, loop through quad faces
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

                    outputWriter.writef(connectivityFormat, fp1PntIdx, fp2PntIdx, sp1PntIdx, sp2PntIdx);
                  }
                }

                // Sub-cells on corners, loop through quad vertices
                for cellNodeIdx in 1..4
                {
                  // Get the global index of the corner node
                  var meshNodeIdx : int = frMesh.cellList[cellIdx].nodes[cellNodeIdx];
                  var vertIdx : int = nodeVertMap[meshNodeIdx];

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

                  var nodePntIdx : int = spCnt + fpCnt + vertIdx;
                  var fp1PntIdx  : int = spCnt + meshFP1Idx;
                  var spPntIdx   : int = spIdx;
                  var fp2PntIdx  : int = spCnt + meshFP2Idx;

                  outputWriter.writef(connectivityFormat, nodePntIdx, fp1PntIdx, spPntIdx, fp2PntIdx);
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
    catch {
      writeln("Failed to write solution to Tecplot file");
    }
  }

  proc output_fr_avg_tecplot_dat(outputDir : string, fileRoot : string, fileSulfix : string,
      frMesh : fr_mesh_c,
      flagPressure : bool = true, flagTemperature : bool = false, flagMach : bool = true, flagEntropy : bool = false,
      flagDouble : bool = false)
  {
    use FileSystem;
    use IO;
    use Path;
    use Flux;
    use Quadrature;
    use Parameters.ParamMesh;
    import Math.log10;
    import Input.fGamma;
    import Input.fR;
    import Dimensional.scales;

    //param fileRoot : string = "sol_gnuplt";
    param fileExt  : string = ".dat";
    param nameSep  : string = "-";
    var realFormat : string = "  %{ 14.7er}\n";

    var fileName   : string = fileRoot + nameSep + fileSulfix + fileExt;
    var outputFile : file;

    try {
      outputFile = open(outputDir + pathSep + fileName , ioMode.cw);
    } catch {
      try! stdout.writeln("Unknown Error creating/opening output file.");
      try! stderr.writeln("Unknown Error creating/opening output file.");
    }

    var vertCnt : int = 0;
    var nodeVertMap : [frMesh.nodeList.domain] int = 0;
    for thisCell in frMesh.cellList
    {
      // Loop though this cell's vertices
      for cellNodeIdx in 1..4//elem_vertices(thisCell.elemTopo())
      {
        // If this vertex has not been indexed yet
        if nodeVertMap[thisCell.nodes[cellNodeIdx]] == 0
        {
          vertCnt += 1;
          nodeVertMap[thisCell.nodes[cellNodeIdx]] = vertCnt;
        }
      }
    }

    try {
      var outputWriter = outputFile.writer();

      // File Header
      {
        var tecplotTitle : string = "FREI CELL_AVERAGED_OUTPUT";
        var tecplotFileType : string = "FULL";
        var tecplotVariables : string = "\"NodeIdx\"";

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

        tecplotVariables += ", \"CharLength\"";

        for varIdx in 1..frMesh.nVars do
          tecplotVariables += ", "+varName[varIdx];

        if flagPressure    then tecplotVariables += ", \"Pressure\"";
        if flagTemperature then tecplotVariables += ", \"Temperature\"";
        if flagMach        then tecplotVariables += ", \"Mach\"";
        if flagEntropy     then tecplotVariables += ", \"Entropy\"";

        outputWriter.writef("TITLE = %\"S\n", tecplotTitle);
        outputWriter.writef("FILETYPE = %\"S\n", tecplotFileType);
        outputWriter.writef("VARIABLES = %s\n", tecplotVariables);
      }

      // Initially every family is combined into one single zone, but in the future more zones should be added
      // correnpoding to each mesh family.
      //
      // Ex 1, Ringleb Flow: Zone 1 = Inlet
      //                     Zone 2 = Outlet
      //                     Zone 3 = Inner Wall
      //                     Zone 4 = Outer Wall
      //                     Zone 5 = Mesh Interior / Flow
      //
      // Ex 2, Airfoil Flow: Zone 1 = Farfield
      //                     Zone 2 = Airfoil
      //                     Zone 3 = Mesh Interior / Flow

      // Finite Element Zone Record (Mesh Interior / Flow)
      {
        // Control Line
        outputWriter.writef("ZONE\n");

        // Zone Header
        {
          var tecplotZoneTitle : string = "Flow";

          outputWriter.writef("T = %\"S", tecplotZoneTitle);
          outputWriter.writef(", ZONETYPE = FEQUADRILATERAL");
          outputWriter.writef(", DATAPACKING = BLOCK"); // Default is BLOCK
          outputWriter.writef(", VARLOCATION = ([4-10] = CELLCENTERED)"); // Default is NODAL
          outputWriter.writef(", NODES = %i, ELEMENTS = %i", vertCnt, frMesh.nCells);
          outputWriter.writef(", NV = 1"); // Position of the variable with the node index, assumed to be ordered
          outputWriter.writef("\n");

          if flagDouble // Default is SINGLE
          {
            outputWriter.writef("DT = (DOUBLE, DOUBLE, DOUBLE ...)");
            realFormat = "  %{ 22.15er}\n";
            outputWriter.writef("\n");
          }
          //outputWriter.writef("\n");
        }

        // Zone Data
        {
          var pointIdxFormat : string = "  %" + ceil(log10(vertCnt+1)):int:string + "i\n";

          // Loop through spatial coordinates printing locations of all nodes
          for nodeIdx in frMesh.nodeList.domain do
            if nodeVertMap[nodeIdx] != 0 then
              outputWriter.writef(pointIdxFormat, nodeVertMap[nodeIdx]);

          // Loop through spatial coordinates printing locations of all nodes
          for dimIdx in 1..frMesh.nDims do
            for nodeIdx in frMesh.nodeList.domain do
              if nodeVertMap[nodeIdx] != 0 then
                outputWriter.writef(realFormat, frMesh.nodeList[nodeIdx].xyz[dimIdx]);

          // Cell characteristic length
          for cellIdx in 1..frMesh.nCells do
            outputWriter.writef(realFormat, frMesh.cellCharLeng[cellIdx]);

          // Loop though all cells calculating averages of the conserved variables
          var solAvg : [1..frMesh.nCells, 1..frMesh.nVars] real;
          var dimSolAvg : [1..frMesh.nCells, 1..frMesh.nVars] real;

          for cellIdx in frMesh.cellList.domain
          {
            ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
            ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

            for varIdx in 1..frMesh.nVars do
              solAvg[cellIdx, varIdx] =
                average( frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                         frMesh.jacSP[        cellSPini.. #cellSPcnt],
                         frMesh.cellList[cellIdx].elemTopo(),
                         frMesh.solOrder
                );

            dimSolAvg[cellIdx, ..] = scales!.non2dim_cv(solAvg[cellIdx, ..]);
          }

          // Loop through conserved variables
          for varIdx in 1..frMesh.nVars do
          {
            for cellIdx in 1..frMesh.nCells do
              outputWriter.writef(realFormat, dimSolAvg[cellIdx, varIdx]);

            outputWriter.writef("\n");
          }

          if flagPressure
          {
            for cellIdx in 1..frMesh.nCells do
              outputWriter.writef(realFormat, pressure_cv(dimSolAvg[cellIdx, ..], fGamma));

            outputWriter.writef("\n");
          }

          if flagTemperature then
          {
            for cellIdx in 1..frMesh.nCells do
              outputWriter.writef(realFormat, temperature_cv(dimSolAvg[cellIdx, ..], fGamma, fR));

            outputWriter.writef("\n");
          }

          if flagMach then
          {
            for cellIdx in 1..frMesh.nCells do
              outputWriter.writef(realFormat, mach_cv(dimSolAvg[cellIdx, ..], fGamma));

            outputWriter.writef("\n");
          }

          if flagEntropy then
          {
            for cellIdx in 1..frMesh.nCells do
              outputWriter.writef(realFormat, entropy_cv(dimSolAvg[cellIdx, ..], fGamma));

            outputWriter.writef("\n");
          }

          outputWriter.writef("\n");
        }

        // Connectivity Data
        {
          for thisCell in frMesh.cellList do
            select thisCell.elemTopo()
            {
              when TOPO_QUAD
              {
                var connectivityFormat : string = "  %" + ceil(log10(vertCnt+1)):int:string + "i"
                                                 +"  %" + ceil(log10(vertCnt+1)):int:string + "i"
                                                 +"  %" + ceil(log10(vertCnt+1)):int:string + "i"
                                                 +"  %" + ceil(log10(vertCnt+1)):int:string + "i"
                                                 +"\n";

                // Get the global index of the corner node
                var node1Idx : int = nodeVertMap[thisCell.nodes[1]];
                var node2Idx : int = nodeVertMap[thisCell.nodes[2]];
                var node3Idx : int = nodeVertMap[thisCell.nodes[3]];
                var node4Idx : int = nodeVertMap[thisCell.nodes[4]];

                outputWriter.writef(connectivityFormat, node1Idx, node2Idx, node3Idx, node4Idx);
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
    catch {
      writeln("Failed to write solution to Tecplot file");
    }
  }


  proc output_tecplot_plt(outputDir : string, fileRoot : string, fileSulfix : string, xyz : [] real, vars : [] real,
      flagPressure : bool = false, flagMach : bool = false)
  {}

  proc output_tecplot_szplt(outputDir : string, fileRoot : string, fileSulfix : string, xyz : [] real, vars : [] real,
      flagPressure : bool = false, flagMach : bool = false)
  {}
}
