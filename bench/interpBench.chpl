module InterpBench
{
  use Time;
  use Random;
  use LinearAlgebra;
  param TIMEUNIT = TimeUnits.milliseconds;

  // Input variables
  config var nDOF  : int = 131072;
  config var nReps : int = 1;
  config var nDims : int = 2;
  config var nIters : int = 1;
  config var nCells : int = 0;
  config var minDegree : int = 1;
  config var maxDegree : int = 8;
  config var modeAuto : bool = true;
  config var testManual : bool = true;
  config var testLinAlg : bool = true;

  proc main()
  {
    use IO;

    var testTime : real = 0;
    var compTime : real = 0;
    var minTest : int =  0;
    var maxTest : int = 16;

    if testManual == false then minTest = 9;
    if testLinAlg == false then maxTest = 8;

    var cnt  : [minDegree..maxDegree, minTest..maxTest] int  = 0;
    var sum1 : [minDegree..maxDegree, minTest..maxTest] real = 0;
    var sum2 : [minDegree..maxDegree, minTest..maxTest] real = 0;

    // Run tests
    if modeAuto
    {
      writeln("Running on automatic mode:");
      var outputTimer : Timer;
      outputTimer.start();

      while modeAuto
      {
        // Find the test we spent less time on
        var (minTime, minTimeTest) = minloc reduce zip(sum1, sum1.domain);

        var interpDegree = minTimeTest(0);
        var cellCnt : int;
        if nDOF > 0 then
          cellCnt = nDOF/(interpDegree+1)**nDims;
        else
          cellCnt = nCells;

        select minTimeTest(1)
        {
          when  0 do (compTime, testTime) = manual_dense_0(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  1 do (compTime, testTime) = manual_dense_1(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  2 do (compTime, testTime) = manual_dense_2(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  3 do (compTime, testTime) = manual_dense_3(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  4 do (compTime, testTime) = manual_dense_4(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  5 do (compTime, testTime) = manual_dense_5(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  6 do (compTime, testTime) = manual_dense_6(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  7 do (compTime, testTime) = manual_dense_7(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  8 do (compTime, testTime) = manual_dense_8(   nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when  9 do (compTime, testTime) = LinAlgDot_dense_1(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when 10 do (compTime, testTime) = LinAlgDot_dense_2(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when 11 do (compTime, testTime) = LinAlgDot_dense_3(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when 12 do (compTime, testTime) = LinAlgDot_dense_4(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when 13 do (compTime, testTime) = LinAlgDot_dense_5(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when 14 do (compTime, testTime) = LinAlgDot_dense_6(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when 15 do (compTime, testTime) = LinAlgDot_dense_7(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
          when 16 do (compTime, testTime) = LinAlgDot_dense_8(nDims, nIters, cellCnt, interpDegree, seed=cnt[minTimeTest]+1);
        }

        cnt[ minTimeTest(0), minTimeTest(1)] += 1;
        sum1[minTimeTest(0), minTimeTest(1)] += testTime;
        sum2[minTimeTest(0), minTimeTest(1)] += testTime**2;
        writef("    - Test %2i, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
            minTimeTest(1), cellCnt, minTimeTest(0), testTime, compTime);

        if outputTimer.elapsed(TimeUnits.seconds) > 60
        {
          dump_stats(sum1=sum1, sum2=sum2, cnt=cnt);
          outputTimer.clear();
        }
      }
    }
    else
    {
      for rep in 1..nReps
      {
        writeln();
        writeln("Starting test rep ", rep);
        for interpDegree in minDegree..maxDegree
        {
          writeln("   Degree: ", interpDegree);
          var cellCnt : int;
          if nDOF > 0 then
            cellCnt = nDOF/(interpDegree+1)**nDims;
          else
            cellCnt = nCells;

          // Manual implementations
          if testManual
          {
            (compTime, testTime) = manual_dense_1(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  1, interpDegree] += 1;
            sum1[ 1, interpDegree] += testTime;
            sum2[ 1, interpDegree] += testTime**2;
            writef("      - Manual    Dense 1: %{ 9.3er} ms\n", testTime);
            writef("    - Test  1, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = manual_dense_2(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  2, interpDegree] += 1;
            sum1[ 2, interpDegree] += testTime;
            sum2[ 2, interpDegree] += testTime**2;
            writef("    - Test  2, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = manual_dense_3(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  3, interpDegree] += 1;
            sum1[ 3, interpDegree] += testTime;
            sum2[ 3, interpDegree] += testTime**2;
            writef("    - Test  3, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = manual_dense_4(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  4, interpDegree] += 1;
            sum1[ 4, interpDegree] += testTime;
            sum2[ 4, interpDegree] += testTime**2;
            writef("    - Test  4, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = manual_dense_5(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  5, interpDegree] += 1;
            sum1[ 5, interpDegree] += testTime;
            sum2[ 5, interpDegree] += testTime**2;
            writef("    - Test  5, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = manual_dense_6(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  6, interpDegree] += 1;
            sum1[ 6, interpDegree] += testTime;
            sum2[ 6, interpDegree] += testTime**2;
            writef("    - Test  6, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = manual_dense_7(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  7, interpDegree] += 1;
            sum1[ 7, interpDegree] += testTime;
            sum2[ 7, interpDegree] += testTime**2;
            writef("    - Test  7, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = manual_dense_8(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  8, interpDegree] += 1;
            sum1[ 8, interpDegree] += testTime;
            sum2[ 8, interpDegree] += testTime**2;
            writef("    - Test  8, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);
          }

          // LinearAlgebra implementations
          if testLinAlg
          {
            (compTime, testTime) = LinAlgDot_dense_1(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[  9, interpDegree] += 1;
            sum1[ 9, interpDegree] += testTime;
            sum2[ 9, interpDegree] += testTime**2;
            writef("    - Test  9, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = LinAlgDot_dense_2(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[ 10, interpDegree] += 1;
            sum1[10, interpDegree] += testTime;
            sum2[10, interpDegree] += testTime**2;
            writef("    - Test 10, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = LinAlgDot_dense_3(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[ 11, interpDegree] += 1;
            sum1[11, interpDegree] += testTime;
            sum2[11, interpDegree] += testTime**2;
            writef("    - Test 11, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = LinAlgDot_dense_4(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[ 12, interpDegree] += 1;
            sum1[12, interpDegree] += testTime;
            sum2[12, interpDegree] += testTime**2;
            writef("    - Test 12, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = LinAlgDot_dense_5(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[ 13, interpDegree] += 1;
            sum1[13, interpDegree] += testTime;
            sum2[13, interpDegree] += testTime**2;
            writef("    - Test 13, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = LinAlgDot_dense_6(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[ 14, interpDegree] += 1;
            sum1[14, interpDegree] += testTime;
            sum2[14, interpDegree] += testTime**2;
            writef("    - Test 14, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = LinAlgDot_dense_7(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[ 15, interpDegree] += 1;
            sum1[15, interpDegree] += testTime;
            sum2[15, interpDegree] += testTime**2;
            writef("    - Test 15, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);

            (compTime, testTime) = LinAlgDot_dense_8(nDims, nIters, cellCnt, interpDegree, seed=rep);
            cnt[ 16, interpDegree] += 1;
            sum1[16, interpDegree] += testTime;
            sum2[16, interpDegree] += testTime**2;
            writef("    - Test 16, %8i cells @ degree %2i ran in: %{ 9.3er} ms total with %{ 9.3er} ms of computation\n",
                cellCnt, interpDegree, testTime, compTime);
          }
        }
      }
    }

    // Print times
    dump_stats(sum1=sum1, sum2=sum2, cnt=cnt);
  }

  proc dump_stats(fileName : string = "interpBench.out", sum1 : [] real, sum2 : [] real, cnt : [] int,
      dumpPerIter : bool = true)
  {
    use FileSystem;
    use Path;
    use IO;

    var cellCnt : int;

    // Open a file
    var outputFile : file;
    var outputChan : channel;
    try! outputFile = open(fileName, iomode.cw);
    try! outputChan = outputFile.writer();

    // Print input parameters
    try! outputChan.writeln("Input Parameters");
    try! outputChan.writeln("    - nDOF      : ", nDOF);
    try! outputChan.writeln("    - nReps     : ", nReps);
    try! outputChan.writeln("    - nDims     : ", nDims);
    try! outputChan.writeln("    - nIters    : ", nIters);
    try! outputChan.writeln("    - nCells    : ", nCells);
    try! outputChan.writeln("    - minDegree : ", minDegree);
    try! outputChan.writeln("    - maxDegree : ", maxDegree);
    try! outputChan.writeln("    - modeAuto  : ", modeAuto);
    try! outputChan.writeln("    - testManual: ", testManual);
    try! outputChan.writeln("    - testLinAlg: ", testLinAlg);

    // Print the number of reps run of each test
    try! outputChan.writeln();
    try! outputChan.writeln("Test reps count:");
    for interpDegree in minDegree..maxDegree
    {
      if nDOF > 0 then
        cellCnt = nDOF/(interpDegree+1)**nDims;
      else
        cellCnt = nCells;
      try! outputChan.writef("    Degree %2i, %8i Cells", interpDegree, cellCnt);
      for testIdx in cnt.dim(1) do
        try! outputChan.writef(", %9i", cnt[interpDegree, testIdx]);
      try! outputChan.writef("\n");
    }

    // Print the average result of each test
    try! outputChan.writeln();
    try! outputChan.writeln("Test reps average:");
    for interpDegree in minDegree..maxDegree
    {
      if nDOF > 0 then
        cellCnt = nDOF/(interpDegree+1)**nDims;
      else
        cellCnt = nCells;
      try! outputChan.writef("    Degree %2i, %8i Cells", interpDegree, cellCnt);
      for testIdx in cnt.dim(1) do
        try! outputChan.writef(", %{9.3er}", sum1[interpDegree, testIdx]/cnt[interpDegree, testIdx]);
      try! outputChan.writef("\n");
    }

    // Print the standard deviation of each test
    try! outputChan.writeln();
    try! outputChan.writeln("Test reps Std-Dev:");
    for interpDegree in minDegree..maxDegree
    {
      if nDOF > 0 then
        cellCnt = nDOF/(interpDegree+1)**nDims;
      else
        cellCnt = nCells;
      try! outputChan.writef("    Degree %2i, %8i Cells", interpDegree, cellCnt);
      for testIdx in cnt.dim(1)
      {
        var avg : real = sum1[interpDegree, testIdx]/cnt[interpDegree, testIdx];
        var stdDev : real = sqrt((sum2[interpDegree, testIdx] - 2*avg*sum1[interpDegree, testIdx])/cnt[interpDegree, testIdx] + avg**2);
        try! outputChan.writef(", %{9.3er}", stdDev);
      }
      try! outputChan.writef("\n");
    }

    // Print the ration of the standard deviation by the average duration of each test
    try! outputChan.writeln();
    try! outputChan.writeln("Test reps Std-Dev / Average:");
    for interpDegree in minDegree..maxDegree
    {
      if nDOF > 0 then
        cellCnt = nDOF/(interpDegree+1)**nDims;
      else
        cellCnt = nCells;
      try! outputChan.writef("    Degree %2i, %8i Cells", interpDegree, cellCnt);
      for testIdx in cnt.dim(1)
      {
        var avg : real = sum1[interpDegree, testIdx]/cnt[interpDegree, testIdx];
        var stdDev : real = sqrt((sum2[interpDegree, testIdx] - 2*avg*sum1[interpDegree, testIdx])/cnt[interpDegree, testIdx] + avg**2);
        try! outputChan.writef(", %{8.2dr}%%", 100*stdDev/avg);
      }
      try! outputChan.writef("\n");
    }

    if dumpPerIter
    {
      // Print the number of reps run of each test
        try! outputChan.writeln();
        try! outputChan.writeln("Test reps time/iter average:");
        for interpDegree in minDegree..maxDegree
        {
        if nDOF > 0 then
          cellCnt = nDOF/(interpDegree+1)**nDims;
        else
          cellCnt = nCells;
        try! outputChan.writef("    Degree %2i, %8i Cells", interpDegree, cellCnt);
        for testIdx in cnt.dim(1) do
          try! outputChan.writef(", %{9.3er}", sum1[interpDegree, testIdx]/cnt[interpDegree, testIdx]/nIters);
        try! outputChan.writef("\n");
      }

      // Print the number of reps run of each test
      try! outputChan.writeln();
      try! outputChan.writeln("Test reps time/iter Std-Dev:");
      for interpDegree in minDegree..maxDegree
      {
        if nDOF > 0 then
          cellCnt = nDOF/(interpDegree+1)**nDims;
        else
          cellCnt = nCells;
        try! outputChan.writef("    Degree %2i, %8i Cells", interpDegree, cellCnt);
        for testIdx in cnt.dim(1)
        {
          var avg : real = sum1[interpDegree, testIdx]/cnt[interpDegree, testIdx];
          var stdDev : real = sqrt((sum2[interpDegree, testIdx] - 2*avg*sum1[interpDegree, testIdx])/cnt[interpDegree, testIdx] + avg**2);
          try! outputChan.writef(", %{9.3er}", stdDev/nIters);
        }
        try! outputChan.writef("\n");
      }
    }

    // Close file
    try! outputChan.close();
    try! outputFile.close();
  }

  /////////////////////////////////
  // Manual Loop Implementations //
  /////////////////////////////////

  // Loops: vfs,  F[v,f] += S[v,s]*M[f,s]
 proc manual_dense_0(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for varIdx in 1..nDims+2 do
        for cellIdx in 1..nCells
        {
          var spFirst : int = (cellIdx-1)*cellSPcnt;

            for fpIdx in 1..cellFPcnt do
              for spIdx in 1..cellSPcnt do
                varFP[varIdx, (cellIdx-1)*cellFPcnt+fpIdx] += varSP[varIdx, spFirst+spIdx] * sp2fp[fpIdx, spIdx];
        }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: vfs,  F[v,f] += S[v,s]*M[f,s]
  proc manual_dense_2(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for varIdx in 1..nDims+2 do
          for fpIdx in 1..cellFPcnt do
            for spIdx in 1..cellSPcnt do
              varFP[varIdx, (cellIdx-1)*cellFPcnt+fpIdx] += varSP[varIdx, spFirst+spIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: vfs,  F[f,v] += S[v,s]*M[f,s]
  proc manual_dense_6(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for varIdx in 1..nDims+2 do
          for fpIdx in 1..cellFPcnt do
            for spIdx in 1..cellSPcnt do
              varFP[(cellIdx-1)*cellFPcnt+fpIdx, varIdx] += varSP[varIdx, spFirst+spIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: fvs,  F[f,v] += S[v,s]*M[f,s]
  proc manual_dense_5(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          for varIdx in 1..nDims+2 do
            for spIdx in 1..cellSPcnt do
              varFP[(cellIdx-1)*cellFPcnt+fpIdx, varIdx] += varSP[varIdx, spFirst+spIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: fvs,  F[v,f] += S[v,s]*M[f,s]
  proc manual_dense_1(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          for varIdx in 1..nDims+2 do
            for spIdx in 1..cellSPcnt do
              varFP[varIdx, (cellIdx-1)*cellFPcnt+fpIdx] += varSP[varIdx, spFirst+spIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: vfs,  F[v,f] += S[s,v]*M[f,s]
  proc manual_dense_4(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for varIdx in 1..nDims+2 do
          for fpIdx in 1..cellFPcnt do
            for spIdx in 1..cellSPcnt do
              varFP[varIdx, (cellIdx-1)*cellFPcnt+fpIdx] += varSP[spFirst+spIdx, varIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: vfs,  F[f,v] += S[s,v]*M[f,s]
  proc manual_dense_8(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    //  FP[f,v] = SP[s,v] * M[f,s]
    //  Loop: fvs
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for varIdx in 1..nDims+2 do
          for fpIdx in 1..cellFPcnt do
            for spIdx in 1..cellSPcnt do
              varFP[(cellIdx-1)*cellFPcnt+fpIdx, varIdx] += varSP[spFirst+spIdx, varIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: fvs,  F[f,v] += S[s,v]*M[f,s]
  proc manual_dense_7(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          for varIdx in 1..nDims+2 do
            for spIdx in 1..cellSPcnt do
              varFP[(cellIdx-1)*cellFPcnt+fpIdx, varIdx] += varSP[spFirst+spIdx, varIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Loops: fvs,  F[v,f] += S[s,v]*M[f,s]
  proc manual_dense_3(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          for varIdx in 1..nDims+2 do
            for spIdx in 1..cellSPcnt do
              varFP[varIdx, (cellIdx-1)*cellFPcnt+fpIdx] += varSP[spFirst+spIdx, varIdx]*sp2fp[fpIdx, spIdx];
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  ////////////////////////////////////
  // Linear Algebra Implementations //
  ////////////////////////////////////

  proc LinAlgDot_dense_1(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          varFP[.., (cellIdx-1)*cellFPcnt+fpIdx] = dot(varSP[.., spFirst..spLast], sp2fp[fpIdx, ..]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  proc LinAlgDot_dense_3(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          varFP[.., (cellIdx-1)*cellFPcnt+fpIdx] = dot(sp2fp[fpIdx, ..], varSP[spFirst..spLast, ..]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  proc LinAlgDot_dense_5(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          varFP[(cellIdx-1)*cellFPcnt+fpIdx, ..] = dot(varSP[.., spFirst..spLast], sp2fp[fpIdx, ..]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  proc LinAlgDot_dense_7(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx*cellSPcnt;

        for fpIdx in 1..cellFPcnt do
          varFP[(cellIdx-1)*cellFPcnt+fpIdx, ..] = dot(sp2fp[fpIdx, ..], varSP[spFirst..spLast, ..]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  // Transposed varSP alternates

  proc LinAlgDot_dense_2(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx    *cellSPcnt;

        var fpFirst : int = (cellIdx-1)*cellFPcnt+1;
        var fpLast  : int =  cellIdx   *cellFPcnt;

        for varIdx in 1..nDims+2 do
          varFP[varIdx, fpFirst..fpLast] = dot(sp2fp, varSP[varIdx, spFirst..spLast]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  proc LinAlgDot_dense_4(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nDims+2, 1..nCells*cellFPcnt] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx*cellSPcnt;

        var fpFirst : int = (cellIdx-1)*cellFPcnt+1;
        var fpLast  : int =  cellIdx   *cellFPcnt;

        for varIdx in 1..nDims+2 do
          varFP[varIdx, fpFirst..fpLast] = dot(sp2fp, varSP[spFirst..spLast, varIdx]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  proc LinAlgDot_dense_6(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nDims+2, 1..nCells*cellSPcnt] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx*cellSPcnt;

        var fpFirst : int = (cellIdx-1)*cellFPcnt+1;
        var fpLast  : int =  cellIdx   *cellFPcnt;

        for varIdx in 1..nDims+2 do
          varFP[fpFirst..fpLast, varIdx] = dot(sp2fp, varSP[varIdx, spFirst..spLast]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  proc LinAlgDot_dense_8(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  {
    var interpTimer : Timer;
    interpTimer.start();

    // Hard coded for quadrilaterals
    var nFaces    : int = 2*nDims;
    var cellSPcnt : int =        (interpDegree+1)**(nDims);
    var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

    // Allocate calculation arrays
    var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
    var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
    var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

    fillRandom(sp2fp, seed);
    fillRandom(varSP, seed);
    fillRandom(varFP, seed);

    var initTime : real = interpTimer.elapsed(TIMEUNIT);
    // Do lots of dot products
    for iteration in 1..nIters do
      for cellIdx in 1..nCells
      {
        var spFirst : int = (cellIdx-1)*cellSPcnt+1;
        var spLast  : int = cellIdx*cellSPcnt;

        var fpFirst : int = (cellIdx-1)*cellFPcnt+1;
        var fpLast  : int =  cellIdx   *cellFPcnt;

        for varIdx in 1..nDims+2 do
          varFP[fpFirst..fpLast, varIdx] = dot(sp2fp, varSP[spFirst..spLast, varIdx]);
      }
    var compTime : real = interpTimer.elapsed(TIMEUNIT);

    return (compTime-initTime, compTime);
  }

  //proc LinAlgDot_dense_9(nDims, nIters, nCells, interpDegree, seed) : (2*real)
  //{
  //  var interpTimer : Timer;
  //  interpTimer.start();

  //  // Hard coded for quadrilaterals
  //  var nFaces    : int = 2*nDims;
  //  var cellSPcnt : int =        (interpDegree+1)**(nDims);
  //  var cellFPcnt : int = nFaces*(interpDegree+1)**(nDims-1);

  //  // Allocate calculation arrays
  //  var sp2fp : [1..cellFPcnt, 1..cellSPcnt]      real;
  //  var varSP : [1..nCells*cellSPcnt, 1..nDims+2] real;
  //  var varFP : [1..nCells*cellFPcnt, 1..nDims+2] real;

  //  fillRandom(sp2fp, seed);
  //  fillRandom(varSP, seed);
  //  fillRandom(varFP, seed);

  //  var initTime : real = interpTimer.elapsed(TIMEUNIT);
  //  // Do lots of dot products
  //  for iteration in 1..nIters do
  //    for cellIdx in 1..nCells
  //    {
  //      var spFirst : int = (cellIdx-1)*cellSPcnt+1;
  //      var spLast  : int = cellIdx*cellSPcnt;

  //      var fpFirst : int = (cellIdx-1)*cellFPcnt+1;
  //      var fpLast  : int =  cellIdx   *cellFPcnt;

  //      for fpIdx in 1..cellFPcnt do
  //        varFP[(cellIdx-1)*cellFPcnt+fpIdx, ..] = dot(sp2fp[, ..], varSP[spFirst..spLast, ..]);
  //    }
  //  var compTime : real = interpTimer.elapsed(TIMEUNIT);

  //  return (compTime-initTime, compTime);
  //}
}
