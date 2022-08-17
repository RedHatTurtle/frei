module interpBenchMod
    use iso_fortran_env
    implicit none

    ! Input variables
    integer :: nDOF  = 131072
    integer :: nReps = 1
    integer :: nDims = 2
    integer :: nIters = 1
    integer :: nCells = 0
    integer :: minDegree = 1
    integer :: maxDegree = 8
    logical :: modeAuto  = .true.
    logical :: testManual = .true.
    logical :: testMATMUL = .true.
    logical :: testLAPACK = .true.

    ! Indirect configuration variables
    integer :: minTest = 1
    integer :: maxTest = 1

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Output Test statistics !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine dump_stats(fileName, sum1, sum2, cnt, dumpPerIter)
        use iso_fortran_env
        implicit none

        ! Input variables
        character (len=*) :: fileName
        integer, dimension(:,:) :: cnt(minDegree:maxDegree, minTest:maxTest)
        real(real64), dimension(:,:) :: sum1(minDegree:maxDegree, minTest:maxTest)
        real(real64), dimension(:,:) :: sum2(minDegree:maxDegree, minTest:maxTest)
        logical :: dumpPerIter

        ! Auxiliary varaibles
        integer :: interpDegree, testIdx, cellCnt
        real(real64) :: avg, stdDev

        ! Open output file
        open(unit = 100, file = fileName)

        ! Print the input parameters
        write(100,*) "Input Parameters"
        write(100,*) "    - nDOF      : ", nDOF
        write(100,*) "    - nReps     : ", nReps
        write(100,*) "    - nDims     : ", nDims
        write(100,*) "    - nIters    : ", nIters
        write(100,*) "    - nCells    : ", nCells
        write(100,*) "    - minDegree : ", minDegree
        write(100,*) "    - maxDegree : ", maxDegree
        write(100,*) "    - modeAuto  : ", modeAuto
        write(100,*) "    - testManual: ", testManual
        write(100,*) "    - testLinAlg: ", testMATMUL
        write(100,*) "    - testLinAlg: ", testLAPACK

        ! Print the number of reps run of each test
        write(100,*)
        write(100,*) "Test reps count:"
        do interpDegree = minDegree, maxDegree
            if (nDOF > 0) then
                cellCnt = nDOF/((interpDegree+1)**nDims)
            endif
            write(100, "(A, I2, A, I8, A)", advance="no") "    Degree ", interpDegree, ", ", cellCnt, " Cells"
            do testIdx = minTest, maxTest
                write(100, "(A, I9)", advance="no") ", ", cnt(interpDegree, testIdx)
            enddo
            write(100,*)
        enddo

        ! Print the average result of each test
        write(100,*)
        write(100,*) "Test reps average:"
        do interpDegree = minDegree, maxDegree
            if (nDOF > 0) then
                cellCnt = nDOF/((interpDegree+1)**nDims)
            endif
            write(100, "(A, I2, A, I8, A)", advance="no") "    Degree ", interpDegree, ", ", cellCnt, " Cells"
            do testIdx = minTest, maxTest
                avg = sum1(interpDegree, testIdx)/cnt(interpDegree, testIdx)
                write(100, "(A, E9.3E2)", advance="no") ", ", avg
            enddo
            write(100,*)
        enddo

        ! Print the standard deviation of each test
        write(100,*)
        write(100,*) "Test reps Std-Dev:"
        do interpDegree = minDegree, maxDegree
            if (nDOF > 0) then
                cellCnt = nDOF/((interpDegree+1)**nDims)
            endif
            write(100, "(A, I2, A, I8, A)", advance="no") "    Degree ", interpDegree, ", ", cellCnt, " Cells"
            do testIdx = minTest, maxTest
                ! stddev = sqrt((sum(x**2)-sum(x)**2/size)/(size-1))
                avg = sum1(interpDegree, testIdx)/cnt(interpDegree, testIdx)
                stdDev = sqrt(   (sum2(interpDegree, testIdx) - 2*avg*sum1(interpDegree, testIdx))                      &
                    /(cnt(interpDegree, testIdx)-1.0) + avg**2)
                write(100, "(A, E9.3E2)", advance="no") ", ", stdDev
            enddo
            write(100,*)
        enddo

        ! Print the ration of the standard deviation by the average duration of each test
        write(100,*)
        write(100,*) "Test reps Std-Dev / Average:"
        do interpDegree = minDegree, maxDegree
            if (nDOF > 0) then
                cellCnt = nDOF/((interpDegree+1)**nDims)
            endif
            write(100, "(A, I2, A, I8, A)", advance="no") "    Degree ", interpDegree, ", ", cellCnt, " Cells"
            do testIdx = minTest, maxTest
                avg = sum1(interpDegree, testIdx)/cnt(interpDegree, testIdx)
                stdDev = sqrt(   (sum2(interpDegree, testIdx) - 2*avg*sum1(interpDegree, testIdx))                      &
                    /(cnt(interpDegree, testIdx)-1.0) + avg**2)
                write(100, "(A, F8.2, A)", advance="no") ", ", 100.0_real64*stdDev/avg, "%"
            enddo
            write(100,*)
        enddo

        ! Close file
        close(100)
    end subroutine dump_stats

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Manual Loop Implementations !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Loops: vfs,  F[v,f] += S[v,s]*M[f,s]
    real(real64) function manual_dense_0(nDims, nIters, nCells, interpDegree, seed) result(exec_time)
        use iso_fortran_env
        implicit none

        integer, intent(in) :: nDims, nIters, nCells, interpDegree, seed

        integer :: seed_size
        integer, allocatable :: seed_array(:)

        real(real64), allocatable :: sp2fp(:,:), varSP(:,:), varFP(:,:)

        integer :: nFaces, CellSPcnt, cellFPcnt
        integer :: iter, varIdx, cellIdx, fpIdx, spIdx
        integer :: spFirst, fpFirst
        integer(int64) :: t_start, t_init, t_end

        call system_clock(t_start)

        !! Hard coded for quadrilaterals
        nFaces    = 2*nDims
        cellSPcnt =        (interpDegree+1)**(nDims)
        cellFPcnt = nFaces*(interpDegree+1)**(nDims-1)

        !! Allocate calculation arrays
        allocate(sp2fp(cellSPcnt, cellFPcnt))
        allocate(varSP(nCells*cellSPcnt, nDims+2))
        allocate(varFP(nCells*cellFPCnt, nDims+2))

        ! Set up repeatable seed
        call random_seed(size=seed_size)
        allocate( seed_array(1:seed_size) )
        seed_array = seed
        call random_seed(put=seed_array)

        call random_number(sp2fp)
        call random_number(varSP)
        call random_number(varFP)

        call system_clock(t_init)
        !! Do lots of dot products
        do iter=1, nIters
            do varIdx=1, nDims+2
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    fpFirst = (cellIdx-1)*cellFPcnt+1

                    do fpIdx=1, cellFPcnt
                        do spIdx=1, cellSPcnt
                            varFP(fpFirst+fpIdx-1, varIdx) = sp2fp(spIdx, fpIdx) * varSP(spFirst+spIdx-1, varIdx)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call system_clock(t_end)

        exec_time = real(t_end-t_init, real64)/(1000000.0_real64)
    end function manual_dense_0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! DOT_PRODUCT Loop Implementations !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Loops: vfs,  F[v,f] += S[v,s]*M[f,s]
    real(real64) function dotprod_dense_0(nDims, nIters, nCells, interpDegree, seed) result(exec_time)
        use iso_fortran_env
        implicit none

        integer, intent(in) :: nDims, nIters, nCells, interpDegree, seed

        integer :: seed_size
        integer, allocatable :: seed_array(:)

        real(real64), allocatable :: sp2fp(:,:), varSP(:,:), varFP(:,:)

        integer :: nFaces, CellSPcnt, cellFPcnt
        integer :: iter, varIdx, cellIdx, fpIdx, spIdx
        integer :: spFirst, spLast
        integer(int64) :: t_start, t_init, t_end

        call system_clock(t_start)

        !! Hard coded for quadrilaterals
        nFaces    = 2*nDims
        cellSPcnt =        (interpDegree+1)**(nDims)
        cellFPcnt = nFaces*(interpDegree+1)**(nDims-1)

        !! Allocate calculation arrays
        allocate(sp2fp(cellSPcnt, cellFPcnt))
        allocate(varSP(nCells*cellSPcnt, nDims+2))
        allocate(varFP(nCells*cellFPCnt, nDims+2))

        ! Set up repeatable seed
        call random_seed(size=seed_size)
        allocate( seed_array(1:seed_size) )
        seed_array = seed
        call random_seed(put=seed_array)

        call random_number(sp2fp)
        call random_number(varSP)
        call random_number(varFP)

        call system_clock(t_init)
        !! Do lots of dot products
        do iter=1, nIters
            do varIdx=1, nDims+2
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    spLast  =  cellIdx   *cellSPcnt

                    do fpIdx=1, cellFPcnt
                        varFP((cellIdx-1)*cellFPcnt+fpIdx, varIdx) = dot_product(sp2fp(:, fpIdx), varSP(spFirst:spLast, varIdx))
                    enddo
                enddo
            enddo
        enddo
        call system_clock(t_end)

        exec_time = real(t_end-t_init, real64)/(1000000.0_real64)
    end function dotprod_dense_0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! MATMUL Loop Implementations !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Loops: vfs,  F[v,f] += S[v,s]*M[f,s]
    real(real64) function matmul_dense_0(nDims, nIters, nCells, interpDegree, seed) result(exec_time)
        use iso_fortran_env
        implicit none

        integer, intent(in) :: nDims, nIters, nCells, interpDegree, seed

        integer :: seed_size
        integer, allocatable :: seed_array(:)

        real(real64), allocatable :: sp2fp(:,:), varSP(:,:), varFP(:,:)

        integer :: nFaces, CellSPcnt, cellFPcnt
        integer :: iter, varIdx, cellIdx, fpIdx, spIdx
        integer :: spFirst, spLast
        integer :: fpFirst, fpLast
        integer(int64) :: t_start, t_init, t_end

        call system_clock(t_start)

        !! Hard coded for quadrilaterals
        nFaces    = 2*nDims
        cellSPcnt =        (interpDegree+1)**(nDims)
        cellFPcnt = nFaces*(interpDegree+1)**(nDims-1)

        !! Allocate calculation arrays
        allocate(sp2fp(cellSPcnt, cellFPcnt))
        allocate(varSP(nCells*cellSPcnt, nDims+2))
        allocate(varFP(nCells*cellFPCnt, nDims+2))

        ! Set up repeatable seed
        call random_seed(size=seed_size)
        allocate( seed_array(1:seed_size) )
        seed_array = seed
        call random_seed(put=seed_array)

        call random_number(sp2fp)
        call random_number(varSP)
        call random_number(varFP)

        call system_clock(t_init)
        !! Do lots of dot products
        do iter=1, nIters
            do varIdx=1, nDims+2
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    spLast  =  cellIdx   *cellSPcnt

                    fpFirst = (cellIdx-1)*cellFPcnt+1
                    fpLast  =  cellIdx   *cellFPcnt

                    varFP(fpFirst:fpLast, varIdx) = matmul(transpose(sp2fp(:, :)), varSP(spFirst:spLast, varIdx))
                enddo
            enddo
        enddo
        call system_clock(t_end)

        exec_time = real(t_end-t_init, real64)/(1000000.0_real64)
    end function matmul_dense_0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! LAPACK Loop Implementations !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Loops: vfs,  F[v,f] += S[v,s]*M[f,s]
    real(real64) function dgemmT_dense_0(nDims, nIters, nCells, interpDegree, seed) result(exec_time)
        use iso_fortran_env
        implicit none

        integer, intent(in) :: nDims, nIters, nCells, interpDegree, seed

        integer :: seed_size
        integer, allocatable :: seed_array(:)

        real(real64), allocatable :: sp2fp(:,:), varSP(:,:), varFP(:,:)

        integer :: nFaces, CellSPcnt, cellFPcnt
        integer :: iter, varIdx, cellIdx, fpIdx, spIdx
        integer :: spFirst, spLast
        integer :: fpFirst, fpLast
        integer(int64) :: t_start, t_init, t_end

        call system_clock(t_start)

        !! Hard coded for quadrilaterals
        nFaces    = 2*nDims
        cellSPcnt =        (interpDegree+1)**(nDims)
        cellFPcnt = nFaces*(interpDegree+1)**(nDims-1)

        !! Allocate calculation arrays
        allocate(sp2fp(cellSPcnt, cellFPcnt))
        allocate(varSP(nCells*cellSPcnt, nDims+2))
        allocate(varFP(nCells*cellFPCnt, nDims+2))

        ! Set up repeatable seed
        call random_seed(size=seed_size)
        allocate( seed_array(1:seed_size) )
        seed_array = seed
        call random_seed(put=seed_array)

        call random_number(sp2fp)
        call random_number(varSP)
        call random_number(varFP)

        call system_clock(t_init)
        !! Do lots of dot products
        do iter=1, nIters
            do varIdx=1, nDims+2
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    spLast  =  cellIdx   *cellSPcnt

                    fpFirst = (cellIdx-1)*cellFPcnt+1
                    fpLast  =  cellIdx   *cellFPcnt

                    call dgemm('T', 'N', cellFPcnt, 1, cellSPcnt, 0, sp2fp(:,:), cellSPcnt,         &
                        varSP(spFirst:spLast, varIdx), cellSPcnt, 0, varFP(fpFirst:fpLast, varIdx), cellFPcnt)
                enddo
            enddo
        enddo
        call system_clock(t_end)

        exec_time = real(t_end-t_init, real64)/(1000000.0_real64)
    end function dgemmT_dense_0
end module interpBenchMod

program InterpBench
    use iso_fortran_env
    use interpBenchMod
    implicit none

    ! Old computation variables
    real(real64), allocatable :: sp2fp(:,:), varSP(:,:), varFP(:,:)
    integer :: cellSPcnt, cellFPcnt
    integer :: cellIdx, fpIdx
    integer :: spFirst, spLast
    integer :: interpDegree
    integer :: nFaces
    integer :: nTest

    integer :: rep, seed

    integer(int64) :: t_start, t_end, iter

    ! Old results array
    real(real64), allocatable :: results(:,:,:) ! Degree, Test, Repetition

    ! Automatic mode varaibles
    real(real64) :: minTime, testTime
    integer :: testIdx, testDegree, minTimeTest, testCells

    ! Results array
    real(real64), dimension(:,:), allocatable :: sum1, sum2
    integer, dimension(:,:), allocatable :: cnt

    write(*,*) "Number of Dimensions"
    read(*,*)  nDims
    write(*,*) nDims, " dimensions"
    write(*,*) "Number of Iterations"
    read(*,*)  nIters
    write(*,*) nIters, "iterations"
    write(*,*) "Number of DOFs (0 to define number of cells)"
    read(*,*)  nDOF
    write(*,*) nDOF, " DOFs"
    if (nDOF == 0) then
        write(*,*) "Undefined number of DOFs"
        modeAuto = .false.
        write(*,*) "How many Cells then?"
        read(*,*)  nCells
        write(*,*) nCells, " cells"
    endif
    write(*,*) "Min Polynomial Degree"
    read(*,*)  minDegree
    write(*,*) "    From", minDegree, " degree"
    write(*,*) "Max Polynomial Degree"
    read(*,*)  maxDegree
    write(*,*) "    To", maxDegree, " degree"
    write(*,*) "Number of test repetitions"
    read(*,*)  nReps
    write(*,*) nReps, " repetitions"

    minTest = 1
    maxTest = 4

    allocate(results(minDegree:maxDegree, 1:4, 1:nReps))
    allocate(sum1(minDegree:maxDegree, minTest:maxTest))
    allocate(sum2(minDegree:maxDegree, minTest:maxTest))
    allocate( cnt(minDegree:maxDegree, minTest:maxTest))

    results = 0.0_real64
    sum1 = 0.0_real64
    sum2 = 0.0_real64
    cnt  = 0

    if (modeAuto) then
        write(*,*)
        write(*,*) "Running on automatic mode:"

        do while (modeAuto)
            ! Find the test we spent less time on
            minTime = huge(minTime)
            do testDegree = minDegree, maxDegree
                do testIdx = minTest, maxTest
                    if ( sum1(testDegree, testIdx) < minTime ) then
                        minTime = sum1(testDegree, testIdx)
                        interpDegree = testDegree
                        minTimeTest = testIdx
                    endif
                enddo
            enddo

            ! Calculate the number of cells if constant DOFs
            if (nDOF > 0) then
                testCells = nDOF/((interpDegree+1)**2)
            endif

            select case (minTimeTest)
              case (1)
                testTime = manual_dense_0(nDims, nIters, testCells, interpDegree, cnt(interpDegree, minTimeTest))
              case (2)
                testTime = dotprod_dense_0(nDims, nIters, testCells, interpDegree, cnt(interpDegree, minTimeTest))
              case (3)
                testTime = matmul_dense_0(nDims, nIters, testCells, interpDegree, cnt(interpDegree, minTimeTest))
              case (4)
                testTime = dgemmT_dense_0(nDims, nIters, testCells, interpDegree, cnt(interpDegree, minTimeTest))
            endselect

            sum1(interpDegree, minTimeTest) = sum1(interpDegree, minTimeTest) + testTime
            sum2(interpDegree, minTimeTest) = sum2(interpDegree, minTimeTest) + testTime*testTime
            cnt(interpDegree, minTimeTest)  =  cnt(interpdegree, mintimetest) + 1

            write(*, "(A, I8, A, I2, A, I2, A, E9.3E2, A)") "Repetition ", cnt(interpDegree, minTimeTest), " of test ", &
                minTimeTest, " @ degree ", interpDegree, " in ", testTime, " ms"

            call dump_stats(fileName="interpBench.fortran.out", sum1=sum1, sum2=sum2, cnt=cnt, dumpPerIter=.true.)
        enddo
    else
        write(*,*) "Starting test with ", nDims, " dimensions"
        do rep=1, nReps
            write(*,*)
            write(*,*) "Starting test rep ", rep
            do interpDegree=minDegree, maxDegree
                write(*,*) "    Degree: ", interpDegree

                nFaces    = 2*nDims
                cellSPcnt =        (interpDegree+1)**(nDims)
                cellFPcnt = nFaces*(interpDegree+1)**(nDims-1)

                allocate(sp2fp(cellSPcnt, cellFPcnt))
                allocate(varSP(nCells*cellSPcnt, nDims+2))
                allocate(varFP(nCells*cellFPCnt, nDims+2))

                call random_number(sp2fp)
                call random_number(varSP)
                call random_number(varFP)

                call system_clock(t_start)
                do iter=1, nIters
                    do cellIdx=1, nCells
                        spFirst = (cellIdx-1)*cellSPcnt+1
                        spLast  = cellIdx*cellSPcnt

                        do fpIdx=1, cellFPcnt
                            varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(sp2fp(:, fpIdx), varSP(spFirst:spLast, :))
                        end do
                    end do
                enddo
                call system_clock(t_end)
                results(interpDegree, 1, rep) = real(t_end-t_start, real64)/(1000000.0_real64*real(nIters, real64))
                write(*,"(A, ES10.4E1, A)") "        Test 1: ", results(interpDegree, 1, rep), " ms avg per iteration"

                deallocate(sp2fp)
                deallocate(varSP)
                deallocate(varFP)

                allocate(sp2fp(cellSPcnt, cellFPcnt))
                allocate(varSP(nCells*cellSPcnt, nDims+2))
                allocate(varFP(nCells*cellFPCnt, nDims+2))

                call random_number(sp2fp)
                call random_number(varSP)
                call random_number(varFP)

                call system_clock(t_start)
                do iter=1, nIters
                    do cellIdx=1, nCells
                        spFirst = (cellIdx-1)*cellSPcnt+1
                        spLast  = cellIdx*cellSPcnt

                        do fpIdx=1, cellFPcnt
                            varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(transpose(varSP(spFirst:spLast, :)), sp2fp(:, fpIdx))
                        end do
                    end do
                enddo
                call system_clock(t_end)
                results(interpDegree, 2, rep) = real(t_end-t_start, real64)/(1000000.0_real64*real(nIters, real64))
                write(*,"(A, ES10.4E1, A)") "        Test 2: ", results(interpDegree, 2, rep), " ms avg per iteration"

                deallocate(sp2fp)
                deallocate(varSP)
                deallocate(varFP)

                allocate(sp2fp(cellSPcnt, cellFPcnt))
                allocate(varSP(nDims+2, nCells*cellSPcnt))
                allocate(varFP(nCells*cellFPCnt, nDims+2))

                call random_number(sp2fp)
                call random_number(varSP)
                call random_number(varFP)

                call system_clock(t_start)
                do iter=1, nIters
                    do cellIdx=1, nCells
                        spFirst = (cellIdx-1)*cellSPcnt+1
                        spLast  = cellIdx*cellSPcnt

                        do fpIdx=1, cellFPcnt
                            varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(varSP(:, spFirst:spLast), sp2fp(:, fpIdx))
                        end do
                    end do
                enddo
                call system_clock(t_end)
                results(interpDegree, 3, rep) = real(t_end-t_start, real64)/(1000000.0_real64*real(nIters, real64))
                write(*,"(A, ES10.4E1, A)") "        Test 3: ", results(interpDegree, 3, rep), " ms avg per iteration"

                deallocate(sp2fp)
                deallocate(varSP)
                deallocate(varFP)

                allocate(sp2fp(cellSPcnt, cellFPcnt))
                allocate(varSP(nDims+2, nCells*cellSPcnt))
                allocate(varFP(nCells*cellFPCnt, nDims+2))

                call random_number(sp2fp)
                call random_number(varSP)
                call random_number(varFP)

                ! Test 4
                ! This seems to scale a bit better with larger nDims+2 than test 1 and scales a lot better than 1 with larger interpDegree
                call system_clock(t_start)
                do iter=1, nIters
                    do cellIdx=1, nCells
                        spFirst = (cellIdx-1)*cellSPcnt+1
                        spLast  = cellIdx*cellSPcnt

                        do fpIdx=1, cellFPcnt
                            varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(sp2fp(:, fpIdx), transpose(varSP(:, spFirst:spLast)))
                        end do
                    end do
                enddo
                call system_clock(t_end)
                results(interpDegree, 4, rep) = real(t_end-t_start, real64)/(1000000.0_real64*real(nIters, real64))
                write(*,"(A, ES10.4E1, A)") "        Test 4: ", results(interpDegree, 4, rep), " ms avg per iteration"

                deallocate(sp2fp)
                deallocate(varSP)
                deallocate(varFP)
            enddo
        enddo
    endif

    call dump_stats(fileName="interpBench.fortran.out", sum1=sum1, sum2=sum2, cnt=cnt, dumpPerIter=.true.)

end program InterpBench
