program InterpBench
    use iso_fortran_env
    implicit none

    integer(int32) :: nDims
    integer(int32) :: nFaces
    integer(int32) :: nIter
    integer(int32) :: nCells
    integer(int32) :: interpDegree

    real(real64), allocatable :: sp2fp(:,:), varSP(:,:), varFP(:,:)
    real(real128), allocatable :: results(:,:,:) ! Degree, Test, Repetition

    integer(int32) :: cellSPcnt
    integer(int32) :: cellFPcnt

    integer(int32) :: cellIdx, fpIdx
    integer(int32) :: spFirst, spLast
    integer(int32) :: nTest

    integer(int64) :: t_start, t_end, iter

    integer(int32) :: rep, nReps
    integer(int32) :: minDegree, maxDegree

    write(*,*) "Number of Dimensions"
    read(*,*) nDims
    write(*,*) "Number of Iterations"
    read(*,*) nIter
    write(*,*) "Number of Cells"
    read(*,*) nCells
    write(*,*) "Min Polynomial Degree"
    read(*,*) minDegree
    write(*,*) "Max Polynomial Degree"
    read(*,*) maxDegree
    write(*,*) "Number of test repetitions"
    read(*,*) nReps

    allocate(results(minDegree:maxDegree, 1:4, 1:nReps))
    results = 0.0_real64

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
            do iter=1, nIter
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    spLast  = cellIdx*cellSPcnt

                    do fpIdx=1, cellFPcnt
                        varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(sp2fp(:, fpIdx), varSP(spFirst:spLast, :))
                    end do
                end do
            enddo
            call system_clock(t_end)
            results(interpDegree, 1, rep) = real(t_end-t_start, real128)/(1000000.0_real128*real(nIter, real128))
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
            do iter=1, nIter
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    spLast  = cellIdx*cellSPcnt

                    do fpIdx=1, cellFPcnt
                        varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(transpose(varSP(spFirst:spLast, :)), sp2fp(:, fpIdx))
                    end do
                end do
            enddo
            call system_clock(t_end)
            results(interpDegree, 2, rep) = real(t_end-t_start, real128)/(1000000.0_real128*real(nIter, real128))
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
            do iter=1, nIter
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    spLast  = cellIdx*cellSPcnt

                    do fpIdx=1, cellFPcnt
                        varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(varSP(:, spFirst:spLast), sp2fp(:, fpIdx))
                    end do
                end do
            enddo
            call system_clock(t_end)
            results(interpDegree, 3, rep) = real(t_end-t_start, real128)/(1000000.0_real128*real(nIter, real128))
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
            do iter=1, nIter
                do cellIdx=1, nCells
                    spFirst = (cellIdx-1)*cellSPcnt+1
                    spLast  = cellIdx*cellSPcnt

                    do fpIdx=1, cellFPcnt
                        varFP((cellIdx-1)*cellFPcnt+fpIdx, :) = matmul(sp2fp(:, fpIdx), transpose(varSP(:, spFirst:spLast)))
                    end do
                end do
            enddo
            call system_clock(t_end)
            results(interpDegree, 4, rep) = real(t_end-t_start, real128)/(1000000.0_real128*real(nIter, real128))
            write(*,"(A, ES10.4E1, A)") "        Test 4: ", results(interpDegree, 4, rep), " ms avg per iteration"

            deallocate(sp2fp)
            deallocate(varSP)
            deallocate(varFP)
        enddo
    enddo

    write(*,*)
    do rep=1, nReps
        write(*,"(A, I2, A)") "Rep ", rep, " results"
        do interpDegree=minDegree, maxDegree
            write(*,"(A, I2, A, 4ES10.4E1)") "  Degree ", interpDegree, ": ", results(interpDegree, :, rep)
        end do
    end do
    write(*,*)
    write(*,*) "Rep average"
    do interpDegree=minDegree, maxDegree
        write(*,"(A, I2, A, 4ES10.4E1)") "    Degree ", interpDegree, ": ", sum(results(interpDegree, 1, :))/real(nReps, real128), &
                                                                            sum(results(interpDegree, 2, :))/real(nReps, real128), &
                                                                            sum(results(interpDegree, 3, :))/real(nReps, real128), &
                                                                            sum(results(interpDegree, 4, :))/real(nReps, real128)
    end do
    write(*,*)
    write(*,*) "Rep Std-Dev"
    do interpDegree=minDegree, maxDegree
        write(*,"(A, I2, A, 4ES10.4E1)") "    Degree ", interpDegree, ": ", sqrt(sum(results(interpDegree, 1, :)**2)-sum(results(interpDegree, 1, :))**2/real(nReps, real128))/real(nReps-1, real128), &
                                                                            sqrt(sum(results(interpDegree, 2, :)**2)-sum(results(interpDegree, 2, :))**2/real(nReps, real128))/real(nReps-1, real128), &
                                                                            sqrt(sum(results(interpDegree, 3, :)**2)-sum(results(interpDegree, 3, :))**2/real(nReps, real128))/real(nReps-1, real128), &
                                                                            sqrt(sum(results(interpDegree, 4, :)**2)-sum(results(interpDegree, 4, :))**2/real(nReps, real128))/real(nReps-1, real128)
    end do
    ! stddev = sqrt((sum(x**2)-sum(x)**2/size)/(size-1))

end program InterpBench
