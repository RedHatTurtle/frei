prototype module Output
{
  proc iterOutput(nIter)
  {
    integer(isp), intent(in) :: nIter
    character(len=20) :: stringIter, dirName

    write(stringIter,'(i)') nIter
    dirName  = 'iter' // trim(adjustL(stringIter))

    call execute_command_line('mkdir -p '//trim(dirName) )
    call stdOut(dirName)
  }

  proc timeOutput(curTime)
  {
    real(rdp), intent(in) :: curTime
    character(len=20) :: stringTime, dirName

    ! Define the string format to write the current time on the file/dir name string
    ! Build string with the file/directory name
    write(stringTime,'(f0.4)') curTime
    dirName = 'time' // trim(stringTime)

    call execute_command_line('mkdir -p '//trim(dirName) )
    call stdOut(dirName)
  }

  proc output_solution_gnuplot_1d(in dirName : string)
  {
    fileName : string;

    fileName = trim(dirName) // '/stdOut.out'
    stdOutfile = 101; open( unit=stdOutFile, file=trim(fileName), status='replace' )

    for cell in 0..nCells
      for dof in 0..nSPs
        write(stdOutFile,'(i7,i3,7(es12.3E3))') cID, dofID, mesh%cellList(cID)%DOF(dofID)%x, mesh%cellList(cID)%DOF(dofID)%u(:), mesh%cellList(cID)%DOF(dofID)%p(), mesh%cellList(cID)%DOF(dofID)%T(), mesh%cellList(cID)%DOF(dofID)%u(2)/mesh%cellList(cID)%DOF(dofID)%u(1)
      enddo
    enddo

    close( stdOutFile )
  }
}
