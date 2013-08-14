program main
  use shaker
  implicit none
  integer :: i, linecnt, readerr = 0, failcode = 0, ierr = 0
  character(len=60) :: current_line, jobfile
  character(len=1) :: userinput

  write(6,'(A)') '***************************************************************************'
  write(6,'(A)') '**                                                                       **'
  write(6,'(A)') '**                              SHAKER v2.1                              **'
  write(6,'(A)') '**                                                                       **'
  write(6,'(A)') '**                     enhanced molecule mixing tool                     **'
  write(6,'(A)') '**                                                                       **'
  write(6,'(A)') '***************************************************************************'
  write(6,'(A)')
  write(6,'(A)')

  call getarg(1,jobfile)

  open(110,file='config.txt',status='old',iostat=ierr)
  if (ierr .eq. 0) then
    linecnt = 0
    do
      linecnt = linecnt + 1
      read(110,'(A60)',iostat=ierr) current_line
      if (ierr .eq. 0 .and. state .ne. 9) then
        call shaker_read_recipe(current_line,readerr)
        if (readerr .ne. 0) then
          write(6,'(A,I3.1,A)') 'Error reading line',linecnt,' of config.txt; aborting read'
          exit
        end if
      else
        write(6,'(A)') 'config.txt processed'
        exit
      end if
      state = 0
    end do
    close(110)
  else
    write(6,'(A)') 'Configuration file not found (config.txt)'
  end if

  open(111,file=trim(jobfile),status='old',iostat=ierr)
  if (ierr .eq. 0) then
    linecnt = 0
    do
      linecnt = linecnt + 1
      read(111,'(A60)',iostat=ierr) current_line
      if (ierr .eq. 0) then
        call shaker_read_recipe(current_line,readerr)
        if (readerr .ne. 0) then
          write(6,'(A,I3.1,A)') 'Error reading line',linecnt,' of "'//trim(jobfile)//'"; aborting read'
          exit
        end if
        if (struc_r_err) then
          state = 9
          exit
        end if
      else
        write(6,'(A)') '"'//trim(jobfile)//'" processed'
        exit
      end if
    end do
    close(111)
  else
    write(6,'(A)') 'Jobfile not found ('//trim(jobfile)//')'
  end if

  if (state .ne. 1) then
    write(6,'(A,I1)') 'STOP: irregular termination of config procedure; state = ', state
    stop
  end if

  call shaker_write_job()
  
  write(6,'(A)')
  write(6,'(A)') 'Do you wish to proceed? (*/n)'
  read(5,*) userinput
  write(6,'(A)')
  
  if (userinput .eq. 'n' .or. userinput .eq. 'N') then
    call shaker_clr()
    stop
  end if
  
  call shaker_power_up()
  
  do i=1,runs
    call shaker_shake(i,failcode)
    if (failcode .ne. 0) exit
  end do

  call shaker_clr()
  
  if (failcode .ne. 0) then
    write(6,'(A,I1)') 'STOP: unable to shake; reason = ', failcode
  else
    write(6,'(A)')
    write(6,'(A,F5.1,A)') ' Building efficiency: ', 1.d2 - totalofails*1.d2/dble(runs + totalofails),&
      &'% of created structures were viable'
    write(6,'(A)')
    write(6,'(A)')
    write(6,'(A)') '**                      ... Shaken, not stirred ...                      **'
  end if
end program
