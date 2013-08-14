module shaker
  use structure_c
  use geo_data_c
  use rng
  implicit none
  save

  double precision, parameter :: pi = 3.141592653589793d0

  type ingredient
    type(structure), pointer :: s => null()
    type(geo_data), pointer :: g => null()
    type(ingredient), pointer :: next => null()
  end type

  !~public
  integer :: state = 0, totalnatm = 0, totalofails = 0
  integer :: runs = 10, tfail = 10, ofail = 10
  double precision :: tthresh = .7d0, othresh = 1.1d0, gd = 4.d0, totalenergy = 0.d0
  character(len=72) :: indir = '', outdir = 'out/', cmdstr = 'command string', oname = 'struc'
  character(len=3) :: ofmt = 'xyz', ifmt = 'xyz'
  logical :: rng_seeded = .false., struc_r_err = .false.
  !~private
  integer :: rotcnt = 0, naxs = 0
  integer, dimension(:,:), allocatable :: taxspv
  double precision :: rmin = 0.d0, rmax = 0.d0, gdtmp
  double precision, dimension(3) :: vmove = (/0.d0,0.d0,0.d0/)
  character(len=71) :: inpath
  logical :: glue_option = .false., ch_ori = .false.

  type(ingredient), target :: head

contains

  subroutine shaker_read_recipe(line,rerr)
    integer :: i, rerr
    character(*) :: line
    character(len=1) :: dummy

    select case(state)
      case(0)
        if (line(1:3) .eq. 'def') then
          read(line,'(A4,A71)',iostat=rerr,err=666,end=666) dummy, inpath
          gdtmp = gd
          state = 2
        else if (line(1:3) .eq. 'run') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, runs
          if (runs .ge. 1000000) then
            write(6,'(A)') 'Maximal number of runs is limited to 999999 currently'
            runs = 999999
          end if
        else if (line(1:2) .eq. 'id') then
          read(line,'(A3,A72)',iostat=rerr,err=666,end=666) dummy, indir
        else if (line(1:2) .eq. 'od') then
          read(line,'(A3,A72)',iostat=rerr,err=666,end=666) dummy, outdir
        else if (line(1:2) .eq. 'of') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, ofmt
        else if (line(1:2) .eq. 'if') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, ifmt
        else if (line(1:3) .eq. 'cmd') then
          read(line,'(A4,A71)',iostat=rerr,err=666,end=666) dummy, cmdstr
        else if (line(1:2) .eq. 'on') then
          read(line,'(A3,A72)',iostat=rerr,err=666,end=666) dummy, oname
        else if (line(1:2) .eq. 'gd') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, gd
        else if (line(1:3) .eq. 'tth') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, tthresh
        else if (line(1:3) .eq. 'pth') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, othresh
        else if (line(1:3) .eq. 'tfl') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, tfail
          if (tfail .ge. 100000) then
            write(6,'(A)') 'tfail is limited to 99999 currently'
            tfail = 99999
          end if
        else if (line(1:3) .eq. 'pfl') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, ofail
          if (ofail .ge. 100000) then
            write(6,'(A)') 'pfail is limited to 99999 currently'
            ofail = 99999
          end if
        else if (line(1:3) .eq. 'srn') then
          rng_seeded = .true.
        end if
      case(1)
        if (line(1:3) .eq. 'def') then
          read(line,'(A4,A71)',iostat=rerr,err=666,end=666) dummy, inpath
          gdtmp = gd
          state = 2
        end if
      case(2)
        if (line(1:2) .eq. 'mv') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, vmove(:)
        else if (line(1:2) .eq. 'os') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, rmin, rmax
          ch_ori = .true.
        else if (line(1:3) .eq. 'gts') then
          read(line,*,iostat=rerr) dummy, gdtmp
          glue_option = .true.
          rerr = 0
        else if (line(1:2) .eq. 'rb') then
          read(line,*,iostat=rerr,err=666,end=666) dummy, naxs
          if (naxs .gt. 0) then
            allocate(taxspv(2,naxs))
            state = 3
          else
            naxs = 0
          end if
        else if (line(1:3) .eq. 'end') then
          if (.not. allocated(taxspv)) then
            allocate(taxspv(2,naxs))
          end if
          call shaker_add_ing()
          deallocate(taxspv)
          ch_ori = .false.
          glue_option = .false.
          naxs = 0
          vmove(:) = 0.d0
          rmin = 0.d0
          rmax = 0.d0
          state = 1
        end if
      case(3)
        rotcnt = rotcnt + 1
        read(line,*,iostat=rerr,err=666,end=666) dummy, taxspv(:,rotcnt)
        if (rotcnt .eq. naxs) then
          rotcnt = 0
          state = 2
        end if
      case default
        write(*,*) 'State unkown / fail'
    end select
666 return
  end subroutine

  subroutine shaker_add_ing()
    type(structure), pointer :: struc => null()
    type(geo_data), pointer :: geo => null()
    type(ingredient), pointer :: current => null(), putpos => null()
    integer :: cnt
    
    cnt = 0
    putpos => head
    
    do while (associated(putpos%next))
      cnt = cnt + 1
      putpos => putpos%next
    end do
    
    allocate(current,struc,geo)
    
    call struc_read_file(struc,trim(join_fpath(indir,inpath)),ifmt,struc_r_err)
    call struc_set(struc,t=cnt+1)
    
    call geo_def(geo,vmove,ch_ori,rmin,rmax,glue_option,gdtmp,naxs,taxspv,struc%bonds)
    
    current%s => struc
    current%g => geo
    
    totalnatm = totalnatm + struc_get_natoms(current%s)
    totalenergy = totalenergy + struc_get_energy(current%s)
    
    putpos%next => current
    
    nullify(struc,geo,putpos)
  end subroutine

  subroutine shaker_power_up()
    call rng_init_seed(rng_seeded)
  end subroutine

  subroutine shaker_shake(filenum,retcode)
    integer :: filenum, retcode
    type(ingredient), pointer :: current => null()
    integer :: i, failcnt, natmtmp, naxs, loopcnt
    character(len=20) :: failcntc
    logical :: fail1, fail2
    type(structure) :: multimer
    character(len=80) :: fpath
    character(len=2), dimension(:), allocatable :: atmtype
    double precision, dimension(:,:), allocatable :: coords
    
    failcnt = 0
    current => head%next
    
    do while (associated(current))
      
      naxs = size(current%g%axes,2)
      
      if (naxs .gt. 0) then
        
        do i=1,naxs
          call struc_rot_bond(current%s,rng_get_double(0.d0,2*pi),current%g%axes(:,i),current%g%rotors(i)%points(:))
        end do
        
        call struc_test(current%s,tthresh,fail1)
        
        if (fail1) then
          
          failcnt = failcnt + 1
          
          if (failcnt .lt. tfail) then
            cycle
          else
            write(failcntc,'(I6)') failcnt
            write(6,'(A,I2)') 'Rotating bonds failed '//trim(adjustl(failcntc))//' times in a row for ingredient ',&
              &struc_get_tag(current%s)
            retcode = 1
            return
          end if
          
        else
          
          failcnt = 0
          
        end if
        
      end if
      
      current => current%next
      
    end do
    
    do
      current => head%next
      
      do while (associated(current))
        
        call struc_ctrans(current%s)
        call struc_move(current%s,current%g%vt)
        
        if (current%g%ch_ori) then
          call struc_move(current%s,rng_get_vec_ss(current%g%rmin,current%g%rmax))
          call struc_rot(current%s,rng_get_rotmat())
        end if
        
        current => current%next
      end do
      
      call shaker_test_glue(fail2)
      call shaker_test_shake(fail1)
      
      if (fail1 .or. fail2) then
        
        failcnt = failcnt + 1
        totalofails = totalofails + 1
        
        if (failcnt .lt. ofail) then
          cycle
        else
          write(failcntc,'(I6)') failcnt
          write(6,'(A)') 'Positioning failed '//trim(adjustl(failcntc))//' times in a row;'//&
            &' the setup is probably too restrictive'
          retcode = 2
          return
        end if
        
      else
        exit
      end if
      
    end do
    
    write(fpath,'(A,I6.6,A)') trim(join_fpath(outdir,oname)),filenum,'.'//ofmt
    
    allocate(atmtype(totalnatm), coords(3,totalnatm))
    current => head%next
    
    natmtmp = 0
    loopcnt = 0
    
    do while (associated(current))
      
      do i = 1, struc_get_natoms(current%s)
        atmtype(i+natmtmp) = current%s%atmtype(i)
        coords(:,i+natmtmp) = current%s%coords(:,i)
        loopcnt = loopcnt + 1
      end do
      
      natmtmp = loopcnt
      current => current%next
      
    end do
    
    
    call struc_alloc(multimer,totalnatm)
    call struc_set(multimer,atmtype,coords)
    
    deallocate(atmtype,coords)
    
    call struc_write_file(multimer,trim(fpath),ofmt,fail1,trim(cmdstr))
    
    if (fail1) then
      retcode = 3
      return
    else
      write(6,'(A)') 'Writing '//trim(fpath)
      retcode = 0
    end if
  end subroutine

  subroutine shaker_test_shake(rslt)
    type(ingredient), pointer :: ing1 => null(), ing2 => null()
    integer :: i, j
    logical :: rslt
    
    rslt = .false.
    ing1 => head
    
    do
      
      if (associated(ing1%next)) then
        
        ing1 => ing1%next
        if (associated(ing1%next)) then
          ing2 => ing1%next
        else
          exit
        end if
        
        do
          
          rslt = shaker_chk_dist(ing1%s,ing2%s,othresh)
          if (rslt) return
          
          if (associated(ing2%next)) then
            ing2 => ing2%next
          else
            exit
          end if
          
        end do
        
      else
        exit
      end if
      
    end do

    nullify(ing1,ing2)
  end subroutine
  
  subroutine shaker_test_glue(rslt)
    type(ingredient), pointer :: ing1 => null(), ing2 => null()
    integer :: i, j
    logical :: rslt, tmprslt
    double precision :: gluedist
    
    rslt = .false.
    
    ing1 => head
    
    do
      
      if (associated(ing1%next)) then
        
        ing1 => ing1%next
        if (.not. ing1%g%glue) cycle
        gluedist = ing1%g%gd
        
        ing2 => head
        
        do
          
          if (associated(ing2%next)) then
            ing2 => ing2%next
          else
            exit
          end if
          
          if (struc_get_tag(ing2%s) .eq. struc_get_tag(ing1%s)) cycle
          
          tmprslt = shaker_chk_dist(ing1%s,ing2%s,gluedist)
          
          if (tmprslt) exit
          
        end do
        
        if (.not. tmprslt) then
          rslt = .true.
          return
        end if
        
      else
        exit
      end if
      
    end do
  end subroutine
  
  function shaker_chk_dist(struc1,struc2,dist) result(rslt)
    type(structure) :: struc1, struc2
    double precision :: dist
    logical :: rslt
    integer :: i, j
    
    rslt = .false.
    
    do i = 1, struc_get_natoms(struc1)
      do j = 1, struc_get_natoms(struc2)
        if (get_dist(struc1%coords(:,i),struc2%coords(:,j)) .le. dist) then
          rslt = .true.
          return
        end if
      end do
    end do
  end function

  subroutine shaker_clr()
    type(ingredient), pointer :: current => null(), next => null()
    
    current => head%next
    nullify(head%next)
    
    do while (associated(current))
      next => current%next
      call struc_clr(current%s)
      call geo_clr(current%g)
      deallocate(current)
      nullify(current)
      current => next
    end do
    
    totalnatm = 0
    totalenergy = 0.d0
  end subroutine

  subroutine shaker_write_job()
    type(ingredient), pointer :: current => null()
    character(len=6) :: intch
    character(len=8) :: rngseedstat
    
    write(intch,'(I6)') runs
    
    if (rng_seeded) then
      rngseedstat = 'seeded'
    else
      rngseedstat = 'unseeded'
    end if
    write(6,'(A)')
    write(6,'(A)') ' Input directory: '//trim(adjustl(indir))
    write(6,'(A)') ' Output directory: '//trim(adjustl(outdir))
    write(6,'(A)') ' File formats: '//ifmt//' => '//ofmt
    write(6,'(A)') ' Output filenames: '//trim(oname)//'####'
    write(6,'(A)')
    write(6,'(A)') ' Number of runs: '//adjustl(intch)
    write(6,'(A)') ' Random number generator is '//trim(rngseedstat)
    write(6,'(A)')
    write(6,'(A,F5.2)') ' Min intramol distances: ', tthresh
    write(6,'(A,F5.2)') ' Min intermol distances: ', othresh
    write(6,'(A)')

    current => head%next

    do while (associated(current))
      write(6,'(A)')
      write(6,'(A,I2)') 'INGREDIENT ', struc_get_tag(current%s)
      write(6,'(A)')
      call struc_print_bs(current%s,6)
      write(6,'(A)')
      call geo_write_all(current%g,6)
      current => current%next
    end do
    
    write(intch,'(I5)') totalnatm
    
    write(6,'(A)')
    write(6,'(A)') 'TOTAL'
    write(6,'(A)')
    write(6,'(A)') ' Number of atoms: '//adjustl(intch)
    write(6,'(A,F15.9)') ' Energy: ', totalenergy
    write(6,'(A)')
  end subroutine
  
  function join_fpath(dir,fname) result(string)
    character(*) :: dir, fname
    character(len=99) :: string
    integer :: lpos
    
    if (dir .ne. '') then
      lpos = len_trim(dir)
      
      if (dir(lpos:lpos) .ne. '/') then
        string = trim(adjustl(dir))//'/'//trim(adjustl(fname))
      else
        string = trim(adjustl(dir))//trim(adjustl(fname))
      end if
    else
      string = trim(adjustl(fname))
    end if
  end function
end module

