module geo_data_c

  type intray
    integer, dimension(:), pointer :: points => null()
  end type

  type geo_data
    logical :: ch_ori, glue
    double precision :: rmin, rmax, gd
    double precision, dimension(3) :: vt
    integer, dimension(:,:), pointer :: axes => null()
    type(intray), dimension(:), pointer :: rotors => null()
  end type
  
contains

  subroutine geo_def(this,vt,ch_ori,rmin,rmax,glue,gd,naxs,taxspv,bm)
    implicit none
    type(geo_data) :: this
    integer :: naxs
    logical :: ch_ori, glue
    double precision :: rmin, rmax, gd
    double precision, dimension(3) :: vt
    integer, dimension(2,naxs) :: taxspv
    integer, dimension(:,:) :: bm
    
    integer :: i
    
    this%vt = vt
    this%ch_ori = ch_ori
    this%glue = glue
    this%rmin = rmin
    this%rmax = rmax
    this%gd = gd
    
    allocate(this%axes(2,naxs),this%rotors(naxs))
    do i=1,naxs
      this%axes(:,i) = taxspv(:,i)
    end do
    call geo_detect_rotors(this,naxs,bm)
    
  end subroutine
  
  subroutine geo_detect_rotors(this,naxs,bm_in)
    type(geo_data) :: this
    integer, dimension(:,:) :: bm_in
    integer, dimension(:,:), allocatable :: bm
    integer, dimension(:), allocatable :: rotor
    integer :: i, j, n, naxs, nrot
    
    n = size(bm_in,1)
    allocate(bm(n,n),rotor(n))
    
    do j = 1, naxs
      bm = bm_in
      rotor = 0
      nrot = 0
      
      call geo_detect_rotor_r(bm,this%axes(2,j),this%axes(1,j),rotor,nrot)
      
      allocate(this%rotors(j)%points(nrot))
      
      do i = 1, nrot
        this%rotors(j)%points(i) = rotor(i)
      end do
    end do
    
    deallocate(bm,rotor)
  end subroutine
  
  recursive subroutine geo_detect_rotor_r(bm,c,l,rotor,nrot)
    integer, dimension(:,:) :: bm
    integer, dimension(:) :: rotor
    integer :: c, l, n, i, nrot
    
    n = size(bm,1)
    
    bm(l,:) = 0
    
    do i = 1, n
      if (bm(i,c) .eq. 1) then
        nrot = nrot + 1
        rotor(nrot) = i
        call geo_detect_rotor_r(bm,i,c,rotor,nrot)
      end if
    end do
    
  end subroutine
  
  subroutine geo_write_all(this,unt)
    type(geo_data) :: this
    integer :: unt, naxs, i
    character(len=20) :: wfmt
    
    write(unt,'(A,3F7.3)') 'Position of geometric center: ', this%vt(:)
    write(unt,'(A)')
    
    if (this%ch_ori) then
      write(unt,'(A)') 'Random orientation in a spherical shell'
      write(unt,'(A,F7.3,A,F7.3)') 'R_min = ', this%rmin, ' ,    R_max = ', this%rmax
      write(unt,'(A)')
    end if
    
    if (this%glue) then
      write(unt,'(A,F6.2)') 'Maximum distance from another ingredient: ', this%gd
      write(unt,'(A)')
    end if
    
    naxs = size(this%axes,2)
    
    if (naxs .gt. 0) then
      write(unt,'(A)') 'All torsion axes with corresponding rotor'
      do i = 1, naxs
        write(wfmt,'(A,I3.3,A)') '(2I3,A,', size(this%rotors(i)%points), 'I3)'
        write(unt,trim(wfmt)) this%axes(:,i), ' >',this%rotors(i)%points(:)
      end do
      write(unt,'(A)')
    end if
  end subroutine
  
  subroutine geo_clr(this)
    implicit none
    type(geo_data) :: this
    integer :: i
    
    do i=1,size(this%axes,2)
      deallocate(this%rotors(i)%points)
      nullify(this%rotors(i)%points)
    end do
    deallocate(this%axes,this%rotors)
    nullify(this%axes,this%rotors)
    
    this%vt = 0.d0
    this%ch_ori = .false.
    this%glue = .false.
    this%rmin = 0.d0
    this%rmax = 0.d0
    this%gd = 0.d0
    
  end subroutine

end module
