module transform_tools  

  !******************************************************************
  ! Contains procedures transforming grid<->wave and for getting
  ! energy and enstrophy conserving jacobian via use of the transform.  
  ! Uses staggered grid to remove aliasing via method of Orszag 71.
  !
  ! Routines:  Init_transform, Spec2grid, Grid2spec, Jacob, 
  !            ir_prod, ir_pwr
  !
  ! Dependencies: fft_mod, par_tools, op_rules, io_tools
  !
  ! Algorithm by Geoff Vallis (c1991) and refitted for F90 by Shafer Smith
  ! (c1998).  Redone again in 2000 by SS. Parallelized by Guillaume Lapeyre
  ! to use mpi in 2001. Totally rewritten and thoroughly tested in 2005 by SS.
  !
  ! Variable names:
  !
  ! kmax :  largest wavenumber allowed in each horizontal direction
  ! ngrid:  size of square grid-space array: nx = ny = 2*(kmax+1)
  ! kx...:  proper values of x-wavenumber
  ! ky...:  proper values of y-wavenumber
  ! ..vec:  vector (1d array) of values
  ! nkx  :  extent of kx-dimension of LOCAL spectral fields 
  ! ny   :  extent of y-dimension (dim 1) of LOCAL grid field
  !    --- NOTE THAT nkx = ny IS ENFORCED BY FFT_FFTW ---
  ! nky  :  extent of ky-dimension of spectral fields
  !******************************************************************

  implicit none
  private
  save

  real,parameter                          :: pi = 3.14159265358979324
  integer,parameter                       :: forward = 1, backward = -1
  complex,parameter                       :: i = (0.,1.)
  complex,dimension(:,:),allocatable      :: one_plus_ialpha,one_minus_ialpha
  complex,dimension(:,:),allocatable      :: one_plus_icalpha,one_minus_icalpha
  real,dimension(:,:),allocatable         :: sgn
  integer,dimension(:),allocatable        :: kxvec,kyvec
  integer                                 :: kmax, nz, ngrid, nky      
  integer                                 :: kx_start, kx_end, nkx, y_start, y_end, ny
  integer,dimension(:),allocatable        :: cpu
  logical                                 :: first_lhp_pe, first_rhp_pe

  public :: Init_transform, Spec2grid, Grid2spec, Jacob, ir_prod, ir_pwr
  
  interface Spec2grid
     module procedure Spec2grid_cc2, Spec2grid_cc3
     module procedure Spec2grid_cr2, Spec2grid_cr3
  end interface
  interface Grid2spec
     module procedure  Grid2spec_cc3, Grid2spec_cc2
     module procedure  Grid2spec_rc2, Grid2spec_rc3
  end interface
  interface ir_prod
     module procedure ir_prod2, ir_prod3
  end interface
  interface ir_pwr
     module procedure ir_pwr2, ir_pwr3
  end interface
  interface jacob
     module procedure jacob2, jacob3
  end interface

contains

  !*************************************************************************
  
  subroutine Init_transform(kmaxi, nzi, kx_starto, kx_endo, nkxo, y_starto, y_endo, nyo)

    ! Set up index arrays for transform routines and initialize fft

    use fft_mod,   only: Init_fft
    use par_tools, only: processor_id, num_processors, par_allgather, par_sync
    use io_tools,  only: message

    integer,intent(in)                 :: kmaxi, nzi
    integer,intent(out)                :: kx_starto, kx_endo, nkxo, y_starto, y_endo, nyo
    integer                            :: kx, ky, n, m, p, first_kx, last_kx
    integer,dimension(:),allocatable   :: first_and_last_kx, indeces, num_ind
    complex,dimension(:,:),allocatable :: alpha

    kmax = kmaxi
    nz = nzi
    ngrid = 2*(kmax+1)
    nky = kmax+1

    call par_sync

    ! Sends kmax and nz, and returns kx_start, nkx, y_start, ny
    call Init_fft(kmax,nz, kx_start, kx_end, nkx, y_start, y_end, ny)      

    ! Set dummy args for return to caller
    kx_starto = kx_start;  kx_endo = kx_end; nkxo = nkx;  
    y_starto = y_start;    y_endo = y_end;   nyo = ny

    ! Set up array 'cpu' that maps kx to processor_ids.  This requires info obtained by fftw init 
    ! above.  cpu will be stored identically on every processor, so need to broadcast each pe's 
    ! info to all other processors.  

    allocate(cpu(-kmax-1:kmax))      
    allocate(first_and_last_kx(2), indeces(2*num_processors), num_ind(num_processors))

    first_and_last_kx(1) = kx_start               ! kx_start starts at -kmax-1
    first_and_last_kx(2) = kx_end

    call par_allgather(first_and_last_kx,indeces) ! First arg on each pe is concatenated 
                                                  ! into second argument, existing on each pe
    do p = 0 , num_processors-1                   ! Now unpack all_indeces and make cpu
       first_kx = indeces(2*p+1)
       last_kx = indeces(2*p+2)
       cpu(first_kx:last_kx) = p
    enddo

    ! print*, 'mype', processor_id, 'kx_start', kx_start
    ! if (processor_id==0) print*, 'cpu',cpu

    ! Set up some arrays for packing and unpacking fields for transform
    
    allocate(kxvec(nkx),kyvec(nky));      kxvec=0; kyvec=0

    kxvec  = (/ (kx, kx = kx_start, kx_end) /)  ! True wavenumber vectors     
    kyvec  = (/ (ky, ky = 0,        kmax)   /)

    ! flags to know whether we are on first LHP or first RHP processor, for which first 
    ! column of ky values is simply coppied to lower half plane.  
    if (kx_start == -kmax-1)  then
       first_lhp_pe = .true.
    else
       first_lhp_pe = .false.
    endif
    if (kx_start == 0)  then
       first_rhp_pe = .true.
    else
       first_rhp_pe = .false.
    endif

    ! Get packing/gridshift factor for Orszag transform.  F*exp(i*pi*K) is equivalent
    ! to F on a staggered grid....  See (nonexistant) notes.

    allocate(alpha(nkx,nky))
    allocate(one_plus_ialpha(nkx,nky),one_minus_ialpha(nkx,nky))
    allocate(one_plus_icalpha(nkx,nky),one_minus_icalpha(nkx,nky))

    alpha =  cexp(i*pi*(spread(kxvec,2,nky)+spread(kyvec,1,nkx))/ngrid) ! = e^(i*pi*(kx+ky)/ngrid)
    one_plus_ialpha   = 1. + i*alpha
    one_minus_ialpha  = 1. - i*alpha
    one_plus_icalpha  = 0.25*(1. + i*conjg(alpha))
    one_minus_icalpha = 0.25*(1. - i*conjg(alpha))

    deallocate(alpha)

    ! sgn multiplies phys space field to compensate for having 0 frequency in 
    ! middle of spectral domain.  Old trick in fft-ology.
    ! sgn is transposed here because it is a grid array, assuming FFTW_TRANSPOSE.
    allocate(sgn(ny,ngrid));  sgn = 0;
    do n=1,ny
       sgn(n,:) = (/ ((-1)**(m+n), m=1,ngrid) /)
    enddo

    call par_sync
    call Message('Transform tools initialized')

  end subroutine Init_transform

  !*************************************************************************

  subroutine Jacob2(fk,gk,jack) 

    ! This routine returns the jacobian of arrarys f() and g().

    use op_rules,  only: operator(*)

    complex,dimension(nkx,0:kmax),intent(in)  :: fk,gk
    complex,dimension(nkx,0:kmax),intent(out) :: jack
    complex,dimension(ny,ngrid)               :: f,g,temp

    call Spec2grid_cc2(i*spread(kxvec,2,nky)*fk,f)
    call Spec2grid_cc2(i*spread(kyvec,1,nkx)*gk,g)
    temp = ir_prod2(f,g)
    
    call Spec2grid_cc2(i*spread(kyvec,1,nkx)*fk,f)
    call Spec2grid_cc2(i*spread(kxvec,2,nky)*gk,g)
    temp = temp - ir_prod2(f,g)
    
    call Grid2spec_cc2(temp, jack)

  end subroutine Jacob2

  !*************************************************************************

  subroutine Jacob3(fk,gk,jack)

    ! This routine returns the jacobian of arrarys f() and g().

    use op_rules,  only: operator(*)

    complex,dimension(nz,nkx,0:kmax),intent(in)  :: fk,gk
    complex,dimension(nz,nkx,0:kmax),intent(out) :: jack
    complex,dimension(nz,ny,ngrid)               :: f,g,temp

    call Spec2grid_cc3(i*spread(kxvec,2,nky)*fk, f)
    call Spec2grid_cc3(i*spread(kyvec,1,nkx)*gk, g)
    temp = ir_prod3(f,g)
    
    call Spec2grid_cc3(i*spread(kyvec,1,nkx)*fk, f)
    call Spec2grid_cc3(i*spread(kxvec,2,nky)*gk, g)
    temp = temp - ir_prod3(f,g)
    
    call Grid2spec_cc3(temp,jack)

  end subroutine Jacob3

  !*******************************************************************

  subroutine Spec2grid_cc2(wavefield_uhp,gridfield)

    ! Spectral to grid transform, passing straight grid field
    ! to real part of complex output variable (gridfield), and same
    ! physical field *on staggered grid* to *imaginary* component
    ! of result variable  (for use in Orszag transform method).
    ! Input:  upper half plane of spectral field (wavefield_uhp)
    ! Output: resulting grid field (gridfield)

    use fft_mod,   only: fft

    complex,dimension(nkx,0:kmax),intent(in)  :: wavefield_uhp
    complex,dimension(ny,ngrid ),intent(out)  :: gridfield
    complex,dimension(nkx,-kmax-1:kmax)       :: wavefield
    complex,dimension(nkx,-kmax:0)            :: wavefield_lhp

    ! Lower half plane requires swaps between pairs of processors that are symmetric in kx, 
    ! but first kx (column of ky values) on each pe is special.  pes containing kx=-kmax-1 
    ! (dummy column) and kx=0 trade this column directly to lower half plane.  All others 
    ! need to swap it, but not to same processor that the rest of the kx columns on current 
    ! pe swap with.  

    ! Copy packed (w/ staggered and straight grids) upper half plane into full field
    wavefield = 0.
    wavefield(:,0:kmax) = wavefield_uhp*one_plus_ialpha
    wavefield_lhp       = conj_symm2(wavefield_uhp*one_minus_ialpha)

    ! Fill LHP 
    wavefield(:,-kmax:-1) = wavefield_lhp(:,-kmax:-1)

    ! if on left LHP, do  modify ky=0, otherwise don't
    where (kxvec<=0) wavefield(:,0) = wavefield_lhp(:,0)

    if (first_lhp_pe)  wavefield(1,:) = 0. ! dummy column

    gridfield = fft(wavefield,backward)
    gridfield = sgn*gridfield

  end subroutine Spec2grid_cc2

  !*******************************************************************

  subroutine Spec2grid_cc3(wavefield_uhp,gridfield)

    ! Same as Spec2grid_cc2, but for 3d fields

    use fft_mod,   only: fft
    use op_rules,  only: operator(*)

    complex,dimension(nz,nkx,0:kmax),intent(in)  :: wavefield_uhp
    complex,dimension(nz,ny,ngrid),intent(out)   :: gridfield
    complex,dimension(nz,nkx,-kmax-1:kmax)       :: wavefield
    complex,dimension(nz,nkx,-kmax:0)            :: wavefield_lhp
    integer                                      :: n

    wavefield = 0.
    wavefield(:,:,0:kmax) = wavefield_uhp*one_plus_ialpha 
    wavefield_lhp = conj_symm3(wavefield_uhp*one_minus_ialpha)

    ! Fill LHP 
    wavefield(:,:,-kmax:-1) = wavefield_lhp(:,:,-kmax:-1)
    do n=1,nz
       where (kxvec<=0) wavefield(n,:,0) = wavefield_lhp(n,:,0)
    enddo

    if (first_lhp_pe) wavefield(:,1,:) = 0. ! dummy column

    gridfield = fft(wavefield,backward)
    gridfield = sgn*gridfield

  end subroutine Spec2grid_cc3

  !*******************************************************************

  subroutine Spec2grid_cr2(wavefield_uhp,gridfield)

    ! Spectral to grid transform, returning real result (only for 
    ! diagnostic use)

    use fft_mod,   only: fft

    complex,dimension(nkx,0:kmax),intent(in)  :: wavefield_uhp
    real,   dimension(ny,ngrid),intent(out)   :: gridfield
    complex,dimension(nkx,-kmax-1:kmax)       :: wavefield
    complex,dimension(nkx,-kmax:0)            :: wavefield_lhp

    ! Lower half plane requires swaps between pairs of processors that are symmetric in kx, 
    ! but first kx (column of ky values) on each pe is special.  pes containing kx=-kmax-1 
    ! (dummy column) and kx=0 trade this column directly to lower half plane.  All others 
    ! need to swap it, but not to same processor that the rest of the kx columns on current 
    ! pe swap with.  

    ! Copy packed (w/ staggered and straight grids) upper half plane into full field
    wavefield = 0.
    wavefield(:,0:kmax) = wavefield_uhp
    wavefield_lhp       = conj_symm2(wavefield_uhp)

    ! Fill LHP 
    wavefield(:,-kmax:-1) = wavefield_lhp(:,-kmax:-1)
    where (kxvec<=0) wavefield(:,0) = wavefield_lhp(:,0)

    if (first_lhp_pe)  wavefield(1,:) = 0. ! dummy column

    wavefield = fft(wavefield,backward)
    wavefield = sgn*wavefield
    gridfield = real(wavefield)

  end subroutine Spec2grid_cr2

  !*******************************************************************

  subroutine Spec2grid_cr3(wavefield_uhp,gridfield)

    ! Same as Spec2grid_cr2, but for 3d fields

    use fft_mod,   only: fft
    use op_rules,  only: operator(*)

    complex,dimension(nz,nkx,0:kmax),intent(in)  :: wavefield_uhp
    real,   dimension(nz,ny,ngrid),intent(out)   :: gridfield
    complex,dimension(nz,nkx,-kmax-1:kmax)       :: wavefield
    complex,dimension(nz,nkx,-kmax:0)            :: wavefield_lhp
    integer                                      :: n

    wavefield = 0.
    wavefield(:,:,kmax+2:ngrid) = wavefield_uhp 
    wavefield_lhp = conj_symm3(wavefield_uhp)

    ! Fill LHP (if on left LHP, do  modify ky=0, otherwise don't)
    wavefield(:,:,-kmax:-1) = wavefield_lhp(:,:,-kmax:-1)
    do n=1,nz
       where (kxvec<=0) wavefield(n,:,0) = wavefield_lhp(n,:,0)
    enddo

    if (first_lhp_pe) wavefield(:,1,:) = 0. ! dummy column

    wavefield = fft(wavefield,backward)
    gridfield = real(wavefield)
    wavefield = sgn*wavefield

  end subroutine Spec2grid_cr3

  !*************************************************************************

  subroutine Grid2spec_cc2(gridfield,wavefield_uhp)

    ! Physical to spectral transform:  assumes complex input variable 
    ! contains field as real part and same field on staggered grid
    ! in imaginary part.  

    use fft_mod,  only: fft

    complex,dimension(ny,ngrid),intent(in)    :: gridfield
    complex,dimension(nkx,nky),intent(out)    :: wavefield_uhp
    complex,dimension(nkx,-kmax-1:kmax)       :: wavefield
    complex,dimension(nkx,nky)                :: wavefield_lhp

    wavefield     = fft(sgn*gridfield,forward)
    wavefield_uhp = wavefield(:,0:kmax)
    wavefield_lhp = conj_symm2(wavefield(:,-kmax:0))
    wavefield_uhp = wavefield_uhp*one_minus_icalpha + wavefield_lhp*one_plus_icalpha

    if (first_lhp_pe)    wavefield_uhp(1,:) = 0. ! dummy column
    where (kxvec<=0)     wavefield_uhp(:,1) = 0. ! LHP ky=0 row

  end subroutine Grid2spec_cc2

  !*************************************************************************

  subroutine Grid2spec_cc3(gridfield,wavefield_uhp)

    ! Physical to spectral transform:  assumes complex input variable 
    ! contains field as real part and same field on staggered grid
    ! in imaginary part.  

    use fft_mod,  only: fft
    use op_rules, only: operator(*)

    complex,dimension(nz,ny,ngrid),intent(in)    :: gridfield
    complex,dimension(nz,nkx,nky),intent(out)    :: wavefield_uhp
    complex,dimension(nz,nkx,-kmax-1:kmax)       :: wavefield
    complex,dimension(nz,nkx,nky)                :: wavefield_lhp
    integer                                      :: n

    wavefield     = fft(sgn*gridfield,forward)
    wavefield_uhp = wavefield(:,:,0:kmax)
    wavefield_lhp = conj_symm3(wavefield(:,:,-kmax:0))
    wavefield_uhp = wavefield_uhp*one_minus_icalpha + wavefield_lhp*one_plus_icalpha

    if (first_lhp_pe)    wavefield_uhp(:,1,:) = 0. ! dummy column
    do n=1,nz
       where (kxvec<=0)  wavefield_uhp(n,:,1) = 0. ! LHP ky=0 row
    enddo

  end subroutine Grid2spec_cc3

  !*************************************************************************

  subroutine Grid2spec_rc2(gridfield,wavefield_uhp)

    ! Physical to spectral transform - just for diagnostic use, as 
    ! product is not dealiased.

    use fft_mod,  only: fft

    real,   dimension(ny,ngrid), intent(in)   :: gridfield
    complex,dimension(nkx,nky),intent(out)    :: wavefield_uhp
    complex,dimension(nkx,-kmax-1:kmax)       :: wavefield

    wavefield   = sgn*gridfield*cmplx(1.,0.)
    wavefield   = fft(wavefield,forward)
    wavefield_uhp = wavefield(:,0:kmax)

    if (first_lhp_pe)    wavefield_uhp(1,:) = 0. ! dummy column
    where (kxvec<=0)     wavefield_uhp(:,1) = 0. ! LHP ky=0 row

  end subroutine Grid2spec_rc2

  !*************************************************************************

  subroutine Grid2spec_rc3(gridfield,wavefield_uhp)

    ! Physical to spectral transform - just for diagnostic use, as 
    ! product is not dealiased.

    use fft_mod,  only: fft
    use op_rules, only: operator(*)

    real,   dimension(nz,ny,ngrid), intent(in)   :: gridfield
    complex,dimension(nz,nkx,nky),intent(out)    :: wavefield_uhp
    complex,dimension(nz,nkx,-kmax-1:kmax)       :: wavefield
    integer                                      :: n

    wavefield   = sgn*gridfield*cmplx(1.,0.)
    wavefield   = fft(wavefield,forward)
    wavefield_uhp = wavefield(:,:,0:kmax)

    if (first_lhp_pe)    wavefield_uhp(:,1,:) = 0. ! dummy column
    do n=1,nz
       where (kxvec<=0)  wavefield_uhp(n,:,1) = 0. ! LHP ky=0 row
    enddo

  end subroutine Grid2spec_rc3

  !*******************************************************************

  function ir_prod2(f,g) result(prod)

    ! This is for doing the special multiply of physical fields,
    ! assuming that fields on the staggered grid are stored as the
    ! imaginary component.

    complex,dimension(:,:),intent(in)       :: f,g
    complex,dimension(size(f,1),size(f,2))  :: prod
    
    prod = cmplx(real(f)*real(g),aimag(f)*aimag(g))

  end function ir_prod2

  !*******************************************************************

  function ir_prod3(f,g) result(prod)

    ! This is for doing the special multiply of physical fields,
    ! assuming that fields on the staggered grid are stored as the
    ! imaginary component.

    complex,dimension(:,:,:),intent(in)               :: f,g
    complex,dimension(size(f,1),size(f,2),size(f,3))  :: prod
    
    prod = cmplx(real(f)*real(g),aimag(f)*aimag(g))

  end function ir_prod3

  !*******************************************************************

  function ir_pwr2(f,pwr) result(prod)

    ! This is for doing the special power raise of physical fields,
    ! assuming that fields on the staggered grid are stored as the
    ! imaginary component.

    complex,dimension(:,:),intent(in)       :: f
    real,intent(in)                         :: pwr
    complex,dimension(size(f,1),size(f,2))  :: prod
    
    if (pwr==2.) then
       prod = cmplx(real(f)**2.,aimag(f)**2.)
    elseif (pwr==.5) then
       prod = cmplx(real(f)**.5,aimag(f)**.5)
    endif
    
  end function ir_pwr2

  !*******************************************************************

  function ir_pwr3(f,pwr) result(prod)

    ! This is for doing the special power raise of physical fields,
    ! assuming that fields on the staggered grid are stored as the
    ! imaginary component.

    complex,dimension(:,:,:),intent(in)               :: f
    real,intent(in)                                   :: pwr
    complex,dimension(size(f,1),size(f,2),size(f,3))  :: prod
    
    if (pwr==2.) then
       prod = cmplx(real(f)**2.,aimag(f)**2.)
    elseif (pwr==.5) then
       prod = cmplx(real(f)**.5,aimag(f)**.5)
    endif

  end function ir_pwr3


  !*******************************************************************

  function conj_symm2(f) result(fout)

    ! For input f_[k,l], get conj(f_[-k,-l]).  This requires 
    ! interprocessor communication.

    use par_tools, only: par_swapdata

    complex,dimension(:,:),intent(in)          :: f
    complex,dimension(size(f,1),size(f,2))     :: fout
    complex,dimension(size(f,1)-1,size(f,2))   :: tempp
    complex,dimension(size(f,2))               :: tempcol

    ! Prepare flipped and conjugated LHP temp for swapping (in place).
    tempp   = conjg(f(ubound(f,1):2:-1,ubound(f,2):1:-1))
    tempcol = conjg(f(1               ,ubound(f,2):1:-1))

    ! exchange with symmetric processor.  kx_start starts at -kmax-1 on pe 0, 
    ! and kx_start+1 is first kx value on section tempp
    call par_swapdata(tempp, cpu( -(kx_start+1) )) 

    ! exchange first column with assymetric processor
    if (.not.(first_lhp_pe.or.first_rhp_pe)) then 
       call par_swapdata(tempcol, cpu(-kx_start))
    endif

    fout = 0.
    fout(1,:)                = tempcol
    fout(2:ubound(fout,1),:) = tempp

  end function conj_symm2

  !*******************************************************************

  function conj_symm3(f) result(fout)

    ! For input f_[k,l], get conj(f_[-k,-l]).  This requires 
    ! interprocessor communication.

    use par_tools, only: par_swapdata

    complex,dimension(:,:,:),intent(in)                :: f
    complex,dimension(size(f,1),size(f,2),size(f,3))   :: fout
    complex,dimension(size(f,1),size(f,2)-1,size(f,3)) :: tempp
    complex,dimension(size(f,1),size(f,3))             :: tempcol

    ! Prepare flipped and conjugated LHP temp for swapping (in place).
    tempp   = conjg(f(:,ubound(f,2):2:-1,ubound(f,3):1:-1))
    tempcol = conjg(f(:,1               ,ubound(f,3):1:-1))

    ! exchange with symmetric processor  (kx_start starts at zero)
    call par_swapdata(tempp, cpu( -(kx_start+1) ))

    ! exchange first_lhp_pe column with assymetric processor
    if (.not.(first_lhp_pe.or.first_rhp_pe)) then 
       call par_swapdata(tempcol, cpu(-kx_start))
    endif

    fout = 0.
    fout(:,1,:)                = tempcol
    fout(:,2:ubound(fout,2),:) = tempp

  end function conj_symm3

  !*******************************************************************

end module transform_tools
