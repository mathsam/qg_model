module fft_mod                    !-*-f90-*-  <- tells emacs: use f90-mode

  !********************************************************************
  ! Calls complex to complex 2d forward and backward ffts.  
  ! Uses FFTW package 2.1.5, parallelized for MPI, double precision
  !
  ! Routines: init_fft, fft
  !
  ! Dependencies:  You must link with FFTW libraries.
  !
  !********************************************************************

  implicit none
  private
  save

  ! for FFTW
  integer,parameter        :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
  integer,parameter        :: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1
  integer,parameter        :: FFTW_ESTIMATE=1, FFTW_MEASURE=0
  integer,parameter        :: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8
  integer,parameter        :: FFTW_USE_WISDOM=16
  integer,parameter        :: FFTW_THREADSAFE=128

  ! for FFTW  MPI
  integer,parameter        :: FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0
  integer,parameter        :: FFTW_SCRAMBLED_INPUT=8192
  integer,parameter        :: FFTW_SCRAMBLED_OUTPUT=16384
  integer,parameter        :: USE_WORK=1

  ! For my transform
  integer,parameter        :: forward=1, backward=-1  
  integer*8                :: plan_f, plan_b
  integer,dimension(2)     :: n_dim
  integer                  :: total_local_size, ngrid, nkx, nz
  real                     :: scale_f

  complex, dimension(:), allocatable  :: work

  public :: Init_fft, fft

  interface fft
     module procedure fft2, fft3
  end interface

contains

  !********************************************************************

  subroutine Init_fft(kmax, nzi, kx_start, kx_end, nkxout, y_start, y_end, ny)

    ! Initialize fft routine

    use par_tools, only: MPI_COMM_WORLD, processor_id

    integer,intent(in)  :: kmax, nzi
    integer,intent(out) :: kx_start, kx_end, nkxout, y_start, y_end, ny

    ngrid = 2*(kmax+1)
    n_dim = (/ ngrid, ngrid /)
    nz = nzi
    scale_f = 1./(n_dim(1)*n_dim(2))

    ! create FFTW plan for 2D complex to complex forward transform
    call fftwnd_f77_mpi_create_plan(plan_f,MPI_COMM_WORLD,2,n_dim, &
         FFTW_FORWARD,FFTW_ESTIMATE)

    ! local data size - this call returns ny, y_start, nkx, kx_start, 
    ! and total_local_size
    call fftwnd_f77_mpi_local_sizes(plan_f,ny,y_start, &
         nkx,kx_start,total_local_size)

        ! Check local congruence
    if (nkx.ne.ny) then
       print*, 'FATAL: nkx not equal to ny on cpu ',processor_id
       stop
    endif

    y_start = y_start + 1            ! Make it start at 1
    y_end = y_start + ny - 1
    kx_start = kx_start - kmax - 1   ! Make it start at -kmax-1
    kx_end = kx_start + nkx -1
    nkxout = nkx  ! Set dummy variable for return to caller

    ! create FFTW plan for 2D complex to complex backward transform
    call fftwnd_f77_mpi_create_plan(plan_b,MPI_COMM_WORLD,2,n_dim, &
         FFTW_BACKWARD,FFTW_ESTIMATE)

    ! work array allows faster performance.  must be same size as local data.
    if (USE_WORK == 1) then
       allocate( work(nz*total_local_size) ); work=0. 
    else
       allocate( work(1) ); work=0. 
    endif

  end subroutine Init_fft

  !********************************************************************

  function fft2(f,dirn) result(fr)

    ! Calculate 2d complex to complex fft.  
    ! dirn = -1  ==>  backward fft:  spectral --> grid
    ! dirn = +1  ==>  forward  fft:  grid --> spectral

    complex,dimension(nkx,ngrid)  :: f
    complex,dimension(nkx,ngrid)  :: fr
    complex,dimension(ngrid,nkx)  :: temp1
    integer,intent(in)            :: dirn
    real                          :: scale=1.0

    temp1 = transpose(f)
     
    if (dirn==backward) then
       scale=1.0
       call fftwnd_f77_mpi(plan_b,1,temp1,work,USE_WORK,FFTW_TRANSPOSED_ORDER)
    elseif (dirn==forward)  then
       scale=scale_f
       call fftwnd_f77_mpi(plan_f,1,temp1,work,USE_WORK,FFTW_TRANSPOSED_ORDER)
    endif

    fr = scale*transpose(temp1)

  end function fft2

  !********************************************************************

  function fft3(f,dirn) result(fr)
 
    ! Calculate 2d complex to complex fft for each vertical level
    ! dirn = -1  ==>  backward fft:  spectral --> grid
    ! dirn = +1  ==>  forward  fft:  grid --> spectral

    complex,dimension(:,:,:)                          :: f
    complex,dimension(size(f,1),size(f,2),size(f,3))  :: fr
    complex,dimension(size(f,1),size(f,3),size(f,2))  :: temp1
    integer,intent(in)               :: dirn
    real                             :: scale=1.0
    integer                          :: iz, n_z
    
    n_z = size(f,1)
    do iz = 1,n_z
       temp1(iz,:,:) = transpose(f(iz,:,:))
    enddo
    if (dirn==backward) then
       scale=1.0
       call fftwnd_f77_mpi(plan_b,n_z,temp1,work,USE_WORK,FFTW_TRANSPOSED_ORDER)
    elseif (dirn==forward)  then
       scale=scale_f
       call fftwnd_f77_mpi(plan_f,n_z,temp1,work,USE_WORK,FFTW_TRANSPOSED_ORDER)
    endif

    do iz=1,n_z
       fr(iz,:,:) = scale*transpose(temp1(iz,:,:))
    enddo

  end function fft3

  !********************************************************************

end module fft_mod
