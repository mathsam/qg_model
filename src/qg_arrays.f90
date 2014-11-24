module qg_arrays                !-*-f90-*-

  !************************************************************************
  ! This module contains all of the dynamically allocated field variables.
  !
  ! Routines:  init_qg_arrays
  !
  ! Dependencies:  io_tools, qg_params, par_tools
  !
  !************************************************************************

  implicit none
  public
  save

  ! Main arrays
  complex,dimension(:,:,:),allocatable   :: q, q_o, psi, psi_o
  complex,dimension(:,:,:),allocatable   :: rhs
  complex,dimension(:,:,:),allocatable   :: b, b_o, rhs_b

  ! Quad drag
  complex,dimension(:,:,:),allocatable   :: ug, vg, qxg, qyg
  complex,dimension(:,:),  allocatable   :: qdrag, qdrag2, unormbg

  ! Filter
  real,   dimension(:,:),  allocatable   :: filter

  ! Wavenumber arrays
  real,   dimension(:,:),  allocatable   :: ksqd_, kx_, ky_
  integer,dimension(:),    allocatable   :: kxv, kyv, lin2kx, lin2ky

contains

  !************************************************************************

  subroutine init_qg_arrays

    use qg_params,       only: kx_start, kx_end, kmax, nkx, nky, nz, ny, ngrid, &
                               quad_drag
    use io_tools,        only: message
    use par_tools,       only: par_sync

    integer                                 :: kx, ky

    ! The first column of spectral arrays is a padding column at 
    ! kx = -kmax-1. 

    ! Allocate primary fields (velocities, pv, streamfunc, etc)
    allocate(ug     (nz,ny,ngrid));                         ug = 0.
    allocate(vg     (nz,ny,ngrid));                         vg = 0.
    allocate(qxg    (nz,ny,ngrid));                         qxg = 0.
    allocate(qyg    (nz,ny,ngrid));                         qyg = 0.
    allocate(q      (1:nz,kx_start:kx_end,0:kmax));         q = 0. 
    allocate(q_o    (1:nz,kx_start:kx_end,0:kmax));         q_o = 0.
    allocate(psi    (1:nz,kx_start:kx_end,0:kmax));         psi = 0.
    allocate(psi_o  (1:nz,kx_start:kx_end,0:kmax));         psi_o = 0. 
    allocate(rhs    (1:nz,kx_start:kx_end,0:kmax));         rhs = 0.
    allocate(filter (     kx_start:kx_end,0:kmax));         filter = 1.

    ! Initialize quadratic drag arrays
    if (quad_drag/=0) then
       allocate(qdrag2(kx_start:kx_end,0:kmax));            qdrag2 = 0.
       allocate(qdrag(kx_start:kx_end,0:kmax));             qdrag = 0.
       allocate(unormbg(ny,ngrid));                         unormbg = 0.
    endif

    call Message('Primary fields allocated')

    ! Store values of kx, ky and k^2 in 2d arrays.
    allocate(ksqd_  (kx_start:kx_end,0:kmax));              ksqd_= 0.
    allocate(kx_    (kx_start:kx_end,0:kmax));              kx_ = 0.
    allocate(ky_    (kx_start:kx_end,0:kmax));              ky_= 0.
    allocate(kxv(nkx), kyv(nky));              kxv = 0.; kyv = 0.

    kxv = (/ (kx, kx = kx_start, kx_end) /)
    if (kx_start == -kmax-1) kxv(1) = 0         ! dummy column
    kyv = (/ (ky,ky=0,kmax) /)
    kx_ = float(spread(kxv,2,nky))
    ky_ = float(spread(kyv,1,nkx))
    ksqd_ = kx_**2 + ky_**2
    
    ! 0,0 never used - this way can divide by K^2 w/o worry
    where (kx_==0.and.ky_==0) ksqd_ = 0.1

    call par_sync
    call Message('Wavenumber arrays initialized')

  end subroutine init_qg_arrays


end module qg_arrays
