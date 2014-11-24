module qg_topo

  !************************************************************************
  ! Topography initialization and storage.
  !
  ! Routines:  init_topo
  !
  ! Dependencies:  io_tools, transform_tools, op_rules, qg_params, 
  !                qg_arrays, numerics_lib, par_tools
  !
  !************************************************************************

  implicit none
  private
  save

  complex,dimension(:,:),  allocatable   :: hb, toposhift

  public :: init_topo, hb, toposhift

contains

  !************************************************************************

  subroutine init_topo

    !************************************************************************
    ! Read in or create bottom topography in a manner specified by 
    ! 'topo_type', which can have the values:
    !   
    !   gaussbump : Gaussian bump in center of domain with width 'del_topo'.
    !               Note that in physical space, -pi<=x,y<=+pi
    !   
    !   spectral  : Random topography whose representation in spectral
    !               space is isotropic and centered at 'k_o_topo' with 
    !               width 'del_topo'
    !
    !   xslope    : Slope in zonal direction 
    !
    !   yslope    : Slope in meridional direction (equiv to beta in bottom 
    !               layer only).
    !
    !   read      : Read in initial spectral topography from 'hb_in_file'
    !
    ! All parameters read in via input namelist file.  Also, magnitude
    ! of topography is set by 'toposcale', regardless of topo type.
    ! If not set, its default value is 1.
    !************************************************************************
  
    use io_tools,        only: Message, Write_field, Read_field
    use qg_params,       only: kmax, kx_start, kx_end, ngrid,nky, nkx, ny, nz, &
                               toposcale, hb_in_file,                           &
                               del_topo,k_o_topo,pi,i,idum,hb_file,             &
                               parameters_ok,cr,ci,io_root,                     &
                               topo_type, restarting
    use qg_strat_and_shear, only: ubar, vbar
    use qg_arrays,       only: ksqd_, kx_, ky_
    use transform_tools, only: Grid2spec
    use numerics_lib,    only: Ran
    use par_tools,       only: par_scatter, par_gather, par_sync


    complex,dimension(:,:), allocatable      :: hb_global
    real,dimension(:,:),allocatable          :: hbg
    integer                                  :: midx, midy, ix, iy
    real                                     :: radius2
   

    ! Allocate permanent arrays
    allocate(hb       (kx_start:kx_end,0:kmax))
    allocate(toposhift(kx_start:kx_end,0:kmax))

    if (restarting) then

       allocate(hb_global(-kmax-1:kmax,0:kmax));    hb_global = 0.
       call Read_field(hb_global(-kmax:kmax,:),hb_file)
       call par_scatter(hb_global,hb,io_root)
       deallocate(hb_global)

    else 

       call Message('Using bottom topography')
       midx = kmax+1; midy = kmax+1

       select case (trim(topo_type))

       case('gaussbump')
          
          call Message('Gaussian bump topography selected')

          if (del_topo==cr) then
             call Message('Error: bump width del_topo must be set to make &
                           &gaussian bump topography')
             Parameters_ok=.false.
          endif

          if (del_topo<=0) then
             call Message('Error: require del_topo>=0 - yours is:',&
                           r_tag=del_topo)
             Parameters_ok=.false.
          endif
          
          allocate(hbg(ny,ngrid)); hbg = 0.
          do iy = 1,ny
             do ix = 1,ngrid
                radius2 = ((ix-midx)**2+(iy-midy)**2)*((2*pi)**2/ngrid**2)
                hbg(iy,ix) = exp(-radius2/del_topo**2)
             enddo
          enddo
          call Grid2spec(hbg,hb)  
          deallocate(hbg)
          
       case('spectral')           ! Make a spectral peak at k_o_topo
          
          call Message('Spectral topography selected')
          if (del_topo==cr) then
             call Message('Error: spectral width del_topo must be set for&
                  & spectral topography')
             Parameters_ok=.false.
          endif
          if (del_topo<=0) then
             call Message('Error: require del_topo>=0 - yours is:',&
                           r_tag=del_topo)
             Parameters_ok=.false.
          endif
          if (k_o_topo==cr) then
             call Message('Error: spectral centroid k_o_topo must be set to &
                  &make spectral topography')
             Parameters_ok=.false.
          endif
          
          hb = exp(-(sqrt(ksqd_)-k_o_topo)**2/del_topo**2)*cexp(i*2*pi*Ran(idum,nkx,nky))
          
       case('xslope')
          
          call Message('Linear zonal slope topography selected')
          
          allocate(hbg(ny,ngrid)); hbg = 0.
          do iy = 1,ny
             do ix = 1,ngrid
                hbg(iy,ix) = (ix-1)*(2*pi/ngrid)
             enddo
          enddo
          call Grid2spec(hbg,hb)  
          deallocate(hbg)
          
       case('yslope')
          
          call Message('Linear meridional slope topography selected')
          
          allocate(hbg(ny,ngrid)); hbg = 0.
          do iy = 1,ny
             do ix = 1,ngrid
                hbg(iy,ix) = (iy-1)*(2*pi/ngrid)
             enddo
          enddo
          call Grid2spec(hbg,hb)  
          deallocate(hbg)
          
       case('read')
          
          if (trim(hb_file)=='') then
             call Message('Error: no input file for topography given')
             Parameters_ok=.false.
          endif
          call Message('Will read spectral topography from: '//trim(hb_in_file))

          allocate(hb_global(-kmax-1:kmax,0:kmax));    hb_global = 0.
          call Read_field(hb_global(-kmax:kmax,:),hb_in_file,exclude_dd=1)
          call par_scatter(hb_global,hb,io_root)
          deallocate(hb_global)

       case default
          
          call Message('Error: must select topo_type = &
                       &spectral|gaussbump|xslope|yslope|read&
                       & - yours is:'//trim(topo_type))
          Parameters_ok=.false.
          
       end select
       
       if (toposcale==cr) then
          call Message('Warning: toposcale not set - setting to 1')
          toposcale = 1.
       endif
       
       hb = toposcale*hb
       
       allocate(hb_global(-kmax-1:kmax,0:kmax));    hb_global = 0.
       call par_gather(hb,hb_global,io_root)
       call Write_field(hb_global(-kmax:kmax,:),hb_file)
       deallocate(hb_global)

    endif

    toposhift = - i*ubar(nz)*(kx_*hb)- i*vbar(nz)*(ky_*hb)

    call par_sync

  end subroutine init_topo


end module qg_topo
