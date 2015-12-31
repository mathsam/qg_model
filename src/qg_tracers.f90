module qg_tracers            !-*-f90-*-

  !************************************************************************
  ! This module contains all of the tracer arrays and routines
  !
  ! Routines:     init_tracers, step_tracers
  !
  ! Dependencies: op_rules, io_tools, qg_arrays, qg_params, 
  !               qg_strat_and_shear, transform_tools, numerics_lib, 
  !               strat_tools, par_tools
  !
  !************************************************************************

  implicit none
  private
  save

  complex,dimension(:,:,:),allocatable   :: tracer_y, tracer_y_o
  complex,dimension(:,:,:),allocatable   :: tracer_x, tracer_x_o
  complex,dimension(:,:,:),allocatable   :: tracer_bt,tracer_bt_o !tracer advected by barotropic flow
  complex,dimension(:,:,:),allocatable   :: rhs_tx, rhs_ty, rhs_tbt
  complex,dimension(:,:,:),allocatable   :: psim, psi_stir, psi_bt_stir
  real,   dimension(:,:),  allocatable   :: vmodet, filter_t, op_tracer
  real,   dimension(:),    allocatable   :: dzt, ubart, vbart, deltaz
  real,   dimension(:),    allocatable   :: ubar_proj, vbar_proj !ubar projected to modes 1:maxmode
  real,   dimension(:),    allocatable   :: vmode_x_dz

  public :: tracer_x, tracer_y, tracer_bt, filter_t, dzt, psi_stir, psi_bt_stir, &
            init_tracers, step_tracers
            
contains

  !************************************************************************

  subroutine Init_tracers
    
    use op_rules,           only: operator(+), operator(-), operator(*)
    use io_tools,           only: message, read_field
    use qg_filter_tools,    only: init_filter
    use qg_strat_and_shear, only: dz, ubar, vbar, vmode
    use strat_tools,        only: Layer2mode, Mode2layer
    use qg_params,          only: kx_start, kx_end, kmax, nz, nzt,    &
                                  cr, ci, parameters_ok, maxmode,     &
                                  uscale, vscale,                     &
                                  kappa_v,                            &
                                  use_tracer_x, use_tracer_y,         &
                                  use_tracer_bt,                      &
                                  tracer_x_init_file,                 & 
                                  tracer_y_init_file,                 &
                                  tracer_bt_init_file,                &
                                  filter_type_t, filter_exp_t,        &
                                  filt_tune_t, dealiasing_t,          &
                                  k_cut_t,                            &
                                  dzt_in_file,                        &
                                  ubart_in_file,                      & 
                                  vbart_in_file,                      & 
                                  vmodet_in_file,                     &
                                  io_root,                            &
                                  reset_tracer
    use par_tools,          only: par_bcast
    real,   dimension(:),    allocatable   :: ubarm
    

    call Message('Tracers on')
    call Message('Any following messages regarding filter parameters')
    call Message('  refer to tracer filter: append _t to variable names')
    
    if (nzt==ci) then
       nzt = nz
       dzt = dz
       call Message('nz_tracer = ',nzt)
    endif
    if (maxmode==ci) then
       maxmode = nz-1
    elseif (maxmode>=nz.or.maxmode<0) then
       call Message('Error: need 0<=maxmode<nz')
       parameters_ok=.false.
    endif
    call Message('maxmode = ',maxmode)
    if (kappa_v/=0) call Message('Vertical tracer diffusivity kappa_v = ',r_tag=kappa_v)

    ! Set up tracer filter

    allocate(filter_t(kx_start:kx_end,0:kmax))
    filter_t = Init_filter(filter_type_t,filter_exp_t,k_cut_t, &
                           dealiasing_t,filt_tune_t,1.)
    call Message('Tracer filter created')
     
    ! Init stirring fields
    allocate(psim    (1:nz ,kx_start:kx_end,0:kmax));         psim = 0.
    allocate(psi_stir(1:nzt,kx_start:kx_end,0:kmax));         psi_stir = 0.
    allocate(psi_bt_stir (1,kx_start:kx_end,0:kmax));         psi_bt_stir = 0.

    if (use_tracer_x .OR. use_tracer_y) then
       allocate(ubar_proj(nz))
       allocate(vbar_proj(nz))
       allocate(ubarm(nz))
       ubarm     = layer2mode(ubar, vmode, dz, maxmode)
       ubar_proj = mode2layer(ubarm, vmode, maxmode)
       ubarm     = layer2mode(vbar, vmode, dz, maxmode)
       vbar_proj = mode2layer(ubarm, vmode, maxmode)
       deallocate(ubarm)
    endif    

    if (use_tracer_bt) then
       allocate(vmode_x_dz(nz))
       vmode_x_dz(:) =  vmode(:,1)*dz(:)
       allocate(tracer_bt  (1,kx_start:kx_end,0:kmax));    tracer_bt   = 0.
       allocate(tracer_bt_o(1,kx_start:kx_end,0:kmax));    tracer_bt_o = 0.
       allocate(rhs_tbt    (1,kx_start:kx_end,0:kmax));    rhs_tbt     = 0.
       tracer_bt = filter_t* &
            make_tracer('',tracer_bt_init_file,"tracer_bt", 1)
       call Message('tracer_bt fields allocated')
    endif

    if (use_tracer_x) then
       allocate(tracer_x  (1:nzt,kx_start:kx_end,0:kmax));    tracer_x = 0.
       allocate(tracer_x_o(1:nzt,kx_start:kx_end,0:kmax));    tracer_x_o = 0.
       allocate(rhs_tx    (1:nzt,kx_start:kx_end,0:kmax));    rhs_tx = 0.
       tracer_x = filter_t*  &
            make_tracer('',tracer_x_init_file,"tracer_x", nzt)
       call Message('Tracer_x fields allocated')
    endif

    if (use_tracer_y) then
       allocate(tracer_y  (1:nzt,kx_start:kx_end,0:kmax));     tracer_y = 0.
       allocate(tracer_y_o(1:nzt,kx_start:kx_end,0:kmax));     tracer_y_o = 0.
       allocate(rhs_ty    (1:nzt,kx_start:kx_end,0:kmax));     rhs_ty = 0.
       tracer_y = filter_t*  &
            make_tracer('',tracer_y_init_file,"tracer_y", nzt)
       call Message('Tracer_y fields allocated')
    endif

    reset_tracer = .false.

    ! Get strat data for extrapolation to higher number of layers
    allocate(dzt(1:nzt));                                     dzt = 0.
    allocate(ubart(1:nzt));                                   ubart = 0.
    allocate(vbart(1:nzt));                                   vbart = 0.
    allocate(vmodet(1:nzt,1:nz));                             vmodet = 0.
    if (nzt>nz) then
       call Message('Will read dzt from '//trim(dzt_in_file))
       call Read_field(dzt,dzt_in_file,exclude_dd=1)
       call par_bcast(dzt,io_root)
       dzt = dzt/sum(dzt)
       call Message('Will read vmodet from '//trim(vmodet_in_file))
       call Read_field(vmodet,vmodet_in_file,exclude_dd=1)
       call par_bcast(vmodet,io_root)
       if (uscale/=0) then
          call Message('Will read ubart from '//trim(ubart_in_file))
          call Read_field(ubart,ubart_in_file,exclude_dd=1)
          call par_bcast(ubart,io_root)
       else
          ubart = 0.
       endif
       if (vscale/=0) then
          call Message('Will read vbart from '//trim(vbart_in_file))
          call Read_field(vbart,vbart_in_file,exclude_dd=1)
          call par_bcast(vbart,io_root)
       else
          vbart = 0.
       endif
    else if (nzt==nz) then
       dzt = dz
       vmodet = vmode
       ubart = ubar
       vbart = vbar
    else
       call Message('nzt < nz: take caution of what you are doing!')
       dzt = 1.0
       vmodet = 1.0
    end if

    ! Initialize operator for implicit method if vertical diffusion turned on
    if (kappa_v/=0) then

       if (nzt==1) then 
          call Message('Error: init_tracers: set kappa_v=0 if nzt=1')
          parameters_ok=.false.
       endif

       allocate(deltaz(1:nzt-1),op_tracer(1:nzt,-1:1))
       op_tracer = 0.
       deltaz    = 0.5*(dzt(1:nzt-1) + dzt(2:nzt))

       op_tracer(2:nzt  ,-1) =  kappa_v/(dzt(2:nzt)  *deltaz(1:nzt-1))
       op_tracer(1:nzt-1, 1) =  kappa_v/(dzt(1:nzt-1)*deltaz(1:nzt-1)) 
       op_tracer(:      , 0) = -( op_tracer(:,-1) + op_tracer(:,1) )

    endif

    tracer_x_o  = tracer_x
    tracer_y_o  = tracer_y
    tracer_bt_o = tracer_bt

    call get_rhs_tracer

    call Message('Tracer(s) initialized')
    
  end subroutine init_tracers

  !************************************************************************

  function make_tracer(tracer_file, tracer_init_file, tracer_name, num_layers) &
       result(tracer)

    !************************************************************************
    ! Read in or create initial tracer distribution in manner specified 
    ! by 'tracer_init_type', which can have the values:
    !
    !   spectral_m         : Gaussian spread of initial energy about isotropic horiz.
    !                        wavenumber 'k_o' and with width 'delk', and all energy
    !                        in vertical mode 'm_o'
    !   read               : read in inital spectral tracer distribution 
    !                        from file 'tracer_init_file'
    !
    ! In all but 'read' case, total initial variance is set by 'tvar_o'.
    ! Default is tvar_o = 0.
    !************************************************************************

    use op_rules,        only: operator(+), operator(-), operator(*)
    use io_tools,        only: Message, Read_field
    use qg_params,       only: kx_start, kx_end, kmax, ngrid, nkx, nky, nz,      &
                               parameters_ok, cr, ci, io_root,                   &
                               frame, start_frame, start_frame_t, rewindfrm,     &
                               tracer_init_type, restarting,                     &
                               nc_restartfile, reset_tracer
    use qg_arrays,       only: ksqd_
    use qg_strat_and_shear, only: dz, vmode
    use transform_tools, only: Grid2spec
    use strat_tools,     only: Layer2mode, Mode2layer
    use numerics_lib,    only: ran
    use par_tools,       only: par_scatter, par_sync
    use nc_io_tools,     only: read_nc

    integer, intent(in)                          :: num_layers
    complex,dimension(1:num_layers,1:nkx,0:kmax) :: tracer
    character(*), intent(in)                     :: tracer_file
    character(*), intent(in)                     :: tracer_init_file
    character(*), intent(in)                     :: tracer_name
    complex,dimension(:,:,:),allocatable         :: tracer_global
    real                                         :: tv

    if (restarting .and. (.not. reset_tracer)) then

       allocate(tracer_global(1:num_layers,-kmax-1:kmax,0:kmax)); tracer_global=0.
       call read_nc("./INPUT/"//trim(nc_restartfile), tracer_name, &
                       tracer_global(:,-kmax:kmax,:)) 
       call par_scatter(tracer_global,tracer,io_root)
       deallocate(tracer_global)

    else 

       select case (trim(tracer_init_type))    
       case ('')
          tracer = 0.

       case ('read')
          if (trim(tracer_file)=='') then
             call Message('Error: no input file for tracer given')
             Parameters_ok=.false.
          endif

          if (start_frame_t==ci) then
             call Message('Warning: start_frame_t not initialized-setting to 1')
             start_frame_t = 1
          elseif (start_frame_t<=0) then
             call Message('Error: require start_frame_t>=0')
             Parameters_ok=.false.
          endif
          call Message('Initial tracer will be read from: '&
               &//trim(tracer_init_file)//', frame:', tag=start_frame_t)


          allocate(tracer_global(1:num_layers,-kmax-1:kmax,0:kmax)); tracer_global=0.
          call Message('Will read initial tracer field from: '&
               &//trim(tracer_init_file))
          call Read_field(tracer_global(:,-kmax:kmax,:),tracer_init_file, &
               frame=start_frame_t,exclude_dd=1,zfirst=1)
          call par_scatter(tracer_global,tracer,io_root)
          deallocate(tracer_global)

       case default
          call Message('Error: legal tracer init types are &
               &tracer_init_type=read|spatially_centered|spatially_constant&
               & - yours is:'//trim(tracer_init_type))
          parameters_ok=.false.
       end select

    end if

    call par_sync

  end function  make_tracer

  !*********************************************************************

  subroutine step_tracers

    use op_rules,      only: operator(+), operator(-), operator(*)
    use numerics_lib,  only: march, march_implicit
    use qg_params,     only: dt, call_tx, call_ty, call_tbt, robert, &
                             use_tracer_x, use_tracer_y, use_tracer_bt, &
                             kappa_v

    if (kappa_v==0) then

       if (use_tracer_x) tracer_x = filter_t*March(tracer_x,tracer_x_o,rhs_tx,dt,robert,call_tx)
       if (use_tracer_y) tracer_y = filter_t*March(tracer_y,tracer_y_o,rhs_ty,dt,robert,call_ty)

    else  ! Use implicit method to calculate vertical diffusion and step tracers
       
       if (use_tracer_x) tracer_x = filter_t*March_implicit(tracer_x,tracer_x_o,rhs_tx, &
                                                            op_tracer,dt,robert,call_tx)
       if (use_tracer_y) tracer_y = filter_t*March_implicit(tracer_y,tracer_y_o,rhs_ty, &
                                                            op_tracer,dt,robert,call_ty)
    endif

    if (use_tracer_bt) tracer_bt = filter_t*March(tracer_bt,tracer_bt_o,rhs_tbt,dt,robert,call_tbt)

    call get_rhs_tracer

  end subroutine step_tracers

  !*********************************************************************

  subroutine  get_rhs_tracer 

    use op_rules,           only: operator(+), operator(-), operator(*)
    use strat_tools,        only: Layer2mode, Mode2layer
    use qg_arrays,          only: psi, kx_, ky_, ksqd_
    use qg_strat_and_shear, only: dz, vmode, ubar, vbar
    use qg_params,          only: maxmode, nz, nzt, i, kmax, nkx, nky, kappa_v, &
                                  kappa_h, use_tracer_x, use_tracer_y, &
                                  use_tracer_bt, use_mean_grad_t
    use transform_tools,    only: jacob
    use numerics_lib,       only: diffz2

    integer :: iz, ix, iy

    if (maxmode<nz-1) then
       ! Stir tracer with modes 0 through maxmode
       if (nzt /= 1) then
           psim = layer2mode(psi,vmode,dz,maxmode)
           psi_stir = mode2layer(psim,vmodet,maxmode)
       else
           do iy = 1, size(psi, 3)
             do ix = 1, size(psi, 2)
             psi_stir(1,ix,iy) = dot_product(vmode(:,1)*dz, psi(:,ix,iy))
             enddo
           enddo
       endif
    else
       psi_stir = psi
    endif

    if (use_tracer_bt) then
        do iy = 1, size(psi, 3)
          do ix = 1, size(psi, 2)
          psi_bt_stir(1,ix,iy) = dot_product(vmode_x_dz, psi(:,ix,iy))
          enddo
        enddo
        call Jacob(tracer_bt(1,:,:), psi_bt_stir(1,:,:), rhs_tbt(1,:,:))
        if (use_mean_grad_t) rhs_tbt = rhs_tbt - i*(kx_*psi_bt_stir)
    endif

    if (use_tracer_x) then
       if (nzt /= nz) then
          do iz=1,nzt
             call Jacob(tracer_x(iz,:,:),psi_stir(iz,:,:),rhs_tx(iz,:,:))
          enddo
       else
          call Jacob(tracer_x,psi_stir,rhs_tx)
       endif
!       if (any(ubar/=0)) rhs_tx = rhs_tx - i*(kx_*(ubart*tracer_x))
!       if (any(vbar/=0)) rhs_tx = rhs_tx - i*(ky_*(vbart*tracer_x))
       ! Note: above form does not take into account time varying mean state, when
       ! used, but below form does not allow higher res tracer.
       if (any(ubar_proj/=0)) rhs_tx = rhs_tx - i*(kx_*(ubar_proj*tracer_x))
       if (any(vbar_proj/=0)) rhs_tx = rhs_tx - i*(ky_*(vbar_proj*tracer_x))
       if (use_mean_grad_t) rhs_tx = rhs_tx + i*(ky_*psi_stir)
!       if (kappa_v/=0) rhs_tx = rhs_tx + kappa_v*diffz2(tracer_x_o,dzt)
       if (kappa_h/=0) rhs_tx = rhs_tx - kappa_h*ksqd_*tracer_x_o

    endif

    if (use_tracer_y) then
       if (nzt /= nz) then
          do iz=1,nzt
             call Jacob(tracer_y(iz,:,:),psi_stir(iz,:,:),rhs_ty(iz,:,:))
          enddo
       else
          call Jacob(tracer_y,psi_stir,rhs_ty)
       endif
!       if (any(ubar/=0)) rhs_ty = rhs_ty - i*(kx_*(ubart*tracer_y))
!       if (any(vbar/=0)) rhs_ty = rhs_ty - i*(ky_*(vbart*tracer_y))
       ! Note: above form does not take into account time varying mean state, when
       ! used, but below form does not allow higher res tracer.
       if (any(ubar_proj/=0)) rhs_ty = rhs_ty - i*(kx_*(ubar_proj*tracer_y))
       if (any(vbar_proj/=0)) rhs_ty = rhs_ty - i*(ky_*(vbar_proj*tracer_y))
       if (use_mean_grad_t) rhs_ty = rhs_ty - i*(kx_*psi_stir)
!       if (kappa_v/=0) rhs_ty = rhs_ty + kappa_v*diffz2(tracer_y_o,dzt)
       if (kappa_h/=0) rhs_ty = rhs_ty - kappa_h*ksqd_*tracer_y_o

    endif
    
  end subroutine  get_rhs_tracer

  !*********************************************************************

end module qg_tracers
