module qg_diagnostics

  !************************************************************************
  ! Energetics and other diagnostics
  ! 
  ! Routines:  get_energetics, get_spectra, energy, enstrophy
  !
  ! Dependencies: qg_arrays, qg_params, io_tools, op_rules, transform_tools,
  !               qg_diag_tools, rmf_forcing, qg_tracers, par_tools,
  !               qg_strat_and_shear, strat_tools, numerics_lib
  !
  !************************************************************************

  implicit none
  private
  integer :: diag1_time_var_id, ke_var_id, ens_var_id, meanke_var_id,    &
             filter_rate_var_id, bd_rate_var_id, qd_rate_var_id,         &
             td_rate_var_id, gen_rmf_rate_var_id, ape_var_id,            &
             thd_rate_var_id, gen_bci_rate_var_id, filter_rate_tx_var_id,&
             tvarx_var_id, tvary_var_id,                                 &
             gen_tg_tx_var_id, filter_rate_ty_var_id, gen_tg_ty_var_id,  &
             eddy_time_var_id
  integer :: diagspec_time_var_id, kes_var_id, kesx_var_id, kesy_var_id, &
             gens_var_id, gens_rmf_var_id, kems_var_id, kemsx_var_id,    &
             kemsy_var_id, apems_var_id, apes_var_id, genms_var_id,      &
             thdms_var_id, bdms_var_id, tdms_var_id, qdms_var_id,        &
             filterms_var_id, xferms_var_id, energy_xfers_var_id,        &
             enstrophy_xfers_var_id, txvars_var_id, txfluxs_var_id,      &
             tyvars_var_id, tyfluxs_var_id, uv_avg_x_var_id,             &
             vq_avg_x_var_id
  save

  public :: Get_energetics, Get_spectra, energy, enstrophy, get_corrs,   &
            init_get_energetics, init_get_spectra

  real                                      :: tvary0, tvarx0
  real,dimension(:),allocatable             :: stcor0,shcor0
  complex,dimension(:,:,:),allocatable      :: psi0

contains

  !*************************************************************************

  real function energy(psi)

    ! Get total energy
    
    use op_rules,  only:  operator(*)
    use qg_params, only:  nz, F
    use qg_arrays, only:  ksqd_ 
    use qg_strat_and_shear, only: dz, drho
    use par_tools, only:  par_sum, processor_id

    complex,dimension(:,:,:),intent(in) :: psi
    real                                :: ke=0., ape=0.

    ke = sum(dz*(ksqd_*(psi*conjg(psi))))
    call par_sum(ke)

    if (nz>1) then
       ape = F*sum((conjg(psi(2:nz,:,:)-psi(1:nz-1,:,:))* &
                         (psi(2:nz,:,:)-psi(1:nz-1,:,:)))*(1/drho))
       call par_sum(ape)
    else
       ape = F*sum(psi*conjg(psi))
       call par_sum(ape)
    endif

    energy = ke + ape

  end function energy
 
  !*************************************************************************

  real function enstrophy(q)

    ! Get total enstrophy
    
    use op_rules,  only: operator(*)
    use par_tools, only: par_sum
    use qg_strat_and_shear, only: dz

    complex,dimension(:,:,:),intent(in) :: q
    
    enstrophy = sum( dz*q*conjg(q) )
    call par_sum(enstrophy)

  end function enstrophy

  !*************************************************************************

  subroutine Get_corrs(dframe)

    !************************************************************************
    ! Compute horizontally averaged time correlation of strain and shear 
    !************************************************************************

    use op_rules,           only: operator(+), operator(-), operator(*)
    use io_tools,           only: Write_field, Message
    use qg_params,          only: nz,kx_start,kx_end,kmax
    use qg_arrays,          only: ksqd_, psi
    use qg_strat_and_shear, only: dz
    use numerics_lib,       only: diffz
    use par_tools,          only: par_sum

    integer,intent(in)                        :: dframe
    real,dimension(:),allocatable             :: stcor,shcor
    logical,save                              :: called_yet=.false.

    if (.not.called_yet) then

       call Message('Initializing correlation calculator')
       ! Allocate and store initial streamfunction
       allocate(stcor0(nz),shcor0(nz)); stcor0=0.; shcor0=0.
       allocate(psi0(nz,kx_start:kx_end,0:kmax))
       psi0 = psi
       stcor0 = sum(sum((ksqd_*ksqd_) * real( psi0*conjg(psi0) ),2),2)      
       call par_sum(stcor0)
       shcor0 = sum(sum(real( diffz(psi0,dz) * conjg(diffz(psi0,dz))),2),2)
       call par_sum(shcor0)

       called_yet=.true.

    end if
       
    allocate(stcor(nz),shcor(nz));  stcor=0.; shcor=0.

    stcor = sum(sum((ksqd_*ksqd_) * real( psi*conjg(psi0) ),2),2)
    call par_sum(stcor)
    stcor = stcor/stcor0
    call write_field(stcor,'straincorr',dframe)

    shcor  = sum(sum(real( diffz(psi,dz)  * conjg(diffz(psi0,dz))),2),2)
    call par_sum(shcor)
    shcor = shcor/shcor0
    call write_field(shcor,'shearcorr',dframe)

  end subroutine Get_corrs

  !*************************************************************************

  subroutine init_get_energetics()

    !************************************************************************
    ! Prepare energetics.nc file and variables for Get_energetics function
    ! below to use
    !************************************************************************

    use io_tools,    only: Message
    use qg_params,   only: bot_drag,top_drag,therm_drag,                     &
                           use_forcing, use_tracer_x,use_tracer_y,           &
                           quad_drag,surface_bc,filter_type,filter_type_t,   &
                           time_varying_mean,use_mean_grad_t,                &
                           nz, F
    use nc_io_tools, only: create_file, enddef_file, register_variable,      &
                           create_axis_time

    integer :: energetics_file_id, axis_time_id

    energetics_file_id = create_file("energetics.nc")
    axis_time_id       = create_axis_time(energetics_file_id)

    diag1_time_var_id = register_variable(energetics_file_id, "time", &
                                    (/axis_time_id/), .false.)

    ke_var_id   = register_variable(energetics_file_id, "ke", &
                                    (/axis_time_id/), .false.)

    ens_var_id  = register_variable(energetics_file_id, "ens", &
                                    (/axis_time_id/), .false.)

    if (time_varying_mean) then
        meanke_var_id = register_variable(energetics_file_id, "meanke", &
                                    (/axis_time_id/), .false.)
    endif

    if (trim(filter_type)/='none') then
        filter_rate_var_id = register_variable(energetics_file_id, "filter_rate", &
                                    (/axis_time_id/), .false.)
    endif

    if (bot_drag/=0) then
        bd_rate_var_id = register_variable(energetics_file_id, "bottom_drag_rate", &
                                    (/axis_time_id/), .false.)
    endif

    if (quad_drag/=0) then
        qd_rate_var_id = register_variable(energetics_file_id, "quad_drag_rate", &
                                    (/axis_time_id/), .false.)
    endif

    if (top_drag/=0) then
        td_rate_var_id = register_variable(energetics_file_id, "top_drag_rate", &
                                    (/axis_time_id/), .false.)
    endif

    if (use_forcing) then
        gen_rmf_rate_var_id = register_variable(energetics_file_id, "gen_rmf_rate", &
                                    (/axis_time_id/), .false.)
    endif

    if (nz>1 .or. F/=0) then
        ape_var_id = register_variable(energetics_file_id, "ape", &
                                    (/axis_time_id/), .false.)
        if (therm_drag/=0) then
            thd_rate_var_id = register_variable(energetics_file_id, "therm_drag_rate", &
                                    (/axis_time_id/), .false.)
        endif

        if (nz>1) then
            gen_bci_rate_var_id = register_variable(energetics_file_id, "gen_bci_rate", &
                                    (/axis_time_id/), .false.)
        endif
    endif

    if (use_tracer_x) then
        tvarx_var_id = register_variable(energetics_file_id, "tvarx", &
                                    (/axis_time_id/), .false.)

        if(trim(filter_type_t)/='none') then
            filter_rate_tx_var_id = register_variable(energetics_file_id, "filter_rate_tx", &
                                    (/axis_time_id/), .false.)
        endif

        if (use_mean_grad_t) then
            gen_tg_tx_var_id = register_variable(energetics_file_id, "gen_tg_tx", &
                                    (/axis_time_id/), .false.)
        endif
    endif

    if (use_tracer_y) then
        tvary_var_id = register_variable(energetics_file_id, "tvary", &
                                    (/axis_time_id/), .false.)

        if(trim(filter_type_t)/='none') then
            filter_rate_ty_var_id = register_variable(energetics_file_id, "filter_rate_ty", &
                                    (/axis_time_id/), .false.)
        endif

        if (use_mean_grad_t) then
            gen_tg_ty_var_id = register_variable(energetics_file_id, "gen_tg_ty", &
                                    (/axis_time_id/), .false.)
        endif
    endif

    eddy_time_var_id = register_variable(energetics_file_id, "eddy_time", &
                                    (/axis_time_id/), .false.)

    call enddef_file(energetics_file_id)
    call Message("energetics.nc initialized")

  end subroutine init_get_energetics


  function Get_energetics(framein) result(dframe)

    !************************************************************************
    ! Calculate KE, APE, and ENS, as well as all rates of generation
    ! and dissipation.  Also calculates current eddy rotation period.
    ! All energetics factors multiplied by 2 since we are only using
    ! upper-half plane spectral fields
    !************************************************************************
    
    use op_rules,    only: operator(+), operator(-), operator(*)
    use io_tools,    only: Message, Write_field
    use qg_arrays,   only: ksqd_, kx_, ky_, psi, q, qdrag, filter
    use rmf_forcing, only: force_o, get_gen_rmf_rate
    use qg_strat_and_shear, only: dz,drho,shearu,shearv,ubar,vbar
    use qg_tracers,  only: tracer_x, tracer_y, psi_stir, dzt, filter_t
    use qg_params,   only: dt,i,pi,F,Fe,kmax,nkx,nky,nz,nzt,                 &
                           bot_drag,top_drag,therm_drag,                     &
                           time,cntr,uscale,vscale,use_forcing,              &
                           use_tracer_x,use_tracer_y,                        &
                           quad_drag,surface_bc,filter_type,filter_type_t,   &
                           time_varying_mean,use_mean_grad_t
    use qg_filter_tools, only: get_filter_rate
    use par_tools,    only: par_sum
    use nc_io_tools,  only: write_nc
    use, intrinsic :: ieee_arithmetic


    integer,intent(in) :: framein
    integer            :: dframe
    real               :: ke=0., ens=0., ape=0., dedt=0.
    real               :: gen_bci_rate=0., gen_rmf_rate=0., thd_rate=0.
    real               :: bd_rate=0., qd_rate=0., filter_rate=0.
    real               :: filter_rate_ty=0., filter_rate_tx=0., filter_rate_m=0.
    real               :: td_rate=0.,eddy_time=0.,tvary=0.,tvarx=0.,zeta_rms=0.
    real               :: gen_tg_ty=0., gen_tg_tx=0.
    real               :: gen_tg_m=0.,tmoi=0.,gen_div_m=0.,foo=0.
    real               :: meanke=0.
    logical,save       :: called_yet=.false.

    dframe = framein + 1            ! Update diagnostics frame counter

!    call Write_field(time,'diag1_time',dframe) 
                                     ! Track diagnostic-writes
    call write_nc(diag1_time_var_id, time) 
    
    ke = sum(ksqd_*dz*psi*conjg(psi))
    call par_sum(ke)
!    call Write_field(ke,'ke',dframe)
    call write_nc(ke_var_id, ke)

    ens = sum( dz*(q*conjg(q)) )
    call par_sum(ens)
!    call Write_field(ens,'ens',dframe) 
    call write_nc(ens_var_id, ens)

    if (time_varying_mean) then
       meanke = sum(dz*(ubar**2 + vbar**2))/2
!       call Write_field(meanke,'mean_ke',dframe)
       call write_nc(meanke_var_id,meanke)
    end if

    if (trim(filter_type)/='none') then
       filter_rate = get_filter_rate(psi,q,filter,dz,dt)
!       call Write_field(filter_rate,'filter_rate',dframe)
       call write_nc(filter_rate_var_id, filter_rate)
    endif

    if (bot_drag/=0) then
       bd_rate = -2*sum(bot_drag*ksqd_*dz(nz)*conjg(psi(nz,:,:))*psi(nz,:,:))
       call par_sum(bd_rate)
!       call Write_field(bd_rate,'bd_rate',dframe)
       call write_nc(bd_rate_var_id, bd_rate)
    endif
    if (quad_drag/=0) then
       qd_rate = 2*sum(dz(nz)*conjg(psi(nz,:,:))*qdrag)
       call par_sum(qd_rate)
!       call Write_field(qd_rate,'qd_rate',dframe)
       call write_nc(qd_rate_var_id, qd_rate)
    endif
    if (top_drag/=0) then
       td_rate = -2*sum(top_drag*ksqd_*dz(1)*conjg(psi(1,:,:))*psi(1,:,:))
       call par_sum(td_rate)
!       call Write_field(td_rate,'td_rate',dframe)
       call write_nc(td_rate_var_id, td_rate)
    endif
    if (use_forcing) then
       gen_rmf_rate = get_gen_rmf_rate(psi)
!       call Write_field(gen_rmf_rate,'gen_rmf_rate',dframe)
       call write_nc(gen_rmf_rate_var_id, gen_rmf_rate)
    endif
    if (nz>1) then
       ape = sum(F*conjg(psi(2:nz,:,:)-psi(1:nz-1,:,:))* &
            (psi(2:nz,:,:)-psi(1:nz-1,:,:))*(1/drho(1:nz-1)))
       call par_sum(ape)
!       call Write_field(ape,'ape',dframe)
       call write_nc(ape_var_id, ape)

       gen_bci_rate = -2*sum(real(i*(kx_*(dz*conjg(psi)*(shearu*psi-ubar*q)))))&
                      -2*sum(real(i*(ky_*(dz*conjg(psi)*(shearv*psi-vbar*q)))))
       call par_sum(gen_bci_rate)
!       call Write_field(gen_bci_rate,'gen_bci_rate',dframe)
       call write_nc(gen_bci_rate_var_id, gen_bci_rate)
       if (therm_drag/=0) then
!          if (trim(surface_bc)=='surf_buoy') then
!             ! extrapolate psi(1) to surface
!             thd_rate = 2*F*therm_drag*sum(conjg(psi(:,:,1)+b*dz(1)/2)*b)
          thd_rate = 2*sum(therm_drag*conjg(psi)*(q + (ksqd_*psi(1:nz,:,:))))
!          elseif (Fe==0) then
!             thd_rate = 2*F*sum(real(therm_drag* &
!                  (-conjg(psi(1,:,:))*(psi(1,:,:)-psi(2,:,:)) &
!                   -conjg(psi(2,:,:))*(psi(2,:,:)-psi(1,:,:)))))
!          endif
          call par_sum(thd_rate)
!          call Write_field(thd_rate,'thd_rate',dframe)
          call write_nc(thd_rate_var_id, thd_rate)
       endif
    elseif (nz==1.and.F/=0) then     ! checked
       ape = F*sum(psi*conjg(psi))
       call par_sum(ape)
!       call Write_field(ape,'ape',dframe)
       call write_nc(ape_var_id, ape)
       if (therm_drag/=0) then
          thd_rate = -2*sum(therm_drag*F*conjg(psi)*psi)
          call par_sum(thd_rate)
!          call Write_field(thd_rate,'thd_rate',dframe)
          call write_nc(thd_rate_var_id, thd_rate)
       endif
    endif

    if (use_tracer_x) then
       tvarx = sum(dzt*tracer_x*conjg(tracer_x))
       call par_sum(tvarx)
!       call Write_field(tvarx,'tvarx',dframe)
       call write_nc(tvarx_var_id, tvarx)
       if (trim(filter_type_t)/='none') then
          filter_rate_tx = -get_filter_rate(tracer_x,tracer_x,filter_t,dzt,dt)
!          call Write_field(filter_rate_tx,'filter_rate_tx',dframe)
          call write_nc(filter_rate_tx_var_id, filter_rate_tx)
       endif
       if (use_mean_grad_t) then
          gen_tg_tx = - 2*sum(real(dzt*tracer_x*conjg(i*ky_*psi_stir)))
          call par_sum(gen_tg_tx)
!          call Write_field(gen_tg_tx,'gen_tg_tx',dframe)
          call write_nc(gen_tg_tx_var_id, gen_tg_tx)
       endif
    endif
    if (use_tracer_y) then
       tvary = sum(dzt*tracer_y*conjg(tracer_y))
       call par_sum(tvary)
!       call Write_field(tvary,'tvary',dframe)
       call write_nc(tvary_var_id, tvary)
       if (trim(filter_type_t)/='none') then
          filter_rate_ty = -get_filter_rate(tracer_y,tracer_y,filter_t,dzt,dt)
!          call Write_field(filter_rate_ty,'filter_rate_ty',dframe)
          call write_nc(filter_rate_ty_var_id, filter_rate_ty)
       endif
       if (use_mean_grad_t) then
          gen_tg_ty = - 2*sum(real(dzt*tracer_y*conjg(i*kx_*psi_stir)))
          call par_sum(gen_tg_ty)
!          call Write_field(gen_tg_ty,'gen_tg_ty',dframe)
          call write_nc(gen_tg_ty_var_id, gen_tg_ty)
       endif
    endif
       
    zeta_rms = 2*sum( dz*(ksqd_**2*psi*conjg(psi)) )
    call par_sum(zeta_rms)
    eddy_time = 2*pi/sqrt(zeta_rms)
!    call Write_field(eddy_time,'eddy_time',dframe)
    call write_nc(eddy_time_var_id, eddy_time)
    
    ! Check for floating point exceptions 

    if (.not.ieee_is_finite(ke)) call Message('INF or NAN in ke - quitting!',fatal='y')

    ! Write some information to screen and to log file
    call Message('')
    call Message('time step     =',tag=cntr)
    call Message('energy        =',r_tag=ke+ape)
    call Message('enstrophy     =',r_tag=ens)
    if (nz>1.and.(uscale/=0..or.vscale/=0.)) call Message('gen_bci_rate  =',r_tag=gen_bci_rate)
    if (time_varying_mean) call Message('Mean KE       =',r_tag=meanke)
    if (use_forcing)  call Message('gen_rmf_rate  =',r_tag=gen_rmf_rate)
    if (bot_drag/=0)  call Message('botdrag_rate  =',r_tag=bd_rate)
    if (quad_drag/=0) call Message('quaddrag_rate =',r_tag=qd_rate)
    if (therm_drag/=0)call Message('thermdrag_rate=',r_tag=thd_rate)
    if (trim(filter_type)/='none') call Message('filter_rate   =',r_tag=filter_rate) 
    if (use_tracer_x) then
       call Message('txvariance     =',r_tag=tvarx)
       call Message('txvar_gen      =',r_tag=gen_tg_tx)
       call Message('txvar_dissip   =',r_tag=filter_rate_tx)
       if (.not.called_yet) tvarx0 = tvarx
       if (.not.ieee_is_finite(tvarx)) call Message('INF or NAN in tvarx - quitting!',fatal='y')
       tvarx0 = tvarx
    endif
    if (use_tracer_y) then
       call Message('tyvariance     =',r_tag=tvary)
       call Message('tyvar_gen      =',r_tag=gen_tg_ty)
       call Message('tyvar_dissip   =',r_tag=filter_rate_ty)
       if (.not.called_yet) tvary0 = tvary
       if (.not.ieee_is_finite(tvary)) call Message('INF or NAN in tvary - quitting!',fatal='y')
       tvary0 = tvary
    endif

    called_yet = .true.

  end function Get_energetics

  subroutine init_get_spectra()
    !************************************************************************
    ! Prepare spectra.nc file and variables for Get_spectra function
    ! below to use
    !************************************************************************
    use io_tools,    only: Message
    use nc_io_tools, only: create_file, enddef_file, register_variable,      &
                           create_axis_time, create_axis_ktotal,             &
                           create_axis_z, create_axis_custom,                &
                           create_axis_lat
    use qg_params,   only: do_aniso_spectra, uscale, vscale, use_forcing,    &
                           nz, nzt, do_genm_spectra, therm_drag, bot_drag,   &
                           top_drag, quad_drag, filter_type, do_xfer_spectra,&
                           use_tracer_x, use_tracer_y, do_x_avgs,            &
                           use_mean_grad_t
                           

    integer :: spectra_file_id, axis_time_id, axis_ktotal_id, &
               axis_z_id, axis_zt_id, axis_m3_id, axis_lat_id

    spectra_file_id = create_file('spectra.nc')
    axis_time_id    = create_axis_time(spectra_file_id)
    axis_z_id       = create_axis_z(spectra_file_id)
    axis_zt_id      = create_axis_custom(spectra_file_id, 'ztracer', nzt)
    axis_ktotal_id  = create_axis_ktotal(spectra_file_id)
    axis_m3_id      = create_axis_custom(spectra_file_id, 'm3', nz**3)
    axis_lat_id     = create_axis_lat(spectra_file_id)
    
    diagspec_time_var_id = register_variable(spectra_file_id, 'time', &
                           (/axis_time_id/), .false.)
    kes_var_id           = register_variable(spectra_file_id, 'kes',  &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
    
    if (do_aniso_spectra) then
        kesx_var_id      = register_variable(spectra_file_id, 'kesx', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
        kesy_var_id      = register_variable(spectra_file_id, 'kesy', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
    endif

    if (uscale/=0 .OR. vscale/=0) then
        gens_var_id      = register_variable(spectra_file_id, 'gens', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
    endif 

    if (use_forcing) then
        gens_rmf_var_id  = register_variable(spectra_file_id, 'gens_rmf', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
    endif

!for either single layer or multiple layers
    if (bot_drag /= 0) then
        bdms_var_id = register_variable(spectra_file_id, 'bdms', &
                       (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
    endif

    if (quad_drag /= 0) then
        qdms_var_id = register_variable(spectra_file_id, 'qdms', &
                       (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
    endif

    if (trim(filter_type)/='none') then
        filterms_var_id= register_variable(spectra_file_id, 'filterms', &
                       (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
    endif

    multilayer : if (nz>1) then
        kems_var_id      = register_variable(spectra_file_id, 'kems', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)

        if (do_aniso_spectra) then
            kemsx_var_id = register_variable(spectra_file_id, 'kemsx', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
            kemsy_var_id = register_variable(spectra_file_id, 'kemsy', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
        endif

        apems_var_id    = register_variable(spectra_file_id, 'apems', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
        apes_var_id     = register_variable(spectra_file_id, 'apes', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)

        if (do_genm_spectra .AND. uscale/=0) then
            genms_var_id= register_variable(spectra_file_id, 'genms', &
                           (/axis_ktotal_id, axis_m3_id, axis_time_id/), .false.)
        endif

        if (therm_drag /= 0) then
            thdms_var_id= register_variable(spectra_file_id, 'thdms', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
        endif

        if (top_drag /= 0) then
            tdms_var_id = register_variable(spectra_file_id, 'tdms', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
        endif

        if (do_xfer_spectra) then
            xferms_var_id= register_variable(spectra_file_id, 'xferms', &
                           (/axis_ktotal_id, axis_m3_id, axis_time_id/), .false.)
        endif

    elseif (nz==1) then

        if (do_xfer_spectra) then
            energy_xfers_var_id= register_variable(spectra_file_id, 'energy_xfers', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
            enstrophy_xfers_var_id= register_variable(spectra_file_id, 'enstrophy_xfers', &
                           (/axis_ktotal_id, axis_z_id, axis_time_id/), .false.)
        endif

    endif multilayer

    if (use_tracer_x) then
        txvars_var_id   = register_variable(spectra_file_id, 'txvars', &
                           (/axis_ktotal_id, axis_zt_id, axis_time_id/), .false.)
        if (use_mean_grad_t) then
            txfluxs_var_id = register_variable(spectra_file_id, 'txfluxs', &
                           (/axis_ktotal_id, axis_zt_id, axis_time_id/), .false.)
        endif
    endif

    if (use_tracer_y) then
        tyvars_var_id   = register_variable(spectra_file_id, 'tyvars', &
                           (/axis_ktotal_id, axis_zt_id, axis_time_id/), .false.)
        if (use_mean_grad_t) then
            tyfluxs_var_id = register_variable(spectra_file_id, 'tyfluxs', &
                           (/axis_ktotal_id, axis_zt_id, axis_time_id/), .false.)
        endif
    endif

    if (do_x_avgs) then 
        uv_avg_x_var_id = register_variable(spectra_file_id, 'uv_avg_x', &
                           (/axis_z_id, axis_lat_id, axis_time_id/), .false.)
        vq_avg_x_var_id = register_variable(spectra_file_id, 'vq_avg_x', &
                           (/axis_z_id, axis_lat_id, axis_time_id/), .false.)
    endif

    call enddef_file(spectra_file_id)
    call Message('spectra.nc initialized')
  end subroutine init_get_spectra

  function Get_spectra(framein) result(dframe)

    !************************************************************************
    ! Calculate the isotropic horizontal wavenumber vs vertical
    ! wavenumber spectra of modal and layered energetics
    !************************************************************************

    use op_rules,        only: operator(+), operator(-), operator(*)
    use io_tools,        only: Write_field,Message
    use transform_tools, only: Jacob, Spec2grid
    use numerics_lib,    only: Ring_integral, sub2ind
    use strat_tools,     only: Layer2mode
    use qg_arrays,       only: ksqd_,kx_,ky_,kxv,kyv,                                     &
                               psi, q, filter, ug, vg, qdrag                               
    use rmf_forcing,     only: force_o
    use qg_strat_and_shear, only: dz,drho,shearu,shearv,ubar,vbar,                        &
                                  tripint, um, vm, vmode, kz
    use qg_tracers,      only: tracer_x, tracer_y, psi_stir, dzt, filter_t
    use qg_params,       only: i,pi,F,Fe,                                                 &
                               ngrid,kmax,nkx,kx_start,kx_end,ny,nky,ngrid,nz,nzt,        &
                               time,dt,cntr,bot_drag,top_drag,therm_drag,quad_drag,       &
                               uscale,vscale,use_forcing,                                 &
                               filter_type_t,use_mean_grad_t,use_tracer_x,use_tracer_y,   &
                               do_xfer_spectra,do_x_avgs,do_genm_spectra,filter_type,     &
                               do_aniso_spectra, time_varying_mean
    use par_tools,       only: par_sum
    use nc_io_tools,     only: write_nc

    integer,intent(in)                     :: framein
    integer                                :: dframe
    real,dimension(:,:,:),allocatable      :: field, fieldt
    real,dimension(:,:),allocatable        :: spec, spec_m3, xavg, field2d, spect
    real,dimension(:),allocatable          :: uq, vq
    complex,dimension(:,:,:),allocatable   :: psim, jack, qg
    integer                                :: m, j, k, iz
    logical,save                           :: called_yet=.false.

    dframe = framein+1 
    call Write_field(time,'diag2_time',dframe)     ! Track diagnostic-writes
    call write_nc(diagspec_time_var_id, time)

    ! Allocate fields to collect spectral integration results and write to file 
    allocate(spec(1:kmax,1:nz));                                  spec=0.
    allocate(field(nz,kx_start:kx_end,0:kmax));                   field=0.
    allocate(field2d(kx_start:kx_end,nky));                       field2d=0.
    allocate(uq(nz),vq(nz));                                      vq=0.

    if (.not.called_yet) then
       if (do_xfer_spectra) then
          call Message('Will calculate modal transfer spectra - very expensive for nz>2')
       endif
       if (do_genm_spectra) then
          call Message('Will calculate modal generation spectra - very expensive for nz>2.')
          call Message('Not ammended for nonzonal flow yet...')
       endif
       called_yet=.true.
    endif

    ! KE spectra
    field = real(dz*(ksqd_*psi*conjg(psi)))
    spec = Ring_integral(field,kxv,kyv,kmax)
    call par_sum(spec)
    call Write_field(spec,'kes',dframe)
    call write_nc(kes_var_id, spec)

    ! Spectra of energy along kx and ky axes if anisotropy expected
    if (do_aniso_spectra) then         ! KE(kx>0,0)
       spec = 0.
       if (kx_start>0) then
          spec(kx_start:kx_end,:) = transpose(field(:,:,0))
       elseif (kx_start==0) then  ! Skip kx=0 
          spec(1:kx_end,:) = transpose(field(:,1:kx_end,0))
       endif
       call par_sum(spec)
       call Write_field(spec(1:kmax,1:nz),'kesx',dframe)
       call write_nc(kesx_var_id, spec)
       spec=0.
       if (kx_start==0) spec = transpose(field(1:nz,0,1:kmax))
       call par_sum(spec)            ! Adds zeros from all other pes
       call Write_field(spec,'kesy',dframe)
       call write_nc(kesy_var_id, spec)
    endif   

    ! Calculate eddy generation spectra

    spec=0.;  field=0.
    if (uscale/=0.or.vscale/=0) then     ! From mean shear forcing/BCI
       if (uscale/=0) field = -real(2*i*(dz*(shearu*psi-ubar*q)*(kx_*conjg(psi))))
       if (vscale/=0) field = field - real(2*i*(dz*(shearv*psi-vbar*q)*(ky_*conjg(psi))))
       spec = Ring_integral(2*field,kxv,kyv,kmax)
       call Write_field(spec,'gens',dframe)
       call write_nc(gens_var_id, spec)
    endif
    field = 0.; spec = 0.
    if (use_forcing) then        ! From random Markovian forcing
       field(1,:,:) = real(conjg(psi(1,:,:))*force_o)
       if (nz>1) field(2:nz,:,:) = 0.
       spec = Ring_integral(2*field,kxv,kyv,kmax)
       call par_sum(spec)
       call Write_field(spec,'gens_rmf',dframe)
       call write_nc(gens_rmf_var_id, spec)
    endif

    ! Calculate BC diagnostics

    multilayer: if (nz>1) then  

       ! checked
       allocate(psim(1:nz,1:nkx,0:kmax)); psim=0.
       psim = Layer2mode(psi,vmode,dz)
       field = real(ksqd_*psim*conjg(psim)) ! Modal KE
       spec = Ring_integral(field,kxv,kyv,kmax)
       call par_sum(spec)
       call Write_field(spec,'kems',dframe)
       call write_nc(kems_var_id, spec)

       ! Spectra of energy along kx and ky axes if anisotropy expected
       if (do_aniso_spectra) then         ! KE(kx>0,0)
          spec = 0.
          if (kx_start>0) then
             spec(kx_start:kx_end,:) = transpose(field(:,:,0))
          elseif (kx_start==0) then   ! Skip kx=0 
             spec(1:kx_end,:) = transpose(field(:,1:kx_end,0))
          endif
          call par_sum(spec)
          call Write_field(spec(1:kmax,1:nz),'kemsx',dframe)
          call write_nc(kemsx_var_id, spec)
          spec=0.
          if (kx_start==0) spec = transpose(field(1:nz,0,1:kmax))
          call par_sum(spec)            ! Adds zeros from all other pes
          call Write_field(spec,'kemsy',dframe)
          call write_nc(kemsy_var_id, spec)
       endif

       field = real((kz**2)*psim*conjg(psim))    ! Modal APE
       spec = Ring_integral(field,kxv,kyv,kmax)
       call par_sum(spec)
       call Write_field(spec,'apems',dframe)
       call write_nc(apems_var_id, spec)

       field=0.
       field(1:nz-1,:,:) = &                  ! Interface APE
            real(F*conjg(psi(2:nz,:,:)-psi(1:nz-1,:,:))* &
            (psi(2:nz,:,:)-psi(1:nz-1,:,:))*(1/drho(1:nz-1)))
       spec=0.
       spec(:,1:nz-1) = Ring_integral(field(1:nz-1,:,:),kxv,kyv,kmax)
       call par_sum(spec)
       call Write_field(spec,'apes',dframe)
       call write_nc(apes_var_id, spec)

! Need to add vscale and check and calc
       if (do_genm_spectra.and.uscale/=0) then   ! Modal eddy generation
          allocate(spec_m3(kmax,nz**3));           spec_m3=0.
          do m = 1,nz
             do j = 1,nz
                do k = 1,nz
                   field(1,:,:) = 2*real(tripint(j,k,m)*um(j)* &
                        conjg(psim(m,:,:))*i*kx_* &
                        (kz(j)**2-kz(k)**2-ksqd_)*psim(k,:,:))
                   spec_m3(:,sub2ind((/ k,j,m /),nz)) = &
                        Ring_integral(field(1,:,:),kxv,kyv,kmax)
                enddo
             enddo
          enddo
          call par_sum(spec_m3)
          call write_field(spec_m3,'genms',dframe)
          call write_nc(genms_var_id, spec_m3)
          deallocate(spec_m3)
       endif

       if (therm_drag/=0) then   ! Modal thermal drag dissipation 
          do m = 1,nz
             field(m,:,:) = -2*real(therm_drag*sum((vmode(1,:)-vmode(2,:))*psim,1) &
                  *(conjg(psim(m,:,:))*(vmode(1,m)*dz(1)-vmode(2,m)*dz(2))))
          enddo
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'thdms',dframe)
          call write_nc(thdms_var_id, spec)
       endif

       ! checked          
       if (bot_drag/=0) then     ! Modal bottom drag dissipation 
          do m = 1,nz
             field(m,:,:) = -2*real(bot_drag*dz(nz)*ksqd_*sum(vmode(nz,:)*psim,1) &
                  *(vmode(nz,m)*conjg(psim(m,:,:))) )
          enddo
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'bdms',dframe)
          call write_nc(bdms_var_id, spec)
       endif
    
       if (top_drag/=0) then     ! Modal top drag dissipation 
          do m = 1,nz
             field(m,:,:) = -2*real(top_drag*dz(1)*ksqd_*sum(vmode(1,:)*psim,1)  &
                            *(vmode(1,m)*conjg(psim(m,:,:))))
          enddo
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'tdms',dframe)
          call write_nc(tdms_var_id, spec)
       endif

       if (quad_drag/=0) then     ! Modal quadratic drag dissipation 
          do m = 1,nz
             field(m,:,:) = -2*real(dz(nz)*sum(vmode(nz,:)*qdrag,1) &
                  *(vmode(nz,m)*conjg(psim(m,:,:))) )
          enddo
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'qdms',dframe)
          call write_nc(qdms_var_id, spec)
       endif

       ! checked
       if (trim(filter_type)/='none') then   ! Modal filter dissipation
          field2d = 0.
          where (filter/=0) field2d = (filter**(-1)-1)/(2*dt)
          field = -2*real((field2d*(ksqd_+kz**2))*psim*conjg(psim))
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'filterms',dframe)
          call write_nc(filterms_var_id, spec)
       endif

       if (do_xfer_spectra) then         ! Internal transfer terms
          allocate(spec_m3(kmax,nz**3));      spec_m3=0.
          allocate(jack(nz,nkx,0:kmax));      jack=0.
          do m = 1,nz           
             do j = 1,nz
                do k = 1,nz
                   call jacob(psim(j,:,:),-(kz(k)**2+ksqd_)*psim(k,:,:),jack(k,:,:))
                   if (abs(tripint(m,j,k))>1000*epsilon(tripint(m,j,k))) then
                      field(1,:,:) = 2*real(tripint(m,j,k)*conjg(psim(m,:,:))*jack(k,:,:) )
                      spec_m3(:,sub2ind((/ k,j,m /),nz)) = Ring_integral(field(1,:,:),kxv,kyv,kmax)
                   else
                      spec_m3(:,sub2ind((/ k,j,m /),nz)) = 0.
                   endif
                enddo
             enddo
          enddo
          call par_sum(spec_m3)
          call Write_field(spec_m3,'xferms',dframe)
          call write_nc(xferms_var_id, spec_m3)
          deallocate(jack,spec_m3)
       endif
       deallocate(psim)

       ! Area averaged eddy PV flux as function of depth

       uq = 2*sum(sum(real((-i*ky_*psi)*conjg(q)),2),2) 
       call par_sum(uq)
       call write_field(uq,'uq',dframe)
       vq = 2*sum(sum(real((i*kx_*psi)*conjg(q)),2),2)
       call par_sum(vq)
       call write_field(vq,'vq',dframe)

       if (time_varying_mean) then
          call write_field(ubar,'ubar',dframe)
          call write_field(sum( dz*ubar*vq ),'int_ubar_vq_dz',dframe)
          call write_field(vbar,'vbar',dframe)
          call write_field(sum( -dz*vbar*uq ),'int_vbar_uq_dz',dframe)
       end if

    elseif (nz==1) then

       if (F/=0) then             ! Calculate APE spectrum
          field = F*real(conjg(psi)*psi)
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'apes',dframe)
          call write_nc(apes_var_id, spec)
       endif

       if (do_xfer_spectra) then
          allocate(jack(nz,nkx,0:kmax));      jack=0.
          call jacob(psi,q,jack)
          field = -2*real(conjg(psi)*jack)
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'energy_xfers',dframe)
          call write_nc(energy_xfers_var_id, spec)
          field = -2*real(conjg(q)*jack)
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'enstrophy_xfers',dframe)
          call write_nc(enstrophy_xfers_var_id, spec)
       endif

       if (bot_drag/=0) then     ! Bottom drag dissipation 
          field = -2*real(bot_drag*ksqd_*psi*conjg(psi))
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'bds',dframe)
          call write_nc(bdms_var_id, spec)
       endif

       if (trim(filter_type)/='none') then    ! Filter dissipation
          field2d=0.
          where (filter/=0) field2d = (filter**(-1)-1)/(2*dt)
          field = -2*real(field2d*(ksqd_*psi*conjg(psi)))
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'filters',dframe)        
          call write_nc(filterms_var_id, spec)
       endif

       if (quad_drag/=0) then      ! Quadratic drag dissipation
          field = 2*real(dz(nz)*conjg(psi)*qdrag)
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'qds',dframe)
          call write_nc(qdms_var_id, spec)
       endif

    endif multilayer

    ! Write tracer spectra

    if (use_tracer_x.or.use_tracer_y) then
       allocate(spect(kmax,nzt));                      spect=0.
       allocate(fieldt(nzt,kx_start:kx_end,0:kmax));   fieldt=0.;
   
       if (use_tracer_x) then
          fieldt = real(dzt*tracer_x*conjg(tracer_x))     ! <t't'>
          spect = Ring_integral(fieldt,kxv,kyv,kmax)
          call par_sum(spect)
          call Write_field(spect,'txvars',dframe)
          call write_nc(txvars_var_id, spect)
          if (use_mean_grad_t) then                       ! <tx'u'> vs. k
             fieldt = -2*real(dzt*tracer_x*conjg(i*ky_*psi_stir))
             spect = Ring_integral(fieldt,kxv,kyv,kmax)
             call par_sum(spect)
             call Write_field(spect,'txfluxs',dframe)
             call write_nc(txfluxs_var_id, spect)
          endif
       endif
       if (use_tracer_y) then
          fieldt = real(dzt*tracer_y*conjg(tracer_y))     ! <t't'>
          spect = Ring_integral(fieldt,kxv,kyv,kmax)
          call par_sum(spect)
          call Write_field(spect,'tyvars',dframe)
          call write_nc(tyvars_var_id, spect)
          if (use_mean_grad_t) then                       ! <ty'v'> vs. k
             fieldt = 2*real(dzt*tracer_y*conjg(i*kx_*psi_stir))
             spect = Ring_integral(fieldt,kxv,kyv,kmax)
             call par_sum(spect)
             call Write_field(spect,'tyfluxs',dframe)
             call write_nc(tyfluxs_var_id, spect)
          endif
       endif
       deallocate(spect,fieldt)
    endif

    ! Calculate zonal average momentum and PV fluxes

    ! checked
    if (do_x_avgs) then    ! Zonal averages 
       allocate(xavg(nz,ngrid));             xavg=0.
       allocate(qg(nz,ny,ngrid));            qg=0.
       xavg=0.
       xavg(1:nz,1:ngrid)=sum(real(ug)*real(vg),2)/ngrid
       call par_sum(xavg)
       call Write_field(xavg,'uv_avg_x',dframe)
       call write_nc(uv_avg_x_var_id, xavg)
       xavg=0.
       call Spec2grid(q,qg)
       xavg(1:nz,1:ngrid)= sum(real(vg)*real(qg),2)/ngrid ! <v'q'>
       call par_sum(xavg)
       call Write_field(xavg,'vq_avg_x',dframe)
       call write_nc(vq_avg_x_var_id, xavg)
       deallocate(xavg,qg)
    endif

    deallocate(field2d,field,spec,vq)

    call Message('Finished spectral diagnostics')

  end function Get_spectra

   elemental logical function ieee_is_nan(x)
     real,intent(in):: x
     ieee_is_nan = isnan(x)
   end function ieee_is_nan
   elemental logical function ieee_is_finite(x)
     real,intent(in):: x
     ieee_is_finite = .not. (isnan(x) .or. abs(x) > huge(x))
   end function ieee_is_finite

  !*********************************************************************

end module qg_diagnostics
