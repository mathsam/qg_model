program qg_driver 

  ! Stratified or barotropic spectral, homogeneous QG model, using MPI and FFTW
  !
  ! Routines: main, Get_rhs, Get_rhs_t
  !
  ! Dependencies: everything.

  use qg_arrays
  use qg_params
  use qg_strat_and_shear,     only: init_strat_and_shear
  use qg_tracers,             only: init_tracers, step_tracers
  use qg_topo,                only: init_topo
  use qg_init_streamfunction, only: init_streamfunction
  use qg_output,              only: init_counters, write_restarts, write_snapshots
  use qg_diagnostics,         only: get_energetics, get_spectra, enstrophy, get_corrs, &
                                    init_get_energetics, init_get_spectra
  use qg_filter_tools,        only: init_filter
  use rmf_forcing,            only: init_rmf_forcing, markovian
  use transform_tools,        only: init_transform
  use numerics_lib,           only: march
  use io_tools,               only: Message, read_field
  use par_tools,              only: init_par, end_par, par_sync, processor_id
  use op_rules,               only: operator(+), operator(-), operator(*)
  use nc_io_tools,            only: close_all_files
  use qg_sponge,              only: init_qg_sponge

  implicit none

  ! *********** Model initialization *********************

  call Init_par
  call Init_parameters           ! Read command line & namelist, init io_tools
  call Init_transform(kmax, nz, kx_start, kx_end, nkx, y_start, y_end, ny)
  call Init_counters
  call Init_qg_arrays
  filter = Init_filter(filter_type,filter_exp,k_cut,dealiasing,filt_tune,1.)
  call Init_pv_inversion         ! below
  call Init_strat_and_shear
  call par_sync
  if (use_forcing) call Init_rmf_forcing
  if (use_topo)    call Init_topo
  call Init_streamfunction
  if (use_tracer_x.or.use_tracer_y) call Init_tracers  ! requires psi
  call Get_pv                    ! below

  if (sponge_rate > 0.) then
     call init_qg_sponge()
  endif

  call Get_rhs                   ! below
  if (do_energetics) call init_get_energetics()
  if (do_spectra)    call init_get_spectra()

  if (.not.parameters_ok) then
     call Message('The listed errors pertain to values set in your input')
     call Message('namelist file - correct the entries and try again.',fatal='y')
  endif

  ! *********** Main time loop *************************

  call Message('Beginning calculation')

  do cntr = cnt, total_counts         

     call par_sync         ! Make sure we're synched at each time step

     start = (cntr==1)     ! Start flag

     ! Calculate diagnostics, write output
     if (mod(cntr,diag1_step)==0) then
        if (do_energetics) d1frame = Get_energetics(d1frame)
     endif
     if (do_spectra.and.(mod(cntr,diag2_step)==0)) then
        d2frame = Get_spectra(d2frame)
        if (do_corrs) call Get_corrs(d2frame)
     endif
     if (mod(cntr,write_step)==0) then
        frame = Write_snapshots(frame)
     endif
     ! Adapt dt - enstrophy ensured non-zero in initialization.
     if (adapt_dt.and.(mod(cntr,dt_step)==0.or.start)) then
        dt = dt_tune*pi/(kmax*sqrt(max(enstrophy(q),beta,1.)))
     endif

     ! Do the work
     q = filter*March(q,q_o,rhs,dt,robert,call_q)
     psi_o = psi
     call Invert_pv
     call Get_rhs
     if (time_varying_mean) call Step_mean_UV

     if (use_tracer_x.or.use_tracer_y) call step_tracers

     ! Update clock
     time = time + dt 

  enddo  ! End of main time loop

  call Message('Calculation done')
  restart_frame = Write_restarts(restart_frame)
  call close_all_files()  ! close all the NetCDF files
  call End_par
  
!*********************************************************************

contains

  !*********************************************************************

  subroutine Get_rhs

    ! Get physical space velocity and pv gradient terms for
    ! use in calculation of advection (jacobian) and quadratic drag
    
    use qg_strat_and_shear, only: qbarx, qbary, ubar, vbar, dz
    use qg_topo,            only: hb, toposhift
    use transform_tools,    only: grid2spec, spec2grid, ir_pwr, ir_prod
    use qg_sponge,          only: apply_sponge
    complex,dimension(nz,ny,ngrid) :: q_pgrid ! PV on physical grid

    if (.not.linear) then
       call Spec2grid(-i*ky_*psi,ug)
       call Spec2grid(i*kx_*psi,vg) 
       if (use_topo) q(nz,:,:) = q(nz,:,:) + hb
       call Spec2grid(i*kx_*q,qxg) 
       call Spec2grid(i*ky_*q,qyg)
       if (use_topo) q(nz,:,:) = q(nz,:,:) - hb
       if (sponge_rate > 0.) then
           call spec2grid(q, q_pgrid) 
           call Grid2spec(-ir_prod(ug,qxg) - ir_prod(vg,qyg) &
                          + apply_sponge(q_pgrid), rhs)
       else
           call Grid2spec(-ir_prod(ug,qxg) - ir_prod(vg,qyg),rhs)
       endif
         
      if (quad_drag/=0) then  
         unormbg = ir_pwr(ir_pwr(ug(nz,:,:),2.) + ir_pwr(vg(nz,:,:),2.),.5)
         call Grid2spec( ir_prod(unormbg,&
               (ug(nz,:,:)*sin(qd_angle)+vg(nz,:,:)*cos(qd_angle))), qdrag  )
         call Grid2spec( ir_prod(unormbg, &
               (ug(nz,:,:)*cos(qd_angle)-vg(nz,:,:)*sin(qd_angle))), qdrag2 )
         qdrag = quad_drag*filter*(i*kx_*qdrag - i*ky_*qdrag2)
         rhs(nz,:,:) = rhs(nz,:,:) - qdrag
      endif
    endif

    if (any(qbarx/=0)) rhs = rhs - i*(ky_*(-qbarx*psi + vbar*q))
    if (any(qbary/=0)) rhs = rhs - i*(kx_*( qbary*psi + ubar*q))
    if (use_topo) then
       if (time_varying_mean) toposhift = - i*ubar(nz)*(kx_*hb)- i*vbar(nz)*(ky_*hb)
       rhs(nz,:,:) = rhs(nz,:,:) + toposhift
    end if
    if (top_drag/=0)   rhs(1,:,:) = rhs(1,:,:) + top_drag*ksqd_*psi_o(1,:,:)
    if (bot_drag/=0) then
       if (limit_bot_drag) then
          where (ksqd_<kb_max**2)
             rhs(nz,:,:) = rhs(nz,:,:) + bot_drag*ksqd_*psi_o(nz,:,:)
          end where
       else
          rhs(nz,:,:) = rhs(nz,:,:) + bot_drag*ksqd_*psi_o(nz,:,:)
       endif
    endif
    if (use_forcing) rhs(z_force,:,:)  = rhs(z_force,:,:) +  &
         Markovian(kf_min,kf_max,forc_coef,forc_corr,norm_forcing,norm_diss, &
         filter,filter_type,psi,q,dz,dt)

    if (therm_drag/=0) then
       if (nz>1) then
          rhs = rhs - therm_drag*(q + (ksqd_*psi(1:nz,:,:)))
       elseif (nz==1.and.F/=0) then
          rhs = rhs + therm_drag*F*psi_o
       endif
    endif

    rhs = filter*rhs

  end subroutine Get_rhs

  !*******************************************************************

  subroutine Step_mean_UV

    !**************************************************************
    ! Forward step mean velocity via
    !  dU/dt =  <vq> - T_restore^(-1) * (U - U_o)
    !  dV/dt = -<uq> - T_restore^(-1) * (V - V_o)
    !*************************************************************

    use par_tools,              only: par_sum
    use qg_strat_and_shear,     only: ubar, ubar_o, ubar_target, &
                                      vbar, vbar_o, vbar_target, &
                                      qbarx, qbary, shearu, shearv, &
                                      dz, psiq, um, vm, vmode
    use qg_topo,                only: hb
    use numerics_lib,           only: tri2mat

    real, dimension(nz)             :: uq, vq, rhsU, rhsV

    if (use_topo) q(nz,:,:) = q(nz,:,:) + hb
    vq = 2*sum(sum(real((i*kx_*psi)*conjg(q)),2),2)
    call par_sum(vq)
    
    rhsU = vq - inv_restore_time*(ubar-ubar_target)
    ubar = March(ubar,ubar_o,rhsU,dt,robert,call_U)
    shearu = matmul(tri2mat(psiq),ubar)
    um = matmul(transpose(vmode),dz*ubar)    ! Project shear onto modes
    qbary = -shearu + beta  

    uq = 2*sum(sum(real((-i*ky_*psi)*conjg(q)),2),2)
    call par_sum(uq)

    rhsV = -uq - inv_restore_time*(vbar-vbar_target)
    vbar = March(vbar,vbar_o,rhsV,dt,robert,call_V)
    shearv = matmul(tri2mat(psiq),vbar)
    vm = matmul(transpose(vmode),dz*vbar)    ! Project shear onto modes
    qbarx = shearv
    if (use_topo) q(nz,:,:) = q(nz,:,:) - hb

  end subroutine Step_mean_UV

  !*******************************************************************

  subroutine Get_pv

    !**************************************************************
    ! Calculate PV field from streamfunction
    !**************************************************************

    use qg_strat_and_shear, only: psiq

    integer   :: top

    select case (trim(surface_bc))
    case ('rigid_lid');    top = 1
    case ('periodic');     top = nz
    end select

    if (nz==2) then
       q(1,:,:) = &
            psiq(1,-1)*psi(top,:,:) &
            + psiq(1,0)*psi(1,:,:) &
            + psiq(1,1)*psi(2,:,:) 
       q(nz,:,:) = &
            psiq(nz,-1)*psi(nz-1,:,:) &
            + psiq(nz,0)*psi(nz,:,:) &
            + psiq(nz,1)*psi(1,:,:)
       
    elseif (nz>1) then    
       q(1,:,:) = &
            psiq(1,-1)*psi(top,:,:) &
            + psiq(1,0)*psi(1,:,:) &
            + psiq(1,1)*psi(2,:,:) 
       q(2:nz-1,:,:) = &
            psiq(2:nz-1,-1)*psi(1:nz-2,:,:) &
            + psiq(2:nz-1,0)*psi(2:nz-1,:,:) &
            + psiq(2:nz-1,1)*psi(3:nz,:,:)
       q(nz,:,:) = &
            psiq(nz,-1)*psi(nz-1,:,:) &
            + psiq(nz,0)*psi(nz,:,:) &
            + psiq(nz,1)*psi(1,:,:)
    else
       q = -F*psi
    endif

    q = -ksqd_*psi + q
    
  end subroutine Get_pv

  !*******************************************************************

  subroutine init_pv_inversion

    ! Initialize PV inversion

    ! First check upper BC

    if (nz>1) then
       select case (trim(surface_bc))
       case ('rigid_lid')
          call Message('Rigid lid surface BC selected')
!       case ('surf_buoy')
!          call Message('Surface B term selected - DOES NOTHING')
!          surf_buoy = .true.
       case ('periodic')
          call Message('Periodic vertical BC selected')
          if (trim(strat_type)/='linear') &
               call Message('Warning: Non-uniform stratification is not& 
               & well-posed with periodic vertical BC')
       case default
          call Message('Info:  No surface BC selected -- selecting rigid_lid &
               &by default.')
          call Message('Legal surface BC choices: &
               &rigid_lid|surf_buoy|periodic')
          surface_bc = 'rigid_lid'
       end select
    endif
    
    ! Get # of true elements in de-aliasing mask and set up index maps so 
    ! that we can choose only to invert dynamic part of psi

    nmask = sum(ceiling(filter))   ! Get non-zero total of filter 
    ! print*, 'my pe', processor_id, 'nmask ', nmask

    allocate(lin2kx(nmask),lin2ky(nmask));     lin2kx=0; lin2ky=0

    lin2kx = pack(kx_,MASK=(filter>0))
    lin2ky = pack(ky_,MASK=(filter>0))

    call par_sync

    call Message('PV inversion initialized')

  end subroutine init_pv_inversion

  !*******************************************************************

  subroutine Invert_pv

    !**************************************************************
    ! Invert PV (q) in manner depending on surface_bc, which can be 
    ! either 'rigid_lid', 'free_surface', or 'periodic'.  These
    ! types are checked in qg_params/Check_parameters right now.
    !**************************************************************

    use numerics_lib,       only: tridiag, tridiag_cyc
    use qg_strat_and_shear, only: psiq, dz, twolayfac

    real,dimension(1:nz,-1:1)                       :: psiq_m_ksqd
    complex,dimension(nz)                           :: qvec
    integer                                         :: kx, ky, kl, knum, iz

    psi = 0.
    psiq_m_ksqd = psiq

    if (nz==1) then       ! Barotropic model inversion

       where (ksqd_/=0)
          psi(1,:,:) = -q(1,:,:)/(F+ksqd_)
       endwhere

    elseif(nz==2.and.(trim(surface_bc)=='rigid_lid')) then
       
       ! 2-layer inversion for rigid lid case - see init_strat for denom
       psi(1,:,:)=-((F/dz(2) + ksqd_)*q(1,:,:) + (F/dz(1))*q(2,:,:))*twolayfac
       psi(2,:,:)=-((F/dz(1) + ksqd_)*q(2,:,:) + (F/dz(2))*q(1,:,:))*twolayfac

    else    ! Multi-layer inversion

       select case (trim(surface_bc))
       case ('rigid_lid')

          do kl = 1,nmask

             kx = lin2kx(kl)
             ky = lin2ky(kl)
             psiq_m_ksqd(:,0) = psiq(:,0) - ksqd_(kx,ky) 
             psi(:,kx,ky) = tridiag(q(:,kx,ky),psiq_m_ksqd)

          enddo

       case ('periodic')

          do kl = 1,nmask
             kx = lin2kx(kl)
             ky = lin2ky(kl)
             psiq_m_ksqd(:,0) = psiq(:,0) - ksqd_(kx,ky) 
             psi(:,kx,ky) = tridiag_cyc(q(:,kx,ky),psiq_m_ksqd)
          enddo          

       end select

    endif

  end subroutine Invert_pv

  !******************* END OF PROGRAM ********************************* 

end program qg_driver
