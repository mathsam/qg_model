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
  save

  public :: Get_energetics, Get_spectra, energy, enstrophy, get_corrs

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
                           time,cntr,uscale,vscale,use_forcing,                     &
                           use_tracer_x,use_tracer_y,                        &
                           quad_drag,surface_bc,filter_type,filter_type_t,   &
                           time_varying_mean,use_mean_grad_t, do_energetics
    use qg_filter_tools, only: get_filter_rate
    use par_tools,    only: par_sum
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

    call Write_field(time,'diag1_time',dframe) 
                                     ! Track diagnostic-writes
    
    ke = sum(ksqd_*dz*psi*conjg(psi))
    call par_sum(ke)
    call Write_field(ke,'ke',dframe)

    ens = sum( dz*(q*conjg(q)) )
    call par_sum(ens)
    call Write_field(ens,'ens',dframe) 

    if (time_varying_mean) then
       meanke = sum(dz*(ubar**2 + vbar**2))/2
       call Write_field(meanke,'mean_ke',dframe)
    end if

    if (trim(filter_type)/='none') then
       filter_rate = get_filter_rate(psi,q,filter,dz,dt)
       call Write_field(filter_rate,'filter_rate',dframe)
    endif

    if (bot_drag/=0) then
       bd_rate = -2*sum(bot_drag*ksqd_*dz(nz)*conjg(psi(nz,:,:))*psi(nz,:,:))
       call par_sum(bd_rate)
       call Write_field(bd_rate,'bd_rate',dframe)
    endif
    if (quad_drag/=0) then
       qd_rate = 2*sum(dz(nz)*conjg(psi(nz,:,:))*qdrag)
       call par_sum(qd_rate)
       call Write_field(qd_rate,'qd_rate',dframe)
    endif
    if (top_drag/=0) then
       td_rate = -2*sum(top_drag*ksqd_*dz(1)*conjg(psi(1,:,:))*psi(1,:,:))
       call par_sum(td_rate)
       call Write_field(td_rate,'td_rate',dframe)
    endif
    if (use_forcing) then
       gen_rmf_rate = get_gen_rmf_rate(psi)
       call Write_field(gen_rmf_rate,'gen_rmf_rate',dframe)
    endif
    if (nz>1) then
       ape = sum(F*conjg(psi(2:nz,:,:)-psi(1:nz-1,:,:))* &
            (psi(2:nz,:,:)-psi(1:nz-1,:,:))*(1/drho(1:nz-1)))
       call par_sum(ape)
       call Write_field(ape,'ape',dframe)

       gen_bci_rate = -2*sum(real(i*(kx_*(dz*conjg(psi)*(shearu*psi-ubar*q)))))&
                      -2*sum(real(i*(ky_*(dz*conjg(psi)*(shearv*psi-vbar*q)))))
       call par_sum(gen_bci_rate)
       call Write_field(gen_bci_rate,'gen_bci_rate',dframe)
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
          call Write_field(thd_rate,'thd_rate',dframe)
       endif
    elseif (nz==1.and.F/=0) then     ! checked
       ape = F*sum(psi*conjg(psi))
       call par_sum(ape)
       call Write_field(ape,'ape',dframe)
       if (therm_drag/=0) then
          thd_rate = -2*sum(therm_drag*F*conjg(psi)*psi)
          call par_sum(thd_rate)
          call Write_field(thd_rate,'thd_rate',dframe)
       endif
    endif

    if (use_tracer_x) then
       tvarx = sum(dzt*tracer_x*conjg(tracer_x))
       call par_sum(tvarx)
       call Write_field(tvarx,'tvarx',dframe)
       if (trim(filter_type_t)/='none') then
          filter_rate_tx = -get_filter_rate(tracer_x,tracer_x,filter_t,dzt,dt)
          call Write_field(filter_rate_tx,'filter_rate_tx',dframe)
       endif
       if (use_mean_grad_t) then
          gen_tg_tx = - 2*sum(real(dzt*tracer_x*conjg(i*ky_*psi_stir)))
          call par_sum(gen_tg_tx)
          call Write_field(gen_tg_tx,'gen_tg_tx',dframe)
       endif
    endif
    if (use_tracer_y) then
       tvary = sum(dzt*tracer_y*conjg(tracer_y))
       call par_sum(tvary)
       call Write_field(tvary,'tvary',dframe)
       if (trim(filter_type_t)/='none') then
          filter_rate_ty = -get_filter_rate(tracer_y,tracer_y,filter_t,dzt,dt)
          call Write_field(filter_rate_ty,'filter_rate_ty',dframe)
       endif
       if (use_mean_grad_t) then
          gen_tg_ty = - 2*sum(real(dzt*tracer_y*conjg(i*kx_*psi_stir)))
          call par_sum(gen_tg_ty)
          call Write_field(gen_tg_ty,'gen_tg_ty',dframe)
       endif
    endif
       
    zeta_rms = 2*sum( dz*(ksqd_**2*psi*conjg(psi)) )
    call par_sum(zeta_rms)
    eddy_time = 2*pi/sqrt(zeta_rms)
    call Write_field(eddy_time,'eddy_time',dframe)
    
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

  !*************************************************************************

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
       spec=0.
       if (kx_start==0) spec = transpose(field(1:nz,0,1:kmax))
       call par_sum(spec)            ! Adds zeros from all other pes
       call Write_field(spec,'kesy',dframe)
    endif   

    ! Calculate eddy generation spectra

    spec=0.;  field=0.
    if (uscale/=0.or.vscale/=0) then     ! From mean shear forcing/BCI
       if (uscale/=0) field = -real(2*i*(dz*(shearu*psi-ubar*q)*(kx_*conjg(psi))))
       if (vscale/=0) field = field - real(2*i*(dz*(shearv*psi-vbar*q)*(ky_*conjg(psi))))
       spec = Ring_integral(2*field,kxv,kyv,kmax)
       call Write_field(spec,'gens',dframe)
    endif
    field = 0.; spec = 0.
    if (use_forcing) then        ! From random Markovian forcing
       field(1,:,:) = real(conjg(psi(1,:,:))*force_o)
       if (nz>1) field(2:nz,:,:) = 0.
       spec = Ring_integral(2*field,kxv,kyv,kmax)
       call par_sum(spec)
       call Write_field(spec,'gens_rmf',dframe)
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
          spec=0.
          if (kx_start==0) spec = transpose(field(1:nz,0,1:kmax))
          call par_sum(spec)            ! Adds zeros from all other pes
          call Write_field(spec,'kemsy',dframe)
       endif

       field = real((kz**2)*psim*conjg(psim))    ! Modal APE
       spec = Ring_integral(field,kxv,kyv,kmax)
       call par_sum(spec)
       call Write_field(spec,'apems',dframe)

       field=0.
       field(1:nz-1,:,:) = &                  ! Interface APE
            real(F*conjg(psi(2:nz,:,:)-psi(1:nz-1,:,:))* &
            (psi(2:nz,:,:)-psi(1:nz-1,:,:))*(1/drho(1:nz-1)))
       spec=0.
       spec(:,1:nz-1) = Ring_integral(field(1:nz-1,:,:),kxv,kyv,kmax)
       call par_sum(spec)
       call Write_field(spec,'apes',dframe)

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
       endif
    
       if (top_drag/=0) then     ! Modal top drag dissipation 
          do m = 1,nz
             field(m,:,:) = -2*real(top_drag*dz(1)*ksqd_*sum(vmode(1,:)*psim,1)  &
                            *(vmode(1,m)*conjg(psim(m,:,:))))
          enddo
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'tdms',dframe)
       endif

       if (quad_drag/=0) then     ! Modal quadratic drag dissipation 
          do m = 1,nz
             field(m,:,:) = -2*real(dz(nz)*sum(vmode(nz,:)*qdrag,1) &
                  *(vmode(nz,m)*conjg(psim(m,:,:))) )
          enddo
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'qdms',dframe)
       endif

       ! checked
       if (trim(filter_type)/='none') then   ! Modal filter dissipation
          field2d = 0.
          where (filter/=0) field2d = (filter**(-1)-1)/(2*dt)
          field = -2*real((field2d*(ksqd_+kz**2))*psim*conjg(psim))
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'filterms',dframe)
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
       endif

       if (do_xfer_spectra) then
          allocate(jack(nz,nkx,0:kmax));      jack=0.
          call jacob(psi,q,jack)
          field = -2*real(conjg(psi)*jack)
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'energy_xfers',dframe)
          field = -2*real(conjg(q)*jack)
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'enstrophy_xfers',dframe)
       endif

       if (bot_drag/=0) then     ! Bottom drag dissipation 
          field = -2*real(bot_drag*ksqd_*psi*conjg(psi))
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'bds',dframe)
       endif

       if (trim(filter_type)/='none') then    ! Filter dissipation
          field2d=0.
          where (filter/=0) field2d = (filter**(-1)-1)/(2*dt)
          field = -2*real(field2d*(ksqd_*psi*conjg(psi)))
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'filters',dframe)        
       endif

       if (quad_drag/=0) then      ! Quadratic drag dissipation
          field = 2*real(dz(nz)*conjg(psi)*qdrag)
          spec = Ring_integral(field,kxv,kyv,kmax)
          call par_sum(spec)
          call Write_field(spec,'qds',dframe)
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
          if (use_mean_grad_t) then                       ! <tx'u'> vs. k
             fieldt = -2*real(dzt*tracer_x*conjg(i*ky_*psi_stir))
             spect = Ring_integral(fieldt,kxv,kyv,kmax)
             call par_sum(spect)
             call Write_field(spect,'txfluxs',dframe)
          endif
       endif
       if (use_tracer_y) then
          fieldt = real(dzt*tracer_y*conjg(tracer_y))     ! <t't'>
          spect = Ring_integral(fieldt,kxv,kyv,kmax)
          call par_sum(spect)
          call Write_field(spect,'tyvars',dframe)
          if (use_mean_grad_t) then                       ! <ty'v'> vs. k
             fieldt = 2*real(dzt*tracer_y*conjg(i*kx_*psi_stir))
             spect = Ring_integral(fieldt,kxv,kyv,kmax)
             call par_sum(spect)
             call Write_field(spect,'tyfluxs',dframe)
          endif
       endif
       deallocate(spect,fieldt)
    endif

    ! Calculate zonal average momentum and PV fluxes

    ! checked
    if (do_x_avgs) then    ! Zonal averages 
       allocate(xavg(ngrid,nz));             xavg=0.
       allocate(qg(nz,ny,ngrid));            qg=0.
       xavg=0.
       xavg(1:nz,1:ngrid)=sum(real(ug)*real(vg),2)/ngrid
       call par_sum(xavg)
       call Write_field(xavg,'uv_avg_x',dframe)
       xavg=0.
       call Spec2grid(q,qg)
       xavg(1:nz,1:ngrid)= sum(real(vg)*real(qg),2)/ngrid ! <v'q'>
       call par_sum(xavg)
       call Write_field(xavg,'vq_avg_x',dframe)
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
