module qg_strat_and_shear

  !************************************************************************
  ! Stratification and shear initialization and storage.
  !
  ! Routines:      init_strat_and_shear
  !
  ! Dependencies:  io_tools, op_rules, par_tools, strat_tools, qg_params, 
  !                qg_arrays
  !
  !************************************************************************

  implicit none
  private
  save

  real,   dimension(:),    allocatable   :: dz, rho, drho, kz
  real,   dimension(:),    allocatable   :: shearu, shearv, qbarx, qbary
  real,   dimension(:),    allocatable   :: ubar, vbar, um, vm
  real,   dimension(:),    allocatable   :: ubar_o, ubar_target
  real,   dimension(:),    allocatable   :: vbar_o, vbar_target
  real,   dimension(:,:),  allocatable   :: vmode, psiq, twolayfac
  real,   dimension(:,:,:),allocatable   :: tripint

  public ::  vmode, psiq, twolayfac, dz, rho, drho, kz, shearu, shearv,   &
             qbarx, qbary, ubar, vbar, um, vm, ubar_o, ubar_target, &
             vbar_o, vbar_target, tripint, init_strat_and_shear

contains

  !************************************************************************

  subroutine init_strat_and_shear

    use io_tools,     only: message, read_field, write_field
    use qg_arrays,    only: ksqd_
    use qg_params,    only: nz, beta, umode, cr, io_root, F,         &
                            uscale,ubar_type,ubar_in_file,ubar_file, &
                            vscale,vbar_type,vbar_in_file,vbar_file, &
                            kmax, kx_start, kx_end, parameters_ok,   &
                            inv_restore_time, time_varying_mean,     &
                            restarting, d2frame
    use numerics_lib, only: tri2mat
    use par_tools,    only: par_bcast, processor_id

    ! Allocate strat arrays

    allocate(dz(nz));                         dz    = 1.
    allocate(rho(nz));                        rho   = 1.
    allocate(drho(1:nz-1));                   drho  = 1.
    allocate(psiq(1:nz,-1:1));                psiq  = 0.
    allocate(vmode(1:nz,1:nz));               vmode = 1.
    allocate(kz(nz));                         kz    = 0.
    allocate(tripint(nz,nz,nz));              tripint = 0.
    allocate(ubar(nz));                       ubar   = 0.
    allocate(vbar(nz));                       vbar   = 0.
    allocate(shearu(nz));                     shearu = 0.
    allocate(shearv(nz));                     shearv = 0.
    allocate(um(nz));                         um     = 0.
    allocate(vm(nz));                         vm     = 0.
    allocate(qbarx(nz));                      qbarx  = 0.
    allocate(qbary(nz));                      qbary  = 0.

    if (inv_restore_time/=0) then
       call Message('You have selected to use a time-varying mean velocity&
            & with inverse restore time =',r_tag=inv_restore_time)
       allocate(ubar_o(nz),ubar_target(nz),vbar_o(nz),vbar_target(nz))
       ubar_o = 0.; ubar_target = 0.; vbar_o = 0.; vbar_target = 0.
       time_varying_mean=.true.
    endif

    if (uscale==cr) then
       call Message('Info: uscale not set - setting to 0')
       uscale = 0.
    elseif (uscale/=cr) then
       call Message('Zonal mean velocity: uscale =',r_tag=uscale)
    endif
    
    if (vscale==cr) then
       call Message('Info: vscale not set - setting to 0')
       vscale = 0.
    elseif (vscale/=cr) then
       call Message('Meridional mean velocity: uscale =',r_tag=vscale)
    endif

    if (nz==1) then

       ubar = uscale
       vbar = vscale

    elseif (nz>1) then

       if (processor_id==io_root) then  ! Do all vertical calcs on io_root

          call init_strat  ! Sets dz, rho, drho, psiq, vmode, kz, tripint
          
          if (uscale/=0.and.ubar_type/='none') then
             if (restarting.and.inv_restore_time/=0) then
                call Read_field(ubar,ubar_file,frame=d2frame)
                call Read_field(ubar_target,ubar_file,frame=1)
             else
                ubar = Init_ubar(ubar_type,ubar_in_file,uscale,umode,ubar_file)
                ubar_target = ubar
             endif
             shearu = matmul(tri2mat(psiq),ubar)
             um = matmul(transpose(vmode),dz*ubar)    ! Project shear onto modes
          endif
          
          if (vscale/=0.and.vbar_type/='none') then
             if (restarting.and.inv_restore_time/=0) then
                call Read_field(vbar,vbar_file,frame=d2frame)
                call Read_field(vbar_target,vbar_file,frame=1)
             else
                vbar = Init_ubar(vbar_type,vbar_in_file,vscale,umode,vbar_file)
                vbar_target = vbar
             endif 
             shearv = matmul(tri2mat(psiq),vbar)
             vm = matmul(transpose(vmode),dz*vbar)    ! Project shear onto modes
          endif
          
       endif

       ! Now broadcast the results
       call par_bcast(dz,io_root)
       call par_bcast(rho,io_root)
       call par_bcast(drho,io_root)
       call par_bcast(psiq,io_root)
       call par_bcast(vmode,io_root)
       call par_bcast(kz,io_root)
       call par_bcast(tripint,io_root)
       call par_bcast(ubar,io_root)
       call par_bcast(vbar,io_root)
       call par_bcast(shearu,io_root)
       call par_bcast(shearv,io_root)
       call par_bcast(um,io_root)
       call par_bcast(vm,io_root)
       if (time_varying_mean) then
          call par_bcast(ubar_target,io_root)
          call par_bcast(vbar_target,io_root)
       end if

       if (nz==2) then
          ! Set up array for nz=2 analytic PV inversion
          allocate(twolayfac(kx_start:kx_end,0:kmax))
          twolayfac = (ksqd_**2 + ksqd_*F*(1/dz(1)+1/dz(2)))**(-1)
       endif

    endif

    qbary = -shearu + beta                      ! Mean PV gradient
    qbarx = shearv

    call Message('Stratification and shear initialized')

  end subroutine init_strat_and_shear

  !************************************************************************

  subroutine Init_strat

    !************************************************************************
    ! If its a multi-level run, read in or create density profile (rho),
    ! layer thicknesses (dz).  Also set inversion matrix psiq,
    ! vertical modes (vmode), defn wavenumbers (kz) and triple interaction
    ! coefficients (tripint matrix)
    !
    ! Legal values of 'strat_type' are:
    !
    !  linear   : uniform stratification
    !
    !  twolayer : two unequal layers - in this case 'deltc' = dz(1)
    !             and vmode, kz and tripint calculated analytically
    !
    !  exp      : exponential stratification (see qg_strat_tools/Make_strat)
    !             with decay scale 'deltc'
    !
    !  stc      : tanh profile for stratification ( " ") with decay scale
    !             'deltc'
    !
    !  read     : read from 'dz_file' and 'rho_file'
    !************************************************************************

    use io_tools,    only: Message, Write_field, Read_field
    use qg_params,   only: dz_file,rho_file,F,Fe,nz,deltc,hf,drt,drb,        &
                           rho_slope,psiq_file,vmode_file,kz_file,           &
                           tripint_file,dz_in_file,rho_in_file,              &
                           parameters_ok,cr,ci,tripint_in_file,read_tripint, &
                           kmax, kx_start, kx_end, strat_type,surface_bc
    use strat_tools

    real,dimension(1:nz)    :: zl

    ! Create or read density (rho) and layer thickness (dz) profiles

    select case (trim(strat_type))
    case ('linear')
          
       call Message('Linear stratification selected')

       deltc = 0.
       dz = 1./float(nz)
       zl = Get_z(dz)
       rho = 1+((zl-zl(1))/(zl(nz)-zl(1))-.5)

    case ('twolayer')
          
       call Message('Two (unequal) layer stratification selected')
       if (nz/=2) then
          call Message('Error: this setting invalid for nz /= 2')
          Parameters_ok=.false.
       endif
       if (deltc==cr) then
          call Message('Error: deltc must be set to give top layer thickness')
          Parameters_ok=.false.
       else
          call Message('Upper layer thickness ratio deltc =',r_tag=deltc)
       endif

       dz(1) = deltc; dz(2) = 1-deltc
       zl = Get_z(dz)
       rho = 1+((zl-zl(1))/(zl(nz)-zl(1))-.5)

    case('exp')
          
       call Message('Exponential stratification selected')
       if (nz==2) call Message('Caution: consider using strat_type = twolayer &
            &which interprets deltc as upper layer thickness fraction')
       if (deltc==cr) then
          call Message('Error: deltc must be set to make EXP rho')
          Parameters_ok=.false.
       else
          call Message('Scale depth deltc =',r_tag=deltc)
       endif
          
       call Make_strat(dz,rho,deltc,strat_type,rho_slope)

    case ('stc')
          
       call Message('Surface intensified stratification selected')
       if (nz==2) call Message('Caution: consider using strat_type = twolayer &
            &which interprets deltc as upper layer thickness fraction')
       if (deltc==cr) then
          call Message('Error: deltc must be set to make STC rho')
          Parameters_ok=.false.
       else
          call Message('Scale depth deltc =',r_tag=deltc)
       endif

       call Make_strat(dz,rho,deltc,strat_type,rho_slope)

    case ('read')
          
       if (trim(rho_file)=='') then
          call Message('Error: no input file for rho given')
          Parameters_ok=.false.
       endif
       if (trim(dz_file)=='') then
          call Message('Error: no input file for dz given')
          Parameters_ok=.false.
       endif
       call Message('Will read density from: '//trim(rho_in_file))
       call Message('Will read layer thicknesses from: '//trim(dz_in_file))

       call Read_field(dz,dz_in_file,exclude_dd=1)
       call Read_field(rho,rho_in_file,exclude_dd=1)

       dz = dz/sum(dz) 
       rho = rho/sum(rho*dz) 

    case default
          
       call Message('Error: must select strat_type = &
            &linear|stc|exp|read &
            &with nz>1 - yours is:'//trim(strat_type))
       Parameters_ok=.false.
       
    end select

    drho = rho(2:nz)-rho(1:nz-1)       ! Get delta_rho
    drho = drho/(sum(drho)/size(drho))         ! Normalize drho

    call Get_vmodes(dz,drho,F,Fe,kz,vmode,surface_bc)

    call Message('First internal deformation wavenumber is: ', r_tag=kz(2))

    if (read_tripint) then
       call Read_field(tripint,tripint_in_file,exclude_dd=1)
    else
       tripint = Get_tripint(dz,rho,surface_bc,hf,drt,drb)
    endif

    psiq = strat_params(dz,drho,F,Fe,surface_bc)

    ! Save strat fields for post-processing
    call Write_field(psiq,psiq_file)
    call Write_field(vmode,vmode_file)
    call Write_field(kz,kz_file)            
    call Write_field(dz,dz_file)
    call Write_field(rho,rho_file)
    call Write_field(tripint,tripint_file)

  end subroutine Init_strat

  !************************************************************************

  function Init_ubar(ubar_type,ubar_in_file,uscale,umode,ubar_file) result(ubar)

    !************************************************************************
    ! Make mean zonal velocity profile. Legal 'ubar_type' values are
    !
    !  stc     : half gaussian profile with decay scale 'delu'
    !
    !  exp     : exponential profile with decay scale 'delu'
    !  
    !  linear  : linearly varying profile (Eady profile)
    !
    !  modal   : ubar projects exactly onto mode 'umode'
    !
    !  read    : read from file 'ubar_in_file'
    !
    ! Result of any init type is multiplied by uscale.  All init types
    ! except 'read' and 'modal' are normalized by Normu.
    !************************************************************************

    use io_tools,    only: Message, Write_field, Read_field
    use qg_params,   only: nz,delu,parameters_ok,cr,ci,normalize_ubar
    use strat_tools, only: Get_z

    character(*),intent(in)  :: ubar_type,ubar_in_file,ubar_file
    integer,intent(in)       :: umode
    real,intent(inout)       :: uscale
    real,dimension(1:nz)     :: ubar
    real,dimension(1:nz)     :: zl

    zl = Get_z(dz)
    select case (trim(ubar_type))
    case ('stc')

       call Message('Surface intensified mean zonal velocity selected')
       if (delu==ci) then
          call Message('Error: delu must be set for making ubar profile')
          Parameters_ok=.false.
       elseif (delu<0) then
          call Message('Error: require delu>=0 - yours is:',r_tag=delu)
          Parameters_ok=.false.
       endif

       ubar = Normu(exp(-(zl/delu)**2),dz)

    case ('exp')

       call Message('Exponential mean zonal velocity selected')
       if (delu==ci) then
          call Message('Error: delu must be set for making ubar profile')
          Parameters_ok=.false.
       elseif (delu<0) then
          call Message('Error: require delu>=0 - yours is:',r_tag=delu)
          Parameters_ok=.false.
       endif

       ubar = Normu(exp(-zl/delu),dz)

    case ('linear')

       call Message('Linear mean zonal velocity selected')
       delu = 0.

       ubar = Normu(2*(1-(zl-zl(1))/(zl(nz)-zl(1)))-1,dz)

    case ('modal') 

       call Message('Modal mean zonal velocity selected')
       if (umode==ci) then
          call Message('Error: umode must be set for modal ubar')
          Parameters_ok=.false.
       elseif ((umode<0).or.(umode>nz-1)) then
          call Message('Error: require 0<umode<=nz; yours is:',tag=umode)
          Parameters_ok=.false.
       endif

       ubar = vmode(:,umode+1)

    case ('read')

       if (trim(ubar_file)=='') then
          call Message('Error: no input file for ubar given')
          Parameters_ok=.false.
       endif
       call Message('Will read mean U from: '//trim(ubar_in_file))

       call Read_field(ubar,ubar_in_file,exclude_dd=1)

       if (normalize_ubar) ubar = Normu(ubar,dz)

    case default

       call Message('Error: must select ubar_type = linear|stc|exp|modal|&
            &read with nz>1 - yours is:'//trim(ubar_type))
       Parameters_ok=.false.

    end select
    
    ubar = ubar*uscale
    
    call Write_field(ubar,ubar_file)

  end function Init_ubar

  !*********************************************************************

  function Normu(ubarin,dz) result(ubarout)

    !************************************************************************
    ! Normalize initial mean velocities such that (1) there is no barotropic
    ! component, (2) ubar = ubar/Int(abs(ubar)dz).
    !************************************************************************
    
    real,dimension(:),intent(in)   :: ubarin, dz
    real,dimension(size(ubarin,1)) :: ubarout

    ubarout = ubarin - sum(ubarin*dz)
    if (maxval(abs(ubarout))/=0) ubarout = ubarout/maxval(abs(ubarout))
    
  end function Normu

  !*********************************************************************

end module qg_strat_and_shear
