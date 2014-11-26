module qg_params                   !-*-f90-*-

  !************************************************************************
  ! Contains all global parameters for program, and routines for
  ! processing their I/O
  !
  ! Routines: Initialize_parameters, Write_parameters, Check_parameters
  !
  ! Dependencies: IO_tools, par_tools
  !************************************************************************

  implicit none
  public
  save

  integer,parameter       :: ci=-9999      ! For checking init status of variables
  real,parameter          :: cr=-9999.

  ! Parameters in namelist input -- initialize some parameterss to
  ! negative values (cr and ci) in order to check they are intentionally
  ! initialized in Parameters_ok (below)

  ! Spatial parameters

  integer                 :: nz          = ci          ! Vertical resolution
  integer                 :: kmax        = ci          ! Spectral resolution

  ! Temporal parameters

  real                    :: dt              = cr      ! Time step
  logical                 :: adapt_dt        = .true.  ! Adaptive timestep
  real                    :: time            = 0
  logical                 :: linear          = .false. ! Omit non-linear terms

  ! Counters, flags, etc

  integer                 :: write_step      = ci      ! Frame snapshot step
  integer                 :: diag1_step      = ci      ! Diagnostics 1 step
  integer                 :: diag2_step      = ci      ! Diagnostics 2 step
  integer                 :: restart_step    = ci      ! Restart save frequency
  integer                 :: total_counts    = ci      ! Total timesteps to do
  integer                 :: start_frame     = ci      ! Frame to start from
  integer                 :: start_frame_t   = ci      ! Frame to start trc frm
  integer                 :: cntr            = 1       ! Timestep counter value
  integer                 :: frame           = 0       ! Field snapshot frame
  integer                 :: restart_frame   = 0       ! Restart frame
  integer                 :: d1frame         = 0       ! Diagnostics 1 frame
  integer                 :: d2frame         = 0       ! Diagnostics 2 frame
!  integer                 :: d3frame         = 0       ! Diagnostics 3 frame
  character(8)            :: end_date        = ''      ! End date of sim
  character(10)           :: end_time        = ''      ! End time of sim
  logical                 :: restarting  = .false.     ! Is this a restart?
  logical                 :: do_energetics   = .false.

  ! Scale parameters

  real                    :: beta        = cr          ! beta_0*L^2/[(2pi)^2 U]
  real                    :: F           = cr          ! f^2*L^2/[(2pi)^2g'H_0]
  real                    :: Fe          = 0.          ! F*drho(1)/drho_ext

  ! Stratification

  character(20)           :: surface_bc     = ''       ! Surface BC type
  character(20)           :: strat_type     = ''       ! Stratification type
  real                    :: deltc       = cr          ! ~thermocline thickness
  character(85)           :: dz_in_file  = ''          ! Layer thicknesses
  character(85)           :: rho_in_file = ''          ! Density profile
  logical                 :: read_tripint = .false.    ! Override tripint calc
  character(85)           :: tripint_in_file= ''       ! tripint file profile

  ! Mean velocity

  character(20)           :: ubar_type   = ''          ! Mean u type
  character(20)           :: vbar_type   = ''          ! Mean u type
  integer                 :: umode       = ci          ! Mode of Ubar ('mode')
  real                    :: delu        = cr          ! ~Ubar surf intens'n
  real                    :: uscale      = 0.          ! Scale factor for Ubar
  real                    :: vscale      = 0.          ! Scale factor for Vbar
  character(85)           :: ubar_in_file= ''          ! Ubar profile
  character(85)           :: vbar_in_file= ''          ! Ubar profile
  logical                 :: normalize_ubar=.true.     ! Remove BT & norm by max

  ! Time-varying mean velocity

  real                    :: inv_restore_time= 0      ! Inverse mean restore time

  ! Initial streamfunction

  character(20)           :: psi_init_type  = ''       ! Init streamfn type
  logical                 :: initialize_energy= .true. ! Set init energy to e_o
  real                    :: e_o         = cr          ! Initial energy
  real                    :: k_o         = cr          ! Initial peak k
  real                    :: delk        = cr          ! Initial spread in k
  real                    :: aspect_vort = cr          ! Initl vort aspct ratio
  real                    :: del_vort    = cr          ! Initial vortex width
  integer                 :: m_o         = ci          ! Initial modal peak
  integer                 :: z_o         = ci          ! Initial level
  character(85)           :: psi_init_file = ''        ! Psi input field

  ! Dissipation

  character(20)           :: filter_type = ''          ! Spatial filter 
  real                    :: filter_exp  = cr          ! Filter exponent
  real                    :: filt_tune   = 1.           ! Tuning for ens filter
  real                    :: k_cut       = cr          ! Exp cutoff scale
  logical                 :: limit_bot_drag = .false.  ! Apply drg to k<kb_min
  real                    :: kb_max      = cr          ! max k for bot drag  
  real                    :: bot_drag    = 0.          ! Bottom drag
  real                    :: therm_drag  = 0.          ! Thermal drag
  real                    :: top_drag    = 0.          ! Top drag
  real                    :: quad_drag   = 0.          ! Quadratic bottom drag
  real                    :: qd_angle    = 0.          ! Quad drag turn angle

  ! Topography

  logical                 :: use_topo    = .false.     ! Include topography
  character(20)           :: topo_type   = ''          ! Topography type
  real                    :: toposcale   = cr          ! Scale factor for topo
  real                    :: del_topo    = cr          ! Bump width in k or x
  real                    :: k_o_topo    = cr          ! Peak posn for kspc hb
  character(85)           :: hb_in_file  = ''          ! Bottom topo (spec)

  ! Random Forcing

  logical                 :: use_forcing = .false.     ! Use rndm markov frcing
  logical                 :: norm_forcing= .false.     ! Norm gen rate from RMF
  logical                 :: norm_diss   = .false.     ! Norm diss rate in RMF
  real                    :: forc_coef   = cr          ! BT forcing coefficient
  real                    :: forc_corr   = cr          ! BT forcing correlation
  real                    :: kf_min      = cr          ! min k^2 for BT frc
  real                    :: kf_max      = cr          ! max k^2 for BT frc
  integer                 :: z_force     = ci          ! level to force at

  ! Tracer

  logical                 :: use_tracer_x  = .false.   ! Calc tracers
  logical                 :: use_tracer_y  = .false.   ! Calc tracers
  logical                 :: use_mean_grad_t=.false.   ! Use mean trcr gradient
  character(20)           :: tracer_init_type = ''     ! Initial tracer type
  integer                 :: nzt           = ci        ! Vertical tracer res
  integer                 :: maxmode       = ci        ! highest mode for stirring
  character(20)           :: filter_type_t    = ''     ! Spatial filter tracer 
  real                    :: k_cut_t       = cr        ! Exp cutoff wavenumber
  real                    :: filter_exp_t  = cr        ! Filter_t exponent
  real                    :: filt_tune_t   = 1.        ! Tuning for tvar filt
  real                    :: kappa_v       = 0.        ! Vertical diffusivity
  real                    :: kappa_h       = 0.        ! Laplacian diffusivity
  character(85)           :: tracer_y_init_file = ''   ! Initial tracer field
  character(85)           :: tracer_x_init_file = ''   ! Initial tracer field
  character(85)           :: ubart_in_file  = ''       ! For hi z res tracer 
  character(85)           :: vbart_in_file  = ''       ! For hi z res tracer 
  character(85)           :: dzt_in_file    = ''       ! "
  character(85)           :: vmodet_in_file = ''       ! "

  ! Spectral and other diagnostics

  logical                 :: do_spectra       = .true. ! Calc/save spectra
  logical                 :: do_aniso_spectra = .false.! Calc/save spectra
  logical                 :: do_xfer_spectra  = .false.! Calc transfer specs
  logical                 :: do_genm_spectra  = .false.! Calc modal generation
  logical                 :: do_x_avgs        = .false.! Calc zonal averages
  logical                 :: do_corrs         = .false.! Calc shear/strain corrs


  ! DIP Switches and tuning factors - factory settings.
  ! All of these are included in namelist input as well, but they
  ! are not checked for initialization in Parameters_ok so that you
  ! dont have to include them in the input file (but can if you need to).

  integer                 :: recunit     = 8            ! For direct access IO
  integer                 :: idum        = -7           ! Random generator seed

  ! Numerical stability tuning parameters

  real                    :: robert      = 0.01         ! Robert filter value
  real                    :: dt_tune     = 0.3          ! Tuning for adaptv dt
  integer                 :: dt_step     = 10           ! dt re-calc interval
  real                    :: rmf_norm_min= 1.e-5        ! RMF genn min 4 normn
                                                        !   for get_tripint
  real                    :: drt         = -9000.       ! 1st derivs at bnds for
  real                    :: drb         = -9000.       ! spline in get_tripint
  real                    :: rho_slope   = .00005       ! Lin slope * exp prof
  integer                 :: hf          = 10           ! High res interp factr

  character(20)           :: dealiasing  = 'orszag'     ! Set de-aliasing form
  character(20)           :: dealiasing_t= 'orszag'     ! Set de-aliasing form

  ! **************End of namelist input parameters*************
  !
  ! Parameters for global internal use - NOT included in any namelist input
  !

  character(32)           :: version                = '4.0/MPI'

  ! Output file and directory names - diagnostics outputs are defined in 
  ! respective diagnostics module.

  character(80)           :: datadir                = '.'
  character(20)           :: inputfile              = 'input.nml'
  character(20)           :: restartfile            = 'restart.nml'
  character(20)           :: nc_restartfile         = 'restart.nc' !NetCDF restart
  character(20)           :: psi_file               = 'psi'         ! Psi snapshots 
  character(20)           :: psi_restart_file       = 'psi_restart' ! Psi snapshots 
  character(20)           :: force_o_file           = 'force_o'     ! Forcing 
  character(20)           :: tracer_y_file          = 'tracer_y'    ! Tracer snapshots 
  character(20)           :: tracer_x_file          = 'tracer_x'    ! Tracer snapshots 
  character(20)           :: tracer_y_restart_file  = 'tracer_y_restart' ! Tracer snapshots 
  character(20)           :: tracer_x_restart_file  = 'tracer_x_restart' ! Tracer snapshots 
  character(20)           :: psiq_file              = 'psiq'
  character(20)           :: vmode_file             = 'vmode'
  character(20)           :: kz_file                = 'kz'
  character(20)           :: tripint_file           = 'tripint'
  character(20)           :: write_time_file        = 'write_time'
  character(20)           :: restart_time_file      = 'restart_time'
  character(20)           :: ubar_file              = 'ubar'        ! Ubar profile
  character(20)           :: vbar_file              = 'vbar'        ! Vbar profile
  character(20)           :: dz_file                = 'dz'          ! Layer thicknesses
  character(20)           :: rho_file               = 'rho'         ! Density profile
  character(20)           :: hb_file                = 'hb'          ! Bottom topo (spec)
  character(20)           :: stop_file              = 'stop'        ! Exist of file=>stop

  ! Resolution and tuning parameters (set as functions of kmax ) 

  integer                 :: numk, nky, ngrid, nmask

  ! Internal flags and counters

  logical                 :: surf_buoy         = .false.
  logical                 :: time_varying_mean = .false. 
  logical                 :: parameters_ok     = .true.
  logical                 :: rewindfrm         = .false.
  logical                 :: start
  integer                 :: cnt = 1, call_q=0, call_b=0
  integer                 :: call_ty=0, call_tx=0, call_U=0, call_V=0

  ! Parameters that are set by FFTW / par_tools if MPI is used

  integer                 :: kx_start, kx_end, nkx, y_start, y_end, ny
  integer                 :: io_root = 0

  ! Cabalistic numbers 

  real,parameter          :: pi          = 3.14159265358979
  complex,parameter       :: i           = (0.,1.)


  ! Namelist declarations

  ! Spatial parameters
  namelist/run_params/nz, kmax

  ! Temporal parameters
  namelist/run_params/dt,adapt_dt,time,linear

  ! Counters, flags, etc
  namelist/run_params/write_step,diag1_step,diag2_step,restart_step
  namelist/run_params/total_counts,start_frame,start_frame_t
  namelist/run_params/cntr,frame,restart_frame,d1frame,d2frame !,d3frame
  namelist/run_params/end_date,end_time,restarting

  ! Scale parameters
  namelist/run_params/F,Fe,beta

  ! Stratification
  namelist/run_params/strat_type,surface_bc
  namelist/run_params/deltc,read_tripint
  namelist/run_params/dz_in_file,rho_in_file,tripint_in_file

  ! Mean velocity
  namelist/run_params/ubar_type,vbar_type
  namelist/run_params/delu,umode,uscale,vscale
  namelist/run_params/ubar_in_file,vbar_in_file
  namelist/run_params/normalize_ubar

  ! Time-varying mean velocity
  namelist/run_params/inv_restore_time

  ! Initial streamfunction
  namelist/run_params/initialize_energy,psi_init_type
  namelist/run_params/e_o,k_o,delk,aspect_vort,del_vort,m_o,z_o
  namelist/run_params/psi_init_file

  ! Dissipation
  namelist/run_params/filter_type
  namelist/run_params/filter_exp,filt_tune,k_cut
  namelist/run_params/limit_bot_drag,kb_max,bot_drag,therm_drag,top_drag
  namelist/run_params/quad_drag,qd_angle

  ! Topography
  namelist/run_params/use_topo,topo_type
  namelist/run_params/toposcale,del_topo,k_o_topo
  namelist/run_params/hb_in_file
  
  ! Random Forcing
  namelist/run_params/use_forcing,norm_forcing,norm_diss
  namelist/run_params/forc_coef,forc_corr,kf_min,kf_max,z_force

  ! Tracer
  namelist/run_params/use_tracer_x,use_tracer_y,tracer_init_type
  namelist/run_params/use_mean_grad_t,kappa_v,kappa_h
  namelist/run_params/nzt,maxmode
  namelist/run_params/filter_type_t,filter_exp_t,k_cut_t,filt_tune_t
  namelist/run_params/tracer_y_init_file,tracer_x_init_file
  namelist/run_params/ubart_in_file,vbart_in_file,dzt_in_file,vmodet_in_file

  ! Spectral diagnostics
  namelist/run_params/do_spectra,do_xfer_spectra,do_genm_spectra
  namelist/run_params/do_aniso_spectra,do_x_avgs,do_corrs

  ! Tuning factors, DIP switches
  namelist/run_params/recunit, idum
  namelist/run_params/robert,dt_tune,dt_step
  namelist/run_params/rmf_norm_min
  namelist/run_params/drt,drb,rho_slope,hf
  namelist/run_params/dealiasing,dealiasing_t

  ! switchs: if do energetics
  namelist/run_params/do_energetics

  ! End of Namelist definition

contains

  !************** Routines for parameter I/O*********************

  subroutine Init_parameters

    !************************************************************************
    ! Read in command line arguments 'datadir' and 'inputfile'
    ! then read namelist in 'inputfile' for input parameters,
    ! and pass some parameters to 'io_tools' which are needed for I/O.
    ! If no first arg supplied, then datadir = ./, and if no
    ! 2nd arg, then inputfile = input.nml (in datadir)
    !
    ! MPI must already be initialized before calling
    !
    !************************************************************************

    use io_tools,  only: Message, Pass_params
    use par_tools, only: processor_id, par_bcast
    use nc_io_tools, only : pass_params_nc_io

    character(80)         :: progname='',datadirin='',inputfilein='',fnamein=''
    integer               :: fin=7, iock, nchars
    logical               :: restart_exists=.false.
    logical               :: restart_dir_exists = .false.

    if (processor_id == io_root) then
       ! Get executable name
       call getarg(0,progname)

       ! Get data directory
       call getarg(1,datadirin)
       if (trim(datadirin)/='') datadir = datadirin
       nchars = len_trim(datadir)
       if (datadir(nchars:nchars)/='/') datadir=trim(datadir)//'/'
       
       ! Get input file.  If restart exists, use this, but override if
       ! a second command line arg is specified.
       inquire(file=trim(datadir)//trim(restartfile),exist=restart_exists)
       if (restart_exists) inputfile = restartfile
       call getarg(2,inputfilein)
       if (trim(inputfilein)/='') inputfile = inputfilein

       ! Check if ./INPUT and ./RESTART directories exists; otherwise stop
       inquire(DIRECTORY="./INPUT",   EXIST = restart_dir_exists)
       if (.not. restart_dir_exists) then
           call Message('INPUT dir does not exist. Exiting!')
           Stop "Stopping"
       endif

       inquire(DIRECTORY="./RESTART", EXIST = restart_dir_exists)
       if (.not. restart_dir_exists) then
           call Message('RESTART dir does not exist. Exiting!')
           Stop "Stopping"
       endif

    endif
    call par_bcast(progname,io_root)
    call par_bcast(datadir,io_root)
    call par_bcast(inputfile,io_root)

    ! Send info to io_tools 

    call Pass_params(datadir,recunit,processor_id,io_root) 

    call Message('')
    call Message('SQG model version: '//trim(version))
    call Message('')
    call Message('Executable:        '//trim(progname))
    call Message('Data directory:    '//trim(datadir))
    call Message('Input file:        '//trim(datadir)//inputfile)
    call Message('')

    ! Read in namelist on every pe

    call Message('Reading namelist')
    fnamein = trim(datadir)//trim(inputfile)
    open(unit=fin,file=fnamein,status='old',delim="apostrophe",iostat=iock)
    if (iock/=0) call Message('Open input namelist error; file:'//fnamein,&
         tag=iock,fatal='y')
    read (unit=fin,nml=run_params,iostat=iock)
    if (iock/=0) call Message('Read input namelist error; file, iock ='//fnamein,&
         tag=iock,fatal='y')
    close(fin)

    call Message('Checking parameters for consistency')
    call Message('')
    parameters_ok = Check_parameters()  ! See below - check and set counters

    ngrid = 2*(kmax+1)                  ! Physical resolution in x
    nky = kmax+1                        ! Number of meridional wavenumbers

    ! Make sure random num gen is set to start right
    idum = -abs(idum)

    call pass_params_nc_io(processor_id,io_root,kmax,nz)

  end subroutine Init_parameters

  !**************************************************************

  subroutine Write_parameters

    !************************************************************************
    ! Write parameters to restart namelist file and selected
    ! params to a special file for use by SQG Matlab routines.
    !************************************************************************

    use io_tools,  only: Message, Open_file
    use par_tools, only: processor_id

    character(85) :: restartname, exitcheckname
    integer       :: outf = 40,iock

    restartname = trim(datadir)//trim(restartfile)

    ! Set certain parameters the way they should be set for restarting,
    ! regardless of what goes on in the program.  Store temp vals as nec.

    restarting = .true.
    start_frame = frame
    start_frame_t = frame
    call date_and_time(DATE=end_date,TIME=end_time)

    if (processor_id==io_root) then

       ! Write full namelist to restartname (restart.nml)

       open(unit=outf,file=restartname,status='unknown',delim="apostrophe", iostat=iock)
       if (iock /= 0) call Message('write params nml open error; &
            &iock =',tag=iock,fatal='y')
       write (unit=outf,nml=run_params,iostat=iock)
       if (iock /= 0) call Message('write params nml write error; &
            &iock =',tag=iock,fatal='y')
       close(outf)

    endif
 
  end subroutine Write_parameters

  !**************************************************************

  logical function Stop_run()

    !**************************************************************
    ! Check for existence of file stopfile in datadir
    !**************************************************************
    use io_tools,  only: Message

    inquire(file=trim(datadir)//trim(stop_file),exist=stop_run)
    if (stop_run) call message('Stop file found')

  end function Stop_run

  !**************************************************************

  logical function Check_parameters()

    !**************************************************************
    ! This function will test that certain input namelist params are consistent
    ! and print some messages with basic info for run
    !**************************************************************

    use io_tools,  only: Message
    use par_tools, only: num_processors

    Check_parameters=.true.

    ! Resolution

    if (kmax==ci) then
       call Message('Error: kmax not set')
       Check_parameters=.false.
    elseif (mod(kmax+1,2)/=0) then
       call Message('Error: kmax must be odd - yours is: ',tag=kmax)
       Check_parameters=.false.
    else
       call Message('Horizontal spatial resolution =', tag=2*(kmax+1))
    endif

    ! Check that resolution requested is ok for MPI
    if( modulo(2*(kmax+1),num_processors) /= 0  ) then
       call Message('Error: num_processors must divide ngrid.  Yours is ', tag=2*(kmax+1))
       Check_parameters=.false.
    endif

    if (nz==ci) then
       call Message('Error: nz not set')
       Check_parameters=.false.
    elseif (nz<1) then
       call Message('Error: nz must be positive - yours is:',tag=nz)
       Check_parameters=.false.
    else
       call Message('Vertical resolution =', tag=nz)
    endif

    if (adapt_dt) then
       call Message('Using adaptive time-stepping')
    else
       call Message('Using fixed time-step')
       if (dt==cr) then
          call Message('Error: dt must be set explicitly in non-adaptive mode')
          Check_parameters=.false.
       endif
    endif

    ! Check basic scales
    
    if (F==cr.and.(nz==1)) then
       call Message('Info: F not set - setting to 0')
       F = 0.
    elseif ((F<=0.).and.(nz>1)) then
       call Message('Error: require F>0 with nz>1')
       Check_parameters=.false.
    else
       call Message('F parameter =',r_tag=F)
    endif

    if (Fe/=0) then
       call Message('Free surface selected:  Fe=',r_tag=Fe)
       if (nz==2) then
          call Message('For nz=2, need equal layer depths.  Note that tripint&
                        & not calculated correctly for this case yet')
       endif
    endif

    if (beta==cr) then
       call Message('Error: beta not set')
       Check_parameters=.false.
    else
       call Message('Beta-plane: beta =',r_tag=beta)
    endif

    call Message('Linear bottom drag = ',r_tag=bot_drag)
    if (bot_drag/=0) then
       if (limit_bot_drag) then
          if (kb_max==cr) then
             call Message('Error: Must set kb_max if limit_bot_drag is true')
             Check_parameters = .false.
          else
             call Message('Linear bottom drag limited to K < kb_max =',r_tag=kb_max)
          endif
       endif
    endif

    call Message('Thermal drag = ',r_tag=therm_drag)
    call Message('Quadratic bottom drag =', r_tag=quad_drag)

  end function Check_parameters

!*************************************************************************

end module qg_params
