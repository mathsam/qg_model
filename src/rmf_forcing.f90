module rmf_forcing

  !************************************************************************
  ! Tools for random markovian forcing.
  !
  ! Routines: init_rmf_forcing, markovian, get_gen_rmf_rate
  !
  ! Dependencies: io_tools, par_tools, op_rules, qg_params, qg_arrays,
  !               qg_filter_tools, numerics_lib
  !
  !************************************************************************

  implicit none
  private
  save

  complex, dimension(:,:), allocatable   :: force_o
  
  public :: init_rmf_forcing, markovian, get_gen_rmf_rate, force_o

contains

  !************************************************************************

  subroutine Init_rmf_forcing

    !************************************************************************
    ! If random Markovian forcing is to be used, check input params and,
    ! if its a restart run, read in forcing from previous timestep
    !************************************************************************
 
    use io_tools,  only: Message, Read_field
    use qg_params, only: kmax, kx_start, kx_end, nkx, nky,       &
                         parameters_ok,cr,ci,datadir,io_root,    &
                         forc_coef, forc_corr, kf_min, kf_max,   &
                         z_force, norm_forcing, force_o_file,    &
                         nc_restartfile
    use par_tools, only: par_scatter, par_sync
    use nc_io_tools, only: read_nc

    complex,dimension(:,:),allocatable      :: force_global
    logical                                 :: file_exists

    allocate(force_o(kx_start:kx_end,0:kmax))

    call Message('Random Markovian forcing on')
    if (forc_coef==cr) then
       call Message('Error: forc_coef must be set for RM forcing')
       parameters_ok=.false.
    else
       call Message('Forcing with coefficient forc_coef =',r_tag=forc_coef)
    endif
    if (forc_corr==cr) then
       call Message('Info: forc_corr not set - setting to .5')
       forc_corr = .5
    endif
    if ((forc_corr<0).or.(forc_corr>1)) then
       call Message('Error: require 0<=forc_corr<=1 - &
            &yours is:',r_tag=forc_corr)
       parameters_ok=.false.
    endif
    if (kf_min==cr)  then
       call Message('Error: kf_min must be set for RM forcing')
       parameters_ok=.false.
    else
       call Message('Minimum forcing wavenumber kf_min =',tag=int(kf_min))
    endif
    if (kf_max==cr)  then
       call Message('Error: kf_max must be set for RM forcing')
       parameters_ok=.false.
    else
       call Message('Maximum forcing wavenumber kf_max =',tag=int(kf_max))
    endif
    if (norm_forcing) then
       call Message('Info: generation rate set to forc_coef')
    endif
    if (z_force==ci) then
       call Message('Forcing level not set -- setting to 1.')
       z_force = 1
    endif

    force_o = 0.
    inquire(file="./INPUT/"//trim(nc_restartfile), exist=file_exists)
    if (file_exists) then
       call Message('Found force_o file -- reading in')
       allocate(force_global(-kmax-1:kmax,0:kmax))
       force_global = 0.
!       call Read_field(force_global(-kmax:kmax,:),force_o_file)
       call read_nc("./INPUT/"//trim(nc_restartfile), "forcing", &
                    force_global(-kmax:kmax,:))
       call Message("Read in forcing restart")
       call par_scatter(force_global,force_o,io_root)
       deallocate(force_global)
    endif

    call par_sync

    call Message('Forcing initialized')

  end subroutine Init_rmf_forcing

  !************************************************************************

  function Markovian(kf_min,kf_max,amp,lambda,norm_forcing,norm_diss, &
                     filter,filter_type,psi,q,dz,dt) &
       result(forc)
    
    !**************************************************************
    ! Random Markovian forcing function.  If norm_forcing = T, function
    ! will normalize the forcing such that the total generation = amp
    !**************************************************************
    
    use op_rules,     only: operator(+), operator(-), operator(*)
    use qg_arrays,    only: ksqd_, kx_, ky_
    use qg_params,    only: nkx,nky,kmax,nz,idum,i,pi,rmf_norm_min
    use qg_filter_tools, only: get_filter_rate
    use numerics_lib, only: Ran
    use io_tools,     only: Message

    complex,dimension(nkx,nky)                :: forc
    real,intent(in)                           :: kf_min,kf_max,amp,lambda
    logical,intent(in),optional               :: norm_forcing,norm_diss
    complex,dimension(:,:,:),intent(in)       :: psi, q
    real,dimension(:,:),intent(in)            :: filter
    character(*), intent(in)                  :: filter_type
    real,dimension(:),intent(in)              :: dz
    real,intent(in)                           :: dt
    ! Local
    real                                      :: gamma=1.,gr=0.,fr=0.

    where((ksqd_ > kf_min**2).and.(ksqd_ <= kf_max**2))
       force_o = lambda*force_o &
             + amp*sqrt(1-lambda**2)*cexp((i*2*pi)*Ran(idum,nkx,nky))
    endwhere
    forc = force_o
    if (norm_forcing) then

       gr = get_gen_rmf_rate(psi)
       if (gr<0.) then
          gr = -gr
          forc = -forc
       endif
       if (trim(filter_type)/='none'.and.norm_diss) then
          fr = get_filter_rate(psi,q,filter,dz,dt)
       else
          fr = 0.
       endif
       gamma = gr/(amp-fr)
       if (gamma>rmf_norm_min) then
           forc = forc/gamma
       else
          call Message('RandMarkForc: norm factor too small, =',r_tag=gamma)
       endif
    endif
 
    where (kx_<=0.and.ky_==0) forc = 0.
    force_o = forc
    
  end function Markovian

  !*************************************************************************

  real function get_gen_rmf_rate(psi)

    ! Get energy generation rate due to random Markovian forcing
    
    use op_rules,  only: operator(*)
    use par_tools, only: par_sum

    complex,dimension(:,:,:),intent(in) :: psi

    get_gen_rmf_rate = - 2*sum(real(conjg(psi)*force_o))
    call par_sum(get_gen_rmf_rate)

  end function get_gen_rmf_rate
 
  !*************************************************************************


end module rmf_forcing
