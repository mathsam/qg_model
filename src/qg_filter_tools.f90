module qg_filter_tools

  !************************************************************************
  ! Routines to set up filter and find its dissipation
  !
  ! Routine:  init_filter, get_filter_rate
  !
  ! Dependencies: io_tools, qg_params, qg_arrays, par_tools
  !
  !************************************************************************

  implicit none
  private
  save

  public :: init_filter, get_filter_rate

contains

  !************************************************************************

  function init_filter(filter_type,filter_exp,k_cut,dealiasing,ft,Nexp) &
       result(filter)
    
    !************************************************************************
    ! Set dealiasing mask for isotropic truncation (semicircle is just tangent
    ! to line in '4/3' rule) and combine it with small scale spatial filter.
    ! Parameter 'filtdec' below is set so that filter decays to 
    ! (1+4*pi/nx)**(-1) at max K allowed by de-aliasing form
    ! Initialize small scale filter/de-aliasing mask.  Legal filter types are:
    !
    !   hyperviscous : equivalent to RHS dissipation 
    !                  nu*del^(2*filter_exp)*field, with nu set optimally.
    !                  Requires 'filter_exp'
    !
    !   exp_cutoff   : exponential cutoff filter (see code below)
    !                  Requires 'filter_exp' and 'k_cut'
    !   
    !   none         : none
    !
    ! Character switch 'dealiasing' sets the de-aliasing type:
    !
    !   orszag       : clips corners of spectral fields at 
    !                  |k|+|l| = (4/3)*(kmax+1) as per Orszag '72
    !
    !   isotropic    : clips spectral fields at all points ouside circle
    !                  which inscribes limits of orszag filter (circle
    !                  tangent to |k|+|l| = (4/3)*(kmax+1) ==>
    !                  Kda = sqrt(8/9)*(kmax+1) ).  In this case, the 
    !                  *effective* Kmax = Kda.
    !
    !   none         : none
    !
    ! ft is tuning factor for max filter value, and Nexp sets
    ! exponent on resolution in (1+4*pi/nx**Nexp) (since SUQG requires 
    ! Nexp = 2).
    !
    !************************************************************************
 
    use io_tools,  only: Message
    use qg_params, only: nkx,kx_start,kx_end,ngrid,kmax,nky,pi,parameters_ok,cr,ci
    use qg_arrays, only: kx_, ky_, ksqd_
    use par_tools, only: par_sync

    character(*),intent(in)                 :: filter_type,dealiasing
    real,intent(in)                         :: filter_exp,k_cut,ft,Nexp
    real,dimension(kx_start:kx_end,0:kmax)  :: filter       ! Local
    real                                    :: filtdec
    integer                                 :: kmax_da

    filter = 1.  
    where (kx_<=0.and.ky_==0) filter = 0.            ! 0 since given by conj sym
    if (kx_start==-kmax-1) filter(kx_start,:) = 0.   ! Dummy column

    select case (trim(dealiasing))  
    case ('orszag')

       call Message('Using Orszag (non-isotropic) de-aliasing')
       kmax_da = sqrt(8./9.)*(kmax+1)
       where ( abs(kx_)+abs(ky_) >= (4./3.)*(kmax+1) ) filter = 0.

    case ('isotropic')
       
       call Message('Using isotropic de-aliasing')
       kmax_da = sqrt(8./9.)*(kmax+1)
       where ( ksqd_ >= (8./9.)*(kmax+1)**2 ) filter = 0.

    case ('none')

       call Message('Spectral convolutions will not be de-aliased')
       kmax_da = kmax

    case default

       call Message('Error: dealiasing must be one of orszag|isotropic|none; &
                     &yours is'//trim(dealiasing))
       parameters_ok = .false.

    end select

    ! Now set enstrophy filtering **only in region where filter is not zeroed
    ! from de-aliasing above**

    select case (trim(filter_type))  
    case ('hyperviscous') 

       call Message('Using hyperviscous filter')
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('...with (del**2)**',r_tag=filter_exp)
       endif

       where (filter>0.) 
          filter = 1/(1+ft*(4*pi/ngrid**Nexp)*(ksqd_/kmax_da**2)**filter_exp)
       endwhere

    case ('exp_cutoff') 

       call Message('Using exponential cutoff filter')
       if (k_cut==cr) then
          call Message('Error: cutoff scale k_cut not set')
          parameters_ok=.false.
       else
          call Message('Cutoff scale k_cut =',r_tag=k_cut)
       endif
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('Filter exponent filter_exp =',r_tag=filter_exp)
       endif

       filtdec = -log(1+ft*4*pi/ngrid**Nexp)/(kmax_da-k_cut)**filter_exp
       where ((ksqd_>k_cut**2).and.(filter>0.)) 
          filter = exp(filtdec*(sqrt(ksqd_)-k_cut)**filter_exp)
       end where

    case ('none')

       call Message('No spatial filtering')

    case default ! or filter_type = 'none' -- make that the default in decln.

       call Message('Error: Must select filter_type.  Legal choices are &
                     &filter_type = hyperviscous|exp_cutoff|none')
       parameters_ok=.false.

    end select

    call par_sync

  end function init_filter

  !*************************************************************************

  real function get_filter_rate(psi,q,filter,dz,dt)

    ! Get filter dissipation rate
    
    use op_rules,  only: operator(*)
    use qg_params, only: nky, nkx
    use par_tools, only: par_sum

    complex,dimension(:,:,:),intent(in) :: psi,q
    real, dimension(:,:), intent(in)    :: filter
    real,dimension(:), intent(in)       :: dz
    real,intent(in)                     :: dt
    real,dimension(:,:),allocatable     :: temp

    allocate(temp(nkx,nky)); temp = 0.
    where (filter/=0) temp = (filter**(-1)-1)/(2*dt)
    get_filter_rate = 2*sum(real( temp*(dz*(conjg(psi)*q)) ))
    call par_sum(get_filter_rate)
    deallocate(temp)

  end function get_filter_rate
 


end module qg_filter_tools
