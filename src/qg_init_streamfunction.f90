module qg_init_streamfunction       !-*-f90-*- 

  !************************************************************************
  ! Streamfunction initialization routin.
  !
  ! Routines:  init_streamfunction
  !
  ! Dependencies: op_rules, io_tools, par_tools, qg_params, qg_arrays, 
  !               qg_strat_and_shear, strat_tools, numerics_lib,
  !               qg_diagnostics, transform_tools
  !************************************************************************

  implicit none
  private
  save

  public :: init_streamfunction

contains

  !************************************************************************

  subroutine init_streamfunction

    !************************************************************************
    ! Read in or create initial streamfunction field in manner specified
    !   by variable 'psi_init_type', which can have the values:
    !
    ! spectral_m :  Gaussian spread of initial energy about isotropic horiz.
    !               wavenumber 'k_o' and with width 'delk', and all energy
    !               in vertical mode 'm_o'
    ! all_scales:   All scales given energy
    ! spectral_z :  Same as above in horizontal plane, but with all initial
    !               energy in level 'z_o'
    ! elliptical_vortex :  Elliptical gaussian bump in initial grid vorticity
    !               field, aspect ratio 'aspect_vort' and width 'del_vort',
    !               and contained in mode 'm_o'
    ! read :        Read in from 'psi_init_file' at frame 'start_frame'
    ! 
    ! All of the values in quotes can be set in input namelist
    !************************************************************************

    use op_rules,        only: operator(+), operator(-), operator(*)
    use io_tools,        only: Message, Read_field
    use qg_params,       only: kmax, nky, ngrid, nz, nkx, y_start, y_end, ny,  &
                               z_o, k_o, m_o, e_o, delk, aspect_vort, del_vort,&
                               i, pi, idum, F,                                 &
                               psi_init_file,                                  &
                               frame, start_frame, rewindfrm,                  &
                               parameters_ok, cr, ci, io_root,                 &
                               initialize_energy, psi_init_type, restarting,   &
                               nc_restartfile
    use qg_arrays,       only: ksqd_, kx_, ky_, psi
    use qg_strat_and_shear, only: kz, vmode
    use qg_diagnostics,  only: energy
    use numerics_lib,    only: Ran
    use transform_tools, only: Grid2spec
    use strat_tools,     only: Layer2mode, Mode2layer
    use par_tools,       only: par_scatter, par_sync, processor_id
    use nc_io_tools,     only: read_nc

    real,dimension(:,:),allocatable            :: espec,mu
    real,dimension(:,:,:),allocatable          :: zetag
    complex,dimension(:,:,:),allocatable       :: zeta, psim
    complex,dimension(:,:,:),allocatable       :: psi_global
    real                                       :: e, radius2, total_energy
    integer                                    :: ix, iy, y, midx, midy, iz


    if (trim(psi_init_type)=='spectral') psi_init_type = 'spectral_z'

    restart: if (restarting) then

       allocate(psi_global(1:nz,-kmax-1:kmax,0:kmax)); psi_global=0.
       if (.not.rewindfrm) then
          call read_nc("./INPUT/"//trim(nc_restartfile),"psi", psi_global(:,-kmax:kmax,:))
       else
          call read_nc("./INPUT/"//trim(nc_restartfile),"psi", psi_global(:,-kmax:kmax,:))
       endif
       call par_scatter(psi_global,psi,io_root)
       deallocate(psi_global)

    else 

       psi = 0.
       select case (trim(psi_init_type))

       case('all_scales')

          call Message('Initial energy will be spectrally global')

          do iz=1,nz
             psi(iz,:,:) = 1./sqrt(ksqd_+4.*F)*cexp(i*2*pi*Ran(idum,nkx,nky))
          enddo

       case('spectral_m')

          call Message('Initial streamfunction will be spectrally local')
          if (k_o==cr) then
             call Message('Error: k_o must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (k_o<=0) then
             call Message('Error: require k_o > 0 - yours is:', r_tag=k_o)
             Parameters_ok=.false.
          else
             call Message('Initial energy centroid at isotropic wavenumber &
                  &k_o =', r_tag=k_o)
          endif
          if (delk==cr) then
             call Message('Error: delk must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (delk<=0) then
             call Message('Error: need delk>0 - yours is:',r_tag=delk)
             Parameters_ok=.false.
          else
             call Message('Initial energy peak wavenumber width delk =',&
                  r_tag=delk)
          endif
          if (nz>1) then
             if (m_o==ci) then
                call Message('Error: need m_o for modal streamfunction')
                Parameters_ok=.false.
             elseif ((m_o<0).or.(m_o>nz-1)) then
                call Message('Error: require 0<=m_o<= nz-1 - yours is:',&
                              tag=m_o)
                Parameters_ok=.false.
             else
                call Message('Streamfunction will be initd in mode m_o =',&
                     tag=m_o)
             endif
          else
             m_o = 0
          endif

          allocate(espec(nkx,nky),mu(nkx,nky)); espec = 1.; mu = 0.
          allocate(psim(nz,nkx,nky)); psim = 0.
          if (delk/=0) espec = exp(-(sqrt(ksqd_)-k_o)**2/delk**2) 
          mu = sqrt(ksqd_+kz(m_o+1)**2)        ! Total geostrophic wavenumber
          psim(m_o+1,:,:)= sqrt(espec)/mu*cexp(i*2*pi*Ran(idum,nkx,nky))
          psi(1:nz,:,:) = Mode2layer(psim,vmode)
          deallocate(espec,mu,psim)
          
       case('spectral_z')
          
          call Message('Initial streamfunction will be spectral by layer')
          if (k_o==cr) then
             call Message('Error: k_o must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (k_o<=0) then
             call Message('Error: require k_o > 0 - yours is:', r_tag=k_o)
             Parameters_ok=.false.
          else
             call Message('Initial energy centroid at isotropic wavenumber &
                  &k_o =', r_tag=k_o)
          endif
          if (delk==cr) then
             call Message('Error: delk must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (delk<0) then
             call Message('Error: require delk>=0 - yours is:',r_tag=delk)
             Parameters_ok=.false.
          else
             call Message('Initial energy peak wavenumber width delk =',&
                  r_tag=delk)
          endif
         if (nz>1) then
             if (z_o==ci) then
                call Message('Error: z_o must be set for streamfunction')
                Parameters_ok=.false.
             elseif ((z_o<=0).or.(z_o>nz)) then
                call Message('Error: need 0<z_o<=nz - yours is:',tag=z_o)
                Parameters_ok=.false.
             else
                call Message('Streamfunction will be initd in layer z_o =',&
                           tag=z_o)
             endif
          else
             z_o = 1
          endif
       
          allocate(espec(nkx,nky)); espec = 1.
          if (delk/=0) espec = exp(-(sqrt(ksqd_)-k_o)**2/delk) 
          psi(z_o,:,:) = sqrt(espec)/sqrt(ksqd_)*cexp(i*2*pi*Ran(idum,nkx,nky))
          deallocate(espec)

       case('elliptical_vortex')

          call Message('Initial vorticity field will be an elliptical vortex')
          if (del_vort==cr) then
             call Message('Error: del_vort must be set to make streamfunction')
             Parameters_ok=.false.
          else
             call Message('Initial vortex width del_vort =',r_tag=del_vort)
          endif
          if (aspect_vort==cr) then
             call Message('Error: aspect_vort must be set to make streamfuncn')
             Parameters_ok=.false.
          else
             call Message('Initial vortex aspect ratio aspect_vort =', &
                  r_tag=aspect_vort)
          endif
          if (nz>1) then
             if (m_o==ci) then
                call Message('Error: m_o must be set for modal streamfunction')
                Parameters_ok=.false.
             elseif ((m_o<0).or.(m_o>nz-1)) then
                call Message('Error: require 0<= m_o<=nz-1: yours is:',tag=m_o)
                Parameters_ok=.false.
             else
                call Message('Streamfunction will be initd in mode m_o =',&
                     tag=m_o+1)
             endif
          else
             m_o = 0
          endif

          allocate(psim(nz,nkx,nky),zeta(1,nkx,nky))
          psim=0.; zeta=0.
          allocate(zetag(1,ny,ngrid)); zetag = 0.
          midx = kmax+1; midy = midx
             
          do y = y_start, y_end
             iy = y - y_start + 1  ! goes from 1 to ny on local proc
             do  ix = 1, ngrid
                radius2 = ((ix-midx)**2+aspect_vort*(y-midy)**2)
                radius2 = radius2*((2*pi)**2/ngrid**2)
                zetag(1,iy,ix) = exp(-radius2/del_vort**2)
             enddo
          enddo
          
          call Grid2spec(zetag(1,:,:),zeta(1,:,:))
          psi=0.
          psim(m_o+1,:,:) = -(1./ksqd_)*zeta(1,:,:)
 !         psi(1:nz,:,:) = Mode2layer(psim,vmode)
          deallocate(zetag,zeta,psim)
          
       case('read')
          
          if (trim(psi_init_file)=='') then
             call Message('Error: no input file for psi given')
             Parameters_ok=.false.
          endif
          if (start_frame==ci) then
             call Message('Warning: start_frame not initialized-setting to 1')
             start_frame = 1
          elseif (start_frame<=0) then
             call Message('Error: require start_frame>=0')
             Parameters_ok=.false.
          endif
          call Message('Initial streamfunction will be read from: '&
               &//trim(psi_init_file)//', frame:', tag=start_frame)

          allocate(psi_global(1:nz,-kmax-1:kmax,0:kmax)); psi_global=0.
          call Read_field(psi_global(:,-kmax:kmax,:),psi_init_file,frame=start_frame, &
               exclude_dd=1,zfirst=1)  ! only reads on proc 0 
          call par_scatter(psi_global,psi,io_root)  ! scatters from proc 0
          deallocate(psi_global)

       case default
          
          call Message('Error: must select psi_init_type = &
               &read|spectral_m|spectral_z|spectral|elliptical_vortex &
               &|all_scales - yours is:'//trim(psi_init_type))
          Parameters_ok=.false.
          
       end select

       if (initialize_energy) then
          if (e_o==cr) then
             call Message('Info: e_o not set - setting to 1.')
          elseif (e_o<=0.) then
             call Message('Error: must have e_o > 0')
             parameters_ok = .false.
          else
             call Message('Initial energy will be set to e_o =', r_tag=e_o)
          endif

          total_energy = energy(psi)
          
          if (total_energy>0) then
             psi = sqrt(e_o/total_energy)*psi
          elseif (total_energy<=epsilon(e_o)) then
             call Message('Error: no amplitude in initial psi field:',&
                  r_tag=total_energy)
             parameters_ok = .false.
          endif
       endif
    endif restart 
    
    do iz=1,nz
       where (kx_<=0.and.ky_==0) psi(iz,:,:) = 0.
    enddo

    call par_sync

    call Message('Streamfunction initialized')

  end subroutine  init_streamfunction

  !*********************************************************************

end module qg_init_streamfunction
