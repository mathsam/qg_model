module qg_output

  !**********************************************************************
  ! Main output routines for restarts and snapshots.  Also initializes counters
  !
  ! Routines: init_counters, write_restarts, write_snapshots
  !
  ! Dependencies: io_tools, par_tools, qg_params, qg_arrays, rmf_forcing, qg_tracers 
  !
  !**********************************************************************

  implicit none
  private
  integer :: history_file_id, restart_file_id
  integer :: psi_var_id, tracerx_var_id, tracery_var_id, time_var_id
  integer :: restart_psi_var_id, restart_tracerx_var_id,   &
             restart_tracery_var_id, restart_force_var_id
  save

  public :: init_counters, write_restarts, write_snapshots

contains

  !************************************************************************

  subroutine Init_counters

    !************************************************************************
    ! Initialize diagnostic and snapshot frame counters
    !************************************************************************

    use io_tools,  only: Message
    use qg_params, only: ci, cr, parameters_ok, cntr, cnt, total_counts,   &
                         start_frame, frame, d1frame, d2frame, rewindfrm,  & 
                         write_step, diag1_step, diag2_step,               &
                         do_spectra, restarting, restart_step, psi_file,   &
                         use_tracer_x, use_tracer_y
    use nc_io_tools, only : create_file, enddef_file, register_variable,     &
                            create_axis_time, create_axis_kx, create_axis_ky,&
                            create_axis_z, create_axis_real_and_imag

    integer :: axis_kx_id, axis_ky_id, axis_z_id, axis_time_id, axis_compl_id

    restart: if (restarting) then

       call Message('This is a restart run')
       if (start_frame<0.or.start_frame>frame) then
          call Message('Error: require 0 <= start_frame <= frame')
          parameters_ok = .false.
       elseif (start_frame/=frame) then      ! Set counters for restart
          rewindfrm = .true.
          frame = start_frame 
          if (frame==0) then
             cntr = 1
             d1frame = 0
             d2frame = 0
          else
             cntr = frame*write_step
             d1frame = (cntr-1-mod(cntr-1,diag1_step))/diag1_step + 1 
             d2frame = (cntr-1-mod(cntr-1,diag2_step))/diag2_step + 1
          endif
          call Message('Will read streamfunction from:'// &
                trim(psi_file)//', frame:',tag=start_frame+1)
       endif

    elseif (.not.restarting) then

       frame = 0; d1frame = 0; d2frame = 0; cntr = 1

       call Message('Checking counters...')  ! Check internal counters
       
       if (total_counts==ci) then
          call Message('Error: total_counts not set')
          parameters_ok=.false.
       else
          call Message('Total number of timesteps model will run =', &
                        tag=total_counts)
       endif
       if (diag1_step==ci) then
          call Message('Info: timeseries write interval not set - &
               &setting diag1_step = 50')
          diag1_step = 100
       endif
       call Message('Timeseries will be written at timestep interval =', &
                     tag=diag1_step)
       if (diag2_step==ci.and.do_spectra) then
          call Message('Info: spectra write interval not set - setting diag2_step = 100')
          diag2_step = 100
       endif
       call Message('Spectra will be written at timestep interval =', &
                     tag=diag2_step)
       if (write_step==ci) then
          call Message('Error: Full field snapshot write interval not set - &
               &choose a value for write_step')
          parameters_ok=.false.
       else
          call Message('Snapshots will be written at timestep interval =', &
                     tag=write_step)
       endif
       if (restart_step==ci) then
          call Message('Restart_step not set - setting to diag1_step')
          restart_step=diag1_step
       else
          call Message('Restart files written at timestep interval =', &
                                tag=restart_step)
       endif

    endif restart

    cnt = 1 

    call Message('Counters initialized')

    history_file_id = create_file("history.nc")
    axis_time_id    = create_axis_time(history_file_id)
    axis_kx_id      = create_axis_kx(history_file_id)
    axis_ky_id      = create_axis_ky(history_file_id)
    axis_z_id       = create_axis_z(history_file_id)
    axis_compl_id   = create_axis_real_and_imag(history_file_id)

    psi_var_id  = register_variable(history_file_id, "psi",  &
        (/axis_z_id, axis_kx_id, axis_ky_id, axis_compl_id, axis_time_id/), .true.)

    if (use_tracer_x) tracerx_var_id = &
                  register_variable(history_file_id, "tracer_x",  &
        (/axis_z_id, axis_kx_id, axis_ky_id, axis_compl_id, axis_time_id/), .true.)

    if (use_tracer_y) tracery_var_id = &
                  register_variable(history_file_id, "tracer_y",  &
        (/axis_z_id, axis_kx_id, axis_ky_id, axis_compl_id, axis_time_id/), .true.)

    time_var_id = register_variable(history_file_id, "time", &
        (/axis_time_id/), .false.)

    call enddef_file(history_file_id)

    call Message("history.nc initialized")

  end subroutine Init_counters
   
  !*********************************************************************

  function Write_snapshots(framein) result(frameout)

    !**************************************************************
    ! Write full snapshot of dynamic field and make restart files
    !**************************************************************

    use io_tools,    only: Message, Write_field
    use qg_params,   only: kmax, nz, nzt, psi_file,       &
                           use_tracer_x, use_tracer_y,    &
                           tracer_x_file, tracer_y_file,  &
                           time, write_time_file, io_root
    use qg_arrays,   only: psi
    use qg_tracers,  only: tracer_x, tracer_y
    use par_tools,   only: par_gather, processor_id, par_sync
    use nc_io_tools, only: write_nc

    integer,intent(in)                    :: framein
    integer                               :: frameout
    complex,dimension(:,:,:),allocatable  :: psi_global, tracer_global

    frameout = framein + 1                ! Update field frame counter

    allocate(psi_global(1:nz,-kmax-1:kmax,0:kmax))
    psi_global=0.
    call par_gather(psi,psi_global,io_root)
    call Write_field(psi_global(:,-kmax:kmax,:),psi_file,frame=frameout,zfirst=1)
    call write_nc(psi_var_id, psi_global(:,-kmax:kmax,:))
    deallocate(psi_global)
   
    if (use_tracer_x.or.use_tracer_y) then
       allocate(tracer_global(1:nzt,-kmax-1:kmax,0:kmax))
       tracer_global = 0.
       if (use_tracer_x) then
          call par_gather(tracer_x,tracer_global,io_root)
          call Write_field(tracer_global(:,-kmax:kmax,:),tracer_x_file, &
                           frame=frameout,zfirst=1) 
          call write_nc(tracerx_var_id, tracer_global(:,-kmax:kmax,:))
       endif
       if (use_tracer_y) then
          call par_gather(tracer_y,tracer_global,io_root)
          call Write_field(tracer_global(:,-kmax:kmax,:),tracer_y_file, &
                           frame=frameout,zfirst=1) 
          call write_nc(tracery_var_id, tracer_global(:,-kmax:kmax,:))
       endif
       deallocate(tracer_global)
   endif

   call Write_field(time,write_time_file,frameout)  
   call write_nc(time_var_id, time)
   call Message('Wrote snapshots, frame: ',tag=frameout)

  end function Write_snapshots

  !*********************************************************************

  function Write_restarts(framein) result(frameout)

    !**************************************************************
    ! Write full snapshot of dynamic field and make restart files
    !**************************************************************

    use io_tools,    only: Message, Write_field
    use qg_params,   only: kmax,nz,nzt,                             &
                           psi_restart_file,time,restart_time_file, &
                           use_tracer_x,tracer_x_restart_file,      &
                           use_tracer_y,tracer_y_restart_file,      &
                           use_forcing,force_o_file,                &
                           Write_parameters, io_root
    use qg_arrays,   only: psi
    use qg_tracers,  only: tracer_x, tracer_y
    use rmf_forcing, only: force_o
    use par_tools,   only: par_gather

    integer,intent(in)                      :: framein
    integer                                 :: frameout
    complex,dimension(:,:),allocatable      :: force_o_global
    complex,dimension(:,:,:),allocatable    :: psi_global, tracer_global

    frameout = framein + 1                ! Update field frame counter

    allocate(psi_global(1:nz,-kmax-1:kmax,0:kmax))
    psi_global=0.
    call par_gather(psi,psi_global,io_root)
    call Write_field(psi_global(:,-kmax:kmax,:),psi_restart_file,zfirst=1)
    deallocate(psi_global)

    if (use_tracer_x.or.use_tracer_y) then
       allocate(tracer_global(1:nzt,-kmax-1:kmax,0:kmax))
       tracer_global = 0.
       if (use_tracer_x) then
          call par_gather(tracer_x,tracer_global,io_root)
          call Write_field(tracer_global(:,-kmax:kmax,:),tracer_x_restart_file,zfirst=1) 
       endif
       if (use_tracer_y) then
          call par_gather(tracer_y,tracer_global,io_root)
          call Write_field(tracer_global(:,-kmax:kmax,:),tracer_y_restart_file,zfirst=1) 
       endif
       deallocate(tracer_global)
    endif

    if (use_forcing) then
       allocate(force_o_global(-kmax-1:kmax,0:kmax))
       force_o_global = 0.
       call par_gather(force_o,force_o_global,io_root)
       call Write_field(force_o_global,force_o_file)
       deallocate(force_o_global)
    endif

    call Write_parameters            ! Write params for restart.nml
!    call Write_field(time,restart_time_file,frame=frameout) ! no need for this 
    call Message('Wrote restarts')

  end function Write_restarts

  !*********************************************************************

                      
end module qg_output
