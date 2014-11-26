module nc_io_tools

  !*************************************************************************
  ! Output model results/Read restart files in NetCDF format
  ! Intents to store most information inside this module
  !
  ! Intents for serial I/O, therefore, only root node can do I/O. And I/O
  ! should be done for global domain!
  ! 
  ! Avialable dimensions are
  !   kx, ky, z, time, real_and_imag, lon, lat
  ! The last dimension must be time. And if variable is complex, the second
  ! to last dimension must be real_and_imag
  ! In the example below, psi in the model takes 5 dimensions to output
  ! For another example, total energy in qg_diagnostics takes 2 dimensions
  !
  ! Sample usage:
  ! only need to store var_id in the calling modules, and other variables 
  ! can/should be kept local
  ! in this paralleled version, fields are stored differently compared with
  ! the serial version
  ! for example, psi(1:nz, kx_start:kx_end, 0:kmax) in this version
  !
  ! Note that this paralleled version cannot do surface qg. And psi in
  ! calculation is always 3 dimensional array in the code
  !
  !   integer :: file_id, var_psi_id
  !   integer :: axis_kx_id, axis_ky_id, axis_z_id, axis_time_id, axis_compl_id
  !
  !   file_id = create_file("history.nc")
  !   axis_time_id = create_axis_time(file_id)
  !   axis_kx_id   = create_axis_kx(file_id)
  !   axis_ky_id   = create_axis_ky(file_id)
  !   axis_z_id    = create_axis_z(file_id)
  !   axis_compl_id= create_axis_real_and_imag(file_id)
  ! 
  !   var_psi_id = register_variable(file_id,"psi", &
  !       (/axis_z_id, axis_kx_id, axis_ky_id, axis_compl_id, axis_time_id/), .true.)
  !   call enddef_file(file_id)
  !  
  !   call write_nc(var_psi_id, psi_global)
  !   
  !   call close_all_files()
  !
  ! Algorithm:
  !  Use file_info_list to keep track of all the NetCDF files that has been
  !  opened. current_file_num is initially set to 1, and is increased by 1
  !  each time a new file is opened. current_file_num -1 is the total number 
  !  of files that has been opened. 
  !  Use var_info_list to keep track of all variable information in a
  !  similar way.
  !  nc_file_id is returned by nf90_create.
  !*************************************************************************

  use netcdf
  implicit none
  private

  type t_var_info
    integer:: framenum   = 1
    integer:: nc_var_id     = -1
    integer:: nc_file_id    = -1
    integer:: num_dim    = 1 ! e.g. (x,y,z,real_and_imag,t), dimension is 5
    logical:: is_complex = .FALSE.
    integer, allocatable :: var_start_real(:)
    integer, allocatable :: var_start_imag(:)
    integer, allocatable :: var_count(:)
  end type t_var_info

  integer, parameter :: MAX_NUM_VARS = 50
  type(t_var_info)   :: var_info_list(MAX_NUM_VARS)
  integer            :: current_var_num = 1

  integer, parameter :: NAN = -9999
  integer, parameter :: MAX_NUM_DIMS = 5 
  type t_file_info
    integer                          :: nc_file_id
    integer, dimension(MAX_NUM_DIMS) :: axis_id_list = NAN
    integer, dimension(MAX_NUM_DIMS) :: output_length_list = NAN
    integer                          :: current_num_dims = 1
  end type t_file_info

  integer, parameter :: MAX_FILES = 10
  type(t_file_info)  :: file_info_list(MAX_FILES)
  integer            :: current_file_num = 1

  integer :: my_pe, io_root
  integer :: kmax, nz ! same as that in qg_params
  save

  public :: pass_params_nc_io, create_file, enddef_file,             &
            close_all_files,                                         &
            register_variable, write_nc,                             &
            create_axis_time, create_axis_kx, create_axis_ky,        &
            create_axis_z, create_axis_real_and_imag,                &
            create_axis_lon, create_axis_lat

  interface write_nc
    module procedure write_nc_0d, write_nc_0d_complex, &
                     write_nc_1d, write_nc_1d_complex, &
                     write_nc_2d, write_nc_2d_complex, &
                     write_nc_3d, write_nc_3d_complex
  end interface

  interface read_nc
    module procedure read_nc_0d, read_nc_3d_complex
  end interface

contains

  subroutine pass_params_nc_io(processor_id_in, io_root_in, kmax_in, nz_in)
    integer, intent(in):: processor_id_in, io_root_in, kmax_in, nz_in
    my_pe   = processor_id_in
    io_root = io_root_in
    kmax    = kmax_in
    nz      = nz_in
  end subroutine pass_params_nc_io

  subroutine handle_err(status, err_message)
     integer, intent ( in) :: status
     character(*), intent (in), optional :: err_message

     if(present(err_message)) print *, err_message
 
     if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       print *, "Netcdf error"
       stop "Stopped"
     end if
  end subroutine handle_err

!*************************************************************************
! Create a file for output and returns an integer identifier
! Notice that the identifier is not the id returned from nf90_create
! But it is the index of the file_info_list, which stores the file id
! together with the axis id returned from nf90_...
!*************************************************************************
  integer function create_file(file_name)
    character(*), intent(in) :: file_name
    integer :: file_id, status

    if (my_pe /= io_root) then
        create_file = -1
        return
    endif

    if (current_file_num > MAX_FILES) then
        print *, "Number of output files exceed maximum allowed"
        stop "Stopped"
    endif

    status = nf90_create(path = file_name, cmode = nf90_clobber,&
                         ncid = file_id)
    if (status /= nf90_noerr) call handle_err(status, "Create file error at "//file_name)
    file_info_list(current_file_num)%nc_file_id = file_id
    create_file = current_file_num
    current_file_num = current_file_num + 1
    return
  end function create_file

  subroutine enddef_file(file_id)
    integer, intent(in) :: file_id
    integer :: status, nc_file_id

    if (my_pe /= io_root) then
         return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id
    status = nf90_enddef(nc_file_id)
    if (status /= nf90_noerr) call handle_err(status, "End def error")
  end subroutine enddef_file

  subroutine close_all_files()
    integer :: i
    integer :: status, nc_file_id

    if (my_pe /= io_root) then
        return
    endif

    do i = 1, current_file_num - 1
        nc_file_id = file_info_list(i)%nc_file_id
        status = nf90_close(nc_file_id)
        if (status /= nf90_noerr) call handle_err(status)
    enddo

    do i = 1, current_var_num - 1
        deallocate(var_info_list(i)%var_start_real)
        deallocate(var_info_list(i)%var_start_imag)
        deallocate(var_info_list(i)%var_count)
    enddo

    !after these operations, file_info_list and var_info_list can
    !be used to keep track of files/vars again
    current_file_num = 1
    current_var_num  = 1
  end subroutine close_all_files

  integer function register_variable(file_id_in,var_name,axis_id_list,is_complex)
    integer, intent(in)               :: file_id_in
    character(*), intent(in)          :: var_name
    integer, dimension(:), intent(in) :: axis_id_list
    logical, intent(in)               :: is_complex
    integer, dimension(size(axis_id_list))  :: nc_axis_id_list
    integer :: i, status
    integer :: nc_file_id, nc_var_id
    integer :: var_num_dims

    if (my_pe /= io_root) then
        register_variable = -1
        return
    endif

    if (current_var_num > MAX_NUM_VARS) then
        print *, "Number of variables exceed maximum"
        stop "Stopped"
    endif 

    nc_file_id = file_info_list(file_id_in)%nc_file_id
    do i = 1, size(axis_id_list)
        nc_axis_id_list(i) = file_info_list(file_id_in)%axis_id_list(axis_id_list(i))
    enddo
    status = nf90_def_var(nc_file_id, trim(var_name), nf90_double, &
                          nc_axis_id_list, nc_var_id)
    if (status /= nf90_noerr) call handle_err(status,"error in register "//var_name)

    var_info_list(current_var_num)%nc_var_id  = nc_var_id
    var_num_dims = size(axis_id_list)
    var_info_list(current_var_num)%num_dim = var_num_dims
    var_info_list(current_var_num)%nc_file_id = nc_file_id

    allocate(var_info_list(current_var_num)%var_start_real(var_num_dims))
    allocate(var_info_list(current_var_num)%var_start_imag(var_num_dims))
    allocate(var_info_list(current_var_num)%var_count(var_num_dims))
    var_info_list(current_var_num)%var_start_real(:) = 1
    var_info_list(current_var_num)%var_start_imag(:) = 1
    if(is_complex) then
      var_info_list(current_var_num)%var_start_imag(var_num_dims - 1) = 2
      var_info_list(current_var_num)%is_complex = is_complex
    endif
   
    do i = 1, size(axis_id_list)
        var_info_list(current_var_num)%var_count(i) = &
               file_info_list(file_id_in)%output_length_list(axis_id_list(i))
    enddo

    register_variable = current_var_num 
    current_var_num   = current_var_num + 1
    return
  end function register_variable

  subroutine write_nc_0d(var_id_in, field0d)
    integer, intent(in) :: var_id_in
    real, intent(in)    :: field0d ! scalar
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(1) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, field0d, &
                          start = var_info_list(var_id_in)%var_start_real)
    if (status /= nf90_noerr) call handle_err(status)
    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_0d

  subroutine write_nc_0d_complex(var_id_in, field0d)
    integer, intent(in) :: var_id_in
    complex, intent(in)    :: field0d ! scalar
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(2) = var_info_list(var_id_in)%framenum
    var_info_list(var_id_in)%var_start_imag(2) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, real(field0d), &
                          start = var_info_list(var_id_in)%var_start_real)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(nc_file_id, nc_var_id, aimag(field0d), &
                          start = var_info_list(var_id_in)%var_start_imag)
    if (status /= nf90_noerr) call handle_err(status)

    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_0d_complex

  subroutine write_nc_1d(var_id_in, field1d)
    integer, intent(in)               :: var_id_in
    real, dimension(:), intent(in)    :: field1d ! 1d array, e.g., PV(lat) output everytime 
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(2) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, field1d, &
                          start = var_info_list(var_id_in)%var_start_real, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_1d

  subroutine write_nc_1d_complex(var_id_in, field1d)
    integer, intent(in)                  :: var_id_in
    complex, dimension(:), intent(in)    :: field1d
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(3) = var_info_list(var_id_in)%framenum
    var_info_list(var_id_in)%var_start_imag(3) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, real(field1d), &
                          start = var_info_list(var_id_in)%var_start_real, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(nc_file_id, nc_var_id, aimag(field1d), &
                          start = var_info_list(var_id_in)%var_start_imag, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_1d_complex

  subroutine write_nc_2d(var_id_in, field2d)
    integer, intent(in)               :: var_id_in
    real, dimension(:,:), intent(in)    :: field2d ! 2d array, e.g., PV(lat,lon) output everytime 
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(3) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, field2d, &
                          start = var_info_list(var_id_in)%var_start_real, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_2d

  subroutine write_nc_2d_complex(var_id_in, field2d)
    integer, intent(in)                  :: var_id_in
    complex, dimension(:,:), intent(in)    :: field2d
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(4) = var_info_list(var_id_in)%framenum
    var_info_list(var_id_in)%var_start_imag(4) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, real(field2d), &
                          start = var_info_list(var_id_in)%var_start_real, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(nc_file_id, nc_var_id, aimag(field2d), &
                          start = var_info_list(var_id_in)%var_start_imag, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_2d_complex

  subroutine write_nc_3d(var_id_in, field3d)
    integer, intent(in)               :: var_id_in
    real, dimension(:,:,:), intent(in)    :: field3d ! 3d array, e.g., PV(lat,lon,z) output everytime 
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(4) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, field3d, &
                          start = var_info_list(var_id_in)%var_start_real, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_3d

  subroutine write_nc_3d_complex(var_id_in, field3d)
    integer, intent(in)                    :: var_id_in
    complex, dimension(:,:,:), intent(in)  :: field3d
    integer :: nc_file_id, nc_var_id, status

    if (my_pe /= io_root) then
        return
    endif

    nc_file_id = var_info_list(var_id_in)%nc_file_id
    nc_var_id  = var_info_list(var_id_in)%nc_var_id 

    var_info_list(var_id_in)%var_start_real(5) = var_info_list(var_id_in)%framenum
    var_info_list(var_id_in)%var_start_imag(5) = var_info_list(var_id_in)%framenum

    status = nf90_put_var(nc_file_id, nc_var_id, real(field3d), &
                          start = var_info_list(var_id_in)%var_start_real, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(nc_file_id, nc_var_id, aimag(field3d), &
                          start = var_info_list(var_id_in)%var_start_imag, &
                          count = var_info_list(var_id_in)%var_count)
    if (status /= nf90_noerr) call handle_err(status)
    var_info_list(var_id_in)%framenum = var_info_list(var_id_in)%framenum + 1
  end subroutine write_nc_3d_complex

  subroutine axis_dim_exceed_err()
    print *, "axis dimensions exceed maxmium"
    stop "Stopped"
  end subroutine axis_dim_exceed_err

  integer function create_axis_time(file_id)
    integer, intent(in):: file_id
    integer:: status, time_dimid
    integer:: nc_file_id, current_num_dims

    if (my_pe /= io_root) then
        create_axis_time = -1
        return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id
    status = nf90_def_dim(nc_file_id, "time_step", nf90_unlimited, time_dimid)
    if (status /= nf90_noerr) call handle_err(status,"create time axis error")

    current_num_dims = file_info_list(file_id)%current_num_dims
    if(current_num_dims > MAX_NUM_DIMS) call axis_dim_exceed_err
    file_info_list(file_id)%axis_id_list(current_num_dims) = time_dimid
    file_info_list(file_id)%output_length_list(current_num_dims) = 1
    create_axis_time = current_num_dims
    file_info_list(file_id)%current_num_dims = current_num_dims + 1
    return
  end function create_axis_time

  integer function create_axis_kx(file_id)
    integer, intent(in):: file_id
    integer:: status, kx_dimid, nkx
    integer:: nc_file_id, current_num_dims 

    if (my_pe /= io_root) then
        create_axis_kx = -1
        return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id
    nkx = 2*kmax + 1
    status = nf90_def_dim(nc_file_id, "kx", nkx, kx_dimid)
    if (status /= nf90_noerr) call handle_err(status,"create kx error")

    current_num_dims = file_info_list(file_id)%current_num_dims
    if(current_num_dims > MAX_NUM_DIMS) call axis_dim_exceed_err
    file_info_list(file_id)%axis_id_list(current_num_dims) = kx_dimid
    file_info_list(file_id)%output_length_list(current_num_dims) = nkx
    create_axis_kx = current_num_dims
    file_info_list(file_id)%current_num_dims = current_num_dims + 1
    return
  end function create_axis_kx

  integer function create_axis_ky(file_id)
    integer, intent(in):: file_id
    integer:: status, ky_dimid, nky
    integer:: nc_file_id, current_num_dims 

    if (my_pe /= io_root) then
        create_axis_ky = -1
        return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id
 
    nky = kmax + 1
    status = nf90_def_dim(nc_file_id, "ky", nky, ky_dimid)
    if (status /= nf90_noerr) call handle_err(status,"create ky error")

    current_num_dims = file_info_list(file_id)%current_num_dims
    if(current_num_dims > MAX_NUM_DIMS) call axis_dim_exceed_err
    file_info_list(file_id)%axis_id_list(current_num_dims) = ky_dimid
    file_info_list(file_id)%output_length_list(current_num_dims) = nky
    create_axis_ky = current_num_dims
    file_info_list(file_id)%current_num_dims = current_num_dims + 1
    return
  end function create_axis_ky

  integer function create_axis_z(file_id)
    integer, intent(in):: file_id
    integer:: status, z_dimid
    integer:: nc_file_id, current_num_dims 

    if (my_pe /= io_root) then
        create_axis_z = -1
        return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id
 
    status = nf90_def_dim(nc_file_id, "z", nz, z_dimid)
    if (status /= nf90_noerr) call handle_err(status,"create z error")

    current_num_dims = file_info_list(file_id)%current_num_dims
    if(current_num_dims > MAX_NUM_DIMS) call axis_dim_exceed_err
    file_info_list(file_id)%axis_id_list(current_num_dims) = z_dimid
    file_info_list(file_id)%output_length_list(current_num_dims) = nz
    create_axis_z = current_num_dims
    file_info_list(file_id)%current_num_dims = current_num_dims + 1
    return
  end function create_axis_z

  integer function create_axis_real_and_imag(file_id)
    integer, intent(in):: file_id
    integer:: status, real_and_imag_dimid
    integer:: nc_file_id, current_num_dims 

    if (my_pe /= io_root) then
        create_axis_real_and_imag = -1
        return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id

    status = nf90_def_dim(nc_file_id, "real_and_imag", 2, real_and_imag_dimid)
    if (status /= nf90_noerr) call handle_err(status,"create real_and imag error")

    current_num_dims = file_info_list(file_id)%current_num_dims
    if(current_num_dims > MAX_NUM_DIMS) call axis_dim_exceed_err
    file_info_list(file_id)%axis_id_list(current_num_dims) = real_and_imag_dimid 
    file_info_list(file_id)%output_length_list(current_num_dims) = 1
    create_axis_real_and_imag = current_num_dims
    file_info_list(file_id)%current_num_dims = current_num_dims + 1
    return
  end function create_axis_real_and_imag

  integer function create_axis_lon(file_id)
    integer, intent(in):: file_id
    integer:: status, nx, lon_dimid
    integer:: nc_file_id, current_num_dims 

    if (my_pe /= io_root) then
        create_axis_lon = -1
        return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id

    nx = 2*(kmax + 1)
    status = nf90_def_dim(nc_file_id, "lon", nx, lon_dimid)
    if (status /= nf90_noerr) call handle_err(status,"create lon error")

    current_num_dims = file_info_list(file_id)%current_num_dims
    if(current_num_dims > MAX_NUM_DIMS) call axis_dim_exceed_err
    file_info_list(file_id)%axis_id_list(current_num_dims) = lon_dimid 
    file_info_list(file_id)%output_length_list(current_num_dims) = nx
    create_axis_lon = current_num_dims
    file_info_list(file_id)%current_num_dims = current_num_dims + 1
    return
  end function create_axis_lon

  integer function create_axis_lat(file_id)
    integer, intent(in):: file_id
    integer:: status, ny, lat_dimid
    integer:: nc_file_id, current_num_dims 

    if (my_pe /= io_root) then
        create_axis_lat = -1
        return
    endif

    nc_file_id = file_info_list(file_id)%nc_file_id

    ny = 2*(kmax + 1)
    status = nf90_def_dim(nc_file_id, "lat", ny, lat_dimid)
    if (status /= nf90_noerr) call handle_err(status,"create lat error")

    current_num_dims = file_info_list(file_id)%current_num_dims
    if(current_num_dims > MAX_NUM_DIMS) call axis_dim_exceed_err
    file_info_list(file_id)%axis_id_list(current_num_dims) = lat_dimid 
    file_info_list(file_id)%output_length_list(current_num_dims) = ny
    create_axis_lat = current_num_dims
    file_info_list(file_id)%current_num_dims = current_num_dims + 1
    return
  end function create_axis_lat

  subroutine read_nc_0d(filename, varname, var_out)
    character(*), intent(in) :: filename, varname
    real, intent(inout)      :: var_out
    integer                  :: status, file_id, var_id

    if (my_pe /= io_root) then
        return
    endif

    status = nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = file_id)
    if (status /= nf90_noerr) call handle_err(status,"read file"//filename//" error")

    status = nf90_inq_varid(file_id, varname, var_id)
    if (status /= nf90_noerr) call handle_err(status,"inquire variable " & 
                                          //varname//" in file "//filename//" error")

    status = nf90_get_var(file_id, var_id, var_out)
    if (status /= nf90_noerr) call handle_err(status,"read variable " & 
                                         //varname//" in file "//filename//" error")
  end subroutine read_nc_0d

  subroutine read_nc_3d_complex(filename, varname, var_out)
    character(*), intent(in)                  :: filename, varname
    complex, dimension(:,:,:), intent(inout)  :: var_out
    integer                                   :: status, file_id, var_id
    real, dimension(size(var_out,1), size(var_out,2), size(var_out,3),2) &
                                              :: field_real_and_imag

    if (my_pe /= io_root) then
        return
    endif

    status = nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = file_id)
    if (status /= nf90_noerr) call handle_err(status,"read file"//filename//" error")

    status = nf90_inq_varid(file_id, varname, var_id)
    if (status /= nf90_noerr) call handle_err(status,"inquire variable " &
                                          //varname//" in file "//filename//" error")

    status = nf90_get_var(file_id, var_id, field_real_and_imag)
    if (status /= nf90_noerr) call handle_err(status,"read variable " &
                                         //varname//" in file "//filename//" error")

    var_out = field_real_and_imag(:,:,:,1) + (0.0, 1.0)*field_real_and_imag(:,:,:,2)
  end subroutine read_nc_3d_complex

end module nc_io_tools
