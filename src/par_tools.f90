module par_tools                   !-*-f90-*-

  !********************************************************************
  ! Contains global parameters and routines for using MPI
  !
  ! Routines:  par_sync, par_bcast, par_allgather, par_swapdata, par_sum, 
  !            par_scatter, par_gather, init_par, end_par
  !
  ! Dependencies: mpif.h, and must be linked with mpi libraries  
  !
  !********************************************************************

  implicit none
  private
  save

  integer                 :: processor_id = 0    ! Default values for 1 pe
  integer                 :: num_processors = 1 
  integer                 :: ierr

  include "mpif.h" 

  public :: par_bcast, par_allgather, par_swapdata, par_sum, par_scatter, par_gather, &
            par_sync, init_par, end_par, processor_id, num_processors, MPI_COMM_WORLD

  interface par_bcast
     module procedure bcast_int0, bcast_int1, bcast_char
     module procedure bcast_real0, bcast_real1, bcast_real2, bcast_real3
  end interface
  interface par_allgather
     module procedure allgather_int1
  end interface
  interface par_swapdata
     module procedure swapdata_complex1, swapdata_complex2, swapdata_complex3
  end interface
  interface par_sum
     module procedure sum_real0, sum_real1, sum_real2
     module procedure sum_int0
  end interface
  interface par_scatter
     module procedure scatter_complex2, scatter_complex3
  end interface
  interface par_gather
     module procedure gather_complex2, gather_complex3
  end interface

 contains

   !*************************************************************************

   subroutine init_par

     ! Initialize MPI, get processor id for current processor, and
     ! get total number of processors.

     call mpi_init(ierr)
     call mpi_comm_rank(MPI_COMM_WORLD,processor_id,ierr)
     call mpi_comm_size(MPI_COMM_WORLD,num_processors,ierr)

   end subroutine init_par

   !*************************************************************************

   subroutine end_par

     ! Finalize MPI

     call mpi_finalize(ierr)

   end subroutine end_par

   !*************************************************************************

   subroutine par_sync

     ! Synchronize processes
     
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

   end subroutine par_sync

   !*************************************************************************

   subroutine bcast_char(c,p)
     
     ! Broadcast a char string from processor p to all other processors.

     character(*),intent(in)  :: c
     integer, intent(in)      :: p     
    
     call MPI_Bcast(c,len(c),MPI_CHARACTER,p,MPI_COMM_WORLD,ierr)

   end subroutine bcast_char

   !*************************************************************************

   subroutine bcast_int0(val,p)

     ! Broadcast a scalar integer from processor p to all other processors.
     
     integer,intent(in)     :: val
     integer, intent(in)    :: p     
    
     call MPI_Bcast(val,1,MPI_INTEGER,p,MPI_COMM_WORLD,ierr)

   end subroutine bcast_int0

   !*************************************************************************

   subroutine bcast_int1(arr,p)
     
     ! Broadcast a 1d integer array from processor p to all other processors.
     
     integer,dimension(:),intent(in)  :: arr
     integer, intent(in)              :: p     
    
     call MPI_Bcast(arr,size(arr),MPI_INTEGER,p,MPI_COMM_WORLD,ierr)

   end subroutine bcast_int1

   !*************************************************************************

   subroutine bcast_real0(val,p)
     
     ! Broadcast a scalar real var from processor p to all other processors.
     
     real,intent(in)        :: val
     integer, intent(in)    :: p     
    
     call MPI_Bcast(val,1,MPI_REAL8,p,MPI_COMM_WORLD,ierr)

   end subroutine bcast_real0

   !*************************************************************************

   subroutine bcast_real1(arr,p)
     
     ! Broadcast a 1d real array from processor p to all other processors.
     
     real,dimension(:),intent(in)  :: arr
     integer, intent(in)           :: p     
    
     call MPI_Bcast(arr,size(arr),MPI_REAL8,p,MPI_COMM_WORLD,ierr)

   end subroutine bcast_real1

   !*************************************************************************

   subroutine bcast_real2(arr,p)
     
     ! Broadcast a 2d real array from processor p to all other processors.
     
     real,dimension(:,:),intent(in)  :: arr
     integer, intent(in)             :: p     
    
     call MPI_Bcast(arr,size(arr),MPI_REAL8,p,MPI_COMM_WORLD,ierr)

   end subroutine bcast_real2

   !*************************************************************************

   subroutine bcast_real3(arr,p)
     
     ! Broadcast a 3d real array from processor p to all other processors.
     
     real,dimension(:,:,:),intent(in)  :: arr
     integer, intent(in)               :: p     
    
     call MPI_Bcast(arr,size(arr),MPI_REAL8,p,MPI_COMM_WORLD,ierr)

   end subroutine bcast_real3

   !*************************************************************************

   subroutine allgather_int1(fin,fout)

     ! Concatenate data from send bufferes on all processors into receive buffers
     ! of size num_processors*size(send buffer) on each processor

     integer,dimension(:),intent(in)   :: fin
     integer,dimension(size(fin)*num_processors),intent(out)  :: fout
     integer                           :: count

     count = size(fin)
     call mpi_allgather(fin,  count, MPI_INTEGER, &
                        fout, count, MPI_INTEGER, &
                        MPI_COMM_WORLD, ierr)

   end subroutine allgather_int1

   !*************************************************************************
   
   subroutine swapdata_complex1(f,p)

     ! Implementation of MPI_SENDRECV_REPLACE.  Trade 1d complex array f 
     ! between current processor and processor p.  

     complex,dimension(:),intent(inout)  :: f
     integer, intent(in)                 :: p     
     integer                             :: tag = 1
     integer, dimension(MPI_STATUS_SIZE) :: status

     call mpi_sendrecv_replace(f, size(f), MPI_DOUBLE_COMPLEX, p, tag, p, tag, &
                               MPI_COMM_WORLD, status, ierr)


   end subroutine swapdata_complex1

   !*************************************************************************
   
   subroutine swapdata_complex2(f,p)

     ! Implementation of MPI_SENDRECV_REPLACE.  Trade 2d complex array f
     ! between current processor and processor p

     complex,dimension(:,:),intent(inout)  :: f
     integer, intent(in)                   :: p     
     integer                               :: tag = 1
     integer, dimension(MPI_STATUS_SIZE)   :: status

     call mpi_sendrecv_replace(f, size(f), MPI_DOUBLE_COMPLEX, p, tag, p, tag, &
                               MPI_COMM_WORLD, status, ierr)


   end subroutine swapdata_complex2

   !*************************************************************************
   
   subroutine swapdata_complex3(f,p)

     ! Implementation of MPI_SENDRECV_REPLACE.  Trade  3d complex array f
     ! between current processor and processor p

     complex,dimension(:,:,:),intent(inout)  :: f
     integer, intent(in)                     :: p     
     integer                                 :: tag = 1
     integer, dimension(MPI_STATUS_SIZE)     :: status

     call mpi_sendrecv_replace(f, size(f), MPI_DOUBLE_COMPLEX, p, tag, p, tag, &
                               MPI_COMM_WORLD, status, ierr)


   end subroutine swapdata_complex3

   !*************************************************************************

   subroutine sum_int0(f)

     ! Sum integer scalar f over all processors.

     integer, intent(inout) :: f
     integer                :: f_summed
     
     call mpi_allreduce(f, f_summed, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr) 

     f = f_summed

   end subroutine sum_int0
   
   !*************************************************************************

   subroutine sum_real0(f)

     ! Sum real scalar f over all processors.

     real, intent(inout) :: f
     real                :: f_summed
     
     call mpi_allreduce(f, f_summed, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr) 

     f = f_summed

   end subroutine sum_real0
   
   !*************************************************************************

   subroutine sum_real1(f)

     ! Sum 1d real array f over all processors.

     real, dimension(:), intent(inout) :: f
     real, dimension(size(f))          :: f_summed
     
     call mpi_allreduce(f, f_summed, size(f), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr) 

     f = f_summed

   end subroutine sum_real1
   
   !*************************************************************************

   subroutine sum_real2(f)

     ! Sum 2d real array f over all processors.

     real, dimension(:,:), intent(inout)  :: f
     real, dimension(size(f,1),size(f,2)) :: f_summed
     
     call mpi_allreduce(f, f_summed, size(f), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr) 

     f = f_summed

   end subroutine sum_real2
   
   !*************************************************************************
   
   subroutine gather_complex2(fin,fout,root)

     ! Like allgather, but data in buffers fin are only concatenated on processor 
     ! 'root'.  For 2d complex array

     complex,dimension(:,:),intent(in)             :: fin
     complex,dimension(:,:),intent(inout)          :: fout
     integer,intent(in)                            :: root
     complex,dimension(size(fin,2),size(fin,1))    :: temp_in
     complex,dimension(size(fout,2),size(fout,1))  :: temp_out

     ! We need y to be contiguous in memory, so transpose to make x the
     ! 2nd dimension.
     temp_in = transpose(fin)

     CALL mpi_gather(temp_in,  size(temp_in), MPI_DOUBLE_COMPLEX, &
                     temp_out, size(temp_in), MPI_DOUBLE_COMPLEX,  &
                     root, MPI_COMM_WORLD , ierr)

     ! On root proc (0), undo transpose
     if (processor_id==root) then
        fout = transpose(temp_out)
     endif

   end subroutine gather_complex2

   !*************************************************************************

    subroutine gather_complex3(fin,fout,root)

     ! Like allgather, but data in buffers fin are only concatenated on processor 
     ! 'root'.  For 3d complex array

     complex,dimension(:,:,:),intent(in)                       :: fin
     complex,dimension(:,:,:),intent(inout)                    :: fout
     integer,intent(in)                                        :: root
     complex,dimension(size(fin,1),size(fin,3),size(fin,2))    :: temp_in
     complex,dimension(size(fout,1),size(fout,3),size(fout,2)) :: temp_out
     integer                                                   :: iz

     ! We need z and y to be contiguous in memory, so transpose to make x the
     ! 3rd dimension.

     do iz = 1,size(fin,1)
        temp_in(iz,:,:) = transpose(fin(iz,:,:))
     enddo

     CALL mpi_gather(temp_in,  size(temp_in), MPI_DOUBLE_COMPLEX, &
                     temp_out, size(temp_in), MPI_DOUBLE_COMPLEX,  &
                     root,MPI_COMM_WORLD , ierr)

     ! On root proc, undo transpose
     if (processor_id==root) then
        do iz = 1,size(fout,1)
           fout(iz,:,:) = transpose(temp_out(iz,:,:))
        enddo
     endif

   end subroutine gather_complex3

   !*************************************************************************

   subroutine scatter_complex2(fin,fout,root)
     
     ! Scatter data in array fin on processor root to all processors.
     ! Not the same as broadcast because data in fin is broken up into
     ! segments of size fout, and each part is sent to only one particular
     ! processor.  For 2d complex array

     complex,dimension(:,:),intent(in)              :: fin
     complex,dimension(:,:),intent(inout)           :: fout
     integer,intent(in)                             :: root
     complex,dimension(size(fin,2),size(fin,1))     :: temp_in
     complex,dimension(size(fout,2),size(fout,1))   :: temp_out
 
     if (processor_id==0) then
        temp_in = transpose(fin)
     endif

     CALL mpi_scatter(temp_in,  size(temp_out), MPI_DOUBLE_COMPLEX, &
                      temp_out, size(temp_out), MPI_DOUBLE_COMPLEX, &
                      root, MPI_COMM_WORLD, ierr)

     fout = transpose(temp_out)

   end subroutine scatter_complex2

   !*************************************************************************

   subroutine scatter_complex3(fin,fout,root)
     
     ! Scatter data in array fin on processor root to all processors.
     ! Not the same as broadcast because data in fin is broken up into
     ! segments of size fout, and each part is sent to only one particular
     ! processor.  For 3d complex array

      complex,dimension(:,:,:),intent(in)                          :: fin
     complex,dimension(:,:,:),intent(inout)                       :: fout
     integer,intent(in)                                           :: root
     complex,dimension(size(fin,1),size(fin,3),size(fin,2))       :: temp_in
     complex,dimension(size(fout,1),size(fout,3),size(fout,2))    :: temp_out
     integer                                                      :: iz
 
     if (processor_id==0) then
        do iz = 1,size(fin,1)
           temp_in(iz,:,:) = transpose(fin(iz,:,:))
        enddo
     endif

     CALL mpi_scatter(temp_in,  size(temp_out), MPI_DOUBLE_COMPLEX, &
                      temp_out, size(temp_out), MPI_DOUBLE_COMPLEX, &
                      root, MPI_COMM_WORLD, ierr)

     do iz = 1,size(fout,1)
        fout(iz,:,:) = transpose(temp_out(iz,:,:))
     enddo

   end subroutine scatter_complex3

   !*************************************************************************

 end module par_tools
