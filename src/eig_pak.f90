module eig_pak                     !-*-f90-*-  <- tells emacs: use f90-mode

  !********************************************************************
  !
  ! F90 wrapper for routine eigrg1.f, which computes eigenvalues
  ! of general complex nonsymmetric matrix.  
  !
  ! USAGE:  CALL eig(array, eval, evec_l=evec_l, evec_r=evec_r, &
  !                  info=info)
  !
  ! {note that you must use explicit argument calls for last three 
  ! optional arguments of form: arg = <your choice of argument>}
  !
  ! ARGUMENTS: 
  !
  ! array:   (input) complex 2D array of size NxN - values are 
  !          preserved.
  !
  ! eval:    (output) 1D complex array of size N containing eigenvalues.
  !
  ! evec:    (output,optional) 2D complex array of size NxN containing 
  !          eigenvectors as columns: v(j) = evec(:,j) is the 
  !          j-th eigenvector, corresponding to eval(j).
  !
  ! info:    (output,optional) INTEGER 
  !          = 0:  successful EXIT
  !          32+j: IF INFO = 32+j, the j-th eigenvalue could not be computed
  !                after 30 iterations
  !
  !**********************************************************************

  implicit none

contains

  subroutine eig(array,eval,evec,info)

    real,dimension(:,:),intent(in)               :: array
    real,dimension(:),intent(out),optional       :: eval
    real,dimension(:,:),intent(out),optional     :: evec
    integer,intent(out),optional                 :: info

    ! Local

    real,dimension(:,:),allocatable :: work1,evecd
    real,dimension(:),allocatable   :: work2
    real,dimension(:),allocatable   :: wr, wi
    integer                         :: n, ierr

    n = size(array,1)
    allocate(work1(n,n),work2(2*n+2),wr(n),wi(n),evecd(n,n))

    call eigrg1(n,n,array,n,wr,wi,evecd,work1,work2,ierr)

    if (present(eval)) eval = wr
    if (present(info)) info = ierr
    if (present(evec)) evec = evecd

    deallocate(work1,work2,wr,wi,evecd)

  end subroutine eig

end module eig_pak
