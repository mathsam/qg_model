module qg_advect

 !******************************************************************************
 ! Calculating -J(psi, zeta) using full nonlinear dynamics or quasi-linear 
 ! dynamics without eddy-eddy interactions
 ! The quasi-linear dynamics is
 ! \partial\zeta/\partial t=-\left[v'\partial\overline{\zeta}/\partial y
 !                          +\overline{u}\partial\zeta'/\partial x
 !                          +\overline{J(\psi',\zeta')}\right]
 !                          +...
 !******************************************************************************

 use qg_params, only: nz, ny, ngrid, turn_off_eddy_eddy
 implicit none
 private

 complex,dimension(:,:,:),allocatable :: u_eddy, qyg_eddy
 complex,dimension(:,:),allocatable   :: u_mean, qyg_mean !shape is (nz,ny)
 save

 public:: init_qg_advect, get_advection

contains

  subroutine init_qg_advect()
    use io_tools,    only: Message

    if(turn_off_eddy_eddy) then
      allocate(u_eddy(nz, ny, ngrid))
      allocate(qyg_eddy(nz, ny, ngrid))
      allocate(u_mean(nz,ny))
      allocate(qyg_mean(nz,ny))
      call Message('use quasi-linear dynamics: eddy-eddy interactions removed')
    endif
  end subroutine init_qg_advect


  function get_advection(ug, vg, qxg, qyg) result (minus_j_psi_vor)
    use transform_tools, only: ir_prod, ir_prod_eddy_mean, ir_prod_mean_eddy
    complex,dimension(:,:,:),intent(in) :: ug, vg, qxg, qyg
    complex,dimension(size(ug,1),size(ug,2),size(ug,3)) :: minus_j_psi_vor
    integer :: ix 

    if(.NOT. turn_off_eddy_eddy) then
        minus_j_psi_vor = -ir_prod(ug,qxg) - ir_prod(vg,qyg)
    else
        u_mean = SUM(ug, 3)/ngrid
        qyg_mean = SUM(qyg, 3)/ngrid
        do ix = 1, ngrid
            u_eddy(:,:,ix) = ug(:,:,ix) - u_mean
            qyg_eddy(:,:,ix) = qyg(:,:,ix) - qyg_mean
        enddo

        minus_j_psi_vor(:,:,1) = SUM(-ir_prod(u_eddy,qxg) &
                                     -ir_prod(vg,qyg_eddy), 3)/ngrid
        do ix = 2, ngrid
            minus_j_psi_vor(:,:,ix) = minus_j_psi_vor(:,:,1)
        enddo

        minus_j_psi_vor = minus_j_psi_vor - ir_prod_eddy_mean(vg,qyg_mean) &
                                          - ir_prod_mean_eddy(u_mean,qxg)
    endif 

  end function get_advection

end module qg_advect
