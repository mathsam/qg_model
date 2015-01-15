module qg_sponge

 !******************************************************************************
 ! Add sponges at the walls of the channel (boundaries in y-direction)
 !
 !******************************************************************************

 use qg_params, only: nz, ny, ngrid
 implicit none
 private

 ! mask that contains coefficients to damp potential vorticity
 ! qg: potential vorticity on physical grid
 ! Dqg/Dt = -sponge_mask*qg
 ! real and imaginary for non-staggered and staggered parts
 complex,dimension(:,:,:),allocatable :: sponge_mask
 save

 public:: init_qg_sponge, apply_sponge

contains

  subroutine init_qg_sponge()
    use qg_params,    only: sponge_rate, sponge_width,     &
                            y_start, y_end,                &
                            pi
    use io_tools,     only: Message
    integer :: y, iy, ix
    real    :: lower_sponge, upper_sponge, lat  ! within [0, 2pi]


    ! note: ny is extent of y-dimension (dim 1) of LOCAL grid field
    ! ngrid = 2*(kmax+1) Physical resolution in x

    allocate(sponge_mask(nz, ny, ngrid)) 
    sponge_mask = 0.0

    do y = y_start, y_end
        iy = y - y_start + 1
        do ix = 1, ngrid
           ! non-staggered
           lat = (DBLE(y)-1.)/DBLE(ngrid)*2.*pi   
           lower_sponge = sponge_rate*(1. - lat/sponge_width)
           upper_sponge = sponge_rate*(lat - 2.*pi + sponge_width)/sponge_width
           sponge_mask(:,iy,ix) = max(0., lower_sponge, upper_sponge)
           ! staggered
           lat = (DBLE(y)-0.5)/DBLE(ngrid)*2.*pi             
           lower_sponge = sponge_rate*(1. - lat/sponge_width)
           upper_sponge = sponge_rate*(lat - 2.*pi + sponge_width)/sponge_width
           sponge_mask(:,iy,ix) = sponge_mask(:,iy,ix) + &
                                  cmplx(0., 1.) * max(0., lower_sponge, upper_sponge)
        enddo
    enddo
    call Message('qg_sponge initialized')
  end subroutine init_qg_sponge


  function apply_sponge(q_in) result (rhs_out)
    ! q is on physical grid. rhs_out is also on physical grid
    ! return -q_in*sponge_mask
    use transform_tools, only: spec2grid, ir_prod
    complex,dimension(:,:,:),intent(in)                       :: q_in
    complex,dimension(size(q_in,1),size(q_in,2),size(q_in,3)) :: rhs_out

    rhs_out = - ir_prod(sponge_mask, q_in)
  end function apply_sponge

end module qg_sponge
