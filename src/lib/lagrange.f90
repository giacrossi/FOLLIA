program lagrange

  implicit none

  real(R_P), allocatable :: coef(:,:)
  real(R_P), allocatable :: x(:)
  real(R_P) :: x_tar, prod
  integer(I_P) :: S, i, j, k

  write(*,*) 'Insert stencil number'
  read(*,*) S

  allocate(coef(0:S-1,0:S-1))
  coef(:,:) = 0._R_P
  allocate(x(-S+1:S-1))

  write(*,*) 'Insert target coordinate'
  read(*,*) x_tar
  write(*,*) 'Insert stencil coordinates'
  do i=-S+1,S-1
    read(*,*) x(i)
  enddo

  do k=0,S-1 !stencils loop
    do j=0,S-1 !values loop
      prod = 1._R_P
      do i=0,S-1
        if (i==j) cycle
        prod = prod * ((x_tar - x(1-S+k+i)) / (x(1-S+k+j) - x(1-S+k+i)))
      enddo
      coef(j,k) = prod
    enddo
  enddo

  print *, coef(:,:)

endprogram lagrange
