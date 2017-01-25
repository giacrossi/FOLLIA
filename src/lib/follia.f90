module follia
!< Module providing the Lagrange Inerpolator Object

use penf, only : I_P, R_P

type :: lagrange_interpolator
  !< Object providing the Lagrange Interpolator.
  integer(I_P)           :: S         !< Stencil number and dimension(s).
  real(R_P), allocatable :: x(:)      !< Abscissas of the interpolation points [-S+1:S-1].
  real(R_P), allocatable :: y(:)      !< Function values at the interpolation points [-S+1:S-1].
  real(R_P), allocatable :: coef(:,:) !< Lagrange coefficients of the interpolation [1:S,1:S].
  real(R_P)              :: x_target  !< Abscissa where the function must be interpolated.
  real(R_P),allocatable  :: interp(:) !< Interpolated value(s) [1:S].
  contains
    procedure, pass(self) :: allocate_interpolator     !< Allocate the interpolator object.
    procedure, pass(self) :: deallocate_interpolator   !< Deallocate the interpolator object.
    procedure, pass(self) :: init_uniform_grid         !< Initialize a uniform grid.
    procedure, pass(self) :: compute_coefficients      !< Compute the interpolation coefficients.
    procedure, pass(self) :: compute_interpolations    !< Compute the interpolation.
endtype lagrange_interpolator

contains

  subroutine allocate_interpolator(self, S)
  !< Allocate the interpolator using the stencil number.
  class(lagrange_interpolator), intent(inout) :: self   !< Lagrange interpolator.
  integer(I_P)                , intent(in)    :: S      !< Stencil number and dimension(s).

  call self%deallocate_interpolator
  self%S = S
  allocate(self%x     (-S+1:S-1    ))
  allocate(self%y     (-S+1:S-1    ))
  allocate(self%coef  (     1:S,1:S))
  allocate(self%interp(         1:S))
  self%x        = 0.0_R_P
  self%y        = 0.0_R_P
  self%coef     = 0.0_R_P
  self%interp   = 0.0_R_P
  self%x_target = 0.0_R_P
  endsubroutine allocate_interpolator

  subroutine deallocate_interpolator(self)
  !< Dellocate the interpolator.
  class(lagrange_interpolator), intent(inout) :: self   !< Lagrange interpolator.

  if (allocated(self%x     )) deallocate(self%x     )
  if (allocated(self%y     )) deallocate(self%y     )
  if (allocated(self%coef  )) deallocate(self%coef  )
  if (allocated(self%interp)) deallocate(self%interp)
  endsubroutine deallocate_interpolator

  subroutine init_uniform_grid(self)
  !< Initialize a uniform grid for interpolation stencils.
  class(lagrange_interpolator), intent(inout) :: self   !< Lagrange interpolator.
  integer(I_P)                                :: i      !< Counter.

  associate(S=>self%S,x=>self%x)
    do i=1,2*S-1  !< Whole interpolation domain loop
      x(i) = -S + i
    enddo
  endassociate
  endsubroutine init_uniform_grid

  subroutine compute_coefficients(self)
  !< Compute the coefficients of the Lagrange polynomial(s).
  class(lagrange_interpolator), intent(inout) :: self    !< Lagrange interpolator.
  real(R_P)                                   :: prod    !< Temporary variable.
  integer(I_P)                                :: i, j ,k !< Counters.

  associate(S=>self%S,x=>self%x,coef=>self%coef,x_tar=>self%x_target)
    do k=1,S  !stencils loop
      do j=1,S  !values loop
        prod = 1._R_P
        do i=0,S-1
          if (i==j) cycle
          prod = prod * ((x_tar - x(-S+k+i)) / (x(-S+k+j-1) - x(-S+k+i)))
        enddo
        coef(j,k) = prod
      enddo
    enddo
  endassociate
  endsubroutine compute_coefficients

  subroutine compute_interpolations(self)
  !< Compute the Lagrange interpolation(s).
  class(lagrange_interpolator), intent(inout) :: self    !< Lagrange interpolator.
  real(R_P)                                   :: inter   !< Temporary variable.
  integer(I_P)                                :: i, j ,k !< Counters.

  associate(S=>self%S,y=>self%y,coef=>self%coef,interp=>self%interp)
    do j=1,S  !stencils loop
      inter = 0.0_R_P
      do i=1,S  !values loop
        inter = inter + coef(i,j) * y(-S+i+j-1)
      enddo
      interp(j) = inter
    enddo
  endassociate
  endsubroutine compute_interpolations
endmodule follia
