module follia
!< Module providing the Lagrange Inerpolator Object

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif

public :: lagrange_interpolator

type :: lagrange_interpolator
  !< Object providing the Lagrange Interpolator.
  integer(I_P)           :: S         !< Stencil number and dimension(s).
  real(RPP), allocatable :: x(:)      !< Abscissas of the interpolation points [-S+1:S-1].
  real(RPP), allocatable :: y(:)      !< Function values at the interpolation points [-S+1:S-1].
  real(RPP), allocatable :: coef(:,:) !< Lagrange coefficients of the interpolation [1:S,1:S].
  real(RPP)              :: x_target  !< Abscissa where the function must be interpolated.
  real(RPP),allocatable  :: interp(:) !< Interpolated value(s) [1:S].
  contains
    procedure, pass(self) :: allocate_interpolator     !< Allocate the interpolator object.
    procedure, pass(self) :: deallocate_interpolator   !< Deallocate the interpolator object.
    procedure, pass(self) :: initialize_interpolator   !< Initialize the interpolator.
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
  self%x        = 0.0_RPP
  self%y        = 0.0_RPP
  self%coef     = 0.0_RPP
  self%interp   = 0.0_RPP
  self%x_target = 0.0_RPP
  endsubroutine allocate_interpolator

  subroutine deallocate_interpolator(self)
  !< Dellocate the interpolator.
  class(lagrange_interpolator), intent(inout) :: self   !< Lagrange interpolator.

  if (allocated(self%x     )) deallocate(self%x     )
  if (allocated(self%y     )) deallocate(self%y     )
  if (allocated(self%coef  )) deallocate(self%coef  )
  if (allocated(self%interp)) deallocate(self%interp)
  endsubroutine deallocate_interpolator

  subroutine initialize_interpolator(self, S, x_target, absc)
  !< Initialize the interpolator
  class(lagrange_interpolator), intent(inout) :: self      !< Lagrange interpolator.
  integer(I_P)                , intent(in)    :: S      !< Stencil number and dimension(s).
  real(RPP)                   , intent(in)    :: x_target  !< Abscissa where interpolation takes place.
  real(RPP), optional         , intent(in)    :: absc(:)   !< Abscissas of the interpolation points[-S+1:S-1].
  real(RPP)                                   :: ubound_x  !< Upper bound for x_target [(x(0)+x(+1))/2].
  real(RPP)                                   :: lbound_x  !< Lower bound for x_target [(x(0)+x(-1))/2].
  integer(I_P)                                :: i         !< Counter.

  call self%allocate_interpolator(S)
  associate(x=>self%x,x_tar=>self%x_target)
    if (present(absc)) then
      x = absc
    else
      do i=1,2*S-1  !< Whole interpolation domain loop
        x(-S+i) = -S + i
      enddo
    endif
    ubound_x = (x(0) + x( 1)) / 2._RPP
    lbound_x = (x(0) + x(-1)) / 2._RPP
    if (x_target>ubound_x) then
      error stop 'The abscissa of the target point is bigger than the upper bound of the cell'
    elseif (x_target<lbound_x) then
      error stop 'The abscissa of the target point is lower than the upper bound of the cell'
    elseif (x_target==x(0)) then
      error stop 'The abscissa of the target point corresponds to the center of the cell'
    else
      x_tar = x_target
    endif
  endassociate
  endsubroutine initialize_interpolator

  subroutine compute_coefficients(self)
  !< Compute the coefficients of the Lagrange polynomial(s).
  class(lagrange_interpolator), intent(inout) :: self    !< Lagrange interpolator.
  real(RPP)                                   :: prod    !< Temporary variable.
  integer(I_P)                                :: i, j ,k !< Counters.

  associate(S=>self%S,x=>self%x,coef=>self%coef,x_tar=>self%x_target)
    do k=1,S  !stencils loop
      do j=1,S  !values loop
        prod = 1._RPP
        do i=0,S-1
          if (-S+k+j-1==-S+k+i) cycle
          prod = prod * ((x_tar - x(-S+k+i)) / (x(-S+k+j-1) - x(-S+k+i)))
        enddo
        coef(j,k) = prod
        if (prod==0._RPP) then
          print *, i, j
          do i=0,S-1
            if (-S+k+j-1==-S+k+i) cycle
            print *, x(-S+k+i), x(-S+k+j-1)
          enddo
        endif
      enddo
    enddo
  endassociate
  endsubroutine compute_coefficients

  subroutine compute_interpolations(self)
  !< Compute the Lagrange interpolation(s).
  class(lagrange_interpolator), intent(inout) :: self   !< Lagrange interpolator.
  real(RPP)                                   :: inter  !< Temporary variable.
  integer(I_P)                                :: i, j   !< Counters.

  associate(S=>self%S,y=>self%y,coef=>self%coef,interp=>self%interp)
    do j=1,S  !stencils loop
      inter = 0.0_RPP
      do i=1,S  !values loop
        inter = inter + coef(i,j) * y(-S+i+j-1)
      enddo
      interp(j) = inter
    enddo
  endassociate
  endsubroutine compute_interpolations
endmodule follia
