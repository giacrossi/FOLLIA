module interpolation
use flap, only : command_line_interface
use penf, only : I_P, R_P, FR_P, str, strz
use pyplot_module, only :  pyplot
use follia, only : lagrange_interpolator

implicit none
private
public :: interpolate

type :: interpolate
  !< Class to handle interpolation(s).
  !<
  !< Interpolation is driven by the Command Line Interface (CLI) options.
  !<
  !< Interpolation has only 1 public method `execute`: it executes intrpolation(s) accordingly to cli options.
  private
  type(command_line_interface)     :: cli                     !< Command line interface handler.
  integer(I_P)                     :: error=0                 !< Error handler.
  integer(I_P)                     :: S_number                !< Number of different stencils used.
  character(99)                    :: interf='r'              !< Interface evaluation.
  logical                          :: unif_grid=.false.       !< Flag for activating uniform grid.
  logical                          :: coeff=.false.           !< Flag for activating Lagrange coefficients saving.
  logical                          :: inter=.false.           !< Flag for activating Lagrange inverpolation evalutation.
  contains
    ! public methods
    procedure, pass(self) :: execute !< Execute interpolation(s).
    ! private methods
    procedure, pass(self), private :: initialize                 !< Initialize interpolation(s).
    procedure, pass(self), private :: perform                    !< Perform interpolation(s).
    procedure, pass(self), private :: save_results               !< Save Lagrange interpolation result(s).
endtype interpolate

contains

  !< Public methods
  subroutine execute(self, interpolator)
  !< Execute the interpolation.
  class(interpolate),           intent(inout) :: self !< Lagrange interpolator.
  class(lagrange_interpolator), intent(inout) :: interpolator

  call self%initialize
  call self%perform(interpolator)
  endsubroutine execute

  !< Private methods
  subroutine initialize(self)
  !< Initialize interpolation, set Command Line Interface, parse it and check its validity.
  class(interpolate), intent(inout) :: self !< Lagrange interpolator.

  call set_cli
  call parse_cli
  contains
    subroutine set_cli()
    !< Set Command Line Interface.

    associate(cli => self%cli)
      call cli%init(progname    = 'Lagrange interpolator',                                          &
                    authors     = 'Fortran-FOSS-Programmers',                                       &
                    license     = 'GNU GPLv3',                                                      &
                    description = 'Evaluate coefficients and value(s) of a Lagrange interpolation', &
                    examples    = ["lagrange_interpolation -i            ",                         &
                                   "lagrange_interpolation -i  -c  -int r",                         &
                                   "lagrange_interpolation -i  -p -r     "])
      call cli%add(switch='--interpolate', switch_ab='-i', help='Perform interpolation', required=.false., &
                   def='.false.', act='store_true')
      call cli%add(switch='--interface', switch_ab='-int', help='Interface chosen', required=.false., &
                   def='r', act='store')
      call cli%add(switch='--stencils', switch_ab='-s', help='Stencils used for interpolation', &
                   required=.true., act='store')
      call cli%add(switch='--coeff', switch_ab='-c', help='Save Lagrange coefficients', required=.false., &
                   act='store', def='.false')
      call cli%add(switch='--uniform_grid', switch_ab='-u', help='Activate uniform grid', required=.false., &
                   act='store', def='.false')
    endassociate
    endsubroutine set_cli

    subroutine parse_cli()
    !< Parse Command Line Interface and check its validity.

    call self%cli%parse(error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-i', val=self%inter, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-int', val=self%interf, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-s', val=self%S_number, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-c', val=self%coeff, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-u', val=self%unif_grid, error=self%error) ; if (self%error/=0) stop
    endsubroutine parse_cli
  endsubroutine initialize

  subroutine save_results(self, interpolator)
  !< Save Lagrange interpolation results.

  class(interpolate),           intent(inout) :: self
  class(lagrange_interpolator), intent(in)    :: interpolator
  character(len=:), allocatable       :: buffer       !< Buffer string.
  integer(I_P)                        :: file_unit
  integer(I_P)                        :: i, j         !< Counters.

  if (self%inter) then
    open(newunit=file_unit, file='interpolations.dat')
    write(file_unit, "(A)") 'Interpolation results'
    do i=1,interpolator%S
      write(file_unit, "(A, FR_I, A, FR_P)") 'interpolation(', i, ')= ', interpolator%interp(i)
    enddo
    close(unit=file_unit)
  endif
  if (self%coeff) then
    open(newunit=file_unit, file='coefficients.dat')
    write(file_unit, "(A)") 'Coefficients of Lagrange interpolation'
    do j=1,interpolator%S
      do i=1,interpolator%S
        write(file_unit, "(A, FR_I, FR_I, A, FR_P)") 'coef(', i, j, ')= ', interpolator%coef(i,j)
      enddo
    enddo
    close(unit=file_unit)
  endif
  endsubroutine save_results

  subroutine perform(self, interpolator)
  !< Perform the test.
  class(interpolate),           intent(inout) :: self          !< Lagrange interpolator.
  class(lagrange_interpolator), intent(inout) :: interpolator

  call interpolator%allocate_interpolator(self%S_number)
  select case(self%interf)
  case('r')
    interpolator%x_target = 0.5_R_P
  case('l')
    interpolator%x_target = -0.5_R_P
  endselect
  if (self%unif_grid) then
    call interpolator%init_uniform_grid
  endif
  if (self%inter) then
    call interpolator%compute_interpolations
  endif
  call self%save_results(interpolator)
  endsubroutine perform
endmodule interpolation

program lagrange_interpolation
!< Program for computing Lagrange coefficients and interpolation(s).

use interpolation

implicit none
type(interpolate) :: interp

call interp%execute
endprogram lagrange_interpolation
