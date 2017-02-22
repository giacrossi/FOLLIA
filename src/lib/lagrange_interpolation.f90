recursive subroutine symplify(num, den, symplified)
use PENF, only : I_P
implicit none

integer(I_P), intent(inout) :: num, den
logical,      intent(inout) :: symplified

if (.not.(symplified)) then
  if (mod(num,8_I_P)/=0_I_P) then
    if (mod(num,4_I_P)/=0_I_P) then
      if (mod(num,2_I_P)/=0_I_P) then
        symplified=.true.
      else
        if ((abs(num)>=2_I_P).and.(abs(den)>=2_I_P)) then
          num = num / 2_I_P
          den = den / 2_I_P
          symplified=.true.
        endif
      endif
    else
      if (mod(den,4_I_P)==0_I_P.and.(abs(num)>=4_I_P).and.(abs(den)>=4_I_P)) then
        num = num / 4_I_P
        den = den / 4_I_P
      elseif (mod(den,2_I_P)==0_I_P.and.(abs(num)>=2_I_P).and.(abs(den)>=2_I_P)) then
        num = num / 2_I_P
        den = den / 2_I_P
        symplified=.true.
      else
        symplified=.true.
      endif
    endif
  else
    if (mod(den,8_I_P)==0_I_P.and.(abs(num)>=8_I_P).and.(abs(den)>=8_I_P)) then
      num = num / 8_I_P
      den = den / 8_I_P
    elseif (mod(den,4_I_P)/=0_I_P) then
      if (mod(num,2_I_P)==0_I_P.and.(abs(num)>=2_I_P).and.(abs(den)>=2_I_P)) then
        num = num / 2_I_P
        den = den / 2_I_P
        symplified=.true.
      else
        symplified=.true.
      endif
    else
      if ((abs(num)>=4_I_P).and.(abs(den)>=4_I_P)) then
        num = num / 4_I_P
        den = den / 4_I_P
      endif
    endif
  endif
  call symplify(num, den, symplified)
endif
endsubroutine symplify

program lagrange_interpolation
!< Program for computing Lagrange coefficients and interpolation(s).

use follia
#ifdef r16p
use penf, only: I_P, str, RPP=>R16P
#else
use penf, only: I_P, str, RPP=>R8P
#endif

implicit none

interface
  subroutine symplify(num, den, symplified)
    use PENF, only : I_P
    integer(I_P), intent(inout) :: num, den
    logical,      intent(inout) :: symplified
  endsubroutine
endinterface

type(lagrange_interpolator) :: interp
integer(I_P), allocatable   :: den(:,:)
integer(I_P)                :: i, j, S
integer(I_P)                :: file_unit
logical                     :: symp
integer(I_P), allocatable   :: int_coef(:,:)

print *, 'Insert stencil dimension and number'
read *, S

call interp%initialize_interpolator(S,-0.5_RPP)

call interp%compute_coefficients
call interp%compute_interpolations

interp%coef(:,:) = interp%coef(:,:) * (2._RPP**20_RPP)

allocate(int_coef(1:interp%S,1:interp%S))
allocate(den(1:interp%S,1:interp%S))

int_coef=nint(interp%coef)

den(:,:) = 2_I_P**20_I_P

symp=.false.

associate(coef=>int_coef)
  do j=1,interp%S
    do i=1,interp%S
      call symplify(coef(i,j),den(i,j),symp)
      symp=.false.
    enddo
  enddo
endassociate

open(newunit=file_unit, file='coefficients_S-'//trim(str(interp%S))//'.dat')
write(file_unit,"(A)") 'LEFT INTERFACE'

do j=1,interp%S
  do i=1,interp%S
    write(file_unit,"(A,I1,A,I1,A,I7,A,I7,A)") 'c(1,', i-1_I_P,',', j-1_I_P,')= ', int_coef(i,j),'._RPP /', den(i,j),'._RPP'
  enddo
enddo

call interp%initialize_interpolator(S,0.5_RPP)

call interp%compute_coefficients
call interp%compute_interpolations

interp%coef(:,:) = interp%coef(:,:) * (2._RPP**20_RPP)

int_coef=nint(interp%coef)

den(:,:) = 2_I_P**20_I_P

symp=.false.

associate(coef=>int_coef)
  do j=1,interp%S
    do i=1,interp%S
      call symplify(coef(i,j),den(i,j),symp)
      symp=.false.
    enddo
  enddo
endassociate

write(file_unit,"(A)") 'RIGHT INTERFACE'

do j=1,interp%S
  do i=1,interp%S
    write(file_unit,"(A,I1,A,I1,A,I7,A,I7,A)") 'c(2,', i-1_I_P,',', j-1_I_P,')= ', int_coef(i,j),'._RPP /', den(i,j),'._RPP'
  enddo
enddo

close(unit=file_unit)
endprogram lagrange_interpolation
