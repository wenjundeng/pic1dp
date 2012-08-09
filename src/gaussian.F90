! generating Gaussian distribution random numbers

! if this module is used in a PETSc program,
! use -D__PETSc option when compiling, or uncomment the following line
! #define __PETSc

module gaussian
implicit none

! real kind # used in this module
#ifdef __PETSc
#include "finclude/petscdef.h"
PetscReal, parameter :: gaussian_testkind = 0d0
PetscInt, parameter :: grk = kind(gaussian_testkind)
PetscInt, parameter :: gik = kind(grk)
#else
integer, parameter :: grk = selected_real_kind(10)
integer, parameter :: gik = selected_int_kind(10)
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize random seeds !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gaussian_init(seed_type, mype)
implicit none
! seed type: 1: random seeds (using system_clock); 2: constant seeds
integer(kind = gik), intent(in) :: seed_type

! MPI process rank, used in MPI programs to avoid all processes generating the
! same random number sequence
integer(kind = gik), intent(in), optional :: mype

integer(kind = gik) :: clock, iseed, nseed
integer(kind = gik), dimension(:), allocatable :: seeds

! initialize random number seeds
call random_seed(size = nseed)
allocate (seeds(nseed))
if (seed_type == 2) then
  clock = 1111
else
  call system_clock(count = clock)
end if
if (present(mype)) then
  seeds = clock + 991 * mype + 997 * (/ (iseed - 1, iseed = 1, nseed) /)
else
  seeds = clock + 997 * (/ (iseed - 1, iseed = 1, nseed) /)
end if
call random_seed(put = seeds)
deallocate (seeds)

end subroutine gaussian_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate Gaussian random numbers and put in array !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gaussian_generate(array)
implicit none

real(kind = grk), dimension(:), intent(out) :: array

integer(kind = gik) :: iarray_low, iarray_high, iarray
real(kind = grk), dimension(2) :: unirand1, unirand2 ! uniform random numbers
real(kind = grk) :: w

iarray_low = lbound(array, 1)
iarray_high = ubound(array, 1)

do iarray = iarray_low, iarray_high, 2
  do
    call random_number(unirand1)
    unirand2 = 2.0_grk * unirand1 - 1.0_grk
    w = unirand2(1) * unirand2(1) + unirand2(2) * unirand2(2)
    if (w < 1.0_grk) exit
  end do
  w = sqrt((-2.0_grk * log(w)) / w)
  array(iarray) = unirand2(1) * w
  if (iarray < iarray_high) array(iarray + 1) = unirand2(2) * w
end do

end subroutine gaussian_generate

end module gaussian

