! Copyright 2012 Wenjun Deng <wdeng@wdeng.info>
!
! This file is part of PIC1D-PETSc
!
! PIC1D-PETSc is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PIC1D-PETSc is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with PIC1D-PETSc.  If not, see <http://www.gnu.org/licenses/>.


! a Gaussian random number generator using polar form of Box-Muller transform

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
integer, parameter :: gik = selected_int_kind(5)
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

integer, parameter :: nprime = 100
integer, dimension(0 : nprime - 1), parameter :: &
  primes1 = (/ &
    15484219, 15484223, 15484243, 15484247, 15484279, &
    15484333, 15484363, 15484387, 15484393, 15484409, &
    15484421, 15484453, 15484457, 15484459, 15484471, &
    15484489, 15484517, 15484519, 15484549, 15484559, &
    15484591, 15484627, 15484631, 15484643, 15484661, &
    15484697, 15484709, 15484723, 15484769, 15484771, &
    15484783, 15484817, 15484823, 15484873, 15484877, &
    15484879, 15484901, 15484919, 15484939, 15484951, &
    15484961, 15484999, 15485039, 15485053, 15485059, &
    15485077, 15485083, 15485143, 15485161, 15485179, &
    15485191, 15485221, 15485243, 15485251, 15485257, &
    15485273, 15485287, 15485291, 15485293, 15485299, &
    15485311, 15485321, 15485339, 15485341, 15485357, &
    15485363, 15485383, 15485389, 15485401, 15485411, &
    15485429, 15485441, 15485447, 15485471, 15485473, &
    15485497, 15485537, 15485539, 15485543, 15485549, &
    15485557, 15485567, 15485581, 15485609, 15485611, &
    15485621, 15485651, 15485653, 15485669, 15485677, &
    15485689, 15485711, 15485737, 15485747, 15485761, &
    15485773, 15485783, 15485801, 15485807, 15485837 /), &
  primes2 = (/ &
    7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
    7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
    7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
    7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
    7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
    7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
    7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
    7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
    7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
    7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

integer :: clock, iseed, nseed
integer, dimension(:), allocatable :: seeds

real(kind = grk), dimension(:), allocatable :: rseeds

! initialize random number seeds
call random_seed(size = nseed)
allocate (seeds(nseed), rseeds(nseed))
if (seed_type == 2) then
  clock = primes1(1)
else
  call system_clock(count = clock)
end if
seeds(:) = clock
if (present(mype)) then
  seeds(:) = seeds(:) + primes1( &
    mod(clock + primes2(mod(clock, nprime)) * mype, nprime) &
  ) * mype
end if
do iseed = 1, nseed
  seeds(iseed) = seeds(iseed) + primes2( &
    mod(seeds(iseed) + primes1(mod(clock, nprime)) * iseed, nprime) &
  ) * iseed
end do
call random_seed(put = seeds)

! use random numbers to make seeds more random
call random_number(rseeds)
seeds(:) = (rseeds(:) - 0.5_grk) * 2.0_grk * huge(0)
call random_seed(put = seeds)

deallocate (seeds, rseeds)

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

