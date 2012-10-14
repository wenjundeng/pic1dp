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


! a 64-bit puedo-random number generator with multiple choices of algorithms
!   and multiple choices of distributions
! this module uses some Fortran 2003 features, make sure that your compiler
!   supports Fortran 2003 standard
module multirand
implicit none

! kind number used in this module
! 64-bit integers and real numbers
integer, parameter :: i64 = selected_int_kind(10)
integer, parameter :: r64 = selected_real_kind(10)
! 32-bit integers and real numbers
integer, parameter :: i32 = selected_int_kind(5)
integer, parameter :: r32 = selected_real_kind(5)

! real number of 2**63-1, the maximum value for signed 64-bit integer
real(kind = r64), parameter :: multirand_max64 = 9223372036854775807.0_r64
! real number of 2**64-1, the maximum value for unsigned 64-bit integer
real(kind = r64), parameter :: multirand_maxu64 = 18446744073709551615.0_r64

! real number of 2**31-1, the maximum value for signed 32-bit integer
real(kind = r32), parameter :: multirand_max32 = 2147483647.0_r32
! real number of 2**32-1, the maximum value for unsigned 32-bit integer
real(kind = r32), parameter :: multirand_maxu32 = 4294967295.0_r32
! max # of seeds needed for all algorithms
integer, parameter :: multirand_nseed = 312

! # of algorithms for generating random integers
!integer, parameter :: multirand_nal_int = 2

! # of seeds for each algorithm
!integer, dimension(multirand_nal_int), parameter :: &
!  multirand_nseed = (/ 4, 312 /)

! multirand_int is a wrapper for uniform random integer generation
! choose different generation algorithms using al_int argument 
!   of multirand_init()
procedure(multirand_kiss), pointer :: multirand_int => multirand_kiss

! seeds
integer(kind = i64), dimension(0 : multirand_nseed - 1) :: multirand_seeds

! seed index
integer :: multirand_iseed

! Gaussian random number buffer
real(kind = r64) :: multirand_gbuf
logical :: multirand_gbuf_filled = .false.

contains

!!!!!!!!!!!!!!!!!!!!!!!!
! initialize generator !
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_init(al_int, seed_type, mype, warmup, selftest)
implicit none

! specify algorithm for generating uniform random integers
! 1: George Marsaglia's 64-bit KISS
! 2: 64-bit Mersenne Twister 19937
! (more to be implemented)
integer, intent(in), optional :: al_int

! seed type:
! 1: constant seeds
! 2: random seeds from system_clock (and mype if provided)
! 3: random seeds from /dev/urandom (default)
integer, intent(in), optional :: seed_type

! MPI process rank, used in MPI programs to avoid all processes generating the
!   same random number sequence
integer, intent(in), optional :: mype

! discard (warmup * multirand_nseed) of random numbers as warm-up
! if not specified, warmup = 5
integer, intent(in), optional :: warmup

! specify whether to run self-test
! it is recommended to run self-test at least once before production usage
!   and rerun self-test whenever hardware, operating system, or compiler is
!   is changed
logical, intent(in), optional :: selftest

integer, parameter :: nprime = 100
integer(kind = i64), dimension(0 : nprime - 1), parameter :: &
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
integer, parameter :: furandom = 89

integer :: seed_type_act, warmup_act, flag
integer(kind = i64) :: clock, iseed, nseed
integer(kind = i64), dimension(0 : multirand_nseed - 1) :: tmpseeds
logical :: selftest_act, exist_urandom
character(len = 20) :: msg, msg2

if (present(al_int)) then
  if (al_int == 2) multirand_int => multirand_mt19937
end if
if (associated(multirand_int, multirand_kiss)) then
  nseed = 4
else ! implying multirand_mt19937
  nseed = 312
  multirand_iseed = nseed
end if

! test if the random number generator is working as expected
if (present(selftest)) then
  selftest_act = selftest
else
  selftest_act = .true.
end if
if (selftest_act) then
  if (present(mype)) then
    call multirand_selftest(mype)
  else
    call multirand_selftest()
  end if
end if

if (present(seed_type)) then
  seed_type_act = seed_type
else
  seed_type_act = 3
end if
if (present(mype)) then
  write (msg2, *) mype
  write (msg, *) '[', trim(adjustl(msg2)), ']'
else
  write (msg, *) ''
end if
if (seed_type_act == 3) then
  inquire (file = '/dev/urandom', exist = exist_urandom)
  if (.not. exist_urandom) then
    write (*, '(4a)') trim(adjustl(msg)), '[multirand_init] Warning: ', &
      '/dev/urandom does not exist, ', &
      'falling back to use system_clock for random seeds.'
    seed_type_act = 2
  end if
end if
if (seed_type_act == 3) then
  open (furandom, file = '/dev/urandom', status = 'old', action = 'read', &
    iostat = flag, access = 'stream', form = 'unformatted')
  if (flag /= 0) then
    write (*, '(3a)') trim(adjustl(msg)), '[multirand_init] Warning: ', &
      'error trying to open /dev/urandom, ', &
      'falling back to use system_clock for random seeds.'
    seed_type_act = 2
  end if
end if
if (seed_type_act == 3) then
  read (furandom) multirand_seeds(0 : nseed - 1)
  ! make corrections to satisfy certain seed requirements
  if (associated(multirand_int, multirand_kiss)) then
    do while (multirand_seeds(1) == 0)
      read (furandom) multirand_seeds(1)
    end do
    do while (multirand_seeds(0) == 0 .and. multirand_seeds(3) == 0)
      read (furandom) multirand_seeds(0), multirand_seeds(3)
    end do
  end if
  close (furandom)
  !write (*, *) multirand_seeds
else
  ! first use clock and mype to initialize KISS seeds
  if (seed_type_act == 2) then
    call system_clock(count = clock)
  else ! implying seed_type_act == 1
    clock = primes1(1)
  end if
  multirand_seeds(0 : 3) = clock
  if (present(mype)) then
    multirand_seeds(0 : 3) = multirand_seeds(0 : 3) + primes1( &
      mod(abs(clock + primes2(mod(abs(clock), nprime)) * mype), nprime) &
    ) * mype
  end if
  do iseed = 0, 3
    multirand_seeds(iseed) = multirand_seeds(iseed) + primes2( &
      mod(abs(multirand_seeds(iseed) &
      + primes1(mod(abs(clock), nprime)) * iseed), nprime) &
    ) * iseed
  end do
  ! then use KISS to randomize seeds
  do iseed = 0, nseed - 1
    tmpseeds(iseed) = multirand_kiss()
  end do
  ! make corrections to satisfy certain seed requirements
  if (associated(multirand_int, multirand_kiss)) then
    do while (tmpseeds(1) == 0)
      tmpseeds(1) = multirand_kiss()
    end do
    do while (tmpseeds(0) == 0 .and. tmpseeds(3) == 0)
      tmpseeds(0) = multirand_kiss()
      tmpseeds(3) = multirand_kiss()
    end do
  end if
  multirand_seeds = tmpseeds
end if ! else of if (seed_type_act == 3)

! warm up generator
if (present(warmup)) then
  warmup_act = warmup
else
  warmup_act = 5
end if
do iseed = 1, warmup_act * multirand_nseed
  flag = multirand_int()
end do

if (selftest_act) then
  write (*, '(3a, 5z17)') trim(adjustl(msg)), '[multirand_init] Info: ', &
    'first 5 seeds:', multirand_seeds(0 : 4)
end if

end subroutine multirand_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test if the random number generator is working correctly !
! after calling this you should re-initialized seeds       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_selftest(mype)
implicit none
! MPI process rank
integer, intent(in), optional :: mype

integer, parameter :: ntest = 10
integer(kind = i64), dimension(ntest), target :: test_kiss = (/ &
  8932985056925012148_i64,  5710300428094272059_i64, &
  -104233206776033023_i64, -4143107803135683366_i64, &
   542381058189297533_i64, -4244931820854714191_i64, &
  6853720724624422285_i64,  -767542866500872268_i64, &
  -257204313086867125_i64,  8128797625455304420_i64 /), &
test_mt19937 = (/ &
 -3932459287431434586_i64, 4620546740167642908_i64, &
 -5337173792191653896_i64, -983805426561117294_i64, &
   355488278567739596_i64, 7469126240319926998_i64, &
  4635995468481642529_i64,  418970542659199878_i64, &
 -8842573084457035060_i64, 6358044926049913402_i64 /)

integer(kind =i64), dimension(:), pointer :: ptest
integer(kind = i64) :: i, x
real(kind = r64) :: r
character(len = 20) :: msg, msg2
logical :: match

if (present(mype)) then
  write (msg2, *) mype
  write (msg, *) '[', trim(adjustl(msg2)), ']'
else
  write (msg, *) ''
end if

i = huge(i)
if (i /= 9223372036854775807_i64) then
  write (*, '(3a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    'largest integer is not 2**63-1; multirand is not working correctly.'
end if
r = i / multirand_maxu64 + 0.5_r64
if (r > 1.0_r64) then
  write (*, '(3a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    '[0, 1] real random number upper boundary is larger than 1.0.'
elseif (r < 1.0_r64) then
  write (*, '(3a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    '[0, 1] real random number upper boundary is smaller than 1.0.'
end if
i = i + 1
if (i /= -9223372036854775807_i64 - 1_i64) then
  write (*, '(3a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    'smallest integer is not -2**63; multirand is not working correctly.'
end if
r = i / multirand_maxu64 + 0.5_r64
if (r < 0.0_r64) then
  write (*, *) 'r = ', r
  write (*, '(3a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    '[0, 1] real random number lower boundary is smaller than 0.0.'
elseif (r > 0.0_r64) then
  write (*, '(3a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    '[0, 1] real random number lower boundary is larger than 0.0.'
end if

! test default seeds
if (associated(multirand_int, multirand_kiss)) then
  multirand_seeds(0 : 3) = (/ 1234567890987654321_i64, &
    362436362436362436_i64, 1066149217761810_i64, 123456123456123456_i64/)
  ptest => test_kiss
else ! implying multirand_mt19937
  do i = 1, 312 - 1
    multirand_seeds(i) = 6364136223846793005_i64 &
      * ieor(multirand_seeds(i - 1), ishft(multirand_seeds(i - 1), -62)) &
      + i
  end do
  ptest => test_mt19937
end if
match = .true.
do i = 1, ntest
  if (multirand_int() /= ptest(i)) then
    match = .false.
    exit
  end if
end do
if (.not. match) then
  write (*, '(4a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    'random number generator with default seeds is giving unexpected ', &
    'sequence.'
end if
end subroutine multirand_selftest


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate uniform integer random numbers to fill an array !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_int_array(array)
implicit none

! array to fill in random numbers
integer(kind = i64), dimension(:), intent(out) :: array

integer(kind = i64) :: iarray

do iarray = 1, size(array, kind = i64)
  array(iarray) = multirand_int()
end do
end subroutine multirand_int_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate a uniform [0, 1] real random number !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind = r64) function multirand_real()
implicit none
multirand_real = multirand_int() / multirand_maxu64 + 0.5_r64
! uncomment the following two lines to avoid overflow/underflow
!   due to round-off error
! if (multirand_real > 1.0_r64) multirand_real = 1.0_r64
! if (multirand_real < 0.0_r64) multirand_real = 0.0_r64
! doing this boundary check may sacrifice about 7% efficiency
! use multirand_selftest() to test if overflow/underflow may occur
end function multirand_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate uniform [0, 1] real random numbers to fill an array !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_real_array(array)
implicit none

! array to fill in random numbers
real(kind = r64), dimension(:), intent(out) :: array

integer(kind = i64) :: iarray

do iarray = 1, size(array, kind = i64)
  array(iarray) = multirand_real()
end do
end subroutine multirand_real_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate uniform [0, 1] 32-bit real random numbers to fill an array !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_real_array32(array)
implicit none

! array to fill in random numbers
real(kind = r32), dimension(:), intent(out) :: array

integer(kind = i64) :: iarray, iarray_high, irand
integer(kind = i32) :: irand32_1, irand32_2

iarray_high = size(array, kind = i64)
do iarray = 1, iarray_high, 2
  irand = multirand_int()
  irand32_1 = iand(irand, 4294967295_i64) ! 4294967295 = 0xFFFFFFFF
  array(iarray) = irand32_1 / multirand_maxu32 + 0.5_r32
  if (iarray < iarray_high) then
    irand32_2 = ishft(irand, -32)
    array(iarray + 1) = irand32_2 / multirand_maxu32 + 0.5_r32
  end if
end do
end subroutine multirand_real_array32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! standard Gaussian random number using George Marsaglia's polar method !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind = r64) function multirand_gaussian()
implicit none

real(kind = r64) :: x, y, w
if (multirand_gbuf_filled) then
  multirand_gaussian = multirand_gbuf
  multirand_gbuf_filled = .false.
  return
end if
do
  x = multirand_int() / multirand_max64
  y = multirand_int() / multirand_max64
  w = x * x + y * y
  if (w > 0.0_r64 .and. w < 1.0_r64) exit
end do
w = sqrt((-2.0_r64 * log(w)) / w)
multirand_gaussian = w * x
multirand_gbuf = w * y
multirand_gbuf_filled = .true.
end function multirand_gaussian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill an array with standard Gaussian random numbers !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_gaussian_array(array)
implicit none

! optional array to fill in random numbers
! if omitted, then generate one random number and return
real(kind = r64), dimension(:), intent(out) :: array

integer(kind = i64) :: iarray_low, iarray_high, iarray
real(kind = r64) :: x, y, w

iarray_low = 1
iarray_high = size(array, kind = i64)

if (multirand_gbuf_filled) then
  array(iarray_low) = multirand_gbuf
  iarray_low = iarray_low + 1
  multirand_gbuf_filled = .false.
end if
do iarray = iarray_low, iarray_high, 2
  do
    x = multirand_int() / multirand_max64
    y = multirand_int() / multirand_max64
    w = x * x + y * y
    if (w > 0.0_r64 .and. w < 1.0_r64) exit
  end do
  w = sqrt((-2.0_r64 * log(w)) / w)
  array(iarray) = x * w
  if (iarray < iarray_high) then
    array(iarray + 1) = y * w
  else
    multirand_gbuf = y * w
    multirand_gbuf_filled = .true.
  end if
end do

end subroutine multirand_gaussian_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill an 32-bit array with standard Gaussian random numbers !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_gaussian_array32(array)
implicit none

! optional array to fill in random numbers
! if omitted, then generate one random number and return
real(kind = r32), dimension(:), intent(out) :: array

integer(kind = i64) :: iarray_low, iarray_high, iarray, i
real(kind = r32) :: x, y, w

iarray_low = 1
iarray_high = size(array, kind = i64)

if (multirand_gbuf_filled) then
  array(iarray_low) = multirand_gbuf
  iarray_low = iarray_low + 1
  multirand_gbuf_filled = .false.
end if
do iarray = iarray_low, iarray_high, 2
  do
    i = multirand_int()
    x = iand(i, 4294967295_i64) / multirand_max32
    y = ishft(i, -32) / multirand_max32
    w = x * x + y * y
    if (w > 0.0_r32 .and. w < 1.0_r32) exit
  end do
  w = sqrt((-2.0_r32 * log(w)) / w)
  array(iarray) = x * w
  if (iarray < iarray_high) then
    array(iarray + 1) = y * w
  end if
end do

end subroutine multirand_gaussian_array32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! George Marsaglia's 64-bit KISS                                             !
! https://groups.google.com/d/msg/comp.lang.fortran/qFv18ql_WlU/IK8KGZZFJx4J !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind = i64) function multirand_kiss()
implicit none

integer(kind = i64) :: t

t = ishft(multirand_seeds(0), 58) + multirand_seeds(3)
if (s(multirand_seeds(0)) == s(t)) then
  multirand_seeds(3) = ishft(multirand_seeds(0), -6) &
    + s(multirand_seeds(0))
else
  multirand_seeds(3) = ishft(multirand_seeds(0), -6) &
    - s(multirand_seeds(0) + t) + 1
end if
multirand_seeds(0) = multirand_seeds(0) + t
multirand_seeds(1) = m(m(m(multirand_seeds(1), 13_i64), -17_i64), 43_i64)
multirand_seeds(2) = 6906969069_i64 * multirand_seeds(2) + 1234567_i64
multirand_kiss = multirand_seeds(0) + multirand_seeds(1) + multirand_seeds(2)

contains

  integer(kind = i64) function m(x, k)
  implicit none
  integer(kind = i64), intent(in) :: x, k
  m = ieor(x, ishft(x, k))
  end function m

  integer(kind = i64) function s(x)
  implicit none
  integer(kind = i64), intent(in) :: x
  s = ishft(x, -63)
  end function s

end function multirand_kiss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 64-bit Mersenne Twister 19937                              !
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind = i64) function multirand_mt19937()
implicit none

integer, parameter :: &
  NN = 312, &
  MM = 156
integer(kind = i64), parameter :: &
  UM = -2147483648_i64, & ! 0xFFFFFFFF80000000ULL
  LM = 2147483647_i64 ! 0x7FFFFFFFULL
integer(kind = i64), dimension(0 : 1), parameter :: mag01 = &
  (/ 0_i64, -5403634167711393303_i64 & ! MATRIX_A 0xB5026F5AA96619E9ULL
  /)

integer :: i
integer(kind = i64) :: x

if (multirand_iseed >= NN) then
  do i = 0, NN - MM - 1
    x = ior(iand(multirand_seeds(i), UM), iand(multirand_seeds(i + 1), LM))
    multirand_seeds(i) = ieor( &
      ieor(multirand_seeds(i + MM), ishft(x, -1)), &
      mag01(iand(x, 1_i64)))
  end do
  do i = NN - MM, NN - 2
    x = ior(iand(multirand_seeds(i), UM), iand(multirand_seeds(i + 1), LM))
    multirand_seeds(i) = ieor( &
      ieor(multirand_seeds(i + (MM - NN)), ishft(x, -1)), &
      mag01(iand(x, 1_i64)))
  end do
  x = ior(iand(multirand_seeds(NN - 1), UM), iand(multirand_seeds(0), LM))
  multirand_seeds(NN - 1) = ieor( &
    ieor(multirand_seeds(MM - 1), ishft(x, -1)), &
    mag01(iand(x, 1_i64)))
  multirand_iseed = 0
end if

x = multirand_seeds(multirand_iseed)
multirand_iseed = multirand_iseed + 1
x = ieor(x, iand(ishft(x, -29), 6148914691236517205_i64))
x = ieor(x, iand(ishft(x, 17), 8202884508482404352_i64))
x = ieor(x, iand(ishft(x, 37), -2270628950310912_i64))
! the 3 long integers in the above 3 lines are:
! 0x5555555555555555ULL, 0x71D67FFFEDA60000ULL, 0xFFF7EEE000000000ULL
x = ieor(x, ishft(x, -43))

multirand_mt19937 = x
end function multirand_mt19937

end module multirand

