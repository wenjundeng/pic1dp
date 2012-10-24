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


! multirand is a 64-bit puedo-random number generator with multiple choices
!   of algorithms and multiple choices of distributions
! (interfaces for 32-bit random numbers are also provided)
! (32-bit interfaces are not fully tested yet)
! multirand uses some Fortran 2003 features, make sure that your compiler
!   supports Fortran 2003 standard (GNU Fortran 4.6 suffices)
! if your compiler does not support procedure pointer well, then
!   use -DNO_PROC_POINTER option when compiling or uncomment the
!   following line:
! #define NO_PROC_POINTER
!
! usage of multirand:
!   use multirand_init(...) to initialize multirand
!   then you can use these functions to generate a single random number:
!     multirand_int64(), multirand_int32(),
!     multirand_real64(), multirand_real32(),
!     multirand_gaussian64(), multirand_gaussian32()
!   you can use these subroutines to fill an array with random numbers:
!     multirand_int_array64(array), multirand_int_array32(array),
!     multirand_real_array64(array), multirand_real_array32(array),
!     multirand_gaussian_array64(array), multirand_gaussian_array32(array)
!   you can also use these generic interfaces to fill an array:
!     multirand_int_array(array),
!     multirand_real_array(array),
!     multirand_gaussian_array(array)
! note that these functions and subroutines are NOT thread-safe,
!   which needs to be addressed in future
module multirand
implicit none

! macros for random integer to [0, 1] real conversion
#define INT2REAL64(i) i / multirand_maxu64 + 0.5_mrkr64
#define INT2REAL32(i) i / multirand_maxu32 + 0.5_mrkr32

! macros for extracting two 32-bit integers from a 64-bit integer
! 4294967295 = 0xFFFFFFFF
#define INT64TO32_1(i) iand(i, 4294967295_mrki64)
#define INT64TO32_2(i) ishft(i, -32)

! macro for xor-shift operation
#define XORSHFT(x, shft) ieor(x, ishft(x, shft))
! macro for xor-and-shift operation
#define XORANDSHFT(x, shft, and) ieor(x, iand(ishft(x, shft), and))

! kind numbers used in this module
! 64-bit integers and real numbers
integer, parameter :: mrki64 = selected_int_kind(10)
integer, parameter :: mrkr64 = selected_real_kind(10)
! 32-bit integers and real numbers
integer, parameter :: mrki32 = selected_int_kind(5)
integer, parameter :: mrkr32 = selected_real_kind(5)

! real number of 2**63-1, the maximum value for signed 64-bit integer
real(kind = mrkr64), parameter :: multirand_max64 &
  = 9223372036854775807.0_mrkr64
! real number of 2**64-1, the maximum value for unsigned 64-bit integer
real(kind = mrkr64), parameter :: multirand_maxu64 &
  = 18446744073709551615.0_mrkr64

! real number of 2**31-1, the maximum value for signed 32-bit integer
real(kind = mrkr32), parameter :: multirand_max32 = 2147483647.0_mrkr32
! real number of 2**32-1, the maximum value for unsigned 32-bit integer
real(kind = mrkr32), parameter :: multirand_maxu32 = 4294967295.0_mrkr32

! max # of seeds needed for all algorithms
integer, parameter :: multirand_nseed = 20635

! seeds
integer(kind = mrki64), dimension(0 : multirand_nseed - 1) :: multirand_seeds

! seed index
integer :: multirand_iseed

#ifdef NO_PROC_POINTER
! algorithm index
integer :: multirand_al_int
#else
! multirand_int64 is a wrapper for 64-bit uniform random integer generation
! choose different generation algorithms using al_int argument 
!   of multirand_init()
procedure(multirand_kiss64), pointer :: multirand_int64 => null()
#endif

! 32-bit random integer buffer
integer(kind = mrki32) :: multirand_int32buf
logical :: multirand_int32buf_filled = .false.

! Gaussian random number buffers
real(kind = mrkr64) :: multirand_gaussian64buf
logical :: multirand_gaussian64buf_filled = .false.
real(kind = mrkr32) :: multirand_gaussian32buf
logical :: multirand_gaussian32buf_filled = .false.

! generic interfaces for filling an array with random numbers
interface multirand_int_array
  module procedure multirand_int_array64
  module procedure multirand_int_array32
end interface

interface multirand_real_array
  module procedure multirand_real_array64
  module procedure multirand_real_array32
end interface

interface multirand_gaussian_array
  module procedure multirand_gaussian_array64
  module procedure multirand_gaussian_array32
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!
! initialize generator !
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_init(al_int, seed_type, mype, warmup, selftest)
implicit none

! specify algorithm for generating uniform random integers
! 1: George Marsaglia's 64-bit KISS (period ~ 2**(247.42) ~ 10**(74.48))
! 2: 64-bit Mersenne Twister 19937 (period = 2**19937 - 1 ~ 10**6001)
! 3: George Marsaglia's 64-bit SuperKISS
!   (period = 5*2**1320480*(2**64-1) ~ 10**397524)
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
integer(kind = mrki64), dimension(0 : nprime - 1), parameter :: &
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

integer :: al_int_act, seed_type_act, warmup_act, flag
integer(kind = mrki64) :: clock, iseed, nseed
integer(kind = mrki64), dimension(0 : multirand_nseed - 1) :: tmpseeds
logical :: selftest_act, exist_urandom
character(len = 20) :: msg, msg2

! choose algorithm for 64-bit uniform random integer
if (present(al_int)) then
  al_int_act = al_int
else
  al_int_act = 3
end if
#ifdef NO_PROC_POINTER
multirand_al_int = al_int_act
#endif
if (al_int_act == 2) then
#ifndef NO_PROC_POINTER
  multirand_int64 => multirand_mt19937_64
#endif
  nseed = 312
elseif (al_int_act == 3) then
#ifndef NO_PROC_POINTER
  multirand_int64 => multirand_superkiss64
#endif
  nseed = 20635
else
#ifndef NO_PROC_POINTER
  multirand_int64 => multirand_kiss64
#endif
  nseed = 4
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

! initialize seeds
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
#ifdef NO_PROC_POINTER
  if (multirand_al_int == 1) then
#else
  if (associated(multirand_int64, multirand_kiss64)) then
#endif
    do while (multirand_seeds(1) == 0)
      read (furandom) multirand_seeds(1)
    end do
    do while (multirand_seeds(0) == 0 .and. multirand_seeds(3) == 0)
      read (furandom) multirand_seeds(0), multirand_seeds(3)
    end do
#ifdef NO_PROC_POINTER
  elseif (multirand_al_int == 3) then
#else
  elseif (associated(multirand_int64, multirand_superkiss64)) then
#endif
    do while (multirand_seeds(20634) == 0)
      read (furandom) multirand_seeds(20634)
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
      mod(abs(clock + primes2(mod(abs(clock), int(nprime, mrki64))) * mype), &
      int(nprime, mrki64)) &
    ) * mype
  end if
  do iseed = 0, 3
    multirand_seeds(iseed) = multirand_seeds(iseed) + primes2( &
      mod(abs(multirand_seeds(iseed) &
      + primes1(mod(abs(clock), int(nprime, mrki64))) * iseed), &
      int(nprime, mrki64)) &
    ) * iseed
  end do
  ! then use KISS to randomize seeds
  do iseed = 1, 20 ! warm up KISS generator
    tmpseeds(0) = multirand_kiss64()
  end do
  do iseed = 1, nseed - 1
    tmpseeds(iseed) = multirand_kiss64()
  end do
  ! make corrections to satisfy certain seed requirements
#ifdef NO_PROC_POINTER
  if (multirand_al_int == 1) then
#else
  if (associated(multirand_int64, multirand_kiss64)) then
#endif
    do while (tmpseeds(1) == 0)
      tmpseeds(1) = multirand_kiss64()
    end do
    do while (tmpseeds(0) == 0 .and. tmpseeds(3) == 0)
      tmpseeds(0) = multirand_kiss64()
      tmpseeds(3) = multirand_kiss64()
    end do
#ifdef NO_PROC_POINTER
  elseif (multirand_al_int == 3) then
#else
  elseif (associated(multirand_int64, multirand_superkiss64)) then
#endif
    do while (multirand_seeds(20634) == 0)
      tmpseeds(20634) = multirand_kiss64()
    end do
  end if
  multirand_seeds(:) = tmpseeds(:)
end if ! else of if (seed_type_act == 3)

! set seed array index
#ifdef NO_PROC_POINTER
if (multirand_al_int == 2) then
#else
if (associated(multirand_int64, multirand_mt19937_64)) then
#endif
  multirand_iseed = 312
#ifdef NO_PROC_POINTER
elseif (multirand_al_int == 3) then
#else
elseif (associated(multirand_int64, multirand_superkiss64)) then
#endif
  multirand_iseed = 20632
end if

!if (selftest_act) then
!  write (*, '(3a, 5z17)') trim(adjustl(msg)), '[multirand_init] Info: ', &
!    'first 5 seeds:', multirand_seeds(0 : 4)
!end if

! warm up the generator
if (present(warmup)) then
  warmup_act = warmup
else
  warmup_act = 5
end if
do iseed = 1, warmup_act * nseed
  flag = multirand_int64()
end do

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
integer(kind = mrki64), dimension(ntest), target :: test_kiss64 = (/ &
  8932985056925012148_mrki64,  5710300428094272059_mrki64, &
  -104233206776033023_mrki64, -4143107803135683366_mrki64, &
   542381058189297533_mrki64, -4244931820854714191_mrki64, &
  6853720724624422285_mrki64,  -767542866500872268_mrki64, &
  -257204313086867125_mrki64,  8128797625455304420_mrki64 /), &
test_mt19937_64_head = (/ &
 -3932459287431434586_mrki64, 4620546740167642908_mrki64, &
 -5337173792191653896_mrki64, -983805426561117294_mrki64, &
   355488278567739596_mrki64, 7469126240319926998_mrki64, &
  4635995468481642529_mrki64,  418970542659199878_mrki64, &
 -8842573084457035060_mrki64, 6358044926049913402_mrki64 /), &
test_mt19937_64_tail = (/ &
-7948593974297132281_mrki64,  1921007855220546564_mrki64, &
 7643484074408755248_mrki64, -7128315020423208677_mrki64, &
 1370093900783164344_mrki64,  6776537281339823025_mrki64, &
 3450492372588984223_mrki64, -9045729527952115285_mrki64, &
 7896519943553875907_mrki64, -4143300141377237606_mrki64 /), &
test_superkiss64_head = (/ &
  6140839658375754198_mrki64,   -95225469143006167_mrki64, &
 -9148462456964506707_mrki64,  3912874252778582253_mrki64, &
  6801212277726928591_mrki64,  -809575511391043410_mrki64, &
  -397286769868273005_mrki64,  4963780769400405858_mrki64, &
  2406624640673457322_mrki64,  1246843699883922102_mrki64 /), &
test_superkiss64_tail = (/ &
 -1387224431860786161_mrki64, -8846516422183390713_mrki64, &
  8111357788999165247_mrki64,   444070776306226770_mrki64, &
 -7730678117654887867_mrki64,  -296399128303442035_mrki64, &
 -1658509282659454084_mrki64, -8190332265239255687_mrki64, &
 -1492517620356299342_mrki64, -5016179395587873849_mrki64 /)

integer(kind = mrki64), dimension(:), pointer :: ptest_head, ptest_tail
integer :: i, itail
integer(kind = mrki64) :: k
real(kind = mrkr64) :: r
character(len = 20) :: msg, msg2
logical :: match

if (present(mype)) then
  write (msg2, *) mype
  write (msg, *) '[', trim(adjustl(msg2)), ']'
else
  write (msg, *) ''
end if

! test boundaries
k = huge(k)
if (k /= 9223372036854775807_mrki64) then
  write (msg2, *) k
  write (*, '(5a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    'largest integer is ', trim(adjustl(msg2)), &
    ' instead of 2**63-1=9223372036854775807; ', &
    'multirand is not working correctly.'
end if
r = k / multirand_maxu64 + 0.5_mrkr64
if (r /= 1.0_mrkr64) then
  write (msg2, *) r
  write (*, '(5a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    '[0, 1] real random number upper boundary is ', trim(adjustl(msg2)), &
    ' instead of 1.0.'
end if
k = k + 1
if (k /= -9223372036854775807_mrki64 - 1_mrki64) then
  write (msg2, *) k
  write (*, '(5a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    'smallest integer is ', trim(adjustl(msg2)), &
    ' instead of -2**63=-9223372036854775808; ', &
    'multirand is not working correctly.'
end if
r = k / multirand_maxu64 + 0.5_mrkr64
if (r /= 0.0_mrkr64) then
  write (msg2, *) r
  write (*, '(5a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    '[0, 1] real random number lower boundary is ', trim(adjustl(msg2)), &
    ' instead of 0.0.'
end if

! pending to add test of 32-bit numbers

! test default seeds
#ifdef NO_PROC_POINTER
if (multirand_al_int == 2) then
#else
if (associated(multirand_int64, multirand_mt19937_64)) then
#endif
  multirand_seeds(0) = 5489_mrki64
  do multirand_iseed = 1, 312 - 1
    multirand_seeds(multirand_iseed) = 6364136223846793005_mrki64 &
      * XORSHFT(multirand_seeds(multirand_iseed - 1), -62) &
      + multirand_iseed
  end do
  multirand_iseed = 312
  ptest_head => test_mt19937_64_head
  ptest_tail => test_mt19937_64_tail
  itail = 312 - ntest / 2
#ifdef NO_PROC_POINTER
elseif (multirand_al_int == 3) then
#else
elseif (associated(multirand_int64, multirand_superkiss64)) then
#endif
  multirand_seeds(20632 : 20634) = (/ 36243678541_mrki64, &
    12367890123456_mrki64, 521288629546311_mrki64 /)
  do multirand_iseed = 0, 20631
    multirand_seeds(20633) = multirand_seeds(20633) * 6906969069_mrki64 + 123
    multirand_seeds(20634) = XORSHFT(multirand_seeds(20634), 13)
    multirand_seeds(20634) = XORSHFT(multirand_seeds(20634), -17)
    multirand_seeds(20634) = XORSHFT(multirand_seeds(20634), 43)
    multirand_seeds(multirand_iseed) &
      = multirand_seeds(20633) + multirand_seeds(20634)
  end do
  !write (*, *) multirand_seeds(20632 : 20634), &
  !  multirand_seeds(0), multirand_seeds(20631)
  multirand_iseed = 20632
  ptest_head => test_superkiss64_head
  ptest_tail => test_superkiss64_tail
  itail = 20632 - ntest / 2
else ! implying multirand_kiss64
  multirand_seeds(0 : 3) = (/ 1234567890987654321_mrki64, &
    362436362436362436_mrki64, 1066149217761810_mrki64, &
    123456123456123456_mrki64/)
  ptest_head => test_kiss64
  nullify (ptest_tail)
end if
match = .true.
do i = 1, ntest
!  k = multirand_int64()
!  if (k /= ptest_head(i)) then
!    write (*, *) k, ptest_head(i)
  if (multirand_int64() /= ptest_head(i)) then
    match = .false.
    exit
  end if
end do
if (.not. match) then
  write (*, '(4a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
    'random number generator with default seeds is giving unexpected ', &
    'head sequence.'
elseif (associated(ptest_tail)) then
  ! test tail sequence
  do i = ntest + 1, itail
    k = multirand_int64()
  end do
  do i = itail + 1, itail + ntest
!    k = multirand_int64()
!    if (k /= ptest_tail(i)) then
!      write (*, *) i, k, ptest_tail(i)
    if (multirand_int64() /= ptest_tail(i - itail)) then
      match = .false.
      exit
    end if
  end do
  if (.not. match) then
    write (*, '(4a)') trim(adjustl(msg)), '[multirand_selftest] Warning: ', &
      'random number generator with default seeds is giving unexpected ', &
      'tail sequence.'
  end if
end if ! elseif (associated(ptest_tail))
end subroutine multirand_selftest


#ifdef NO_PROC_POINTER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate a 64-bit uniform random integer !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind = mrki64) function multirand_int64()
implicit none
if (multirand_al_int == 2) then
  multirand_int64 = multirand_mt19937_64()
elseif (multirand_al_int == 3) then
  multirand_int64 = multirand_superkiss64()
else
  multirand_int64 = multirand_kiss64()
end if
end function multirand_int64
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate a 32-bit uniform random integer !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind = mrki32) function multirand_int32()
implicit none
integer(kind = mrki64) :: irand64

if (multirand_int32buf_filled) then
  multirand_int32 = multirand_int32buf
  multirand_int32buf_filled = .false.
  return
end if
irand64 = multirand_int64()
multirand_int32 = INT64TO32_1(irand64)
multirand_int32buf = INT64TO32_2(irand64)
multirand_int32buf_filled = .true.
end function multirand_int32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill array with 64-bit uniform random integers !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_int_array64(array)
implicit none

! array to fill in random numbers
integer(kind = mrki64), dimension(:), intent(out) :: array

integer(kind = mrki64) :: iarray

do iarray = 1, size(array, kind = mrki64)
  array(iarray) = multirand_int64()
end do
end subroutine multirand_int_array64


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill array with 32-bit uniform random integers !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_int_array32(array)
implicit none

! array to fill in random numbers
integer(kind = mrkr32), dimension(:), intent(out) :: array

integer(kind = mrki64) :: iarray, iarray_low, iarray_high, irand64

iarray_low = 1
iarray_high = size(array, kind = mrki64)
if (multirand_int32buf_filled) then
  array(iarray_low) = multirand_int32buf
  multirand_int32buf_filled = .false.
  iarray_low = iarray_low + 1
end if
do iarray = iarray_low, iarray_high, 2
  irand64 = multirand_int64()
  array(iarray) = INT64TO32_1(irand64)
  if (iarray < iarray_high) then
    array(iarray + 1) = INT64TO32_2(irand64)
  else
    multirand_int32buf = INT64TO32_2(irand64)
    multirand_int32buf_filled = .true.
  end if
end do
end subroutine multirand_int_array32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate a 64-bit uniform [0, 1] real random number !
! use multirand_selftest() for safe boundary check    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind = mrkr64) function multirand_real64()
implicit none
multirand_real64 = INT2REAL64(multirand_int64())
end function multirand_real64


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate a 32-bit uniform [0, 1] real random number !
! use multirand_selftest() for safe boundary check    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind = mrkr32) function multirand_real32()
implicit none
multirand_real32 = INT2REAL32(multirand_int32())
end function multirand_real32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill array with 64-bit uniform [0, 1] real random numbers !
! use multirand_selftest() for safe boundary check          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_real_array64(array)
implicit none

! array to fill in random numbers
real(kind = mrkr64), dimension(:), intent(out) :: array

integer(kind = mrki64) :: iarray

do iarray = 1, size(array, kind = mrki64)
  array(iarray) = INT2REAL64(multirand_int64())
end do
end subroutine multirand_real_array64


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill array with 32-bit uniform [0, 1] real random numbers !
! use multirand_selftest() for safe boundary check          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_real_array32(array)
implicit none

! array to fill in random numbers
real(kind = mrkr32), dimension(:), intent(out) :: array

integer(kind = mrki64) :: iarray, iarray_low, iarray_high, irand64

iarray_low = 1
iarray_high = size(array, kind = mrki64)
if (multirand_int32buf_filled) then
  array(iarray_low) = INT2REAL32(multirand_int32buf)
  multirand_int32buf_filled = .false.
  iarray_low = iarray_low + 1
end if
do iarray = iarray_low, iarray_high, 2
  irand64 = multirand_int64()
  array(iarray) = INT2REAL32(INT64TO32_1(irand64))
  if (iarray < iarray_high) then
    array(iarray + 1) = INT2REAL32(INT64TO32_2(irand64))
  else
    multirand_int32buf = INT64TO32_2(irand64)
    multirand_int32buf_filled = .true.
  end if
end do
end subroutine multirand_real_array32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate a 64-bit standard Gaussian random number !
! using George Marsaglia's polar method             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind = mrkr64) function multirand_gaussian64()
implicit none

real(kind = mrkr64) :: x, y, w
if (multirand_gaussian64buf_filled) then
  multirand_gaussian64 = multirand_gaussian64buf
  multirand_gaussian64buf_filled = .false.
  return
end if
do
  x = multirand_int64() / multirand_max64
  y = multirand_int64() / multirand_max64
  w = x * x + y * y
  if (w > 0.0_mrkr64 .and. w < 1.0_mrkr64) exit
end do
w = sqrt((-2.0_mrkr64 * log(w)) / w)
multirand_gaussian64 = w * x
multirand_gaussian64buf = w * y
multirand_gaussian64buf_filled = .true.
end function multirand_gaussian64


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate a 32-bit standard Gaussian random number !
! using George Marsaglia's polar method             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind = mrkr32) function multirand_gaussian32()
implicit none

integer(kind = mrki64) :: irand64
real(kind = mrkr32) :: x, y, w
if (multirand_gaussian32buf_filled) then
  multirand_gaussian32 = multirand_gaussian32buf
  multirand_gaussian32buf_filled = .false.
  return
end if
do
  irand64 = multirand_int64()
  x = INT64TO32_1(irand64) / multirand_max32
  y = INT64TO32_2(irand64) / multirand_max32
  w = x * x + y * y
  if (w > 0.0_mrkr32 .and. w < 1.0_mrkr32) exit
end do
w = sqrt((-2.0_mrkr32 * log(w)) / w)
multirand_gaussian32 = w * x
multirand_gaussian32buf = w * y
multirand_gaussian32buf_filled = .true.
end function multirand_gaussian32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill array with 64-bit standard Gaussian random numbers !
! using George Marsaglia's polar method                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_gaussian_array64(array)
implicit none

! array to fill in random numbers
real(kind = mrkr64), dimension(:), intent(out) :: array

integer(kind = mrki64) :: iarray_low, iarray_high, iarray
real(kind = mrkr64) :: x, y, w

iarray_low = 1
iarray_high = size(array, kind = mrki64)

if (multirand_gaussian64buf_filled) then
  array(iarray_low) = multirand_gaussian64buf
  iarray_low = iarray_low + 1
  multirand_gaussian64buf_filled = .false.
end if
do iarray = iarray_low, iarray_high, 2
  do
    x = multirand_int64() / multirand_max64
    y = multirand_int64() / multirand_max64
    w = x * x + y * y
    if (w > 0.0_mrkr64 .and. w < 1.0_mrkr64) exit
  end do
  w = sqrt((-2.0_mrkr64 * log(w)) / w)
  array(iarray) = x * w
  if (iarray < iarray_high) then
    array(iarray + 1) = y * w
  else
    multirand_gaussian64buf = y * w
    multirand_gaussian64buf_filled = .true.
  end if
end do

end subroutine multirand_gaussian_array64


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill array with 32-bit standard Gaussian random numbers !
! using George Marsaglia's polar method                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multirand_gaussian_array32(array)
implicit none

! array to fill in random numbers
real(kind = mrkr32), dimension(:), intent(out) :: array

integer(kind = mrki64) :: iarray_low, iarray_high, iarray, irand64
real(kind = mrkr32) :: x, y, w

iarray_low = 1
iarray_high = size(array, kind = mrki64)

if (multirand_gaussian32buf_filled) then
  array(iarray_low) = multirand_gaussian32buf
  iarray_low = iarray_low + 1
  multirand_gaussian32buf_filled = .false.
end if
do iarray = iarray_low, iarray_high, 2
  do
    irand64 = multirand_int64()
    x = INT64TO32_1(irand64) / multirand_max32
    y = INT64TO32_2(irand64) / multirand_max32
    w = x * x + y * y
    if (w > 0.0_mrkr32 .and. w < 1.0_mrkr32) exit
  end do
  w = sqrt((-2.0_mrkr32 * log(w)) / w)
  array(iarray) = x * w
  if (iarray < iarray_high) then
    array(iarray + 1) = y * w
  else
    multirand_gaussian32buf = y * w
    multirand_gaussian32buf_filled = .true.
  end if
end do

end subroutine multirand_gaussian_array32


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! George Marsaglia's 64-bit KISS                                             !
! https://groups.google.com/d/msg/comp.lang.fortran/qFv18ql_WlU/IK8KGZZFJx4J !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind = mrki64) function multirand_kiss64()
implicit none

integer(kind = mrki64) :: t

#define S(x) ishft(x, -63)

t = ishft(multirand_seeds(0), 58) + multirand_seeds(3)
if (S(multirand_seeds(0)) == S(t)) then
  multirand_seeds(3) = ishft(multirand_seeds(0), -6) &
    + S(multirand_seeds(0))
else
  multirand_seeds(3) = ishft(multirand_seeds(0), -6) &
    - S(multirand_seeds(0) + t) + 1
end if
multirand_seeds(0) = multirand_seeds(0) + t
multirand_seeds(1) = XORSHFT(multirand_seeds(1), 13)
multirand_seeds(1) = XORSHFT(multirand_seeds(1), -17)
multirand_seeds(1) = XORSHFT(multirand_seeds(1), 43)
multirand_seeds(2) = 6906969069_mrki64 * multirand_seeds(2) + 1234567
multirand_kiss64 = multirand_seeds(0) + multirand_seeds(1) + multirand_seeds(2)

#undef S

end function multirand_kiss64


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 64-bit Mersenne Twister 19937                              !
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind = mrki64) function multirand_mt19937_64() result(x)
implicit none

integer, parameter :: &
  nn = 312, &
  mm = 156
integer(kind = mrki64), parameter :: &
  um = -2147483648_mrki64, & ! 0xFFFFFFFF80000000
  lm = 2147483647_mrki64, &  ! 0x000000007FFFFFFF
  and1 = 6148914691236517205_mrki64, & ! 0x5555555555555555
  and2 = 8202884508482404352_mrki64, & ! 0x71D67FFFEDA60000
  and3 = -2270628950310912_mrki64      ! 0xFFF7EEE000000000
integer(kind = mrki64), dimension(0 : 1), parameter :: mag01 = &
  (/ 0_mrki64, -5403634167711393303_mrki64 & ! MATRIX_A 0xB5026F5AA96619E9ULL
  /)

! refill seed array
if (multirand_iseed >= nn) then
  do multirand_iseed = 0, nn - mm - 1
    x = ior(iand(multirand_seeds(multirand_iseed), um), &
      iand(multirand_seeds(multirand_iseed + 1), lm))
    multirand_seeds(multirand_iseed) = ieor( &
      ieor(multirand_seeds(multirand_iseed + mm), ishft(x, -1)), &
      mag01(iand(x, 1_mrki64)))
  end do
  do multirand_iseed = nn - mm, nn - 2
    x = ior(iand(multirand_seeds(multirand_iseed), um), &
      iand(multirand_seeds(multirand_iseed + 1), lm))
    multirand_seeds(multirand_iseed) = ieor( &
      ieor(multirand_seeds(multirand_iseed + (mm - nn)), ishft(x, -1)), &
      mag01(iand(x, 1_mrki64)))
  end do
  x = ior(iand(multirand_seeds(nn - 1), um), iand(multirand_seeds(0), lm))
  multirand_seeds(nn - 1) = ieor( &
    ieor(multirand_seeds(mm - 1), ishft(x, -1)), &
    mag01(iand(x, 1_mrki64)))
  multirand_iseed = 0
end if

x = XORANDSHFT(multirand_seeds(multirand_iseed), -29, and1)
x = XORANDSHFT(x, 17, and2)
x = XORANDSHFT(x, 37, and3)
x = XORSHFT(x, -43)
multirand_iseed = multirand_iseed + 1

end function multirand_mt19937_64


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! George Marsaglia's 64-bit SuperKISS                    !
! http://mathforum.org/kb/message.jspa?messageID=6914945 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind = mrki64) function multirand_superkiss64()
implicit none

integer, parameter :: nn = 20632, &
  icarry = nn, &
  ixcng = nn + 1, &
  ixs = nn + 2
!integer(kind = mrki64)

integer(kind = mrki64) :: s, z, h

! refill seed array
if (multirand_iseed >= nn) then
  do multirand_iseed = 0, nn - 1
    h = iand(multirand_seeds(icarry), 1_mrki64)
    z = ishft(ishft(multirand_seeds(multirand_iseed), 41), -1) &
      + ishft(ishft(multirand_seeds(multirand_iseed), 39), -1) &
      + ishft(multirand_seeds(icarry), -1)
    multirand_seeds(icarry) = ishft(multirand_seeds(multirand_iseed), -23) &
      + ishft(multirand_seeds(multirand_iseed), -25) &
      + ishft(z, -63)
    multirand_seeds(multirand_iseed) = not(ishft(z, 1) + h)
  end do
  multirand_iseed = 0
end if

multirand_seeds(ixcng) = multirand_seeds(ixcng) * 6906969069_mrki64 + 123
multirand_seeds(ixs) = XORSHFT(multirand_seeds(ixs), 13)
multirand_seeds(ixs) = XORSHFT(multirand_seeds(ixs), -17)
multirand_seeds(ixs) = XORSHFT(multirand_seeds(ixs), 43)

multirand_superkiss64 = multirand_seeds(multirand_iseed) &
  + multirand_seeds(ixcng) + multirand_seeds(ixs)
multirand_iseed = multirand_iseed + 1

end function multirand_superkiss64

end module multirand

