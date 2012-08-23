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


! a wall clock timer using MPI_Wtime()

! if this module is used in a PETSc program,
! use -D__PETSc option when compiling, or uncomment the following line
! #define __PETSc

module wtimer
implicit none

! real kind # used in this module
#ifdef __PETSc
#include "finclude/petscdef.h"
PetscReal, parameter :: wtimer_testkind = 0d0
PetscInt, parameter :: wtrk = kind(wtimer_testkind)
PetscInt, parameter :: wtik = kind(wtrk)
#else
integer, parameter :: wtrk = selected_real_kind(10)
integer, parameter :: wtik = selected_int_kind(10)
#endif

! maximum # of timers
integer(kind = wtik), parameter :: wtimer_nwt_max = 40

! timers and start time of timers
real(kind = wtrk), dimension(wtimer_nwt_max) :: &
  wtimer_wt = 0.0_wtrk, wtimer_wt_start = 0.0_wtrk

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start timer with index iwt !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wtimer_start(iwt)
implicit none
#ifdef __PETSc
#include "finclude/petsc.h90"
#else
#include "mpif.h"
#endif

! timer index
integer(kind = wtik), intent(in) :: iwt

! iwt out of range
if (iwt < 1 .or. iwt > wtimer_nwt_max) return

wtimer_wt_start(iwt) = MPI_Wtime()

end subroutine wtimer_start


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! stop timer with index iwt !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wtimer_stop(iwt)
implicit none
#ifdef __PETSc
#include "finclude/petsc.h90"
#else
#include "mpif.h"
#endif

! timer index
integer(kind = wtik), intent(in) :: iwt

! iwt out of range
if (iwt < 1 .or. iwt > wtimer_nwt_max) return

wtimer_wt(iwt) = wtimer_wt(iwt) + MPI_Wtime() - wtimer_wt_start(iwt)

end subroutine wtimer_stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print timer result to a string !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wtimer_print(iwt, string, iwt_percent, escape_percent)
implicit none

! timer index
integer(kind = wtik), intent(in) :: iwt

! output string
character(len = *), intent(out) :: string

! print the percentage of time of iwt in time of iwt_percent
! if not specified, then disable printing percentage
integer(kind = wtik), intent(in), optional :: iwt_percent

! whether to escape the % symbol by %%
! useful if the string will be passed to PetscPrintf
! default to false
logical, intent(in), optional :: escape_percent

character(len = 3) :: padding

! iwt out of range
if (iwt < 1 .or. iwt > wtimer_nwt_max) return

if (present(iwt_percent)) then
  ! iwt_percent out of range
  if (iwt_percent < 1 .or. iwt_percent > wtimer_nwt_max) return

  if (present(escape_percent)) then
    if (escape_percent) then
      padding = '%%)'
    else
      padding = '%)'
    end if
  else
    padding = '%)'
  end if

  write (string, '(2a, f5.1, a)') wtimer_sec2text(wtimer_wt(iwt)), '(', &
    wtimer_wt(iwt) / wtimer_wt(iwt_percent) * 100.0_wtrk, trim(padding)
else
  write (string, '(a)') wtimer_sec2text(wtimer_wt(iwt))
end if

end subroutine wtimer_print


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! converts # of seconds to human readable text !
! e.g., 5.2min, 3.2hr, 4.1day                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len = 9) function wtimer_sec2text(sec)
implicit none

real(kind = wtrk), intent(in) :: sec 
real(kind = wtrk) :: numtime

character(len = 3) :: suffix
character(len = 9) :: text

suffix = 'sec'
numtime = sec
if (numtime > 60.0_wtrk) then
  numtime = numtime / 60.0_wtrk
  suffix = 'min'
  if (numtime > 60.0_wtrk) then
    numtime = numtime / 60.0_wtrk
    suffix = 'hr'
    if (numtime > 24.0_wtrk) then
      numtime = numtime / 24.0_wtrk
      suffix = 'day'
    end if
  end if
end if
write (text, '(f6.1, a)') numtime, suffix
wtimer_sec2text = adjustr(text)

end function wtimer_sec2text

end module wtimer

