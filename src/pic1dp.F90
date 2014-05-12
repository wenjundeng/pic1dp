! Copyright 2012-2014 Wenjun Deng <wdeng@wdeng.info>
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


! PIC1D-PETSc
program pic1dp
use wtimer
use pic1dp_global
use pic1dp_input
use pic1dp_particle
use pic1dp_field
use pic1dp_interaction
use pic1dp_output
implicit none

! the following line is to work around a bug in PETSc 3.3-p2 and before
!#include "finclude/petsctsdef.h"

#include "finclude/petsc.h90"

character(len = 25), parameter :: version = '2014-05-12 12:03:27-04:00'

! status of termination condition: 0: not to terminate; 1: to terminate
PetscInt :: itermination

logical :: flag_optimized


call PetscInitialize(PETSC_NULL_CHARACTER, global_ierr)
CHKERRQ(global_ierr)

call wtimer_start(global_iwt_total) ! start recording total time

! initialization
call wtimer_start(global_iwt_init)
call MPI_Comm_rank(MPI_COMM_WORLD, global_mype, global_ierr)
CHKERRQ(global_ierr)
call MPI_Comm_size(MPI_COMM_WORLD, global_npe, global_ierr)
CHKERRQ(global_ierr)
if (input_verbosity >= 1) &
  call global_pp("PIC1D-PETSc version " // version // "\n")

call input_init
call particle_init
call field_init
call output_init
call wtimer_stop(global_iwt_init)

! load particles
call particle_load
if (input_iptclshape < 4) call particle_compute_shape_x

global_itime = 0
global_time = 0.0_kpr

! solve initial field
call interaction_collect_charge ! collect charges
call field_solve_electric ! solve electric field
call output_progress(0) ! output progress header
call output_all ! output for 0th time step
!call field_test

itermination = check_termination()
do while (itermination == 0) ! main time evolution loop
  do global_irk = 1, 2
    call interaction_push_particle ! push particles

    call particle_optimize(flag_optimized) ! particle optimization
    if (flag_optimized) call output_progress(2)

    ! update particle shape matrix
    if (input_iptclshape < 4) call particle_compute_shape_x

    call interaction_collect_charge ! collect charges
    call field_solve_electric ! solve electric field
  end do
  ! update time
  global_itime = global_itime + 1
  global_time = global_time + input_dt

  itermination = check_termination() ! check if needs to terminate

  ! output
  if ( &
    mod( &
      global_time + PETSC_SQRT_MACHINE_EPSILON, input_output_interval &
    ) < mod( &
      global_time + PETSC_SQRT_MACHINE_EPSILON - input_dt, &
      input_output_interval &
    ) & ! time just passed a full interval
    .or. itermination == 1 & ! final time step
  ) then
    call output_all
  end if
end do

! finalization
call wtimer_start(global_iwt_final)
call particle_final
call field_final
call output_final
call wtimer_stop(global_iwt_final)

call wtimer_stop(global_iwt_total) ! stop recording total time

! print out wall clock timers
! no problem for being after output_final
call output_wtimer

call PetscFinalize(global_ierr)
CHKERRQ(global_ierr)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check termination condition !
! return termination flag     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PetscInt function check_termination()
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

if ( &
  global_itime >= input_ntime_max &
  .or. global_time + PETSC_SQRT_MACHINE_EPSILON >= input_time_max &
) then
  check_termination = 1
else
  check_termination = 0
end if

end function check_termination

end program pic1dp

