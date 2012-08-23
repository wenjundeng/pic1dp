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
#include "finclude/petsc.h90"

character(len = 25), parameter :: version = '2012-07-22 19:33:37-04:00'

! wall clock timer indexes
PetscInt, parameter :: &
  iwt_total = 1, &
  iwt_init = 2, &
  iwt_particle_load = 3, &
  iwt_push_particle = 4, &
  iwt_particle_shape = 5, &
  iwt_collect_charge = 6, &
  iwt_field_electric = 7, &
  iwt_output = 8, &
  iwt_final = 9

PetscInt :: nrk ! # of Runge-Kutta sub-steps
PetscInt :: irk ! indexing Runge-Kutta sub-steps

! termination related variables
PetscInt :: itermination ! status of termination condition: 0: not to terminate; 1: to terminate

! for wall clock timers
character(len = 18), dimension(30) :: string_wt


call PetscInitialize(PETSC_NULL_CHARACTER, global_ierr)
CHKERRQ(global_ierr)

call wtimer_start(iwt_total) ! start recording total time

! initialization
call wtimer_start(iwt_init)
call MPI_Comm_rank(MPI_COMM_WORLD, global_mype, global_ierr)
CHKERRQ(global_ierr)
if (input_verbosity >= 1) call global_pp("PIC1D-PETSc version " // version // "\n")
call particle_init
call field_init
call output_init
call wtimer_stop(iwt_init)

! load particles
call wtimer_start(iwt_particle_load)
call particle_load
call wtimer_stop(iwt_particle_load)
if (input_iptclshape < 4) then
  call wtimer_start(iwt_particle_shape)
  call particle_compute_shape_x
  call wtimer_stop(iwt_particle_shape)
end if

global_itime = 0
global_time = 0.0_kpr

! solve initial field
! collect charges
call wtimer_start(iwt_collect_charge)
call interaction_collect_charge
call wtimer_stop(iwt_collect_charge)
! solve electric field
call wtimer_start(iwt_field_electric)
call field_solve_electric
call wtimer_stop(iwt_field_electric)

if (input_verbosity == 1) then
  call global_pp("Info: progress:\n")
  call global_pp("progrss  itime     time  int E^2 dx\n")
end if

! output for 0th time step
call wtimer_start(iwt_output)
call output_all
call wtimer_stop(iwt_output)

!call field_test

nrk = 2 ! 2-step 2nd-order Runge-Kutta
itermination = check_termination()
do while (itermination == 0)
  do irk = 1, nrk
    ! push particles
    call wtimer_start(iwt_push_particle)
    call interaction_push_particle(irk)
    call wtimer_stop(iwt_push_particle)

    if (input_iptclshape < 4) then
      !call wtimer_start(21)
      !call particle_zeroshape_x
      !call wtimer_stop(21)
      ! update particle shape matrix
      call wtimer_start(iwt_particle_shape)
      call particle_compute_shape_x
      call wtimer_stop(iwt_particle_shape)
    end if

    ! collect charges
    call wtimer_start(iwt_collect_charge)
    call interaction_collect_charge
    call wtimer_stop(iwt_collect_charge)

    ! solve electric field
    call wtimer_start(iwt_field_electric)
    call field_solve_electric
    call wtimer_stop(iwt_field_electric)
  end do
  ! update time
  global_itime = global_itime + 1
  global_time = global_time + input_dt

  ! check if needs to terminate
  itermination = check_termination()

  ! output
  if ( &
    mod( &
      global_time + PETSC_SQRT_MACHINE_EPSILON, input_output_interval &
    ) < mod( &
      global_time + PETSc_SQRT_MACHINE_EPSILON - input_dt, input_output_interval &
    ) & ! time just passed a full interval
    .or. itermination == 1 &
  ) then
    call wtimer_start(iwt_output)
    call output_all
    call wtimer_stop(iwt_output)
  end if
end do

! finalization
call wtimer_start(iwt_final)
call particle_final
call field_final
call output_final
call wtimer_stop(iwt_final)

call wtimer_stop(iwt_total) ! stop recording total time

! print timers
if (input_verbosity >= 1) then
  call global_pp("Info: timers:\n")
  call global_pp( &
    "            total   initialization    particle load    push particle\n")
  call wtimer_print(iwt_total, string_wt(iwt_total), iwt_total, .true.)
  call wtimer_print(iwt_init, string_wt(iwt_init), iwt_total, .true.)
  call wtimer_print(iwt_particle_load, string_wt(iwt_particle_load), &
    iwt_total, .true.)
  call wtimer_print(iwt_push_particle, string_wt(iwt_push_particle), iwt_total, .true.)
  call global_pp(string_wt(iwt_total) // string_wt(iwt_init) &
    // string_wt(iwt_particle_load) // string_wt(iwt_push_particle) // "\n")

  call global_pp( &
    "   particle shape   collect charge   electric field           output\n")
  call wtimer_print(iwt_particle_shape, string_wt(iwt_particle_shape), &
    iwt_total, .true.)
  call wtimer_print(iwt_collect_charge, string_wt(iwt_collect_charge), &
    iwt_total, .true.)
  call wtimer_print(iwt_field_electric, string_wt(iwt_field_electric), &
    iwt_total, .true.)
  call wtimer_print(iwt_output, string_wt(iwt_output), iwt_total, .true.)
  call global_pp(string_wt(iwt_particle_shape) // string_wt(iwt_collect_charge) &
    // string_wt(iwt_field_electric) // string_wt(iwt_output) // "\n")

  call global_pp( &
    "     finalization    MPI_Allreduce          scatter                .\n")
  call wtimer_print(iwt_final, string_wt(iwt_final), iwt_total, .true.)
  call wtimer_print(21, string_wt(21), iwt_total, .true.)
  call wtimer_print(22, string_wt(22), iwt_total, .true.)
  call global_pp(string_wt(iwt_final) // &
    string_wt(21) // string_wt(22) // "\n")
!  call global_pp(string_wt(iwt_final) // "\n")
end if

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
  .or. global_time >= input_time_max &
) then
  check_termination = 1
else
  check_termination = 0
end if

end function check_termination

end program pic1dp

