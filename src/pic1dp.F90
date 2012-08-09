! PIC1D-PETSc
program pic1dp
use wtimer
use pic1dp_global
use pic1dp_input
use pic1dp_particle
use pic1dp_output
implicit none
#include "finclude/petsc.h90"

! wall clock timer indexes
PetscInt, parameter :: &
  iwt_total = 1, &
  iwt_particle_load = 2, &
  iwt_particle_push = 3, &
  iwt_field = 4, &
  iwt_output = 5

PetscInt :: nrk ! # of Runge-Kutta sub-steps
PetscInt :: irk ! indexing Runge-Kutta sub-steps

PetscInt :: iconvergence ! status of convergence: 0: not convergent; 1: convergent

PetscInt :: itime
PetscReal :: time
PetscInt :: nimit, nimit_max ! # and max # of implicit iteration per time step
PetscInt :: nimit_dt ! # of implicit iteration per trial of time step
PetscReal :: nimit_avg ! average # of implicit iteration per time step

! termination related variables
PetscInt :: itermination ! status of termination condition: 0: not to terminate; 1: to terminate

character(len = 18), dimension(4) :: string_wt


call PetscInitialize(PETSC_NULL_CHARACTER, global_ierr)
CHKERRQ(global_ierr)

call wtimer_start(iwt_total) ! start recording total time

call MPI_Comm_rank(MPI_COMM_WORLD, global_mype, global_ierr)
CHKERRQ(global_ierr)

if (input_verbosity >= 1) call global_pp("PIC1D-PETSc\n")

call particle_init
call output_init

call wtimer_start(iwt_particle_load)
call particle_load
call wtimer_stop(iwt_particle_load)

call wtimer_start(iwt_output)
call output_xv_ptcldist
call wtimer_stop(iwt_output)

call particle_final
call output_final

call wtimer_stop(iwt_total) ! stop recording total time

if (input_verbosity >= 1) then
  call global_pp("            total    particle load           output\n")
  call wtimer_print(iwt_total, string_wt(1), iwt_total, .true.)
  call wtimer_print(iwt_particle_load, string_wt(2), iwt_total, .true.)
  call wtimer_print(iwt_output, string_wt(3), iwt_total, .true.)
  call global_pp(string_wt(1) // string_wt(2) // string_wt(3) // "\n")
end if

call PetscFinalize(global_ierr)
CHKERRQ(global_ierr)

end program pic1dp

