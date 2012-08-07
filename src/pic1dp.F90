program pic1dp
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

! # of wall clock timers
! 1: total
PetscInt, parameter :: nwtimer = 1

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

! wall clock timers and related variables
PetscInt :: iwtimer
PetscReal :: wt1, wt2, wt3, wt4
PetscReal, dimension(nwtimer) :: wtimer = 0.0_kpr

end program pic1dp

