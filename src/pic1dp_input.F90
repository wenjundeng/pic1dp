! input parameters and functions
module pic1dp_input
use pic1dp_global, only: kpr
implicit none
#include "finclude/petscdef.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!
! termination conditions !
!!!!!!!!!!!!!!!!!!!!!!!!!!
! whichever of the following conditions is reached, the code stops

! maximum # of time steps
PetscInt, parameter :: input_ntime_max = 900000

! maximum physical time (normalized by 1 / omega_pe)
PetscReal, parameter :: input_time_max = 100.0_kpr


!!!!!!!!!!!!!!!!!!!!!!!
! physical parameters !
!!!!!!!!!!!!!!!!!!!!!!!

! linear or nonlinear run. 0: nonlinear; 1: linear.
PetscInt, parameter :: input_linear = 1

! # of particle species
PetscInt, parameter :: input_nspecies = 1

! particle charge (normalized by proton charge e),
! and mass (normalized by electron mass)
PetscReal, dimension(input_nspecies), parameter :: &
  input_charge = (/ -1.0_kpr /), &
  input_mass = (/ 1.0_kpr /)


!!!!!!!!!!!!!!!!!!!!!!!!
! numerical parameters !
!!!!!!!!!!!!!!!!!!!!!!!!

! initial time step (normalized by 1 / omega_pe)
PetscReal, parameter :: input_dt = 0.1_kpr

! # of marker particles per species
PetscInt, parameter :: input_nparticle = 1000000

! length in real space (normalized by electron Debye length)
PetscReal, parameter :: input_lx = 16.0_kpr

! # of grid points in real space
PetscInt, parameter :: input_nx = 64

! # of modes kept
PetscInt, parameter :: input_nmode = 1

! modes kept
! for each mode, the number here gives the # of mode periods in real space
PetscInt, dimension(input_nmode), parameter :: input_mode = (/ 1 /)

! random seed type
! 1: random seeds (using system_clock)
! 2: constant seeds
PetscInt, parameter :: input_seed_type = 1

!!!!!!!!!!!!!!!!!!!!!
! output parameters !
!!!!!!!!!!!!!!!!!!!!!

! verbose mode, controlling how  much information is shown to screen (stdout)
! 0: not showing any information except warnings and errors
! 1: besides warnings and errors, show essential information about the program
! 2: show additional progress information (only for debugging)
! 3: show additional diagnostic information (only for debugging, printing a
!   lot of variables)
PetscInt, parameter :: input_verbosity = 1

! time interval between data output (normalized by 1 / omega_pe)
PetscReal, parameter :: input_output_interval = 0.5_kpr

! # of velocity grid for output
PetscInt, parameter :: input_output_nv = 128

! maximum velocity in output
PetscReal, parameter :: input_output_v_max = 5.0_kpr

end module pic1dp_input

