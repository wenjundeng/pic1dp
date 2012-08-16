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
PetscReal, parameter :: input_time_max = 50.0_kpr


!!!!!!!!!!!!!!!!!!!!!!!
! physical parameters !
!!!!!!!!!!!!!!!!!!!!!!!

! linear or nonlinear run. 0: nonlinear; 1: linear.
PetscInt, parameter :: input_linear = 1

! length in real space (normalized by electron Debye length)
PetscReal, parameter :: input_lx = 16.0_kpr

! # of particle species
PetscInt, parameter :: input_nspecies = 1

! particle charge (normalized by proton charge e),
! mass (normalized by electron mass),
! and temperature (normalized by electron temperature)
PetscReal, dimension(input_nspecies), parameter :: &
  input_charge = (/ -1.0_kpr /), &
  input_mass = (/ 1.0_kpr /), &
  input_temperature = (/ 1.0_kpr /)

! # of modes kept
PetscInt, parameter :: input_nmode = 2

! modes kept
! for each mode, the number here gives the # of mode periods in real space
PetscInt, dimension(0 : input_nmode - 1), parameter :: &
  input_mode = (/ 1, 2 /)


!!!!!!!!!!!!!!!!!!!!!
! initial condition !
!!!!!!!!!!!!!!!!!!!!!
! # of modes
PetscInt, parameter :: input_init_nmode = 1

! initial modes
PetscInt, dimension(0 : input_init_nmode - 1), parameter :: &
  input_init_mode = (/ 1 /)

! initial amplitude of each mode (cos and sin part)
PetscScalar, dimension(0 : input_init_nmode - 1), parameter :: &
  input_init_mode_cos = (/ 0.0_kpr /), &
  input_init_mode_sin = (/ 1e-1_kpr /)


!!!!!!!!!!!!!!!!!!!!!!!!
! numerical parameters !
!!!!!!!!!!!!!!!!!!!!!!!!

! initial time step (normalized by 1 / omega_pe)
PetscReal, parameter :: input_dt = 0.1_kpr

! # of marker particles per species
PetscInt, parameter :: input_nparticle = 6000000

! # of grid points in real space
PetscInt, parameter :: input_nx = 64

! particle shape calculation
! 1: use PETSc matrix for shape matrix, destroy and re-create matrix every time
! 2: use PETSc matrix for shape matrix, use MatZeroEntries to reset matrix
! 3: use regular array for shape matrix
! 4: calculate shape function as needed, no matrix
PetscInt, parameter :: input_iptclshape = 4

! random seed type
! 1: random seeds (using system_clock)
! 2: constant seeds
PetscInt, parameter :: input_seed_type = 2


!!!!!!!!!!!!!!!!!!!!!
! output parameters !
!!!!!!!!!!!!!!!!!!!!!

! verbose mode, controlling how much information is shown to screen (stdout)
! 0: not showing any information except warnings and errors
! 1: besides warnings and errors, show essential information about the program
! 2: show additional progress information (only for debugging)
! 3: show additional diagnostic information (only for debugging, printing a
!   lot of variables)
PetscInt, parameter :: input_verbosity = 1

! time interval between data output (normalized by 1 / omega_pe)
PetscReal, parameter :: input_output_interval = 1.0_kpr

! # of velocity grid for output
PetscInt, parameter :: input_output_nv = 64

! maximum velocity in output
PetscReal, parameter :: input_output_v_max = 5.0_kpr

end module pic1dp_input

