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


! module for managing input parameters and functions
! change input parameters in this file
module pic1dp_input
use pic1dp_global
implicit none
#include "finclude/petscdef.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!
! termination conditions !
!!!!!!!!!!!!!!!!!!!!!!!!!!
! whichever of the following conditions is reached, the code stops

! maximum # of time steps
PetscInt, parameter :: input_ntime_max = 900000

! maximum physical time (normalized by 1 / omega_pe)
PetscReal, parameter :: input_time_max = 500.0_kpr


!!!!!!!!!!!!!!!!!!!!!!!
! physical parameters !
!!!!!!!!!!!!!!!!!!!!!!!

! linear or nonlinear run. 0: nonlinear; 1: linear.
PetscInt, parameter :: input_linear = 0

! length in real space (normalized by electron Debye length)
PetscReal, parameter :: &
  input_lx = 2.0_kpr * 3.1415926535897932384626_kpr / 0.36_kpr

! equilibrium particle velocity distribution.
! 0: (shifted) Maxwellian; 1: two-stream1; 2: two-stream2; 3: bump-on-tail
! two-stream1: f(v) = 1/sqrt(2 pi) * v**2 * exp(-v**2 / 2)
! two-stream2: f(v) = 1/2 * (fm(v - v0) + fm(v + v0)); fm: Maxwellian
! bump-on-tail: f(v) = density * fm(v; T) + (1 - density) * fm(v - v0; T2)
PetscInt, parameter :: input_iptcldist = 3

! # of particle species
PetscInt, parameter :: input_nspecies = 1

! particle charge (normalized by proton charge e),
! mass (normalized by electron mass),
! temperature (normalized by electron temperature),
! density (normalized by electron equilibrium density),
! and equilibrium flow (normalized by electron thermal velocity)
! temperature and equilibrium flow are not used for input_iptcldist = 1
! temperature2 is the beam temperature for input_iptcldist = 3
PetscReal, dimension(input_nspecies), parameter :: &
  input_species_charge = (/ -1.0_kpr /), &
  input_species_mass = (/ 1.0_kpr /), &
  input_species_temperature = (/ 1.0_kpr /), &
  input_species_temperature2 = (/ 1.0_kpr /), &
  input_species_density = (/ 0.9_kpr /), &
  input_species_v0 = (/ 5.0_kpr /)

! # of modes kept
PetscInt, parameter :: input_nmode = 1

! modes kept
! for each mode, the number here gives the # of mode periods in real space
PetscInt, dimension(0 : input_nmode - 1), parameter :: &
  input_modes = (/ 1 /)


!!!!!!!!!!!!!!!!!!!!!
! initial condition !
!!!!!!!!!!!!!!!!!!!!!
! # of modes
PetscInt, parameter :: input_init_nmode = 1

! initial modes
PetscInt, dimension(0 : input_init_nmode - 1), parameter :: &
  input_init_mode = (/ 1 /)

! initial amplitude of each mode (cos and sin part)
! note that the initial distribution in velocity space depends on
! initial perturbation shape function input_pertb_shape(v)
PetscScalar, dimension(0 : input_init_nmode - 1), parameter :: &
  input_init_mode_cos = (/ 0.00_kpr /), &
  input_init_mode_sin = (/ 1e-5_kpr /)


!!!!!!!!!!!!!!!!!!!!!!!!
! numerical parameters !
!!!!!!!!!!!!!!!!!!!!!!!!

! delta f or full-f. 0: full-f; 1: delta f.
PetscInt, parameter :: input_deltaf = 1

! initial time step (normalized by 1 / omega_pe)
PetscReal, parameter :: input_dt = 0.05_kpr

! allocation for # of marker particles per species
! this is also the maximum allowed # of marker particles during simulation
PetscInt, parameter :: input_nparticle_max = 6400000

! # of initial loaded particles for each species
PetscInt, dimension(input_nspecies), parameter :: &
  input_species_nparticle_init = (/ 6400000 /)

! marker distribution in velocity space:
! 1: same as physical distribution; 2: uniform
! if set to 1, now only input_iptcldist = 0 (shifted Maxwellian) is implemented
PetscInt, parameter :: input_imarker = 2

! maximum velocity for uniform loading and for output
PetscReal, parameter :: input_v_max = 8.0_kpr

! # of grid points in real space
PetscInt, parameter :: input_nx = 192

! # of grid points in velocity space for resonant detection
PetscInt, parameter :: input_nv = 128

! particle shape calculation
! 1: use PETSc matrix for shape matrix, destroy and re-create matrix every time
! 2: use PETSc matrix for shape matrix, use MatZeroEntries to reset matrix
! 3: use regular array for shape matrix
! 4: calculate shape function as needed, no matrix
PetscInt, parameter :: input_iptclshape = 4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! marker particle optimizations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! # of times of merging particles
! (set to 0 to disable merging)
PetscInt, parameter :: input_nmerge = 0

! list of times to merge particles, must be in ascending order
PetscReal, dimension(input_nmerge), parameter :: &
  input_tmerge = (/ (50.0_kpr + global_itime * 0.5_kpr, &
    global_itime = 1, input_nmerge) /)

! list of merging thresholds in terms of
! fraction of absolute value of distribution in v
PetscReal, dimension(input_nmerge), parameter :: &
  input_thshmerge = (/ &
    (0.1_kpr / max(input_nmerge, 1) * real(global_itime, kpr), &
    global_itime = 1, input_nmerge) /)

! # of times of throwing away particles
! (set to 0 to disable throwing away)
PetscInt, parameter :: input_nthrowaway = 0

! list of times to throw away particles, must be in ascending order
PetscReal, dimension(input_nthrowaway), parameter :: &
  input_tthrowaway = (/ (50.0_kpr + global_itime * 0.5_kpr, &
    global_itime = 1, input_nthrowaway) /)

! type of throwing away
! 1: throw away based on threshold given by input_thshthrowaway
! 2: throw away based on profile of int |delta f| dx
PetscInt, parameter :: input_typethrowaway = 2

! list of throwing away thresholds in terms of
! fraction of absolute value of distribution in v
! only useful when input_typethrowaway == 1
PetscReal, dimension(input_nthrowaway), parameter :: &
  input_thshthrowaway = (/ &
    (0.1_kpr / max(input_nthrowaway, 1) * real(global_itime, kpr), &
    global_itime = 1, input_nthrowaway) /)

! fraction of not important particles to be thrown away
! only useful when input_typethrowaway == 1
PetscReal, parameter :: input_throwaway_frac = 0.9_kpr

! # of times of splitting particles
! (set to 0 to disable splitting)
PetscInt, parameter :: input_nsplit = 0

! list of times to split particles, must be in ascending order
PetscReal, dimension(input_nsplit), parameter :: &
  input_tsplit = (/ (50.0_kpr + global_itime * 0.5_kpr, &
    global_itime = 1, input_nsplit) /)

! list of splitting thresholds in terms of
! fraction of absolute value of distribution in v
PetscReal, dimension(input_nsplit), parameter :: &
  input_thshsplit = (/ &
    (1.0_kpr - 0.9_kpr / max(input_nsplit, 1) * real(global_itime, kpr), &
    global_itime = 1, input_nsplit) /)

! # of groups for splitting particles
PetscInt, parameter :: input_split_ngroup = 5

! variation of change in v in terms of fraction of grid size
PetscReal, parameter :: input_split_dv_sig_frac = 0.1_kpr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! random number generator (multirand) parameters !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! algorithm for random integer (engine for all types of random numbers)
! 1: George Marsaglia's 64-bit KISS (period ~ 2**(247.42) ~ 10**(74.48))
! 2: 64-bit Mersenne Twister 19937 (period = 2**19937 - 1 ~ 10**6001)
! 3: George Marsaglia's 64-bit SuperKISS
!   (period = 5*2**1320480*(2**64-1) ~ 10**397524)
PetscInt, parameter :: input_multirand_al_int = 3

! random seed type
! 1: constant seeds
! 2: random seeds using system_clock
! 3: random seeds using /dev/urandom (if fail, fall back to system_clock)
PetscInt, parameter :: input_multirand_seed_type = 3

! # of warm up rounds
PetscInt, parameter :: input_multirand_warmup = 5

! whether to run self test during initialization of multirand
! after a few runs if you don't see warnings from multirand_selftest
!   in stdout, then you can turn this off
! once CPU, OS or operating system is changed, you should turn this back on
!   for a few runs to test if the random number generator is working correctly
logical, parameter :: input_multirand_selftest = .true.


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
PetscReal, parameter :: input_output_interval = 0.5_kpr

! # of x grids for output of particle distribution
PetscInt, parameter :: input_nx_opd = 64

! # of v grids for output of particle distribution
PetscInt, parameter :: input_nv_opd = 64

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initial perturbation shpae in velocity space !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PetscScalar function input_pertb_shape(v, ispecies)
!use multirand
implicit none
PetscScalar, intent(in) :: v
PetscInt, intent(in) :: ispecies

! for constant 1.0, the perturbation shape is the same as the marker particle
! distribution in velocity space
input_pertb_shape = 1.0_kpr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! below are other testing shapes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! random noise
!input_pertb_shape = multirand_real64() - 0.5_kpr

!input_pertb_shape = v**30 * exp(-v * v) * 1e-8_kpr

end function input_pertb_shape


!!!!!!!!!!!!!!!!!!!!!!!!
! check input validity !
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine input_init
implicit none
#include "finclude/petsc.h90"

! check input parameters
if (&
  input_iptcldist >= 1 .and. input_imarker == 1 &
) then
  call global_pp("Error: case of input_iptcldist >= 1 and input_")
  call global_pp("imarker = 1 not implemented yet.\n")
  call PetscFinalize(global_ierr)
  CHKERRQ(global_ierr)
  stop 1
end if
if (input_linear == 1 .and. input_deltaf == 0) then
  call global_pp("Error: case of input_linear = 1 and input_deltaf = 0 no")
  call global_pp("t implemented yet.\n")
  call PetscFinalize(global_ierr)
  CHKERRQ(global_ierr)
  stop 1
end if
end subroutine input_init

end module pic1dp_input

