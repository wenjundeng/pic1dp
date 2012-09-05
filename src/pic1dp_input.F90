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
PetscInt, parameter :: input_linear = 0

! length in real space (normalized by electron Debye length)
!PetscReal, parameter :: input_lx = 15.708  ! k = 0.4
PetscReal, parameter :: input_lx = 4.0_kpr * 3.141592653589793238_kpr

! equilibrium particle velocity distribution.
! 0: (shifted) Maxwellian; 1: two-stream1; 2: two-stream2
! two-stream1: f(v) = 1/sqrt(2 pi) * v**2 * exp(-v**2 / 2)
! two-stream2: f(v) = 1/2 * (fm(v - v0) + fm(v + v0)); fm: Maxwellian
PetscInt, parameter :: input_iptcldist = 2

! # of particle species
PetscInt, parameter :: input_nspecies = 1

! particle charge (normalized by proton charge e),
! mass (normalized by electron mass),
! temperature (normalized by electron temperature),
! density (normalized by electron equilibrium density),
! and equilibrium flow (normalized by electron thermal velocity)
! temperature and equilibrium flow are not used for input_iptcldist = 1
PetscReal, dimension(input_nspecies), parameter :: &
  input_species_charge = (/ -1.0_kpr /), &
  input_species_mass = (/ 1.0_kpr /), &
  input_species_temperature = (/ 0.5_kpr /), &
  input_species_density = (/ 1.0_kpr /), &
  input_species_v0 = (/ sqrt(2.0_kpr) /)

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
! marker distribution input_imarker and initial perturbation shape function
! input_pertb_shape(v)
PetscScalar, dimension(0 : input_init_nmode - 1), parameter :: &
  input_init_mode_cos = (/ 0e-5_kpr /), &
  input_init_mode_sin = (/ 1e-4_kpr /)


!!!!!!!!!!!!!!!!!!!!!!!!
! numerical parameters !
!!!!!!!!!!!!!!!!!!!!!!!!

! delta f or full-f. 0: full-f; 1: delta f.
PetscInt, parameter :: input_deltaf = 1

! initial time step (normalized by 1 / omega_pe)
PetscReal, parameter :: input_dt = 0.1_kpr

! # of marker particles per species
PetscInt, parameter :: input_nparticle = 3200000

! marker distribution in velocity space:
! 1: same as physical distribution; 2: uniform
! if set to 1, now only input_iptcldist = 0 (shifted Maxwellian) is implemented
PetscInt, parameter :: input_imarker = 2

! maximum velocity for uniform loading and for output
PetscReal, parameter :: input_v_max = 6.0_kpr

! # of grid points in real space
PetscInt, parameter :: input_nx = 64

! # of grid points in velocity space for output and resonant detection
PetscInt, parameter :: input_nv = 128

! # of times of merging particles (set to 0 to disable merging)
PetscInt, parameter :: input_nmerge = 30

! list of times to merge particles, must be in ascending order
PetscReal, dimension(input_nmerge), parameter :: &
!  input_tmerge = 0.0_kpr ! if input_nmerge = 0, use this line
!  input_tmerge = (/ 30.0_kpr /)
!  input_tmerge = (/ 30.0_kpr, 32.0_kpr, 34.0_kpr, 36.0_kpr /)
  input_tmerge = (/ 20.0_kpr, 21.0_kpr, 22.0_kpr, 23.0_kpr, 24.0_kpr, &
    25.0_kpr, 26.0_kpr, 27.0_kpr, 28.0_kpr, 29.0_kpr, &
    30.0_kpr, 31.0_kpr, 32.0_kpr, 33.0_kpr, 34.0_kpr, &
    35.0_kpr, 36.0_kpr, 37.0_kpr, 38.0_kpr, 39.0_kpr, &
    40.0_kpr, 41.0_kpr, 42.0_kpr, 43.0_kpr, 44.0_kpr, &
    45.0_kpr, 46.0_kpr, 47.0_kpr, 48.0_kpr, 49.0_kpr /)

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
PetscReal, parameter :: input_output_interval = 0.5_kpr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initial perturbation shpae in velocity space !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PetscScalar function input_pertb_shape(v, ispecies)
implicit none
PetscScalar, intent(in) :: v
PetscInt, intent(in) :: ispecies

PetscScalar :: omega_r, omega_i, k

! for constant 1.0, the perturbation shape is the same as the marker particle
! distribution in velocity space
input_pertb_shape = 1.0_kpr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! below are other testing shapes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! random noise
!call random_number(input_pertb_shape)
!input_pertb_shape = input_pertb_shape - 0.5

!input_pertb_shape = v**30 * exp(-v * v) * 1e-8_kpr

!k = 0.4_kpr
!omega_r = 1.10625_kpr
!omega_i = 0.0136088_kpr

!input_pertb_shape = (omega_r - k * v) / ((omega_r - k * v)**2 + omega_i**2) &
!  * (v - input_species_v0(ispecies)) / input_species_temperature(ispecies)
end function input_pertb_shape

end module pic1dp_input

