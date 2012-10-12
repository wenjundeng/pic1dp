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


! module for managing marker particles
module pic1dp_particle
use pic1dp_input
implicit none
#include "finclude/petscdef.h"

! indexing particle optimization operations
PetscInt :: particle_imerge, particle_ithrowaway, particle_isplit

! particle x coordinate, velocity
! equilibrium weight: p = f / g (nonlinear); p = f_0 / g (linear)
! perturbed weight: w = delta f / g
! f is total distribution, delta f is perturbed distribution
! g is marker distribution
! 1: resonant, to be split; 2: non-resonant, to be merged
Vec, dimension(input_nspecies) :: &
  particle_x, particle_v, particle_p, particle_w, &
  particle_x_bak, particle_v_bak, particle_w_bak

! electric field at particle position, need to be used in pic1dp_interaction
Vec :: particle_electric

! temporary particle vectors, need to be used in pic1dp_output
Vec :: particle_tmp1, particle_tmp2

Mat, dimension(input_nspecies) :: particle_shape_x

! shape arrays
PetscInt, dimension(:, :), allocatable :: particle_shape_x_indexes
PetscReal, dimension(:, :), allocatable :: particle_shape_x_values

! local index range of particle
PetscInt :: particle_ip_low, particle_ip_high

! # of valid particles
PetscInt :: particle_np(input_nspecies)

! absolute value of perturbed particle distribution in v
PetscScalar, dimension(input_nspecies, 0 : input_nv - 1) :: &
  particle_dist_pertb_abs_v

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects and other variables !
! (for particle loading, see particle_load)    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_init
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

if (input_nmerge > 0) then
  particle_imerge = 1 
else
  particle_imerge = 0 
end if
if (input_nthrowaway > 0) then
  particle_ithrowaway = 1 
else
  particle_ithrowaway = 0 
end if
if (input_nsplit > 0) then
  particle_isplit = 1 
else
  particle_isplit = 0 
end if

call VecCreate(MPI_COMM_WORLD, particle_x(1), global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(particle_x(1), PETSC_DECIDE, input_nparticle_max, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(particle_x(1), global_ierr)
CHKERRQ(global_ierr)

do ispecies = 2, input_nspecies
  call VecDuplicate(particle_x(1), particle_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do

do ispecies = 1, input_nspecies
  call VecDuplicate(particle_x(1), particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_p(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_w(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call VecDuplicate(particle_x(1), particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_v_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_w_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  if (input_iptclshape <= 2) then
    call global_matcreate(particle_shape_x(ispecies), &
      input_nparticle_max, input_nx, 2, 2)
  end if
end do

call VecDuplicate(particle_x(1), particle_electric, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(particle_x(1), particle_tmp1, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(particle_x(1), particle_tmp2, global_ierr)
CHKERRQ(global_ierr)

call VecGetOwnershipRange(particle_x(1), particle_ip_low, particle_ip_high, global_ierr)
CHKERRQ(global_ierr)
if (input_iptclshape == 3) then
  allocate (particle_shape_x_indexes( &
    input_nspecies, particle_ip_high - particle_ip_low &
  ), particle_shape_x_values( &
    input_nspecies, particle_ip_high - particle_ip_low &
  ))
end if

end subroutine particle_init


!!!!!!!!!!!!!!!!!!
! load particles !
!!!!!!!!!!!!!!!!!!
subroutine particle_load
use wtimer
use gaussian
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, ip, imode
PetscInt :: nparticle_unload
PetscScalar, dimension(:), pointer :: px, pv, pp, pw

call wtimer_start(global_iwt_particle_load)

call gaussian_init(input_seed_type, global_mype)

do ispecies = 1, input_nspecies
  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  if (input_imarker == 1) then ! load markers same as physical distribution
    ! currently only supports input_iptcldist == 0: (shifted) Maxwellian
    call gaussian_generate(pv)
    pv(:) = pv(:) * sqrt(input_species_temperature(ispecies) &
      / input_species_mass(ispecies)) + input_species_v0(ispecies)
    pp(:) = input_species_density(ispecies) * input_lx &
    / input_species_nparticle_init(ispecies)
  else ! input_imarker == 2, uniform in velocity space
    call random_number(pv)
    pv(:) = (pv(:) - 0.5_kpr) * 2.0_kpr * input_v_max
    if (input_iptcldist == 1) then ! two-stream1
      pp(:) = input_species_density(ispecies) * input_lx &
        * 2.0_kpr * input_v_max &
        / input_species_nparticle_init(ispecies) &
        * pv(:)**2 * exp(-pv(:)**2 / 2.0_kpr) / sqrt(2.0_kpr * PETSC_PI)
    elseif (input_iptcldist == 2) then ! two-stream2
      pp(:) = input_species_density(ispecies) * input_lx &
        * 2.0_kpr * input_v_max &
        / input_species_nparticle_init(ispecies) &
        * (exp(-(pv(:) + input_species_v0(ispecies))**2 / (2.0_kpr &
        * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
        + exp(-(pv(:) - input_species_v0(ispecies))**2 / (2.0_kpr &
        * input_species_temperature(ispecies) / input_species_mass(ispecies)))) &
        / sqrt(8.0_kpr * PETSC_PI &
        * input_species_temperature(ispecies) / input_species_mass(ispecies))
    elseif (input_iptcldist == 3) then ! bump-on-tail
      pp(:) = 1.0_kpr * input_lx &
        * 2.0_kpr * input_v_max &
        / input_species_nparticle_init(ispecies) &
        * (input_species_density(ispecies) * exp(-pv(:)**2 / (2.0_kpr &
        * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
        / sqrt(2.0_kpr * PETSC_PI &
        * input_species_temperature(ispecies) / input_species_mass(ispecies)) &
        + (1.0_kpr - input_species_density(ispecies)) &
        * exp(-(pv(:) - input_species_v0(ispecies))**2 / (2.0_kpr &
        * input_species_temperature2(ispecies) / input_species_mass(ispecies))) &
        / sqrt(2.0_kpr * PETSC_PI &
        * input_species_temperature2(ispecies) / input_species_mass(ispecies)))
    else ! (shifted) Maxwellian
      pp(:) = input_species_density(ispecies) * input_lx &
        * 2.0_kpr * input_v_max &
        / input_species_nparticle_init(ispecies) &
        * exp(-(pv(:) - input_species_v0(ispecies))**2 / (2.0_kpr &
        * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
        / sqrt(2.0_kpr * PETSC_PI &
        * input_species_temperature(ispecies) / input_species_mass(ispecies))
    end if
  end if

  ! uniform in x
  call random_number(px)
  px(:) = px(:) * input_lx

  pw(:) = 0.0_kpr
  do imode = 0, input_init_nmode - 1
    pw(:) = pw(:) &
      + input_init_mode_cos(imode) * cos(2.0_kpr * PETSC_PI / input_lx &
        * real(input_init_mode(imode), kpr) * px(:)) &
      + input_init_mode_sin(imode) * sin(2.0_kpr * PETSC_PI / input_lx &
        * real(input_init_mode(imode), kpr) * px(:))
  end do
  ! apply perturbation shape in velocity space
  do ip = 1, particle_ip_high - particle_ip_low
    pw(ip) = pw(ip) * pp(ip) &
      * input_pertb_shape(pv(ip), ispecies)
  end do

  ! unload not used particles
  nparticle_unload &
    = (input_nparticle_max - input_species_nparticle_init(ispecies)) &
    / global_npe
  if (global_mype == 0) then
    nparticle_unload = nparticle_unload &
      + mod(input_nparticle_max - input_species_nparticle_init(ispecies), &
      global_npe)
  end if
  particle_np(ispecies) = particle_ip_high - particle_ip_low - nparticle_unload

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  ! for linear, p = f_0 / g; for nonlinear, p = f / g = f_0 / g + delta f / g
  if (input_linear == 0) then
    call VecAXPY(particle_p(ispecies), &
      1.0_kpr, particle_w(ispecies), global_ierr)
    CHKERRQ(global_ierr)
  end if
end do

call wtimer_stop(global_iwt_particle_load)

end subroutine particle_load


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute shape matrix in x !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_shape_x
use wtimer
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip
PetscInt :: ix1, ix2
PetscScalar :: sx
PetscScalar, dimension(:), pointer :: px
PetscInt, dimension(0 : 1) :: indexes
PetscScalar, dimension(0 : 1) :: values

call wtimer_start(global_iwt_particle_shape)

do ispecies = 1, input_nspecies
  if (input_iptclshape == 1) then
    call MatDestroy(particle_shape_x(ispecies), global_ierr)
    CHKERRQ(global_ierr)
    call global_matcreate(particle_shape_x(ispecies), &
      input_nparticle_max, input_nx, 2, 2)
  elseif (input_iptclshape == 2) then
    call MatZeroEntries(particle_shape_x(ispecies), global_ierr)
    CHKERRQ(global_ierr)
  end if

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)

  do ip = 1, particle_np(ispecies)
    ! enforce periodic boundary condition
    px(ip) = mod(px(ip), input_lx)
    ! if x is negative, mod gives negative result, so shift it to positive
    if (px(ip) < 0.0_kpr) px(ip) = px(ip) + input_lx

    sx = px(ip) / input_lx * input_nx
    ix1 = floor(sx)
    sx = sx - real(ix1, kpr)
    if (input_iptclshape <= 2) then
      ix2 = ix1 + 1
      if (ix2 == input_nx) ix2 = 0

      nindex = 2
      indexes(0) = ix1
      indexes(1) = ix2
      values(0) = (1.0_kpr - sx)
      values(1) = sx
      call MatSetValues( &
        particle_shape_x(ispecies), 1, ip + particle_ip_low - 1, &
        nindex, indexes, values, INSERT_VALUES, global_ierr &
      )
      CHKERRQ(global_ierr)
    else ! implying input_iptclshape == 3
      particle_shape_x_indexes(ispecies, ip) = ix1
      particle_shape_x_values(ispecies, ip) = (1.0_kpr - sx)
    end if
  end do ! ip = 1, particle_np(ispecies)
  if (input_iptclshape <= 2) then
    call MatAssemblyBegin(particle_shape_x(ispecies), &
      MAT_FINAL_ASSEMBLY, global_ierr)
    CHKERRQ(global_ierr)
    call MatAssemblyEnd(particle_shape_x(ispecies), &
      MAT_FINAL_ASSEMBLY, global_ierr)
    CHKERRQ(global_ierr)
  end if

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
end do ! ispecies = 1, input_nspecies

call wtimer_stop(global_iwt_particle_shape)

end subroutine particle_compute_shape_x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute absolute value of perturbed distribution in v !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_dist_pertb_abs_v
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscScalar, dimension(:), pointer :: pv, pw

PetscScalar, dimension(0 : input_nv - 1) :: dist_pertb_abs_v

PetscInt :: ispecies, ip, iv
PetscScalar :: sv, df

do ispecies = 1, input_nspecies
  ! collect particles to v grids
  ! (v is equilibrium constant of motion)
  dist_pertb_abs_v(:) = 0.0_kpr

  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  do ip = 1, particle_np(ispecies)
    ! ignore too fast particle
    if (abs(pv(ip)) >= input_v_max) cycle

    sv = (pv(ip) + input_v_max) &
      / (input_v_max * 2.0_kpr) * (input_nv - 1)
    iv = floor(sv)
    sv = 1.0_kpr - (sv - real(iv, kpr))

    dist_pertb_abs_v(iv) = dist_pertb_abs_v(iv) + sv * abs(pw(ip))
    dist_pertb_abs_v(iv + 1) = dist_pertb_abs_v(iv + 1) &
      + (1.0_kpr - sv) * abs(pw(ip))
  end do ! ip = 1, particle_np(ispecies)

  call MPI_Allreduce(dist_pertb_abs_v(:), &
    particle_dist_pertb_abs_v(ispecies, :), &
    input_nv, MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)
end do ! ispecies = 1, input_nspecies

end subroutine particle_compute_dist_pertb_abs_v


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! merge not important particles               !
! using particle_dist_pertb_abs_v computed by !
! particle_compute_dist_pertb_abs_v           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_merge(thsh_frac_dist_pertb_abs_v)
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

! merging threshold in terms of fraction of absolute value of
! distribtuion in v
PetscReal, intent(in) :: thsh_frac_dist_pertb_abs_v

PetscInt, dimension(:, :, :), allocatable :: ipbin_top
PetscInt, dimension(:, :, :, :), allocatable :: ipbin
PetscInt :: ispecies, ip, ip1, ix, iv, iw

PetscScalar, dimension(:), pointer :: px, pv, pp, pw
PetscScalar :: sx, sv, df, df_thsh

allocate (ipbin(0 : input_nx - 1, 0 : input_nv - 1, 2, 1))
allocate (ipbin_top(0 : input_nx - 1, 0 : input_nv - 1, 2))

do ispecies = 1, input_nspecies
  ipbin_top(:, :, :) = 1
  df_thsh = maxval(particle_dist_pertb_abs_v(ispecies, :)) &
    * thsh_frac_dist_pertb_abs_v

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  ip = 0
  do
    ip = ip + 1
    if (ip > particle_np(ispecies)) exit

    ! ignore too fast particle
    !if (abs(pv(ip)) >= input_v_max) cycle

    sv = (pv(ip) + input_v_max) / (input_v_max * 2.0_kpr) * (input_nv - 1)
    iv = floor(sv)
    if (iv < 0) then
      iv = 0
      sv = 1.0_kpr
      df = particle_dist_pertb_abs_v(ispecies, iv)
    elseif (iv >= input_nv - 1) then
      iv = input_nv - 1
      sv = 1.0_kpr
      df = particle_dist_pertb_abs_v(ispecies, iv)
    else
      sv = 1.0_kpr - (sv - real(iv, kpr))
      df = particle_dist_pertb_abs_v(ispecies, iv) * sv &
        + particle_dist_pertb_abs_v(ispecies, iv + 1) * (1.0_kpr - sv)
    end if ! else of if (iv < 0) elseif (iv > input_nv - 1)
    ! ignore important particle
    if (df >= df_thsh) cycle

    ! enforce periodic boundary condition
    px(ip) = mod(px(ip), input_lx)
    ! if x is negative, mod gives negative result, so shift it to positive
    if (px(ip) < 0.0_kpr) px(ip) = px(ip) + input_lx

    sx = px(ip) / input_lx * input_nx
    ix = floor(sx)
    if (pw(ip) > 0.0_kpr) then
      iw = 2
    else
      iw = 1
    end if
    if (ipbin_top(ix, iv, iw) < 2) then ! bin not full yet
      ipbin(ix, iv, iw, ipbin_top(ix, iv, iw)) = ip
      ipbin_top(ix, iv, iw) = ipbin_top(ix, iv, iw) + 1
    else ! bin full, merge particles
      ! calculate merged particle and put in slot of index ip1
      ip1 = ipbin(ix, iv, iw, 1)
      px(ip1) = (pw(ip1) * px(ip1) + pw(ip) * px(ip)) / (pw(ip1) + pw(ip))
      ! no need to enforce boundary condition as it will be enforced
      ! in particle_compute_shape_x or interaction_collect_charge
      pv(ip1) = (pw(ip1) * pv(ip1) + pw(ip) * pv(ip)) / (pw(ip1) + pw(ip))
      pp(ip1) = pp(ip1) + pp(ip)
      pw(ip1) = pw(ip1) + pw(ip)

      if (ip < particle_np(ispecies)) then
        ! move the last particle to current index
        px(ip) = px(particle_np(ispecies))
        pv(ip) = pv(particle_np(ispecies))
        pp(ip) = pp(particle_np(ispecies))
        pw(ip) = pw(particle_np(ispecies))
        ip = ip - 1 ! make the loop to go over the current index again
      end if
      particle_np(ispecies) = particle_np(ispecies) - 1

      ! reset bin
      ipbin_top(ix, iv, iw) = 1
    end if
  end do ! ip = 1, particle_np(ispecies)

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)
end do ! ispecies = 1, input_nspecies

deallocate (ipbin, ipbin_top)

end subroutine particle_merge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! throw away some not important particles     !
! using particle_dist_pertb_abs_v computed by !
! particle_compute_dist_pertb_abs_v           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_throwaway(thsh_frac_dist_pertb_abs_v)
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

! throwing away threshold in terms of fraction of absolute value of
! distribtuion in v
PetscReal, intent(in) :: thsh_frac_dist_pertb_abs_v

PetscInt :: ispecies, ip, iv

PetscScalar, dimension(:), pointer :: px, pv, pp, pw
PetscScalar :: sv, df, df_thsh
PetscReal :: dice

do ispecies = 1, input_nspecies
  df_thsh = maxval(particle_dist_pertb_abs_v(ispecies, :)) &
    * thsh_frac_dist_pertb_abs_v

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  ip = 0
  do
    ip = ip + 1
    if (ip > particle_np(ispecies)) exit

    ! ignore too fast particle
    !if (abs(pv(ip)) >= input_v_max) cycle

    sv = (pv(ip) + input_v_max) / (input_v_max * 2.0_kpr) * (input_nv - 1)
    iv = floor(sv)
    if (iv < 0) then
      iv = 0
      sv = 1.0_kpr
      df = particle_dist_pertb_abs_v(ispecies, iv)
    elseif (iv >= input_nv - 1) then
      iv = input_nv - 1
      sv = 1.0_kpr
      df = particle_dist_pertb_abs_v(ispecies, iv)
    else
      sv = 1.0_kpr - (sv - real(iv, kpr))
      df = particle_dist_pertb_abs_v(ispecies, iv) * sv &
        + particle_dist_pertb_abs_v(ispecies, iv + 1) * (1.0_kpr - sv)
    end if ! else of if (iv < 0) elseif (iv > input_nv - 1)
    if (input_typethrowaway == 1) then
      ! ignore important particle
      if (df >= df_thsh) cycle
    end if

    df = df / maxval(particle_dist_pertb_abs_v(ispecies, :))
    call random_number(dice)

    if ((input_typethrowaway == 1 .and. dice < input_throwaway_frac) &
      .or. (input_typethrowaway == 2 .and. dice > df)) then
      ! throw away particle, and move the last particle to current index
      if (ip < particle_np(ispecies)) then
        px(ip) = px(particle_np(ispecies))
        pv(ip) = pv(particle_np(ispecies))
        pp(ip) = pp(particle_np(ispecies))
        pw(ip) = pw(particle_np(ispecies))
        ip = ip - 1 ! make the loop to go over current index again
      end if
      particle_np(ispecies) = particle_np(ispecies) - 1
    else
      ! keep particle, but scale up weight
      if (input_typethrowaway == 1) then
        pp(ip) = pp(ip) / (1.0_kpr - input_throwaway_frac)
        pw(ip) = pw(ip) / (1.0_kpr - input_throwaway_frac)
      else
        pp(ip) = pp(ip) / df
        pw(ip) = pw(ip) / df
      end if
    end if
  end do ! ip = 1, particle_np(ispecies)

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)
end do ! ispecies = 1, input_nspecies

end subroutine particle_throwaway


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! split resonant particles                    !
! using particle_dist_pertb_abs_v computed by !
! particle_compute_dist_pertb_abs_v           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_split(thsh_frac_dist_pertb_abs_v)
use pic1dp_global
use gaussian
implicit none
#include "finclude/petsc.h90"

! splitting threshold in terms of fraction of absolute value of
! distribtuion in v
PetscReal, intent(in) :: thsh_frac_dist_pertb_abs_v

PetscInt :: ispecies, ip, ip1, iv, igroup, np_inc
PetscScalar, dimension(:), pointer :: px, pv, pp, pw
PetscReal, dimension(input_split_ngroup) :: grand
PetscScalar :: sv
PetscScalar :: df, df_thsh
!PetscBool :: bout

!bout = PETSC_TRUE

do ispecies = 1, input_nspecies
  ! if particle array for this species is full, no splitting can be done
  if (particle_ip_high - particle_ip_low - particle_np(ispecies) &
    < 2 * input_split_ngroup - 1) cycle

  np_inc = 0
  df_thsh = maxval(particle_dist_pertb_abs_v(ispecies, :)) &
    * thsh_frac_dist_pertb_abs_v

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  do ip = 1, particle_np(ispecies)
    ! if particle array for this species is full, no more splitting can be done
    if (particle_ip_high - particle_ip_low - (particle_np(ispecies) + np_inc) &
      < 2 * input_split_ngroup - 1) exit
    ! ignore too fast particle
    !if (abs(pv(ip)) >= input_v_max) cycle

    sv = (pv(ip) + input_v_max) / (input_v_max * 2.0_kpr) * (input_nv - 1)
    iv = floor(sv)
    if (iv < 0) then
      iv = 0
      sv = 1.0_kpr
      df = particle_dist_pertb_abs_v(ispecies, iv)
    elseif (iv >= input_nv - 1) then
      iv = input_nv - 1
      sv = 1.0_kpr
      df = particle_dist_pertb_abs_v(ispecies, iv)
    else
      sv = 1.0_kpr - (sv - real(iv, kpr))
      df = particle_dist_pertb_abs_v(ispecies, iv) * sv &
        + particle_dist_pertb_abs_v(ispecies, iv + 1) * (1.0_kpr - sv)
    end if ! else of if (iv < 0) elseif (iv > nv - 1)
    ! ignore important particle
    if (df <= df_thsh) cycle

    call gaussian_generate(grand)
    !grand(:) = grand(:) * input_lx / input_nx * dx_sig_frac
    grand(:) = grand(:) * 2.0_kpr * input_v_max / input_nv &
      * input_split_dv_sig_frac

!    if (bout) then
!      write (*, *) 'grand:', grand
!      write (*, *) 'x, v, p, w:', px(ip), pv(ip), pp(ip), pw(ip)
!    end if
    do igroup = 1, input_split_ngroup
      ip1 = particle_np(ispecies) + np_inc + igroup * 2 - 1
      px(ip1) = px(ip)
      pv(ip1) = pv(ip) + grand(igroup)
      pp(ip1) = pp(ip) / (input_split_ngroup * 2.0_kpr)
      if (input_deltaf == 1) pw(ip1) = pw(ip) / (input_split_ngroup * 2.0_kpr)
!      if (bout) then
!        write (*, *) 'x, v, p, w:', px(ip1), pv(ip1), pp(ip1), pw(ip1)
!      end if

      if (igroup == input_split_ngroup) then
        ip1 = ip
      else
        ip1 = particle_np(ispecies) + np_inc + igroup * 2
      end if
      px(ip1) = px(ip)
      pv(ip1) = pv(ip) - grand(igroup)
      pp(ip1) = pp(ip) / (input_split_ngroup * 2.0_kpr)
      if (input_deltaf == 1) pw(ip1) = pw(ip) / (input_split_ngroup * 2.0_kpr)
!      if (bout) then
!        write (*, *) 'x, v, p, w:', px(ip1), pv(ip1), pp(ip1), pw(ip1)
!      end if
    end do

    np_inc = np_inc + (2 * input_split_ngroup - 1)
!    bout = PETSC_FALSE
  end do ! ip = 1, particle_np(ispecies)

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  particle_np(ispecies) = particle_np(ispecies) + np_inc
end do ! ispecies = 1, input_nspecies

end subroutine particle_split


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! particle optimization, manage calling of merge/throwaway/split !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_optimize(flag_optimized)
use wtimer
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

logical, intent(out) :: flag_optimized ! whether optimization is performed

flag_optimized = .false.

if (input_deltaf == 0) return ! now only support optimization for delta f

call wtimer_start(global_iwt_particle_optimize)

! merge not important particles
if (particle_imerge > 0 .and. particle_imerge <= input_nmerge) then
  if ( &
    global_time + input_dt >= input_tmerge(particle_imerge) &
    .and. global_irk == 2 &
  ) then
    ! calculate absolute value of perturbed distribution in v
    call particle_compute_dist_pertb_abs_v
    ! perform merging
    call particle_merge(input_thshmerge(particle_imerge))
    particle_imerge = particle_imerge + 1
    flag_optimized = .true.
  end if
end if

! throw away some not important particles
if (particle_ithrowaway > 0 .and. particle_ithrowaway <= input_nthrowaway) then
  if ( &
    global_time + input_dt >= input_tthrowaway(particle_ithrowaway) &
    .and. global_irk == 2 &
  ) then
    ! calculate absolute value of perturbed distribution in v
    call particle_compute_dist_pertb_abs_v
    ! perform throwing away
    call particle_throwaway(input_thshthrowaway(particle_ithrowaway))
    particle_ithrowaway = particle_ithrowaway + 1
    flag_optimized = .true.
  end if
end if

! split important particles
if (particle_isplit > 0 .and. particle_isplit <= input_nsplit) then
  if ( &
    global_time + input_dt >= input_tsplit(particle_isplit) &
    .and. global_irk == 2 &
  ) then
    ! calculate absolute value of perturbed distribution in v
    call particle_compute_dist_pertb_abs_v
    ! perform splitting
    call particle_split(input_thshsplit(particle_isplit))
    particle_isplit = particle_isplit + 1
    flag_optimized = .true.
  end if
end if

call wtimer_stop(global_iwt_particle_optimize)

end subroutine particle_optimize


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_final
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

do ispecies = 1, input_nspecies
  call VecDestroy(particle_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_p(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_w(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call VecDestroy(particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_v_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_w_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call MatDestroy(particle_shape_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do

call VecDestroy(particle_electric, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(particle_tmp1, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(particle_tmp2, global_ierr)
CHKERRQ(global_ierr)

if (input_iptclshape == 3) then
  deallocate (particle_shape_x_indexes, particle_shape_x_values)
end if

end subroutine particle_final

end module pic1dp_particle

