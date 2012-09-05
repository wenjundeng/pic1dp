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

! particle x coordinate, velocity
! equilibrium weight: p = f / g (nonlinear); p = f_0 / g (linear)
! perturbed weight: w = delta f / g
! f is total distribution, delta f is perturbed distribution
! g is marker distribution
! s is particle status: 0: normal/non-resonant; -1: invalid; 1: resonant
Vec, dimension(input_nspecies) :: &
  particle_x, particle_v, particle_p, particle_w, &
  particle_x_bak, particle_v_bak, particle_w_bak, &
  particle_s

! electric field at particle position, need to be used in pic1dp_interaction
Vec :: particle_electric

! temporary particle vectors, need to be used 
! in pic1dp_interaction and pic1dp_output
Vec :: particle_tmp1, particle_tmp2

Mat, dimension(input_nspecies) :: particle_shape_x

! shape arrays
PetscInt, dimension(:, :), allocatable :: particle_shape_x_indexes
PetscReal, dimension(:, :), allocatable :: particle_shape_x_values

! local index range of particle
PetscInt :: particle_ip_low, particle_ip_high

! delta f threshold (fraction of max(abs(delta f))) for resonant particle
PetscReal :: particle_df_thsh_res_frac

! # of reduced particles due to merge
PetscInt :: particle_nredu(input_nspecies) = 0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_init
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

call VecCreate(MPI_COMM_WORLD, particle_x(1), global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(particle_x(1), PETSC_DECIDE, input_nparticle, global_ierr)
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
  call VecDuplicate(particle_x(1), particle_s(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecSet(particle_s(ispecies), 0.0_kpr, global_ierr)
  CHKERRQ(global_ierr)

  call VecDuplicate(particle_x(1), particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_v_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_w_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  if (input_iptclshape <= 2) then
    call global_matcreate(particle_shape_x(ispecies), &
      input_nparticle, input_nx, 2, 2)
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
use pic1dp_global
use pic1dp_input
use gaussian
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

PetscInt :: ip, imode
PetscScalar, dimension(:), pointer :: px, pv, pp, pw

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
    ! currently only supports input_iptcldist == 1: (shifted) Maxwellian
    call gaussian_generate(pv)
    pv(:) = pv(:) * sqrt(input_species_temperature(ispecies) &
      / input_species_mass(ispecies)) + input_species_v0(ispecies)
    pp(:) = input_species_density(ispecies) * input_lx / input_nparticle
  else ! input_imarker == 2, uniform in velocity space
    call random_number(pv)
    pv(:) = (pv(:) - 0.5_kpr) * 2.0_kpr * input_v_max
    if (input_iptcldist == 1) then ! two-stream1
      pp(:) = input_species_density(ispecies) * input_lx / input_nparticle &
        * pv(:)**2 * exp(-pv(:)**2 / 2.0_kpr) / sqrt(2.0_kpr * PETSC_PI) &
        * 2.0_kpr * input_v_max
    elseif (input_iptcldist == 2) then ! two-stream2
      pp(:) =  input_species_density(ispecies) * input_lx / input_nparticle &
        * (exp(-(pv(:) + input_species_v0(ispecies))**2 / (2.0_kpr &
        * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
        + exp(-(pv(:) - input_species_v0(ispecies))**2 / (2.0_kpr &
        * input_species_temperature(ispecies) / input_species_mass(ispecies)))) &
        / sqrt(8.0_kpr * PETSC_PI &
        * input_species_temperature(ispecies) / input_species_mass(ispecies)) &
        * 2.0_kpr * input_v_max
    else ! (shifted) Maxwellian
      pp(:) = input_species_density(ispecies) * input_lx / input_nparticle &
        * exp(-(pv(:) - input_species_v0(ispecies))**2 / (2.0_kpr &
        * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
        / sqrt(2.0_kpr * PETSC_PI &
        * input_species_temperature(ispecies) / input_species_mass(ispecies)) &
        * 2.0_kpr * input_v_max
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
    pw(ip) = pw(ip) &
      * input_species_density(ispecies) * input_lx / input_nparticle &
      * input_pertb_shape(pv(ip), ispecies)
  end do

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

end subroutine particle_load


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute shape matrix in x !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_shape_x
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip
PetscInt :: ix1, ix2
PetscScalar :: sx
PetscScalar, dimension(:), pointer :: px, ps
PetscInt, dimension(0 : 1) :: indexes
PetscScalar, dimension(0 : 1) :: values

do ispecies = 1, input_nspecies
  if (input_iptclshape == 1) then
    call MatDestroy(particle_shape_x(ispecies), global_ierr)
    CHKERRQ(global_ierr)
    call global_matcreate(particle_shape_x(ispecies), &
      input_nparticle, input_nx, 2, 2)
  elseif (input_iptclshape == 2) then
    call MatZeroEntries(particle_shape_x(ispecies), global_ierr)
    CHKERRQ(global_ierr)
  end if

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)

  do ip = 1, particle_ip_high - particle_ip_low
    ! ignore invalid particles
    if (ps(ip) < -0.5_kpr) cycle
    ! enforce periodic boundary condition
    px(ip) = mod(px(ip), input_lx)
    ! if x is negative, mod gives negative result, so shift it to positive
    if (px(ip) < 0.0_kpr) px(ip) = px(ip) + input_lx

    sx = px(ip) / input_lx * input_nx
    ix1 = floor(sx)
    sx = sx - real(ix1, kpr)
!    if (global_mype == 0 .and. ip == 1) write (*, *) ix1, 1.0_kpr - sx
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
    else ! input_iptclshape == 3
      particle_shape_x_indexes(ispecies, ip) = ix1
      particle_shape_x_values(ispecies, ip) = (1.0_kpr - sx)
    end if
  end do ! ip = 1, particle_ip_high - particle_ip_low
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
  call VecRestoreArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine particle_compute_shape_x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! detect resonant particles and set resonant status !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_detect_resonant
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscScalar, dimension(0 : input_nv - 1) :: &
  ptcldist_pertb_v, ptcldist_pertb_v_redu

PetscScalar, dimension(:), pointer :: pv, pw, ps

PetscInt :: ispecies, ip, iv
PetscScalar :: sv, df, df_thsh_res

do ispecies = 1, input_nspecies
  ! first step, collect particles to v grids
  ! (v is equilibrium constant of motion)
  ptcldist_pertb_v(:) = 0.0_kpr

  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)

  do ip = 1, particle_ip_high - particle_ip_low
    ! ignore invalid particle
    if (ps(ip) < -0.5_kpr) cycle
    ! ignore too fast particle
    if (abs(pv(ip)) >= input_v_max) cycle

    sv = (pv(ip) + input_v_max) &
      / (input_v_max * 2.0_kpr) * (input_nv - 1)
    iv = floor(sv)
    sv = 1.0_kpr - (sv - real(iv, kpr))

    ! note that abs() is applied to w, otherwise positive and negative
    ! resonant regions may get canceled out
    ptcldist_pertb_v(iv) = ptcldist_pertb_v(iv) + sv * abs(pw(ip))
    ptcldist_pertb_v(iv + 1) = ptcldist_pertb_v(iv + 1) &
      + (1.0_kpr - sv) * abs(pw(ip))
  end do ! ip = 1, particle_ip_high - particle_ip_low

  call MPI_Allreduce(ptcldist_pertb_v, ptcldist_pertb_v_redu, &
    input_nv, MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)

  df_thsh_res = maxval(ptcldist_pertb_v_redu) * particle_df_thsh_res_frac

  ! second step, mark resonant particles
  do ip = 1, particle_ip_high - particle_ip_low
    ! ignore invalid particle
    if (ps(ip) < -0.5_kpr) cycle
    ! ignore too fast particle
    if (abs(pv(ip)) >= input_v_max) cycle

    sv = (pv(ip) + input_v_max) &
      / (input_v_max * 2.0_kpr) * (input_nv - 1)
    iv = floor(sv)
    sv = 1.0_kpr - (sv - real(iv, kpr))

    df = ptcldist_pertb_v_redu(iv) * sv &
      + ptcldist_pertb_v_redu(iv + 1) * (1.0_kpr - sv)

    if (df > df_thsh_res) then
      ps(ip) = 1.0_kpr ! resonant
    else
      ps(ip) = 0.0_kpr ! non-resonant
    end if
  end do ! ip = 1, particle_ip_high - particle_ip_low

  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)

end do ! ispecies = 1, input_nspecies

end subroutine particle_detect_resonant


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! merge non-resonant particles                             !
! using resonant status marked by particle_detect_resonant !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_merge
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt, dimension(:, :, :), allocatable :: ipbin_top
PetscInt, dimension(:, :, :, :), allocatable :: ipbin
PetscInt :: ispecies, ip, ip1, ix, iv, iw, nx, nv
PetscInt :: nredu1, nredu2

PetscScalar, dimension(:), pointer :: px, pv, pp, pw, ps
PetscScalar :: sx, sv

nx = input_nx * 2
nv = input_nv * 2

allocate (ipbin(0 : nx - 1, 0 : nv - 1, 2, 1))
allocate (ipbin_top(0 : nx - 1, 0 : nv - 1, 2))

do ispecies = 1, input_nspecies
  ipbin_top(:, :, :) = 1
  nredu1 = 0

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)

  do ip = 1, particle_ip_high - particle_ip_low
    ! ignore invalid and resonant particle
    if (ps(ip) < -0.5_kpr .or. ps(ip) > 0.5_kpr) cycle
    ! ignore too fast particle
    if (abs(pv(ip)) >= input_v_max) cycle

    ! enforce periodic boundary condition
    px(ip) = mod(px(ip), input_lx)
    ! if x is negative, mod gives negative result, so shift it to positive
    if (px(ip) < 0.0_kpr) px(ip) = px(ip) + input_lx

    sx = px(ip) / input_lx * nx
    ix = floor(sx)
    sv = (pv(ip) + input_v_max) / (input_v_max * 2.0_kpr) * (nv - 1)
    iv = floor(sv)

    if (pw(ip) > 0.0_kpr) then
      iw = 2
    else
      iw = 1
    end if
    if (ipbin_top(ix, iv, iw) < 2) then ! bin not full yet
      ipbin(ix, iv, iw, ipbin_top(ix, iv, iw)) = ip
      ipbin_top(ix, iv, iw) = ipbin_top(ix, iv, iw) + 1
    else ! bin full, merge particles
      ! calculate merged particle
      ip1 = ipbin(ix, iv, iw, 1)
      px(ip1) = (pw(ip1) * px(ip1) + pw(ip) * px(ip)) / (pw(ip1) + pw(ip))
      ! no need to enforce boundary condition as it will be enforced
      ! in particle_compute_shape_x or interaction_collect_charge
      pv(ip1) = (pw(ip1) * pv(ip1) + pw(ip) * pv(ip)) / (pw(ip1) + pw(ip))
      pp(ip1) = pp(ip1) + pp(ip)
      pw(ip1) = pw(ip1) + pw(ip)

      ! set obsolete particle status to be invalid and reset other properties
      ps(ip) = -1.0_kpr
      pv(ip) = 0.0_kpr
      pp(ip) = 0.0_kpr
      pw(ip) = 0.0_kpr
      nredu1 = nredu1 + 1

      ! reset bin
      ipbin_top(ix, iv, iw) = 1
    end if
  end do ! ip = 1, particle_ip_high - particle_ip_low

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)

  ! combine # of reduced particles
  call MPI_Allreduce(nredu1, nredu2, 1, MPIU_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  particle_nredu(ispecies) = particle_nredu(ispecies) + nredu2
end do ! ispecies = 1, input_nspecies

deallocate (ipbin, ipbin_top)

end subroutine particle_merge


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_final
use pic1dp_global
use pic1dp_input
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
  call VecDestroy(particle_s(ispecies), global_ierr)
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

