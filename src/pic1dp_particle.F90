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
Vec, dimension(input_nspecies) :: &
  particle_x, particle_v, particle_p, particle_w, &
  particle_x_bak, particle_v_bak, particle_w_bak

! temporary particle vectors, need to be used 
! in pic1dp_interaction and pic1dp_output
Vec :: particle_tmp1, particle_tmp2

Mat, dimension(input_nspecies) :: &
  particle_shape_xv, particle_shape_x, particle_shape_v

! shape arrays
PetscInt, dimension(:, :), allocatable :: particle_shape_x_indexes
PetscReal, dimension(:, :), allocatable :: particle_shape_x_values

! local index range of particle
PetscInt :: particle_ip_low, particle_ip_high

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

  call VecDuplicate(particle_x(1), particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_v_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_w_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call global_matcreate(particle_shape_xv(ispecies), &
    input_nparticle, input_nx * input_output_nv, 4, 4)
  call global_matcreate(particle_shape_x(ispecies), &
    input_nparticle, input_nx, 2, 2)
  call global_matcreate(particle_shape_v(ispecies), &
    input_nparticle, input_output_nv, 2, 2)
end do

call VecDuplicate(particle_x(1), particle_tmp1, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(particle_x(1), particle_tmp2, global_ierr)
CHKERRQ(global_ierr)

call VecGetOwnershipRange(particle_x(1), particle_ip_low, particle_ip_high, global_ierr)
CHKERRQ(global_ierr)
allocate (particle_shape_x_indexes( &
  input_nspecies, particle_ip_low : particle_ip_high - 1 &
), particle_shape_x_values( &
  input_nspecies, particle_ip_low : particle_ip_high - 1 &
))

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

  if (input_imarker == 1) then ! Maxwellian in velocity space
    call gaussian_generate(pv)
    pv(:) = pv(:) * sqrt( &
      input_species_temperature(ispecies) / input_species_mass(ispecies))
    pp(:) = input_lx / input_nparticle
  else ! input_imarker == 2, uniform in velocity space
    call random_number(pv)
    pv(:) = (pv(:) - 0.5_kpr) * 2.0_kpr * input_v_max
    pp(:) = input_lx / input_nparticle * exp(-pv(:)**2 / (2.0_kpr &
      * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
      / sqrt(2.0_kpr * PETSC_PI &
      * input_species_temperature(ispecies) / input_species_mass(ispecies)) &
      * 2.0_kpr * input_v_max
  end if

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
  do ip = particle_ip_low, particle_ip_high - 1
    pw(ip - particle_ip_low + 1) = pw(ip - particle_ip_low + 1) &
      * input_lx / input_nparticle &
      * input_pertb_shape(pv(ip - particle_ip_low + 1))
  end do

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  ! for linear, p = f_0 / g; for nonlinear, p = f / g
  if (input_linear == 0) then
    call VecAXPY(particle_p(ispecies), &
      1.0_kpr, particle_w(ispecies), global_ierr)
    CHKERRQ(global_ierr)
  end if
end do

!  call VecView(particle_x(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_p(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_w(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)

end subroutine particle_load


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute shape matrix in x-v plane !
! obsolete                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_shape_xv
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip
PetscInt :: ix1, ix2, iv, ixv
PetscScalar :: sx, sv
PetscScalar, dimension(:), pointer :: px, pv
PetscInt, dimension(0 : 3) :: indexes
PetscScalar, dimension(0 : 3) :: values


do ispecies = 1, input_nspecies
!  call MatDestroy(particle_shape_xv(ispecies), global_ierr)
!  CHKERRQ(global_ierr)
  call MatZeroEntries(particle_shape_xv(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)

  do ip = particle_ip_low, particle_ip_high - 1
    if (abs(pv(ip - particle_ip_low + 1)) >= input_v_max) cycle

    sx = px(ip - particle_ip_low + 1) / input_lx * input_nx
    ix1 = floor(sx)
    sx = sx - real(ix1, kpr)
    ix2 = ix1 + 1
    if (ix2 == input_nx) ix2 = 0

    sv = (pv(ip - particle_ip_low + 1) + input_v_max) &
      / (input_v_max * 2.0_kpr) * (input_output_nv - 1)
    iv = floor(sv)
    sv = sv - real(iv, kpr)

    nindex = 4
    indexes(0) = iv * input_nx + ix1
    indexes(1) = iv * input_nx + ix2
    indexes(2) = (iv + 1) * input_nx + ix1
    indexes(3) = (iv + 1) * input_nx + ix2
    values(0) = (1.0_kpr - sx) * (1.0_kpr - sv)
    values(1) = sx * (1.0_kpr - sv)
    values(2) = (1.0_kpr - sx) * sv
    values(3) = sx * sv
    call MatSetValues( &
      particle_shape_xv(ispecies), 1, ip, &
      nindex, indexes, values, INSERT_VALUES, global_ierr &
    )
    CHKERRQ(global_ierr)
  end do
  call MatAssemblyBegin(particle_shape_xv(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)
  call MatAssemblyEnd(particle_shape_xv(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)


!  call VecView(particle_x(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call MatView(particle_shape_xv(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
end do

end subroutine particle_compute_shape_xv


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
PetscScalar, dimension(:), pointer :: px
PetscInt, dimension(0 : 1) :: indexes
PetscScalar, dimension(0 : 1) :: values

!call VecGetOwnershipRange(particle_x(1), ip_low, ip_high, global_ierr)
!CHKERRQ(global_ierr)

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

  do ip = particle_ip_low, particle_ip_high - 1
    ! enforce periodic boundary condition
    px(ip - particle_ip_low + 1) = mod(px(ip - particle_ip_low + 1), input_lx)
    ! if x is negative, mod gives negative result, so shift it to positive
    if (px(ip - particle_ip_low + 1) < 0.0_kpr) &
      px(ip - particle_ip_low + 1) = px(ip - particle_ip_low + 1) + input_lx

    sx = px(ip - particle_ip_low + 1) / input_lx * input_nx
    ix1 = floor(sx)
    sx = sx - real(ix1, kpr)
!    if (global_mype == 0 .and. ip == particle_ip_low) write (*, *) ix1, 1.0_kpr - sx
    if (input_iptclshape < 3) then
      ix2 = ix1 + 1
      if (ix2 == input_nx) ix2 = 0

      nindex = 2
      indexes(0) = ix1
      indexes(1) = ix2
      values(0) = (1.0_kpr - sx)
      values(1) = sx
      call MatSetValues( &
        particle_shape_x(ispecies), 1, ip, &
        nindex, indexes, values, INSERT_VALUES, global_ierr &
      )
      CHKERRQ(global_ierr)
    else ! input_iptclshape == 3
      particle_shape_x_indexes(ispecies, ip) = ix1
      particle_shape_x_values(ispecies, ip) = (1.0_kpr - sx)
    end if
  end do
  if (input_iptclshape < 3) then
    call MatAssemblyBegin(particle_shape_x(ispecies), &
      MAT_FINAL_ASSEMBLY, global_ierr)
    CHKERRQ(global_ierr)
    call MatAssemblyEnd(particle_shape_x(ispecies), &
      MAT_FINAL_ASSEMBLY, global_ierr)
    CHKERRQ(global_ierr)
  end if

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)

!  call VecView(particle_x(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call MatView(particle_shape_xv(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
end do

end subroutine particle_compute_shape_x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute shape matrix in v !
! obsolete                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_shape_v
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip
PetscInt :: iv
PetscScalar :: sv
PetscScalar, dimension(:), pointer :: pv
PetscInt, dimension(0 : 1) :: indexes
PetscScalar, dimension(0 : 1) :: values

!call VecGetOwnershipRange(particle_v(1), ip_low, ip_high, global_ierr)
!CHKERRQ(global_ierr)

do ispecies = 1, input_nspecies
!  call MatDestroy(particle_shape_v(ispecies), global_ierr)
!  CHKERRQ(global_ierr)
!  call global_matcreate(particle_shape_v(ispecies), &
!    input_nparticle, input_nx, 2, 2)
  call MatZeroEntries(particle_shape_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)

  do ip = particle_ip_low, particle_ip_high - 1
    ! if particle velocity out of v_max, ignore this particle
    if (abs(pv(ip - particle_ip_low + 1)) >= input_v_max) cycle

    sv = (pv(ip - particle_ip_low + 1) + input_v_max) &
      / (input_v_max * 2.0_kpr) * (input_output_nv - 1)
    iv = floor(sv)
    sv = sv - real(iv, kpr)

    nindex = 2
    indexes(0) = iv
    indexes(1) = iv + 1
    values(0) = (1.0_kpr - sv)
    values(1) = sv
    call MatSetValues( &
      particle_shape_v(ispecies), 1, ip, &
      nindex, indexes, values, INSERT_VALUES, global_ierr &
    )
    CHKERRQ(global_ierr)
  end do
  call MatAssemblyBegin(particle_shape_v(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)
  call MatAssemblyEnd(particle_shape_v(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)

!  call VecView(particle_x(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call MatView(particle_shape_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
end do

end subroutine particle_compute_shape_v


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

  call MatDestroy(particle_shape_xv(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call MatDestroy(particle_shape_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call MatDestroy(particle_shape_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do

call VecDestroy(particle_tmp1, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(particle_tmp2, global_ierr)
CHKERRQ(global_ierr)

deallocate (particle_shape_x_indexes, particle_shape_x_values)

end subroutine particle_final

end module pic1dp_particle

