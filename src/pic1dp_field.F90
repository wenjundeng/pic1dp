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


! module for managing field quantities
module pic1dp_field
!use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petscdef.h"

! electric field and charge density
Vec :: field_electric, field_electric_seq, field_chargeden
Vec :: field_tmp

! for collecting charge and calculate electric field in Fourier space
Vec :: field_mode_re, field_mode_im

! inverse grad operator (scaled by i) in Fourier space
Vec :: field_mode_grad_inv

! for partial DFT and inverse DFT
Mat :: field_fourier_re, field_fourier_im

IS :: field_is_electric, field_is_electric_seq
VecScatter :: field_vs_electric

PetscScalar, dimension(0 : input_nx - 1) :: field_arr_charge1, field_arr_charge2

! local range of index of x
PetscInt :: field_ix_low, field_ix_high

! local range of index of mode
PetscInt :: field_imode_low, field_imode_high

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_init
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: ix, imode !, ix_low, ix_high
PetscInt, dimension(1) :: arrix
PetscInt :: nindex
PetscInt, dimension(0 : input_nmode - 1) :: indexes
PetscScalar, dimension(0 : input_nmode - 1) :: values

call VecCreate(MPI_COMM_WORLD, field_electric, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(field_electric, PETSC_DECIDE, input_nx, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(field_electric, global_ierr)
CHKERRQ(global_ierr)

call VecDuplicate(field_electric, field_chargeden, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(field_electric, field_tmp, global_ierr)
CHKERRQ(global_ierr)

call VecCreate(MPI_COMM_SELF, field_electric_seq, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(field_electric_seq, PETSC_DECIDE, input_nx, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(field_electric_seq, global_ierr)
CHKERRQ(global_ierr)

call VecCreate(MPI_COMM_WORLD, field_mode_re, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(field_mode_re, PETSC_DECIDE, input_nmode, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(field_mode_re, global_ierr)
CHKERRQ(global_ierr)

call VecDuplicate(field_mode_re, field_mode_im, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(field_mode_re, field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)

!call MatCreate(MPI_COMM_WORLD, field_fourier_re, global_ierr)
!CHKERRQ(global_ierr)
!call MatSetType(field_fourier_re, MATDENSE, global_ierr)
!CHKERRQ(global_ierr)
!call MatSetSizes( &
!  field_fourier_re, PETSC_DECIDE, PETSC_DECIDE, &
!  input_nx, input_nmode, global_ierr &
!)
!CHKERRQ(global_ierr)
!call MatCreateDense( &
!  MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, &
!  input_nx, input_nmode, PETSC_NULL_SCALAR, &
!  field_fourier_re, global_ierr &
!)
!CHKERRQ(global_ierr)
!call MatSetUp(field_fourier_re, global_ierr)
!CHKERRQ(global_ierr)
!call MatSetFromOptions(field_fourier_re, global_ierr)
!CHKERRQ(global_ierr)

!call MatCreate(MPI_COMM_WORLD, field_fourier_im, global_ierr)
!CHKERRQ(global_ierr)
!call MatSetType(field_fourier_im, MATDENSE, global_ierr)
!CHKERRQ(global_ierr)
!call MatSetSizes( &
!  field_fourier_im, PETSC_DECIDE, PETSC_DECIDE, &
!  input_nx, input_nmode, global_ierr &
!)
!CHKERRQ(global_ierr)
!call MatCreateDense( &
!  MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, &
!  input_nx, input_nmode, PETSC_NULL_SCALAR, &
!  field_fourier_im, global_ierr &
!)
!CHKERRQ(global_ierr)
!call MatSetUp(field_fourier_re, global_ierr)
!CHKERRQ(global_ierr)
!call MatSetFromOptions(field_fourier_im, global_ierr)
!CHKERRQ(global_ierr)

! I failed to make the PETSc dense matrix work,
! so here I use sparse matrix to work around
call global_matcreate(field_fourier_re, input_nx, input_nmode, &
  input_nmode, input_nmode)
call global_matcreate(field_fourier_im, input_nx, input_nmode, &
  input_nmode, input_nmode)

call ISCreateStride(MPI_COMM_WORLD, input_nx, 0, 1, &
  field_is_electric, global_ierr)
CHKERRQ(global_ierr)
call ISCreateStride(MPI_COMM_SELF, input_nx, 0, 1, &
  field_is_electric_seq, global_ierr)
CHKERRQ(global_ierr)
call VecScatterCreate(field_electric, field_is_electric, &
  field_electric_seq, field_is_electric_seq, field_vs_electric, global_ierr)
CHKERRQ(global_ierr)

call VecGetOwnershipRange(field_electric, field_ix_low, field_ix_high, global_ierr)
CHKERRQ(global_ierr)

! initialize inverse grad operator in partial Fourier space
call VecGetOwnershipRange(field_mode_grad_inv, &
  field_imode_low, field_imode_high, global_ierr)
CHKERRQ(global_ierr)
nindex = field_imode_high - field_imode_low
do imode = field_imode_low, field_imode_high - 1
  indexes(imode - field_imode_low) = imode
  values(imode - field_imode_low) &
    = 1.0_kpr / (2.0_kpr * PETSC_PI / input_lx * input_modes(imode))
end do
call VecSetValues(field_mode_grad_inv, nindex, indexes, values, &
  INSERT_VALUES, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyBegin(field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyEnd(field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)

! initialize Fourier matrixes
!call MatGetOwnershipRange(field_fourier_re, ix_low, ix_high, global_ierr)
!CHKERRQ(global_ierr)

nindex = input_nmode
do imode = 0, input_nmode - 1
  indexes(imode) = imode
end do
do ix = field_ix_low, field_ix_high - 1
  arrix(1) = ix
  do imode = 0, input_nmode - 1
    values(imode) = &
      cos(2.0_kpr * PETSC_PI / input_nx * real(input_modes(imode), kpr) * ix)
  end do
  call MatSetValues( &
    field_fourier_re, 1, arrix, &
    nindex, indexes, values, INSERT_VALUES, global_ierr &
  )
  do imode = 0, input_nmode - 1
    values(imode) = &
      -sin(2.0_kpr * PETSC_PI / input_nx * real(input_modes(imode), kpr) * ix)
  end do
  call MatSetValues( &
    field_fourier_im, 1, arrix, &
    nindex, indexes, values, INSERT_VALUES, global_ierr &
  )
end do
call MatAssemblyBegin(field_fourier_re, MAT_FINAL_ASSEMBLY, global_ierr)
CHKERRQ(global_ierr)
call MatAssemblyBegin(field_fourier_im, MAT_FINAL_ASSEMBLY, global_ierr)
CHKERRQ(global_ierr)
call MatAssemblyEnd(field_fourier_re, MAT_FINAL_ASSEMBLY, global_ierr)
CHKERRQ(global_ierr)
call MatAssemblyEnd(field_fourier_im, MAT_FINAL_ASSEMBLY, global_ierr)
CHKERRQ(global_ierr)

end subroutine field_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve electric field from charge !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_solve_electric
use wtimer
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: m, n
PetscReal :: norm2

call wtimer_start(global_iwt_field_electric)

! transform charge to partial Fourier space and scale by -i
call MatMultTranspose(field_fourier_re, field_chargeden, &
  field_mode_im, global_ierr)
CHKERRQ(global_ierr)
call VecScale(field_mode_im, -1.0_kpr / input_nx, global_ierr)
CHKERRQ(global_ierr)
call MatMultTranspose(field_fourier_im, field_chargeden, &
  field_mode_re, global_ierr)
CHKERRQ(global_ierr)
call VecScale(field_mode_re, 1.0_kpr / input_nx, global_ierr)
CHKERRQ(global_ierr)

! apply inverse grad to charge to get electric field in partial Fourier space
call VecPointwiseMult(field_mode_re, field_mode_re, &
  field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)
call VecPointwiseMult(field_mode_im, field_mode_im, &
  field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)

! inverse transform to real space
call MatMult(field_fourier_re, field_mode_re, field_electric, global_ierr)
CHKERRQ(global_ierr)
call MatMultAdd(field_fourier_im, field_mode_im, &
  field_electric, field_electric, global_ierr)
CHKERRQ(global_ierr)
call VecScale(field_electric, 2.0_kpr, global_ierr)
CHKERRQ(global_ierr)

!call VecNorm(field_chargeden, NORM_2, norm2, global_ierr)
!CHKERRQ(global_ierr)
!write (global_msg, *) "norm2 of charge=", norm2, "\n"
!call global_pp(global_msg)
!call VecNorm(field_electric, NORM_2, norm2, global_ierr)
!CHKERRQ(global_ierr)
!write (global_msg, *) "norm2 of electric=", norm2, "\n"
!call global_pp(global_msg)

call wtimer_start(global_iwt_field_electric)

end subroutine field_solve_electric


!!!!!!!!!!!!!!!!!!!!!
! test field solver !
!!!!!!!!!!!!!!!!!!!!!
subroutine field_test
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: ix, nindex
PetscInt, dimension(0 : input_nx - 1) :: indexes
PetscScalar, dimension(0 : input_nx - 1) :: values

nindex = field_ix_high - field_ix_low
do ix = field_ix_low, field_ix_high - 1
  indexes(ix - field_ix_low) = ix
  values(ix - field_ix_low) = cos(2.0_kpr * PETSC_PI * ix / real(input_nx, kpr))
end do
call VecSetValues(field_chargeden, nindex, indexes, values, &
  INSERT_VALUES, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyBegin(field_chargeden, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyEnd(field_chargeden, global_ierr)
CHKERRQ(global_ierr)

!call VecView(field_chargeden, PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!CHKERRQ(global_ierr)
!call MatView(field_fourier_re, PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!CHKERRQ(global_ierr)

call field_solve_electric

call VecView(field_electric, PETSC_VIEWER_STDOUT_WORLD, global_ierr)
CHKERRQ(global_ierr)

end subroutine field_test


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_final
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

call VecDestroy(field_electric, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_chargeden, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_tmp, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_electric_seq, global_ierr)
CHKERRQ(global_ierr)

call VecDestroy(field_mode_re, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_mode_im, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)

call MatDestroy(field_fourier_re, global_ierr)
CHKERRQ(global_ierr)
call MatDestroy(field_fourier_im, global_ierr)
CHKERRQ(global_ierr)

call ISDestroy(field_is_electric, global_ierr)
CHKERRQ(global_ierr)
call ISDestroy(field_is_electric_seq, global_ierr)
CHKERRQ(global_ierr)
call VecScatterDestroy(field_vs_electric, global_ierr)
CHKERRQ(global_ierr)

end subroutine field_final

end module pic1dp_field

