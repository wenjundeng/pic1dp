! manage field quantities
module pic1dp_field
use pic1dp_global
implicit none
#include "finclude/petscdef.h"

Vec :: field_electric, field_charge
Vec :: field_tmp

! for collecting charge and calculate electric field in Fourier space
Vec :: field_mode_re, field_mode_im

! inverse grad operator (scaled by i) in Fourier space
Vec :: field_mode_grad_inv

! for particle DFT and inverse DFT
Mat :: field_fourier_re, field_fourier_im

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_init
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: ix, ix_low, ix_high
PetscInt :: imode, imode_low, imode_high
PetscInt :: nindex
PetscInt, dimension(0 : input_nmode - 1) :: indexes
PetscScalar, dimension(0 : input_nmode - 1) :: values

call VecCreate(MPI_COMM_WORLD, field_electric, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(field_electric, PETSC_DECIDE, input_nx, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(field_electric, global_ierr)
CHKERRQ(global_ierr)

call VecDuplicate(field_electric, field_charge, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(field_electric, field_tmp, global_ierr)
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

! initialize Fourier matrixes
call MatCreate(MPI_COMM_WORLD, field_fourier_re, global_ierr)
CHKERRQ(global_ierr)
call MatSetType(field_fourier_re, MATDENSE, global_ierr)
CHKERRQ(global_ierr)
call MatSetSizes( &
  field_fourier_re, PETSC_DECIDE, PETSC_DECIDE, &
  input_nx, input_nmode, global_ierr &
)
CHKERRQ(global_ierr)
call MatSetFromOptions(field_fourier_re, global_ierr)
CHKERRQ(global_ierr)

call MatCreate(MPI_COMM_WORLD, field_fourier_im, global_ierr)
CHKERRQ(global_ierr)
call MatSetType(field_fourier_im, MATDENSE, global_ierr)
CHKERRQ(global_ierr)
call MatSetSizes( &
  field_fourier_im, PETSC_DECIDE, PETSC_DECIDE, &
  input_nx, input_nmode, global_ierr &
)
CHKERRQ(global_ierr)
call MatSetFromOptions(field_fourier_im, global_ierr)
CHKERRQ(global_ierr)
!call global_matcreate(field_fourier_re, input_nx, input_nmode, &
!  input_nmode, input_nmode)
!call global_matcreate(field_fourier_im, input_nx, input_nmode, &
!  input_nmode, input_nmode)

call MatGetOwnershipRange(field_fourier_re, ix_low, ix_high, global_ierr)
CHKERRQ(global_ierr)

nindex = input_nmode
indexes = (/ (imode, imode = 0, input_nmode - 1) /)
do ix = ix_low, ix_high - 1
  values = (/ ( &
    cos(2.0_kpr * PETSC_PI / input_nx * input_mode(imode) * ix), &
    imode = 0, input_nmode - 1 &
  ) /)
  call MatSetValues( &
    field_fourier_re, 1, ix, &
    nindex, indexes, values, INSERT_VALUES, global_ierr &
  )
  values = (/ ( &
    -sin(2.0_kpr * PETSC_PI / input_nx * input_mode(imode) * ix), &
    imode = 0, input_nmode - 1 &
  ) /)
  call MatSetValues( &
    field_fourier_im, 1, ix, &
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

! initialize inverse grad operator in partial Fourier space
call VecGetOwnershipRange(field_mode_grad_inv, imode_low, imode_high, global_ierr)
CHKERRQ(global_ierr)
nindex = imode_high - imode_low
do imode = imode_low, imode_high - 1
  indexes(imode - imode_low) = imode
  values(imode - imode_low) = 1.0_kpr / (2.0_kpr * PETSC_PI / input_lx * input_mode(imode))
end do
call VecSetValues(field_mode_grad_inv, nindex, indexes, values, &
  INSERT_VALUES, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyBegin(field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyEnd(field_mode_grad_inv, global_ierr)
CHKERRQ(global_ierr)

end subroutine field_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve electric field from charge !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_solve_electric
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: m, n

! transform charge to partial Fourier space and scale by -i
call MatMultTranspose(field_fourier_re, field_charge, &
  field_mode_im, global_ierr)
CHKERRQ(global_ierr)
call VecScale(field_mode_im, -1.0_kpr / input_nx, global_ierr)
CHKERRQ(global_ierr)
call MatMultTranspose(field_fourier_im, field_charge, &
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

end subroutine field_solve_electric


!!!!!!!!!!!!!!!!!!!!!
! test field solver !
!!!!!!!!!!!!!!!!!!!!!
subroutine field_test
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: ix, ix_low, ix_high, nindex
PetscInt, dimension(0 : input_nx - 1) :: indexes
PetscScalar, dimension(0 : input_nx - 1) :: values

call VecGetOwnershipRange(field_charge, ix_low, ix_high, global_ierr)
CHKERRQ(global_ierr)
nindex = ix_high - ix_low
do ix = ix_low, ix_high - 1
  indexes(ix - ix_low) = ix
  values(ix - ix_low) = cos(2.0_kpr * PETSC_PI * ix / real(input_nx, kpr))
end do
call VecSetValues(field_charge, nindex, indexes, values, &
  INSERT_VALUES, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyBegin(field_charge, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyEnd(field_charge, global_ierr)
CHKERRQ(global_ierr)

!call VecView(field_charge, PETSC_VIEWER_STDOUT_WORLD, global_ierr)
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
call VecDestroy(field_charge, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_tmp, global_ierr)
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

end subroutine field_final

end module pic1dp_field

