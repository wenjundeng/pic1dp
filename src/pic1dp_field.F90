! manage field quantities
module pic1dp_field
use pic1dp_global
implicit none
#include "finclude/petscdef.h"

Vec :: field_phi, field_charge
Vec :: field_tmp

! for collecting charge and calculate phi in Fourier space
Vec :: field_mode_re, field_mode_im

! inverse Laplacian operator in Fourier space
Vec :: field_mode_laplacian_inv

! for particle DFT and inverse DFT
Mat :: field_fourier_re, field_fourier_im

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_init
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ix, ix_low, ix_high
PetscInt :: imode, imode_low, imode_high
PetscInt :: nindex
PetscInt, dimension(0 : input_nmode - 1) :: indexes
PetscScalar, dimension(0 : input_nmode - 1) :: values

call VecCreate(MPI_COMM_WORLD, field_phi, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(field_phi, PETSC_DECIDE, input_nx, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(field_phi, global_ierr)
CHKERRQ(global_ierr)

call VecDuplicate(field_phi, field_charge, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(field_phi, field_tmp, global_ierr)
CHKERRQ(global_ierr)

call VecCreate(MPI_COMM_SELF, field_mode_re, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(field_mode_re, PETSC_DECIDE, input_nmode, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(field_mode_re, global_ierr)
CHKERRQ(global_ierr)

call VecDuplicate(field_mode_re, field_mode_im, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(field_mode_re, field_mode_laplacian_inv, global_ierr)
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

call MatGetOwnershipRange(field_fourier_re, ix_low, ix_high, global_ierr)
CHKERRQ(global_ierr)

nindex = input_nmode

do ix = ix_low, ix_high - 1
  indexes = (/ (imode, imode = 0, input_nmode - 1) /)
  values = (/ ( &
    cos(2.0_kpr * PETSC_PI / input_nx * input_mode(imode) * ix, &
    imode = 0, input_nmode - 1 &
  ) /)
  call MatSetValues( &
    field_fourier_re, 1, ix, &
    nindex, indexes, values, INSERT_VALUES, global_ierr &
  )
  values = (/ ( &
    sin(2.0_kpr * PETSC_PI / input_nx * input_mode(imode) * ix, &
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

call VecGetOwnershipRange(field_mode_laplacian_inv, imode_low, imode_high, global_ierr)
CHKERRQ(global_ierr)
nindex = imode_high - imode_low
indexes = (/ (imode, imode = imode_low, imode_high - 1) /)
values = (/ (1.0_kpr / (2.0_kpr * PETSC_PI / input_lx * input_mode(imode))**2, &
  imode = imode_low, imode_high - 1) /)
call VecSetValues(field_mode_laplacian_inv, nindex, indexes, values, &
  INSERT_VALUES, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyBegin(field_mode_laplacian_inv, global_ierr)
CHKERRQ(global_ierr)
call VecAssemblyEnd(field_mode_laplacian_inv, global_ierr)
CHKERRQ(global_ierr)

end subroutine field_init


!!!!!!!!!!!!!!!!!!!!!!!!!
! solve phi from charge !
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_solve_phi
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

! transform charge to partial Fourier space
call MatMultTranspose(field_fourier_re, field_charge, &
  field_mode_re, global_ierr)
CHKERRQ(global_ierr)
call MatMultTranspose(field_fourier_im, field_charge, &
  field_mode_im, global_ierr)
CHKERRQ(global_ierr)

! apply inverse Laplacian to charge to get phi in partial Fourier space
call VecPointwiseMult(field_fourier_re, field_fourier_re, &
  field_mode_laplacian_inv, global_ierr)
CHKERRQ(global_ierr)
call VecPointwiseMult(field_fourier_im, field_fourier_im, &
  field_mode_laplacian_inv, global_ierr)
CHKERRQ(global_ierr)

! inverse transform to real space

end subroutine field_solve_phi


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine field_final
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

call VecDestroy(field_phi, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_charge, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_tmp, global_ierr)
CHKERRQ(global_ierr)

call VecDestroy(field_mode_re, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_mode_im, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(field_mode_laplacian_inv, global_ierr)
CHKERRQ(global_ierr)

call MatDestroy(field_fourier_re, global_ierr)
CHKERRQ(global_ierr)
call MatDestroy(field_fourier_im, global_ierr)
CHKERRQ(global_ierr)

end subroutine field_final

end module pic1dp_field

