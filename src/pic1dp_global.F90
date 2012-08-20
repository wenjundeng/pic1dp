! global constants, variables, and subroutines
module pic1dp_global
implicit none
#include "finclude/petscdef.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the kind # of PetscReal !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a real constant to test the kind # of PetscReal
PetscReal, parameter :: testkindPetscReal = 0d0 
! kind # of PetscReal
PetscInt, parameter :: kpr = kind(testkindPetscReal)


!!!!!!!!!!!!!!!!!!!!!
! program variables !
!!!!!!!!!!!!!!!!!!!!!

PetscInt :: global_mype ! rank of current MPI process
PetscErrorCode :: global_ierr ! for storing MPI and PETSc error code
character(len = 5000) :: global_msg ! for messages and temporary strings

PetscInt :: global_itime ! indexing time step
PetscReal :: global_time ! physical time

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!
! wrapper of PetscPrintf !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine global_pp(string, adjl)
implicit none
#include "finclude/petsc.h90"

character(len = *), intent(in) :: string ! string to print
logical, intent(in), optional :: adjl ! whether to apply adjustl to string

logical :: adjl_act

adjl_act = .false.
if (present(adjl)) adjl_act = adjl

if (adjl_act) then
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(string)), global_ierr)
else
  call PetscPrintf(MPI_COMM_WORLD, trim(string), global_ierr)
end if
CHKERRQ(global_ierr)

end subroutine global_pp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create and initialize matrix !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine global_matcreate(mattocreate, nrow, ncol, d_nz, o_nz)
implicit none
#include "finclude/petsc.h90"

Mat, intent(in) :: mattocreate ! matrix to create
PetscInt, intent(in) :: &
  nrow, ncol, & ! # of rows and columns of the matrix
  d_nz, & ! # of nonzeros per row
  o_nz ! # of nonzeros per row in off-diagonal portion

call MatCreate(MPI_COMM_WORLD, mattocreate, global_ierr)
CHKERRQ(global_ierr)
call MatSetType(mattocreate, MATAIJ, global_ierr)
CHKERRQ(global_ierr)
call MatSetSizes( &
  mattocreate, PETSC_DECIDE, PETSC_DECIDE, &
  nrow, ncol, global_ierr &
)
CHKERRQ(global_ierr)
call MatSeqAIJSetPreallocation( &
  mattocreate, d_nz, PETSC_NULL_INTEGER, global_ierr &
)
CHKERRQ(global_ierr)
call MatMPIAIJSetPreallocation( &
  mattocreate, d_nz, PETSC_NULL_INTEGER, &
  o_nz, PETSC_NULL_INTEGER, global_ierr &
)
CHKERRQ(global_ierr)
call MatSetFromOptions(mattocreate, global_ierr)
CHKERRQ(global_ierr)

end subroutine global_matcreate

end module pic1dp_global

