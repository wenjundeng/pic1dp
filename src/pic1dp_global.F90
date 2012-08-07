! global constants and variables
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

end module pic1dp_global

