! Copyright 2012, 2013 Wenjun Deng <wdeng@wdeng.info>
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


! module for managing global constants, variables, and subroutines
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
PetscInt, parameter :: kpi = kind(kpr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! wall clock timer parameters and variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! wall clock timer indexes
PetscInt, parameter :: &
  global_iwt_total = 1, &
  global_iwt_init = 2, &
  global_iwt_particle_load = 3, &
  global_iwt_push_particle = 4, &
  global_iwt_particle_shape = 5, &
  global_iwt_collect_charge = 6, &
  global_iwt_field_electric = 7, &
  global_iwt_particle_optimize = 8, &
  global_iwt_output = 9, &
  global_iwt_final = 10, &
  global_iwt_mpiallredu = 21, &
  global_iwt_scatter = 22


!!!!!!!!!!!!!!!!!!!!!
! program variables !
!!!!!!!!!!!!!!!!!!!!!

PetscInt :: global_mype ! rank of current MPI process
PetscInt :: global_npe ! # of MPI processes
PetscErrorCode :: global_ierr ! for storing MPI and PETSc error code
character(len = 5000) :: global_msg ! for messages and temporary strings

PetscInt :: global_itime ! indexing time step
PetscReal :: global_time ! physical time
PetscInt :: global_irk ! indexing Runge-Kutta sub-step

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

!call MatSeqAIJSetPreallocation( &
!  mattocreate, d_nz, PETSC_NULL_INTEGER, global_ierr &
!)
!CHKERRQ(global_ierr)
!call MatMPIAIJSetPreallocation( &
!  mattocreate, d_nz, PETSC_NULL_INTEGER, &
!  o_nz, PETSC_NULL_INTEGER, global_ierr &
!)
!CHKERRQ(global_ierr)

! use default preallocation to avoid problem with manual preallocation
call MatSetUp(mattocreate, global_ierr)
CHKERRQ(global_ierr)

call MatSetFromOptions(mattocreate, global_ierr)
CHKERRQ(global_ierr)

end subroutine global_matcreate

end module pic1dp_global

