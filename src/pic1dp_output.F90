module pic1dp_output
implicit none
#include "finclude/petscdef.h"

Vec :: output_ptcldist

PetscViewer :: output_viewer

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_init
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

call VecCreate(MPI_COMM_WORLD, output_ptcldist, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(output_ptcldist, PETSC_DECIDE, &
  input_nx * input_output_nv, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(output_ptcldist, global_ierr)
CHKERRQ(global_ierr)

call PetscViewerBinaryOpen( &
  MPI_COMM_WORLD, 'pic1dp.out', FILE_MODE_WRITE, &
  output_viewer, global_ierr &
)
CHKERRQ(global_ierr)

end subroutine output_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output marker particle distribution on x-v plane !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_xv_ptcldist
use pic1dp_global
use pic1dp_input
use pic1dp_particle
implicit none
#include "finclude/petsc.h90"

PetscInt, dimension(2) :: intbuf
PetscReal, dimension(2) :: realbuf

call particle_compute_shape_xv
call MatMultTranspose(particle_shape_xv, particle_p, &
  output_ptcldist, global_ierr)
CHKERRQ(global_ierr)

call VecView(output_ptcldist, PETSC_VIEWER_STDOUT_WORLD, global_ierr)
CHKERRQ(global_ierr)

realbuf(1) = input_lx
realbuf(2) = input_output_v_max
call PetscViewerBinaryWriteReal(output_viewer, realbuf, 2, &
  PETSC_TRUE, global_ierr)
CHKERRQ(global_ierr)

intbuf(1) = input_nx
intbuf(2) = input_output_nv
call PetscViewerBinaryWriteInt(output_viewer, intbuf, 2, &
  PETSC_TRUE, global_ierr)
CHKERRQ(global_ierr)

call VecView(output_ptcldist, output_viewer, global_ierr)
CHKERRQ(global_ierr)

end subroutine output_xv_ptcldist


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

call VecDestroy(output_ptcldist, global_ierr)
CHKERRQ(global_ierr)
call PetscViewerDestroy(output_viewer, global_ierr)
CHKERRQ(global_ierr)

end subroutine output_final

end module pic1dp_output

