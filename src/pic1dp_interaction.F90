! manage field-particle interactions
module pic1dp_interaction
use pic1dp_field
use pic1dp_particle
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collect charge from particle to field grid              !
! particle_compute_shape_x needs to be called before this !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaction_collect_charge
use pic1d_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

call VecSet(field_charge, 0.0_kpr, global_ierr)
CHKERRQ(global_ierr)

do ispecies = 1, input_nspecies
  call MatMultTranspose( &
    particle_shape_x(ispecies), particle_w(ispecies), &
    field_tmp, global_ierr &
  )
  CHKERRQ(global_ierr)
  call VecAXPY(field_charge, input_charge(ispecies), field_tmp, global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine interaction_collect_charge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! using field to push particle !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaction_push_particle
use pic1d_global
implicit none
#include "finclude/petsc.h90"

end subroutine interaction_push_particle

end module pic1dp_interaction

