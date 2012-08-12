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
use pic1dp_global
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
subroutine interaction_push_particle(irk)
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(in) :: irk ! index of Runge-Kutta sub-step

PetscInt :: ispecies

PetscScalar :: dt ! time step

if (irk == 1) then
  dt = 0.5_kpr * input_dt
  do ispecies = 1, input_nspecies
    call VecCopy(particle_x(ispecies), particle_x_bak(ispecies), global_ierr)
    CHKERRQ(global_ierr)
    call VecCopy(particle_v(ispecies), particle_v_bak(ispecies), global_ierr)
    CHKERRQ(global_ierr)
    call VecCopy(particle_w(ispecies), particle_w_bak(ispecies), global_ierr)
    CHKERRQ(global_ierr)
  end do
else
  ! 2nd Runge-Kutta sub-step
  dt = input_dt
end if

do ispecies = 1, input_nspecies
  ! push x
  call VecWAXPY(particle_x(ispecies), dt, particle_v(ispecies), &
    particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  ! periodic boundary condition is enforced in particle_compute_shape_x

  ! calculate electric field at particle
  call MatMult(particle_shape_x(ispecies), field_electric, &
    particle_tmp1, global_ierr)
  CHKERRQ(global_ierr)
  ! at this point particle_tmp1 stores electric field

  ! push w
  if (input_linear == 1) then
    call VecCopy(particle_p(ispecies), particle_tmp2, global_ierr)
    CHKERRQ(global_ierr)
  else
    call VecWAXPY(particle_tmp2, -1.0_kpr, particle_w(ispecies), &
      particle_p(ispecies), global_ierr)
    CHKERRQ(global_ierr)
  end if
  call VecPointwiseMult(particle_tmp2, &
    particle_tmp2, particle_tmp1, global_ierr)
  CHKERRQ(global_ierr)
  call VecPointwiseMult(particle_tmp2, &
    particle_tmp2, particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecWAXPY( &
    particle_w(ispecies), &
    dt * input_charge(ispecies) / input_temperature(ispecies), &
    particle_tmp2, particle_w_bak(ispecies), global_ierr &
  )
  CHKERRQ(global_ierr)

  ! push v
  ! (this must be after push x and w because dx/dt and dw/dt depend on v)
  if (input_linear /= 1) then
    call VecWAXPY( &
      particle_v(ispecies), &
      dt * input_charge(ispecies) / input_mass(ispecies), &
      particle_tmp1, particle_v_bak(ispecies), global_ierr &
    )
  end if
end do

end subroutine interaction_push_particle

end module pic1dp_interaction

