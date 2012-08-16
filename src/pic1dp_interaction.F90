! manage field-particle interactions
module pic1dp_interaction
use pic1dp_field
use pic1dp_particle
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collect charge from particle to field grid !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaction_collect_charge
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, ip, ix
PetscScalar, dimension(:), pointer :: pw, pc, px
PetscScalar :: sx

if (input_iptclshape < 3) then
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
else
  ! first collect local portion of particle charges to field_arr_charge2
  field_arr_charge1 = 0.0_kpr
  field_arr_charge2 = 0.0_kpr
  do ispecies = 1, input_nspecies
    call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
    CHKERRQ(global_ierr)
    if (input_iptclshape == 4) then
      call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
      CHKERRQ(global_ierr)
    end if
    do ip = particle_ip_low, particle_ip_high - 1
      if (input_iptclshape == 3) then
        ix = particle_shape_x_indexes(ispecies, ip)
        sx = particle_shape_x_values(ispecies, ip)
      else ! input_iptclshape == 4
        ! enforce periodic boundary condition
        px(ip - particle_ip_low + 1) &
          = mod(px(ip - particle_ip_low + 1), input_lx)
        ! if x is negative, mod gives negative result, so shift it to positive
        if (px(ip - particle_ip_low + 1) < 0.0_kpr) then
          px(ip - particle_ip_low + 1) &
            = px(ip - particle_ip_low + 1) + input_lx
        end if

        sx = px(ip - particle_ip_low + 1) / input_lx * input_nx
        ix = floor(sx)
        sx = 1.0_kpr - (sx - real(ix, kpr))
      end if
      field_arr_charge1(ix) = field_arr_charge1(ix) &
        + sx * pw(ip - particle_ip_low + 1)
      ix = ix + 1
      if (ix > input_nx - 1) ix = 0
      field_arr_charge1(ix) = field_arr_charge1(ix) &
        + (1.0_kpr - sx) * pw(ip - particle_ip_low + 1)
    end do ! ip = particle_ip_low, particle_ip_high - 1
    if (input_iptclshape == 4) then
      call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
      CHKERRQ(global_ierr)
    end if
    call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
    CHKERRQ(global_ierr)
    field_arr_charge2 = field_arr_charge2 &
      + field_arr_charge1 * input_charge(ispecies)
  end do ! ispecies = 1, input_nspecies

  ! then use MPI_Allreduce to get charges of all particles in field_arr_charge1
  call MPI_Allreduce(field_arr_charge2, field_arr_charge1, input_nx, &
    MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD, global_ierr)

  ! copy values from field_arr_charge1 to field_charge
  call VecGetArrayF90(field_charge, pc, global_ierr)
  CHKERRQ(global_ierr)
  do ix = field_ix_low, field_ix_high - 1
    pc(ix - field_ix_low + 1) = field_arr_charge1(ix)
  end do
  call VecRestoreArrayF90(field_charge, pc, global_ierr)
  CHKERRQ(global_ierr)
end if

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

PetscInt :: ispecies, ix
PetscInt :: ip
PetscScalar :: sx
PetscScalar, dimension(:), pointer :: pe, pp, px

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

if (input_iptclshape == 3 .or. input_iptclshape == 4) then
  call VecScatterBegin( &
    field_vs_electric, field_electric, field_electric_seq, &
    INSERT_VALUES, SCATTER_FORWARD, global_ierr &
  )
  CHKERRQ(global_ierr)
  call VecScatterEnd( &
    field_vs_electric, field_electric, field_electric_seq, &
    INSERT_VALUES, SCATTER_FORWARD, global_ierr &
  )
  CHKERRQ(global_ierr)
  call VecGetArrayF90(field_electric_seq, pe, global_ierr)
  CHKERRQ(global_ierr)
end if

do ispecies = 1, input_nspecies
  ! calculate electric field at particle
  if (input_iptclshape < 3) then
    call MatMult(particle_shape_x(ispecies), field_electric, &
      particle_tmp1, global_ierr)
    CHKERRQ(global_ierr)
  else
    call VecGetArrayF90(particle_tmp1, pp, global_ierr)
    CHKERRQ(global_ierr)
    if (input_iptclshape == 4) then
      call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
      CHKERRQ(global_ierr)
    end if
    do ip = particle_ip_low, particle_ip_high - 1
      if (input_iptclshape == 3) then
        ix = particle_shape_x_indexes(ispecies, ip)
        sx = particle_shape_x_values(ispecies, ip)
      else
        ! periodic boundary condition is enforced in interaction_collect_charge
        ! no need to enforce boundary condition here

        sx = px(ip - particle_ip_low + 1) / input_lx * input_nx
        ix = floor(sx)
        sx = 1.0_kpr - (sx - real(ix, kpr))
      end if
!      if (global_mype == 0 .and. ip == particle_ip_low) write (*, *) ix, sx
      pp(ip - particle_ip_low + 1) = pe(ix + 1) * sx
      ix = ix + 1
      if (ix > input_nx - 1) ix = 0
      pp(ip - particle_ip_low + 1) = &
        pp(ip - particle_ip_low + 1) + pe(ix + 1) * (1.0_kpr - sx)
    end do
    if (input_iptclshape == 4) then
      call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
      CHKERRQ(global_ierr)
    end if
    call VecRestoreArrayF90(particle_tmp1, pp, global_ierr)
    CHKERRQ(global_ierr)
  end if
  ! at this point particle_tmp1 stores electric field at particle positions

  ! push x
  call VecWAXPY(particle_x(ispecies), dt, particle_v(ispecies), &
    particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  ! periodic boundary condition is enforced in particle_compute_shape_x
  ! or interaction_collect_charge, depending on input_iptclshape

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
end do ! ispecies = 1, input_nspecies

if (input_iptclshape == 3 .or. input_iptclshape == 4) then
  call VecRestoreArrayF90(field_electric_seq, pe, global_ierr)
  CHKERRQ(global_ierr)
end if

end subroutine interaction_push_particle

end module pic1dp_interaction

