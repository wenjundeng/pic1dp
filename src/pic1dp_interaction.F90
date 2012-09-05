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


! module for managing field-particle interactions
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
use wtimer
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, ip, ix
PetscScalar, dimension(:), pointer :: pw, pc, px, ps
PetscScalar :: sx

if (input_iptclshape < 3) then
  call VecSet(field_chargeden, 0.0_kpr, global_ierr)
  CHKERRQ(global_ierr)

  do ispecies = 1, input_nspecies
    if (input_deltaf == 1) then
      call MatMultTranspose( &
        particle_shape_x(ispecies), particle_w(ispecies), &
        field_tmp, global_ierr &
      )
      CHKERRQ(global_ierr)
    else
      ! full-f case, first use p to calculate total density
      call MatMultTranspose( &
        particle_shape_x(ispecies), particle_p(ispecies), &
        field_tmp, global_ierr &
      )
      CHKERRQ(global_ierr)
      ! then subtract equilibrium density to get perturbed density
      call VecGetArrayF90(field_tmp, pc, global_ierr)
      CHKERRQ(global_ierr)
      pc(:) = pc(:) - input_species_density(ispecies) * input_lx / input_nx
      call VecRestoreArrayF90(field_tmp, pc, global_ierr)
      CHKERRQ(global_ierr)
    end if
    call VecAXPY(field_chargeden, &
      input_species_charge(ispecies), field_tmp, global_ierr)
    CHKERRQ(global_ierr)
  end do
  ! at this point field_chargeden stores charge on a grid
  ! now divide by grid size to get charge density
  call VecScale(field_chargeden, input_nx / input_lx, global_ierr)
  CHKERRQ(global_ierr)
else
  ! first collect local portion of particle charges to field_arr_charge2
  field_arr_charge2 = 0.0_kpr
  do ispecies = 1, input_nspecies
    field_arr_charge1 = 0.0_kpr
    if (input_deltaf == 1) then
      call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
      CHKERRQ(global_ierr)
    else
      ! full-f case, first use p to get total density
      call VecGetArrayF90(particle_p(ispecies), pw, global_ierr)
      CHKERRQ(global_ierr)
    end if
    if (input_iptclshape == 4) then
      call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
      CHKERRQ(global_ierr)
    end if
    call VecGetArrayF90(particle_s(ispecies), ps, global_ierr)
    CHKERRQ(global_ierr)
    do ip = 1, particle_ip_high - particle_ip_low
      ! ignore invalid particles
      if (ps(ip) < -0.5_kpr) cycle

      if (input_iptclshape == 3) then
        ix = particle_shape_x_indexes(ispecies, ip)
        sx = particle_shape_x_values(ispecies, ip)
      else ! input_iptclshape == 4
        ! enforce periodic boundary condition
        px(ip) = mod(px(ip), input_lx)
        ! if x is negative, mod gives negative result, so shift it to positive
        if (px(ip) < 0.0_kpr) px(ip) = px(ip) + input_lx

        sx = px(ip) / input_lx * input_nx
        ix = floor(sx)
        sx = 1.0_kpr - (sx - real(ix, kpr))
      end if
      if (ix > input_nx - 1) then
        write (*, *) 'global_mype, ip, x, sx, ix = ', global_mype, ip, px(ip), sx, ix
      end if
      field_arr_charge1(ix) = field_arr_charge1(ix) + sx * pw(ip)
      ix = ix + 1
      if (ix > input_nx - 1) ix = 0
      field_arr_charge1(ix) = field_arr_charge1(ix) + (1.0_kpr - sx) * pw(ip)
    end do ! ip = 1, particle_ip_high - particle_ip_low
    if (input_iptclshape == 4) then
      call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
      CHKERRQ(global_ierr)
    end if
    if (input_deltaf == 1) then
      call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
      CHKERRQ(global_ierr)
    else
      call VecRestoreArrayF90(particle_p(ispecies), pw, global_ierr)
      CHKERRQ(global_ierr)
    end if
    call VecRestoreArrayF90(particle_s(ispecies), ps, global_ierr)
    CHKERRQ(global_ierr)
    field_arr_charge2(:) = field_arr_charge2(:) &
      + field_arr_charge1(:) * input_species_charge(ispecies)
  end do ! ispecies = 1, input_nspecies

  call wtimer_start(21)
  ! then use MPI_Allreduce to get charges of all particles in field_arr_charge1
  call MPI_Allreduce(field_arr_charge2, field_arr_charge1, input_nx, &
    MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  call wtimer_stop(21)

  ! copy values from field_arr_charge1 to field_chargeden
  call VecGetArrayF90(field_chargeden, pc, global_ierr)
  CHKERRQ(global_ierr)
  pc(1 : field_ix_high - field_ix_low) &
    = field_arr_charge1(field_ix_low : field_ix_high - 1) * input_nx / input_lx
  if (input_deltaf == 0) then
    ! full-f case, subtract equilibrium charge density
    do ispecies = 1, input_nspecies
      pc(:) = pc(:) &
        - input_species_charge(ispecies) * input_species_density(ispecies)
    end do
  end if
  call VecRestoreArrayF90(field_chargeden, pc, global_ierr)
  CHKERRQ(global_ierr)
end if

end subroutine interaction_collect_charge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! using field to push particle !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaction_push_particle(irk)
use pic1dp_global
use pic1dp_input
use wtimer
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(in) :: irk ! index of Runge-Kutta sub-step

PetscInt :: ispecies, ix
PetscInt :: ip
PetscScalar :: sx, electric, tmp1, tmp2
PetscScalar, dimension(:), pointer :: pe, pptcl, ppe, px, pv, pp, pw, ps
PetscScalar, dimension(:), pointer :: pxb, pvb, pwb

PetscScalar :: dt ! time step

if (irk == 1) then
  dt = 0.5_kpr * input_dt
  do ispecies = 1, input_nspecies
    call VecCopy(particle_x(ispecies), particle_x_bak(ispecies), global_ierr)
    CHKERRQ(global_ierr)
    call VecCopy(particle_v(ispecies), particle_v_bak(ispecies), global_ierr)
    CHKERRQ(global_ierr)
    if (input_deltaf == 1) then
      call VecCopy(particle_w(ispecies), particle_w_bak(ispecies), global_ierr)
      CHKERRQ(global_ierr)
    end if
  end do
else
  ! 2nd Runge-Kutta sub-step
  dt = input_dt
end if

call wtimer_start(22)
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
call wtimer_stop(22)

do ispecies = 1, input_nspecies
  if (input_iptclshape < 3) then
    ! calculate electric field at particle if using PETSc for shape matrix
    call MatMult(particle_shape_x(ispecies), field_electric, &
      particle_electric, global_ierr)
    CHKERRQ(global_ierr)
    call VecGetArrayF90(particle_electric, ppe, global_ierr)
    CHKERRQ(global_ierr)
  end if
  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)

  call VecGetArrayF90(particle_x_bak(ispecies), pxb, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v_bak(ispecies), pvb, global_ierr)
  CHKERRQ(global_ierr)
  if (input_deltaf == 1) then
    call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
    CHKERRQ(global_ierr)
    call VecGetArrayF90(particle_w_bak(ispecies), pwb, global_ierr)
    CHKERRQ(global_ierr)
  end if
  do ip = 1, particle_ip_high - particle_ip_low
    ! ignore invalid particles
    if (ps(ip) < -0.5_kpr) cycle

    ! get electric field at particle
    if (input_iptclshape <= 2) then
      electric = ppe(ip)
    else
      if (input_iptclshape == 3) then
        ix = particle_shape_x_indexes(ispecies, ip)
        sx = particle_shape_x_values(ispecies, ip)
      else ! indicating (input_iptclshape == 4)
        ! periodic boundary condition is enforced in interaction_collect_charge
        ! no need to enforce boundary condition here

        sx = px(ip) / input_lx * input_nx
        ix = floor(sx)
        sx = 1.0_kpr - (sx - real(ix, kpr))
      end if
      electric = pe(ix + 1) * sx
      ix = ix + 1
      if (ix > input_nx - 1) ix = 0
      electric = electric + pe(ix + 1) * (1.0_kpr - sx)
    end if

    ! push x
    px(ip) = pxb(ip) + dt * pv(ip)
    ! periodic boundary condition is enforced in particle_compute_shape_x
    ! or interaction_collect_charge, depending on input_iptclshape

    ! push w
    if (input_deltaf == 1) then
      ! calculate (p - w) * E (or p * E for linear) and store in tmp1
      if (input_linear == 1) then
        tmp1 = pp(ip) * electric
      else
        tmp1 = (pp(ip) - pw(ip)) * electric
      end if

      ! calculate -d f_0 / d v / f_0 and store in tmp2
      if (input_iptcldist == 1) then ! two-stream1
        tmp2 = pv(ip) - 2.0_kpr / pv(ip)
      elseif (input_iptcldist == 2) then ! two-stream2
        tmp2 = ((pv(ip) + input_species_v0(ispecies)) &
          * exp(-(pv(ip) + input_species_v0(ispecies))**2 / (2.0_kpr &
          * input_species_temperature(ispecies) &
          / input_species_mass(ispecies))) &
          + (pv(ip) - input_species_v0(ispecies)) &
          * exp(-(pv(ip) - input_species_v0(ispecies))**2 / (2.0_kpr &
          * input_species_temperature(ispecies) &
          / input_species_mass(ispecies)))) &
          / (exp(-(pv(ip) + input_species_v0(ispecies))**2 / (2.0_kpr &
          * input_species_temperature(ispecies) &
          / input_species_mass(ispecies))) &
          + exp(-(pv(ip) - input_species_v0(ispecies))**2 / (2.0_kpr &
          * input_species_temperature(ispecies) &
          / input_species_mass(ispecies)))) &
          * input_species_mass(ispecies) / input_species_temperature(ispecies)
      else ! (shifted) Maxwellian
        tmp2 = (pv(ip) - input_species_v0(ispecies)) &
          * input_species_temperature(ispecies) &
          / input_species_mass(ispecies)
      end if
      ! end of calculating -d f_0 / d v / f_0 and storing in tmp2

      pw(ip) = pwb(ip) + dt * tmp1 * tmp2 * input_species_charge(ispecies) &
        / input_species_mass(ispecies)
    end if ! (input_deltaf == 1)

    ! push v
    ! (this must be after push x and w because dx/dt and dw/dt depend on v)
    if (input_linear == 0) then
      pv(ip) = pvb(ip) + dt * electric * input_species_charge(ispecies) &
        / input_species_mass(ispecies)
    end if 
  end do ! ip = 1, particle_ip_high - particle_ip_high
  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_x_bak(ispecies), pxb, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v_bak(ispecies), pvb, global_ierr)
  CHKERRQ(global_ierr)
  if (input_deltaf == 1) then
    call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
    CHKERRQ(global_ierr)
    call VecRestoreArrayF90(particle_w_bak(ispecies), pwb, global_ierr)
    CHKERRQ(global_ierr)
  end if
  if (input_iptclshape < 3) then
    call VecRestoreArrayF90(particle_electric, ppe, global_ierr)
    CHKERRQ(global_ierr)
  end if
end do ! ispecies = 1, input_nspecies

if (input_iptclshape == 3 .or. input_iptclshape == 4) then
  call VecRestoreArrayF90(field_electric_seq, pe, global_ierr)
  CHKERRQ(global_ierr)
end if

end subroutine interaction_push_particle

end module pic1dp_interaction

