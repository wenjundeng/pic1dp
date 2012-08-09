! manage particles
module pic1dp_particle
use pic1dp_input
implicit none
#include "finclude/petscdef.h"

! particle x coordinate, velocity, p = f / g, w = delta f / g
! f is total distribution, delta f is perturbed distribution
! g is marker distribution
Vec, dimension(input_nspecies) :: &
  particle_x, particle_v, particle_p, particle_w

Mat, dimension(input_nspecies) :: particle_mat_shape_xv

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_init
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

call VecCreate(MPI_COMM_WORLD, particle_x(1), global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(particle_x(1), PETSC_DECIDE, input_nparticle, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(particle_x(1), global_ierr)
CHKERRQ(global_ierr)

do ispecies = 2, input_nspecies
  call VecDuplicate(particle_x(1), particle_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do

do ispecies = 1, input_nspecies
  call VecDuplicate(particle_x(1), particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_p(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_w(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine particle_init


!!!!!!!!!!!!!!!!!!
! load particles !
!!!!!!!!!!!!!!!!!!
subroutine particle_load
use pic1dp_global
use pic1dp_input
use gaussian
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

PetscInt :: ip_low, ip_high, ip, nindex
PetscInt, dimension(:), allocatable :: indexes
PetscScalar, dimension(:), allocatable :: values

call gaussian_init(input_seed_type, global_mype)

call VecGetOwnershipRange(particle_x(1), ip_low, ip_high, global_ierr)
CHKERRQ(global_ierr)
allocate (indexes(0 : ip_high - ip_low), values(0 : ip_high - ip_low))
nindex = ip_high - ip_low
indexes = (/ (ip, ip = ip_low, ip_high - 1) /)
do ispecies = 1, input_nspecies
  call gaussian_generate(values)

  call VecSetValues(particle_v(ispecies), nindex, indexes, values, &
    INSERT_VALUES, global_ierr)
  CHKERRQ(global_ierr)
  call VecAssemblyBegin(particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call random_number(values)
  values = values * input_lx
  call VecSetValues(particle_x(ispecies), nindex, indexes, values, &
    INSERT_VALUES, global_ierr)
  CHKERRQ(global_ierr)
  call VecAssemblyBegin(particle_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call VecAssemblyEnd(particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecAssemblyEnd(particle_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do
deallocate (indexes, values)

end subroutine particle_load


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! construct shape matrix in x-v plane !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_construct_mat_shape_xv
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip_low, ip_high, ip
PetscInt :: ix1, ix2, iv, ixv
PetscScalar :: sx, sv
PetscScalar, dimension(:), pointer :: px, pv
PetscInt, dimension(0 : 3) :: indexes
PetscScalar, dimension(0 : 3) :: values

call VecGetOwnershipRange(particle_x(1), ip_low, ip_high, global_ierr)
CHKERRQ(global_ierr)

do ispecies = 1, input_nspecies
  call MatDestroy(particle_mat_shape_xv(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call global_matcreate(particle_mat_shape_xv(ispecies), &
    input_nparticle, input_nx * input_output_nv, 4, 4)

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)

  do ip = ip_low, ip_high - 1
    if ( &
      pv(ip - ip_low + 1) <= -input_output_v_max &
      .or. pv(ip - ip_low + 1) >= input_output_v_max &
    ) cycle

    sx = px(ip - ip_low + 1) / input_lx * input_nx
    ix1 = floor(sx)
    sx = sx - real(ix1, kpr)
    ix2 = ix1 + 1
    if (ix2 == input_nx) ix2 = 0

    sv = (pv(ip - ip_low + 1) + input_output_v_max) &
      / (input_output_v_max * 2.0_kpr) * (input_output_nv - 1)
    iv = floor(sv)
    sv = sv - real(iv, kpr)

    nindex = 4
    indexes(0) = iv * input_nx + ix1
    indexes(1) = iv * input_nx + ix2
    indexes(2) = (iv + 1) * input_nx + ix1
    indexes(3) = (iv + 1) * input_nx + ix2
    values(0) = sx * sv
    values(1) = (1.0_kpr - sx) * sv
    values(2) = sx * (1.0_kpr - sv)
    values(3) = (1.0_kpr - sx) * (1.0_kpr - sv)
    call MatSetValues( &
      particle_mat_shape_xv(ispecies), 1, ip, &
      nindex, indexes, values, INSERT_VALUES, global_ierr &
    )
    CHKERRQ(global_ierr)
  end do
  call MatAssemblyBegin(particle_mat_shape_xv(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)

  call MatAssemblyEnd(particle_mat_shape_xv(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine particle_construct_mat_shape_xv


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_final
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

do ispecies = 1, input_nspecies
  call VecDestroy(particle_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_p(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_w(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call MatDestroy(particle_mat_shape_xv(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine particle_final

end module pic1dp_particle

