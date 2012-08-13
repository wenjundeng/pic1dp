! manage particles
module pic1dp_particle
use pic1dp_input
implicit none
#include "finclude/petscdef.h"

! particle x coordinate, velocity, p = f / g, w = delta f / g
! f is total distribution, delta f is perturbed distribution
! g is marker distribution
Vec, dimension(input_nspecies) :: &
  particle_x, particle_v, particle_p, particle_w, &
  particle_x_bak, particle_v_bak, particle_w_bak

! temporary particle vectors, need to be used 
! in pic1dp_interaction and pic1dp_output
Vec :: particle_tmp1, particle_tmp2

Mat, dimension(input_nspecies) :: &
  particle_shape_xv, particle_shape_x, particle_shape_v

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

  call VecDuplicate(particle_x(1), particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_v_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDuplicate(particle_x(1), particle_w_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call global_matcreate(particle_shape_xv(ispecies), &
    input_nparticle, input_nx * input_output_nv, 4, 4)
  call global_matcreate(particle_shape_x(ispecies), &
    input_nparticle, input_nx, 2, 2)
  call global_matcreate(particle_shape_v(ispecies), &
    input_nparticle, input_output_nv, 2, 2)
end do

call VecDuplicate(particle_x(1), particle_tmp1, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(particle_x(1), particle_tmp2, global_ierr)
CHKERRQ(global_ierr)

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
allocate (indexes(0 : ip_high - ip_low - 1), values(0 : ip_high - ip_low - 1))
nindex = ip_high - ip_low
indexes = (/ (ip, ip = ip_low, ip_high - 1) /)
do ispecies = 1, input_nspecies
  call gaussian_generate(values)
  values = values * sqrt(input_temperature(ispecies) / input_mass(ispecies))

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

  call VecSet(particle_p(ispecies), input_lx / input_nparticle, global_ierr)
  CHKERRQ(global_ierr)
  call VecSet(particle_w(ispecies), 0.0_kpr, global_ierr)
  CHKERRQ(global_ierr)
end do
deallocate (indexes, values)

!  call VecView(particle_x(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_p(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_w(1), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)

end subroutine particle_load


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute shape matrix in x-v plane !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_shape_xv
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip_low, ip_high, ip
PetscInt :: ip_low1, ip_high1
PetscInt :: ix1, ix2, iv, ixv
PetscScalar :: sx, sv
PetscScalar, dimension(:), pointer :: px, pv
PetscInt, dimension(0 : 3) :: indexes
PetscScalar, dimension(0 : 3) :: values


call VecGetOwnershipRange(particle_x(1), ip_low, ip_high, global_ierr)
CHKERRQ(global_ierr)

do ispecies = 1, input_nspecies
!  call MatDestroy(particle_shape_xv(ispecies), global_ierr)
!  CHKERRQ(global_ierr)
  call MatZeroEntries(particle_shape_xv(ispecies), global_ierr)
  CHKERRQ(global_ierr)

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
    values(0) = (1.0_kpr - sx) * (1.0_kpr - sv)
    values(1) = sx * (1.0_kpr - sv)
    values(2) = (1.0_kpr - sx) * sv
    values(3) = sx * sv
    call MatSetValues( &
      particle_shape_xv(ispecies), 1, ip, &
      nindex, indexes, values, INSERT_VALUES, global_ierr &
    )
    CHKERRQ(global_ierr)
  end do
  call MatAssemblyBegin(particle_shape_xv(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)
  call MatAssemblyEnd(particle_shape_xv(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)


!  call VecView(particle_x(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call MatView(particle_shape_xv(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
end do

end subroutine particle_compute_shape_xv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute shape matrix in x !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_shape_x
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip_low, ip_high, ip
PetscInt :: ix1, ix2
PetscScalar :: sx
PetscScalar, dimension(:), pointer :: px
PetscInt, dimension(0 : 1) :: indexes
PetscScalar, dimension(0 : 1) :: values

call VecGetOwnershipRange(particle_x(1), ip_low, ip_high, global_ierr)
CHKERRQ(global_ierr)

do ispecies = 1, input_nspecies
!  call MatDestroy(particle_shape_x(ispecies), global_ierr)
!  CHKERRQ(global_ierr)
!  call global_matcreate(particle_shape_x(ispecies), &
!    input_nparticle, input_nx, 2, 2)
  call MatZeroEntries(particle_shape_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)

  do ip = ip_low, ip_high - 1
    ! enforce periodic boundary condition
    px(ip - ip_low + 1) = mod(px(ip - ip_low + 1), input_lx)
    ! if x is negative, mod gives negative result, so shift it to positive
    if (px(ip - ip_low + 1) < 0.0_kpr) &
      px(ip - ip_low + 1) = px(ip - ip_low + 1) + input_lx

    sx = px(ip - ip_low + 1) / input_lx * input_nx
    ix1 = floor(sx)
    sx = sx - real(ix1, kpr)
    ix2 = ix1 + 1
    if (ix2 == input_nx) ix2 = 0

    nindex = 2
    indexes(0) = ix1
    indexes(1) = ix2
    values(0) = (1.0_kpr - sx)
    values(1) = sx
    call MatSetValues( &
      particle_shape_x(ispecies), 1, ip, &
      nindex, indexes, values, INSERT_VALUES, global_ierr &
    )
    CHKERRQ(global_ierr)
  end do
  call MatAssemblyBegin(particle_shape_x(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)
  call MatAssemblyEnd(particle_shape_x(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)

!  call VecView(particle_x(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call MatView(particle_shape_xv(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
end do

end subroutine particle_compute_shape_x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute shape matrix in v !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine particle_compute_shape_v
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, nindex
PetscInt :: ip_low, ip_high, ip
PetscInt :: iv
PetscScalar :: sv
PetscScalar, dimension(:), pointer :: pv
PetscInt, dimension(0 : 1) :: indexes
PetscScalar, dimension(0 : 1) :: values

call VecGetOwnershipRange(particle_v(1), ip_low, ip_high, global_ierr)
CHKERRQ(global_ierr)

do ispecies = 1, input_nspecies
!  call MatDestroy(particle_shape_v(ispecies), global_ierr)
!  CHKERRQ(global_ierr)
!  call global_matcreate(particle_shape_v(ispecies), &
!    input_nparticle, input_nx, 2, 2)
  call MatZeroEntries(particle_shape_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)

  do ip = ip_low, ip_high - 1
    ! if particle velocity out of output_v_max, ignore this particle
    if ( &
      pv(ip - ip_low + 1) <= -input_output_v_max &
      .or. pv(ip - ip_low + 1) >= input_output_v_max &
    ) cycle

    sv = (pv(ip - ip_low + 1) + input_output_v_max) &
      / (input_output_v_max * 2.0_kpr) * (input_output_nv - 1)
    iv = floor(sv)
    sv = sv - real(iv, kpr)

    nindex = 2
    indexes(0) = iv
    indexes(1) = iv + 1
    values(0) = (1.0_kpr - sv)
    values(1) = sv
    call MatSetValues( &
      particle_shape_v(ispecies), 1, ip, &
      nindex, indexes, values, INSERT_VALUES, global_ierr &
    )
    CHKERRQ(global_ierr)
  end do
  call MatAssemblyBegin(particle_shape_v(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)
  call MatAssemblyEnd(particle_shape_v(ispecies), &
    MAT_FINAL_ASSEMBLY, global_ierr)
  CHKERRQ(global_ierr)

  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)

!  call VecView(particle_x(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(particle_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
!  call MatView(particle_shape_v(ispecies), PETSC_VIEWER_STDOUT_WORLD, global_ierr)
!  CHKERRQ(global_ierr)
end do

end subroutine particle_compute_shape_v


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

  call VecDestroy(particle_x_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_v_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecDestroy(particle_w_bak(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  call MatDestroy(particle_shape_xv(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call MatDestroy(particle_shape_x(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call MatDestroy(particle_shape_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)
end do

call VecDestroy(particle_tmp1, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(particle_tmp2, global_ierr)
CHKERRQ(global_ierr)

end subroutine particle_final

end module pic1dp_particle

