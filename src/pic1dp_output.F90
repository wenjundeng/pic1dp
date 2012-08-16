module pic1dp_output
implicit none
#include "finclude/petscdef.h"

! for particle distribution
Vec :: output_vec_ptcldist_xv, output_vec_ptcldist_v

PetscViewer :: output_viewer ! for output file writing

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_init
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt, parameter :: &
  nintparameter = 4 + input_nmode, &
  nrealparameter = 2
PetscInt, dimension(nintparameter) :: intbuf
PetscReal, dimension(nrealparameter) :: realbuf

PetscInt :: iparameter

call VecCreate(MPI_COMM_WORLD, output_vec_ptcldist_xv, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(output_vec_ptcldist_xv, PETSC_DECIDE, &
  input_nx * input_output_nv, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(output_vec_ptcldist_xv, global_ierr)
CHKERRQ(global_ierr)

call VecCreate(MPI_COMM_WORLD, output_vec_ptcldist_v, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(output_vec_ptcldist_v, PETSC_DECIDE, &
  input_output_nv, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(output_vec_ptcldist_v, global_ierr)
CHKERRQ(global_ierr)

call PetscViewerBinaryOpen( &
  MPI_COMM_WORLD, 'pic1dp.out', FILE_MODE_WRITE, &
  output_viewer, global_ierr &
)
CHKERRQ(global_ierr)

! output parameters
intbuf(1) = input_nspecies
intbuf(2) = input_nmode
intbuf(3) = input_nx
intbuf(4) = input_output_nv
do iparameter = 5, 4 + input_nmode
  intbuf(iparameter) = input_mode(iparameter - 5)
end do
call PetscViewerBinaryWriteInt(output_viewer, intbuf, nintparameter, &
  PETSC_TRUE, global_ierr)
CHKERRQ(global_ierr)

realbuf(1) = input_lx
realbuf(2) = input_output_v_max
call PetscViewerBinaryWriteReal(output_viewer, realbuf, nrealparameter, &
  PETSC_TRUE, global_ierr)
CHKERRQ(global_ierr)

end subroutine output_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output field quantities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_field(electric_energe)
use pic1dp_global
use pic1dp_field
use pic1dp_particle
implicit none
#include "finclude/petsc.h90"

PetscReal, intent(out) :: electric_energe

! # of scalars for output
! electric energe + total and perturbed kinetic energe for each species
PetscInt, parameter :: nscalar = 1 + input_nspecies * 2

PetscInt :: ispecies
PetscReal :: energe
PetscReal, dimension(nscalar) :: realbuf

! output E^2 (electric energe) and particle kinetic energe
call VecNorm(field_electric, NORM_2, energe, global_ierr)
CHKERRQ(global_ierr)
energe = energe * energe
realbuf(1) = energe
electric_energe = energe

do ispecies = 1, input_nspecies
  ! put v*v in particle_tmp1
  call VecPointwiseMult(particle_tmp1, &
    particle_v(ispecies), particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  ! calculate sum(v*v*p) to get total energe
  call VecPointwiseMult(particle_tmp2, &
    particle_tmp1, particle_p(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecSum(particle_tmp2, energe, global_ierr)
  CHKERRQ(global_ierr)
  realbuf(ispecies * 2) = energe

  ! calculate sum(v*v*w) to get perturbed energe
  call VecPointwiseMult(particle_tmp2, &
    particle_tmp1, particle_w(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecSum(particle_tmp2, energe, global_ierr)
  CHKERRQ(global_ierr)
  realbuf(ispecies * 2 + 1) = energe
end do
call PetscViewerBinaryWriteReal(output_viewer, realbuf, nscalar, &
  PETSC_TRUE, global_ierr)
CHKERRQ(global_ierr)

! output electric field Fourier components
call VecView(field_mode_re, output_viewer, global_ierr)
CHKERRQ(global_ierr)
call VecView(field_mode_im, output_viewer, global_ierr)
CHKERRQ(global_ierr)

! output electric field and charge in x space
call VecView(field_electric, output_viewer, global_ierr)
CHKERRQ(global_ierr)
call VecView(field_charge, output_viewer, global_ierr)
CHKERRQ(global_ierr)

end subroutine output_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output physical particle distribution function f and delta f on x-v plane !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_ptcldist_xv
use pic1dp_global
use pic1dp_particle
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

! calculate shape matrix in x-v plane
call particle_compute_shape_xv

do ispecies = 1, input_nspecies
  ! calculate and output f
  call MatMultTranspose(particle_shape_xv, particle_p(ispecies), &
    output_vec_ptcldist_xv, global_ierr)
  CHKERRQ(global_ierr)
  call VecView(output_vec_ptcldist_xv, output_viewer, global_ierr)
  CHKERRQ(global_ierr)

  ! calculate and output delta f
  call MatMultTranspose(particle_shape_xv, particle_w(ispecies), &
    output_vec_ptcldist_xv, global_ierr)
  CHKERRQ(global_ierr)
  call VecView(output_vec_ptcldist_xv, output_viewer, global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine output_ptcldist_xv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output physical particle distribution function f and delta f on v space !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_ptcldist_v
use pic1dp_global
use pic1dp_particle
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies

! calculate shape matrix in v space
call particle_compute_shape_v

do ispecies = 1, input_nspecies
  ! calculate and output f
  call MatMultTranspose(particle_shape_v, particle_p(ispecies), &
    output_vec_ptcldist_v, global_ierr)
  CHKERRQ(global_ierr)
  call VecView(output_vec_ptcldist_v, output_viewer, global_ierr)
  CHKERRQ(global_ierr)

  ! calculate and output delta f
  call MatMultTranspose(particle_shape_v, particle_w(ispecies), &
    output_vec_ptcldist_v, global_ierr)
  CHKERRQ(global_ierr)
  call VecView(output_vec_ptcldist_v, output_viewer, global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine output_ptcldist_v


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output progress information to stdout !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_progress(itime, time, electric_energe)
use pic1dp_global
use pic1dp_input
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(in) :: itime
PetscReal, intent(in) :: time, electric_energe

PetscReal, dimension(2) :: progress ! 1: itime, 2: time
PetscInt :: iprogress
character :: cprogress

if (input_verbosity == 1) then
  progress(1) = 1e2_kpr * real(itime, kpr) / input_ntime_max
  progress(2) = 1e2_kpr * time / input_time_max
  iprogress = maxloc(progress, 1)
  if (iprogress == 1) then
    cprogress = 'i'
  else
    cprogress = 't'
  end if
  write (global_msg, '(a, f5.1, a, i7, f9.3, es11.2e3, a)') &
    cprogress, progress(iprogress), "%%", itime, time, electric_energe, "\n"
  call global_pp(global_msg)
elseif (input_verbosity >= 2) then
  write (global_msg, '(a, i7, a, f9.3, a)') 'Info: finished itime = ', itime, &
    ', time = ', time, "\n"
  call global_pp(global_msg)
end if
end subroutine output_progress


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output all needed quantities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_all(itime, time)
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(in) :: itime
PetscReal, intent(in) :: time

PetscReal :: electric_energe

call output_field(electric_energe)
!call output_ptcldist_xv
!call output_ptcldist_v
call output_progress(itime, time, electric_energe)

end subroutine output_all


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

call VecDestroy(output_vec_ptcldist_xv, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(output_vec_ptcldist_v, global_ierr)
CHKERRQ(global_ierr)
call PetscViewerDestroy(output_viewer, global_ierr)
CHKERRQ(global_ierr)

end subroutine output_final

end module pic1dp_output

