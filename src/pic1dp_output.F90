module pic1dp_output
implicit none
#include "finclude/petscdef.h"

! for particle distribution
Vec :: output_vec_ptcldist_xv1, output_vec_ptcldist_xv2
Vec :: output_vec_ptcldist_v

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

call VecCreate(MPI_COMM_WORLD, output_vec_ptcldist_xv1, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(output_vec_ptcldist_xv1, PETSC_DECIDE, &
  input_nx * input_output_nv, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(output_vec_ptcldist_xv1, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(output_vec_ptcldist_xv1, &
  output_vec_ptcldist_xv2, global_ierr)
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output physical particle distribution function f and delta f !
! on x-v plane and in v space                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_ptcldist
use pic1dp_global
use pic1dp_particle
implicit none
#include "finclude/petsc.h90"

PetscInt :: ispecies, ip, ix, iv
PetscScalar :: sx, sv
PetscScalar, dimension(:), pointer :: px, pv, pp, pw
PetscScalar, dimension(0 : input_nx * input_output_nv - 1) :: &
  ptcldist_total_xv, ptcldist_pertb_xv, &
  ptcldist_total_xv_redu, ptcldist_pertb_xv_redu
PetscScalar, dimension(0 : input_output_nv - 1) :: &
  ptcldist_total_v, ptcldist_pertb_v, &
  ptcldist_total_v_redu, ptcldist_pertb_v_redu

! calculate shape matrix in x-v plane (too slow, no longer used)
!call particle_compute_shape_xv

do ispecies = 1, input_nspecies
  ! calculate and output f
!  call MatMultTranspose(particle_shape_xv, particle_p(ispecies), &
!    output_vec_ptcldist_xv, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(output_vec_ptcldist_xv1, output_viewer, global_ierr)
!  CHKERRQ(global_ierr)

  ! calculate and output delta f
!  call MatMultTranspose(particle_shape_xv, particle_w(ispecies), &
!    output_vec_ptcldist_xv, global_ierr)
!  CHKERRQ(global_ierr)
!  call VecView(output_vec_ptcldist_xv2, output_viewer, global_ierr)
!  CHKERRQ(global_ierr)


  ptcldist_total_xv = 0.0_kpr
  ptcldist_pertb_xv = 0.0_kpr
  ptcldist_total_v = 0.0_kpr
  ptcldist_pertb_v = 0.0_kpr

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  do ip = particle_ip_low, particle_ip_high
    if (pv(ip - particle_ip_low + 1) <= -input_output_v_max &
      .or. pv(ip - particle_ip_low + 1) >= input_output_v_max &
    ) cycle

    sx = px(ip - particle_ip_low + 1) / input_lx * input_nx
    ix = floor(sx)
    sx = 1.0_kpr - (sx - real(ix, kpr))

    sv = (pv(ip - particle_ip_low + 1) + input_output_v_max) &
      / (input_output_v_max * 2.0_kpr) * (input_output_nv - 1)
    iv = floor(sv)
    sv = 1.0_kpr - (sv - real(iv, kpr))

    ptcldist_total_xv(iv * input_nx + ix) &
      = ptcldist_total_xv(iv * input_nx + ix) &
      + sx * sv * pp(ip - particle_ip_low + 1)
    ptcldist_pertb_xv(iv * input_nx + ix) &
      = ptcldist_pertb_xv(iv * input_nx + ix) &
      + sx * sv * pw(ip - particle_ip_low + 1)
    ptcldist_total_xv((iv + 1) * input_nx + ix) &
      = ptcldist_total_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv) * pp(ip - particle_ip_low + 1)
    ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
      = ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv) * pw(ip - particle_ip_low + 1)
    ix = ix + 1
    if (ix > input_nx - 1) ix = 0
    sx = 1.0_kpr - sx
    ptcldist_total_xv(iv * input_nx + ix) &
      = ptcldist_total_xv(iv * input_nx + ix) &
      + sx * sv * pp(ip - particle_ip_low + 1)
    ptcldist_pertb_xv(iv * input_nx + ix) &
      = ptcldist_pertb_xv(iv * input_nx + ix) &
      + sx * sv * pw(ip - particle_ip_low + 1)
    ptcldist_total_xv((iv + 1) * input_nx + ix) &
      = ptcldist_total_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv) * pp(ip - particle_ip_low + 1)
    ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
      = ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv) * pw(ip - particle_ip_low + 1)

    ptcldist_total_v(iv) = ptcldist_total_v(iv) &
      + sv * pp(ip - particle_ip_low + 1)
    ptcldist_pertb_v(iv) = ptcldist_pertb_v(iv) &
      + sv * pw(ip - particle_ip_low + 1)
    ptcldist_total_v(iv + 1) = ptcldist_total_v(iv + 1) &
      + (1.0_kpr - sv) * pp(ip - particle_ip_low + 1)
    ptcldist_pertb_v(iv + 1) = ptcldist_pertb_v(iv + 1) &
      + (1.0_kpr - sv) * pw(ip - particle_ip_low + 1)
  end do
  
  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
  CHKERRQ(global_ierr)

  call MPI_Reduce(ptcldist_total_xv, ptcldist_total_xv_redu, &
    input_nx * input_output_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  call MPI_Reduce(ptcldist_pertb_xv, ptcldist_pertb_xv_redu, &
    input_nx * input_output_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  call MPI_Reduce(ptcldist_total_v, ptcldist_total_v_redu, &
    input_output_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  call MPI_Reduce(ptcldist_pertb_v, ptcldist_pertb_v_redu, &
    input_output_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)

  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_total_xv_redu, &
    input_nx * input_output_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_pertb_xv_redu, &
    input_nx * input_output_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_total_v_redu, &
    input_output_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_pertb_v_redu, &
    input_output_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
end do

end subroutine output_ptcldist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output physical particle distribution function f and delta f on v space !
! obsolete                                                                !
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
call output_ptcldist
call output_progress(itime, time, electric_energe)

end subroutine output_all


!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalize PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

call VecDestroy(output_vec_ptcldist_xv1, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(output_vec_ptcldist_xv2, global_ierr)
CHKERRQ(global_ierr)
call VecDestroy(output_vec_ptcldist_v, global_ierr)
CHKERRQ(global_ierr)
call PetscViewerDestroy(output_viewer, global_ierr)
CHKERRQ(global_ierr)

end subroutine output_final

end module pic1dp_output

