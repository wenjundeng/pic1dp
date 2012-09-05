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


! module for managing output
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
  input_nx * input_nv, global_ierr)
CHKERRQ(global_ierr)
call VecSetFromOptions(output_vec_ptcldist_xv1, global_ierr)
CHKERRQ(global_ierr)
call VecDuplicate(output_vec_ptcldist_xv1, &
  output_vec_ptcldist_xv2, global_ierr)
CHKERRQ(global_ierr)

call VecCreate(MPI_COMM_WORLD, output_vec_ptcldist_v, global_ierr)
CHKERRQ(global_ierr)
call VecSetSizes(output_vec_ptcldist_v, PETSC_DECIDE, &
  input_nv, global_ierr)
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
intbuf(4) = input_nv
do iparameter = 5, 4 + input_nmode
  intbuf(iparameter) = input_modes(iparameter - 5)
end do
call PetscViewerBinaryWriteInt(output_viewer, intbuf, nintparameter, &
  PETSC_TRUE, global_ierr)
CHKERRQ(global_ierr)

realbuf(1) = input_lx
realbuf(2) = input_v_max
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
! electric energe + marker, total and perturbed kinetic energe for each species
PetscInt, parameter :: nscalar = 2 + input_nspecies * 3

PetscInt :: ispecies
PetscReal :: energe
PetscReal, dimension(nscalar) :: realbuf

realbuf(1) = global_time

! output E^2 (electric energe) and particle kinetic energe
call VecNorm(field_electric, NORM_2, energe, global_ierr)
CHKERRQ(global_ierr)
energe = energe * energe * input_lx / input_nx
realbuf(2) = energe
electric_energe = energe

do ispecies = 1, input_nspecies
  ! put v*v in particle_tmp1
  call VecPointwiseMult(particle_tmp1, &
    particle_v(ispecies), particle_v(ispecies), global_ierr)
  CHKERRQ(global_ierr)

  ! calculate sum(v*v) to get marker energe
  call VecSum(particle_tmp1, energe, global_ierr)
  CHKERRQ(global_ierr)
  realbuf(ispecies * 3) = energe

  ! calculate sum(v*v*p) to get total energe
  call VecPointwiseMult(particle_tmp2, &
    particle_tmp1, particle_p(ispecies), global_ierr)
  CHKERRQ(global_ierr)
  call VecSum(particle_tmp2, energe, global_ierr)
  CHKERRQ(global_ierr)
  realbuf(ispecies * 3 + 1) = energe

  ! calculate sum(v*v*w) to get perturbed energe
  if (input_deltaf == 1) then
    call VecPointwiseMult(particle_tmp2, &
      particle_tmp1, particle_w(ispecies), global_ierr)
    CHKERRQ(global_ierr)
    call VecSum(particle_tmp2, energe, global_ierr)
    CHKERRQ(global_ierr)
    if (input_linear == 1) then
      ! linear case, add perturbed energe to get total energe
      realbuf(ispecies * 3 + 1) = realbuf(ispecies * 3 + 1) + energe
    end if
  else
    ! at this point energe is total energe
    ! subtract equilibrium energe to get perturbed energe
    if (input_iptcldist == 1) then ! two-stream1
      energe = energe - 3.0_kpr * input_species_density(ispecies) * input_lx
    elseif (input_iptcldist == 2) then ! two-stream2
      ! need to calculate equilibrium energe for this case
    else ! (shifted) Maxwellian
      energe = energe - input_species_temperature(ispecies) &
        / input_species_mass(ispecies) * input_species_density(ispecies) &
        * input_lx
    end if
  end if
  realbuf(ispecies * 3 + 2) = energe
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
call VecView(field_chargeden, output_viewer, global_ierr)
CHKERRQ(global_ierr)

end subroutine output_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output physical particle distribution function f and delta f !
! on x-v plane and in v space                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_ptcldist
use pic1dp_global
use pic1dp_input
use pic1dp_particle
implicit none
#include "finclude/petsc.h90"

PetscScalar, parameter :: delv_inv &
  = (input_nv - 1) / (2.0_kpr * input_v_max)
PetscScalar, parameter :: delx_inv = input_nx / input_lx

PetscInt :: ispecies, ip, ix, iv
PetscScalar :: sx, sv
PetscScalar, dimension(:), pointer :: px, pv, pp, pw, ps
PetscScalar, dimension(0 : input_nx * input_nv - 1) :: &
  ptcldist_markr_xv, ptcldist_total_xv, ptcldist_pertb_xv, &
  ptcldist_markr_xv_redu, ptcldist_total_xv_redu, ptcldist_pertb_xv_redu
PetscScalar, dimension(0 : input_nv - 1) :: &
  ptcldist_markr_v, ptcldist_total_v, ptcldist_pertb_v, &
  ptcldist_markr_v_redu, ptcldist_total_v_redu, ptcldist_pertb_v_redu

! calculate shape matrix in x-v plane (too slow, no longer used)
!call particle_compute_shape_xv

do ispecies = 1, input_nspecies
  ptcldist_markr_xv = 0.0_kpr
  ptcldist_total_xv = 0.0_kpr
  ptcldist_pertb_xv = 0.0_kpr
  ptcldist_markr_v = 0.0_kpr
  ptcldist_total_v = 0.0_kpr
  ptcldist_pertb_v = 0.0_kpr

  call VecGetArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecGetArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)
  if (input_deltaf == 1) then
    call VecGetArrayF90(particle_w(ispecies), pw, global_ierr)
    CHKERRQ(global_ierr)
  end if

  do ip = 1, particle_ip_high - particle_ip_low
    ! ignore invalid particle
    if (ps(ip) < -0.5_kpr) cycle
    ! if particle speed is out of v_max, ignore this particle
    if (abs(pv(ip)) >= input_v_max) cycle

    sx = px(ip) / input_lx * input_nx
    ix = floor(sx)
    sx = 1.0_kpr - (sx - real(ix, kpr))

    sv = (pv(ip) + input_v_max) &
      / (input_v_max * 2.0_kpr) * (input_nv - 1)
    iv = floor(sv)
    sv = 1.0_kpr - (sv - real(iv, kpr))

    ptcldist_markr_xv(iv * input_nx + ix) &
      = ptcldist_markr_xv(iv * input_nx + ix) &
      + sx * sv
    ptcldist_total_xv(iv * input_nx + ix) &
      = ptcldist_total_xv(iv * input_nx + ix) &
      + sx * sv * pp(ip)
    if (input_deltaf == 1) then
      ptcldist_pertb_xv(iv * input_nx + ix) &
        = ptcldist_pertb_xv(iv * input_nx + ix) &
        + sx * sv * pw(ip)
    end if
    ptcldist_markr_xv((iv + 1) * input_nx + ix) &
      = ptcldist_markr_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv)
    ptcldist_total_xv((iv + 1) * input_nx + ix) &
      = ptcldist_total_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv) * pp(ip)
    if (input_deltaf == 1) then
      ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
        = ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
        + sx * (1.0_kpr - sv) * pw(ip)
    end if
    ix = ix + 1
    if (ix > input_nx - 1) ix = 0
    sx = 1.0_kpr - sx
    ptcldist_markr_xv(iv * input_nx + ix) &
      = ptcldist_markr_xv(iv * input_nx + ix) &
      + sx * sv
    ptcldist_total_xv(iv * input_nx + ix) &
      = ptcldist_total_xv(iv * input_nx + ix) &
      + sx * sv * pp(ip)
    if (input_deltaf == 1) then
      ptcldist_pertb_xv(iv * input_nx + ix) &
        = ptcldist_pertb_xv(iv * input_nx + ix) &
        + sx * sv * pw(ip)
    end if
    ptcldist_markr_xv((iv + 1) * input_nx + ix) &
      = ptcldist_markr_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv)
    ptcldist_total_xv((iv + 1) * input_nx + ix) &
      = ptcldist_total_xv((iv + 1) * input_nx + ix) &
      + sx * (1.0_kpr - sv) * pp(ip)
    if (input_deltaf == 1) then
      ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
        = ptcldist_pertb_xv((iv + 1) * input_nx + ix) &
        + sx * (1.0_kpr - sv) * pw(ip)
    end if

    ptcldist_markr_v(iv) = ptcldist_markr_v(iv) + sv
    ptcldist_total_v(iv) = ptcldist_total_v(iv) &
      + sv * pp(ip)
    if (input_deltaf == 1) then
      ptcldist_pertb_v(iv) = ptcldist_pertb_v(iv) &
        + sv * pw(ip)
    end if
    ptcldist_markr_v(iv + 1) = ptcldist_markr_v(iv + 1) &
      + (1.0_kpr - sv)
    ptcldist_total_v(iv + 1) = ptcldist_total_v(iv + 1) &
      + (1.0_kpr - sv) * pp(ip)
    if (input_deltaf == 1) then
      ptcldist_pertb_v(iv + 1) = ptcldist_pertb_v(iv + 1) &
        + (1.0_kpr - sv) * pw(ip)
    end if
  end do ! ip = 1, particle_ip_high - particle_ip_low
  call VecRestoreArrayF90(particle_x(ispecies), px, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_v(ispecies), pv, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_p(ispecies), pp, global_ierr)
  CHKERRQ(global_ierr)
  call VecRestoreArrayF90(particle_s(ispecies), ps, global_ierr)
  CHKERRQ(global_ierr)
  if (input_deltaf == 1) then
    call VecRestoreArrayF90(particle_w(ispecies), pw, global_ierr)
    CHKERRQ(global_ierr)
  end if

  ! if linear, p = f_0 / g, needs to add perturbed dist. to get total dist.
  if (input_linear == 1) then
    ptcldist_total_xv(:) = ptcldist_total_xv(:) + ptcldist_pertb_xv(:)
    ptcldist_total_v(:) = ptcldist_total_v(:) + ptcldist_pertb_v(:)
  end if

  call MPI_Reduce(ptcldist_markr_xv, ptcldist_markr_xv_redu, &
    input_nx * input_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  call MPI_Reduce(ptcldist_total_xv, ptcldist_total_xv_redu, &
    input_nx * input_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  call MPI_Reduce(ptcldist_markr_v, ptcldist_markr_v_redu, &
    input_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  call MPI_Reduce(ptcldist_total_v, ptcldist_total_v_redu, &
    input_nv, MPIU_SCALAR, MPI_SUM, 0, &
    MPI_COMM_WORLD, global_ierr)
  CHKERRQ(global_ierr)
  if (input_deltaf == 1) then
    call MPI_Reduce(ptcldist_pertb_xv, ptcldist_pertb_xv_redu, &
      input_nx * input_nv, MPIU_SCALAR, MPI_SUM, 0, &
      MPI_COMM_WORLD, global_ierr)
    CHKERRQ(global_ierr)
    call MPI_Reduce(ptcldist_pertb_v, ptcldist_pertb_v_redu, &
      input_nv, MPIU_SCALAR, MPI_SUM, 0, &
      MPI_COMM_WORLD, global_ierr)
    CHKERRQ(global_ierr)
  end if

  ! divide by grid size to get actual distribution function
  if (global_mype == 0) then ! only root MPI process needs to do this
    ptcldist_markr_xv_redu(:) = ptcldist_markr_xv_redu(:) * delx_inv * delv_inv
    ptcldist_total_xv_redu(:) = ptcldist_total_xv_redu(:) * delx_inv * delv_inv
    ptcldist_markr_v_redu(:) = ptcldist_markr_v_redu(:) * delv_inv
    ptcldist_total_v_redu(:) = ptcldist_total_v_redu(:) * delv_inv
    if (input_deltaf == 1) then
      ptcldist_pertb_xv_redu(:) = ptcldist_pertb_xv_redu(:) &
        * delx_inv * delv_inv
      ptcldist_pertb_v_redu(:) = ptcldist_pertb_v_redu(:) * delv_inv
    else
      do iv = 0, input_nv - 1
        ! reuse sv for velocity value
        sv = (real(iv, kpr) / (input_nv - 1) * 2.0_kpr - 1.0_kpr) &
          * input_v_max
        if (input_iptcldist == 1) then ! two-stream1
          ptcldist_pertb_xv_redu(iv * input_nx : (iv + 1) * input_nx - 1) &
            = ptcldist_total_xv_redu(iv * input_nx : (iv + 1) * input_nx - 1) &
            - input_species_density(ispecies) * sv**2 * exp(-sv**2 / 2.0_kpr) &
            / sqrt(2.0_kpr * PETSC_PI)
          ptcldist_pertb_v_redu(iv) = ptcldist_total_v_redu(iv) - input_lx &
            * input_species_density(ispecies) * sv**2 * exp(-sv**2 / 2.0_kpr) &
            / sqrt(2.0_kpr * PETSC_PI)
        elseif (input_iptcldist == 2) then ! two-stream2
          ptcldist_pertb_xv_redu(iv * input_nx : (iv + 1) * input_nx - 1) &
            = ptcldist_total_xv_redu(iv * input_nx : (iv + 1) * input_nx - 1) &
            - input_species_density(ispecies) &
            * (exp(-(sv + input_species_v0(ispecies))**2 / (2.0_kpr &
            * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
            + exp(-(sv - input_species_v0(ispecies))**2 / (2.0_kpr &
            * input_species_temperature(ispecies) / input_species_mass(ispecies)))) &
            / (sqrt(8.0_kpr * PETSC_PI) * input_species_temperature(ispecies) &
            / input_species_mass(ispecies))
          ptcldist_pertb_v_redu(iv) = ptcldist_total_v_redu(iv) - input_lx &
            * input_species_density(ispecies) &
            * (exp(-(sv + input_species_v0(ispecies))**2 / (2.0_kpr &
            * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
            + exp(-(sv - input_species_v0(ispecies))**2 / (2.0_kpr &
            * input_species_temperature(ispecies) / input_species_mass(ispecies)))) &
            / (sqrt(8.0_kpr * PETSC_PI) * input_species_temperature(ispecies) &
            / input_species_mass(ispecies))
        else ! (shifted) Maxwellian
          ptcldist_pertb_xv_redu(iv * input_nx : (iv + 1) * input_nx - 1) &
            = ptcldist_total_xv_redu(iv * input_nx : (iv + 1) * input_nx - 1) &
            - input_species_density(ispecies) &
            * exp(-(sv - input_species_v0(ispecies))**2 / (2.0_kpr &
            * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
            / (sqrt(2.0_kpr * PETSC_PI) * input_species_temperature(ispecies) &
            / input_species_mass(ispecies))
          ptcldist_pertb_v_redu(iv) = ptcldist_total_v_redu(iv) - input_lx &
            * input_species_density(ispecies) &
            * exp(-(sv - input_species_v0(ispecies))**2 / (2.0_kpr &
            * input_species_temperature(ispecies) / input_species_mass(ispecies))) &
            / (sqrt(2.0_kpr * PETSC_PI) * input_species_temperature(ispecies) &
            / input_species_mass(ispecies))
        end if
      end do ! iv = 0, input_nv - 1
    end if ! else of if (input_deltaf == 1)
  end if ! (global_mype == 0)

  ! only root MPI process will write to file
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_markr_xv_redu, &
    input_nx * input_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_total_xv_redu, &
    input_nx * input_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_pertb_xv_redu, &
    input_nx * input_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_markr_v_redu, &
    input_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_total_v_redu, &
    input_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
  call PetscViewerBinaryWriteScalar(output_viewer, ptcldist_pertb_v_redu, &
    input_nv, PETSC_TRUE, global_ierr)
  CHKERRQ(global_ierr)
end do ! ispecies = 1, input_nspecies

end subroutine output_ptcldist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output progress information to stdout !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_progress(progress_type, electric_energe)
use pic1dp_global
use pic1dp_input
use pic1dp_particle
implicit none
#include "finclude/petsc.h90"

! progress type. 0: regular; 1: merge particle output
PetscInt, intent(in) :: progress_type

! electric energe, has to be provided if progress_type = 0
PetscReal, intent(in), optional :: electric_energe

PetscReal, dimension(2) :: progress ! 1: itime, 2: time
PetscInt :: iprogress
character :: cprogress

if (input_verbosity == 1) then
  progress(1) = 1e2_kpr * real(global_itime, kpr) / input_ntime_max
  progress(2) = 1e2_kpr * global_time / input_time_max
  iprogress = maxloc(progress, 1)
  if (iprogress == 1) then
    cprogress = 'i'
  else
    cprogress = 't'
  end if
  if (progress_type == 1) then
    write (global_msg, '(a, f5.1, a, i7, f9.3, a, i9, a, i9, a)') &
      cprogress, progress(iprogress), "%%", global_itime, &
      global_time + input_dt, &
      ' : merged ', sum(particle_nredu), '; left: ', &
      input_nparticle * input_nspecies - sum(particle_nredu), "\n"
  else ! regular progress type
    write (global_msg, '(a, f5.1, a, i7, f9.3, es12.3e3, a)') &
      cprogress, progress(iprogress), "%%", global_itime, global_time, &
      electric_energe, "\n"
  end if ! else of if (progress_type == 1)
  call global_pp(global_msg)
elseif (input_verbosity >= 2) then
  if (progress_type == 1) then
    write (global_msg, '(a, i9, a, i9, a)') &
      'Info: particle_merge reduced # of particles:', &
      sum(particle_nredu), '; left:', &
      input_nparticle * input_nspecies - sum(particle_nredu), "\n"
  else ! regular progress type
    write (global_msg, '(a, i7, a, f9.3, a)') 'Info: finished itime = ', &
      global_itime, ', time = ', global_time, "\n"
  end if ! else of if (progress_type == 1)
  call global_pp(global_msg)
end if

end subroutine output_progress


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output all needed quantities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_all
use pic1dp_global
implicit none
#include "finclude/petsc.h90"

PetscReal :: electric_energe

call output_field(electric_energe)
call output_ptcldist
call output_progress(0, electric_energe)

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

