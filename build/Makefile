# Copyright 2012 Wenjun Deng <wdeng@wdeng.info>
#
# This file is part of PIC1D-PETSc
#
# PIC1D-PETSc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIC1D-PETSc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PIC1D-PETSc.  If not, see <http://www.gnu.org/licenses/>.


# Makefile for compiling PIC1D-PETSc
#
# MPI Fortran 95 compiler
ifndef MPIF90
	MPIF90 := mpif90
endif

# compiling options
ifndef FFLAGS
#	FFLAGS := -g -fbounds-check -ffpe-trap=invalid,zero,overflow # for debugging
	FFLAGS := -O3
endif

# if your Fortran compiler does not support procedure pointer well,
#   uncomment the following line
# NO_PROC_POINTER := -DNO_PROC_POINTER

# source file path
SRCDIR := ../src

.PHONY : build clean

build : pic1dp

# PETSc variables
-include $(PETSC_DIR)/conf/variables

pic1dp : pic1dp.o multirand.o wtimer.o pic1dp_global.o pic1dp_input.o \
 pic1dp_particle.o pic1dp_field.o pic1dp_interaction.o pic1dp_output.o
	$(MPIF90) $(FFLAGS) -o $@ $^ $(PETSC_LIB)

pic1dp.o : $(SRCDIR)/pic1dp.F90 multirand.o wtimer.o pic1dp_global.o \
 pic1dp_input.o pic1dp_particle.o pic1dp_field.o pic1dp_interaction.o \
 pic1dp_output.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

multirand.o : $(SRCDIR)/multirand.F90
	$(MPIF90) $(FFLAGS) $(NO_PROC_POINTER) -c -o $@ $< $(FCPPFLAGS)

wtimer.o : $(SRCDIR)/wtimer.F90
	$(MPIF90) $(FFLAGS) -D__PETSc -c -o $@ $< $(FCPPFLAGS)

pic1dp_global.o : $(SRCDIR)/pic1dp_global.F90
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

pic1dp_input.o : $(SRCDIR)/pic1dp_input.F90 pic1dp_global.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

pic1dp_particle.o : $(SRCDIR)/pic1dp_particle.F90 \
 multirand.o wtimer.o pic1dp_global.o pic1dp_input.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

pic1dp_field.o : $(SRCDIR)/pic1dp_field.F90 \
 wtimer.o pic1dp_global.o pic1dp_input.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

pic1dp_interaction.o : $(SRCDIR)/pic1dp_interaction.F90 wtimer.o \
 pic1dp_global.o pic1dp_input.o pic1dp_field.o pic1dp_particle.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

pic1dp_output.o : $(SRCDIR)/pic1dp_output.F90 \
 wtimer.o pic1dp_global.o pic1dp_input.o pic1dp_field.o \
 pic1dp_particle.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

clean :
	rm -f pic1dp *.o *.mod

