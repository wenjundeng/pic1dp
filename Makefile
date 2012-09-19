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


# Makefile for compiling and running PIC1D-PETSc

# MPI Fortran 90 compiler
MPIF90 := mpif90

# compiling options
FFLAGS := -O3

# MPI executor
MPIEXEC := mpiexec

# number of processes to build (compile) PIC1D-PETSc
NPE_BUILD := 4

# number of MPI processes to run PIC1D-PETSc
NPE_RUN := 4

export

.PHONY : build run visual doc cleanbuild cleanrun

build : build/pic1dp

build/pic1dp : $(wildcard src/*)
	$(MAKE) -C ./build -j $(NPE_BUILD) build

run : build/pic1dp
	#cp build/pic1dp run/pic1dp
	$(MAKE) -C ./run run

visual :
	./tools/visual.py ./run

doc :
	$(MAKE) -C ./doc

cleanbuild :
	$(MAKE) -C ./build clean

cleanrun :
	$(MAKE) -C ./run clean

