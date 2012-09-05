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


# Makefile for launching PIC1D-PETSc

# MPI executor
ifndef MPIEXEC
	MPIEXEC := mpiexec
endif

# number of processes to build (compile) PIC1D-PETSc
ifndef NPE_BUILD
	NPE_BUILD := 4
endif

# number of MPI processes to run PIC1D-PETSc
ifndef NPE_RUN
	NPE_RUN := 4
endif

.PHONY : run clean

run : pic1dp clean
	$(MPIEXEC) -n $(NPE_RUN) ./$< 2>&1 | tee $<.log

pic1dp : ../build/pic1dp
	cp ../build/pic1dp ./

../build/pic1dp : $(wildcard ../src/*) ../build/Makefile
	$(MAKE) -C ../build -j $(NPE_BUILD) build

clean :
	rm -f pic1dp.log pic1dp.out pic1dp.out.info

