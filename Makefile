# MPI Fortran compiler
MPIF90 := mpif90

# compiling options
FFLAGS := -O2

# MPI executor
MPIEXEC := mpiexec

# number of processes to build (compile) PIC1D-PETSc
NPE_BUILD := 4

# number of MPI processes to run PIC1D-PETSc
NPE_RUN := 20

export

.PHONY : build run restart cleanbuild cleanrun

build : build/pic1dp

build/pic1dp : $(wildcard src/*)
	${MAKE} -C ./build -j ${NPE_BUILD} build

run : build/pic1dp
	cp build/pic1dp run/pic1dp
	${MAKE} -C ./run run

restart : build/pic1dp
	cp build/pic1dp run/pic1dp
	${MAKE} -C ./run restart

cleanbuild :
	${MAKE} -C ./build clean

cleanrun :
	${MAKE} -C ./run clean

