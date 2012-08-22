PIC1D-PETSc
===========
version 2012-08-16 00:38:26-04:00
---------------------------------

PIC1D-PETSc is a code simulating 1D electrostatic plasma by solving
Vlasov-Poisson equation using particle-in-cell (PIC) method.  PIC1D-PETSc is a
reimplementation of Zhihong Lin's
[PIC1D](http://phoenix.ps.uci.edu/zlin/pic1d/).  This reimplementation
demonstrates the possibility of formulating PIC method in a vector-matrix form
and tests the practicability of implementing the vector-matrix PIC formulation
by PETSc.


Copying
-------

System requirements
-------------------

PIC1D-PETSc requires MPI-2 and [PETSc](http://www.mcs.anl.gov/petsc/) with real
scalars.  PIC1D-PETSc has been tested with
[MPICH2](http://www.mcs.anl.gov/research/projects/mpich2/) 1.4,
[OpenMPI](http://www.open-mpi.org/) 1.5 and PETSc 3.2.

PIC1D-PETSc also requires a Fortran 90 compiler for compilation.  PIC1D-PETSc
has been tested with [GNU Fortran](http://gcc.gnu.org/fortran/) 4.6.


Usage
-----

### Know about the directory tree

+ `README.md` -- This instruction.
+ `Makefile` -- For using GNU make to compile and run PIC1D-PETSc.
+ `build/` -- Place to compile PIC1D-PETSc.  Object files generated during
compilation will be put here.
	+ `build/Makefile` -- For using GNU make to compile PIC1D-PETSc.
+ `run/` -- Place to run PIC1D-PETSc.  Output data will be put here.
	+ `run/Makefile` -- For using GNU make to launch PIC1D-PETSc.
+ `src/` -- Place to store source files.
	+ `src/pic1dp.F90` -- The main program file.


