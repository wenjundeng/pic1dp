# MPI Fortran 90 compiler
ifndef MPIF90
	MPIF90 := mpif90
endif

# compiling options
ifndef FFLAGS
#	FFLAGS := -g -fbounds-check
	FFLAGS := -O2
endif

SRCDIR := ../src

.PHONY : build clean

build : pic1dp

-include ${PETSC_DIR}/conf/variables

pic1dp : pic1dp.o pic1dp_global.o pic1dp_input.o pic1dp_particle.o \
 pic1dp_field.o pic1dp_interaction.o pic1dp_output.o wtimer.o gaussian.o
	${MPIF90} ${FFLAGS} -o $@ $^ ${PETSC_LIB}

pic1dp.o : ${SRCDIR}/pic1dp.F90 pic1dp_global.o pic1dp_input.o \
 pic1dp_particle.o pic1dp_field.o pic1dp_interaction.o pic1dp_output.o \
 wtimer.o gaussian.o
	${MPIF90} ${FFLAGS} -c -o $@ $< ${FCPPFLAGS}

pic1dp_global.o : ${SRCDIR}/pic1dp_global.F90
	${MPIF90} ${FFLAGS} -c -o $@ $< ${FCPPFLAGS}

pic1dp_input.o : ${SRCDIR}/pic1dp_input.F90 pic1dp_global.o
	${MPIF90} ${FFLAGS} -c -o $@ $< ${FCPPFLAGS}

pic1dp_particle.o : ${SRCDIR}/pic1dp_particle.F90 pic1dp_global.o \
 pic1dp_input.o gaussian.o
	${MPIF90} ${FFLAGS} -c -o $@ $< ${FCPPFLAGS}

pic1dp_field.o : ${SRCDIR}/pic1dp_field.F90 pic1dp_global.o pic1dp_input.o
	${MPIF90} ${FFLAGS} -c -o $@ $< ${FCPPFLAGS}

pic1dp_interaction.o : ${SRCDIR}/pic1dp_interaction.F90 pic1dp_global.o \
 pic1dp_input.o pic1dp_field.o pic1dp_particle.o
	${MPIF90} ${FFLAGS} -c -o $@ $< ${FCPPFLAGS}

pic1dp_output.o : ${SRCDIR}/pic1dp_output.F90 pic1dp_global.o pic1dp_input.o \
 pic1dp_field.o pic1dp_particle.o
	${MPIF90} ${FFLAGS} -c -o $@ $< ${FCPPFLAGS}

wtimer.o : ${SRCDIR}/wtimer.F90
	${MPIF90} ${FFLAGS} -D__PETSc -c -o $@ $< ${FCPPFLAGS}

gaussian.o : ${SRCDIR}/gaussian.F90
	${MPIF90} ${FFLAGS} -D__PETSc -c -o $@ $< ${FCPPFLAGS}

clean :
	rm -f pic1dp *.o *.mod

