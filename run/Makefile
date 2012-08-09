# MPI executor
ifndef MPIEXEC
	MPIEXEC := mpiexec
endif

# number of processes to build (compile) PIC1D-PETSc
ifndef NPE_BUILD
	NPE_BUILD := 2
endif

# number of MPI processes to run PIC1D-PETSc
ifndef NPE_RUN
	NPE_RUN := 2
endif

.PHONY : run clean

run : pic1dp clean
	${MPIEXEC} -n ${NPE_RUN} ./$< 2>&1 | tee $<.log

pic1dp : ../build/pic1dp
	cp ../build/pic1dp ./

../build/pic1dp : $(wildcard ../src/*)
	${MAKE} -C ../build -j ${NPE_BUILD} build

clean :
	rm -f pic1dp.log pic1dp.out

