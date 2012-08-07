! input parameters and functions
module pic1dp_input
use pic1dp_global, only: kpr
implicit none
#include "finclude/petscdef.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!
! termination conditions !
!!!!!!!!!!!!!!!!!!!!!!!!!!
! whichever of the following conditions is reached, the code stops

! maximum # of time steps
PetscInt, parameter :: input_ntime_max = 900000

! maximum physical time (1 / omega_pe)
PetscReal, parameter :: input_time_max = 100.0_kpr

end module pic1dp_input

