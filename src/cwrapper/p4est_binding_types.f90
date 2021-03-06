#include "hopest_f.h"

MODULE MOD_P4estBindingTypes
!===================================================================================================================================
! Fortran <-> C++ wrapper routine for the P4est Routines
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING  
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

TYPE :: t_p4est_ptr
  TYPE(C_PTR)                :: p4est
  TYPE(C_PTR)                :: mesh
END TYPE t_p4est_ptr


END MODULE MOD_P4estBindingTypes
