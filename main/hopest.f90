#include "hopest_f.h"

PROGRAM Hopest
!===================================================================================================================================
! Control program of the Flexi code. Initialization of the computation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,  ONLY:IgnoredStrings
USE MOD_Mesh,         ONLY:InitMesh,FinalizeMesh
USE MOD_IO_HDF5,      ONLY:InitIO

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL :: Time
!===================================================================================================================================
CALL InitIO()
! hier fehlt der HEADER

! Measure init duration
StartTime=RUNTIME()

! Initialization
CALL InitMesh()
CALL IgnoredStrings()

!Finalize
CALL FinalizeMesh()
! Measure simulation duration
Time=RUNTIME()
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A,F8.2,A)') ' HOPEST FINISHED! [',Time-StartTime,' sec ]'
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM Hopest
