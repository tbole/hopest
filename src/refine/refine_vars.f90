#include "hopest_f.h"
MODULE MODH_Refine_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------


TYPE tRefinePtr 
  TYPE(tRefine),POINTER ::rp
END TYPE tRefinePtr 


TYPE tRefine 
  INTEGER                   :: refineType
  INTEGER                   :: refineLevel
  TYPE(tRefineGeom),POINTER :: rGeo
  TYPE(tRefineBC),POINTER   :: rBC
END TYPE tRefine


TYPE tRefineBC 
  CHARACTER(LEN=255)        :: RefineBoundaryName
  INTEGER                   :: refineBCindex
  INTEGER,ALLOCATABLE       :: TreeSidesToRefine(:,:)
  REAL                      :: targetLength
END TYPE tRefineBC

TYPE tRefineGeom 
  INTEGER                   :: refineGeomType
  REAL                      :: sphereCenter(3)
  REAL                      :: sphereRadius
  REAL                      :: cylinderCenter(3)
  REAL                      :: cylinderAxis(3)
  REAL                      :: cylinderRadius
  REAL                      :: shellCenter(3)
  REAL                      :: shellRadius_inner
  REAL                      :: shellRadius_outer
  REAL                      :: boxBoundary(6)
  REAL                      :: boxSurf_x0(1:3,0:5)
  REAL                      :: boxSurf_n(1:3,0:5)
END TYPE tRefineGeom

!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                          :: Nsuper
INTEGER                          :: MaxrefineLevel
REAL,ALLOCATABLE                 :: xi_Nsuper(:) 
INTEGER                          :: nRefines
TYPE(tRefinePtr),POINTER         :: refines(:)
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE GETNEWREFINE
  MODULE PROCEDURE GETNEWREFINE
END INTERFACE
CONTAINS

FUNCTION GETNEWREFINE()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tRefine),POINTER :: getNewRefine
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNode,iLocSide
!===================================================================================================================================
ALLOCATE(getNewRefine)
ALLOCATE(getNewRefine%rBC)
ALLOCATE(getNewRefine%rGeo)
getNewRefine%RefineType=0
getNewRefine%RefineLevel=0
END FUNCTION GETNEWREFINE

END MODULE MODH_Refine_Vars
