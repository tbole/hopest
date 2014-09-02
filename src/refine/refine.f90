#include "hopest_f.h"

MODULE MODH_Refine
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE RefineMesh
  MODULE PROCEDURE RefineMesh
END INTERFACE

PUBLIC::RefineMesh
!===================================================================================================================================

CONTAINS


SUBROUTINE RefineMesh()
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Refine_Vars
USE MODH_Mesh_Vars,     ONLY: Ngeo
USE MODH_Refine_Binding,ONLY: p4_refine_mesh
USE MODH_P4EST_Vars,    ONLY: p4est,mesh,geom
USE MODH_Readintools,   ONLY: GETINT,CNTSTR
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_FUNPTR)              :: refineFunc
INTEGER                     :: iRefine,nRefines,i
INTEGER,ALLOCATABLE         :: refineType(:)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'BUILD P4EST MESH AND REFINE ...'
SWRITE(UNIT_StdOut,'(132("-"))')

!ALLOCATE(RefineList(nTrees))
!RefineList=0

nRefines=CNTSTR('refineType',1)
SWRITE(UNIT_stdOut,'(A,I0)')'Number of refines : ',nRefines
ALLOCATE(refineType(nRefines))

NSuper =GETINT('Nsuper','2') ! default conform refinement
ALLOCATE(Xi_Nsuper(0:Nsuper))
DO i=0,Nsuper
  Xi_Nsuper(i)=-1.+2.*REAL(i)/REAL(Nsuper)
END DO
  

DO iRefine=1,nRefines
  refineType =GETINT('refineType','1') ! default conform refinement
  refineLevel=GETINT('refineLevel','1')
  ! Transform input mesh to adapted mesh
  ! Do refinement and save p4est refine
  SELECT CASE(refineType(iRefine))
  CASE(1)  ! conform refinement of the whole mesh
    refineFunc=C_FUNLOC(RefineAll)
  CASE(2)  ! Boundary Element refinement
    refineBCIndex=GETINT('refineBCIndex','1')
    CALL InitRefineBoundaryElems()
    refineFunc=C_FUNLOC(RefineByList)
  CASE(3) ! Refine by Geometry
    CALL InitRefineGeom()
    refineFunc=C_FUNLOC(RefineByGeom)
  CASE(11)
    refineFunc=C_FUNLOC(RefineFirst)
  CASE DEFAULT
    STOP 'refineType is not defined'
  END SELECT
   CALL p4_refine_mesh(p4est,refineFunc,refineLevel,geom, & !IN
                    mesh)                               !OUT
   SDEALLOCATE(TreeSidesToRefine)
END DO

!SDEALLOCATE(RefineList)
END SUBROUTINE RefineMesh


SUBROUTINE InitRefineBoundaryElems()
!===================================================================================================================================
! init the refinment list
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Refine_Vars,ONLY: TreeSidesToRefine,refineBCIndex
USE MODH_Mesh_Vars,  ONLY: Trees,nTrees
USE MODH_P4EST_Vars,  ONLY: P2H_FaceMap
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iTree,iSide
!===================================================================================================================================
! These are the refinement functions which are called by p4est
ALLOCATE(TreeSidesToRefine(0:5,1:nTrees))
TreeSidesToRefine=0
DO iTree=1,nTrees
  DO iSide=0,5
    IF (Trees(iTree)%ep%Side(P2H_FaceMap(iSide))%sp%BCIndex.EQ.refineBCIndex) THEN
      TreeSidesToRefine(iSide,iTree)=1
    END IF
  END DO
END DO


END SUBROUTINE InitRefineBoundaryElems


SUBROUTINE InitRefineGeom()
!===================================================================================================================================
! init the geometric refinment
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Refine_Vars, ONLY: refineGeomType,sphereCenter,sphereRadius,boxBoundary 
USE MODH_Refine_Vars, ONLY: shellRadius_inner,shellRadius_outer,shellCenter
USE MODH_Refine_Vars, ONLY: cylinderCenter,CylinderAxis,CylinderRadius 
USE MODH_Refine_Vars, ONLY: BoxSurf_x0,BoxSurf_n
USE MODH_Readintools, ONLY:GETREALARRAY,GETINT,GETREAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: boxC(1:3,0:1,0:1,0:1)
!===================================================================================================================================
! These are the refinement functions which are called by p4est
refineGeomType =GETINT('refineGeomType') ! 
SELECT CASE(refineGeomType)
  CASE(1) ! Sphere
    sphereCenter =GETREALARRAY('sphereCenter',3)
    sphereRadius =GETREAL('sphereRadius')
  CASE(11) ! Sphere Shell
    shellCenter =GETREALARRAY('shellCenter',3)
    shellRadius_inner =GETREAL('shellRadius_inner')
    shellRadius_outer =GETREAL('shellRadius_outer')
    IF(shellRadius_inner .GE. shellRadius_outer) STOP 'inner radius is greater or equal to outer radius!'
  CASE(2) ! Box 
    boxBoundary=GETREALARRAY('boxBoundary',6)
  CASE(3) ! CYLINDER
    cylinderCenter =GETREALARRAY('cylinderCenter',3)
    CylinderAxis =GETREALARRAY('CylinderAxis',3)
    CylinderAxis=CylinderAxis/VECNORM(CylinderAxis)
    CylinderRadius =GETREAL('CylinderRadius')
  CASE(4) ! 8 corner Box (p4est standard) only planar surfaces!
    boxC(:,0,0,0)=GETREALARRAY('boxCorner000',3) 
    boxC(:,1,0,0)=GETREALARRAY('boxCorner100',3) 
    boxC(:,0,1,0)=GETREALARRAY('boxCorner010',3) 
    boxC(:,1,1,0)=GETREALARRAY('boxCorner110',3) 
    boxC(:,0,0,1)=GETREALARRAY('boxCorner001',3) 
    boxC(:,1,0,1)=GETREALARRAY('boxCorner101',3) 
    boxC(:,0,1,1)=GETREALARRAY('boxCorner011',3) 
    boxC(:,1,1,1)=GETREALARRAY('boxCorner111',3) 
    boxSurf_x0(:,0)=boxC(:,0,0,0)
    boxSurf_x0(:,2)=boxC(:,0,0,0)
    boxSurf_x0(:,4)=boxC(:,0,0,0)
    boxSurf_n(:, 0)=CROSS((boxC(:,0,0,1)-boxC(:,0,0,0)),boxC(:,0,1,0)-boxC(:,0,0,0))
    boxSurf_n(:, 2)=CROSS((boxC(:,1,0,0)-boxC(:,0,0,0)),boxC(:,0,0,1)-boxC(:,0,0,0))
    boxSurf_n(:, 4)=CROSS((boxC(:,0,1,0)-boxC(:,0,0,0)),boxC(:,1,0,0)-boxC(:,0,0,0))
    boxSurf_x0(:,1)=boxC(:,1,1,1)
    boxSurf_x0(:,3)=boxC(:,1,1,1)
    boxSurf_x0(:,5)=boxC(:,1,1,1)
    boxSurf_n(:, 1)=CROSS((boxC(:,1,0,1)-boxC(:,1,1,1)),boxC(:,1,1,0)-boxC(:,1,1,1))
    boxSurf_n(:, 3)=CROSS((boxC(:,1,1,0)-boxC(:,1,1,1)),boxC(:,0,1,1)-boxC(:,1,1,1))
    boxSurf_n(:, 5)=CROSS((boxC(:,0,1,1)-boxC(:,1,1,1)),boxC(:,1,0,1)-boxC(:,1,1,1))
  CASE DEFAULT
    STOP 'no refineGeomType defined'
END SELECT

END SUBROUTINE InitRefineGeom

FUNCTION RefineAll(x,y,z,tree,level,childID) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Refine_Vars, ONLY: RefineLevel
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
P4EST_F90_QCOORD,    INTENT(IN),VALUE :: x,y,z
P4EST_F90_TOPIDX,    INTENT(IN),VALUE :: tree
P4EST_F90_QLEVEL,    INTENT(IN),VALUE :: level
INTEGER(KIND=C_INT ),INTENT(IN),VALUE :: childID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT) :: refineAll
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
RefineAll=0
IF(level.GE.refineLevel) RETURN
refineAll=1
END FUNCTION RefineAll

FUNCTION RefineByList(x,y,z,tree,level,childID) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Refine_Vars, ONLY: refineLevel,TreeSidesToRefine
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
P4EST_F90_QCOORD,    INTENT(IN),VALUE :: x,y,z
P4EST_F90_TOPIDX,    INTENT(IN),VALUE :: tree
P4EST_F90_QLEVEL,    INTENT(IN),VALUE :: level
INTEGER(KIND=C_INT ),INTENT(IN),VALUE :: childID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                   :: refineByList
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
P4EST_F90_QCOORD                      ::  IntSize,length
!-----------------------------------------------------------------------------------------------------------------------------------
RefineByList=0
IF(level.GE.refineLevel) RETURN
!IF (level.EQ.0) RefineByList=SUM(TreeSidesToRefine(:,tree))
!IF (level.GE.1) RefineByList=TreeSidesToRefine(childID+1,tree)

IntSize = 2**19
length  = 2**(19-level)
IF(x.EQ.0) THEN
  RefineByList=RefineBylist+TreeSidesToRefine(0,tree)
END IF
IF(y.EQ.0) THEN
   RefineByList=RefineBylist+TreeSidesToRefine(2,tree)
END IF
IF(z.EQ.0) THEN
  RefineByList=RefineBylist+TreeSidesToRefine(4,tree)
END IF
IF(x.EQ.intSize-length) THEN
  RefineByList=RefineBylist+TreeSidesToRefine(1,tree)
END IF
IF(y.EQ.intSize-length) THEN
  RefineByList=RefineBylist+TreeSidesToRefine(3,tree)
END IF
IF(z.EQ.intSize-length) THEN
  RefineByList=RefineBylist+TreeSidesToRefine(5,tree)
END IF

END FUNCTION RefineByList


FUNCTION RefineByGeom(x,y,z,tree,level,childID) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Refine_Vars, ONLY: refineLevel
USE MODH_Refine_Vars, ONLY: NSuper,Xi_Nsuper
USE MODH_Mesh_Vars,   ONLY: XGeo,Ngeo
USE MODH_Mesh_Vars,   ONLY: wBary_Ngeo,xi_Ngeo
USE MODH_Mesh_Vars,   ONLY: doSplineInterpolation
USE MODH_Spline1D,    ONLY: GetSplineVandermonde
USE MODH_Basis,       ONLY: LagrangeInterpolationPolys 
USE MODH_ChangeBasis, ONLY: ChangeBasis3D_XYZ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
P4EST_F90_QCOORD,    INTENT(IN),VALUE :: x,y,z
P4EST_F90_TOPIDX,    INTENT(IN),VALUE :: tree
P4EST_F90_QLEVEL,    INTENT(IN),VALUE :: level
INTEGER(KIND=C_INT ),INTENT(IN),VALUE :: childID
REAL                                  :: xi0(3)
REAL                                  :: xiBary(3)
REAL                                  :: dxi,length
REAL,DIMENSION(0:Nsuper,0:Ngeo)       :: Vdm_xi,Vdm_eta,Vdm_zeta
REAL                                  :: XCorner(3), XBaryElem(3),test
INTEGER                               :: IntSize
REAL                                  :: sIntSize
REAL                                  :: XGeoSuper(3,0:Nsuper,0:Nsuper,0:Nsuper)
REAL                                  :: l_xi(0:NGeo),l_eta(0:NGeo),l_zeta(0:NGeo),l_etazeta
INTEGER                               :: i,j,k
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT) :: refineByGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! TODO :
! Sichereres Kriterium: kleinste Abstand zwischen Baryzentrum und geom. Objekt muss kleiner sein als Element-Diagonale!
! Beispiel Kugel/Shell: Länge des Vektors von Elementbaryzentrum zum Kugelmittelpunkt minus Kugelradius kleiner als Element-Diagonale
refineByGeom = 0
IF(level.GE.refineLevel) RETURN

IntSize=2**19 !TODO: use sIntSize from mesh_vars. not initialized at this stage yet
sIntSize=1./REAL(IntSize)

! transform p4est first corner coordinates (integer from 0... intsize) to [-1,1] reference element of its tree
xi0(1)=-1.+2.*REAL(x)*sIntSize
xi0(2)=-1.+2.*REAL(y)*sIntSize
xi0(3)=-1.+2.*REAL(z)*sIntSize
! length of the element / quadrant in reference coordinates of its tree [-1,1]
length=2./REAL(2**level)
! Build Vandermonde matrices for each parameter range in xi, eta,zeta
IF(doSplineInterpolation)THEN
  CALL getSplineVandermonde(Ngeo+1,Nsuper+1,Vdm_xi(:,:)  ,xi0(1)+0.5*(xi_Nsuper(:)+1.)*Length)
  CALL getSplineVandermonde(Ngeo+1,Nsuper+1,Vdm_eta(:,:) ,xi0(2)+0.5*(xi_Nsuper(:)+1.)*Length)
  CALL getSplineVandermonde(Ngeo+1,Nsuper+1,Vdm_zeta(:,:),xi0(3)+0.5*(xi_Nsuper(:)+1.)*Length)
ELSE !polynomials
  DO i=0,Nsuper 
    dxi=0.5*(xi_Nsuper(i)+1.)*Length
    CALL LagrangeInterpolationPolys(xi0(1) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_xi(i,:)) 
    CALL LagrangeInterpolationPolys(xi0(2) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_eta(i,:)) 
    CALL LagrangeInterpolationPolys(xi0(3) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_zeta(i,:)) 
  END DO
END IF
!interpolate tree HO mapping to quadrant HO mapping
CALL ChangeBasis3D_XYZ(3,Ngeo,Nsuper,Vdm_xi,Vdm_eta,Vdm_zeta,XGeo(:,:,:,:,tree),XGeoSuper(:,:,:,:))


! Barycenter
xiBary=xi0+0.5*length
! Build Vandermonde matrices for each parameter range in xi, eta,zeta
IF(doSplineInterpolation)THEN
  CALL getSplineVandermonde(Ngeo+1,1,l_xi(:)  ,xiBary(1))
  CALL getSplineVandermonde(Ngeo+1,1,l_eta(:) ,xiBary(2))
  CALL getSplineVandermonde(Ngeo+1,1,l_zeta(:),xiBary(3))
ELSE !polynomials
  CALL LagrangeInterpolationPolys(xiBary(1),Ngeo,xi_Ngeo,wBary_Ngeo,l_xi(:)) 
  CALL LagrangeInterpolationPolys(xiBary(2),Ngeo,xi_Ngeo,wBary_Ngeo,l_eta(:)) 
  CALL LagrangeInterpolationPolys(xiBary(3),Ngeo,xi_Ngeo,wBary_Ngeo,l_zeta(:)) 
END IF
!interpolate tree HO mapping to quadrant HO mapping
XBaryElem(:)=0.
DO k=0,NGeo
  DO j=0,NGeo
    l_etazeta=l_eta(j)*l_zeta(k)
    DO i=0,NGeo
      XBaryElem(:)=XBaryElem(:)+Xgeo(:,i,j,k,tree)*l_xi(i)*l_etazeta
    END DO
  END DO
END DO
! check barycenter
IF (CheckGeom(XBaryElem)) THEN
  refineByGeom = 1
  RETURN
END IF
! check Nsuper Nodes
DO k=0,Nsuper
  DO j=0,Nsuper
    DO i=0,Nsuper
      IF (CheckGeom(XGeoSuper(:,i,j,k))) THEN
        refineByGeom = 1
        RETURN
      END IF
    END DO ! i
  END DO ! j 
END DO ! k
END FUNCTION RefineByGeom


FUNCTION CheckGeom(x_in)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MODH_Globals,     ONLY: VECNORM
USE MODH_Refine_Vars, ONLY: refineGeomType
USE MODH_Refine_Vars, ONLY: sphereCenter,sphereRadius,boxBoundary
USE MODH_Refine_Vars, ONLY: shellRadius_inner,shellRadius_outer,shellCenter
USE MODH_Refine_Vars, ONLY: cylinderCenter,CylinderAxis,CylinderRadius 
USE MODH_Refine_Vars, ONLY: BoxSurf_x0,BoxSurf_n
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) ::x_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL :: checkGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iSide
REAL            :: test
!-----------------------------------------------------------------------------------------------------------------------------------
CheckGeom=.FALSE.
SELECT CASE (refineGeomType)
CASE(1)   ! SPHERE
  test=VECNORM(x_in-SphereCenter)
  IF (test.LE.sphereRadius) THEN
    CheckGeom= .TRUE.
  END IF
CASE(11) ! sphere shell
  test=VECNORM(x_in-shellCenter)
  IF ((test.LE.shellRadius_outer) .AND. (test.GE.shellRadius_inner)) THEN
    CheckGeom= .TRUE.
  END IF
CASE(2) ! box
  IF (x_in(1) .GE. boxBoundary(1) .AND. x_in(1) .LE. boxBoundary(2) .AND. & 
      x_in(2) .GE. boxBoundary(3) .AND. x_in(2) .LE. boxBoundary(4) .AND. &
      x_in(3) .GE. boxBoundary(5) .AND. x_in(3) .LE. boxBoundary(6)) THEN
    CheckGeom= .TRUE.
  END IF
CASE(3)   ! CYLINDER
  test=VECNORM((x_in-CylinderCenter)-CylinderAxis*SUM(CylinderAxis*(x_in-CylinderCenter)))
  IF (test.LE.cylinderRadius) THEN
    CheckGeom= .TRUE.
  END IF
CASE(4)   ! CORNER BOX
  CheckGeom= .TRUE.
  DO iSide=0,5
    test=SUM((x_in(:)-BoxSurf_x0(:,iSide))*boxSurf_n(:,iSide))
    IF(test.GT.0.)THEN
      CheckGeom= .FALSE.
      EXIT
    END IF
  END DO
      
END SELECT
END FUNCTION CheckGeom

FUNCTION RefineFirst(x,y,z,tree,level,childID) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Refine_Vars, ONLY: refineLevel
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
P4EST_F90_QCOORD,INTENT(IN),VALUE :: x,y,z
P4EST_F90_TOPIDX,INTENT(IN),VALUE :: tree
P4EST_F90_QLEVEL,INTENT(IN),VALUE :: level
P4EST_F90_LOCIDX,INTENT(IN),VALUE :: childID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                      :: refineFirst
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
refineFirst = 0
IF(level.GE.refineLevel) RETURN
IF(tree.EQ.1)THEN
  refineFirst=1
ELSE
  refineFirst=0
END IF
END FUNCTION RefineFirst

END MODULE MODH_Refine
