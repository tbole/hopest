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
INTERFACE InitRefine
  MODULE PROCEDURE InitRefine
END INTERFACE

INTERFACE RefineMesh
  MODULE PROCEDURE RefineMesh
END INTERFACE

PUBLIC::InitRefine
PUBLIC::RefineMesh
!===================================================================================================================================

CONTAINS


SUBROUTINE InitRefine()
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Refine_Vars
USE MODH_Mesh_Vars,     ONLY: Ngeo
USE MODH_Readintools,   ONLY: GETINT,CNTSTR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iRefine,i
TYPE(tRefine),POINTER :: aRef
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'INIT REFINE ...'
SWRITE(UNIT_StdOut,'(132("-"))')


NSuper =GETINT('Nsuper','2') ! default conform refinement
ALLOCATE(Xi_Nsuper(0:Nsuper))
DO i=0,Nsuper
  Xi_Nsuper(i)=-1.+2.*REAL(i)/REAL(Nsuper)
END DO
  
nRefines=CNTSTR('refineType',1)
SWRITE(UNIT_stdOut,'(A,I0)')'Number of refine Types : ',nRefines
MaxRefineLevel=0
ALLOCATE(Refines(nRefines))
DO iRefine=1,nRefines
  Refines(iRefine)%rp=>GETNEWREFINE()
  aRef=>Refines(iRefine)%rp
  aRef%refineType =GETINT('refineType','1') ! default conform refinement
  aRef%refineLevel=GETINT('refineLevel','1')
  MaxRefineLevel=MAX(MaxRefineLevel,aRef%refineLevel)
  ! Transform input mesh to adapted mesh
  ! Do refinement and save p4est refine
  SELECT CASE(aRef%refineType)
  CASE(1,11)  ! conform refinement of the whole mesh / first element, do nothing
  CASE(2)  ! Boundary Element refinement
    CALL InitRefineBC(aRef%rBC,aRef%refineLevel)
  CASE(3) ! Refine by Geometry
    CALL InitRefineGeom(aRef%rGeo)
  CASE DEFAULT
    STOP 'refineType is not defined'
  END SELECT
END DO

SWRITE(UNIT_stdOut,'(A)')'INIT REFINE DONE.'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRefine


SUBROUTINE InitRefineBC(rBC,refineLevel)
!===================================================================================================================================
! init the refinment list
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,   ONLY: Trees,nTrees
USE MODH_Mesh_Vars,   ONLY: BoundaryName,nBCs
USE MODH_Refine_Vars, ONLY: tRefineBC
USE MODH_P4EST_Vars,  ONLY: P2H_FaceMap
USE MODH_Readintools, ONLY: GETSTR,GETREAL
USE MODH_p4est_Vars,  ONLY: sIntSize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: refineLevel 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tRefineBC),POINTER :: rBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iTree,iSide,iBC
!===================================================================================================================================
! These are the refinement functions which are called by p4est
rBC%RefineBoundaryName=GETSTR('refineBoundaryName')
rBC%refineBCIndex=0
DO iBC=1,nBCs
  IF(INDEX(TRIM(rBC%RefineBoundaryName),TRIM(BoundaryName(iBC))) .NE.0)THEN
    WRITE(Unit_StdOut,'(A,A)')    ' |  Boundary for refinement found | ',TRIM(BoundaryName(iBC))
    rBC%refineBCIndex=iBC
    EXIT
  END IF
END DO
IF(rBC%refineBCindex.EQ.0) THEN
  WRITE(Unit_stdOut,'(A,A,A)')'Boundary ',TRIM(rBC%RefineBoundaryName),' not found!'
  STOP
END IF

ALLOCATE(rBC%TreeSidesToRefine(0:5,1:nTrees))
rBC%TreeSidesToRefine=0
DO iTree=1,nTrees
  DO iSide=0,5
    IF (Trees(iTree)%ep%Side(P2H_FaceMap(iSide))%sp%BCIndex.EQ.rBC%refineBCIndex) THEN
      rBC%TreeSidesToRefine(iSide,iTree)=1
    END IF
  END DO
END DO

! set the target length in p4est qcoords
rBC%targetLength=GETREAL('refineBoundaryThickness','0')
IF(rBC%targetLength.EQ.0.) THEN  ! only refine in direct neighborhood of the BC
  rBC%targetLength=2./REAL(2**(refineLevel))
END IF
END SUBROUTINE InitRefineBC


SUBROUTINE InitRefineGeom(rGeo)
!===================================================================================================================================
! init the geometric refinment
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Refine_Vars, ONLY: tRefineGeom
USE MODH_Readintools, ONLY: GETREALARRAY,GETINT,GETREAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tRefineGeom),POINTER :: rGeo
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: boxC(1:3,0:1,0:1,0:1)
!===================================================================================================================================
rGeo%refineGeomType =GETINT('refineGeomType') ! 
! These are the refinement functions which are called by p4est
SELECT CASE(rGeo%refineGeomType)
  CASE(1) ! Sphere
    rGeo%sphereCenter =GETREALARRAY('sphereCenter',3)
    rGeo%sphereRadius =GETREAL('sphereRadius')
  CASE(11) ! Sphere Shell
    rGeo%shellCenter =GETREALARRAY('shellCenter',3)
    rGeo%shellRadius_inner =GETREAL('shellRadius_inner')
    rGeo%shellRadius_outer =GETREAL('shellRadius_outer')
    IF(rGeo%shellRadius_inner .GE. rGeo%shellRadius_outer) STOP 'inner radius is greater or equal to outer radius!'
  CASE(2) ! Box 
    rGeo%boxBoundary=GETREALARRAY('boxBoundary',6)
  CASE(3) ! CYLINDER
    rGeo%cylinderCenter =GETREALARRAY('cylinderCenter',3)
    rGeo%CylinderAxis   =GETREALARRAY('CylinderAxis',3)
    rGeo%CylinderAxis   =rGeo%CylinderAxis/VECNORM(rGeo%CylinderAxis)
    rGeo%CylinderRadius =GETREAL('CylinderRadius')
  CASE(4) ! 8 corner Box (p4est standard) only planar surfaces!
    boxC(:,0,0,0)=GETREALARRAY('boxCorner000',3) 
    boxC(:,1,0,0)=GETREALARRAY('boxCorner100',3) 
    boxC(:,0,1,0)=GETREALARRAY('boxCorner010',3) 
    boxC(:,1,1,0)=GETREALARRAY('boxCorner110',3) 
    boxC(:,0,0,1)=GETREALARRAY('boxCorner001',3) 
    boxC(:,1,0,1)=GETREALARRAY('boxCorner101',3) 
    boxC(:,0,1,1)=GETREALARRAY('boxCorner011',3) 
    boxC(:,1,1,1)=GETREALARRAY('boxCorner111',3) 
    rGeo%boxSurf_x0(:,0)=boxC(:,0,0,0)
    rGeo%boxSurf_x0(:,2)=boxC(:,0,0,0)
    rGeo%boxSurf_x0(:,4)=boxC(:,0,0,0)
    rGeo%boxSurf_n(:, 0)=CROSS((boxC(:,0,0,1)-boxC(:,0,0,0)),boxC(:,0,1,0)-boxC(:,0,0,0))
    rGeo%boxSurf_n(:, 2)=CROSS((boxC(:,1,0,0)-boxC(:,0,0,0)),boxC(:,0,0,1)-boxC(:,0,0,0))
    rGeo%boxSurf_n(:, 4)=CROSS((boxC(:,0,1,0)-boxC(:,0,0,0)),boxC(:,1,0,0)-boxC(:,0,0,0))
    rGeo%boxSurf_x0(:,1)=boxC(:,1,1,1)
    rGeo%boxSurf_x0(:,3)=boxC(:,1,1,1)
    rGeo%boxSurf_x0(:,5)=boxC(:,1,1,1)
    rGeo%boxSurf_n(:, 1)=CROSS((boxC(:,1,0,1)-boxC(:,1,1,1)),boxC(:,1,1,0)-boxC(:,1,1,1))
    rGeo%boxSurf_n(:, 3)=CROSS((boxC(:,1,1,0)-boxC(:,1,1,1)),boxC(:,0,1,1)-boxC(:,1,1,1))
    rGeo%boxSurf_n(:, 5)=CROSS((boxC(:,0,1,1)-boxC(:,1,1,1)),boxC(:,1,0,1)-boxC(:,1,1,1))
  CASE DEFAULT
    STOP 'no refineGeomType defined'
END SELECT

END SUBROUTINE InitRefineGeom


SUBROUTINE RefineMesh()
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Refine_Vars,   ONLY: MaxRefineLevel
USE MODH_Refine_Binding,ONLY: p4_refine_mesh
USE MODH_P4EST_Vars,    ONLY: p4est,mesh,geom
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_FUNPTR)        :: refineFunction
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'BUILD P4EST MESH AND REFINE ...'
SWRITE(UNIT_StdOut,'(132("-"))')
refineFunction=C_FUNLOC(RefineFunc)
CALL p4_refine_mesh(p4est,refineFunction,maxRefineLevel,geom, & !IN
                    mesh)                               !OUT
END SUBROUTINE RefineMesh


FUNCTION RefineFunc(x,y,z,tree,level,childID) BIND(C)
!===================================================================================================================================
! Function used by p4est for each quadrant
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Refine_Vars, ONLY: Refines,tRefine,nRefines
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
INTEGER(KIND=C_INT) :: refineFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tRefine),POINTER :: aRef
INTEGER               :: iRefine
!===================================================================================================================================
refineFunc=0
DO iRefine=1,nRefines
  aRef=>Refines(iRefine)%rp
  IF(level.GE. aRef%refineLevel) CYCLE
  SELECT CASE(aRef%refineType)
  CASE(1) !Refineall
    refineFunc=refineFunc+1
  CASE(11) !RefineFirst (Debugging)
    IF(tree.EQ.1)THEN
      refineFunc=refineFunc+1
    END IF
  CASE(2)
    refineFunc=refineFunc+RefineBC(x,y,z,tree,level,aRef%rBC) 
  CASE(3)
    refineFunc=refineFunc+RefineByGeom(x,y,z,tree,level,aRef%rGeo) 
  END SELECT
  IF(refineFunc.GT.0) RETURN
END DO
END FUNCTION refineFunc


FUNCTION RefineBC(x,y,z,tree,level,rBC) 
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_p4est_Vars, ONLY:IntSize,sIntSize
USE MODH_Refine_vars,ONLY:tRefineBC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
P4EST_F90_QCOORD,    INTENT(IN),VALUE :: x,y,z
P4EST_F90_TOPIDX,    INTENT(IN),VALUE :: tree
P4EST_F90_QLEVEL,    INTENT(IN),VALUE :: level
TYPE(tRefineBC),POINTER               :: rBC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)      :: refineBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: dir
REAL                     :: xi(3),length
!P4EST_F90_QCOORD         :: length
!-----------------------------------------------------------------------------------------------------------------------------------
RefineBC=0
length=1./REAL(2**level)
xi(1)=REAL(x)*sIntSize
xi(2)=REAL(y)*sIntSize
xi(3)=REAL(z)*sIntSize
DO dir=1,3
  IF(xi(dir).LT.rBC%targetLength) THEN
    RefineBC=RefineBC+rBC%TreeSidesToRefine(2*(dir-1),tree)
  END IF
  IF((1.-xi(dir)-length).LT.rBC%targetLength) THEN
    RefineBC=RefineBC+rBC%TreeSidesToRefine(2*(dir-1)+1,tree)
  END IF
END DO
END FUNCTION RefineBC


FUNCTION RefineByGeom(x,y,z,tree,level,rGeo)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_p4est_Vars,  ONLY: sIntSize
USE MODH_Refine_Vars, ONLY: NSuper,Xi_Nsuper
USE MODH_Refine_Vars, ONLY: tRefineGeom 
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
TYPE(tRefineGeom),POINTER             :: rGeo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT) :: refineByGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                  :: xi0(3)
REAL                                  :: xiBary(3)
REAL                                  :: dxi,length
REAL,DIMENSION(0:Nsuper,0:Ngeo)       :: Vdm_xi,Vdm_eta,Vdm_zeta
REAL                                  :: XCorner(3), XBaryElem(3),test
REAL                                  :: XGeoSuper(3,0:Nsuper,0:Nsuper,0:Nsuper)
REAL                                  :: l_xi(0:NGeo),l_eta(0:NGeo),l_zeta(0:NGeo),l_etazeta
INTEGER                               :: i,j,k
!===================================================================================================================================
! TODO :
! Sichereres Kriterium: kleinste Abstand zwischen Baryzentrum und geom. Objekt muss kleiner sein als Element-Diagonale!
! Beispiel Kugel/Shell: Länge des Vektors von Elementbaryzentrum zum Kugelmittelpunkt minus Kugelradius kleiner als Element-Diagonale
refineByGeom = 0


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
IF (CheckGeom(XBaryElem,rGeo)) THEN
  refineByGeom = 1
  RETURN
END IF
! check Nsuper Nodes
DO k=0,Nsuper
  DO j=0,Nsuper
    DO i=0,Nsuper
      IF (CheckGeom(XGeoSuper(:,i,j,k),rGeo)) THEN
        refineByGeom = 1
        RETURN
      END IF
    END DO ! i
  END DO ! j 
END DO ! k
END FUNCTION RefineByGeom


FUNCTION CheckGeom(x_in,rGeo)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MODH_Globals,     ONLY: VECNORM
USE MODH_Refine_Vars, ONLY: tRefineGeom 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) ::x_in(3)
TYPE(tRefineGeom),POINTER     :: rGeo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL :: checkGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iSide
REAL            :: test
!===================================================================================================================================
CheckGeom=.FALSE.
SELECT CASE (rGeo%refineGeomType)
CASE(1)   ! SPHERE
  test=VECNORM(x_in-rGeo%SphereCenter)
  IF (test.LE.rGeo%sphereRadius) THEN
    CheckGeom= .TRUE.
  END IF
CASE(11) ! sphere shell
  test=VECNORM(x_in-rGeo%shellCenter)
  IF ((test.LE.rGeo%shellRadius_outer) .AND. (test.GE.rGeo%shellRadius_inner)) THEN
    CheckGeom= .TRUE.
  END IF
CASE(2) ! box
  IF (x_in(1) .GE. rGeo%boxBoundary(1) .AND. x_in(1) .LE. rGeo%boxBoundary(2) .AND. & 
      x_in(2) .GE. rGeo%boxBoundary(3) .AND. x_in(2) .LE. rGeo%boxBoundary(4) .AND. &
      x_in(3) .GE. rGeo%boxBoundary(5) .AND. x_in(3) .LE. rGeo%boxBoundary(6)) THEN
    CheckGeom= .TRUE.
  END IF
CASE(3)   ! CYLINDER
  test=VECNORM((x_in-rGeo%CylinderCenter)-rGeo%CylinderAxis*SUM(rGeo%CylinderAxis*(x_in-rGeo%CylinderCenter)))
  IF (test.LE.rGeo%cylinderRadius) THEN
    CheckGeom= .TRUE.
  END IF
CASE(4)   ! CORNER BOX
  CheckGeom= .TRUE.
  DO iSide=0,5
    test=SUM((x_in(:)-rGeo%BoxSurf_x0(:,iSide))*rGeo%boxSurf_n(:,iSide))
    IF(test.GT.0.)THEN
      CheckGeom= .FALSE.
      EXIT
    END IF
  END DO
END SELECT
END FUNCTION CheckGeom

END MODULE MODH_Refine
