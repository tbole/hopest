#include "hopest_f.h"

MODULE MOD_MeshFromP4EST
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE BuildMeshFromP4EST
  MODULE PROCEDURE BuildMeshFromP4EST
END INTERFACE

INTERFACE BuildHOMeh
  MODULE PROCEDURE BuildHOMesh
END INTERFACE

PUBLIC::BuildMeshFromP4EST
PUBLIC::BuildHOMesh
!===================================================================================================================================

CONTAINS


SUBROUTINE BuildMeshFromP4EST()
!===================================================================================================================================
! Subroutine to translate p4est mesh datastructure to HOPR datastructure
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_p4estBinding
USE MOD_Output_Vars, ONLY:Projectname
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                 :: QT,QQ,QF,QH
TYPE(tElem),POINTER         :: aQuad,nbQuad,Tree
TYPE(tSide),POINTER         :: aSide
INTEGER                     :: iQuad,iMortar,iElem
INTEGER                     :: PSide,PnbSide,nbSide
INTEGER                     :: nbQuadInd
INTEGER                     :: PMortar,PFlip,HFlip,QHInd
INTEGER                     :: iLocSide
INTEGER                     :: StartQuad,EndQuad
INTEGER                     :: BClocSide,BCindex
!===================================================================================================================================
IF(MESHInitIsDone) RETURN
SWRITE(UNIT_stdOut,'(A)')'BUILD P4EST MESH AND REFINE ...'
SWRITE(UNIT_StdOut,'(132("-"))')

! Get arrays from p4est: use pointers for c arrays (QT,QQ,..), duplicate data for QuadCoords,Level
CALL p4est_get_mesh_info(p4est_ptr%p4est,p4est_ptr%mesh,nQuadrants,nHalfFaces)

ALLOCATE(QuadCoords(3,nQuadrants),QuadLevel(nQuadrants)) ! big to small flip
QuadCoords=0
QuadLevel=0

CALL p4est_get_quadrants(p4est_ptr%p4est,p4est_ptr%mesh,nQuadrants,nHalfFaces,& !IN
                         intsize,QT,QQ,QF,QH,QuadCoords,QuadLevel)              !OUT

CALL C_F_POINTER(QT,QuadToTree,(/nQuadrants/))
CALL C_F_POINTER(QQ,QuadToQuad,(/6,nQuadrants/))
CALL C_F_POINTER(QF,QuadToFace,(/6,nQuadrants/))
IF(nHalfFaces.GT.0) CALL C_F_POINTER(QH,QuadToHalf,(/4,nHalfFaces/))

ALLOCATE(TreeToQuad(2,nElems))
TreeToQuad(1,1)=0
TreeToQuad(2,nElems)=nQuadrants
StartQuad=0
EndQuad=0
DO iElem=1,nElems
  TreeToQuad(1,iElem)=StartQuad
  DO iQuad=StartQuad+1,nQuadrants
    IF(QuadToTree(iQuad)+1.EQ.iElem)THEN
      EndQuad=EndQuad+1
    ELSE
      TreeToQuad(2,iElem)=EndQuad
      StartQuad=EndQuad
      EXIT 
    END IF
  END DO !iQuad
END DO !iElem
!----------------------------------------------------------------------------------------------------------------------------
!             Start to build p4est datastructure in HOPEST
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
ALLOCATE(Quads(1:nQuadrants))
DO iQuad=1,nQuadrants
  Quads(iQuad)%ep=>GETNEWELEM()
  aQuad=>Quads(iQuad)%ep
  aQuad%Ind    = iQuad
  CALL CreateSides(aQuad)
  DO iLocSide=1,6
    aQuad%Side(iLocSide)%sp%flip=-999
  END DO
END DO

DO iQuad=1,nQuadrants
  aQuad=>Quads(iQuad)%ep
  Tree=>Elems(QuadToTree(iQuad)+1)%ep
  aQuad%type=Tree%type
  DO iLocSide=1,6
    aSide=>aQuad%Side(iLocSide)%sp
    ! Get P4est local side
    PSide=H2P_FaceMap(iLocSide)
    ! Get P4est neighbour side/flip/morter
    CALL EvalP4ESTConnectivity(QuadToFace(PSide+1,iQuad),PnbSide,PFlip,PMortar)
    ! transform p4est orientation to HOPR flip (magic)
    HFlip=GetHFlip(PSide,PnbSide,PFlip)  !Hflip of neighbor side!!!
    IF(PMortar.EQ.4)THEN
      ! Neighbour side is mortar (4 sides), all neighbour element sides have same orientation and local side ind
      QHInd=QuadToQuad(PSide+1,iQuad)+1
      aSide%nMortars=4
      aSide%MortarType=1             ! 1->4 case
      ALLOCATE(aSide%MortarSide(4))
      DO iMortar=1,4
        nbQuadInd=QuadToHalf(iMortar,QHInd)+1
        nbQuad=>Quads(nbQuadInd)%ep
        nbSide=P2H_FaceMap(PnbSide)
        aSide%MortarSide(iMortar)%sp=>nbQuad%side(nbSide)%sp
        aSide%MortarSide(iMortar)%sp%flip=HFlip
      END DO ! iMortar
    ELSE
      nbQuadInd=QuadToQuad(PSide+1,iQuad)+1
      nbQuad=>Quads(nbQuadInd)%ep
      nbSide=P2H_FaceMap(PnbSide)
      IF((nbQuadInd.EQ.iQuad).AND.(nbSide.EQ.iLocSide))THEN
        ! this is a boundary side: 
        BCindex=Tree%Side(iLocSide)%sp%BCindex 
        IF(BCIndex.EQ.0) STOP 'Problem in Boundary assignment'
        aSide%BCIndex=BCIndex
        NULLIFY(aSide%connection)
        aSide%Flip=0
        
      ELSE
        aSide%connection=>nbQuad%side(nbSide)%sp
        aSide%connection%flip=HFlip
      END IF !BC side
      IF(PMortar.NE.-1) aSide%MortarType= - (PMortar+1)  ! Pmortar 0...3, small side belonging to  mortar group -> -1..-4
    END IF ! PMortar
  END DO
END DO

! set master slave,  element with lower element ID is master (flip=0)
DO iQuad=1,nQuadrants
  aQuad=>Quads(iQuad)%ep
  DO iLocSide=1,6
    aSide=>aQuad%Side(iLocSide)%sp
    IF(ASSOCIATED(aSide%connection))THEN
      IF(aSide%connection%elem%ind.GT.iQuad)THEN
        aSide%flip=0
      END IF
    END IF
  END DO
END DO


!sanity check
DO iQuad=1,nQuadrants
  aQuad=>Quads(iQuad)%ep
  DO iLocSide=1,6
    IF(aQuad%Side(iLocSide)%sp%flip.LT.0) THEN
      WRITE(*,*) 'flip assignmenti failed, iQuad= ',iQuad,', iLocSide= ',iLocSide 
      STOP
    END IF
  END DO
END DO
END SUBROUTINE BuildMeshFromP4EST


SUBROUTINE EvalP4ESTConnectivity(Conn,nbSide,Flip,Mortar)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=1),INTENT(IN)   :: Conn   ! p4est Side,Flip,Mortar encoding
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)          :: nbSide ! Neighbour side in p4est convention: 0..5
INTEGER,INTENT(OUT)          :: Flip   ! Flip in p4est convention: 0..3
INTEGER,INTENT(OUT)          :: Mortar ! Mortar in p4est convention: 0..3,
                                       ! -1 if conformal, 4 if half-size neighbour
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: tmp
!-----------------------------------------------------------------------------------------------------------------------------------
! The quad_to_quad list stores one value for each local quadrant's face.
! This value is in 0..local_num_quadrants-1 for local quadrants, or in
! local_num_quadrants + (0..ghost_num_quadrants-1) for ghost quadrants.
! The quad_to_face list has equally many entries which are either:
! 1. A value of v = 0..23 indicates one same-size neighbor.
!    This value is decoded as v = r * 6 + nf, where nf = 0..5 is the
!    neigbbor's connecting face number and r = 0..3 is the relative
!    orientation of the neighbor's face, see p8est_connectivity.h.
! 2. A value of v = 24..119 indicates a double-size neighbor.
!    This value is decoded as v = 24 + h * 24 + r * 6 + nf, where
!    r and nf are as above and h = 0..3 is the number of the subface.
! 3. A value of v = -24..-1 indicates four half-size neighbors.
!    In this case the corresponding quad_to_quad index points into the
!    quad_to_half array which stores four quadrant numbers per index,
!    and the orientation of the smaller faces follows from 24 + v.
!    The entries of quad_to_half encode between local and ghost quadrant
!    in the same way as the quad_to_quad values described above.
! A quadrant on the boundary of the forest sees itself and its face number.

SELECT CASE(Conn)
CASE(0:23)   ! 1. conformal neighbour
  nbSide = MOD(Conn,6)       ! 0..5
  Flip   = (Conn-nbSide)/6   ! 0..3
  Mortar = -1
CASE(24:119) ! 2. double-size neighbour
  tmp    = MOD(Conn,24)      ! 0..3
  nbSide = MOD(tmp,6)        ! 0..5  
  Flip   = (tmp-nbSide)/6    ! 0..3
  Mortar = (Conn-tmp-24)/24  ! 0..3 
CASE(-24:-1) ! 3. half-size neighbour
  tmp    = Conn+24
  nbSide = MOD(tmp,6)       ! 0..5
  Flip   = (tmp-nbSide)/6   ! 0..3
  Mortar = 4
CASE DEFAULT
  STOP 'This type of face connectivity does not exist, has to be -24<f<23'
END SELECT

END SUBROUTINE EvalP4ESTConnectivity


FUNCTION GetHFlip(PSide0,PSide1,PFlip)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:P2H_FaceMap,H2P_FaceNodeMap,P2H_FaceNodeMap,P4R,P4Q,P4P
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: PSide0  ! local side ID of HOPEST
INTEGER,INTENT(IN)   :: PSide1  ! P4EST neighbour local side id
INTEGER,INTENT(IN)   :: PFlip   ! Neighbour side in p4est convention: 0..5
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER              :: GetHFlip   ! Neighbour side in p4est convention: 0..5
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: HSide0,PNode0,PNode1 ! p4est Side,Flip,Mortar encoding
!-----------------------------------------------------------------------------------------------------------------------------------
!1. Get CGNS side from P4 side
HSide0=P2H_FaceMap(PSide0)
!2. First node CGNS -> p4est
PNode0=H2P_FaceNodeMap(1,HSide0)
!3. Get oriented node on neighbour side, formula and matrices see paper Burstedde p4est, 2011
PNode1=P4P(P4Q(P4R(PSide0,PSide1),PFlip),PNode0)
!4. P4EST node -> CGNS
GetHFlip=P2H_FaceNodeMap(PNode1,PSide1)

END FUNCTION GetHFlip

SUBROUTINE BuildHOMesh()
!===================================================================================================================================
! uses XGeo High order data from trees and interpolates it to the quadrants 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ALLOCATE(XgeoQuads(3,0:Ngeo,0:Ngeo,0:Ngeo,nQuadrants))
IF(refineLevel.EQ.0)THEN
  XgeoQuads=Xgeo
ELSE
  STOP 'interpolation only for refinelevel 0'
END IF

END SUBROUTINE BuildHOMesh

END MODULE MOD_MeshFromP4EST
