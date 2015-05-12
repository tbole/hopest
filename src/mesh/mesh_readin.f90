#include "hopest_f.h"

MODULE MODH_Mesh_ReadIn
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MODH_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMeshFromHDF5
  MODULE PROCEDURE ReadMeshFromHDF5
END INTERFACE

INTERFACE ReadGeoFromHDF5
  MODULE PROCEDURE ReadGeoFromHDF5
END INTERFACE

INTERFACE ReadMeshHeader
  MODULE PROCEDURE ReadMeshHeader
END INTERFACE

PUBLIC::ReadMeshFromHDF5
PUBLIC::ReadGeoFromHDF5
PUBLIC::ReadMeshHeader
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadBCs()
!===================================================================================================================================
! Read boundary conditions from data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,ONLY:BoundaryName,BoundaryType,nBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iBC,Offset
!===================================================================================================================================
! Get number of boundary condtions
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
nBCs=HSize(1)
DEALLOCATE(HSize)
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF(HSize(2).NE.nBCs) STOP 'Problem in readBC'
DEALLOCATE(HSize)

ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(BC_SIZE,nBCs))
offset=0
CALL ReadArray('BCNames',1,(/nBCs/),        Offset,1,StrArray    =BoundaryName)
CALL ReadArray('BCType' ,2,(/BC_SIZE,nBCs/),Offset,1,IntegerArray=BoundaryType)

SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(UNIT_StdOut,'(A,A16,A20,A9,A9,A9,A9)')'BOUNDARY CONDITIONS','|','Name','Type','Curved','State','Alpha'
DO iBC=1,nBCs
  SWRITE(UNIT_StdOut,'(A,A33,A20,I9,I9,I9,I9)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(:,iBC)
END DO
SWRITE(UNIT_StdOut,'(132("."))')
END SUBROUTINE ReadBCs


SUBROUTINE SetUserBCs()
!===================================================================================================================================
! The user can redefine boundaries in the ini file. We create the mappings for the boundaries.
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,  ONLY: BoundaryName,BoundaryType,nBCs,nUserBCs
USE MODH_ReadinTools,ONLY: CNTSTR,GETSTR,GETINTARRAY
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: BCMapping(nBCs)
CHARACTER(LEN=255),ALLOCATABLE :: BoundaryNameUser(:)
INTEGER,ALLOCATABLE            :: BoundaryTypeUser(:,:)
INTEGER                        :: iBC,iUserBC
!===================================================================================================================================
! read in boundary conditions, will overwrite BCs from meshfile!
nUserBCs = CNTSTR('BoundaryName',0)
IF(nUserBCs.EQ.0) RETURN

! Read user BC
ALLOCATE(BoundaryNameUser(nUserBCs))
ALLOCATE(BoundaryTypeUser(nUserBCs,2))
DO iBC=1,nUserBCs
  BoundaryNameUser(iBC)   = GETSTR('BoundaryName')
  BoundaryTypeUser(iBC,:) = GETINTARRAY('BoundaryType',2) !(/Type,State/)
END DO

! Override BCs
BCMapping=0
DO iBC=1,nBCs
  DO iUserBC=1,nUserBCs
    IF(INDEX(TRIM(BoundaryNameUser(iUserBC)),TRIM(BoundaryName(iBC))) .NE.0)THEN
      SWRITE(Unit_StdOut,'(A,A)')    ' |     Boundary in HDF file found | ',TRIM(BoundaryName(iBC))
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                            was | ',BoundaryType(1,iBC),BoundaryType(3,iBC)
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                      is set to | ',BoundaryTypeUser(iUserBC,1:2)
      BoundaryType(BC_TYPE ,iBC) = BoundaryTypeUser(iUserBC,1)
      BoundaryType(BC_STATE,iBC) = BoundaryTypeUser(iUserBC,2)
    END IF
  END DO
END DO

SWRITE(UNIT_StdOut,'(132("."))')
DEALLOCATE(BoundaryNameUser,BoundaryTypeUser)
END SUBROUTINE SetUserBCs


SUBROUTINE ReadMeshHeader(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,ONLY: NGeo,useCurveds,nGlobalTrees
USE MODH_Mesh,     ONLY: SetCurvedInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: isMesh
!===================================================================================================================================
CALL CheckIfMesh(FileString,isMesh)
IF(.NOT.isMesh) CALL abort(__STAMP__,&
                           'ERROR: Given file '//TRIM(FileString)//' is no valid mesh file.')

CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)
CALL ReadAttribute(File_ID,'nElems',1,IntegerScalar=nGlobalTrees) !global number of elements
CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
CALL SetCurvedInfo()
useCurveds=.TRUE. ! TODO: maybe implement as optional parameter

CALL GetDataSize(File_ID,'NodeCoords',nDims,HSize)
IF(HSize(2).NE.(NGeo+1)**3*nGlobalTrees) CALL abort(__STAMP__,&
                      'ERROR: Number of nodes in NodeCoords is not consistent with nTrees and NGeo.')
DEALLOCATE(HSize)

CALL readBCs()
CALL setUserBCs()
CALL CloseDataFile() 

END SUBROUTINE ReadMeshHeader


SUBROUTINE ReadMeshFromHDF5(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file and build p4est_connectivity
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars
USE MODH_P4EST_Vars,    ONLY: connectivity,H2P_VertexMap,H2P_FaceMap
USE MODH_P4EST_Binding, ONLY: p4_connectivity_treevertex,p4_build_p4est
USE MODH_P4EST,         ONLY: getHFlip
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: NodeCoordsTmp(:,:,:,:,:)
INTEGER                        :: i,j,k,l
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: FirstSideInd,LastSideInd
INTEGER                        :: nbLocSide,offsetSideID
INTEGER                        :: BCindex
INTEGER                        :: iElem,iTree,ElemID
INTEGER                        :: iNode,jNode,NodeID,SideID
INTEGER                        :: iLocSide,jLocSide
INTEGER                        :: iSide
LOGICAL                        :: oriented
INTEGER                        :: nPeriodicSides 
LOGICAL                        :: fileExists
LOGICAL                        :: doConnection
TYPE(tElem),POINTER            :: aElem
TYPE(tSide),POINTER            :: aSide,bSide
TYPE(tNode),POINTER            :: aNode
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:)
INTEGER                        :: nNodeIDs,nSideIDs
INTEGER,ALLOCATABLE            :: GlobalNodeIDs(:,:,:,:)
INTEGER                        :: C2V(3,8) ! Volume to 1D corner node mapping
! p4est interface
INTEGER                        :: num_vertices
INTEGER                        :: num_trees
INTEGER                        :: num_periodics,iPeriodic,PFlip,HFlip,HFlip_test
INTEGER,ALLOCATABLE            :: tree_to_vertex(:,:)
REAL,ALLOCATABLE               :: vertices(:,:)
INTEGER,ALLOCATABLE            :: JoinFaces(:,:)
!===================================================================================================================================
INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
IF(.NOT.FileExists)  &
    CALL abort(__STAMP__, &
       'readMesh from data file "'//TRIM(FileString)//'" does not exist')

SWRITE(UNIT_stdOut,'(A)')'READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_StdOut,'(132("-"))')

! map variables to use as much of Flexi codebase as possible
offsetElem=0
nElems=nGlobalTrees   !local number of Elements 
nTrees=nGlobalTrees   !local number of Elements 
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------
! Define corner to volume mapping: C2V in CGNS ordering
C2V(:,1)=(/0,   0,   0   /) 
C2V(:,2)=(/NGeo,0,   0   /)
C2V(:,3)=(/NGeo,NGeo,0   /)
C2V(:,4)=(/0,   NGeo,0   /)
C2V(:,5)=(/0,   0,   NGeo/)
C2V(:,6)=(/NGeo,0,   NGeo/)
C2V(:,7)=(/NGeo,NGeo,NGeo/)
C2V(:,8)=(/0,   NGeo,NGeo/)

! get unique node mapping
nNodes=nElems*(NGeo+1)**3
ALLOCATE(GlobalNodeIDs(0:NGeo,0:NGeo,0:NGeo,nElems))
ALLOCATE(XGeo(       3,0:NGeo,0:NGeo,0:NGeo,nElems))

CALL ReadArray('GlobalNodeIDs',1,(/nNodes/),offsetElem*(NGeo+1)**3,1,IntegerArray=GlobalNodeIDs)
CALL ReadArray('NodeCoords',2, (/3,nNodes/),offsetElem*(NGeo+1)**3,2,RealArray=XGeo)

! Unique nodes for elements
CALL ReadAttribute(File_ID,'nUniqueNodes',1,IntegerScalar=nUniqueNodes)
ALLOCATE(UniqueNodes(1:nUniqueNodes))
DO iNode=1,nUniqueNodes
  NULLIFY(UniqueNodes(iNode)%np)
END DO

num_vertices=0
DO iElem=1,nElems
  DO iNode=1,8
    i=C2V(1,iNode);j=C2V(2,iNode);k=C2V(3,iNode);
    l=GlobalNodeIDs(i,j,k,iElem)
    IF(.NOT.ASSOCIATED(UniqueNodes(l)%np))THEN
      ALLOCATE(UniqueNodes(l)%np)
      num_vertices=num_vertices+1
    END IF
    UniqueNodes(l)%np%x  =XGeo(:,i,j,k,iElem)
    UniqueNodes(l)%np%ind=l
    UniqueNodes(l)%np%tmp=-1
  END DO
END DO

! remap to linear mesh, use only corner nodes
!IF(.NOT.useCurveds)THEN
!  ALLOCATE(XGeoTmp(3,0:NGeo,0:NGeo,0:NGeo,nElems))
!  XGeoTmp=XGeo
!  DEALLOCATE(XGeo)
!  ALLOCATE(XGeo(3,0:1,0:1,0:1,nElems))
!  XGeo(:,0,0,0,:)=XGeoTmp(:,0,   0,   0,   :)
!  XGeo(:,1,0,0,:)=XGeoTmp(:,NGeo,0,   0,   :)
!  XGeo(:,0,1,0,:)=XGeoTmp(:,0,   NGeo,0,   :)
!  XGeo(:,1,1,0,:)=XGeoTmp(:,NGeo,NGeo,0,   :)
!  XGeo(:,0,0,1,:)=XGeoTmp(:,0,   0,   NGeo,:)
!  XGeo(:,1,0,1,:)=XGeoTmp(:,NGeo,0,   NGeo,:)
!  XGeo(:,0,1,1,:)=XGeoTmp(:,0,   NGeo,NGeo,:)
!  XGeo(:,1,1,1,:)=XGeoTmp(:,NGeo,NGeo,NGeo,:)
!  DEALLOCATE(XGeoTmp)
!  NGeo=1
!  nNodes=nElems*(NGeo+1)**3
!ENDIF

!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
FirstElemInd=offsetElem+1
LastElemInd=offsetElem+nElems
ALLOCATE(Trees(                 FirstElemInd:LastElemInd))
ALLOCATE(ElemInfo(ElemInfoSize,FirstElemInd:LastElemInd))
CALL ReadArray('ElemInfo',2,(/ElemInfoSize,nElems/),offsetElem,2,IntegerArray=ElemInfo)

DO iElem=FirstElemInd,LastElemInd
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  Trees(iElem)%ep=>GETNEWELEM()
  aElem=>Trees(iElem)%ep
  aElem%Ind    = iElem
  aElem%Type   = ElemInfo(ELEM_Type,iElem)
  aElem%Zone   = ElemInfo(ELEM_Zone,iElem)
  ! WARNING: THIS IS CARTESIAN ORDERING NOT CGNS ORDERING AS IN FLEXI !!!
  DO iNode=1,8
    i=C2V(1,iNode);j=C2V(2,iNode);k=C2V(3,iNode);
    aElem%Node(iNode)%np=>UniqueNodes(GlobalNodeIDs(i,j,k,iElem))%np
  END DO
  CALL createSides(aElem)
END DO

DEALLOCATE(GlobalNodeIDs)

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------

offsetSideID=ElemInfo(ELEM_FirstSideInd,FirstElemInd) ! hdf5 array starts at 0-> -1  
nSideIDs=ElemInfo(ELEM_LastSideInd,LastElemInd)-ElemInfo(ELEM_FirstSideInd,FirstElemInd)
!read local SideInfo from data file 
FirstSideInd=offsetSideID+1
LastSideInd=offsetSideID+nSideIDs
ALLOCATE(SideInfo(SideInfoSize,FirstSideInd:LastSideInd))
CALL ReadArray('SideInfo',2,(/SideInfoSize,nSideIDs/),offsetSideID,2,IntegerArray=SideInfo)


DO iElem=FirstElemInd,LastElemInd
  aElem=>Trees(iElem)%ep
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  !build up sides of the element according to CGNS standard
  ! assign flip
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    aSide%Elem=>aElem
    oriented=(Sideinfo(SIDE_ID,iSide).GT.0)
    aSide%Ind=ABS(SideInfo(SIDE_ID,iSide))
    IF(oriented)THEN !oriented side
      aSide%flip=0
    ELSE !not oriented
      aSide%flip=MOD(Sideinfo(SIDE_Flip,iSide),10)
      IF((aSide%flip.LT.0).OR.(aSide%flip.GT.4)) STOP 'NodeID doesnt belong to side'
    END IF

    ! Check if mortar element
    ElemID=SideInfo(SIDE_nbElemID,iSide) !IF nbElemID <0, this marks a mortar master side. 
                                         ! The number (-1,-2,-3) is the Type of mortar
    IF(ElemID.LT.0)THEN ! mortar Sides attached!
      CALL abort(__STAMP__, &
           'Only conforming meshes in readin.')
    END IF
  END DO !i=1,locnSides
END DO !iElem

 
! build up side connection 
DO iElem=FirstElemInd,LastElemInd
  aElem=>Trees(iElem)%ep
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    elemID  = SideInfo(SIDE_nbElemID,iSide)
    BCindex = SideInfo(SIDE_BCID,iSide)
    doConnection=.TRUE. ! for periodic sides if BC is reassigned as non periodic
    IF(BCindex.NE.0)THEN !BC
      aSide%BCindex = BCindex
      IF((BoundaryType(BC_TYPE,aSide%BCindex).NE.1).AND.&
         (BoundaryType(BC_TYPE,aSide%BCindex).NE.100))THEN ! Reassignement from periodic to non-periodic
        doConnection=.FALSE.
        aSide%flip  =0
      END IF
    ELSE
      aSide%BCindex = 0
    END IF

    IF(.NOT.doConnection) CYCLE
    IF(ASSOCIATED(aSide%connection)) CYCLE

    ! check if neighbor on local proc or MPI connection
    IF((elemID.LE.LastElemInd).AND.(elemID.GE.FirstElemInd))THEN !local
      nbLocSide=Sideinfo(SIDE_Flip,iSide)/10
      IF((nbLocSide.LT.1).OR.(nbLocSide.GT.6))&
        CALL abort(__STAMP__,'SideInfo: Index of local side must be between 1 and 6!')
      bSide=>Trees(elemID)%ep%Side(nbLocSide)%sp
      aSide%connection=>bSide
      bSide%connection=>aSide
      IF(bSide%ind.NE.aSide%ind)&
        CALL abort(__STAMP__,&
        'SideInfo: Index of side and neighbor side have to be identical!')
    ELSE !MPI
#ifdef MPI
      aSide%connection=>GETNEWSIDE()            
      aSide%connection%flip=aSide%flip
      aSide%connection%Elem=>GETNEWELEM()
      aSide%NbProc = ELEMIPROC(elemID)
#else
      CALL abort(__STAMP__, &
        ' elemID of neighbor not in global Elem list ')
#endif
    END IF
  END DO !iLocSide 
END DO !iElem

DEALLOCATE(ElemInfo,SideInfo)
CALL CloseDataFile()


!----------------------------------------------------------------------------------------------------------------------------
!                  P4EST MESH CONNECTIVITY
!----------------------------------------------------------------------------------------------------------------------------
! should be replaced by connectivity information ?
! needs unique corner nodes for mesh connectivity

num_vertices=0
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iNode=1,8
    aNode=>aElem%Node(iNode)%np
    IF(aNode%tmp.EQ.-1)THEN
      num_vertices=num_vertices+1
      aElem%Node(iNode)%np%tmp=num_vertices
    END IF
  END DO
END DO !iTree

ALLOCATE(Vertices(3,num_vertices))
DO iNode=1,nUniqueNodes
  aNode=>UniqueNodes(iNode)%np
  IF(ASSOCIATED(aNode))THEN ! only corner nodes are associated
    IF(aNode%tmp.GT.0) Vertices(:,aNode%tmp)=aNode%x
  END IF
END DO

num_trees=nTrees
ALLOCATE(tree_to_vertex(8,num_trees))
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iNode=1,8
    tree_to_vertex(iNode,iTree)=aElem%Node(H2P_VertexMap(iNode)+1)%np%tmp-1
  END DO
END DO

!periodic Boundaries
num_periodics=0
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%BCIndex.EQ.0) CYCLE ! NO Boundary Condition
    IF((BoundaryType(BC_TYPE,aSide%BCIndex).EQ.1).AND.(aSide%flip.EQ.0))THEN !periodic side 
      num_periodics=num_periodics+1
    END IF
  END DO !iLocSide
END DO !iTree


IF(num_periodics.GT.0) THEN
  ALLOCATE(JoinFaces(5,num_periodics))
  DO iTree=1,nTrees
    aElem=>Trees(iTree)%ep
    DO iLocSide=1,6
      aElem%Side(iLocSide)%sp%tmp=H2P_FaceMap(iLocSide)  !local Face ID in p4est
    END DO
  END DO
  
  iperiodic=0
  DO iTree=1,nTrees
    aElem=>Trees(iTree)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      IF(aSide%BCIndex.EQ.0) CYCLE ! NO Boundary Condition
      IF((BoundaryType(BC_TYPE,aSide%BCIndex).EQ.1).AND.(aSide%flip.EQ.0))THEN !periodic side 
        HFlip=aSide%connection%flip
        iperiodic=iperiodic+1
        bSide=>aSide%connection
        IF(aSide%tmp.GT.bSide%tmp)THEN
          aSide=>aSide%connection
        END IF
        bSide=>aSide%connection
        JoinFaces(1,iPeriodic)=aSide%Elem%ind-1  !treeID of face with smaller p4est locfaceID
        JoinFaces(2,iPeriodic)=bSide%Elem%ind-1  ! neighbor tree id
        JoinFaces(3,iPeriodic)=aSide%tmp         ! p4est locSideID
        JoinFaces(4,iPeriodic)=bSide%tmp         ! p4est neighbor locSideID
        DO PFlip=0,3
          Hflip_test=getHflip(aSide%tmp,bSide%tmp,PFlip)
          IF(HFlip_test.EQ.HFlip) EXIT
        END DO
        JoinFaces(5,iPeriodic)=PFlip
      END IF
    END DO !iLocSide
  END DO !iTree
END IF !num_periodics>0

CALL p4_connectivity_treevertex(num_vertices,num_trees,vertices,tree_to_vertex, &
                                   num_periodics,JoinFaces,connectivity)

DEALLOCATE(Vertices,tree_to_vertex)
IF(num_periodics.GT.0) DEALLOCATE(JoinFaces) 
 

! COUNT SIDES

nBCSides=0
nSides=0
nPeriodicSides=0
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    aSide%tmp=0 
  END DO !iLocSide
END DO !iTree
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp

    IF(aSide%tmp.EQ.0)THEN
      nSides=nSides+1
      aSide%tmp=-1 !used as marker
      IF(ASSOCIATED(aSide%connection)) aSide%connection%tmp=-1
      IF(aSide%BCindex.NE.0)THEN !side is BC or periodic side
        IF(ASSOCIATED(aSide%connection))THEN
          nPeriodicSides=nPeriodicSides+1
        ELSE
          nBCSides=nBCSides+1
        END IF
      END IF
    END IF
  END DO !iLocSide
END DO !iTree


WRITE(*,*)'-------------------------------------------------------'
WRITE(*,'(A22,I8)' )'NGeo:',NGeo
WRITE(*,'(A22,X7L)')'useCurveds:',useCurveds
WRITE(*,'(A22,I8)' )'nTrees:',nTrees
WRITE(*,'(A22,I8)' )'nNodes:',nNodes
WRITE(*,'(A22,I8)' )'nSides:',nSides
WRITE(*,'(A22,I8)' )'nBCSides:',nBCSides
WRITE(*,'(A22,I8)' )'nPeriodicSides:',nPeriodicSides
WRITE(*,*)'-------------------------------------------------------'

END SUBROUTINE ReadMeshFromHDF5



SUBROUTINE ReadGeoFromHDF5(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars, ONLY: nTrees,XGeo,Ngeo,nNodes,HexMap
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k
INTEGER                        :: iTree
INTEGER                        :: nNodeIDs
INTEGER                        :: FirstNodeInd,LastNodeInd
LOGICAL                        :: fileExists
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE               :: NodeCoords(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
IF(.NOT.FileExists)  &
    CALL abort(__STAMP__, &
       'readMesh from data file "'//TRIM(FileString)//'" does not exist')

SWRITE(UNIT_stdOut,'(A)')'READ GEOMETRY DATA FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_StdOut,'(132("-"))')
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)

CALL GetDataSize(File_ID,'NodeCoords',nDims,HSize)
nNodes=HSize(1) !global number of unique nodes
DEALLOCATE(HSize)

!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
ALLOCATE(ElemInfo(1:nTrees,ElemInfoSize))
CALL ReadArray('ElemInfo',2,(/nTrees,ElemInfoSize/),0,1,IntegerArray=ElemInfo)

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from data file 
nNodeIDs=ElemInfo(nTrees,ELEM_LastNodeInd)-ElemInfo(1,ELEM_FirstNodeInd)
ALLOCATE(NodeInfo(1:nNodeIDs))
CALL ReadArray('NodeInfo',1,(/nNodeIDs/),0,1,IntegerArray=NodeInfo)

! get physical coordinates
ALLOCATE(NodeCoords(nNodes,3))

CALL ReadArray('NodeCoords',2,(/nNodes,3/),0,1,RealArray=NodeCoords)

CALL CloseDataFile() 

ALLOCATE(Xgeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nTrees))
IF(Ngeo.EQ.1)THEN !use the corner nodes
  DO iTree=1,nTrees
    firstNodeInd=ElemInfo(iTree,ELEM_FirstNodeInd) !first index -1 in NodeInfo
    Xgeo(:,0,0,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+1),:)
    Xgeo(:,1,0,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+2),:)
    Xgeo(:,1,1,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+3),:)
    Xgeo(:,0,1,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+4),:)
    Xgeo(:,0,0,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+5),:)
    Xgeo(:,1,0,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+6),:)
    Xgeo(:,1,1,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+7),:)
    Xgeo(:,0,1,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+8),:)
 END DO !iTree=1,nTrees
ELSE
  DO iTree=1,nTrees
    firstNodeInd=ElemInfo(iTree,ELEM_FirstNodeInd) !first index -1 in NodeInfo
    lastNodeInd=ElemInfo(iTree,ELEM_LastNodeInd)
    IF(LastNodeInd-firstNodeInd-14.NE.(Ngeo+1)**3) STOP 'Problem with curved'
    firstNodeInd=firstNodeInd +8
    DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
      Xgeo(:,i,j,k,iTree)=NodeCoords(NodeInfo(firstNodeInd+HexMap(i,j,k)),:)
    END DO ; END DO ; END DO 
 END DO !iTree=1,nTrees
END IF

DEALLOCATE(ElemInfo,NodeInfo,NodeCoords)

END SUBROUTINE ReadGeoFromHDF5



SUBROUTINE CheckIfMesh(MeshFile_in,isMeshFile)
!===================================================================================================================================
! Check if the file is a mesh file
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
USE MODH_HDF5_Input
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),INTENT(IN) :: MeshFile_in
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: isMeshFile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: FileVersion
INTEGER            :: i,nDataSets
CHARACTER(LEN=255) :: DataSets(9)
!===================================================================================================================================
CALL OpenDataFile(MeshFile_in,create=.FALSE.,single=.FALSE.)
nDataSets=9
DataSets(1)= 'BCNames'
DataSets(2)='BCType'
DataSets(3)='ElemBarycenters'
DataSets(4)='ElemCounter'
DataSets(5)='ElemInfo'
DataSets(6)='ElemWeight'
DataSets(7)='NodeCoords'
DataSets(8)='GlobalNodeIDs'
DataSets(9)='SideInfo'
CALL DatasetExists(File_ID,'Version',isMeshFile,attrib=.TRUE.)
IF(isMeshFile)THEN
  CALL ReadAttribute(File_ID,'Version',1,RealScalar=FileVersion)
  IF(FileVersion.LT.0.999) isMeshFile=.FALSE.
END IF
DO i=1,nDataSets
  IF(.NOT.isMeshFile) EXIT
  CALL DatasetExists(File_ID,DataSets(i),isMeshFile)
END DO! i=1,nDataSets
CALL CloseDataFile()
END SUBROUTINE CheckIfMesh



END MODULE MODH_Mesh_ReadIn
