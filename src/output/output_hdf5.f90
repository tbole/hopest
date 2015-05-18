#include "hopest_f.h"
MODULE MODH_Output_HDF5
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER(HID_T)                 :: File_ID
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteMeshToHDF5
  MODULE PROCEDURE WriteMeshToHDF5
END INTERFACE

PUBLIC::WriteMeshToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteMeshToHDF5(FileString)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MODH_Mesh_Vars
USE MODH_IO_HDF5
USE MODH_HDF5_output
USE MODH_ChangeBasis, ONLY:ChangeBasis3D
USE MODH_P4EST_Vars,  ONLY:QuadToTree,QuadCoords,QuadLevel,sIntSize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

TYPE(tElem),POINTER            :: Elem
TYPE(tElem),POINTER            :: master
TYPE(tSide),POINTER            :: Side
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:),GlobalNodeIDs(:)
REAL,ALLOCATABLE               :: ElemBary(:,:)
REAL,ALLOCATABLE               :: ElemWeight(:)
REAL,ALLOCATABLE               :: xiMinMax(:,:,:)
REAL                           :: length
INTEGER                        :: iElem,i
INTEGER                        :: NodeID,iNode,iType
INTEGER                        :: iSide,SideID,iLocSide,iMortar
INTEGER                        :: ElemCounter(2,11)
INTEGER                        :: nSideIDs,nNodeIDs
INTEGER                        :: nTotalSides
INTEGER                        :: locnSides,locnNodes,offsetID
!===================================================================================================================================
WRITE(*,'(132("~"))')
WRITE(*,'(A)')' WRITE MESH TO HDF5 FILE... ' // TRIM(FileString) 

!set all side indices =0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    Side%ind=0
  END DO !iLocSide=1,6
END DO !iElem=1,nElems

! count Elements , unique sides 
nSideIDs=0 !number of unique side IDs (side and side%connection have the same sideID)


DO iElem=1,nElems
  Elem=>Elems(iElem)%ep

  ! Count sides
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF(Side%ind.EQ.0) THEN
      IF(Side%MortarType.EQ.0)THEN
        nSideIDs=nSideIDs+1
        Side%ind=-88888
        IF(ASSOCIATED(Side%connection))THEN      
          IF(Side%connection%ind.EQ.0) nSideIDs=nSideIDs-1 ! count inner and periodic sides only once 
        END IF
      ELSEIF(Side%MortarType.GT.0)THEN
        nSideIDs=nSideIDs+1
        Side%ind=-88888
        DO iMortar=1,Side%nMortars
          IF(Side%MortarSide(iMortar)%sp%ind.EQ.0)THEN
            nSideIDs=nSideIDs+1
            Side%MortarSide(iMortar)%sp%ind=-88888
          END IF 
        END DO !iMortar
      ELSE
        nSideIDs=nSideIDs+1
        Side%ind=-88888
      END IF
    END IF
  END DO
END DO


!set unique nodes and Side Indices
SideID=0
NodeID=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  Elem%ind=iElem
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF(side%ind.EQ.-88888) THEN  ! assign side ID only for non MPI sides and lower MPI sides
      IF(Side%MortarType.EQ.0)THEN
        SideID=SideID+1
        Side%ind=SideID
        IF(ASSOCIATED(Side%connection))THEN      
          IF(Side%connection%ind.EQ.-88888) Side%connection%ind=SideID ! count inner and periodic sides only once 
        END IF
      ELSEIF(Side%MortarType.GT.0)THEN
        SideID=SideID+1
        Side%ind=SideID
        DO iMortar=1,Side%nMortars
          IF(Side%MortarSide(iMortar)%sp%ind.EQ.-88888)THEN
            SideID=SideID+1
            Side%MortarSide(iMortar)%sp%ind=SideID
          END IF 
        END DO !iMortar
      ELSE
        SideID=SideID+1
        Side%ind=SideID
      END IF
    END IF
  END DO !iLocSide
END DO !Elem

DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF((Side%ind.LE.0).OR.(Side%ind.GT.nSideIds)) STOP 'Problem with sideID assigment'
  END DO !iLocSide
END DO !Elem

IF(SideID.NE.nSideIDs) STOP' problem: SideID <> nSideIDs'

! start output
CALL OpenHDF5File(FileString,create=.TRUE.,single=.TRUE.)  

!-----------------------------------------------------------------
!attributes 
!-----------------------------------------------------------------
IF(Ngeo_out.EQ.1) THEN
  useCurveds=.FALSE.
  nCurvedNodes=0
ELSE
  nCurvedNodes=(Ngeo_out+1)**3
END IF

CALL WriteAttributeToHDF5(File_ID,'Ngeo',1,IntScalar=Ngeo_out)

!-----------------------------------------------------------------
! WRITE BC 
!-----------------------------------------------------------------
CALL WriteArrayToHDF5(File_ID,'BCNames',1,(/nBCs/),(/nBCs/),(/0/),&
                      collective=.FALSE.,StrArray=BoundaryName)
CALL WriteArrayToHDF5(File_ID,'BCType' ,2,(/BC_SIZE,nBCs/),(/BC_SIZE,nBCs/),(/0,0/),&
                      collective=.FALSE.,IntArray=BoundaryType)

!-----------------------------------------------------------------
! Barycenters
!-----------------------------------------------------------------

ALLOCATE(ElemBary(3,nElems))
DO iElem=1,nElems
  ElemBary(1,iElem)=SUM(XGeoElem(1,:,:,:,iElem))
  ElemBary(2,iElem)=SUM(XGeoElem(2,:,:,:,iElem))
  ElemBary(3,iElem)=SUM(XGeoElem(3,:,:,:,iElem))
END DO !iElem=1,nElem
ElemBary=ElemBary*(1./(Ngeo_out+1.)**3.)
CALL WriteArrayToHDF5(File_ID,'ElemBarycenters',2,(/3,nElems/),(/3,nElems/),(/0,0/),&
                      collective=.FALSE.,RealArray=ElemBary)
DEALLOCATE(ElemBary)

!-----------------------------------------------------------------
! WRITE NodeCoords  for each element !!!! (multiple nodes!!!)
!-----------------------------------------------------------------
!transform to equidistant nodes (overwrite!!!):
DO iElem=1,nElems
  CALL ChangeBasis3D(3,Ngeo_out,Ngeo_out,Vdm_CL_EQ_out,XgeoElem(:,:,:,:,iElem),XgeoElem(:,:,:,:,iElem))
END DO
nNodeIDs=(Ngeo_out+1)**3*nElems
CALL WriteArrayToHDF5(File_ID,'NodeCoords',2,(/3,nNodeIDs/),(/3,nNodeIDs/),(/0,0/),&
                      collective=.FALSE.,RealArray=XGeoElem)

!-----------------------------------------------------------------
!fill ElementInfo. 
!-----------------------------------------------------------------

ALLOCATE(ElemInfo(ElemInfoSize,nElems))
ElemInfo=0
Elemcounter=0
Elemcounter(1,:)=(/104,204,105,115,205,106,116,206,108,118,208/)
iNode  = 0 
iSide  = 0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  locnSides=0
  ! for element sides
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    locnSides = locnSides+1

    !Mortar
    SELECT CASE(Side%MortarType)
    CASE(0) !do nothing
    CASE(1)
      locnSides=locnSides+4
    CASE(2,3)
      locnSides=locnSides+2
    END SELECT

  END DO !iLocSide=1,6
  ElemInfo(ELEM_Type        ,iElem) = Elem%Type          ! Element Type
  ElemInfo(ELEM_Zone        ,iElem) = Elem%Zone          ! Zone Number
  ElemInfo(ELEM_FirstSideInd,iElem) = iSide              ! first index -1 in SideInfo
  ElemInfo(ELEM_LastSideInd ,iElem) = iSide+locnSides    ! last index in SideInfo
  ElemInfo(ELEM_FirstNodeInd,iElem) = iNode              ! first index -1 in NodeInfo
  ElemInfo(ELEM_LastNodeInd ,iElem) = iNode+nCurvedNodes ! last index in NodeInfo
  iNode = iNode + nCurvedNodes
  iSide = iSide + locnSides
  !approximate weight: locnNodes
  CALL AddToElemCounter(Elem%type,ElemCounter)
END DO !iElem
nTotalSides = iSide

CALL WriteArrayToHDF5(File_ID,'ElemInfo',2,(/ElemInfoSize,nElems/),(/ElemInfoSize,nElems/),&
                              (/0,0/),collective=.FALSE.,IntArray=ElemInfo)
DEALLOCATE(ElemInfo)

CALL WriteArrayToHDF5(File_ID,'ElemCounter',2,(/2,11/),(/2,11/),(/0,0/),&
                      collective=.FALSE.,IntArray=ElemCounter)
WRITE(*,*)'Mesh statistics:'
WRITE(*,*)'Element Type | number of elements'
DO i=1,11
  WRITE(*,'(I4,A,I8)') Elemcounter(1,i),'        | ',Elemcounter(2,i)
END DO

!-----------------------------------------------------------------
! element weights
!-----------------------------------------------------------------
!WRITE ElemWeight,into (1,nElems)  
ALLOCATE(ElemWeight(1:nElems))
ElemWeight=1.
CALL WriteArrayToHDF5(File_ID,'ElemWeight',1,(/nElems/),(/nElems/),(/0/),&
                              collective=.FALSE.,RealArray=ElemWeight)
DEALLOCATE(ElemWeight)

!-----------------------------------------------------------------
!fill SideInfo
!-----------------------------------------------------------------

ALLOCATE(SideInfo(SideInfoSize,nTotalSides)) 
SideInfo=0 
iSide=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    iSide=iSide+1

    i=MERGE(104,204,Ngeo_out.EQ.1)                              ! 104:bilinear,204:curved
    SideInfo(SIDE_Type,iSide)=MERGE(i,-i,Side%MortarType.GE.0)  ! mark side with mortar neighbour
    SideInfo(SIDE_ID,  iSide)=MERGE(Side%ind,-Side%ind,Side%flip.EQ.0) ! neg. sideID for slaves
    SideInfo(SIDE_Flip,iSide)=Side%flip
    SideInfo(SIDE_BCID,iSide)=Side%BCIndex
    IF(ASSOCIATED(Side%Connection))THEN
      SideInfo(SIDE_nbElemID,iSide)=Side%Connection%Elem%ind    ! neighbour Element ID
    END IF

    !MORTAR
    IF(Side%MortarType.GT.0)THEN
      IF(ASSOCIATED(Side%Connection)) CALL abort(__STAMP__,&
                                                 'Mortar master with connection is not allowed')
      IF(Side%flip.NE.0) STOP 'Problem with flip on mortar'
      SideInfo(SIDE_nbElemID,iSide)=-Side%MortarType !marker for attached mortar sides
      DO iMortar=1,Side%nMortars
        iSide=iSide+1
        SideInfo(SIDE_Type,    iSide)= MERGE(104,204,Ngeo_out.EQ.1)
        SideInfo(SIDE_ID,      iSide)= Side%MortarSide(iMortar)%sp%ind       ! small are always master
        SideInfo(SIDE_nbElemID,iSide)= Side%MortarSide(iMortar)%sp%Elem%ind  ! neighbour Element ID
        SideInfo(SIDE_BCID,    iSide)= Side%MortarSide(iMortar)%sp%BCIndex
      END DO
    ENDIF
  END DO !iLocSide=1,6
END DO !iElem=1,nElems

CALL WriteArrayToHDF5(File_ID,'SideInfo',2,(/SideInfoSize,nTotalSides/),(/SideInfoSize,nTotalSides/),&
                      (/0,0/),collective=.FALSE.,IntArray=SideInfo)
DEALLOCATE(SideInfo)

!-----------------------------------------------------------------
!fill GlobalNodeIDs
!-----------------------------------------------------------------

! TODO: this is only dummy GlobalNodeID array, without unique nodes
ALLOCATE(GlobalNodeIDs(1:nNodeIDs))
DO iNode=1,nNodeIDs
  GlobalNodeIDs(iNode)=iNode
END DO

CALL WriteArrayToHDF5(File_ID,'GlobalNodeIDs',1,(/nNodeIDs/),(/nNodeIDs/),(/0/),&
                      collective=.FALSE.,IntArray=GlobalNodeIDs)
DEALLOCATE(GlobalNodeIDs)

!-----------------------------------------------------------------
! Additional information which is only defined for mortar meshes
!-----------------------------------------------------------------
nNodeIDs=(Ngeo+1)**3*nTrees
CALL WriteAttributeToHDF5(File_ID,'isMortarMesh',1,IntScalar=1)
CALL WriteAttributeToHDF5(File_ID,'NgeoTree',1,IntScalar=Ngeo)
CALL WriteAttributeToHDF5(File_ID,'nTrees',1,IntScalar=nTrees)
CALL WriteAttributeToHDF5(File_ID,'nNodesTree',1,IntScalar=nNodeIDs)
!-----------------------------------------------------------------
! WRITE TreeCoords for each tree !!!! (multiple nodes!!!)
!-----------------------------------------------------------------
CALL WriteArrayToHDF5(File_ID,'TreeCoords',2,(/3,nNodeIDs/),(/3,nNodeIDs/),(/0,0/),&
                      collective=.FALSE.,RealArray=XGeo)
!-----------------------------------------------------------------
! WRITE QuadToTree for each quadrant !!!!
!-----------------------------------------------------------------
QuadToTree=QuadToTree+1 ! 0-based to 1-based
CALL WriteArrayToHDF5(File_ID,'ElemToTree',1,(/nElems/),(/nElems/),(/0/),&
                      collective=.FALSE.,IntArray=QuadToTree)
QuadToTree=QuadToTree-1 ! revert
!-----------------------------------------------------------------
! WRITE xiMinMax for each quadrant !!!!
!-----------------------------------------------------------------
ALLOCATE(xiMinMax(3,2,nElems))
DO iElem=1,nElems
  xiMinMax(:,1,iElem)=-1.+2.*REAL(QuadCoords(:,iElem))*sIntSize
  length=2./REAL(2**QuadLevel(iElem))
  xiMinMax(:,2,iElem)=xiMinMax(:,1,iElem)+length
END DO
CALL WriteArrayToHDF5(File_ID,'xiMinMax',3,(/3,2,nElems/),(/3,2,nElems/),(/0,0,0/),&
                      collective=.FALSE.,RealArray=xiMinMax)
DEALLOCATE(xiMinMax)
!-----------------------------------------------------------------
CALL CloseHDF5File()

WRITE(*,'(A)')' DONE WRITING MESH.'
WRITE(*,'(132("~"))')
END SUBROUTINE WriteMeshToHDF5



SUBROUTINE AddtoElemCounter(elemtype,elemCounter) 
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: ElemType
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)    :: ElemCounter(2,11)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(ElemType)
    CASE(104) !linear tet
      elemcounter(2,1)=elemcounter(2,1)+1
    CASE(204) !spli tet
      elemcounter(2,2)=elemcounter(2,2)+1
    CASE(105) !line pyr
      elemcounter(2,3)=elemcounter(2,3)+1
    CASE(115) !non-near pyr
      elemcounter(2,4)=elemcounter(2,4)+1
    CASE(205) !spli pyr
      elemcounter(2,5)=elemcounter(2,5)+1
    CASE(106) !line prism
      elemcounter(2,6)=elemcounter(2,6)+1
    CASE(116) !non-near prism
      elemcounter(2,7)=elemcounter(2,7)+1
    CASE(206) !spli prism
      elemcounter(2,8)=elemcounter(2,8)+1
    CASE(108) !line hex
      elemcounter(2,9)=elemcounter(2,9)+1
    CASE(118) !non-linear hex
      elemcounter(2,10)=elemcounter(2,10)+1
    CASE(208) !splinhex
      elemcounter(2,11)=elemcounter(2,11)+1
    CASE DEFAULT
      STOP 'elem type not defined in elemcounter'
  END SELECT
END SUBROUTINE AddToElemCounter


END MODULE MODH_Output_HDF5
