#include "hopest_f.h"

MODULE MODH_VTK
!===================================================================================================================================
! Module for generic data output in vtk xml fromat
!
! WARNING: WriteDataToVTK works only for POSTPROCESSING
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDataToVTK3D
  MODULE PROCEDURE WriteDataToVTK3D
END INTERFACE

PUBLIC::WriteDataToVTK3D
PUBLIC::LinkVTKFiles
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteDataToVTK3D(NPlot,nElems,nVal,VarNames,Coord,Value,FileString)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MODH_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)            :: nElems                  ! Number of output elements
REAL,INTENT(IN)               :: Coord(3,0:NPlot,0:NPlot,0:NPlot,nElems)      ! CoordsVector 
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Value(nVal,0:NPlot,0:NPlot,0:NPlot,nElems)   ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iVal,iElem,Offset,nBytes,nVTKElems,nVTKCells,ivtk=44
INTEGER            :: INT
INTEGER            :: Vertex(8,(NPlot+1)**3*nElems)
INTEGER            :: NPlot_p1_3,NPlot_p1_2,NodeID,NodeIDElem,ElemType
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
REAL(KIND=4)       :: Float
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY (VTU) FILE..."

NPlot_p1_3=(NPlot+1)*(NPlot+1)*(NPlot+1)
NPlot_p1_2=(NPlot+1)*(NPlot+1)

! Line feed character
lf = char(10)

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify file type
nVTKElems=NPlot_p1_3*nElems
nVTKCells=NPlot*NPlot*NPlot*nElems
Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKElems
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=0
WRITE(StrOffset,'(I16)')Offset
DO iVal=1,nVal
  Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'" &
format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF(INT)+nVTKElems*SIZEOF(FLOAT)
  WRITE(StrOffset,'(I16)')Offset
END DO
Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify cell data
Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(INT)+3*nVTKElems*SIZEOF(FLOAT)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(INT)+8*nVTKElems*SIZEOF(INT)
WRITE(StrOffset,'(I16)')Offset
! Offsets
Buffer='        <DataArray type="Int32" Name="offsets" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(INT)+nVTKElems*SIZEOF(INT)
WRITE(StrOffset,'(I16)')Offset
! Elem types
Buffer='        <DataArray type="Int32" Name="types" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
! Prepare append section
Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
! Write leading data underscore
Buffer='_';WRITE(ivtk) TRIM(Buffer)

! Write binary raw data into append section
! Point data
nBytes = nVTKElems*SIZEOF(FLOAT)
DO iVal=1,nVal
  WRITE(ivtk) nBytes,REAL(Value(iVal,:,:,:,:),4)
END DO
! Points
nBytes = nBytes * 3
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coord(:,:,:,:,:),4)
! Connectivity
NodeID=0
NodeIDElem=0
DO iElem=1,nElems
  DO k=1,NPlot
    DO j=1,NPlot
      DO i=1,NPlot
        NodeID=NodeID+1
        !
        Vertex(:,NodeID)=(/                                       &
          NodeIDElem+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P4(CGNS=tecplot standard)
          NodeIDElem+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P1
          NodeIDElem+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P2
          NodeIDElem+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P3
          NodeIDElem+i+   j   *(NPlot+1)+ k   *NPlot_p1_2-1,      & !P8
          NodeIDElem+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P5
          NodeIDElem+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P6
          NodeIDElem+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2-1      /) !P7
      END DO
    END DO
  END DO
  !
  NodeIDElem=NodeIDElem+NPlot_p1_3
END DO
nBytes = 8*nVTKElems*SIZEOF(INT)
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
! Offset
nBytes = nVTKElems*SIZEOF(INT)
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=8,8*nVTKElems,8)
! Elem type
ElemType = 12 ! VTK_HEXAHEDRON
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKElems)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK3D



SUBROUTINE LinkVTKFiles(FileString,OrigFileString,nVal,VarNames)
!===================================================================================================================================
! Linkes VTK data- und mesh-files together
!===================================================================================================================================
! MODULES
USE MODH_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString              ! Output file name
CHARACTER(LEN=*),INTENT(IN)  :: OrigFileString          ! Filename of CGNS data files to be linked
INTEGER,INTENT(IN)           :: nVal                    ! Number of nodal output variables
CHARACTER(LEN=*),INTENT(IN)  :: VarNames(nVal)          ! Names of all variables that will be written out
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVal,iProc,iExt,ivtk=44
CHARACTER(LEN=255) :: OrigFileName,Continuous
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   LINKING VTK FILES..."
! Check if continuous solution is to be linked
iExt=INDEX(FileString,'.continuous',BACK = .TRUE.)
IF(iExt.NE.0)THEN
  Continuous='.continuous'
ELSE
  Continuous=''
END IF
! Line feed character
lf = char(10)

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify file type
Buffer='  <PUnstructuredGrid GhostLevel="0">'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    <PPointData>'//lf;WRITE(ivtk) TRIM(Buffer)
DO iVal=1,nVal
  Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
END DO
Buffer='    </PPointData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    <PCellData> </PCellData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    <PPoints>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    </PPoints>'//lf;WRITE(ivtk) TRIM(Buffer)
! Link files
DO iProc=0,nProcessors-1
  OrigFileName=TRIM(INTSTAMP(OrigFileString,iProc))//TRIM(Continuous)//'.vtu'
  Buffer='    <Piece Source="'//TRIM(OrigFileName)//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
END DO
! Write footer
Buffer='  </PUnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE LinkVTKFiles

END MODULE MODH_VTK
