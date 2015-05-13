#include "hopest_f.h"

MODULE MODH_HDF5_output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MODH_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteHDF5Header
  MODULE PROCEDURE WriteHDF5Header
END INTERFACE

!INTERFACE WriteArrayToHDF5
!  MODULE PROCEDURE WriteArrayToHDF5
!END INTERFACE

INTERFACE WriteAttributeToHDF5
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

PUBLIC :: WriteArrayToHDF5
PUBLIC :: WriteHDF5Header
PUBLIC :: WriteAttributeToHDF5
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteHDF5Header(FileType_in,File_ID)
!===================================================================================================================================
! Subroutine to write a distinct file header to each HDF5 file
!===================================================================================================================================
! MODULES
USE MODH_Output_Vars,ONLY:ProgramName,FileVersion,ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)              :: FileType_in
INTEGER(HID_T),INTENT(IN)                :: File_ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files
! Attributes are program name, file type identifier, project name and version number
CALL WriteAttributeToHDF5(File_ID,'Program'     ,1,StrScalar=(/TRIM(ProgramName)/))
CALL WriteAttributeToHDF5(File_ID,'File_Type'   ,1,StrScalar=(/TRIM(FileType_in)/))
CALL WriteAttributeToHDF5(File_ID,'Project_Name',1,StrScalar=(/TRIM(ProjectName)/))
CALL WriteAttributeToHDF5(File_ID,'File_Version',1,RealScalar=FileVersion)
END SUBROUTINE WriteHDF5Header


SUBROUTINE WriteArrayToHDF5(File_ID,DataSetName,rank,nValGlobal,nVal,offset,&
                            collective,resizeDim,chunkSize,&
                            RealArray,IntArray,StrArray)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MODH_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T),INTENT(IN)     :: File_ID
CHARACTER(LEN=*),INTENT(IN)   :: DataSetName
INTEGER,INTENT(IN)            :: rank             ! number of dimensions of the array
INTEGER,INTENT(IN)            :: nValGlobal(rank) ! max size of array in offset dimension
INTEGER,INTENT(IN)            :: nVal(rank)       ! size of complete (local) array to write
INTEGER,INTENT(IN)            :: offset(rank)     ! offset =0, start at beginning of the array
LOGICAL,INTENT(IN)            :: collective       ! use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL   :: resizeDim(rank)  ! specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL   :: chunkSize(rank)  ! specify chunksize
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray(PRODUCT(nVal))
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,Type_ID,dsetparams
INTEGER(HSIZE_T)               :: Dimsf(Rank),OffsetHDF(Rank),nValMax(Rank)
INTEGER(SIZE_T)                :: SizeSet=255
LOGICAL                        :: chunky
TYPE(C_PTR)                    :: buf
!INTEGER(SIZE_T)                :: tmp1,tmp2,typeoffset,typesize
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'
! Get global array size, always first dimension!!
nValMax=nValGlobal
chunky=.FALSE.
CALL H5PCREATE_F(H5P_DATASET_CREATE_F,dsetparams,iError)
IF(PRESENT(chunkSize))THEN
  chunky=.TRUE.
  Dimsf=chunkSize
  CALL H5PSET_CHUNK_F(dsetparams,rank,dimsf,iError)
END IF

IF(PRESENT(resizeDim))THEN
  IF(.NOT.PRESENT(chunkSize))&
    CALL abort(__STAMP__,&
               'Chunk size has to be specified when using resizable arrays.')
  nValMax = MERGE(H5S_UNLIMITED_F,nValMax,resizeDim)
END IF

! Create the dataset with default properties.
IF(PRESENT(RealArray)) Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, SizeSet, iError)
END IF

Dimsf = nValGlobal ! we need the global array size
CALL H5ESET_AUTO_F(0,iError)
CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
IF(iError.NE.0)THEN ! does not exist
  ! Create the data space for the  dataset.
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError, nValMax)
  CALL H5DCREATE_F(File_ID, TRIM(DataSetName), Type_ID, FileSpace, DSet_ID,iError,dsetparams)
CALL H5SCLOSE_F(FileSpace, iError)
END IF
CALL H5ESET_AUTO_F(1,iError)
IF(chunky)THEN
  CALL H5DSET_EXTENT_F(DSet_ID,Dimsf,iError) ! if resizable then dataset may need to be extended
END IF

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
OffsetHDF = Offset
! Create the data space in the memory
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
ELSE
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SSELECT_NONE_F(FileSpace,iError)
ELSE
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, OffsetHDF, Dimsf, iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#ifdef MPI
IF(collective)THEN
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F,  iError)
ELSE
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)
END IF
#endif
!Write the dataset collectively.

#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(IntArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,IntArray, Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(RealArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,RealArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,StrArray, Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
#else
IF(PRESENT(IntArray))  buf=C_LOC(IntArray)
IF(PRESENT(RealArray)) buf=C_LOC(RealArray)
IF(PRESENT(StrArray))  buf=C_LOC(StrArray(1))
CALL H5DWRITE_F(DSet_ID,Type_ID,buf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
#endif /* HDF5_F90 */

IF(PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close the property list.
CALL H5PCLOSE_F(dsetparams, iError)
CALL H5PCLOSE_F(PList_ID, iError)
! Close dataspaces.
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
! Close the dataset.
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteArrayToHDF5


SUBROUTINE WriteAttributeToHDF5(Loc_ID_in,AttribName,nVal,DataSetname,&
                                RealScalar,IntScalar,StrScalar,LogicalScalar, &
                                RealArray,IntArray,StrArray)
!===================================================================================================================================
! Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,goupID. This must be opened
! outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T),INTENT(IN)              :: Loc_ID_in
CHARACTER(LEN=*), INTENT(IN)           :: AttribName
INTEGER,INTENT(IN)                     :: nVal
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL  :: DatasetName
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealScalar
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntScalar
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrScalar(1)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(nVal)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray(nVal)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(nVal)
LOGICAL           ,INTENT(IN),OPTIONAL        :: LogicalScalar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: DataSpace,Attr_ID,Loc_ID,Type_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER,TARGET                 :: logtoint
TYPE(C_PTR)                    :: buf
!===================================================================================================================================
LOGWRITE(*,*)' WRITE ATTRIBUTE "',TRIM(AttribName),'" TO HDF5 FILE...'
IF(PRESENT(DataSetName))THEN
 ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
Rank=1
Dimsf(:)=0
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
! Write the attribute data.
IF(PRESENT(RealScalar)) Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealArray))  Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntScalar))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntArray))   Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(LogicalScalar))THEN
  LogToInt=MERGE(1,0,LogicalScalar)
  Type_ID=H5T_NATIVE_INTEGER
  END IF
IF(PRESENT(StrScalar).OR.PRESENT(StrArray))THEN
  ! Create character string datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  IF(PRESENT(StrScalar))THEN
    AttrLen=LEN(StrScalar(1))
  ELSE
  AttrLen=255
END IF
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
ENDIF

CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), Type_ID, DataSpace, Attr_ID, iError)
! Write the attribute data.
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(RealArray))     CALL H5AWRITE_F(Attr_ID, Type_ID, RealArray,  Dimsf, iError)
IF(PRESENT(RealScalar))    CALL H5AWRITE_F(Attr_ID, Type_ID, RealScalar, Dimsf, iError)
IF(PRESENT(IntArray))      CALL H5AWRITE_F(Attr_ID, Type_ID, IntArray,   Dimsf, iError)
IF(PRESENT(IntScalar))     CALL H5AWRITE_F(Attr_ID, Type_ID, IntScalar,  Dimsf, iError)
IF(PRESENT(LogicalScalar)) CALL H5AWRITE_F(Attr_ID, Type_ID, LogToInt,   Dimsf, iError)
IF(PRESENT(StrScalar))     CALL H5AWRITE_F(Attr_ID, Type_ID, StrScalar,  Dimsf, iError)
IF(PRESENT(StrArray))      CALL H5AWRITE_F(Attr_ID, Type_ID, StrArray,   Dimsf, iError)
#else /* HDF5_F90 */
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntArray))      buf=C_LOC(IntArray)
IF(PRESENT(IntScalar))     buf=C_LOC(IntScalar)
IF(PRESENT(LogicalScalar)) buf=C_LOC(LogToInt)
IF(PRESENT(StrScalar))     buf=C_LOC(StrScalar(1))
IF(PRESENT(StrArray))      buf=C_LOC(StrArray(1))
CALL H5AWRITE_F(Attr_ID, Type_ID, buf, iError)
#endif /* HDF5_F90 */

! Close datatype
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close dataspace
CALL H5SCLOSE_F(DataSpace, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteAttributeToHDF5


END MODULE MODH_HDF5_output
