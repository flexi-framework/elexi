!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
!> Module for the Particle HDF5 Output
!==================================================================================================================================
MODULE MOD_Particle_HDF5_output
! MODULES
USE MOD_IO_HDF5
USE MOD_HDF5_WriteArray
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE WriteParticleToHDF5
  MODULE PROCEDURE WriteParticleToHDF5
END INTERFACE

INTERFACE WriteAttributeToHDF5
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

INTERFACE WriteHDF5Header
  MODULE PROCEDURE WriteHDF5Header
END INTERFACE

!INTERFACE WriteArrayToHDF5
!  MODULE PROCEDURE WriteArrayToHDF5
!END INTERFACE

!INTERFACE DistributedWriteArray
!  MODULE PROCEDURE DistributedWriteArray
!END INTERFACE

PUBLIC :: WriteParticleToHDF5
PUBLIC :: WriteAttributeToHDF5
PUBLIC :: WriteHDF5Header
PUBLIC :: WriteArrayToHDF5
#if USE_MPI
PUBLIC :: DistributedWriteArray
#endif
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
SUBROUTINE WriteParticleToHDF5(FileName)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Mesh_Vars,             ONLY: nGlobalElems, offsetElem
USE MOD_Part_Tools,            ONLY: UpdateNextFreePosition
USE MOD_Particle_Erosion_Vars, ONLY: PartTrackReflection
USE MOD_Particle_Vars,         ONLY: PDM, PEM, PartState, PartSpecies,PartReflCount
USE MOD_Particle_Vars,         ONLY: useLinkedList
#if USE_MPI
USE MOD_Particle_MPI_Vars,     ONLY: PartMPI
#endif /*MPI*/
#if CODE_ANALYZE
USE MOD_Particle_Tracking_Vars,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVar
#if USE_MPI
INTEGER                        :: sendbuf(2),recvbuf(2)
INTEGER                        :: nParticles(0:nProcessors-1)
#endif
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
INTEGER                        :: locnPart,offsetnPart
INTEGER                        :: iPart,nPart_glob, iElem_glob, iElem_loc
INTEGER,ALLOCATABLE            :: PartInt(:,:)
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER,PARAMETER              :: PartIntSize=2      !number of entries in each line of PartInt
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER                        :: locnPart_max
#if USE_SM
INTEGER                        :: i,j,k,m,iPartition
REAL                           :: L_xi(3,0:PP_N), L_eta_zeta
#endif
!=============================================
! Write properties -----------------------------------------------------------------------------------------------------------------

  ! Check if we are counting reflections
  IF (PartTrackReflection) THEN
    PartDataSize=8
  ELSE
    PartDataSize=7
  END IF

  locnPart =   0
  DO pcount = 1,PDM%ParticleVecLength
    IF(PDM%ParticleInside(pcount)) THEN
      locnPart = locnPart + 1
    END IF
  END DO

#if USE_MPI
  sendbuf(1)=locnPart
  recvbuf=0
  CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
  offsetnPart=recvbuf(1)
  sendbuf(1)=recvbuf(1)+locnPart
  CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
  !global numbers
  nPart_glob=sendbuf(1)
  CALL MPI_GATHER(locnPart,1,MPI_INTEGER,nParticles,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
  LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
  CALL MPI_REDUCE(locnPart, locnPart_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)
#else
  offsetnPart=0
  nPart_glob=locnPart
  locnPart_max=locnPart
#endif

  ALLOCATE(PartInt( PartIntSize ,offsetElem+1 :offsetElem+PP_nElems))
  ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))

!!! Kleiner Hack von JN (Teil 1/2):

  ALLOCATE(PEM%pStart(1:PP_nElems)           , &
           PEM%pNumber(1:PP_nElems)          , &
           PEM%pNext(1:PDM%maxParticleNumber), &
           PEM%pEnd(1:PP_nElems) )
  useLinkedList=.TRUE.
  CALL UpdateNextFreePosition()

!!! Ende kleiner Hack von JN (Teil 1/2)
  iPart=offsetnPart
  DO iElem_loc=1,PP_nElems
    iElem_glob = iElem_loc + offsetElem
    PartInt(1,iElem_glob)=iPart
    IF (ALLOCATED(PEM%pNumber)) THEN
      PartInt(2,iElem_glob) = PartInt(1,iElem_glob) + PEM%pNumber(iElem_loc)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(1,iElem_glob)+1,PartInt(2,iElem_glob)
        PartData(1,iPart)=PartState(1,pcount)
        PartData(2,iPart)=PartState(2,pcount)
        PartData(3,iPart)=PartState(3,pcount)
        PartData(4,iPart)=PartState(4,pcount)
        PartData(5,iPart)=PartState(5,pcount)
        PartData(6,iPart)=PartState(6,pcount)
        PartData(7,iPart)=REAL(PartSpecies(pcount))
        IF (PartTrackReflection) THEN
            PartData(8,iPart)=REAL(PartReflCount(pcount))
        END IF

#if CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(pcount.EQ.PARTOUT)THEN
            PartData(7,iPart)=-PartData(7,iPart)
          END IF
        END IF
#endif /*CODE_ANALYZE*/

        pcount = PEM%pNext(pcount)
      END DO
      iPart = PartInt(2,iElem_glob)
    ELSE
      CALL abort(&
      __STAMP__&
      , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
    END IF
    PartInt(2,iElem_glob)=iPart
  END DO

  nVar=2
  ALLOCATE(StrVarNames(nVar))
  StrVarNames(1)='FirstPartID'
  StrVarNames(2)='LastPartID'

  IF(MPIRoot)THEN
#if USE_MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif
    CALL WriteAttributeToHDF5(File_ID,'VarNamesPartInt',nVar,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

  reSwitch=.FALSE.
  IF(gatheredWrite)THEN
    ! gatheredwrite not working with distributed particles
    ! particles require own routine for which the communicator has to be build each time
    reSwitch=.TRUE.
    gatheredWrite=.FALSE.
  END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems)    ,&
      nVar            => INT(nVar)            ,&
      PP_nElems       => INT(PP_nElems)       ,&
      offsetElem      => INT(offsetElem)      ,&
      PartDataSize    => INT(PartDataSize)    )

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='PartInt',rank=2        ,&
                          nValGlobal=(/nVar    ,nGlobalElems/),&
                          nVal=      (/nVar    ,PP_nElems   /),&
                          offset=    (/0       ,offsetElem  /),&
                          collective=.TRUE.    ,IntArray=PartInt)
  DEALLOCATE(StrVarNames)

  ALLOCATE(StrVarNames(PartDataSize))
  StrVarNames(1)='ParticlePositionX'
  StrVarNames(2)='ParticlePositionY'
  StrVarNames(3)='ParticlePositionZ'
  StrVarNames(4)='VelocityX'
  StrVarNames(5)='VelocityY'
  StrVarNames(6)='VelocityZ'
  StrVarNames(7)='Species'
  IF (PartTrackReflection) THEN
      StrVarNames(8)='ReflectionCount'
  END IF

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',PartDataSize,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

#if USE_MPI
 CALL DistributedWriteArray(FileName,&
                            DataSetName='PartData'   ,rank=2       ,&
                            nValGlobal=(/PartDataSize,nPart_glob/) ,&
                            nVal=      (/PartDataSize,locnPart  /) ,&
                            offset=    (/0           ,offsetnPart/),&
                            collective=.FALSE.,offSetDim=2         ,&
                            communicator=PartMPI%COMM,RealArray=PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName='PartData'   ,rank=2          ,&
                        nValGlobal=(/PartDataSize,nPart_glob/)    ,&
                        nVal=      (/PartDataSize,locnPart/)      ,&
                        offset=    (/0           ,offsetnPart/)   ,&
                        collective=.TRUE., RealArray=PartData)
  CALL CloseDataFile()
#endif /*MPI*/


  END ASSOCIATE
  ! reswitch
  IF(reSwitch) gatheredWrite=.TRUE.

  DEALLOCATE(StrVarNames)
  DEALLOCATE(PartInt)
  DEALLOCATE(PartData)

!!! Kleiner Hack von JN (Teil 2/2):
  useLinkedList=.FALSE.
  DEALLOCATE(PEM%pStart , &
             PEM%pNumber, &
             PEM%pNext  , &
             PEM%pEnd   )
!!! Ende kleiner Hack von JN (Teil 2/2)


END SUBROUTINE WriteParticleToHDF5


SUBROUTINE WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,offset,&
                            collective,resizeDim,chunkSize,&
                            RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: DataSetName
INTEGER,INTENT(IN)            :: rank             ! number of dimensions of the array
INTEGER,INTENT(IN)            :: nValGlobal(rank) ! max size of array in offset dimension
INTEGER,INTENT(IN)            :: nVal(rank)       ! size of complete (local) array to write
INTEGER,INTENT(IN)            :: offset(rank)     ! offset =0, start at beginning of the array
LOGICAL,INTENT(IN)            :: collective       ! use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL   :: resizeDim(rank)  ! specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL   :: chunkSize(rank)  ! specify chunksize
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(rank)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray(rank)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(rank)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,Type_ID,dsetparams
INTEGER(HSIZE_T)               :: Dimsf(Rank),OffsetHDF(Rank),nValMax(Rank)
INTEGER(SIZE_T)                :: SizeSet=255
LOGICAL                        :: chunky
#ifndef HDF5_F90 /* HDF5 compiled with fortran2003 flag */
TYPE(C_PTR)                    :: buf
#endif
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'

! specify chunk size if desired
nValMax=nValGlobal
chunky=.FALSE.
CALL H5PCREATE_F(H5P_DATASET_CREATE_F,dsetparams,iError)
IF(PRESENT(chunkSize))THEN
  chunky=.TRUE.
  Dimsf=chunkSize
  CALL H5PSET_CHUNK_F(dsetparams,rank,dimsf,iError)
END IF
! make array extendable in case you want to append something
IF(PRESENT(resizeDim))THEN
  IF(.NOT.PRESENT(chunkSize))&
    CALL abort(&
    __STAMP__&
    ,'Chunk size has to be specified when using resizable arrays.')
  nValMax = MERGE(H5S_UNLIMITED_F,nValMax,resizeDim)
END IF

! Create the dataset with default properties.
IF(PRESENT(RealArray))     Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  SizeSet=255
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
#if USE_MPI
IF(collective)THEN
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F,  iError)
ELSE
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)
END IF
#endif

!Write the dataset collectively.
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(IntegerArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,IntegerArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(RealArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,RealArray   ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,StrArray    ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
#else
IF(PRESENT(IntegerArray)) buf=C_LOC(IntegerArray)
IF(PRESENT(RealArray))    buf=C_LOC(RealArray)
IF(PRESENT(StrArray))     buf=C_LOC(StrArray(1))
!IF(ANY(Dimsf.EQ.0)) buf =NULL()
IF(ANY(Dimsf.EQ.0)) THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,C_NULL_PTR,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
ELSE
  CALL H5DWRITE_F(DSet_ID,Type_ID,buf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
#endif /* HDF5_F90 */

IF(PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close the property list, dataspaces and dataset.
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
                                RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,groupID. This must be opened
! outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)           :: Loc_ID_in
CHARACTER(LEN=*)  ,INTENT(IN)           :: AttribName
INTEGER           ,INTENT(IN)           :: nVal
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL  :: DatasetName
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealScalar
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerScalar
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrScalar(1)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(nVal)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray(nVal)
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
#ifndef HDF5_F90 /* HDF5 compiled with fortran2003 flag */
TYPE(C_PTR)                    :: buf
#endif
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
Dimsf(:)=0 !???
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
IF(PRESENT(RealScalar))    Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealArray))     Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerScalar)) Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntegerArray))  Type_ID=H5T_NATIVE_INTEGER
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
IF(PRESENT(RealArray))     CALL H5AWRITE_F(Attr_ID, Type_ID, RealArray,     Dimsf, iError)
IF(PRESENT(RealScalar))    CALL H5AWRITE_F(Attr_ID, Type_ID, RealScalar,    Dimsf, iError)
IF(PRESENT(IntegerArray))  CALL H5AWRITE_F(Attr_ID, Type_ID, IntegerArray,  Dimsf, iError)
IF(PRESENT(IntegerScalar)) CALL H5AWRITE_F(Attr_ID, Type_ID, IntegerScalar, Dimsf, iError)
IF(PRESENT(LogicalScalar)) CALL H5AWRITE_F(Attr_ID, Type_ID, LogToInt,      Dimsf, iError)
IF(PRESENT(StrScalar))     CALL H5AWRITE_F(Attr_ID, Type_ID, StrScalar,     Dimsf, iError)
IF(PRESENT(StrArray))      CALL H5AWRITE_F(Attr_ID, Type_ID, StrArray,      Dimsf, iError)
#else /* HDF5_F90 */
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntegerArray))  buf=C_LOC(IntegerArray)
IF(PRESENT(IntegerScalar)) buf=C_LOC(IntegerScalar)
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

SUBROUTINE WriteHDF5Header(FileType_in,File_ID)
!===================================================================================================================================
! Subroutine to write a distinct file header to each HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,                      ONLY:ProgramName,FileVersion,ProjectName
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
CHARACTER(LEN=255)                       :: tmp255
!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files
! Attributes are program name, file type identifier, project name and version number

!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files

! First write program name
tmp255=TRIM(ProgramName)
CALL WriteAttributeToHDF5(File_ID,'Program'     ,1,StrScalar=(/tmp255/))
tmp255=TRIM(FileType_in)
CALL WriteAttributeToHDF5(File_ID,'File_Type'   ,1,StrScalar=(/tmp255/))
tmp255=TRIM(ProjectName)
CALL WriteAttributeToHDF5(File_ID,'Project_Name',1,StrScalar=(/tmp255/))
CALL WriteAttributeToHDF5(File_ID,'File_Version',1,RealScalar=FileVersion)
END SUBROUTINE WriteHDF5Header

#if USE_MPI
SUBROUTINE DistributedWriteArray(FileName,DataSetName,rank,nValGlobal,nVal,offset,collective,&
                                 offSetDim,communicator,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Write distributed data to proc, e.g. particles which are not hosted by each proc
! a new output-communicator is build and afterwards killed
! offset is in the last dimension
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5,                ONLY: OpenDataFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName,DataSetName
LOGICAL,INTENT(IN)             :: collective
INTEGER,INTENT(IN)             :: offSetDim,communicator
INTEGER,INTENT(IN)             :: rank,nVal(rank),nValGlobal(rank),offset(rank)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray( PRODUCT(nVal))
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                        :: Color, OutPutCOMM,nOutPutProcs,MyOutputRank
LOGICAL                        :: DataOnProc, DoNotSplit
!===================================================================================================================================

DataOnProc=.FALSE.
IF(nVal(offSetDim).GT.0) DataOnProc=.TRUE.
CALL MPI_ALLREDUCE(DataOnProc,DoNotSplit, 1, MPI_LOGICAL, MPI_LAND, COMMUNICATOR, IERROR)


IF(.NOT.DoNotSplit)THEN
  color=MPI_UNDEFINED
  IF(DataOnProc) color=87
  MyOutputRank=0

  CALL MPI_COMM_SPLIT(COMMUNICATOR, color, MyOutputRank, OutputCOMM,iError)
  IF(DataOnProc) THEN
    CALL MPI_COMM_SIZE(OutputCOMM, nOutPutProcs,iError)
    IF(nOutPutProcs.EQ.1)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective=.FALSE.,RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective=.FALSE.,IntegerArray =IntegerArray)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective=.FALSE.,StrArray =StrArray)
      CALL CloseDataFile()
    ELSE
      CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective,RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective,IntegerArray =IntegerArray)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective,StrArray =StrArray)
      CALL CloseDataFile()
    END IF
    CALL MPI_BARRIER(OutputCOMM,IERROR)
    CALL MPI_COMM_FREE(OutputCOMM,iERROR)
  END IF
  ! MPI Barrier is requried, that the other procs don't open the datafile while this procs are still writring
  CALL MPI_BARRIER(COMMUNICATOR,IERROR)
  OutputCOMM=MPI_UNDEFINED
ELSE
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#else
!  CALL OpenDataFile(FileName,create=.FALSE.,readOnly=.FALSE.)
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif
  IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,RealArray=RealArray)
  IF(PRESENT(IntegerArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,IntegerArray =IntegerArray)
  IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,StrArray =StrArray)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

END SUBROUTINE DistributedWriteArray
#endif /*MPI*/

END MODULE MOD_Particle_HDF5_output
