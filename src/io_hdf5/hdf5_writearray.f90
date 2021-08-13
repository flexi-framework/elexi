!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
!> Module providing IO routines for parallel output in HDF5 format: solution, time averaged files, baseflow, record points,...
!==================================================================================================================================
MODULE MOD_HDF5_WriteArray
! MODULES
USE MOD_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
!INTERFACE WriteArray
!  MODULE PROCEDURE WriteArray
!END INTERFACE

!INTERFACE GatheredWriteArray
!  MODULE PROCEDURE GatheredWriteArray
!END INTERFACE

PUBLIC :: WriteArray
PUBLIC :: GatheredWriteArray
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine is a wrapper routine for WriteArray and first gathers all output arrays of an MPI sub group,
!> then only the master will write the data. Significantly reduces IO overhead for a large number of processes!
!==================================================================================================================================
SUBROUTINE GatheredWriteArray(FileName,create,DataSetName,rank,nValGlobal,nVal,offset,collective,RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName          !< Name of the file to write to
CHARACTER(LEN=*),INTENT(IN)    :: DataSetName       !< Name of the dataset to write
LOGICAL,INTENT(IN)             :: create            !< Should the file be created or not
LOGICAL,INTENT(IN)             :: collective        !< Collective write or not
INTEGER,INTENT(IN)             :: rank              !< Rank of array
INTEGER,INTENT(IN)             :: nVal(rank)        !< Local number of variables in each rank
INTEGER,INTENT(IN)             :: nValGlobal(rank)  !< Global number of variables in each rank
INTEGER,INTENT(IN)             :: offset(rank)      !< Offset in each rank
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal)) !< Real array to write
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray( PRODUCT(nVal)) !< Integer array to write
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal)) !< String array to write
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
REAL,              ALLOCATABLE :: UReal(:)
CHARACTER(LEN=255),ALLOCATABLE :: UStr(:)
INTEGER,           ALLOCATABLE :: UInt(:)
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nLocalProcs) :: nDOFPerNode,offsetNode
!==================================================================================================================================
! HDF5 with MPI can only write max. (32 bit integer / 8) elements
IF(REAL(PRODUCT(nVal)).GT.((2**28-1)/8.))  CALL Abort(__STAMP__, &
    'Total size of HDF5 array "'//TRIM(DataSetName)//'" is too big! Reduce number of entries per rank or compile without MPI!')

IF(gatheredWrite)THEN
  IF(ANY(offset(1:rank-1).NE.0)) &
    CALL abort(__STAMP__,'Offset only allowed in last dimension for gathered IO.')

  ! Get last dim of each array on IO nodes
  nDOFLocal=PRODUCT(nVal)
  CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerNode,1,MPI_INTEGER,0,MPI_COMM_NODE,iError)

  ! Allocate big array and compute offsets of small arrs inside big
  offsetNode=0
  IF(MPILocalRoot)THEN
    nValGather=nVal
    nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
    DO i=2,nLocalProcs
      offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
    END DO
    IF(PRESENT(RealArray)) ALLOCATE(UReal(PRODUCT(nValGather)))
    IF(PRESENT(IntArray))  ALLOCATE(UInt( PRODUCT(nValGather)))
    IF(PRESENT(StrArray))  ALLOCATE(UStr( PRODUCT(nValGather)))
  ELSE
    IF(PRESENT(RealArray)) ALLOCATE(UReal(1))
    IF(PRESENT(IntArray))  ALLOCATE(UInt( 1))
    IF(PRESENT(StrArray))  ALLOCATE(UStr( 1))
  ENDIF

  ! Gather small arrays on IO nodes
  IF(PRESENT(RealArray)) CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
                                          UReal,nDOFPerNode,offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_NODE,iError)
  IF(PRESENT(IntArray))  CALL MPI_GATHERV(IntArray, nDOFLocal,MPI_INTEGER,&
                                          UInt, nDOFPerNode,offsetNode,MPI_INTEGER,0,MPI_COMM_NODE,iError)
  !IF(PRESENT(StrArray))  CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
  !                                        UReal,nDOFPerNode, offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_NODE,iError)

  IF(MPILocalRoot)THEN
    ! Reopen file and write DG solution (only IO nodes)
    CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS)
    IF(PRESENT(RealArray)) CALL WriteArray(DataSetName,rank,nValGlobal,nValGather,&
                                                 offset,collective=collective,RealArray=UReal)
    IF(PRESENT(IntArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nValGather,&
                                                 offset,collective=collective,IntArray =UInt)
    !IF(PRESENT(StrArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nValGather,&
    !                                             offset,collective=collective,StrArr =UStr)
    CALL CloseDataFile()
  END IF

  SDEALLOCATE(UReal)
  SDEALLOCATE(UInt)
  SDEALLOCATE(UStr)
ELSE
#endif
  CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.)
  IF(PRESENT(RealArray)) CALL WriteArray(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,RealArray=RealArray)
  IF(PRESENT(IntArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,IntArray =IntArray)
  IF(PRESENT(StrArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,StrArray =StrArray)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

END SUBROUTINE GatheredWriteArray

!==================================================================================================================================
!> Low-level subroutine to actually write data to HDF5 format
!==================================================================================================================================
SUBROUTINE WriteArray(DataSetName,rank,nValGlobal,nVal,offset,&
                            collective,resizeDim,chunkSize,&
                            RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: DataSetName      !< name of the dataset to write the data into
INTEGER,INTENT(IN)            :: rank             !< number of dimensions of the array
INTEGER,INTENT(IN)            :: nValGlobal(rank) !< max size of array in offset dimension
INTEGER,INTENT(IN)            :: nVal(rank)       !< size of complete (local) array to write
INTEGER,INTENT(IN)            :: offset(rank)     !< offset =0, start at beginning of the array
LOGICAL,INTENT(IN)            :: collective       !< use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL   :: resizeDim(rank)  !< specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL   :: chunkSize(rank)  !< specify chunksize
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal)) !< number of array entries
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray(PRODUCT(nVal))  !< number of array entries
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(PRODUCT(nVal))  !< number of array entries (length 255)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,Type_ID,dsetparams
INTEGER(HSIZE_T)               :: Dimsf(Rank),OffsetHDF(Rank),nValMax(Rank)
INTEGER(SIZE_T)                :: SizeSet=255
LOGICAL                        :: chunky,exists
TYPE(C_PTR)                    :: buf
!==================================================================================================================================
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
! Check data set. Data sets can be checked by determining the existence of the corresponding link
CALL H5LEXISTS_F(File_ID, TRIM(DataSetName), exists, iError)
IF (.NOT. exists) THEN
  ! Create the data space for the  dataset.
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError, nValMax)
  CALL H5DCREATE_F(File_ID, TRIM(DataSetName), Type_ID, FileSpace, DSet_ID,iError,dsetparams)
  CALL H5SCLOSE_F(FileSpace, iError)
ELSE
  CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
END IF
IF(chunky)THEN
  CALL H5DSET_EXTENT_F(DSet_ID,Dimsf,iError) ! if resizable then dataset may need to be extended
END IF

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf     = nVal  ! Now we need the local array size
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
IF(PRESENT(IntArray))  buf=C_LOC(IntArray)
IF(PRESENT(RealArray)) buf=C_LOC(RealArray)
IF(PRESENT(StrArray))  buf=C_LOC(StrArray(1))
CALL H5DWRITE_F(DSet_ID,Type_ID,buf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)

IF(PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close the property list, dataspaces and dataset.
CALL H5PCLOSE_F(dsetparams, iError)
CALL H5PCLOSE_F(PList_ID, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(*,*)'...DONE!'

END SUBROUTINE WriteArray

END MODULE MOD_HDF5_WriteArray
