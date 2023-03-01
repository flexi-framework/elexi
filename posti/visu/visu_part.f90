!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving ListInEs with discontinuous Galerkin methods.
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
!=================================================================================================================================
! Routines needed to visualize particles
!=================================================================================================================================
MODULE MOD_Posti_Part_Tools
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

#if USE_PARTICLES
INTERFACE InitParticleOutput
  MODULE PROCEDURE InitParticleOutput
END INTERFACE

INTERFACE InitPartStatistics
  MODULE PROCEDURE InitPartStatistics
END INTERFACE

INTERFACE InitPartState
  MODULE PROCEDURE InitPartState
END INTERFACE

INTERFACE ReadPartStatistics
  MODULE PROCEDURE ReadPartStatistics
END INTERFACE

INTERFACE ReadPartStateFile
  MODULE PROCEDURE ReadPartStateFile
END INTERFACE

INTERFACE PartitionPartMPI
  MODULE PROCEDURE PartitionPartMPI
END INTERFACE

INTERFACE FinalizeReadPartStateFile
  MODULE PROCEDURE FinalizeReadPartStateFile
END INTERFACE

PUBLIC :: InitParticleOutput
PUBLIC :: InitPartState
PUBLIC :: InitPartStatistics
PUBLIC :: ReadPartStateFile
PUBLIC :: ReadPartStatistics
PUBLIC :: PartitionPartMPI
PUBLIC :: FinalizeReadPartStateFile
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleOutput(ListIn)
!===================================================================================================================================
! Initialize the visualization and map the variable names to classify these in conservative and derived quantities.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Visu_Vars                     ,ONLY: tVisuParticle
USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tVisuParticle),INTENT(INOUT)     :: ListIn
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: i, iVar
CHARACTER(LEN=255)                    :: tmp, tmp2
!===================================================================================================================================

IF(ListIn%nPartVar_Visu.GT.0)THEN
  SDEALLOCATE(ListIn%VarNamePartCombine)
  SDEALLOCATE(ListIn%VarNamePartCombineLen)
  ALLOCATE (ListIn%VarNamePartCombineLen(ListIn%nPartVar_Visu) &
           ,ListIn%VarNamePartCombine(   ListIn%nPartVar_Visu))
  ListIn%VarNamePartCombineLen = 0
  ListIn%VarNamePartCombine    = 0

  DO iVar=2,ListIn%nPartVar_Visu
    i = LEN(TRIM(ListIn%VarNamePartVisu(iVar)))
    tmp = ListIn%VarNamePartVisu(iVar)
    tmp2 = ListIn%VarNamePartVisu(iVar-1)
    IF (TRIM(tmp(:i-1)) .EQ. TRIM(tmp2(:i-1))) THEN
      IF (ListIn%VarNamePartCombine(iVar-1) .EQ. 0) ListIn%VarNamePartCombine(iVar-1) = 1
      ListIn%VarNamePartCombine(iVar) = ListIn%VarNamePartCombine(iVar-1) + 1
    END IF
  END DO

  ListIn%VarNamePartCombineLen(ListIn%nPartVar_visu) = ListIn%VarNamePartCombine(ListIn%nPartVar_visu)

  DO iVar=ListIn%nPartVar_visu-1,1,-1
    IF (ListIn%VarNamePartCombine(iVar).GT.0) THEN
      ListIn%VarNamePartCombineLen(iVar) = MAX(ListIn%VarNamePartCombine(iVar), ListIn%VarNamePartCombineLen(iVar+1))
    END IF
  END DO
END IF

END SUBROUTINE InitParticleOutput


SUBROUTINE InitPartStatistics(AttribName,statefile,VarNames,AttribExists)
!===================================================================================================================================
! Read in main attributes from given HDF5 state file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5                ,ONLY: File_ID
USE MOD_HDF5_Input             ,ONLY: HSize,ReadArray,GetVarNames,GetDataSize,GetAttributeSize
USE MOD_HDF5_Input             ,ONLY: DatasetExists,ReadAttribute
USE MOD_HDF5_Input             ,ONLY: OpenDataFile,CloseDataFile
USE MOD_Visu_Vars              ,ONLY: nSpecies,PartSpeciesIndex
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=*),INTENT(IN)                :: AttribName
CHARACTER(LEN=255),INTENT(IN)              :: statefile
CHARACTER(LEN=255),ALLOCATABLE,INTENT(OUT) :: VarNames(:)
LOGICAL,INTENT(OUT)                        :: AttribExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Parameters
INTEGER,PARAMETER                  :: ELEM_FirstPartInd = 1
INTEGER,PARAMETER                  :: ELEM_LastPartInd  = 2
INTEGER,PARAMETER                  :: PartIntSize       = 2           ! number of entries in each line of PartInt
! Counters
INTEGER                            :: dims,nVal,iVar
! HDF5
INTEGER                            :: PartDim                         ! dummy for rank of partData
INTEGER                            :: PartDataSize
CHARACTER(LEN=255),ALLOCATABLE     :: VarNamesPart(:)
! Formatting
! CHARACTER(LEN=10)                  :: formatstr
CHARACTER(LEN=40)                  :: tmpStr
! Particle variables
INTEGER                            :: nGlobalParts
INTEGER                            :: nGlobalElems
INTEGER,ALLOCATABLE                :: PartInt( :,:)
REAL,ALLOCATABLE                   :: PartData(:,:)
!===================================================================================================================================

SDEALLOCATE(VarNames)
CALL DatasetExists(File_ID,AttribName,AttribExists,attrib=.TRUE.)

IF (.NOT.AttribExists) RETURN

! get size of array
CALL GetAttributeSize(File_ID,AttribName,dims,HSize)
nVal = INT(HSize(1))
DEALLOCATE(HSize)
ALLOCATE(VarNamesPart(nVal))

! read variable names
CALL ReadAttribute(File_ID,TRIM(AttribName),nVal,StrArray=VarNamesPart)

! Find the index of the PartSpecies
PartSpeciesIndex = 0
DO iVar = 1,nVal
  IF (TRIM(VarNamesPart(iVar)).EQ.'Species') PartSpeciesIndex = iVar
END DO ! iVar = 1,nVal

IF (iVar.EQ.0) RETURN

! Read number of local particles and their offset from HDF5
nSpecies = 0

CALL CloseDataFile()
IF (MPIRoot) THEN
  CALL OpenDataFile(statefile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  CALL GetDataSize(File_ID,'PartInt',PartDim,HSize)
  CHECKSAFEINT(HSize(2),4)
  nGlobalElems = INT(HSize(2))
  ALLOCATE(PartInt(PartIntSize,nGlobalElems))
  CALL ReadArray('PartInt',2,(/PartIntSize,nGlobalElems/),0,2,IntArray=PartInt)

  nGlobalParts = PartInt(ELEM_LastPartInd,nGlobalElems)
  DEALLOCATE(PartInt)

  ! Get the number of species
  CALL GetDataSize(File_ID,'PartData',PartDim,HSize)
  CHECKSAFEINT(HSize(2),4)
  PartDataSize = INT(HSize(1))

  ! Get PartData
  ALLOCATE(PartData(PartDataSize,nGlobalParts))
  CALL ReadArray('PartData',2,(/PartDataSize,nGlobalParts/),0,2,RealArray=PartData)

  ! Find the number of species
  nSpecies = MAXVAL(INT(PartData(PartSpeciesIndex,:)))
  CALL CloseDataFile()
END IF
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

#if USE_MPI
CALL MPI_BCAST(nSpecies,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

! Create dummy varnames for the PartStatistics
ALLOCATE(VarNames(nSpecies))
DO iVar = 1,nSpecies
  WRITE(tmpStr,'(A,I0.3)') 'PartDensity',iVar
  VarNames(iVar) = TRIM(tmpStr)
END DO ! iVar = 1,nSpecies

END SUBROUTINE InitPartStatistics


SUBROUTINE ReadPartStatistics(VariableName,ElemData)
!===================================================================================================================================
! Read in particle state from given HDF5 state file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars           ,ONLY: ElemVol,wGPVol
USE MOD_HDF5_Input             ,ONLY: HSize,ReadArray,GetDataSize
USE MOD_HDF5_Input             ,ONLY: ReadArray, OpenDataFile, CloseDataFile, DatasetExists
USE MOD_Interpolation_Vars     ,ONLY: wGP
USE MOD_IO_HDF5                ,ONLY: File_ID
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nElems,sJ
USE MOD_Visu_Vars              ,ONLY: PartSpeciesIndex
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)      :: VariableName
REAL,INTENT(INOUT)                 :: ElemData(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Parameters
INTEGER,PARAMETER                  :: ELEM_FirstPartInd = 1
INTEGER,PARAMETER                  :: ELEM_LastPartInd  = 2
INTEGER,PARAMETER                  :: PartIntSize       = 2           ! number of entries in each line of PartInt
! Counters
INTEGER                            :: i,j,k
INTEGER                            :: locnPart,offsetnPart
INTEGER                            :: FirstElemInd,LastElemInd,iElem
! HDF5
INTEGER                            :: PartDim                         ! dummy for rank of partData
INTEGER                            :: PartDataSize
CHARACTER(LEN=255)                 :: dummyStr
! Mesh
! REAL                               :: ElemVol
! Particle variables
INTEGER                            :: SpeciesInd,iPart
INTEGER,ALLOCATABLE                :: PartInt( :,:)
REAL,ALLOCATABLE                   :: PartData(:,:)
!===================================================================================================================================

FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems

! Get PartInt
ALLOCATE(PartInt(PartIntSize,FirstElemInd:LastElemInd))
CALL ReadArray('PartInt',2,(/PartIntSize,nElems/),offsetElem,2,IntArray=PartInt)

locnPart    = PartInt(ELEM_LastPartInd ,LastElemInd )-PartInt(ELEM_FirstPartInd,FirstElemInd)
offsetnPart = PartInt(ELEM_FirstPartInd,FirstElemInd)

! Get the number of particle variables
CALL GetDataSize(File_ID,'PartData',PartDim,HSize)
CHECKSAFEINT(HSize(2),4)
PartDataSize = INT(HSize(1))

! Get PartData
ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
CALL ReadArray('PartData',2,(/PartDataSize,locnPart/),offsetnPart,2,RealArray=PartData)

! Count the number of particles per element
READ(VariableName,'(A11,I3)') dummyStr,SpeciesInd
ElemData = 0
DO iElem = 1,nElems
  locnPart    = PartInt(ELEM_LastPartInd ,offsetElem+iElem)-PartInt(ELEM_FirstPartInd,offsetElem+iElem)
  offsetnPart = PartInt(ELEM_FirstPartInd,offsetElem+iElem)

  DO iPart = offsetnPart+1,offsetnPart+locnPart
    IF (PartData(PartSpeciesIndex,iPart).EQ.SpeciesInd) &
      ElemData(iElem) = ElemData(iElem) + 1
  END DO ! iPart
END DO ! iElem

! Divide by element volume
IF (.NOT. ALLOCATED(ElemVol)) THEN
  ALLOCATE(ElemVol(nElems))
  ALLOCATE(wGPVol(0:PP_N,0:PP_N,0:PP_NZ))
  ElemVol = 0.

  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    wGPVol(i,j,k) = wGP(i)*wGP(j)*wGP(k)
  END DO; END DO; END DO

  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      ElemVol(iElem) = ElemVol(iElem) + wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
END IF

ElemData = ElemData/ElemVol

END SUBROUTINE ReadPartStatistics


SUBROUTINE InitPartState(datasetNames,ListIn)
!===================================================================================================================================
! Read in main attributes from given HDF5 state file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Visu_Vars              ,ONLY: tVisuParticle
USE MOD_IO_HDF5                ,ONLY: File_ID
USE MOD_HDF5_Input             ,ONLY: HSize,ReadArray,GetVarNames,GetDataSize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),INTENT(IN)           :: datasetNames
TYPE(tVisuParticle),INTENT(INOUT)       :: ListIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: dims, j
LOGICAL                                 :: VarNamesExist
CHARACTER(LEN=255),ALLOCATABLE          :: varnames(:)
!===================================================================================================================================
SDEALLOCATE(ListIn%VarNamesPart_HDF5)
SDEALLOCATE(ListIn%VarNamePartDummy)
SDEALLOCATE(ListIn%mapAllVarsToVisuVars)

SELECT CASE(TRIM(datasetNames))
  CASE('PartData')
    CALL GetVarNames("VarNamesParticles"     ,varnames,VarNamesExist)
  CASE('ImpactData')
    CALL GetVarNames("VarNamesImpactTracking",varnames,VarNamesExist)
END SELECT

IF(VarNamesExist)THEN
  CALL GetDataSize(File_ID,TRIM(datasetNames),dims,HSize)
  ListIn%nPartVar_HDF5 = INT(HSize(1))-3
  ListIn%nGlobalParts  = INT(HSize(2))

  ALLOCATE(ListIn%VarNamesPart_HDF5(   1:ListIn%nPartVar_HDF5) &
          ,ListIn%VarNamePartDummy(    1:ListIn%nPartVar_HDF5)  &
          ,ListIn%mapAllVarsToVisuVars(1:ListIn%nPartVar_HDF5))

  ListIn%VarNamesPart_HDF5    = varnames(4:)
  ListIn%VarNamePartDummy     = ''
  ListIn%mapAllVarsToVisuVars = 0

  DO j=1,ListIn%nPartVar_HDF5
    ListIn%VarNamesPart_HDF5(j)=TRIM(datasetNames)//":"//TRIM(ListIn%VarNamesPart_HDF5(j))
  END DO
END IF

SDEALLOCATE(varnames)

END SUBROUTINE InitPartState


SUBROUTINE ReadPartStateFile(InputFile, DataArrayIn, ListIn)
!===================================================================================================================================
! Read in particle state from given HDF5 state file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Visu_Vars              ,ONLY: tVisuParticle, VisuPart
USE MOD_IO_HDF5                ,ONLY: File_ID
USE MOD_HDF5_Input             ,ONLY: ReadArray, OpenDataFile, CloseDataFile, DatasetExists
USE MOD_StringTools            ,ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),INTENT(IN)           :: InputFile
CHARACTER(LEN=255),INTENT(IN),OPTIONAL  :: DataArrayIn
TYPE(tVisuParticle),INTENT(INOUT)       :: ListIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)                      :: DataArray
INTEGER                                 :: iPart,iVar,iVar2,i
INTEGER                                 :: nInvalidPartsLoc,nInvalidPartsGlob
REAL,ALLOCATABLE                        :: PartData(:,:),PartDataCheck(:,:)
CHARACTER(LEN=255)                      :: tmp, tmp2
LOGICAL                                 :: datasetFound
!===================================================================================================================================
SDEALLOCATE(ListIn%PartData_HDF5)

IF (ListIn%nPartVar_HDF5.EQ.0) THEN
  ListIn%nPartVar_Visu = 0
  RETURN
END IF

IF(PRESENT(DataArrayIn))THEN
  DataArray = DataArrayIn
ELSE
  DataArray = 'PartData'
END IF

SWRITE(Unit_stdOut,'(A,A,A)') ' ',TRIM(DataArray),' DECOMPOSITION'

CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,TRIM(DataArray),datasetFound)
IF(.NOT.datasetFound) RETURN

! Assign particles to corresponding procs
CALL PartitionPartMPI(ListIn)

ALLOCATE(PartData     (1:ListIn%nPartVar_HDF5+3,1:ListIn%nLocalParts) &
        ,PartDataCheck(1:ListIn%nPartVar_HDF5+3,1:ListIn%nLocalParts))
CALL ReadArray(DataArray,2,(/ListIn%nPartVar_HDF5+3, ListIn%nLocalParts/),ListIn%offsetPart,2,RealArray=PartData)

ListIn%nPart_Visu = 0
nInvalidPartsLoc  = 0

DO iPart = 1,ListIn%nLocalParts
  ! Remove invalid particles
  IF (ANY(IEEE_IS_NAN(PartData(:,iPart)))) THEN
    nInvalidPartsLoc  = nInvalidPartsLoc  + 1
  ELSE
    ListIn%nPart_Visu = ListIn%nPart_Visu + 1
    PartDataCheck(:,ListIn%nPart_Visu) = PartData(:,iPart)
  END IF
END DO

nInvalidPartsGlob = nInvalidPartsLoc
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,nInvalidPartsGlob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif

! Write data back to PartData
IF (nInvalidPartsGlob.GT.0 .AND. MPIRoot) WRITE(Unit_stdOut,'(A,I0,A)') ' ',nInvalidPartsGlob,' invalid particles detected. Ignoring ...'
IF (nInvalidPartsLoc.GT.0) THEN
  DEALLOCATE(PartData)
  ALLOCATE(  PartData   (1:ListIn%nPartVar_HDF5+3,1:ListIn%nPart_Visu))
  PartData(:,1:ListIn%nPart_Visu) = PartDataCheck(:,1:ListIn%nPart_Visu)
  ListIn%nLocalParts = ListIn%nPart_Visu
END IF
DEALLOCATE(PartDataCheck)

! Visualize all particle variables
IF(VisuPart)THEN
  ListIn%nPartVar_Visu = ListIn%nPartVar_HDF5

  SDEALLOCATE(ListIn%VarNamePartVisu)
  ALLOCATE(ListIn%PartData_HDF5(  1:ListIn%nPartVar_Visu+3,1:ListIn%nLocalParts) &
          ,ListIn%VarNamePartVisu(1:ListIn%nPartVar_Visu))
  ListIn%PartData_HDF5        = PartData
  ListIn%VarNamePartVisu      = ListIn%VarNamesPart_HDF5
  ListIn%mapAllVarsToVisuVars = 1
! Only visualize explicitely requested variables
ELSE
  ALLOCATE(ListIn%PartData_HDF5(1:ListIn%nPartVar_Visu+3,1:ListIn%nLocalParts))
  ListIn%PartData_HDF5(1:3,:) = PartData(1:3,:)
  DO iVar = 1,ListIn%nPartVar_Visu
    DO iVar2 = 1,ListIn%nPartVar_HDF5
      IF(STRICMP(ListIn%VarNamePartVisu(iVar),ListIn%VarNamesPart_HDF5(iVar2)))THEN
        ListIn%PartData_HDF5(iVar+3,:) = PartData(iVar2+3,:)
      END IF
    END DO
  END DO
END IF

DEALLOCATE(PartData)
CALL CloseDataFile()

! Needed for output
IF(ListIn%nPartVar_Visu.GT.0)THEN
  SDEALLOCATE(ListIn%VarNamePartCombine)
  SDEALLOCATE(ListIn%VarNamePartCombineLen)
  ALLOCATE (ListIn%VarNamePartCombine(   ListIn%nPartVar_Visu) &
           ,ListIn%VarNamePartCombineLen(ListIn%nPartVar_Visu))
  ListIn%VarNamePartCombineLen = 0
  ListIn%VarNamePartCombine    = 0

  DO iVar=2,ListIn%nPartVar_Visu
    i = LEN(TRIM(ListIn%VarNamePartVisu(iVar)))
    tmp = ListIn%VarNamePartVisu(iVar)
    tmp2 = ListIn%VarNamePartVisu(iVar-1)
    IF (TRIM(tmp(:i-1)) .EQ. TRIM(tmp2(:i-1))) THEN
      IF (ListIn%VarNamePartCombine(iVar-1) .EQ. 0) ListIn%VarNamePartCombine(iVar-1) = 1
      ListIn%VarNamePartCombine(iVar) = ListIn%VarNamePartCombine(iVar-1) + 1
    END IF
  END DO

  ListIn%VarNamePartCombineLen(ListIn%nPartVar_visu) = ListIn%VarNamePartCombine(ListIn%nPartVar_visu)
  DO iVar=ListIn%nPartVar_visu-1,1,-1
    IF (ListIn%VarNamePartCombine(iVar).GT.0) THEN
      ListIn%VarNamePartCombineLen(iVar) = MAX(ListIn%VarNamePartCombine(iVar), ListIn%VarNamePartCombineLen(iVar+1))
    END IF
  END DO
END IF

END SUBROUTINE ReadPartStateFile


SUBROUTINE FinalizeReadPartStateFile(ListIn)
!==================================================================================================================================
! Deallocates solution vector and reset main attributes
!==================================================================================================================================
! MODULES
USE MOD_Visu_Vars,                       ONLY: tVisuParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tVisuParticle),INTENT(INOUT)       :: ListIn
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(ListIn%VarNamePartVisu)
SDEALLOCATE(ListIn%VarNamePartDummy)
SDEALLOCATE(ListIn%VarNamePartCombine)
SDEALLOCATE(ListIn%VarNamePartCombineLen)
SDEALLOCATE(ListIn%VarNamesPart_HDF5)
SDEALLOCATE(ListIn%PartData_HDF5)
SDEALLOCATE(ListIn%mapAllVarsToVisuVars)
END SUBROUTINE FinalizeReadPartStateFile


SUBROUTINE PartitionPartMPI(ListIn)
!==================================================================================================================================
! Basic MPI initialization for visualization of particles
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Visu_Vars                ,ONLY: tVisuParticle
USE MOD_Particle_MPI_Vars        ,ONLY: PartitionPartIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tVisuParticle),INTENT(INOUT)       :: ListIn
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER,ALLOCATABLE                     :: offsetPartMPI(:)
INTEGER,ALLOCATABLE                     :: nPartsMPI(:)
INTEGER                                 :: iProc,iPart
#endif /*USE_MPI*/
!==================================================================================================================================
#if USE_MPI
ALLOCATE(offsetPartMPI(0:nProcessors) &
        ,nPartsMPI(    0:nProcessors-1))
offsetPartMPI=0

ListIn%nLocalParts = ListIn%nGlobalParts/nProcessors
iPart = ListIn%nGlobalParts-ListIn%nLocalParts*nProcessors

DO iProc = 0,nProcessors-1
  offsetPartMPI(iProc) = ListIn%nLocalParts*iProc+MIN(iProc,iPart)
END DO

offsetPartMPI(nProcessors)   = ListIn%nGlobalParts
nPartsMPI(  0:nProcessors-1) = offsetPartMPI(1:nProcessors)-offsetPartMPI(0:nProcessors-1)
ListIn%nLocalParts = nPartsMPI(myRank)
ListIn%offsetPart  = offsetPartMPI(myRank)
#else /* MPI */
ListIn%nLocalParts = ListIn%nGlobalParts
ListIn%offsetPart  = 0
#endif /* MPI */

PartitionPartIsDone = .TRUE.

#if USE_MPI
DEALLOCATE(offsetPartMPI)
DEALLOCATE(nPartsMPI)
#endif /* MPI */

END SUBROUTINE PartitionPartMPI
#endif /*PARTICLES*/

END MODULE MOD_Posti_Part_Tools

