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

INTERFACE InitPartState
  MODULE PROCEDURE InitPartState
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

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitPartState, ReadPartStateFile, FinalizeReadPartStateFile
PUBLIC :: InitParticleOutput,PartitionPartMPI
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
  ALLOCATE (ListIn%VarNamePartCombine(ListIn%nPartVar_Visu))
  SDEALLOCATE(ListIn%VarNamePartCombineLen)
  ALLOCATE (ListIn%VarNamePartCombineLen(ListIn%nPartVar_Visu))
  ListIn%VarNamePartCombine = 0
  DO iVar=2,ListIn%nPartVar_Visu
    i = LEN(TRIM(ListIn%VarNamePartVisu(iVar)))
    tmp = ListIn%VarNamePartVisu(iVar)
    tmp2 = ListIn%VarNamePartVisu(iVar-1)
    IF (TRIM(tmp(:i-1)) .EQ. TRIM(tmp2(:i-1))) THEN
      IF (ListIn%VarNamePartCombine(iVar-1) .EQ. 0) ListIn%VarNamePartCombine(iVar-1) = 1
      ListIn%VarNamePartCombine(iVar) = ListIn%VarNamePartCombine(iVar-1) + 1
    END IF
  END DO
  ListIn%VarNamePartCombineLen = 0
  ListIn%VarNamePartCombineLen(ListIn%nPartVar_visu) = ListIn%VarNamePartCombine(ListIn%nPartVar_visu)
  DO iVar=ListIn%nPartVar_visu-1,1,-1
    IF (ListIn%VarNamePartCombine(iVar).GT.0) THEN
      ListIn%VarNamePartCombineLen(iVar) = MAX(ListIn%VarNamePartCombine(iVar), ListIn%VarNamePartCombineLen(iVar+1)) 
    END IF
  END DO
END IF

END SUBROUTINE InitParticleOutput


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

IF(datasetNames.EQ.'PartData')THEN
  CALL GetVarNames("VarNamesParticles",varnames,VarNamesExist)
ELSE IF(datasetNames.EQ.'ErosionData')THEN
  CALL GetVarNames("VarNamesErosion",varnames,VarNamesExist)
END IF

IF(VarNamesExist)THEN
  CALL GetDataSize(File_ID,TRIM(datasetNames),dims,HSize)
  ListIn%nPartVar_HDF5=INT(HSize(1))-3
  ListIn%nGlobalParts=INT(HSize(2))
  CALL PartitionPartMPI(ListIn)
  ALLOCATE(ListIn%VarNamesPart_HDF5(1:ListIn%nPartVar_HDF5))
  ListIn%VarNamesPart_HDF5=varnames(4:)
  ALLOCATE(ListIn%VarNamePartDummy(1:ListIn%nPartVar_HDF5))
  ListIn%VarNamePartDummy=''
  ALLOCATE(ListIn%mapAllVarsToVisuVars(1:ListIn%nPartVar_HDF5))
  ListIn%mapAllVarsToVisuVars=0
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
INTEGER                                 :: iPart,iVar,iVar2, i
REAL,ALLOCATABLE                        :: PartData(:,:)
CHARACTER(LEN=255)                      :: tmp, tmp2
LOGICAL                                 :: datasetFound
!===================================================================================================================================
SDEALLOCATE(ListIn%PartData_HDF5)

IF(PRESENT(DataArrayIn))THEN
  DataArray = DataArrayIn
ELSE
  DataArray = 'PartData'
END IF

CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

CALL DatasetExists(File_ID,TRIM(DataArray),datasetFound)

IF(.NOT.datasetFound) RETURN

ALLOCATE(PartData(1:ListIn%nPartVar_HDF5+3,1:ListIn%nLocalParts))
CALL ReadArray(DataArray,2,(/ListIn%nPartVar_HDF5+3, ListIn%nLocalParts/),ListIn%offsetPart,2,RealArray=PartData)

ListIn%nPart_Visu=0

DO iPart=1,ListIn%nLocalParts
  ListIn%nPart_Visu=ListIn%nPart_Visu+1
END DO

IF(VisuPart)THEN
  ListIn%nPartVar_Visu=ListIn%nPartVar_HDF5
  ALLOCATE(ListIn%PartData_HDF5(1:ListIn%nPartVar_Visu+3,1:ListIn%nLocalParts))
  ListIn%PartData_HDF5=PartData
  SDEALLOCATE(ListIn%VarNamePartVisu)
  ALLOCATE(ListIn%VarNamePartVisu(1:ListIn%nPartVar_Visu))
  ListIn%VarNamePartVisu=ListIn%VarNamesPart_HDF5
  ListIn%mapAllVarsToVisuVars=1
ELSE
  ALLOCATE(ListIn%PartData_HDF5(1:ListIn%nPartVar_Visu+3,1:ListIn%nLocalParts))
  ListIn%PartData_HDF5(1:3,:)=PartData(1:3,:)
  DO iVar=1,ListIn%nPartVar_Visu
    DO iVar2=1,ListIn%nPartVar_HDF5
      IF(STRICMP(ListIn%VarNamePartVisu(iVar),ListIn%VarNamesPart_HDF5(iVar2)))THEN
        ListIn%PartData_HDF5(iVar+3,:)=PartData(iVar2+3,:)
        ListIn%mapAllVarsToVisuVars(iVar2)=1
      END IF
    END DO
  END DO
END IF

SDEALLOCATE(PartData)
CALL CloseDataFile()

! Needed for output
IF(ListIn%nPartVar_Visu.GT.0)THEN
  SDEALLOCATE(ListIn%VarNamePartCombine)
  ALLOCATE (ListIn%VarNamePartCombine(ListIn%nPartVar_Visu))
  SDEALLOCATE(ListIn%VarNamePartCombineLen)
  ALLOCATE (ListIn%VarNamePartCombineLen(ListIn%nPartVar_Visu))
  ListIn%VarNamePartCombine = 0
  DO iVar=2,ListIn%nPartVar_Visu
    i = LEN(TRIM(ListIn%VarNamePartVisu(iVar)))
    tmp = ListIn%VarNamePartVisu(iVar)
    tmp2 = ListIn%VarNamePartVisu(iVar-1)
    IF (TRIM(tmp(:i-1)) .EQ. TRIM(tmp2(:i-1))) THEN
      IF (ListIn%VarNamePartCombine(iVar-1) .EQ. 0) ListIn%VarNamePartCombine(iVar-1) = 1
      ListIn%VarNamePartCombine(iVar) = ListIn%VarNamePartCombine(iVar-1) + 1
    END IF
  END DO
  ListIn%VarNamePartCombineLen = 0
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
USE MOD_Particle_MPI_Vars        ,ONLY: PartitionPartIsDone, nPartsMPI, offsetPartMPI
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
INTEGER                                 :: iProc, iPart
#endif /* MPI */
!==================================================================================================================================
#if USE_MPI
ALLOCATE(offsetPartMPI(0:nProcessors),nPartsMPI(0:nProcessors-1))
offsetPartMPI=0
ListIn%nLocalParts=ListIn%nGlobalParts/nProcessors
iPart=ListIn%nGlobalParts-ListIn%nLocalParts*nProcessors
DO iProc=0,nProcessors-1  
  offsetPartMPI(iProc)=ListIn%nLocalParts*iProc+MIN(iProc,iPart)
END DO
offsetPartMPI(nProcessors)=ListIn%nGlobalParts
nPartsMPI(0:nProcessors-1)=offsetPartMPI(1:nProcessors)-offsetPartMPI(0:nProcessors-1)
ListIn%nLocalParts=nPartsMPI(myRank)
ListIn%offsetPart=offsetPartMPI(myRank)
#else /* MPI */
ListIn%nLocalParts=ListIn%nGlobalParts
ListIn%offsetPart=0
#endif /* MPI */
PartitionPartIsDone=.TRUE.
SDEALLOCATE(offsetPartMPI)
SDEALLOCATE(nPartsMPI)
END SUBROUTINE PartitionPartMPI
#endif

END MODULE MOD_Posti_Part_Tools

