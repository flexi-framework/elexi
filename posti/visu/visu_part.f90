!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
!=================================================================================================================================
! Routines needed to visualize particles
!=================================================================================================================================
MODULE MOD_Posti_Part_Tools
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

#if USE_PARTICLES
INTERFACE InitParticle
  MODULE PROCEDURE InitParticle
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
PUBLIC :: InitParticle,PartitionPartMPI
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticle(ListIn,notRequestInformation)
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
LOGICAL,INTENT(IN),OPTIONAL           :: notRequestInformation
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: i, iVar
character(len=255)                    :: tmp
character(len=255)                    :: tmp2
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ...'

SWRITE(UNIT_StdOut,*) 'Preparing paricle variables from state file...'

! In case state variables are missing: use varnames from hDF5 file instead (used for timeavg-files)
SWRITE(UNIT_StdOut,*) 'WARNING: Not all Conservative Variables available in State File,'
SWRITE(UNIT_StdOut,*) '         cannot calculate any derived quantities.'
SWRITE(UNIT_StdOut,*) '         Visualizing State File Variables instead!'

! visualize all state file variables
ListIn%nPartVar_visu   = ListIn%nPartVar_HDF5
SDEALLOCATE(ListIn%VarNamePartVisu)
ALLOCATE(ListIn%VarNamePartVisu(ListIn%nPartVar_visu))

! Create state variables
ListIn%VarNamePartVisu(1:ListIn%nPartVar_HDF5) = ListIn%VarNamesPart_HDF5

SDEALLOCATE(ListIn%VarNamePartCombine)
ALLOCATE (ListIn%VarNamePartCombine(ListIn%nPartVar_visu))
SDEALLOCATE(ListIn%VarNamePartCombineLen)
ALLOCATE (ListIn%VarNamePartCombineLen(ListIn%nPartVar_visu))
ListIn%VarNamePartCombine = 0
DO iVar=2,ListIn%nPartVar_visu
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

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticle


SUBROUTINE InitPartState(InputFile,DataArrayIn,ListIn)
!===================================================================================================================================
! Read in main attributes from given HDF5 state file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Visu_Vars              ,ONLY: tVisuParticle, SpeciesID
USE MOD_IO_HDF5                ,ONLY: File_ID
USE MOD_HDF5_Input             ,ONLY: DatasetExists, GetDataSize, OpenDataFile, CloseDataFile, HSize, ReadAttribute
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
INTEGER                                 :: dims
LOGICAL                                 :: PartDataFound
CHARACTER(LEN=255)                      :: DataArray, Dataset
CHARACTER(LEN=255),ALLOCATABLE          :: tmpArray(:)
!===================================================================================================================================

IF(PRESENT(DataArrayIn))THEN
  DataArray = DataArrayIn
ELSE
  DataArray = 'PartData'
END IF

Dataset=''

! Get number and names of variables
CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,DataArray,PartDataFound)
IF(PartDataFound)THEN
  CALL GetDataSize(File_ID,DataArray,dims,HSize) 
  ListIn%nPartVar_HDF5=HSize(1)-3 ! remove position x,y,z 
  ListIn%nTotalParts_HDF5=HSize(2) 
  ALLOCATE(ListIn%VarNamesPart_HDF5(ListIn%nPartVar_HDF5))
  ALLOCATE(tmpArray(ListIn%nPartVar_HDF5+3))
   
  ! Read erosion names here
  SELECT CASE(DataArray)
      CASE('ErosionData')
          CALL ReadAttribute(File_ID,'VarNamesErosion',ListIn%nPartVar_HDF5+3,TRIM(DataSet),StrArray=tmpArray)
      CASE DEFAULT
          CALL ReadAttribute(File_ID,'VarNamesParticles',ListIn%nPartVar_HDF5+3,TRIM(DataSet),StrArray=tmpArray)
  END SELECT

  ListIn%VarNamesPart_HDF5(1:ListIn%nPartVar_HDF5)=tmpArray(4:ListIn%nPartVar_HDF5+3)
  DEALLOCATE(tmpArray)
  ListIn%nGlobalParts=ListIn%nTotalParts_HDF5
  CALL PartitionPartMPI(ListIn)
ELSE
  ListIn%nTotalParts_HDF5=0
  ListIn%nGlobalParts=0
  ListIn%nPartVar_HDF5=0
  ListIn%nGlobalParts=ListIn%nTotalParts_HDF5
END IF

CALL CloseDataFile()

END SUBROUTINE InitPartState


SUBROUTINE ReadPartStateFile(InputFile, DataArrayIn, ListIn)
!===================================================================================================================================
! Read in particle state from given HDF5 state file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Visu_Vars              ,ONLY: tVisuParticle
USE MOD_IO_HDF5                ,ONLY: File_ID
USE MOD_HDF5_Input             ,ONLY: DatasetExists, ReadArray, OpenDataFile, CloseDataFile, HSize
USE MOD_Particle_Vars          ,ONLY: PDM
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
LOGICAL                                 :: PartDataFound
CHARACTER(LEN=255)                      :: DataArray
INTEGER                                 :: iPart,iVar
!===================================================================================================================================
CALL FinalizeReadPartStateFile(ListIn)

IF(PRESENT(DataArrayIn))THEN
  DataArray = DataArrayIn
ELSE
  DataArray = 'PartData'
END IF

CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,DataArray,PartDataFound)
IF(.NOT.PartDataFound)THEN
  ListIn%nPart_Visu=0
  RETURN
END IF

ALLOCATE(ListIn%PartData_HDF5(1:ListIn%nPartVar_HDF5+3,1:ListIn%nLocalParts))
CALL ReadArray(DataArray,2,(/ListIn%nPartVar_HDF5+3, ListIn%nLocalParts/),ListIn%offsetPart,2,RealArray=ListIn%PartData_HDF5)

ALLOCATE(ListIn%VisualizePart(1:ListIn%nLocalParts))
ListIn%VisualizePart=.FALSE.
ListIn%nPart_Visu=0

DO iPart=1,ListIn%nLocalParts
  ListIn%nPart_Visu=ListIn%nPart_Visu+1
  ListIn%VisualizePart(iPart)=.TRUE.
END DO

CALL CloseDataFile()

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
SDEALLOCATE(ListIn%VarNamesPart_HDF5)
SDEALLOCATE(ListIn%PartData_HDF5)
SDEALLOCATE(ListIn%VisualizePart)
END SUBROUTINE FinalizeReadPartStateFile

!SUBROUTINE CalcPartEquation()
!!==================================================================================================================================
!!
!!==================================================================================================================================
!! MODULES
!USE MOD_Globals
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!==================================================================================================================================
!
!SWRITE(UNIT_StdOut,'(132("-"))')
!SWRITE(UNIT_stdOut,'(A)')" CONVERT PARTICLE DERIVED QUANTITIES..."
!
!ALLOCATE(Part_calc(PartQ%nCalc))
!
!END SUBROUTINE CalcPartEquation


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
!SDEALLOCATE(offsetPartMPI)
!SDEALLOCATE(nPartsMPI)
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

