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
#include "particle.h"

!===================================================================================================================================
!> Module contains the tools for load_balancing
!===================================================================================================================================
MODULE MOD_LoadBalance_Tools
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE

#if USE_MPI
INTERFACE DomainDecomposition
  MODULE PROCEDURE DomainDecomposition
END INTERFACE

PUBLIC::DomainDecomposition
#endif /*USE_MPI*/
!==================================================================================================================================

CONTAINS


#if USE_MPI
SUBROUTINE DomainDecomposition()
!===================================================================================================================================
!> Read ElemTime from .h5 container and compute domain decomposition
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars         ,ONLY: DomainDecompositionWallTime
USE MOD_LoadDistribution     ,ONLY: ApplyWeightDistributionMethod
USE MOD_LoadBalance_Vars     ,ONLY: NewImbalance,MaxWeight,MinWeight,ElemGlobalTime,LoadDistri,PartDistri,TargetWeight,ElemTime
USE MOD_HDF5_Input           ,ONLY: OpenDataFile
USE MOD_IO_HDF5              ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars            ,ONLY: MeshFile,offsetElem,nElems,nGlobalElems
USE MOD_MPI_Vars             ,ONLY: offsetElemMPI
USE MOD_Restart_Vars         ,ONLY: DoRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
LOGICAL                        :: ElemTimeExists
REAL,ALLOCATABLE               :: WeightSum_proc(:)
INTEGER                        :: iProc
INTEGER                        :: iElem
REAL                           :: StartT,EndT
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(UNIT_stdOut,'(A)')' DOMAIN DECOMPOSITION ...'

#if USE_MPI
StartT=MPI_WTIME()
#else
CALL CPU_TIME(StartT)
#endif

!simple partition: nGlobalelems/nprocs, do this on proc 0
SDEALLOCATE(offsetElemMPI)
ALLOCATE(offsetElemMPI(0:nProcessors))
offsetElemMPI  = 0
SDEALLOCATE(LoadDistri)
ALLOCATE(LoadDistri(   0:nProcessors-1))
LoadDistri(:)  = 0.
SDEALLOCATE(PartDistri)
ALLOCATE(PartDistri(   0:nProcessors-1))
PartDistri(:)  = 0
ElemTimeExists = .FALSE.

IF (DoRestart) THEN
  !--------------------------------------------------------------------------------------------------------------------------------!
  ! Readin of ElemTime: Read in only by MPIRoot in single mode, only communicate logical ElemTimeExists
  ! because the root performs the distribution of elements (domain decomposition) due to the load distribution scheme
  SDEALLOCATE(ElemGlobalTime)
  ALLOCATE(ElemGlobalTime(1:nGlobalElems)) ! Allocate ElemGlobalTime for all MPI ranks
  ElemGlobalTime = 0.

  ! 1) Only MPIRoot does readin of ElemTime
  IF (MPIRoot) THEN
    ! read ElemTime by root only
    CALL ReadElemTime(single=.TRUE.)

    ! if the elemtime is 0.0, the value must be changed in order to prevent a division by zero
    IF (MAXVAL(ElemGlobalTime).LE.0.) THEN
      ElemGlobalTime = 1.0
      ElemTimeExists = .FALSE.
    ELSE
      ElemTimeExists = .TRUE.
    END IF
  END IF

  ! 2) Distribute logical information ElemTimeExists
  CALL MPI_BCAST(ElemTimeExists,1,MPI_LOGICAL,0,MPI_COMM_FLEXI,iError)

  ! Distribute the elements according to the selected distribution method
  CALL ApplyWeightDistributionMethod(ElemTimeExists)
ELSE
  ! Simple partition: nGlobalelems/nProcessors
  nElems = nGlobalElems/nProcessors
  iElem  = nGlobalElems-nElems*nProcessors
  DO iProc = 0,nProcessors-1
    offsetElemMPI(iProc) = nElems*iProc + MIN(iProc,iElem)
  END DO
  offsetElemMPI(nProcessors) = nGlobalElems
END IF ! IF(DoRestart)

! Set local number of elements
nElems     = offsetElemMPI(myRank+1) - offsetElemMPI(myRank)
offsetElem = offsetElemMPI(myRank)
LOGWRITE(*,*)'offset,nElems',offsetElem,nElems

! Sanity check: local nElems and offset
IF (nElems.LE.0) CALL ABORT(__STAMP__,' Process did not receive any elements/load! ')

! Read the ElemTime again, but this time with every proc, depending on the domain decomposition in order to write the data
! to the state file (keep ElemTime on restart, if no new ElemTime is calculated during the run or replace with newly measured values
! if LoadBalance is on)
IF (ElemTimeExists) CALL ReadElemTime(single=.FALSE.)

! Set new ElemTime depending on new load distribution
SDEALLOCATE(ElemTime)
ALLOCATE(ElemTime(1:nElems))
ElemTime = 0.
CALL AddToElemData(ElementOut,'ElemTime',ElemTime)
!CALL AddToElemData(ElementOut,'ElemTime',RealArray=ElemTime(1:nElems))

! Calculate new (theoretical) imbalance with offsetElemMPI information
IF (ElemTimeExists.AND.MPIRoot) THEN
  ALLOCATE(WeightSum_proc(0:nProcessors-1))
  DO iProc=0,nProcessors-1
    WeightSum_proc(iProc) = SUM(ElemGlobalTime(1+offsetElemMPI(iProc):offsetElemMPI(iProc+1)))
  END DO
  SDEALLOCATE(ElemGlobalTime)
  MaxWeight = MAXVAL(WeightSum_proc)
  MinWeight = MINVAL(WeightSum_proc)
  ! WeightSum (Mesh global value) is already set in BalanceMethod scheme

  ! new computation of current imbalance
  TargetWeight = SUM(WeightSum_proc)/nProcessors
  NewImbalance = (MaxWeight-TargetWeight )/TargetWeight

  ! abort if invalid targetWeight
  IF(TargetWeight.LE.0.0) CALL abort(__STAMP__,' LoadBalance: TargetWeight = ',RealInfo=TargetWeight)

  ! valid decomposition, output result
  SWRITE(UNIT_stdOut,'(A)') ' Calculated new (theoretical) imbalance with offsetElemMPI information'
  SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | MaxWeight:                                     ', MaxWeight
  SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | MinWeight:                                     ', MinWeight
  SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | TargetWeight:                                  ', TargetWeight
  SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | NewImbalance:                                  ', NewImbalance
  DEALLOCATE(WeightSum_proc)
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' No ElemTime found in restart file'
  NewImbalance = -1.
  MaxWeight    = -1.
  MinWeight    = -1.
END IF

! Re-open mesh file to continue readin. Meshfile is not set if this routine is called from posti
IF (INDEX(MeshFile,'h5').NE.0)  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

EndT                        = FLEXITIME()
DomainDecompositionWallTime = EndT-StartT
SWRITE(UNIT_stdOut,'(A,F0.3,A)')' DOMAIN DECOMPOSITION ... DONE  [',DomainDecompositionWallTime,'s]'

END SUBROUTINE DomainDecomposition


SUBROUTINE ReadElemTime(single)
!===================================================================================================================================
!> Read ElemTime from .h5 container either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Input       ,ONLY: ReadArray,DatasetExists
USE MOD_LoadBalance_Vars ,ONLY: ElemGlobalTime
USE MOD_LoadBalance_Vars ,ONLY: ElemTime_tmp
USE MOD_Mesh_Vars        ,ONLY: offsetElem,nElems,nGlobalElems
USE MOD_Restart_Vars     ,ONLY: RestartFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)  :: single !< read data file either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
LOGICAL             :: ElemTimeExists
!===================================================================================================================================

! Read data file either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
IF (single) THEN
  nElems         = nGlobalElems   ! Temporarily set nElems as nGlobalElems for GetArrayAndName
  offsetElem     = 0              ! Offset is the index of first entry, hdf5 array starts at 0-.GT. -1

  ! NEW method
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  CALL DatasetExists(File_ID,'ElemTime',ElemTimeExists)
  IF (ElemTimeExists) THEN
    CALL ReadArray('ElemTime',2,(/1,nGlobalElems/),0,2,RealArray=ElemGlobalTime)
    WRITE(UNIT_stdOut,*) "Read ElemTime from restart file: "//TRIM(RestartFile)
  END IF ! ElemTimeExists
  CALL CloseDataFile()

ELSE
  SDEALLOCATE(ElemTime_tmp)
  ALLOCATE(ElemTime_tmp(1:nElems))
  ElemTime_tmp  = 0.
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_FLEXI)
  CALL ReadArray('ElemTime',2,(/1,nElems/),OffsetElem,2,RealArray=ElemTime_tmp)
  CALL CloseDataFile()
END IF ! single

END SUBROUTINE ReadElemTime
#endif /*USE_MPI*/


END MODULE MOD_LoadBalance_Tools
