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

#if USE_LOADBALANCE && USE_MPI
INTERFACE DomainDecomposition
  MODULE PROCEDURE DomainDecomposition
END INTERFACE

PUBLIC::DomainDecomposition
#endif /*USE_LOADBALANCE && USE_MPI*/
!==================================================================================================================================

CONTAINS


#if USE_LOADBALANCE && USE_MPI
SUBROUTINE DomainDecomposition()
!===================================================================================================================================
!> Read ElemTime from .h5 container and compute domain decomposition
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars               ,ONLY: DomainDecompositionWallTime
USE MOD_HDF5_Input                 ,ONLY: OpenDataFile,CloseDataFile
USE MOD_IO_HDF5                    ,ONLY: AddToElemData,ElementOut
USE MOD_LoadDistribution           ,ONLY: ApplyWeightDistributionMethod,WeightDistribution_Equal
USE MOD_LoadBalance_Vars           ,ONLY: NewImbalance,MaxWeight,MinWeight,ElemGlobalTime,LoadDistri,TargetWeight
USE MOD_LoadBalance_Vars           ,ONLY: ElemTime,ProcTime,PerformLoadBalance
USE MOD_Mesh_Vars                  ,ONLY: offsetElem,nElems,nGlobalElems!,MeshFile
USE MOD_MPI_Vars                   ,ONLY: offsetElemMPI
USE MOD_Restart_Vars               ,ONLY: DoRestart
USE MOD_ReadInTools                ,ONLY: PrintOption
! Element Redistribution
USE MOD_LoadBalance_Vars           ,ONLY: MPInElemSend,MPIoffsetElemSend,MPInElemRecv,MPIoffsetElemRecv
USE MOD_LoadBalance_Vars           ,ONLY: nElemsOld,offsetElemOld
USE MOD_LoadBalance_Vars           ,ONLY: ElemInfoRank_Shared,ElemInfoRank_Shared_Win
USE MOD_Particle_Mesh_Vars         ,ONLY: ElemInfo_Shared,ElemInfo_Shared_Win
USE MOD_Particle_MPI_Shared        ,ONLY: Allocate_Shared,BARRIER_AND_SYNC
USE MOD_Particle_MPI_Shared_Vars   ,ONLY: myComputeNodeRank,MPI_COMM_SHARED
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
LOGICAL                        :: ElemTimeExists
REAL,ALLOCATABLE               :: WeightSum_proc(:)
INTEGER                        :: iProc
! Element Redistribution
INTEGER                        :: iElem,ElemRank,nElemsProc
INTEGER                        :: offsetElemSend,offsetElemRecv
! Timers
REAL                           :: StartT,EndT
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' DOMAIN DECOMPOSITION...'

GETTIME(StartT)

IF (PerformLoadBalance) THEN
  nElemsOld     = nElems
  offsetElemOld = offsetElem
  IF (myComputeNodeRank.EQ.0) &
    ElemInfoRank_Shared  = ElemInfo_Shared(ELEM_RANK,:)
  CALL BARRIER_AND_SYNC(ElemInfoRank_Shared_Win,MPI_COMM_SHARED)
END IF

SDEALLOCATE(LoadDistri)
ALLOCATE(LoadDistri(0:nProcessors-1))
LoadDistri(:)  = 0.
ElemTimeExists = .FALSE.

IF (DoRestart .OR. PerformLoadBalance) THEN
  !--------------------------------------------------------------------------------------------------------------------------------!
  ! Readin of ElemTime: Read in only by MPIRoot in single mode, only communicate logical ElemTimeExists
  ! because the root performs the distribution of elements (domain decomposition) due to the load distribution scheme
  SDEALLOCATE(ElemGlobalTime)
  ALLOCATE(ElemGlobalTime(1:nGlobalElems)) ! Allocate ElemGlobalTime for all MPI ranks
  ElemGlobalTime = 0.

  IF (PerformLoadBalance) THEN
    CALL ReadElemTime(single=.TRUE.)
  ! 1) Only MPIRoot does readin of ElemTime during restart
  ELSEIF (MPIRoot) THEN
    CALL ReadElemTime(single=.TRUE.)
  END IF

  IF (MPIRoot) THEN
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
  CALL WeightDistribution_Equal(nProcessors,nGlobalElems,offsetElemMPI)

  ! Send the load distribution to all other procs
  CALL MPI_BCAST(offsetElemMPI,nProcessors+1,MPI_INTEGER,0,MPI_COMM_FLEXI,iERROR)
END IF ! IF(DoRestart)

! Set local number of elements
nElems     = offsetElemMPI(myRank+1) - offsetElemMPI(myRank)
offsetElem = offsetElemMPI(myRank)
LOGWRITE(*,*)'offset,nElems',offsetElem,nElems

IF (PerformLoadBalance) THEN
  ! Only update the mapping of element to rank
  IF (myComputeNodeRank.EQ.0) THEN
    ! Shared array is allocated on compute-node level, compute-node root must update the mapping
    DO iProc = 0,nProcessors-1
      nElemsProc = offsetElemMPI(iProc+1) - offsetElemMPI(iProc)
      ElemInfo_Shared(ELEM_RANK,offsetElemMPI(iProc)+1:offsetElemMPI(iProc)+nElemsProc) = iProc
    END DO ! iProc = 0,nProcessors-1
  END IF ! myComputeNodeRank.EQ.0
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)

  ! Calculate the elements to send
  MPInElemSend      = 0
  MPIoffsetElemSend = 0
  ! Loop with the old element over the new elem distribution
  DO iElem = 1,nElemsOld
    ElemRank               = ElemInfo_Shared(ELEM_RANK,offsetElemOld+iElem)+1
    MPInElemSend(ElemRank) = MPInElemSend(ElemRank) + 1
  END DO

  offsetElemSend = 0
  DO iProc = 2,nProcessors
    MPIoffsetElemSend(iProc) = SUM(MPInElemSend(1:iProc-1))
  END DO

  ! Calculate the elements to recv
  MPInElemRecv      = 0
  MPIoffsetElemRecv = 0
  ! Loop with the new element over the old elem distribution
  DO iElem = 1,nElems
    ElemRank               = ElemInfoRank_Shared(offsetElem+iElem)+1
    MPInElemRecv(ElemRank) = MPInElemRecv(ElemRank) + 1
  END DO

  offsetElemRecv = 0
  DO iProc = 2,nProcessors
    MPIoffsetElemRecv(iProc) = SUM(MPInElemRecv(1:iProc-1))
  END DO
END IF

! Sanity check: local nElems and offset
IF (nElems.LE.0) CALL Abort(__STAMP__,' Process did not receive any elements/load! ')

! Read the ElemTime again, but this time with every proc, depending on the domain decomposition in order to write the data
! to the state file (keep ElemTime on restart, if no new ElemTime is calculated during the run or replace with newly measured values
! if LoadBalance is on)
IF (ElemTimeExists) CALL ReadElemTime(single=.FALSE.)

! Set new ElemTime depending on new load distribution
SDEALLOCATE(ElemTime)
SDEALLOCATE(ProcTime)
ALLOCATE(ElemTime(1:nElems),&
         ProcTime(1:nElems))
ElemTime = 0.
ProcTime = 0.
CALL AddToElemData(ElementOut,'ElemTime',ElemTime)
CALL AddToElemData(ElementOut,'ProcTime',ProcTime)

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

  ! Abort if invalid targetWeight
  IF(TargetWeight.LE.0.0) CALL Abort(__STAMP__,' LoadBalance: TargetWeight = ',RealInfo=TargetWeight)

  ! valid decomposition, output result
  SWRITE(UNIT_stdOut,'(132("."))')
  SWRITE(UNIT_stdOut,'(A)') ' Calculated new (theoretical) imbalance with offsetElemMPI information'
  CALL PrintOption('MaxWeight'   ,'CALC',RealOpt=MaxWeight)
  CALL PrintOption('MinWeight'   ,'CALC',RealOpt=MinWeight)
  CALL PrintOption('TargetWeight','CALC',RealOpt=TargetWeight)
  CALL PrintOption('NewImbalance','CALC',RealOpt=NewImbalance)
  SWRITE(UNIT_stdOut,'(132("."))')
  DEALLOCATE(WeightSum_proc)
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' No ElemTime found in restart file'
  NewImbalance = -1.
  MaxWeight    = -1.
  MinWeight    = -1.
END IF

! Re-open mesh file to continue readin. Meshfile is not set if this routine is called from posti
! IF (INDEX(MeshFile,'h5').NE.0)  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

EndT                        = FLEXITIME()
DomainDecompositionWallTime = EndT-StartT
SWRITE(UNIT_stdOut,'(A,F0.3,A)')' DOMAIN DECOMPOSITION... DONE! [',DomainDecompositionWallTime,'s]'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE DomainDecomposition


SUBROUTINE ReadElemTime(single)
!===================================================================================================================================
!> Read ElemTime from .h5 container either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Input             ,ONLY: ReadArray,DatasetExists
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime,ElemGlobalTime
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime_tmp
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_LoadBalance_Vars       ,ONLY: offsetElemMPIOld
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nElems,nGlobalElems
USE MOD_Restart_Vars           ,ONLY: RestartFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)  :: single !< read data file either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
LOGICAL             :: ElemTimeExists
INTEGER             :: iProc
REAL                :: StartT,EndT
REAL,ALLOCATABLE    :: ElemTimeTmp(:)
INTEGER             :: ElemPerProc(0:nProcessors-1)
!===================================================================================================================================

IF (PerformLoadBalance) THEN
  IF (single) THEN
    DO iProc = 0,nProcessors-1
      ElemPerProc(iProc) = offsetElemMPIOld(iProc+1) - offsetElemMPIOld(iProc)
    END DO
    CALL MPI_GATHERV(ElemTime,nElems,MPI_DOUBLE_PRECISION,ElemGlobalTime,ElemPerProc,offsetElemMPIOld(0:nProcessors-1),MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
  ELSE
    ALLOCATE(ElemTimeTmp(1:nElems))

    ASSOCIATE (&
            counts_send  => (MPInElemSend     ) ,&
            disp_send    => (MPIoffsetElemSend) ,&
            counts_recv  => (MPInElemRecv     ) ,&
            disp_recv    => (MPIoffsetElemRecv))
      ! Communicate PartInt over MPI
      CALL MPI_ALLTOALLV(ElemTime,counts_send,disp_send,MPI_DOUBLE_PRECISION,ElemTimeTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
    END ASSOCIATE

    DEALLOCATE(ElemTime)
    ALLOCATE(ElemTime(1:nElems))
    ElemTime = ElemTimeTmp
    DEALLOCATE(ElemTimeTmp)
  END IF
ELSE
  ! Read data file either single=.TRUE. (only MPI root) or single=.FALSE. (all ranks)
  IF (single) THEN
    WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')   ' | Reading ElemTime from restart file (single-mode): ',TRIM(RestartFile),' ...'
    GETTIME(StartT)

    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
    CALL DatasetExists(File_ID,'ElemTime',ElemTimeExists)
    IF (ElemTimeExists) &
      CALL ReadArray('ElemTime',2,(/1,nGlobalElems/),0,2,RealArray=ElemGlobalTime)
    CALL CloseDataFile()

    GETTIME(EndT)
    WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE! [',EndT-StartT,'s]'
  ELSE
    IF(MPIRoot)THEN
      WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO') ' | Reading ElemTime from restart file (MPI-mode)   : ',TRIM(RestartFile),' ...'
      GETTIME(StartT)
    END IF

    SDEALLOCATE(ElemTime_tmp)
    ALLOCATE(ElemTime_tmp(1:nElems))
    ElemTime_tmp  = 0.
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_FLEXI)
    CALL ReadArray('ElemTime',2,(/1,nElems/),OffsetElem,2,RealArray=ElemTime_tmp)
    CALL CloseDataFile()

    IF(MPIRoot)THEN
      GETTIME(EndT)
      WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE! [',EndT-StartT,'s]'
    END IF
  END IF ! single
END IF

END SUBROUTINE ReadElemTime
#endif /*USE_LOADBALANCE && USE_MPI*/

END MODULE MOD_LoadBalance_Tools
