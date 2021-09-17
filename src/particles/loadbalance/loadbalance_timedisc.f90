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
!> Module contains the routines for load balancing
!===================================================================================================================================
MODULE MOD_LoadBalance_TimeDisc
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_LOADBALANCE
INTERFACE InitialAutoRestart
  MODULE PROCEDURE InitialAutoRestart
END INTERFACE

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

PUBLIC :: InitialAutoRestart
PUBLIC :: LoadBalance
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

CONTAINS

#if USE_LOADBALANCE

SUBROUTINE InitialAutoRestart(t,dt,dt_min,tend,tAnalyze,Analyze_dt,tmp_LoadBalanceSample,tmp_DoLoadBalance)
!===================================================================================================================================
! routine to adjust analyze parameters to handle initial auto restart
!===================================================================================================================================
! USED MODULES
USE MOD_LoadBalance_Vars
USE MOD_Restart_Vars               ,ONLY: RestartTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)              :: t
REAL,INTENT(INOUT)           :: dt
REAL,INTENT(INOUT)           :: dt_min
REAL,INTENT(IN)              :: tend
REAL,INTENT(INOUT)           :: Analyze_dt
REAL,INTENT(INOUT)           :: tAnalyze
INTEGER,INTENT(INOUT)        :: tmp_LoadBalanceSample    !> loadbalance sample saved until initial autorestart ist finished
LOGICAL,INTENT(INOUT)        :: tmp_DoLoadBalance        !> loadbalance flag saved until initial autorestart ist finished
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: tEndDiff
!===================================================================================================================================

IF (.NOT.DoInitialAutoRestart) RETURN

tmp_DoLoadBalance     = DoLoadBalance
DoLoadBalance         = .TRUE.
tmp_LoadbalanceSample = LoadBalanceSample
LoadBalanceSample     = InitialAutoRestartSample

! correct initialautrestartSample if partweight_initialautorestart is enabled so tAnalyze is calculated correctly
! LoadBalanceSample still needs to be zero
IF (IAR_PerformPartWeightLB) InitialAutoRestartSample=1

! correction for first analyzetime due to auto initial restart
tEndDiff = tend - t
IF (MIN(RestartTime+Analyze_dt,tEnd,RestartTime+InitialAutoRestartSample*dt).LT.tAnalyze) THEN
  tAnalyze     = MIN(RestartTime*Analyze_dt,tEnd,RestartTime+InitialAutoRestartSample*dt)
  Analyze_dt   = tAnalyze-t
  dt           = MINVAL((/dt_Min,Analyze_dt,tEndDiff/))
END IF

END SUBROUTINE InitialAutoRestart


SUBROUTINE LoadBalance(OutputTime)
!===================================================================================================================================
! routine perfoming the load balancing
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_Globals_Vars               ,ONLY: InitializationWallTime
USE MOD_Preproc
USE MOD_Analyze                    ,ONLY: InitAnalyze,FinalizeAnalyze
USE MOD_DG                         ,ONLY: InitDG,FinalizeDG
USE MOD_Equation                   ,ONLY: InitEquation,FinalizeEquation
USE MOD_Filter                     ,ONLY: InitFilter,FinalizeFilter
USE MOD_Lifting                    ,ONLY: InitLifting,FinalizeLifting
USE MOD_LoadBalance_Vars           ,ONLY: ElemTime,nLoadBalanceSteps,NewImbalance,MinWeight,MaxWeight
USE MOD_LoadBalance_Vars           ,ONLY: CurrentImbalance,MaxWeight,MinWeight
USE MOD_LoadBalance_Vars           ,ONLY: PerformLoadBalance
USE MOD_Indicator                  ,ONLY: InitIndicator,FinalizeIndicator
USE MOD_IO_HDF5                    ,ONLY: ElementOut,FieldOut
USE MOD_Mesh                       ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars                  ,ONLY: nElems
USE MOD_MPI                        ,ONLY: InitMPIVars,FinalizeMPI
USE MOD_Output_Vars                ,ONLY: ProjectName
USE MOD_Overintegration            ,ONLY: InitOverintegration,FinalizeOverintegration
USE MOD_Particle_Init              ,ONLY: InitParticles,FinalizeParticles
USE MOD_Particle_MPI               ,ONLY: InitParticleMPI,FinalizeParticleMPI
USE MOD_RecordPoints               ,ONLY: InitRecordPoints,FinalizeRecordPoints
USE MOD_ReadInTools                ,ONLY: prms
USE MOD_Restart                    ,ONLY: InitRestart,FinalizeRestart,Restart
USE MOD_Restart_Vars               ,ONLY: RestartFile,doRestart
USE MOD_Sponge                     ,ONLY: InitSponge,FinalizeSponge
!USE MOD_TimeDisc                   ,ONLY: InitTimeDisc,FinalizeTimeDisc
USE MOD_TimeDisc_Vars              ,ONLY: dtElem
#if PARABOLIC
USE MOD_Lifting                    ,ONLY: InitLifting
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV                         ,ONLY: InitFV,FinalizeFV
#endif /* FV_ENABLED */
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
 REAL,INTENT(IN)                :: OutputTime     !< simulation time when output is performed
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: LB_Time,LB_StartTime
!===================================================================================================================================
! only do load-balance if necessary
IF (.NOT.PerformLoadBalance) THEN
  ElemTime = 0.
  InitializationWallTime = 0.
  RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE...'

! Record time spent for load balance restart
LB_StartTime=FLEXITIME()

nLoadBalanceSteps = nLoadBalanceSteps + 1

! finalize all arrays
!CALL FinalizeFlexi(IsLoadBalance=.TRUE.)
! reallocate. FLEXI determines new imbalance in InitMesh() -> ReadMesh()
!CALL InitFlexi    (IsLoadBalance=.TRUE.)

! Either have to move InitFlexi & FinalizeFlexi into separate file or call the subroutines here manually. Chose the latter option to
! keep compatibility with master branch
!-- Save the last output file name as restart file name
RestartFile = TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//'State',OutputTime))//'.h5'
doRestart   = .TRUE.

!-- Finalize every mesh dependent routine
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
#if PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/
#if FV_ENABLED
CALL FinalizeFV()
#endif /*PARABOLIC*/
CALL FinalizeDG()
CALL FinalizeEquation()
! Calling timedisc causes circular depends. We only need to reallocate dtElem
!CALL FinalizeTimeDisc()
CALL FinalizeRestart()
SDEALLOCATE(dtElem)

CALL FinalizeMesh()
CALL FinalizeSponge()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
CALL FinalizeIndicator()
CALL FinalizeParticleMPI()
CALL FinalizeParticles()
CALL FinalizeMPI()

!-- Set removed flag to false
CALL prms%finalize()
ElementOut   => NULL() !< linked list of output pointers
FieldOut     => NULL() !< linked list of output pointers

!-- Restart every mesh dependent routine
CALL InitMesh(meshMode=2)
CALL InitRestart(RestartFile)
CALL InitFilter()
CALL InitOverintegration()
CALL InitIndicator()
CALL InitMPIvars()
CALL InitEquation()
CALL InitDG()
#if FV_ENABLED
CALL InitFV()
#endif /*FV_ENABLED*/
#if PARABOLIC
CALL InitLifting()
#endif /*PARABOLIC*/
CALL InitSponge()

! Calling timedisc causes circular depends. We only need to reallocate dtElem
!CALL InitTimeDisc()
ALLOCATE(dtElem(nElems))
dtElem = 0.

CALL InitAnalyze()
CALL InitRecordpoints()
CALL Restart()
CALL InitParticleMPI()
CALL InitParticles(doLoadBalance_opt=.TRUE.)

! zero ElemTime, the measurement starts again
ElemTime = 0.

IF(NewImbalance.GT.CurrentImbalance) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful!'
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' LoadBalance successful!'
END IF
SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | OldImbalance:                                  ', CurrentImbalance
SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | NewImbalance:                                  ', NewImbalance
SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | MaxWeight:                                     ', MaxWeight
SWRITE(UNIT_stdOut,'(A,ES15.7)') ' | MinWeight:                                     ', MinWeight

! Calculate time spent for load balance restart
LB_Time=FLEXITIME()
InitializationWallTime = LB_Time - LB_StartTime
!SWRITE(UNIT_stdOut,'(A,F14.2,A)') '  INITIALIZATION DONE IN [',InitializationWallTime,' sec ]'
!SWRITE(UNIT_stdOut,'(A)')         ' LOAD BALANCE DONE!'
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' LOAD BALANCE DONE IN [',InitializationWallTime,' sec ]'
SWRITE(UNIT_StdOut,'(132("-"))')

PerformLoadBalance = .FALSE.

END SUBROUTINE LoadBalance
#endif /*USE_LOADBALANCE*/

END MODULE MOD_LoadBalance_TimeDisc
