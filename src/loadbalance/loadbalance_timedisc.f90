!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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

!===================================================================================================================================
!> Module contains the routines for load balancing
!> This file is needed only in FLEXI to avoid a circular depencendy
!===================================================================================================================================
MODULE MOD_LoadBalance_TimeDisc
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

PUBLIC :: LoadBalance
!===================================================================================================================================

CONTAINS

SUBROUTINE LoadBalance(OutputTime)
!===================================================================================================================================
! Routine performing the load balancing
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_Globals_Vars               ,ONLY: InitializationWallTime
USE MOD_Preproc
USE MOD_Analyze                    ,ONLY: InitAnalyze,FinalizeAnalyze
USE MOD_BaseFlow                   ,ONLY: InitBaseFlow,FinalizeBaseFlow
USE MOD_DG                         ,ONLY: InitDG,FinalizeDG
USE MOD_Equation                   ,ONLY: InitEquation,FinalizeEquation
USE MOD_Filter                     ,ONLY: InitFilter,FinalizeFilter
USE MOD_LoadBalance_BaseFlow       ,ONLY: BaseFlowRestart
USE MOD_LoadBalance_Restart        ,ONLY: FieldRestart
USE MOD_LoadBalance_Vars           ,ONLY: ElemTime,ElemTimeField
USE MOD_LoadBalance_Vars           ,ONLY: nLoadBalanceSteps,LoadBalanceMaxSteps,NewImbalance,MinWeight,MaxWeight
USE MOD_LoadBalance_Vars           ,ONLY: CurrentImbalance,MaxWeight,MinWeight
USE MOD_LoadBalance_Vars           ,ONLY: PerformLoadBalance
USE MOD_IO_HDF5                    ,ONLY: ElementOut,FieldOut,FinalizeIOHDF5
USE MOD_Mesh                       ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars                  ,ONLY: nElems
USE MOD_MPI                        ,ONLY: InitMPIVars,FinalizeMPI
USE MOD_Output_Vars                ,ONLY: ProjectName
USE MOD_Overintegration            ,ONLY: InitOverintegration,FinalizeOverintegration
USE MOD_Predictor                  ,ONLY: InitPredictor,FinalizePredictor
USE MOD_RecordPoints               ,ONLY: InitRecordPoints,FinalizeRecordPoints
USE MOD_ReadInTools                ,ONLY: prms
USE MOD_Restart                    ,ONLY: InitRestart,FinalizeRestart,Restart
USE MOD_Restart_Vars               ,ONLY: RestartFile,doRestart,RestartMode
USE MOD_Sponge                     ,ONLY: InitSponge,FinalizeSponge
USE MOD_StringTools                ,ONLY: set_formatting,clear_formatting
!USE MOD_TimeDisc                   ,ONLY: InitTimeDisc,FinalizeTimeDisc
USE MOD_TimeDisc_Vars              ,ONLY: dtElem
USE MOD_TimeDisc_Vars              ,ONLY: Ut_tmp,UPrev,S2
USE MOD_TimeDisc_Vars              ,ONLY: TimeDiscType,TimeDiscMethod
#if PARABOLIC
USE MOD_Lifting                    ,ONLY: InitLifting,FinalizeLifting
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV                         ,ONLY: InitFV,FinalizeFV
USE MOD_Indicator                  ,ONLY: InitIndicator,FinalizeIndicator
#endif /* FV_ENABLED */
#if USE_PARTICLES
USE MOD_LoadBalance_Vars           ,ONLY: ElemTimePart
USE MOD_Particle_Init              ,ONLY: InitParticles,FinalizeParticles
USE MOD_Particle_MPI               ,ONLY: InitParticleMPI,FinalizeParticleMPI
#if PARTICLES_COUPLING >= 2
USE MOD_LoadBalance_Vars           ,ONLY: ElemTimePartDepo
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
USE MOD_LoadBalance_Vars           ,ONLY: ElemTimePartColl
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/
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
  ElemTime               = 0.
  ElemTimeField          = 0.
  InitializationWallTime = 0.
#if USE_PARTICLES
  ElemTimePart           = 0.
#endif /*USE_PARTICLES*/
#if PARTICLES_COUPLING >= 2
  ElemTimePartDepo       = 0.
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
  ElemTimePartColl       = 0.
#endif /*PARTICLES_COUPLING*/
  RETURN
END IF

nLoadBalanceSteps  = nLoadBalanceSteps+1
SWRITE(UNIT_stdOut,'(132("-"))')
CALL set_formatting('green')
IF (LoadBalanceMaxSteps.GT.0) THEN
  SWRITE(UNIT_stdOut,'(A,I0,A,I0,A)') ' PERFORMING LOAD BALANCE ',nLoadBalanceSteps,' of ',LoadBalanceMaxSteps,' ...'
ELSE
  SWRITE(UNIT_stdOut,'(A,I0,A)'     ) ' PERFORMING LOAD BALANCE ',nLoadBalanceSteps,' ...'
END IF
CALL clear_formatting()
! Measure init duration
LB_StartTime=FLEXITIME()

! finalize all arrays
!CALL FinalizeFlexi(IsLoadBalance=.TRUE.)
! reallocate. FLEXI determines new imbalance in InitMesh() -> ReadMesh()
!CALL InitFlexi    (IsLoadBalance=.TRUE.)

! Either have to move InitFlexi & FinalizeFlexi into separate file or call the subroutines here manually. Chose the latter option to
! keep compatibility with master branch
!-- Save the last output file name as restart file name
RestartFile = TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//'State',OutputTime))//'.h5'
doRestart   = .TRUE.
RestartMode = 1

!-- Finalize every mesh dependent routine
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
#if PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/
#if FV_ENABLED
CALL FinalizeFV()
CALL FinalizeIndicator()
#endif /*PARABOLIC*/
CALL FinalizeDG()
CALL FinalizeEquation()
! Calling timedisc causes circular depends. We only need to reallocate arrays
!CALL FinalizeTimeDisc()
SDEALLOCATE(Ut_tmp)
SDEALLOCATE(S2)
SDEALLOCATE(UPrev)
SDEALLOCATE(dtElem)

CALL FinalizeMesh()
CALL FinalizeSponge()
CALL FinalizeBaseFlow()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
#if USE_PARTICLES
CALL FinalizeParticleMPI()
CALL FinalizeParticles()
#endif /*USE_PARTICLES*/
CALL FinalizeMPI()
CALL FinalizeIOHDF5()

!-- Set removed flag to false
CALL prms%finalize(.TRUE.)
ElementOut   => NULL() !< linked list of output pointers
FieldOut     => NULL() !< linked list of output pointers

!-- Restart every mesh dependent routine
CALL InitMesh(meshMode=2)
CALL InitFilter()
CALL InitOverintegration()
CALL InitMPIVars()
CALL InitEquation()
CALL InitBaseFlow()
CALL InitDG()
#if FV_ENABLED
CALL InitIndicator()
CALL InitFV()
#endif /*FV_ENABLED*/
#if PARABOLIC
CALL InitLifting()
#endif /*PARABOLIC*/

! Calling timedisc causes circular depends. We only need to reallocate arrays
!CALL InitTimeDisc()
SELECT CASE(TimeDiscType)
  CASE('LSERKW2')
    ALLOCATE(Ut_tmp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE('LSERKK3')
    ALLOCATE(S2   (1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) &
            ,UPrev(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
    ! Implicit time integration
  CASE('ESDIRK')
    ! Predictor for Newton
    CALL FinalizePredictor()
    CALL InitPredictor(TimeDiscMethod)
END SELECT
ALLOCATE(dtElem(nElems))
dtElem = 0.

CALL InitAnalyze()
CALL InitRecordpoints()
CALL BaseFlowRestart()
CALL FieldRestart()
CALL InitSponge()
#if USE_PARTICLES
CALL InitParticleMPI()
CALL InitParticles(doLoadBalance_opt=.TRUE.)
#endif /*USE_PARTICLES*/

! zero ElemTime, the measurement starts again
ElemTime         = 0.
#if USE_PARTICLES
ElemTimePart     = 0.
#if PARTICLES_COUPLING >= 2
ElemTimePartDepo = 0.
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
ElemTimePartColl = 0.
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/
ElemTimeField    = 0.

IF(NewImbalance.GT.CurrentImbalance) THEN
  CALL set_formatting('red')
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful! Determined new (theoretical) imbalance from elem distribution'
  CALL clear_formatting()
ELSE
  CALL set_formatting('green')
  SWRITE(UNIT_stdOut,'(A)') ' LoadBalance successful! Determined new (theoretical) imbalance from elem distribution'
  CALL clear_formatting()
END IF
SWRITE(UNIT_stdOut,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)')&
  ' MinWeight: ', MinWeight, '    MaxWeight: ', MaxWeight, '    OldImbalance: ', CurrentImbalance,'    NewImbalance: ', NewImbalance

! Calculate time spent for load balance restart
LB_Time                = FLEXITIME()
InitializationWallTime = LB_Time - LB_StartTime
CALL DisplayMessageAndTime(InitializationWallTime,'LOAD BALANCE DONE!',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')

END SUBROUTINE LoadBalance

END MODULE MOD_LoadBalance_TimeDisc
