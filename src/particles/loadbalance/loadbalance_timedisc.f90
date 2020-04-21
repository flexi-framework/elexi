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

#if USE_MPI
INTERFACE InitialAutoRestart
  MODULE PROCEDURE InitialAutoRestart
END INTERFACE

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

PUBLIC :: InitialAutoRestart
PUBLIC :: LoadBalance
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

#if USE_MPI

SUBROUTINE InitialAutoRestart(t,dt,dt_min,tend,tAnalyze,Analyze_dt,tmp_LoadBalanceSample,tmp_DoLoadBalance)
!===================================================================================================================================
! routine to adjust analyze parameters to handle initial auto restart
!===================================================================================================================================
! USED MODULES
USE MOD_LoadBalance_Vars
USE MOD_Restart_Vars           ,ONLY: RestartTime
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


SUBROUTINE LoadBalance()
!===================================================================================================================================
! routine perfoming the load balancing
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime,nLoadBalanceSteps,NewImbalance,MinWeight,MaxWeight
USE MOD_LoadBalance_Vars       ,ONLY: CurrentImbalance,MaxWeight,MinWeight
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
USE MOD_Particle_Globals       ,ONLY: InitializationWallTime
!USE MOD_Restart                ,ONLY: Restart
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: LB_Time,LB_StartTime
!===================================================================================================================================
! only do load-balance if necessary
IF (.NOT.PerformLoadBalance) THEN
  ElemTime=0.
  InitializationWallTime = 0.
  RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE ...'

! Record time spent for load balance restart
LB_StartTime=FLEXITIME()

nLoadBalanceSteps = nLoadBalanceSteps + 1

! finalize all arrays
!CALL FinalizeFlexi(IsLoadBalance=.TRUE.)
! reallocate. FLEXI determines new imbalance in InitMesh() -> ReadMesh()
!CALL InitFlexi    (IsLoadBalance=.TRUE.)

! Either have to move InitFlexi & FinalizeFlexi into separate file or call the subroutines here manually. Chose the latter option to
! keep compatibility with master branch

! restart
!CALL Restart()

! zero ElemTime, the measurement starts again
ElemTime     = 0.

IF(NewImbalance.GT.CurrentImbalance) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful!'
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' LoadBalance successful!'
END IF
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' OldImbalance: ', CurrentImbalance
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' NewImbalance: ', NewImbalance
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MaxWeight:    ', MaxWeight
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MinWeight:    ', MinWeight

! Calculate time spent for load balance restart
LB_Time=FLEXITIME()
InitializationWallTime = LB_Time - LB_StartTime
SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' INITIALIZATION DONE! [',InitializationWallTime,' sec ]'
SWRITE(UNIT_stdOut,'(A)')         ' LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE LoadBalance
#endif /*USE_MPI*/

END MODULE MOD_LoadBalance_TimeDisc