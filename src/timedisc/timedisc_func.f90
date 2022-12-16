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
!> Module for the GTS Temporal discretization
!==================================================================================================================================
MODULE MOD_TimeDisc_Functions
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersTimeDisc
  MODULE PROCEDURE DefineParametersTimeDisc
END INTERFACE

INTERFACE InitTimeDisc
  MODULE PROCEDURE InitTimeDisc
END INTERFACE

INTERFACE InitTimeStep
  MODULE PROCEDURE InitTimeStep
END INTERFACE

INTERFACE UpdateTimeStep
  MODULE PROCEDURE UpdateTimeStep
END INTERFACE

INTERFACE AnalyzeTimeStep
  MODULE PROCEDURE AnalyzeTimeStep
END INTERFACE

INTERFACE TimeDisc_Info
  MODULE PROCEDURE TimeDisc_Info
END INTERFACE

INTERFACE FinalizeTimeDisc
  MODULE PROCEDURE FinalizeTimeDisc
END INTERFACE

PUBLIC :: DefineParametersTimeDisc
PUBLIC :: InitTimeDisc
PUBLIC :: InitTimeStep
PUBLIC :: UpdateTimeStep
PUBLIC :: AnalyzeTimeStep
PUBLIC :: TimeDisc_Info
PUBLIC :: FinalizeTimeDisc
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersTimeDisc()
! MODULES
USE MOD_ReadInTools         ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("TimeDisc")
CALL prms%CreateStringOption('TimeDiscMethod', "Specifies the type of time-discretization to be used, e.g. the name of&
                                               & a specific Runge-Kutta scheme. Possible values:\n"//&
                                               "  * standardrk3-3\n  * carpenterrk4-5\n  * niegemannrk4-14\n"//&
                                               "  * toulorgerk4-8c\n  * toulorgerk3-7c\n  * toulorgerk4-8f\n"//&
                                               "  * ketchesonrk4-20\n  * ketchesonrk4-18\n  * eulerimplicit\n"//&
                                               "  * cranknicolson2-2\n  * esdirk2-3\n  * esdirk3-4\n"//&
                                               "  * esdirk4-6" , value='CarpenterRK4-5')
CALL prms%CreateRealOption(  'TEnd',           "End time of the simulation (mandatory).")
CALL prms%CreateRealOption(  'TStart',         "Start time of the simulation (optional, conflicts with restart).","0.0")
CALL prms%CreateRealOption(  'CFLScale',       "Scaling factor for the theoretical CFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateRealOption(  'DFLScale',       "Scaling factor for the theoretical DFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateRealOption(  'dtmin',          "Minimal allowed timestep (optional)","-1.0")
CALL prms%CreateRealOption(  'dtkill',         "Kill FLEXI if dt gets below this value (optional)","-1.0")
CALL prms%CreateIntOption(   'maxIter',        "Stop simulation when specified number of timesteps has been performed.", value='-1')
CALL prms%CreateIntOption(   'nCalcTimeStepMax',"Compute dt at least after every Nth timestep.", value='1')
#if USE_PARTICLES
CALL prms%CreateStringOption('ParticleTimeDiscMethod', "Specifies the type of particle time-discretization to be used."//&
                                                       "Possible values:\n:"//&
                                                       "* Runge-Kutta: Use the standard Runge-Kutta from the DG solver\n"//&
                                                       "* Euler: Use the low order Euler scheme.")
#endif
END SUBROUTINE DefineParametersTimeDisc


!==================================================================================================================================
!> Get information for end time and max time steps from ini file
!==================================================================================================================================
SUBROUTINE InitTimeDisc()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars         ,ONLY: NFilter,FilterType
USE MOD_IO_HDF5             ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Overintegration_Vars,ONLY: NUnder
USE MOD_Predictor           ,ONLY: InitPredictor
USE MOD_ReadInTools         ,ONLY: GETREAL,GETINT,GETSTR
USE MOD_StringTools         ,ONLY: LowCase,StripSpaces
USE MOD_TimeDisc_Vars       ,ONLY: b_dt,CFLScale,dtElem,dt,tend
USE MOD_TimeDisc_Vars       ,ONLY: tStart,dt_dynmin,dt_kill
USE MOD_TimeDisc_Vars       ,ONLY: Ut_tmp,UPrev,S2
USE MOD_TimeDisc_Vars       ,ONLY: maxIter,nCalcTimeStepMax
USE MOD_TimeDisc_Vars       ,ONLY: SetTimeDiscCoefs,TimeStep,TimeDiscName,TimeDiscType,TimeDiscInitIsDone,nRKStages
USE MOD_TimeDisc_Vars       ,ONLY: TimeDiscMethod
USE MOD_TimeStep            ,ONLY: TimeStepByLSERKW2,TimeStepByLSERKK3,TimeStepByESDIRK
#if PARABOLIC
USE MOD_TimeDisc_Vars       ,ONLY: DFLScale
#endif /*PARABOLIC*/
#if USE_PARTICLES
USE MOD_Particle_Boundary_Vars,ONLY: WriteMacroSurfaceValues,MacroValSampTime
USE MOD_TimeDisc_Vars       ,ONLY: t
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: NEff
!==================================================================================================================================
IF(TimeDiscInitIsDone)THEN
  SWRITE(UNIT_stdOut,'(A)') "InitTimeDisc already called."
  RETURN
END IF

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TIMEDISC...'

! Read max number of iterations to perform
maxIter = GETINT('maxIter')

! Get nCalcTimeStepMax: check if advanced settings should be used!
nCalcTimeStepMax = GETINT('nCalcTimeStepMax')

! Get TimeDisc Method
TimeDiscMethod = GETSTR('TimeDiscMethod')
CALL StripSpaces(TimeDiscMethod)
CALL LowCase(TimeDiscMethod)
CALL SetTimeDiscCoefs(TimeDiscMethod)

SELECT CASE(TimeDiscType)
  CASE('LSERKW2')
    TimeStep=>TimeStepByLSERKW2
    ALLOCATE(Ut_tmp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE('LSERKK3')
    TimeStep=>TimeStepByLSERKK3
    ALLOCATE(S2   (1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) &
            ,UPrev(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
  CASE('ESDIRK')
    ! Implicit time integration
    TimeStep=>TimeStepByESDIRK
    ! Predictor for Newton
    CALL InitPredictor(TimeDiscMethod)
END SELECT

ALLOCATE(b_dt(1:nRKStages))

! Read the end time TEnd from ini file
TEnd     = GETREAL('TEnd')

! Read the end time TEnd from ini file
TStart   = GETREAL('TStart')

! Read the normalized CFL number
CFLScale = GETREAL('CFLScale')
#if PARABOLIC
! Read the normalized DFL number
DFLScale = GETREAL('DFLScale')
#endif /*PARABOLIC*/
NEff     = MIN(PP_N,NFilter,NUnder)
IF(FilterType.GT.2) NEff = PP_N!LAF,HESTHAVEN no timestep effect
CALL fillCFL_DFL(NEff,PP_N)
! Read in minimal timestep
dt_dynmin = GETREAL("dtmin")
! Read in kill timestep
dt_kill   = GETREAL("dtkill")

! Set timestep to a large number
dt=HUGE(1.)

! Allocate and Initialize elementwise dt array
ALLOCATE(dtElem(nElems))
dtElem=0.

CALL AddToElemData(ElementOut,'dt',dtElem)

#if USE_PARTICLES
IF(WriteMacroSurfaceValues) MacroValSampTime = t
#endif

SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: '//TRIM(TimeDiscName)

TimeDiscInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TIMEDISC DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitTimeDisc


!===================================================================================================================================
!> Initial time step calculation for new timedisc-loop
!===================================================================================================================================
SUBROUTINE InitTimeStep()
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY: t,tAnalyze,tEnd,dt,dt_min,dt_minOld
USE MOD_TimeDisc_Vars       ,ONLY: ViscousTimeStep,CalcTimeStart,nCalcTimeStep
USE MOD_TimeDisc_Vars       ,ONLY: doAnalyze,doFinalize
#if USE_PARTICLES
USE MOD_Particle_Vars       ,ONLY: nSpecies
USE MOD_Particle_TimeDisc_Vars,ONLY: UseManualTimeStep,ManualTimeStep
#endif /*USE_PARTICLES*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars    ,ONLY: DoLoadBalanceBackup,LoadBalanceSampleBackup,DoLoadBalance
USE MOD_LoadBalance_Vars    ,ONLY: LoadBalanceSample,PerformLBSample
USE MOD_LoadBalance_Vars    ,ONLY: DoInitialAutoRestart,InitialAutoRestartSample
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: errType
!===================================================================================================================================
#if USE_PARTICLES
! Check if we are running in SteadyState mode
IF (UseManualTimeStep) dt = ManualTimeStep

! Get time step if needed. For DG time stepping, this is already calculated in BuildBGMAndIdentifyHaloRegion
IF ((UseManualTimeStep .AND. dt.EQ.0.) .OR. nSpecies.EQ.0) THEN
#endif /*USE_PARTICLES*/
dt      = EvalInitialTimeStep(errType)
#if USE_PARTICLES
ELSE
errType = 0
END IF
#endif /*USE_PARTICLES*/
dt_min(DT_MIN)     = dt
dt_min(DT_ANALYZE) = tAnalyze-t             ! Time to next analysis, put in extra variable so number does not change due to numerical errors
dt_min(DT_END)     = tEnd    -t             ! Do the same for end time
dt                 = MINVAL(dt_min)

IF (dt.EQ.dt_min(DT_ANALYZE))       doAnalyze  = .TRUE.
IF (dt.EQ.dt_min(DT_END    )) THEN; doAnalyze  = .TRUE.; doFinalize = .TRUE.; END IF
dt                 = MINVAL(dt_min,MASK=dt_min.GT.0)

nCalcTimeStep = 0
dt_minOld     = -999.
IF (errType.NE.0) CALL Abort(__STAMP__,&
#if EQNSYSNR == 3
  'Error: (1) density, (2) convective / (3) viscous timestep / muTilde (4) is NaN. Type/time:',errType,t)
#else
  'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)
#endif

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,ES16.7)') ' Initial Timestep  : ', dt
IF(ViscousTimeStep.AND.MPIRoot) WRITE(UNIT_stdOut,'(A)') ' Viscous timestep dominates! '

#if USE_LOADBALANCE
IF (DoInitialAutoRestart) THEN

  ! Set general load balance flag ON
  DoLoadBalanceBackup   = DoLoadBalance ! Backup
  DoLoadBalance         = .TRUE.        ! Force TRUE (during automatic restart, original variable might be false)

  ! Backup number of samples required for each load balance
  LoadBalanceSampleBackup = LoadBalanceSample        ! Backup: this is zero when PerformPartWeightLB=.TRUE.
  LoadBalanceSample       = InitialAutoRestartSample ! this is zero when InitialAutoRestartPartWeight=.TRUE.

  ! Activate sampling in first time step
  PerformLBSample       = .TRUE.

  ! Sanity check: initial automatic restart must happen before tAnalyze is reached (tAnalyze < LoadBalanceSample*dt not implemented)
  DO WHILE(LoadBalanceSample*dt.GT.dt_Min(DT_ANALYZE).AND.(LoadBalanceSample.GT.1))
   LoadBalanceSample = LoadBalanceSample-1
  END DO

END IF
#endif /*USE_LOADBALANCE*/

CalcTimeStart = FLEXITIME()

END SUBROUTINE InitTimeStep


!===================================================================================================================================
!> Update time step at the beginning of each timedisc loop
!===================================================================================================================================
SUBROUTINE UpdateTimeStep()
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars        ,ONLY: tWriteData
USE MOD_HDF5_Output_State   ,ONLY: WriteState
USE MOD_Mesh_Vars           ,ONLY: MeshFile
USE MOD_TimeDisc_Vars       ,ONLY: t,tAnalyze,tEnd,dt,dt_min,dt_minOld
USE MOD_TimeDisc_Vars       ,ONLY: nCalcTimeStep,nCalcTimeStepMax
USE MOD_TimeDisc_Vars       ,ONLY: doAnalyze,doFinalize
#if USE_PARTICLES
USE MOD_Particle_Output     ,ONLY: FillParticleData
USE MOD_Particle_TimeDisc_Vars,ONLY: UseManualTimeStep
#endif /*USE_PARTICLES*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars    ,ONLY: DoLoadBalance,LoadBalanceSample,PerformLBSample
USE MOD_Particle_Globals    ,ONLY: LESSEQUALTOLERANCE
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: errType
!===================================================================================================================================

! Return if no timestep update requested in this iteration
IF (nCalcTimeStep.GE.1) THEN
  nCalcTimeStep = nCalcTimeStep-1
  RETURN
END IF

#if USE_PARTICLES
IF (UseManualTimeStep) THEN
  dt_min(DT_ANALYZE) = tAnalyze-t             ! Time to next analysis, put in extra variable so number does not change due to numerical errors
  dt_min(DT_END)     = tEnd    -t             ! Do the same for end time
  dt                 = MINVAL(dt_min)

  IF (dt.EQ.dt_min(DT_ANALYZE))       doAnalyze  = .TRUE.
  IF (dt.EQ.dt_min(DT_END    )) THEN; doAnalyze  = .TRUE.; doFinalize = .TRUE.; END IF
  dt                 = MINVAL(dt_min,MASK=dt_min.GT.0)

#if USE_LOADBALANCE
  ! Activate normal load balancing (NOT initial restart load balancing)
  ! 1.) Catch all iterations within sampling interval (make sure to get the first iteration in interval): LESSEQUALTOLERANCE(a,b,tol)
  ! 2.)             Load balancing is activated: DoLoadBalance=T
  IF( (LESSEQUALTOLERANCE(dt_Min(DT_ANALYZE), LoadBalanceSample*dt, 1e-5) &
      .OR. LESSEQUALTOLERANCE(dt_Min(DT_END), LoadBalanceSample*dt, 1e-5))&
      .AND. DoLoadBalance) PerformLBSample=.TRUE. ! Activate load balancing in this time step
#endif /*USE_LOADBALANCE*/

  ! This is strictly not a manual time step but we need to compensate for small errors along the way ...
  ! Increase time step if the NEXT time step would be smaller than dt/100
  IF(dt_min(DT_ANALYZE)-dt.LT.dt/100.0 .AND. dt_min(DT_ANALYZE).GT.0) THEN; dt = dt_min(DT_ANALYZE); doAnalyze  = .TRUE.; END IF
  ! Increase time step if the LAST time step would be smaller than dt/100
  IF(    dt_min(DT_END)-dt.LT.dt/100.0 .AND. dt_min(DT_END    ).GT.0) THEN; dt = dt_min(DT_END)    ; doAnalyze  = .TRUE.; doFinalize = .TRUE.; END IF

  RETURN
END IF
#endif /*USE_PARTICLES*/

dt_min(DT_MIN)     = EvalTimeStep(errType)
dt_min(DT_ANALYZE) = tAnalyze-t             ! Time to next analysis, put in extra variable so number does not change due to numerical errors
dt_min(DT_END)     = tEnd    -t             ! Do the same for end time
dt                 = MINVAL(dt_min)

IF (dt.EQ.dt_min(DT_ANALYZE))       doAnalyze  = .TRUE.
IF (dt.EQ.dt_min(DT_END    )) THEN; doAnalyze  = .TRUE.; doFinalize = .TRUE.; END IF
dt                 = MINVAL(dt_min,MASK=dt_min.GT.0)

nCalcTimeStep = MIN(FLOOR(ABS(LOG10(ABS(dt_minOld/dt-1.)**2.*100.+EPSILON(0.)))),nCalcTimeStepMax) - 1
dt_minOld     = dt
IF (errType.NE.0) THEN
#if USE_PARTICLES
  ! Fill the SFC-ordered particle arrays
  CALL FillParticleData()
#endif /*USE_PARTICLES*/
  CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,FutureTime=tWriteData,isErrorFile=.TRUE.)
  CALL Abort(__STAMP__,&
#if EQNSYSNR == 3
  'Error: (1) density, (2) convective / (3) viscous timestep / muTilde (4) is NaN. Type/time:',errType,t)
#else
  'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)
#endif
END IF

#if USE_LOADBALANCE
! Activate normal load balancing (NOT initial restart load balancing)
! 1.) Catch all iterations within sampling interval (make sure to get the first iteration in interval): LESSEQUALTOLERANCE(a,b,tol)
! 2.)             Load balancing is activated: DoLoadBalance=T
IF( (LESSEQUALTOLERANCE(dt_Min(DT_ANALYZE), LoadBalanceSample*dt, 1e-5) &
    .OR. LESSEQUALTOLERANCE(dt_Min(DT_END), LoadBalanceSample*dt, 1e-5))&
    .AND. DoLoadBalance) PerformLBSample=.TRUE. ! Activate load balancing in this time step
#endif /*USE_LOADBALANCE*/

! Increase time step if the NEXT time step would be smaller than dt/100
IF(dt_min(DT_ANALYZE)-dt.LT.dt/100.0 .AND. dt_min(DT_ANALYZE).GT.0) THEN; dt = dt_min(DT_ANALYZE); doAnalyze  = .TRUE.; END IF
! Increase time step if the LAST time step would be smaller than dt/100
IF(    dt_min(DT_END)-dt.LT.dt/100.0 .AND. dt_min(DT_END    ).GT.0) THEN; dt = dt_min(DT_END)    ; doAnalyze  = .TRUE.; doFinalize = .TRUE.; END IF

END SUBROUTINE UpdateTimeStep


!===================================================================================================================================
!> Update time step at the beginning of each timedisc loop
!===================================================================================================================================
SUBROUTINE AnalyzeTimeStep()
! MODULES
USE MOD_Globals
USE MOD_Analyze             ,ONLY: Analyze
USE MOD_Analyze_Vars        ,ONLY: analyze_dt,WriteData_dt,tWriteData,nWriteData
USE MOD_AnalyzeEquation_Vars,ONLY: doCalcTimeAverage
USE MOD_DG                  ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars             ,ONLY: U
USE MOD_Equation_Vars       ,ONLY: StrVarNames
USE MOD_HDF5_Output         ,ONLY: WriteBaseFlow
USE MOD_HDF5_Output_State   ,ONLY: WriteState
USE MOD_Mesh_Vars           ,ONLY: MeshFile
USE MOD_Output              ,ONLY: Visualize,PrintAnalyze,PrintStatusLine
USE MOD_PruettDamping       ,ONLY: TempFilterTimeDeriv
USE MOD_RecordPoints        ,ONLY: RecordPoints,WriteRP
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc
USE MOD_Sponge_Vars         ,ONLY: CalcPruettDamping
USE MOD_TestCase            ,ONLY: AnalyzeTestCase
USE MOD_TestCase_Vars       ,ONLY: nAnalyzeTestCase
USE MOD_TimeAverage         ,ONLY: CalcTimeAverage
USE MOD_TimeDisc_Vars       ,ONLY: t,dt,dt_min,tAnalyze,tEnd,CalcTimeStart
USE MOD_TimeDisc_Vars       ,ONLY: iter,iter_analyze,maxIter
USE MOD_TimeDisc_Vars       ,ONLY: doAnalyze,doFinalize,writeCounter
#if FV_ENABLED == 1
USE MOD_FV_Switching        ,ONLY: FV_Info
#elif FV_ENABLED == 2
USE MOD_FV_Blending         ,ONLY: FV_Info
#endif
#if PP_LIMITER
USE MOD_PPLimiter           ,ONLY: PPLimiter_Info,PPLimiter
#endif /*PP_LIMITER*/
#if USE_PARTICLES
USE MOD_Particle_Output     ,ONLY: FillParticleData
USE MOD_Particle_TimeDisc_Vars,ONLY: PreviousTime,UseManualTimeStep
#endif /*USE_PARTICLES*/
#if USE_LOADBALANCE
USE MOD_HDF5_output         ,ONLY: RemoveHDF5
USE MOD_LoadBalance         ,ONLY: ComputeElemLoad
USE MOD_LoadBalance_TimeDisc,ONLY: LoadBalance
USE MOD_LoadBalance_Vars    ,ONLY: DoLoadBalance,ElemTime,ElemTimePart,ElemTimeField,DoLoadBalanceBackup,LoadBalanceSampleBackup
USE MOD_LoadBalance_Vars    ,ONLY: LoadBalanceSample,PerformLBSample,PerformLoadBalance,LoadBalanceMaxSteps,nLoadBalanceSteps
USE MOD_LoadBalance_Vars    ,ONLY: DoInitialAutoRestart,ForceInitialLoadBalance
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimeField,RestartTimeBackup!,RestartWallTime
USE MOD_Particle_Globals    ,ONLY: ALMOSTEQUAL
USE MOD_Particle_Localization,ONLY:CountPartsPerElem
USE MOD_Restart_Vars        ,ONLY: RestartTime,RestartFile
USE MOD_TimeDisc_Vars       ,ONLY: maxIter,time_start
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(iter.EQ.maxIter) THEN
  tEnd=t; tAnalyze=t; tWriteData=t
  doAnalyze=.TRUE.; doFinalize=.TRUE.
END IF

! Call DG operator to fill face data, fluxes, gradients for analyze
IF(doAnalyze) THEN
#if USE_PARTICLES
  ! Skip the call, otherwise particles get incremented twice
  IF (.NOT.UseManualTimeStep) THEN
  PreviousTime = t
#endif /*USE_PARTICLES*/
  CALL DGTimeDerivative_weakForm(t)
#if USE_PARTICLES
  PreviousTime = -1
  END IF
#endif /*USE_PARTICLES*/
END IF

! Analysis (possible PerformAnalyze+WriteStateToHDF5 and/or LoadBalance)
#if USE_LOADBALANCE
! For automatic initial restart, check if the number of sampling steps has been achieved and force a load balance step, but skip
! this procedure in the final iteration after which the simulation if finished
!      DoInitialAutoRestart: user-activated load balance restart in first time step (could already be a restart)
! iter.GE.LoadBalanceSample: as soon as the number of time steps for sampling is reached, perform the load balance restart
!                   maxIter: prevent removal of last state file even though no load balance restart was performed
ForceInitialLoadBalance = .FALSE. ! Initialize
IF(DoInitialAutoRestart.AND.(iter.GE.LoadBalanceSample).AND.(iter.NE.maxIter)) ForceInitialLoadBalance=.TRUE.

IF(   ALMOSTEQUAL(dt,dt_Min(DT_ANALYZE)) &
  .OR.ALMOSTEQUAL(dt,dt_Min(DT_END))     &
  .OR.ForceInitialLoadBalance) THEN
  CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB

  ! Check if loadbalancing is enabled with partweight and set PerformLBSample true to calculate elemtimes with partweight
  ! LoadBalanceSample is 0 if partweightLB or IAR_partweighlb are enabled. If only one of them is set Loadbalancesample switches
  ! during time loop
  IF (LoadBalanceSample.EQ.0 .AND. DoLoadBalance) PerformLBSample=.TRUE.
  ! Routine calculates imbalance and if greater than threshold sets PerformLoadBalance=.TRUE.
  CALL ComputeElemLoad()
  ! Force load balance step after elem time has been calculated when doing an initial load balance step at iter=0
  IF (ForceInitialLoadBalance) PerformLoadBalance=.TRUE.
  ! Do not perform a load balance restart when the last timestep is performed
  IF (iter.EQ.maxIter .OR. ALMOSTEQUAL(dt,dt_Min(DT_END))) PerformLoadBalance = .FALSE.
END IF
#endif /*USE_LOADBALANCE*/

! Call your analysis routine for your testcase here.
IF((MOD(iter,INT(nAnalyzeTestCase,KIND=8)).EQ.0).OR.doAnalyze) CALL AnalyzeTestCase(t,doFinalize)
! Evaluate recordpoints
IF(RP_onProc) CALL RecordPoints(PP_nVar,StrVarNames,iter,t,doAnalyze)
! Update Pruett filter base flow
IF(CalcPruettDamping) CALL TempFilterTimeDeriv(U,dt)

! Analyze and output now
#if USE_LOADBALANCE
IF(doAnalyze .OR. PerformLoadBalance) THEN
#else
IF(doAnalyze)THEN
#endif /*USE_LOADBALANCE*/
  CALL PrintAnalyze(dt_Min(DT_MIN))
  CALL TimeDisc_Info(iter_analyze+1)
#if FV_ENABLED
  ! Summation has one more iter step
  CALL FV_Info(iter_analyze+1)
#endif
#if PP_LIMITER
  CALL PPLimiter_Info(iter_analyze+1)
#endif

  writeCounter = writeCounter+1
#if USE_PARTICLES
  ! Fill the SFC-ordered particle arrays
#if !USE_LOADBALANCE
  ! Do not call without load balance and outside nWriteData
  IF ((writeCounter.EQ.nWriteData).OR.doFinalize) &
#endif /*!USE_LOADBALANCE*/
  CALL FillParticleData()
#endif /*USE_PARTICLES*/
  ! Visualize data and write solution
  IF((writeCounter.EQ.nWriteData).OR.doFinalize)THEN
    ! Write various derived data
    IF(doCalcTimeAverage) CALL CalcTimeAverage(.TRUE.,dt,t)
    IF(RP_onProc)         CALL WriteRP(PP_nVar,StrVarNames,t,.TRUE.)
    IF(CalcPruettDamping) CALL WriteBaseFlow(TRIM(MeshFile),t,tWriteData)
    ! Write state file
    ! NOTE: this should be last in the series, so we know all previous data
    ! has been written correctly when the state file is present
    tWriteData = MIN(tAnalyze+WriteData_dt,tEnd)
    CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,FutureTime=tWriteData,isErrorFile=.FALSE.)
    ! Visualize data
    CALL Visualize(t,U)
  END IF

  ! do analysis
  CALL Analyze(t,iter)
  ! Overwrite WriteCounter after Analyze to keep particle analysis synchronous
  IF((writeCounter.EQ.nWriteData).OR.doFinalize) writeCounter=0

  iter_analyze  = 0
  CalcTimeStart = FLEXITIME()
  tAnalyze      = MIN(tAnalyze+analyze_dt,tEnd)

  ! Disable analyze for next time step
  doAnalyze     = .FALSE.

#if USE_LOADBALANCE
  IF((DoLoadBalance.AND.PerformLBSample.AND.(LoadBalanceMaxSteps.GT.nLoadBalanceSteps).OR.(LoadBalanceMaxSteps.EQ.0)).OR.ForceInitialLoadBalance)THEN
    IF(PerformLoadBalance) THEN
      ! DO NOT DELETE THIS: ONLY recalculate the timestep when the mesh is changed!
      !CALL InitTimeStep() ! re-calculate time step after load balance is performed
      RestartTimeBackup = RestartTime! make backup of original restart time
      RestartTime       = t          ! Set restart simulation time to current simulation time because the time is not read from
                                     ! the state file
    END IF

    CALL LoadBalance(OutputTime=t)

    IF(PerformLoadBalance) THEN
      ! RestartWallTime = FLEXITIME()  ! Set restart wall time if a load balance step is performed
      CALL CPU_TIME(time_start)      ! Set the start CPU time to the time after loadbalance init
    END IF
  ELSE
    ElemTime      = 0. ! nullify ElemTime before measuring the time in the next cycle
    ElemTimePart  = 0.
    ElemTimeField = 0.
  END IF

  ! Switch off Initial Auto Restart (initial load balance) after the restart was performed
  IF (DoInitialAutoRestart) THEN
    ! Remove the extra state file written for load balance (only when load balance restart was performed)
    IF(PerformLoadBalance) CALL RemoveHDF5(RestartFile)
    ! Get original settings from backup variables
    DoInitialAutoRestart    = .FALSE.
    ForceInitialLoadBalance = .FALSE.
    DoLoadBalance           = DoLoadBalanceBackup
    LoadBalanceSample       = LoadBalanceSampleBackup
    ! Set to iter_analyze zero so that this first analysis is not counted and the next analysis is the first one,
    ! but only if the initial load balance restart and dt_Analyze did not coincide
    IF(.NOT.ALMOSTEQUALRELATIVE(dt, dt_Min(DT_ANALYZE), 1E-5)) iter_analyze = 0
    ! Set time of the state file that was created before automatic initial restart (to be written in the next state file)
    ! tPreviousAnalyze = RestartTimeBackup
  END IF

  PerformLBSample    = .FALSE. ! Deactivate load balance sampling
  PerformLoadBalance = .FALSE. ! Unset      load balance flag
#endif /*USE_LOADBALANCE*/
END IF

END SUBROUTINE AnalyzeTimeStep


!==================================================================================================================================
!> Evaluates the initial time step for the current update of U
!==================================================================================================================================
FUNCTION EvalInitialTimeStep(errType) RESULT(dt)
! MODULES
USE MOD_Globals
USE MOD_CalcTimeStep        ,ONLY: CalcTimeStep
USE MOD_TimeDisc_Vars       ,ONLY: dt_kill,dt_dynmin,dt_analyzemin,dtElem,nDtLimited
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: dt
INTEGER,INTENT(OUT)          :: errType
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Initially set dt_analyzemin
dt_analyzemin = HUGE(1.)
! Initially set nDtLimited
nDtLimited    = 0

dt            = CalcTimeStep(errType)
dt_analyzemin = MIN(dt_analyzemin,dt)
IF (dt.LT.dt_dynmin) THEN
  CALL PrintWarning("TimeDisc INFO - Timestep dropped below predefined minimum! - LIMITING!")
  dt = dt_dynmin;   dtElem = dt_dynmin
  nDtLimited  = nDtLimited + 1
END IF
IF (dt.LT.dt_kill) &
  CALL Abort(__STAMP__,"TimeDisc ERROR - Initial timestep below critical kill timestep!")
END FUNCTION EvalInitialTimeStep

!==================================================================================================================================
!> Evaluates the time step for the current update of U
!==================================================================================================================================
FUNCTION EvalTimeStep(errType) RESULT(dt_Min)
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars        ,ONLY: tWriteData
USE MOD_CalcTimeStep        ,ONLY: CalcTimeStep
USE MOD_HDF5_Output_State   ,ONLY: WriteState
USE MOD_Mesh_Vars           ,ONLY: MeshFile
USE MOD_TimeDisc_Vars       ,ONLY: dt_kill,dt_dynmin,t,dt_analyzemin,dtElem,nDtLimited
#if USE_PARTICLES
USE MOD_Particle_Output     ,ONLY: FillParticleData
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: dt_Min
INTEGER,INTENT(OUT)          :: errType
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
dt_Min        = CalcTimeStep(errType)
dt_analyzemin = MIN(dt_analyzemin,dt_Min)
IF (dt_Min.LT.dt_dynmin) THEN
  dt_Min = dt_dynmin;   dtElem = dt_dynmin
  nDtLimited  = nDtLimited + 1
END IF
IF (dt_Min.LT.dt_kill) THEN
#if USE_PARTICLES
  CALL FillParticleData()
#endif /*USE_PARTICLES*/
  CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
                  FutureTime=tWriteData,isErrorFile=.TRUE.)
  CALL Abort(__STAMP__,&
    'TimeDisc ERROR - Critical Kill timestep reached! Time: ',RealInfo=t)
END IF
END FUNCTION EvalTimeStep


!===================================================================================================================================
!> Scaling of the CFL number, from paper GASSNER, KOPRIVA, "A comparision of the Gauss and Gauss-Lobatto
!> Discontinuous Galerkin Spectral Element Method for Wave Propagation Problems" .
!> For N=1-10 input CFLscale can now be (nearly) 1. and will be scaled adequately depending on
!> polynomial degree N, NodeType and TimeDisc method.
!===================================================================================================================================
SUBROUTINE FillCFL_DFL(Nin_CFL,Nin_DFL)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY:CFLScale,CFLScale_Readin,CFLScaleAlpha
#if PARABOLIC
USE MOD_TimeDisc_Vars       ,ONLY:DFLScale,DFLScale_Readin,DFLScaleAlpha,RelativeDFL
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_TimeDisc_Vars       ,ONLY:CFLScaleFV
USE MOD_FV_Vars             ,ONLY:FV_w
#endif /*FV*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nin_CFL !< input polynomial degree for advection terms
INTEGER,INTENT(IN) :: Nin_DFL !< input polynomial degree for viscous terms
                              !< for overintegration via Filtering, the gradients and thus the visscous flux remains of order N+1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: alpha
#if !(PARABOLIC)
INTEGER            :: dummy
#endif
!===================================================================================================================================
! CFL in DG depends on the polynomial degree
! Runge-Kutta methods
alpha       = CFLScaleAlpha(MIN(15,Nin_CFL))
CFLScale(0) = CFLScale(0)*alpha
#if FV_ENABLED
CFLScale(1) = CFLScale(1)*CFLScaleFV*MINVAL(FV_w)/2. ! Scale with smallest FV subcell
#endif
IF((Nin_CFL.GT.15).OR.(CFLScale(0).GT.alpha))THEN
  SWRITE(UNIT_stdOut,'(132("!"))')
  SWRITE(UNIT_stdOut,'(A)')'Warning: The chosen CFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_stdOut,'(132("!"))')
END IF
!scale with 2N+1
CFLScale(0) = CFLScale(0)/(2.*Nin_CFL+1.)
SWRITE(UNIT_stdOut,'(A,2ES16.7)') '   CFL (DG/FV):',CFLScale
CFLScale_Readin = CFLScale

#if PARABOLIC
!########################### DFL ########################################
! DFL in DG depends on the polynomial degree
! since DFl is only on real axis, stability numbers are defined for RK3 and then scaled for RK4

alpha       = DFLScaleAlpha(MIN(10,Nin_DFL))*RelativeDFL
DFLScale(0) = DFLScale(0)*alpha
#if FV_ENABLED
DFLScale(1) = DFLScale(1)*DFLScaleAlpha(1)*RelativeDFL*(MINVAL(FV_w)/2.)**2 ! Scale with smallest FV subcell
#endif
IF((Nin_DFL.GT.10).OR.(DFLScale(0).GT.alpha))THEN
  SWRITE(UNIT_stdOut,'(132("!"))')
  SWRITE(UNIT_stdOut,'(A)')'Warning: The chosen DFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_stdOut,'(132("!"))')
END IF
DFLScale(0) = DFLScale(0)/(2.*Nin_DFL+1.)**2
SWRITE(UNIT_stdOut,'(A,2ES16.7)') '   DFL (DG/FV):',DFLScale
DFLScale_Readin = DFLScale
#else
dummy = Nin_DFL ! prevent compile warning
#endif /*PARABOLIC*/
END SUBROUTINE fillCFL_DFL


!==================================================================================================================================
!> Print information on the timestep
!==================================================================================================================================
SUBROUTINE TimeDisc_Info(iter)
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars       ,ONLY: dt_analyzemin,dt_dynmin,nDtLimited
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN) :: iter !< number of iterations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Output Data
IF ((dt_analyzemin.LT.dt_dynmin) .AND. nDtLimited.GT.0) THEN
    CALL PrintWarning("TimeDisc INFO - Timestep dropped below predefined minimum! - LIMITING!")
    SWRITE(UNIT_stdOut,'(A,ES16.7,A)')' TimeDisc   : physical timestep acc. to CFL/DFL:', dt_analyzemin
    SWRITE(UNIT_stdOut,'(A,F8.3,A,I0,A)')' TimeDisc   : limited timestep amount %: ',REAL(nDtLimited)/iter*100,', ',nDtLimited,' timesteps'
END IF

! reset dt_analyzemin
dt_analyzemin = HUGE(1.)
! reset nDtLimited
nDtLimited = 0

END SUBROUTINE TimeDisc_Info


!==================================================================================================================================
!> Finalizes variables necessary for timedisc subroutines
!==================================================================================================================================
SUBROUTINE FinalizeTimeDisc()
! MODULES
USE MOD_TimeDisc_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
TimeDiscInitIsDone = .FALSE.
SDEALLOCATE(b_dt)
SDEALLOCATE(dtElem)
SDEALLOCATE(Ut_tmp)
SDEALLOCATE(S2)
SDEALLOCATE(UPrev)
SDEALLOCATE(RKA)
SDEALLOCATE(RKb)
SDEALLOCATE(RKc)
SDEALLOCATE(RKg1)
SDEALLOCATE(RKg2)
SDEALLOCATE(RKg3)
SDEALLOCATE(RKdelta)
NULLIFY(TimeStep)
END SUBROUTINE FinalizeTimeDisc

END MODULE MOD_TimeDisc_Functions
