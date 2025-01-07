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
#if USE_PARTICLES
#include "particle.h"
#endif /*USE_PARTICLES*/

!===================================================================================================================================
!> Module contains the routines for load balancing
!===================================================================================================================================
MODULE MOD_LoadBalance
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: DefineParametersLoadBalance
#if USE_MPI
PUBLIC:: InitLoadBalance
PUBLIC:: FinalizeLoadBalance
#if USE_LOADBALANCE
#if USE_PARTICLES
PUBLIC:: InitLoadBalanceTracking
#endif /*USE_PARTICLES*/
PUBLIC:: ComputeElemLoad
PUBLIC:: PrintImbalance
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Define parameters
!===================================================================================================================================
SUBROUTINE DefineParametersLoadBalance()
! MODULES
USE MOD_ReadInTools            ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('LoadBalance')
CALL prms%CreateLogicalOption( 'DoLoadBalance'                ,  "Set flag for doing dynamic LoadBalance."                        &
                                                              ,  '.FALSE.')
CALL prms%CreateLogicalOption( 'UseH5IOLoadBalance'           , 'Use HDF5 IO for dynamic load balancing instead of MPI_ALLGATHERV'&
                                                              ,  '.FALSE.')
CALL prms%CreateIntOption(     'LoadBalanceSample'            ,  "Define number of iterations (before analyze_dt)"              //&
                                                                 " that are used for calculation of elemtime information"         &
                                                              ,  '1')
CALL prms%CreateIntOption(     'LoadBalanceMaxSteps'          ,  'Define number of maximum load balacing steps that are allowed.' &
                                                              ,  '0')
CALL prms%CreateRealOption(    'Load-DeviationThreshold'      ,  "Define threshold for dynamic load-balancing.\n"                //&
                                                                 "Restart performed if (Maxweight-Targetweight)/Targetweight >"  //&
                                                                 " defined value."                                                 &
                                                              ,  '0.10')
CALL prms%CreateIntOption(     'WeightDistributionMethod'     ,  "Method for distributing the elem to procs.\n"                  //&
                                                                 "DEFAULT: 1 if Elemtime exits, else -1\n"                       //&
                                                                 "-1   : elements are equally distributed\n"                     //&
                                                                 " 0   : distribute to procs using elemloads\n"                  //&
                                                                 " 1   : distribute to procs using elemloads, last proc receives"//&
                                                                 " least\n"                                                      //&
                                                                 " 3/4 : parallel single step distribution\n"                    //&
                                                                 " 5/6 : iterative smoothing of loads towards last proc\n")

CALL prms%SetSection("LoadBalance Restart")
CALL prms%CreateLogicalOption( 'DoInitialAutoRestart',           "Set Flag for doing automatic initial restart with"             //&
                                                                 " loadbalancing routines after first 'InitialAutoRestartSample'"//&
                                                                 "-number of iterations.\n"                                      //&
                                                                 "Restart is done if Imbalance > 'Load-DeviationThreshold'."       &
                                                              ,  '.FALSE.')
CALL prms%CreateIntOption(     'InitialAutoRestartSample',       "Define number of iterations at simulation start used for"      //&
                                                                 " elemtime sampling before performing automatic initial"        //&
                                                                 " restart.\n"                                                   //&
                                                                 "IF 0 than one iteration is sampled and statefile written has"  //&
                                                                 " zero timeflag.\n"                                             //&
                                                                 " DEFAULT: LoadBalanceSample.")
#if USE_PARTICLES
CALL prms%CreateLogicalOption( 'PartWeightLoadBalance'        ,  'Set flag for doing LoadBalance with partMPIWeight instead of '//&
                                                                 'elemtimes. Elemtime array in state file is filled with '      //&
                                                                 'nParts*PartMPIWeight for each Elem. '                         //&
                                                                 ' If Flag [TRUE] LoadBalanceSample is set to 0 and vice versa.'  &
                                                              ,  '.FALSE.')
CALL prms%CreateRealOption(    'Part-MPIWeight'               ,  "Define weight of particles for elem loads."                      &
                                                              ,  '0.02')
CALL prms%CreateLogicalOption( 'InitialAutoRestart-PartWeightLoadBalance', "Set flag for doing initial auto restart with"        //&
                                                                 " partMPIWeight instead of  ElemTimes. ElemTime array in state" //&
                                                                 " file is filled with nParts*PartMPIWeight for each Elem. "     //&
                                                                 " If Flag [TRUE] InitialAutoRestartSample is set to 0 and vice" //&
                                                                 "versa.", '.FALSE.')
#endif /*USE_PARTICLES*/


END SUBROUTINE DefineParametersLoadBalance


#if USE_MPI
SUBROUTINE InitLoadBalance()
!===================================================================================================================================
! Init load balancing, new initialization of variables for load balancing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars       ,ONLY: InitLoadBalanceIsDone,DoLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: PerformLBSample,LoadBalanceSample
USE MOD_LoadBalance_Vars       ,ONLY: nLoadBalance,nLoadBalanceSteps,DeviationThreshold
USE MOD_ReadInTools            ,ONLY: GETLOGICAL,GETREAL,GETINT
#if USE_LOADBALANCE
USE MOD_Analyze_Vars           ,ONLY: nWriteData
USE MOD_LoadBalance_Vars       ,ONLY: LoadBalanceMaxSteps
USE MOD_LoadBalance_Vars       ,ONLY: tCurrent
USE MOD_LoadBalance_Vars       ,ONLY: MPInElemSend,MPIoffsetElemSend,MPInElemRecv,MPIoffsetElemRecv
USE MOD_LoadBalance_Vars       ,ONLY: ElemInfoRank_Shared,ElemInfoRank_Shared_Win
USE MOD_LoadBalance_Vars       ,ONLY: DoInitialAutoRestart,InitialAutoRestartSample
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems
USE MOD_MPI_Shared
USE MOD_ReadInTools            ,ONLY: PrintOption
#endif /*USE_LOADBALANCE*/
#if USE_PARTICLES
USE MOD_LoadBalance_Vars       ,ONLY: PerformPartWeightLB
USE MOD_LoadBalance_Vars       ,ONLY: ParticleMPIWeight
USE MOD_LoadBalance_Vars       ,ONLY: MPInPartSend,MPIoffsetPartSend,MPInPartRecv,MPIoffsetPartRecv
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=5)  :: tmpStr
#if USE_PARTICLES
LOGICAL           :: InitialAutoRestartPartWeight
#endif /*USE_PARTICLES*/
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LOAD BALANCE...'

#if USE_LOADBALANCE
IF(nProcessors.EQ.1)THEN
  ! deactivate load balance for single core computations
  DoLoadBalance        = .FALSE.
  UseH5IOLoadBalance   = .FALSE.
  SWRITE(UNIT_stdOut,'(A)') ' | No LoadBalance (nProcessors=1)'
  DeviationThreshold   = HUGE(1.0)
#if USE_PARTICLES
  PerformPartWeightLB  = .FALSE.
#endif /*USE_PARTICLES*/
ELSE
  WRITE(tmpStr,'(I5.5)') nWriteData
  DoLoadBalance        = GETLOGICAL('DoLoadBalance')
  UseH5IOLoadBalance   = GETLOGICAL('UseH5IOLoadBalance')
  LoadBalanceSample    = GETINT    ('LoadBalanceSample')
  LoadBalanceMaxSteps  = GETINT    ('LoadBalanceMaxSteps')
  IF (LoadBalanceMaxSteps.LT.0) LoadBalanceMaxSteps = 0
  DeviationThreshold   = GETREAL   ('Load-DeviationThreshold')
#if USE_PARTICLES
  PerformPartWeightLB  = GETLOGICAL('PartWeightLoadBalance')
#endif /*USE_PARTICLES*/
END IF

#if USE_PARTICLES
IF (PerformPartWeightLB) THEN
  ! Read particle MPI weight
  ParticleMPIWeight   = GETREAL('Part-MPIWeight')
  ! deactivate loadbalance sampling of elemtimes if balancing with partweight is enabled
  LoadBalanceSample = 0
  ! CALL PrintOption('PartWeightLoadBalance = T : LoadBalanceSample','INFO',IntOpt=LoadBalanceSample)
ELSE IF (LoadBalanceSample.EQ.0) THEN
  ! loadbalance (elemtimes) is done with ParticleMPIWeight if LoadBalanceSample is set to zero
  PerformPartWeightLB = .TRUE.
  ! CALL PrintOption('LoadbalanceSample = 0 : PartWeightLoadBalance','INFO',LogOpt=PerformPartWeightLB)
END IF
#endif /*USE_PARTICLES*/
#else  /*USE_LOADBALANCE*/
DoLoadBalance          = .FALSE. ! deactivate loadbalance if no preproc flag is set
DeviationThreshold     = HUGE(1.0)
LoadBalanceSample      = 0
#if USE_PARTICLES
PerformPartWeightLB    = .FALSE.
#endif /*USE_PARTICLES*/
#endif /*USE_LOADBALANCE*/
nLoadBalance           = 0
nLoadBalanceSteps      = 0
PerformLBSample        = .FALSE.

#if USE_LOADBALANCE
! Initialize and nullify variables
ALLOCATE(tCurrent(1:LB_NTIMES))
! Allocation length (1:number of loadbalance times)
! look into piclas.h for more info about time names
tcurrent               = 0.
ALLOCATE(MPInElemSend(nProcessors),MPIoffsetElemSend(nProcessors),MPInElemRecv(nProcessors),MPIoffsetElemRecv(nProcessors))
#if USE_PARTICLES
ALLOCATE(MPInPartSend(nProcessors),MPIoffsetPartSend(nProcessors),MPInPartRecv(nProcessors),MPIoffsetPartRecv(nProcessors))
#endif /*USE_PARTICLES*/
CALL Allocate_Shared((/nGlobalElems/),ElemInfoRank_Shared_Win,ElemInfoRank_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemInfoRank_Shared_Win,iError)

! Automatically do a load balance step at the beginning of a new simulation or a user-restarted simulation
DoInitialAutoRestart = GETLOGICAL('DoInitialAutoRestart')
IF(nProcessors.LT.2) DoInitialAutoRestart = .FALSE.

WRITE(UNIT=tmpStr,FMT='(I0)') LoadBalanceSample
InitialAutoRestartSample     = GETINT('InitialAutoRestartSample',TRIM(tmpStr))
#if USE_PARTICLES
InitialAutoRestartPartWeight = GETLOGICAL('InitialAutoRestart-PartWeightLoadBalance','F')
IF (InitialAutoRestartPartWeight) THEN
  InitialAutoRestartSample = 0 ! deactivate loadbalance sampling of ElemTimes if balancing with partweight is enabled
  CALL PrintOption('InitialAutoRestart-PartWeightLoadBalance = T : InitialAutoRestartSample','INFO',IntOpt=InitialAutoRestartSample)
ELSE IF (InitialAutoRestartSample.EQ.0) THEN
  InitialAutoRestartPartWeight = .TRUE. ! loadbalance (ElemTimes) is done with partmpiweight if loadbalancesampling is set to zero
  CALL PrintOption('InitialAutoRestart-PartWeightLoadBalance','INFO',LogOpt=InitialAutoRestartPartWeight)
END IF
#endif /*USE_PARTICLES*/
#endif /*USE_LOADBALANCE*/

InitLoadBalanceIsDone  = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitLoadBalance


#if USE_LOADBALANCE
#if USE_PARTICLES
SUBROUTINE InitLoadBalanceTracking
!===================================================================================================================================
! Re-allocate nPartsPerElem depending on new number of elements
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nTracksPerElem,nSurfacefluxPerElem
#if PARTICLES_COUPLING >= 2
USE MOD_LoadBalance_Vars       ,ONLY: nDeposPerElem
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
USE MOD_LoadBalance_Vars       ,ONLY: nCollsPerElem
#endif /*PARTICLES_COULING*/
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! allocate element counters
! Maybe this could happen in shared memory
SDEALLOCATE(nPartsPerElem)
ALLOCATE(nPartsPerElem(      1:nElems))
SDEALLOCATE(nTracksPerElem)
ALLOCATE(nTracksPerElem(     1:nElems))
SDEALLOCATE(nSurfacefluxPerElem)
ALLOCATE(nSurfacefluxPerElem(1:nElems))
#if PARTICLES_COUPLING >= 2
SDEALLOCATE(nDeposPerElem)
ALLOCATE(nDeposPerElem(      1:nElems))
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
SDEALLOCATE(nCollsPerElem)
ALLOCATE(nCollsPerElem(      1:nElems))
#endif /*PARTICLES_COUPLING == 4*/

CALL AddToElemData(ElementOut,'nPartsPerElem',IntArray=nPartsPerElem(:))
nPartsPerElem       = 0
nTracksPerElem      = 0
nSurfacefluxPerElem = 0
#if PARTICLES_COUPLING >= 2
nDeposPerElem       = 0
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
nCollsPerElem       = 0
#endif /*PARTICLES_COUPLING*/

END SUBROUTINE InitLoadBalanceTracking
#endif /*USE_PARTICLES*/


SUBROUTINE ComputeElemLoad()
!----------------------------------------------------------------------------------------------------------------------------------!
! compute the element load
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime,ProcTime,tCurrent,nLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: DeviationThreshold,PerformLoadBalance,LoadBalanceSample
USE MOD_LoadBalance_Vars       ,ONLY: CurrentImbalance,PerformLBSample
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimeFieldTot,ElemTimeField
USE MOD_LoadDistribution       ,ONLY: WriteElemTimeStatistics
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_TimeDisc_Vars          ,ONLY: t
#if FV_ENABLED
USE MOD_FV_Vars                ,ONLY: FV_Elems
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimeFVTot,ElemTimeFV
#endif /*FV_ENABLED*/
#if USE_PARTICLES
! USE MOD_Particle_Globals
USE MOD_LoadBalance_Vars       ,ONLY: ParticleMPIWeight
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nTracksPerElem,nSurfacefluxPerElem
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimePartTot,ElemTimePart
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
#if PARTICLES_COUPLING >= 2
USE MOD_LoadBalance_Vars       ,ONLY: nDeposPerElem
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimePartDepoTot,ElemTimePartDepo
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
USE MOD_LoadBalance_Vars       ,ONLY: nCollsPerElem
USE MOD_LoadBalance_Vars       ,ONLY: ElemTimePartCollTot,ElemTimePartColl
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/
!USE MOD_TimeDisc_Vars          ,ONLY: nRKStages
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem
REAL                  :: ElemTimeFieldElem
#if FV_ENABLED
REAL                  :: nElemsFV
REAL                  :: ElemTimeFVElem
#endif /*FV_ENABLED*/
#if USE_PARTICLES
INTEGER(KIND=8)       :: helpSum
REAL                  :: ElemTimePartElem
REAL                  :: stotalParts,sTotalTracks
REAL                  :: sTotalSurfaceFluxes
#endif /*USE_PARTICLES*/
!===================================================================================================================================
! Initialize
ElemTimeFieldTot    = 0.
#if FV_ENABLED
ElemTimeFVTot       = 0.
#endif /*FV_ENABLED*/
#if USE_PARTICLES
ElemTimePartTot     = 0.
#if PARTICLES_COUPLING >= 2
ElemTimePartDepoTot = 0.
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
ElemTimePartCollTot = 0.
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/

! If elem times are calculated by time measurement (PerformLBSample) and no Partweight Loadbalance is enabled
IF(PerformLBSample .AND. LoadBalanceSample.GT.0) THEN
  ! number of load balance calls to Compute Elem Load
  nLoadBalance = nLoadBalance + 1

#if FV_ENABLED
  nElemsFV = REAL(COUNT(FV_Elems.EQ.1))
#endif /*FV_ENABLED*/

#if USE_PARTICLES
  ! Calculate weightings, these are the denominators
  sTotalTracks        = 1.
  sTotalSurfaceFluxes = 1.

  ! Recount particles
  CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.)

  ! Calculate and weight particle number per element
  helpSum = SUM(nPartsPerElem)
  IF(helpSum.GT.0) THEN
    stotalParts   = 1.0/REAL(helpSum)
  ! No particle count, distribute load equally on all elements
  ELSE
    stotalParts   = 1.0/REAL(nElems)
    nPartsPerElem = 1
  END IF

  ! Calculate and weight tracks per element
  IF (TrackingMethod.EQ.REFMAPPING) THEN
    helpSum = SUM(nTracksPerElem)
    IF (SUM(nTracksPerElem).GT.0) THEN
      sTotalTracks = 1.0/REAL(helpSum)
    END IF
  END IF
  ! calculate and weight number of surface fluxes
  helpSum = SUM(nSurfacefluxPerElem)
  IF (helpSum.GT.0) THEN
    sTotalSurfaceFluxes = 1.0/REAL(helpSum)
  END IF
  ! Calculate weight number of particles used in sampling of bc element for adaptive particle bcs
  ! helpSum = SUM(nPartsPerBCElem)
  ! IF (helpSum.GT.0) THEN
  !   stotalBCParts=1.0/REAL(helpSum)
  ! ELSE
  !   stotalBCParts=1.0/REAL(nElems)
  !   nPartsPerBCElem=1
  ! END IF
  ! ----------------------------------------------
#endif /*USE_PARTICLES*/

  ! distribute times of different routines on elements with respective weightings
  DO iElem = 1,nElems
    ! Add LB times to elements with respective weightings
    ! ElemTimeFieldElem = (tCurrent(LB_DG) + tCurrent(LB_DGANALYZE))/REAL(nElems)
    ElemTimeFieldElem = tCurrent(LB_DG)/REAL(nElems)
    ElemTime(iElem)   = ElemTime(iElem)  + ElemTimeFieldElem
    ElemTimeField     = ElemTimeField    + ElemTimeFieldElem

#if FV_ENABLED
    ! Add the time for FV
    IF (FV_Elems(iElem).EQ.1) THEN
      ElemTimeFVElem  = tCurrent(LB_FV)/nElemsFV
      ElemTime(iElem) = ElemTime(iElem)  + ElemTimeFVElem
      ElemTimeField   = ElemTimeField    + ElemTimeFVElem
      ElemTimeFV      = ElemTimeFV       + ElemTimeFVElem
    END IF ! FV_Elems
#endif /*FV_ENABLED*/

#if USE_PARTICLES
    ElemTimePartElem =                                                             &
        + tCurrent(LB_INTERPOLATION)  * nPartsPerElem(iElem)       * sTotalParts   &
        + tCurrent(LB_PUSH)           * nPartsPerElem(iElem)       * sTotalParts   &
        + tCurrent(LB_TRACK)          * nTracksPerElem(iElem)      * sTotalTracks  &
        + tCurrent(LB_SURFFLUX)       * nSurfacefluxPerElem(iElem) * stotalSurfacefluxes

    ElemTime(iElem) = ElemTime(iElem) + ElemTimePartElem
    ElemTimePart    = ElemTimePart    + ElemTimePartElem

#if PARTICLES_COUPLING >= 2
    ElemTimePartElem = 0.
    ElemTimePartElem =                                                             &
        + tCurrent(LB_DEPOSITION)     * nDeposPerElem(iElem)       * sTotalParts

    ElemTime(iElem)  = ElemTime(iElem) + ElemTimePartElem
    ElemTimePart     = ElemTimePart    + ElemTimePartElem

    ! Also count the deposition time separately
    ElemTimePartDepo = ElemTimePartDepo + ElemTimePartElem
#endif /*PARTICLES_COUPLING*/

#if PARTICLES_COUPLING == 4
    ElemTimePartElem = 0.
    ElemTimePartElem =                                                             &
        + tCurrent(LB_COLLISION)     * nCollsPerElem(iElem)       * sTotalParts

    ElemTime(iElem)  = ElemTime(iElem) + ElemTimePartElem
    ElemTimePart     = ElemTimePart    + ElemTimePartElem

    ! Also count the collision time separately
    ElemTimePartColl = ElemTimePartColl + ElemTimePartElem
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/
  END DO ! iElem=1,nElems

! If no Elem times are calculated but Partweight Loadbalance is enabled
ELSE IF(PerformLBSample .AND. LoadBalanceSample.EQ.0) THEN
  ! number of load balance calls to ComputeElemLoad
  nLoadBalance = nLoadBalance+1
#if USE_PARTICLES
  ! no time measurement and particles are present: simply add the ParticleMPIWeight times the number of particles present
  DO iElem = 1,nElems
    ElemTimePartElem = nPartsPerElem(iElem) * ParticleMPIWeight + 1.0
    ElemTime(iElem)  = ElemTime(iElem) + ElemTimePartElem
    ElemTimePart     = ElemTimePart    + ElemTimePartElem
  END DO ! iElem=1,nElems

  ! Sanity check ElemTime
  IF((MAXVAL(nPartsPerElem).GT.0).AND.(MAXVAL(ElemTime).LE.1.0)) THEN
    IPWRITE(UNIT_stdOut,'(I0,A,I0,F12.6)') 'parts, time =', MAXVAL(nPartsPerElem),MAXVAL(ElemTime)
    CALL Abort(__STAMP__&
        ,' ERROR: MAXVAL(nPartsPerElem).GT.0 but MAXVAL(ElemTime).LE.1.0 with ParticleMPIWeight=',RealInfo=ParticleMPIWeight)
  END IF
#endif /*USE_PARTICLES*/
END IF

! Determine load on complete proc
ProcTime(1:nElems) = SUM(ElemTime(1:nElems))

! Determine sum of balance and calculate target balanced weight and communicate via MPI_ALLREDUCE
CALL ComputeImbalance()

! Fill .csv file for performance analysis and load balance: write data line
CALL WriteElemTimeStatistics(WriteHeader=.FALSE.,time=t)

! only check if imbalance is > a given threshold
PerformLoadBalance = MERGE(.TRUE.,.FALSE.,CurrentImbalance.GT.DeviationThreshold)

! Reset counters
#if USE_PARTICLES
nPartsPerElem      = 0
nTracksPerElem     = 0
! nSurfacePartsPerElem = 0
#endif /*USE_PARTICLES*/
tCurrent           = 0.

END SUBROUTINE ComputeElemLoad


SUBROUTINE ComputeImbalance()
!----------------------------------------------------------------------------------------------------------------------------------!
! subroutine to compute the imbalance
! Maxweight:        maximum weight of all processes
! Minweight:        minimum weight of all processes
! targetweight:     current target weight
! CurrentImbalance: rel. deviation of current imbalance
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_LoadBalance_Vars    ,ONLY: WeightSum,TargetWeight,CurrentImbalance,MaxWeight,MinWeight
USE MOD_LoadBalance_Vars    ,ONLY: ElemTime,PerformLBSample
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimeFieldTot,ElemTimeField
USE MOD_Utils               ,ONLY: ALMOSTZERO
#if FV_ENABLED
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimeFVTot,ElemTimeFV
#endif /*FV_ENABLED*/
#if USE_PARTICLES
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimePartTot,ElemTimePart
USE MOD_Particle_Tools      ,ONLY: GetOffsetAndGlobalNumberOfParts
USE MOD_Particle_Output_Vars,ONLY: offsetnPart,locnPart
#if PARTICLES_COUPLING >= 2
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimePartDepoTot,ElemTimePartDepo
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimePartCollTot,ElemTimePartColl
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: WeightSum_loc
!===================================================================================================================================

IF(.NOT.PerformLBSample                               ) THEN
  WeightSum        = 0.
  TargetWeight     = 0.
  CurrentImbalance = -1.
ELSE
  ! Collect ElemTime for particles and field separately (only on root process)
  CALL MPI_REDUCE(ElemTimeField,ElemTimeFieldTot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
  WeightSum = ElemTimeFieldTot            ! only correct on MPI root

#if FV_ENABLED
  CALL MPI_REDUCE(ElemTimeFV   ,ElemTimeFVTot   ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
#endif /*FV_ENABLED*/

#if USE_PARTICLES
  ! Collect total number of particles for output
  CALL GetOffsetAndGlobalNumberOfParts('ComputeImbalance',offsetnPart,nGlobalNbrOfParticles,locnPart,.FALSE.)

  ! Collect ElemTime for particles and field separately (only on root process)
  CALL MPI_REDUCE(ElemTimePart ,ElemTimePartTot ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
  WeightSum = WeightSum + ElemTimePartTot ! only correct on MPI root

#if PARTICLES_COUPLING >= 2
  ! Also count the deposition time separately
  CALL MPI_REDUCE(ElemTimePartDepo,ElemTimePartDepoTot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
#endif /*PARTICLES_COUPLING*/

#if PARTICLES_COUPLING == 4
  ! Also count the collision time separately
  CALL MPI_REDUCE(ElemTimePartColl,ElemTimePartCollTot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/
  ! send WeightSum from MPI root to all other procs
  CALL MPI_BCAST(WeightSum,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)

  WeightSum_loc = SUM(ElemTime)

  IF(ALMOSTZERO(WeightSum_loc))THEN
    IPWRITE(UNIT_stdOut,'(I0,A,F12.6)') ' Info: The measured time of all elems is zero. ALMOSTZERO(WeightSum)=.TRUE., SUM(ElemTime)=',WeightSum_loc
  END IF

  !CALL MPI_ALLREDUCE(WeightSum_loc,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
  CALL MPI_ALLREDUCE(WeightSum_loc,MaxWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_FLEXI,iError)
  CALL MPI_ALLREDUCE(WeightSum_loc,MinWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_FLEXI,iError)

  ! Calculate the average value that is supposed to be the optimally distributed weight
  TargetWeight = WeightSum/nProcessors

  IF((TargetWeight.GT.MaxWeight).OR.(TargetWeight.LT.MinWeight))THEN
    SWRITE (UNIT_stdOut,'(A)') " ERROR: after ALLREDUCE, TargetWeight is either smaller than MinWeight or larger than MaxWeight!"
  END IF ! (TargetWeight.GT.MaxWeight).OR.(TargetWeight.LT.MinWeight)

  ! Computation of current imbalance
  IF(ABS(TargetWeight).LE.0.)THEN
    CurrentImbalance = 0.
    SWRITE(UNIT_stdOut,'(A,F14.2,A1)')&
        ' ERROR: after ALLREDUCE, TargetWeight (=WeightSum/nProcessors) cannot be zero! TargetWeight=[',TargetWeight,']'
  ELSE
    CurrentImbalance =  (MaxWeight-TargetWeight)/TargetWeight
  END IF
END IF

END SUBROUTINE ComputeImbalance


!===================================================================================================================================
! Routine to print loadbalance stats
!===================================================================================================================================
SUBROUTINE PrintImbalance()
! MODULES
USE MOD_Globals
USE MOD_LoadBalance_Vars    ,ONLY: DoLoadBalance
USE MOD_LoadBalance_Vars    ,ONLY: TargetWeight,CurrentImbalance,MaxWeight,MinWeight,DeviationThreshold
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.MPIRoot      ) RETURN
IF(.NOT.doLoadBalance) RETURN
IF(MaxWeight.EQ.-1   ) RETURN

WRITE(UNIT_stdOut,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.2,A)')&
   ' MinWeight: ', MinWeight, '    MaxWeight: ', MaxWeight, '    TargetWeight: ', TargetWeight,'    CurrentImbalance: ',&
     CurrentImbalance, '    (Threshold: ', DeviationThreshold, ')'

END SUBROUTINE PrintImbalance


#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
! Deallocate arrays
!===================================================================================================================================
SUBROUTINE FinalizeLoadBalance()
! MODULES
USE MOD_Globals
USE MOD_LoadBalance_Vars
#if USE_LOADBALANCE
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars   ,ONLY: MPI_COMM_SHARED
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Deallocate arrays
SDEALLOCATE(tCurrent)
SDEALLOCATE(ElemTime)
SDEALLOCATE(ProcTime)
SDEALLOCATE(LoadDistri)

SDEALLOCATE(MPInElemSend)
SDEALLOCATE(MPIoffsetElemSend)
SDEALLOCATE(MPInElemRecv)
SDEALLOCATE(MPIoffsetElemRecv)
#if USE_PARTICLES
SDEALLOCATE(MPInPartSend)
SDEALLOCATE(MPIoffsetPartSend)
SDEALLOCATE(MPInPartRecv)
SDEALLOCATE(MPIoffsetPartRecv)
#endif /*USE_PARTICLES*/

#if USE_LOADBALANCE
! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
CALL MPI_WIN_UNLOCK_ALL(ElemInfoRank_Shared_Win,iError)
CALL MPI_WIN_FREE(ElemInfoRank_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

! Then, free the pointers or arrays
NULLIFY(ElemInfoRank_Shared)
#endif /*USE_LOADBALANCE*/

InitLoadBalanceIsDone = .FALSE.

END SUBROUTINE FinalizeLoadBalance
#endif /*MPI*/

END MODULE MOD_LoadBalance
