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
MODULE MOD_LoadBalance
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersLoadBalance
  MODULE PROCEDURE DefineParametersLoadBalance
END INTERFACE

#if USE_MPI
INTERFACE InitLoadBalance
  MODULE PROCEDURE InitLoadBalance
END INTERFACE

INTERFACE FinalizeLoadBalance
  MODULE PROCEDURE FinalizeLoadBalance
END INTERFACE

#if USE_LOADBALANCE
INTERFACE InitLoadBalanceTracking
  MODULE PROCEDURE InitLoadBalanceTracking
END INTERFACE

INTERFACE ComputeElemLoad
  MODULE PROCEDURE ComputeElemLoad
END INTERFACE

INTERFACE PrintImbalance
  MODULE PROCEDURE PrintImbalance
END INTERFACE
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

PUBLIC :: DefineParametersLoadBalance
#if USE_MPI
PUBLIC :: InitLoadBalance
PUBLIC :: FinalizeLoadBalance
#if USE_LOADBALANCE
PUBLIC :: InitLoadBalanceTracking
PUBLIC :: ComputeElemLoad
PUBLIC :: PrintImbalance
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
CALL prms%SetSection("LoadBalance")
CALL prms%CreateLogicalOption( 'DoLoadBalance'                ,  "Set flag for doing dynamic LoadBalance."                        &
                                                              ,  '.FALSE.')
CALL prms%CreateIntOption(     'LoadBalanceSample'            ,  "Define number of iterations (before analyze_dt)"              //&
                                                                 " that are used for calculation of elemtime information"         &
                                                              ,  '1')
CALL prms%CreateIntOption(    'LoadBalanceMaxSteps'           ,  'Define number of maximum load balacing steps that are allowed.' &
                                                              ,  '0')
CALL prms%CreateLogicalOption('PartWeightLoadBalance'         ,  'Set flag for doing LoadBalance with partMPIWeight instead of '//&
                                                                 'elemtimes. Elemtime array in state file is filled with '      //&
                                                                 'nParts*PartMPIWeight for each Elem. '                         //&
                                                                 ' If Flag [TRUE] LoadBalanceSample is set to 0 and vice versa.'  &
                                                              ,  '.FALSE.')
CALL prms%CreateRealOption(    'Load-DeviationThreshold'      ,  "Define threshold for dynamic load-balancing.\n"                //&
                                                                 "Restart performed if (Maxweight-Targetweight)/Targetweight >"  //&
                                                                 " defined value."                                                 &
                                                              ,  '0.10')
CALL prms%CreateRealOption(    'Part-MPIWeight'               ,  "Define weight of particles for elem loads."                      &
                                                              ,  '0.02')
CALL prms%CreateIntOption(     'WeightDistributionMethod'     ,  "Method for distributing the elem to procs.\n"                  //&
                                                                 "DEFAULT: 1 if Elemtime exits, else -1\n"                       //&
                                                                 "-1   : elements are equally distributed\n"                     //&
                                                                 " 0   : distribute to procs using elemloads\n"                  //&
                                                                 " 1   : distribute to procs using elemloads, last proc receives"//&
                                                                 " least\n"                                                      //&
                                                                 " 3/4 : parallel single step distribution\n"                    //&
                                                                 " 5/6 : iterative smoothing of loads towards last proc\n")

CALL prms%SetSection("Restart")
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

END SUBROUTINE DefineParametersLoadBalance


#if USE_MPI
SUBROUTINE InitLoadBalance()
!===================================================================================================================================
! Init load balancing, new initialization of variables for load balancing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars
USE MOD_ReadInTools            ,ONLY: GETLOGICAL, GETREAL, GETINT
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LOAD BALANCE...'

#if USE_LOADBALANCE
IF(nProcessors.EQ.1)THEN
  ! deactivate load balance for single core computations
  DoLoadBalance        = .FALSE.
  SWRITE(UNIT_stdOut,'(A)') ' | No LoadBalance (nProcessors=1)'
  DeviationThreshold   = HUGE(1.0)
  PerformPartWeightLB  = .FALSE.
ELSE
  DoLoadBalance        = GETLOGICAL('DoLoadBalance')
  LoadBalanceSample    = GETINT    ('LoadBalanceSample')
  LoadBalanceMaxSteps  = GETINT    ('LoadBalanceMaxSteps')
  DeviationThreshold   = GETREAL   ('Load-DeviationThreshold')
  PerformPartWeightLB  = GETLOGICAL('PartWeightLoadBalance')
END IF

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
#else
DoLoadBalance          = .FALSE. ! deactivate loadbalance if no preproc flag is set
DeviationThreshold     = HUGE(1.0)
LoadBalanceSample      = 0
PerformPartWeightLB    = .FALSE.
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
#endif /*USE_LOADBALANCE*/

InitLoadBalanceIsDone  = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
!SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitLoadBalance


#if USE_LOADBALANCE
SUBROUTINE InitLoadBalanceTracking
!===================================================================================================================================
! Re-allocate nPartsPerElem depending on new number of elements
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nTracksPerElem,nSurfacefluxPerElem
USE MOD_Mesh_Vars              ,ONLY: nElems
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
ALLOCATE(nPartsPerElem(1:nElems))
SDEALLOCATE(nTracksPerElem)
ALLOCATE(nTracksPerElem(1:nElems))
SDEALLOCATE(nSurfacefluxPerElem)
ALLOCATE(nSurfacefluxPerElem(1:nElems))

CALL AddToElemData(ElementOut,'nPartsPerElem',IntArray=nPartsPerElem(:))
nPartsPerElem       = 0
nTracksPerElem      = 0
nSurfacefluxPerElem = 0

END SUBROUTINE InitLoadBalanceTracking


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
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nTracksPerElem,nSurfacefluxPerElem
USE MOD_LoadBalance_Vars       ,ONLY: ParticleMPIWeight,ElemTimePartTot,ElemTimePart
USE MOD_LoadBalance_Vars       ,ONLY: CurrentImbalance,PerformLBSample,ElemTimeFieldTot,ElemTimeField
USE MOD_LoadDistribution       ,ONLY: WriteElemTimeStatistics
USE MOD_Particle_Globals
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_TimeDisc_Vars          ,ONLY: t
!USE MOD_TimeDisc_Vars          ,ONLY: nRKStages
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem
INTEGER(KIND=8)       :: HelpSum
REAL                  :: stotalParts,sTotalTracks
REAL                  :: sTotalSurfaceFluxes
REAL                  :: ElemTimeFieldElem,ElemTimePartElem
!===================================================================================================================================
! Initialize
ElemTimeFieldTot = 0.
ElemTimePartTot  = 0.

! If elem times are calculated by time measurement (PerformLBSample) and no Partweight Loadbalance is enabled
IF(PerformLBSample .AND. LoadBalanceSample.GT.0) THEN
  ! number of load balance calls to Compute Elem Load
  nLoadBalance = nLoadBalance + 1

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
    stotalParts   = 1.0/REAL(PP_nElems)
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
  !   stotalBCParts=1.0/REAL(PP_nElems)
  !   nPartsPerBCElem=1
  ! END IF
  ! ----------------------------------------------
! #endif /*PARTICLES*/

  ! distribute times of different routines on elements with respective weightings
  DO iElem = 1,PP_nElems
    ! Add particle LB times to elements with respective weightings
    ! ElemTimeFieldElem = (tCurrent(LB_DG) + tCurrent(LB_DGANALYZE))/REAL(PP_nElems)
    ElemTimeFieldElem = tCurrent(LB_DG)/REAL(PP_nElems)
    ElemTime(iElem) = ElemTime(iElem)  + ElemTimeFieldElem
    ElemTimeField   = ElemTimeField    + ElemTimeFieldElem

    ElemTimePartElem =                                                            &
        + tCurrent(LB_INTERPOLATION)   * nPartsPerElem(iElem)*sTotalParts         &
        + tCurrent(LB_PUSH)            * nPartsPerElem(iElem)*sTotalParts         &
        + tCurrent(LB_TRACK)          * nTracksPerElem(iElem)*sTotalTracks        &
        + tCurrent(LB_SURFFLUX)  * nSurfacefluxPerElem(iElem)*stotalSurfacefluxes

    ElemTime(iElem) = ElemTime(iElem) + ElemTimePartElem
    ElemTimePart    = ElemTimePart    + ElemTimePartElem
  END DO ! iElem=1,PP_nElems

! If no Elem times are calculated but Partweight Loadbalance is enabled
ELSE IF(PerformLBSample .AND. LoadBalanceSample.EQ.0) THEN
  ! number of load balance calls to ComputeElemLoad
  nLoadBalance = nLoadBalance+1
  ! no time measurement and particles are present: simply add the ParticleMPIWeight times the number of particles present
  DO iElem = 1,PP_nElems
    ElemTimePartElem = nPartsPerElem(iElem)*ParticleMPIWeight + 1.0
    ElemTime(iElem)  = ElemTime(iElem) + ElemTimePartElem
    ElemTimePart     = ElemTimePart    + ElemTimePartElem
  END DO ! iElem=1,PP_nElems

  ! Sanity check ElemTime
  IF((MAXVAL(nPartsPerElem).GT.0).AND.(MAXVAL(ElemTime).LE.1.0)) THEN
    IPWRITE (*,*) "parts, time =", MAXVAL(nPartsPerElem),MAXVAL(ElemTime)
    CALL ABORT(__STAMP__&
        ,' ERROR: MAXVAL(nPartsPerElem).GT.0 but MAXVAL(ElemTime).LE.1.0 with ParticleMPIWeight=',RealInfo=ParticleMPIWeight)
  END IF
END IF

! Determine load on complete proc
ProcTime(1:PP_nElems) = SUM(ElemTime(1:PP_nElems))

! Determine sum of balance and calculate target balanced weight and communicate via MPI_ALLREDUCE
CALL ComputeImbalance()

! Fill .csv file for performance analysis and load balance: write data line
CALL WriteElemTimeStatistics(WriteHeader=.FALSE.,time=t)

! only check if imbalance is > a given threshold
PerformLoadBalance = MERGE(.TRUE.,.FALSE.,CurrentImbalance.GT.DeviationThreshold)

! Reset counters
nTracksPerElem       = 0
nPartsPerElem        = 0
!nSurfacePartsPerElem = 0
tCurrent             = 0.

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
USE MOD_Particle_Globals
USE MOD_LoadBalance_Vars    ,ONLY: WeightSum,TargetWeight,CurrentImbalance,MaxWeight,MinWeight
USE MOD_LoadBalance_Vars    ,ONLY: ElemTime,PerformLBSample,PerformPartWeightLB
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimeFieldTot,ElemTimeField
USE MOD_LoadBalance_Vars    ,ONLY: ElemTimePartTot,ElemTimePart
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: WeightSum_loc
!===================================================================================================================================

IF(.NOT.PerformLBSample .AND. .NOT.PerformPartWeightLB) THEN
  WeightSum        = 0.
  TargetWeight     = 0.
  CurrentImbalance = -1.0
ELSE

  ! Collect ElemTime for particles and field separately (only on root process)
  CALL MPI_REDUCE(ElemTimeField , ElemTimeFieldTot , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_WORLD , IERROR)
  CALL MPI_REDUCE(ElemTimePart  , ElemTimePartTot  , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_WORLD , IERROR)
  WeightSum = ElemTimeFieldTot            ! only correct on MPI root
  WeightSum = WeightSum + ElemTimePartTot ! only correct on MPI root
  ! send WeightSum from MPI root to all other procs
  CALL MPI_BCAST(WeightSum,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)

  WeightSum_loc = SUM(ElemTime)

  IF(ALMOSTZERO(WeightSum_loc))THEN
    IPWRITE(*,*) 'Info: The measured time of all elems is zero. ALMOSTZERO(WeightSum)=.TRUE., SUM(ElemTime)=',WeightSum_loc
  END IF

  !CALL MPI_ALLREDUCE(WeightSum_loc,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
  CALL MPI_ALLREDUCE(WeightSum_loc,MaxWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
  CALL MPI_ALLREDUCE(WeightSum_loc,MinWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)

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
USE MOD_LoadBalance_Vars    ,ONLY: TargetWeight,CurrentImbalance,MaxWeight,MinWeight,DeviationThreshold
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.MPIRoot) RETURN

  SWRITE(UNIT_stdOut,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.2,A)')&
      '   MinWeight: ', MinWeight, '   MaxWeight: ', MaxWeight, '   TargetWeight: ', TargetWeight,'    CurrentImbalance: ',&
        CurrentImbalance, '    (Threshold: ', DeviationThreshold, ')'

END SUBROUTINE PrintImbalance


#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
! Deallocate arrays
!===================================================================================================================================
SUBROUTINE FinalizeLoadBalance()
! MODULES
USE MOD_LoadBalance_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Deallocate arrays
SDEALLOCATE(tCurrent)

InitLoadBalanceIsDone = .FALSE.

END SUBROUTINE FinalizeLoadBalance
#endif /*MPI*/

END MODULE MOD_LoadBalance
