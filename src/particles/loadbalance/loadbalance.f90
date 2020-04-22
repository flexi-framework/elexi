!=================================================================================================================================
! Copyright (c) 2010-2019  Prof. Claus-Dieter Munz
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

INTERFACE AnalyzeLoadBalance
  MODULE PROCEDURE AnalyzeLoadBalance
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
PUBLIC :: AnalyzeLoadBalance
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
CALL prms%CreateLogicalOption( 'DoLoadBalance'                ,  "Set flag for doing dynamic LoadBalance.", '.FALSE.')
CALL prms%CreateLogicalOption( 'DoLoadBalanceTimeBased'       ,  "Set flag to calculate load (im)balance based on elapsed time"  //&
                                                                 " instead of DG and particle weights", '.TRUE.')
CALL prms%CreateIntOption(     'LoadBalanceSample'            ,  "Define number of iterations (before analyze_dt)"               //&
                                                                 " that are used for calculation of elemtime information",         &
                                                                 value='1')

CALL prms%CreateRealOption(    'Load-DeviationThreshold'      ,  "Define threshold for dynamic load-balancing.\n"                //&
                                                                 "Restart performed if (Maxweight-Targetweight)/Targetweight >"  //&
                                                                 " defined value.",                                                &
                                                                 value='0.10')
CALL prms%CreateRealOption(    'Particles-MPIWeight'          ,  "Define weight of particles for elem loads.",                     &
                                                                 value='0.02')
CALL prms%CreateIntOption(     'WeightDistributionMethod'     ,  "Method for distributing the elem to procs.\n"                  //&
                                                                 "DEFAULT: 1 if Elemtime exits else -1\n"                        //&
                                                                 "-1: elements are equally distributed\n"                        //&
                                                                 " 0: distribute to procs using elemloads\n"                     //&
                                                                 " 1: distribute to procs using elemloads, last proc recieves"   //&
                                                                 " least\n"                                                      //&
                                                                 " 2: TODO DEFINE\n"//&
                                                                 " 3: TODO DEFINE\n"//&
                                                                 " 4/5: iterative smoothing of loads towards last proc\n")

CALL prms%SetSection("Restart")
CALL prms%CreateLogicalOption( 'DoInitialAutoRestart',           "Set Flag for doing automatic initial restart with"             //&
                                                                 " loadbalancing routines after first 'InitialAutoRestartSample'"//&
                                                                 "-number of iterations.\n"                                      //&
                                                                 "Restart is done if Imbalance > 'Load-DeviationThreshold'.",      &
                                                                 '.FALSE.')
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
#if USE_LOADBALANCE
USE MOD_Particle_Tracking_Vars ,ONLY: MeasureTrackTime
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LOAD BALANCE ...'

#if USE_LOADBALANCE
IF(nProcessors.EQ.1)THEN
  ! deactivate load balance for single core computations
  DoLoadBalance        = .FALSE.
  SWRITE(UNIT_stdOut,'(A)') 'No LoadBalance (nProcessors=1): DoLoadBalance=', DoLoadBalance
  DeviationThreshold   = HUGE(1.0)
ELSE
  DoLoadBalance        = GETLOGICAL('DoLoadBalance'          ,'F')
  LoadBalanceTimeBased = GETLOGICAL('DOLoadBalanceTimeBased' ,'T')
  LoadBalanceSample    = GETINT    ('LoadBalanceSample'      ,'1')
  DeviationThreshold   = GETREAL   ('Load-DeviationThreshold','0.10')
  MeasureTrackTime     = GETLOGICAL('MeasureTrackTime'       ,'F')
END IF

IF (DoLoadBalance.AND.(.NOT.LoadBalanceTimeBased)) THEN
  ! Read particle MPI weight
  ParticleMPIWeight   = GETREAL('Particles-MPIWeight','0.02')
  LoadBalanceSample   = 0.
  IF (ParticleMPIWeight.LT.0.0) CALL abort(__STAMP__,' ERROR: Particle weight cannot be negative!')

  ! Rescale ParticleMPIWeight with PP_N^3 to account for ElemTime(DG) = 1
  !ParticleMPIWeight   = ParticleMPIWeight * PP_N**3
END IF
#else
DoLoadBalance          = .FALSE. ! deactivate loadbalance if no preproc flag is set
DeviationThreshold     = HUGE(1.0)
LoadBalanceSample      = 0

#endif /*USE_LOADBALANCE*/

! Initialize and nullify variables
ALLOCATE(tCurrent(1:LB_NTIMES))
tcurrent               = 0.
nLoadBalance           = 0
nLoadBalanceSteps      = 0

InitLoadBalanceIsDone  = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLoadBalance


#if USE_LOADBALANCE
SUBROUTINE InitLoadBalanceTracking
!===================================================================================================================================
! Re-allocate nPartsPerElem depending on new number of elements
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5                ,ONLY: AddToElemData,ElementOut
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nTracksPerElem
USE MOD_Mesh_Vars              ,ONLY: nElems
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Maybe this could happen in shared memory
IF (.NOT.ALLOCATED(nPartsPerElem)) THEN
  ALLOCATE(nPartsPerElem(1:nElems))
ELSE
  SDEALLOCATE(nPartsPerElem)
  ALLOCATE(nPartsPerElem(1:nElems))
END IF
CALL AddToElemData(ElementOut,'nPartsPerElem',IntArray=nPartsPerElem(:))
SDEALLOCATE(nTracksPerElem)
ALLOCATE(nTracksPerElem(1:nElems))
nPartsPerElem          = 0
nTracksPerElem         = 0

END SUBROUTINE InitLoadBalanceTracking


SUBROUTINE ComputeElemLoad()
!----------------------------------------------------------------------------------------------------------------------------------!
! compute the element load
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars       ,ONLY: LoadBalanceTimeBased,ElemTime,tCurrent,nLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: DeviationThreshold,PerformLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nTracksPerElem
USE MOD_LoadBalance_Vars       ,ONLY: ParticleMPIWeight
USE MOD_LoadBalance_Vars       ,ONLY: CurrentImbalance
USE MOD_LoadDistribution       ,ONLY: WriteElemTimeStatistics
USE MOD_Particle_Globals
USE MOD_Particle_Localization  ,ONLY: CountPartsPerElem
USE MOD_Particle_Tracking_vars ,ONLY: DoRefMapping
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
!===================================================================================================================================

! Calculate Load Balance based on elapsed time, linearly distributed to each element
IF (LoadBalanceTimeBased) THEN
  ! Increment number of load balance calls to Compute Elem Load
  nLoadBalance = nLoadBalance + 1

  ! Calculate weightings, these are the denominators
  sTotalTracks = 1.0

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
  IF (DoRefMapping) THEN
    helpSum = SUM(nTracksPerElem)
    IF(SUM(nTracksPerElem).GT.0) THEN
      sTotalTracks=1.0/REAL(helpSum)
    END IF
  END IF

  ! distribute times of different routines on elements with respective weightings
  DO iElem = 1,PP_nElems
    ! Time used in DG routines, assume equal load per element
!    !> Missed the first call to DGTimeDerivative_WeakForm, rescale the time accordingly
!    ElemTime(iElem) = ElemTime(iElem)            * ((nRKStages+1)/nRKStages)           &
!                    + tCurrent(LB_DG)            / REAL(PP_nElems)
    ElemTime(iElem) = ElemTime(iElem) + tCurrent(LB_DG)/ REAL(PP_nElems)
    ! Add particle LB times to elements with respective weightings
    ElemTime(iElem) = ElemTime(iElem)                                                  &
                    + tCurrent(LB_INTERPOLATION) * nPartsPerElem(iElem) *sTotalParts   &
                    + tCurrent(LB_PUSH)          * nPartsPerElem(iElem) *sTotalParts   &
                    + tCurrent(LB_TRACK)         * nTracksPerElem(iElem)*sTotalTracks
  END DO ! iElem=1,PP_nElems

! Calculate Load Balance based on the amount of particles per cell
ELSE
  ! simply add the ParticleMPIWeight times the number of particles present
  DO iElem = 1,PP_nElems
    ElemTime(iElem) = 1. + nPartsPerElem(iElem)*ParticleMPIWeight
  END DO ! iElem=1,PP_nElems

  ! Sanity check ElemTime
  IF((MAXVAL(nPartsPerElem).GT.0).AND.(MAXVAL(ElemTime).LE.1.0))THEN
    IPWRITE (*,*) "parts, time =", MAXVAL(nPartsPerElem),MAXVAL(ElemTime)
    CALL abort(__STAMP__&
        ,' ERROR: MAXVAL(nPartsPerElem).GT.0 but MAXVAL(ElemTime).LE.1.0 with ParticleMPIWeight=',RealInfo=ParticleMPIWeight)
  END IF
END IF ! LoadBalanceTimeBased

! Determine sum of balance and calculate target balanced weight and communicate via MPI_ALLREDUCE
CALL ComputeImbalance()

! Fill .csv file for performance analysis and load balance: write data line
CALL WriteElemTimeStatistics(WriteHeader=.FALSE.,time=t)

! only check if imbalance is > a given threshold
PerformLoadBalance = .FALSE.
IF(CurrentImbalance.GT.DeviationThreshold) PerformLoadBalance=.TRUE.

! Reset counters
nTracksPerElem = 0
nPartsPerElem  = 0
tCurrent       = 0.

END SUBROUTINE ComputeElemLoad
#endif /*USE_LOADBALANCE*/


#if USE_LOADBALANCE
SUBROUTINE ComputeImbalance(output_opt)
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
USE MOD_LoadBalance_Vars,    ONLY:WeightSum,TargetWeight,CurrentImbalance,MaxWeight,MinWeight
USE MOD_LoadBalance_Vars,    ONLY:ElemTime
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL  :: output_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Total load on proc
WeightSum=SUM(ElemTime)

IF(ALMOSTZERO(WeightSum))THEN
    IPWRITE(*,*) 'Info: The calculated time of all elems is zero. ALMOSTZERO(WeightSum)=.TRUE., WeightSum=',WeightSum
END IF

CALL MPI_ALLREDUCE(WeightSum,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
CALL MPI_ALLREDUCE(WeightSum,MaxWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_FLEXI,iError)
CALL MPI_ALLREDUCE(WeightSum,MinWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_FLEXI,iError)

WeightSum    = TargetWeight             ! Set total weight for writing to file
TargetWeight = TargetWeight/nProcessors ! Calculate the average value that is supposed to be the optimally distributed weight

! new computation of current imbalance on this proc
IF(ABS(TargetWeight).EQ.0.) THEN
    CurrentImbalance = 0.
ELSEIF(ABS(TargetWeight).LT.0.) THEN
    SWRITE(UNIT_stdOut,'(A,F14.2,A1)')                                                                                             &
        ' ERROR: after ALLREDUCE, WeightSum/TargetWeight cannot be zero! TargetWeight=[',TargetWeight,']'
    CurrentImbalance = HUGE(1.0)
ELSE
    CurrentImbalance = (MaxWeight-TargetWeight ) / TargetWeight
END IF

! Only output if at analyze step
IF (PRESENT(output_opt)) THEN
!    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MaxWeight:        ', MaxWeight
!    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MinWeight:        ', MinWeight
!    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' TargetWeight:     ', TargetWeight
!    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' CurrentImbalance: ', CurrentImbalance
!    WRITE(UNIT_StdOut,'(132("."))')
    SWRITE(UNIT_stdOut,'(A14,ES16.7,A20,ES16.7,A20,ES16.7,A19,ES11.5)') ' MaxWeight  :      ', MaxWeight,                          &
                                                                        '     MinWeight   : ', MinWeight,                          &
                                                                        '     TargetWeight: ', TargetWeight,                       &
                                                                        ' CurrentImbalance: ', CurrentImbalance
END IF

END SUBROUTINE ComputeImbalance


!===================================================================================================================================
! Routine to print loadbalance stats
!===================================================================================================================================
SUBROUTINE AnalyzeLoadBalance()
! MODULES
USE MOD_Globals
USE MOD_LoadBalance_Vars
USE MOD_Particle_Globals    ,ONLY: ALMOSTZERO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: ElemSum,ElemTimeMin,ElemTimeMax
!===================================================================================================================================
! Calculate particle imbalance
ElemSum = SUM(ElemTime)

! Procs without particles would cause divide by zero
CALL MPI_REDUCE(ElemSum,ElemTimeMin,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(ElemSum,ElemTimeMax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)

IF (MPIROOT) THEN
  IF (ElemTimeMin.GT.0) THEN
!    WRITE(UNIT_StdOut,'(132("."))')
    WRITE(UNIT_StdOut,'(A,F5.2,A)')  ' PART. LOAD INBALANCE (MAX/MIN): [',ElemTimeMax/ElemTimeMin,' ]'
  ELSE
    WRITE(UNIT_StdOut,'(A,F5.2,A)')  ' PART. LOAD INBALANCE (MAX/MIN): [',0.,' ]'
  END IF
END IF

CALL ComputeImbalance(output_opt=.TRUE.)

IF(MPIROOT) THEN
    WRITE(UNIT_StdOut,'(132("."))')
END IF

END SUBROUTINE AnalyzeLoadBalance


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
