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

!===================================================================================================================================
! Variables needed for load balancing
!===================================================================================================================================
MODULE MOD_LoadBalance_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                             :: DoLoadBalance                              !> Use dynamic load balancing
INTEGER                             :: WeightDistributionMethod                   !> Load distribution method
! INTEGER                             :: LoadBalanceSampleBackup                    !> Loadbalance sample saved until initial autorestart ist finished
! LOGICAL                             :: DoLoadBalanceBackup                        !> Loadbalance flag saved until initial autorestart ist finished
LOGICAL                             :: PerformLoadBalance=.FALSE.                 !> Flag if load balance is performed in current time step iteration
INTEGER                             :: LoadBalanceSample                          !> Number of samples for loadbalance
LOGICAL                             :: PerformLBSample                            !> Flag for enabling time measurement in current
                                                                                  !> Time step (automatically set depending on LB
                                                                                  !> sampling method)
LOGICAL                             :: PerformPartWeightLB                        !> Flag for performing LB with partMPIWeight
                                                                                  !> instead of summed ElemTimes
                                                                                  !> -> nParts*PartWeight written into elemtime array
LOGICAL                             :: InitLoadBalanceIsDone                      !> Switch for checking

! time measurement
REAL,ALLOCATABLE                    :: tCurrent(:)                                !> Time measurement over one step
                                                                                  !> measured elem-independent and later weighted
                                                                                  !> for indices look into piclas.h
#if USE_LOADBALANCE
REAL,ALLOCATABLE                    :: tCurrent_LB_DG(:)                          !> Time measurement over one step
#endif /*USE_LOADBALANCE*/

! counter
INTEGER                             :: nLoadBalance                               !> Number of load balances calculations (calls of ComputeElemLoad)
INTEGER                             :: nLoadBalanceSteps                          !> Number of performed load balances steps
#if USE_LOADBALANCE
INTEGER                             :: LoadBalanceMaxSteps                        !> Number of maximum allowed performed load balances steps
REAL,ALLOCATABLE                    :: LoadDistri(:)                              !> Weighted load distribution of all procs
INTEGER,ALLOCATABLE                 :: PartDistri(:)                              !> Part distribution of all procs
REAL                                :: MaxWeight                                  !> Maximum Weight of proc on domain
REAL                                :: MinWeight                                  !> Minimum Weight of proc on domain
REAL                                :: CurrentImbalance
REAL                                :: NewImbalance                               ! >Imbalance after rebalance step

! Variables moved over from mesh_reading/timedisc/restart
LOGICAL                             :: ElemTimeExists
REAL                                :: RestartWallTime                            ! wall time at the beginning of a simulation OR
                                                                                  ! when a restart is performed via Load Balance
REAL                                :: RestartTimeBackup                          !
LOGICAl                             :: DoInitialAutoRestart= .FALSE.
INTEGER                             :: InitialAutoRestartSample
LOGICAL                             :: ForceInitialLoadBalance  !> Set true when initial load balance steps are completed and force the load balance

TYPE tData
  INTEGER, ALLOCATABLE :: offsetElemMPI(:)
  INTEGER              :: numOfCalls
  TYPE(tData), POINTER :: nextData => null()
END TYPE tData
TYPE(tData), POINTER :: firstData => null() !linked-list of old offsetElemMPI for WeightDistributionMethod 5 and 6
#endif /*USE_LOADBALANCE*/

!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                                :: ParticleMPIWeight
REAL                                :: DeviationThreshold                         ! threshold for load-balancing
LOGICAL                             :: writePartitionInfo                         ! write partitioninfo file
REAL                                :: WeightSum                                  ! global sum of all weights
REAL                                :: targetWeight                               ! optimal weight for each proc

!-----------------------------------------------------------------------------------------------------------------------------------
! Element Local measurement
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                    :: ElemTime(:)
REAL,ALLOCATABLE                    :: ProcTime(:)
#if USE_LOADBALANCE
REAL,ALLOCATABLE                    :: ElemTime_tmp(:)                          ! Additional container for restarting and keeping the old ElemTime values in
                                                                                ! the state.h5 file
REAL                                :: ElemTimePartTot                          ! Total time spent for particle routines (all procs)
REAL                                :: ElemTimeFieldTot                         ! Total time spent for field routines (all procs)
REAL                                :: ElemTimePart                             ! Time spent for particle routines
REAL                                :: ElemTimeField                            ! Time spent for field routines
REAL,ALLOCATABLE                    :: ElemGlobalTime(:)
INTEGER(KIND=4),ALLOCATABLE         :: nPartsPerElem(:)
INTEGER(KIND=4),ALLOCATABLE         :: nTracksPerElem(:)
INTEGER(KIND=4),ALLOCATABLE         :: nSurfacefluxPerElem(:)
INTEGER(KIND=4),ALLOCATABLE         :: nPartsPerBCElem(:)
#endif /*USE_LOADBALANCE*/

END MODULE MOD_LoadBalance_Vars
