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
LOGICAL                             :: DoLoadBalance                              ! DoLoadBalance
LOGICAL                             :: LoadBalanceTimeBased                       ! Flag if loadbalance is performed based on
                                                                                  ! elapsed time
LOGICAL                             :: PerformLoadBalance=.FALSE.                 ! Flag if loadbalance is performed in current iter
INTEGER                             :: LoadBalanceSample                          ! Number of samples for loadbalance
LOGICAL                             :: InitLoadBalanceIsDone                      ! switch for checking

! counter
INTEGER                             :: nLoadBalance                               ! Number of load balances calculations (calls of ComputeElemLoad)
INTEGER                             :: nLoadBalanceSteps                          ! Number of performed load balances steps
INTEGER                             :: LoadBalanceMaxSteps                        ! Number of maximum allowed performed load balances steps
REAL,ALLOCATABLE                    :: LoadDistri(:)                              ! Weighted load distribution of all procs
INTEGER,ALLOCATABLE                 :: PartDistri(:)                              ! Part distribution of all procs
REAL                                :: MaxWeight                                  ! Maximum Weight of proc on domain
REAL                                :: MinWeight                                  ! Minimum Weight of proc on domain
REAL                                :: CurrentImbalance
REAL                                :: NewImbalance                               ! Imbalance after rebalance step

! Variables moved over from mesh_reading/timedisc/restart
LOGICAL                             :: ElemTimeExists
REAL                                :: RestartWallTime                            ! wall time at the beginning of a simulation OR
                                                                                  ! when a restart is performed via Load Balance
LOGICAl                             :: DoInitialAutoRestart= .FALSE.
INTEGER                             :: InitialAutoRestartSample
LOGICAl                             :: IAR_PerformPartWeightLB= .FALSE.

TYPE tData
  INTEGER, ALLOCATABLE :: offsetElemMPI(:)
  INTEGER              :: numOfCalls
  TYPE(tData), POINTER :: nextData => null()
END TYPE tData
TYPE(tData), POINTER :: firstData => null() !linked-list of old offsetElemMPI for WeightDistributionMethod 5 and 6

!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER(KIND=8)                     :: nSkipAnalyze                               ! Skip Analyze-Dt
REAL                                :: ParticleMPIWeight
REAL                                :: DeviationThreshold                         ! threshold for load-balancing
LOGICAL                             :: writePartitionInfo                         ! write partitioninfo file
REAL                                :: WeightSum                                  ! global sum of all weights
REAL                                :: targetWeight                               ! optimal weight for each proc

!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                    :: tCurrent(:)                                ! Time measured over one step
                                                                                  ! measured for all elements on one proc and later
                                                                                  ! weighted onto each element

!-----------------------------------------------------------------------------------------------------------------------------------
! Element Local measurement
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                    :: ElemTime(:)
REAL,ALLOCATABLE                    :: ElemTime_tmp(:)                          ! Additional container for restarting and keeping the old ElemTime values in
                                                                                ! the state.h5 file
REAL,ALLOCATABLE                    :: ElemGlobalTime(:)
INTEGER(KIND=4),ALLOCATABLE         :: nPartsPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nSurfacefluxPerElem(:)
INTEGER(KIND=4),ALLOCATABLE         :: nTracksPerElem(:)


END MODULE MOD_LoadBalance_Vars
