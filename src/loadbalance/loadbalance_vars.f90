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

!===================================================================================================================================
! Variables needed for load balancing
!===================================================================================================================================
MODULE MOD_LoadBalance_Vars
! MODULES
! USE MOD_Globals
USE ISO_C_BINDING ! < Directly include to avoid circular dependency
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: HP => INT16,  & ! half precision (only defined for integer)
                                         SP => REAL32, & ! single precision
                                         DP => REAL64, & ! double precision
                                         QP => REAL128   ! quadruple precision
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                             :: DoLoadBalance                              !> Use dynamic load balancing
LOGICAL                             :: UseH5IOLoadBalance                         !> Use hdf5 IO for dynamic load balancing instead of MPI_ALLGATHERV
INTEGER                             :: LoadBalanceSampleBackup                    !> Loadbalance sample saved until initial autorestart ist finished
LOGICAL                             :: DoLoadBalanceBackup                        !> Loadbalance flag saved until initial autorestart ist finished
LOGICAL                             :: PerformLoadBalance=.FALSE.                 !> Flag if load balance is performed in current time step iteration
INTEGER                             :: LoadBalanceSample                          !> Number of samples for loadbalance
LOGICAL                             :: PerformLBSample                            !> Flag for enabling time measurement in current
                                                                                  !> Time step (automatically set depending on LB
                                                                                  !> sampling method)
LOGICAL                             :: PerformPartWeightLB                        !> Flag for performing LB with partMPIWeight
                                                                                  !> instead of summed ElemTimes
                                                                                  !> -> nParts*PartWeight written into elemtime array
LOGICAL                             :: InitLoadBalanceIsDone                      !> Switch for checking
INTEGER                             :: WeightDistributionMethod   !> Method used for distributing the elements among the available processors

! time measurement
REAL,ALLOCATABLE                    :: tCurrent(:)                                !> Time measurement over one step
                                                                                  !> measured elem-independent and later weighted
                                                                                  !> for indices look into piclas.h
REAL,ALLOCATABLE                    :: tCurrent_LB_DG(:)                          !> Time measurement over one step

! counter
INTEGER                             :: nLoadBalance                               !> Number of load balances calculations (calls of ComputeElemLoad)
INTEGER                             :: nLoadBalanceSteps                          !> Number of performed load balances steps
INTEGER                             :: LoadBalanceMaxSteps                        !> Number of maximum allowed performed load balances steps
REAL,ALLOCATABLE                    :: LoadDistri(:)                              !> Weighted load distribution of all procs
REAL                                :: MaxWeight                                  !> Maximum Weight of proc on domain
REAL                                :: MinWeight                                  !> Minimum Weight of proc on domain
REAL                                :: CurrentImbalance
REAL                                :: NewImbalance                               ! >Imbalance after rebalance step

! Variables moved over from mesh_reading/timedisc/restart
LOGICAL                             :: ElemTimeExists
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

!-----------------------------------------------------------------------------------------------------------------------------------
! element load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                             :: nElemsOld
INTEGER                             :: offsetElemOld
INTEGER,ALLOCATABLE                 :: offsetElemMPIOld(:)
INTEGER,ALLOCATABLE                 :: MPInElemSend(:)
INTEGER,ALLOCATABLE                 :: MPIoffsetElemSend(:)
INTEGER,ALLOCATABLE                 :: MPInElemRecv(:)
INTEGER,ALLOCATABLE                 :: MPIoffsetElemRecv(:)
INTEGER,ALLOCATABLE                 :: MPInPartSend(:)
INTEGER,ALLOCATABLE                 :: MPIoffsetPartSend(:)
INTEGER,ALLOCATABLE                 :: MPInPartRecv(:)
INTEGER,ALLOCATABLE                 :: MPIoffsetPartRecv(:)
INTEGER,POINTER                     :: ElemInfoRank_Shared(:) => NULL()
INTEGER                             :: ElemInfoRank_Shared_Win

!-----------------------------------------------------------------------------------------------------------------------------------
! general load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                                :: DeviationThreshold                         ! threshold for load-balancing
LOGICAL                             :: writePartitionInfo                         ! write partitioninfo file
REAL                                :: WeightSum                                  ! global sum of all weights
REAL                                :: targetWeight                               ! optimal weight for each proc
#if USE_PARTICLES
REAL                                :: ParticleMPIWeight
#endif /*USE_PARTICLES*/

!-----------------------------------------------------------------------------------------------------------------------------------
! Element Local measurement
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                    :: ElemTime(:)
REAL,ALLOCATABLE                    :: ProcTime(:)
REAL,ALLOCATABLE                    :: ElemTime_tmp(:)                          ! Additional container for restarting and keeping the old ElemTime values in
                                                                                ! the state.h5 file
REAL                                :: ElemTimeFieldTot                         ! Total time spent for field routines (all procs)
REAL                                :: ElemTimeField                            ! Time spent for field routines
REAL,ALLOCATABLE                    :: ElemGlobalTime(:)
#if FV_ENABLED
REAL                                :: ElemTimeFVTot                            ! Total time spent for FV routines (all procs)
REAL                                :: ElemTimeFV                               ! Time spent for FV routines
#endif /*FV_ENABLED*/
#if USE_PARTICLES
REAL                                :: ElemTimePartTot                          ! Total time spent for particle routines (all procs)
REAL                                :: ElemTimePart                             ! Time spent for particle routines
INTEGER(KIND=SP),ALLOCATABLE        :: nPartsPerElem(:)
INTEGER(KIND=SP),ALLOCATABLE        :: nTracksPerElem(:)
INTEGER(KIND=SP),ALLOCATABLE        :: nSurfacefluxPerElem(:)
INTEGER(KIND=SP),ALLOCATABLE        :: nPartsPerBCElem(:)
#if PARTICLES_COUPLING >= 2
REAL                                :: ElemTimePartDepoTot                      ! Total time spent for particle deposition (all procs)
REAL                                :: ElemTimePartDepo                         ! Time spent for particle deposition
INTEGER(KIND=SP),ALLOCATABLE        :: nDeposPerElem(:)
#endif /*PARTICLES_COUPLING*/
#if PARTICLES_COUPLING == 4
REAL                                :: ElemTimePartCollTot                      ! Total time spent for particle collision (all procs)
REAL                                :: ElemTimePartColl                         ! Time spent for particle collision
INTEGER(KIND=SP),ALLOCATABLE        :: nCollsPerElem(:)
#endif /*PARTICLES_COUPLING*/
#endif /*USE_PARTICLES*/

END MODULE MOD_LoadBalance_Vars
