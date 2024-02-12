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
!> Contains variables to exchange data using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_MPI_Shared_Vars
! MODULES
#if USE_MPI
USE __MPI__
#endif

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
#if USE_MPI
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL            :: MPISharedInitIsDone=.FALSE.

! Communication
INTEGER            :: ComputeNodeRootRank                      !> Rank of compute-node root in global comm
INTEGER            :: myComputeNodeRank                        !> Rank of current proc on current compute-node
INTEGER            :: myLeaderGroupRank                        !> Rank of compute-node root in compute-node-root comm
INTEGER,ALLOCATABLE:: MPIRankGlobal(:)                         !> Array of size nProcessors holding the global rank of each proc
INTEGER,ALLOCATABLE:: MPIRankShared(:)                         !> Array of size nProcessors holding the shared rank of each proc
INTEGER,ALLOCATABLE:: MPIRankLeader(:)                         !> Array of size nLeaderGroupProcs holding the global rank of each proc
INTEGER            :: nComputeNodeProcessors                   !> Number of procs on current compute-node
INTEGER            :: nComputeNodes                            !> Number of nodes
INTEGER            :: nLeaderGroupProcs                        !> Number of nodes (only valid on CN roots)
#if ! (CORE_SPLIT==0)
! When core-level splitting is used, it is not clear how many cores are on the same physical compute node.
INTEGER            :: NbrOfPhysicalNodes                       !> Number of physical nodes (as opposed to virtual nodes) on which the simulation is executed
#endif /*! (CORE_SPLIT==0)*/
INTEGER            :: nProcessors_Global                       !> Number of total procs
MPI_TYPE_COMM      :: MPI_COMM_SHARED          = MPI_COMM_NULL !> Communicator on current compute-node
MPI_TYPE_COMM      :: MPI_COMM_LEADERS_SHARED  = MPI_COMM_NULL !> Communicator compute-node roots (my_rank_shared=0)
MPI_TYPE_REQUEST,ALLOCATABLE:: MPI_COMM_LEADERS_REQUEST(:)     !> Request handle for non-blocking communication
INTEGER            :: MPI_COMM_LEADERS_REQUEST_SIZE            !> Size of request handle for non-blocking communication
MPI_TYPE_GROUP     :: MPI_GROUP_SHARED         = MPI_GROUP_NULL
MPI_TYPE_GROUP     :: MPI_GROUP_LEADERS_SHARED = MPI_GROUP_NULL

! Mesh
! !> Counters
! INTEGER            :: nComputeNodeElems                      !> Number of elems on current compute-node
! INTEGER            :: nComputeNodeSides                      !> Number of sides on current compute-node
! INTEGER            :: nComputeNodeNodes                      !> Number of nodes on current compute-node
! INTEGER            :: nComputeNodeTrees                      !> Number of trees on current compute-node
! INTEGER            :: offsetComputeNodeElem                  !> elem offset of compute-node root
! INTEGER            :: offsetComputeNodeSide                  !> side offset of compute-node root
! INTEGER            :: offsetComputeNodeNode                  !> node offset of compute-node root
! INTEGER            :: offsetComputeNodeTree                  !> tree offset of compute-node root
#if USE_PARTICLES
INTEGER            :: nComputeNodeTotalElems                   !> Number of elems on current compute-node (including halo region)
INTEGER            :: nComputeNodeTotalSides                   !> Number of sides on current compute-node (including halo region)
INTEGER            :: nComputeNodeTotalNodes                   !> Number of nodes on current compute-node (including halo region)
#endif /*USE_PARTICLES*/

! Offsets for MPI_ALLGATHERV
INTEGER,ALLOCATABLE:: displsElem(:),recvcountElem(:)
INTEGER,ALLOCATABLE:: displsSide(:),recvcountSide(:)
INTEGER,ALLOCATABLE:: displsNode(:),recvcountNode(:)
INTEGER,ALLOCATABLE:: displsTree(:),recvcountTree(:)

#if USE_PARTICLES
! Surface sampling
INTEGER,ALLOCATABLE:: MPIRankSharedLeader(:)                   !> Array of size nLeaderGroupProcs holding the leader rank of each proc
INTEGER,ALLOCATABLE:: MPIRankSurfLeader(:)                     !> Array of size nLeaderGroupProcs holding the surf rank of each proc
MPI_TYPE_COMM      :: MPI_COMM_LEADERS_SURF = MPI_COMM_NULL    !> Communicator compute-node roots on surface communicator (my_rank_shared=0)
INTEGER            :: mySurfRank            = -888             !> rank on MPI_COMM_LEADERS_SURF
INTEGER            :: nSurfLeaders                             !> compute-node leaders on MPI_COMM_LEADERS_SURF
! INTEGER            :: nSurfCommProc                            !> compute-nodes which send or receive sides from us

INTEGER,ALLOCATABLE,DIMENSION(:,:):: nSurfSidesLeader          !> number of surf sides per leader proc
                                                               !> 1 - sides from local leader to other leader
                                                               !> 2 - sides from other leader to local leader
#endif /*USE_PARTICLES*/

!> Solution
REAL,POINTER       :: U_Shared(:,:,:,:,:)                      !> DG solution on current node
MPI_TYPE_WIN       :: U_Shared_Win                             !> Pointer to shared memory window

MPI_TYPE_INFO      :: MPI_INFO_SHARED_LOOSE                    !> MPI_INFO object allowing for re-ordering of same origin atomic RMA operations
! INTEGER            :: MPI_INFO_SHARED_STRICT                   !> MPI_INFO object not allowing for re-ordering of same origin atomic RMA operations

!> Other variables in particle_mesh_vars.f90
#endif /* USE_MPI */
END MODULE MOD_MPI_Shared_Vars
