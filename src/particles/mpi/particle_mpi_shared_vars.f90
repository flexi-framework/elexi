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
!===================================================================================================================================
!> Contains variables to exchange data using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_Particle_MPI_Shared_Vars
! MODULES
#if USE_MPI
USE mpi
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
INTEGER            :: ComputeNodeRootRank                   !> Rank of compute-node root in global comm
INTEGER            :: myComputeNodeRank                     !> Rank of current proc on current compute-node
INTEGER            :: myLeaderGroupRank                     !> Rank of compute-node root in compute-node-root comm
INTEGER,ALLOCATABLE:: MPIRankGlobal(:)                      !> Array of size nProcessors holding the global rank of each proc
INTEGER,ALLOCATABLE:: MPIRankShared(:)                      !> Array of size nProcessors holding the shared rank of each proc
INTEGER,ALLOCATABLE:: MPIRankLeader(:)                      !> Array of size nLeaderGroupProcs holding the global rank of each proc
INTEGER            :: nComputeNodeProcessors                !> Number of procs on current compute-node
INTEGER            :: nLeaderGroupProcs                     !> Number of nodes
INTEGER            :: nProcessors_Global                    !> Number of total procs
INTEGER            :: MPI_COMM_SHARED                       !> Communicator on current compute-node
INTEGER            :: MPI_COMM_LEADERS_SHARED               !> Communicator compute-node roots (my_rank_shared=0)

! Mesh
!> Counters
!INTEGER            :: nNonUniqueGlobalSides                 !> total nb. of non-unique sides of mesh (hexahedral: 6*nElems)
!INTEGER            :: nNonUniqueGlobalNodes                 !> total nb. of non-unique nodes of mesh (hexahedral: 8**NGeo * nElems)
!INTEGER            :: nNonUniqueGlobalTrees                 !> total nb. of trees
!INTEGER            :: nUniqueMasterMortarSides              !> total nb. of master mortar sides in the mesh
!INTEGER            :: nComputeNodeElems                     !> Number of elems on current compute-node
!INTEGER            :: nComputeNodeSides                     !> Number of sides on current compute-node
!INTEGER            :: nComputeNodeNodes                     !> Number of nodes on current compute-node
!INTEGER            :: nComputeNodeTrees                     !> Number of trees on current compute-node
INTEGER            :: nComputeNodeTotalElems                !> Number of elems on current compute-node (including halo region)
INTEGER            :: nComputeNodeTotalSides                !> Number of sides on current compute-node (including halo region)
INTEGER            :: nComputeNodeTotalNodes                !> Number of nodes on current compute-node (including halo region)
!INTEGER            :: offsetComputeNodeElem                 !> elem offset of compute-node root
!INTEGER            :: offsetComputeNodeSide                 !> side offset of compute-node root
!INTEGER            :: offsetComputeNodeNode                 !> node offset of compute-node root
!INTEGER            :: offsetComputeNodeTree                 !> tree offset of compute-node root

! Offsets for MPI_ALLGATHERV
INTEGER,ALLOCATABLE:: displsElem(:),recvcountElem(:)
INTEGER,ALLOCATABLE:: displsSide(:),recvcountSide(:)
INTEGER,ALLOCATABLE:: displsNode(:),recvcountNode(:)
INTEGER,ALLOCATABLE:: displsTree(:),recvcountTree(:)

! Surface sampling
INTEGER,ALLOCATABLE:: MPIRankSharedLeader(:)                !> Array of size nLeaderGroupProcs holding the leader rank of each proc
INTEGER,ALLOCATABLE:: MPIRankSurfLeader(:)                  !> Array of size nLeaderGroupProcs holding the surf rank of each proc
INTEGER            :: MPI_COMM_LEADERS_SURF=MPI_COMM_NULL   !> Communicator compute-node roots on surface communicator (my_rank_shared=0)
INTEGER            :: mySurfRank           =-888            !> rank on MPI_COMM_LEADERS_SURF
INTEGER            :: nSurfLeaders                          !> compute-node leaders on MPI_COMM_LEADERS_SURF
!INTEGER            :: nSurfCommProc                         !> compute-nodes which send or receive sides from us

INTEGER,ALLOCATABLE,DIMENSION(:,:):: nSurfSidesLeader       !> number of surf sides per leader proc
                                                            !> 1 - sides from local leader to other leader
                                                            !> 2 - sides from other leader to local leader

INTEGER, ALLOCATABLE :: CNTotalElem2GlobalElem(:)           !> Compute Nodes mapping 1:nTotal -> 1:nGlobal
INTEGER, ALLOCATABLE :: GlobalElem2CNTotalElem(:)           !> Reverse Mapping

!> Solution
REAL,POINTER       :: U_Shared(:,:,:,:,:)             !> DG solution on current node
INTEGER            :: U_Shared_Win                    !> Pointer to shared memory window

!> Other variables in particle_mesh_vars.f90
#endif /* USE_MPI */
END MODULE
