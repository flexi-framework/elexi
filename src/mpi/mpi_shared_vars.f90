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
MODULE MOD_MPI_Shared_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
#if USE_MPI_SHARED
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL            :: MPISharedInitIsDone=.FALSE.

! Communication
INTEGER            :: myRank_Shared                   !> Rank of current proc on current node
INTEGER,ALLOCATABLE:: MPIRankGlobal(:)                !> Array of size nProcessors holding the global rank of each proc
INTEGER,ALLOCATABLE:: MPIRankShared(:)                !> Array of size nProcessors holding the shared rank of each proc
INTEGER            :: nProcessors_Shared              !> Number of procs on current node
INTEGER            :: nProcessors_Global              !> Number of total procs
INTEGER            :: MPI_COMM_SHARED                 !> Communicator on current node

! Mesh
!> Counters
INTEGER            :: nElems_Shared                   !> Number of elems on current node
INTEGER            :: nSides_Shared                   !> Number of sides on current node
INTEGER            :: nTotalElems_Shared              !> Number of elems on current node (including halo region)
INTEGER            :: nTotalSides_Shared              !> Number of sides on current node (including halo region)
INTEGER            :: nElemsHalo_Shared               !> Number of elems on current node (only halo region)
INTEGER            :: nSidesHalo_Shared               !> Number of sides on current node (only halo region)
INTEGER            :: OffsetElem_Shared_Root          !> offsetElem of root on current node
INTEGER            :: OffsetSide_Shared_Root          !> offsetSide of root on current node
INTEGER            :: OffsetElemHalo_Shared           !> offsetElem on current node (only halo region)
INTEGER            :: OffsetSideHalo_Shared           !> offsetSide on current node (only halo region)

!> Coordinates and relationships
REAL,POINTER       :: Elem_xGP_Shared(:,:,:,:,:)      !> Gauss points on current node
INTEGER            :: Elem_xGP_Shared_Win             !> Pointer to shared memory window
INTEGER,POINTER    :: ElemToSide_Shared(:,:,:)        !> ElemToSide mapping on current node
INTEGER            :: ElemToSide_Shared_Win           !> Pointer to shared memory window
INTEGER,POINTER    :: SideToElem_Shared(:,:)          !> SideToElem mapping on current node
INTEGER            :: SideToElem_Shared_Win           !> Pointer to shared memory window
REAL,POINTER       :: NodeCoords_Shared(:,:,:,:,:)    !> NodeCoords on current node
INTEGER            :: NodeCoords_Shared_Win           !> Pointer to shared memory window

!> Solution
REAL,POINTER       :: U_Shared(:,:,:,:,:)             !> DG solution on current node
INTEGER            :: U_Shared_Win                    !> Pointer to shared memory window


#endif /* MPI_SHARED */
END MODULE
