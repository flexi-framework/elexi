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
! Builds the mesh for particle tracking, separate from the DG mesh
!===================================================================================================================================
MODULE MOD_Particle_Mesh_Readin
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables
INTEGER,ALLOCATABLE       :: ElemInfo_Shared_tmp(:)
INTEGER,ALLOCATABLE       :: SideInfo_Shared_tmp(:)

!-----------------------------------------------------------------------------------------------------------------------------------
! Interfaces

INTERFACE ReadMeshBasics
  MODULE PROCEDURE ReadMeshBasics
END INTERFACE

INTERFACE ReadMeshElems
  MODULE PROCEDURE ReadMeshElems
END INTERFACE

INTERFACE ReadMeshSides
  MODULE PROCEDURE ReadMeshSides
END INTERFACE

INTERFACE ReadMeshSideNeighbors
  MODULE PROCEDURE ReadMeshSideNeighbors
END INTERFACE

INTERFACE ReadMeshNodes
  MODULE PROCEDURE ReadMeshNodes
END INTERFACE

INTERFACE ReadMeshTrees
  MODULE PROCEDURE ReadMeshTrees
END INTERFACE

INTERFACE CommunicateMeshReadin
  MODULE PROCEDURE CommunicateMeshReadin
END INTERFACE

PUBLIC :: ReadMeshBasics
PUBLIC :: ReadMeshElems
PUBLIC :: ReadMeshSides
PUBLIC :: ReadMeshSideNeighbors
PUBLIC :: ReadMeshNodes
PUBLIC :: ReadMeshTrees
PUBLIC :: CommunicateMeshReadin
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadMeshBasics()
!===================================================================================================================================
! Read basic global counters from mesh file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input                ,ONLY: File_ID,ReadAttribute
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: nNonUniqueGlobalSides,nNonUniqueGlobalNodes
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

!CALL ReadAttribute(File_ID,'nUniqueSides',1,IntScalar=nGlobalUniqueSidesFromMesh)
CALL ReadAttribute(File_ID,'nSides'      ,1,IntScalar=nNonUniqueGlobalSides)
CALL ReadAttribute(File_ID,'nNodes'      ,1,IntScalar=nNonUniqueGlobalNodes)

END SUBROUTINE ReadMeshBasics


SUBROUTINE ReadMeshElems()
!===================================================================================================================================
! Create particle mesh arrays for elems
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Vars                  ,ONLY: offsetElemMPI
USE MOD_Particle_MPI_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iProc
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
!===================================================================================================================================
#if USE_MPI
! allocate shared array for ElemInfo
CALL MPI_ALLREDUCE(nElems,nComputeNodeElems,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
MPISharedSize = INT((ELEM_HALOFLAG)*nGlobalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/ELEMINFOSIZE,nGlobalElems/),ElemInfo_Shared_Win,ElemInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemInfo_Shared_Win,IERROR)

ElemInfo_Shared(1:ELEMINFOSIZE_H5,offsetElem+1:offsetElem+nElems) = ElemInfo(:,:)
ElemInfo_Shared(ELEM_RANK        ,offsetElem+1:offsetElem+nElems) = myRank
CALL MPI_WIN_SYNC(ElemInfo_Shared_Win,IERROR)

! allocate temporary array to hold processor rank for each elem
ALLOCATE(ElemInfo_Shared_tmp(offsetElem+1:offsetElem+nElems))
ElemInfo_Shared_tmp(offsetElem+1:offsetElem+nElems) = myRank

! broadcast elem offset of compute-node root
offsetComputeNodeElem=offsetElem
CALL MPI_BCAST(offsetComputeNodeElem,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

#else
! allocate local array for ElemInfo
ALLOCATE(ElemInfo_Shared(1:ELEMINFOSIZE,1:nElems))
ElemInfo_Shared(1:ELEMINFOSIZE_H5,1:nElems) = ElemInfo(:,:)
#endif  /*USE_MPI*/

! create global element index to global processor index mapping
#if USE_MPI
MPISharedSize = INT(nGlobalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nGlobalElems/),ElemToProcID_Shared_Win,ElemToProcID_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToProcID_Shared_Win,IERROR)

! build mapping from elems to procs. This contains all procs, not just on the compute-node
IF (myComputeNodeRank.EQ.0) THEN
  DO iProc = 1, nProcessors
    ElemToProcID_Shared(offsetElemMPI(iProc-1)+1:offsetElemMPI(iProc)) = iProc - 1
  END DO ! iProc = 1, nProcessors
END IF
CALL MPI_WIN_SYNC(ElemToProcID_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif  /*USE_MPI*/

END SUBROUTINE ReadMeshElems


SUBROUTINE ReadMeshSides()
!===================================================================================================================================
! Create particle mesh arrays for sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,FirstElemInd,LastElemInd,NbElemID
INTEGER                        :: iSide,sideCount
INTEGER                        :: nSideIDs,offsetSideID
INTEGER                        :: iLocSide,jLocSide,nlocSides,nlocSidesNb,NbSideID
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
!===================================================================================================================================
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo(ELEM_LASTSIDEIND ,LastElemInd)-ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd)

#if USE_MPI
! all procs on my compute-node communicate the number of non-unique sides
CALL MPI_ALLREDUCE(nSideIDs,nComputeNodeSides,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
MPISharedSize = INT((SIDEINFOSIZE+1)*nNonUniqueGlobalSides,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/SIDEINFOSIZE+1,nNonUniqueGlobalSides/),SideInfo_Shared_Win,SideInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideInfo_Shared_Win,IERROR)
SideInfo_Shared(1                :SIDEINFOSIZE_H5,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE_H5+1:SIDEINFOSIZE+1 ,offsetSideID+1:offsetSideID+nSideIDs) = 0
CALL MPI_WIN_SYNC(SideInfo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

DO iSide=1,nSideIDs
  END DO

! fill the SIDE_LOCALID. Basically, this array contains the 1:6 local sides of an element. ! If an element has hanging nodes (i.e.
! has a big mortar side), the big side has negative index (-1,-2 or -3) and the next 2 (-2, -3) or 4 (-1) sides are the subsides.
! Consequently, a hexahedral element can have more than 6 non-unique sides. If we find a small mortar side in the small element,
! the SIDE_LOCALID points to the global ID of the big mortar side, indicated by negative sign
DO iElem=FirstElemInd,LastElemInd
  iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  SideInfo_Shared(SIDE_ELEMID,iSide+1:ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)) = iElem
  sideCount = 0
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  DO iLocSide = 1,nlocSides
    iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
    ! Big mortar side
    IF (MOD(SideInfo_Shared(SIDE_TYPE,iSide),10).GT.4) THEN
      ! Check all sides on the small element side to find the small mortar side pointing back
      NbElemID    = SideInfo_Shared(SIDE_NBELEMID,iSide)
      nlocSidesNb = ElemInfo_Shared(ELEM_LASTSIDEIND,NbElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID)
      DO jLocSide = 1,nlocSidesNb
        NbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID) + jLocSide
        IF (ABS(SideInfo_Shared(SIDE_ID,iSide)).EQ.ABS(SideInfo_Shared(SIDE_ID,NbSideID))) THEN
          SideInfo_Shared(SIDE_LOCALID,iSide) = -NbSideID
          EXIT
        END IF
      END DO
    ! Regular side
    ELSE
      sideCount = sideCount + 1
      SideInfo_Shared(SIDE_LOCALID,iSide) = sideCount
    END IF
  END DO
END DO
CALL MPI_WIN_SYNC(SideInfo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

#else
ALLOCATE(SideInfo_Shared(1:SIDEINFOSIZE,1:nSideIDs))
SideInfo_Shared(1:SIDEINFOSIZE,1:nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE+1,1:nSideIDs) = 0
#endif  /*USE_MPI*/
ALLOCATE(SideInfo_Shared_tmp(offsetSideID+1:offsetSideID+nSideIDs))

END SUBROUTINE ReadMeshSides


SUBROUTINE ReadMeshSideNeighbors(ElemID,SideID)
!===================================================================================================================================
! Fills temporary array to add side neighbors to SideInfo(_Shared)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: ElemID
INTEGER,INTENT(IN)             :: SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (ElemID.LE.offsetComputeNodeElem+1 .OR. ElemID.GT.offsetComputeNodeElem+nComputeNodeElems) THEN
  ! neighbour element is outside of compute-node
  SideInfo_Shared_tmp(SideID) = 2
ELSE
  SideInfo_Shared_tmp(SideID) = 1
END IF

END SUBROUTINE ReadMeshSideNeighbors


SUBROUTINE ReadMeshNodes()
!===================================================================================================================================
! Create particle mesh arrays for nodes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input                ,ONLY: ReadArray
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: FirstNodeInd,LastNodeInd
INTEGER                        :: nNodeIDs,offsetNodeID
INTEGER,ALLOCATABLE            :: NodeInfo(:)
REAL,ALLOCATABLE               :: NodeCoords_indx(:,:)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
!===================================================================================================================================
! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)
FirstNodeInd = offsetNodeID+1
LastNodeInd  = offsetNodeID+nNodeIDs

!read local Node Info from data file
ALLOCATE(NodeInfo(FirstNodeInd:LastNodeInd))
CALL ReadArray('GlobalNodeIDs',1,(/nNodeIDs/),offsetNodeID,1,IntArray=NodeInfo)
ALLOCATE(NodeCoords_indx(3,nNodeIDs))
CALL ReadArray('NodeCoords',2,(/3,nNodeIDs/),offsetNodeID,2,RealArray=NodeCoords_indx)

#if USE_MPI
! allocate shared array for NodeInfo
CALL MPI_ALLREDUCE(nNodeIDs,nComputeNodeNodes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
MPISharedSize = INT(nNonUniqueGlobalNodes,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalNodes/),NodeInfo_Shared_Win,NodeInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeInfo_Shared_Win,IERROR)
NodeInfo_Shared(offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeInfo(:)
CALL MPI_WIN_SYNC(NodeInfo_Shared_Win,IERROR)

MPISharedSize = INT(3*nNonUniqueGlobalNodes,MPI_ADDRESS_KIND)*MPI_DOUBLE
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalNodes/),NodeCoords_Shared_Win,NodeCoords_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
NodeCoords_Shared(:,offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeCoords_indx(:,:)
CALL MPI_WIN_SYNC(NodeCoords_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#else
ALLOCATE(NodeInfo_Shared(1:nNodeIDs))
NodeInfo_Shared(1:nNodeIDs) = NodeInfo(:)
ALLOCATE(NodeCoords_Shared(3,nNodeIDs))
NodeCoords_Shared(:,:) = NodeCoords_indx(:,:)
#endif  /*USE_MPI*/

DEALLOCATE(NodeInfo,NodeCoords_indx)

END SUBROUTINE ReadMeshNodes


SUBROUTINE ReadMeshTrees()
!===================================================================================================================================
! Create particle mesh arrays for trees
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
!===================================================================================================================================
#if USE_MPI
! allocate shared array for TreeCoords
CALL MPI_ALLREDUCE(nTrees,nNonUniqueGlobalTrees,1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,IERROR)
MPISharedSize = INT((NGeoTree+1)**3*nNonUniqueGlobalTrees,MPI_ADDRESS_KIND)*MPI_DOUBLE
CALL Allocate_Shared(MPISharedSize,(/3,nGeoTree+1,nGeoTree+1,nGeoTree+1,nNonUniqueGlobalTrees/),TreeCoords_Shared_Win,TreeCoords_Shared)
CALL MPI_WIN_LOCK_ALL(0,TreeCoords_Shared_Win,IERROR)
TreeCoords_Shared(:,:,:,:,offsetTree:offsetTree+nTrees) = TreeCoords(:,:,:,:,:)
CALL MPI_WIN_SYNC(TreeCoords_Shared_Win,IERROR)
#endif  /*USE_MPI*/

END SUBROUTINE ReadMeshTrees


SUBROUTINE CommunicateMeshReadin()
!===================================================================================================================================
! Create particle mesh arrays for nodes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iProc
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: nSideIDs,offsetSideID
INTEGER                        :: nNodeIDs,offsetNodeID
INTEGER,ALLOCATABLE            :: displsCN(:),recvcountCN(:)
INTEGER,ALLOCATABLE            :: displsElem(:),recvcountElem(:)
INTEGER,ALLOCATABLE            :: displsSide(:),recvcountSide(:)
INTEGER,ALLOCATABLE            :: displsNode(:),recvcountNode(:)
INTEGER,ALLOCATABLE            :: displsTree(:),recvcountTree(:)
!===================================================================================================================================
#if USE_MPI
SWRITE(UNIT_stdOut,'(A)') ' COMMUNICATING MESH ON SHARED MEMORY...'

CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd)
offsetNodeID = ElemInfo(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo(ELEM_LASTNODEIND,LastElemInd)-ElemInfo(ELEM_FIRSTNODEIND,FirstElemind)

IF(myComputeNodeRank.EQ.0)THEN
  ! Arrays for the compute-node to communicate their offsets
  ALLOCATE(displsCN   (0:nLeaderGroupProcs-1))
  ALLOCATE(recvcountCN(0:nLeaderGroupProcs-1))
  DO iProc=0,nLeaderGroupProcs-1
    displsCN(iProc) = iProc
  END DO
  recvcountCN(:) = 1
  ! Arrays for the compute node to hold the elem offsets
  ALLOCATE(displsElem   (0:nLeaderGroupProcs-1))
  ALLOCATE(recvcountElem(0:nLeaderGroupProcs-1))
  displsElem(myLeaderGroupRank) = offsetComputeNodeElem
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsElem,recvcountCN,displsCN,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountElem(iProc-1) = displsElem(iProc)-displsElem(iProc-1)
  END DO
  recvcountElem(nLeaderGroupProcs-1) = nGlobalElems - displsElem(nLeaderGroupProcs-1)
END IF

! Broadcast compute node side offset on node
offsetComputeNodeSide=offsetSideID
CALL MPI_BCAST(offsetComputeNodeSide,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF(myComputeNodeRank.EQ.0)THEN
  ! Arrays for the compute node to hold the side offsets
  ALLOCATE(displsSide   (0:nLeaderGroupProcs-1))
  ALLOCATE(recvcountSide(0:nLeaderGroupProcs-1))
  displsSide(myLeaderGroupRank) = offsetComputeNodeSide
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsSide,recvcountCN,displsCN,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountSide(iProc-1) = displsSide(iProc)-displsSide(iProc-1)
  END DO
  recvcountSide(nLeaderGroupProcs-1) = nNonUniqueGlobalSides - displsSide(nLeaderGroupProcs-1)
END IF

! Broadcast compute node node offset on node
offsetComputeNodeNode=offsetNodeID
CALL MPI_BCAST(offsetComputeNodeNode,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF(myComputeNodeRank.EQ.0)THEN
  ! Arrays for the compute node to hold the node offsets
  ALLOCATE(displsNode(0:nLeaderGroupProcs-1))
  ALLOCATE(recvcountNode(0:nLeaderGroupProcs-1))
  displsNode(myLeaderGroupRank) = offsetComputeNodeNode
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsNode,recvcountCN,displsCN,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountNode(iProc-1) = displsNode(iProc)-displsNode(iProc-1)
  END DO
  recvcountNode(nLeaderGroupProcs-1) = nNonUniqueGlobalNodes - displsNode(nLeaderGroupProcs-1)
END IF

! Broadcast compute node tree offset on node
offsetComputeNodeTree=offsetTree
CALL MPI_BCAST(offsetComputeNodeTree,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF(myComputeNodeRank.EQ.0)THEN
  ! Arrays for the compute node to hold the node offsets
  ALLOCATE(displsTree(0:nLeaderGroupProcs-1))
  ALLOCATE(recvcountTree(0:nLeaderGroupProcs-1))
  displsTree(myLeaderGroupRank) = offsetComputeNodeTree
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsTree,recvcountCN,displsCN,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountTree(iProc-1) = displsTree(iProc)-displsTree(iProc-1)
  END DO
  recvcountTree(nLeaderGroupProcs-1) = nNonUniqueGlobalTrees - displsTree(nLeaderGroupProcs-1)

  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemInfo_Shared,ELEMINFOSIZE       *recvcountElem  &
      ,ELEMINFOSIZE*displsElem    ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)   *recvcountSide  &
      ,(SIDEINFOSIZE+1)*displsSide,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeInfo_Shared,                    recvcountNode  &
      ,displsNode                 ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeCoords_Shared,3                *recvcountNode  &
      ,3*displsNode               ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,IERROR)
  IF(isMortarMesh)THEN
    CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,xiMinMax_Shared,3*2              *recvcountElem  &
        ,3*2*displsElem           ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,IERROR)
    CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemToTree_Shared ,               recvcountElem  &
        ,displsElem               ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
    CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,TreeCoords_Shared ,(NGeoTree+1)**3*recvcountTree &
        ,displsTree               ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,IERROR)
  END IF
END IF
#endif  /*USE_MPI*/
ElemInfo_Shared(ELEM_RANK,     offsetElem+1  :offsetElem+nElems)     = ElemInfo_Shared_tmp
SideInfo_Shared(SIDEINFOSIZE+1,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo_Shared_tmp
DEALLOCATE(ElemInfo_Shared_tmp,SideInfo_Shared_tmp)

! final sync of all mesh shared arrays
CALL MPI_WIN_SYNC(ElemInfo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideInfo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(NodeInfo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(NodeCoords_Shared_Win,IERROR)
IF (isMortarMesh) THEN
  CALL MPI_WIN_SYNC(xiMinMax_Shared_Win,IERROR)
  CALL MPI_WIN_SYNC(ElemToTree_Shared_Win,IERROR)
  CALL MPI_WIN_SYNC(TreeCoords_Shared_Win,IERROR)
END IF

CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

END SUBROUTINE CommunicateMeshReadin


END MODULE MOD_Particle_Mesh_Readin
