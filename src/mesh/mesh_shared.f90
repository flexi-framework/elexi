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
! Builds the mesh for particle tracking, separate from the DG mesh
!===================================================================================================================================
MODULE MOD_Mesh_Shared
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables
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

INTERFACE StartCommunicateMeshReadin
  MODULE PROCEDURE StartCommunicateMeshReadin
END INTERFACE

INTERFACE FinishCommunicateMeshReadin
  MODULE PROCEDURE FinishCommunicateMeshReadin
END INTERFACE

INTERFACE FinalizeMeshShared
  MODULE PROCEDURE FinalizeMeshShared
END INTERFACE

PUBLIC :: ReadMeshBasics
PUBLIC :: ReadMeshElems
PUBLIC :: ReadMeshSides
PUBLIC :: ReadMeshSideNeighbors
PUBLIC :: ReadMeshNodes
PUBLIC :: ReadMeshTrees
PUBLIC :: StartCommunicateMeshReadin
PUBLIC :: FinishCommunicateMeshReadin
PUBLIC :: FinalizeMeshShared
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadMeshBasics(FileString)
!===================================================================================================================================
! Read basic global counters from mesh file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input                ,ONLY: File_ID,ReadAttribute
USE MOD_HDF5_Input                ,ONLY: OpenDataFile,CloseDataFile
USE MOD_Mesh_Vars                 ,ONLY: NGeo,nGlobalElems
USE MOD_Mesh_Vars                 ,ONLY: nNonUniqueGlobalSides,nNonUniqueGlobalNodes
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString !< (IN) mesh filename
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

! INTEGER KIND=4 check for number of nodes
CHECKSAFEINT((NGeo+1)**3.*INT(nGlobalElems,8),4)

CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
!CALL ReadAttribute(File_ID,'nUniqueSides',1,IntScalar=nGlobalUniqueSidesFromMesh)
CALL ReadAttribute(File_ID,'nSides'      ,1,IntScalar=nNonUniqueGlobalSides)
CALL ReadAttribute(File_ID,'nNodes'      ,1,IntScalar=nNonUniqueGlobalNodes)
CALL CloseDataFile()

END SUBROUTINE ReadMeshBasics


SUBROUTINE ReadMeshElems()
!===================================================================================================================================
! Create shared mesh arrays for elems
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Vars                  ,ONLY: offsetElemMPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

#if USE_MPI
! allocate shared array for ElemInfo
nComputeNodeElems = offsetElemMPI(ComputeNodeRootRank+nComputeNodeProcessors) - offsetElemMPI(ComputeNodeRootRank)

#if USE_LOADBALANCE
! Only update the mapping of element to rank, done in loadbalance_tools.f90
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  ! allocate shared array for ElemInfo
  CALL Allocate_Shared((/ELEMINFOSIZE,nGlobalElems/),ElemInfo_Shared_Win,ElemInfo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemInfo_Shared_Win,IERROR)

  ElemInfo_Shared(1:ELEMINFOSIZE_H5,offsetElem+1:offsetElem+nElems) = ElemInfo(:,:)
  ElemInfo_Shared(ELEM_RANK        ,offsetElem+1:offsetElem+nElems) = myRank
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/
#endif  /*USE_MPI*/

#if USE_MPI
! broadcast elem offset of compute-node root
offsetComputeNodeElem=offsetElem
CALL MPI_BCAST(offsetComputeNodeElem,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

#else
! allocate local array for ElemInfo
nComputeNodeElems = nElems
ALLOCATE(ElemInfo_Shared(1:ELEMINFOSIZE,1:nElems))
ElemInfo_Shared(1:ELEMINFOSIZE_H5,1:nElems) = ElemInfo(:,:)
#endif  /*USE_MPI*/

END SUBROUTINE ReadMeshElems


SUBROUTINE ReadMeshSides()
!===================================================================================================================================
! Create shared mesh arrays for sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: FirstElemInd,LastElemInd,iElem
INTEGER                        :: nSideIDs,offsetSideID,iSide
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo(ELEM_LASTSIDEIND ,LastElemInd)-ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd)

ALLOCATE(SideInfo_Shared_tmp(offsetSideID+1:offsetSideID+nSideIDs))
SideInfo_Shared_tmp = 0

#if USE_MPI
! all procs on my compute-node communicate the number of non-unique sides
CALL MPI_ALLREDUCE(nSideIDs,nComputeNodeSides,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

CALL Allocate_Shared((/SIDEINFOSIZE+1,nNonUniqueGlobalSides/),SideInfo_Shared_Win,SideInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideInfo_Shared_Win,IERROR)
SideInfo_Shared(1                :SIDEINFOSIZE_H5,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE_H5+1:SIDEINFOSIZE+1 ,offsetSideID+1:offsetSideID+nSideIDs) = 0
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
#else
nComputeNodeSides = nSideIDs
ALLOCATE(SideInfo_Shared(1:SIDEINFOSIZE+1,1:nSideIDs))
SideInfo_Shared(1                :SIDEINFOSIZE_H5,1:nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE_H5+1:SIDEINFOSIZE+1 ,1:nSideIDs) = 0
#endif /*USE_MPI*/

! Identify all elements containing a mortar side
ElemInfo_Shared(ELEM_HASMORTAR,offsetElem+1:offsetElem+nElems) = 0
DO iElem = FirstElemInd, LastElemInd
  ! Loop over all sides and check for mortar sides
  DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)
    IF (SideInfo_Shared(SIDE_NBELEMID,iSide).LT.0) THEN
      ElemInfo_Shared(ELEM_HASMORTAR,iElem) = 1
      EXIT
    END IF
  END DO
END DO

END SUBROUTINE ReadMeshSides


SUBROUTINE ReadMeshSideNeighbors(ElemID,SideID)
!===================================================================================================================================
! Fills temporary array to add side neighbors to SideInfo(_Shared)
!===================================================================================================================================
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_Mesh_Vars                 ,ONLY: nComputeNodeElems,offsetComputeNodeElem
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: ElemID
INTEGER,INTENT(IN)             :: SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

IF (ElemID.EQ.0) THEN
  ! no connection
  SideInfo_Shared_tmp(SideID) = 0
ELSE
#if USE_MPI
  IF (ElemID.LE.offsetComputeNodeElem+1 .OR. ElemID.GT.offsetComputeNodeElem+nComputeNodeElems) THEN
    ! neighbour element is outside of compute-node
    SideInfo_Shared_tmp(SideID) = 2
  ELSE
#endif /*USE_MPI*/
    SideInfo_Shared_tmp(SideID) = 1
#if USE_MPI
  END IF
#endif /*USE_MPI*/
END IF

END SUBROUTINE ReadMeshSideNeighbors


SUBROUTINE ReadMeshNodes()
!===================================================================================================================================
! Create shared mesh arrays for nodes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ChangeBasisByDim          ,ONLY: ChangeBasisVolume
USE MOD_HDF5_Input                ,ONLY: ReadArray,OpenDataFile
USE MOD_IO_HDF5                   ,ONLY: CloseDataFile
USE MOD_Interpolation             ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars        ,ONLY: NodeType,NodeTypeVISU
USE MOD_Mesh_Vars                 ,ONLY: MeshFile,NGeo,NGeoOverride,useCurveds
USE MOD_Mesh_Vars                 ,ONLY: nElems,offsetElem,nGlobalElems
USE MOD_Mesh_Vars                 ,ONLY: nComputeNodeElems,nComputeNodeNodes
USE MOD_Mesh_Vars                 ,ONLY: meshScale
USE MOD_Mesh_Vars                 ,ONLY: ElemInfo_Shared
USE MOD_Mesh_Vars                 ,ONLY: NodeCoords_Shared
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
USE MOD_Mesh_Vars                 ,ONLY: nNonUniqueGlobalNodes
USE MOD_Mesh_Vars                 ,ONLY: NodeCoords_Shared_Win
USE MOD_Mesh_Vars                 ,ONLY: ElemInfo_Shared_Win
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iNode,i,j,k
INTEGER                        :: FirstElemInd,LastElemInd,LocElemID
INTEGER                        :: FirstNodeInd,LastNodeInd
INTEGER                        :: nNodeIDs,offsetNodeID
INTEGER,ALLOCATABLE            :: NodeInfo(:)
REAL,ALLOCATABLE               :: NodeCoords_indx(:,:)
INTEGER                        :: CornerNodeIDswitch(8)
REAL,ALLOCATABLE               :: Vdm_EQNGeo_EQNGeoOverride(:,:)
REAL,ALLOCATABLE               :: NodeCoordsTmp(:,:,:,:),NodeCoordsNew(:,:,:,:)
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  CALL MPI_ALLREDUCE(nNodeIDs,nComputeNodeNodes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
  RETURN
END IF
#endif /*USE_LOADBALANCE*/

FirstNodeInd = offsetNodeID+1
LastNodeInd  = offsetNodeID+nNodeIDs

! read local Node Info from data file
ALLOCATE(NodeInfo(FirstNodeInd:LastNodeInd))
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_FLEXI)
CALL ReadArray('GlobalNodeIDs',1,(/nNodeIDs/),offsetNodeID,1,IntArray=NodeInfo)
ALLOCATE(NodeCoords_indx(3,nNodeIDs))
CALL ReadArray('NodeCoords',2,(/3,nNodeIDs/),offsetNodeID,2,RealArray=NodeCoords_indx)
CALL CloseDataFile()

! Keep all nodes if elements are curved and not interpolated
IF (NGeoOverride.LE.0 .AND. (useCurveds.OR.NGeo.EQ.1) .OR. NGeoOverride.EQ.NGeo) THEN

#if USE_MPI
!  ! allocate shared array for NodeInfo
  CALL MPI_ALLREDUCE(nNodeIDs,nComputeNodeNodes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
!  CALL Allocate_Shared((/nNonUniqueGlobalNodes/),NodeInfo_Shared_Win,NodeInfo_Shared)
!  CALL MPI_WIN_LOCK_ALL(0,NodeInfo_Shared_Win,IERROR)
!  NodeInfo_Shared(offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeInfo(:)
!  CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win,MPI_COMM_SHARED)

  CALL Allocate_Shared((/3,nNonUniqueGlobalNodes/),NodeCoords_Shared_Win,NodeCoords_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
  NodeCoords_Shared(:,offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeCoords_indx(:,:)
#else
  nComputeNodeNodes = nNodeIDs
!  ALLOCATE(NodeInfo_Shared(1:nNodeIDs))
!  NodeInfo_Shared(1:nNodeIDs) = NodeInfo(:)
  ALLOCATE(NodeCoords_Shared(3,nNodeIDs))
  NodeCoords_Shared(:,:) = NodeCoords_indx(:,:)
#endif  /*USE_MPI*/

! Reduce NodeCoords if no curved elements are to be used
ELSE IF (NGeoOverride.LE.0 .AND. .NOT.useCurveds .AND. NGeo.GT.1) THEN
  ! the cornernodes are not the first 8 entries (for Ngeo>1) of NodeInfo array so mapping is built
  CornerNodeIDswitch(1)=1
  CornerNodeIDswitch(2)=(Ngeo+1)
  CornerNodeIDswitch(3)=(Ngeo+1)*Ngeo+1
  CornerNodeIDswitch(4)=(Ngeo+1)**2
  CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
  CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
  CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1
  CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2

  ASSOCIATE(CNS => CornerNodeIDswitch)

  ! Only the 8 corner nodes count for nodes. (NGeo+1)**2 = 8
  nComputeNodeNodes = 8*nComputeNodeElems

#if USE_MPI
  CALL Allocate_Shared((/3,8*nGlobalElems/),NodeCoords_Shared_Win,NodeCoords_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
#else
  ALLOCATE(NodeCoords_Shared(3,8*nGlobalElems))
#endif  /*USE_MPI*/

#if USE_MPI
  ! throw away all nodes except the 8 corner nodes of each hexa
  nNonUniqueGlobalNodes = 8*nGlobalElems
#endif  /*USE_MPI*/

  DO iElem = FirstElemInd,LastElemInd
    FirstNodeInd = ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) - offsetNodeID
    ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) = 8*(iElem-1)
    ElemInfo_Shared(ELEM_LASTNODEIND ,iElem) = 8* iElem
    DO iNode = 1,8
      NodeCoords_Shared(:,8*(iElem-1) + iNode) = NodeCoords_indx(:,FirstNodeInd+CNS(iNode))
    END DO
  END DO

  END ASSOCIATE

! Interpolate NodeCoords if NGeo != NGeoOverride
ELSE
  ALLOCATE(Vdm_EQNGeo_EQNGeoOverride(0:NGeoOverride,0:NGeo))
  CALL GetVandermonde(NGeo,NodeTypeVISU,NGeoOverride,NodeTypeVISU,Vdm_EQNGeo_EQNGeoOverride,modal=.FALSE.)

  nComputeNodeNodes = (NGeoOverride+1)**3 * nComputeNodeElems

#if USE_MPI
  CALL Allocate_Shared((/3,(NGeoOverride+1)**3*nGlobalElems/),NodeCoords_Shared_Win,NodeCoords_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
#else
  ALLOCATE(NodeCoords_Shared(3,(NGeoOverride+1)**3*nGlobalElems))
#endif  /*USE_MPI*/

  ALLOCATE( NodeCoordsTmp(1:3,0:NGeo        ,0:NGeo        ,0:NGeo)                                                                &
          , NodeCoordsNew(1:3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride))

#if USE_MPI
  nNonUniqueGlobalNodes = (NGeoOverride+1)**3*nGlobalElems
#endif  /*USE_MPI*/

  DO iElem=FirstElemInd,LastElemInd
    LocElemID = iElem - offsetElem
    ! change ElemInfo_Shared to reflect new NodeCoords
    ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) = (NGeoOverride+1)**3*(iElem-1)
    ElemInfo_Shared(ELEM_LASTNODEIND ,iElem) = (NGeoOverride+1)**3*(iElem)

    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      NodeCoordsTmp(:,i,j,k) = NodeCoords_indx(:,(LocElemID-1)*(NGeo+1)**3 + k*(NGeo+1)**2 + j*(NGeo+1) + i + 1)
    END DO; END DO; END DO

    ! interpolate NodeCoords to new NGeo
    CALL ChangeBasisVolume(3,NGeo,NGeoOverride,Vdm_EQNGeo_EQNGeoOverride                                                           &
                          ,NodeCoordsTmp                                                                                           &
                          ,NodeCoordsNew)

    DO k = 0,NGeoOverride; DO j = 0,NGeoOverride; DO i = 0,NGeoOverride
      NodeCoords_Shared(:,ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) + k*(NGeoOverride+1)**2 + j*(NGeoOverride+1) + i + 1) = NodeCoordsNew(:,i,j,k)
    END DO; END DO; END DO
  END DO
END IF

! Update node counters
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)
FirstNodeInd = offsetNodeID+1
LastNodeInd  = offsetNodeID+nNodeIDs

! scale mesh if desired. Mesh deformation currently not supported!
IF (ABS(meshScale-1.).GT.1e-14) THEN
  NodeCoords_Shared(:,FirstNodeInd:LastNodeInd) = NodeCoords_Shared(:,FirstNodeInd:LastNodeInd) * meshScale
END IF

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeCoords_Shared_Win,MPI_COMM_SHARED)
#endif  /*USE_MPI*/

DEALLOCATE(NodeInfo,NodeCoords_indx)

END SUBROUTINE ReadMeshNodes


SUBROUTINE ReadMeshTrees()
!===================================================================================================================================
! Create shared mesh arrays for trees
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars                 ,ONLY: nElems
USE MOD_Mesh_Vars                 ,ONLY: nNonUniqueGlobalTrees
USE MOD_Mesh_Vars                 ,ONLY: nTrees,NGeoTree
USE MOD_Mesh_Vars                 ,ONLY: xiMinMax,xiMinMax_Shared
USE MOD_Mesh_Vars                 ,ONLY: ElemToTree,ElemToTree_Shared
USE MOD_Mesh_Vars                 ,ONLY: TreeCoords,TreeCoords_Shared
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars           ,ONLY: MPI_COMM_SHARED
USE MOD_Mesh_Vars                 ,ONLY: offsetElem,nGlobalElems
USE MOD_Mesh_Vars                 ,ONLY: offsetTree
USE MOD_Mesh_Vars                 ,ONLY: xiMinMax_Shared_Win
USE MOD_Mesh_Vars                 ,ONLY: ElemToTree_Shared_Win
USE MOD_Mesh_Vars                 ,ONLY: TreeCoords_Shared_Win
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

#if USE_MPI
CALL Allocate_Shared((/3,2,nGlobalElems/),xiMinMax_Shared_Win,xiMinMax_Shared)
CALL MPI_WIN_LOCK_ALL(0,xiMinMax_Shared_Win,IERROR)
xiMinMax_Shared(:,:,offsetElem+1:offsetElem+nElems) = xiMinMax(:,:,:)
CALL Allocate_Shared((/nGlobalElems/),ElemToTree_Shared_Win,ElemToTree_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToTree_Shared_Win,IERROR)
ElemToTree_Shared(offsetElem+1:offsetElem+nElems) = ElemToTree(:)
! allocate shared array for TreeCoords
CALL MPI_ALLREDUCE(nTrees,nNonUniqueGlobalTrees,1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,IERROR)
CALL Allocate_Shared((/3,nGeoTree+1,nGeoTree+1,nGeoTree+1,nNonUniqueGlobalTrees/),TreeCoords_Shared_Win,TreeCoords_Shared)
CALL MPI_WIN_LOCK_ALL(0,TreeCoords_Shared_Win,IERROR)
TreeCoords_Shared(:,:,:,:,offsetTree:offsetTree+nTrees) = TreeCoords(:,:,:,:,:)

CALL BARRIER_AND_SYNC(xiMinMax_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemToTree_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(TreeCoords_Shared_Win,MPI_COMM_SHARED)
#else
nNonUniqueGlobalTrees = nTrees
ALLOCATE(xiMinMax_Shared(3,2,nElems))
xiMinMax_Shared(:,:,:) = xiMinMax(:,:,:)
ALLOCATE(ElemToTree_Shared(nElems))
ElemToTree_Shared(:) = ElemToTree(:)
ALLOCATE(TreeCoords_Shared(3,nGeoTree+1,nGeoTree+1,nGeoTree+1,nNonUniqueGlobalTrees))
TreeCoords_Shared(:,:,:,:,:) = TreeCoords(:,:,:,:,:)
#endif  /*USE_MPI*/

END SUBROUTINE ReadMeshTrees


SUBROUTINE StartCommunicateMeshReadin()
!===================================================================================================================================
! Communicates the readin mesh between MPI leaders
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: StartT
USE MOD_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: nSideIDs,offsetSideID
#if USE_MPI
INTEGER                        :: iProc
INTEGER                        :: offsetNodeID!,nNodeIDs
#endif /*USE_MPI*/
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

! Start timer: finished in FinishCommunicateMeshReadin()
GETTIME(StartT)

#if USE_MPI
! MPISharedInitIsDone is not set if this routine is called from posti
IF (.NOT.MPISharedInitIsDone) RETURN

CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND ,LastElemInd) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
! nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd) - ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)
#else
FirstElemInd = 1
LastElemInd  = nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  ! Update SideInfo with new information
  SideInfo_Shared(SIDE_NBELEMTYPE,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo_Shared_tmp
  DEALLOCATE(SideInfo_Shared_tmp)
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  IF (myComputeNodeRank.EQ.0) THEN
    SWRITE(UNIT_stdOut,'(A)',ADVANCE="NO") ' | Updating mesh on shared memory ...'

    ! Arrays for the compute node to hold the elem offsets
    ALLOCATE(displsElem(   0:nLeaderGroupProcs-1),&
             recvcountElem(0:nLeaderGroupProcs-1))
    displsElem(myLeaderGroupRank) = offsetComputeNodeElem
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsElem,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountElem(iProc-1) = displsElem(iProc)-displsElem(iProc-1)
    END DO
    recvcountElem(nLeaderGroupProcs-1) = nGlobalElems - displsElem(nLeaderGroupProcs-1)
  END IF

  ! Broadcast compute node side offset on node
  offsetComputeNodeSide=offsetSideID
  CALL MPI_BCAST(offsetComputeNodeSide,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

  IF (myComputeNodeRank.EQ.0) THEN
    ! Arrays for the compute node to hold the side offsets
    ALLOCATE(displsSide(   0:nLeaderGroupProcs-1),&
             recvcountSide(0:nLeaderGroupProcs-1))
    displsSide(myLeaderGroupRank) = offsetComputeNodeSide
    CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsSide,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
    DO iProc=1,nLeaderGroupProcs-1
      recvcountSide(iProc-1) = displsSide(iProc)-displsSide(iProc-1)
    END DO
    recvcountSide(nLeaderGroupProcs-1) = nNonUniqueGlobalSides - displsSide(nLeaderGroupProcs-1)

    ! Gather mesh information in a non-blocking way
    ALLOCATE(MPI_COMM_LEADERS_REQUEST(1))
    ! ElemInfo_Shared only needs ELEM_RANK updated, performed in loadbalance_tools.f90
    ! CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemInfo_Shared,ELEMINFOSIZE       *recvcountElem  &
    !     ,ELEMINFOSIZE*displsElem     ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)
    CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)   *recvcountSide  &
        ,(SIDEINFOSIZE+1)*displsSide ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)
  END IF

  ! Broadcast compute node node offset on node
  offsetComputeNodeNode=offsetNodeID
  CALL MPI_BCAST(offsetComputeNodeNode,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

  RETURN
END IF
#endif /*USE_LOADBALANCE*/

SWRITE(UNIT_stdOut,'(A)',ADVANCE="NO") ' | Communicating mesh on shared memory...'

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the elem offsets
  ALLOCATE(displsElem(   0:nLeaderGroupProcs-1),&
           recvcountElem(0:nLeaderGroupProcs-1))
  displsElem(myLeaderGroupRank) = offsetComputeNodeElem
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsElem,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountElem(iProc-1) = displsElem(iProc)-displsElem(iProc-1)
  END DO
  recvcountElem(nLeaderGroupProcs-1) = nGlobalElems - displsElem(nLeaderGroupProcs-1)
END IF

! Broadcast compute node side offset on node
offsetComputeNodeSide=offsetSideID
CALL MPI_BCAST(offsetComputeNodeSide,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the side offsets
  ALLOCATE(displsSide(   0:nLeaderGroupProcs-1),&
           recvcountSide(0:nLeaderGroupProcs-1))
  displsSide(myLeaderGroupRank) = offsetComputeNodeSide
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsSide,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountSide(iProc-1) = displsSide(iProc)-displsSide(iProc-1)
  END DO
  recvcountSide(nLeaderGroupProcs-1) = nNonUniqueGlobalSides - displsSide(nLeaderGroupProcs-1)
END IF

! Broadcast compute node node offset on node
offsetComputeNodeNode=offsetNodeID
CALL MPI_BCAST(offsetComputeNodeNode,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the node offsets
  ALLOCATE(displsNode(   0:nLeaderGroupProcs-1),&
           recvcountNode(0:nLeaderGroupProcs-1))
  displsNode(myLeaderGroupRank) = offsetComputeNodeNode
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsNode,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountNode(iProc-1) = displsNode(iProc)-displsNode(iProc-1)
  END DO
  recvcountNode(nLeaderGroupProcs-1) = nNonUniqueGlobalNodes - displsNode(nLeaderGroupProcs-1)
END IF

! Broadcast compute node tree offset on node
offsetComputeNodeTree=offsetTree
CALL MPI_BCAST(offsetComputeNodeTree,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the node offsets
  ALLOCATE(displsTree(   0:nLeaderGroupProcs-1),&
           recvcountTree(0:nLeaderGroupProcs-1))
  displsTree(myLeaderGroupRank) = offsetComputeNodeTree
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsTree,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountTree(iProc-1) = displsTree(iProc)-displsTree(iProc-1)
  END DO
  recvcountTree(nLeaderGroupProcs-1) = nNonUniqueGlobalTrees - displsTree(nLeaderGroupProcs-1)

  ! Gather mesh information in a non-blocking way
  MPI_COMM_LEADERS_REQUEST_SIZE = MERGE(6,3,isMortarMesh)
  ALLOCATE(MPI_COMM_LEADERS_REQUEST(1:MPI_COMM_LEADERS_REQUEST_SIZE))
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemInfo_Shared,ELEMINFOSIZE       *recvcountElem  &
      ,ELEMINFOSIZE*displsElem     ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)   *recvcountSide  &
      ,(SIDEINFOSIZE+1)*displsSide ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(2),IERROR)
!  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeInfo_Shared,                    recvcountNode  &
!      ,displsNode                  ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(3),IERROR)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NodeCoords_Shared,3                *recvcountNode  &
      ,3*displsNode                ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(3),IERROR)
  IF(isMortarMesh)THEN
    CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,xiMinMax_Shared,3*2              *recvcountElem  &
        ,3*2*displsElem            ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(4),IERROR)
    CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemToTree_Shared ,               recvcountElem  &
        ,displsElem                ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(5),IERROR)
    CALL MPI_IALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,TreeCoords_Shared ,(NGeoTree+1)**3*recvcountTree &
        ,displsTree                ,MPI_DOUBLE_PRECISION,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(6),IERROR)
  END IF
END IF
#endif /*USE_MPI*/

END SUBROUTINE StartCommunicateMeshReadin


SUBROUTINE FinishCommunicateMeshReadin()
!===================================================================================================================================
! Communicates the readin mesh between MPI leaders
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: CommMeshReadinWallTime,StartT
USE MOD_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: FirstElemInd,LastElemInd
INTEGER :: iElem,NbElemID
INTEGER :: nSideIDs,offsetSideID
INTEGER :: iSide,sideCount
INTEGER :: iLocSide,jLocSide,nlocSides,nlocSidesNb,NbSideID
REAL    :: EndT
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  ! Finish non-blocking mesh communication
  IF (myComputeNodeRank.EQ.0) THEN
    CALL MPI_WAITALL(1,MPI_COMM_LEADERS_REQUEST,MPI_STATUSES_IGNORE,IERROR)
    DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
  END IF

  ! final sync of all mesh shared arrays
  CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

  CALL DisplayMessageAndTime(CommMeshReadinWallTime,'DONE!',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
  RETURN
END IF
#endif /*USE_LOADBALANCE*/

#if USE_MPI
! Finish non-blocking mesh communication
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_WAITALL(MPI_COMM_LEADERS_REQUEST_SIZE,MPI_COMM_LEADERS_REQUEST,MPI_STATUSES_IGNORE,IERROR)
  DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
END IF

! Ensure communication for determination of SIDE_LOCALID
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND ,LastElemInd) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#else
FirstElemInd = 1
LastElemInd  = nElems
offsetSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
#endif /*USE_MPI*/

! fill the SIDE_LOCALID. Basically, this array contains the 1:6 local sides of an element. ! If an element has hanging nodes (i.e.
! has a big mortar side), the big side has negative index (-1,-2 or -3) and the next 2 (-2, -3) or 4 (-1) sides are the subsides.
! Consequently, a hexahedral element can have more than 6 non-unique sides. If we find a small mortar side in the small element,
! the SIDE_LOCALID points to the global ID of the big mortar side, indicated by negative sign
!
! This step has to be done after ElemInfo and SideInfo are communicated as some side might be missing on the current node otherwise
DO iElem = FirstElemInd,LastElemInd
  iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  SideInfo_Shared(SIDE_ELEMID,iSide+1:ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)) = iElem
  sideCount = 0
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  DO iLocSide = 1,nlocSides
    iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
    ! Big mortar side
    IF (SideInfo_Shared(SIDE_TYPE,iSide).LE.100) THEN
      sideCount = sideCount + 1
      SideInfo_Shared(SIDE_LOCALID,iSide) = sideCount
    ELSE
      ! Mortar case
      SideInfo_Shared(SIDE_LOCALID,iSide) = -1
    END IF
      ! Check all sides on the small element side to find the small mortar side pointing back
      NbElemID    = SideInfo_Shared(SIDE_NBELEMID,iSide)
    IF (NbElemID.EQ.0) THEN
      SideInfo_Shared(SIDE_NBSIDEID,iSide) = 0
    ELSE IF (NbElemID.LE.-1) THEN
      SideInfo_Shared(SIDE_NBSIDEID,iSide) = -1
    ELSE
      nlocSidesNb = ElemInfo_Shared(ELEM_LASTSIDEIND,NbElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID)
      DO jLocSide = 1,nlocSidesNb
        NbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID) + jLocSide
        IF (ABS(SideInfo_Shared(SIDE_ID,iSide)).EQ.ABS(SideInfo_Shared(SIDE_ID,NbSideID))) THEN
          SideInfo_Shared(SIDE_NBSIDEID,iSide) = NbSideID
          EXIT
        END IF
      END DO
    END IF
  END DO
END DO

#if USE_MPI
! Perform second communication step to distribute updated SIDE_LOCALID
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)

IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,SideInfo_Shared,(SIDEINFOSIZE+1)   *recvcountSide  &
      ,(SIDEINFOSIZE+1)*displsSide,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
END IF

! Write compute-node local SIDE_NBELEMTYPE
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
#endif  /*USE_MPI*/

SideInfo_Shared(SIDE_NBELEMTYPE,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo_Shared_tmp
DEALLOCATE(SideInfo_Shared_tmp)

#if USE_MPI
! final sync of all mesh shared arrays
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win  ,MPI_COMM_SHARED)
! CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win  ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(NodeCoords_Shared_Win,MPI_COMM_SHARED)
IF (isMortarMesh) THEN
  CALL BARRIER_AND_SYNC(xiMinMax_Shared_Win  ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(ElemToTree_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(TreeCoords_Shared_Win,MPI_COMM_SHARED)
END IF
#endif  /*USE_MPI*/

EndT                   = FLEXITIME()
CommMeshReadinWallTime = EndT-StartT
CALL DisplayMessageAndTime(CommMeshReadinWallTime,'DONE',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("."))')

END SUBROUTINE FinishCommunicateMeshReadin


SUBROUTINE FinalizeMeshShared()
!===================================================================================================================================
! Finalizes the shared mesh readin
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! do not build shared mesh information in posti mode
IF (postiMode) RETURN

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! Free communication arrays
SDEALLOCATE(displsElem)
SDEALLOCATE(recvcountElem)
SDEALLOCATE(displsSide)
SDEALLOCATE(recvcountSide)

#if USE_LOADBALANCE
! Keep pure mesh geometry available during loadbalance, distributed arrays will be restored in mesh_readin
IF (PerformLoadBalance) THEN
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  RETURN
END IF
#endif /*USE_LOADBALANCE*/

! elems
CALL MPI_WIN_UNLOCK_ALL(ElemInfo_Shared_Win,iError)
CALL MPI_WIN_FREE(ElemInfo_Shared_Win,iError)

! sides
CALL MPI_WIN_UNLOCK_ALL(SideInfo_Shared_Win,iError)
CALL MPI_WIN_FREE(SideInfo_Shared_Win,iError)

! nodes
!CALL MPI_WIN_UNLOCK_ALL(NodeInfo_Shared_Win,iError)
!CALL MPI_WIN_FREE(NodeInfo_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(NodeCoords_Shared_Win,iError)
CALL MPI_WIN_FREE(NodeCoords_Shared_Win,iError)

! trees
IF (ASSOCIATED(TreeCoords_Shared)) THEN
  CALL MPI_WIN_UNLOCK_ALL(TreeCoords_Shared_Win,iError)
  CALL MPI_WIN_FREE(TreeCoords_Shared_Win,iError)
END IF
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

! Then, free the pointers or arrays
MDEALLOCATE(ElemInfo_Shared)
MDEALLOCATE(SideInfo_Shared)
!MDEALLOCATE(NodeInfo_Shared)
MDEALLOCATE(NodeCoords_Shared)
MDEALLOCATE(TreeCoords_Shared)

! Free communication arrays
#if USE_MPI
SDEALLOCATE(displsNode)
SDEALLOCATE(recvcountNode)
SDEALLOCATE(displsTree)
SDEALLOCATE(recvcountTree)
#endif /*USE_MPI*/

END SUBROUTINE FinalizeMeshShared

END MODULE MOD_Mesh_Shared
