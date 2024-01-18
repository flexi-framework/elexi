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

MODULE MOD_Particle_Mesh_Tools
!===================================================================================================================================
! Contains subroutines for mappings of particle mesh
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetGlobalElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalElemIDInterface),POINTER :: GetGlobalElemID    !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetCNElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetCNElemIDInterface),POINTER     :: GetCNElemID        !< pointer defining the mapping: global element ID -> compute-node element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetGlobalSideIDInterface(iSide)
    INTEGER,INTENT(IN) :: iSide
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalSideIDInterface),POINTER :: GetGlobalSideID    !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetCNSideIDInterface(iSide)
    INTEGER,INTENT(IN) :: iSide
  END FUNCTION
END INTERFACE

PROCEDURE(GetCNSideIDInterface),POINTER     :: GetCNSideID        !< pointer defining the mapping: global element ID -> compute-node element ID

! Initialization routines
INTERFACE InitGetGlobalElemID
  MODULE PROCEDURE InitGetGlobalElemID
END INTERFACE

INTERFACE InitGetCNElemID
  MODULE PROCEDURE InitGetCNElemID
END INTERFACE

INTERFACE InitGetGlobalSideID
  MODULE PROCEDURE InitGetGlobalSideID
END INTERFACE

INTERFACE InitGetCNSideID
  MODULE PROCEDURE InitGetCNSideID
END INTERFACE

INTERFACE GetGlobalElemID
  PROCEDURE GetGlobalElemID
END INTERFACE

INTERFACE GetCNElemID
  PROCEDURE GetCNElemID
END INTERFACE

INTERFACE GetGlobalSideID
  PROCEDURE GetGlobalSideID
END INTERFACE

INTERFACE GetCNSideID
  PROCEDURE GetCNSideID
END INTERFACE

INTERFACE GetGlobalNonUniqueSideID
  MODULE PROCEDURE GetGlobalNonUniqueSideID
END INTERFACE

INTERFACE GetSideBoundingBoxTria
  MODULE PROCEDURE GetSideBoundingBoxTria
END INTERFACE

PUBLIC :: InitGetGlobalElemID
PUBLIC :: InitGetCNElemID
PUBLIC :: InitGetGlobalSideID
PUBLIC :: InitGetCNSideID
PUBLIC :: GetGlobalElemID
PUBLIC :: GetCNElemID
PUBLIC :: GetGlobalSideID
PUBLIC :: GetCNSideID
PUBLIC :: GetGlobalNonUniqueSideID
PUBLIC :: GetSideBoundingBoxTria
PUBLIC :: GetCornerNodes
PUBLIC :: GetCornerNodeMapCGNS
!===================================================================================================================================

CONTAINS


!==================================================================================================================================!
!> Initialize GetGlobalElemID function (mapping of compute-node element ID to global element ID)
!==================================================================================================================================!
SUBROUTINE InitGetGlobalElemID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetGlobalElemID => GetGlobalElemID_iElem
ELSE
  GetGlobalElemID => GetGlobalElemID_fromTotalElem
END IF
#else
GetGlobalElemID => GetGlobalElemID_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy = GetGlobalElemID_fromTotalElem(1)
#endif
dummy = GetGlobalElemID_iElem(1)
END SUBROUTINE InitGetGlobalElemID


!==================================================================================================================================!
!> Initialize GetGlobalSideID function (mapping of compute-node side ID to global side ID)
!==================================================================================================================================!
SUBROUTINE InitGetGlobalSideID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetGlobalSideID => GetGlobalSideID_iSide
ELSE
  GetGlobalSideID => GetGlobalSideID_fromTotalSide
END IF
#else
GetGlobalSideID => GetGlobalSideID_iSide
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalSideID_fromTotalSide(1)
#endif
dummy=GetGlobalSideID_iSide(1)
END SUBROUTINE InitGetGlobalSideID


!==================================================================================================================================!
!> Get the compute-node element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION GetGlobalElemID_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: GetGlobalElemID_iElem
!===================================================================================================================================
GetGlobalElemID_iElem = iElem
END FUNCTION GetGlobalElemID_iElem


!==================================================================================================================================!
!> Get the compute-node element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION GetGlobalSideID_iSide(iSide)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: GetGlobalSideID_iSide
!===================================================================================================================================
GetGlobalSideID_iSide = iSide
END FUNCTION GetGlobalSideID_iSide


#if USE_MPI
!==================================================================================================================================!
!> Get the global element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION GetGlobalElemID_fromTotalElem(iElem)
! MODULES
USE MOD_Particle_Mesh_Vars, ONLY:CNTotalElem2GlobalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalElemID_fromTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalElemID_fromTotalElem = CNTotalElem2GlobalElem(iElem)
END FUNCTION GetGlobalElemID_fromTotalElem
#endif /*USE_MPI*/


#if USE_MPI
!==================================================================================================================================!
!> Get the global element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION GetGlobalSideID_fromTotalSide(iSide)
! MODULES
USE MOD_Particle_Mesh_Vars, ONLY:CNTotalSide2GlobalSide
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalSideID_fromTotalSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalSideID_fromTotalSide = CNTotalSide2GlobalSide(iSide)
END FUNCTION GetGlobalSideID_fromTotalSide
#endif /*USE_MPI*/


!==================================================================================================================================!
!> Initialize GetCNElemID function (mapping of global element ID to compute-node element ID)
!==================================================================================================================================!
SUBROUTINE InitGetCNElemID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetCNElemID => CNElemID_is_iElem
ELSE
  GetCNElemID => GetGlobalElem2CNTotalElem
END IF
#else
GetCNElemID => CNElemID_is_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalElem2CNTotalElem(1)
#endif
dummy=CNElemID_is_iElem(1)
END SUBROUTINE InitGetCNElemID


!==================================================================================================================================!
!> Initialize GetCNSideID function (mapping of global element ID to compute-node element ID)
!==================================================================================================================================!
SUBROUTINE InitGetCNSideID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetCNSideID => CNSideID_is_iSide
ELSE
  GetCNSideID => GetGlobalSide2CNTotalSide
END IF
#else
GetCNSideID => CNSideID_is_iSide
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalSide2CNTotalSide(1)
#endif
dummy=CNSideID_is_iSide(1)
END SUBROUTINE InitGetCNSideID


!==================================================================================================================================!
!> Get the CN element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION CNElemID_is_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem ! Global and local element ID are the same
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: CNElemID_is_iElem
!===================================================================================================================================
CNElemID_is_iElem = iElem
END FUNCTION CNElemID_is_iElem


!==================================================================================================================================!
!> Get the CN element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION CNSideID_is_iSide(iSide)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide ! Global and local element ID are the same
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: CNSideID_is_iSide
!===================================================================================================================================
CNSideID_is_iSide = iSide
END FUNCTION CNSideID_is_iSide


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION GetGlobalElem2CNTotalElem(iElem)
! MODULES
USE MOD_Particle_Mesh_Vars, ONLY:GlobalElem2CNTotalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem ! Global element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalElem2CNTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalElem2CNTotalElem = GlobalElem2CNTotalElem(iElem)
END FUNCTION GetGlobalElem2CNTotalElem
#endif /*USE_MPI*/


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE INTEGER FUNCTION GetGlobalSide2CNTotalSide(iSide)
! MODULES
USE MOD_Particle_Mesh_Vars, ONLY:GlobalSide2CNTotalSide
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide ! Global element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalSide2CNTotalSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalSide2CNTotalSide = GlobalSide2CNTotalSide(iSide)
END FUNCTION GetGlobalSide2CNTotalSide
#endif /*USE_MPI*/


PURE INTEGER FUNCTION GetGlobalNonUniqueSideID(ElemID,localSideID)
!===================================================================================================================================
!> Determines the non-unique global side ID of the local side in global element ElemID
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars   ,ONLY: ElemInfo_Shared, SideInfo_Shared
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID                              !< global element ID
INTEGER,INTENT(IN) :: localSideID                         !< local side id of an element (1:6)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!INTEGER :: GetGlobalNonUniqueSideID
INTEGER :: iSide,firstSide,lastSide
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + 1
lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND, ElemID)

! Default to -1
GetGlobalNonUniqueSideID = -1

! Small mortar sides are added after
DO iSide = firstSide,lastSide
  IF (SideInfo_Shared(SIDE_LOCALID,iSide).EQ.localSideID) THEN
    GetGlobalNonUniqueSideID = iSide
    RETURN
  END IF
END DO

END FUNCTION GetGlobalNonUniqueSideID


!==================================================================================================================================!
!> Determines the bounding box of a TriaSurfaceSide
!==================================================================================================================================!
SUBROUTINE GetSideBoundingBoxTria(SideID, BoundingBox)
! MODULES
USE MOD_Mesh_Vars               ,ONLY: SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemSideNodeID_Shared
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: SideID
REAL, INTENT(OUT)             :: BoundingBox(1:3,1:8)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iLocSide,CNElemID,iNode
REAL                      :: NodePoints(1:3,1:4)
REAL                      :: xMin,xMax,yMin,yMax,zMin,zMax
!==================================================================================================================================

CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)

DO iNode = 1, 4
  NodePoints(1:3,iNode) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)
END DO

xMin = MINVAL(NodePoints(1,:))
yMin = MINVAL(NodePoints(2,:))
zMin = MINVAL(NodePoints(3,:))
xMax = MAXVAL(NodePoints(1,:))
yMax = MAXVAL(NodePoints(2,:))
zMax = MAXVAL(NodePoints(3,:))

BoundingBox(1:3,1) = (/xMin,yMin,zMin/)
BoundingBox(1:3,2) = (/xMax,yMin,zMin/)
BoundingBox(1:3,3) = (/xMax,yMax,zMin/)
BoundingBox(1:3,4) = (/xMin,yMax,zMin/)
BoundingBox(1:3,5) = (/xMin,yMin,zMax/)
BoundingBox(1:3,6) = (/xMax,yMin,zMax/)
BoundingBox(1:3,7) = (/xMax,yMax,zMax/)
BoundingBox(1:3,8) = (/xMin,yMax,zMax/)

END SUBROUTINE GetSideBoundingBoxTria


PPURE FUNCTION GetCornerNodes(NGeoLoc)
!===================================================================================================================================
!> Get the corner nodes, depending on the given NGeo, for high-order elements
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: NGeoLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                        :: GetCornerNodes(8)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

GetCornerNodes(1) = 1
GetCornerNodes(2) = (NGeoLoc+1)
GetCornerNodes(3) = (NGeoLoc+1)*NGeoLoc+1
GetCornerNodes(4) = (NGeoLoc+1)**2
GetCornerNodes(5) = (NGeoLoc+1)**2*NGeoLoc+1
GetCornerNodes(6) = (NGeoLoc+1)**2*NGeoLoc+(NGeoLoc+1)
GetCornerNodes(7) = (NGeoLoc+1)**2*NGeoLoc+(NGeoLoc+1)*NGeoLoc+1
GetCornerNodes(8) = (NGeoLoc+1)**3

END FUNCTION GetCornerNodes


PPURE SUBROUTINE GetCornerNodeMapCGNS(NGeoLoc,CornerNodesCGNS,NodeMapCGNS)
!===================================================================================================================================
!> Get the corner nodes and node map when converting from CGNS format, depending on the given NGeo.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: NGeoLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT),OPTIONAL   :: CornerNodesCGNS(8),NodeMapCGNS(4,6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: CNS(8)
!===================================================================================================================================

! CornerNodeSwitch: a mapping is required for Ngeo > 1 as the corner nodes are not the first 8 entries of NodeInfo array
! For NGeo = 1, CNS(3) = 4 and vice versa, and CNS(7) = 8 and vice versa
CNS(1) = 1
CNS(2) = (NGeoLoc+1)
CNS(3) = (NGeoLoc+1)**2
CNS(4) = (NGeoLoc+1)*NGeoLoc+1
CNS(5) = (NGeoLoc+1)**2*NGeoLoc+1
CNS(6) = (NGeoLoc+1)**2*NGeoLoc+(NGeoLoc+1)
CNS(7) = (NGeoLoc+1)**3
CNS(8) = (NGeoLoc+1)**2*NGeoLoc+(NGeoLoc+1)*NGeoLoc+1

IF(PRESENT(CornerNodesCGNS)) CornerNodesCGNS = CNS

IF(PRESENT(NodeMapCGNS)) THEN
  ! Set the corresponding mapping for HOPR coordinates in CGNS format
  NodeMapCGNS(:,1) = (/CNS(1),CNS(4),CNS(3),CNS(2)/)
  NodeMapCGNS(:,2) = (/CNS(1),CNS(2),CNS(6),CNS(5)/)
  NodeMapCGNS(:,3) = (/CNS(2),CNS(3),CNS(7),CNS(6)/)
  NodeMapCGNS(:,4) = (/CNS(3),CNS(4),CNS(8),CNS(7)/)
  NodeMapCGNS(:,5) = (/CNS(1),CNS(5),CNS(8),CNS(4)/)
  NodeMapCGNS(:,6) = (/CNS(5),CNS(6),CNS(7),CNS(8)/)
END IF

END SUBROUTINE GetCornerNodeMapCGNS

END MODULE MOD_Particle_Mesh_Tools
