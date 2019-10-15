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
#include "flexi.h"

!===================================================================================================================================
! Contains routines to build the halo region using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_Particle_MPI_Shared
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
#if USE_MPI_SHARED
INTERFACE InitMeshShared
  MODULE PROCEDURE InitMeshShared
END INTERFACE

INTERFACE InitParticleMeshShared
  MODULE PROCEDURE InitParticleMeshShared
END INTERFACE

INTERFACE FinalizeMeshShared
  MODULE PROCEDURE FinalizeMeshShared
END INTERFACE

INTERFACE FinalizeParticleMeshShared
  MODULE PROCEDURE FinalizeParticleMeshShared
END INTERFACE

PUBLIC:: InitMeshShared
PUBLIC:: InitParticleMeshShared
PUBLIC:: FinalizeMeshShared
PUBLIC:: FinalizeParticleMeshShared
!===================================================================================================================================

CONTAINS

SUBROUTINE InitMeshShared()
!===================================================================================================================================
! Loads the mesh on the current node into a MPI-3 shared memory window
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc,            ONLY:N
USE MOD_Mesh_Vars,          ONLY:nGeo,nElems,nSides,nMPISides_YOUR,offsetElem
USE MOD_Mesh_Vars,          ONLY:Elem_xGP,ElemToSide,SideToElem,NodeCoords
USE MOD_Mesh_Vars,          ONLY:useCurveds
USE MOD_MPI_Shared,         ONLY:Allocate_Shared
USE MOD_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars, ONLY:offsetSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=MPI_ADDRESS_KIND)  :: MPISharedSize
INTEGER                         :: nElems_Shared_Glob
INTEGER                         :: FirstElemShared,LastElemShared
INTEGER                         :: FirstSideShared,LastSideShared
INTEGER                         :: iElem,iSide,ilocSide
!=================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,I1,A)') ' INIT SHARED MESH...'

!> First communicate total number of elems on the node
CALL MPI_ALLREDUCE(nElems,nElems_Shared,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)

!> Sanity check number of cells per node
CALL MPI_REDUCE(nElems_Shared,nElems_Shared_Glob,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
SWRITE(UNIT_stdOut,'(A,F14.1,A)') ' | Copying information for ',REAL(nElems_Shared_Glob)/REAL(nProcessors_Shared),' cells/node in shared memory'

!> Send offsetElem of node root to all other procs on node
IF (myRank_shared.EQ.0) offsetElem_shared_root = offsetElem
CALL MPI_BCAST(offsetElem_shared_root,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

!> Calculate the local offset relative to the node MPI root
FirstElemShared = offsetElem-offsetElem_shared_root+1
LastElemShared  = offsetElem-offsetElem_shared_root+nElems

!> Then communicate total number of sides on the node
CALL MPI_ALLREDUCE(nSides-nMPISides_YOUR,nSides_Shared,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)

!> Send offsetSide of node root to all other procs on node
IF (myRank_shared.EQ.0) offsetSide_shared_root = offsetSide

!> Calculate the local offset relative to the node MPI root. MPI sides on neighbor processors are ignored
FirstSideShared = offsetSide-offsetSide_shared_root+1
LastSideShared  = offsetSide-offsetSide_shared_root+nSides-nMPISides_YOUR

!==== NodeCoords ================================================================================================================
!> DataSizeLength for NodeCoords
IF (useCurveds) THEN
  MPISharedSize = INT(3*(NGeo+1)*(NGeo+1)*(NGeo+1)*nElems_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,NGeo+1,NGeo+1,NGeo+1,nElems_Shared/),NodeCoords_Shared_Win,NodeCoords_Shared)
ELSE
  MPISharedSize = INT(3*(1   +1)*(1   +1)*(1   +1)*nElems_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,1   +1,   1+1,   1+1,nElems_Shared/),NodeCoords_Shared_Win,NodeCoords_Shared)
END IF

!> NodeCoords from each proc
!>>> Caution: NGeo starts from 0 in NodeCoords, but from 1 in NodeCoords_Shared
CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
IF (useCurveds) THEN
  NodeCoords_Shared(:,1:NGeo+1,1:Ngeo+1,1:NGeo+1,FirstElemShared:LastElemShared) = NodeCoords(:,:,:,:,:)
ELSE
  NodeCoords_Shared(:,1:1   +1,1:1   +1,1:1   +1,FirstElemShared:LastElemShared) = NodeCoords(:,:,:,:,:)
END IF

! Deallocate the leftover NodeCoords array
DEALLOCATE(NodeCoords)

!==== Elem_xGP ==================================================================================================================
!> DataSizeLength for Elem_xGP
MPISharedSize = INT(3*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*nElems_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,PP_N+1,PP_N+1,PP_NZ+1,nElems_Shared/),Elem_xGP_Shared_Win,Elem_xGP_Shared)

!> Elem_xGP from each proc
!>>> Caution: PP starts from 0 in Elem_xGP, but from 1 in Elem_xGP_Shared
CALL MPI_WIN_LOCK_ALL(0,Elem_xGP_Shared_Win,IERROR)
Elem_xGP_Shared(:,1:PP_N+1,1:PP_N+1,1:PP_NZ+1,FirstElemShared:LastElemShared) = Elem_xGP(:,:,:,:,:)

!==== ElemToSide ================================================================================================================
!> DataSizeLength for ElemToSide
MPISharedSize = INT(5*nElems_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,6,nElems_Shared/),ElemToSide_Shared_Win,ElemToSide_Shared)

!> ElemToSide from each proc
CALL MPI_WIN_LOCK_ALL(0,ElemToSide_Shared_Win,IERROR)

!> Shift shared sides by FirstSideShared
DO iElem = 1,nElems
  DO ilocSide = 1,6
    ElemToSide_Shared(E2S_SIDE_ID,ilocSide,iElem+FirstElemShared-1) = ElemToSide(E2S_SIDE_ID,ilocSide,iElem) + FirstSideShared - 1
    ElemToSide_Shared(E2S_FLIP   ,ilocSide,iElem+FirstElemShared-1) = ElemToSide(E2S_FLIP   ,ilocSide,iElem)
  END DO
END DO

!==== SideToElem ================================================================================================================
!> DataSizeLength for ElemToSide
MPISharedSize = INT(5*nSides_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/5,nSides_Shared/),SideToElem_Shared_Win,SideToElem_Shared)

!> ElemToSide from each proc
CALL MPI_WIN_LOCK_ALL(0,SideToElem_Shared_Win,IERROR)

!> Shift shared sides by FirstSideShared
DO iSide = 1,nSides-nMPISides_YOUR
  SideToElem_Shared(S2E_ELEM_ID       ,iSide+FirstSideShared-1) = SideToElem(S2E_ELEM_ID       ,iSide) + FirstElemShared - 1
  SideToElem_Shared(S2E_NB_ELEM_ID    ,iSide+FirstSideShared-1) = SideToElem(S2E_NB_ELEM_ID    ,iSide) + FirstElemShared - 1
  SideToElem_Shared(S2E_LOC_SIDE_ID   ,iSide+FirstSideShared-1) = SideToElem(S2E_LOC_SIDE_ID   ,iSide) + FirstSideShared - 1
  SideToElem_Shared(S2E_NB_LOC_SIDE_ID,iSide+FirstSideShared-1) = SideToElem(S2E_NB_LOC_SIDE_ID,iSide) + FirstSideShared - 1
  SideToElem_Shared(S2E_FLIP          ,iSide+FirstSideShared-1) = SideToElem(S2E_FLIP          ,iSide)
END DO

! Synchronize all RMA communication
CALL MPI_WIN_SYNC(Elem_xGP_Shared_Win,  IERROR)
CALL MPI_WIN_SYNC(NodeCoords_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemToSide_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideToElem_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

SWRITE(UNIT_stdOut,'(A)')' INIT SHARED MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMeshShared


SUBROUTINE InitParticleMeshShared
!===================================================================================================================================
! Loads the particle mesh on the current node into a MPI-3 shared memory window
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,          ONLY:nGeo,nElems,nSides,nMPISides_YOUR,offsetElem
USE MOD_Mesh_Vars,          ONLY:useCurveds
USE MOD_MPI_Shared,         ONLY:Allocate_Shared
USE MOD_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars, ONLY:offsetSide,XiEtaZetaBasis,slenXiEtaZetaBasis
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=MPI_ADDRESS_KIND)  :: MPISharedSize
INTEGER                         :: FirstElemShared,LastElemShared
INTEGER                         :: FirstSideShared,LastSideShared
!=================================================================================================================================

!> Calculate the local offset relative to the node MPI root
FirstElemShared = offsetElem-offsetElem_shared_root+1
LastElemShared  = offsetElem-offsetElem_shared_root+nElems

!> Calculate the local offset relative to the node MPI root. MPI sides on neighbor processors are ignored
FirstSideShared = offsetSide-offsetSide_shared_root+1
LastSideShared  = offsetSide-offsetSide_shared_root+nSides-nMPISides_YOUR

!==== BezierControlPoint3D ======================================================================================================
!> DataSizeLength for BezierControlPoint3D
IF (useCurveds) THEN
  MPISharedSize = INT(3*(NGeo+1)*(NGeo+1)*nSides_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,NGeo+1,NGeo+1,nSides_Shared/),BezierControlPoints3D_Shared_Win,BezierControlPoints3D_Shared)
ELSE
  MPISharedSize = INT(3*(1   +1)*(1   +1)*nSides_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,   1+1,   1+1,nSides_Shared/),BezierControlPoints3D_Shared_Win,BezierControlPoints3D_Shared)
END IF

!> BezierControlPoints3D from each proc
!>>> Caution: NGeo starts from 0 in BezierControlPoints3D, but from 1 in BezierControlPoints3D_Shared
CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3D_Shared_Win,IERROR)
IF (useCurveds) THEN
  BezierControlPoints3D_Shared(:,1:NGeo+1,1:Ngeo+1,FirstSideShared:LastSideShared) = &
         BezierControlPoints3D(:,:,:,1:nSides-nMPISides_YOUR)
ELSE
  BezierControlPoints3D_Shared(:,1:1   +1,1:1   +1,FirstSideShared:LastSideShared) = &
         BezierControlPoints3D(:,:,:,1:nSides-nMPISides_YOUR)
END IF

! !==== XiEtaZetaBasis ============================================================================================================
! !> DataSizeLength for XiEtaZetaBasis
! MPISharedSize = INT(3*6*nElems_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
! CALL Allocate_Shared(MPISharedSize,(/3,6,nElems_Shared/),XiEtaZetaBasis_Shared_Win,XiEtaZetaBasis_Shared)
!
! !> XiEtaZetaBasis each proc
! CALL MPI_WIN_LOCK_ALL(0,XiEtaZetaBasis_Shared_Win,IERROR)
! XiEtaZetaBasis_Shared(:,:,FirstElemShared:LastElemShared) = XiEtaZetaBasis(:,:,:)
!
! !==== slenXiEtaZetaBasis ========================================================================================================
! !> DataSizeLength for slenXiEtaZetaBasis
! MPISharedSize = INT(6*nElems_Shared,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
! CALL Allocate_Shared(MPISharedSize,(/6,nElems_Shared/),slenXiEtaZetaBasis_Shared_Win,slenXiEtaZetaBasis_Shared)
!
! !> slenXiEtaZetaBasis each proc
! CALL MPI_WIN_LOCK_ALL(0,slenXiEtaZetaBasis_Shared_Win,IERROR)
! slenXiEtaZetaBasis_Shared(:,FirstElemShared:LastElemShared) = slenXiEtaZetaBasis(:,:)


! Synchronize all RMA communication
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
! CALL MPI_WIN_SYNC(XiEtaZetaBasis_Shared_Win       ,IERROR)
! CALL MPI_WIN_SYNC(slenXiEtaZetaBasis_Shared_Win   ,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

END SUBROUTINE InitParticleMeshShared


SUBROUTINE FinalizeMeshShared
!===================================================================================================================================
! Finalizes the mesh on the current node into a MPI-3 shared memory window
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MPI_Shared_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================

! Free RMA windows
CALL MPI_WIN_FREE(Elem_xGP_Shared_Win,  IERROR)
CALL MPI_WIN_FREE(NodeCoords_Shared_Win,IERROR)
CALL MPI_WIN_FREE(ElemToSide_Shared_Win,IERROR)
CALL MPI_WIN_FREE(SideToElem_Shared_Win,IERROR)

END SUBROUTINE FinalizeMeshShared


SUBROUTINE FinalizeParticleMeshShared
!===================================================================================================================================
! Finalizes the particle mesh on the current node into a MPI-3 shared memory window
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_MPI_Shared_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================

! Free RMA windows
CALL MPI_WIN_FREE(BezierControlPoints3D_Shared_Win,IERROR)
! CALL MPI_WIN_FREE(XiEtaZetaBasis_Shared_Win       ,IERROR)
! CALL MPI_WIN_FREE(slenXiEtaZetaBasis_Shared_Win   ,IERROR)

END SUBROUTINE FinalizeParticleMeshShared

#endif

END MODULE MOD_Particle_MPI_Shared
