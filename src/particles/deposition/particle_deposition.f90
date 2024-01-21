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
#if PART_TWO_WAY
#include "flexi.h"
#include "particle.h"

!===================================================================================================================================
! Module containing the different deposition methods (NGP, linear (inter-cell) weighting, shape function
!===================================================================================================================================
MODULE MOD_Particle_Deposition
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitializeDeposition
  MODULE PROCEDURE InitializeDeposition
END INTERFACE

INTERFACE FinalizeDeposition
  MODULE PROCEDURE FinalizeDeposition
END INTERFACE

PUBLIC :: InitializeDeposition
PUBLIC :: FinalizeDeposition
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the deposition variables first
!===================================================================================================================================
SUBROUTINE InitializeDeposition()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input                ,ONLY: OpenDataFile,DatasetExists,ReadArray,CloseDataFile
USE MOD_Interpolation_Vars        ,ONLY: xGP
USE MOD_IO_HDF5                   ,ONLY: File_ID
USE MOD_Mesh_Vars
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Deposition_Method
USE MOD_Particle_Mesh_Vars
USE MOD_ReadInTools               ,ONLY: GETINTFROMSTR
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
LOGICAL             :: exists
INTEGER             :: FirstElemInd,LastElemInd
INTEGER             :: FirstNodeInd,LastNodeInd
INTEGER             :: nNodes,offsetNode
#if USE_MPI
INTEGER,ALLOCATABLE :: FEMElemInfo(:,:)
INTEGER,ALLOCATABLE :: VertexInfo(:,:)
#else
INTEGER,PARAMETER   :: myComputeNodeRank = 0
#endif /*USE_MPI*/
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE DEPOSITION...'

! General deposition variables
ALLOCATE(PartSource(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(Ut_src(    1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

DepositionType          = GETINTFROMSTR('Part-DepositionType')

! Associate the DepositionMethod pointer
CALL InitDepositionMethod()

SELECT CASE(DepositionType)
  CASE(DEPO_CV,DEPO_STEP)
    ALLOCATE(PartSource_tmp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ))

  CASE(DEPO_CVLM,DEPO_SF_GAUSS,DEPO_SF_POLY)
    ! Only FEM meshes supported
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL DatasetExists(File_ID,'FEMconnect',exists,attrib=.TRUE.)
    IF (.NOT.exists) CALL Abort(__STAMP__,'Only FEM meshes supported for deposition!')

    ALLOCATE(CellVolWeightFac(  0:PP_N))
    CellVolWeightFac(0:PP_N) = xGP(0:PP_N)
    CellVolWeightFac(0:PP_N) = (CellVolWeightFac(0:PP_N)+1.0)/2.0

    FirstElemInd =  offsetElem+1
    LastElemInd  =  offsetElem+nElems
    FirstNodeInd = (offsetElem*8)+1
    LastNodeInd  = (offsetElem+nElems)*8
    nNodes       =  LastNodeInd - FirstNodeInd + 1
    offsetNode   =  offsetElem*8

#if USE_MPI
#if USE_LOADBALANCE
    ! Only update the mapping of element to rank, done in loadbalance_tools.f90
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ! Read in the HDF5 data
      ALLOCATE(FEMElemInfo(FEMELEMINFOSIZE,FirstElemInd:LastElemInd))
      ALLOCATE(VertexInfo (VERTEXINFOSIZE ,FirstNodeInd:LastNodeInd))
      CALL ReadArray('FEMElemInfo',2,(/FEMELEMINFOSIZE,nElems/),offsetElem,2,IntArray=FEMElemInfo)
      CALL ReadArray('VertexInfo ',2,(/VERTEXINFOSIZE ,nNodes/),offsetNode,2,IntArray=VertexInfo)

      ! allocate shared array for FEMElemInfo
      CALL Allocate_Shared((/FEMELEMINFOSIZE,nGlobalElems  /),FEMElemInfo_Shared_Win,FEMElemInfo_Shared)
      CALL Allocate_Shared((/VERTEXINFOSIZE ,nGlobalElems*8/),VertexInfo_Shared_Win ,VertexInfo_Shared)
      CALL MPI_WIN_LOCK_ALL(0,FEMElemInfo_Shared_Win,iError)
      CALL MPI_WIN_LOCK_ALL(0,VertexInfo_Shared_Win ,iError)

      FEMElemInfo_Shared(1:FEMELEMINFOSIZE,FirstElemInd:LastElemInd) = FEMElemInfo(:,:)
      VertexInfo_Shared (1:VERTEXINFOSIZE ,FirstNodeInd:LastNodeInd) = VertexInfo( :,:)
      CALL BARRIER_AND_SYNC(FEMElemInfo_Shared_Win,MPI_COMM_SHARED)
      CALL BARRIER_AND_SYNC(VertexInfo_Shared_Win ,MPI_COMM_SHARED)
      DEALLOCATE(FEMElemInfo)
      DEALLOCATE(VertexInfo)

      IF (myComputeNodeRank.EQ.0) THEN
        ! Gather mesh information in a blocking way
        CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,FEMElemInfo_Shared,FEMELEMINFOSIZE*recvcountElem    &
            ,FEMELEMINFOSIZE*displsElem     ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
        CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,VertexInfo_Shared ,VERTEXINFOSIZE *recvcountElem*8  &
            ,VERTEXINFOSIZE*displsElem*8    ,MPI_INTEGER         ,MPI_COMM_LEADERS_SHARED,IERROR)
      END IF

      ! Ensure communication finished
      CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
      CALL BARRIER_AND_SYNC(FEMElemInfo_Shared_Win,MPI_COMM_SHARED)
      CALL BARRIER_AND_SYNC(VertexInfo_Shared_Win ,MPI_COMM_SHARED)
#if USE_LOADBALANCE
    END IF
#endif /*USE_LOADBALANCE*/
#else
    ! allocate local array for FEM information
    ALLOCATE(FEMElemInfo_Shared(1:FEMELEMINFOSIZE,FirstElemInd:LastElemInd))
    ALLOCATE(VertexInfo_Shared (1:VERTEXINFOSIZE ,FirstNodeInd:LastNodeInd))
    CALL ReadArray('FEMElemInfo',2,(/FEMELEMINFOSIZE,nElems/),offsetElem,2,IntArray=FEMElemInfo_Shared)
    CALL ReadArray('VertexInfo ',2,(/VERTEXINFOSIZE ,nNodes/),offsetNode,2,IntArray=VertexInfo_Shared)
#endif  /*USE_MPI*/

    CALL CloseDataFile()

    IF (myComputeNodeRank.EQ.0) nUniqueFEMNodes = MAXVAL(VertexInfo_Shared(1,:))
#if USE_MPI
    CALL MPI_BCAST(nUniqueFEMNodes,1,MPI_INTEGER,0,MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

    ALLOCATE(FEMNodeSource(PP_nVar,nUniqueFEMNodes))

END SELECT

LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DEPOSITION DONE!'

END SUBROUTINE InitializeDeposition


!===================================================================================================================================
!> Finalize the deposition variables
!===================================================================================================================================
SUBROUTINE FinalizeDeposition()
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars        ,ONLY: FEMElemInfo_Shared
USE MOD_Particle_Mesh_Vars        ,ONLY: VertexInfo_Shared
USE MOD_Particle_Deposition_Vars
#if USE_MPI
USE MOD_MPI_Shared_Vars           ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars        ,ONLY: FEMElemInfo_Shared_Win
USE MOD_Particle_Mesh_Vars        ,ONLY: VertexInfo_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(PartSource)
SDEALLOCATE(PartSource_tmp)
SDEALLOCATE(Ut_src)

SELECT CASE(DepositionType)
  CASE(DEPO_CVLM)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

    CALL MPI_WIN_UNLOCK_ALL(FEMElemInfo_Shared_Win           ,iError)
    CALL MPI_WIN_FREE(      FEMElemInfo_Shared_Win           ,iError)
    CALL MPI_WIN_UNLOCK_ALL(VertexInfo_Shared_Win            ,iError)
    CALL MPI_WIN_FREE(      VertexInfo_Shared_Win            ,iError)

    CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

    ! Then, free the pointers or arrays
    MDEALLOCATE(FEMElemInfo_Shared)
    MDEALLOCATE(VertexInfo_Shared)
END SELECT

! CALL Abort(__STAMP__,'TODO!')

END SUBROUTINE FinalizeDeposition

END MODULE MOD_Particle_Deposition
#endif /*PART_TWO_WAY*/
