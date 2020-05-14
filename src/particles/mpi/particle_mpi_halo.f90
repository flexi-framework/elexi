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
#include "particle.h"

!===================================================================================================================================
! Contains routines to build the halo exchange
!===================================================================================================================================
MODULE MOD_Particle_MPI_Halo
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE IdentifyPartExchangeProcs
  MODULE PROCEDURE IdentifyPartExchangeProcs
END INTERFACE

INTERFACE FinalizePartExchangeProcs
  MODULE PROCEDURE FinalizePartExchangeProcs
END INTERFACE

PUBLIC :: IdentifyPartExchangeProcs
PUBLIC :: FinalizePartExchangeProcs
!===================================================================================================================================

CONTAINS

SUBROUTINE IdentifyPartExchangeProcs
!===================================================================================================================================
! Identifies processors in physical range for particle exchange communication. This communication has to occur at every RK step and
! would be too costly if done as an all-to-all communication
!> Procs on the same compute node are always assumed to be in range since communication is handled on proc
!> Procs in the compute node halo region are only considered in range if they lie within halo_eps of the mesh on the current proc
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_MPI_Vars                ,ONLY: offsetElemMPI
USE MOD_Particle_Globals        ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeElems
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID,GetGlobalNonUniqueSideID
USE MOD_Particle_MPI_Vars       ,ONLY: halo_eps
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputenodeProcessors,nProcessors_Global,ComputeNodeRootRank
USE MOD_Particle_MPI_Vars       ,ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Partner identification
INTEGER                        :: iProc
INTEGER                        :: iPeriodicVector,jPeriodicVector
INTEGER                        :: iPeriodicDir,jPeriodicDir,kPeriodicDir,CoordDir
INTEGER,DIMENSION(2)           :: DirPeriodicVector = [-1,1]
INTEGER                        :: firstElem,lastElem
REAL,DIMENSION(6)              :: xCoordsProc,xCoordsTest,xCoordsOrigin
! Communication test
INTEGER                        :: nExchangeProcessorsGlobal
INTEGER,ALLOCATABLE            :: TestRecv(:),TestSend(:),SendRequest(:),RecvRequest(:)
!=================================================================================================================================

! Keep everything in sync here
CALL MPI_BARRIER(MPI_COMM_FLEXI,IERROR)

!SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors...'

! Allocate arrays
ALLOCATE(GlobalProcToExchangeProc(EXCHANGE_PROC_SIZE,0:nProcessors_Global-1))
GlobalProcToExchangeProc(:,:) = -1

! Identify all procs on same node
nExchangeProcessors = 0
DO iProc = ComputeNodeRootRank,ComputeNodeRootRank+nComputeNodeProcessors-1
  ! Do not attempt to communicate with myself
  IF (iProc.EQ.myRank) CYCLE

  ! Build mapping global to compute-node
  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
  nExchangeProcessors = nExchangeProcessors + 1
END DO

! Identify all procs with elements in range. If all elements are on the current node, every proc will be added. Otherwise, compare
! Cartesian bounding boxes and add procs within halo_eps
IF (nComputeNodeElems.NE.nComputeNodeTotalElems) THEN
!  ! Simple method, use the Cartesian bounding box around the proc. Add halo_eps while testing
  xCoordsProc(1) = GEO%xmin
  xCoordsProc(2) = GEO%xmax
  xCoordsProc(3) = GEO%ymin
  xCoordsProc(4) = GEO%ymax
  xCoordsProc(5) = GEO%zmin
  xCoordsProc(6) = GEO%zmax

  ! Test against the Cartesian bounding box of every other proc
  ProcLoop: DO iProc = 0,nProcessors_Global-1
    ! Ignore procs on the same node
    IF (iProc.GE.ComputeNodeRootRank .AND. iProc.LE.ComputeNodeRootRank+nComputeNodeProcessors-1) CYCLE

    firstElem = offsetElemMPI(iProc) + 1
    lastElem  = offsetElemMPI(iProc+1)

    xCoordsOrigin(1) = MINVAL(NodeCoords_Shared(1,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                 :ElemInfo_Shared(ELEM_LASTNODEIND,lastElem)))
    xCoordsOrigin(2) = MAXVAL(NodeCoords_Shared(1,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                 :ElemInfo_Shared(ELEM_LASTNODEIND,lastElem)))
    xCoordsOrigin(3) = MINVAL(NodeCoords_Shared(2,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                 :ElemInfo_Shared(ELEM_LASTNODEIND,lastElem)))
    xCoordsOrigin(4) = MAXVAL(NodeCoords_Shared(2,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                 :ElemInfo_Shared(ELEM_LASTNODEIND,lastElem)))
    xCoordsOrigin(5) = MINVAL(NodeCoords_Shared(3,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                 :ElemInfo_Shared(ELEM_LASTNODEIND,lastElem)))
    xCoordsOrigin(6) = MAXVAL(NodeCoords_Shared(3,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                 :ElemInfo_Shared(ELEM_LASTNODEIND,lastElem)))

    ! Check if proc is in range
    IF (HaloBoxInProc(xCoordsOrigin,xCoordsProc,halo_eps)) THEN
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
      GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
      nExchangeProcessors = nExchangeProcessors + 1
      CYCLE ProcLoop
    END IF

    ! Proc out of range, also check periodic displacements
    SELECT CASE(GEO%nPeriodicVectors)
      ! No periodic vectors, proc out of range
      CASE(0)
        ! Do nothing, just don't throw errors

      ! One periodic vector
      CASE(1)
        ! Check both directions
        DO iPeriodicDir = 1,2
          DO CoordDir = 1,3
            xCoordsTest(2*CoordDir-1:2*CoordDir) = xCoordsOrigin(2*CoordDir-1:2*CoordDir)                                        &
                                                 + GEO%PeriodicVectors(CoordDir,1) * DirPeriodicVector(iPeriodicDir)
          END DO

          ! Check if proc is in range
          IF (HaloBoxInProc(xCoordsTest,xCoordsProc,halo_eps)) THEN
            GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
            GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
            nExchangeProcessors = nExchangeProcessors + 1
            CYCLE ProcLoop
          END IF
        END DO

      ! Two periodic vectors
      CASE(2)
        ! Check both directions
        DO iPeriodicVector = 1,2; DO iPeriodicDir = 1,2
          DO CoordDir = 1,3
            xCoordsTest(2*CoordDir-1:2*CoordDir) = xCoordsOrigin(2*CoordDir-1:2*CoordDir)                                        &
                                             + GEO%PeriodicVectors(CoordDir,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)
          END DO
          ! Check if proc is in range
          IF (HaloBoxInProc(xCoordsTest,xCoordsProc,halo_eps)) THEN
            GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
            GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
            nExchangeProcessors = nExchangeProcessors + 1
            CYCLE ProcLoop
          END IF

          ! Also check linear combination, see particle_bgm.f90
          DO jPeriodicVector = 1,2; DO jPeriodicDir = 1,2
            DO CoordDir = 1,3
              xCoordsTest(2*CoordDir-1:2*CoordDir) = xCoordsOrigin(2*CoordDir-1:2*CoordDir)                                      &
                                               + GEO%PeriodicVectors(CoordDir,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                                               + GEO%PeriodicVectors(CoordDir,jPeriodicVector) * DirPeriodicVector(jPeriodicDir)
            END DO
            ! Check if proc is in range
            IF (HaloBoxInProc(xCoordsTest,xCoordsProc,halo_eps)) THEN
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ProcLoop
            END IF
          END DO; END DO ! jPeriodicDir, jPeriodicVector
        END DO; END DO ! iPeriodicDir, iPeriodicVector

      ! Three periodic vectors
      CASE(3)
        ! Check both directions
        DO iPeriodicVector = 1,3; DO iPeriodicDir = 1,2
          DO CoordDir = 1,3
            xCoordsTest(2*CoordDir-1:2*CoordDir) = xCoordsOrigin(2*CoordDir-1:2*CoordDir)                                        &
                                             + GEO%PeriodicVectors(CoordDir,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)
          END DO
          ! Check if proc is in range
          IF (HaloBoxInProc(xCoordsTest,xCoordsProc,halo_eps)) THEN
            GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
            GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
            nExchangeProcessors = nExchangeProcessors + 1
            CYCLE ProcLoop
          END IF

          ! Also check linear combination, see particle_bgm.f90
          DO jPeriodicVector = 1,3; DO jPeriodicDir = 1,2
            ! Only check the lower triangular entries
            IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

            DO CoordDir = 1,3
              xCoordsTest(2*CoordDir-1:2*CoordDir) = xCoordsOrigin(2*CoordDir-1:2*CoordDir)                                      &
                                               + GEO%PeriodicVectors(CoordDir,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                                               + GEO%PeriodicVectors(CoordDir,jPeriodicVector) * DirPeriodicVector(jPeriodicDir)
            END DO
            ! Check if proc is in range
            IF (HaloBoxInProc(xCoordsTest,xCoordsProc,halo_eps)) THEN
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ProcLoop
            END IF
          END DO; END DO ! jPeriodicDir, jPeriodicVector
        END DO; END DO ! iPeriodicDir, iPeriodicVector

        ! Finally, check the linear combination of all three vectors
        DO iPeriodicDir = 1,2; DO jPeriodicDir = 1,2; DO kPeriodicDir = 1,2
            ! Only check the lower triangular entries
            IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

            DO CoordDir = 1,3
              xCoordsTest(2*CoordDir-1:2*CoordDir) = xCoordsOrigin(2*CoordDir-1:2*CoordDir)                                      &
                                               + GEO%PeriodicVectors(CoordDir,1) * DirPeriodicVector(iPeriodicDir)               &
                                               + GEO%PeriodicVectors(CoordDir,2) * DirPeriodicVector(jPeriodicDir)               &
                                               + GEO%PeriodicVectors(CoordDir,3) * DirPeriodicVector(kPeriodicDir)
            END DO
            ! Check if proc is in range
            IF (HaloBoxInProc(xCoordsTest,xCoordsProc,halo_eps)) THEN
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ProcLoop
            END IF
          END DO; END DO; END DO ! iPeriodicDir, jPeriodicDir, kPeriodicDir

      CASE DEFAULT
          CALL ABORT(__STAMP__,'Invalid number of periodic vectors in particle_mpi_halo.f90')

    END SELECT
  END DO ProcLoop
! Node contains the entire mesh, all-to-all communication is required
ELSE
  DO iProc = 0,nProcessors_Global-1
    ! Ignore procs on the same node
    IF (iProc.GE.ComputeNodeRootRank .AND. iProc.LE.ComputeNodeRootRank+nComputeNodeProcessors-1) CYCLE

    ! Add the other procs to the mapping
    GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
    GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
    nExchangeProcessors = nExchangeProcessors + 1
  END DO
END IF

!
ALLOCATE(ExchangeProcToGlobalProc(2,0:nExchangeProcessors-1))

! Loop through all procs and build reverse mapping
nExchangeProcessors = 0

DO iProc = ComputeNodeRootRank,ComputeNodeRootRank + nComputeNodeProcessors - 1
  ! Do not attempt to communicate with myself
  IF (iProc.EQ.myRank) CYCLE

  ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,nExchangeProcessors) = 1
  ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,nExchangeProcessors) = iProc

  nExchangeProcessors = nExchangeProcessors + 1
END DO

DO iProc = 0,nProcessors_Global-1
  ! Only consider procs in halo region
  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).NE.2) CYCLE

  ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,nExchangeProcessors) = 2
  ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,nExchangeProcessors) = iProc

  nExchangeProcessors = nExchangeProcessors +1
END DO

! Test the communication
ALLOCATE( TestRecv   (0:nExchangeProcessors-1) &
        , TestSend   (0:nExchangeProcessors-1) &
        , SendRequest(0:nExchangeProcessors-1) &
        , RecvRequest(0:nExchangeProcessors-1))

TestSend = 1
TestRecv = 0

!--- Asynchronous communication, just recv here and check for success later.
DO iProc = 0,nExchangeProcessors-1
  CALL MPI_IRECV( TestRecv(iProc)                                            &
                , 1                                                          &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , MPI_COMM_FLEXI                                             &
                , RecvRequest(iProc)                                         &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

!--- Asynchronous communication, just send here and check for success later.
DO iProc = 0,nExchangeProcessors-1
  CALL MPI_ISEND( TestSend(iProc)                                            &
                , 1                                                          &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , MPI_COMM_FLEXI                                             &
                , SendRequest(iProc)                                         &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

DO iProc = 0,nExchangeProcessors-1
  CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

!--- Check if any communication has no reverse
IF (ANY(TestRecv.EQ.0)) &
  CALL ABORT(__STAMP__, 'Initialization of particle exchange failed!')

!--- Stop here if any exchange failed
CALL MPI_BARRIER(MPI_COMM_FLEXI,IERROR)

! -- Average number of exchange processors
CALL MPI_REDUCE(nExchangeProcessors,nExchangeProcessorsGlobal,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
SWRITE(UNIT_stdOut,'(A,I0,A)') ' | Started particle exchange communication with average ', &
                                 nExchangeProcessorsGlobal/nProcessors_Global            , &
                                 ' partners per proc'

SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE IdentifyPartExchangeProcs


!===================================================================================================================================
! Deallocates arrays for halo exchange
!===================================================================================================================================
SUBROUTINE FinalizePartExchangeProcs()
! MODULES
USE MOD_Particle_MPI_Vars       ,ONLY: ExchangeProcToGlobalProc,GlobalProcToExchangeProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(ExchangeProcToGlobalProc)
SDEALLOCATE(GlobalProcToExchangeProc)

END SUBROUTINE FinalizePartExchangeProcs


PURE FUNCTION HaloBoxInProc(CartNodes,CartProc,halo_eps)
!===================================================================================================================================
! check if bounding box is on proc
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: CartNodes(6)
REAL,INTENT(IN)       :: CartProc( 6)
REAL,INTENT(IN)       :: halo_eps
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL               :: HaloBoxInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iNode
REAL,DIMENSION(1:3,8) :: xCoordsProc,xCordsTest
!===================================================================================================================================

HaloBoxInProc = .FALSE.

! Reconstruct the eight corner nodes of the cuboid
xCoordsProc(1:3,1) = (/CartProc(1), CartProc(3), CartProc(5)/)
xCoordsProc(1:3,2) = (/CartProc(2), CartProc(3), CartProc(5)/)
xCoordsProc(1:3,3) = (/CartProc(2), CartProc(4), CartProc(5)/)
xCoordsProc(1:3,4) = (/CartProc(1), CartProc(4), CartProc(5)/)
xCoordsProc(1:3,5) = (/CartProc(1), CartProc(3), CartProc(6)/)
xCoordsProc(1:3,6) = (/CartProc(2), CartProc(3), CartProc(6)/)
xCoordsProc(1:3,7) = (/CartProc(2), CartProc(4), CartProc(6)/)
xCoordsProc(1:3,8) = (/CartProc(1), CartProc(4), CartProc(6)/)

xCordsTest(1:3,1) = (/CartNodes(1), CartNodes(3), CartNodes(5)/)
xCordsTest(1:3,2) = (/CartNodes(2), CartNodes(3), CartNodes(5)/)
xCordsTest(1:3,3) = (/CartNodes(2), CartNodes(4), CartNodes(5)/)
xCordsTest(1:3,4) = (/CartNodes(1), CartNodes(4), CartNodes(5)/)
xCordsTest(1:3,5) = (/CartNodes(1), CartNodes(3), CartNodes(6)/)
xCordsTest(1:3,6) = (/CartNodes(2), CartNodes(3), CartNodes(6)/)
xCordsTest(1:3,7) = (/CartNodes(2), CartNodes(4), CartNodes(6)/)
xCordsTest(1:3,8) = (/CartNodes(1), CartNodes(4), CartNodes(6)/)

! Check if any of the eight test corner nodes is within the current proc
DO iNode = 1,8
  IF (   ((xCordsTest (1,iNode).LE.CartProc (2)+halo_eps).AND.(xCordsTest (1,iNode).GE.CartProc (1)-halo_eps))  &
    .AND.((xCordsTest (2,iNode).LE.CartProc (4)+halo_eps).AND.(xCordsTest (2,iNode).GE.CartProc (3)-halo_eps))  &
    .AND.((xCordsTest (3,iNode).LE.CartProc (6)+halo_eps).AND.(xCordsTest (3,iNode).GE.CartProc (5)-halo_eps))) THEN
    HaloBoxInProc = .TRUE.
  END IF
END DO

! Reverse check if any of the proc corner nodes is within the test proc
DO iNode = 1,8
  IF (   ((xCoordsProc(1,iNode).LE.CartNodes(2)+halo_eps).AND.(xCoordsProc(1,iNode).GE.CartNodes(1)-halo_eps))  &
    .AND.((xCoordsProc(2,iNode).LE.CartNodes(4)+halo_eps).AND.(xCoordsProc(2,iNode).GE.CartNodes(3)-halo_eps))  &
    .AND.((xCoordsProc(3,iNode).LE.CartNodes(6)+halo_eps).AND.(xCoordsProc(3,iNode).GE.CartNodes(5)-halo_eps))) THEN
    HaloBoxInProc = .TRUE.
  END IF
END DO

END FUNCTION HaloBoxInProc
#endif /*MPI*/

END MODULE MOD_Particle_MPI_Halo
