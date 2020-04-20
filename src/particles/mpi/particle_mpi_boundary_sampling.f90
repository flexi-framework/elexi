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

MODULE MOD_Particle_MPI_Boundary_Sampling
!===================================================================================================================================
! module for MPI communication of particle surface sampling
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE InitSurfCommunication
  MODULE PROCEDURE InitSurfCommunication
END INTERFACE

INTERFACE ExchangeSurfData
  MODULE PROCEDURE ExchangeSurfData
END INTERFACE


PUBLIC :: InitSurfCommunication
PUBLIC :: ExchangeSurfData
!===================================================================================================================================

CONTAINS


SUBROUTINE InitSurfCommunication()
!----------------------------------------------------------------------------------------------------------------------------------!
!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_SURF
!USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeProcessors,myComputeNoderank
USE MOD_Particle_MPI_Shared_Vars,ONLY: myLeaderGroupRank,nLeaderGroupProcs
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPIRankSharedLeader,MPIRankSurfLeader
USE MOD_Particle_MPI_Shared_Vars,ONLY: mySurfRank,nSurfLeaders!,nSurfCommProc
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides,nComputeNodeSurfTotalSides,offsetComputeNodeSurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode,SurfSampSize,nSurfSample
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
!USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2GlobalSide
USE MOD_Particle_MPI_Vars       ,ONLY: SurfSendBuf,SurfRecvBuf
!USE MOD_Particle_Vars           ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: msg_status(1:MPI_STATUS_SIZE)
INTEGER                       :: iProc,color
INTEGER                       :: leadersGroup,LeaderID,surfGroup
INTEGER                       :: iSide
INTEGER                       :: sendbuf,recvbuf
INTEGER                       :: nSendSurfSidesTmp(0:nLeaderGroupProcs-1)
INTEGER                       :: nRecvSurfSidesTmp(0:nLeaderGroupProcs-1)
!INTEGER                       :: nSurfSidesLeader(1:2,0:nLeaderGroupProcs-1)
INTEGER                       :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1)
INTEGER                       :: SendSurfGlobalID(0:nLeaderGroupProcs-1,1:nComputeNodeSurfTotalSides)
!INTEGER                       :: SampSizeAllocate
!===================================================================================================================================
!--- Open receive buffer (number of sampling surfaces in other node's halo region)
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_IRECV( nRecvSurfSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SHARED                                     &
                , RecvRequest(iProc)                                          &
                , IERROR)
END DO

!--- count all surf sides per other compute-node which get sampling data from current leader
nRecvSurfSidesTmp = 0
nSendSurfSidesTmp = 0

DO iSide = 1,nComputeNodeSurfTotalSides
  ! count surf sides per compute node
  LeaderID = SurfSide2GlobalSide(SURF_LEADER,iSide)
  nSendSurfSidesTmp(LeaderID) = nSendSurfSidesTmp(LeaderID) + 1
  SendSurfGlobalID(LeaderID,nSendSurfSidesTmp(LeaderID)) = SurfSide2GlobalSide(SURF_SIDEID,iSide)
END DO

!--- send all other leaders the number of sampling sides coming from current node
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_ISEND( nSendSurfSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SHARED                                     &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_WAIT(SendRequest(iProc),msg_status(:),IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(RecvRequest(iProc),msg_status(:),IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

!--- Split communicator from MPI_COMM_LEADER_SHARED
color = MPI_UNDEFINED
IF (SurfOnNode) color = 1201

! create new SurfMesh communicator for SurfMesh communication. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_LEADERS_SHARED, color, MPI_INFO_NULL, MPI_COMM_LEADERS_SURF, IERROR)

! Find my rank on the shared communicator, comm size and proc name
CALL MPI_COMM_RANK(MPI_COMM_LEADERS_SURF, mySurfRank  , IERROR)
CALL MPI_COMM_SIZE(MPI_COMM_LEADERS_SURF, nSurfLeaders, IERROR)

! Map global rank number into shared rank number. Returns MPI_UNDEFINED if not on the same communicator
ALLOCATE(MPIRankSharedLeader(0:nLeaderGroupProcs-1))
ALLOCATE(MPIRankSurfLeader  (0:nLeaderGroupProcs-1))
DO iProc=0,nLeaderGroupProcs-1
  MPIRankSharedLeader(iProc) = iProc
END DO

! Get handles for each group
CALL MPI_COMM_GROUP(MPI_COMM_LEADERS_SHARED,leadersGroup,IERROR)
CALL MPI_COMM_GROUP(MPI_COMM_LEADERS_SURF  ,surfGroup   ,IERROR)

! Finally translate global rank to local rank
CALL MPI_GROUP_TRANSLATE_RANKS(leadersGroup,nLeaderGroupProcs,MPIRankSharedLeader,surfGroup,MPIRankSurfLeader,IERROR)
SWRITE(UNIT_stdOUt,'(A,I3,A)') ' Starting surface communication between ', nSurfLeaders, ' compute nodes'

!!--- Count all communicated sides and build mapping for other leaders
!ALLOCATE(nSurfSidesLeader(1:2,0:nSurfLeaders-1))
!
!nSurfCommProc = 0
!DO iProc = 0,nLeaderGroupProcs-1
!  ! a leader defines itself as if it has surf sides within its local domain. However, there might be procs which neither send nor
!  ! receive sides from us. We can reduce nSurfLeaders to nSurfCommProc
!  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
!!  IF ((nRecvSurfSidesTmp(iProc).EQ.0) .AND. (nSendSurfSidesTmp(iProc).EQ.0)) CYCLE
!
!  ! MPI ranks, start at 0
!  nSurfSidesLeader(1,nSurfCommProc) = nSendSurfSidesTmp(iProc)
!  nSurfSidesLeader(2,nSurfCommProc) = nRecvSurfSidesTmp(iProc)
!  nSurfCommProc = nSurfCommProc + 1
!END DO

!--- Open receive buffer (mapping from message surface ID to global side ID)
ALLOCATE(SurfMapping(0:nSurfLeaders-1))
DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  ! Save number of send and recv sides
  SurfMapping(iProc)%nRecvSurfSides = nSendSurfSidesTmp(MPIRankSurfLeader(iProc))
  SurfMapping(iProc)%nSendSurfSides = nRecvSurfSidesTmp(MPIRankSurfLeader(iProc))

  ! Only open recv buffer if we are expecting sides from this leader node
  IF (nRecvSurfSidesTmp(MPIRankSurfLeader(iProc)).EQ.0) CYCLE

  ALLOCATE(SurfMapping(iProc)%RecvSurfGlobalID(1:nRecvSurfSidesTmp(MPIRankSurfLeader(iProc))))

  CALL MPI_IRECV( SurfMapping(iProc)%RecvSurfGlobalID                         &
                , nRecvSurfSidesTmp(MPIRankSurfLeader(iProc))                 &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , RecvRequest(iProc)                                          &
                , IERROR)
END DO

DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  ! Only open send buffer if we are expecting sides from this leader node
  IF (nSendSurfSidesTmp(MPIRankSurfLeader(iProc)).EQ.0) CYCLE

  ALLOCATE(SurfMapping(iProc)%SendSurfGlobalID(1:nSendSurfSidesTmp(MPIRankSurfLeader(iProc))))

  SurfMapping(iProc)%SendSurfGlobalID = SendSurfGlobalID(MPIRankSurfLeader(iProc),1:nSendSurfSidesTmp(MPIRankSurfLeader(iProc)))

  CALL MPI_ISEND( SurfMapping(iProc)%SendSurfGlobalID                         &
                , nSendSurfSidesTmp(MPIRankSurfLeader(iProc))                 &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  IF (nSendSurfSidesTmp(MPIRankSurfLeader(iProc)).NE.0) THEN
    CALL MPI_WAIT(SendRequest(iProc),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF

  IF (nRecvSurfSidesTmp(MPIRankSurfLeader(iProc)).NE.0) THEN
    CALL MPI_WAIT(RecvRequest(iProc),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

!--- Allocate send and recv buffer for each surf leader
ALLOCATE(SurfSendBuf(0:nSurfLeaders-1))
ALLOCATE(SurfRecvBuf(0:nSurfLeaders-1))

DO iProc = 0,nSurfLeaders-1
  ! Only allocate send buffer if we are expecting sides from this leader node
  IF (SurfMapping(iProc)%nSendSurfSides.GT.0) THEN
    ALLOCATE(SurfSendBuf(iProc)%content(SurfSampSize*(nSurfSample**2)*SurfMapping(iProc)%nSendSurfSides))
    SurfSendBuf(iProc)%content = 0.
  END IF

  ! Only allocate recv buffer if we are expecting sides from this leader node
  IF (SurfMapping(iProc)%nRecvSurfSides.GT.0) THEN
    ALLOCATE(SurfRecvBuf(iProc)%content(SurfSampSize*(nSurfSample**2)*SurfMapping(iProc)%nRecvSurfSides))
    SurfRecvBuf(iProc)%content = 0.
  END IF
END DO ! iProc


!--- Save number of total surf sides
IF (surfOnNode) THEN
  IF (nSurfLeaders.EQ.1) THEN
    offsetComputeNodeSurfSide = 0
    nSurfTotalSides           = nComputeNodeSurfSides
  ELSE
    sendbuf = nComputeNodeSurfSides
    recvbuf = 0
    CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_LEADERS_SURF,iError)
    offsetComputeNodeSurfSide = recvbuf
    ! last proc knows CN total number of BC elems
    sendbuf = offsetComputeNodeSurfSide + nSurfTotalSides
    CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nSurfLeaders-1,MPI_COMM_LEADERS_SURF,iError)
    nSurfTotalSides = sendbuf
  END IF
END IF


END SUBROUTINE InitSurfCommunication


SUBROUTINE ExchangeSurfData()
!===================================================================================================================================
! exchange the surface data
!> 1) collect the information on the local compute-node
!> 2) compute-node leaders with sampling sides in their halo region and the original node communicate the sampling information
!> 3) compute-node leaders ensure synchronization of shared arrays on their node
!!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSampSize,nSurfSample
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState,SampWallState_Shared,SampWallState_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: nErosionVars
USE MOD_Particle_MPI_Vars       ,ONLY: SurfSendBuf,SurfRecvBuf
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_Particle_MPI_Shared_Vars,ONLY: nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_Particle_Vars           ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iProc,SideID
INTEGER                         :: iPos,p,q
INTEGER                         :: iSpec,nShift
INTEGER                         :: MessageSize,iSurfSide,SurfSideID
INTEGER                         :: nValues
INTEGER                         :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
REAL                            :: SampWallTmp(1:SurfSampSize)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfOnNode) RETURN

! collect the information from the proc-local shadow arrays in the compute-node shared array
MessageSize = SurfSampSize*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides
CALL MPI_REDUCE(SampWallState       ,SampWallState_Shared       ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
CALL MPI_WIN_SYNC(SampWallState_Shared_Win       ,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! now the shadow arrays can be nullified
SampWallState        = 0.

! prepare buffers for surf leader communication
IF (myComputeNodeRank.EQ.0) THEN
  nValues = SurfSampSize*nSurfSample**2

! open receive buffer
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nRecvSurfSides * nValues
    CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , RecvRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! build message
  DO iProc = 0,nSurfLeaders-1
    ! Ignore myself
    IF (iProc .EQ. mySurfRank) CYCLE

    ! Only assemble message if we are expecting sides to send to this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0
    SurfSendBuf(iProc)%content = 0.

    DO iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
      SideID     = SurfMapping(iProc)%SendSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

      ! Assemble message
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
          SurfSendBuf(iProc)%content(iPos+1:iPos+SurfSampSize) = SampWallState_Shared(:,p,q,SurfSideID)
          iPos = iPos + SurfSampSize
        END DO ! p=0,nSurfSample
      END DO ! q=0,nSurfSample

      SampWallState_Shared(:,:,:,SurfSideID)=0.
    END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nSendSurfSides * nValues
    CALL MPI_ISEND( SurfSendBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , SendRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! Finish received number of sampling surfaces
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    IF (SurfMapping(iProc)%nSendSurfSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nRecvSurfSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
      SideID     = SurfMapping(iProc)%RecvSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

      DO q=1,nSurfSample
        DO p=1,nSurfSample
          ! Treat each variable individual, depending on if they can be added
          SampWallTmp(:)                    = SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfSampSize)

          ! All Variables are saved DOUBLE. First Total, then per SPECIES
          !-- 1. - .. / Impact Counter
          SampWallState_Shared(1,p,q,SurfSideID) = SampWallState_Shared(1,p,q,SurfSideID) + SampWallTmp(1)
          !-- 2. - 6. / Kinetic energy on impact (mean, min, max, M2, variance)
          SampWallState_Shared(2,p,q,SurfSideID) = ( SampWallState_Shared(2,p,q,SurfSideID)                                        &
                                                   *(SampWallState_Shared(1,p,q,SurfSideID)                                        &
                                                   - SampWallTmp(1)) + SampWallTmp(2)*SampWallTmp(1))                              &
                                                   / SampWallState_Shared(1,p,q,SurfSideID)
          IF (SampWallTmp(3).LT.SampWallState_Shared(3,p,q,SurfSideID))                                                            &
            SampWallState_Shared(3,p,q,SurfSideID) = SampWallTmp(3)
          IF (SampWallTmp(4).GT.SampWallState_Shared(4,p,q,SurfSideID))                                                            &
            SampWallState_Shared(4,p,q,SurfSideID) = SampWallTmp(4)
          !-- 5-6 are M2 and variance. Nothing to be done about those as we need data after each particle impact
          !>>
          !-- 7. - 11 / Impact angle (mean, min, max, M2, variance)
          SampWallState_Shared(7,p,q,SurfSideID) = ( SampWallState_Shared(7,p,q,SurfSideID)                                        &
                                                   *(SampWallState_Shared(1,p,q,SurfSideID)                                        &
                                                   - SampWallTmp(1)) + SampWallTmp(7)*SampWallTmp(1))                              &
                                                   / SampWallState_Shared(1,p,q,SurfSideID)
          IF (SampWallTmp(8).LT.SampWallState_Shared(8,p,q,SurfSideID))                                                            &
            SampWallState_Shared(3,p,q,SurfSideID) = SampWallTmp(8)
          IF (SampWallTmp(9).GT.SampWallState_Shared(9,p,q,SurfSideID))                                                            &
            SampWallState_Shared(4,p,q,SurfSideID) = SampWallTmp(9)
          !-- 10-11 are M2 and variance. Nothing to be done about those as we need data after each particle impact
          !>>
          !-- 12 - 14 / Sampling Current Forces at walls - we can simply add those
          SampWallState_Shared(12:14,p,q,SurfSideID) = SampWallState_Shared(12:14,p,q,SurfSideID) + SampWallTmp(12:14)
          !-- 15 - 17 / Sampling Average Forces at walls  - we can simply add those
          SampWallState_Shared(15:17,p,q,SurfSideID) = SampWallState_Shared(15:17,p,q,SurfSideID) + SampWallTmp(15:17)
          !<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          IF (nSpecies.GT.1) THEN
            DO iSpec = 1,nSpecies
              nShift = iSpec * nErosionVars
              !-- 1. - .. / Impact Counter
              SampWallState_Shared(1+nShift,p,q,SurfSideID) = SampWallState_Shared(1+nShift,p,q,SurfSideID) + SampWallTmp(1+nShift)
              !-- 2. - 6. / Kinetic energy on impact (mean, min, max, M2, variance)
              SampWallState_Shared(2+nShift,p,q,SurfSideID) = ( SampWallState_Shared(2+nShift,p,q,SurfSideID)                      &
                                                              *(SampWallState_Shared(1+nShift,p,q,SurfSideID)                      &
                                                              - SampWallTmp(1+nShift))                                             &
                                                              + SampWallTmp(2+nShift)*SampWallTmp(1+nShift))                       &
                                                              / SampWallState_Shared(1+nShift,p,q,SurfSideID)
              IF (SampWallTmp(3+nShift).LT.SampWallState_Shared(3+nShift,p,q,SurfSideID))                                          &
                  SampWallState_Shared(3+nShift,p,q,SurfSideID) = SampWallTmp(3+nShift)
              IF (SampWallTmp(4+nShift).GT.SampWallState_Shared(4+nShift,p,q,SurfSideID))                                          &
                  SampWallState_Shared(4+nShift,p,q,SurfSideID) = SampWallTmp(4+nShift)
              !-- 5-6 are M2 and variance. Nothing to be done about those as we need data after each particle impact
              !>>
              !-- 7. - 11 / Impact angle (mean, min, max, M2, variance)
              SampWallState_Shared(7+nShift,p,q,SurfSideID) = (SampWallState_Shared(7+nShift,p,q,SurfSideID)                       &
                                                             *(SampWallState_Shared(1+nShift,p,q,SurfSideID)                       &
                                                             - SampWallTmp(1+nShift))                                              &
                                                             + SampWallTmp(7+nShift)*SampWallTmp(1+nShift))                        &
                                                             / SampWallState_Shared(1+nShift,p,q,SurfSideID)
              IF (SampWallTmp(8+nShift).LT.SampWallState_Shared(8+nShift,p,q,SurfSideID))                                          &
                  SampWallState_Shared(8+nShift,p,q,SurfSideID) = SampWallTmp(8+nShift)
              IF (SampWallTmp(9+nShift).GT.SampWallState_Shared(9+nShift,p,q,SurfSideID))                                          &
                  SampWallState_Shared(9+nShift,p,q,SurfSideID) = SampWallTmp(9+nShift)
              !-- 10-11 are M2 and variance. Nothing to be done about those as we need data after each particle impact
              !>>
              !-- 12 - 14 / Sampling Current Forces at walls - we can simply add those
              SampWallState_Shared(12+nShift:14+nShift,p,q,SurfSideID) = SampWallState_Shared(12+nShift:14+nShift,p,q,SurfSideID)  &
                                                                       + SampWallTmp         (12+nShift:14+nShift)
              !-- 15 - 17 / Sampling Average Forces at walls  - we can simply add those
              SampWallState_Shared(15+nShift:17+nShift,p,q,SurfSideID) = SampWallState_Shared(15+nShift:17+nShift,p,q,SurfSideID)  &
                                                                       + SampWallTmp         (15+nShift:17+nShift)
            END DO ! iSpec = 1,nSpecies
          END IF ! nSpecies.GT.1
          iPos = iPos + SurfSampSize
        END DO ! p = 0,nSurfSample
      END DO ! q = 0,nSurfSample
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides

     ! Nullify buffer
    SurfRecvBuf(iProc)%content = 0.
  END DO ! iProc
END IF

! ensure synchronization on compute node
CALL MPI_WIN_SYNC(SampWallState_Shared_Win       ,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

END SUBROUTINE ExchangeSurfData
#endif /*USE_MPI*/


END MODULE MOD_Particle_MPI_Boundary_Sampling
