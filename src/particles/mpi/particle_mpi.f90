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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_MPI
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleMPI
  MODULE PROCEDURE InitParticleMPI
END INTERFACE

#if USE_MPI
INTERFACE IRecvNbOfParticles
  MODULE PROCEDURE IRecvNbOfParticles
END INTERFACE

INTERFACE SendNbOfParticles
  MODULE PROCEDURE SendNbOfParticles
END INTERFACE

INTERFACE FinalizeParticleMPI
  MODULE PROCEDURE FinalizeParticleMPI
END INTERFACE

INTERFACE MPIParticleSend
  MODULE PROCEDURE MPIParticleSend
END INTERFACE

INTERFACE MPIParticleRecv
  MODULE PROCEDURE MPIParticleRecv
END INTERFACE

INTERFACE InitHaloMesh
  MODULE PROCEDURE InitHaloMesh
END INTERFACE

INTERFACE InitParticleCommSize
  MODULE PROCEDURE InitParticleCommSize
END INTERFACE

INTERFACE InitEmissionComm
  MODULE PROCEDURE InitEmissionComm
END INTERFACE

INTERFACE ExchangeBezierControlPoints3D
  MODULE PROCEDURE ExchangeBezierControlPoints3D
END INTERFACE

PUBLIC :: InitParticleMPI
PUBLIC :: InitHaloMesh
PUBLIC :: InitParticleCommSize
PUBLIC :: InitEmissionComm
PUBLIC :: FinalizeParticleMPI
PUBLIC :: ExchangeBezierControlPoints3D
PUBLIC :: IRecvNbOfParticles
PUBLIC :: SendNbofParticles
PUBLIC :: MPIParticleSend
PUBLIC :: MPIParticleRecv
#else
PUBLIC :: InitParticleMPI
#endif /*MPI*/
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
#if USE_MPI
USE MOD_MPI_Vars,               ONLY: nNbProcs
#endif
USE MOD_Particle_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                         :: color
#endif
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI ... '
IF(ParticleMPIInitIsDone) &
  CALL abort(&
    __STAMP__&
  ,' Particle MPI already initialized!')

#if USE_MPI
PartMPI%myRank = myRank
color = 999
! Get particle MPI communicator and size
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,PartMPI%MyRank,PartMPI%COMM,iERROR)
CALL MPI_COMM_SIZE (PartMPI%COMM,PartMPI%nProcs ,iError)

! Check if particle COMM has the same size as DG comm
IF(PartMPI%nProcs.NE.nProcessors) CALL abort(&
    __STAMP__&
    ,' MPI Communicater-size does not match!', IERROR)
PartCommSize   = 0

! Check if this proc is root on PartMPI
IF(PartMPI%MyRank.EQ.0) THEN
  PartMPI%MPIRoot=.TRUE.
ELSE
  PartMPI%MPIRoot=.FALSE.
END IF
iMessage=0

ALLOCATE(SendRequest_Flux(nNbProcs)  )
ALLOCATE(RecRequest_Flux(nNbProcs))
#else
PartMPI%myRank = 0
PartMPI%nProcs = 1
PartMPI%MPIRoot=.TRUE.
#endif  /*MPI*/

ParticleMPIInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMPI


#if USE_MPI
SUBROUTINE InitParticleCommSize()
!===================================================================================================================================
! get size of Particle-MPI-Message. Unfortunately, this subroutine have to be called after particle_init because
! all required features have to be read from the ini-File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Erosion_Vars,    ONLY:PartTrackReflection
USE MOD_Particle_MPI_Vars
USE MOD_Particle_SGS_Vars,        ONLY:nSGSVars,SGSinUse
USE MOD_Particle_Tracking_Vars,   ONLY:DoRefMapping
USE MOD_Particle_Vars,            ONLY:PDM
#if USE_RW
USE MOD_Particle_RandomWalk_Vars, ONLY:nRWVars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: ALLOCSTAT
!===================================================================================================================================

PartCommSize   = 0
! PartState: position and velocity
PartCommSize   = PartCommSize + 6
! TurbPartState: SGS turbulent velocity and random draw
PartCommSize   = PartCommSize + nSGSVars
#if USE_RW
! TurbPartState: RW turbulent velocity, interaction time and random draw
PartCommSize   = PartCommSize + nRWVars
#endif
! Tracking: Include Reference coordinates
IF(DoRefMapping) PartCommSize=PartCommSize+3
! Species-ID
PartCommSize   = PartCommSize + 1
! id of element
PartCommSize   = PartCommSize + 1

! time integration
! communication after each Runge-Kutta stage, so send time derivative must be communicated to the new proc
! Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
PartCommSize   = PartCommSize + 6
! TurbPt_tmp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
IF (SGSinUse)    PartCommSize = PartCommSize + 3
! IsNewPart for RK-Reconstruction
PartCommSize   = PartCommSize + 1

! Include reflection counter
IF(PartTrackReflection) PartCommSize   = PartCommSize + 1

! additional stuff for full RK schemes, e.g. implicit and imex RK
! if iStage=0, then the PartStateN is not communicated
PartCommSize0  = PartCommSize

ALLOCATE( PartMPIExchange%nPartsSend(2,PartMPI%nMPINeighbors)  &
        , PartMPIExchange%nPartsRecv(2,PartMPI%nMPINeighbors)  &
        , PartRecvBuf(1:PartMPI%nMPINeighbors)                 &
        , PartSendBuf(1:PartMPI%nMPINeighbors)                 &
        , PartMPIExchange%SendRequest(2,PartMPI%nMPINeighbors) &
        , PartMPIExchange%RecvRequest(2,PartMPI%nMPINeighbors) &
        , PartTargetProc(1:PDM%MaxParticleNumber)              &
        , STAT=ALLOCSTAT                                       )

IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)

PartMPIExchange%nPartsSend=0
PartMPIExchange%nPartsRecv=0

END SUBROUTINE InitParticleCommSize


SUBROUTINE IRecvNbOfParticles()
!===================================================================================================================================
! Open Recv-Buffer for number of received particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartMPIExchange
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iProc
!===================================================================================================================================

! Asynchronous communication. Open receive buffer to all neighboring procs to get the number of particles THEY want to send
PartMPIExchange%nPartsRecv=0
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_IRECV( PartMPIExchange%nPartsRecv(:,iProc)                        &
                , 2                                                          &
                , MPI_INTEGER                                                &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(1,iProc)                       &
                , IERROR )
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


SUBROUTINE SendNbOfParticles(doParticle_In)
!===================================================================================================================================
! this routine sends the number of send particles. Following steps are performed
! 1) Compute number of Send Particles
! 2) Perform MPI_ISEND with number of particles
! Rest is performed in SendParticles
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartHaloElemToProc, PartTargetProc
USE MOD_Particle_Vars,            ONLY:PEM,PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: doParticle(1:PDM%ParticleVecLength)
INTEGER                       :: iPart,ElemID,iProc
!===================================================================================================================================

IF(PRESENT(DoParticle_IN))THEN
  DoParticle=PDM%ParticleInside(1:PDM%ParticleVecLength).AND.DoParticle_In
ELSE
  DoParticle(1:PDM%ParticleVecLength)=PDM%ParticleInside(1:PDM%ParticleVecLength)
END IF

! 1) get number of send particles
!--- Count number of particles in cells in the halo region and add them to the message
!--- CAUTION: using local indices for halo elem -> proc association with PartHaloElemToProc. PartMPI%Neighbor contains the inverse
!--- mapping
PartMPIExchange%nPartsSend=0
PartTargetProc=-1
DO iPart=1,PDM%ParticleVecLength
  !IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
  IF(.NOT.DoParticle(iPart)) CYCLE
  ElemID=PEM%Element(iPart)
  IF(ElemID.GT.PP_nElems) THEN
    PartMPIExchange%nPartsSend(1,PartHaloElemToProc(LOCAL_PROC_ID,ElemID))=             &
                                        PartMPIExchange%nPartsSend(1,PartHaloElemToProc(LOCAL_PROC_ID,ElemID))+1
    PartTargetProc(iPart)=PartHaloElemToProc(LOCAL_PROC_ID,ElemID)
  END IF
END DO ! iPart


! 2) send number of send particles
!--- Loop over all neighboring procs. Map local proc ID to global through PartMPI%Neighbor.
!--- Asynchronous communication, just send here and check for success later.
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_ISEND( PartMPIExchange%nPartsSend(:,iProc)                      &
                , 2                                                          &
                , MPI_INTEGER                                                &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(1,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE SendNbOfParticles


SUBROUTINE MPIParticleSend()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! first steps are perforemd in SendNbOfParticles
! 1) Compute number of Send Particles
! 2) Performe MPI_ISEND with number of particles
! Starting Here:
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
! CAUTION: If particles are send for deposition, PartTargetProc has the information, if a particle is send
!          and after the buld and wait for number of particles reused to build array with external parts
!          informations in PartState,.. can be reusded, because they are not overwritten
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartHaloElemToProc,PartCommSize,PartSendBuf, PartRecvBuf &
                                      ,PartTargetProc
USE MOD_Particle_Vars,            ONLY:PartSpecies,PEM,PDM,PartPosRef
USE MOD_Particle_Vars,            ONLY:PartState,Pt_temp
USE MOD_Particle_Vars,            ONLY:TurbPartState,TurbPt_temp
! Variables for erosion tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Erosion_Vars
! Variables for SGS model
USE MOD_Particle_SGS_Vars,        ONLY:SGSinUse,nSGSVars
#if USE_RW
! Variables for RW model
USE MOD_Particle_RandomWalk_Vars, ONLY:nRWVars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID,iPos,iProc,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize, nRecvParticles, nSendParticles
INTEGER                       :: ALLOCSTAT
!===================================================================================================================================

! 3) Build Message
DO iProc=1, PartMPI%nMPINeighbors
  ! find number of particles to send
  nSendParticles=PartMPIExchange%nPartsSend(1,iProc)
  iPos=0

  ! only build message if we have particles to send
  IF(nSendParticles.EQ.0) CYCLE

  ! allocate SendBuff of required size
  MessageSize=nSendParticles*PartCommSize
  ALLOCATE(PartSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) &
    CALL abort(__STAMP__,'  Cannot allocate PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))

  ! fill message
  !--- Loop over all particles
  DO iPart=1,PDM%ParticleVecLength
    IF(nSendParticles.EQ.0) EXIT

    ! find particles to send to the current neighboring proc
    IF(PartTargetProc(iPart).EQ.iProc) THEN
      ! fill content
      ElemID=PEM%Element(iPart)

      ! position and velocity in physical space
      PartSendBuf(iProc)%content(1+iPos:6+iPos) = PartState(1:6,iPart)
      jpos=iPos+6

      ! SGS turbulent velocity and random draw
      PartSendBuf(iProc)%content(1+jPos:nSGSVars+jPos) = TurbPartState(1:nSGSVars,iPart)
      jpos=jpos+nSGSVars
#if USE_RW
      ! RW turbulent velocity, interaction time and random draw
      PartSendBuf(iProc)%content(1+jPos:nRWVars+jPos) = TurbPartState(1:nRWVars,iPart)
      jpos=jpos+nRWVars
#endif

      ! position in reference space (if required)
      IF(DoRefMapping) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = PartPosRef(1:3,iPart)
        jPos=jPos+3
      END IF

      ! reflection counter
      IF (PartTrackReflection) THEN
        PartSendBuf(iProc)%content(1+jPos) = REAL(PartReflCount(iPart),KIND=8)
        jPos=jPos+1
      END IF

      ! particles species
      PartSendBuf(iProc)%content(       1+jPos) = REAL(PartSpecies(iPart),KIND=8)
      jPos=jPos+1

      ! Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
      PartSendBuf(iProc)%content(1+jPos:6+jPos) = Pt_temp(1:6,iPart)
      jPos=jPos+6

      ! TurbPt_tmp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
      IF (SGSinUse) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = TurbPt_temp(1:3,iPart)
        jpos=jpos+3
      END IF

      ! IsNewPart for RK-Reconstruction
      IF (PDM%IsNewPart(iPart)) THEN
        PartSendBuf(iProc)%content(1+jPos) = 1.
      ELSE
        PartSendBuf(iProc)%content(1+jPos) = 0.
      END IF
      jPos=jPos+1

      ! native ElemID of particle position on receiving proc
      PartSendBuf(iProc)%content(  1+jPos) = REAL(PartHaloElemToProc(NATIVE_ELEM_ID,ElemID),KIND=8)
      jPos=jPos+1

      ! PartCommSize must be a multiple of particles to send
      IF(MOD(jPos,PartCommSize).NE.0) THEN
        IPWRITE(UNIT_stdOut,*)  'PartCommSize',PartCommSize
        IPWRITE(UNIT_stdOut,*)  'jPos',jPos
        CALL Abort(__STAMP__,' Particle-wrong sending message size!')
      END IF

      ! Move index counter to beginning of new part variables
      iPos=iPos+PartCommSize

      ! particle is ready for send, now it can deleted
      PDM%ParticleInside(iPart) = .FALSE.
    END IF ! Particle is particle with target proc-id equals local proc id
  END DO  ! iPart

  IF(iPos.NE.(MessageSize)) IPWRITE(*,*) ' error message size', iPos,MessageSize
END DO ! iProc

! 4) Finish Received number of particles
!--- Wait for all neighboring procs to acknowlege both our send and our recv request. Then every neighbor proc knows the number of
!--- particles to communicate
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_WAIT(PartMPIExchange%SendRequest(1,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

! total number of received particles
PartMPIExchange%nMPIParticles=SUM(PartMPIExchange%nPartsRecv(1,:))


DO iPart=1,PDM%ParticleVecLength
  IF(PartTargetProc(iPart).EQ.-1) CYCLE
  PartState(1:6,iPart) = 0.
  PartSpecies(iPart)   = 0
  Pt_temp(1:6,iPart)   = 0.
END DO ! iPart=1,PDM%ParticleVecLength


! 5) Allocate received buffer and open MPI_IRECV
DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE

  ! count number of particles from each proc and determine size of message
  nRecvParticles = PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize    = nRecvParticles * PartCommSize

  ! allocate recv buffer with the correct size
  ALLOCATE(PartRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    IPWRITE(*,*) 'sum of total received particles            ', SUM(PartMPIExchange%nPartsRecv(1,:))
    IPWRITE(*,*) 'sum of total received deposition particles ', SUM(PartMPIExchange%nPartsRecv(2,:))
    CALL abort(__STAMP__,'  Cannot allocate PartRecvBuf, local source ProcId, Allocstat',iProc,REAL(ALLOCSTAT))
  END IF

  ! start asynchronous MPI recv for each proc
  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) &
    CALL abort(__STAMP__,' MPI Communication error', IERROR)

END DO ! iProc

! 6) Send Particles
DO iProc=1,PartMPI%nMPINeighbors
  ! only consider procs where we have particles to send
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE

  ! count number of particles for each proc and determine size of message
  nSendParticles = PartMPIExchange%nPartsSend(1,iProc)
  MessageSize    = nSendParticles*PartCommSize

  ! start asychronous MPI send for each proc
  CALL MPI_ISEND( PartSendBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) &
    CALL abort(__STAMP__,' MPI Communication error', IERROR)

END DO ! iProc

END SUBROUTINE MPIParticleSend


SUBROUTINE MPIParticleRecv()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! 1) Compute number of Send Particles
! 2) Performe MPI_ISEND with number of particles
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartCommSize, PartRecvBuf,PartSendBuf
USE MOD_Particle_Vars,            ONLY:PartSpecies,PEM,PDM,PartPosRef
USE MOD_Particle_Vars,            ONLY:PartState,Pt_temp
USE MOD_Particle_Vars,            ONLY:TurbPartState,TurbPt_temp
USE MOD_Particle_Tracking_Vars,   ONLY:DoRefMapping
! variables for erosion tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Erosion_Vars
! Variables for SGS model
USE MOD_Particle_SGS_Vars,        ONLY:SGSinUse,nSGSVars
#if USE_RW
! Variables for RW model
USE MOD_Particle_RandomWalk_Vars, ONLY:nRWVars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc, iPos, nRecv, PartID,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize, nRecvParticles
!===================================================================================================================================

! wait for all neighboring procs to acknowledge our MPI send
DO iProc=1,PartMPI%nMPINeighbors
  ! ignore procs we did not send anything
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE

  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) &
    CALL abort(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc


nRecv=0
! decompose the recv message from each proc and set local variables accordingly
DO iProc=1,PartMPI%nMPINeighbors
  ! ignore procs we did not receive anything from
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE

  ! count number of recv particles and determine message sage
  nRecvParticles = PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize    = nRecvParticles*PartCommSize

  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
  ! correct loop shape
  ! DO iPart=1,nRecvParticles
  ! nParts 1 Pos=1..17
  ! nPart2 2 Pos=1..17,18..34
  DO iPos=0,MessageSize-1,PartCommSize
    IF(nRecvParticles.EQ.0) EXIT

    ! increment counter of received particles
    nRecv  = nRecv+1

    ! particles get a local ID on each proc, therefore put it at the next free position
    PartID = PDM%nextFreePosition(nRecv+PDM%CurrentNextFreePosition)
    IF(PartID.EQ.0) &
      CALL abort(__STAMP__,' Error in ParticleExchange_parallel. Corrupted list: PIC%nextFreePosition', nRecv)

    ! position and velocity in physical space
    PartState(1:6,PartID)   = PartRecvBuf(iProc)%content(1+iPos:6+iPos)
    jpos=iPos+6

    TurbPartState(1:nSGSVars,PartID) = PartRecvBuf(iProc)%content(1+jpos:nSGSVars+jpos)
    jpos=jpos+nSGSVars
#if USE_RW
      ! RW turbulent velocity, interaction time and random draw
      TurbPartState(1:nRWVars,PartID) = PartRecvBuf(iProc)%content(1+jPos:nRWVars+jPos)
      jpos=jpos+nRWVars
#endif

    ! position in reference space (if required)
    IF(DoRefMapping) THEN
      PartPosRef(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos:3+jPos)
      jPos=jPos+3
    END IF

    ! reflection counter
    IF (PartTrackReflection) THEN
      PartReflCount(PartID)   = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
      jpos=jpos+1
    END IF

    ! particles species
    PartSpecies(PartID)     = INT(PartRecvBuf(iProc)%content( 1+jPos),KIND=4)
    jPos=jPos+1

    ! Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
    Pt_temp(1:6,PartID)     = PartRecvBuf(iProc)%content( 1+jPos:6+jPos)
    jpos=jpos+6

    ! TurbPt_tmp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
    IF (SGSinUse) THEN
      TurbPt_temp(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos:3+jPos)
      jpos=jpos+3
    END IF

    ! IsNewPart for RK-Reconstruction
    IF      ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 1) THEN
      PDM%IsNewPart(PartID)=.TRUE.
    ELSE IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 0) THEN
      PDM%IsNewPart(PartID)=.FALSE.
    ELSE
      CALL Abort(__STAMP__,'Error with IsNewPart in MPIParticleRecv!')
    END IF
    jPos=jPos+1

    ! native ElemID of particle position on my proc
    PEM%Element(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1

    ! PartCommSize must be a multiple of particles to receive
    IF(MOD(jPos,PartCommSize).NE.0)THEN
      IPWRITE(UNIT_stdOut,*)  'jPos',jPos
      CALL Abort(__STAMP__,' Particle-wrong receiving message size!')
    END IF

    ! Set Flag for received parts in order to localize them later
    PDM%ParticleInside(PartID) = .TRUE.
    PEM%lastElement(PartID)    = -888
  END DO
END DO ! iProc


PDM%ParticleVecLength       = PDM%ParticleVecLength + PartMPIExchange%nMPIParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles
IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) &
  CALL abort(__STAMP__,' ParticleVecLegnth>MaxParticleNumber due to MPI-communication!')

! deallocate send,receive buffer
DO iProc=1,PartMPI%nMPINeighbors
  SDEALLOCATE(PartRecvBuf(iProc)%content)
  SDEALLOCATE(PartSendBuf(iProc)%content)
END DO ! iProc

! finally nullify the send and recv counters
PartMPIExchange%nPartsRecv=0
PartMPIExchange%nPartsSend=0


END SUBROUTINE MPIParticleRecv


SUBROUTINE FinalizeParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_MPI_Vars
USE MOD_Particle_Vars,            ONLY:Species,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nInitRegions,iInitRegions,iSpec
!===================================================================================================================================

nInitRegions=0
DO iSpec=1,nSpecies
  nInitRegions=nInitRegions+Species(iSpec)%NumberOfInits+(1-Species(iSpec)%StartnumberOfInits)
END DO ! iSpec
IF(nInitRegions.GT.0) THEN
  DO iInitRegions=1,nInitRegions
    IF(PartMPI%InitGroup(iInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      CALL MPI_COMM_FREE(PartMPI%InitGroup(iInitRegions)%Comm,iERROR)
    END IF
  END DO ! iInitRegions
END IF
CALL MPI_COMM_FREE(PartMPI%COMM,iERROR)


SDEALLOCATE( PartHaloElemToProc)
SDEALLOCATE( PartMPI%isMPINeighbor)
SDEALLOCATE( PartMPI%MPINeighbor )
SDEALLOCATE( PartMPI%GlobalToLocal )
SDEALLOCATE( PartMPIExchange%nPartsSend)
SDEALLOCATE( PartMPIExchange%nPartsRecv)
SDEALLOCATE( PartMPIExchange%RecvRequest)
SDEALLOCATE( PartMPIExchange%SendRequest)
SDEALLOCATE( PartMPIExchange%Send_message)
SDEALLOCATE( PartMPI%isMPINeighbor)
SDEALLOCATE( PartMPI%MPINeighbor)
SDEALLOCATE( PartMPI%InitGroup)
SDEALLOCATE( PartSendBuf)
SDEALLOCATE( PartRecvBuf)

! and for communication
SDEALLOCATE( PartTargetProc )

ParticleMPIInitIsDone=.FALSE.
END SUBROUTINE FinalizeParticleMPI


SUBROUTINE ExchangeBezierControlPoints3D()
!===================================================================================================================================
! exchange all beziercontrolpoints at MPI interfaces
! maybe extended to periodic sides, to be tested
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_Particle_MPI_Vars
USE MOD_Mesh_Vars,                  ONLY:NGeo,nSides,firstMPISide_YOUR&
                                        ,firstMortarMPISide,lastMortarMPISide,lastMPISide_YOUR,SideToElem
USE MOD_Particle_Mesh_Vars,         ONLY:NGeoElevated,MortarSlave2MasterInfo
USE MOD_Particle_Surfaces,          ONLY:GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_vars,     ONLY:BezierControlPoints3D,SideSlabIntervals,BezierControlPoints3DElevated &
                                        ,SideSlabIntervals,SideSlabNormals,BoundingBoxIsEmpty
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 ::BezierSideSize,SendID, iSide,SideID,ElemID
!===================================================================================================================================

! funny: should not be required, as sides are built for master and slave sides??
! communicate the MPI Master Sides to Slaves
! all processes have now filled sides and can compute the particles inside the proc region

SendID=1
BezierSideSize=3*(NGeo+1)*(NGeo+1)
DO iNbProc=1,nNbProcs
  ! Start receive face data
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =BezierSideSize*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(BezierControlPoints3D(:,:,:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,RecRequest_Flux(iNbProc),iError)
  END IF
  ! Start send face data
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =BezierSideSize*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(BezierControlPoints3D(:,:,:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,SendRequest_Flux(iNbProc),iError)
  END IF
END DO !iProc=1,nNBProcs

DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0) CALL MPI_WAIT(RecRequest_Flux(iNbProc) ,MPIStatus,iError)
  IF(iERROR.NE.0) CALL abort(&
  __STAMP__&
  ,' MPI-Error during BezierControlPoint-exchange. iError', iERROR)
END DO !iProc=1,nNBProcs
! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0) CALL MPI_WAIT(SendRequest_Flux(iNbProc),MPIStatus,iError)
  IF(iERROR.NE.0) CALL abort(&
  __STAMP__&
  ,' MPI-Error during BezierControlPoint-exchange. iError', iERROR)
END DO !iProc=1,nNBProcs

! build the bounding box for YOUR-MPI-sides without mortar sides
DO iSide=firstMPISide_YOUR,lastMPISide_YOUR
  ! elevation occurs within this routine
  CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)                         &
                                     ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide) &
                                     ,SideSlabNormals(1:3,1:3,iSide)                                         &
                                     ,SideSlabInterVals(1:6,iSide)                                           &
                                     ,BoundingBoxIsEmpty(iSide)                                              )
END DO

! build the bounding box for missing MPI-mortar sides, or YOUR mortar sides
! actually, I do not know, if this is requried
DO iSide=firstMortarMPISide,lastMortarMPISide
  ElemID=SideToElem(S2E_ELEM_ID,iSide)
  SideID=MortarSlave2MasterInfo(iSide)
  IF(ElemID.NE.-1) CYCLE
  IF(SideID.EQ.-1) CYCLE
  ! elevation occurs within this routine
  CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)                         &
                                     ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide) &
                                     ,SideSlabNormals(1:3,1:3,iSide)                                         &
                                     ,SideSlabInterVals(1:6,iSide)                                           &
                                     ,BoundingBoxIsEmpty(iSide)                                              )
END DO

DO iSide=1,nSides
  IF(SUM(ABS(SideSlabIntervals(:,iSide))).EQ.0)THEN
    CALL abort(&
    __STAMP__&
    ,'  Zero bounding box found!, iSide',iSide)
  END IF
END DO

END SUBROUTINE ExchangeBezierControlPoints3D


SUBROUTINE InitHaloMesh()
!===================================================================================================================================
! communicate all direct neighbor sides from master to slave
! has to be called after GetSideType and MPI_INIT of DG solver
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_MPI_Vars
USE MOD_PreProc
USE MOD_Particle_MPI_Vars,          ONLY:PartMPI,PartHaloElemToProc,printMPINeighborWarnings
USE MOD_Particle_MPI_Halo,          ONLY:IdentifyHaloMPINeighborhood,ExchangeHaloGeometry
USE MOD_Particle_Mesh_Vars,         ONLY:nTotalElems,nPartSides
#if CODE_ANALYZE
USE MOD_Particle_Tracking_vars,     ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,         ONLY:nTotalSides,nTotalBCSides
#endif
#if USE_MPI_SHARED
USE MOD_MPI_Shared_Vars,            ONLY:MPIRankShared
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 ::iElem
INTEGER                 ::iProc,ALLOCSTAT,iMPINeighbor
LOGICAL                 ::TmpNeigh
INTEGER,ALLOCATABLE     ::SideIndex(:),ElemIndex(:)
!===================================================================================================================================

ALLOCATE(SideIndex(1:nPartSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'  Cannot allocate SideIndex!')
SideIndex=0
ALLOCATE(ElemIndex(1:PP_nElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'  Cannot allocate ElemIndex!')
ElemIndex=0

! check epsilondistance
DO iProc=0,PartMPI%nProcs-1
  IF(iProc.EQ.PartMPI%MyRank) CYCLE
#if USE_MPI_SHARED
  IF(MPIRankShared(iProc).NE.MPI_UNDEFINED) CYCLE
#endif
  LOGWRITE(*,*)'  - Identify non-immediate MPI-Neighborhood...'
  !--- AS: identifies which of my node have to be sent to iProc w.r.t. to
  !        eps vicinity region.
  CALL IdentifyHaloMPINeighborhood(iProc,SideIndex,ElemIndex)
  LOGWRITE(*,*)'    ...Done'

  LOGWRITE(*,*)'  - Exchange Geometry of MPI-Neighborhood...'
  CALL ExchangeHaloGeometry(iProc,ElemIndex)
  LOGWRITE(*,*)'    ...Done'
  SideIndex(:)=0
  ElemIndex(:)=0
END DO
DEALLOCATE(SideIndex,STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
,'Could not deallocate SideIndex')
END IF

#if CODE_ANALYZE
IF(DoRefMapping) CALL CheckArrays(nTotalSides,nTotalElems,nTotalBCSides)
#endif /*CODE_ANALYZE*/

! Make sure PMPIVAR%MPINeighbor is consistent
DO iProc=0,PartMPI%nProcs-1
  IF (PartMPI%MyRank.EQ.iProc) CYCLE
#if USE_MPI_SHARED
  IF(MPIRankShared(iProc).NE.MPI_UNDEFINED) CYCLE
#endif
  IF (PartMPI%MyRank.LT.iProc) THEN
    CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
    CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  ELSE IF (PartMPI%MyRank.GT.iProc) THEN
    CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
    CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
  END IF

  IF (TmpNeigh.NEQV.PartMPI%isMPINeighbor(iProc)) THEN
    IF(printMPINeighborWarnings)THEN
      WRITE(*,*) 'WARNING: MPINeighbor set to TRUE',PartMPI%MyRank,iProc
    END IF
    IF(.NOT.PartMPI%isMPINeighbor(iProc))THEN
      PartMPI%isMPINeighbor(iProc) = .TRUE.
      PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
    END IF
  END IF
END DO


! fill list with neighbor proc id and add local neighbor id to PartHaloElemToProc
ALLOCATE( PartMPI%MPINeighbor(PartMPI%nMPINeighbors) &
        , PartMPI%GlobalToLocal(0:PartMPI%nProcs-1)  )
iMPINeighbor=0
PartMPI%GlobalToLocal=-1
DO iProc=0,PartMPI%nProcs-1
  IF(PartMPI%isMPINeighbor(iProc))THEN
    iMPINeighbor=iMPINeighbor+1
    PartMPI%MPINeighbor(iMPINeighbor)=iProc
    PartMPI%GlobalToLocal(iProc)     =iMPINeighbor
    DO iElem=PP_nElems+1,nTotalElems
      IF(iProc.EQ.PartHaloElemToProc(NATIVE_PROC_ID,iElem)) PartHaloElemToProc(LOCAL_PROC_ID,iElem)=iMPINeighbor
    END DO ! iElem
  END IF
END DO

IF(iMPINeighbor.NE.PartMPI%nMPINeighbors) CALL abort(&
  __STAMP__&
  , ' Found number of mpi neighbors does not match! ', iMPINeighbor,REAL(PartMPI%nMPINeighbors))

IF(PartMPI%nMPINeighbors.GT.0)THEN
  IF(ANY(PartHaloElemToProc(LOCAL_PROC_ID,:).EQ.-1)) IPWRITE(UNIT_stdOut,*) ' Local proc id not found'
  IF(MAXVAL(PartHaloElemToProc(LOCAL_PROC_ID,:)).GT.PartMPI%nMPINeighbors) IPWRITE(UNIT_stdOut,*) ' Local proc id too high.'
  IF(MINVAL(PartHaloElemToProc(NATIVE_ELEM_ID,:)).LT.1) IPWRITE(UNIT_stdOut,*) ' native elem id too low'
  IF(MINVAL(PartHaloElemToProc(NATIVE_PROC_ID,:)).LT.0) IPWRITE(UNIT_stdOut,*) ' native proc id not found'
  IF(MAXVAL(PartHaloElemToProc(NATIVE_PROC_ID,:)).GT.PartMPI%nProcs-1) IPWRITE(UNIT_stdOut,*) ' native proc id too high.'
END IF

END SUBROUTINE InitHaloMesh


SUBROUTINE InitEmissionComm()
!===================================================================================================================================
! build emission communicators for particle emission regions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
USE MOD_Particle_Vars,          ONLY:Species,nSpecies
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
#ifndef PP_HDG
USE MOD_CalcTimeStep,           ONLY:CalcTimeStep
#endif /*PP_HDG*/
USE MOD_Particle_MPI_Vars,      ONLY:halo_eps
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec,iInit,iNode,iRank
INTEGER                         :: nInitRegions
LOGICAL                         :: RegionOnProc
REAL                            :: xCoords(3,8),lineVector(3),radius,height
REAL                            :: xlen,ylen,zlen
INTEGER                         :: color,iProc
INTEGER                         :: noInitRank,InitRank
!INTEGER,ALLOCATABLE             :: DummyRank(:)
LOGICAL                         :: hasRegion
!===================================================================================================================================

! get number of total init regions
nInitRegions=0
DO iSpec=1,nSpecies
  nInitRegions=nInitRegions+Species(iSpec)%NumberOfInits+(1-Species(iSpec)%StartnumberOfInits)
END DO ! iSpec
IF(nInitRegions.EQ.0) RETURN

! allocate communicators
ALLOCATE( PartMPI%InitGroup(1:nInitRegions))

nInitRegions=0
DO iSpec=1,nSpecies
  RegionOnProc=.FALSE.
  DO iInit=Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
    nInitRegions=nInitRegions+1
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE ('point')
       xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
       RegionOnProc=PointInProc(xCoords(1:3,1))
    CASE ('line_with_equidistant_distribution')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE ('line')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE ('Gaussian')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('disc')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('circle')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('circle_equidistant')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('cuboid')
      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL abort(&
         __STAMP__&
         ,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC&
                                                           +Species(iSpec)%Init(iInit)%BaseVector2IC

      IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
        height = halo_eps
      ELSE
        height= Species(iSpec)%Init(iInit)%CuboidHeightIC
      END IF
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cylinder')
      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL abort(&
         __STAMP__&
         ,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      radius = Species(iSpec)%Init(iInit)%RadiusIC
      ! here no radius, already inclueded
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC-Species(iSpec)%Init(iInit)%BaseVector1IC &
                                                           -Species(iSpec)%Init(iInit)%BaseVector2IC

      xCoords(1:3,2)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector1IC&
                                   +2.0*Species(iSpec)%Init(iInit)%BaseVector2IC

      IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
        height = halo_eps
      ELSE
        height= Species(iSpec)%Init(iInit)%CylinderHeightIC
      END IF
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cuboid_equal')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL abort(&
          __STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)

     !~j CALL abort(&
     !~j __STAMP__&
     !~j ,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
    CASE ('cuboid_with_equidistant_distribution')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL abort(&
          __STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE('sin_deviation')
       IF(Species(iSpec)%Init(iInit)%initialParticleNumber.NE. &
            (Species(iSpec)%Init(iInit)%maxParticleNumberX * Species(iSpec)%Init(iInit)%maxParticleNumberY &
            * Species(iSpec)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',iSpec,' does not match number of particles in each direction!'
         CALL abort(&
         __STAMP__&
         ,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = abs(GEO%xmaxglob  - GEO%xminglob)
       ylen = abs(GEO%ymaxglob  - GEO%yminglob)
       zlen = abs(GEO%zmaxglob  - GEO%zminglob)
       xCoords(1:3,1) = (/GEO%xminglob,GEO%yminglob,GEO%zminglob/)
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,'not implemented')
    END SELECT
    ! create new communicator
    color=MPI_UNDEFINED
    IF(RegionOnProc) color=nInitRegions!+1
    ! set communicator id
    Species(iSpec)%Init(iInit)%InitComm=nInitRegions

    ! create ranks for RP communicator
    IF(PartMPI%MPIRoot) THEN
      InitRank=-1
      noInitRank=-1
      iRank=0
      PartMPI%InitGroup(nInitRegions)%MyRank=0
      IF(RegionOnProc) THEN
        InitRank=0
      ELSE
        noInitRank=0
      END IF
      DO iProc=1,PartMPI%nProcs-1
        CALL MPI_RECV(hasRegion,1,MPI_LOGICAL,iProc,0,PartMPI%COMM,MPIstatus,iError)
        IF(hasRegion) THEN
          InitRank=InitRank+1
          CALL MPI_SEND(InitRank,1,MPI_INTEGER,iProc,0,PartMPI%COMM,iError)
        ELSE
          noInitRank=noInitRank+1
          CALL MPI_SEND(noInitRank,1,MPI_INTEGER,iProc,0,PartMPI%COMM,iError)
        END IF
      END DO
    ELSE
      CALL MPI_SEND(RegionOnProc,1,MPI_LOGICAL,0,0,PartMPI%COMM,iError)
      CALL MPI_RECV(PartMPI%InitGroup(nInitRegions)%MyRank,1,MPI_INTEGER,0,0,PartMPI%COMM,MPIstatus,iError)
    END IF

    ! create new emission communicator
    CALL MPI_COMM_SPLIT(PartMPI%COMM, color, PartMPI%InitGroup(nInitRegions)%MyRank, PartMPI%InitGroup(nInitRegions)%COMM,iError)
    IF(RegionOnProc) CALL MPI_COMM_SIZE(PartMPI%InitGroup(nInitRegions)%COMM,PartMPI%InitGroup(nInitRegions)%nProcs ,iError)
    IF(PartMPI%InitGroup(nInitRegions)%MyRank.EQ.0 .AND. RegionOnProc) &
    WRITE(UNIT_StdOut,*) 'Emission-Region,Emission-Communicator:',nInitRegions,PartMPI%InitGroup(nInitRegions)%nProcs,' procs'
    IF(PartMPI%InitGroup(nInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      IF(PartMPI%InitGroup(nInitRegions)%MyRank.EQ.0) THEN
        PartMPI%InitGroup(nInitRegions)%MPIRoot=.TRUE.
      ELSE
        PartMPI%InitGroup(nInitRegions)%MPIRoot=.FALSE.
      END IF
      ALLOCATE(PartMPI%InitGroup(nInitRegions)%GroupToComm(0:PartMPI%InitGroup(nInitRegions)%nProcs-1))
      PartMPI%InitGroup(nInitRegions)%GroupToComm(PartMPI%InitGroup(nInitRegions)%MyRank)= PartMPI%MyRank
      CALL MPI_ALLGATHER(PartMPI%MyRank,1,MPI_INTEGER&
                        ,PartMPI%InitGroup(nInitRegions)%GroupToComm(0:PartMPI%InitGroup(nInitRegions)%nProcs-1)&
                       ,1,MPI_INTEGER,PartMPI%InitGroup(nInitRegions)%COMM,iERROR)
      ALLOCATE(PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1))
      PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1)=-1
      DO iRank=0,PartMPI%InitGroup(nInitRegions)%nProcs-1
        PartMPI%InitGroup(nInitRegions)%CommToGroup(PartMPI%InitGroup(nInitRegions)%GroupToComm(iRank))=iRank
      END DO ! iRank
    END IF
  END DO ! iniT
END DO ! iSpec

END SUBROUTINE InitEmissionComm


FUNCTION BoxInProc(CartNodes,nNodes)
!===================================================================================================================================
! check if bounding box is on proc
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: CartNodes(1:3,1:nNodes)
INTEGER,INTENT(IN):: nNodes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: BoxInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
!===================================================================================================================================

BoxInProc=.FALSE.
! get background of nodes
xmin=HUGE(1)
xmax=-HUGE(1)
ymin=HUGE(1)
ymax=-HUGE(1)
zmin=HUGE(1)
zmax=-HUGE(1)
testval = CEILING((MINVAL(CartNodes(1,:))-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmin    = MIN(xmin,testval)
testval = CEILING((MAXVAL(CartNodes(1,:))-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmax    = MAX(xmax,testval)
testval = CEILING((MINVAL(CartNodes(2,:))-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymin    = MIN(ymin,testval)
testval = CEILING((MAXVAL(CartNodes(2,:))-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymax    = MAX(ymax,testval)
testval = CEILING((MINVAL(CartNodes(3,:))-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmin    = MIN(zmin,testval)
testval = CEILING((MAXVAL(CartNodes(3,:))-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmax    = MAX(zmax,testval)

IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) BoxInProc=.TRUE.

END FUNCTION BoxInProc


FUNCTION PointInProc(CartNode)
!===================================================================================================================================
! check if point is on proc
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: CartNode(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: PointInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
!===================================================================================================================================

PointInProc=.FALSE.
! get background of nodes
xmin=HUGE(1)
xmax=-HUGE(1)
ymin=HUGE(1)
ymax=-HUGE(1)
zmin=HUGE(1)
zmax=-HUGE(1)
testval = CEILING((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmin    = MIN(xmin,testval)
testval = CEILING((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmax    = MAX(xmax,testval)
testval = CEILING((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymin    = MIN(ymin,testval)
testval = CEILING((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymax    = MAX(ymax,testval)
testval = CEILING((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmin    = MIN(zmin,testval)
testval = CEILING((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmax    = MAX(zmax,testval)

IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) PointInProc=.TRUE.

END FUNCTION PointInProc


#if CODE_ANALYZE
SUBROUTINE CheckArrays(nTotalSides,nTotalElems,nTotalBCSides)
!===================================================================================================================================
! check if any entry of the checked arrays exists and if the entry is not NAN
! instead of using IEEE standard, the infamous nan-check a(i).NE.a(i) is used
! Sanity check for refmapping and mpi-communication.
! PO: not sure if it is required any more.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartHaloElemToProc
USE MOD_Mesh_Vars,              ONLY:BC,nGeo
USE MOD_Particle_Mesh_Vars,     ONLY:XCL_NGeo,DXCL_NGEO
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicType,PartBCSideList
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: nTotalSides,nTotalElems
INTEGER,INTENT(IN),OPTIONAL        :: nTotalBCSides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iElem,iVar,iVar2,i,j,k
INTEGER                            :: ilocSide,iSide
!===================================================================================================================================

DO iElem=1,nTotalElems
  ! PartElemToSide
  DO ilocSide=1,6
    DO iVar=1,2
      IF(PartElemToSide(iVar,ilocSide,iElem).NE.PartElemToSide(iVar,ilocSide,iElem)) CALL abort(&
__STAMP__&
, ' Error in PartElemToSide')
    END DO ! iVar=1,2
  END DO ! ilocSide=1,6
  IF(DoRefMapping)THEN
    ! XCL_NGeo & dXCL_NGeo
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          DO iVar=1,3
            IF(XCL_NGeo(iVar,i,j,k,iElem).NE.XCL_NGeo(iVar,i,j,k,iElem)) CALL abort(&
__STAMP__&
, ' Error in XCL_NGeo')
            DO iVar2=1,3
              IF(dXCL_NGeo(iVar2,iVar,i,j,k,iElem).NE.dXCL_NGeo(iVar2,iVar,i,j,k,iElem)) CALL abort(&
__STAMP__&
, ' Error in dXCL_NGeo')
            END DO ! iVar2=1,3
          END DO ! iVar=1,3
        END DO ! i=0,NGeo
      END DO ! j=0,NGeo
    END DO ! k=0,NGeo
  END IF ! DoRefMapping

END DO ! iElem=1,nTotalElems
IF(DoRefMapping)THEN
  ! PartBCSideList
  DO iSide = 1,nTotalSides
    IF(PartBCSideList(iSide).NE.PartBCSideList(iSide)) CALL abort(&
__STAMP__&
        , ' Error in dXCL_NGeo')
  END DO ! iSide=1,nTotalSides
  ! BezierControlPoints3D
  DO iSide=1,nTotalBCSides
    DO k=0,NGeo
      DO j=0,NGeo
        DO iVar=1,3
          IF(BezierControlPoints3D(iVar,j,k,iSide) &
         .NE.BezierControlPoints3D(iVar,j,k,iSide)) CALL abort(&
__STAMP__&
, ' Error in dXCL_NGeo')
        END DO ! iVar=1,3
      END DO ! j=0,nGeo
    END DO ! k=0,nGeo
    ! Slabnormals & SideSlabIntervals
    DO iVar=1,3
      DO iVar2=1,3
        IF(SideSlabNormals(iVar2,iVar,iSide).NE.SideSlabNormals(iVar2,iVar,iSide)) CALL abort(&
__STAMP__&
, ' Error in PartHaloElemToProc')
      END DO ! iVar2=1,PP_nVar
    END DO ! iVar=1,PP_nVar
    DO ilocSide=1,6
      IF(SideSlabIntervals(ilocSide,iSide).NE.SideSlabIntervals(ilocSide,iSide)) CALL abort(&
__STAMP__&
, ' Error in SlabInvervalls')
    END DO ! ilocSide=1,6
    IF(BoundingBoxIsEmpty(iSide).NEQV.BoundingBoxIsEmpty(iSide)) CALL abort(&
__STAMP__&
, ' Error in BoundingBoxIsEmpty')
  END DO ! iSide=1,nTotalBCSides
END IF ! DoRefMapping
! PartHaloElemToProc
DO iElem=PP_nElems+1,nTotalElems
  DO iVar=1,3
    IF(PartHaloElemToProc(iVar,iElem).NE.PartHaloElemToProc(iVar,iElem)) CALL abort(&
__STAMP__&
, ' Error in PartHaloElemToProc')
  END DO ! iVar=1,3
END DO ! iElem=PP_nElems+1,nTotalElems
DO iSide = 1,nTotalSides
  DO iVar=1,5
    IF(PartSideToElem(iVar,iSide).NE.PartSideToElem(iVar,iSide)) CALL abort(&
        __STAMP__&
        , ' Error in PartSideToElem')
  END DO ! iVar=1,5
  IF(SidePeriodicType(iSide).NE.SidePeriodicType(iSide)) CALL abort(&
      __STAMP__&
      , ' Error in BCSideType')
  IF(BC(iSide).NE.BC(iSide)) CALL abort(&
      __STAMP__&
      , ' Error in BC')
END DO ! iSide=1,nTotalSides


END SUBROUTINE CheckArrays
#endif /*CODE_ANALYZE*/
#endif /*MPI*/

END MODULE MOD_Particle_MPI
