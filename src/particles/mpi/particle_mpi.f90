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

INTERFACE InitParticleCommSize
  MODULE PROCEDURE InitParticleCommSize
END INTERFACE

INTERFACE InitEmissionComm
  MODULE PROCEDURE InitEmissionComm
END INTERFACE

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

!INTERFACE InitHaloMesh
!  MODULE PROCEDURE InitHaloMesh
!END INTERFACE

!INTERFACE ExchangeBezierControlPoints3D
!  MODULE PROCEDURE ExchangeBezierControlPoints3D
!END INTERFACE

PUBLIC :: InitParticleMPI
PUBLIC :: InitEmissionComm
PUBLIC :: InitParticleCommSize
PUBLIC :: SendNbOfParticles
PUBLIC :: IRecvNbOfParticles
PUBLIC :: MPIParticleSend
PUBLIC :: MPIParticleRecv
PUBLIC :: FinalizeParticleMPI
!PUBLIC :: InitHaloMesh
!PUBLIC :: ExchangeBezierControlPoints3D
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
USE MOD_Particle_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                         :: color
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI ... '
IF(ParticleMPIInitIsDone) &
  CALL ABORT(__STAMP__,' Particle MPI already initialized!')

#if USE_MPI
PartMPI%myRank = myRank
color = 999
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,PartMPI%MyRank,PartMPI%COMM,iERROR)
CALL MPI_COMM_SIZE (PartMPI%COMM,PartMPI%nProcs ,iError)
!PartMPI%COMM   = MPI_COMM_WORLD
IF(PartMPI%nProcs.NE.nProcessors) CALL ABORT(__STAMP__,' MPI Communicater-size does not match!', IERROR)
PartCommSize   = 0
IF(PartMPI%MyRank.EQ.0) THEN
  PartMPI%MPIRoot=.TRUE.
ELSE
  PartMPI%MPIRoot=.FALSE.
END IF
iMessage=0
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
USE MOD_Particle_Analyze_Vars,    ONLY:doParticleDispersionTrack
USE MOD_Particle_Boundary_Vars,    ONLY:doParticleReflectionTrack
USE MOD_Particle_MPI_Vars
USE MOD_Particle_SGS_Vars,        ONLY:nSGSVars!,SGSinUse
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
! IF (SGSinUse)    PartCommSize = PartCommSize + 3
! IsNewPart for RK-Reconstruction
PartCommSize   = PartCommSize + 1

! Include reflection counter
IF(doParticleReflectionTrack) PartCommSize = PartCommSize + 1

! Include dispersion path
IF(doParticleDispersionTrack) PartCommSize = PartCommSize + 3

! additional stuff for full RK schemes, e.g. implicit and imex RK
! if iStage=0, then the PartStateN is not communicated
PartCommSize0  = PartCommSize

ALLOCATE( PartMPIExchange%nPartsSend(2,0:nExchangeProcessors-1)  &
        , PartMPIExchange%nPartsRecv(2,0:nExchangeProcessors-1)  &
        , PartRecvBuf(0:nExchangeProcessors-1)                   &
        , PartSendBuf(0:nExchangeProcessors-1)                   &
        , PartMPIExchange%SendRequest(2,0:nExchangeProcessors-1) &
        , PartMPIExchange%RecvRequest(2,0:nExchangeProcessors-1) &
        , PartTargetProc(1:PDM%MaxParticleNumber)              &
        , STAT=ALLOCSTAT                                       )

IF (ALLOCSTAT.NE.0) &
  CALL ABORT(__STAMP__,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)

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
USE MOD_Particle_MPI_Vars,      ONLY:nExchangeProcessors,ExchangeProcToGlobalProc
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
DO iProc=0,nExchangeProcessors-1
  CALL MPI_IRECV( PartMPIExchange%nPartsRecv(:,iProc)                        &
                , 2                                                          &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(1,iProc)                       &
                , IERROR )
! IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__&
!         ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


SUBROUTINE SendNbOfParticles(doParticle_In)
!===================================================================================================================================
! this routine sends the number of send particles. The following steps are performed:
! 1) Compute number of Send Particles
! 2) Perform MPI_ISEND with number of particles
! The remaining steps are performed in SendParticles
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI,PartMPIExchange,PartTargetProc
USE MOD_Particle_MPI_Vars      ,ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc
USE MOD_Particle_Vars          ,ONLY: PEM,PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: doPartInExists
INTEGER                       :: iPart,ElemID
INTEGER                       :: iProc,ProcID
!===================================================================================================================================
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.

! 1) get number of send particles
!--- Count number of particles in cells in the halo region and add them to the message
!--- CAUTION: using local indices for halo elem -> proc association with PartHaloElemToProc. PartMPI%Neighbor contains the inverse
!--- mapping
PartMPIExchange%nPartsSend=0
PartTargetProc=-1
DO iPart=1,PDM%ParticleVecLength
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  
  ElemID=PEM%Element(iPart)
  ProcID = ElemInfo_Shared(ELEM_RANK,ElemID)

  ! Particle on local proc, do nothing
  IF (ProcID.EQ.myRank) CYCLE

  ! Add particle to target proc count
    PartMPIExchange%nPartsSend(1,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) =  &
      PartMPIExchange%nPartsSend(1,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)) + 1
    PartTargetProc(iPart) = GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,ProcID)
END DO ! iPart


! 2) send number of send particles
!--- Loop over all neighboring procs. Map local proc ID to global through ExchangeProcToGlobalProc.
!--- Asynchronous communication, just send here and check for success later.
DO iProc=0,nExchangeProcessors-1
  CALL MPI_ISEND( PartMPIExchange%nPartsSend(:,iProc)                        &
                , 2                                                          &
                , MPI_INTEGER                                                &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(1,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE SendNbOfParticles


SUBROUTINE MPIParticleSend()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! first steps are performed in SendNbOfParticles
! 1) Compute number of Send Particles
! 2) Perform MPI_ISEND with number of particles
! Starting Here:
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
! CAUTION: If particles are sent for deposition, PartTargetProc has the information, if a particle is sent
!          and after the build and wait for number of particles reused to build array with external parts
!          information in PartState,.. can be reused, because they are not overwritten
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars,    ONLY:PartPath,doParticleDispersionTrack
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartCommSize,PartSendBuf,PartRecvBuf,PartTargetProc
USE MOD_Particle_MPI_Vars,        ONLY:nExchangeProcessors,ExchangeProcToGlobalProc
USE MOD_Particle_Vars,            ONLY:PartSpecies,PEM,PDM,PartPosRef
USE MOD_Particle_Vars,            ONLY:PartState,Pt_temp
USE MOD_Particle_Vars,            ONLY:TurbPartState!,TurbPt_temp
USE MOD_Particle_Tracking_Vars,   ONLY:DoRefMapping
! Variables for erosion tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Boundary_Vars
! Variables for SGS model
USE MOD_Particle_SGS_Vars,        ONLY:nSGSVars
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
INTEGER                       :: iPart,iPos,iProc,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,0:nExchangeProcessors-1)
INTEGER                       :: MessageSize, nRecvParticles, nSendParticles
INTEGER                       :: ALLOCSTAT
!===================================================================================================================================

! 3) Build Message
DO iProc=0,nExchangeProcessors-1
  ! allocate SendBuf and prepare to build message
  nSendParticles = PartMPIExchange%nPartsSend(1,iProc)
  iPos=0

  ! only build message if we have particles to send
  IF(nSendParticles.EQ.0) CYCLE

  ! allocate SendBuff of required size
  MessageSize=nSendParticles*PartCommSize
  ALLOCATE(PartSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) &
    CALL ABORT(__STAMP__,'  Cannot allocate PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))

  ! fill message
  !--- Loop over all particles
  DO iPart=1,PDM%ParticleVecLength
    IF(nSendParticles.EQ.0) EXIT

    ! particle belongs on target proc
    IF(PartTargetProc(iPart).EQ.iProc) THEN
      !>> particle position in physical space
      PartSendBuf(iProc)%content(1+iPos:6+iPos) = PartState(1:6,iPart)
      jpos=iPos+6
      IF (ALLOCATED(TurbPartState)) THEN
        !>> SGS turbulent velocity and random draw
        PartSendBuf(iProc)%content(1+jPos:nSGSVars+jPos) = TurbPartState(1:nSGSVars,iPart)
        jpos=jpos+nSGSVars
#if USE_RW
        !>> RW turbulent velocity, interaction time and random draw
        PartSendBuf(iProc)%content(1+jPos:nRWVars+jPos) = TurbPartState(1:nRWVars,iPart)
        jpos=jpos+nRWVars
#endif
      END IF
      !>> particle position in reference space
      IF(DoRefMapping) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = PartPosRef(1:3,iPart)
        jPos=jPos+3
      END IF
      !>> reflection counter
      IF (doParticleReflectionTrack) THEN
        PartSendBuf(iProc)%content(1+jPos) = REAL(PartReflCount(iPart),KIND=8)
        jPos=jPos+1
      END IF
      !>> absolute particle path
      IF (doParticleDispersionTrack) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = PartPath(1:3,iPart)
        jPos=jPos+3
      END IF
      !>> particles species
      PartSendBuf(iProc)%content(1+jPos) = REAL(PartSpecies(iPart),KIND=8)
      jPos=jPos+1
      !>> Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
      PartSendBuf(iProc)%content(1+jPos:6+jPos) = Pt_temp(1:6,iPart)
      jPos=jPos+6
      !>> TurbPt_tmp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
!      IF (SGSinUse) THEN
!        PartSendBuf(iProc)%content(1+jPos:3+jPos) = TurbPt_temp(1:3,iPart)
!        jpos=jpos+3
!      END IF
      !>> IsNewPart for RK-Reconstruction
      IF (PDM%IsNewPart(iPart)) THEN
        PartSendBuf(iProc)%content(1+jPos) = 1.
      ELSE
        PartSendBuf(iProc)%content(1+jPos) = 0.
      END IF
      jPos=jPos+1
      !>> particle element
      PartSendBuf(iProc)%content(    1+jPos) = REAL(PEM%Element(iPart),KIND=8)
      jPos=jPos+1

      ! sanity check the message length. PartCommSize must be a multiple of particles to send
      IF(MOD(jPos,PartCommSize).NE.0) THEN
        IPWRITE(UNIT_stdOut,*)  'PartCommSize',PartCommSize
        IPWRITE(UNIT_stdOut,*)  'jPos',jPos
        CALL ABORT( __STAMP__,' Particle-wrong sending message size!')
      END IF

      ! increment message position to next element, PartCommSize.EQ.jPos
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
DO iProc=0,nExchangeProcessors-1
  CALL MPI_WAIT(PartMPIExchange%SendRequest(1,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

! total number of received particles
PartMPIExchange%nMPIParticles=SUM(PartMPIExchange%nPartsRecv(1,:))

! nullify data on old particle position for safety
DO iPart=1,PDM%ParticleVecLength
  IF(PartTargetProc(iPart).EQ.-1) CYCLE
  PartState(1:6,iPart) = 0.
  PartSpecies(iPart)   = 0
  Pt_temp(1:6,iPart)   = 0.
  IF (doParticleReflectionTrack) PartReflCount(  iPart) = 0
  IF (doParticleDispersionTrack) PartPath(    :,iPart) = 0.
  IF (ALLOCATED(TurbPartState))  TurbPartState(:,iPart) = 0.
END DO ! iPart=1,PDM%ParticleVecLength

! 5) Allocate received buffer and open MPI_IRECV
DO iProc=0,nExchangeProcessors-1

  ! skip proc if no particles are to be received
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE

  ! count number of particles from each proc and determine size of message
  nRecvParticles = PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize    = nRecvParticles * PartCommSize

  ! allocate recv buffer with the correct size
  ALLOCATE(PartRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    IPWRITE(*,*) 'sum of total received particles            ', SUM(PartMPIExchange%nPartsRecv(1,:))
    IPWRITE(*,*) 'sum of total received deposition particles ', SUM(PartMPIExchange%nPartsRecv(2,:))
    CALL ABORT(__STAMP__,'  Cannot allocate PartRecvBuf, local source ProcId, Allocstat',iProc,REAL(ALLOCSTAT))
  END IF

  ! start asynchronous MPI recv for each proc
  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

! 6) Send Particles
DO iProc=0,nExchangeProcessors-1

  ! skip proc if no particles are to be sent
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE

  ! count number of particles for each proc and determine size of message
  nSendParticles = PartMPIExchange%nPartsSend(1,iProc)
  MessageSize    = nSendParticles*PartCommSize

  ! start asychronous MPI send for each proc
  CALL MPI_ISEND( PartSendBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc)         &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)

  ! Deallocate sendBuffer after send was successful, see MPIParticleRecv
END DO ! iProc

END SUBROUTINE MPIParticleSend


SUBROUTINE MPIParticleRecv()
!===================================================================================================================================
! this routine finishes the communication and places the particle information in the correct arrays. Following steps are performed
! 1) Finish all send requests -> MPI_WAIT
! 2) Finish all recv requests -> MPI_WAIT
! 3) Place particle information in correct arrays
! 4) Deallocate send and recv buffers
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars,    ONLY:PartPath,doParticleDispersionTrack
USE MOD_Particle_MPI_Vars,        ONLY:PartMPIExchange,PartCommSize,PartRecvBuf,PartSendBuf!,PartMPI
USE MOD_Particle_MPI_Vars,        ONLY:nExchangeProcessors
USE MOD_Particle_Vars,            ONLY:PartSpecies,PEM,PDM,PartPosRef
USE MOD_Particle_Vars,            ONLY:PartState,Pt_temp
USE MOD_Particle_Vars,            ONLY:TurbPartState!,TurbPt_temp
USE MOD_Particle_Tracking_Vars,   ONLY:DoRefMapping
! variables for erosion tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Boundary_Vars
! Variables for SGS model
USE MOD_Particle_SGS_Vars,        ONLY:nSGSVars
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
INTEGER                       :: iProc,iPos,nRecv,PartID,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,0:nExchangeProcessors-1)
INTEGER                       :: MessageSize,nRecvParticles
!===================================================================================================================================

! wait for all send requests to be successful
DO iProc=0,nExchangeProcessors-1
  ! skip proc if no particles are to be sent
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE
  
  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc


nRecv=0
! decompose the recv message from each proc and set local variables accordingly
DO iProc=0,nExchangeProcessors-1
  ! skip proc if no particles are to be received
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE

  ! count number of recv particles and determine message sage
  nRecvParticles = PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize    = nRecvParticles*PartCommSize

  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
  
  ! place particle information in correct arrays
  DO iPos=0,MessageSize-1,PartCommSize
    ! find free position in particle array
    nRecv  = nRecv+1

    ! particles get a local ID on each proc, therefore put it at the next free position
    PartID = PDM%nextFreePosition(nRecv+PDM%CurrentNextFreePosition)
    IF(PartID.EQ.0) &
      CALL ABORT(__STAMP__,' Error in ParticleExchange_parallel. Corrupted list: PIC%nextFreePosition', nRecv)

    !>> position and velocity in physical space
    PartState(1:6,PartID)   = PartRecvBuf(iProc)%content(1+iPos:6+iPos)
    jpos=iPos+6
    IF (ALLOCATED(TurbPartState)) THEN
      !>> SGS turbulent velocity and random draw
      TurbPartState(1:nSGSVars,PartID) = PartRecvBuf(iProc)%content(1+jpos:nSGSVars+jpos)
      jpos=jpos+nSGSVars
#if USE_RW
      !>>! RW turbulent velocity, interaction time and random draw
      TurbPartState(1:nRWVars,PartID) = PartRecvBuf(iProc)%content(1+jPos:nRWVars+jPos)
      jpos=jpos+nRWVars
#endif
    END IF
    !>> particle position in reference space
    IF(DoRefMapping) THEN
      PartPosRef(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos:3+jPos)
      jPos=jPos+3
    END IF
    !>> reflection counter
    IF (doParticleReflectionTrack) THEN
      PartReflCount(PartID)  = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
      jpos=jpos+1
    END IF
    !>> absolute particle path
    IF (doParticleDispersionTrack) THEN
      PartPath(1:3,PartID)   = PartRecvBuf(iProc)%content(1+jPos:3+jPos)
      jPos=jPos+3
    END IF
    !>> particles species
    PartSpecies(PartID)      = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1
    !>> Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
    Pt_temp(1:6,PartID)      = PartRecvBuf(iProc)%content(1+jPos:6+jPos)
    jpos=jpos+6
    !>> TurbPt_tmp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
!    IF (SGSinUse) THEN
!      TurbPt_temp(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos:3+jPos)
!      jpos=jpos+3
!    END IF
    !>> IsNewPart for RK-Reconstruction
    IF      ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 1) THEN
      PDM%IsNewPart(PartID)=.TRUE.
    ELSE IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 0) THEN
      PDM%IsNewPart(PartID)=.FALSE.
    ELSE
      CALL ABORT(__STAMP__,'Error with IsNewPart in MPIParticleRecv!',1,PartRecvBuf(iProc)%content(1+jPos))
    END IF
    jPos=jPos+1
    !>> particle elment
    PEM%Element(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1

    ! PartCommSize must be a multiple of particles to receive
    IF(MOD(jPos,PartCommSize).NE.0)THEN
      IPWRITE(UNIT_stdOut,*)  'jPos',jPos
      CALL ABORT(__STAMP__,' Particle-wrong receiving message size!')
    END IF

    ! Set Flag for received parts in order to localize them later
    PDM%ParticleInside(PartID) = .TRUE.
    PEM%lastElement(PartID)    = -888
  END DO
END DO ! iProc


PDM%ParticleVecLength       = PDM%ParticleVecLength + PartMPIExchange%nMPIParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles
IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) &
  CALL ABORT(__STAMP__,' ParticleVecLegnth>MaxParticleNumber due to MPI-communication!')

! deallocate send,receive buffer
DO iProc=0,nExchangeProcessors-1
  SDEALLOCATE(PartRecvBuf(iProc)%content)
  SDEALLOCATE(PartSendBuf(iProc)%content)
END DO ! iProc

! last step, nullify number of sent and received particles
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
IF(PartMPI%COMM.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(PartMPI%COMM,iERROR)

!SDEALLOCATE( PartHaloElemToProc)
!SDEALLOCATE( PartHaloNodeToProc)
!SDEALLOCATE( PartMPI%isMPINeighbor)
!SDEALLOCATE( PartMPI%MPINeighbor )
!SDEALLOCATE( PartMPI%GlobalToLocal )
SDEALLOCATE( PartMPIExchange%nPartsSend)
SDEALLOCATE( PartMPIExchange%nPartsRecv)
SDEALLOCATE( PartMPIExchange%RecvRequest)
SDEALLOCATE( PartMPIExchange%SendRequest)
SDEALLOCATE( PartMPIExchange%Send_message)
!SDEALLOCATE( PartMPI%isMPINodeNeighbor)
!SDEALLOCATE( PartMPI%MPINodeNeighbor)
SDEALLOCATE( PartMPI%InitGroup)
SDEALLOCATE( PartSendBuf)
SDEALLOCATE( PartRecvBuf)

SDEALLOCATE( ExchangeProcToGlobalProc)
SDEALLOCATE( GlobalProcToExchangeProc)

! and for communication
SDEALLOCATE( PartTargetProc )

ParticleMPIInitIsDone=.FALSE.
END SUBROUTINE FinalizeParticleMPI


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
         CALL ABORT(&
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
         CALL ABORT(&
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
          CALL ABORT(&
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

     !~j CALL ABORT(&
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
          CALL ABORT(&
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
         CALL ABORT(&
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
      CALL ABORT(&
      __STAMP__&
      ,'ERROR: Given SpaceIC is not implemented!')
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
#endif /*MPI*/

END MODULE MOD_Particle_MPI
