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

!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_MPI
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE DefineParticleMPI
  MODULE PROCEDURE DefineParticleMPI
END INTERFACE

INTERFACE InitParticleMPI
  MODULE PROCEDURE InitParticleMPI
END INTERFACE

INTERFACE InitParticleCommSize
  MODULE PROCEDURE InitParticleCommSize
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

PUBLIC :: DefineParticleMPI
PUBLIC :: InitParticleMPI
PUBLIC :: InitParticleCommSize
PUBLIC :: SendNbOfParticles
PUBLIC :: IRecvNbOfParticles
PUBLIC :: MPIParticleSend
PUBLIC :: MPIParticleRecv
PUBLIC :: FinalizeParticleMPI
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

#if USE_MPI
SUBROUTINE DefineParticleMPI
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools              ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL prms%CreateLogicalOption('CheckExchangeProcs' , 'Check if proc communication of particle info is non-symmetric', '.TRUE.')

END SUBROUTINE DefineParticleMPI


SUBROUTINE InitParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars
#if USE_MPI
USE MOD_ReadInTools              ,ONLY: GETLOGICAL
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!#if USE_MPI
!INTEGER                         :: color
!#endif /*USE_MPI*/
!===================================================================================================================================

!SWRITE(UNIT_stdOut,'(132("-"))')
!SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI... '
IF(ParticleMPIInitIsDone) CALL Abort(__STAMP__,' Particle MPI already initialized!')

#if USE_MPI
! Get flag for ignoring the abort if the number of global exchange procs is non-symmetric
CheckExchangeProcs = GETLOGICAL('CheckExchangeProcs')

! Duplicate particle communicator from MPI_COMM_FLEXI
CALL MPI_COMM_DUP (MPI_COMM_FLEXI,PartMPI%COMM,iError)
CALL MPI_COMM_RANK(PartMPI%COMM,PartMPI%myRank,iError)
CALL MPI_COMM_SIZE(PartMPI%COMM,PartMPI%nProcs,iError)

PartCommSize    = 0
PartMPI%MPIRoot = MERGE(.TRUE.,.FALSE.,PartMPI%MyRank.EQ.0)
#else
PartMPI%myRank  = 0
PartMPI%nProcs  = 1
PartMPI%MPIRoot=.TRUE.
#endif  /*MPI*/

ParticleMPIInitIsDone=.TRUE.
!SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI DONE!'
!SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitParticleMPI


SUBROUTINE InitParticleCommSize()
!===================================================================================================================================
! get size of Particle-MPI-Message. Unfortunately, this subroutine have to be called after particle_init because
! all required features have to be read from the ini-File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars,    ONLY:doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars,   ONLY:doParticleReflectionTrack
USE MOD_Particle_MPI_Vars
USE MOD_Particle_SGS_Vars,        ONLY:nSGSVars!,SGSinUse
USE MOD_Particle_Tracking_Vars,   ONLY:TrackingMethod
USE MOD_Particle_Vars,            ONLY:PDM,doPartIndex
#if USE_RW
USE MOD_Particle_RandomWalk_Vars, ONLY:nRWVars
#endif
#if USE_BASSETFORCE
USE MOD_Particle_Vars,            ONLY:nBassetVars,N_Basset!,FbCoeffm
#endif /* USE_BASSETFORCE */
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
PartCommSize   = PartCommSize + PP_nVarPart
! TurbPartState: SGS turbulent velocity and random draw
PartCommSize   = PartCommSize + nSGSVars
#if USE_RW
! TurbPartState: RW turbulent velocity, interaction time and random draw
PartCommSize   = PartCommSize + nRWVars
#endif
! Tracking: Include Reference coordinates
IF(TrackingMethod.EQ.REFMAPPING) PartCommSize=PartCommSize+3
! Species-ID
PartCommSize   = PartCommSize + 1
! Index
IF(doPartIndex)  PartCommSize = PartCommSize + 1
#if USE_BASSETFORCE
! Extended RHS: Basset force
PartCommSize   = PartCommSize + nBassetVars
PartCommSize   = PartCommSize + 1
PartCommSize   = PartCommSize + N_Basset+2
!PartCommSize   = PartCommSize + 3*FbCoeffm
#endif /* USE_BASSETFORCE */
! id of element
PartCommSize   = PartCommSize + 1

! time integration
! communication after each Runge-Kutta stage, so send time derivative must be communicated to the new proc
! Pt_temp for pushing: Runge-Kutta derivative of position and velocity
PartCommSize   = PartCommSize + (PP_nVarPart-1)
! TurbPt_tmp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
! IF (SGSinUse)    PartCommSize = PartCommSize + 3
! IsNewPart for RK-Reconstruction
PartCommSize   = PartCommSize + 1

! Include reflection counter
IF(doParticleReflectionTrack)                        PartCommSize = PartCommSize + 1

! Include dispersion path
IF(doParticleDispersionTrack.OR.doParticlePathTrack) PartCommSize = PartCommSize + 3

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
  CALL Abort(__STAMP__,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)

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
USE MOD_Particle_MPI_Vars,      ONLY: PartMPI,PartMPIExchange
USE MOD_Particle_MPI_Vars,      ONLY: nExchangeProcessors,ExchangeProcToGlobalProc
USE MOD_Particle_Vars,          ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iProc
!===================================================================================================================================

! Return if running particle code without any species
IF (nSpecies.LE.0) RETURN

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
! IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__&
!         ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


SUBROUTINE SendNbOfParticles()
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
USE MOD_Particle_Vars          ,ONLY: PEM,PDM,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID
INTEGER                       :: iProc,ProcID
!===================================================================================================================================

! Return if running particle code without any species
IF (nSpecies.LE.0) RETURN

! 1) get number of send particles
!--- Count number of particles in cells in the halo region and add them to the message
PartMPIExchange%nPartsSend = 0
PartTargetProc = -1

DO iPart=1,PDM%ParticleVecLength
  ! Cycle if particle does not exist on proc
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ElemID = PEM%Element(iPart)
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
  IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)
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
USE MOD_Particle_Analyze_Vars,    ONLY:PartPath,doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartCommSize,PartSendBuf,PartRecvBuf,PartTargetProc
USE MOD_Particle_MPI_Vars,        ONLY:nExchangeProcessors,ExchangeProcToGlobalProc
USE MOD_Particle_Vars,            ONLY:PartSpecies,PEM,PDM,PartPosRef,PartIndex,doPartIndex
USE MOD_Particle_Vars,            ONLY:PartState,Pt_temp,nSpecies
USE MOD_Particle_Vars,            ONLY:TurbPartState!,TurbPt_temp
USE MOD_Particle_Tracking_Vars,   ONLY:TrackingMethod
! Variables for impact tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Boundary_Vars
! Variables for SGS model
USE MOD_Particle_SGS_Vars,        ONLY:nSGSVars
#if USE_RW
! Variables for RW model
USE MOD_Particle_RandomWalk_Vars, ONLY:nRWVars
#endif
#if USE_BASSETFORCE
USE MOD_Particle_Vars,            ONLY:durdt,nBassetVars,bIter,N_Basset,Fbdt!,FbCoeffm,Fbi
#endif /* USE_BASSETFORCE */
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

! Return if running particle code without any species
IF (nSpecies.LE.0) RETURN

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
    CALL Abort(__STAMP__,'  Cannot allocate PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))

  ! fill message
  !--- Loop over all particles
  DO iPart=1,PDM%ParticleVecLength
    IF(nSendParticles.EQ.0) EXIT

    ! particle belongs on target proc
    IF(PartTargetProc(iPart).EQ.iProc) THEN
      !>> particle position in physical space
      PartSendBuf(iProc)%content(1+iPos:PP_nVarPart+iPos) = PartState(1:PP_nVarPart,iPart)
      jpos=iPos+PP_nVarPart
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
      IF(TrackingMethod.EQ.REFMAPPING) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = PartPosRef(1:3,iPart)
        jPos=jPos+3
      END IF
      !>> reflection counter
      IF (doParticleReflectionTrack) THEN
        PartSendBuf(iProc)%content(1+jPos) = REAL(PartReflCount(iPart),KIND=8)
        jPos=jPos+1
      END IF
      !>> absolute particle path
      IF (doParticleDispersionTrack.OR.doParticlePathTrack) THEN
        PartSendBuf(iProc)%content(1+jPos:3+jPos) = PartPath(1:3,iPart)
        jPos=jPos+3
      END IF
      !>> particles species
      PartSendBuf(iProc)%content(1+jPos) = REAL(PartSpecies(iPart),KIND=8)
      jPos=jPos+1
      !>> particles index
      IF (doPartIndex) THEN
        PartSendBuf(iProc)%content(1+jPos) = REAL(PartIndex(iPart),KIND=8)
        jPos=jPos+1
      END IF
#if USE_BASSETFORCE
      !>> Basset force
      PartSendBuf(iProc)%content(1+jPos:nBassetVars+jPos) = durdt(1:nBassetVars,iPart)
      jPos=jPos+nBassetVars
      PartSendBuf(iProc)%content(1+jPos) = REAL(bIter(iPart),KIND=8)
      jPos=jPos+1
      PartSendBuf(iProc)%content(1+jPos:N_Basset+2+jPos) = Fbdt(1:N_Basset+2,iPart)
      jPos=jPos+N_Basset+2
      !PartSendBuf(iProc)%content(1+jPos:3*FbCoeffm+jPos) = RESHAPE(Fbi(1:3,1:FbCoeffm,iPart),(/3*FbCoeffm/))
      !jPos=jPos+3*FbCoeffm
#endif /* USE_BASSETFORCE */
      !>> Pt_temp for pushing: Runge-Kutta derivative of position and velocity
      PartSendBuf(iProc)%content(1+jPos:PP_nVarPart-1+jPos) = Pt_temp(1:PP_nVarPart-1,iPart)
      jPos=jPos+PP_nVarPart-1
      !>> TurbPt_temp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
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
        IPWRITE(UNIT_stdOut,'(I0,A,I0)')  'PartCommSize',PartCommSize
        IPWRITE(UNIT_stdOut,'(I0,A,I0)')  'jPos',jPos
        CALL Abort( __STAMP__,' Particle-wrong sending message size!')
      END IF

      ! increment message position to next element, PartCommSize.EQ.jPos
      iPos=iPos+PartCommSize
      ! particle is ready for send, now it can deleted
      PDM%ParticleInside(iPart) = .FALSE.
    END IF ! Particle is particle with target proc-id equals local proc id
  END DO  ! iPart

  IF(iPos.NE.(MessageSize)) IPWRITE(UNIT_stdOut,'(I0,A,I0,I0)') ' error message size', iPos,MessageSize
END DO ! iProc

! 4) Finish Received number of particles
!--- Wait for all neighboring procs to acknowlege both our send and our recv request. Then every neighbor proc knows the number of
!--- particles to communicate
DO iProc=0,nExchangeProcessors-1
  CALL MPI_WAIT(PartMPIExchange%SendRequest(1,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)
END DO ! iProc

! total number of received particles
PartMPIExchange%nMPIParticles=SUM(PartMPIExchange%nPartsRecv(1,:))

! nullify data on old particle position for safety
DO iPart=1,PDM%ParticleVecLength
  IF(PartTargetProc(iPart).EQ.-1) CYCLE
  PartState(:,iPart) = 0.
  PartSpecies(iPart) = 0
  IF (doPartIndex) PartIndex(iPart) = 0
#if USE_BASSETFORCE
  durdt(:,iPart) = 0.
  bIter(iPart) = 0.
  Fbdt(:,iPart) = 0.
  !Fbi(:,:,iPart) = 0.
#endif
  Pt_temp(:,iPart) = 0.
  IF (doParticleReflectionTrack)                        PartReflCount(  iPart) = 0
  IF (doParticleDispersionTrack.OR.doParticlePathTrack) PartPath(     :,iPart) = 0.
  IF (ALLOCATED(TurbPartState))                         TurbPartState(:,iPart) = 0.
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
    IPWRITE(UNIT_stdOut,'(I0,A,I0)') 'sum of total received particles            ', SUM(PartMPIExchange%nPartsRecv(1,:))
    IPWRITE(UNIT_stdOut,'(I0,A,I0)') 'sum of total received deposition particles ', SUM(PartMPIExchange%nPartsRecv(2,:))
    CALL Abort(__STAMP__,'  Cannot allocate PartRecvBuf, local source ProcId, Allocstat',iProc,REAL(ALLOCSTAT))
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
  IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)
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
  IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)

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
USE MOD_Particle_Analyze_Vars,    ONLY:PartPath,doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_MPI_Vars,        ONLY:PartMPIExchange,PartCommSize,PartRecvBuf,PartSendBuf!,PartMPI
USE MOD_Particle_MPI_Vars,        ONLY:nExchangeProcessors
USE MOD_Particle_Vars,            ONLY:PartSpecies,PEM,PDM,PartPosRef,PartIndex,doPartIndex
USE MOD_Particle_Vars,            ONLY:PartState,Pt_temp,nSpecies
USE MOD_Particle_Vars,            ONLY:TurbPartState!,TurbPt_temp
USE MOD_Particle_Tracking_Vars,   ONLY:TrackingMethod
! variables for impact tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Boundary_Vars
! Variables for SGS model
USE MOD_Particle_SGS_Vars,        ONLY:nSGSVars
#if USE_RW
! Variables for RW model
USE MOD_Particle_RandomWalk_Vars, ONLY:nRWVars
#endif
#if USE_BASSETFORCE
USE MOD_Particle_Vars,            ONLY:durdt,nBassetVars,N_Basset,bIter,Fbdt!,FbCoeffm,Fbi
#endif /* USE_BASSETFORCE */
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

! Return if running particle code without any species
IF (nSpecies.LE.0) RETURN

! wait for all send requests to be successful
DO iProc=0,nExchangeProcessors-1
  ! skip proc if no particles are to be sent
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE

  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)
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
      CALL Abort(__STAMP__,' Error in ParticleExchange_parallel. Corrupted list: PIC%nextFreePosition', nRecv)

    !>> position and velocity in physical space
    PartState(1:PP_nVarPart,PartID)   = PartRecvBuf(iProc)%content(1+iPos:PP_nVarPart+iPos)
    jpos=iPos+PP_nVarPart
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
    IF(TrackingMethod.EQ.REFMAPPING) THEN
      PartPosRef(1:3,PartID) = PartRecvBuf(iProc)%content(1+jPos:3+jPos)
      jPos=jPos+3
    END IF
    !>> reflection counter
    IF (doParticleReflectionTrack) THEN
      PartReflCount(PartID)  = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
      jpos=jpos+1
    END IF
    !>> absolute particle path
    IF (doParticleDispersionTrack.OR.doParticlePathTrack) THEN
      PartPath(1:3,PartID)   = PartRecvBuf(iProc)%content(1+jPos:3+jPos)
      jPos=jPos+3
    END IF
    !>> particles species
    PartSpecies(PartID)      = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1
    !>> particles index
    IF (doPartIndex) THEN
      PartIndex(PartID)        = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
      jPos=jPos+1
    END IF
#if USE_BASSETFORCE
    !>> Basset force
    durdt(1:nBassetVars,PartID) = PartRecvBuf(iProc)%content(1+jPos:nBassetVars+jPos)
    jPos=jPos+nBassetVars
    bIter(PartID)               = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1
    Fbdt(1:N_Basset+2,PartID)   = PartRecvBuf(iProc)%content(1+jPos:N_Basset+2+jPos)
    jPos=jPos+N_Basset+2
    !Fbi(1:3,1:FbCoeffm,PartID)  = RESHAPE(PartRecvBuf(iProc)%content(1+jPos:3*FbCoeffm+jPos),(/3,FbCoeffm/))
    !jPos=jPos+3*FbCoeffm
#endif /* USE_BASSETFORCE */
    !>> Pt_temp for pushing: Runge-Kutta derivative of position and velocity
    Pt_temp(1:PP_nVarPart-1,PartID) = PartRecvBuf(iProc)%content(1+jPos:PP_nVarPart-1+jPos)
    jpos=jpos+PP_nVarPart-1
    !>> TurbPt_temp for pushing: Runge-Kutta derivative of turbulent velocity fluctuation
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
      CALL Abort(__STAMP__,'Error with IsNewPart in MPIParticleRecv!',1,PartRecvBuf(iProc)%content(1+jPos))
    END IF
    jPos=jPos+1
    !>> particle elment
    PEM%Element(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1

    ! PartCommSize must be a multiple of particles to receive
    IF(MOD(jPos,PartCommSize).NE.0)THEN
      IPWRITE(UNIT_stdOut,'(A,I0)')  'jPos',jPos
      CALL Abort(__STAMP__,' Particle-wrong receiving message size!')
    END IF

    ! Set Flag for received parts in order to localize them later
    PDM%ParticleInside(PartID) = .TRUE.
    PEM%lastElement(PartID)    = -888
  END DO
END DO ! iProc


PDM%ParticleVecLength       = PDM%ParticleVecLength + PartMPIExchange%nMPIParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles
IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber)THEN
  WRITE(  UNIT_stdOut,'(A)')         ""
  IPWRITE(UNIT_stdOut,'(I0,A,I0)')   " PDM%ParticleVecLength = ", PDM%ParticleVecLength
  IPWRITE(UNIT_stdOut,'(I0,A,I0,A)') " PDM%MaxParticleNumber = ", PDM%MaxParticleNumber," for current processor"
  IPWRITE(UNIT_stdOut,'(I0,A)')      " Increase value for [Part-maxParticleNumber]!"
  CALL Abort(__STAMP__,' ParticleVecLegnth > MaxParticleNumber due to MPI-communication!')
END IF

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
  nInitRegions = nInitRegions + Species(iSpec)%NumberOfInits !+ (1-Species(iSpec)%StartnumberOfInits)
  ! old style parameters has been defined for inits/emissions but might have no particles
  IF (Species(iSpec)%Init(0)%UseForEmission) nInitRegions = nInitRegions + 1
END DO ! iSpec
IF(nInitRegions.GT.0) THEN
  DO iInitRegions=1,nInitRegions
    IF(PartMPI%InitGroup(iInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      IF(PartMPI%InitGroup(iInitRegions)%Comm.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(PartMPI%InitGroup(iInitRegions)%Comm,iERROR)
    END IF
  END DO ! iInitRegions
END IF
IF(PartMPI%COMM.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(PartMPI%COMM,iERROR)

SDEALLOCATE( PartMPIExchange%nPartsSend)
SDEALLOCATE( PartMPIExchange%nPartsRecv)
SDEALLOCATE( PartMPIExchange%RecvRequest)
SDEALLOCATE( PartMPIExchange%SendRequest)
SDEALLOCATE( PartMPIExchange%Send_message)
SDEALLOCATE( PartSendBuf)
SDEALLOCATE( PartRecvBuf)
SDEALLOCATE( ExchangeProcToGlobalProc)
SDEALLOCATE( GlobalProcToExchangeProc)

! and for communication
SDEALLOCATE( PartTargetProc )

ParticleMPIInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMPI
#endif /*MPI*/

END MODULE MOD_Particle_MPI
