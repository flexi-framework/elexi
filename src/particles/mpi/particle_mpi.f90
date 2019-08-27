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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

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

INTERFACE BoxInProc
  MODULE PROCEDURE BoxInProc
END INTERFACE

INTERFACE ExchangeBezierControlPoints3D
  MODULE PROCEDURE ExchangeBezierControlPoints3D
END INTERFACE

PUBLIC :: InitParticleMPI,FinalizeParticleMPI,InitHaloMesh, InitParticleCommSize, IRecvNbOfParticles, MPIParticleSend
PUBLIC :: MPIParticleRecv
PUBLIC :: InitEmissionComm
PUBLIC :: ExchangeBezierControlPoints3D
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
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                            :: myRealTestValue
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
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,PartMPI%MyRank,PartMPI%COMM,iERROR)
CALL MPI_COMM_SIZE (PartMPI%COMM,PartMPI%nProcs ,iError)

IF(PartMPI%nProcs.NE.nProcessors) CALL abort(&
    __STAMP__&
    ,' MPI Communicater-size does not match!', IERROR)
PartCommSize   = 0  
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
!! determine datatype length for variables to be sent
!myRealKind = KIND(myRealTestValue)
!IF (myRealKind.EQ.4) THEN
!  myRealKind = MPI_REAL
!ELSE IF (myRealKind.EQ.8) THEN
!  myRealKind = MPI_DOUBLE_PRECISION
!ELSE
!  myRealKind = MPI_REAL
!END IF

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
USE MOD_Particle_MPI_Vars
USE MOD_Particle_Vars,          ONLY:PDM
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Erosion_Vars,  ONLY:PartTrackReflection
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: ALLOCSTAT
!===================================================================================================================================

PartCommSize   = 0  
! PartState: position and velocity
PartCommSize   = PartCommSize + 6
! Tracking: Include Reference coordinates
IF(DoRefMapping) PartCommSize=PartCommSize+3
! Species-ID
PartCommSize   = PartCommSize + 1
! id of element
PartCommSize   = PartCommSize + 1

! time integration
! Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
PartCommSize   = PartCommSize + 6
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
        , PartMPIDepoSend(1:PDM%MaxParticleNumber)             &
        , STAT=ALLOCSTAT                                       )

IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)

PartMPIExchange%nPartsSend=0
PartMPIExchange%nPartsRecv=0

IF(DoExternalParts)THEN
  ExtPartCommSize=7
END IF

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
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iProc
!===================================================================================================================================

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
 ! IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__&
 !         ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


SUBROUTINE SendNbOfParticles(doParticle_In)
!===================================================================================================================================
! this routine sends the number of send particles. Following steps are performed
! 1) Compute number of Send Particles
! 2) Performe MPI_ISEND with number of particles
! Rest is perforemd in SendParticles
! 3) Build Message 
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartHaloElemToProc, PartTargetProc
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,PEM,PDM,Species,PartPosRef
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! variables for parallel deposition
USE MOD_Particle_MPI_Vars,        ONLY:DoExternalParts,PartMPIDepoSend
USE MOD_Particle_MPI_Vars,        ONLY:PartShiftVector
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
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
! shape function 
INTEGER                       :: CellX,CellY,CellZ
INTEGER                       :: PartDepoProcs(1:PartMPI%nProcs+1), nDepoProcs, ProcID,LocalProcID
INTEGER                       :: nPartShape
REAL                          :: ShiftedPart(1:3)
LOGICAL                       :: PartInBGM
!===================================================================================================================================

IF(PRESENT(DoParticle_IN))THEN
  DoParticle=PDM%ParticleInside(1:PDM%ParticleVecLength).AND.DoParticle_In
ELSE
  DoParticle(1:PDM%ParticleVecLength)=PDM%ParticleInside(1:PDM%ParticleVecLength)
END IF

! 1) get number of send particles
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

! external particles for communication
IF(DoExternalParts)THEN
  PartMPIDepoSend=.FALSE.
  nPartShape=0
  DO iPart=1,PDM%ParticleVecLength
    !IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    IF(.NOT.DoParticle(iPart)) CYCLE
    IF (Species(PartSpecies(iPart))%ChargeIC.EQ.0) CYCLE        ! Don't deposite neutral particles!
    CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
    CellX = MIN(GEO%FIBGMimax,CellX)
    CellX = MAX(GEO%FIBGMimin,CellX)
    CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
    CellY = MIN(GEO%FIBGMjmax,CellY)
    CellY = MAX(GEO%FIBGMjmin,CellY)
    CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
    CellZ = MIN(GEO%FIBGMkmax,CellZ)
    CellZ = MAX(GEO%FIBGMkmin,CellZ)
    IF(ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
      IF(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1) .GT. 0) THEN
        nPartShape=nPartShape+1
        PartMPIDepoSend(iPart) = .TRUE.
      END IF
    END IF
  END DO ! iPart=1,PDM%ParticleVecLength
  ! now, get correct BGM cell for particle 
  ! including periodic displacement or BGM element without mpi neighbors
  ! shape-padding could be modified for all other deposition methods? reuse?
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartMPIDepoSend(iPart)) CYCLE
    IF (Species(PartSpecies(iPart))%ChargeIC.EQ.0) CYCLE ! get BMG cell
    CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
    CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
    CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
    PartInBGM = .TRUE.
    ! check if particle is in range of my FIBGM 
    ! first check is outside
    IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
        (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
        (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
      PartInBGM = .FALSE.
    ELSE
      ! particle inside, check if particle is not moved by periodic BC
      ! if this is the case, then the ShapeProcs is not allocated!
      IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
      END IF
    END IF
    IF (.NOT.PartInBGM) THEN
      ! it is possible that the particle has been moved over a periodic side 
      IF (GEO%nPeriodicVectors.GT.0) THEN
        ShiftedPart(1:3) = PartState(iPart,1:3) + partShiftVector(1:3,iPart)
        CellX = INT((ShiftedPart(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
        CellY = INT((ShiftedPart(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
        CellZ = INT((ShiftedPart(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
        IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
            (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
            (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
          CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          CellX = MIN(GEO%FIBGMimax,CellX)
          CellX = MAX(GEO%FIBGMimin,CellX)
          CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          CellY = MIN(GEO%FIBGMjmax,CellY)
          CellY = MAX(GEO%FIBGMjmin,CellY)
          CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          CellZ = MIN(GEO%FIBGMkmax,CellZ)
          CellZ = MAX(GEO%FIBGMkmin,CellZ)
        ELSE
          IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
            IPWRITE(UNIT_errOut,*)'ERROR in SendNbOfParticles: Particle outside BGM! Err2'
            IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',DoParticle(iPart)
            IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
            IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
            IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
            IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',ShiftedPart(1),',PartY=',ShiftedPart(2),',PartZ=',&
                    ShiftedPart(3)
            IF(DoRefMapping)THEN
              IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartXi=',PartPosRef(1,iPart)   &
                                                    ,',PartEta=',PartPosRef(2,iPart) &
                                                    ,',PartZeta=',PartPosRef(3,iPart)
            END IF
            CALL Abort(&
            __STAMP__&
            ,'Particle outside BGM! Err2')
          END IF
        END IF
      ELSE
        IPWRITE(UNIT_errOut,*)'Warning in SendNbOfParticles: Particle outside BGM!'
        IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',DoParticle(iPart)
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
        IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',PartState(iPart,1),',PartY=',PartState(iPart,2),',PartZ=',&
                PartState(iPart,3)
        IF(DoRefMapping)THEN
          IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartXi=',PartPosRef(1,iPart),',PartEta=',PartPosRef(2,iPart),',PartZeta=',&
                  PartPosRef(3,iPart)
        END IF
        IPWRITE(UNIT_errOut,*)'Remap particle!'

        CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
        CellX = MIN(GEO%FIBGMimax,CellX)
        CellX = MAX(GEO%FIBGMimin,CellX)
        CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
        CellY = MIN(GEO%FIBGMjmax,CellY)
        CellY = MAX(GEO%FIBGMjmin,CellY)
        CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
        CellZ = MIN(GEO%FIBGMkmax,CellZ)
        CellZ = MAX(GEO%FIBGMkmin,CellZ)
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'New-CellX=',CellX,',New-CellY=',CellY,',New-CellZ=',CellZ
        ! nothing to do, because of tolerance, particle could be outside
        !IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
        !    (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
        !    (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
 
        !  CALL Abort(&
        !       __STAMP__&
        !      'Particle outside BGM!')
        !END IF
      END IF ! GEO%nPeriodicVectors
    END IF ! PartInBGM
    nDepoProcs=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)
    PartDepoProcs(1:nDepoProcs)=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(2:nDepoProcs+1)
    DO iProc=1,nDepoProcs
      ProcID=PartDepoProcs(iProc)
      ! particle shall not be send to MyRank, is fixed without MPI communication in MPIParticleSend
      IF(ProcID.EQ.PartMPI%MyRank) CYCLE
      ! if shapeproc is target proc, to net send
      ! short version
      LocalProcID=PartMPI%GlobalToLocal(ProcID)
      IF(PartTargetProc(iPart).EQ.LocalProcID) CYCLE
      PartMPIExchange%nPartsSend(2,LocalProcID)= PartMPIExchange%nPartsSend(2,LocalProcID) +1
    END DO ! iProc=1,nDepoProcs
  END DO ! iPart=1,PDM%ParticleVecLength
END IF


! 2) send number of send particles
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
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,PEM,PDM,Species,PartPosRef
USE MOD_Particle_Vars,            ONLY:Pt_temp
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! variables for parallel deposition
USE MOD_Particle_MPI_Vars,        ONLY:DoExternalParts,PartMPIDepoSend,PartShiftVector, ExtPartCommSize, PartMPIDepoSend
USE MOD_Particle_MPI_Vars,        ONLY:ExtPartState,ExtPartSpecies,NbrOfExtParticles
! variables for erosion tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Erosion_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID,iPos,iProc,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize, nRecvParticles, nSendParticles, nSendExtParticles, nRecvExtParticles
INTEGER                       :: ALLOCSTAT
! shape function 
INTEGER                       :: CellX,CellY,CellZ
INTEGER                       :: PartDepoProcs(1:PartMPI%nProcs+1), nDepoProcs, ProcID, jProc, iExtPart, LocalProcID
REAL                          :: ShiftedPart(1:3)
LOGICAL                       :: PartInBGM
! Polyatomic Molecules
INTEGER                       :: MsgLengthPoly(1:PartMPI%nMPINeighbors)
!===================================================================================================================================

!--- Determining the number of additional variables due to VibQuantsPar of polyatomic particles
!--- (size varies depending on the species of particle)
MsgLengthPoly(:) = 0

! 3) Build Message
DO iProc=1, PartMPI%nMPINeighbors
  ! allocate SendBuf
  nSendParticles=PartMPIExchange%nPartsSend(1,iProc)
  iPos=0
  IF(DoExternalParts)THEN
    nSendExtParticles=PartMPIExchange%nPartsSend(2,iProc)
    IF((nSendExtParticles.EQ.0).AND.(nSendParticles.EQ.0)) CYCLE
    MessageSize=nSendParticles*PartCommSize &
               +nSendExtParticles*ExtPartCommSize
  ELSE
    IF(nSendParticles.EQ.0) CYCLE
    MessageSize=nSendParticles*PartCommSize
  END IF
  
  ALLOCATE(PartSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'  Cannot allocate PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))
  ! fill message
  DO iPart=1,PDM%ParticleVecLength
    IF(nSendParticles.EQ.0) EXIT
    IF(PartTargetProc(iPart).EQ.iProc) THEN
      !iPos=iPos+1
        
      ! fill content
      ElemID=PEM%Element(iPart)
      PartSendBuf(iProc)%content(1+iPos:6+iPos) = PartState(iPart,1:6)
      IF(DoRefMapping) THEN ! + deposition type....
        PartSendBuf(iProc)%content(7+iPos:9+iPos) = PartPosRef(1:3,iPart)
        jPos=iPos+9
      ELSE
        jPos=iPos+6
      END IF
      
      ! Send reflection counter
      IF (PartTrackReflection) THEN
          IF(DoRefMapping) THEN
              PartSendBuf(iProc)%content(10+iPos) = REAL(PartReflCount(iPart),KIND=8)
          ELSE
              PartSendBuf(iProc)%content(7+iPos)  = REAL(PartReflCount(iPart),KIND=8)
          ENDIF
      ! Now update send buffer position
          IF(DoRefMapping) THEN
              jPos=iPos+10
          ELSE
              jPos=iPos+7
          END IF
      END IF
      
      PartSendBuf(iProc)%content(       1+jPos) = REAL(PartSpecies(iPart),KIND=8)
      jPos=jPos+1
      PartSendBuf(iProc)%content(1+jPos:6+jPos) = Pt_temp(iPart,1:6)
      IF (PDM%IsNewPart(iPart)) THEN
        PartSendBuf(iProc)%content(7+jPos) = 1.
      ELSE
        PartSendBuf(iProc)%content(7+jPos) = 0.
      END IF
      jPos=jPos+7
      PartSendBuf(iProc)%content(    1+jPos) = REAL(PartHaloElemToProc(NATIVE_ELEM_ID,ElemID),KIND=8)
      jPos=jPos+1
      !--- put the polyatomic vibquants per particle at the end of the message
      IF(MOD(jPos,PartCommSize).NE.0) THEN
        IPWRITE(UNIT_stdOut,*)  'PartCommSize',PartCommSize
        IPWRITE(UNIT_stdOut,*)  'jPos',jPos
        CALL Abort(&
            __STAMP__&
            ,' Particle-wrong sending message size!')
      END IF
      ! endif timedisc
      ! here iPos because PartCommSize contains DoRefMapping
      iPos=iPos+PartCommSize
      ! particle is ready for send, now it can deleted
      PDM%ParticleInside(iPart) = .FALSE.  
    END IF ! Particle is particle with target proc-id equals local proc id
  END DO  ! iPart
  ! next, external particles has to be handled for deposition
  IF(DoExternalParts)THEN
    IF(nSendExtParticles.EQ.0) CYCLE
    iPos=nSendParticles*PartCommSize
    DO iPart=1,PDM%ParticleVecLength
      IF(PartTargetProc(iPart).EQ.iProc) CYCLE
      IF(.NOT.PartMPIDepoSend(iPart)) CYCLE
      IF (Species(PartSpecies(iPart))%ChargeIC.EQ.0) CYCLE
      ! get BMG cell
      CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
      CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
      CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
      PartInBGM = .TRUE.
      ! check if particle is in range of my FIBGM 
      ! first check is outside
      IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
          (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
          (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
        PartInBGM = .FALSE.
      ELSE
        ! particle inside, check if particle is not moved by periodic BC
        ! if this is the case, then the ShapeProcs is not allocated!
        IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
          PartInBGM = .FALSE.
        END IF
      END IF
      IF (.NOT.PartInBGM) THEN
        ! it is possible that the particle has been moved over a periodic side 
        IF (GEO%nPeriodicVectors.GT.0) THEN
          ShiftedPart(1:3) = PartState(iPart,1:3) + partShiftVector(1:3,iPart)
          CellX = INT((ShiftedPart(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          CellY = INT((ShiftedPart(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          CellZ = INT((ShiftedPart(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
              (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
              (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
            CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
            CellX = MIN(GEO%FIBGMimax,CellX)
            CellX = MAX(GEO%FIBGMimin,CellX)
            CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
            CellY = MIN(GEO%FIBGMjmax,CellY)
            CellY = MAX(GEO%FIBGMjmin,CellY)
            CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
            CellZ = MIN(GEO%FIBGMkmax,CellZ)
            CellZ = MAX(GEO%FIBGMkmin,CellZ)
          ELSE
            IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
              IPWRITE(UNIT_errOut,*)'ERROR in SendNbOfParticles: Particle outside BGM! Err2'
              IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',PDM%ParticleInside(iPart)
              IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
              IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
              IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
              IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',ShiftedPart(1),',PartY=',ShiftedPart(2),',PartZ=',&
                      ShiftedPart(3)
              CALL Abort(&
              __STAMP__&
              ,'Particle outside BGM! Err2')
            END IF
          END IF
        ELSE
          IPWRITE(UNIT_errOut,*)'Remap particle!'

          CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          CellX = MIN(GEO%FIBGMimax,CellX)
          CellX = MAX(GEO%FIBGMimin,CellX)
          CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          CellY = MIN(GEO%FIBGMjmax,CellY)
          CellY = MAX(GEO%FIBGMjmin,CellY)
          CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          CellZ = MIN(GEO%FIBGMkmax,CellZ)
          CellZ = MAX(GEO%FIBGMkmin,CellZ)
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'New-CellX=',CellX,',New-CellY=',CellY,',New-CellZ=',CellZ
          IPWRITE(UNIT_errOut,*)'ERROR in SendNbOfParticles: Particle outside BGM!'
          IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',PDM%ParticleInside(iPart)
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
          IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',PartState(iPart,1),',PartY=',PartState(iPart,2),',PartZ=',&
                  PartState(iPart,3)
          !CALL Abort(&
          !     __STAMP__&
          !    'Particle outside BGM!')
        END IF ! GEO%nPeriodicVectors
      END IF ! PartInBGM
      nDepoProcs=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)
      PartDepoProcs(1:nDepoProcs)=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(2:nDepoProcs+1)
      DO jProc=1,nDepoProcs
        ProcID=PartDepoProcs(jProc)
        LocalProcID=PartMPI%GlobalToLocal(ProcID)
        IF(PartTargetProc(iPart).EQ.LocalProcID) CYCLE
        IF(LocalProcID.NE.iProc) CYCLE
        PartSendBuf(iProc)%content(1+iPos:6+iPos) = PartState(iPart,1:6)
        PartSendBuf(iProc)%content(       7+iPos) = REAL(PartSpecies(iPart),KIND=8)
        ! count only, if particle is sent
        iPos=iPos+ExtPartCommSize
      END DO ! jProc=1,nDepoProcs
    END DO ! iPart=1,PDM%ParticleVecLength 
  END IF ! DoExternalParts
  IF(iPos.NE.(MessageSize-MsgLengthPoly(iProc))) IPWRITE(*,*) ' error message size', iPos,(MessageSize-MsgLengthPoly(iProc))
END DO ! iProc

! 4) Finish Received number of particles
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


! caution, fancy trick, particles are sent, but information is not deleted
! temporary storage
IF(DoExternalParts) THEN
  NbrOfExtParticles =SUM(PartMPIExchange%nPartsSend(1,:))+SUM(PartMPIExchange%nPartsRecv(2,:))
  ALLOCATE(ExtPartState  (1:NbrOfExtParticles,1:6) &
          ,ExtPartSpecies(1:NbrOfExtParticles)     &
          !,ExtPartToFIBGM(1:6,1:NbrOfExtParticles) &
          ,STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'  Cannot allocate ExtPartState on Rank',PartMPI%MyRank,REAL(ALLOCSTAT))
  ! map alt state to ext
  iExtPart=0
  DO iPart=1,PDM%ParticleVecLength
    IF(PartTargetProc(iPart).EQ.-1) CYCLE
    iExtPart=iExtPart+1
    ExtPartState(iExtPart,1:6)        = PartState(iPart,1:6)
    ExtPartSpecies(iExtPart)          = PartSpecies(iPart)
  END DO ! iPart=1,PDM%ParticleVecLength 
END IF

DO iPart=1,PDM%ParticleVecLength
  IF(PartTargetProc(iPart).EQ.-1) CYCLE
  PartState(iPart,1:6)=0.
  PartSpecies(iPart)=0
  Pt_temp(iPart,1:6)=0.
END DO ! iPart=1,PDM%ParticleVecLength 


! 5) Allocate received buffer and open MPI_IRECV
DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE
  nRecvParticles=PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize=nRecvParticles*PartCommSize
  IF(DoExternalParts) THEN
    nRecvExtParticles=PartMPIExchange%nPartsRecv(2,iProc)
    MessageSize=MessageSize   &
               +nRecvExtParticles*ExtPartCommSize
  END IF
  ! determine the maximal possible polyatomic addition to the regular recv message
  ALLOCATE(PartRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    IPWRITE(*,*) 'sum of total received particles            ', SUM(PartMPIExchange%nPartsRecv(1,:))
    IPWRITE(*,*) 'sum of total received deposition particles ', SUM(PartMPIExchange%nPartsRecv(2,:))
    CALL abort(&
    __STAMP__&
    ,'  Cannot allocate PartRecvBuf, local source ProcId, Allocstat',iProc,REAL(ALLOCSTAT))
  END IF
  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)

END DO ! iProc

! 6) Send Particles
DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE
  nSendParticles=PartMPIExchange%nPartsSend(1,iProc)
  MessageSize=nSendParticles*PartCommSize
  IF(DoExternalParts)THEN
    nSendExtParticles=PartMPIExchange%nPartsSend(2,iProc)
    MessageSize=MessageSize &
               +nSendExtParticles*ExtPartCommSize
  END IF

  CALL MPI_ISEND( PartSendBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)

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
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartCommSize, PartRecvBuf,PartSendBuf
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,PEM,PDM,PartPosRef
USE MOD_Particle_Vars,            ONLY:Pt_temp
! variables for parallel deposition
USE MOD_Particle_MPI_Vars,        ONLY:DoExternalParts,ExtPartCommSize
USE MOD_Particle_MPI_Vars,        ONLY:ExtPartState,ExtPartSpecies
! variables for erosion tracking
USE MOD_Particle_Vars,            ONLY:PartReflCount
USE MOD_Particle_Erosion_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc, iPos, nRecv, PartID,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize, nRecvParticles, nRecvExtParticles
! shape function 
INTEGER                       :: iExtPart
! Polyatomic Molecules
INTEGER                       :: MsgLengthPoly
!===================================================================================================================================

DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE
  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

! old number of already filled ExtParticles
IF(DoExternalParts) iExtPart=SUM(PartMPIExchange%nPartsSend(1,:))

nRecv=0
DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE
  nRecvParticles=PartMPIExchange%nPartsRecv(1,iProc)
  IF(DoExternalParts) THEN
    nRecvExtParticles=PartMPIExchange%nPartsRecv(2,iProc)
  END IF
  ! determine the maximal possible polyatomic addition to the regular message
  MsgLengthPoly = 0
  !IF(nRecvParticles.EQ.0) CYCLE
  MessageSize=nRecvParticles*PartCommSize
  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
  ! correct loop shape
  ! DO iPart=1,nRecvParticles
  ! nParts 1 Pos=1..17 
  ! nPart2 2 Pos=1..17,18..34
  DO iPos=0,MessageSize-1-MsgLengthPoly,PartCommSize
    IF(nRecvParticles.EQ.0) EXIT
    nRecv=nRecv+1
    PartID = PDM%nextFreePosition(nRecv+PDM%CurrentNextFreePosition)
    IF(PartID.EQ.0)  CALL abort(&
      __STAMP__&
      ,' Error in ParticleExchange_parallel. Corrupted list: PIC%nextFreePosition', nRecv)
    PartState(PartID,1:6)   = PartRecvBuf(iProc)%content( 1+iPos: 6+iPos)
    IF(DoRefMapping) THEN
      PartPosRef(1:3,PartID) = PartRecvBuf(iProc)%content(7+iPos: 9+iPos)
      jPos=iPos+9
    ELSE
      jPos=iPos+6
    END IF
      
    ! Receive reflection counter
      IF (PartTrackReflection) THEN
          IF(DoRefMapping)THEN
              PartReflCount(PartID)   = INT(PartRecvBuf(iProc)%content(10+iPos),KIND=4)
          ELSE
              PartReflCount(PartID)   = INT(PartRecvBuf(iProc)%content(7+iPos),KIND=4)
          END IF
    ! Now update receive buffer position
          IF(DoRefMapping) THEN
              jPos=iPos+10
          ELSE
              jPos=iPos+7
          END IF
      END IF
      
    PartSpecies(PartID)     = INT(PartRecvBuf(iProc)%content( 1+jPos),KIND=4)
    jPos=jPos+1
    
    Pt_temp(PartID,1:6)     = PartRecvBuf(iProc)%content( 1+jPos:6+jPos)
    IF ( INT(PartRecvBuf(iProc)%content( 7+jPos)) .EQ. 1) THEN
      PDM%IsNewPart(PartID)=.TRUE.
    ELSE IF ( INT(PartRecvBuf(iProc)%content( 7+jPos)) .EQ. 0) THEN
      PDM%IsNewPart(PartID)=.FALSE.
    ELSE
      CALL Abort(&
        __STAMP__&
        ,'Error with IsNewPart in MPIParticleRecv!')
    END IF
    jPos=jPos+7

    PEM%Element(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1
    IF(MOD(jPos,PartCommSize).NE.0)THEN
      IPWRITE(UNIT_stdOut,*)  'jPos',jPos
      CALL Abort(&
          __STAMP__&
          ,' Particle-wrong receiving message size!')
    END IF
    ! Set Flag for received parts in order to localize them later
    PDM%ParticleInside(PartID) = .TRUE.
    PEM%lastElement(PartID) = -888 
  END DO
  IF(DoExternalParts)THEN
    jPos=MessageSize
    IF(nRecvExtParticles.EQ.0) CYCLE
    MessageSize=nRecvExtParticles*ExtPartCommSize+jPos
    DO iPos=jPos,MessageSize-1,ExtPartCommSize
      iExtPart=iExtPart+1
      ExtPartState(iExtPart,1:6) = PartRecvBuf(iProc)%content( 1+iPos: 6+iPos)
      ExtPartSpecies(iExtPart)   = INT(PartRecvBuf(iProc)%content( 7+iPos),KIND=4)
    END DO ! iPos
  END IF ! DoExternalParts
  ! be nice: deallocate the receive buffer
  ! deallocate non used array
END DO ! iProc


PDM%ParticleVecLength       = PDM%ParticleVecLength + PartMPIExchange%nMPIParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles
IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) CALL abort(&
    __STAMP__&
    ,' ParticleVecLegnth>MaxParticleNumber due to MPI-communication!')

! deallocate send,receive buffer
DO iProc=1,PartMPI%nMPINeighbors
  SDEALLOCATE(PartRecvBuf(iProc)%content)
  SDEALLOCATE(PartSendBuf(iProc)%content)
END DO ! iProc


! final
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
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
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
SDEALLOCATE( ExtPartState)
SDEALLOCATE( ExtPartSpecies)

! and for communication
SDEALLOCATE( PartTargetProc )
SDEALLOCATE( PartMPIDepoSend )

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

!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
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
  IF (PartMPI%MyRank.LT.iProc) THEN
    CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
    CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  ELSE IF (PartMPI%MyRank.GT.iProc) THEN
    CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
    CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
  END IF
  !IPWRITE(UNIT_stdOut,*) 'check',tmpneigh,PartMPI%isMPINeighbor(iProc)
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
!USE MOD_Particle_Mesh,          ONLY:BoxInProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec,iInit,iNode,iRank
INTEGER                         :: nInitRegions
LOGICAL                         :: RegionOnProc
REAL                            :: xCoords(3,8),lineVector(3),radius,height
REAL                            :: xlen,ylen,zlen
REAL                            :: dt
INTEGER                         :: color,iProc
INTEGER                         :: noInitRank,InitRank
!INTEGER,ALLOCATABLE             :: DummyRank(:)
LOGICAL                         :: hasRegion

! partns
INTEGER                         :: dummy_errtype
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
    CASE('gyrotron_circle')
      Radius=Species(iSpec)%Init(iInit)%RadiusIC+Species(iSpec)%Init(iInit)%RadiusICGyro
      xlen=Radius
      ylen=Radius
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      IF(Species(iSpec)%Init(iInit)%initialParticleNumber.NE.0)THEN
        lineVector(1:3)=(/0.,0.,Species(iSpec)%Init(iInit)%CuboidHeightIC/)
      ELSE
!#ifndef PP_HDG
        dt = CALCTIMESTEP(dummy_errtype)
!#endif /*PP_HDG*/
        lineVector(1:3)= dt* Species(iSpec)%Init(iInit)%VeloIC/Species(iSpec)%Init(iInit)%alpha 
        zlen=0.
      END IF
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,+ylen,+zlen/)
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
    CASE('cuboid_vpi')

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

      height = halo_eps
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cylinder_vpi')
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

      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC-radius*Species(iSpec)%Init(iInit)%BaseVector1IC &
                                                           -radius*Species(iSpec)%Init(iInit)%BaseVector2IC

      xCoords(1:3,2)=xCoords(1:3,1)+2.0*radius*Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=xCoords(1:3,1)+2.0*radius*Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=xCoords(1:3,1)+2.0*radius*Species(iSpec)%Init(iInit)%BaseVector1IC&
                                   +2.0*radius*Species(iSpec)%Init(iInit)%BaseVector2IC

      height = halo_eps
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)


    CASE('LD_insert')
      RegionOnProc=.TRUE.
    CASE('cell_local')
      RegionOnProc=.TRUE.
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
    CASE ('IMD')
       RegionOnProc=.TRUE.
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
    WRITE(UNIT_StdOut,*) ' Emission-Region,Emission-Communicator:',nInitRegions,PartMPI%InitGroup(nInitRegions)%nProcs,' procs'
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
#endif /*MPI*/

END MODULE MOD_Particle_MPI