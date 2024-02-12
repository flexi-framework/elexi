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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_MPI_Vars
! MODULES
#if USE_MPI
USE __MPI__
#endif /*USE_MPI*/
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!TYPE tPartExchange
INTEGER                                  :: nExchangeProcessors              ! number of MPI processes for particles exchange
INTEGER,ALLOCATABLE                      :: ExchangeProcToGlobalProc(:,:)    ! mapping from exchange proc ID to global proc ID
INTEGER,ALLOCATABLE                      :: GlobalProcToExchangeProc(:,:)    ! mapping from global proc ID to exchange proc ID

LOGICAL                                  :: ParticleMPIInitIsDone=.FALSE.
LOGICAL                                  :: CheckExchangeProcs               ! On default, check if proc communication is symmetric

TYPE tPartMPIGROUP
  MPI_TYPE_COMM                          :: COMM                             ! MPI communicator for PIC GTS region
  MPI_TYPE_REQUEST                       :: Request                          ! MPI request for asynchronous communication
  MPI_TYPE_REQUEST                       :: RequestIndex                     ! MPI request for asynchronous communication
  INTEGER                                :: nProcs                           ! number of MPI processes for particles
  INTEGER                                :: MyRank                           ! MyRank of PartMPIVAR%COMM
  LOGICAL                                :: MPIRoot                          ! Root, MPIRank=0
  INTEGER,ALLOCATABLE                    :: GroupToComm(:)                   ! list containing the rank in PartMPI%%InitGroup%COMM
  INTEGER,ALLOCATABLE                    :: CommToGroup(:)                   ! list containing the rank in PartMPI%%InitGroup%COMM
END TYPE

TYPE tPartMPIVAR
  TYPE(tPartMPIGROUP),ALLOCATABLE        :: InitGroup(:)                     ! small communicator for initialization
END TYPE

TYPE (tPartMPIVAR)                       :: PartMPI

REAL                                     :: SafetyFactor                     ! Factor to scale the halo region with MPI
REAL                                     :: halo_eps_velo                    ! halo_eps_velo
REAL                                     :: halo_eps                         ! length of halo-region
REAL                                     :: halo_eps2                        ! length of halo-region^2

LOGICAL                                  :: PartitionPartIsDone
#if USE_MPI
INTEGER                                  :: PartCommSize                     ! Number of REAL entries for particle communication
INTEGER                                  :: PartCommSize0                    ! Number of REAL entries for particle communication
                                                                             ! should think about own MPI-Data-Typ
TYPE tMPIMessage
  REAL,ALLOCATABLE                      :: content(:)                        ! message buffer real
END TYPE

TYPE(tMPIMessage),ALLOCATABLE  :: PartRecvBuf(:)                             ! PartRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: PartSendBuf(:)                             ! PartSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: SurfRecvBuf(:)                             ! SurfRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: SurfSendBuf(:)                             ! SurfSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: NodeRecvBuf(:)                             ! NodeRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: NodeSendBuf(:)                             ! NodeSendBuf with all requried types


TYPE(tMPIMessage),ALLOCATABLE  :: EmissionRecvBuf(:)                         ! EmissionRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: EmissionSendBuf(:)                         ! EmissionSendBuf with all requried types

TYPE tParticleMPIExchange
  INTEGER,ALLOCATABLE            :: nPartsSend(:,:)                          ! Only MPI neighbors
  INTEGER,ALLOCATABLE            :: nPartsRecv(:,:)                          ! Only MPI neighbors
  INTEGER                        :: nMPIParticles                            ! Number of all received particles
  MPI_TYPE_REQUEST,ALLOCATABLE   :: SendRequest(:,:)                         ! Send request message handle 1 - Number, 2-Message
  MPI_TYPE_REQUEST,ALLOCATABLE   :: RecvRequest(:,:)                         ! Receive request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)                          ! Message, required for particle emission
END TYPE
 TYPE (tParticleMPIExchange)     :: PartMPIExchange

TYPE tParticleMPIExchange2
  INTEGER,ALLOCATABLE            :: nPartsSend(:,:)                          ! Only same init communicator
  INTEGER,ALLOCATABLE            :: nPartsRecv(:,:)                          ! Only same init communicator
  INTEGER                        :: nMPIParticles                            ! Number of all received particles
  MPI_TYPE_REQUEST,ALLOCATABLE   :: SendRequest(:,:)                         ! Send requires message handle 1 - Number, 2-Message
  MPI_TYPE_REQUEST,ALLOCATABLE   :: RecvRequest(:,:)                         ! Receive request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)                          ! Message, required for particle emission
END TYPE
TYPE (tParticleMPIExchange2)     :: PartMPIInsert

TYPE tParticleMPIExchange3
  INTEGER,ALLOCATABLE            :: nPartsSend(:,:)                          ! Only same init communicator
  INTEGER,ALLOCATABLE            :: nPartsRecv(:,:)                          ! Only same init communicator
  INTEGER                        :: nMPIParticles                            ! Number of all received particles
  MPI_TYPE_REQUEST,ALLOCATABLE   :: SendRequest(:,:)                         ! Send requires message handle 1 - Number, 2-Message
  MPI_TYPE_REQUEST,ALLOCATABLE   :: RecvRequest(:,:)                         ! Receive request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)                          ! Message, required for particle emission
END TYPE
TYPE (tParticleMPIExchange3)     :: PartMPILocate

INTEGER,ALLOCATABLE              :: PartTargetProc(:)                        ! local proc id for communication
#endif /*USE_MPI*/
!===================================================================================================================================

END MODULE MOD_Particle_MPI_Vars
