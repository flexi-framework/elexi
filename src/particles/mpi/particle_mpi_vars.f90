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
MODULE MOD_Particle_MPI_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
INTEGER             :: iNbProc
INTEGER,ALLOCATABLE :: RecRequest_Flux(:),SendRequest_Flux(:)
INTEGER,ALLOCATABLE :: PartHaloElemToProc(:,:)                               ! containing native elemid and native proc id
                                                                             ! 1 - Native_Elem_ID
                                                                             ! 2 - Rank of Proc
                                                                             ! 3 - local neighbor id

INTEGER,ALLOCATABLE :: PartHaloSideToProc(:,:)                               ! containing native sideid and native proc id
                                                                             ! 1 - Native_Side_ID
                                                                             ! 3 - Rank of Proc
                                                                             ! 4 - local neighbor id


INTEGER             :: myRealKind
LOGICAL                                  :: ParticleMPIInitIsDone=.FALSE.
LOGICAL                                  :: printMPINeighborWarnings         ! print warning messages or not
LOGICAL                                  :: printBezierControlPointsWarnings ! print warning messages or not
INTEGER                                  :: iMessage                         ! Number of MPI-Messages for Debug purpose

TYPE tPartMPIGROUP
  INTEGER                                :: COMM                             ! MPI communicator for PIC GTS region
  INTEGER                                :: nProcs                           ! number of MPI processes for particles
  INTEGER                                :: MyRank                           ! MyRank of PartMPIVAR%COMM
  LOGICAL                                :: MPIRoot                          ! Root, MPIRank=0
  INTEGER,ALLOCATABLE                    :: GroupToComm(:)                   ! list containing the rank in PartMPI%COMM
  INTEGER,ALLOCATABLE                    :: CommToGroup(:)                   ! list containing the rank in PartMPI%COMM
END TYPE

TYPE tPeriodicPtr
  INTEGER                  , ALLOCATABLE  :: BGMPeriodicBorder(:,:)          ! indices of periodic border nodes
END TYPE

#if USE_MPI
TYPE tPartMPIConnect
  TYPE(tPeriodicPtr)       , ALLOCATABLE :: Periodic(:)                     ! data for different periodic borders for process
  LOGICAL                                :: isBGMNeighbor                    ! Flag: which process is neighber wrt. bckgrnd mesh
  LOGICAL                                :: isBGMPeriodicNeighbor            ! Flag: which process is neighber wrt. bckgrnd mesh
  INTEGER                  , ALLOCATABLE :: BGMBorder(:,:)            ! indices of border nodes (1=min 2=max,xyz)
  INTEGER                                :: BGMPeriodicBorderCount            ! Number(#) of overlapping areas due to periodic bc
END TYPE
#endif /*MPI*/


TYPE tPartMPIVAR
#if USE_MPI
  TYPE(tPartMPIConnect)        , ALLOCATABLE :: DepoBGMConnect(:)            ! MPI connect for each process
#endif /*MPI*/
  TYPE(tPartMPIGROUP),ALLOCATABLE        :: InitGroup(:)                     ! small communicator for initialization
  INTEGER                                :: COMM                             ! MPI communicator for PIC GTS region
  INTEGER                                :: nProcs                           ! number of MPI processes for particles
  INTEGER                                :: MyRank                           ! MyRank of PartMPIVAR%COMM
  LOGICAL                                :: MPIRoot                          ! Root, MPIRank=0
  INTEGER                                :: nMPINeighbors                    ! number of MPI-Neighbors with HALO
  LOGICAL,ALLOCATABLE                    :: isMPINeighbor(:)                 ! list of possible neighbors
  INTEGER,ALLOCATABLE                    :: MPINeighbor(:)                   ! list containing the rank of MPI-neighbors
  INTEGER,ALLOCATABLE                    :: GlobalToLocal(:)                 ! map from global proc id to local
END TYPE

TYPE (tPartMPIVAR)                       :: PartMPI

REAL                                     :: SafetyFactor                     ! Factor to scale the halo region with MPI
REAL                                     :: halo_eps_velo                    ! halo_eps_velo
REAL                                     :: halo_eps                         ! length of halo-region
REAL                                     :: halo_eps2                        ! length of halo-region^2

#if USE_MPI
INTEGER                                  :: PartCommSize                     ! Number of REAL entries for particle communication
INTEGER                                  :: PartCommSize0                    ! Number of REAL entries for particle communication
                                                                             ! should think about own MPI-Data-Typ

TYPE tMPIMessage
  REAL,ALLOCATABLE                      :: content(:)                        ! message buffer real
  LOGICAL,ALLOCATABLE                   :: content_log(:)                    ! message buffer logical for BGM
  INTEGER,ALLOCATABLE                   :: content_int(:)                    ! message buffer for integer for adsorption
END TYPE

TYPE(tMPIMessage),ALLOCATABLE  :: PartRecvBuf(:)                             ! PartRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: PartSendBuf(:)                             ! PartSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: SurfRecvBuf(:)                             ! SurfRecvBuf with all required types
TYPE(tMPIMessage),ALLOCATABLE  :: SurfSendBuf(:)                             ! SurfSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: SurfDistRecvBuf(:)                         ! SurfDistRecvBuf with all requried types
TYPE(tMPIMessage),ALLOCATABLE  :: SurfDistSendBuf(:)                         ! SurfDistSendBuf with all requried types

TYPE(tMPIMessage),ALLOCATABLE  :: AdsorbRecvBuf(:)
TYPE(tMPIMessage),ALLOCATABLE  :: AdsorbSendBuf(:)

TYPE(tMPIMessage),ALLOCATABLE  :: SurfCoverageRecvBuf(:)
TYPE(tMPIMessage),ALLOCATABLE  :: SurfCoverageSendBuf(:)

TYPE(tMPIMessage),ALLOCATABLE  :: CondensRecvBuf(:)
TYPE(tMPIMessage),ALLOCATABLE  :: CondensSendBuf(:)

TYPE tParticleMPIExchange
  INTEGER,ALLOCATABLE            :: nPartsSend(:,:)   ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: nPartsRecv(:,:)   ! only mpi neighbors
  INTEGER                        :: nMPIParticles     ! number of all received particles
  INTEGER,ALLOCATABLE            :: SendRequest(:,:)  ! send requirest message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: RecvRequest(:,:)  ! recv request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)   ! message, required for particle emission
END TYPE
 TYPE (tParticleMPIExchange)     :: PartMPIExchange

TYPE tParticleMPIExchange2
  INTEGER,ALLOCATABLE            :: nPartsSend(:)     ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: nPartsRecv(:)     ! only mpi neighbors
  INTEGER                        :: nMPIParticles     ! number of all received particles
  INTEGER,ALLOCATABLE            :: SendRequest(:,:)  ! send requirest message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: RecvRequest(:,:)  ! recv request message handle,  1 - Number, 2-Message
  TYPE(tMPIMessage),ALLOCATABLE  :: send_message(:)   ! message, required for particle emission
END TYPE
TYPE (tParticleMPIExchange2)     :: PartMPIInsert

TYPE tDistNbrComm
  INTEGER,ALLOCATABLE            :: nPosSend(:)   ! number of distribution site to send for surface (nSidesSend,nCoordination)
  INTEGER,ALLOCATABLE            :: nPosRecv(:)   ! number of distribution sites to receive for surface (nSidesRecv,nCoordination)
END TYPE

TYPE tSurfMPIExchange
  INTEGER,ALLOCATABLE            :: nSidesSend(:)     ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: nSidesRecv(:)     ! only mpi neighbors
  INTEGER,ALLOCATABLE            :: SendRequest(:)   ! send requirest message handle 1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: RecvRequest(:)   ! recv request message handle,  1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: nSurfDistSidesSend(:)        ! number of mpi sides to send (nProcs)
  INTEGER,ALLOCATABLE            :: nSurfDistSidesRecv(:)        ! number of sides received from mpi (nProcs)
  INTEGER,ALLOCATABLE            :: nCoverageSidesSend(:)        ! number of mpi sides to send (nProcs)
  INTEGER,ALLOCATABLE            :: nCoverageSidesRecv(:)        ! number of sides received from mpi (nProcs)
  TYPE(tDistNbrComm),ALLOCATABLE :: NbrOfPos(:)          ! array for number of distribution sites sending per proc
  INTEGER,ALLOCATABLE            :: SurfDistSendRequest(:,:)     ! send request message handle,  1 - Number, 2-Message
  INTEGER,ALLOCATABLE            :: SurfDistRecvRequest(:,:)     ! recv request message handle,  1 - Number, 2-Message
END TYPE
TYPE (tSurfMPIExchange)          :: SurfExchange


INTEGER,ALLOCATABLE                      :: PartTargetProc(:)                ! local proc id for communication

REAL, ALLOCATABLE                        :: PartShiftVector(:,:)             ! store particle periodic map
#endif /*MPI*/
!===================================================================================================================================

END MODULE MOD_Particle_MPI_Vars
