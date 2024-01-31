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
#include "particle.h"

!===================================================================================================================================
! Module containing general particles collisions routines
!===================================================================================================================================
MODULE MOD_Particle_Collision
! MODULES
IMPLICIT NONE
PRIVATE

#if PARTICLES_COUPLING == 4
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersCollission
  MODULE PROCEDURE DefineParametersCollission
END INTERFACE

INTERFACE InitializeCollision
  MODULE PROCEDURE InitializeCollision
END INTERFACE

INTERFACE UpdateParticleShared
  MODULE PROCEDURE UpdateParticleShared
END INTERFACE

INTERFACE FinalizeCollision
  MODULE PROCEDURE FinalizeCollision
END INTERFACE

PUBLIC :: DefineParametersCollission
PUBLIC :: InitializeCollision
PUBLIC :: UpdateParticleShared
PUBLIC :: FinalizeCollision
!==================================================================================================================================

CONTAINS

SUBROUTINE DefineParametersCollission
!===================================================================================================================================
!> Initialize particle collision module
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools                ,ONLY: prms
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL prms%SetSection('Particle Collisions')

CALL prms%CreateLogicalOption('Part-CalcCollision'        , 'Flag to enable four-way coupling'                               &
                                                          , '.FALSE.')
CALL prms%CreateRealOption(   'Part-Collisions-RestCoeff' , 'Resitution coefficient (in [0;1]).'                             &
                                                          , '1.' )
CALL prms%CreateLogicalOption('Part-Collisions-Friction'  , 'Enable or disable friction tangential.'                         &
                                                          , '.FALSE.' )
CALL prms%CreateRealOption(   'Part-Collisions-FricCoeff' , 'Friction coefficient.'                                           &
                                                          , '0.3' )

END SUBROUTINE DefineParametersCollission


SUBROUTINE InitializeCollision()
!===================================================================================================================================
!> Initialize particle collisions module
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: nGlobalElems
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_ReadInTools             ,ONLY: GETLOGICAL,GETREAL
USE MOD_Particle_Collision_Vars
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: PartIntSize=2
REAL                           :: StartT,EndT
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE COLLISIONS...'
GETTIME(StartT)

! initialize the parameters of the model particle collision
PartCollisionModel%e        = GETREAL   ('Part-Collisions-RestCoeff')
IF (PartCollisionModel%e.LT.0. .OR. PartCollisionModel%e.GT.1.) CALL Abort(__STAMP__,'Part-Collisions-RestCoeff is not between 0 and 1!')
PartCollisionModel%Friction = GETLOGICAL('Part-Collisions-Friction')
PartCollisionModel%f        = GETREAL   ('Part-Collisions-FricCoeff')

! determine neighboring elements
CALL IdentifyNeighboringElements()

! allocate shared particle arrays
CALL Allocate_Shared((/PartIntSize,nGlobalElems/),PartInt_Shared_Win ,PartInt_Shared)
CALL MPI_WIN_LOCK_ALL(0,PartInt_Shared_Win,iError)
IF (myComputeNodeRank.EQ.0) PartInt_Shared = -1
CALL BARRIER_AND_SYNC(PartInt_Shared_Win,MPI_COMM_SHARED)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'INIT PARTICLE COLLISIONS DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.TRUE.)

END SUBROUTINE InitializeCollision


SUBROUTINE IdentifyNeighboringElements()
!===================================================================================================================================
!> Initialize mesh neighboring elements for neighbor search
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals        ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars      ,ONLY: BoundsOfElem_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: FIBGM_nElems,FIBGM_Element,FIBGM_offsetElem
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_MPI_Vars       ,ONLY: halo_eps
USE MOD_Particle_Vars           ,ONLY: doCalcPartCollision
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeTotalElems
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem,nGlobalElems
USE MOD_Mesh_Vars               ,ONLY: nComputeNodeElems
USE MOD_Utils                   ,ONLY: InsertionSort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: firstElem,lastElem
INTEGER                        :: iElem,iLocElem,ElemID
! Distances
REAL                           :: BoundsOfElemCenter(     1:4)
REAL                           :: LocalBoundsOfElemCenter(1:4)
REAL                           :: OrigBoundsOfElemCenter( 1:4)
REAL                           :: dist
! Periodic
INTEGER                        :: iDir,jDir,kDir
INTEGER                        :: iPeriodic,jPeriodic,kPeriodic
! Mappings
LOGICAL,ALLOCATABLE            :: currentNeighElem(:,:)             ! Elem: 1) tested, 2)
INTEGER                        :: NeighElemID
INTEGER,ALLOCATABLE            :: nNeighElems(:)
INTEGER                        :: offsetNeighElem
INTEGER                        :: currentOffset,currentCounter
! BGM
INTEGER                        :: iBGM,jBGM,kBGM
INTEGER                        :: iBGMmin,jBGMmin,kBGMmin
INTEGER                        :: iBGMmax,jBGMmax,kBGMmax
! Compute-node
INTEGER                        :: CNNeighElemID
! Processor
INTEGER                        :: iNeighElem
LOGICAL,ALLOCATABLE            :: ElemIsNeighElem(:)
! MPI
INTEGER                        :: iProc
INTEGER                        :: sendint,recvint
! Offsets for MPI_ALLGATHERV
INTEGER,ALLOCATABLE            :: displsNeighElem(:),recvcountNeighElem(:)
! Sorting
INTEGER,ALLOCATABLE            :: Neigh_Sort(:)
! Timers
! REAL                           :: StartT,EndT
!===================================================================================================================================

IF (.NOT.doCalcPartCollision) RETURN

! IF (MPIRoot) THEN
!   ! WRITE(UNIT_StdOut,'(132("-"))')
!   WRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors ...'
!   StartT = MPI_WTIME()
! END IF ! MPIRoot

#if USE_MPI
firstElem = INT(REAL(myComputeNodeRank  )*REAL(nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL(myComputeNodeRank+1)*REAL(nComputeNodeElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

! allocate arrays to store detected neighboring elements for the currently considered element only
! i.e. these lists are reset when we consider a new element
ALLOCATE(currentNeighElem(2,nComputeNodeElems))
ALLOCATE(nNeighElems(       nComputeNodeElems))
currentNeighElem = .FALSE.
nNeighElems      = 0

DO iElem = firstElem,lastElem
  ElemID  = GetGlobalElemID(iElem)

  OrigBoundsOfElemCenter(1:3) = (/    SUM(BoundsOfElem_Shared(1:2,1,iElem)),                                                   &
                                      SUM(BoundsOfElem_Shared(1:2,2,iElem)),                                                   &
                                      SUM(BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
  OrigBoundsOfElemCenter(4) = VECNORM ((/ BoundsOfElem_Shared(2,1,iElem)-BoundsOfElem_Shared(1,1,iElem),                       &
                                          BoundsOfElem_Shared(2,2,iElem)-BoundsOfElem_Shared(1,2,iElem),                       &
                                          BoundsOfElem_Shared(2,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)

  ! Check all periodic linear combinations
  DO iPeriodic = 0,GEO%nPeriodicVectors
    DO jPeriodic = iPeriodic,GEO%nPeriodicVectors
      DO kPeriodic = jPeriodic,GEO%nPeriodicVectors
        ! Avoid duplicates
        IF (iPeriodic.NE.0 .AND. iPeriodic.GE.kPeriodic) CYCLE
        IF (jPeriodic.NE.0 .AND. jPeriodic.GE.kPeriodic) CYCLE
        IF (iPeriodic.NE.0 .AND. iPeriodic.GE.jPeriodic) CYCLE

        DO iDir = -1, 1, 2; DO jDir = -1, 1, 2; DO kDir = -1, 1, 2
          ! Avoid duplicates
          IF (iPeriodic.EQ.0 .AND. iDir.EQ.-1) CYCLE
          IF (jPeriodic.EQ.0 .AND. jDir.EQ.-1) CYCLE
          IF (kPeriodic.EQ.0 .AND. kdir.EQ.-1) CYCLE ! No displacement

          BoundsOfElemCenter(4) = OrigBoundsOfElemCenter(4)
          ! One periodic vector
          IF (    kPeriodic.EQ.0) THEN
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3)
          ELSE IF(jPeriodic.EQ.0) THEN
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,kPeriodic)*REAL(kDir)
          ! Two periodic vectors
          ELSE IF (iPeriodic.EQ.0) THEN
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,kPeriodic)*REAL(kDir) &
                                                                  + GEO%PeriodicVectors(1:3,jPeriodic)*REAL(jDir)
          ! Three periodic vectors
          ELSE
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,kPeriodic)*REAL(kDir) &
                                                                  + GEO%PeriodicVectors(1:3,jPeriodic)*REAL(jDir) &
                                                                  + GEO%PeriodicVectors(1:3,iPeriodic)*REAL(iDir)
          END IF
          ! Hier gibt es keine Llamas

          ! compute the limits of the sphere centered on point BoundsOfElemCenter(1:3) and of radius BoundsOfElemCenter(4) + halo_eps
          iBGMmin = MAX(GEO%FIBGMimin,FLOOR  ((BoundsOfElemCenter(1) - (halo_eps+BoundsOfElemCenter(4)) - GEO%xminglob)/GEO%FIBGMdeltas(1)) - 1)
          iBGMmax = MIN(GEO%FIBGMimax,CEILING((BoundsOfElemCenter(1) + (halo_eps+BoundsOfElemCenter(4)) - GEO%xminglob)/GEO%FIBGMdeltas(1)) + 1)
          jBGMmin = MAX(GEO%FIBGMjmin,FLOOR  ((BoundsOfElemCenter(2) - (halo_eps+BoundsOfElemCenter(4)) - GEO%yminglob)/GEO%FIBGMdeltas(2)) - 1)
          jBGMmax = MIN(GEO%FIBGMjmax,CEILING((BoundsOfElemCenter(2) + (halo_eps+BoundsOfElemCenter(4)) - GEO%yminglob)/GEO%FIBGMdeltas(2)) + 1)
          kBGMmin = MAX(GEO%FIBGMkmin,FLOOR  ((BoundsOfElemCenter(3) - (halo_eps+BoundsOfElemCenter(4)) - GEO%zminglob)/GEO%FIBGMdeltas(3)) - 1)
          kBGMmax = MIN(GEO%FIBGMkmax,CEILING((BoundsOfElemCenter(3) + (halo_eps+BoundsOfElemCenter(4)) - GEO%zminglob)/GEO%FIBGMdeltas(3)) + 1)

          ! loop over cells of the BGM that are in the radius of ElemID
          DO kBGM=kBGMmin,kBGMmax ; DO jBGM=jBGMmin,jBGMmax ; DO iBGM=iBGMmin,iBGMmax
  ElemLoop: DO iLocElem = 1,FIBGM_nElems(iBGM,jBGM,kBGM)

              ! retrieve the (global ID) of the tested element
              NeighElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iLocElem)
              CNNeighElemID = GetCNElemID(NeighElemID)

              ! check if tested element NeighElemID was previously tested/found
              IF (ANY(currentNeighElem(:,CNNeighElemID),1)) CYCLE
              currentNeighElem(1,CNNeighElemID) = .TRUE.

              ! ignore myself
              ! IF (ElemID.EQ.NeighElemID) CYCLE

              LocalBoundsOfElemCenter(1:3) = (/ SUM(   BoundsOfElem_Shared(1:2,1,NeighElemID)),                                         &
                                                SUM(   BoundsOfElem_Shared(1:2,2,NeighElemID)),                                         &
                                                SUM(   BoundsOfElem_Shared(1:2,3,NeighElemID)) /) / 2.
              LocalBoundsOfElemCenter(4) = VECNORM ((/ BoundsOfElem_Shared(2  ,1,NeighElemID)-BoundsOfElem_Shared(1,1,NeighElemID),        &
                                                       BoundsOfElem_Shared(2  ,2,NeighElemID)-BoundsOfElem_Shared(1,2,NeighElemID),        &
                                                       BoundsOfElem_Shared(2  ,3,NeighElemID)-BoundsOfElem_Shared(1,3,NeighElemID) /) / 2.)

              ! compute the distance between the center of the bounding box of both elements
              dist = VECNORM(BoundsOfElemCenter(1:3) - LocalBoundsOfElemCenter(1:3))

              ! Also check the periodic combinations
              ! CALL IdentifyPeriodicNeighboringElements(BoundsOfElemCenter,LocalBoundsOfElemCenter,dist_periodic)
              ! dist = MIN(dist,dist_periodic)

              ! the current tested element is discared if it is not close enough to the current element
              IF (dist .GT. halo_eps + BoundsOfElemCenter(4) + LocalBoundsOfElemCenter(4)) CYCLE ElemLoop
              currentNeighElem(2,CNNeighElemID) = .TRUE.

            END DO ElemLoop

            ! Count the number of total neighbor elements
            nNeighElems(iElem) =  COUNT(currentNeighElem(2,:))
          END DO; END DO; END DO
          ! Reset the tested counter as periodic sides must be checked in all directions
          currentNeighElem(1,:) = .FALSE.
        END DO; END DO; END DO
      END DO ! kPeriodic
    END DO ! jPeriodic
  END DO ! iPeriodic

  currentNeighElem(2,:) = .FALSE.
END DO ! iElem = firstElem,lastElem

sendint = SUM(nNeighElems)
recvint = 0
CALL MPI_EXSCAN(sendint,recvint,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetNeighElem = recvint
! last proc knows CN total number of border elems
sendint = offsetNeighElem + SUM(nNeighElems)
CALL MPI_BCAST(sendint,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeNeighElems = sendint

! Arrays to hold the mapping
ALLOCATE(displsNeighElem(   0:nComputeNodeProcessors-1),&
         recvcountNeighElem(0:nComputeNodeProcessors-1))
displsNeighElem(myComputeNodeRank) = firstElem - 1
CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsNeighElem,1,MPI_INTEGER,MPI_COMM_SHARED,iError)
DO iProc=1,nComputeNodeProcessors-1
  recvcountNeighElem(iProc-1) = displsNeighElem(iProc)-displsNeighElem(iProc-1)
END DO
recvcountNeighElem(nComputeNodeProcessors-1) = nComputeNodeElems - displsNeighElem(nComputeNodeProcessors-1)
CALL MPI_ALLGATHERV(MPI_IN_PLACE,0                                 ,MPI_DATATYPE_NULL &
                   ,nNeighElems ,recvcountNeighElem,displsNeighElem,MPI_INTEGER      ,MPI_COMM_SHARED,iError)
DEALLOCATE(displsNeighElem   ,&
           recvcountNeighElem)

! allocated shared memory for nElems per BGM cell
! MPI shared memory is continuous, beginning from 1. All shared arrays have to
! be shifted to BGM[i]min with pointers
CALL Allocate_Shared((/nComputeNodeElems/),Neigh_nElems_Shared_Win,Neigh_nElems_Shared)
CALL MPI_WIN_LOCK_ALL(0,Neigh_nElems_Shared_Win,iError)
! allocated shared memory for BGM cell offset in 1D array of BGM to element mapping
CALL Allocate_Shared((/nComputeNodeElems/),Neigh_offsetElem_Shared_Win,Neigh_offsetElem_Shared)
CALL MPI_WIN_LOCK_ALL(0,Neigh_offsetElem_Shared_Win,iError)
Neigh_nElems     (1:nComputeNodeElems) => Neigh_nElems_Shared
Neigh_offsetElem (1:nComputeNodeElems) => Neigh_offsetElem_Shared

! last proc of compute-node writes into shared memory to make nElems per BGM accessible for every proc
IF (myComputeNodeRank.EQ.nComputeNodeProcessors-1) THEN
  currentOffset = 0
  DO iElem = 1,nComputeNodeElems
    ! senfbuf and recvbuf have to stay on original position. Shift 1 --> BGMimin
    Neigh_nElems(iElem)     = nNeighElems(iElem)
    Neigh_offsetElem(iElem) = currentOffset
    currentOffset           = currentoffset + nNeighElems(iElem)
  END DO ! iElem
END IF
DEALLOCATE(nNeighElems)

CALL BARRIER_AND_SYNC(Neigh_nElems_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(Neigh_offsetElem_Shared_Win,MPI_COMM_SHARED)

CALL Allocate_Shared((/nComputeNodeNeighElems/),CNElem2CNNeighElem_Win,CNElem2CNNeighElem)
CALL MPI_WIN_LOCK_ALL(0,CNElem2CNNeighElem_Win,iError)

currentNeighElem = .FALSE.
currentCounter   = 0.      !< All neighbor elements on current proc

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)

  OrigBoundsOfElemCenter(1:3) = (/    SUM(BoundsOfElem_Shared(1:2,1,iElem)),                                                   &
                                      SUM(BoundsOfElem_Shared(1:2,2,iElem)),                                                   &
                                      SUM(BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
  OrigBoundsOfElemCenter(4) = VECNORM ((/ BoundsOfElem_Shared(2,1,iElem)-BoundsOfElem_Shared(1,1,iElem),                       &
                                          BoundsOfElem_Shared(2,2,iElem)-BoundsOfElem_Shared(1,2,iElem),                       &
                                          BoundsOfElem_Shared(2,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)

  ! Check all periodic linear combinations
  DO iPeriodic = 0,GEO%nPeriodicVectors
    DO jPeriodic = iPeriodic,GEO%nPeriodicVectors
      DO kPeriodic = jPeriodic,GEO%nPeriodicVectors
        ! Avoid duplicates
        IF (iPeriodic.NE.0 .AND. iPeriodic.GE.kPeriodic) CYCLE
        IF (jPeriodic.NE.0 .AND. jPeriodic.GE.kPeriodic) CYCLE
        IF (iPeriodic.NE.0 .AND. iPeriodic.GE.jPeriodic) CYCLE

        DO iDir = -1, 1, 2; DO jDir = -1, 1, 2; DO kDir = -1, 1, 2
          ! Avoid duplicates
          IF (iPeriodic.EQ.0 .AND. iDir.EQ.-1) CYCLE
          IF (jPeriodic.EQ.0 .AND. jDir.EQ.-1) CYCLE
          IF (kPeriodic.EQ.0 .AND. kdir.EQ.-1) CYCLE ! No displacement

          BoundsOfElemCenter(4) = OrigBoundsOfElemCenter(4)
          ! One periodic vector
          IF (    kPeriodic.EQ.0) THEN
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3)
          ELSE IF(jPeriodic.EQ.0) THEN
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,kPeriodic)*REAL(kDir)
          ! Two periodic vectors
          ELSE IF (iPeriodic.EQ.0) THEN
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,kPeriodic)*REAL(kDir) &
                                                                  + GEO%PeriodicVectors(1:3,jPeriodic)*REAL(jDir)
          ! Three periodic vectors
          ELSE
            BoundsOfElemCenter(1:3) = OrigBoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,kPeriodic)*REAL(kDir) &
                                                                  + GEO%PeriodicVectors(1:3,jPeriodic)*REAL(jDir) &
                                                                  + GEO%PeriodicVectors(1:3,iPeriodic)*REAL(iDir)
          END IF
          ! Hier gibt es keine Llamas

          ! compute the limits of the sphere centered on point BoundsOfElemCenter(1:3) and of radius BoundsOfElemCenter(4) + halo_eps
          iBGMmin = MAX(GEO%FIBGMimin,FLOOR  ((BoundsOfElemCenter(1) - (halo_eps+BoundsOfElemCenter(4)) - GEO%xminglob)/GEO%FIBGMdeltas(1)) - 1)
          iBGMmax = MIN(GEO%FIBGMimax,CEILING((BoundsOfElemCenter(1) + (halo_eps+BoundsOfElemCenter(4)) - GEO%xminglob)/GEO%FIBGMdeltas(1)) + 1)
          jBGMmin = MAX(GEO%FIBGMjmin,FLOOR  ((BoundsOfElemCenter(2) - (halo_eps+BoundsOfElemCenter(4)) - GEO%yminglob)/GEO%FIBGMdeltas(2)) - 1)
          jBGMmax = MIN(GEO%FIBGMjmax,CEILING((BoundsOfElemCenter(2) + (halo_eps+BoundsOfElemCenter(4)) - GEO%yminglob)/GEO%FIBGMdeltas(2)) + 1)
          kBGMmin = MAX(GEO%FIBGMkmin,FLOOR  ((BoundsOfElemCenter(3) - (halo_eps+BoundsOfElemCenter(4)) - GEO%zminglob)/GEO%FIBGMdeltas(3)) - 1)
          kBGMmax = MIN(GEO%FIBGMkmax,CEILING((BoundsOfElemCenter(3) + (halo_eps+BoundsOfElemCenter(4)) - GEO%zminglob)/GEO%FIBGMdeltas(3)) + 1)

          ! loop over cells of the BGM that are in the radius of ElemID
          DO kBGM=kBGMmin,kBGMmax ; DO jBGM=jBGMmin,jBGMmax ; DO iBGM=iBGMmin,iBGMmax
 ElemLoop2: DO iLocElem = 1,FIBGM_nElems(iBGM,jBGM,kBGM)

              ! retrieve the (global ID) of the tested element
              NeighElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iLocElem)
              CNNeighElemID = GetCNElemID(NeighElemID)

              ! check if tested element NeighElemID was previously tested/found
              IF (ANY(currentNeighElem(:,CNNeighElemID),1)) CYCLE
              currentNeighElem(1,CNNeighElemID) = .TRUE.

              ! ignore myself
              ! IF (ElemID.EQ.NeighElemID) CYCLE

              LocalBoundsOfElemCenter(1:3) = (/ SUM(   BoundsOfElem_Shared(1:2,1,NeighElemID)),                                         &
                                                SUM(   BoundsOfElem_Shared(1:2,2,NeighElemID)),                                         &
                                                SUM(   BoundsOfElem_Shared(1:2,3,NeighElemID)) /) / 2.
              LocalBoundsOfElemCenter(4) = VECNORM ((/ BoundsOfElem_Shared(2  ,1,NeighElemID)-BoundsOfElem_Shared(1,1,NeighElemID),        &
                                                       BoundsOfElem_Shared(2  ,2,NeighElemID)-BoundsOfElem_Shared(1,2,NeighElemID),        &
                                                       BoundsOfElem_Shared(2  ,3,NeighElemID)-BoundsOfElem_Shared(1,3,NeighElemID) /) / 2.)

              ! compute the distance between the center of the bounding box of both elements
              dist = VECNORM(BoundsOfElemCenter(1:3) - LocalBoundsOfElemCenter(1:3))

              ! the current tested element is discared if it is not close enough to the current element
              IF (dist .GT. halo_eps + BoundsOfElemCenter(4) + LocalBoundsOfElemCenter(4)) CYCLE ElemLoop2
              currentNeighElem(2,CNNeighElemID) = .TRUE.
              currentCounter                    = currentCounter + 1

              ! Save the mapping
              CNElem2CNNeighElem(offsetNeighElem + currentCounter) = CNNeighElemID
            END DO ElemLoop2
          END DO; END DO; END DO
          ! Reset the tested counter as periodic sides must be checked in all directions
          currentNeighElem(1,:) = .FALSE.
        END DO; END DO; END DO
      END DO ! kPeriodic
    END DO ! jPeriodic
  END DO ! iPeriodic

  currentNeighElem(2,:) = .FALSE.
END DO ! iElem = firstElem,lastElem
CALL BARRIER_AND_SYNC(CNElem2CNNeighElem_Win,MPI_COMM_SHARED)
DEALLOCATE(currentNeighElem)

! Sort each neighbor element in an cell with ascending index
IF (myComputeNodeRank.EQ.0) THEN
  ! Temporary array
  ALLOCATE(Neigh_Sort(MAXVAL(FIBGM_nElems)))

  DO iElem = 1,nComputeNodeElems
    IF(Neigh_nElems(iElem).GT.1) THEN
      Neigh_Sort = CNElem2CNNeighElem(Neigh_offsetElem(iElem)+1:Neigh_offsetElem(iElem)+Neigh_nElems(iElem))
      CALL InsertionSort(Neigh_Sort,len=Neigh_nElems(iElem))
      CNElem2CNNeighElem(Neigh_offsetElem(iElem)+1:Neigh_offsetElem(iElem)+Neigh_nElems(iElem)) = Neigh_Sort
    END IF
  END DO ! iElem

  DEALLOCATE(Neigh_Sort)
END IF
CALL BARRIER_AND_SYNC(CNElem2CNNeighElem_Win,MPI_COMM_SHARED)

! Count number of own and neighbor elements on the current proc
ALLOCATE(ElemIsNeighElem(nComputeNodeTotalElems))
ElemIsNeighElem = .FALSE.

DO iElem = offsetElem+1,offsetElem+nElems
  ElemID = GetCNElemID(iElem)

  DO iNeighElem = Neigh_offsetElem(ElemID)+1,Neigh_offsetElem(ElemID)+Neigh_nElems(ElemID)
    ElemIsNeighElem(CNElem2CNNeighElem(iNeighElem)) = .TRUE.
  END DO ! iNeighElem
END DO ! iELem

nProcNeighElems = COUNT(ElemIsNeighElem)
ElemIsNeighElem = .FALSE.

currentCounter = 0.
ALLOCATE(NeighElemsProc(     nProcNeighElems))
ALLOCATE(offsetNeighElemPart(nGlobalElems))
NeighElemsProc      = -1
offsetNeighElemPart = -1

DO iElem = offsetElem+1,offsetElem+nElems
  ElemID = GetCNElemID(iElem)

  DO iNeighElem = Neigh_offsetElem(ElemID)+1,Neigh_offsetElem(ElemID)+Neigh_nElems(ElemID)
    IF (ElemIsNeighElem(CNElem2CNNeighElem(iNeighElem))) CYCLE

    ElemIsNeighElem(CNElem2CNNeighElem(iNeighElem)) = .TRUE.
    currentCounter = currentCounter + 1

    NeighElemsProc(currentCounter) = GetGlobalElemID(CNElem2CNNeighElem(iNeighElem))
  END DO ! iNeighElem
END DO ! iELem

! Sort each neighbor element in an cell with ascending index
CALL InsertionSort(NeighElemsProc,len=nProcNeighElems)

DEALLOCATE(ElemIsNeighElem)

END SUBROUTINE IdentifyNeighboringElements


SUBROUTINE UpdateParticleShared()
!===================================================================================================================================
!> Update particle arrays in shared array
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: displsElem,recvcountElem
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars         ,ONLY: myLeaderGroupRank,nLeaderGroupProcs
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SHARED
USE MOD_Particle_Output_Vars    ,ONLY: offsetnPart,locnPart
USE MOD_Particle_Tools          ,ONLY: GetOffsetAndGlobalNumberOfParts,UpdateNextFreePosition
USE MOD_Particle_Vars           ,ONLY: PDM,PEM,PartState,PartSpecies
USE MOD_Particle_Vars           ,ONLY: useLinkedList
USE MOD_Particle_Vars           ,ONLY: doCalcPartCollision
USE MOD_Particle_Collision_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: PartIntSize=2
INTEGER                        :: pcount
INTEGER                        :: iPart,iElem
INTEGER                        :: iProc
! Offsets for MPI_ALLGATHERV
INTEGER,ALLOCATABLE            :: MPI_COMM_LEADERS_REQUEST(:)           !> Request handle for non-blocking communication
! INTEGER,ALLOCATABLE            :: displsPartInt(:),recvcountPartInt(:)
!===================================================================================================================================

IF (.NOT.doCalcPartCollision) RETURN

! Determine number of particles in the complete domain
locnPart =   0
!>> Count number of particle on local proc
DO pcount = 1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(pcount)) THEN
    locnPart = locnPart + 1
  END IF
END DO

! Communicate the total number and offset
CALL GetOffsetAndGlobalNumberOfParts('UpdateParticleShared',offsetnPart,nGlobalNbrOfParticles,locnPart,.TRUE.)

! Allocate data arrays for mean particle quantities
IF (.NOT.ASSOCIATED(PartData_Shared)) THEN
  ! Allocate array if not yet associated
  ! > FIXME: This array should be changed to nComputeNodeTotalElems. Use GetGlobalElemID to map CN to global elemID and expose
  ! >        remote memory windows using one-sided RMA appraoch. Get the desired data using MPI_IGET and build a mapping
  CALL Allocate_Shared((/PP_nVarPart+1,INT(nGlobalNbrOfParticles(3)*1.2)/),PartData_Shared_Win,PartData_Shared)
  CALL Allocate_Shared((/              INT(nGlobalNbrOfParticles(3)*1.2)/),PartBC_Shared_Win  ,PartBC_Shared)
  CALL MPI_WIN_LOCK_ALL(0,PartData_Shared_Win,iError)
  CALL MPI_WIN_LOCK_ALL(0,PartBC_Shared_Win  ,iError)
ELSEIF (INT(SIZE(PartData_Shared)/PP_nVarPart).LT.nGlobalNbrOfParticles(3)) THEN
  ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
  CALL MPI_WIN_UNLOCK_ALL(PartData_Shared_Win,iError)
  CALL MPI_WIN_FREE(      PartData_Shared_Win,iError)
  CALL MPI_WIN_UNLOCK_ALL(PartBC_Shared_Win  ,iError)
  CALL MPI_WIN_FREE(      PartBC_Shared_Win  ,iError)
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

  ! Then, free the pointers or arrays
  MDEALLOCATE(PartData_Shared)
  MDEALLOCATE(PartBC_Shared)
  ! Increase array size if needed, 20% margin
  CALL Allocate_Shared((/PP_nVarPart+1,INT(nGlobalNbrOfParticles(3)*1.2)/),PartData_Shared_Win,PartData_Shared)
  CALL Allocate_Shared((/              INT(nGlobalNbrOfParticles(3)*1.2)/),PartBC_Shared_Win  ,PartBC_Shared)
  CALL MPI_WIN_LOCK_ALL(0,PartData_Shared_Win,iError)
  CALL MPI_WIN_LOCK_ALL(0,PartBC_Shared_Win  ,iError)
END IF

! Order the particles along the SFC using a linked list
ALLOCATE(PEM%pStart (offsetElem+1:offsetElem+nElems)    , &
         PEM%pNumber(offsetElem+1:offsetElem+nElems)    , &
         PEM%pNext  (1           :PDM%maxParticleNumber), &
         PEM%pEnd   (offsetElem+1:offsetElem+nElems))
useLinkedList = .TRUE.
CALL UpdateNextFreePosition()

! Walk along the linked list and fill the data arrays
iPart = offsetnPart
! Allocate inverse mapping
ALLOCATE(PEM2PartID(offsetnPart+1:offsetnPart+locnPart))
PEM2PartID = -1
! Walk over all elements on local proc
DO iElem = offsetElem+1,offsetElem+nElems
  ! Set start of particle numbers in current element
  PartInt_Shared(1,iElem) = iPart
  ! Find all particles in current element
  PartInt_Shared(2,iElem) = PartInt_Shared(1,iElem) + PEM%pNumber(iElem)
  ! Sum up particles and add properties to output array
  pcount = PEM%pStart(iElem)
  DO iPart = PartInt_Shared(1,iElem)+1,PartInt_Shared(2,iElem)
    PartData_Shared(1:PP_nVarPart,iPart) = PartState(:,pcount)
    PartData_Shared(PP_nVarPart+1,iPart) = PartSpecies(pcount)
    PEM2PartID(iPart)                    = pcount

    ! Set the index to the next particle
    pcount = PEM%pNext(pcount)
  END DO
  ! Set counter to the end of particle number in the current element
  iPart = PartInt_Shared(2,iElem)
  ! PartInt_Shared(2,iElem) = iPart
END DO ! iElem = offsetElem+1,offsetElem+nElems

! De-allocate linked list and return to normal particle array mode
useLinkedList=.FALSE.
! DEALLOCATE( PEM%pStart   &
!           , PEM%pNumber  &
!           , PEM%pNext    &
!           , PEM%pEnd)

CALL BARRIER_AND_SYNC(PartInt_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PartData_Shared_Win,MPI_COMM_SHARED)

IF (myComputeNodeRank.EQ.0) THEN
  ! Arrays for the compute node to hold the part offsets
  ALLOCATE(displsPartInt(   0:nLeaderGroupProcs-1),&
           recvcountPartInt(0:nLeaderGroupProcs-1))

  displsPartInt(myLeaderGroupRank) = offsetnPart
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsPartInt,1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
  DO iProc=1,nLeaderGroupProcs-1
    recvcountPartInt(iProc-1) = displsPartInt(iProc)-displsPartInt(iProc-1)
  END DO
  recvcountPartInt(nLeaderGroupProcs-1) = nGlobalNbrOfParticles(3) - displsPartInt(nLeaderGroupProcs-1)

  ALLOCATE(MPI_COMM_LEADERS_REQUEST(1:2))
  CALL MPI_IALLGATHERV(MPI_IN_PLACE   ,0                                                         ,MPI_DATATYPE_NULL    &
                      ,PartInt_Shared ,PartIntSize*recvcountElem   ,PartIntSize*displsElem       ,MPI_INTEGER          &
                      ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(1),iError)
  CALL MPI_IALLGATHERV(MPI_IN_PLACE   ,0                                                         ,MPI_DATATYPE_NULL    &
                      ,PartData_Shared,(PP_nVarPart+1)*recvcountPartInt,PP_nVarPart*displsPartInt,MPI_DOUBLE_PRECISION &
                      ,MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_REQUEST(2),iError)
  ! DEALLOCATE(displsPartInt)
  ! DEALLOCATE(recvcountPartInt)

  ! Non-blocking but actually blocking
  CALL MPI_WAITALL(1,MPI_COMM_LEADERS_REQUEST,MPI_STATUSES_IGNORE,iError)
  DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
END IF ! myComputeNodeRank.EQ.0

CALL BARRIER_AND_SYNC(PartInt_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PartData_Shared_Win,MPI_COMM_SHARED)

! IPWRITE(*,*) 'PartInt_Shared(1,:):', PartInt_Shared(1,:)
! IPWRITE(*,*) 'PartInt_Shared(2,:):', PartInt_Shared(2,:)

END SUBROUTINE UpdateParticleShared


! SUBROUTINE IdentifyPeriodicNeighboringElements(BoundsOfElemCenter,LocalBoundsOfElemCenter,dist)
! !===================================================================================================================================
! !> Calculate the minimum distance between two elements considering periodic displacements
! !===================================================================================================================================
! ! MODULES
! USE MOD_Particle_Globals        ,ONLY: VECNORM
! USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT/OUTPUT VARIABLES
! REAL,INTENT(IN)              :: BoundsOfElemCenter(4)
! REAL,INTENT(IN)              :: LocalBoundsOfElemCenter(4)
! REAL,INTENT(OUT)             :: dist
! !----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER                      :: iDir,jDir,kDir
! INTEGER                      :: iPeriodicVector,jPeriodicVector
! !===================================================================================================================================
!
! ! Initialize HUGE distance
! dist = HUGE(1.)
!
! SELECT CASE(GEO%nPeriodicVectors)
!   CASE(1)
!     ! check two directions
!     DO iDir = -1, 1, 2
!       dist = MIN(dist,VECNORM(  BoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,1)*REAL(iDir) - LocalBoundsOfElemCenter(1:3)))
!     END DO
!
!   CASE(2)
!     ! check the two possible periodic vectors. Begin with checking the single periodic vector, followed by the combination of
!     ! the first periodic vector with the other, 1,2,1+2
!     DO iPeriodicVector = 1,2
!       DO iDir = -1, 1, 2
!         dist = MIN(dist,VECNORM(BoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,iPeriodicVector)*REAL(iDir) - LocalBoundsOfElemCenter(1:3)))
!       END DO ! iDir = -1, 1, 2
!
!       ! Check linear combination of two periodic vectors
!       DO jPeriodicVector = 1,2
!         IF (iPeriodicVector.GE.jPeriodicVector) CYCLE
!
!         DO iDir = -1, 1, 2
!           DO jDir = -1, 1, 2
!             dist = MIN(dist,VECNORM( BoundsOfElemCenter(1:3)                                                 &
!                       + GEO%PeriodicVectors(1:3,iPeriodicVector)*REAL(iDir)                                  &
!                       + GEO%PeriodicVectors(1:3,jPeriodicVector)*REAL(jDir) - LocalBoundsOfElemCenter(1:3) ))
!           END DO ! jDir = -1, 1, 2
!         END DO ! iDir = -1, 1, 2
!       END DO ! jPeriodicVector = 1,2
!     END DO ! iPeriodicVector = 1,2
!
!   CASE(3)
!     ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
!     ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
!     DO iPeriodicVector = 1,3
!       DO iDir = -1, 1, 2
!         dist = MIN(dist,VECNORM(BoundsOfElemCenter(1:3) + GEO%PeriodicVectors(1:3,iPeriodicVector)*REAL(iDir) - LocalBoundsOfElemCenter(1:3)))
!       END DO ! iDir = -1, 1, 2
!
!       ! Combination of two periodic vectors
!       DO jPeriodicVector = 1,3
!         IF (iPeriodicVector.GE.jPeriodicVector) CYCLE
!
!         DO iDir = -1, 1, 2
!           DO jDir = -1, 1, 2
!             dist = MIN(dist,VECNORM( BoundsOfElemCenter(1:3)                                                 &
!                       + GEO%PeriodicVectors(1:3,iPeriodicVector)*REAL(iDir)                                  &
!                       + GEO%PeriodicVectors(1:3,jPeriodicVector)*REAL(jDir) - LocalBoundsOfElemCenter(1:3) ))
!           END DO ! jDir = -1, 1, 2
!         END DO ! iDir = -1, 1, 2
!       END DO ! jPeriodicVector = 1,3
!     END DO ! iPeriodicVector = 1,3
!
!     ! Combination of three periodic vectors
!     DO iDir = -1, 1, 2
!       DO jDir = -1, 1, 2
!         DO kDir = -1, 1, 2
!           dist = MIN(dist,VECNORM( BoundsOfElemCenter(1:3)                                                   &
!                       + GEO%PeriodicVectors(1:3,1)*REAL(iDir)                                                &
!                       + GEO%PeriodicVectors(1:3,2)*REAL(jDir)                                                &
!                       + GEO%PeriodicVectors(1:3,3)*REAL(kDir) - LocalBoundsOfElemCenter(1:3)))
!         END DO ! kDir = -1, 1, 2
!       END DO ! jDir = -1, 1, 2
!     END DO ! iDir = -1, 1, 2
!
! END SELECT
!
! END SUBROUTINE IdentifyPeriodicNeighboringElements


SUBROUTINE FinalizeCollision()
!===================================================================================================================================
!> Finalize particle collisions module
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Vars           ,ONLY: doCalcPartCollision
#if USE_MPI
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF (.NOT.doCalcPartCollision) RETURN

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

! IdentifyNeighborElements
CALL MPI_WIN_UNLOCK_ALL(Neigh_nElems_Shared_Win        ,iError)
CALL MPI_WIN_FREE(      Neigh_nElems_Shared_Win        ,iError)
CALL MPI_WIN_UNLOCK_ALL(Neigh_offsetElem_Shared_Win    ,iError)
CALL MPI_WIN_FREE(      Neigh_offsetElem_Shared_Win    ,iError)
CALL MPI_WIN_UNLOCK_ALL(CNElem2CNNeighElem_Win         ,iError)
CALL MPI_WIN_FREE(      CNElem2CNNeighElem_Win         ,iError)
! ComputeParticleCollisions
CALL MPI_WIN_UNLOCK_ALL(PartInt_Shared_Win             ,iError)
CALL MPI_WIN_FREE(      PartInt_Shared_Win             ,iError)
CALL MPI_WIN_UNLOCK_ALL(PartData_Shared_Win            ,iError)
CALL MPI_WIN_FREE(      PartData_Shared_Win            ,iError)
CALL MPI_WIN_UNLOCK_ALL(PartBC_Shared_Win              ,iError)
CALL MPI_WIN_FREE(      PartBC_Shared_Win              ,iError)

CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! Then, free the pointers or arrays
! IdentifyNeighborElements
MDEALLOCATE(Neigh_nElems_Shared)
MDEALLOCATE(Neigh_offsetElem_Shared)
MDEALLOCATE(CNElem2CNNeighElem)
! ComputeParticleCollisions
MDEALLOCATE(PartInt_Shared)
MDEALLOCATE(PartData_Shared)
MDEALLOCATE(PartBC_Shared)

SDEALLOCATE(NeighElemsProc)
SDEALLOCATE(offsetNeighElemPart)

END SUBROUTINE FinalizeCollision

#endif /*PARTICLES_COUPLING == 4*/

END MODULE MOD_Particle_Collision
