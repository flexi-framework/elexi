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

! INTERFACE FinalizeCollision
!   MODULE PROCEDURE FinalizeCollision
! END INTERFACE

PUBLIC :: DefineParametersCollission
PUBLIC :: InitializeCollision
! PUBLIC :: FinalizeCollision
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
USE MOD_ReadInTools                ,ONLY: GETLOGICAL,GETREAL
USE MOD_Particle_Collision_Vars
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: StartT,EndT
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE COLLISSIONS...'
GETTIME(StartT)

! initialize the parameters of the model particle collision
PartCollisionModel%e        = GETREAL   ('Part-Collisions-RestCoeff')
IF (PartCollisionModel%e.LT.0. .OR. PartCollisionModel%e.GT.1.) CALL Abort(__STAMP__,'Part-Collisions-RestCoeff is not between 0 and 1!')
PartCollisionModel%Friction = GETLOGICAL('Part-Collisions-Friction')
PartCollisionModel%f        = GETREAL   ('Part-Collisions-FricCoeff')

! determine neighboring elements
CALL IdentifyNeighboringElements()

GETTIME(EndT)
LBWRITE(UNIT_stdOut,'(A,F0.3,A)') ' INIT PARTICLE COLLISIONS DONE! [',EndT-StartT,'s]'
LBWRITE(UNIT_stdOut,'(132("-"))')

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
USE MOD_MPI_Shared              ,ONLY: Allocate_Shared,BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
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
LOGICAL,ALLOCATABLE            :: currentNeighElem(:)
INTEGER                        :: NeighElemID
INTEGER,ALLOCATABLE            :: nNeighElems(:)
INTEGER                        :: offsetNeighElem
INTEGER                        :: currentOffset
! BGM
INTEGER                        :: iBGM,jBGM,kBGM
INTEGER                        :: iBGMmin,jBGMmin,kBGMmin
INTEGER                        :: iBGMmax,jBGMmax,kBGMmax
! Compute-node
INTEGER                        :: CNNeighElemID
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
ALLOCATE(currentNeighElem(nComputeNodeElems))
ALLOCATE(nNeighElems(     nComputeNodeElems))
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
              NeighElemID = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iLocElem)

              ! ignore myself
              IF (ElemID.EQ.NeighElemID) CYCLE

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

              CNNeighElemID = GetCNElemID(NeighElemID)

              ! check if tested element NeighElemID is not already a detected neighbor of ElemID (to avoid duplicates)
              IF (currentNeighElem(CNNeighElemID)) CYCLE
              currentNeighElem(CNNeighElemID) = .TRUE.
            END DO ElemLoop

            ! Count the number of total neighbor elements
            nNeighElems(iElem) =  COUNT(currentNeighElem)
          END DO; END DO; END DO
        END DO; END DO; END DO
      END DO ! kPeriodic
    END DO ! jPeriodic
  END DO ! iPeriodic

  currentNeighElem   = .FALSE.
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
displsNeighElem(myComputeNodeRank) = offsetNeighElem
CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,displsNeighElem,1,MPI_INTEGER,MPI_COMM_SHARED,iError)
DO iProc=1,nComputeNodeProcessors-1
  recvcountNeighElem(iProc-1) = displsNeighElem(iProc)-displsNeighElem(iProc-1)
END DO
recvcountNeighElem(nComputeNodeProcessors-1) = nComputeNodeNeighElems - displsNeighElem(nComputeNodeProcessors-1)
CALL MPI_ALLGATHERV(MPI_IN_PLACE,0                                 ,MPI_DATATYPE_NULL &
                   ,nNeighElems ,recvcountNeighElem,displsNeighElem,MPI_INTEGER      ,MPI_COMM_SHARED,iError)

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
CALL BARRIER_AND_SYNC(Neigh_nElems_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(Neigh_offsetElem_Shared_Win,MPI_COMM_SHARED)

CALL Allocate_Shared((/nComputeNodeNeighElems/),CNElem2CNNeighElem_Win,CNElem2CNNeighElem)
CALL MPI_WIN_LOCK_ALL(0,CNElem2CNNeighElem_Win,IERROR)

currentNeighElem = .FALSE.

! reset the total number of detected neighboring elements
nNeighElems = 0

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
        ElemLoop2: DO iLocElem = 1,FIBGM_nElems(iBGM,jBGM,kBGM)

              ! retrieve the (global ID) of the tested element
              NeighElemID = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iLocElem)

              ! ignore myself
              IF (ElemID.EQ.NeighElemID) CYCLE

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

              CNNeighElemID = GetCNElemID(NeighElemID)

              ! check if tested element NeighElemID is not already a detected neighbor of ElemID (to avoid duplicates)
              IF (currentNeighElem(CNNeighElemID)) CYCLE
              currentNeighElem(CNNeighElemID) = .TRUE.
              nNeighElems                     = nNeighElems + 1

              ! Save the mapping
              CNElem2CNNeighElem(offsetNeighElem + nNeighElems) = CNNeighElemID
            END DO ElemLoop2
          END DO; END DO; END DO
        END DO; END DO; END DO
      END DO ! kPeriodic
    END DO ! jPeriodic
  END DO ! iPeriodic

  currentNeighElem = .FALSE.
END DO ! iElem = firstElem,lastElem

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
END IF

END SUBROUTINE IdentifyNeighboringElements


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

#endif /*PARTICLES_COUPLING == 4*/

END MODULE MOD_Particle_Collision
