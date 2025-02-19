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
MODULE MOD_Particle_Collision_Method
! MODULES
IMPLICIT NONE
PRIVATE

#if PARTICLES_COUPLING == 4
PUBLIC:: ComputeParticleCollisions
!===================================================================================================================================

CONTAINS

SUBROUTINE ComputeParticleCollisions(dtLoc)
!===================================================================================================================================
!> Compute particle collisions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Readin              ,ONLY: ELEMIPROC
USE MOD_Mesh_Vars                ,ONLY: nElems,offsetElem
USE MOD_Mesh_Vars                ,ONLY: nComputeNodeElems
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Localization    ,ONLY: SinglePointToElement
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemToBCSides
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNElemID,GetGlobalElemID
USE MOD_Particle_Vars            ,ONLY: PEM
USE MOD_Particle_Vars            ,ONLY: doCalcPartCollision
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars,          ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars          ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars,          ONLY: MPI_COMM_SHARED!,MPI_COMM_LEADERS_SHARED
USE MOD_MPI_Shared_Vars          ,ONLY: myComputeNodeRank,nComputeNodeProcessors
! USE MOD_MPI_Shared_Vars          ,ONLY: nLeaderGroupProcs
USE MOD_Particle_Vars            ,ONLY: offsetPartMPI
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBPauseTime
USE MOD_LoadBalance_Vars         ,ONLY: nCollsPerElem
USE MOD_Particle_Globals         ,ONLY: ElementOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                :: dtLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: PartBCdt
INTEGER                        :: iElem,ElemID
INTEGER                        :: firstCNElem,lastCNElem
INTEGER                        :: iPart,iPart1,iPart2
INTEGER                        :: CNElemID,CNNeighElemID
INTEGER                        :: NeighElemID,iCNNeighElem
! Trajectories and position
REAL                           :: P1(3),P2(3)
REAL                           :: pColl1(3),pColl2(3)
REAL                           :: V1(3),V2(3)
! Collisions
REAL                           :: a,b,c,delta,dtColl
! Neighbor particles
LOGICAL,ALLOCATABLE            :: PartColl(:)
! Communication
INTEGER                        :: ElemProc
INTEGER                        :: CNRank,CNRootRank
! INTEGER                        :: MPI_WINDOW(   0:nLeaderGroupProcs)
! Timers
#if USE_LOADBALANCE
REAL                           :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

IF (.NOT.doCalcPartCollision) RETURN

! safety check, dt > 1 currently not supported
IF (dtLoc.GT.1) CALL Abort(__STAMP__,'Time steps greater than 1 are unsupported with particle collisions')

! nullify the BC intersection times
IF (myComputeNodeRank.EQ.0) PartBC_Shared = -1.
CALL BARRIER_AND_SYNC(PartBC_Shared_Win,MPI_COMM_SHARED)

! nullify the collision counter
IF (myComputeNodeRank.EQ.0) PartColl_Shared = 0.
CALL BARRIER_AND_SYNC(PartColl_Shared_Win,MPI_COMM_SHARED)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! loop over internal elements
DO iElem = offsetElem+1,offsetElem+nElems
  ! skip if there is no particles inside the current element
  CNElemID = GetCNElemID(iElem)
  IF (PartInt_Shared(2,CNElemID).EQ.PartInt_Shared(1,CNElemID)) CYCLE

  SELECT CASE(TrackingMethod)
    CASE (TRIATRACKING)
      ! FIXME
      CALL Abort(__STAMP__,'TODO: TriaTracking')
    CASE (TRACING)
      ! loop over all particles in the element
      DO iPart = PartInt_Shared(1,CNElemID)+1,PartInt_Shared(2,CNElemID)
        CALL ComputeParticleTracingIntersection(   PartID    = iPart                , &
                                                   ElemID    = iELem                , &
                                                   CNElemID  = CNElemID             , &
                                                   dtLoc     = dtLoc                , &
                                                   dtBC      = PartBCdt)

        ! Scale intersection dt with time step
        PartBC_Shared(iPart) = PartBCdt * dtLoc
      END DO ! iPart
    CASE (REFMAPPING)
      ! loop over all particles in the element
      DO iPart = PartInt_Shared(1,CNElemID)+1,PartInt_Shared(2,CNElemID)
        ASSOCIATE(offsetSide => ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID)           , &
                  nSides     => ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID))
        CALL ComputeParticleRefmappingIntersection(PartID    = iPart                , &
                                                   CNElemID  = CNElemID             , &
                                                   firstSide = offsetSide + 1       , &
                                                   lastSide  = offsetSide + nSides  , &
                                                   nLocSides = nSides               , &
                                                   dtLoc     = dtLoc                , &
                                                   dtBC      = PartBCdt)
        END ASSOCIATE

        ! Scale intersection dt with time step
        PartBC_Shared(iPart) = PartBCdt * dtLoc
      END DO ! iPart
    END SELECT
END DO ! iElem
CALL BARRIER_AND_SYNC(PartBC_Shared_Win,MPI_COMM_SHARED)

#if USE_LOADBALANCE
CALL LBPauseTime(LB_COLLISION,tLBStart)
#endif /*USE_LOADBALANCE*/

IF (myComputeNodeRank.EQ.0) THEN
  ! Communicate the PartBC_shared
  !> Specify a window of existing memory that is exposed to RMA accesses
  !> A process may elect to expose no memory by specifying size = 0
  ! CALL MPI_WIN_CREATE( PartBC_Shared                                                     &
  !                    , INT(SIZE_REAL*nComputeNodeParts,MPI_ADDRESS_KIND)                 & ! Only local particles are to be sent
  !                    , SIZE_REAL                                                         &
  !                    , MPI_INFO_NULL                                                     &
  !                    , MPI_COMM_LEADERS_SHARED                                           &
  !                    , MPI_WINDOW                                                        &
  !                    , iError)

  !> Start an RMA exposure epoch
  ! MPI_WIN_POST must complete first as per https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node281.htm
  ! > "The call to MPI_WIN_START can block until the matching call to MPI_WIN_POST occurs at all target processes."
  ! No local operations prior to this epoch, so give an assertion
  !> > ACTIVE SYNCHRONIZATION
  ! CALL MPI_WIN_FENCE(    MPI_MODE_NOPRECEDE                                                &
  !                   ,    PartBC_Window                                                     &
  !                   ,    iError)
  !> > PASSIVE SYNCHRONIZATION
  CALL MPI_WIN_LOCK_ALL(0,PartBC_Win,iError)

  ! Loop over all remote elements and fill the PartData_Shared array
  DO iElem = nComputeNodeElems+1,nComputeNodeTotalElems
    ElemID     = GetGlobalElemID(iElem)
    ElemProc   = ELEMIPROC(ElemID)
    CNRank     = INT(ElemProc/nComputeNodeProcessors)
    CNRootRank = CNRank*nComputeNodeProcessors

    ! Position in the RMA window is the ElemID minus the compute node offset
    ASSOCIATE(firstPart       => PartInt_Shared(1,iElem)+1                                 &
             ,lastPart        => PartInt_Shared(2,iElem)                                   &
             ,offsetFirstPart => PartInt_Shared(3,iElem)  -offsetPartMPI(CNRank)           &
             ,offsetLastPart  => PartInt_Shared(4,iElem)  -offsetPartMPI(CNRank)           &
             ,nPart           => PartInt_Shared(4,iElem)  -PartInt_Shared(3,iElem))

    IF (nPart.EQ.0) CYCLE

    CALL MPI_GET(          PartBC_Shared(firstPart)                                        &
                         , nPart                                                           &
                         , MPI_DOUBLE_PRECISION                                            &
                         , CNRank                                                          &
                         , INT(offsetFirstPart,MPI_ADDRESS_KIND)                           &
                         , nPart                                                           &
                         , MPI_DOUBLE_PRECISION                                            &
                         , PartBC_Win                                                      &
                         , iError)
    END ASSOCIATE
  END DO ! iElem

  !> Complete the epoch - this will block until MPI_Get is complete
  !> > ACTIVE SYNCHRONIZATION
  ! CALL MPI_WIN_FENCE(    0                                                                 &
  !                   ,    PartBC_Win                                                     &
  !                   ,    iError)
  ! ! All done with the window - tell MPI there are no more epochs
  ! CALL MPI_WIN_FENCE(    MPI_MODE_NOSUCCEED                                                &
  !                   ,    PartBC_Win                                                     &
  !                   ,    iError)
  !> > PASSIVE SYNCHRONIZATION
  CALL MPI_WIN_UNLOCK_ALL(PartBC_Win,iError)
  !
  ! Free up our window
  ! CALL MPI_WIN_FREE(     MPI_WINDOW                                                        &
  !                   ,    iError)
END IF ! myComputeNodeRank.EQ.0
CALL BARRIER_AND_SYNC(PartBC_Shared_Win,MPI_COMM_SHARED)

! Build a mapping from global particle IDs to neighbor particle IDs
! > Count number of particles in own and neighbor elements
nProcParts = 0
! FIXME: SOME COUNTERS HERE ARE WRONG!
DO iElem = 1,nProcNeighElems
  ElemID   = NeighElemsProc(iElem)
  CNElemID = GetCNElemID(ElemID)

  ! > Map global offset to neighbor offset
  offsetNeighElemPart(ElemID) = PartInt_Shared(1,CNElemID) - nProcParts

  !> Increment the counter by the current number of particles
  nProcParts = nProcParts + PartInt_Shared(2,CNElemID)-PartInt_Shared(1,CNElemID)
END DO

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! > Allocate an array to hold the mapping
ALLOCATE(PartColl(nProcParts))
PartColl = .FALSE.

firstCNElem = GetCNElemID(offsetElem+1)
 lastCNElem = GetCNElemID(offsetElem+nElems)

! > Allocate local collision arrays
! ALLOCATE(PartColl(PartInt_Shared(1,offsetElem+1)+1:PartInt_Shared(2,offsetEle+nElems)))

! loop over internal elements
DO iElem = offsetElem+1,offsetElem+nElems
  CNElemID = GetCNElemID(iElem)

  ! skip if there is no particles inside the current element
  IF (PartInt_Shared(2,CNElemID).EQ.PartInt_Shared(1,CNElemID)) CYCLE

  ! loop over all particles in the element
  DO iPart1 = PartInt_Shared(1,CNElemID)+1,PartInt_Shared(2,CNElemID)
    ! Ignore already collided particles
    IF (PartColl(iPart1 - offsetNeighElemPart(iElem))) CYCLE

    ! loop over neighbor elements (includes own element)
PartLoop: DO iCNNeighElem = Neigh_offsetElem(CNElemID)+1,Neigh_offsetElem(CNElemID)+Neigh_nElems(CNElemID)
      CNNeighElemID = CNElem2CNNeighElem(iCNNeighElem)
      NeighElemID   = GetGlobalElemID(CNNeighElemID)

      ! loop over all particles in the element
      DO iPart2 = PartInt_Shared(1,CNNeighElemID)+1,PartInt_Shared(2,CNNeighElemID)
        ! Ignore myself
        IF (iPart1.EQ.iPart2) CYCLE
        ! Ignore already collided particles
        IF (PartColl(iPart2 - offsetNeighElemPart(NeighElemID))) CYCLE

        ! Reconstruct the particle velocity
        V1 = (PartData_Shared(PART_POSV,iPart1) - PartData_Shared(PART_OLDV,iPart1)) / dtLoc
        V2 = (PartData_Shared(PART_POSV,iPart2) - PartData_Shared(PART_OLDV,iPart2)) / dtLoc

        ! FIXME: What about periodic?

        ! Check if particles will collide
        ! Solve: |x-y|**2 = (r1+r2)**2, x = xp1^n-xp2^n, y = dt_coll*up1^{n+1} - dt_coll*up2^{n+1}; dt_coll \in [0;1]
        ! Triangle inequality: |x+y|**2 <= (|x| + |y|)**2
        ! compute coefficients a,b,c of the previous quadratic equation (a*t^2+b*t+c=0)
        ASSOCIATE(Term => V2 - V1)
        a = DOT_PRODUCT(Term,Term)
        END ASSOCIATE

        ! || particle trajectories: if a.EQ.0 then b.EQ.0 and so the previous equation has no solution: go to the next pair
        IF (a.EQ.0.) CYCLE

        ! Cauchy-Schwarz inequality: <x,y> + <y,x> <= 2 |<x,y>| <= 2|x||y|
        P1 = PartData_Shared(PART_OLDV,iPart1)
        P2 = PartData_Shared(PART_OLDV,iPart2)
        b = 2. * DOT_PRODUCT(ABS(V1 - V2), ABS(P1-P2))
        ! b = 2. * DOT_PRODUCT(P1-P2,PartData_Shared(PART_POSV,iPart1) - PartData_Shared(PART_POSV,iPart2)-(P1 - P2))
        ASSOCIATE(Term => P1 - P2)
        c = DOT_PRODUCT(Term,Term) - 0.25*(PartData_Shared(PART_DIAM,iPart2) + PartData_Shared(PART_DIAM,iPart1))**2.
        END ASSOCIATE

        ! if delta < 0, no collision is possible for this pair
        delta = b * b - 4. * a * c
        IF (delta.LT.0.) CYCLE

        ! a collision is possible, determine when relatively to (tStage+dtloc)
        dtColl = MAX(MIN((-b + SQRT(delta)) / (2. * a),1.),MIN((-b - SQRT(delta)) / (2. * a),1.)) * dtLoc
        ! the collision is only valid if it occurred between tStage and tStage+dtLoc
        IF (dtColl.LT.0. .OR. dtColl.GT.dtLoc) CYCLE

        ! the collision is only valid if it occurred before a BC intersection
        ! > here, the particle radius needs to be considered since the collision is on the particle surface but the BC intersection
        ! > is determined based on the particle center of mass
        IF (dtColl.GT.PartBC_Shared(iPart1)) CYCLE
        IF (dtColl.GT.PartBC_Shared(iPart2)) CYCLE

        !IF (VECNORM(P1 + dtColl * V1 - (P2 + dtColl * V2)) &
          !.GT. 0.5*(PartData_Shared(PART_DIAM,iPart2) + PartData_Shared(PART_DIAM,iPart1))) CYCLE

        ! collision is only valid the particles approach each other
        ! FIXME: NOT TRUE, we can collide on with acute angles
        ASSOCIATE(N1 => V1 &
                 ,N2 => V2)
        pColl1 = PartData_Shared(PART_POSV,iPart1) - (dtLoc-dtColl) * V1
        pColl2 = PartData_Shared(PART_POSV,iPart2) - (dtLoc-dtColl) * V2
        IF (DOT_PRODUCT(pColl1 - P1,N1).LE.0) CYCLE
        IF (DOT_PRODUCT(pColl2 - P2,N2).LE.0) CYCLE
        END ASSOCIATE

        ! collision found
        PartColl(iPart1 - offsetNeighElemPart(iElem      )) = .TRUE.
        PartColl(iPart2 - offsetNeighElemPart(NeighElemID)) = .TRUE.

        PartColl_Shared(iPart1) = iPart2
        IF (iPart2.GE.PartInt_Shared(1,firstCNElem)+1 .AND. &
            iPart2.LE.PartInt_Shared(2, lastCNElem)) &
          PartColl_Shared(iPart2) = iPart1

        EXIT PartLoop
      END DO ! iPart2
    END DO PartLoop ! iCNNeighElem

#if USE_LOADBALANCE
    ! Cell is on current proc, assign load to this cell
    IF (ElementOnProc(iElem)) nCollsPerElem(iElem-offsetElem) = nCollsPerElem(iElem-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
  END DO ! iPart1

#if USE_LOADBALANCE
  CALL LBPauseTime(LB_COLLISION,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO ! iElem
SDEALLOCATE(PartColl)

CALL BARRIER_AND_SYNC(PartColl_Shared_Win,MPI_COMM_SHARED)
IF (myComputeNodeRank.EQ.0) THEN
  !> Start an RMA exposure epoch
  ! MPI_WIN_POST must complete first as per https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node281.htm
  ! > "The call to MPI_WIN_START can block until the matching call to MPI_WIN_POST occurs at all target processes."
  ! No local operations prior to this epoch, so give an assertion
  !> > ACTIVE SYNCHRONIZATION
  ! CALL MPI_WIN_FENCE(    MPI_MODE_NOPRECEDE                                                &
  !                   ,    PartColl_Win                                                      &
  !                   ,    iError)
  !> > PASSIVE SYNCHRONIZATION
  CALL MPI_WIN_LOCK_ALL(0,PartColl_Win,iError)

  ! Loop over all remote elements and fill the PartData_Shared array
  DO iElem = nComputeNodeElems+1,nComputeNodeTotalElems
    ElemID     = GetGlobalElemID(iElem)
    ElemProc   = ELEMIPROC(ElemID)
    CNRank     = INT(ElemProc/nComputeNodeProcessors)
    CNRootRank = CNRank*nComputeNodeProcessors

    ! Position in the RMA window is the ElemID minus the compute node offset
    ASSOCIATE(firstPart       => PartInt_Shared(1,iElem)+1                                 &
             ,lastPart        => PartInt_Shared(2,iElem)                                   &
             ,offsetFirstPart => PartInt_Shared(3,iElem)  -offsetPartMPI(CNRank)           &
             ,offsetLastPart  => PartInt_Shared(4,iElem)  -offsetPartMPI(CNRank)           &
             ,nPart           => PartInt_Shared(4,iElem)  -PartInt_Shared(3,iElem))

    IF (nPart.EQ.0) CYCLE

    CALL MPI_GET(          PartColl_Shared(firstPart)                                      &
                         , nPart                                                           &
                         , MPI_INTEGER                                                     &
                         , CNRank                                                          &
                         , INT(offsetFirstPart,MPI_ADDRESS_KIND)                           &
                         , nPart                                                           &
                         , MPI_INTEGER                                                     &
                         , PartColl_Win                                                    &
                         , iError)
    END ASSOCIATE
  END DO ! iElem

  !> Complete the epoch - this will block until MPI_Get is complete
  !> > ACTIVE SYNCHRONIZATION
  ! CALL MPI_WIN_FENCE(    0                                                                 &
  !                   ,    PartColl_Win                                                      &
  !                   ,    iError)
  ! ! All done with the window - tell MPI there are no more epochs
  ! CALL MPI_WIN_FENCE(    MPI_MODE_NOSUCCEED                                                &
  !                   ,    PartColl_Win                                                      &
  !                   ,    iError)
  !> > PASSIVE SYNCHRONIZATION
  CALL MPI_WIN_UNLOCK_ALL(PartColl_Win,iError)
END IF ! myComputeNodeRank.EQ.0
CALL BARRIER_AND_SYNC(PartColl_Shared_Win,MPI_COMM_SHARED)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! loop over internal elements
DO iElem = offsetElem+1,offsetElem+nElems
  CNElemID = GetCNElemID(iElem)

  ! skip if there is no particles inside the current element
  IF (PartInt_Shared(2,CNElemID).EQ.PartInt_Shared(1,CNElemID)) CYCLE

  ! loop over all particles in the element
  DO iPart1 = PartInt_Shared(1,CNElemID)+1,PartInt_Shared(2,CNElemID)
    ! Ignore already collided particles
    IF (PartColl_Shared(iPart1).EQ.0) CYCLE

    ! get particle index of collision partner
    iPart2 = PartColl_Shared(iPart1)
    ! Ignore myself
    ! IF (iPart1.EQ.iPart2) CYCLE
    ! Ignore already collided particles
    IF (iPart2.EQ.0) CYCLE
    IF (PartColl_Shared(iPart2).EQ.0) CYCLE
    IF (PartColl_Shared(iPart2).NE.iPart1) CYCLE

    ! Reconstruct the particle velocity
    V1 = (PartData_Shared(PART_POSV,iPart1) - PartData_Shared(PART_OLDV,iPart1)) / dtLoc
    V2 = (PartData_Shared(PART_POSV,iPart2) - PartData_Shared(PART_OLDV,iPart2)) / dtLoc

    ! FIXME: What about periodic?

    ! Check if particles will collide
    ! Solve: |x-y|**2 = (r1+r2)**2, x = xp1^n-xp2^n, y = dt_coll*up1^{n+1} - dt_coll*up2^{n+1}; dt_coll \in [0;1]
    ! Triangle inequality: |x+y|**2 <= (|x| + |y|)**2
    ! compute coefficients a,b,c of the previous quadratic equation (a*t^2+b*t+c=0)
    ASSOCIATE(Term => V2 - V1)
    a = DOT_PRODUCT(Term,Term)
    END ASSOCIATE

    ! || particle trajectories: if a.EQ.0 then b.EQ.0 and so the previous equation has no solution: go to the next pair
    IF (a.EQ.0.) CYCLE

    ! Cauchy-Schwarz inequality: <x,y> + <y,x> <= 2 |<x,y>| <= 2|x||y|
    P1 = PartData_Shared(PART_OLDV,iPart1)
    P2 = PartData_Shared(PART_OLDV,iPart2)
    b = 2. * DOT_PRODUCT(ABS(V1 - V2), ABS(P1-P2))
    ! b = 2. * DOT_PRODUCT(P1-P2,PartData_Shared(PART_POSV,iPart1) - PartData_Shared(PART_POSV,iPart2)-(P1 - P2))
    ASSOCIATE(Term => P1 - P2)
    c = DOT_PRODUCT(Term,Term) - 0.25*(PartData_Shared(PART_DIAM,iPart2) + PartData_Shared(PART_DIAM,iPart1))**2.
    END ASSOCIATE

    ! if delta < 0, no collision is possible for this pair
    delta = b * b - 4. * a * c
    IF (delta.LT.0.) CYCLE

    ! a collision is possible, determine when relatively to (tStage+dtloc)
    dtColl = MAX(MIN((-b + SQRT(delta)) / (2. * a),1.),MIN((-b - SQRT(delta)) / (2. * a),1.)) * dtLoc
    ! the collision is only valid if it occurred between tStage and tStage+dtLoc
    IF (dtColl.LT.0. .OR. dtColl.GT.dtLoc) CYCLE

    ! the collision is only valid if it occurred before a BC intersection
    ! > here, the particle radius needs to be considered since the collision is on the particle surface but the BC intersection
    ! > is determined based on the particle center of mass
    IF (dtColl.GT.PartBC_Shared(iPart1)) CYCLE
    IF (dtColl.GT.PartBC_Shared(iPart2)) CYCLE

    !IF (VECNORM(P1 + dtColl * V1 - (P2 + dtColl * V2)) &
      !.GT. 0.5*(PartData_Shared(PART_DIAM,iPart2) + PartData_Shared(PART_DIAM,iPart1))) CYCLE

    ! collision is only valid the particles approach each other
    ! FIXME: NOT TRUE, we can collide on with acute angles
    ASSOCIATE(N1 => V1 &
             ,N2 => V2)
    pColl1 = PartData_Shared(PART_POSV,iPart1) - (dtLoc-dtColl) * V1
    pColl2 = PartData_Shared(PART_POSV,iPart2) - (dtLoc-dtColl) * V2
    IF (DOT_PRODUCT(pColl1 - P1,N1).LE.0) CYCLE
    IF (DOT_PRODUCT(pColl2 - P2,N2).LE.0) CYCLE
    END ASSOCIATE

    ! apply the power
    CALL ComputeHardSphereCollision(iPart1,iPart2,dtLoc-dtColl,dtLoc)
  END DO ! iPart1

#if USE_LOADBALANCE
  CALL LBPauseTime(LB_COLLISION,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO ! iElem

DEALLOCATE( PEM%pStart   &
          , PEM%pNumber  &
          , PEM%pNext    &
          , PEM%pEnd)

! Done with particle collisions. Deallocate mappings
SDEALLOCATE(PEM2PartID)

END SUBROUTINE ComputeParticleCollisions


SUBROUTINE ComputeHardSphereCollision(iPart1,iPart2,mdtColl,dtLoc)
!===================================================================================================================================
!> Apply hard sphere collision model
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
USE MOD_Utils                   ,ONLY: ALMOSTZERO
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals        ,ONLY: UNITVECTOR
USE MOD_Particle_Localization   ,ONLY: SinglePointToElement
USE MOD_Particle_Mesh_Tools     ,ONLY: GetCNElemID
USE MOD_Particle_Vars           ,ONLY: PartState,Species,LastPartPos,PDM,PEM
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: iPart1, iPart2
REAL   ,INTENT(IN)                 :: mdtColl
REAL   ,INTENT(IN)                 :: dtLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: LocPartID1,LocPartID2
INTEGER                            :: firstCNElem,lastCNElem
REAL                               :: n_loc(3),t_loc(3)
REAL                               :: P1_old(3),P2_old(3)
REAL                               :: V1(    3),V2(    3)
REAL                               :: v_rel_contactpoint(3)
REAL                               :: vt_rel_contactpoint(3),vt_rel_contactpoint_norm
REAL                               :: m_red,m1,m2
REAL                               :: Jn,Jt,J(3)
#if USE_PARTROT
REAL                               :: Jxn_loc(3)
#endif /*USE_PARTROT*/
!===================================================================================================================================

! Calculate the CN elems
firstCNElem = GetCNElemID(offsetElem+1)
 lastCNElem = GetCNElemID(offsetElem+nElems)

! Move the particles back to their positions at the collision time
P1_old = PartData_Shared(PART_OLDV,iPart1)
P2_old = PartData_Shared(PART_OLDV,iPart2)

! Reconstruct the particle velocity
V1 = (PartData_Shared(PART_POSV,iPart1) - PartData_Shared(PART_OLDV,iPart1)) / dtLoc
V2 = (PartData_Shared(PART_POSV,iPart2) - PartData_Shared(PART_OLDV,iPart2)) / dtLoc

LocPartID1 = PEM2PartID(iPart1)

! Compute normal vector from iPart1 to iPart2 at the time of contact
n_loc = UNITVECTOR(P2_old - P1_old + mdtColl*(V2 - V1))

! Compute reduced mass of particle pair
m1 = MASS_SPHERE(Species(INT(PartData_Shared(PP_nVarPart+1,iPart1)))%DensityIC,PartData_Shared(PART_DIAM,iPart1))
! FIXME: density of iPart2
m2 = MASS_SPHERE(Species(INT(PartData_Shared(PP_nVarPart+1,iPart2)))%DensityIC,PartData_Shared(PART_DIAM,iPart2))
m_red = m1*m2/(m1+m2)

! Compute normal component Jn of impulsive contact force J
Jn = - m_red * (1+PartCollisionModel%e) * DOT_PRODUCT(V1 - V2,n_loc)
! Alternative formulation with same result
! > https://physics.stackexchange.com/a/256992
! Jn = (1+PartCollisionModel%e) * DOT_PRODUCT(n_loc,V2 - V1)/(1/m1 + 1/m2)

! Compute friction (Jt) if enabled
IF (PartCollisionModel%Friction) THEN
  ! Relative velocity of the contact point just before the collision in iPart2's reference frame
  v_rel_contactpoint = V1 - V2
#if USE_PARTROT
  ! Take particle rotation into account
  v_rel_contactpoint = v_rel_contactpoint + &
                       PartData_Shared(PART_DIAM,iPart1) * 0.5 * CROSSPRODUCT(PartData_Shared(PART_AMOMV,iPart1),n_loc) + &
                       PartData_Shared(PART_DIAM,iPart2) * 0.5 * CROSSPRODUCT(PartData_Shared(PART_AMOMV,iPart2),n_loc)
#endif /*USE_PARTROT*/

  ! Compute tangent unit vector t_loc that is in the same plan as n_loc and v_rel_contactpoint
  vt_rel_contactpoint      = v_rel_contactpoint - DOT_PRODUCT(v_rel_contactpoint,n_loc) * n_loc
  vt_rel_contactpoint_norm = NORM2(vt_rel_contactpoint)
  IF (ALMOSTZERO(vt_rel_contactpoint_norm)) THEN; t_loc = 0.
  ELSE;                                           t_loc = vt_rel_contactpoint / vt_rel_contactpoint_norm
  END IF

  ! as particles can slide on each other, the fact that they can stop sliding during the collision due to friction is taken into account
  Jt = MAX(PartCollisionModel%f * Jn, - 2. / 7. * m_red * vt_rel_contactpoint_norm)

  ! Compute impulsive contact force J
  J = Jn * n_loc + Jt * t_loc
ELSE
  ! Compute impulsive contact force J
  J = Jn * n_loc
END IF

! Apply jump relation to momentum equation, i.e., update velocity of particles
PartState(PART_VELV,LocPartID1) = V1 + J / m1

IF (iPart2.GE.PartInt_Shared(1,firstCNElem)+1 .AND. &
    iPart2.LE.PartInt_Shared(2, lastCNElem))  THEN
  LocPartID2 = PEM2PartID(iPart2)
  PartState(PART_VELV,LocPartID2) = V2 - J / m2
END IF

! Apply jump relation also to angular momentum if particles can rotate
#if USE_PARTROT
Jxn_loc = CROSSPRODUCT(n_loc, J)
PartState(PART_AMOMV,LocPartID1) = PartData_Shared(PART_AMOMV,iPart1) + Jxn_loc * 0.5 / (0.1 * m1 * PartData_Shared(PART_DIAM,iPart1))
IF (iPart2.LE.nComputeNodeParts) THEN
  PartState(PART_AMOMV,LocPartID2) = PartData_Shared(PART_AMOMV,iPart2) + Jxn_loc * 0.5 / (0.1 * m2 * PartData_Shared(PART_DIAM,iPart2))
END IF
#endif /*USE_PARTROT*/

! for LSERK algorithm
PDM%IsNewPart(LocPartID1) = .TRUE.
IF (iPart2.GE.PartInt_Shared(1,firstCNElem)+1 .AND. &
    iPart2.LE.PartInt_Shared(2, lastCNElem))        &
  PDM%IsNewPart(LocPartID2) = .TRUE.

! Update the particle's position at time tStage+dtloc
PartState(  PART_POSV,LocPartID1) = P1_old + mdtColl * PartState(PART_VELV,LocPartID1)
! Update the particle host element
PEM%lastElement(LocPartID1)       = SinglePointToElement(P1_old,doHalo=.TRUE.)
IF (PEM%lastElement(LocPartID1).EQ.-1) CALL Abort(__STAMP__,'LastElemID == -1!')
LastPartPos(PART_POSV,LocPartID1) = P1_old
IF (iPart2.GE.PartInt_Shared(1,firstCNElem)+1 .AND. &
    iPart2.LE.PartInt_Shared(2, lastCNElem))  THEN
  PartState(  PART_POSV,LocPartID2) = P2_old + mdtColl * PartState(PART_VELV,LocPartID2)
  ! Update the particle host element
  PEM%lastElement(LocPartID2)       = SinglePointToElement(P2_old,doHalo=.TRUE.)
  IF (PEM%lastElement(LocPartID1).EQ.-1) CALL Abort(__STAMP__,'LastElemID == -1!')
  LastPartPos(PART_POSV,LocPartID2) = P2_old
END IF

! Count number of collisions, only count twice if both particles are on the local processor
CollisionnLoc = CollisionnLoc + MERGE(2,1,iPart2.GE.PartInt_Shared(1,firstCNElem)+1 .AND. iPart2.LE.PartInt_Shared(2,lastCNElem))

#if USE_PARTTEMP
! TODO
#endif

END SUBROUTINE ComputeHardSphereCollision


SUBROUTINE ComputeParticleTracingIntersection(PartID,ElemID,CNElemID,dtLoc,dtBC)
!===================================================================================================================================
!> Compute particle tracing path until first BC intersection
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars                ,ONLY: SideInfo_Shared
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Localization    ,ONLY: PARTHASMOVED,SinglePointToElement
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNSideID,GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemRadiusNGeo
USE MOD_Particle_Intersection    ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection    ,ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection    ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection    ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Surfaces_Vars   ,ONLY: SideType
! LOCAL VARIABLES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID
INTEGER,INTENT(IN)            :: ElemID
INTEGER,INTENT(IN)            :: CNElemID
REAL,INTENT(IN)               :: dtLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: dtBC
!-----------------------------------------------------------------------------------------------------------------------------------
! Counters
INTEGER                       :: iLocSide
! Sides
INTEGER                       :: SideID,CNSideID,flip
! Particles
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
! Tracking
REAL                          :: locAlpha(1:6),xi,eta
REAL                          :: PartState(1:3),LastPartPos(1:3)
LOGICAL                       :: isHit
LOGICAL                       :: isCriticalParallelInFace
!-----------------------------------------------------------------------------------------------------------------------------------

dtBC = HUGE(1.)

! Calculate particle trajectory
PartTrajectory       = PartData_Shared(PART_VELV,PartID) * dtLoc
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))
! Reconstruct the particle position since PartID is for SHM particles
PartState            = PartData_Shared(PART_POSV,PartID)
LastPartPos          = PartData_Shared(PART_POSV,PartID) - PartData_Shared(PART_VELV,PartID)*dtLoc

! Check if the particle moved at all. If not, tracking is done
IF (.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(CNElemID))) RETURN

PartTrajectory       = PartTrajectory/lengthPartTrajectory
locAlpha             = -1.

DO iLocSide = 1,6
  SideID   = GetGlobalNonUniqueSideID(ElemID,iLocSide)
  CNSideID = GetCNSideID(SideID)
  ! If the side is positive, then the element has the actual side
  ! and neighbour element has the negative one which has to be flipped

  ! BezierControlPoints are now built in cell local system. Hence, sides have always the flip from the shared SideInfo
  flip = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

  isCriticalParallelInFace = .FALSE.

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT)
      CALL ComputePlanarRectIntersection(  isHit                      = isHit                &
                                        ,  PartTrajectory             = PartTrajectory       &
                                        ,  lengthPartTrajectory       = lengthPartTrajectory &
                                        ,  LastPartPos                = LastPartPos          &
                                        ,  alpha                      = locAlpha(iLocSide)   &
                                        ,  xi                         = xi                   &
                                        ,  eta                        = eta                  &
#if CODE_ANALYZE
                                        ,  PartID                     = PartID               &
#endif /*CODE_ANALYZE*/
                                        ,  flip                       = flip                 &
                                        ,  SideID                     = SideID               &
                                        ,  opt_CriticalParallelInSide = isCriticalParallelInFace)
    CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(    isHit                      = isHit                &
                                        ,  PartTrajectory             = PartTrajectory       &
                                        ,  lengthPartTrajectory       = lengthPartTrajectory &
                                        ,  PartState                  = PartState            &
                                        ,  LastPartPos                = LastPartPos          &
                                        ,  alpha                      = locAlpha(iLocSide)   &
                                        ,  xitild                     = xi                   &
                                        ,  etatild                    = eta                  &
                                        ,  PartID                     = PartID               &
                                        ,  flip                       = flip                 &
                                        ,  SideID                     = SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection(isHit                      = isHit                &
                                        ,  PartTrajectory             = PartTrajectory       &
                                        ,  lengthPartTrajectory       = lengthPartTrajectory &
                                        ,  LastPartPos                = LastPartPos          &
                                        ,  alpha                      = locAlpha(iLocSide)   &
                                        ,  xi                         = xi                   &
                                        ,  eta                        = eta                  &
                                        ,  PartID                     = PartID               &
                                        ,  flip                       = flip                 &
                                        ,  SideID                     = SideID               &
                                        ,  opt_CriticalParallelInSide = isCriticalParallelInFace)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(      isHit                      = isHit                &
                                        ,  PartTrajectory             = PartTrajectory       &
                                        ,  lengthPartTrajectory       = lengthPartTrajectory &
                                        ,  PartState                  = PartState            &
                                        ,  LastPartPos                = LastPartPos          &
                                        ,  alpha                      = locAlpha(iLocSide)   &
                                        ,  xi                         = xi                   &
                                        ,  eta                        = eta                  &
                                        ,  PartID                     = PartID               &
                                        ,  flip                       = flip                 &
                                        ,  SideID                     = SideID               &
                                        ,  opt_CriticalParallelInSide = isCriticalParallelInFace)
    CASE DEFAULT
      CALL Abort(__STAMP__,' Missing required side-data. Please increase halo region. ',SideID)
  END SELECT
END DO ! iLocSide = 1,6

IF (MAXVAL(locAlpha).GT.0) dtBC = MINVAL(locAlpha,locAlpha.GT.0)

END SUBROUTINE ComputeParticleTracingIntersection


SUBROUTINE ComputeParticleRefmappingIntersection(PartID,CNElemID,firstSide,lastSide,nLocSides,dtLoc,dtBC)
!===================================================================================================================================
!> Compute particle tracing path until first BC intersection
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars                ,ONLY: SideInfo_Shared
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Localization    ,ONLY: PARTHASMOVED
USE MOD_Particle_Mesh_Vars       ,ONLY: SideBCMetrics
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemRadiusNGeo
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNSideID
USE MOD_Particle_Intersection    ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection    ,ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection    ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection    ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Surfaces_Vars   ,ONLY: SideType
! LOCAL VARIABLES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID
INTEGER,INTENT(IN)            :: CNElemID
INTEGER,INTENT(IN)            :: firstSide
INTEGER,INTENT(IN)            :: lastSide
INTEGER,INTENT(IN)            :: nlocSides
REAL,INTENT(IN)               :: dtLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: dtBC
!-----------------------------------------------------------------------------------------------------------------------------------
! Counters
INTEGER                       :: iLocSide
! Sides
INTEGER                       :: SideID,CNSideID,flip
! Particles
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
! Tracking
REAL                          :: localpha(firstSide:lastSide),xi,eta
REAL                          :: PartState(1:3),LastPartPos(1:3)
LOGICAL                       :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------

dtBC = HUGE(1.)

! check if element is a BC element
IF (nLocSides.EQ.0) RETURN

! Calculate particle trajectory
PartTrajectory       = PartData_Shared(PART_VELV,PartID) * dtLoc
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))
! Reconstruct the particle position since PartID is for SHM particles
PartState            = PartData_Shared(PART_POSV,PartID)
LastPartPos          = PartData_Shared(PART_POSV,PartID) - PartData_Shared(PART_VELV,PartID)*dtLoc

! Check if the particle moved at all. If not, tracking is done
IF (.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(CNElemID))) RETURN

PartTrajectory       = PartTrajectory/lengthPartTrajectory
locAlpha             = -1.

DO iLocSide = firstSide,LastSide

  ! SideBCMetrics is sorted by distance. stop if the first side is out of range
  IF (SideBCMetrics(BCSIDE_DISTANCE,iLocSide).GT.lengthPartTrajectory) CYCLE

  ! side potentially in range (halo_eps)
  SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,iLocSide))
  CNSideID = GetCNSideID(SideID)

  ! BezierControlPoints are now built in cell local system. Hence, sides have always the flip from the shared SideInfo
  flip = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT)
      CALL ComputePlanarRectIntersection(  isHit                = isHit                &
                                        ,  PartTrajectory       = PartTrajectory       &
                                        ,  lengthPartTrajectory = lengthPartTrajectory &
                                        ,  LastPartPos          = LastPartPos          &
                                        ,  alpha                = locAlpha(iLocSide)   &
                                        ,  xi                   = xi                   &
                                        ,  eta                  = eta                  &
#if CODE_ANALYZE
                                        ,  PartID               = PartID               &
#endif /*CODE_ANALYZE*/
                                        ,  flip                 = flip                 &
                                        ,  SideID               = SideID)
    CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(    isHit                = isHit                &
                                        ,  PartTrajectory       = PartTrajectory       &
                                        ,  lengthPartTrajectory = lengthPartTrajectory &
                                        ,  PartState            = PartState            &
                                        ,  LastPartPos          = LastPartPos          &
                                        ,  alpha                = locAlpha(iLocSide)   &
                                        ,  xitild               = xi                   &
                                        ,  etatild              = eta                  &
                                        ,  PartID               = PartID               &
                                        ,  flip                 = flip                 &
                                        ,  SideID               = SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection(isHit                = isHit                &
                                        ,  PartTrajectory       = PartTrajectory       &
                                        ,  lengthPartTrajectory = lengthPartTrajectory &
                                        ,  LastPartPos          = LastPartPos          &
                                        ,  alpha                = locAlpha(iLocSide)   &
                                        ,  xi                   = xi                   &
                                        ,  eta                  = eta                  &
                                        ,  PartID               = PartID               &
                                        ,  flip                 = flip                 &
                                        ,  SideID               = SideID)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(      isHit                = isHit                &
                                        ,  PartTrajectory       = PartTrajectory       &
                                        ,  lengthPartTrajectory = lengthPartTrajectory &
                                        ,  PartState            = PartState            &
                                        ,  LastPartPos          = LastPartPos          &
                                        ,  alpha                = locAlpha(iLocSide)   &
                                        ,  xi                   = xi                   &
                                        ,  eta                  = eta                  &
                                        ,  PartID               = PartID               &
                                        ,  flip                 = flip                 &
                                        ,  SideID               = SideID)
  END SELECT
END DO ! iLocSide = firstSide,LastSide

IF (MAXVAL(locAlpha).GT.0) dtBC = MINVAL(locAlpha,locAlpha.GT.0)

END SUBROUTINE ComputeParticleRefmappingIntersection
#endif /*PARTICLES_COUPLING == 4*/

END MODULE MOD_Particle_Collision_Method
