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
! INTERFACE ComputeHardSphereCollision
!   MODULE PROCEDURE ComputeHardSphereCollision
! END INTERFACE

INTERFACE ComputeParticleCollisions
  MODULE PROCEDURE ComputeParticleCollisions
END INTERFACE

PUBLIC :: ComputeParticleCollisions
!===================================================================================================================================

CONTAINS

SUBROUTINE ComputeParticleCollisions(dtLoc)
!===================================================================================================================================
!> Compute particle collisions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars                ,ONLY: nElems,offsetElem
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemToBCSides
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNElemID,GetGlobalElemID
USE MOD_Particle_Vars            ,ONLY: PEM
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars,          ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars,          ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SHARED
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)               :: dtLoc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: PartBCdt
INTEGER                       :: iElem,ElemID
INTEGER                       :: iPart,iPart1,iPart2
INTEGER                       :: CNElemID,CNNeighElemID
INTEGER                       :: NeighElemID,iCNNeighElem
! Trajectories and position
REAL                          :: P1(3),P2(3)
! Collisions
REAL                          :: a,b,c,delta,dtColl
! Neighbor particles
LOGICAL,ALLOCATABLE           :: PartColl(:)
!===================================================================================================================================

! safety check, dt > 1 currently not supported
IF (dtLoc.GT.1) CALL Abort(__STAMP__,'Time steps greater than 1 are unsupported with particle collisions')

! nullify the BC intersection times
IF (myComputeNodeRank.EQ.0) PartBC_Shared = -1.
CALL BARRIER_AND_SYNC(PartBC_Shared_Win,MPI_COMM_SHARED)

! dtColl,dtBC,PartID

! loop over internal elements
DO iElem = offsetElem+1,offsetElem+nElems
  ! skip if there is no particles inside the current element
  IF (PartInt_Shared(2,iElem).EQ.PartInt_Shared(1,iElem)) CYCLE

  SELECT CASE(TrackingMethod)
    CASE (TRIATRACKING)
      ! FIXME
      CALL Abort(__STAMP__,'TODO: TriaTracking')
    CASE (TRACING)
      ! loop over all particles in the element
      DO iPart = PartInt_Shared(1,iElem)+1,PartInt_Shared(2,iElem)
        CNElemID = GetCNElemID(iElem)

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
      DO iPart = PartInt_Shared(1,iElem)+1,PartInt_Shared(2,iElem)
        CNElemID = GetCNElemID(iElem)

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

IF (myComputeNodeRank.EQ.0) THEN
  ! ALLOCATE(MPI_COMM_LEADERS_REQUEST(1:2))
  CALL MPI_ALLGATHERV(MPI_IN_PLACE   ,0                           ,MPI_DATATYPE_NULL    &
                     ,PartBC_Shared,recvcountPartInt,displsPartInt,MPI_DOUBLE_PRECISION &
                     ,MPI_COMM_LEADERS_SHARED,iError)
  ! DEALLOCATE(MPI_COMM_LEADERS_REQUEST)
  DEALLOCATE(displsPartInt)
  DEALLOCATE(recvcountPartInt)
END IF ! myComputeNodeRank.EQ.0
CALL BARRIER_AND_SYNC(PartBC_Shared_Win,MPI_COMM_SHARED)

! Build a mapping from global particle IDs to neighbor particle IDs
! > Count number of particles in own and neighbor elements
nProcParts = 0
DO iElem = 1,nProcNeighElems
  ElemID = NeighElemsProc(iElem)

  ! > Map global offset to neighbor offset
  offsetNeighElemPart(ElemID) = PartInt_Shared(1,ElemID) - nProcParts

  !> Increment the counter by the current number of particles
  nProcParts = nProcParts + PartInt_Shared(2,ElemID)-PartInt_Shared(1,ElemID)
END DO

! > Allocate an array to hold the mapping
ALLOCATE(PartColl(nProcParts))
PartColl = .FALSE.

! > Allocate local collision arrays
! ALLOCATE(PartColl(PartInt_Shared(1,offsetElem+1)+1:PartInt_Shared(2,offsetEle+nElems)))

! loop over internal elements
DO iElem = offsetElem+1,offsetElem+nElems
  ! skip if there is no particles inside the current element
  IF (PartInt_Shared(2,iElem).EQ.PartInt_Shared(1,iElem)) CYCLE

  ! loop over all particles in the element
  DO iPart1 = PartInt_Shared(1,iElem)+1,PartInt_Shared(2,iElem)
    ! Ignore already collided particles
    IF (PartColl(iPart1 - offsetNeighElemPart(iElem))) CYCLE

    ! loop over neighbor elements (includes own element)
PartLoop: DO iCNNeighElem = Neigh_offsetElem(iElem)+1,Neigh_offsetElem(iElem)+Neigh_nElems(iElem)
      CNNeighElemID = CNElem2CNNeighElem(iCNNeighElem)
      NeighElemID   = GetGlobalElemID(CNNeighElemID)

      ! loop over all particles in the element
      DO iPart2 = PartInt_Shared(1,NeighElemID)+1,PartInt_Shared(2,NeighElemID)
        ! Ignore already collided particles
        IF (PartColl(iPart2 - offsetNeighElemPart(NeighElemID))) CYCLE

        ! FIXME: What about periodic?

        ! Check if particles will collide
        P1 = PartData_Shared(PART_POSV,iPart1) - PartData_Shared(PART_VELV,iPart1)*dtLoc
        P2 = PartData_Shared(PART_POSV,iPart2) - PartData_Shared(PART_VELV,iPart2)*dtLoc
        ! compute coefficients a,b,c of the previous quadratic equation (a*t^2+b*t+c=0)
        ASSOCIATE(Term => PartData_Shared(PART_VELV,iPart2) - PartData_Shared(PART_VELV,iPart1))
        a = DOT_PRODUCT(Term,Term)
        END ASSOCIATE

        ! if a.EQ.0 then b.EQ.0 and so the previous equation has no solution: go to the next pair
        IF (a.EQ.0.) CYCLE

        b = 2. * DOT_PRODUCT(PartData_Shared(PART_VELV,iPart2) - PartData_Shared(PART_VELV,iPart1),P2 - P1)
        ASSOCIATE(Term => P2 - P1 - (PartData_Shared(PART_DIAM,iPart2) + PartData_Shared(PART_DIAM,iPart1)))
        c = DOT_PRODUCT(Term,Term) ** 2. / 4.
        END ASSOCIATE

        ! if delta < 0, no collision is possible for this pair
        delta = b * b - 4. * a * c
        IF (delta.LT.0.) CYCLE

        ! a collision is possible, determine when relatively to (tStage+dtloc)
        dtColl = (b + SQRT(delta)) / (2. * a) * dtLoc
        ! the collision is only valid if it occurred between tStage and tStage+dtLoc
        IF (dtColl.LT.0. .OR. dtColl.GT.dtLoc) CYCLE

        ! the collision is only valid if it occurred before a BC intersection
        IF (dtColl.GT.PartBC_Shared(iPart1)) CYCLE
        IF (dtColl.GT.PartBC_Shared(iPart2)) CYCLE

        ! collision found
        PartColl(iPart1 - offsetNeighElemPart(iElem      )) = .TRUE.
        PartColl(iPart2 - offsetNeighElemPart(NeighElemID)) = .TRUE.

        ! FIXME: apply the power

        EXIT PartLoop
        ! stop 1





      END DO ! iPart2
    END DO PartLoop ! iCNNeighElem
  END DO ! iPart1
END DO ! iElem

SDEALLOCATE(PartColl)






  !     ! FIXME: What about periodic?
  !
  !     ! Check if particles will collide
  !     P1 = PartState(PART_POSV,PartID1) - PartState(PART_VELV,PartID1)*dtLoc
  !     P2 = PartState(PART_POSV,PartID2) - PartState(PART_VELV,PartID2)*dtLoc
  !     ! compute coefficients a,b,c of the previous quadratic equation (a*t^2+b*t+c=0)
  !     a = DOTPRODUCT(PartState(PART_VELV,iPart2) - PartState(PART_VELV,PartID1))
  !
  !     ! if a.EQ.0 then b.EQ.0 and so the previous equation has no solution: go to the next pair
  !     IF (a.EQ.0.) CYCLE
  !
  !     b = 2. * DOT_PRODUCT(PartState(PART_VELV,iPart2) - PartState(PART_VELV,PartID1),P2 - P1)
  !     ASSOCIATE(Term => (P2 - P1 - (PartState(PART_DIAM,iPart2) + PartState(PART_DIAM,PartID1))))
  !     c = DOTPRODUCT(Term,Term) ** 2. / 4.
  !     END ASSOCIATE
  !
  !     ! if delta < 0, no collision is possible for this pair
  !     delta = b * b - 4. * a * c
  !     IF (delta.LT.0.) CYCLE
  !
  !     ! a collision is possible, determine when relatively to (tStage+dtloc)
  !     dtColl = (b + SQRT(delta)) / (2. * a)
  !     ! the collision is only valuable if it occurred between tStage and tStage+dtLoc
  !     IF (dtColl.LT.0. .OR. dtColl.GT.dtLoc) CYCLE
  !
  !     ! Store the collision
  !     PartColl(1,PartID1) = dtColl
  !     PartColl(2,PartID1) = PartID2
  !     PartColl(1,PartID2) = dtColl
  !     PartColl(2,PartID2) = PartID1

  ! FIXBME:
  ! PartBC_Shared





















  ! ! skip if there is no particles inside the current element
  ! IF (PEM%pNumber(iElem).EQ.0) CYCLE
  !
  ! ! loop over all particles in the element, linked list!
  ! DO iPart1 = PEM%pStart(iElem),PEM%pStart(iElem)+PEM%pNumber(iElem)
  !   PartID1 = PEM%pStart(iElem)
  !
  !   ! Check collisions with particles in the same element
  !   DO iPart2 = PEM%pStart(iElem),PEM%pStart(iElem)+PEM%pNumber(iElem)
  !     PartID2 = PEM%pStart(iElem)
  !
  !     ! Particle not checked, set pointer to next particle
  !     IF (PartID1.GE.PartID2) THEN
  !       PartID2 = PEM%pNext(PartID2)
  !       CYCLE
  !     END IF
  !
  !     ! > Check for collisions
  !     ! Ignore already collided particles
  !     IF (PartColl(1,PartID1).GT.0) CYCLE
  !     IF (PartColl(1,PartID2).GT.0) CYCLE
  !
  !     ! FIXME: What about periodic?
  !
  !     ! Check if particles will collide
  !     P1 = PartState(PART_POSV,PartID1) - PartState(PART_VELV,PartID1)*dtLoc
  !     P2 = PartState(PART_POSV,PartID2) - PartState(PART_VELV,PartID2)*dtLoc
  !     ! compute coefficients a,b,c of the previous quadratic equation (a*t^2+b*t+c=0)
  !     a = DOTPRODUCT(PartState(PART_VELV,iPart2) - PartState(PART_VELV,PartID1))
  !
  !     ! if a.EQ.0 then b.EQ.0 and so the previous equation has no solution: go to the next pair
  !     IF (a.EQ.0.) CYCLE
  !
  !     b = 2. * DOT_PRODUCT(PartState(PART_VELV,iPart2) - PartState(PART_VELV,PartID1),P2 - P1)
  !     ASSOCIATE(Term => (P2 - P1 - (PartState(PART_DIAM,iPart2) + PartState(PART_DIAM,PartID1))))
  !     c = DOTPRODUCT(Term,Term) ** 2. / 4.
  !     END ASSOCIATE
  !
  !     ! if delta < 0, no collision is possible for this pair
  !     delta = b * b - 4. * a * c
  !     IF (delta.LT.0.) CYCLE
  !
  !     ! a collision is possible, determine when relatively to (tStage+dtloc)
  !     dtColl = (b + SQRT(delta)) / (2. * a)
  !     ! the collision is only valuable if it occurred between tStage and tStage+dtLoc
  !     IF (dtColl.LT.0. .OR. dtColl.GT.dtLoc) CYCLE
  !
  !     ! Store the collision
  !     PartColl(1,PartID1) = dtColl
  !     PartColl(2,PartID1) = PartID2
  !     PartColl(1,PartID2) = dtColl
  !     PartColl(2,PartID2) = PartID1
  !
  !     ! Move particle index to next particle
  !     PartID2 = PEM%pNext(PartID2)
  !   END DO
  !
  !   ! Check collisions with particles in neighbor elements
  !   DO iNeighElem = Neigh_offsetElem(iElem)+1:Neigh_offsetElem(iElem)+Neigh_nElems(iElem)
  !
  !   END DO
  !   PartID1 = PEM%pNext(PartID1)
  ! END DO ! iPart1 = 1,PEM%pNumber(iElem)







! ! register info about the detected collision
! PartCollisionProc%nCollisions = PartCollisionProc%nCollisions + 1
! IF (PartCollisionProc%nCollisions.GT.MAXPartCollisionPairs) CALL Abort(__STAMP__, "Too many candidate collision pairs to save. Increase Part-Collision-MAXPartCollisionPairs!")
! #if USE_MPI
! PartCollisionProc%type  (PartCollisionProc%nCollisions) = 0
! #endif /*USE_MPI*/
! PartCollisionProc%iPartPeriodicShift(PartCollisionProc%nCollisions) = iPartPeriodicShift
! PartCollisionProc%iPart1            (PartCollisionProc%nCollisions) = iPart1
! PartCollisionProc%iPart2            (PartCollisionProc%nCollisions) = iPart2
! PartCollisionProc%dtColl            (PartCollisionProc%nCollisions) = dtcand

DEALLOCATE( PEM%pStart   &
          , PEM%pNumber  &
          , PEM%pNext    &
          , PEM%pEnd)

! Done with particle collisions. Deallocate mappings
SDEALLOCATE(PEM2PartID)

END SUBROUTINE ComputeParticleCollisions


SUBROUTINE ComputeHardSphereCollision(iPart1,iPart2,dtColl)
!===================================================================================================================================
!> Apply hard sphere collision model
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Utils                   , ONLY: ALMOSTZERO
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals        , ONLY: UNITVECTOR
USE MOD_Particle_Vars           , ONLY: PartState, Species, PartSpecies, LastPartPos!, PDM
!USE MOD_TimeDisc_Vars           , only: CurrentStage, t, iter, dt, RKc
! USE MOD_Particle_Mesh_Vars      , ONLY: GEO
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: iPart1, iPart2
REAL   , INTENT(IN)                :: dtColl
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: n_loc(3),t_loc(3)
REAL                               :: v_rel_contactpoint(3)
REAL                               :: vt_rel_contactpoint(3),vt_rel_contactpoint_norm
REAL                               :: m_red,m1,m2
REAL                               :: Jn, Jt,J(3)
#if USE_PARTROT
REAL                               :: Jxn_loc(3)
#endif /*USE_PARTROT*/
! INTEGER                            :: iPart1Track, iPart2Track
! INTEGER                            :: iPartPeriodicShiftLoc, jPeriodic
!===================================================================================================================================

! ! apply periodic vectors (if there are some)
! IF (GEO%nPeriodicVectors.GT.0) THEN
!   iPartPeriodicShiftLoc = PartCollisionProc%iPartPeriodicShift(iCollision)
!   DO jPeriodic=1,GEO%nPeriodicVectors
!     PartState(PART_POSV,iPart2) = PartState(PART_POSV,iPart2) - REAL(MOD(iPartPeriodicShiftLoc, 3) - 1) * GEO%PeriodicVectors(1:3,jPeriodic)
!     iPartPeriodicShiftLoc = iPartPeriodicShiftLoc / 3
!   END DO
! END IF

! Move the particles back to their positions at the collision time
LastPartPos(PART_POSV,iPart1) = PartState(PART_POSV,iPart1) - dtColl * PartState(PART_VELV,iPart1)
LastPartPos(PART_POSV,iPart2) = PartState(PART_POSV,iPart2) - dtColl * PartState(PART_VELV,iPart2)

! Compute normal vector from iPart1 to iPart2
n_loc = UNITVECTOR(PartState(PART_POSV,iPart2) - PartState(PART_POSV,iPart1))

! Relative velocity of the contact point just before the collision in iPart2's reference frame
v_rel_contactpoint = PartState(PART_VELV,iPart1) - PartState(PART_VELV,iPart2)
#if USE_PARTROT
! Take particle rotation into account
v_rel_contactpoint = v_rel_contactpoint + &
                     PartState(PART_DIAM,iPart1) * 0.5 * CROSSPRODUCT(PartState(PART_AMOMV,iPart1),n_loc) + &
                     PartState(PART_DIAM,iPart2) * 0.5 * CROSSPRODUCT(PartState(PART_AMOMV,iPart2),n_loc)
#endif /*USE_PARTROT*/

! Compute reduced mass of particle pair
m1 = MASS_SPHERE(Species(PartSpecies(iPart1))%DensityIC,PartState(PART_DIAM,iPart1))
m2 = MASS_SPHERE(Species(PartSpecies(iPart2))%DensityIC,PartState(PART_DIAM,iPart2))
m_red = m1*m2/(m1+m2)

! Compute normal component Jn of impulsive contact force J
Jn = - m_red * (1+PartCollisionModel%e) * DOT_PRODUCT(PartState(PART_VELV,iPart1) - PartState(PART_VELV,iPart2),n_loc)

! Compute tangent unit vector t_loc that is in the same plan as n_loc and v_rel_contactpoint
vt_rel_contactpoint = v_rel_contactpoint - DOT_PRODUCT(v_rel_contactpoint,n_loc) * n_loc
vt_rel_contactpoint_norm = NORM2(vt_rel_contactpoint)
IF (ALMOSTZERO(vt_rel_contactpoint_norm)) THEN; t_loc = 0.
ELSE;                                           t_loc = vt_rel_contactpoint / vt_rel_contactpoint_norm
END IF

! Compute friction (Jt) if enabled
IF (PartCollisionModel%Friction .AND. vt_rel_contactpoint_norm.NE.0.) THEN
  ! as particles can slide on each other, the fact that they can stop sliding during the collision due to friction is taken into account
  Jt = MAX(PartCollisionModel%f * Jn, - 2. / 7. * m_red * vt_rel_contactpoint_norm)
ELSE
  Jt = 0.
END IF

! Compute impulsive contact force J
J = Jn * n_loc + Jt * t_loc

! Apply jump relation to momentum equation, i.e., update velocity of particles
PartState(PART_VELV,iPart1) = PartState(PART_VELV,iPart1) + J / m1
PartState(PART_VELV,iPart2) = PartState(PART_VELV,iPart2) - J / m2

! Apply jump relation also to angular momentum if particles can rotate
#if USE_PARTROT
Jxn_loc = CROSSPRODUCT(n_loc, J)
PartState(PART_AMOMV,iPart1) = PartState(PART_AMOMV,iPart1) + Jxn_loc * 0.5 / (0.1 * m1 * PartState(PART_DIAM,iPart1))
PartState(PART_AMOMV,iPart2) = PartState(PART_AMOMV,iPart2) + Jxn_loc * 0.5 / (0.1 * m2 * PartState(PART_DIAM,iPart2))
#endif /*USE_PARTROT*/

! for LSERK algorithm
! PDM%IsNewPart(iPart1) = .TRUE.
! PDM%IsNewPart(iPart2) = .TRUE.

! Update the particle's position at time tStage+dtloc
PartState(PART_POSV,iPart1) = LastPartPos(PART_POSV,iPart1) + dtColl * PartState(PART_VELV,iPart1)
PartState(PART_POSV,iPart2) = LastPartPos(PART_POSV,iPart2) + dtColl * PartState(PART_VELV,iPart2)

! PartState(PART_POSV,iPart1) = LastPartPos(PART_POSV,iPart1) - dtColl * PartState(PART_VELV,iPart1)
! PartState(PART_POSV,iPart2) = LastPartPos(PART_POSV,iPart2) - dtColl * PartState(PART_VELV,iPart2)

! both particles can no longer be considered for collision
! PartCollisionProc%canPartCollide(iPart1Track) = .FALSE.
! PartCollisionProc%canPartCollide(iPart2Track) = .FALSE.

! de-apply periodic vectors (if there are some)
! IF (GEO%nPeriodicVectors.GT.0) THEN
!   iPartPeriodicShiftLoc = PartCollisionProc%iPartPeriodicShift(iCollision)
!   DO jPeriodic=1,GEO%nPeriodicVectors
!     PartState(PART_POSV,iPart2) = PartState(PART_POSV,iPart2) + REAL(MOD(iPartPeriodicShiftLoc, 3) - 1) * GEO%PeriodicVectors(1:3,jPeriodic)
!     iPartPeriodicShiftLoc = iPartPeriodicShiftLoc / 3
!   END DO
! END IF

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
USE MOD_Particle_Localization    ,ONLY: PARTHASMOVED
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
LOGICAL                       :: isHit
LOGICAL                       :: isCriticalParallelInFace
!-----------------------------------------------------------------------------------------------------------------------------------

dtBC = HUGE(1.)

! Calculate particle trajectory
PartTrajectory       = PartData_Shared(PART_VELV,PartID) * dtLoc
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))

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
      CALL ComputePlanarRectIntersection(   isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                           ,xi,eta,PartID,flip,SideID,isCriticalParallelInFace)
    CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(     isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                           ,xi,eta,PartID,flip,SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection( isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                           ,xi,eta,PartID,flip,SideID,isCriticalParallelInFace)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(       isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                           ,xi,eta,PartID,flip,SideID,isCriticalParallelInFace)
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
LOGICAL                       :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------

dtBC = HUGE(1.)

! check if element is a BC element
IF (nLocSides.EQ.0) RETURN

! Calculate particle trajectory
PartTrajectory       = PartData_Shared(PART_VELV,PartID) * dtLoc
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))

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
      CALL ComputePlanarRectIntersection(  isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                        ,  xi,eta,PartID,flip,SideID)
    CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(    isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                      ,    xi,eta,PartID,flip,SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                          ,xi,eta,PartID,flip,SideID)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(      isHit,PartTrajectory,lengthPartTrajectory,locAlpha(iLocSide) &
                                    ,      xi,eta,PartID,flip,SideID)
  END SELECT
END DO ! iLocSide = firstSide,LastSide

IF (MAXVAL(locAlpha).GT.0) dtBC = MINVAL(locAlpha,locAlpha.GT.0)

END SUBROUTINE ComputeParticleRefmappingIntersection
#endif /*PARTICLES_COUPLING == 4*/

END MODULE MOD_Particle_Collision_Method
