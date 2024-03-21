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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_Tracing
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ParticleTracing
  MODULE PROCEDURE ParticleTracing
END INTERFACE

PUBLIC::ParticleTracing

TYPE :: tIntersectLink
  REAL                          :: alpha  = HUGE(1.)
  REAL                          :: alpha2 = HUGE(1.)
  REAL                          :: xi     = -1
  REAL                          :: eta    = -1
  INTEGER                       :: Side   = 0
  INTEGER                       :: IntersectCase = 0
  TYPE(tIntersectLink), POINTER :: prev => null()
  TYPE(tIntersectLink), POINTER :: next => null()
END TYPE tIntersectLink
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleTracing()
!===================================================================================================================================
!> Routine for tracking of moving particles using polynomial description of sides.
!> Routine calculates intersection and boundary interaction for TrackingMethod.EQ.TRACING
!> Time is analyzed for LoadBalancing purposes for each element independently because elements with e.g. surface are more costly
!> ---------------------------------------------------------------------------------------------------------------------------------
!> - Loop over all particles, which are in own proc --> PDM%ParticleInside(1:PDM%ParticleVecLength)
!> -- 1. Initialize particle path and tracking info
!> -- 2. Track particle vector up to final particle position
!> -- 3. special check if some double check has to be performed (only necessary for bilinear sides)
!> -- 4. Check if particle intersected a side and also which side
!>         For each side only one intersection is chosen, but particle might intersect more than one side. Assign pointer list
!> -- 5. Loop over all intersections in pointer list and check intersection type: inner side, BC
!>       and calculate interaction
!> -- 6. Update particle position and decide if double check might be necessary
!> -- 7. Correct intersection list if double check will be performed and leave loop to do double check
!> -- 8. Reset intersection list if no double check is performed
!> -- 9. If tolerance was marked, check if particle is inside of proc volume and try to find it in case it was lost
!> ---------------------------------------------------------------------------------------------------------------------------------
!> - DoubleCheck:
!> -- If a tracked particle hits a bilinear side but the PartTrajectory points inside of the element,
!>    then the second alpha for this side might have been the actual intersection, which has been dropped in intersection routine.
!> -- Consequently, alpha for doublecheck side is saved (moved to the last position in intersectionlist)
!>    and neglected during the second check of for the appropriate sideID.
!> -- This occurs after surfaceflux, reflection, or for periodic particles moving almost in tangential direction to bilinear side.
!> -- The DoubleCheck replaces the need of tolerances
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Eval_xyz                    ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars                   ,ONLY: SideInfo_Shared
USE MOD_Particle_Localization       ,ONLY: LocateParticleInElement,PARTHASMOVED
USE MOD_Particle_Localization       ,ONLY: PartInElemCheck
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Mesh_Tools         ,ONLY: GetGlobalElemID,GetCNElemID,GetCNSideID,GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemRadiusNGeo
USE MOD_Particle_Operations         ,ONLY: RemoveParticle
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Tracking_vars,      ONLY: ntracks,MeasureTrackTime, CountNbOfLostParts , NbrOfLostParticles
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Utils                       ,ONLY: InsertionSort
#if CODE_ANALYZE
USE MOD_Particle_Globals            ,ONLY: ALMOSTEQUAL
USE MOD_Particle_Intersection       ,ONLY: OutputTrajectory
USE MOD_Particle_Tracking_Vars      ,ONLY: PartOut,MPIRankOut
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO,ElemBaryNGeo
USE MOD_TimeDisc_Vars               ,ONLY: currentStage
#endif /*CODE_ANALYZE*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime,LBElemPauseTime,LBElemSplitTime
USE MOD_Mesh_Vars                   ,ONLY: nElems,offsetElem
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart
INTEGER                       :: ElemID,CNElemID,OldElemID,firstElem
INTEGER                       :: ilocSide,SideID,CNSideID,flip
LOGICAL                       :: dolocSide(1:6)
LOGICAL                       :: PartisDone,foundHit,markTol,crossedBC,SwitchedElement,isCriticalParallelInFace
REAL                          :: localpha,xi,eta
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
LOGICAL                       :: moveList, PartDoubleCheck
REAL                          :: alphaDoneRel, oldLengthPartTrajectory
! intersection info list
TYPE(tIntersectLink),POINTER  :: firstIntersect   => NULL()
TYPE(tIntersectLink),POINTER  :: lastIntersect    => NULL()
TYPE(tIntersectLink),POINTER  :: currentIntersect => NULL()
TYPE(tIntersectLink),POINTER  :: tmp => NULL()
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
#if CODE_ANALYZE
REAL                          :: refpos(1:3)
INTEGER                       :: nIntersections
#endif /*CODE_ANALYZE*/
!===================================================================================================================================

! initialize the first and last pointer in intersection info
IF (.NOT.ASSOCIATED(firstIntersect)) THEN
  ALLOCATE(firstIntersect)
  IF (.NOT.ASSOCIATED(firstIntersect%next)) ALLOCATE(firstIntersect%next)
  lastIntersect      => firstIntersect%next
  lastIntersect%prev => firstIntersect
END IF

DO iPart = 1,PDM%ParticleVecLength
  PartDoubleCheck = .FALSE.

  ! Check if the PartID belongs to a particle that needs to be tracked
  IF (PDM%ParticleInside(iPart)) THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
    ! check if particle is inside domain bounding box
    IF(GEO%nPeriodicVectors.EQ.0)THEN
      IF(   (LastPartPos(1,iPart).GT.GEO%xmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,iPart),GEO%xmaxglob) &
        .OR.(LastPartPos(1,iPart).LT.GEO%xminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,iPart),GEO%xminglob) &
        .OR.(LastPartPos(2,iPart).GT.GEO%ymaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,iPart),GEO%ymaxglob) &
        .OR.(LastPartPos(2,iPart).LT.GEO%yminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,iPart),GEO%yminglob) &
        .OR.(LastPartPos(3,iPart).GT.GEO%zmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,iPart),GEO%zmaxglob) &
        .OR.(LastPartPos(3,iPart).LT.GEO%zminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,iPart),GEO%zminglob) ) THEN
        IPWRITE(UNIT_stdOut,'(I0,A18,L)')                            ' ParticleInside ', PDM%ParticleInside(iPart)
        IPWRITE(UNIT_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(iPart)
        IPWRITE(UNIT_stdOut,'(I0,A18,1X,A18,1X,A18)')                '    min ', ' value ', ' max '
        IPWRITE(UNIT_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' x', GEO%xminglob, LastPartPos(1,iPart), GEO%xmaxglob
        IPWRITE(UNIT_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' y', GEO%yminglob, LastPartPos(2,iPart), GEO%ymaxglob
        IPWRITE(UNIT_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' z', GEO%zminglob, LastPartPos(3,iPart), GEO%zmaxglob
        CALL Abort(__STAMP__,' LastPartPos outside of mesh. iPart=, currentStage',iPart,REAL(currentStage))
      END IF
    END IF
    ! caution: reuse of variable, foundHit=TRUE == inside
    ElemID   = PEM%lastElement(iPart)
    CNElemID = GetCNElemID(ElemID)
    CALL GetPositionInRefElem(LastPartPos(1:3,iPart),RefPos,ElemID)
    IF (MAXVAL(ABS(RefPos)).LE.1.0+1e-4) THEN
      foundHit = .TRUE.
    ELSE
      foundHit = .FALSE.
    END IF
    IF(.NOT.foundHit)THEN  ! particle not inside
      IPWRITE(UNIT_stdOut,'(I0,A)')                 ' PartPos not inside of element! '
      IPWRITE(UNIT_stdOut,'(I0,A,I0)')              ' PartID             ', iPart
      IPWRITE(UNIT_stdOut,'(I0,A,I0)')              ' global ElemID      ', ElemID
      IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' ElemBaryNGeo:      ', ElemBaryNGeo(1:3,CNElemID)
      IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' LastPartPos:       ', LastPartPos(1:3,iPart)
      IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' PartPos:           ', PartState(1:3,iPart)
      IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' PartRefPos:        ', RefPos(1:3)
      IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' Velocity:          ', PartState(4:6,iPart)
      IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' PartTrajectory:    ', PartState(1:3,iPart) - LastPartPos(1:3,iPart)
      IPWRITE(UNIT_stdOut,'(I0,A,ES25.14E3)')       ' lengthPT:          ', SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
      CALL Abort(__STAMP__,'ERROR: Lastpartpos in wrong element. PartID:',iPart)
    END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

! -- 1. Initialize particle path and tracking info
    IF (MeasureTrackTime) nTracks = nTracks+1
    PartisDone = .FALSE.
    ElemID     = PEM%lastElement(iPart)
    CNElemID   = GetCNElemID(ElemID)

    ! Calculate particle trajectory
    PartTrajectory          = PartState(1:3,iPart) - LastPartPos(1:3,iPart)
    lengthPartTrajectory    = SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
    oldLengthPartTrajectory = LengthPartTrajectory
    alphaDoneRel            = 0.

    ! Check if the particle moved at all. If not, tracking is done
    IF(.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(CNElemID)) .OR. LengthPartTrajectory.EQ.0) THEN
      PEM%Element(iPart) = ElemID
      PartisDone         = .TRUE.
      CYCLE
    ELSE
      PartTrajectory = PartTrajectory/lengthPartTrajectory
    END IF

#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
    IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN ; IF (iPart.EQ.PARTOUT) THEN
      WRITE(UNIT_stdOut,'(A32)')' ---------------------------------------------------------------'
      WRITE(UNIT_stdOut,'(A)')  '     | Output of Particle information '
      CALL OutputTrajectory(iPart,PartState(1:3,iPart),PartTrajectory,lengthPartTrajectory,LastPartPos(:,iPart))
      WRITE(UNIT_stdOut,'(A,I0)') '     | global ElemID       ', PEM%LastElement(iPart)
    END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

    ! track particle vector until the final particle position is achieved
    dolocSide = .TRUE.
    firstElem = ElemID
    markTol   = .FALSE.

! -- 2. Track particle vector up to the final particle position
    DO WHILE (.NOT.PartisDone)
      ! do not reset markTol after first intersection of for doublecheck.
      ! This prevents particles to get lost unnoticed in case any intersection has marked tolerance.
      ! markTol =.FALSE.

      IF (PartDoubleCheck) THEN
! -- 3. special check if some double check has to be performed (only necessary for bilinear sides)
#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN; IF (iPart.EQ.PARTOUT) THEN
          WRITE(UNIT_stdOut,'(A)')    '     | Calculation of double check: '
        END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        currentIntersect => lastIntersect%prev
        IF (currentIntersect%IntersectCase.EQ.1) THEN
          iLocSide = currentIntersect%Side
          SideID   = GetGlobalNonUniqueSideID(ElemID,iLocSide)
          CNSideID = GetCNSideID(SideID)
          flip     = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

          CALL ComputeBiLinearIntersection(foundHit               &
                                          ,PartTrajectory         &
                                          ,lengthPartTrajectory   &
                                          ,PartState(  1:3,iPart) &
                                          ,LastPartPos(1:3,iPart) &
                                          ,locAlpha               &
                                          ,xi                     &
                                          ,eta                    &
                                          ,iPart                  &
                                          ,flip                   &
                                          ,SideID                 &
                                          ,alpha2=currentIntersect%alpha)
          currentIntersect%alpha         = HUGE(1.)
          currentIntersect%IntersectCase = 0
          IF(foundHit) THEN
            CALL AssignListPosition(currentIntersect,locAlpha,iLocSide,1,xi_IN=xi,eta_IN=eta)
            IF((ABS(xi).GE.0.99).OR.(ABS(eta).GE.0.99)) markTol=.TRUE.
          END IF
        END IF
#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN; IF (iPart.EQ.PARTOUT) THEN
          WRITE(UNIT_stdOut,'(30("-"))')
          WRITE(UNIT_stdOut,'(A)')             '     | Output after compute intersection (tracing double check): '
          IF (currentIntersect%IntersectCase.EQ.1) THEN
            WRITE(UNIT_stdOut,'(2(A,I0),A,L)') '     | SideType: ',SideType(CNSideID),' | SideID: ',SideID,' | Hit: ',foundHit
            WRITE(UNIT_stdOut,'(A,2(1X,G0))')  '     | Intersection xi/eta: ',xi,eta
            WRITE(UNIT_stdOut,'((A,G0))')      '     | RelAlpha: ',locAlpha/lengthpartTrajectory
          END IF
          WRITE(UNIT_stdOut,'(2(A,G0))')       '     | Alpha: ',locAlpha,' | LengthPartTrajectory: ', lengthPartTrajectory
        END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        ! if double check found no intersection reset entry in list and adjust last entry pointer
        IF (.NOT.foundHit) THEN
          currentIntersect%alpha = HUGE(1.)
          currentIntersect%intersectCase = 0
          IF (ASSOCIATED(currentIntersect%prev) .AND. .NOT.ASSOCIATED(currentIntersect%prev,firstIntersect)) THEN
            lastIntersect => currentIntersect
            lastIntersect%prev => currentIntersect%prev
          END IF
        END IF

      ! NOT PartDoubleCheck
      ELSE
! -- 4. Check if particle intersected a side and also which side
!       For each side only one intersection is chosen, but particle might intersect more than one side. Assign pointer list
#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN; IF(iPart.EQ.PARTOUT) THEN
          WRITE(UNIT_stdOut,'(110("="))')
          WRITE(UNIT_stdOut,'(A)')    '     | Calculation of particle intersections: '
        END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        DO ilocSide = 1,6
          locAlpha = -1.
          IF(.NOT.dolocSide(ilocSide)) CYCLE
          SideID   = GetGlobalNonUniqueSideID(ElemID,iLocSide)
          CNSideID = GetCNSideID(SideID)
          ! If the side is positive, then the element has the actual side
          ! and neighbour element has the negative one which has to be flipped

          ! BezierControlPoints are now built in cell local system. Hence, sides have always the flip from the shared SideInfo
          flip = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

          isCriticalParallelInFace = .FALSE.

          SELECT CASE(SideType(CNSideID))
            CASE(PLANAR_RECT)
              CALL ComputePlanarRectIntersection(  foundHit,PartTrajectory,lengthPartTrajectory,                     LastPartPos(:,iPart),locAlpha &
                                                ,  xi,eta,iPart,flip,SideID,isCriticalParallelInFace)
            CASE(BILINEAR,PLANAR_NONRECT)
              CALL ComputeBiLinearIntersection(    foundHit,PartTrajectory,lengthPartTrajectory,PartState(1:3,iPart),LastPartPos(:,iPart),locAlpha &
                                              ,    xi,eta,iPart,flip,SideID)
            CASE(PLANAR_CURVED)
              CALL ComputePlanarCurvedIntersection(foundHit,PartTrajectory,lengthPartTrajectory,                     LastPartPos(:,iPart),locAlpha &
                                                  ,xi,eta,iPart,flip,SideID,isCriticalParallelInFace)
            CASE(CURVED)
              CALL ComputeCurvedIntersection(      foundHit,PartTrajectory,lengthPartTrajectory,PartState(1:3,iPart),LastPartPos(:,iPart),locAlpha &
                                            ,      xi,eta,iPart,flip,SideID,isCriticalParallelInFace)
            CASE DEFAULT
              CALL Abort(__STAMP__,' Missing required side-data. Please increase halo region. ',SideID)
          END SELECT

#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
          IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN; IF (iPart.EQ.PARTOUT) THEN
            WRITE(UNIT_stdOut,'(30("-"))')
            WRITE(UNIT_stdOut,'(A)')           '     | Output after compute intersection (particle tracing): '
            WRITE(UNIT_stdOut,'(2(A,I0),A,L)') '     | SideType: ',SideType(CNSideID),' | SideID: ',SideID,' | Hit: ',foundHit
            WRITE(UNIT_stdOut,'(2(A,G0))')     '     | Alpha: ',locAlpha,' | LengthPartTrajectory: ', lengthPartTrajectory
            WRITE(UNIT_stdOut,'((A,G0))')      '     | RelAlpha: ',locAlpha/lengthpartTrajectory
            WRITE(UNIT_stdOut,'(A,2(1X,G0))')  '     | Intersection xi/eta: ',xi,eta
          END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
          ! Particle detected inside of face and PartTrajectory parallel to face
          IF (isCriticalParallelInFace) THEN
            IPWRITE(UNIT_stdOut,'(I0,A)')    ' Warning: Particle located inside of face and moves parallel to side. Undefined position.'
            IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Removing particle with id: ',iPart
            PartIsDone = .TRUE.
            CALL RemoveParticle(iPart)
            IF(CountNbOfLostParts) NbrOfLostParticles = NbrOfLostParticles + 1
            EXIT
          END IF
          IF (foundHit) THEN
            currentIntersect   => lastIntersect
            CALL AssignListPosition(currentIntersect,locAlpha,iLocSide,1,xi_IN=xi,eta_IN=eta)
            currentIntersect   => lastIntersect
            lastIntersect      => currentIntersect%next
            lastIntersect%prev => currentIntersect
            IF ((ABS(xi).GE.0.99).OR.(ABS(eta).GE.0.99)) markTol=.TRUE.
          END IF
        END DO ! ilocSide
      END IF

! -- 5. Loop over all intersections in pointer list and check intersection type: inner side, BC and calculate interaction
#if CODE_ANALYZE
      nIntersections = 0
#endif /*CODE_ANALYZE*/
      currentIntersect => firstIntersect
      DO WHILE(ASSOCIATED(currentIntersect))
        SwitchedElement = .FALSE.
        crossedBC       = .FALSE.
#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
        nIntersections = nIntersections+1
        IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN; IF (iPart.EQ.PARTOUT) THEN
          WRITE(UNIT_stdOut,'(45(":"))')
          WRITE(UNIT_stdOut,'(A,I0)')  '     -> Check intersection: ', nIntersections
          WRITE(UNIT_stdOut,'(A,I0)')  '     -> Case: '   ,currentIntersect%IntersectCase
          WRITE(UNIT_stdOut,'(A,G0)')  '     -> alpha: '  ,currentIntersect%alpha
          WRITE(UNIT_stdOut,'(A,I0)')  '     -> locSide: ',currentIntersect%Side
          IF (currentIntersect%IntersectCase.EQ.1) THEN
            WRITE(UNIT_stdOut,'(A,I0)') '     -> SideID: ',GetGlobalNonUniqueSideID(ElemID,currentIntersect%Side)
          END IF
        END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
        OldElemID = ElemID
        IF (currentIntersect%IntersectCase.EQ.0) THEN
          ! no intersection
          PEM%Element(iPart ) = ElemID
          PartisDone          = .TRUE.
        ELSE
          SELECT CASE(currentIntersect%IntersectCase)
          !------------------------------------
          CASE(1) ! intersection with cell side
          !------------------------------------
            SideID   = GetGlobalNonUniqueSideID(ElemID,currentIntersect%Side)
            CNSideID = GetCNSideID(SideID)
            flip     = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

            ! missing!!! : mapping from GlobalNonUnique to CNtotalsides
            CALL SelectInterSectionType( PartIsDone                   &
                                       , crossedBC                    &
                                       , doLocSide                    &
                                       , flip                         &
                                       , PartTrajectory               &
                                       , lengthPartTrajectory         &
                                       , currentIntersect%xi          &
                                       , currentIntersect%eta         &
                                       , currentIntersect%alpha       &
                                       , iPart                        &
                                       , SideID                       &
                                       , SideType(CNSideID)           &
                                       , ElemID)

            IF (ElemID.NE.OldElemID) THEN
              IF (.NOT.crossedBC) SwitchedElement=.TRUE.
            END IF
          END SELECT
#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
          IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN; IF (iPart.EQ.PARTOUT) THEN
            IF (crossedBC) THEN
              SELECT CASE(currentIntersect%IntersectCase)
                CASE(1) ! intersection with cell side
                  WRITE(UNIT_stdOut,'(A,L)') '     -> BC was intersected on a side'
              END SELECT
            END IF
          END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

        ! TODO

! -- 6. Update particle position and decide if double check might be necessary
! check what happened with particle (crossedBC or switched element) and set partisdone or double check
#if USE_LOADBALANCE
          IF (OldElemID.GE.offsetElem+1 .AND.OldElemID.LE.offsetElem+nElems) &
            CALL LBElemSplitTime(OldElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
          IF (crossedBC) firstElem = ElemID

          IF (crossedBC .OR. SwitchedElement) THEN
            IF (PartDoubleCheck) THEN
              PartDoubleCheck=.FALSE.
            END IF
            EXIT
          ELSE ! (.NOT.crossedBC .AND. .NOT.SwitchedElement)
            IF (.NOT.PartDoubleCheck) THEN
              PartDoubleCheck = .TRUE.
              PartIsDone      = .FALSE.
            END IF
          END IF
        END IF ! IntersectCase.EQ.0

! -- 7. Correct intersection list if double check will be performed and leave loop to do double check
        ! move current list entry to the end and the total list to the front. exit and check if the last is the correct intersection
        IF (.NOT.crossedBC .AND. .NOT.SwitchedElement .AND. .NOT.PartIsDone .AND. PartDoubleCheck) THEN
          moveList = .FALSE.
          SELECT CASE(currentIntersect%intersectCase)
            CASE(1)
              SideID   = GetGlobalNonUniqueSideID(OldElemID,currentIntersect%Side)
              CNSideID = GetCNSideID(SideID)

              SELECT CASE(SideType(CNSideID))
                CASE(BILINEAR,PLANAR_NONRECT)
                  moveList=.TRUE.
             END SELECT
          END SELECT
          IF (moveList) THEN
            lastIntersect%alpha  = currentIntersect%alpha
            lastIntersect%alpha2 = currentIntersect%alpha2
            lastIntersect%xi     = currentIntersect%xi
            lastIntersect%eta    = currentIntersect%eta
            lastIntersect%Side   = currentIntersect%Side
            lastIntersect%intersectCase = currentIntersect%intersectCase
            tmp => firstIntersect
            DO WHILE (.NOT.ASSOCIATED(tmp,lastIntersect))
              tmp%alpha  = tmp%next%alpha
              tmp%alpha2 = tmp%next%alpha2
              tmp%xi     = tmp%next%xi
              tmp%eta    = tmp%next%eta
              tmp%Side   = tmp%next%Side
              tmp%intersectCase = tmp%next%intersectCase
              tmp => tmp%next
            END DO
            EXIT
          END IF
        END IF

        ! leave loop because particle is found to remain in element (none of the found intersections is valid)
        currentIntersect => currentIntersect%next
        IF (ASSOCIATED(currentIntersect,LastIntersect)) THEN
          PartDoubleCheck = .FALSE.
          PartIsDone      = .TRUE.
          EXIT
        END IF
      END DO ! ASSOCIATED(currentIntersect)
#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank) THEN; IF (iPart.EQ.PARTOUT) THEN
        WRITE(UNIT_stdOut,'(128("="))')
        WRITE(UNIT_stdOut,'(A)')            '     | Output of tracking information after the check of number of intersections: '
        WRITE(UNIT_stdOut,'(4(A,L))')       '     | crossed Side: ',crossedBC,' switched Element: ',SwitchedElement,&
          ' Particle tracking done: ',PartisDone,' Particle is double checked: ',PartDoubleCheck
        IF(SwitchedElement) THEN
          WRITE(UNIT_stdOut,'(A,I0,A,I0)')  '     | First_ElemID: ',PEM%LastElement(iPart),' | new Element: ',ElemID
          WRITE(UNIT_stdOut,'(A,I0)')       '     | first global ElemID     ', PEM%LastElement(iPart)
          WRITE(UNIT_stdOut,'(A,I0)')       '     | new global ElemID       ', ElemID
        END IF
        IF( crossedBC) THEN
          WRITE(UNIT_stdOut,'(A,3(1X,G0))') '     | Last    PartPos:        ',lastPartPos(1:3,iPart)
          WRITE(UNIT_stdOut,'(A,3(1X,G0))') '     | Current PartPos:        ',PartState(1:3,iPart)
          WRITE(UNIT_stdOut,'(A,3(1X,G0))') '     | PartTrajectory:         ',PartTrajectory(1:3)
          WRITE(UNIT_stdOut,'(A,(G0))')     '     | Length PartTrajectory:  ',lengthPartTrajectory
        END IF
        WRITE(UNIT_stdOut,'(128("="))')
      END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
! -- 8. Reset interscetion list if no double check is performed
      ! reset intersection list because no intersections where found or no double check is performed or no interacions occurred
      currentIntersect=>firstIntersect
      IF (currentIntersect%intersectCase.GT.0 .AND. .NOT.PartDoubleCheck)THEN
        DO WHILE (ASSOCIATED(currentIntersect))
          currentIntersect%alpha = HUGE(1.)
          currentIntersect%intersectCase = 0
          IF(ASSOCIATED(currentIntersect,lastIntersect)) THEN
            lastIntersect => firstIntersect%next
            lastIntersect%prev => firstIntersect
            EXIT
        END IF
          currentIntersect => currentIntersect%next
        END DO
      END IF
    END DO ! PartisDone=.FALSE.

#if USE_LOADBALANCE
    IF (PEM%Element(iPart).GE.offsetElem+1 .AND. PEM%Element(iPart).LE.offsetElem+nElems) &
      CALL LBElemPauseTime(PEM%Element(iPart)-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF ! Part inside
END DO ! iPart

#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
! check if particle is still inside of bounding box of domain and in element
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    IF( (PartState(1,iPart).GT.GEO%xmaxglob) &
    .OR.(PartState(1,iPart).LT.GEO%xminglob) &
    .OR.(PartState(2,iPart).GT.GEO%ymaxglob) &
    .OR.(PartState(2,iPart).LT.GEO%yminglob) &
    .OR.(PartState(3,iPart).GT.GEO%zmaxglob) &
    .OR.(PartState(3,iPart).LT.GEO%zminglob) ) THEN
      IPWRITE(UNIT_stdOut,'(I0,A18,L)')                                 ' PDM%IsNewPart ', PDM%IsNewPart(iPart)
      IPWRITE(UNIT_stdOut,'(I0,A18,3(1X,E27.16))')                      ' LastPosition   ', LastPartPos(1:3,iPart)
      IPWRITE(UNIT_stdOut,'(I0,A18,3(1X,E27.16))')                      ' Velocity       ', PartState  (4:6,iPart)
      IPWRITE(UNIT_stdOut,'(I0,A18,1X,A18,1X,A18)')                     '    min ', ' value ', ' max '
      IPWRITE(UNIT_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)')      ' x', GEO%xminglob, PartState(1,iPart), GEO%xmaxglob
      IPWRITE(UNIT_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)')      ' y', GEO%yminglob, PartState(2,iPart), GEO%ymaxglob
      IPWRITE(UNIT_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)')      ' z', GEO%zminglob, PartState(3,iPart), GEO%zmaxglob
      CALL Abort(__STAMP__,' PartPos outside of mesh AFTER tracking. iPart= ,currentStage= ',iPart,REAL(currentStage))
    END IF

    ! caution: reuse of variable, foundHit=TRUE == inside
    ElemID=PEM%Element(iPart)
    CALL GetPositionInRefElem(PartState(1:3,iPart),RefPos,ElemID)
    IF (MAXVAL(ABS(RefPos)).LE.1.0+1e-4) THEN
      foundHit = .TRUE.
    ELSE
      foundHit = .FALSE.
    END IF
    IF (.NOT.foundHit) THEN  ! particle not inside
     IPWRITE(UNIT_stdOut,'(I0,A)') ' PartPos not inside of element! '
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')  ' PartID         ', iPart
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')  ' global ElemID  ', ElemID
     IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' ElemBaryNGeo:      ', ElemBaryNGeo(1:3,ElemID)
     IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' LastPartPos:       ', LastPartPos(1:3,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' PartPos:           ', PartState(1:3,iPart)
     IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' PartRefPos:        ', RefPos(1:3)
     IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.14E3))') ' PartTrajectory:    ', PartTrajectory
     IPWRITE(UNIT_stdOut,'(I0,A,ES25.14E3)')      ' lengthPT:           ', lengthPartTrajectory
     CALL Abort(__STAMP__,'iPart=. ',iPart)
     END IF
  END IF ! Part inside
END  DO ! iPart=1,PDM%ParticleVecLength
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

END SUBROUTINE ParticleTracing


SUBROUTINE AssignListPosition(inLink,alpha_IN,sideID_IN,IntersectCase_IN,xi_IN,eta_IN,alpha2_IN)
!===================================================================================================================================
!> Checks the given intersection linked list starting from the last to the first entry and compares with a given alpha.
!> adds the given alpha at the correct position extending the list if necessary.
!> exits the search if position was assigned.
!> -----------------------
!> first list entry is the smallest found alpha and last is the largest.
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tIntersectLink),POINTER,INTENT(INOUT) :: inLink
REAL,INTENT(IN)              :: alpha_IN
INTEGER,INTENT(IN)           :: sideID_IN
INTEGER,INTENT(IN)           :: IntersectCase_IN
REAL,INTENT(IN),OPTIONAL     :: xi_IN
REAL,INTENT(IN),OPTIONAL     :: eta_IN
REAL,INTENT(IN),OPTIONAL     :: alpha2_IN
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tIntersectLink),POINTER :: tmpLink
!===================================================================================================================================
! start from last intersection entry and place current intersection in correct entry position
tmpLink => inLink
DO WHILE(ASSOCIATED(tmpLink))
  IF (alpha_IN.LE.tmpLink%alpha) THEN
    IF (.NOT. ASSOCIATED(inLink%next)) THEN
      ALLOCATE(inLink%next)
      inLink%next%prev => inLink
    END IF
    tmpLink%next%alpha = tmpLink%alpha
    tmpLink%next%alpha2 = tmpLink%alpha2
    tmpLink%next%xi = tmpLink%xi
    tmpLink%next%eta = tmpLink%eta
    tmpLink%next%Side = tmpLink%Side
    tmpLink%next%intersectCase = tmpLink%intersectCase
    IF (ASSOCIATED(tmpLink%prev)) THEN
      IF (alpha_IN.GT.tmpLink%prev%alpha) THEN
        ! assign new values
        tmpLink%alpha = alpha_IN
        tmpLink%Side = sideID_IN
        tmpLink%intersectCase = IntersectCase_IN
        IF (PRESENT(xi_IN)) tmpLink%xi = xi_IN
        IF (PRESENT(eta_IN)) tmpLink%eta = eta_IN
        IF (PRESENT(alpha2_IN)) tmpLink%alpha2 = alpha2_IN
        EXIT
      END IF
    ELSE
      ! assign new values
      tmpLink%alpha = alpha_IN
      tmpLink%Side = sideID_IN
      tmpLink%intersectCase = IntersectCase_IN
      IF (PRESENT(xi_IN)) tmpLink%xi = xi_IN
      IF (PRESENT(eta_IN)) tmpLink%eta = eta_IN
      IF (PRESENT(alpha2_IN)) tmpLink%alpha2 = alpha2_IN
      EXIT
    END IF
  END IF
  tmpLink => tmpLink%prev
END DO

END SUBROUTINE AssignListPosition


SUBROUTINE SelectInterSectionType(PartIsDone,crossedBC,doLocSide,flip,PartTrajectory,lengthPartTrajectory &
                               ,xi,eta,alpha,PartID,SideID,hitSideType,ElemID)
!===================================================================================================================================
!> Use only for TrackingMethod = TRACING or TracingElement with RefMapping
!> Checks which type of interaction (BC,Periodic,innerSide) has to be applied for the face on the traced particle path
!> - If face is BC-side BoundaryInteraction routine is called
!>   - for triatracking the intersection location of partice trajectory with face is calculated first
!> - If face is innerside switch to respective Element
!>   - for tracing check if path for considered intersection point into current element.
!>     Can happen for particles inserted during surface flux at bilinear faces (double checks filters those intersections after)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars                   ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarNonRectIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Mesh_Tools         ,ONLY: GetGlobalNonUniqueSideID,GetCNSideID
USE MOD_Particle_Surfaces           ,ONLY: CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars,      ONLY: SideNormVec
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Tracking_Vars      ,ONLY: TrackInfo
USE MOD_Particle_Vars,               ONLY: PDM,PartState,LastPartPos
#if CODE_ANALYZE
USE MOD_Mesh_Vars                   ,ONLY: NGeo
USE MOD_Particle_Localization       ,ONLY: SinglePointToElement
USE MOD_Particle_Surfaces_Vars      ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemBaryNGeo
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: PartID                   !< Index of Considered Particle
INTEGER,INTENT(IN)                :: SideID                   !< SideID particle intersected with
INTEGER,INTENT(IN)                :: hitSideType              !< type of SideID (planar,bilinear,...)
INTEGER,INTENT(IN)                :: flip                     !< flip of SideID
REAL,INTENT(INOUT)                :: Xi                       !<
REAL,INTENT(INOUT)                :: Eta                      !<
REAL,INTENT(INOUT)                :: Alpha                    !< portion of PartTrajectory until hit with face
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(INOUT)             :: PartIsDone               !< Flag indicating if tracking of PartID is finished
LOGICAL,INTENT(OUT)               :: crossedBC                !< Flag indicating if BC has been hit
LOGICAL,INTENT(INOUT)             :: DoLocSide(1:6)           !<
INTEGER,INTENT(INOUT)             :: ElemID                   !< Element ID particle is currently in
REAL,INTENT(INOUT),DIMENSION(1:3) :: PartTrajectory           !< normalized particle trajectory (x,y,z)
REAL,INTENT(INOUT)                :: lengthPartTrajectory     !< length of particle trajectory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: isHit
INTEGER                           :: iMortar,nMortarElems
INTEGER                           :: CNSideID,NbElemID,NbSideID,NbCNSideID
INTEGER                           :: iLocalSide
!INTEGER                           :: locFlip
REAL                              :: locAlpha,locXi,locEta
REAL                              :: n_loc(3)
#if CODE_ANALYZE
REAL                              :: v1(3),v2(3)
#endif /* CODE_ANALYZE */
!===================================================================================================================================
! Side is a boundary side
IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
  CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,flip,ElemID,crossedBC)
  TrackInfo%CurrElem = ElemID
  IF(.NOT.PDM%ParticleInside(PartID)) PartisDone = .TRUE.
  dolocSide = .TRUE.

! Side is not a boundary side
ELSE
  ! get side normal vector. Must be calculated with current point if the side is bilinear or curved
  SELECT CASE(hitSideType)
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      CNSideID = GetCNSideID(SideID)
      n_loc    = SideNormVec(1:3,CNSideID)
      ! BezierControlPoints are now built in cell local system. Hence, side always have the flip from the shared SideInfo
      IF (flip.NE.0) n_loc = -n_loc

    ! bilinear sides have no valid side normal vector. Calculate with current intersection point and do not flip
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    ! curved sides have no valid side normal vector. Calculate with current intersection point and do not flip
      CASE(CURVED)
        CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT

#if CODE_ANALYZE
  ! check if normal vector points outwards
  v1 = 0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)  &
           + BezierControlPoints3D(:,NGeo,0   ,SideID)  &
           + BezierControlPoints3D(:,0   ,NGeo,SideID)  &
           + BezierControlPoints3D(:,NGeo,NGeo,SideID))
  v2 = v1  - ElemBaryNGeo(:,GetCNElemID(ElemID))

  IF (DOT_PRODUCT(v2,n_loc).LT.0) THEN
    IPWRITE(UNIT_stdOut,'(I0,A,I0,A,I0,A,I0)') ' Obtained wrong side orientation from flip. SideID:',SideID,'flip:',flip,'PartID:',PartID
    IPWRITE(UNIT_stdOut,'(I0,A,I0,A,3F12.6)')  ' n_loc (flip)', n_loc,'n_loc (estimated):',v2
    CALL Abort(__STAMP__,'SideID',SideID)
  END IF
#endif /* CODE_ANALYZE */

    IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0) RETURN

  ! update particle element
  ! check if the side is a big mortar side
  NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

#if CODE_ANALYZE
  WRITE(UNIT_stdOut,'(30("-"))')
  WRITE(UNIT_stdOut,'(A,I0,A,I0,A,F12.6,A,I0,A,I0)') ' ElemID:',ElemID,'PartID',PartID,'SideID:',SideID,'Move rel. to Side:',DOT_PRODUCT(n_loc,PartTrajectory),'NbElemID:',NbElemID, 'PartElem (w/o refl.)', SinglePointToElement(PartState(1:3,PartID),doHalo=.TRUE.)
  WRITE(UNIT_stdOut,'(A,3F12.6,A,3F12.6)')           ' PartPos',PartState(1:3,PartID), 'PartVel:',PartState(4:6,PartID)
#endif /* CODE_ANALYZE */

  IF (NbElemID.LT.0) THEN ! Mortar side
    nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

    DO iMortar = 1,nMortarElems
      NbSideID   = SideInfo_Shared(SIDE_NBSIDEID,SideID + iMortar)
      NbCNSideID = GetCNSideID(NbSideID)
      ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
      ! no way to recover it during runtime
      IF (NbSideID.LT.1) CALL ABORT(__STAMP__,'Small mortar side not defined!',SideID + iMortar)

      NbElemID = SideInfo_Shared(SIDE_ELEMID,nbSideID)
      ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
      ! no way to recover it during runtime
      IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)

      ! BezierControlPoints are now built in cell local system. We are checking mortar sides, so everything is reversed
      ! locFlip = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,nbSideID),10),SideInfo_Shared(SIDE_ID,nbSideID).GT.0)
      ! Small mortar sides are always slave sides, hence check with flip = 0

      SELECT CASE(SideType(NbCNSideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarRectIntersection(   isHit,PartTrajectory,lengthPartTrajectory,                      LastPartPos(:,PartID),locAlpha &
                                            ,   locXi,locEta,PartID,0      ,NbSideID)
        CASE(PLANAR_NONRECT)
          CALL ComputePlanarNonRectIntersection(isHit,PartTrajectory,lengthPartTrajectory,                      LastPartPos(:,PartID),locAlpha &
                                               ,locXi,locEta,PartID,0      ,NbSideID)
        CASE(BILINEAR)
          CALL ComputeBiLinearIntersection(     isHit,PartTrajectory,lengthPartTrajectory,PartState(1:3,PartID),LastPartPos(:,PartID),locAlpha &
                                          ,     locXi,locEta,PartID,0      ,NbSideID)
        CASE(PLANAR_CURVED)
          CALL ComputePlanarCurvedIntersection( isHit,PartTrajectory,lengthPartTrajectory,                      LastPartPos(:,PartID),locAlpha &
                                          ,     locXi,locEta,PartID,0      ,NbSideID)
        CASE(CURVED)
          CALL ComputeCurvedIntersection(       isHit,PartTrajectory,lengthPartTrajectory,PartState(1:3,PartID),LastPartPos(:,PartID),locAlpha &
                                        ,       locXi,locEta,PartID,0      ,NbSideID)
      END SELECT

      IF (isHit) THEN
        ElemID = SideInfo_Shared(SIDE_ELEMID,NbSideID)
        TrackInfo%CurrElem = ElemID

        ! flag the side the particle passed through, so it does not have to be checked in the next tracing step
  dolocSide = .TRUE.
        DO iLocalSide = 1,6
          IF (NbSideID.EQ.GetGlobalNonUniqueSideID(NbElemID,iLocalSide)) THEN
            dolocSide(iLocalSide) = .FALSE.
            EXIT
          END IF
        END DO
        RETURN
      END IF
    END DO

    ! passed none of the mortar elements. Keep particle inside current element and warn
    IPWRITE(UNIT_stdOut,'(I0,A,I0)') 'Boundary issue with inner mortar element', ElemID

  ! regular side
  ELSE
    ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
    IF (ElemID.LT.1) &
      CALL Abort(__STAMP__,'ERROR in SelectIntersectionType. No Neighbour Elem found!')

    TrackInfo%CurrElem = ElemID

    ! flag the side the particle passed through, so it does not have to be checked in the next tracing step
    dolocSide = .TRUE.
    DO iLocalSide = 1,6
      IF (ABS(SideInfo_Shared(SIDE_ID,SideID)).EQ.ABS(SideInfo_Shared(SIDE_ID,GetGlobalNonUniqueSideID(ElemID,iLocalSide)))) THEN
        dolocSide(iLocalSide) = .FALSE.
        EXIT
      END IF
    END DO

  END IF
END IF

END SUBROUTINE SelectInterSectionType

END MODULE MOD_Particle_Tracing
