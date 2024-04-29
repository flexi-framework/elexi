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
MODULE MOD_Particle_RefTracking
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ParticleRefTracking
  MODULE PROCEDURE ParticleRefTracking
END INTERFACE

PUBLIC::ParticleRefTracking
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleRefTracking()
!===================================================================================================================================
! Reference Tracking for particle without treatment of each inner faces
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars              ,ONLY: useCurveds,NGeo
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_Particle_Localization  ,ONLY: PartInElemCheck
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides
USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGM_nElems,FIBGM_Element,FIBGM_offsetElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_Tracking_Vars ,ONLY: nTracks,Distance,ListDistance,CartesianPeriodic
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef,LastPartPos
USE MOD_Utils                  ,ONLY: ALMOSTEQUAL
USE MOD_Utils                  ,ONLY: InsertionSort
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_LoadBalance_Vars       ,ONLY: nTracksPerElem
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBElemSplitTime,LBElemPauseTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Counters
INTEGER                           :: iPart
! Elements
INTEGER                           :: ElemID,oldElemID,newElemID,LastElemID
INTEGER                           :: CNElemID
! Background mesh
INTEGER                           :: CellX,CellY,CellZ,iBGMElem,nBGMElems
! Particles
REAL                              :: oldXi(3),newXi(3),LastPos(3)
REAL                              :: vec(3),lengthPartTrajectory0
#if USE_MPI
INTEGER                           :: InElem
CHARACTER(LEN=1)                  :: HaloStat
#endif
! Tracking
INTEGER                           :: TestElem,CNTestElem
LOGICAL                           :: PartisDone,PartIsMoved
#if USE_LOADBALANCE
REAL                              :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

DO iPart = 1,PDM%ParticleVecLength
  ! Ignore particles not inside
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  LastElemID = PEM%lastElement(iPart)
  ElemID     = LastElemID
  CNElemID   = GetCNElemID(ElemID)
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  nTracks    = nTracks+1
  ! sanity check
  PartIsDone = .FALSE.

  ! check if element is a BC element. If yes, handle with Tracing instead of RefMapping
  IF (ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID).GT.0) THEN
    lengthPartTrajectory0 = 0.
    CALL ParticleBCTracking(lengthPartTrajectory0                                                                   &
                           ,ElemID                                                                                  &
                           ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + 1                                           &
                           ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID)    &
                           ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNElemID)                                               &
                           ,iPart                                                                                   &
                           ,PartIsDone                                                                              &
                           ,PartIsMoved                                                                             &
                           ,1)
    ! Particle has left domain by a boundary condition
    IF(PartIsDone) CYCLE

    ! Particle is reflected at a wall or has not encountered any boundary condition
    CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)

    ! Position in reference space smaller unity, particle inside
    IF(MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN
       PEM%Element(iPart) = ElemID
#if USE_LOADBALANCE
       ! Particle is on current proc, assign load to new cell
       IF (ElemID.GT.offsetElem+1 .AND. ElemID.LE.offsetElem+nElems) &
         CALL LBElemPauseTime(ElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
      CYCLE
    END IF

  ! No boundary element, therefore, no boundary interaction possible
  ELSE
    ! Account for periodic displacement
    IF (GEO%nPeriodicVectors.GT.0 .AND. CartesianPeriodic) THEN
      LastPos = PartState(1:3,iPart)
      CALL PeriodicMovement(iPart)
      IF (ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID).EQ.-1) THEN
        DO WHILE ( .NOT.ALMOSTEQUAL(LastPos(1),PartState(1,iPart)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(2),PartState(2,iPart)) &
            .OR.   .NOT.ALMOSTEQUAL(LastPos(3),PartState(3,iPart)) )
          LastPos = PartState(1:3,iPart)
          CALL PeriodicMovement(iPart)
        END DO
      END IF
    END IF

    ! Not sure if RefNewtonStart value can be reused at this point. Might be safer to start over new
    ! CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID,DoReUseMap=.TRUE.)
    CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)

    ! Position in reference space smaller unity, particle inside
    IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN
      PEM%Element(iPart) = ElemID
#if USE_LOADBALANCE
       ! Particle is on current proc, assign load to new cell
       IF (ElemID.GT.offsetElem+1 .AND. ElemID.LE.offsetElem+nElems) &
         CALL LBElemPauseTime(ElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/
      CYCLE
    END IF
  END IF ! initial check

#if USE_LOADBALANCE
  ! Particle is on current proc, assign load to new cell
  IF (ElemID.GT.offsetElem+1 .AND. ElemID.LE.offsetElem+nElems) THEN
    CALL LBElemSplitTime(ElemID-offsetElem,tLBStart)
  ELSE
    CALL LBStartTime(tLBStart)
  END IF
#endif /*USE_LOADBALANCE*/

  ! relocate particle
  oldElemID = PEM%lastElement(iPart) ! this is not!  a possible elem

  ! get background mesh cell of particle. MAX might give the wrong element if we are right on the edge, so restrict to valid value
  CellX = MAX(FLOOR((PartState(1,iPart)-GEO%xminglob)/GEO%FIBGMdeltas(1)),0) + 1
  CellX = MIN(GEO%FIBGMimax,CellX)
  CellY = MAX(FLOOR((PartState(2,iPart)-GEO%yminglob)/GEO%FIBGMdeltas(2)),0) + 1
  CellY = MIN(GEO%FIBGMjmax,CellY)
  CellZ = MAX(FLOOR((PartState(3,iPart)-GEO%zminglob)/GEO%FIBGMdeltas(3)),0) + 1
  CellZ = MIN(GEO%FIBGMkmax,CellZ)

  ! check all cells associated with this background mesh cell
  nBGMElems = FIBGM_nElems(CellX,CellY,CellZ)

  ! Multiple potential cells found in BGM. Get closest element barycenter by looping over all elements in BGMcell
  IF (nBGMElems.GT.1) THEN
    Distance     = -HUGE(1.)
    ListDistance = -1
    DO iBGMElem = 1, nBGMElems
      ElemID   = FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+iBGMElem)
      CNElemID = GetCNElemID(ElemID)
      ListDistance(iBGMElem) = ElemID

      ! no element associated with BGM elelemt
      IF (ElemID.EQ.-1) &
        CALL Abort(__STAMP__,'Error during RefMapping: unable to find element associated with BGM element!')

      ! oldElemID was previously checked and particle not found inside
      IF (ElemID.EQ.OldElemID) THEN
        Distance(iBGMElem) = -HUGE(1.)
      ELSE
        Distance(iBGMElem) = SUM((PartState(1:3,iPart)-ElemBaryNGeo(1:3,CNElemID))**2)

        ! Do not consider the element if it is too far away
        IF (Distance(iBGMElem).GT.ElemRadius2NGeo(CNElemID)) THEN
          Distance(iBGMElem) = -HUGE(1.)
        END IF
      END IF
    END DO ! nBGMElems
    CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

  ! Only one potential cell found in BGM. No need to sort
  ELSE IF (nBGMElems.EQ.1) THEN
    Distance(1)     = 0.
    ListDistance(1) = FIBGM_Element(FIBGM_offsetElem(CellX,CellY,CellZ)+1)
  END IF

  oldXi     = PartPosRef(1:3,iPart)
  newXi     = HUGE(1.0)
  newElemID = -1

  ! loop through sorted list and start by closest element
  DO iBGMElem = 1,nBGMElems
    ! ignore old element and elements out of range
    IF (Distance(iBGMELem).EQ.-HUGE(1.)) CYCLE

    ElemID = ListDistance(iBGMElem)
#if USE_LOADBALANCE
    ! Cell is on current proc, assign load to new cell
    IF (ElemID.GT.offsetElem+1 .AND. ElemID.LE.offsetElem+nElems) &
      nTracksPerElem(ElemID-offsetElem) = nTracksPerElem(ElemID-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
    CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),ElemID)

    ! Position in reference space smaller unity, particle inside
    IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) THEN
      PEM%Element(iPart) = ElemID
      PartIsDone         = .TRUE.
      EXIT
    END IF

    ! Position in reference space greater then unity, so particle is not in the current cell. Use it as new guess if it is better than the old guess
    IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.MAXVAL(ABS(newXi))) THEN
      newXi     = PartPosRef(1:3,iPart)
      newElemID = ElemID
    END IF
  END DO ! iBGMElem

  ! Particle not located, continue by using the best xi
  IF (.NOT.PartIsDone) THEN
    ! Old guess was the best, keep it
    IF (MAXVAL(ABS(oldXi)).LT.MAXVAL(ABS(newXi))) THEN
      PartPosRef(1:3,iPart) = oldXi
      PEM%Element(iPart)    = oldElemID
    ! New guess is better, so overwrite the old one
    ELSE
      PartPosRef(1:3,iPart) = newXi
      PEM%Element(iPart)    = NewElemID
      oldElemID             = NewElemID
    END IF

    ! Set test element
    TestElem   = PEM%Element(iPart)
    CNTestElem = GetCNElemID(TestElem)

    IF (CNTestElem.LE.0) THEN
      IPWRITE(UNIT_stdOut,'(I0,A)') ' Particle element not on local compute node!'
      CALL PrintParticleInfo(iPart,CNTestElem,oldXi,newXi)
      CALL Abort(__STAMP__ ,'Particle not inside of element, iPart',iPart)

    ELSE
      ! Position in reference space is outside tolerance of the test element
      IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.ElemEpsOneCell(CNTestElem)) THEN
        PartIsDone = .FALSE.
        ! Element is not a boundary element
        IF (ElemToBCSides(ELEM_NBR_BCSIDES,CNTestElem).LE.0) THEN
          ! Debug output
          IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance issue with internal element!'
          CALL PrintParticleInfo(iPart,CNElemID,oldXi,newXi)
#if USE_MPI
          InElem   = PEM%Element(iPart)
          HaloStat = MERGE('T','F',InElem.GE.offsetElem+1 .AND. InElem.LE.offsetElem+nElems)
          IPWRITE(UNIT_stdOut,'(I0,A,I0,A,A,A)')   ' ElemID       ', InElem,' (halo-elem = ',HaloStat,')'

          InElem = PEM%LastElement(iPart)
          HaloStat = MERGE('T','F',InElem.GE.offsetElem+1 .AND. InElem.LE.offsetElem+nElems)
          IPWRITE(UNIT_stdOut,'(I0,A,I0,A,A,A)')   ' Last-ElemID  ', InElem,' (halo-elem = ',HaloStat,')'
#else
          IPWRITE(UNIT_stdOut,'(I0,A,I0)')         ' ElemID       ', PEM%Element(iPart)
          IPWRITE(UNIT_stdOut,'(I0,A,I0)')         ' Last-ElemID  ', PEM%LastElement(iPart)
#endif
          CALL Abort(__STAMP__,'Particle not inside of element, iPart',iPart)

        ELSE ! BCElem
          IPWRITE(UNIT_stdOut,'(I0,A,1X,I0,A,1X,I0)') ' Fallback for particle ', iPart, ' in element', TestElem
          IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))')   ' ParticlePos           ', PartState(1:3,iPart)
          Vec = PartState(1:3,iPart) - LastPartPos(1:3,iPart)

          ! Tracking on curved meshes with NGeo>1. Check if the particle intersected with face and left the element
          IF (useCurveds .AND. NGeo.GT.1) &
            CALL FallBackFaceIntersection(TestElem                                                                                    &
                                         ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + 1                                             &
                                         ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + ElemToBCSides(ELEM_NBR_BCSIDES,CNTestElem)    &
                                         ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNTestElem)                                                 &
                                         ,iPart)

          ! No fall back algorithm algorithm, try the normal tracking for boundary elements
          lengthPartTrajectory0 = 0.
          CALL ParticleBCTracking(lengthPartTrajectory0                                                                               &
                                 ,TestElem                                                                                            &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + 1                                                     &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNTestElem) + ElemToBCSides(ELEM_NBR_BCSIDES,CNTestElem)            &
                                 ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNTestElem)                                                         &
                                 ,iPart                                                                                               &
                                 ,PartIsDone                                                                                          &
                                 ,PartIsMoved                                                                                         &
                                 ,1)
          IF (PartIsDone) CYCLE

          CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),TestElem)

          ! Position in reference space greater than unity, particle not inside element. Try to relocate
          IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).GT.ElemEpsOneCell(CNTestElem)) THEN
            IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance issue with BC element, relocating!'
            CALL LocateParticleInElement(iPart,doHALO=.TRUE.)

            ! Re-localization to a new element failed
            IF (.NOT.PDM%ParticleInside(iPart)) THEN
              IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance issue with BC element!'
              CALL PrintParticleInfo(iPart,CNTestElem,oldXi,newXi)
              CALL Abort(__STAMP__ ,'Particle not inside of element, iPart',iPart)
            END IF ! .NOT.PDM%ParticleInside(iPart)
          ! Xi.LT.epsOneCell
          ELSE
            PEM%Element(iPart) = TestElem
          END IF ! Xi.LT.epsOneCell
        END IF ! BCElem
      END IF ! Xi.LT.epsOneCell
    END IF ! CNTestElem.GT.1
  END IF ! Particle not located
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO ! iPart

END SUBROUTINE ParticleRefTracking


RECURSIVE SUBROUTINE ParticleBCTracking(lengthPartTrajectory0 &
                                       ,ElemID,firstSide,lastSide,nlocSides,PartId,PartisDone,PartisMoved,iCount)
!===================================================================================================================================
! Calculate intersection with boundary and choose boundary interaction type for reference tracking routine
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars                   ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_Particle_Globals            ,ONLY: VECNORM
USE MOD_Particle_Localization       ,ONLY: PARTHASMOVED
USE MOD_Particle_Mesh_Vars          ,ONLY: SideBCMetrics,ElemToBCSides
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO,ElemRadiusNGeo
USE MOD_Particle_Mesh_Tools         ,ONLY: GetCNElemID,GetCNSideID
USE MOD_Particle_Intersection       ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection       ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Surfaces_Vars      ,ONLY: SideType
USE MOD_Particle_Tracking_Vars      ,ONLY: CartesianPeriodic
USE MOD_Particle_Vars               ,ONLY: PEM,PDM
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
USE MOD_Utils                       ,ONLY: InsertionSort
#if CODE_ANALYZE
USE MOD_Particle_Tracking_Vars,      ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,firstSide,LastSide,nlocSides
INTEGER,INTENT(IN)            :: iCount
LOGICAL,INTENT(INOUT)         :: PartisDone
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES!
INTEGER,INTENT(INOUT)         :: ElemID
LOGICAL,INTENT(INOUT)         :: PartisMoved
REAL,INTENT(INOUT)            :: lengthPartTrajectory0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Counters
INTEGER                       :: ilocSide
INTEGER                       :: nInter
! Elements
INTEGER                       :: OldElemID
INTEGER                       :: CNElemID,CNOldElemID
! Sides
INTEGER                       :: SideID,CNSideID,flip
! Particles
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
LOGICAL                       :: DoTracing,PeriMoved,Reflected
! Tracking
INTEGER                       :: locSideList(firstSide:lastSide),hitlocSide
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
REAL                          :: alphaOld
LOGICAL                       :: isHit,doubleCheck
!===================================================================================================================================

CNElemID   = GetCNElemID(ElemID)

! Calculate particle trajectory
PartTrajectory       = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))

! Check if the particle moved at all. If not, tracking is done
IF (.NOT.PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo(CNElemID))) THEN
  PEM%Element(PartID) = ElemID
  PartisDone = .TRUE.
  RETURN
END IF

PartTrajectory        = PartTrajectory/lengthPartTrajectory
lengthPartTrajectory0 = MAX(lengthPartTrajectory0,lengthPartTrajectory)

! Init variables. Double check if LastPartPos is close to side and first intersection is found for this position (negative alpha)
PartisMoved = .FALSE.
DoTracing   = .TRUE.
doubleCheck = .FALSE.
alphaOld    = -1.

DO WHILE(DoTracing)
  IF (GEO%nPeriodicVectors.GT.0.AND.CartesianPeriodic) THEN
    ! Account for periodic displacement
    CALL PeriodicMovement(PartID,PeriMoved)
    ! Position and trajectory has to be recomputed after periodic displacement
    IF (PeriMoved) THEN
      PartTrajectory       = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
      lengthPartTrajectory = VECNORM(PartTrajectory(1:3))
    END IF
  ELSE
    PeriMoved = .FALSE.
  END IF
  locAlpha = -1.
  nInter   = 0

  ! track particle vector until the final particle position is achieved
  ! check if particle can intersect with current side
  DO iLocSide = firstSide,LastSide

    ! SideBCMetrics is sorted by distance. stop if the first side is out of range
    IF (SideBCMetrics(BCSIDE_DISTANCE,ilocSide).GT.lengthPartTrajectory0) EXIT

    ! side potentially in range (halo_eps)
    SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
    CNSideID = GetCNSideID(SideID)
    locSideList(ilocSide) = ilocSide

    ! BezierControlPoints are now built in cell local system. Hence, sides have always the flip from the shared SideInfo
    flip = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

    ! double check
    IF (doublecheck) THEN
#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.myRank) THEN; IF (PartID.EQ.PARTOUT) THEN
        WRITE(UNIT_stdOut,'(110("="))')
        WRITE(UNIT_stdOut,'(A)')    '     | Particle is double checked: '
      END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

      SELECT CASE(SideType(CNSideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarRectIntersection(  isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xi                   = xi(ilocSide)          &
                                            ,  eta                  = eta(ilocSide)         &
#if CODE_ANALYZE
                                            ,  PartID               = PartID                &
#endif /*CODE_ANALYZE*/
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID)
        CASE(BILINEAR,PLANAR_NONRECT)
          CALL ComputeBiLinearIntersection(    isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  PartState            = PartState(1:3,PartID) &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xitild               = xi(ilocSide)          &
                                            ,  etatild              = eta(ilocSide)         &
                                            ,  PartID               = PartID                &
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID                &
                                            ,  alpha2               = alphaOld)
        CASE(PLANAR_CURVED)
          CALL ComputePlanarCurvedIntersection(isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xi                   = xi(ilocSide)          &
                                            ,  eta                  = eta(ilocSide)         &
                                            ,  PartID               = PartID                &
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID)
        CASE(CURVED)
          CALL ComputeCurvedIntersection(      isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  PartState            = PartState(1:3,PartID) &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xi                   = xi(ilocSide)          &
                                            ,  eta                  = eta(ilocSide)         &
                                            ,  PartID               = PartID                &
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID)
      END SELECT

#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.myRank) THEN; IF (PartID.EQ.PARTOUT) THEN
        WRITE(UNIT_stdOut,'(30("-"))')
        WRITE(UNIT_stdOut,'(A)')           '     | Output after compute intersection (DoubleCheck REFMAPPING): '
        WRITE(UNIT_stdOut,'(2(A,I0),A,L)') '     | SideType: ',SideType(CNSideID),' | SideID: ',SideID,' | Hit: ',isHit
        WRITE(UNIT_stdOut,'(2(A,G0))')     '     | Alpha: ',locAlpha(ilocSide)   ,' | LengthPartTrajectory: ', lengthPartTrajectory
        WRITE(UNIT_stdOut,'((A,G0))')      '     | AlphaOld: ',alphaOld
        WRITE(UNIT_stdOut,'(A,2(1X,G0))')  '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
      END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/

    ! not double check
    ELSE
      SELECT CASE(SideType(CNSideID))
        CASE(PLANAR_RECT)
          CALL ComputePlanarRectIntersection(  isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xi                   = xi(ilocSide)          &
                                            ,  eta                  = eta(ilocSide)         &
#if CODE_ANALYZE
                                            ,   PartID               = PartID                &
#endif /*CODE_ANALYZE*/
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID)
        CASE(BILINEAR,PLANAR_NONRECT)
          CALL ComputeBiLinearIntersection(    isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  PartState            = PartState(1:3,PartID) &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xitild               = xi(ilocSide)          &
                                            ,  etatild              = eta(ilocSide)         &
                                            ,  PartID               = PartID                &
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID                &
                                            ,  alpha2               = alphaOld)
        CASE(PLANAR_CURVED)
          CALL ComputePlanarCurvedIntersection(isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xi                   = xi(ilocSide)          &
                                            ,  eta                  = eta(ilocSide)         &
                                            ,  PartID               = PartID                &
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID)
        CASE(CURVED)
          CALL ComputeCurvedIntersection(      isHit                = isHit                 &
                                            ,  PartTrajectory       = PartTrajectory        &
                                            ,  lengthPartTrajectory = lengthPartTrajectory  &
                                            ,  PartState            = PartState(1:3,PartID) &
                                            ,  LastPartPos          = LastPartPos(:,PartID) &
                                            ,  alpha                = locAlpha(ilocSide)    &
                                            ,  xi                   = xi(ilocSide)          &
                                            ,  eta                  = eta(ilocSide)         &
                                            ,  PartID               = PartID                &
                                            ,  flip                 = flip                  &
                                            ,  SideID               = SideID)
    END SELECT

#if CODE_ANALYZE
!---------------------------------------------CODE_ANALYZE--------------------------------------------------------------------------
      IF (PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.myRank) THEN; IF (PartID.EQ.PARTOUT) THEN
        WRITE(UNIT_stdOut,'(30("-"))')
        WRITE(UNIT_stdOut,'(A)')           '     | Output after compute intersection (REFMAPPING, BCTracing): '
        WRITE(UNIT_stdOut,'(2(A,I0),A,L)') '     | SideType: ',SideType(CNSideID),' | SideID: ',SideID,' | Hit: ',isHit
        WRITE(UNIT_stdOut,'(2(A,G0))')     '     | Alpha: ',locAlpha(ilocSide)   ,' | LengthPartTrajectory: ', lengthPartTrajectory
        WRITE(UNIT_stdOut,'(A,2(1X,G0))')  '     | Intersection xi/eta: ',xi(ilocSide),eta(ilocSide)
      END IF; END IF
!-------------------------------------------END-CODE_ANALYZE------------------------------------------------------------------------
#endif /*CODE_ANALYZE*/
    END IF

    IF (locAlpha(ilocSide).GT.-1.0) nInter = nInter+1
  END DO ! ilocSide

  ! No periodic movement in first intersection
  IF (nInter.EQ.0) THEN
    IF (.NOT.PeriMoved) DoTracing = .FALSE.
  ! Take first possible intersection
  ELSE
    PartIsMoved = .TRUE.
    CALL InsertionSort(locAlpha,locSideList,nlocSides)

    DO ilocSide = firstSide,lastSide
      ! Stop on first intersection
      IF (locAlpha(ilocSide).GT.-1) THEN
        alphaOld = locAlpha(ilocSide)
        EXIT
      END IF
    END DO

    DO ilocSide = firstSide,lastSide
      IF (locAlpha(ilocSide).GT.-1) THEN
        hitlocSide = locSideList(ilocSide)
        SideID     = INT(SideBCMetrics(BCSIDE_SIDEID,hitlocSide))
        flip       = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)
        OldElemID  = ElemID
        CALL GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,locAlpha(ilocSide) &
                                   ,xi(hitlocSide),eta(hitlocSide),PartId,SideID,flip      &
                                   ,ElemID,reflected)

        ! particle moved to a new element in boundary interaction
        IF (ElemID.NE.OldElemID) THEN
          ! Try to recursively calculate the intersection 1000 times. Threshold might be changed...
          IF (iCount.GE.1000 .AND. MOD(iCount,1000).EQ.0) THEN
            IPWRITE(UNIT_stdOut,'(I4,A,I0,A,3(1X,I0))') ' WARNING: proc has called BCTracking ',iCount  &
              ,'x recursively! Part, Side, Elem:',PartId,SideID,ElemID
          END IF

          ! check if a periodic boundary was crossed during boundary interaction
          CNOldElemID = GetCNElemID(OldElemID)

          IF(GEO%nPeriodicVectors.GT.0)THEN
            lengthPartTrajectory0 = MAXVAL(SideBCMetrics(BCSIDE_DISTANCE,                         &
                                           ElemToBCSides(ELEM_FIRST_BCSIDE,CNOldElemID) + 1:      &
                                           ElemToBCSides(ELEM_FIRST_BCSIDE,CNOldElemID) + ElemToBCSides(ELEM_NBR_BCSIDES,CNOldElemID)))
          END IF

          CNElemID = GetCNElemID(ElemID)
          CALL ParticleBCTracking(lengthPartTrajectory0                                                                 &
                                 ,ElemID                                                                                &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + 1                                         &
                                 ,ElemToBCSides(ELEM_FIRST_BCSIDE,CNElemID) + ElemToBCSides(ELEM_NBR_BCSIDES,CNElemID)  &
                                 ,ElemToBCSides(ELEM_NBR_BCSIDES ,CNElemID)                                             &
                                 ,PartID                                                                                &
                                 ,PartIsDone                                                                            &
                                 ,PartIsMoved                                                                           &
                                 ,iCount+1)
          PartisMoved=.TRUE.
          RETURN
        END IF

        ! Particle got reflected and stays in the same element. Done with this particle
        IF (reflected) EXIT
      END IF
    END DO ! ilocSide

    ! Particle left the domain through a boundary
    IF (.NOT.PDM%ParticleInside(PartID)) THEN
      PartisDone = .TRUE.
      RETURN
    END IF

    ! Particle did not leave the domain and did not get reflected
    IF (.NOT.reflected) THEN
      ! Double check the particle if not done already
      IF (.NOT.doubleCheck) THEN
        doubleCheck = .TRUE.
      ! Stop tracing if already double checked
      ELSE
        DoTracing = .FALSE.
      END IF
    END IF
  END IF ! nInter>0
END DO

END SUBROUTINE ParticleBCTracking


SUBROUTINE PeriodicMovement(PartID,isMovedOut)
!----------------------------------------------------------------------------------------------------------------------------------!
! move particle in the periodic direction, if particle is outside of the box
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars,               ONLY:PartState,LastPartPos
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,          ONLY:GEO
USE MOD_Particle_Tracking_Vars,      ONLY:FastPeriodic
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: PartID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL   :: isMovedOut
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPV,iDir
REAL                            :: MoveVector(1:3),GEOLims(2)
LOGICAL                         :: isMoved
CHARACTER(LEN=1)                :: DirChar
!===================================================================================================================================

isMoved = .FALSE.

! Routine if particle domain is rectangular. Periodic displacement can be determined by the bounding box
IF(FastPeriodic)THEN
  DO iDir = 1,3
    IF (.NOT.GEO%directions(iDir)) CYCLE

    SELECT CASE(iDir)
      CASE(1)
        GEOLims = (/GEO%xminglob,GEO%xmaxglob/)
        DirChar = 'x'
      CASE(2)
        GEOLims = (/GEO%yminglob,GEO%ymaxglob/)
        DirChar = 'y'
      CASE(3)
        GEOLims = (/GEO%zminglob,GEO%zmaxglob/)
        DirChar = 'z'
    END SELECT

    ! Check against upper limit
    IF(PartState(iDir,PartID).GT.GEOLims(2)) THEN
      DO iPV = 1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.iDir) EXIT
      END DO
      MoveVector = CEILING(ABS(PartState(iDir,PartID)-GEOLims(2))/ABS(GEO%PeriodicVectors(iDir,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      MoveVector = MoveVector * SIGN(1.,GEO%PeriodicVectors(iDir,iPV))
      PartState(  1:3,PartID) = PartState(  1:3,PartID) - MoveVector
      LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) - MoveVector
      isMoved = .TRUE.
    END IF

    ! Check against lower limit
    IF(PartState(iDir,PartID).LT.GEOLims(1)) THEN
      DO iPV = 1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.iDir) EXIT
      END DO
      MoveVector = CEILING(ABS(PartState(iDir,PartID)-GEOLims(1))/ABS(GEO%PeriodicVectors(iDir,iPV)))*GEO%PeriodicVectors(1:3,iPV)
      MoveVector = MoveVector * SIGN(1.,GEO%PeriodicVectors(iDir,iPV))
      PartState  (1:3,PartID) = PartState(  1:3,PartID) + MoveVector
      LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + MoveVector
      isMoved = .TRUE.
    END IF

    ! Check if particle is inside
    IF(PartState(iDir,PartID).GT.GEOLims(2)) THEN
      IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)') 'PartPos', PartState(:,PartID)
      CALL Abort(__STAMP__,'Particle outside '//DirChar//'+, PartID',PartID)
    END IF
    IF(PartState(iDir,PartID).LT.GEOLims(1)) THEN
      IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)') 'PartPos', PartState(:,PartID)
      CALL Abort(__STAMP__,'Particle outside '//DirChar//'-, PartID',PartID)
    END IF
  END DO

! No fast periodic displacement
ELSE
  DO iDir = 1,3
    IF (.NOT.GEO%directions(iDir)) CYCLE

    SELECT CASE(iDir)
      CASE(1)
        GEOLims = (/GEO%xminglob,GEO%xmaxglob/)
      CASE(2)
        GEOLims = (/GEO%yminglob,GEO%ymaxglob/)
      CASE(3)
        GEOLims = (/GEO%zminglob,GEO%zmaxglob/)
    END SELECT

    ! Check against upper limit
    IF(PartState(iDir,PartID).GT.GEOLims(2)) THEN
      DO iPV = 1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.iDir) EXIT
      END DO
      MoveVector = GEO%PeriodicVectors(1:3,iPV) * SIGN(1.,GEO%PeriodicVectors(iDir,iPV))
      PartState  (1:3,PartID) = PartState  (1:3,PartID) - MoveVector
      LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) - MoveVector
      isMoved = .TRUE.
    END IF

    ! Check against lower limit
    IF(PartState(iDir,PartID).LT.GEOLims(1)) THEN
      DO iPV = 1,GEO%nPeriodicVectors
        IF(GEO%DirPeriodicVectors(iPV).EQ.iDir) EXIT
      END DO
      MoveVector = GEO%PeriodicVectors(1:3,iPV) * SIGN(1.,GEO%PeriodicVectors(iDir,iPV))
      PartState  (1:3,PartID) = PartState  (1:3,PartID) + MoveVector
      LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + MoveVector
      isMoved = .TRUE.
    END IF
  END DO
END IF

IF(PRESENT(isMovedOut)) isMovedOut = isMoved

END SUBROUTINE PeriodicMovement


SUBROUTINE FallBackFaceIntersection(ElemID,firstSide,LastSide,nlocSides,PartID)
!===================================================================================================================================
! Checks if lost particle intersected with face and left Element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,                   ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Condition, ONLY: GetBoundaryInteraction
USE MOD_Particle_Globals,            ONLY: VECNORM
USE MOD_Particle_Localization,       ONLY: LocateParticleInElement
USE MOD_Particle_Intersection,       ONLY: ComputeCurvedIntersection
USE MOD_Particle_Intersection,       ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection,       ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection,       ONLY: ComputePlanarNonRectIntersection
USE MOD_Particle_Intersection,       ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Mesh_Vars,          ONLY: SideBCMetrics
USE MOD_Particle_Mesh_Vars,          ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Tools,         ONLY: GetCNElemID,GetCNSideID
USE MOD_Particle_Surfaces_Vars,      ONLY: SideType
USE MOD_Utils,                       ONLY: InsertionSort
USE MOD_Particle_Vars,               ONLY: PDM,PartState,LastPartPos
! USE MOD_Particle_Vars,               ONLY: PartPosRef
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,ElemID,firstSide,LastSide,nlocSides
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: ilocSide,SideID,CNSideID,locSideList(firstSide:lastSide),hitlocSide
LOGICAL                       :: dolocSide(firstSide:lastSide)
LOGICAL                       :: ishit
REAL                          :: localpha(firstSide:lastSide),xi(firstSide:lastSide),eta(firstSide:lastSide)
INTEGER                       :: nInter,flip
REAL                          :: tmpPos(3),tmpLastPartPos(3),tmpVec(3)
REAL                          :: PartTrajectory(1:3),lengthPartTrajectory
!===================================================================================================================================

! Get particle position and trajectory
tmpPos                  = PartState  (1:3,PartID)
tmpLastPartPos(1:3)     = LastPartPos(1:3,PartID)
tmpVec                  = PartTrajectory
LastPartPos(1:3,PartID) = ElemBaryNGeo(:,GetCNElemID(ElemID))

PartTrajectory       = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))
IF (lengthPartTrajectory.GT.0) PartTrajectory = PartTrajectory/lengthPartTrajectory

locAlpha  = -1.0
nInter    = 0
dolocSide = .TRUE.

! Loop over all sides
DO iLocSide=firstSide,LastSide
  ! track particle vector until the final particle position is achieved
  SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
  CNSideID = GetCNSideID(SideID)
  locSideList(ilocSide) = ilocSide
  flip     = MERGE(0, MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT)
      CALL ComputePlanarRectIntersection(   isHit                = isHit                 &
                                        ,   PartTrajectory       = PartTrajectory        &
                                        ,   lengthPartTrajectory = lengthPartTrajectory  &
                                        ,   LastPartPos          = tmpLastPartPos        &
                                        ,   alpha                = locAlpha(ilocSide)    &
                                        ,   xi                   = xi(ilocSide)          &
                                        ,   eta                  = eta(ilocSide)         &
#if CODE_ANALYZE
                                           ,PartID               = PartID                &
#endif /*CODE_ANALYZE*/
                                        ,   flip                 = flip                  &
                                        ,   SideID               = SideID)
    CASE(PLANAR_NONRECT)
      CALL ComputePlanarNonRectIntersection(isHit                = isHit                 &
                                        ,   PartTrajectory       = PartTrajectory        &
                                        ,   lengthPartTrajectory = lengthPartTrajectory  &
                                        ,   LastPartPos          = tmpLastPartPos        &
                                        ,   alpha                = locAlpha(ilocSide)    &
                                        ,   xi                   = xi(ilocSide)          &
                                        ,   eta                  = eta(ilocSide)         &
#if CODE_ANALYZE
                                        ,   PartID               = PartID                &
#endif /*CODE_ANALYZE*/
                                        ,   flip                 = flip                  &
                                        ,   SideID               = SideID)
    CASE(BILINEAR)
      CALL ComputeBiLinearIntersection(     isHit                = isHit                 &
                                        ,   PartTrajectory       = PartTrajectory        &
                                        ,   lengthPartTrajectory = lengthPartTrajectory  &
                                        ,   PartState            = tmpPos                &
                                        ,   LastPartPos          = tmpLastPartPos        &
                                        ,   alpha                = locAlpha(ilocSide)    &
                                        ,   xitild               = xi(ilocSide)          &
                                        ,   etatild              = eta(ilocSide)         &
                                        ,   PartID               = PartID                &
                                        ,   flip                 = flip                  &
                                        ,   SideID               = SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection( isHit                = isHit                 &
                                        ,   PartTrajectory       = PartTrajectory        &
                                        ,   lengthPartTrajectory = lengthPartTrajectory  &
                                        ,   LastPartPos          = tmpLastPartPos        &
                                        ,   alpha                = locAlpha(ilocSide)    &
                                        ,   xi                   = xi(ilocSide)          &
                                        ,   eta                  = eta(ilocSide)         &
                                        ,   PartID               = PartID                &
                                        ,   flip                 = flip                  &
                                        ,   SideID               = SideID)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(       isHit                = isHit                 &
                                        ,   PartTrajectory       = PartTrajectory        &
                                        ,   lengthPartTrajectory = lengthPartTrajectory  &
                                        ,   PartState            = tmpPos                &
                                        ,   LastPartPos          = tmpLastPartPos        &
                                        ,   alpha                = locAlpha(ilocSide)    &
                                        ,   xi                   = xi(ilocSide)          &
                                        ,   eta                  = eta(ilocSide)         &
                                        ,   PartID               = PartID                &
                                        ,   flip                 = flip                  &
                                        ,   SideID               = SideID)
  END SELECT

  IF (locAlpha(ilocSide).GT.-1.0) THEN
    nInter = nInter+1
  END IF
END DO ! ilocSide

! no intersection found. Try to locate particle manually
IF (nInter.EQ.0) THEN
  PartState(  1:3,PartID) = tmpPos
  LastPartPos(1:3,PartID) = tmpLastPartPos(1:3)
!  IF(PartPosRef(1,PartID).GT. 1.) PartPosRef(1,PartID)= 0.99
!  IF(PartPosRef(1,PartID).LT.-1.) PartPosRef(1,PartID)=-0.99
!  IF(PartPosRef(2,PartID).GT. 1.) PartPosRef(2,PartID)= 0.99
!  IF(PartPosRef(2,PartID).LT.-1.) PartPosRef(2,PartID)=-0.99
!  IF(PartPosRef(3,PartID).GT. 1.) PartPosRef(3,PartID)= 0.99
!  IF(PartPosRef(3,PartID).LT.-1.) PartPosRef(3,PartID)=-0.99
  CALL LocateParticleInElement(PartID,doHalo=.TRUE.)

  ! particle successfully located
  IF (PDM%ParticleInside(PartID)) THEN
    RETURN
  ELSE
    CALL Abort(__STAMP__,'FallBackFaceIntersection failed for particle!')
  END IF

! One or more intersection
ELSE
  ! take first possible intersection and place particle "just a little" further back
  CALL InsertionSort(locAlpha,locSideList,nlocSides)
  DO ilocSide = firstSide,LastSide
    IF (locAlpha(ilocSide).GT.-1.0) THEN
      hitlocSide = locSideList(ilocSide)
      SideID     = INT(SideBCMetrics(BCSIDE_SIDEID,ilocSide))
      LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + 0.97*locAlpha(ilocSide)*PartTrajectory
      PartState  (1:3,PartID) = LastPartPos(1:3,PartID)
    END IF ! locAlpha.GT.-1.0
  END DO ! ilocSide
END IF ! nInter.GT.0

END SUBROUTINE FallBackFaceIntersection


SUBROUTINE PrintParticleInfo(PartID,CNElemID,oldXi,newXi)
!===================================================================================================================================
! Prints information for particle outside reference element during tracking
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,ElemEpsOneCell
USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps2
USE MOD_Particle_Vars          ,ONLY: PartState,TurbPartState,PartPosRef,LastPartPos,PartSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID,CNElemID
REAL,INTENT(IN)               :: oldXi(3),newXi(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: Vec(3)
CHARACTER(LEN=37),PARAMETER   :: FormatPos = '(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)'
!===================================================================================================================================

Vec = PartState(1:3,PartID) - LastPartPos(1:3,PartID)

IPWRITE(UNIT_stdOut,  '(I0,A,3(1X,E15.8))')     ' xi                    ', PartPosRef(   1:3,PartID)
IF (CNElemID.GT.0) &
  IPWRITE(UNIT_stdOut,'(I0,A,1X,E15.8)')        ' epsOneCell            ', ElemEpsOneCell(CNElemID)
IPWRITE(UNIT_stdOut,  '(I0,A,3(1X,E15.8))')     ' oldXi                 ', oldXi
IPWRITE(UNIT_stdOut,  '(I0,A,3(1X,E15.8))')     ' newXi                 ', newXi
IPWRITE(UNIT_stdOut,  '(I0,A)')                 ' PartPos:'
IPWRITE(UNIT_stdOut,  '(I0,A18,1X,A18,1X,A18)') '    min ', ' value ', ' max'
IPWRITE(UNIT_stdOut,  FormatPos)                ' x', GEO%xminglob, PartState(1,PartID), GEO%xmaxglob
IPWRITE(UNIT_stdOut,  FormatPos)                ' y', GEO%yminglob, PartState(2,PartID), GEO%ymaxglob
IPWRITE(UNIT_stdOut,  FormatPos)                ' z', GEO%zminglob, PartState(3,PartID), GEO%zmaxglob
IPWRITE(UNIT_stdOut,  '(I0,A,3(1X,E15.8))')     ' LastPartPos           ', LastPartPos(  1:3,PartID)
IPWRITE(UNIT_stdOut,  '(I0,A,3(1X,E15.8))')     ' Velocity              ', PartState(    4:6,PartID)
IF (ALLOCATED(TurbPartState)) &
  IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))')     ' Velocity (SGS)        ', TurbPartState(1:3,PartID)
IPWRITE(UNIT_stdOut,  '(I0,A,1X,E15.8)')        ' displacement/halo_eps ', DOT_PRODUCT(Vec,Vec)/halo_eps2
IPWRITE(UNIT_stdOut,  '(I0,A,I0)')              ' PartSpecies           ', PartSpecies(PartID)

END SUBROUTINE PrintParticleInfo

END MODULE MOD_Particle_RefTracking
