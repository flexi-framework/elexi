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
MODULE MOD_Particle_TriaTracking
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ParticleTriaTracking
  MODULE PROCEDURE ParticleTriaTracking
END INTERFACE

PUBLIC::ParticleTriaTracking
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleTriaTracking()
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
! 1) Loop over all particles that are still inside
!    2) Perform tracking until the particle is considered "done" (either localized or deleted)
!       2a) Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
!           triangle (ParticleInsideQuad3D)
!       2b) If particle is not within the given element in a), the side through which the particle went is determined by checking
!           each side of the element (ParticleThroughSideCheck3DFast)
!       2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
!           treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
!           the determinants
!    3) In case of a boundary, determine the intersection and perform the appropriate boundary interaction (GetBoundaryInteraction)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Globals            ,ONLY: DOTPRODUCT,VECNORM
USE MOD_Particle_Boundary_Condition ,ONLY: GetBoundaryInteraction
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound
USE MOD_Particle_Localization       ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Intersection       ,ONLY: IntersectionWithWall
USE MOD_Particle_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars          ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_Particle_Tracking_Vars      ,ONLY: TrackInfo !CountNbOfLostParts,NbrOfLostParticles
USE MOD_Particle_Vars               ,ONLY: PEM,PDM,PartSpecies
USE MOD_Particle_Vars               ,ONLY: PartState,LastPartPos
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime, LBElemSplitTime, LBElemPauseTime
USE MOD_Mesh_Vars                   ,ONLY: offsetElem
USE MOD_Particle_Globals            ,ONLY: PP_nElems
USE MOD_Particle_Tracking_Vars      ,ONLY: nTracks,MeasureTrackTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER                   :: eps = 0, xi = -1., eta = -1.
REAL                             :: alpha = -1.
INTEGER                          :: iPart,NblocSideID,NbElemID,ind,nbSideID,nMortarElems,BCType
INTEGER                          :: ElemID,CNElemID,flip,OldElemID,nlocSides
INTEGER                          :: LocalSide
INTEGER                          :: NrOfThroughSides, ind2
INTEGER                          :: SideID,TempSideID,iLocSide, localSideID
INTEGER                          :: TriNum,LocSidesTemp(1:6),TriNumTemp(1:6),GlobSideTemp(1:6)
INTEGER                          :: SecondNrOfThroughSides, indSide
INTEGER                          :: DoneLastElem(1:4,1:6) ! 1:3: 1=Element,2=LocalSide,3=TriNum 1:2: 1=last 2=beforelast
LOGICAL                          :: ThroughSide, InElementCheck,PartisDone
LOGICAL                          :: crossedBC,oldElemIsMortar, isMortarSideTemp(1:6), doCheckSide
REAL                             :: det(6,2),detM,ratio,minRatio,detPartPos
REAL                             :: PartTrajectory(1:3),lengthPartTrajectory
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! -- 1. Loop over all particles that are still inside
DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN

#if USE_LOADBALANCE
    IF (MeasureTrackTime) nTracks = nTracks+1
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

    ! Reset all tracking variables
    PartisDone         = .FALSE.
    ElemID             = PEM%lastElement(iPart)
    TrackInfo%CurrElem = ElemID
    SideID             = 0
    DoneLastElem(:,:)  = 0

! -- 2. Loop tracking until particle is considered "done" (either localized or deleted)
    DO WHILE (.NOT.PartisDone)
! -- 2.a Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
!        triangle (ParticleInsideQuad3D)
      oldElemIsMortar = .FALSE.
      CALL ParticleInsideQuad3D(PartState(1:3,iPart) &
                               ,ElemID               &
                               ,InElementCheck       &
                               ,det)
      ! If it is, set new ElementNumber = ElemID and LocalizeOn = .FALSE. ->PartisDone
      IF (InElementCheck) THEN
        ! If particle is inside the given ElemID, set new PEM%Element and stop tracking this particle ->PartisDone
        PEM%Element(iPart) = ElemID
        PartisDone         = .TRUE.
      ELSE
! -- 2.b If particle is not within the given element in a), the side through which the particle went is determined by checking
!        each side of the element (ParticleThroughSideCheck3DFast)
        NrOfThroughSides = 0
        LocSidesTemp(:)  = 0
        TriNumTemp(:)    = 0
        GlobSideTemp     = 0
        isMortarSideTemp = .FALSE.

        ! Calculate particle trajectory
        PartTrajectory       = PartState(1:3,iPart) - LastPartPos(1:3,iPart)
        lengthPartTrajectory = VECNORM(PartTrajectory)
        IF (ABS(lengthPartTrajectory).GT.0.) PartTrajectory = PartTrajectory/lengthPartTrajectory
        nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)

        DO iLocSide = 1,nlocSides
          TempSideID  = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
          localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)

          ! Side is not one of the 6 local sides
          IF (localSideID.LE.0) CYCLE

          NbElemID = SideInfo_Shared(SIDE_NBELEMID,TempSideID)

          ! Side is a regular side
          IF (NbElemID.GE.0) THEN
            DO TriNum = 1,2
              IF (det(localSideID,TriNum).LE.-eps) THEN
                ThroughSide = .FALSE.
                CALL ParticleThroughSideCheck3DFast(iPart          &
                                                   ,PartTrajectory &
                                                   ,localSideID    &
                                                   ,ElemID         &
                                                   ,ThroughSide    &
                                                   ,TriNum)

                ! Store the information for this side for future checks, if this side was already treated
                IF (ThroughSide) THEN
                  NrOfThroughSides = NrOfThroughSides + 1
                  LocSidesTemp(NrOfThroughSides) = localSideID
                  TriNumTemp(  NrOfThroughSides) = TriNum
                  GlobSideTemp(NrOfThroughSides) = TempSideID
                  SideID    = TempSideID
                  LocalSide = localSideID
                END IF
              END IF
            END DO

          ! Mortar side
          ELSE
            nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,TempSideID).EQ.-1)
            DO ind = 1,nMortarElems
              nbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
              NbElemID = SideInfo_Shared(SIDE_NBELEMID,nbSideID)

              ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
              ! no way to recover it during runtime
              IF (NbElemID.LT.1) CALL Abort(__STAMP__,'Small mortar element not defined!',ElemID)

              ! For small mortar sides, SIDE_NBSIDEID contains the SideID of the corresponding big mortar side
              nbSideID    = SideInfo_Shared(SIDE_NBSIDEID,nbSideID)
              NblocSideID = SideInfo_Shared(SIDE_LOCALID ,nbSideID)

              DO TriNum = 1,2
                ThroughSide = .FALSE.
                CALL ParticleThroughSideCheck3DFast(iPart          &
                                                   ,PartTrajectory &
                                                   ,NblocSideID    &
                                                   ,NbElemID       &
                                                   ,ThroughSide    &
                                                   ,TriNum         &
                                                   ,IsMortar=.TRUE.)
                ! Store the information for this side for future checks, if this side was already treated
                IF (ThroughSide) THEN
                  oldElemIsMortar  = .TRUE.
                  NrOfThroughSides = NrOfThroughSides + 1
                  LocSidesTemp(    NrOfThroughSides) = NblocSideID
                  TriNumTemp(      NrOfThroughSides) = TriNum
                  GlobSideTemp(    NrOfThroughSides) = nbSideID
                  isMortarSideTemp(NrOfThroughSides) = .TRUE.
                  SideID    = nbSideID
                  LocalSide = NblocSideID
                END IF
              END DO
            END DO
          END IF  ! Regular or mortar side
        END DO  ! iLocSide=1,6

        TriNum = TriNumTemp(1)

        ! Additional treatment if particle did not cross any sides or it crossed multiple sides
        IF (NrOfThroughSides.NE.1) THEN
! -- 2.c If no sides are found, the run is aborted (which should not happen). If multiple possible sides are found, additional
!        treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
!        the determinants
          IF (NrOfThroughSides.EQ.0) THEN
            ! Particle appears to have not crossed any of the checked sides. Abort!
            IPWRITE(UNIT_stdOut,'(I0,A,I0,A,I0,A,I0,A)') ' Error in Particle TriaTracking! Particle Number',iPart ,&
                                                                                           'lost. Element:',ElemID,&
                                                                                           '(species:'     ,PartSpecies(iPart),')'
            IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' LastPos: ', LastPartPos(1:3,iPart)
            IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' Pos:     ', PartState(  1:3,iPart)
            IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' Velo:    ', PartState(  4:6,iPart)
            CALL Abort(__STAMP__,'Lost particle detected during TriaTracking!')

          ELSE IF (NrOfThroughSides.GT.1) THEN
            ! Use the slower search method if particle appears to have crossed more than one side (possible for irregular hexagons
            ! and in the case of mortar elements)
            SecondNrOfThroughSides = 0
            minRatio               = 0
            oldElemIsMortar        = .FALSE.

            DO ind2 = 1, NrOfThroughSides
              doCheckSide = .TRUE.

              ! Check if this side was already treated
              DO indSide = 2, 6
                IF((DoneLastElem(1,indSide).EQ.ElemID)             .AND. &
                   (DoneLastElem(4,indSide).EQ.GlobSideTemp(ind2)) .AND. &
                   (DoneLastElem(3,indSide).EQ.TriNumTemp(ind2)))        &
                  doCheckSide = .FALSE.
              END DO

              IF (doCheckSide) THEN
                ! Side is a mortar side
                IF (isMortarSideTemp(ind2)) THEN
                  ! Get the element number of the smaller neighboring element
                  NbElemID = SideInfo_Shared(SIDE_ELEMID,GlobSideTemp(ind2))
                  ! Get the determinant between the old and new particle position and the nodes of the triangle which was crossed
                  CALL ParticleThroughSideLastPosCheck(iPart                  &
                                                      ,LocSidesTemp(ind2)     &
                                                      ,NbElemID               &
                                                      ,InElementCheck         &
                                                      ,TriNumTemp(ind2)       &
                                                      ,detM                   &
                                                      ,isMortarSide=.TRUE.    &
                                                      ,detPartPos  =detPartPos)

                  ! If the particle is inside the neighboring mortar element, it moved through this side
                  IF (InElementCheck) THEN
                    ! Particle moves within side
                    IF ((detM.EQ.0).AND.(detPartPos.EQ.0)) CYCLE

                    ! Determining which side was crossed first by comparing the ratio of the spatial product of PartPos->Tri-Nodes
                    ! and LastPartPos->Tri-Nodes
                    IF ((detM.EQ.0).AND.(minRatio.EQ.0)) THEN
                      SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                      SideID          = GlobSideTemp(ind2)
                      LocalSide       = LocSidesTemp(ind2)
                      TriNum          = TriNumTemp(  ind2)
                      oldElemIsMortar = .TRUE.
                    ! For the extremely unlikely case that the particle landed exactly on the side
                    ELSE
                      IF (detM.EQ.0.0) CYCLE

                      ratio = detPartPos/detM
                      ! Ratio is always negative since detM(=detLastPartPos) is negative or zero, i.e. maximum abs is wanted
                      ! The closer the intersected side is to the last particle position the greater the absolute ratio will be
                      IF (ratio.LT.minRatio) THEN
                        SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                        minRatio        = ratio
                        SideID          = GlobSideTemp(ind2)
                        LocalSide       = LocSidesTemp(ind2)
                        TriNum          = TriNumTemp(  ind2)
                        oldElemIsMortar = .TRUE.
                      END IF
                    END IF
                  END IF  ! InElementCheck
                ELSE  ! Regular side
                  CALL ParticleThroughSideLastPosCheck(iPart              &
                                                      ,LocSidesTemp(ind2) &
                                                      ,ElemID             &
                                                      ,InElementCheck     &
                                                      ,TriNumTemp(ind2)   &
                                                      ,detM)

                  IF (InElementCheck) THEN
                    ! particle moves within side
                    IF ((detM.EQ.0).AND.(det(LocSidesTemp(ind2),TriNumTemp(ind2)).EQ.0)) CYCLE

                    ! Determining which side was crossed first by comparing the ratio of the spatial product of PartPos->Tri-Nodes
                    ! and LastPartPos->Tri-Nodes
                    IF((detM.EQ.0).AND.(minRatio.EQ.0))THEN
                      SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                      SideID          = GlobSideTemp(ind2)
                      LocalSide       = LocSidesTemp(ind2)
                      TriNum          = TriNumTemp(  ind2)
                      oldElemIsMortar = .FALSE.
                    ELSE
                      ! For the extremely unlikely case that the particle landed exactly on the side
                      IF (detM.EQ.0) CYCLE

                      ratio = det(LocSidesTemp(ind2),TriNumTemp(ind2))/detM
                      ! Ratio is always negative since detM(=detLastPartPos) is negative or zero, i.e. maximum abs is wanted
                      ! The closer the intersected side is to the last particle position the greater the absolute ratio will be
                      IF (ratio.LT.minRatio) THEN
                        SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                        minRatio        = ratio
                        SideID          = GlobSideTemp(ind2)
                        LocalSide       = LocSidesTemp(ind2)
                        TriNum          = TriNumTemp(  ind2)
                        oldElemIsMortar = .FALSE.
                      END IF
                    END IF
                  END IF  ! InElementCheck
                END IF  ! isMortarSideTemp = T/F
              END IF  ! doCheckSide
            END DO  ! ind2 = 1, NrOfThroughSides

            ! Particle that went through multiple sides first, but did not cross any sides during the second check -> Deleted!
            IF (SecondNrOfThroughSides.EQ.0) THEN
              ! Particle appears to have not crossed any of the checked sides. Deleted!
              IPWRITE(UNIT_stdOut,'(I0,A,I0,A,I0,A,I0,A)') ' Error in Particle TriaTracking! Particle Number',iPart ,&
                                                                                             'lost. Element:',ElemID,&
                                                                                             '(species:'     ,PartSpecies(iPart),')'
              IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' LastPos: ', LastPartPos(1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' Pos:     ', PartState(  1:3,iPart)
              IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' Velo:    ', PartState(  4:6,iPart)
              CALL Abort(__STAMP__,'Lost particle detected during TriaTracking!')
            END IF
          END IF  ! NrOfThroughSides.EQ.0/.GT.1
        END IF  ! NrOfThroughSides.NE.1

! -- 3. In case of a boundary, perform the appropriate boundary interaction
        crossedBC = .FALSE.
        flip      = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

        ! Side is a boundary side
        IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
          OldElemID = ElemID
          BCType    = PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,SideID))

          IF (BCType.NE.1) CALL IntersectionWithWall(PartTrajectory &
                                                    ,alpha          &
                                                    ,iPart          &
                                                    ,LocalSide      &
                                                    ,ElemID         &
                                                    ,TriNum)

          CALL GetBoundaryInteraction(PartTrajectory       &
                                     ,lengthPartTrajectory &
                                     ,alpha                &
                                     ,xi,eta               &
                                     ,iPart                &
                                     ,SideID               &
                                     ,flip                 &
                                     ,ElemID               &
                                     ,crossedBC            &
                                     ,TriNum=TriNum)

          IF (.NOT.PDM%ParticleInside(iPart)) PartisDone = .TRUE.

#if USE_LOADBALANCE
          IF (OldElemID.GE.offsetElem+1 .AND. OldElemID.LE.offsetElem+PP_nElems) &
            CALL LBElemSplitTime(OldElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/

          SELECT CASE(BCType)
            ! PartBound%ReflectiveBC,PartBound%SymmetryBC
            CASE(2,10)
              DoneLastElem(:,:) = 0
            CASE DEFAULT
              DO ind2 = 5,1,-1
                DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
              END DO
              DoneLastElem(1,1) = OldElemID
              DoneLastElem(2,1) = LocalSide
              DoneLastElem(3,1) = TriNum
              DoneLastElem(4,1) = SideID
          END SELECT

        ! Side is an internal side
        ELSE
          DO ind2 = 5,1,-1
            DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
          END DO
          DoneLastElem(1,1) = ElemID
          DoneLastElem(2,1) = LocalSide
          DoneLastElem(3,1) = TriNum
          DoneLastElem(4,1) = SideID
          IF (oldElemIsMortar) THEN
            ElemID = SideInfo_Shared(SIDE_ELEMID  ,SideID)
          ELSE
            ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
          END IF
        END IF  ! SideInfo_Shared(SIDE_BCID,SideID).GT./.LE. 0
        CNElemID = GetCNElemID(ElemID)

        IF (CNElemID.LT.1) THEN
          IPWRITE(UNIT_stdOut,'(I0,A,F12.6)') 'Particle Velocity: ',VECNORM(PartState(4:6,iPart))
          CALL Abort(__STAMP__,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
        END IF

        TrackInfo%CurrElem = ElemID
      END IF ! InElementCheck = T/F
    END DO  ! .NOT.PartisDone

#if USE_LOADBALANCE
    IF (PEM%Element(iPart).GE.offsetElem+1 .AND. PEM%Element(iPart).LE.offsetElem+PP_nElems) &
      CALL LBElemPauseTime(PEM%Element(iPart)-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/

  END IF
END DO ! iPart = 1,PDM%ParticleVecLength

END SUBROUTINE ParticleTriaTracking


SUBROUTINE ParticleThroughSideCheck3DFast(PartID,PartTrajectory,iLocSide,Element,ThroughSide,TriNum,IsMortar)
!===================================================================================================================================
!> Routine to check whether a particle crossed the given triangle of a side. The determinant between the normalized trajectory
!> vector and the vectors from two of the three nodes to the old particle position is calculated. If the determinants for the three
!> possible combinations are greater than zero, then the particle went through this triangle of the side.
!> Note that if this is a mortar side, the side of the small neighbouring mortar element has to be checked. Thus, the orientation
!> is reversed.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Tools,         ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
INTEGER,INTENT(IN)               :: iLocSide
INTEGER,INTENT(IN)               :: Element
INTEGER,INTENT(IN)               :: TriNum
REAL,   INTENT(IN)               :: PartTrajectory(1:3)
LOGICAL, INTENT(IN), OPTIONAL    :: IsMortar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: ThroughSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER                   :: eps = 0.
INTEGER                          :: CNElemID
INTEGER                          :: n, NodeID
REAL                             :: Px, Py, Pz
REAL                             :: Vx, Vy, Vz
REAL                             :: xNode(3), yNode(3), zNode(3), Ax(3), Ay(3), Az(3)
REAL                             :: det(3)
!===================================================================================================================================

CNElemID = GetCNElemID(Element)

ThroughSide = .FALSE.

Px = lastPartPos(1,PartID)
Py = lastPartPos(2,PartID)
Pz = lastPartPos(3,PartID)

! Normalized particle trajectory (PartPos - lastPartPos)/ABS(PartPos - lastPartPos)
Vx = PartTrajectory(1)
Vy = PartTrajectory(2)
Vz = PartTrajectory(3)

! Get the coordinates of the first node and the vector from the particle position to the node
xNode(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
yNode(1) = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
zNode(1) = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)

Ax(1) = xNode(1) - Px
Ay(1) = yNode(1) - Py
Az(1) = zNode(1) - Pz
! Get the vectors to the other two nodes, depending on the triangle number
IF (PRESENT(IsMortar)) THEN
  ! Note: reverse orientation in the mortar case, as the side is treated from the perspective of the smaller neighbouring element
  !       (TriNum=1: NodeID=3,2; TriNum=2: NodeID=4,3)
  xNode(2) = NodeCoords_Shared(1,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  yNode(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  zNode(2) = NodeCoords_Shared(3,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)

  Ax(2) = xNode(2) - Px
  Ay(2) = yNode(2) - Py
  Az(2) = zNode(2) - Pz

  xNode(3) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
  yNode(3) = NodeCoords_Shared(2,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
  zNode(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)

  Ax(3) = xNode(3) - Px
  Ay(3) = yNode(3) - Py
  Az(3) = zNode(3) - Pz
ELSE
  DO n = 2,3
    NodeID = n+TriNum-1       ! m = true node number of the sides (TriNum=1: NodeID=2,3; TriNum=2: NodeID=3,4)
    xNode(n) = NodeCoords_Shared(1,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
    yNode(n) = NodeCoords_Shared(2,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)
    zNode(n) = NodeCoords_Shared(3,ElemSideNodeID_Shared(NodeID,iLocSide,CNElemID)+1)

    Ax(n) = xNode(n) - Px
    Ay(n) = yNode(n) - Py
    Az(n) = zNode(n) - Pz
  END DO
END IF
!--- check whether v and the vectors from the particle to the two edge nodes build
!--- a right-hand-system. If yes for all edges: vector goes potentially through side
det(1) = ((Ay(1) * Vz - Az(1) * Vy) * Ax(3)  + &
          (Az(1) * Vx - Ax(1) * Vz) * Ay(3)  + &
          (Ax(1) * Vy - Ay(1) * Vx) * Az(3))

det(2) = ((Ay(2) * Vz - Az(2) * Vy) * Ax(1)  + &
          (Az(2) * Vx - Ax(2) * Vz) * Ay(1)  + &
          (Ax(2) * Vy - Ay(2) * Vx) * Az(1))

det(3) = ((Ay(3) * Vz - Az(3) * Vy) * Ax(2)  + &
          (Az(3) * Vx - Ax(3) * Vz) * Ay(2)  + &
          (Ax(3) * Vy - Ay(3) * Vx) * Az(2))

! Comparison of the determinants with against machine precision (eps)
IF ((det(1).GE.-eps).AND.(det(2).GE.-eps).AND.(det(3).GE.-eps)) THEN
 ThroughSide = .TRUE.
END IF

RETURN

END SUBROUTINE ParticleThroughSideCheck3DFast


SUBROUTINE ParticleThroughSideLastPosCheck(i,iLocSide,Element,InElementCheck,TriNum,det,isMortarSide,detPartPos)
!===================================================================================================================================
!> Routine used in the case of a particle crossing multipe sides. Calculates the determinant of the three vectors from the last
!> (and current in the case of mortar sides) particle position to the nodes of the triangle in order to determine if the particle
!> is inside the element. Output of determinant is used to determine which of the sides was crossed first.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Tools,         ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: i, Element, iLocSide, TriNum
LOGICAL, INTENT(IN),OPTIONAL     :: isMortarSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: InElementCheck
REAL   ,INTENT(OUT)              :: det
REAL   ,INTENT(OUT), OPTIONAL    :: detPartPos
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: CNElemID
INTEGER                          :: NodeNum, ind, iNode
REAL                             :: Ax(3),Ay(3),Az(3)
REAL                             :: NodeCoord(1:3,1:3)
!===================================================================================================================================

CNElemID       = GetCNElemID(Element)
InElementCheck = .TRUE.

!--- coords of first node:
DO ind = 1,3
  NodeCoord(ind,1) = NodeCoords_Shared(ind,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
END DO

!--- coords of other two nodes (depending on triangle):
IF(PRESENT(isMortarSide)) THEN
  ! Note: reversed orientation as the triangle is treated from the perspective of the smaller neighbouring mortar element
  NodeCoord(1:3,2) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(2+TriNum,iLocSide,CNElemID)+1)
  NodeCoord(1:3,3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(1+TriNum,iLocSide,CNElemID)+1)
ELSE
  DO iNode = 2,3
    NodeNum = iNode + TriNum - 1
    DO ind = 1,3
      NodeCoord(ind,iNode) = NodeCoords_Shared(ind,ElemSideNodeID_Shared(NodeNum,iLocSide,CNElemID)+1)
    END DO
  END DO
END IF

!--- vector from lastPos(!) to triangle nodes
DO ind = 1,3
  Ax(ind) = NodeCoord(1,ind) - lastPartPos(1,i)
  Ay(ind) = NodeCoord(2,ind) - lastPartPos(2,i)
  Az(ind) = NodeCoord(3,ind) - lastPartPos(3,i)
END DO

!--- determine whether particle is on inner side (rel. to element) of triangle
!--- set corresponding "flag" (see below)
det = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
       (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
       (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))

IF ((det.lt.0).OR.(det.NE.det)) THEN
  InElementCheck = .FALSE.
END IF

IF(PRESENT(isMortarSide).AND.PRESENT(detPartPos)) THEN
  DO ind = 1,3
    Ax(ind) = NodeCoord(1,ind) - PartState(1,i)
    Ay(ind) = NodeCoord(2,ind) - PartState(2,i)
    Az(ind) = NodeCoord(3,ind) - PartState(3,i)
  END DO

  detPartPos = ((Ay(1) * Az(2) - Az(1) * Ay(2)) * Ax(3) +     &
                (Az(1) * Ax(2) - Ax(1) * Az(2)) * Ay(3) +     &
                (Ax(1) * Ay(2) - Ay(1) * Ax(2)) * Az(3))
END IF

END SUBROUTINE ParticleThroughSideLastPosCheck

END MODULE MOD_Particle_TriaTracking
