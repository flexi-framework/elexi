!=================================================================================================================================
! Copyright (c) 2010-2019  Prof. Claus-Dieter Munz
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
! Contains routines to locate the particle on the mesh
!===================================================================================================================================
MODULE MOD_Particle_Localization
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE SinglePointToElement
  MODULE PROCEDURE SinglePointToElement
END INTERFACE

INTERFACE LocateParticleInElement
  MODULE PROCEDURE LocateParticleInElement
END INTERFACE

INTERFACE PartInElemCheck
  MODULE PROCEDURE PartInElemCheck
END INTERFACE

INTERFACE ParticleInsideQuad3D
  MODULE PROCEDURE ParticleInsideQuad3D
END INTERFACE

INTERFACE CountPartsPerElem
  MODULE PROCEDURE CountPartsPerElem
END INTERFACE

PUBLIC:: SinglePointToElement
PUBLIC:: LocateParticleInElement
PUBLIC:: PartInElemCheck
PUBLIC:: ParticleInsideQuad3D
PUBLIC:: CountPartsPerElem
!===================================================================================================================================

CONTAINS


SUBROUTINE LocateParticleInElement(PartID,doHALO)
!----------------------------------------------------------------------------------------------------------------------------------!
! Finds a single particle in its host element
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: PartID
LOGICAL,INTENT(IN)                :: doHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ElemID
!===================================================================================================================================
ElemID = SinglePointToElement(PartState(1:3,PartID),doHALO=doHALO)
PEM%Element(PartID) = ElemID
IF(ElemID.EQ.-1)THEN
  PDM%ParticleInside(PartID)=.FALSE.
ELSE
  PDM%ParticleInside(PartID)=.TRUE.
  IF(TrackingMethod.EQ.REFMAPPING)THEN
    CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef(1:3,PartID),ElemID)
  END IF ! TrackingMethod.EQ.REFMAPPING
END IF ! ElemID.EQ.-1
END SUBROUTINE LocateParticleInElement


!===================================================================================================================================
!> This function maps a 3D point to an element
!> returns elementID or -1 in case no element was found
!===================================================================================================================================
INTEGER FUNCTION SinglePointToElement(Pos3D,doHALO)
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars     ,ONLY: Geo
USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGM_nElems,FIBGM_offsetElem,FIBGM_Element
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_Particle_Tracking_Vars ,ONLY: Distance,ListDistance,TrackingMethod
USE MOD_Particle_Utils         ,ONLY: InsertionSort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)                   :: Pos3D(1:3)
LOGICAL,INTENT(IN)                :: doHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems,ElemID,CNElemID, iBGM,jBGM,kBGM
REAL                              :: Distance2, RefPos(1:3)
REAL                              :: Det(6,2)
LOGICAL                           :: InElementCheck
!===================================================================================================================================

SinglePointToElement = -1

! --- get background mesh cell of point
iBGM = CEILING((Pos3D(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
iBGM = MAX(MIN(GEO%FIBGMimax,iBGM),GEO%FIBGMimin)
jBGM = CEILING((Pos3D(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
jBGM = MAX(MIN(GEO%FIBGMjmax,jBGM),GEO%FIBGMjmin)
kBGM = CEILING((Pos3D(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
kBGM = MAX(MIN(GEO%FIBGMkmax,kBGM),GEO%FIBGMkmin)

!--- check all cells associated with this background mesh cell
nBGMElems = FIBGM_nElems(iBGM,jBGM,kBGM)

! get closest element barycenter
Distance=-1.

ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iBGMElem)
  CNElemID = GetCNElemID(ElemID)
 Distance2 = SUM((Pos3D(1:3)-ElemBaryNGeo(1:3,CNElemID))**2.)

  ! element in range
  Distance(iBGMElem)     = MERGE(Distance2,-1.,Distance2.LE.ElemRadius2NGeo(CNElemID))
  ListDistance(iBGMElem) = ElemID
END DO ! nBGMElems

IF (ALMOSTEQUAL(MAXVAL(Distance),-1.)) THEN
  RETURN
END IF

IF(nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

! loop through sorted list and start by closest element
InElementCheck = .FALSE.

DO iBGMElem = 1,nBGMElems
  ! Element is out of range
  IF (ALMOSTEQUAL(Distance(iBGMElem),-1.)) CYCLE

  ElemID = ListDistance(iBGMElem)

  IF (.NOT.DoHALO) THEN
    IF (ElemID.LT.offsetElem+1 .OR. ElemID.GT.offsetElem+PP_nElems) CYCLE
  END IF

  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING)
      CALL ParticleInsideQuad3D(Pos3D(1:3),ElemID,InElementCheck,Det)

    ! CASE(TRACING)
    !   CALL GetPositionInRefElem(Pos3D(1:3),RefPos,ElemID)
    !   IF (MAXVAL(ABS(RefPos)).LE.1.0) InElementCheck = .TRUE.

    ! FLEXI has ElemEpsOneCell also with TRACING
    ! CASE(REFMAPPING)
    CASE(TRACING,REFMAPPING)
      CALL GetPositionInRefElem(Pos3D(1:3),RefPos,ElemID)
      IF (MAXVAL(ABS(RefPos)).LE.1.0) InElementCheck=.TRUE.

  END SELECT

  IF (InElementCheck) THEN
    SinglePointToElement = ElemID
    RETURN
  END IF

END DO ! iBGMElem

END FUNCTION SinglePointToElement


SUBROUTINE PartInElemCheck(PartPos_In,PartID,ElemID,FoundInElem,IntersectPoint_Opt, &
#if CODE_ANALYZE
                           Sanity_Opt,Tol_Opt,CodeAnalyze_Opt)
#else
                           Tol_Opt)
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
! Checks if particle is in Element
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals       ,ONLY: ALMOSTZERO,VECNORM
USE MOD_Particle_Intersection  ,ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection  ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection  ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Intersection  ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Mesh_Tools    ,ONLY: GetCNElemID,GetCNSideID,GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars ,ONLY: SideType,SideNormVec
USE MOD_Particle_Vars          ,ONLY: LastPartPos
#if CODE_ANALYZE
USE MOD_Globals                ,ONLY: MyRank,UNIT_stdout
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
USE MOD_Particle_Surfaces      ,ONLY: OutputBezierControlPoints
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Intersection  ,ONLY: OutputTrajectory
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: ElemID,PartID
REAL,INTENT(IN)                          :: PartPos_In(1:3)
#if CODE_ANALYZE
LOGICAL,INTENT(IN),OPTIONAL              :: CodeAnalyze_Opt
LOGICAL,INTENT(IN),OPTIONAL              :: Sanity_Opt
#endif /*CODE_ANALYZE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                      :: FoundInElem
REAL,INTENT(OUT),OPTIONAL                :: IntersectPoint_Opt(1:3)
REAL,INTENT(OUT),OPTIONAL                :: Tol_Opt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: CNElemID
INTEGER                                  :: ilocSide,flip,SideID,CNSideID
REAL                                     :: PartTrajectory(1:3),NormVec(1:3)
REAL                                     :: lengthPartTrajectory,PartPos(1:3),LastPosTmp(1:3)
LOGICAL                                  :: isHit
REAL                                     :: alpha,eta,xi,IntersectPoint(1:3)
!===================================================================================================================================

IF(PRESENT(tol_Opt)) tol_Opt=-1.

CNElemID = GetCNElemID(ElemID)

! virtual move to element barycenter
LastPosTmp(1:3)         = LastPartPos(1:3,PartID)
LastPartPos(1:3,PartID) = ElemBaryNGeo(1:3,CNElemID)
PartPos(1:3)            = PartPos_In(1:3)

! get trajectory from element barycenter to current position
PartTrajectory          = PartPos - LastPartPos(1:3,PartID)
lengthPartTrajectory    = VECNORM(PartTrajectory(1:3))

! output the part trajectory
#if CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRankOut.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      IPWRITE(UNIT_stdout,*) ' --------------------------------------------- '
      IPWRITE(UNIT_stdout,*) ' PartInElemCheck '
      CALL OutputTrajectory(PartID,PartPos,PartTrajectory,lengthPartTrajectory)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! we found the particle on the element barycenter
IF(ALMOSTZERO(lengthPartTrajectory))THEN
  FoundInElem =.TRUE.
  LastPartPos(1:3,PartID) = LastPosTmp(1:3)
  RETURN
END IF

! normalize the part trajectory vector
PartTrajectory=PartTrajectory/lengthPartTrajectory

! reset intersection counter and alpha
isHit = .FALSE.
alpha = -1.

DO ilocSide = 1,6
  SideID   = GetGlobalNonUniqueSideID(ElemID,iLocSide)
  CNSideID = GetCNSideID(SideID)
  flip     = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT)
      CALL ComputePlanarRectIntersection(  ishit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,flip,SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,flip,SideID)
    CASE(BILINEAR,PLANAR_NONRECT)
        CALL ComputeBiLinearIntersection(  isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID     ,SideID,ElemCheck_Opt=.TRUE.)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(      isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID     ,SideID,ElemCheck_Opt=.TRUE.)
  END SELECT

#if CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(15("="))')
      WRITE(UNIT_stdout,'(A)')           '     | Output after compute intersection (PartInElemCheck): '
      WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(CNSideID)  ,' | SideID: ',SideID,'| Hit: ',isHit
      WRITE(UNIT_stdout,'(2(A,G0))')     '     | LengthPT: ',LengthPartTrajectory,' | Alpha: ',Alpha
      WRITE(UNIT_stdout,'(A,2(X,G0))')   '     | Intersection xi/eta: ',xi,eta
    END IF
  END IF

  IF(PRESENT(Sanity_Opt))THEN
    IF(Sanity_Opt)THEN
      IF(alpha.GT.-1)THEN
        ! alpha is going from barycenter to point
        ! here, the tolerance for the ratio alpha/LengthPartTrajectory for tracing with element-corners is determined.
        IF(PRESENT(tol_Opt)) tol_Opt=MAX(ABS(1.-alpha/LengthPartTrajectory),tol_Opt)
        ! mark element as trouble element if rel. tol from alpha/LengthPartTrajectory to 1 > 1e-4
        ! tolerance 1e-4 is from ANSA_BOX grid (experimental, arbitrary)
        IF(ALMOSTEQUALRELATIVE(alpha/LengthPartTrajectory,1.,0.002)) THEN
          alpha=-1
        ELSE
          print*,'alpha',alpha,LengthPartTrajectory,alpha/LengthPartTrajectory,ABS(1.-alpha/LengthPartTrajectory),tol_Opt
        END IF
      END IF
    END IF
  END IF

  ! Dirty fix for PartInElemCheck if LastPartPos is almost on side (tolerance issues)
  IF(PRESENT(CodeAnalyze_Opt))THEN
    IF(CodeAnalyze_Opt)THEN
      IF((alpha)/LengthPartTrajectory.GT.0.9)THEN
        alpha = -1.0
      END IF
    END IF
  END IF
#endif /*CODE_ANALYZE*/

  ! found an intersection with a side, but we might be moving into this element from the boundary. Get side normal vector and
  ! intersection point to check.
  IF (alpha.GT.-1) THEN
    SELECT CASE(SideType(CNSideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        CNSideID = GetCNSideID(SideID)
        NormVec  = SideNormVec(1:3,CNSideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    IF(flip.NE.0) NormVec=-NormVec
    IntersectPoint=LastPartPos(1:3,PartID)+alpha*PartTrajectory

#if CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,*) '     | alpha          ',alpha
      WRITE(UNIT_stdout,*) '     | Normal vector  ',NormVec
      WRITE(UNIT_stdout,*) '     | PartTrajectory ',PartTrajectory
      WRITE(UNIT_stdout,*) '     | Dotprod        ',DOT_PRODUCT(NormVec,PartTrajectory)
      WRITE(UNIT_stdout,*) '     | Point 2        ', LastPartPos(1:3,PartID)+alpha*PartTrajectory+NormVec
      WRITE(UNIT_stdout,*) '     | Beziercontrolpoints3d-x'
      CALL OutputBezierControlPoints(BezierControlPoints3D_in=BezierControlPoints3D(1:3,:,:,SideID))
    END IF
  END IF
#endif /*CODE_ANALYZE*/

    ! check if we are moving into this element from the cell boundary. If yes, reset alpha for this side.
    IF(DOT_PRODUCT(NormVec,PartTrajectory).LT.0.)THEN
      alpha=-1.0
    ELSE
      EXIT
    END IF

  END IF
END DO ! ilocSide

! write final decision. We are inside if we did not find an intersection or are moving into this element from the cell boundary
FoundInElem = .TRUE.
IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt = 0.

IF (alpha.GT.-1) THEN
  FoundInElem = .FALSE.
  IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt = IntersectPoint
END IF

! reset the LastParPos to its original value
LastPartPos(1:3,PartID) = LastPosTmp(1:3)

END SUBROUTINE PartInElemCheck


SUBROUTINE ParticleInsideQuad3D(PartStateLoc,ElemID,InElementCheck,Det)
!===================================================================================================================================
!> Checks if particle is inside of a linear element with triangulated faces, compatible with mortars
!> Regular element: The determinant of a 3x3 matrix, where the three vectors point from the particle to the nodes of a triangle, is
!>                  used to determine whether the particle is inside the element. The geometric equivalent is the triple product
!>                  A*(B x C), spanning a signed volume. If the volume/determinant is positive, then the particle is inside.
!> Element with neighbouring mortar elements: Additional checks of the smaller sides are required if the particle is in not in the
!>                                       concave part of the element but in the convex. Analogous procedure using the determinants.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Tools   ,ONLY: GetCNElemID,GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars    ,ONLY: ConcaveElemSide_Shared,ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   ,INTENT(OUT)           :: Det(6,2)
LOGICAL,INTENT(OUT)           :: InElementCheck
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: NbElemID,CNElemID
INTEGER                       :: ilocSide, NodeNum, SideID, SideIDMortar, ind, nNbMortars
LOGICAL                       :: PosCheck, NegCheck, InElementCheckMortar, InElementCheckMortarNb
REAL                          :: A(1:3,1:4), crossP(3)
!===================================================================================================================================

InElementCheck = .TRUE.
InElementCheckMortar = .TRUE.

CNElemID  = GetCNElemID(ElemID)

! loop over all sides attached to an element
DO iLocSide = 1,6
  SideID = GetGlobalNonUniqueSideID(ElemID,iLocSide)

  DO NodeNum = 1,4
    !--- A = vector from particle to node coords
    A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,iLocSide,CNElemID)+1) - PartStateLoc(1:3)
  END DO

  NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
  !--- Treatment of sides which are adjacent to mortar elements
  IF (NbElemID.LT.0) THEN
    PosCheck = .FALSE.
    NegCheck = .FALSE.

    !--- Checking the concave part of the side
    IF (ConcaveElemSide_Shared(iLocSide,CNElemID)) THEN
      ! If the element is actually concave, CalcDetOfTrias determines its determinants
      Det(iLocSide,1:2) = CalcDetOfTrias(A,1)
      IF (Det(iLocSide,1).GE.0) PosCheck = .TRUE.
      IF (Det(iLocSide,2).GE.0) PosCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (.NOT.PosCheck) InElementCheckMortar = .FALSE.
    ELSE
      ! If its a convex element, CalcDetOfTrias determines the concave determinants
      Det(iLocSide,1:2) = CalcDetOfTrias(A,2)
      IF (Det(iLocSide,1).GE.0) PosCheck = .TRUE.
      IF (Det(iLocSide,2).GE.0) PosCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (.NOT.PosCheck) InElementCheckMortar= .FALSE.
    END IF

    !--- Checking the convex part of the side
    IF (.NOT.InElementCheckMortar) THEN
      InElementCheckMortar = .TRUE.
      IF (ConcaveElemSide_Shared(iLocSide,CNElemID)) THEN
        Det(iLocSide,1:2) = CalcDetOfTrias(A,2)
        IF (Det(iLocSide,1).LT.0) NegCheck = .TRUE.
        IF (Det(iLocSide,2).LT.0) NegCheck = .TRUE.
        !--- final determination whether particle is in element
        IF (NegCheck) InElementCheckMortar = .FALSE.
      ELSE
        Det(iLocSide,1:2) = CalcDetOfTrias(A,1)
        IF (Det(iLocSide,1).LT.0) NegCheck = .TRUE.
        IF (Det(iLocSide,2).LT.0) NegCheck = .TRUE.
        !--- final determination whether particle is in element
        IF (NegCheck) InElementCheckMortar= .FALSE.
      END IF
      !--- Particle is in a convex elem but not in concave, checking additionally the mortar neighbors. If particle is not inside
      !    the mortar elements, it has to be in the original element.
      IF (InElementCheckMortar) THEN
        nNbMortars = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)
        DO ind = 1, nNbMortars
          InElementCheckMortarNb = .TRUE.
          SideIDMortar = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
          NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideIDMortar)

          ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
          ! no way to recover it during runtime
          IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)

          CALL ParticleInsideNbMortar(PartStateLoc,NbElemID,InElementCheckMortarNb)
          IF (InElementCheckMortarNb) THEN
            InElementCheck = .FALSE.
            EXIT
          END IF
        END DO
      ELSE
        InElementCheck = .FALSE.
      END IF
    END IF

  ! Treatment of regular elements without mortars
  ELSE
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- compute cross product for vector 1 and 3
    crossP(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
    crossP(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
    crossP(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    !--- negative determinant of triangle 1 (points 1,3,2):
    Det(iLocSide,1) = crossP(1) * A(1,2) + &
                      crossP(2) * A(2,2) + &
                      crossP(3) * A(3,2)
    Det(iLocSide,1) = -det(iLocSide,1)
    !--- determinant of triangle 2 (points 1,3,4):
    Det(iLocSide,2) = crossP(1) * A(1,4) + &
                      crossP(2) * A(2,4) + &
                      crossP(3) * A(3,4)

    IF (Det(iLocSide,1).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF

    IF (Det(iLocSide,2).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF

    !--- final determination whether particle is in element
    IF (ConcaveElemSide_Shared(iLocSide,CNElemID)) THEN
      IF (.NOT.PosCheck) InElementCheck = .FALSE.
    ELSE
      IF (NegCheck) InElementCheck = .FALSE.
    END IF

  END IF ! Mortar element or regular element
END DO ! iLocSide = 1,6

RETURN

END SUBROUTINE ParticleInsideQuad3D


PURE SUBROUTINE ParticleInsideNbMortar(PartStateLoc,ElemID,InElementCheck)
!===================================================================================================================================
!> Routines checks if the particle is inside the neighbouring mortar element. Used for the regular ParticleInsideQuad3D routine
!> after it was determined that the particle is not in the concave part but in the convex part of the element.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Tools   ,ONLY: GetCNElemID,GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars    ,ONLY: SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars    ,ONLY: ConcaveElemSide_Shared,ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
REAL   ,INTENT(IN)            :: PartStateLoc(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: InElementCheck
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: NbElemID,CNElemID
INTEGER                       :: ilocSide, NodeNum, SideID
LOGICAL                       :: PosCheck, NegCheck
REAL                          :: A(1:3,1:4), cross(3)
REAL                          :: Det(2)
!===================================================================================================================================

InElementCheck = .TRUE.
CNElemID  = GetCNElemID(ElemID)

! loop over all sides attached to an element
DO iLocSide = 1,6
  SideID = GetGlobalNonUniqueSideID(ElemID,iLocSide)

  DO NodeNum = 1,4
  !--- A = vector from particle to node coords
    A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,iLocSide,CNElemID)+1) - PartStateLoc(1:3)
  END DO

  NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

  !--- Treatment of sides which are adjacent to mortar elements
  IF (NbElemID.LT.0) THEN
    !--- initialize flags for side checks
    PosCheck = .FALSE.
    NegCheck = .FALSE.

    !--- Check if the particle is inside the convex element. If its outside, it has to be inside the original element
    IF (ConcaveElemSide_Shared(iLocSide,CNElemID)) THEN
      Det(1:2) = CalcDetOfTrias(A,2)
      IF (Det(1).LT.0) NegCheck = .TRUE.
      IF (Det(2).LT.0) NegCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (NegCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    ELSE
      Det(1:2) = CalcDetOfTrias(A,1)
      IF (Det(1).LT.0) NegCheck = .TRUE.
      IF (Det(2).LT.0) NegCheck = .TRUE.
      !--- final determination whether particle is in element
      IF (NegCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    END IF

  ! Regular side
  ELSE
    PosCheck = .FALSE.
    NegCheck = .FALSE.
    !--- compute cross product for vector 1 and 3
    cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
    cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
    cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    !--- negative determinant of triangle 1 (points 1,3,2):
    Det(1) = cross(1) * A(1,2) + &
                      cross(2) * A(2,2) + &
                      cross(3) * A(3,2)
    Det(1) = -det(1)
    !--- determinant of triangle 2 (points 1,3,4):
    Det(2) = cross(1) * A(1,4) + &
                      cross(2) * A(2,4) + &
                      cross(3) * A(3,4)
    IF (Det(1).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF

    IF (Det(2).LT.0) THEN
      NegCheck = .TRUE.
    ELSE
      PosCheck = .TRUE.
    END IF

    !--- final determination whether particle is in element
    IF (ConcaveElemSide_Shared(iLocSide,CNElemID)) THEN
      IF (.NOT.PosCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    ELSE
      IF (NegCheck) THEN
        InElementCheck = .FALSE.
        RETURN
      END IF
    END IF

  END IF  ! Mortar or regular side
END DO  ! iLocSide = 1,6

END SUBROUTINE ParticleInsideNbMortar


SUBROUTINE CountPartsPerElem(ResetNumberOfParticles)
!===================================================================================================================================
! count number of particles in element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_LoadBalance_Vars,        ONLY: nPartsPerElem
USE MOD_Mesh_Vars,               ONLY: offsetElem
USE MOD_Particle_Globals
USE MOD_Particle_Vars,           ONLY: PDM,PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: ResetNumberOfParticles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPart, ElemID
!===================================================================================================================================

! DO NOT NULL this here, if e.g. this routine is called in between RK-stages in which particles are created
IF (ResetNumberOfParticles) THEN
  nPartsPerElem = 0
END IF

! loop over all particles and add them up
DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! LoadBalance operates on local ElemID
    ElemID = PEM%Element(iPart) - offsetElem

    ! Only consider elements currently on the same proc
    IF (ElemID.GE.1 .AND. ElemID.LE.PP_nElems) THEN
      nPartsPerElem(ElemID) = nPartsPerElem(ElemID) + 1
    END IF
  END IF
END DO ! iPart=1,PDM%ParticleVecLength

END SUBROUTINE CountPartsPerElem


PURE FUNCTION CalcDetOfTrias(A,bending)
!================================================================================================================================
!> Calculates the determinant A*(B x C) for both triangles of a side. bending = 1 gives the determinant considering the actual
!> orientation of the side (concave/convex), 2 gives the opposite of the saved form (e.g. a concave side gets the convex analog)
!================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                     :: A(3,4)
INTEGER,INTENT(IN)                  :: bending
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                :: CalcDetOfTrias(2)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: cross(3)
!================================================================================================================================
IF (bending.EQ.1) THEN
  !--- compute cross product for vector 1 and 3
  cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
  cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
  cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

  !--- negative determinant of triangle 1 (points 1,3,2):
  CalcDetOfTrias(1) = cross(1) * A(1,2) + &
                   cross(2) * A(2,2) + &
                   cross(3) * A(3,2)
  CalcDetOfTrias(1)  = -CalcDetOfTrias(1)
  !--- determinant of triangle 2 (points 1,3,4):
  CalcDetOfTrias(2)  = cross(1) * A(1,4) + &
                   cross(2) * A(2,4) + &
                   cross(3) * A(3,4)
ELSE
  !--- compute cross product for vector 2 and 4
  cross(1) = A(2,2) * A(3,4) - A(3,2) * A(2,4)
  cross(2) = A(3,2) * A(1,4) - A(1,2) * A(3,4)
  cross(3) = A(1,2) * A(2,4) - A(2,2) * A(1,4)

  !--- negative determinant of triangle 1 (points 2,4,1):
  CalcDetOfTrias(1) = cross(1) * A(1,1) + &
                   cross(2) * A(2,1) + &
                   cross(3) * A(3,1)
  !--- determinant of triangle 2 (points 2,4,3):
  CalcDetOfTrias(2) = cross(1) * A(1,3) + &
                   cross(2) * A(2,3) + &
                   cross(3) * A(3,3)
  CalcDetOfTrias(2) = -CalcDetOfTrias(2)
END IF

END FUNCTION CalcDetOfTrias

END MODULE MOD_Particle_Localization
