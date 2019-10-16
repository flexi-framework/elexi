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
! Contains routines to locate the particle on the mesh
!===================================================================================================================================
MODULE MOD_Particle_Localization
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
INTERFACE SingleParticleToExactElement
  MODULE PROCEDURE SingleParticleToExactElement
END INTERFACE

INTERFACE SingleParticleToExactElementNoMap
  MODULE PROCEDURE SingleParticleToExactElementNoMap
END INTERFACE

INTERFACE PartInElemCheck
  MODULE PROCEDURE PartInElemCheck
END INTERFACE

INTERFACE ParticleInsideQuad3D
  MODULE PROCEDURE ParticleInsideQuad3D
END INTERFACE

PUBLIC::SingleParticleToExactElement
PUBLIC::SingleParticleToExactElementNoMap
PUBLIC::PartInElemCheck
PUBLIC::ParticleInsideQuad3D
!===================================================================================================================================

CONTAINS


SUBROUTINE SingleParticleToExactElement(iPart,doHalo,initFix,doRelocate)
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,               ONLY:PartState,PEM,PDM,PartPosRef,KeepWallParticles
USE MOD_Particle_Mesh_Vars,          ONLY:Geo
USE MOD_Particle_Tracking_Vars,      ONLY:DoRefMapping,TriaTracking
USE MOD_Particle_Mesh_Vars,          ONLY:epsOneCell,IsTracingBCElem,ElemRadius2NGeo
USE MOD_Eval_xyz,                    ONLY:TensorProductInterpolation
USE MOD_Particle_Utils,              ONLY:InsertionSort
USE MOD_Particle_Tracking_Vars,      ONLY:DoRefMapping,Distance,ListDistance
USE MOD_Particle_Boundary_Condition, ONLY:PARTSWITCHELEMENT
USE MOD_Particle_MPI_Vars,           ONLY:SafetyFactor
#if USE_MPI
USE MOD_Mesh_Vars,                   ONLY:ElemToSide,BC
USE MOD_Particle_MPI_Vars,           ONLY:PartHaloElemToProc
#endif
USE MOD_Particle_Mesh_Vars,          ONLY:ElemBaryNGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: iPart
LOGICAL,INTENT(IN)                :: doHalo
LOGICAL,INTENT(IN)                :: initFix
LOGICAL,INTENT(IN)                :: doRelocate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems,ElemID,CellX,CellY,CellZ
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                           :: InElementCheck,ParticleFound
REAL                              :: xi(1:3),Distance2,Det(6,2)
#if USE_MPI
INTEGER                           :: XiDir,locSideID,flip,SideID
REAL                              :: locXi,locEta,tmpXi
INTEGER                           :: Moved(2)
#endif /*MPI*/
!===================================================================================================================================
ParticleFound = .FALSE.

! --- if the particle is at a wall, leave it in the last element
IF (KeepWallParticles) THEN
  IF (PDM%ParticleAtWall(iPart)) THEN
    PEM%Element(iPart) = PEM%lastElement(iPart)
    ParticleFound = .TRUE.
    RETURN
  END IF
END IF

! --- check if the particle is inside our bounding box (including/excluding the halo region)
IF(DoHALO)THEN
  IF ( (PartState(iPart,1).LT.GEO%xminglob).OR.(PartState(iPart,1).GT.GEO%xmaxglob).OR. &
       (PartState(iPart,2).LT.GEO%yminglob).OR.(PartState(iPart,2).GT.GEO%ymaxglob).OR. &
       (PartState(iPart,3).LT.GEO%zminglob).OR.(PartState(iPart,3).GT.GEO%zmaxglob)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
ELSE
  IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
       (PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
       (PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
END IF

! --- get background mesh cell of particle
CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
CellX = MAX(MIN(GEO%TFIBGMimax,CellX),GEO%TFIBGMimin)
CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MAX(MIN(GEO%TFIBGMjmax,CellY),GEO%TFIBGMjmin)
CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MAX(MIN(GEO%TFIBGMkmax,CellZ),GEO%TFIBGMkmin)


IF (TriaTracking) THEN
  !--- check all cells associated with this background mesh cell
  DO iBGMElem = 1, GEO%FIBGM(CellX,CellY,CellZ)%nElem
    ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
    CALL ParticleInsideQuad3D(PartState(iPart,1:3),ElemID,InElementCheck,Det)
    IF (InElementCheck) THEN
       PEM%Element(iPart) = ElemID
       ParticleFound = .TRUE.
       EXIT
    END IF
  END DO
  IF (.NOT.ParticleFound) THEN
    PDM%ParticleInside(iPart) = .FALSE.
  END IF
   RETURN
END IF

!--- check all cells associated with this background mesh cell
nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem

! get distance to background cell barycenter
Distance=-1.
ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  Distance2=(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
           +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
           +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))
  IF(Distance2.GT.ElemRadius2NGeo(ElemID))THEN
    Distance(iBGMElem)=-1.
  ELSE
    Distance(iBGMElem)=Distance2
  END IF
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

! no candidate in background mesh for our particle position. Can't determine if the current position is valid because we're missing
! the BC information as well. Abort and inform the user to increase the size of the halo mesh.
!>> WARNING: If we did an invalid push, we might abort here as well.
IF(ALMOSTEQUAL(MAXVAL(Distance),-1.))THEN
  PDM%ParticleInside(iPart) = .FALSE.
  IF(DoRelocate) CALL abort(&
__STAMP__&
  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
  RealInfo=SafetyFactor)
  RETURN
END IF

!  multiple cells associated with this background mesh cell. Sort according to distance to the cell barycenter
IF(nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

! loop through sorted list and start by closest element
DO iBGMElem=1,nBGMElems
  IF(ALMOSTEQUAL(Distance(iBGMElem),-1.))CYCLE
  ElemID=ListDistance(iBGMElem)

  ! exclude elements not within our local mesh (without the halo region)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF

  IF(IsTracingBCElem(ElemID))THEN
    CALL PartInElemCheck(PartState(iPart,1:3),iPart,ElemID,InElementCheck)
    IF(.NOT.InElementCheck) CYCLE
  END IF

  ! get position in reference element for element ElemID
  CALL TensorProductInterpolation(PartState(iPart,1:3),xi,ElemID)

  ! check if we are within tolerance for our chosen element
  IF(MAXVAL(ABS(Xi)).LT.epsOneCell(ElemID)) THEN
    IF(.NOT.InitFix) THEN
      InElementCheck=.TRUE.
    ELSE
     ! take extra care in the MPI case. We might encounter xi larger than unity but smaller then epsOneCell [1,epsOneCell], than the
     ! particle is found at least twice.
     InElementCheck=.TRUE.
     ! InElementCheck can only be set to false in the following part
#if USE_MPI
     ! check if xi is larger than unity, than the  particle is found at least twice
     ! particle close to the cell border and possible outside
     IF(MAXVAL(ABS(Xi)).GT.0.99999999) THEN
       ! find the direction closest to the neighbor side
       XiDir = MAXLOC(ABS(Xi),1)
       ! now, get neighbor-side id
       SELECT CASE(XiDir)
       CASE(1) ! Xi
         IF(Xi(XiDir).GT.0)THEN
           ! XI_PLUS
           locSideID=XI_PLUS
           locXi=Xi(3)
           locEta=Xi(2)
         ELSE
           ! XI_MINUS
           locSideID=XI_MINUS
           locXi=Xi(2)
           locEta=Xi(3)
         END IF
       CASE(2) ! Eta
         IF(Xi(XiDir).GT.0)THEN
           locSideID=ETA_PLUS
           locXi=-Xi(1)
           locEta=Xi(3)
         ELSE
           locSideID=ETA_MINUS
           locXi=Xi(1)
           locEta=Xi(3)
         END IF
       CASE(3) ! Zeta
         IF(Xi(XiDir).GT.0)THEN
           locSideID=ZETA_PLUS
           locXi =Xi(1)
           locEta=Xi(2)
         ELSE
           locSideID=ZETA_MINUS
           locXi=Xi(2)
           locEta=Xi(1)
         END IF
       CASE DEFAULT
         CALL abort(&
__STAMP__&
, ' Error in  mesh-connectivity!')
       END SELECT
       ! get flip and rotate xi and eta into side-master system
       flip     =ElemToSide(E2S_FLIP,locSideID,ElemID)
       SideID   =ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
       SELECT CASE(Flip)
       CASE(1) ! slave side, SideID=q,jSide=p
         tmpXi=locEta
         locEta=locXi
         locXi=tmpXi
       CASE(2) ! slave side, SideID=N-p,jSide=q
         locXi=-locXi
         locEta=locEta
       CASE(3) ! slave side, SideID=N-q,jSide=N-p
         tmpXi =-locEta
         locEta=-locXi
         locXi=tmpXi
       CASE(4) ! slave side, SideID=p,jSide=N-q
         locXi =locXi
         locEta=-locEta
       END SELECT
       IF(BC(SideID).GT.0)THEN
         InElementCheck=.FALSE.
       ELSE
        ! check if neighbor element is an mpi-element and if yes, only take the particle if I am the lower rank
         Moved = PARTSWITCHELEMENT(locxi,loceta,locSideID,SideID,ElemID)
         IF(Moved(1).GT.PP_nElems)THEN
           IF(PartHaloElemToProc(NATIVE_PROC_ID,Moved(1)).LT.MyRank)THEN
             InElementCheck=.FALSE.
           END IF
         END IF
       END IF
     END IF
#endif /*MPI*/
    END IF
  ! particle at face,edge or node, check most possible point
  ELSE
    InElementCheck=.FALSE.
  END IF

  ! found the particle on the current proc, set the corresponding ElemID (and save the position in reference space, if required)
  IF (InElementCheck) THEN
    PEM%Element(iPart) = ElemID
    IF(DoRefMapping) PartPosRef(1:3,iPart) = Xi
    ParticleFound = .TRUE.
    EXIT
  END IF
END DO ! iBGMElem

! particle not found
IF (.NOT.ParticleFound) THEN
  ! abort if we are required to find it on the current proc, otherwise just remove it from the list of particles
  IF(DoRelocate) CALL abort(&
__STAMP__&
  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
  RealInfo=SafetyFactor)
  PDM%ParticleInside(iPart) = .FALSE.
END IF

END SUBROUTINE SingleParticleToExactElement


SUBROUTINE SingleParticleToExactElementNoMap(iPart,doHALO,doRelocate)
!===================================================================================================================================
! this subroutine maps each particle to an element
! currently, a background mesh is used to find possible elements. if multiple elements are possible, the element with the smallest
! distance is picked as an initial guess
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,          ONLY:PartState,PEM,PDM
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:Geo
USE MOD_Particle_Utils,         ONLY:InsertionSort
USE MOD_Particle_Tracking_Vars, ONLY:Distance,ListDistance
USE MOD_Particle_Mesh_Vars,     ONLY:ElemBaryNGeo
USE MOD_Particle_MPI_Vars,      ONLY:SafetyFactor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                :: iPart
LOGICAL,INTENT(IN)                :: doHalo
LOGICAL,INTENT(IN)                :: doRelocate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iBGMElem,nBGMElems, ElemID, CellX,CellY,CellZ
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                           :: ParticleFound,InElementCheck
REAL                              :: Distance2
!===================================================================================================================================

ParticleFound = .FALSE.

! --- check if the particle is inside our bounding box (including/excluding the halo region)
IF(DoHALO)THEN
  IF ( (PartState(iPart,1).LT.GEO%xminglob).OR.(PartState(iPart,1).GT.GEO%xmaxglob).OR. &
       (PartState(iPart,2).LT.GEO%yminglob).OR.(PartState(iPart,2).GT.GEO%ymaxglob).OR. &
       (PartState(iPart,3).LT.GEO%zminglob).OR.(PartState(iPart,3).GT.GEO%zmaxglob)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
ELSE
  IF ( (PartState(iPart,1).LT.GEO%xmin).OR.(PartState(iPart,1).GT.GEO%xmax).OR. &
       (PartState(iPart,2).LT.GEO%ymin).OR.(PartState(iPart,2).GT.GEO%ymax).OR. &
       (PartState(iPart,3).LT.GEO%zmin).OR.(PartState(iPart,3).GT.GEO%zmax)) THEN
     PDM%ParticleInside(iPart) = .FALSE.
     RETURN
  END IF
END IF

! --- get background mesh cell of particle
CellX = CEILING((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
CellX = MAX(MIN(GEO%TFIBGMimax,CellX),GEO%TFIBGMimin)
CellY = CEILING((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
CellY = MAX(MIN(GEO%TFIBGMjmax,CellY),GEO%TFIBGMjmin)
CellZ = CEILING((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
CellZ = MAX(MIN(GEO%TFIBGMkmax,CellZ),GEO%TFIBGMkmin)

!--- check all cells associated with this background mesh cell
nBGMElems=GEO%TFIBGM(CellX,CellY,CellZ)%nElem

! get distance to background cell barycenter
Distance=-1.
ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GEO%TFIBGM(CellX,CellY,CellZ)%Element(iBGMElem)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF
  Distance2=(PartState(iPart,1)-ElemBaryNGeo(1,ElemID))*(PartState(iPart,1)-ElemBaryNGeo(1,ElemID)) &
           +(PartState(iPart,2)-ElemBaryNGeo(2,ElemID))*(PartState(iPart,2)-ElemBaryNGeo(2,ElemID)) &
           +(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))*(PartState(iPart,3)-ElemBaryNGeo(3,ElemID))
  IF(Distance2.GT.ElemRadius2NGeo(ElemID))THEN
    Distance(iBGMElem)=-1.
  ELSE
    Distance(iBGMElem)=Distance2
  END IF
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

! no candidate in background mesh for our particle position. Can't determine if the current position is valid because we're missing
! the BC information as well. Abort and inform the user to increase the size of the halo mesh.
!>> WARNING: If we did an invalid push, we might abort here as well.
IF(ALMOSTEQUAL(MAXVAL(Distance),-1.))THEN
  PDM%ParticleInside(iPart) = .FALSE.
  IF(DoRelocate)THEN
    IPWRITE(UNIT_StdOut,*) 'Position',PartState(iPart,1:3)
    CALL abort(&
  __STAMP__&
  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
  RealInfo=SafetyFactor)
  END IF
  RETURN
END IF

!  multiple cells associated with this background mesh cell. Sort according to distance to the cell barycenter
IF(nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

! loop through sorted list and start by closest element
DO iBGMElem=1,nBGMElems
  IF(ALMOSTEQUAL(Distance(iBGMElem),-1.))CYCLE
  ElemID=ListDistance(iBGMElem)

  ! exclude elements not within our local mesh (without the halo region)
  IF(.NOT.DoHALO)THEN
    IF(ElemID.GT.PP_nElems) CYCLE
  END IF

  CALL PartInElemCheck(PartState(iPart,1:3),iPart,ElemID,InElementCheck)

  ! no intersection found and particle is in final element
  IF(InElementCheck)THEN
    PEM%Element(iPart) = ElemID
    ParticleFound=.TRUE.
    EXIT
  END IF
END DO ! iBGMElem

! particle not found
IF (.NOT.ParticleFound) THEN
  IF(DoRelocate) CALL abort(&
  __STAMP__&
  , ' halo mesh too small. increase halo distance by increasing the safety factor. Currently Part-SafetyFactor = ',&
  RealInfo=SafetyFactor)
  PDM%ParticleInside(iPart) = .FALSE.
END IF

END SUBROUTINE SingleParticleToExactElementNoMap


SUBROUTINE PartInElemCheck(PartPos_In,PartID,ElemID,FoundInElem,IntersectPoint_Opt&
#if CODE_ANALYZE
        ,Sanity_Opt&
        ,Tol_Opt&
        ,CodeAnalyze_Opt)
#else
        )
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
! Checks if particle is in Element
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals
USE MOD_Particle_Mesh_Vars,     ONLY:ElemBaryNGeo
USE MOD_Particle_Surfaces_Vars, ONLY:SideType,SideNormVec
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartBCSideList
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Intersection,  ONLY:ComputePlanarRectIntersection
USE MOD_Particle_Intersection,  ONLY:ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection,  ONLY:ComputeBiLinearIntersection
USE MOD_Particle_Intersection,  ONLY:ComputeCurvedIntersection
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
#if CODE_ANALYZE
USE MOD_Globals,                ONLY:MyRank,UNIT_stdout
USE MOD_Mesh_Vars,              ONLY:NGeo
USE MOD_Particle_Tracking_Vars, ONLY:PartOut,MPIRankOut
USE MOD_Particle_Surfaces,      ONLY:OutputBezierControlPoints
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3d
USE MOD_Particle_Intersection,  ONLY:OutputTrajectory
#endif /*CODE_ANALYZE*/
USE MOD_Particle_Vars,          ONLY:LastPartPos
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
#if CODE_ANALYZE
REAL,INTENT(OUT),OPTIONAL                :: Tol_Opt
#endif /*CODE_ANALYZE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if CODE_ANALYZE
INTEGER                                  :: I,J,K
#endif /*CODE_ANALYZE*/
INTEGER                                  :: ilocSide,flip,SideID,BCSideID
REAL                                     :: PartTrajectory(1:3),NormVec(1:3)
REAL                                     :: lengthPartTrajectory,PartPos(1:3),LastPosTmp(1:3)
LOGICAL                                  :: isHit
REAL                                     :: alpha,eta,xi,IntersectPoint(1:3)
!===================================================================================================================================

#if CODE_ANALYZE
IF(PRESENT(tol_Opt)) tol_Opt=-1.
#endif

! virtual move to element barycenter
LastPosTmp(1:3)         = LastPartPos(PartID,1:3)
LastPartPos(PartID,1:3) = ElemBaryNGeo(1:3,ElemID)
PartPos(1:3)            = PartPos_In(1:3)

! get trajectory from element barycenter to current position
PartTrajectory          = PartPos - LastPartPos(PartID,1:3)
lengthPartTrajectory    = SQRT(PartTrajectory(1)*PartTrajectory(1) &
                             + PartTrajectory(2)*PartTrajectory(2) &
                             + PartTrajectory(3)*PartTrajectory(3))

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
  LastPartPos(PartID,1:3) = LastPosTmp(1:3)
  RETURN
END IF

! normalize the part trajectory vector
PartTrajectory=PartTrajectory/lengthPartTrajectory

! reset intersection counter and alpha
isHit = .FALSE.
alpha = -1.

! check for intersections with each side
DO ilocSide=1,6
  SideID = PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID)
  flip   = PartElemToSide(E2S_FLIP,ilocSide,ElemID)
  IF(DoRefMapping)THEN
    IF(SideID.LT.1) CYCLE
    BCSideID=SideID
    SideID=PartBCSideList(BCSideID)
    IF(SideID.LT.1) CYCLE
  END IF

  SELECT CASE(SideType(SideID))
  CASE(PLANAR_RECT)
    CALL ComputePlanarRectIntersection(ishit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta ,PartID,flip,SideID)
  CASE(PLANAR_CURVED)
    CALL ComputePlanarCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,Alpha,xi,eta,PartID,flip,SideID)
  CASE(BILINEAR,PLANAR_NONRECT)
      CALL ComputeBiLinearIntersection(isHit,PartTrajectory,lengthPartTrajectory,Alpha                &
                                                                                       ,xi            &
                                                                                       ,eta           &
                                                                                       ,PartID,SideID &
                                                                                       ,ElemCheck_Opt=.TRUE.)
  CASE(CURVED)
    CALL ComputeCurvedIntersection(isHit,PartTrajectory,lengthPartTrajectory,Alpha,xi,eta,PartID,SideID,ElemCheck_Opt=.TRUE.)
  END SELECT

#if CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(15("="))')
      WRITE(UNIT_stdout,'(A)') '     | Output after compute intersection (PartInElemCheck): '
      WRITE(UNIT_stdout,'(2(A,I0),A,L)') '     | SideType: ',SideType(SideID),' | SideID: ',SideID,'| Hit: ',isHit
      WRITE(UNIT_stdout,'(2(A,G0))')  '     | LengthPT: ',LengthPartTrajectory,' | Alpha: ',Alpha
      WRITE(UNIT_stdout,'(A,2(X,G0))') '     | Intersection xi/eta: ',xi,eta
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
  ! Dirty fix for PartInElemCheck if Lastpartpos is almost on side (tolerance issues)
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
  IF(alpha.GT.-1)THEN
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      NormVec=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    IF(flip.NE.0) NormVec=-NormVec
    IntersectPoint=LastPartPos(PartID,1:3)+alpha*PartTrajectory

#if CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,*) '     | alpha          ',alpha
      WRITE(UNIT_stdout,*) '     | Normal vector  ',NormVec
      WRITE(UNIT_stdout,*) '     | PartTrajectory ',PartTrajectory
      WRITE(UNIT_stdout,*) '     | Dotprod        ',DOT_PRODUCT(NormVec,PartTrajectory)
      WRITE(UNIT_stdout,*) '     | Point 2        ', LastPartPos(PartID,1:3)+alpha*PartTrajectory+NormVec
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
FoundInElem=.TRUE.
IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt=0.
IF(alpha.GT.-1) THEN
  FoundInElem=.FALSE.
  IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt=IntersectPoint
END IF

! reset the LastParPos to its original value
LastPartPos(PartID,1:3) = LastPosTmp(1:3)

END SUBROUTINE PartInElemCheck


SUBROUTINE ParticleInsideQuad3D(PartStateLoc,ElemID,InElementCheck,Det)
!===================================================================================================================================
! checks if particle is inside of linear element with triangulated faces
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,  ONLY : GEO
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
INTEGER                       :: ilocSide, NodeNum
LOGICAL                       :: PosCheck, NegCheck
REAL                          :: A(1:3,1:4), cross(3)
!===================================================================================================================================
InElementCheck = .TRUE.

! for all 6 sides of the element
DO iLocSide = 1,6
   !--- initialize flags for side checks
   PosCheck = .FALSE.
   NegCheck = .FALSE.

   !--- A = vector from particle to node coords
   DO NodeNum = 1,4
     A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,ElemID)) - PartStateLoc(1:3)
   END DO

   !--- compute cross product for vector 1 and 3
   cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
   cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
   cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

   !--- negative determinant of triangle 1 (points 1,3,2):
   Det(iLocSide,1) = cross(1) * A(1,2) + &
                     cross(2) * A(2,2) + &
                     cross(3) * A(3,2)
   Det(iLocSide,1) = -det(iLocSide,1)

   !--- determinant of triangle 2 (points 1,3,4):
   Det(iLocSide,2) = cross(1) * A(1,4) + &
                     cross(2) * A(2,4) + &
                     cross(3) * A(3,4)

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
   IF (GEO%ConcaveElemSide(iLocSide,ElemID)) THEN
     IF (.NOT.PosCheck) InElementCheck = .FALSE.
   ELSE
     IF (NegCheck)      InElementCheck = .FALSE.
   END IF
END DO

END SUBROUTINE ParticleInsideQuad3D

END MODULE MOD_Particle_Localization
