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

!===================================================================================================================================
!> Determines how particles interact with a given boundary condition
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Condition
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE GetBoundaryInteraction
  MODULE PROCEDURE GetBoundaryInteraction
END INTERFACE

INTERFACE GetBoundaryInteractionRef
  MODULE PROCEDURE GetBoundaryInteractionRef
END INTERFACE

INTERFACE GetBoundaryInteractionAuxBC
  MODULE PROCEDURE GetBoundaryInteractionAuxBC
END INTERFACE

INTERFACE PartSwitchElement
  MODULE PROCEDURE PartSwitchElement
END INTERFACE

PUBLIC :: GetBoundaryInteraction
PUBLIC :: GetBoundaryInteractionRef
PUBLIC :: GetBoundaryInteractionAuxBC
PUBLIC :: PartSwitchElement
!===================================================================================================================================

CONTAINS

SUBROUTINE GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,locSideID,ElemID,crossedBC&
                                  ,TriNum)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with a boundary condition
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!  PeriodicBC              = 3
!  MPINeighborhoodBC       = 6
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,                ONLY:Abort
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars, ONLY:PartBound
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
USE MOD_Mesh_Vars,              ONLY:BC,BoundaryName
USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues
USE MOD_StringTools,            ONLY:LowCase
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,flip,locSideID
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: n_loc(1:3)
CHARACTER(20)                        :: BCStringTmp
!===================================================================================================================================

IF (.NOT. ALLOCATED(PartBound%MapToPartBC)) &
  CALL abort(__STAMP__,' ERROR: PartBound not allocated!.',999,999.)

! Reset flag
crossedBC    =.FALSE.

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! If the particle is close to BC, it encounters the BC only if it leaves element/grid
  IF (.NOT.TriaTracking) THEN
    ! Calculate the side normal vector
    !if particle is close to BC, it encounters the BC only if it leaves element/grid
    IF(alpha/lengthPartTrajectory.LE.epsilontol) THEN
      SELECT CASE(SideType(SideID))
        CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
          n_loc=SideNormVec(1:3,SideID)
        CASE(BILINEAR)
          CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
        CASE(CURVED)
          CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
      END SELECT
      ! Flip side orientation if not on the master side
      IF(flip.NE.0) n_loc=-n_loc
      ! Return if particle is not leaving the element
      IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN
    END IF
  END IF

  ! Sample on outlet if requested
  !--> TODO: Move this to an AuxBC
  IF (PDM%ParticleInside(iPart)) THEN
    IF (WriteMacroSurfaceValues .AND. ErosionOutlet) THEN
      ! Match boundary name with 'outlet'
      CALL LowCase(TRIM(BoundaryName(BC(SideID))),BCStringTmp)
      IF (BCStringTmp.EQ.'outlet') THEN
        CALL SideErosion(PartTrajectory,xi,eta,iPart,SideID,flip,alpha,TriNum=TriNum)
      END IF
    END IF
  END IF

  ! Calculate balance of particle count and kinetic energy
  IF(CalcPartBalance) THEN
      nPartOut(   PartSpecies(iPart)) = nPartOut(   PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart)) = PartEkinOut(PartSpecies(iPart)) + CalcEkinPart(iPart)
  END IF

  ! Particle has left the element over an open BC
  PDM%ParticleInside(iPart) = .FALSE.
  crossedBC = .TRUE.
  alpha     = -1.

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (PDM%ParticleInside(iPart)) THEN
    ! simple reflection
    SELECT CASE(PartBound%WallModel(PartBound%MapToPartBC(BC(SideID))))
      ! perfectly reflecting, specular re-emission
      CASE('perfRef')
         CALL  PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,              &
           opt_Symmetry=.FALSE.,opt_Reflected=crossedBC,TriNum=TriNum,opt_LocalSide= locSideID,opt_ElemID=ElemID)
      ! reflection using coefficient of restitution (CoR)
      CASE('coeffRes')
         CALL  DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,              &
           opt_Symmetry=.FALSE.,opt_Reflected=crossedBC,TriNum=TriNum,opt_LocalSide= locSideID,opt_ElemID=ElemID, &
           WallCoeffModel=PartBound%WallCoeffModel(PartBound%MapToPartBC(BC(SideID))))
      CASE DEFAULT
          CALL abort(__STAMP__, ' No or invalid particle wall model given. Please update Part-Boundary[$]-WallModel.')
    END SELECT
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID &
                        ,ElemID,opt_perimoved=crossedBC,TriNum=TriNum)

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(__STAMP__,' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)',999,999.)

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(10) !PartBound%SymmetryBC
!-----------------------------------------------------------------------------------------------------------------------------------
  ! For particles, symmetry equals perfect reflection, also flip forces direction
  CALL  PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
                          opt_Symmetry=.TRUE.,opt_Reflected=crossedBC,TriNum=TriNum,opt_LocalSide= locSideID,opt_ElemID=ElemID)

!CASE(100) !PartBound%AnalyzeBC
!!-----------------------------------------------------------------------------------------------------------------------------------
!  CALL  SideAnalysis(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,ElemID &
!                    ,opt_crossed=crossedBC)

CASE DEFAULT
  CALL abort(__STAMP__,' ERROR: PartBound not associated!. (unknown case)',999,999.)
END SELECT

END SUBROUTINE GetBoundaryInteraction


SUBROUTINE GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,ElemID,crossedBC)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with a boundary condition
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!  PeriodicBC              = 3
!  MPINeighborhoodBC       = 6
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY:BC,nSides,BoundaryName
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
USE MOD_Particle_Boundary_Vars, ONLY:PartBound
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Mesh_Vars,     ONLY:PartBCSideList
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars, ONLY:SideType,SideNormVec,epsilontol
USE MOD_Particle_Tracking_Vars, ONLY:CartesianPeriodic
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,WriteMacroSurfaceValues
USE MOD_StringTools,            ONLY:LowCase
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,flip
REAL,INTENT(IN)                      :: xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
INTEGER,INTENT(INOUT)                :: ElemID
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: n_loc(1:3)
INTEGER                              :: BCSideID
CHARACTER(20)                        :: BCStringTmp
!===================================================================================================================================

IF (.NOT. ALLOCATED(PartBound%MapToPartBC)) &
  CALL abort(__STAMP__,' ERROR: PartBound not allocated!.',999,999.)

! Reset flag
crossedBC    =.FALSE.

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! If the particle is close to BC, it encounters the BC only if it leaves element/grid
  IF(alpha/lengthPartTrajectory.LE.epsilontol) THEN
    ! Calculate the side normal vector
    BCSideID=PartBCSideList(SideID)
    SELECT CASE(SideType(BCSideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        n_loc=SideNormVec(1:3,BCSideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
    END SELECT
    ! Flip side orientation if not on the master side
    IF(flip.NE.0) n_loc=-n_loc
      ! Return if particle is not leaving the element
    IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN
  END IF

  ! Sample on outlet if requested
  !--> TODO: Move this to an AuxBC
  IF (PDM%ParticleInside(iPart)) THEN
    IF (WriteMacroSurfaceValues .AND. ErosionOutlet) THEN
      ! Match boundary name with 'outlet'
      CALL LowCase(TRIM(BoundaryName(BC(SideID))),BCStringTmp)
      IF (BCStringTmp.EQ.'outlet') THEN
        CALL SideErosion(PartTrajectory,xi,eta,iPart,SideID,flip,alpha)
      END IF
    END IF
  END IF

  ! Calculate balance of particle count and kinetic energy
  IF(CalcPartBalance) THEN
      nPartOut(PartSpecies(   iPart)) = nPartOut(PartSpecies(   iPart)) + 1
      PartEkinOut(PartSpecies(iPart)) = PartEkinOut(PartSpecies(iPart)) + CalcEkinPart(iPart)
  END IF

  ! Particle has left the element over an open BC
  PDM%ParticleInside(iPart) = .FALSE.
  crossedBC = .TRUE.
  alpha     = -1.

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  BCSideID=PartBCSideList(SideID)
  ! simple reflection
  IF (PDM%ParticleInside(iPart)) THEN
    SELECT CASE(PartBound%WallModel(PartBound%MapToPartBC(BC(SideID))))
      ! perfectly reflecting, specular re-emission
      CASE('perfRef')
        CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,           &
                               BCSideID=BCSideID,opt_Symmetry=.FALSE.,opt_Reflected=crossedBC)
      CASE('coeffRes')
        CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,           &
                                BCSideID=BCSideID,opt_Symmetry=.FALSE.,opt_Reflected=crossedBC,              &
                                WallCoeffModel=PartBound%WallCoeffModel(PartBound%MapToPartBC(BC(SideID))))
      CASE DEFAULT
          CALL abort(__STAMP__, ' No or invalid particle wall model given. Please update Part-Boundary[$]-WallModel.')
      END SELECT
  ELSE
    ! not inside any-more, removed in last step
    crossedBC=.TRUE.
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! sanity check
  IF (CartesianPeriodic) CALL abort(__STAMP__,' No periodic BCs for CartesianPeriodic!')

  ! move particle periodic distance
  BCSideID=PartBCSideList(SideID)
  CALL PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,ElemID,BCSideID=BCSideID,       &
                  opt_perimoved=crossedBC)

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(__STAMP__,' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)',999,999.)

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(10) !PartBound%SymmetryBC
!-----------------------------------------------------------------------------------------------------------------------------------
  BCSideID=PartBCSideList(SideID)
  CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,             &
                         BCSideID=BCSideID,opt_Symmetry=.TRUE.,opt_reflected=crossedBC)

CASE DEFAULT
  CALL abort(__STAMP__,' ERROR: PartBound not associated!. BC(SideID)',BC(SideID),REAL(SideID/nSides))
END SELECT

END SUBROUTINE GetBoundaryInteractionRef


SUBROUTINE GetBoundaryInteractionAuxBC(PartTrajectory,lengthPartTrajectory,alpha,iPart,AuxBCIdx,crossedBC)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with an auxBC
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies
USE MOD_Particle_Boundary_Vars, ONLY:PartAuxBC
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,AuxBCIdx
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Reset flag
crossedBC    =.FALSE.

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartAuxBC%TargetBoundCond(AuxBCIdx))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartAuxBC%OpenBC
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Calculate balance of particle count and kinetic energy
  IF(CalcPartBalance) THEN
      nPartOut(PartSpecies(   iPart)) = nPartOut(PartSpecies(   iPart)) + 1
      PartEkinOut(PartSpecies(iPart)) = PartEkinOut(PartSpecies(iPart)) + CalcEkinPart(iPart)
  END IF

  ! Particle has left the element over an open BC
  PDM%ParticleInside(iPart) = .FALSE.
  crossedBC = .TRUE.
  alpha     = -1.

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartAuxBC%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (PDM%ParticleInside(iPart)) THEN
    ! simple reflection, auxBC can not have diffuse reflection
    ! perfectly reflecting, specular re-emission
      CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,flip=-1, &
             opt_Symmetry=.FALSE.,opt_Reflected=crossedBC,AuxBCIdx=AuxBCIdx)
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
  CALL abort(__STAMP__,' ERROR: AuxBC bound not associated!. (unknown case)',999,999.)
END SELECT

END SUBROUTINE GetBoundaryInteractionAuxBC


SUBROUTINE PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,flip,BCSideID, &
  opt_Symmetry,opt_Reflected,TriNum,AuxBCIdx,opt_LocalSide,opt_ElemID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,Species,PartSpecies,PartReflCount
USE MOD_Particle_Surfaces_Vars, ONLY:SideNormVec,SideType
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_Particle_Vars,          ONLY:Pt_temp,PDM
USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues
USE MOD_ErosionPoints,          ONLY:RecordErosionPoint
USE MOD_ErosionPoints_Vars,     ONLY:EP_inUse
USE MOD_Particle_Erosion_Vars,  ONLY:PartTrackReflection
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID, flip
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
INTEGER,INTENT(IN),OPTIONAL       :: opt_LocalSide
INTEGER,INTENT(IN),OPTIONAL       :: opt_ElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3),n_loc(1:3),WallVelo(3),intersec(3),r_vec(3),axis(3),cos2inv
REAL                              :: Xitild,EtaTild
INTEGER                           :: p,q,SurfSideID,locBCID
LOGICAL                           :: Symmetry,IsAuxBC
REAL                              :: PartFaceAngle,PartFaceAngle_old
REAL                              :: PartTrajectory_old(3)
REAL                              :: v_magnitude
!===================================================================================================================================

! Check if reflected on AuxBC
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF

! Reflected on AuxBC
IF (IsAuxBC) THEN
  ! Get type of AuxBC and the corresponding normal vector
  SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
    CASE ('plane')
      n_loc    = AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec

    CASE ('cylinder')
      intersec = LastPartPos(1:3,PartID) + alpha*PartTrajectory
      r_vec    = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
      axis     = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
      n_loc    = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis) ) )
      IF (.NOT.AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc

    CASE ('cone')
      intersec = LastPartPos(1:3,PartID) + alpha*PartTrajectory
      r_vec    = AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
      axis     = AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
      cos2inv  = 1./COS(AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle)**2
      n_loc    = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis)*cos2inv ) )
      IF (.NOT.AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc

    CASE ('parabol')
      intersec = LastPartPos(1:3,PartID) + alpha*PartTrajectory
      r_vec    = AuxBC_parabol(AuxBCMap(AuxBCIdx))%r_vec
      axis     = AuxBC_parabol(AuxBCMap(AuxBCIdx))%axis
      n_loc    = UNITVECTOR( intersec - ( r_vec + axis*(DOT_PRODUCT(intersec-r_vec,axis)+0.5*AuxBC_parabol(AuxBCMap(AuxBCIdx))%zfac) ) )
      IF (.NOT.AuxBC_parabol(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc

    CASE DEFAULT
      CALL abort(__STAMP__,'AuxBC does not exist')
  END SELECT

  ! Check if particle is entering the domain through AuxBC
  IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected = .FALSE.
    CALL abort(__STAMP__,'Error in PerfectReflection: Particle coming from outside!')
  ! Particle is leaving the domain through AuxBC
  ELSE IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected = .TRUE.
  ELSE
    CALL abort(__STAMP__,'Error in PerfectReflection: n_vec is perpendicular to PartTrajectory for AuxBC',AuxBCIdx)
  END IF
  WallVelo=PartAuxBC%WallVelo(1:3,AuxBCIdx)

! Normal BC
ELSE
  ! Get wall velo and BCID
  WallVelo  = PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(SideID)))
  locBCID   = PartBound%MapToPartBC(BC(SideID))

  ! Instead of checking DoRefMapping, we pass BCSideID if true
  IF(PRESENT(BCSideID)) THEN
    ! Get side normal vector
    SELECT CASE(SideType(BCSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,BCSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
    END SELECT
  ! not reference mapping
  ELSE
    ! Get side normal vector
    IF (TriaTracking) THEN
      CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,ElemID_opt=opt_ElemID,LocSideID_opt=opt_LocalSide)
    ELSE
      SELECT CASE(SideType(SideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        n_loc=SideNormVec(1:3,SideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
      END SELECT
      ! Flip side orientation if not on the master side
      IF(flip.NE.0) n_loc=-n_loc
    END IF
  END IF

  IF(PRESENT(opt_Symmetry)) THEN
    Symmetry = opt_Symmetry
  ELSE
    Symmetry = .FALSE.
  END IF

  ! Check if particle is encountering the boundary
  IF(DOT_PRODUCT(PartTrajectory,n_loc).LE.0.) THEN
    ! Particle does not encounter the boundary, return
    IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
    RETURN
  ELSE
    IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
  END IF

END IF !IsAuxBC

! Make sure we have to old velocity safe, then update particle velocity
v_old                = PartState(4:6,PartID)
PartState(4:6,PartID)= PartState(4:6,PartID)-2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo

! Ugly hack to catch limited machine accuracy resulting in case DOT_PRODUCT greater than 1
IF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),1.0) .AND. &
  DOT_PRODUCT(PartTrajectory,n_loc) > 1.0) THEN
  PartFaceAngle=ABS(0.5*PI - ACOS(1.))
ELSEIF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),-1.0) .AND. &
  DOT_PRODUCT(PartTrajectory,n_loc) < -1.0) THEN
  PartFaceAngle=ABS(0.5*PI - ACOS(-1.))
ELSE
  PartFaceAngle=ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
ENDIF
! End ugly hack

! Sample macrovalues on boundary
IF ((.NOT.IsAuxBC) .AND. WriteMacroSurfaceValues) THEN
  ! Find correct boundary on SurfMesh
  SurfSideID=SurfMesh%SideIDToSurfID(SideID)
  ! compute p and q for supersampling
  ! correction of xi and eta, can only be applied if xi & eta are not used later!
  IF (TriaTracking) THEN
    p=1 ; q=1
  ELSE
    Xitild =MIN(MAX(-1.,xi ),0.99)
    Etatild=MIN(MAX(-1.,eta),0.99)
    p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
    q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
  END IF

!  ! Ugly hack to catch limited machine accuracy resulting in case DOT_PRODUCT greater than 1
!  IF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),1.0) .AND.              &
!    DOT_PRODUCT(PartTrajectory,n_loc) > 1.0) THEN
!    PartFaceAngle = ABS(0.5*PI - ACOS(1.))
!  ELSEIF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),-1.0) .AND.         &
!    DOT_PRODUCT(PartTrajectory,n_loc) < -1.0) THEN
!    PartFaceAngle = ABS(0.5*PI - ACOS(-1.))
!  ELSE
!    PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
!  ENDIF
!  ! End ugly hack

  ! Record particle impact
  CALL RecordParticleBoundarySampling(PartID,SurfSideID,p,q,v_old,PartFaceAngle)
END IF !.NOT.IsAuxBC

! Make sure we have the old values safe
IF (EP_inUse) THEN
  PartTrajectory_old = PartTrajectory
  PartFaceAngle_old  = PartFaceAngle
END IF

! set particle position on face
!--> first move the particle to the boundary
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
!--> flip trajectory and move the remainder of the particle push
PartTrajectory(1:3)     = PartTrajectory(1:3)-2.*DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc
PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*(lengthPartTrajectory - alpha)

! compute moved particle || rest of movement
PartTrajectory          = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory    = SQRT(PartTrajectory(1)*PartTrajectory(1)          &
                          +    PartTrajectory(2)*PartTrajectory(2)          &
                          +    PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory          = PartTrajectory/lengthPartTrajectory

! Recording of individual particle impacts
IF (EP_inUse) THEN
!  ! Ugly hack to catch limited machine accuracy resulting in case DOT_PRODUCT greater than 1
!  IF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),1.0) .AND.            &
!    DOT_PRODUCT(PartTrajectory,n_loc) > 1.0) THEN
!    PartFaceAngle = ABS(0.5*PI - ACOS(1.))
!  ELSEIF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),-1.0) .AND.       &
!    DOT_PRODUCT(PartTrajectory,n_loc) < -1.0) THEN
!    PartFaceAngle = ABS(0.5*PI - ACOS(-1.))
!  ELSE
!    PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
!  ENDIF
!  ! End ugly hack

  CALL RecordErosionPoint(BCSideID        = BC(SideID),                   &
                          PartID          = PartID,                       &
                          PartFaceAngle   = PartFaceAngle,                &
                          v_old           = v_old,                        &
                          PartFaceAngle_old =PartFaceAngle_old,           &
                          PartReflCount   = PartReflCount(PartID),        &
                          alpha           = alpha)
END IF

! Increase reflection counter
IF (PartTrackReflection) PartReflCount(PartID) = PartReflCount(PartID) + 1

! Correction for Runge-Kutta. Reflect or delete force history / acceleration
IF (.NOT.ALMOSTZERO(DOT_PRODUCT(WallVelo,WallVelo))) THEN
  ! Reconstruction in timedisc during push
  PDM%IsNewPart(PartID) = .TRUE.
ELSE
  Pt_temp(1:3,PartID) = Pt_temp(1:3,PartID)-2.*DOT_PRODUCT(Pt_temp(1:3,PartID),n_loc)*n_loc
  ! Reflect also force history for symmetry
  IF (Symmetry) THEN
    Pt_temp(4:6,PartID) = Pt_temp(4:6,PartID)-2.*DOT_PRODUCT(Pt_temp(4:6,PartID),n_loc)*n_loc
  ! Produces best result compared to analytical solution in place capacitor ...
  ELSE
    Pt_temp(4:6,PartID) = 0.
  END IF
END IF

! Remove sliding low velocity particles
IF (LowVeloRemove) THEN
    v_magnitude   = SQRT(DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID)))

    IF ((Species(PartSpecies(PartID))%LowVeloThreshold.NE.0).AND.(v_magnitude.LT.Species(PartSpecies(PartID))%LowVeloThreshold))THEN
          Species(PartSpecies(PartID))%LowVeloCounter = Species(PartSpecies(PartID))%LowVeloCounter + 1
          PDM%ParticleInside(PartID) = .FALSE.
          IPWRITE(UNIT_stdOut,*) ' Low velocity particle removed after impact. Velocity after reflection:', v_magnitude
    END IF
END IF

END SUBROUTINE PerfectReflection


SUBROUTINE DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,flip,BCSideID, &
  opt_Symmetry,opt_Reflected,TriNum,AuxBCIdx,opt_LocalSide,opt_ElemID,WallCoeffModel)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
! only implemented for DoRefMapping tracking
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_ErosionPoints,          ONLY:RecordErosionPoint
USE MOD_ErosionPoints_Vars,     ONLY:EP_inUse
USE MOD_Mesh_Vars,              ONLY:NGeo,BC
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Erosion_Vars,  ONLY:PartTrackReflection
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars, ONLY:SideNormVec,SideType,BezierControlPoints3D
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,Species,PartSpecies,PartReflCount
USE MOD_Particle_Vars,          ONLY:Pt_temp,PDM
USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues
USE MOD_TimeDisc_Vars,          ONLY:TimeDiscType
! Bons particle rebound model
USE MOD_Particle_Interpolation_Vars,ONLY:FieldAtParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID, flip
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
INTEGER,INTENT(IN),OPTIONAL       :: opt_LocalSide
INTEGER,INTENT(IN),OPTIONAL       :: opt_ElemID
CHARACTER(LEN=255),INTENT(IN)     :: WallCoeffModel
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3),n_loc(1:3),WallVelo(3),intersec(3),r_vec(3),axis(3),cos2inv
REAL                              :: tang1(1:3),tang2(1:3)
REAL                              :: Xitild,EtaTild
INTEGER                           :: p,q,SurfSideID,locBCID
LOGICAL                           :: Symmetry, IsAuxBC
REAL                              :: PartFaceAngle,PartFaceAngleDeg,PartFaceAngle_old
REAL                              :: v_magnitude,v_norm(3),v_tang(3)
REAL                              :: PartTrajectory_old(3)
REAL                              :: PartTrajectoryTang(3),PartTrajectoryNorm(3)
REAL                              :: PartTrajectory_rescale
REAL                              :: eps_n, eps_t, lengthPartTrajectory_old
! Bons particle rebound model
REAL                              :: E_eff
REAL                              :: Vol,r,w,w_crit,sigma_y
!===================================================================================================================================

! Check if reflected on AuxBC
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF

! Reflected on AuxBC
IF (IsAuxBC) THEN
  ! Get type of AuxBC and the corresponding normal vector
  SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
    CASE ('plane')
      n_loc    = AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec

    CASE ('cylinder')
      intersec = LastPartPos(1:3,PartID) + alpha*PartTrajectory
      r_vec    = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
      axis     = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
      n_loc    = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis) ) )
      IF (.NOT.AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc

    CASE ('cone')
      intersec = LastPartPos(1:3,PartID) + alpha*PartTrajectory
      r_vec    = AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
      axis     = AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
      cos2inv  = 1./COS(AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle)**2
      n_loc    = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis)*cos2inv ) )
      IF (.NOT.AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc

    CASE ('parabol')
      intersec = LastPartPos(1:3,PartID) + alpha*PartTrajectory
      r_vec    = AuxBC_parabol(AuxBCMap(AuxBCIdx))%r_vec
      axis     = AuxBC_parabol(AuxBCMap(AuxBCIdx))%axis
      n_loc    = UNITVECTOR( intersec - ( r_vec + axis*(DOT_PRODUCT(intersec-r_vec,axis)+0.5*AuxBC_parabol(AuxBCMap(AuxBCIdx))%zfac) ) )
      IF (.NOT.AuxBC_parabol(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc

    CASE DEFAULT
      CALL abort(__STAMP__,'AuxBC does not exist')
  END SELECT

  ! Check if particle is entering the domain through AuxBC
  IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
    CALL abort(__STAMP__,'Error in DiffuseReflection: Particle coming from outside!')
  ! Particle is leaving the domain through AuxBC
  ELSE IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected = .TRUE.
  ELSE
    CALL abort(__STAMP__,'Error in DiffuseReflection: n_vec is perpendicular to PartTrajectory for AuxBC',AuxBCIdx)
  END IF
  WallVelo=PartAuxBC%WallVelo(1:3,AuxBCIdx)
  IF (n_loc(3).NE.0.) THEN
    tang1(1) = 1.0
    tang1(2) = 1.0
    tang1(3) = -(n_loc(1)+n_loc(2))/n_loc(3)
  ELSE
    IF (n_loc(2).NE.0.) THEN
      tang1(1) = 1.0
      tang1(3) = 1.0
      tang1(2) = -(n_loc(1)+n_loc(3))/n_loc(2)
    ELSE
      IF (n_loc(1).NE.0.) THEN
        tang1(2) = 1.0
        tang1(3) = 1.0
        tang1(1) = -(n_loc(2)+n_loc(3))/n_loc(1)
      ELSE
        CALL abort(&
__STAMP__&
,'Error in DiffuseReflection, n_vec is zero for AuxBC',AuxBCIdx)
      END IF
    END IF
  END IF
  tang1=UNITVECTOR(tang1)
  tang2=CROSSNORM(n_loc,tang1)

! Normal BC
ELSE
  ! Get wall velo and BCID
  WallVelo  = PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(SideID)))
  locBCID   = PartBound%MapToPartBC(BC(SideID))

  ! Instead of checking DoRefMapping, we pass BCSideID if true
  IF(PRESENT(BCSideID))THEN
    ! Get side normal vector
    SELECT CASE(SideType(BCSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,BCSideID)
      tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID))
      tang2=CROSSNORM(n_loc,tang1)
      !tang2=BezierControlPoints3D(:,0,NGeo,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,BCSideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
    END SELECT
  ! not reference mapping
  ELSE
    ! Get side normal vector
    IF (TriaTracking) THEN
      CALL CalcNormAndTangTriangle(nVec=n_loc,tang1=tang1,tang2=tang2, &
          TriNum=TriNum,ElemID_opt=opt_ElemID,LocSideID_opt=opt_LocalSide)
    ELSE
      SELECT CASE(SideType(SideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        n_loc=SideNormVec(1:3,SideID)
        tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,SideID)-BezierControlPoints3D(:,0,0,SideID))
        tang2=CROSSNORM(n_loc,tang1)
        !tang2=BezierControlPoints3D(:,0,NGeo,SideID)-BezierControlPoints3D(:,0,0,SideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,SideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,SideID)
      END SELECT
      ! Flip side orientation if not on the master side
      IF(flip.NE.0) n_loc=-n_loc
    END IF
  END IF

  IF(PRESENT(opt_Symmetry)) THEN
    Symmetry = opt_Symmetry
  ELSE
    Symmetry = .FALSE.
  END IF

  ! Check if particle is encountering the boundary
  IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
    RETURN
  ELSE
    IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
  END IF
END IF !IsAuxBC

! Make sure we have the old velocity safe
v_old   = PartState(4:6,PartID)

! Respect coefficient of restitution
v_norm  = DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc
v_tang  = PartState(4:6,PartID) - v_norm

! Ugly hack to catch limited machine accuracy resulting in case DOT_PRODUCT greater than 1
IF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),1.0) .AND.                    &
  DOT_PRODUCT(PartTrajectory,n_loc) > 1.0) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(1.))
ELSEIF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),-1.0) .AND.               &
  DOT_PRODUCT(PartTrajectory,n_loc) < -1.0) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(-1.))
ELSE
  PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
ENDIF
! End ugly hack

SELECT CASE(WallCoeffModel)
  ! Tabaoff, W.; Wakeman, T.: Basic Erosion Investigation in Small Turbomachinery. / Cincinnati Univ. OH, 1981
  ! >> Ignored variance (random fluctuations in reflection values) and ONLY took mean values!
  CASE('Tabakoff1981')
    ! Transfer from radians to degree
    PartFaceAngleDeg = PartFaceAngle * 180/PI

    ! Compute cubic polynomials for coefficient of restitution
    eps_n = 1.    - 0.0211 * PartFaceAngleDeg  + 0.000228 * PartFaceAngleDeg**2. - 0.000000876 * PartFaceAngleDeg**3.
    eps_t = 0.953                              - 0.00446  * PartFaceAngleDeg**2. + 0.000000648 * PartFaceAngleDeg**3.

  ! Bons, J., Prenter, R., Whitaker, S.: A Simple Physics-Based Model for Particle Rebound and Deposition in Turbomachinery
  ! / J. Turbomach 139(8), 2017
  CASE('Bons2017')
    ! Assume spherical particles for now
    Vol     = Species(PartSpecies(PartID))%MassIC/Species(PartSpecies(PartID))%DensityIC
    r       = (3.*Vol/4./PI)**(1./3.)

    ! Calculate deformation of cylindrical model particle
    w       = 4./3. * r * SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) * (Species(PartSpecies(PartID))%DensityIC &
                                                                       /Species(PartSpecies(PartID))%YoungIC   )**0.5;

    ! Find composite elastic modulus
    E_eff   = ((1. - Species(PartSpecies(PartID))%PoissonIC**2.)/Species(PartSpecies(PartID))%YoungIC +        &
               (1. - PartBound%Poisson(locBCID)            **2.)/PartBound%Young(locBCID)               )**(-1.)

    ! Find critical deformation
    sigma_y = Species(PartSpecies(PartID))%YieldCoeff*SQRT(DOT_PRODUCT(v_old(1:3),v_old(1:3)))
    w_crit  = sigma_y * 4./3. * r / E_eff

    ! Normal coefficient of restitution
    IF (w .LT. w_crit) THEN
      eps_n = 1./SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) *(sigma_y**2. / (Species(PartSpecies(PartID))%DensityIC &
                                                                          *  Species(PartSpecies(PartID))%YoungIC   ))**0.5
    ELSE
      eps_n = 1.;
    END IF

    !> Do not account for adhesion for now <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Tangential coefficient of restitution

    !> Assume change in density from last particle position to wall position to be negligible
    ! Original relation by Barker, B., Casaday, B., Shankara, P., Ameri, A., and Bons, J. P., 2013.
    !> Cosine term added by Bons, J., Prenter, R., Whitaker, S., 2017.
    eps_t   = 1. - FieldAtParticle(PartID,1) / SQRT(DOT_PRODUCT(v_tang(1:3),v_tang(1:3)))  * (eps_n * &
                                               SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3)))) * (1/(eps_n)+1)*COS(PartFaceAngle)**2.

  ! Fong, W.; Amili, O.; Coletti, F.: Velocity and spatial distribution of intertial particles in a turbulent channel flow
  ! / J. FluidMech 872, 2019
  CASE('Fong2019')
    ! Reuse YieldCoeff to modify the normal velocity, keep tangential velocity
    eps_t   = 1.
    eps_n   = PartBound%CoR(locBCID)

  CASE DEFAULT
      CALL abort(__STAMP__, ' No particle wall coefficients given. This should not happen.')

END SELECT

! Sample macrovalues on boundary
IF ((.NOT.IsAuxBC) .AND. WriteMacroSurfaceValues) THEN
  ! Find correct boundary on SurfMesh
  SurfSideID=SurfMesh%SideIDToSurfID(SideID)
  ! compute p and q for supersampling
  ! correction of xi and eta, can only be applied if xi & eta are not used later!
  IF (TriaTracking) THEN
    p=1 ; q=1
  ELSE
    Xitild =MIN(MAX(-1.,xi ),0.99)
    Etatild=MIN(MAX(-1.,eta),0.99)
    p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
    q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
  END IF

  ! Ugly hack to catch limited machine accuracy resulting in case DOT_PRODUCT greater than 1
  IF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),1.0) .AND.                          &
      DOT_PRODUCT(PartTrajectory,n_loc)>1.0) THEN
      PartFaceAngle=ABS(0.5*PI - ACOS(1.))
  ELSEIF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),-1.0) .AND.                     &
      DOT_PRODUCT(PartTrajectory,n_loc)<-1.0) THEN
      PartFaceAngle=ABS(0.5*PI - ACOS(-1.))
  ELSE
  PartFaceAngle=ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
  ENDIF
  ! End ugly hack

  ! Record particle impact
  CALL RecordParticleBoundarySampling(PartID,SurfSideID,p,q,v_old,PartFaceAngle)
END IF !.NOT.IsAuxBC

! Make sure we have the old values safe
IF (EP_inUse) THEN
    PartTrajectory_old = PartTrajectory
    PartFaceAngle_old  = PartFaceAngle
END IF

! set particle position on face
!--> first move the particle to the boundary
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
!--> flip trajectory and move the remainder of the particle push
!===================================================================================================================================
! > modified with coefficients of restitution
! > WARNING: this only works with LSERK. Check and abort if not fulfilled
IF (.NOT.(TimeDiscType.EQ.'LSERKW2').AND..NOT.(TimeDiscType.EQ.'LSWERKW3'))               &
    CALL ABORT(__STAMP__,                                                                 &
    'Time discretization '//TRIM(TimeDiscType)//' is incompatible with current implementation of coefficients of restitution.')

PartTrajectoryTang(1:3) = eps_t*(PartTrajectory(1:3) - DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc)
PartTrajectoryNorm(1:3) = eps_n*(DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc)
PartTrajectory(1:3)     = PartTrajectoryTang(1:3) - PartTrajectoryNorm(1:3)

! Save the old lengthPartTrajectory and rescale to new. We can't compute a complete new lengthPartTrajectory as we do not have the
! current LSERK time step here
lengthPartTrajectory_old= lengthPartTrajectory
PartTrajectory_rescale  = SQRT(PartTrajectory(1)*PartTrajectory(1)                        &
                          +    PartTrajectory(2)*PartTrajectory(2)                        &
                          +    PartTrajectory(3)*PartTrajectory(3) )
lengthPartTrajectory    = lengthPartTrajectory  *PartTrajectory_rescale

! Rescale the PartTrajectory to uniform length, otherwise we apply the coefficients twice
PartTrajectory          = PartTrajectory        /PartTrajectory_rescale

! lengthPartTrajectory is coupled to alpha, so use the old value to non-dimensionalize it
PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1-alpha/lengthPartTrajectory_old) * PartTrajectory(1:3) * lengthPartTrajectory

! compute moved particle || rest of movement
PartTrajectory          = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory    = SQRT(PartTrajectory(1)*PartTrajectory(1)                        &
                          +    PartTrajectory(2)*PartTrajectory(2)                        &
                          +    PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory          = PartTrajectory/lengthPartTrajectory


! compute new particle velocity, modified with coefficents of restitution
PartState(4:6,PartID)= eps_t * v_tang - eps_n * v_norm + WallVelo

IF (EP_inUse) THEN
  ! Ugly hack to catch limited machine accuracy resulting in case DOT_PRODUCT greater than 1
  IF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),1.0) .AND.                            &
    DOT_PRODUCT(PartTrajectory,n_loc) > 1.0) THEN
    PartFaceAngle = ABS(0.5*PI - ACOS(1.))
  ELSEIF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),-1.0) .AND.                       &
    DOT_PRODUCT(PartTrajectory,n_loc) < -1.0) THEN
    PartFaceAngle = ABS(0.5*PI - ACOS(-1.))
  ELSE
    PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
  ENDIF
  ! End ugly hack

  CALL RecordErosionPoint(BCSideID        = BC(SideID),                                   &
                          PartID          = PartID,                                       &
                          PartFaceAngle   = PartFaceAngle,                                &
                          v_old           = v_old,                                        &
                          PartFaceAngle_old =PartFaceAngle_old,                           &
                          PartReflCount   = PartReflCount(PartID),                        &
                          alpha           = alpha)
END IF

! increase reflection counter
IF (PartTrackReflection) PartReflCount(PartID) = PartReflCount(PartID) + 1

! Correction for Runge-Kutta. Reflect or delete force history / acceleration
IF (.NOT.ALMOSTZERO(DOT_PRODUCT(WallVelo,WallVelo))) THEN
  ! Reconstruction in timedisc during push
  PDM%IsNewPart(PartID)=.TRUE.
ELSE
  Pt_temp(1:3,PartID)=Pt_temp(1:3,PartID)-2.*DOT_PRODUCT(Pt_temp(1:3,PartID),n_loc)*n_loc
  ! Reflect also force history for symmetry
  IF (Symmetry) THEN
    Pt_temp(4:6,PartID)=Pt_temp(4:6,PartID)-2.*DOT_PRODUCT(Pt_temp(4:6,PartID),n_loc)*n_loc
  ! Produces best result compared to analytical solution in place capacitor ...
  ELSE
    Pt_temp(4:6,PartID)=0.
  END IF
END IF

! Remove sliding low velocity particles
IF (LowVeloRemove) THEN
    v_magnitude   = SQRT(DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID)))

    IF ((Species(PartSpecies(PartID))%LowVeloThreshold.NE.0).AND.(v_magnitude.LT.Species(PartSpecies(PartID))%LowVeloThreshold))THEN
          Species(PartSpecies(PartID))%LowVeloCounter = Species(PartSpecies(PartID))%LowVeloCounter + 1
          PDM%ParticleInside(PartID) = .FALSE.
          IPWRITE(UNIT_stdOut,*) ' Low velocity particle removed after impact. Velocity after reflection:', v_magnitude
    END IF
END IF

END SUBROUTINE DiffuseReflection


SUBROUTINE PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,ElemID,BCSideID,opt_perimoved,TriNum)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the periodic shift in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking,DoRefMapping
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,SidePeriodicType
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,PEM
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType
USE MOD_Particle_Mesh_Vars,     ONLY:PartSideToElem
#if CODE_ANALYZE
USE MOD_Particle_Tracking_Vars,  ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_perimoved
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: n_loc(1:3)
INTEGER                           :: PVID,moved(2),locSideID
!===================================================================================================================================

! Instead of checking DoRefMapping, we pass opt_BCSideID if true
IF(PRESENT(BCSideID))THEN
  ! Get side normal vector
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  END SELECT
! not reference mapping
ELSE
  ! Get side normal vector
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
  ELSE
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT
  END IF
END IF

! Check if particle is encountering the boundary
IF(DOT_PRODUCT(PartTrajectory,n_loc).LE.0.) THEN
  ! Particle does not encounter the boundary, return
  IF(PRESENT(opt_perimoved)) opt_perimoved=.FALSE.
  RETURN
ELSE
  IF(PRESENT(opt_perimoved)) opt_perimoved=.TRUE.
END IF

! Get periodic vector ID
PVID = SidePeriodicType(SideID)

#if CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' ParticlePosition: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' LastPartPos:      ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/

! set last particle position on face
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
! perform the periodic movement
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
! update particle positon after periodic BC
PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (lengthPartTrajectory-alpha)*PartTrajectory
lengthPartTrajectory    = lengthPartTrajectory - alpha

#if CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' ParticlePosition-pp: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' LastPartPo-pp:       ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/

! refmapping and tracing
! move particle from old element to new element
locSideID = PartSideToElem(S2E_LOC_SIDE_ID,SideID)
Moved     = PARTSWITCHELEMENT(xi,eta,locSideID,SideID,ElemID)
ElemID    = Moved(1)

#if USE_MPI
! New element is neither on proc nor in halo region, abort
IF(ElemID.EQ.-1) &
  CALL abort(__STAMP__,' Halo region to small. Neighbor element is missing!')
#endif /*MPI*/

! Remove history for last element
IF (DoRefMapping) PEM%LastElement(PartID) = 0

END SUBROUTINE PeriodicBC


SUBROUTINE SideErosion(PartTrajectory,xi,eta,PartID,SideID,flip,alpha,BCSideID,TriNum)
!----------------------------------------------------------------------------------------------------------------------------------!
! Tracks erosion on designated sides other than reflective wall
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,PartReflCount
USE MOD_Particle_Surfaces_Vars, ONLY:SideNormVec,SideType
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_ErosionPoints,          ONLY:RecordErosionPoint
USE MOD_ErosionPoints_Vars,     ONLY:EP_inUse
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3)
REAL,INTENT(IN)                   :: xi, eta, alpha
INTEGER,INTENT(IN)                :: PartID, SideID, flip
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3),n_loc(1:3)
REAL                              :: Xitild,EtaTild
INTEGER                           :: p,q,SurfSideID,locBCID
REAL                              :: PartFaceAngle
!===================================================================================================================================
  ! Get BCID
locBCID   = PartBound%MapToPartBC(BC(SideID))

! Instead of checking DoRefMapping, we pass opt_BCSideID if true
IF(PRESENT(BCSideID))THEN
  ! Get side normal vector
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  END SELECT
! not reference mapping
ELSE
  ! Get side normal vector
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
  ELSE
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    IF(flip.NE.0) n_loc=-n_loc
  END IF
END IF

! Wall sampling Macrovalues
! Find correct boundary on SurfMesh
SurfSideID=SurfMesh%SideIDToSurfID(SideID)

! Make sure we have to old velocity safe
v_old = PartState(4:6,PartID)

! compute p and q for supersampling
! correction of xi and eta, can only be applied if xi & eta are not used later!
IF (TriaTracking) THEN
  p=1 ; q=1
ELSE
  Xitild =MIN(MAX(-1.,xi ),0.99)
  Etatild=MIN(MAX(-1.,eta),0.99)
  p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
  q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
END IF

! Ugly hack to catch limited machine accuracy resulting in case DOT_PRODUCT greater than 1
IF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),1.0) .AND.                              &
    DOT_PRODUCT(PartTrajectory,n_loc)>1.0) THEN
    PartFaceAngle=ABS(0.5*PI - ACOS(1.))
ELSEIF (ALMOSTEQUAL(DOT_PRODUCT(PartTrajectory,n_loc),-1.0) .AND.                         &
    DOT_PRODUCT(PartTrajectory,n_loc)<-1.0) THEN
    PartFaceAngle=ABS(0.5*PI - ACOS(-1.))
ELSE
PartFaceAngle=ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
ENDIF
! End ugly hack

CALL RecordParticleBoundarySampling(PartID,SurfSideID,p,q,v_old,PartFaceAngle)

IF (EP_inUse) CALL RecordErosionPoint(BCSideID        = BC(SideID),                       &
                                      PartID          = PartID,                           &
                                      PartFaceAngle   = PartFaceAngle,                    &
                                      v_old           = v_old,                            &
                                      PartFaceAngle_old =PartFaceAngle,                   &
                                      PartReflCount   = PartReflCount(PartID),            &
                                      alpha           = alpha)

END SUBROUTINE SideErosion


SUBROUTINE RecordParticleBoundarySampling(PartID,SurfSideID,p,q,v_old,PartFaceAngle)
!----------------------------------------------------------------------------------------------------------------------------------!
! Combined routine to add calculated erosion variables to tracking array
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Vars,          ONLY:Species,PartSpecies,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: PartFaceAngle, v_old(1:3)
INTEGER,INTENT(IN)                :: PartID, SurfSideID, p, q
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: delta                          ! Reusable variable for variance calculation
REAL                              :: delta2                         ! Reusable variable for variance calculation
INTEGER                           :: nShift                         ! Shift amount for species tracking
REAL                              :: v_magnitude
REAL                              :: e_kin
!===================================================================================================================================
nShift        = PartSpecies(PartID) * nErosionVars

!===================================================================================================================================
! SAMP WALL
!===================================================================================================================================

!----  Sampling kinetic energy at walls
v_magnitude   = SQRT(DOT_PRODUCT(v_old(1:3),v_old(1:3)))
e_kin         = .5*Species(PartSpecies(PartID))%MassIC*v_magnitude**2.

!---- All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
!---- 1. - .. / Impact Counter
    SampWall(SurfSideID)%State(1,p,q)                       = SampWall(SurfSideID)%State(1,p,q) + 1
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN
    SampWall(SurfSideID)%State(1+nShift,p,q)                = SampWall(SurfSideID)%State(1+nShift,p,q) + 1
END IF

!===================================================================================================================================
!---- 2. - 6. / Kinetic energy on impact (mean, min, max, M2, variance)
!<<< Welford's algorithm for variance
!    (count, mean, M2) = existingAggregate
!    count = count + 1
!    delta = newValue - mean
!    mean = mean + delta / count
!    delta2 = newValue - mean
!    M2 = M2 + delta * delta2

!   Record first impact, otherwise min will be frozen at zero
    IF (SampWall(SurfSideID)%State(1,p,q).EQ.1)     THEN
!        SampWall(SurfSideID)%State(2,p,q)                   = e_kin
        SampWall(SurfSideID)%State(3,p,q)                   = e_kin
        SampWall(SurfSideID)%State(4,p,q)                   = e_kin
    END IF

!   All subsequent impacts
    delta                                                   = e_kin - SampWall(SurfSideID)%State(2,p,q)
!   Update mean
    SampWall(SurfSideID)%State(2,p,q)                       = SampWall(SurfSideID)%State(2,p,q) + delta /                          &
                                                              SampWall(SurfSideID)%State(1,p,q)
!    Find min/max of distribution
    IF (e_kin.LT.SampWall(SurfSideID)%State(3,p,q)) THEN
        SampWall(SurfSideID)%State(3,p,q)               = e_kin
    END IF
    IF (e_kin.GT.SampWall(SurfSideID)%State(4,p,q)) THEN
        SampWall(SurfSideID)%State(4,p,q)               = e_kin
    END IF
!   delta2 = newValue - mean
    delta2                                                  = e_kin - SampWall(SurfSideID)%State(2,p,q)
!   M2 = M2 + delta * delta2
    SampWall(SurfSideID)%State(5,p,q)                       = SampWall(SurfSideID)%State(5,p,q) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWall(SurfSideID)%State(6,p,q)                       = SampWall(SurfSideID)%State(5,p,q)/SampWall(SurfSideID)%State(1,p,q)
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN

!   Record first impact, otherwise min will be frozen at zero
    IF (SampWall(SurfSideID)%State(1+nShift,p,q).EQ.1)     THEN
!        SampWall(SurfSideID)%State(2+nShift,p,q)            = e_kin
        SampWall(SurfSideID)%State(3+nShift,p,q)            = e_kin
        SampWall(SurfSideID)%State(4+nShift,p,q)            = e_kin
    END IF

!   All subsequent impacts
    delta                                                   = e_kin - SampWall(SurfSideID)%State(2+nShift,p,q)
!   Update mean
    SampWall(SurfSideID)%State(2+nShift,p,q)                = SampWall(SurfSideID)%State(2+nShift,p,q) + delta /                   &
                                                              SampWall(SurfSideID)%State(1+nShift,p,q)
!   Find min/max of distribution
    IF (e_kin.LT.SampWall(SurfSideID)%State(3+nShift,p,q)) THEN
        SampWall(SurfSideID)%State(3+nShift,p,q)        = e_kin
    END IF
    IF (e_kin.GT.SampWall(SurfSideID)%State(4+nShift,p,q)) THEN
        SampWall(SurfSideID)%State(4+nShift,p,q)        = e_kin
    END IF
!   delta2 = newValue - mean
    delta2                                                  = e_kin - SampWall(SurfSideID)%State(2+nShift,p,q)
!   M2 = M2 + delta * delta2
    SampWall(SurfSideID)%State(5+nShift,p,q)                = SampWall(SurfSideID)%State(5+nShift,p,q) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWall(SurfSideID)%State(6+nShift,p,q)                = SampWall(SurfSideID)%State(5+nShift,p,q) /                           &
                                                              SampWall(SurfSideID)%State(1+nShift,p,q)
END IF

!===================================================================================================================================
!---- 7. - 11 / Impact angle (mean, min, max, M2, variance)
!   Record first impact, otherwise min will be frozen at zero
    IF (SampWall(SurfSideID)%State(1,p,q).EQ.1)     THEN
!        SampWall(SurfSideID)%State(7,p,q)                   = PartFaceAngle
        SampWall(SurfSideID)%State(8,p,q)                   = PartFaceAngle
        SampWall(SurfSideID)%State(9,p,q)                   = PartFaceAngle
    END IF

!   All subsequent impacts
    delta                                                   = PartFaceAngle - SampWall(SurfSideID)%State(7,p,q)
!   Update mean
    SampWall(SurfSideID)%State(7,p,q)                       = SampWall(SurfSideID)%State(7,p,q) + delta /                          &
                                                              SampWall(SurfSideID)%State(1,p,q)
!    Find min/max of distribution
    IF (PartFaceAngle.LT.SampWall(SurfSideID)%State(8,p,q)) THEN
        SampWall(SurfSideID)%State(8,p,q)                = PartFaceAngle
    END IF
    IF (PartFaceAngle.GT.SampWall(SurfSideID)%State(9,p,q)) THEN
        SampWall(SurfSideID)%State(9,p,q)                = PartFaceAngle
    END IF
!   delta2 = newValue - mean
    delta2                                                  = PartFaceAngle - SampWall(SurfSideID)%State(7,p,q)
!   M2 = M2 + delta * delta2
    SampWall(SurfSideID)%State(10,p,q)                      = SampWall(SurfSideID)%State(10,p,q) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWall(SurfSideID)%State(11,p,q)                      = SampWall(SurfSideID)%State(10,p,q)/SampWall(SurfSideID)%State(1,p,q)
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN

!   Record first impact, otherwise min will be frozen at zero
    IF (SampWall(SurfSideID)%State(1+nShift,p,q).EQ.1)     THEN
!        SampWall(SurfSideID)%State(7+nShift,p,q)            = PartFaceAngle
        SampWall(SurfSideID)%State(8+nShift,p,q)            = PartFaceAngle
        SampWall(SurfSideID)%State(9+nShift,p,q)            = PartFaceAngle
    END IF

!   All subsequent impacts
    delta                                                   = PartFaceAngle - SampWall(SurfSideID)%State(7+nShift,p,q)
!   Update mean
    SampWall(SurfSideID)%State(7+nShift,p,q)                = SampWall(SurfSideID)%State(7+nShift,p,q) + delta /                   &
                                                              SampWall(SurfSideID)%State(1+nShift,p,q)
!   Find min/max of distribution
    IF (PartFaceAngle.LT.SampWall(SurfSideID)%State(8+nShift,p,q)) THEN
        SampWall(SurfSideID)%State(8+nShift,p,q)        = PartFaceAngle
    END IF
    IF (PartFaceAngle.GT.SampWall(SurfSideID)%State(9+nShift,p,q)) THEN
        SampWall(SurfSideID)%State(9+nShift,p,q)        = PartFaceAngle
    END IF
!   delta2 = newValue - mean
    delta2                                                  = PartFaceAngle - SampWall(SurfSideID)%State(7+nShift,p,q)
!   M2 = M2 + delta * delta2
    SampWall(SurfSideID)%State(10+nShift,p,q)               = SampWall(SurfSideID)%State(10+nShift,p,q) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWall(SurfSideID)%State(11+nShift,p,q)               = SampWall(SurfSideID)%State(10+nShift,p,q) /                          &
                                                              SampWall(SurfSideID)%State(1+nShift,p,q)
END IF

!===================================================================================================================================
!---- 12 - 14 / Sampling Current Forces at walls
    SampWall(SurfSideID)%State(12,p,q)= SampWall(SurfSideID)%State(12,p,q) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(1))
    SampWall(SurfSideID)%State(13,p,q)= SampWall(SurfSideID)%State(13,p,q) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(2))
    SampWall(SurfSideID)%State(14,p,q)= SampWall(SurfSideID)%State(14,p,q) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(3))
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN
    SampWall(SurfSideID)%State(12+nShift,p,q)= SampWall(SurfSideID)%State(12+nShift,p,q) + Species(PartSpecies(PartID))%MassIC     &
                                               * (v_old(1))
    SampWall(SurfSideID)%State(13+nShift,p,q)= SampWall(SurfSideID)%State(13+nShift,p,q) + Species(PartSpecies(PartID))%MassIC     &
                                               * (v_old(2))
    SampWall(SurfSideID)%State(14+nShift,p,q)= SampWall(SurfSideID)%State(14+nShift,p,q) + Species(PartSpecies(PartID))%MassIC     &
                                               * (v_old(3))
END IF

!===================================================================================================================================
!---- 15 - 17 / Sampling Average Forces at walls
    SampWall(SurfSideID)%State(15,p,q)= SampWall(SurfSideID)%State(15,p,q) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(1))
    SampWall(SurfSideID)%State(16,p,q)= SampWall(SurfSideID)%State(16,p,q) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(2))
    SampWall(SurfSideID)%State(17,p,q)= SampWall(SurfSideID)%State(17,p,q) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(3))
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN
    SampWall(SurfSideID)%State(15+nShift,p,q)= SampWall(SurfSideID)%State(15+nShift,p,q) + Species(PartSpecies(PartID))%MassIC     &
                                        * (v_old(1))
    SampWall(SurfSideID)%State(16+nShift,p,q)= SampWall(SurfSideID)%State(16+nShift,p,q) + Species(PartSpecies(PartID))%MassIC     &
                                        * (v_old(2))
    SampWall(SurfSideID)%State(17+nShift,p,q)= SampWall(SurfSideID)%State(17+nShift,p,q) + Species(PartSpecies(PartID))%MassIC     &
                                        * (v_old(3))
END IF

END SUBROUTINE RecordParticleBoundarySampling


FUNCTION PARTSWITCHELEMENT(xi,eta,locSideID,SideID,ElemID)
!===================================================================================================================================
! particle moves through face and switches element
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToElemAndSide
USE MOD_Mesh_Vars,              ONLY:MortarType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: locSideID, SideID,ElemID
REAL,INTENT(IN)     :: xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: PARTSWITCHELEMENT
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! move particle to new element
!     Type 1               Type 2              Type3
!      eta                  eta                 eta
!       ^                    ^                   ^
!       |                    |                   |
!   +---+---+            +---+---+           +---+---+
!   | 3 | 4 |            |   2   |           |   |   |
!   +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!   | 1 | 2 |            |   1   |           |   |   |
!   +---+---+            +---+---+           +---+---+

SELECT CASE(MortarType(1,SideID))
CASE(1)
  IF(Xi.GT.0.)THEN
    IF(Eta.GT.0.)THEN
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(4  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(4+4,locSideID,ElemID)
    ELSE
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
    END IF
  ELSE
    IF(Eta.GT.0.)THEN
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(3  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(3+4,locSideID,ElemID)
    ELSE
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
    END IF
  END IF
CASE(2)
  IF(Eta.GT.0.)THEN
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
  ELSE
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
  END IF
CASE(3)
  IF(Xi.LE.0.)THEN
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
  ELSE
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
  END IF
CASE DEFAULT ! normal side OR small mortar side
  PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
  PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
END SELECT

END FUNCTION PARTSWITCHELEMENT


END MODULE MOD_Particle_Boundary_Condition
