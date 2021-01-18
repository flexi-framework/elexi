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

INTERFACE GetBoundaryInteractionAuxBC
  MODULE PROCEDURE GetBoundaryInteractionAuxBC
END INTERFACE

PUBLIC :: GetBoundaryInteraction
PUBLIC :: GetBoundaryInteractionAuxBC
!===================================================================================================================================

CONTAINS

SUBROUTINE GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,ElemID,crossedBC&
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
USE MOD_Globals                    ,ONLY: ABORT
USE MOD_Particle_Surfaces          ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Tracking_Vars     ,ONLY: TrackingMethod
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools        ,ONLY: GetCNSideID
USE MOD_Particle_Boundary_Sampling ,ONLY: SideErosion
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Surfaces_vars     ,ONLY: SideNormVec,SideType
USE MOD_Part_Operations            ,ONLY: RemoveParticle
!USE MOD_Mesh_Vars                  ,ONLY: BC
#if CODE_ANALYZE
USE MOD_Globals                    ,ONLY: myRank,UNIT_stdout
USE MOD_Mesh_Vars                  ,ONLY: NGeo
USE MOD_Particle_Mesh_Tools        ,ONLY: GetCNElemID
USE MOD_Particle_Surfaces_Vars     ,ONLY: BezierControlPoints3D
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,flip
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN),OPTIONAL          :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: CNSideID
REAL                                 :: n_loc(1:3)
#if CODE_ANALYZE
REAL                                 :: v1(3),v2(3)
#endif /* CODE_ANALYZE */
!===================================================================================================================================
crossedBC = .FALSE.

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(TrackingMethod)
CASE(REFMAPPING,TRACING)
  ! set BCSideID for normal vector calculation call with (curvi-)linear side description
  ! IF (TrackingMethod.EQ.RefMapping) BCSideID=PartBCSideList(SideID)
  CNSideID = GetCNSideID(SideID)

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc = SideNormVec(1:3,CNSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(  nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
  END SELECT

  ! Flip side orientation if not on the master side
  IF(flip.NE.0) n_loc=-n_loc

#if CODE_ANALYZE
  ! check if normal vector points outwards
  v1 = 0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)  &
           + BezierControlPoints3D(:,NGeo,0   ,SideID)  &
           + BezierControlPoints3D(:,0   ,NGeo,SideID)  &
           + BezierControlPoints3D(:,NGeo,NGeo,SideID))
  v2 = v1  - ElemBaryNGeo(:,GetCNElemID(ElemID))

  IF (DOT_PRODUCT(v2,n_loc).LT.0) THEN
    IPWRITE(UNIT_stdout,*) 'Obtained wrong side orientation from flip. flip:',flip,'PartID:',iPart
    IPWRITE(UNIT_stdout,*) 'n_loc (flip)', n_loc,'n_loc (estimated):',v2
    CALL ABORT(__STAMP__,'SideID',SideID)
  END IF
#endif /* CODE_ANALYZE */

    ! Inserted particles are "pushed" inside the domain and registered as passing through the BC side. If they are very close to the
    ! boundary (first if) than the normal vector is compared with the trajectory. If the particle is entering the domain from outside
    ! it was inserted during surface flux and this routine shall not performed.
    ! Comparing the normal vector with the particle trajectory, if the particle trajectory is pointing inside the domain
    IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN

  CASE(TRIATRACKING)
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
END SELECT

! required for refmapping and tracing, optional for triatracking
crossedBC = .TRUE.

IF (.NOT. ALLOCATED(PartBound%TargetBoundCond)) &
  CALL ABORT(__STAMP__,' ERROR: PartBound not allocated!.')

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,SideID)))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Sample on surface if requested
  CALL SideErosion(PartTrajectory,n_loc,xi,eta,iPart,SideID,alpha)
  CALL RemoveParticle(iPart,alpha=alpha)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
!  IF (PDM%ParticleInside(iPart)) THEN
    ! simple reflection
    SELECT CASE(PartBound%WallModel(SideInfo_Shared(SIDE_BCID,SideID)))
      ! perfectly reflecting, specular re-emission
      CASE('perfRef')
         CALL  PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,n_loc,opt_Symmetry=.FALSE.)
      ! reflection using coefficient of restitution (CoR)
      CASE('coeffRes')
         CALL  DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,n_loc                               &
                                ,WallCoeffModel=PartBound%WallCoeffModel(SideInfo_Shared(SIDE_BCID,SideID)))
      CASE DEFAULT
          CALL ABORT(__STAMP__, ' No or invalid particle wall model given. Please update Part-Boundary[$]-WallModel.')
    END SELECT
!  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,iPart,SideID,ElemID) ! opt_reflected is peri-moved
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)')
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(10) !PartBound%SymmetryBC
!-----------------------------------------------------------------------------------------------------------------------------------
  ! For particles, symmetry equals perfect reflection, also flip forces direction
  CALL  PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,n_loc,opt_Symmetry=.TRUE.)

!CASE(100) !PartBound%AnalyzeBC
!!-----------------------------------------------------------------------------------------------------------------------------------
!  CALL  SideAnalysis(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,ElemID &
!                    ,opt_crossed=crossedBC)

CASE DEFAULT
  CALL abort(__STAMP__,' ERROR: PartBound not associated!. (unknown case)',999,999.)
END SELECT

END SUBROUTINE GetBoundaryInteraction


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
USE MOD_Particle_Vars              ,ONLY: PDM
USE MOD_Particle_Boundary_Sampling ,ONLY: SideErosion
USE MOD_Particle_Boundary_Vars     ,ONLY: PartAuxBC
USE MOD_Particle_Boundary_Vars     ,ONLY: AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol
USE MOD_Particle_Vars              ,ONLY: LastPartPos
USE MOD_Part_Operations            ,ONLY: RemoveParticle
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
REAL                                 :: n_loc(1:3)
REAL                                 :: intersec(3),r_vec(3),axis(3),cos2inv
!===================================================================================================================================

! Reset flag
crossedBC    =.FALSE.
SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
CASE ('plane')
  n_loc = AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec
CASE ('cylinder')
  intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory
  r_vec = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis) ) )
  IF (.NOT.AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE ('cone')
  intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory
  r_vec = AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
  cos2inv = 1./COS(AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle)**2
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis)*cos2inv ) )
  IF (.NOT.AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE ('parabol')
  intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory
  r_vec = AuxBC_parabol(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_parabol(AuxBCMap(AuxBCIdx))%axis
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*(DOT_PRODUCT(intersec-r_vec,axis)+0.5*AuxBC_parabol(AuxBCMap(AuxBCIdx))%zfac) ) )
  IF (.NOT.AuxBC_parabol(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE DEFAULT
  CALL abort(&
    __STAMP__&
    ,'AuxBC does not exist')
END SELECT
IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
  crossedBC=.FALSE.
  !RETURN
  CALL abort(&
    __STAMP__&
    ,'Error in GetBoundaryInteractionAuxBC: Particle coming from outside!')
ELSE IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0.)  THEN
  crossedBC=.TRUE.
ELSE
  CALL abort(&
    __STAMP__&
    ,'Error in GetBoundaryInteractionAuxBC: n_vec is perpendicular to PartTrajectory for AuxBC',AuxBCIdx)
END IF
! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartAuxBC%TargetBoundCond(AuxBCIdx))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartAuxBC%OpenBC
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Sampling on aux surfaces currently not supported
  ! CALL SideErosion(PartTrajectory,n_loc,xi,eta,iPart,SideID,alpha)
  CALL RemoveParticle(iPart,alpha=alpha)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartAuxBC%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (PDM%ParticleInside(iPart)) THEN
    ! simple reflection, auxBC can not have diffuse reflection
    ! perfectly reflecting, specular re-emission
      CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,n_loc=n_loc, &
             opt_Symmetry=.FALSE.,AuxBCIdx=AuxBCIdx)
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
  CALL abort(__STAMP__,' ERROR: AuxBC bound not associated!. (unknown case)')
END SELECT

END SUBROUTINE GetBoundaryInteractionAuxBC


SUBROUTINE PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_Loc, &
                             opt_Symmetry,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_ErosionPoints              ,ONLY: RecordErosionPoint
USE MOD_ErosionPoints_Vars         ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Sampling ,ONLY: SideErosion
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound,PartAuxBC
USE MOD_Particle_Boundary_Vars     ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Boundary_Vars     ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars     ,ONLY: LowVeloRemove
USE MOD_Particle_Mesh_Vars         ,ONLY: SideInfo_Shared
USE MOD_Particle_Globals
USE MOD_Particle_Vars              ,ONLY: PartState,LastPartPos,Species,PartSpecies,PartReflCount
USE MOD_Particle_Vars              ,ONLY: Pt_temp,PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3),lengthPartTrajectory,alpha
REAL,INTENT(IN)                   :: xi, eta
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3),WallVelo(3)
LOGICAL                           :: Symmetry,IsAuxBC
REAL                              :: PartFaceAngle,PartFaceAngle_old
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
  WallVelo=PartAuxBC%WallVelo(1:3,AuxBCIdx)

! Normal BC
ELSE
  ! Get wall velo and BCID
  WallVelo  = PartBound%WallVelo(1:3,SideInfo_Shared(SIDE_BCID,SideID))

  IF(PRESENT(opt_Symmetry)) THEN
    Symmetry = opt_Symmetry
  ELSE
    Symmetry = .FALSE.
  END IF
END IF !IsAuxBC

! Make sure we have the old values safe
v_old                = PartState(4:6,PartID)
IF (doParticleImpactTrack) THEN
  PartFaceAngle_old  = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
END IF

! Sample on boundary
IF ((.NOT.IsAuxBC) .AND. WriteMacroSurfaceValues) THEN
  CALL SideErosion(PartTrajectory,n_loc,xi,eta,PartID,SideID,alpha)
END IF !.NOT.IsAuxBC

! Update particle velocity
PartState(4:6,PartID) = PartState(4:6,PartID)-2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo

! set particle position on face
!--> first move the particle to the boundary
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
!--> flip trajectory and move the remainder of the particle push
PartTrajectory(1:3)     = PartTrajectory(1:3)-2.*DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc
PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*(lengthPartTrajectory - alpha)

! compute moved particle || rest of movement
PartTrajectory          = PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory    = SQRT(PartTrajectory(1)*PartTrajectory(1)            &
                          +    PartTrajectory(2)*PartTrajectory(2)            &
                          +    PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory          = PartTrajectory/lengthPartTrajectory

! Recording of individual particle impacts
IF (doParticleImpactTrack) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

  CALL RecordErosionPoint(BCSideID        = SideInfo_Shared(SIDE_BCID,SideID) &
                         ,PartID          = PartID                        &
                         ,PartFaceAngle   = PartFaceAngle                 &
                         ,v_old           = v_old                         &
                         ,PartFaceAngle_old =PartFaceAngle_old            &
                         ,PartReflCount   = PartReflCount(PartID)         &
                         ,alpha           = alpha                         &
                         ,n_loc           = n_loc)
END IF

! Increase reflection counter
IF (doParticleReflectionTrack) PartReflCount(PartID) = PartReflCount(PartID) + 1

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


SUBROUTINE DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_loc,AuxBCIdx,WallCoeffModel)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_ErosionPoints            ,ONLY: RecordErosionPoint
USE MOD_ErosionPoints_Vars       ,ONLY: doParticleImpactTrack
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound,PartAuxBC
USE MOD_Particle_Boundary_Sampling,ONLY:SideErosion
USE MOD_Particle_Boundary_Vars   ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars   ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Boundary_Vars   ,ONLY: LowVeloRemove
USE MOD_Particle_Mesh_Vars       ,ONLY: SideInfo_Shared
USE MOD_Particle_Surfaces        ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars            ,ONLY: PartState,LastPartPos,Species,PartSpecies,PartReflCount
USE MOD_Particle_Vars            ,ONLY: PDM
! Bons particle rebound model
USE MOD_Particle_Interpolation_Vars,ONLY:FieldAtParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID,SideID
CHARACTER(LEN=255),INTENT(IN)     :: WallCoeffModel
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3),WallVelo(3)
INTEGER                           :: locBCID
LOGICAL                           :: IsAuxBC
REAL                              :: PartFaceAngle,PartFaceAngleDeg,PartFaceAngle_old
REAL                              :: PartTrajectoryTang(3),PartTrajectoryNorm(3)
REAL                              :: v_magnitude,v_norm(3),v_tang(3)
REAL                              :: intersecRemain
REAL                              :: eps_n, eps_t
REAL                              :: tang1(1:3),tang2(1:3)
! Bons particle rebound model
REAL                              :: E_eff
REAL                              :: Vol,d,w,w_crit,sigma_y
!===================================================================================================================================

! check if reflected on AuxBC
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF

! reflected on AuxBC
IF (IsAuxBC) THEN
  WallVelo   = PartAuxBC%WallVelo(1:3,AuxBCIdx)
ELSE
  ! additional states
  locBCID    = SideInfo_Shared(SIDE_BCID,SideID)
  ! get BC values
  WallVelo   = PartBound%WallVelo(1:3,locBCID)
END IF !IsAuxBC
CALL OrthoNormVec(n_loc,tang1,tang2)

! Make sure we have the old velocity safe
v_old   = PartState(4:6,PartID)

! Sample on boundary
IF ((.NOT.IsAuxBC) .AND. WriteMacroSurfaceValues) THEN
  CALL SideErosion(PartTrajectory,n_loc,xi,eta,PartID,SideID,alpha)
END IF

! Calculate wall normal and tangential velocity, impact angle
v_norm  = DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc
v_tang  = PartState(4:6,PartID) - v_norm
PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

SELECT CASE(WallCoeffModel)
  !=================================================================================================================================
  ! Grant, G., Tabakoff, W.. "Erosion prediction in turbomachinery resulting from environmental solid particles."
  ! / Journal of Aircraft 12.5 (1975): 471-478.
  ! >> Ignored variance (random fluctuations in reflection values) and ONLY took mean values!
  !=================================================================================================================================
  CASE('Grant1975')
    ! Transfer from radians to degree
    PartFaceAngleDeg = PartFaceAngle * 180/PI

    ! Compute cubic polynomials for coefficient of restitution
    eps_n = 9.9288e-01 - 3.0708e-02 * PartFaceAngleDeg  + 4.7389e-04 * PartFaceAngleDeg**2. - 2.5845e-06 * PartFaceAngleDeg**3.
    eps_t = 9.8723e-01 - 3.0354e-02 * PartFaceAngleDeg  + 7.1913e-04 * PartFaceAngleDeg**2. - 4.2979e-06 * PartFaceAngleDeg**3.

  !=================================================================================================================================
  ! Tabaoff, W.; Wakeman, T.: Basic Erosion Investigation in Small Turbomachinery. / Cincinnati Univ. OH, 1981
  ! >> Ignored variance (random fluctuations in reflection values) and ONLY took mean values!
  !=================================================================================================================================
  CASE('Tabakoff1981')
    ! Transfer from radians to degree
    PartFaceAngleDeg = PartFaceAngle * 180/PI

    ! Compute cubic polynomials for coefficient of restitution
    eps_n = 1.    - 0.0211 * PartFaceAngleDeg  + 0.000228 * PartFaceAngleDeg**2. - 0.000000876 * PartFaceAngleDeg**3.
    eps_t = 0.953                              - 0.000446 * PartFaceAngleDeg**2. + 0.00000648  * PartFaceAngleDeg**3.

  !===================================================================================================================================
  ! Bons, J., Prenter, R., Whitaker, S.: A Simple Physics-Based Model for Particle Rebound and Deposition in Turbomachinery
  ! / J. Turbomach 139(8), 2017
  !===================================================================================================================================
  CASE('Bons2017')
    ! Assume spherical particles for now
    Vol     = Species(PartSpecies(PartID))%MassIC/Species(PartSpecies(PartID))%DensityIC
    d       = (6.*Vol/PI)**(1./3.)

    ! Find composite elastic modulus
    E_eff   = ((1. - Species(PartSpecies(PartID))%PoissonIC**2.)/Species(PartSpecies(PartID))%YoungIC +        &
               (1. - PartBound%Poisson(SideInfo_Shared(SIDE_BCID,SideID))            **2.)/PartBound%Young(SideInfo_Shared(SIDE_BCID,SideID))               )**(-1.)

    ! Calculate deformation of cylindrical model particle
    w       = SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) * (8*Species(PartSpecies(PartID))%MassIC &
                                                          /(E_eff*3*d)   )**0.5

    ! Find critical deformation
    sigma_y = Species(PartSpecies(PartID))%YieldCoeff
    w_crit  = sigma_y * 2./3. * d / E_eff

    ! Normal coefficient of restitution
    IF (w .GT. w_crit) THEN
      eps_n = 1./SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) *sigma_y*(1/ (Species(PartSpecies(PartID))%DensityIC &
                                                                          *  E_eff   ))**0.5
    ELSE
      eps_n = 1.
    END IF

    ! Tangential coefficient of restitution not considering adhesion
    !> Assume change in density from last particle position to wall position to be negligible
    ! Original relation by Barker, B., Casaday, B., Shankara, P., Ameri, A., and Bons, J. P., 2013.
    !> Cosine term added by Bons, J., Prenter, R., Whitaker, S., 2017.
    eps_t   = 1. - 0.3 / SQRT(DOT_PRODUCT(v_tang(1:3),v_tang(1:3)))  * &
                       SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) * (eps_n+1)*COS(PartFaceAngle)**2.


  !===================================================================================================================================
  ! Whitaker, S., Bons, J.: An improved particle impact model by accounting for rate of strain and stochastic rebound, 2018
  !===================================================================================================================================
  CASE('Whitaker2018')
    ! Assume spherical particles for now
    Vol     = Species(PartSpecies(PartID))%MassIC/Species(PartSpecies(PartID))%DensityIC
    d       = (6.*Vol/PI)**(1./3.)

    ! Find composite elastic modulus
    E_eff   = ((1. - Species(PartSpecies(PartID))%PoissonIC**2.)/Species(PartSpecies(PartID))%YoungIC +        &
               (1. - PartBound%Poisson(SideInfo_Shared(SIDE_BCID,SideID))            **2.)/PartBound%Young(SideInfo_Shared(SIDE_BCID,SideID))               )**(-1.)

    ! Calculate deformation of cylindrical model particle
    w       = SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) * (8*Species(PartSpecies(PartID))%MassIC &
                                                          /(E_eff*3*d)   )**0.5

    ! Find critical deformation
    sigma_y = Species(PartSpecies(PartID))%YieldCoeff*SQRT(DOT_PRODUCT(v_old(1:3),v_old(1:3)))
    w_crit  = sigma_y * 2./3. * d / E_eff

    ! Normal coefficient of restitution
    IF (w .GT. w_crit) THEN
      eps_n = 1./SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) *sigma_y*(1/ (Species(PartSpecies(PartID))%DensityIC &
                                                                          *  E_eff   ))**0.5
    ELSE
      eps_n = 1.
    END IF

    ! Tangential coefficient of restitution not considering adhesion
    !> Assume change in density from last particle position to wall position to be negligible
    ! Original relation by Barker, B., Casaday, B., Shankara, P., Ameri, A., and Bons, J. P., 2013.
    !> Cosine term added by Bons, J., Prenter, R., Whitaker, S., 2017.
    eps_t   = 1. - 0.63 / SQRT(DOT_PRODUCT(v_tang(1:3),v_tang(1:3)))  * &
                       SQRT(DOT_PRODUCT(v_norm(1:3),v_norm(1:3))) * (eps_n+1)*COS(PartFaceAngle)**2.

  !=================================================================================================================================
  ! Fong, W.; Amili, O.; Coletti, F.: Velocity and spatial distribution of intertial particles in a turbulent channel flow
  ! / J. FluidMech 872, 2019
  !=================================================================================================================================
  CASE('Fong2019')
    ! Reuse YieldCoeff to modify the normal velocity, keep tangential velocity
    eps_t   = 1.
    eps_n   = PartBound%CoR(SideInfo_Shared(SIDE_BCID,SideID))

  CASE DEFAULT
      CALL abort(__STAMP__, ' No particle wall coefficients given. This should not happen.')

END SELECT

! Make sure we have the old values safe
IF (doParticleImpactTrack) THEN
    PartFaceAngle_old  = PartFaceAngle
END IF

! set particle position on face
!--> first move the particle to the boundary
#if CODE_ANALYZE
WRITE(UNIT_stdout,'(110("-"))')
WRITE(UNIT_stdout,'(A,I1)') '     | Diffusive reflection on BC: ', SideInfo_Shared(SIDE_BCID,SideID)
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | LastPartPos:                ',LastPartPos(1,PartID),LastPartPos(2,PartID),LastPartPos(3,PartID)
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | PartTrajectory:             ',PartTrajectory(1),PartTrajectory(2),PartTrajectory(3)
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | Velocity:                   ',PartState(4,PartID),PartState(5,PartID),PartState(6,PartID)
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16)')          '     | alpha,lengthPartTrajectory: ',alpha,lengthPartTrajectory
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | Intersection:               ', LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16)')          '     | CoR (normal/tangential):    ',eps_n,eps_t
#endif

LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
!--> flip trajectory and move the remainder of the particle push
!===================================================================================================================================
! > modified with coefficients of restitution
! > WARNING: this only works with LSERK. Check and abort if not fulfilled
!IF (.NOT.(TimeDiscType.EQ.'LSERKW2').AND..NOT.(TimeDiscType.EQ.'LSWERKW3'))               &
!    CALL ABORT(__STAMP__,                                                                 &
!    'Time discretization '//TRIM(TimeDiscType)//' is incompatible with current implementation of coefficients of restitution.')

! Calculate wall normal and tangential velocity components after impact, rescale to uniform length
PartTrajectoryNorm(1:3) = eps_n*(DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc)
PartTrajectoryTang(1:3) = eps_t*(PartTrajectory(1:3) - DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc)
PartTrajectory(1:3)     = PartTrajectoryTang(1:3) - PartTrajectoryNorm(1:3)
PartTrajectory          = PartTrajectory/SQRT(SUM(PartTrajectory**2.))

! Rescale the remainder to the new length
intersecRemain = SQRT(eps_n*eps_n + eps_t*eps_t)/SQRT(2.) * (lengthPartTrajectory - alpha)

! Compute moved particle || rest of movement. PartTrajectory has already been updated
PartState(1:3,PartID) = LastPartPos(1:3,PartID) + intersecRemain*PartTrajectory(1:3)

! Compute new particle velocity, modified with coefficents of restitution
PartState(4:6,PartID)= eps_t * v_tang - eps_n * v_norm + WallVelo

#if CODE_ANALYZE
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | PartTrajectory (CoR)        ',PartTrajectory(1),PartTrajectory(2),PartTrajectory(3)
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | alpha (CoR):                ',interSecRemain
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | NewPartPos:                 ',PartState(1,PartID),PartState(2,PartID),PartState(3,PartID)
WRITE(UNIT_stdout,'(A,E27.16,x,E27.16,x,E27.16)') '     | Velocity (CoR):             ',PartState(4,PartID),PartState(5,PartID),PartState(6,PartID)
#endif

IF (doParticleImpactTrack) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

  CALL RecordErosionPoint(BCSideID        = SideInfo_Shared(SIDE_BCID,SideID),            &
                          PartID          = PartID,                                       &
                          PartFaceAngle   = PartFaceAngle,                                &
                          v_old           = v_old,                                        &
                          PartFaceAngle_old =PartFaceAngle_old,                           &
                          PartReflCount   = PartReflCount(PartID),                        &
                          alpha           = alpha,                                        &
                          n_loc           = n_loc)
END IF

! increase reflection counter
IF (doParticleReflectionTrack) PartReflCount(PartID) = PartReflCount(PartID) + 1

! Correction for Runge-Kutta. Reflect or delete force history / acceleration
! Coefficients of Restitution cause non-differentiable jump in velocity. Always erase force history and reconstruction in timedisc during push
PDM%IsNewPart(PartID)=.TRUE.

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


SUBROUTINE PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,PartID,SideID,ElemID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the periodic shift in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: BoundaryType
!USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Vars          ,ONLY: PartState,LastPartPos!,PEM
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
#if CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
INTEGER,INTENT(IN)                :: PartID,SideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: PVID
!INTEGER                           :: moved(2),locSideID
!===================================================================================================================================

PVID = BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA)

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
ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)


!locSideID = PartSideToElem(S2E_LOC_SIDE_ID,SideID)
!Moved     = PARTSWITCHELEMENT(xi,eta,locSideID,SideID,ElemID)
!ElemID    = Moved(1)
!#if USE_MPI
!IF(ElemID.EQ.-1)THEN
!  CALL abort(&
!__STAMP__&
!,' Halo region to small. Neighbor element is missing!')
!END IF
!#endif /*USE_MPI*/
!IF (TrackingMethod.EQ.REFMAPPING) PEM%LastElement(PartID) = 0

END SUBROUTINE PeriodicBC


!FUNCTION PARTSWITCHELEMENT(xi,eta,locSideID,SideID,ElemID)
!!===================================================================================================================================
!! particle moves through face and switches element
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Particle_Mesh_Vars ,ONLY: PartElemToElemAndSide
!USE MOD_Mesh_Vars          ,ONLY: MortarType
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!INTEGER,INTENT(IN)  :: locSideID, SideID,ElemID
!REAL,INTENT(IN)     :: xi,eta
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!INTEGER,DIMENSION(2) :: PARTSWITCHELEMENT
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!===================================================================================================================================
!
!! move particle to new element
!!     Type 1               Type 2              Type3
!!      eta                  eta                 eta
!!       ^                    ^                   ^
!!       |                    |                   |
!!   +---+---+            +---+---+           +---+---+
!!   | 3 | 4 |            |   2   |           |   |   |
!!   +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!!   | 1 | 2 |            |   1   |           |   |   |
!!   +---+---+            +---+---+           +---+---+
!
!CALL ABORT(__STAMP__,'Not yet implemented for new halo region')
!
!SELECT CASE(MortarType(1,SideID))
!CASE(1)
!  IF(Xi.GT.0.)THEN
!    IF(Eta.GT.0.)THEN
!      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(4  ,locSideID,ElemID)
!      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(4+4,locSideID,ElemID)
!    ELSE
!      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
!      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
!    END IF
!  ELSE
!    IF(Eta.GT.0.)THEN
!      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(3  ,locSideID,ElemID)
!      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(3+4,locSideID,ElemID)
!    ELSE
!      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
!      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
!    END IF
!  END IF
!CASE(2)
!  IF(Eta.GT.0.)THEN
!    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
!    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
!  ELSE
!    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
!    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
!  END IF
!CASE(3)
!  IF(Xi.LE.0.)THEN
!    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
!    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
!  ELSE
!    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
!    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
!  END IF
!CASE DEFAULT ! normal side OR small mortar side
!  PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
!  PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
!END SELECT
!
!END FUNCTION PARTSWITCHELEMENT

END MODULE MOD_Particle_Boundary_Condition
