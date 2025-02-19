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

#define SIGMOID(x)         REAL(1/(1+EXP(-x)),4)
#define SILU(x,beta)       REAL(x/(1+EXP(-beta*x)),4)
#define LOGNORM(x,maxi)    REAL(LOG(ABS(x)+1.)/LOG(maxi+1.),4)
#define LOGNORMINV(x,maxi) REAL((maxi+1.)**x-1.,4)

!===================================================================================================================================
!> Determines how particles interact with a given boundary condition
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Condition
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: GetBoundaryInteraction
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
USE MOD_Globals                    ,ONLY: Abort,FLEXITIME
USE MOD_Globals_Vars               ,ONLY: PI
! USE MOD_Mesh_Vars                  ,ONLY: BC
USE MOD_Mesh_Vars                  ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Sampling ,ONLY: RecordParticleBoundaryImpact
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound,doParticleImpactTrack
USE MOD_Particle_Boundary_Tracking ,ONLY: StoreBoundaryParticleProperties
USE MOD_Particle_Mesh_Tools        ,ONLY: GetCNSideID
USE MOD_Particle_Surfaces          ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars     ,ONLY: SideNormVec,SideType
USE MOD_Particle_Tracking_Vars     ,ONLY: TrackingMethod
USE MOD_Particle_Vars              ,ONLY: PartState,PartReflCount
USE MOD_Particle_Operations        ,ONLY: RemoveParticle
#if CODE_ANALYZE
USE MOD_Globals                    ,ONLY: myRank,UNIT_stdOut
USE MOD_Mesh_Vars                  ,ONLY: NGeo
USE MOD_Particle_Mesh_Tools        ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars         ,ONLY: ElemBaryNGeo
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
REAL                                 :: n_loc(1:3),PartFaceAngle
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
    IPWRITE(UNIT_stdOut,'(I0,A,I1,A,I0)')         'Obtained wrong side orientation from flip. flip: ',flip,', PartID: ',iPart
    IPWRITE(UNIT_stdOut,'(I0,A,3F12.6,A,3F12.6)') 'n_loc (flip)', n_loc,'n_loc (estimated):',v2
    CALL Abort(__STAMP__,'SideID',SideID)
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
  CALL Abort(__STAMP__,' ERROR: PartBound not allocated!.')

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,SideID)))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Sample on surface if requested
  CALL RecordParticleBoundaryImpact(PartTrajectory,n_loc,xi,eta,iPart,SideID)! ,alpha)
  ! Recording of individual particle impacts
  IF (doParticleImpactTrack) THEN
    PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

    CALL StoreBoundaryParticleProperties(BCSideID        = SideInfo_Shared(SIDE_BCID,SideID) &
                                        ,PartID          = iPart                             &
                                        ,PartFaceAngle   = PartFaceAngle                     &
                                        ,v_old           = PartState(PART_VELV,iPart)        &
                                        ,PartFaceAngle_old =PartFaceAngle                    &
                                        ,PartReflCount   = PartReflCount(iPart)              &
                                        ,alpha           = alpha                             &
                                        ,dp_old          = PartState(PART_DIAM,iPart)        &
#if USE_PARTROT
                                        ,rot_old         = PartState(PART_AMOMV,iPart)       &
#endif
                                        )
  END IF
  CALL RemoveParticle(iPart,BCID=SideInfo_Shared(SIDE_BCID,SideID),alpha=alpha)
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
          CALL Abort(__STAMP__, ' No or invalid particle wall model given. Please update Part-Boundary[$]-WallModel.')
    END SELECT
!  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,iPart,SideID,ElemID) ! opt_reflected is peri-moved
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL Abort(__STAMP__,' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)')
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
  CALL Abort(__STAMP__,' ERROR: PartBound not associated!. (unknown case)',999,999.)
END SELECT

END SUBROUTINE GetBoundaryInteraction


SUBROUTINE PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_Loc, &
                             opt_Symmetry)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars               ,ONLY: PI
USE MOD_Mesh_Vars                  ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Sampling ,ONLY: RecordParticleBoundaryImpact
USE MOD_Particle_Boundary_Tracking ,ONLY: StoreBoundaryParticleProperties
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound
USE MOD_Particle_Boundary_Vars     ,ONLY: doParticleReflectionTrack,doParticleImpactTrack
USE MOD_Particle_Boundary_Vars     ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars     ,ONLY: LowVeloRemove
USE MOD_Particle_Globals
USE MOD_Particle_Operations        ,ONLY: RemoveParticle
USE MOD_Particle_Vars              ,ONLY: PartState,LastPartPos,Species,PartSpecies,PartReflCount
USE MOD_Particle_Vars              ,ONLY: Pt_temp,PDM
USE MOD_Utils                      ,ONLY: ALMOSTZERO
#if USE_PARTROT
USE MOD_Mathtools                  ,ONLY: CROSS
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3),lengthPartTrajectory,alpha
REAL,INTENT(IN)                   :: xi, eta
REAL,INTENT(INOUT)                :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3),WallVelo(3)
LOGICAL                           :: Symmetry
REAL                              :: PartFaceAngle,PartFaceAngle_old
REAL                              :: v_magnitude
INTEGER                           :: locBCID
#if USE_PARTROT
REAL                              :: rot_old(1:3)
#endif
!===================================================================================================================================

locBCID   = SideInfo_Shared(SIDE_BCID,SideID)
! Get wall velo and BCID
WallVelo  = PartBound%WallVelo(1:3,locBCID)
IF (PRESENT(opt_Symmetry)) THEN; Symmetry  = opt_Symmetry
ELSE                           ; Symmetry  = .FALSE.
END IF

! Rough wall modelling
IF (PartBound%doRoughWallModelling(locBCID).AND.Species(PartSpecies(PartID))%doRoughWallModelling) THEN
  n_loc = RoughWall(n_loc,locBCID,PartTrajectory)
END IF

! Make sure we have the old values safe
v_old   = PartState(PART_VELV,PartID)
#if USE_PARTROT
rot_old = PartState(PART_AMOMV,PartID)
#endif
IF (doParticleImpactTrack) THEN
  PartFaceAngle_old  = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))
END IF

! Sample on boundary
IF (WriteMacroSurfaceValues) THEN
  CALL RecordParticleBoundaryImpact(PartTrajectory,n_loc,xi,eta,PartID,SideID)! ,alpha)
END IF ! WriteMacroSurfaceValues

! Update particle velocity
PartState(PART_VELV,PartID) = PartState(PART_VELV,PartID) - 2.*DOT_PRODUCT(PartState(PART_VELV,PartID),n_loc)*n_loc + WallVelo

! set particle position on face
!--> first move the particle to the boundary
LastPartPos(1:3,PartID)       = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
!--> flip trajectory and move the remainder of the particle push
PartTrajectory(1:3)           = PartTrajectory(1:3) - 2.*DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc
PartState(PART_POSV,PartID)   = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*(lengthPartTrajectory - alpha)

! compute moved particle || rest of movement
PartTrajectory          = PartState(PART_POSV,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory    = VECNORM(PartTrajectory)
PartTrajectory          = PartTrajectory/lengthPartTrajectory

#if USE_PARTROT
! rotation: I_2*w_2-I_1*w_1 = - r x J , J=m_2*v_2-m_1*v_1, r=d/2*n_loc
! for a constant particle volume: I_2=I_1, m_1=m_2
PartState(PART_AMOMV,PartID) = rot_old(1:3) - 5/PartState(PART_DIAM,PartID)*CROSS(n_loc,(PartState(PART_VELV,PartID)-v_old))
#endif

! Recording of individual particle impacts
IF (doParticleImpactTrack) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

  CALL StoreBoundaryParticleProperties(BCSideID        = SideInfo_Shared(SIDE_BCID,SideID) &
                                      ,PartID          = PartID                            &
                                      ,PartFaceAngle   = PartFaceAngle                     &
                                      ,v_old           = v_old                             &
                                      ,PartFaceAngle_old =PartFaceAngle_old                &
                                      ,PartReflCount   = PartReflCount(PartID)             &
                                      ,alpha           = alpha                             &
                                      ,dp_old          = PartState(PART_DIAM,PartID)       &
#if USE_PARTROT
                                      ,rot_old         = rot_old                           &
#endif
                                      )
END IF

! Increase reflection counter
IF (doParticleReflectionTrack) PartReflCount(PartID) = PartReflCount(PartID) + 1

! Correction for Runge-Kutta. Reflect or delete force history / acceleration
IF (.NOT.ALMOSTZERO(DOT_PRODUCT(WallVelo,WallVelo))) THEN
  ! Reconstruction in timedisc during push
  PDM%IsNewPart(PartID) = .TRUE.
ELSE
  Pt_temp(1:3,PartID) = Pt_temp(1:3,PartID) - 2.*DOT_PRODUCT(Pt_temp(1:3,PartID),n_loc)*n_loc
  ! Reflect also force history for symmetry
  IF (Symmetry) THEN
    Pt_temp(4:6,PartID) = Pt_temp(4:6,PartID) - 2.*DOT_PRODUCT(Pt_temp(4:6,PartID),n_loc)*n_loc
#if USE_PARTROT
    Pt_temp(7:9,PartID) = Pt_temp(7:9,PartID) - 2.*DOT_PRODUCT(Pt_temp(7:9,PartID),n_loc)*n_loc
#endif
  ! Produces best result compared to analytical solution in place capacitor ...
  ELSE
    Pt_temp(4:PP_nVarPartRHS+3,PartID) = 0.
  END IF
END IF

! Remove sliding low velocity particles
IF (LowVeloRemove) THEN
  v_magnitude   = VECNORM(PartState(PART_VELV,PartID))

  IF ((Species(PartSpecies(PartID))%LowVeloThreshold.NE.0).AND.&
    (v_magnitude.LT.Species(PartSpecies(PartID))%LowVeloThreshold)) THEN
    Species(PartSpecies(PartID))%LowVeloCounter = Species(PartSpecies(PartID))%LowVeloCounter + 1
    CALL RemoveParticle(PartID)
    IPWRITE(UNIT_stdOut,'(I0,A,F12.6)') ' Low velocity particle removed after impact. Velocity after reflection:', v_magnitude
  END IF
END IF

END SUBROUTINE PerfectReflection


SUBROUTINE DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_loc,WallCoeffModel)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars               ,ONLY: PI
USE MOD_Mesh_Vars                  ,ONLY: SideInfo_Shared
USE MOD_Particle_Globals           ,ONLY: VECNORM,OrthoNormVec
USE MOD_Particle_Boundary_Sampling ,ONLY: RecordParticleBoundaryImpact
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound,PartBoundANN
USE MOD_Particle_Boundary_Vars     ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars     ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Boundary_Vars     ,ONLY: LowVeloRemove
USE MOD_Particle_Boundary_Vars     ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Tracking ,ONLY: StoreBoundaryParticleProperties
USE MOD_Particle_Fracture          ,ONLY: DiffuseReflectionFracture
USE MOD_Particle_Operations        ,ONLY: RemoveParticle
USE MOD_Particle_Vars              ,ONLY: PartState,LastPartPos,Species,PartSpecies,PartReflCount
USE MOD_Particle_Vars              ,ONLY: PDM
#if USE_PARTROT
USE MOD_Mathtools                  ,ONLY: CROSS
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
REAL,INTENT(INOUT)                :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID,SideID
CHARACTER(LEN=255),INTENT(IN)     :: WallCoeffModel
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3),dp_old,WallVelo(3)
INTEGER                           :: locBCID
REAL                              :: PartFaceAngle,PartFaceAngleDeg,PartFaceAngle_old
REAL                              :: PartTrajectoryTang1(3),PartTrajectoryTang2(3),PartTrajectoryNorm(3)
REAL                              :: v_magnitude,v_norm(3),v_tang1(3),v_tang2(3)
REAL                              :: intersecRemain
REAL                              :: eps_n, eps_t1, eps_t2
REAL                              :: tang1(1:3),tang2(1:3)
#if USE_PARTROT
REAL                              :: rot_old(1:3)
#endif
! Bons particle rebound model
REAL                              :: w,w_crit,sigma_y,E_eff
! RebANN
INTEGER                           :: k
REAL(4)                           :: randnum(6),xin(8)
REAL                              :: dp1,dp2,ekin_1,ekin_2,S
!===================================================================================================================================

! additional states
locBCID    = SideInfo_Shared(SIDE_BCID,SideID)
! get BC values
WallVelo   = PartBound%WallVelo(1:3,locBCID)

! Make sure we have the old velocity/dp/rot safe
v_old   = PartState(PART_VELV,PartID)
dp_old  = PartState(PART_DIAM,PartID)
#if USE_PARTROT
rot_old = PartState(PART_AMOMV,PartID)
#endif

! Compute tangential vectors
CALL OrthoNormVec(n_loc,tang1,tang2)

! Sample on boundary
IF (WriteMacroSurfaceValues) THEN
  CALL RecordParticleBoundaryImpact(PartTrajectory,n_loc,xi,eta,PartID,SideID)! ,alpha)
END IF ! WriteMacroSurfaceValues

! Move particle to wall
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha

! Calculate wall normal and tangential velocities, impact angle
v_norm   = DOT_PRODUCT(PartState(PART_VELV,PartID),n_loc)*n_loc
v_tang1  = DOT_PRODUCT(PartState(PART_VELV,PartID),tang1)*tang1
v_tang2  = DOT_PRODUCT(PartState(PART_VELV,PartID),tang2)*tang2
PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

! Just mirror the velocity in tang2 direction if no eps_t2 is available
eps_t2 = 1.

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
    eps_n  = 9.9288e-01 - 3.0708e-02 * PartFaceAngleDeg  + 4.7389e-04 * PartFaceAngleDeg**2. - 2.5845e-06 * PartFaceAngleDeg**3.
    eps_t1 = 9.8723e-01 - 3.0354e-02 * PartFaceAngleDeg  + 7.1913e-04 * PartFaceAngleDeg**2. - 4.2979e-06 * PartFaceAngleDeg**3.

  !=================================================================================================================================
  ! Tabaoff, W.; Wakeman, T.: Basic Erosion Investigation in Small Turbomachinery. / Cincinnati Univ. OH, 1981
  ! >> Ignored variance (random fluctuations in reflection values) and ONLY took mean values!
  !=================================================================================================================================
  CASE('Tabakoff1981')
    ! Transfer from radians to degree
    PartFaceAngleDeg = PartFaceAngle * 180/PI

    ! Compute cubic polynomials for coefficient of restitution
    eps_n  = 1.    - 0.0211 * PartFaceAngleDeg  + 0.000228 * PartFaceAngleDeg**2. - 0.000000876 * PartFaceAngleDeg**3.
    eps_t1 = 0.953                              - 0.000446 * PartFaceAngleDeg**2. + 0.00000648  * PartFaceAngleDeg**3.

  !===================================================================================================================================
  ! Bons, J., Prenter, R., Whitaker, S.: A Simple Physics-Based Model for Particle Rebound and Deposition in Turbomachinery
  ! / J. Turbomach 139(8), 2017
  !===================================================================================================================================
  CASE('Bons2017')
    ! Find composite elastic modulus
    E_eff   = ((1. - Species(PartSpecies(PartID))%PoissonIC**2.)/Species(PartSpecies(PartID))%YoungIC +        &
               (1. - PartBound%Poisson(SideInfo_Shared(SIDE_BCID,SideID))**2.)/PartBound%Young(SideInfo_Shared(SIDE_BCID,SideID)))**(-1.)

    v_magnitude = VECNORM(v_norm(1:3))

    ! Calculate deformation of cylindrical model particle
    w       = v_magnitude * (8*MASS_SPHERE(Species(PartSpecies(PartID))%DensityIC,PartState(PART_DIAM,PartID)) &
              /(E_eff*3*PartState(PART_DIAM,PartID)))**0.5

    ! Find critical deformation
    sigma_y = Species(PartSpecies(PartID))%YieldCoeff
    w_crit  = sigma_y * 2./3. * PartState(PART_DIAM,PartID) / E_eff

    ! Normal coefficient of restitution
    eps_n = MERGE(1./v_magnitude*sigma_y*(1./(Species(PartSpecies(PartID))%DensityIC * E_eff))**0.5, 1., w .GT. w_crit)

    ! Tangential coefficient of restitution not considering adhesion
    !> Assume change in density from last particle position to wall position to be negligible
    ! Original relation by Barker, B., Casaday, B., Shankara, P., Ameri, A., and Bons, J. P., 2013.
    !> Cosine term added by Bons, J., Prenter, R., Whitaker, S., 2017.
    eps_t1   = 1. - PartBound%FricCoeff(SideInfo_Shared(SIDE_BCID,SideID)) / VECNORM(v_tang1(1:3))  * &
                           v_magnitude * (eps_n+1)*COS(PartFaceAngle)**2.


  !===================================================================================================================================
  ! Whitaker, S., Bons, J.: An improved particle impact model by accounting for rate of strain and stochastic rebound, 2018
  !===================================================================================================================================
  CASE('Whitaker2018')
    ! Find composite elastic modulus
    E_eff   = ((1. - Species(PartSpecies(PartID))%PoissonIC**2.)/Species(PartSpecies(PartID))%YoungIC +        &
               (1. - PartBound%Poisson(SideInfo_Shared(SIDE_BCID,SideID))**2.)/PartBound%Young(SideInfo_Shared(SIDE_BCID,SideID)))**(-1.)

    v_magnitude = VECNORM(v_norm(1:3))

    ! Calculate deformation of cylindrical model particle
    w       = v_magnitude * (8*MASS_SPHERE(Species(PartSpecies(PartID))%DensityIC,PartState(PART_DIAM,PartID)) &
              / (E_eff*3*PartState(PART_DIAM,PartID)))**0.5

    ! Find critical deformation
    sigma_y = Species(PartSpecies(PartID))%Whitaker_a*VECNORM(v_old(1:3))*1e6
    w_crit  = sigma_y * 2./3. * PartState(PART_DIAM,PartID) / E_eff

    ! Normal coefficient of restitution
    eps_n = MERGE(1./v_magnitude*sigma_y*(1./(Species(PartSpecies(PartID))%DensityIC * E_eff))**0.5, 1., w .GT. w_crit)

    ! Tangential coefficient of restitution not considering adhesion
    !> Assume change in density from last particle position to wall position to be negligible
    ! Original relation by Barker, B., Casaday, B., Shankara, P., Ameri, A., and Bons, J. P., 2013.
    !> Cosine term added by Bons, J., Prenter, R., Whitaker, S., 2017.
    eps_t1   = 1. - PartBound%FricCoeff(SideInfo_Shared(SIDE_BCID,SideID)) / VECNORM(v_tang1(1:3)) * &
                           v_magnitude * (eps_n+1)*COS(PartFaceAngle)**2.

  !=================================================================================================================================
  ! Fong, W.; Amili, O.; Coletti, F.: Velocity and spatial distribution of inertial particles in a turbulent channel flow
  ! / J. FluidMech 872, 2019
  !=================================================================================================================================
  CASE('Fong2019')
    ! Reuse YieldCoeff to modify the normal velocity, keep tangential velocity
    eps_t1   = 1.
    eps_n    = PartBound%CoR(SideInfo_Shared(SIDE_BCID,SideID))

  !=================================================================================================================================
  ! Rebound 2D ANN valid for v_{air} \in [150,350] [m/s]
  !=================================================================================================================================
  CASE('RebANN')
    v_magnitude = SQRT(PartState(PART_VEL1,PartID)**2+PartState(PART_VEL3,PartID)**2)
    ! Input with normalization
    xin(1:3) = (/LOGNORM(PartFaceAngle,PartBoundANN%max_in(1)),&
      LOGNORM(v_magnitude,PartBoundANN%max_in(2)), LOGNORM(PartState(PART_DIAM,PartID)*1e6,PartBoundANN%max_in(3))/)
    ! Check if kinetic energy of rebounding particle is lower
    PartBoundANN%output(2) = REAL(v_magnitude + 1,4)
    DO WHILE (PartBoundANN%output(2) > v_magnitude)
      CALL RANDOM_NUMBER(randnum(1:2))
      xin(4) = randnum(1)
      xin(5) = REAL(0.9*randnum(2)+0.1,4)

      CALL CallANN(5,2,xin)

      k = k + 1
      IF (k.GT.5) EXIT
    END DO
    ! Calculate coefficents of restitution
    v_magnitude = PartBoundANN%output(2)
    eps_n  = v_magnitude * SIN(PartBoundANN%output(1)) / DOT_PRODUCT(PartState(PART_VELV,PartID),n_loc)
    eps_t1 = v_magnitude * COS(PartBoundANN%output(1)) / DOT_PRODUCT(PartState(PART_VELV,PartID),tang1)

  !=================================================================================================================================
  ! Fracture 2D ANN valid for v_{air} \in [150,350] [m/s]
  !=================================================================================================================================
  CASE('FracANN')
    v_magnitude = SQRT(PartState(PART_VEL1,PartID)**2+PartState(PART_VEL3,PartID)**2)

    ! Input with normalization
    xin(1:3) = (/LOGNORM(PartFaceAngle,PartBoundANN%max_in(1)),&
      LOGNORM(v_magnitude,PartBoundANN%max_in(2)),&
      LOGNORM(PartState(PART_DIAM,PartID)*1e6,PartBoundANN%max_in(3))/)

    ! Check if kinetic energy of rebounding particle is lower
    ekin_1 = ENERGY_KINETIC(Species(PartSpecies(PartID))%DensityIC,PartState(PART_DIAM,PartID),PartState(PART_VELV,PartID))
    ekin_2 = ekin_1 + 1

    ! fracture probability
    S = 1.-1./(1.+(90.e-6/PartState(PART_DIAM,PartID) - 1.) * ABS(DOT_PRODUCT(PartState(PART_VELV,PartID),n_loc)/350.)**0.45)

    CALL RANDOM_NUMBER(randnum(6))
    DO WHILE (ekin_2 > ekin_1)
      CALL RANDOM_NUMBER(randnum(1:5))

      IF (S .GT. randnum(6)) THEN
        ! particle fractures
        xin(4:8) = (/randnum(1), REAL(0.9*randnum(2)+0.1,4), randnum(3), randnum(4), REAL(0.9*randnum(5)+0.1,4)/)
      ELSE
        xin(4:8) = (/randnum(1), REAL(0.9*randnum(2)+0.1,4), randnum(3), REAL(0.,4), REAL(0.,4)/)
      END IF

      CALL CallANN(8,6,xin)

      IF (S .GT. randnum(6)) THEN
        ! particle fractures
        dp1 = PartBoundANN%output(3)*PartState(PART_DIAM,PartID)
      ELSE
        dp1 = PartState(PART_DIAM,PartID)
      END IF

      IF ((dp1 .LT. PartState(PART_DIAM,PartID))) THEN
        dp2 = VOL_SPHERE_INV((VOL_SPHERE(dp_old)-VOL_SPHERE(dp1)))
        ekin_2 = ENERGY_KINETIC(Species(PartSpecies(PartID))%DensityIC,dp1,(/PartBoundANN%output(2)/)) +&
          ENERGY_KINETIC(Species(PartSpecies(PartID))%DensityIC,dp2,(/PartBoundANN%output(5)/))
      ELSE
        dp2 = 0.
        ekin_2 = ENERGY_KINETIC(Species(PartSpecies(PartID))%DensityIC,dp1,(/PartBoundANN%output(2)/))
      END IF

      k = k + 1
      IF (k.GT.10) THEN; dp2=0.; EXIT; END IF
    END DO

    ! Insert new particle if particle splits
    IF (dp2 .GT. 0.) THEN
      ! Update particle diameter
      PartState(PART_DIAM,PartID) = dp1

      ! Calculate coefficents of restitution
      eps_n  = PartBoundANN%output(5) * SIN(PartBoundANN%output(4)) / NORM2(v_norm(1:3))
      eps_t1 = PartBoundANN%output(5) * COS(PartBoundANN%output(4)) / NORM2(v_tang1(1:3))

      ! Insert new particle
      CALL DiffuseReflectionFracture(PartTrajectory,lengthPartTrajectory,alpha,PartID,SideID,n_loc,tang1,tang2,&
                                     eps_n,eps_t1,eps_t2,dp_old,WallVelo)
    END IF

    ! Calculate coefficents of restitution
    eps_n  = PartBoundANN%output(2) * SIN(PartBoundANN%output(1)) / NORM2(v_norm(1:3))
    eps_t1 = PartBoundANN%output(2) * COS(PartBoundANN%output(1)) / NORM2(v_tang1(1:3))

  CASE DEFAULT
      CALL Abort(__STAMP__, ' No particle wall coefficients given. This should not happen.')

END SELECT

! Make sure we have the old values safe
IF (doParticleImpactTrack) THEN
    PartFaceAngle_old  = PartFaceAngle
END IF

! set particle position on face
!--> first move the particle to the boundary
#if CODE_ANALYZE
WRITE(UNIT_stdOut,'(110("-"))')
WRITE(UNIT_stdOut,'(A,I1)')                       '     | Diffusive reflection on BC: ',SideInfo_Shared(SIDE_BCID,SideID)
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | LastPartPos:                ',LastPartPos(1,PartID),LastPartPos(2,PartID),LastPartPos(3,PartID)
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | PartTrajectory:             ',PartTrajectory(1),PartTrajectory(2),PartTrajectory(3)
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | Velocity:                   ',PartState(4,PartID),PartState(5,PartID),PartState(6,PartID)
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16)')          '     | alpha,lengthPartTrajectory: ',alpha,lengthPartTrajectory
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | Intersection:               ',LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16)')          '     | CoR (normal/tangential):    ',eps_n,eps_t1,eps_t2
#endif

!--> flip trajectory and move the remainder of the particle push
!===================================================================================================================================
! > modified with coefficients of restitution
! > WARNING: this only works with LSERK. Check and Abort if not fulfilled
!IF (.NOT.(TimeDiscType.EQ.'LSERKW2').AND..NOT.(TimeDiscType.EQ.'LSWERKW3'))               &
!    CALL Abort(__STAMP__,                                                                 &
!    'Time discretization '//TRIM(TimeDiscType)//' is incompatible with current implementation of coefficients of restitution.')

! Calculate wall normal and tangential velocity components after impact, rescale to uniform length
PartTrajectoryNorm (1:3) = eps_n  * (DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc)
PartTrajectoryTang1(1:3) = eps_t1 * (DOT_PRODUCT(PartTrajectory(1:3),tang1)*tang1)
PartTrajectoryTang2(1:3) = eps_t2 * (DOT_PRODUCT(PartTrajectory(1:3),tang2)*tang2)
PartTrajectory(1:3)      = PartTrajectoryTang1(1:3) + PartTrajectoryTang2(1:3) - PartTrajectoryNorm(1:3)

! Rescale the remainder to the new length
intersecRemain = (lengthPartTrajectory - alpha)
intersecRemain = SQRT(eps_n*eps_n + eps_t1*eps_t1 + eps_t2*eps_t2)/SQRT(3.) * intersecRemain
! Compute the remainder of the new particle trajectory
PartTrajectory = intersecRemain * PartTrajectory/VECNORM(PartTrajectory)

! Compute moved particle || rest of movement. PartTrajectory has already been updated
PartState(PART_POSV,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)

! Compute new particle velocity, modified with coefficents of restitution
PartState(PART_VELV,PartID) = eps_t1 * v_tang1 + eps_t2 * v_tang2 - eps_n * v_norm + WallVelo

! compute moved particle || rest of movement
lengthPartTrajectory  = VECNORM(PartTrajectory)
PartTrajectory        = PartTrajectory/lengthPartTrajectory

#if USE_PARTROT
! rotation: I_2*w_2-I_1*w_1 = - r x J , J=m_2*v_2-m_1*v_1, r=d/2*n_loc
! for a constant particle volume: I_2=I_1, m_1=m_2
PartState(PART_AMOMV,PartID) = (dp_old/PartState(PART_DIAM,PartID))**5*rot_old - &
                                0.5*dp_old*Species(PartSpecies(PartID))%DensityIC*PI/6*&
                                CROSS(n_loc,(PartState(PART_DIAM,PartID)**3*PartState(PART_VELV,PartID)-v_old*dp_old**3))
#endif

#if CODE_ANALYZE
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | PartTrajectory (CoR)        ',PartTrajectory(1),PartTrajectory(2),PartTrajectory(3)
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | alpha (CoR):                ',interSecRemain
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | NewPartPos:                 ',PartState(1,PartID),PartState(2,PartID),PartState(3,PartID)
WRITE(UNIT_stdOut,'(A,E27.16,x,E27.16,x,E27.16)') '     | Velocity (CoR):             ',PartState(4,PartID),PartState(5,PartID),PartState(6,PartID)
#endif

IF (doParticleImpactTrack) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

  CALL StoreBoundaryParticleProperties(BCSideID           = SideInfo_Shared(SIDE_BCID,SideID)             &
                                      ,PartID             = PartID                                        &
                                      ,PartFaceAngle      = PartFaceAngle                                 &
                                      ,v_old              = v_old                                         &
                                      ,PartFaceAngle_old  = PartFaceAngle_old                             &
                                      ,PartReflCount      = PartReflCount(PartID)                         &
                                      ,alpha              = alpha                                         &
                                      ,dp_old             = dp_old                                        &
#if USE_PARTROT
                                      ,rot_old            = rot_old                                       &
#endif
                                      )
END IF

! increase reflection counter
IF (doParticleReflectionTrack) PartReflCount(PartID) = PartReflCount(PartID) + 1

! Correction for Runge-Kutta. Reflect or delete force history / acceleration
! Coefficients of Restitution cause non-differentiable jump in velocity. Always erase force history and reconstruction in timedisc during push
PDM%IsNewPart(PartID) = .TRUE.

! Remove sliding low velocity particles
IF (LowVeloRemove) THEN
  v_magnitude = VECNORM(PartState(PART_VELV,PartID))

  IF ((Species(PartSpecies(PartID))%LowVeloThreshold.NE.0).AND.&
    (v_magnitude.LT.Species(PartSpecies(PartID))%LowVeloThreshold)) THEN
    Species(PartSpecies(PartID))%LowVeloCounter = Species(PartSpecies(PartID))%LowVeloCounter + 1
    CALL RemoveParticle(PartID)
    IPWRITE(UNIT_stdOut,'(I0,A,F12.6)') ' Low velocity particle removed after impact. Velocity after reflection:', v_magnitude
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
! USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem,TensorProductInterpolation
USE MOD_Mesh_Vars              ,ONLY: BoundaryType
USE MOD_Mesh_Vars              ,ONLY: SideInfo_Shared
! USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Vars          ,ONLY: PartState,LastPartPos
! USE MOD_Particle_Vars          ,ONLY: PartPosRef
! USE MOD_Particle_Mesh_Tools    ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
! USE MOD_Particle_Mesh_Vars     ,ONLY: XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo_Shared,ElemEpsOneCell
! USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod,DoPeriodicCheck,DoPeriodicFix
#if CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
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
! INTEGER                           :: CNElemID
! REAL                              :: PartStateOld(1:3),Displacement(1:3)
!===================================================================================================================================

PVID = GEO%PeriodicMapping(ABS(BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA)))
! Fix the sign for negative BC_ALPHA
PVID = SIGN(PVID,BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA))

#if CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdOut,'(I0,A,I0)')      ' PeriodicBC:      ', SideInfo_Shared(SIDE_BCID,SideID)
    IPWRITE(UNIT_stdOut,'(I0,A,I0)')      ' PartID:          ', PartID
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,G0))') ' ParticlePosition: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,G0))') ' LastPartPos:      ',LastPartPos(1:3,PartID)
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
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,G0))') ' ParticlePosition-pp: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdOut,'(I0,A,3(X,G0))') ' LastPartPo-pp:       ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/

! refmapping and tracing
! move particle from old element to new element
ElemID   = SideInfo_Shared(SIDE_NBELEMID,SideID)

! ! check if particle is correctly inside the new element
! IF (DoPeriodicCheck) THEN
!   SELECT CASE(TrackingMethod)
!     CASE (REFMAPPING,TRACING)
!       CNElemID = GetCNElemID(ElemID)
!       CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef(1:3,PartID),ElemID)
!       ! Position outside of tolerance
!       ! IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).GT.ElemEpsOneCell(CNElemID)) THEN
!       IF(MAXVAL(ABS(PartPosRef(1:3,PartID))).GT.1.) THEN
!         IPWRITE(UNIT_stdOut,'(I0,A)') ' Tolerance Issue with periodic element '
!         IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' PartPos                ', PartState(1:3,PartID)
!         IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' xi                     ', PartPosRef(1:3,PartID)
!         IPWRITE(UNIT_stdOut,'(I0,A,X,E15.8)')    ' epsOneCell             ', ElemEpsOneCell(CNElemID)
!         ! Fix particle position or Abort
!         IF (DoPeriodicFix) THEN
!           ! Save the old PartState
!           PartStateOld = PartState(1:3,PartID)
!           ! Restrict particle position in reference space
!           PartPosRef(1:3,PartID) = MIN(PartPosRef(1:3,PartID), 1.)
!           PartPosRef(1:3,PartID) = MAX(PartPosRef(1:3,PartID),-1.)
!           ! Get the physical coordinates that correspond to the reference coordinates
!           CALL TensorProductInterpolation( PartPosRef(1:3,PartID)                           &
!                                          , 3                                                &
!                                          , NGeo                                             &
!                                          , XiCL_NGeo                                        &
!                                          , wBaryCL_NGeo                                     &
!                                          , XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,ElemID) &
!                                          , PartState(1:3,PartID))
!           ! Calculate required displacement and adjust LastPartPos
!           Displacement = PartState(1:3,PartID) - PartStateOld
!           LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + Displacement
!           IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' Particle relocated by  ', Displacement
!         ELSE
!           CALL Abort(__STAMP__,'Particle not inside of element, PartID'
!       END IF

!     CASE(TRIATRACKING)
!       ! Currently not available

!     END SELECT
! END IF

END SUBROUTINE PeriodicBC


FUNCTION RoughWall(n_in,locBCID,PartTrajectory) RESULT (n_out)
!----------------------------------------------------------------------------------------------------------------------------------!
! Rough wall modelling without multiple rebounds, where the roughness is drawn from a Gaussian distribution with a mean of zero and
! a standard deviation equal to an assumed or experimental wall roughness in radii.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals_Vars               ,ONLY: PI
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBound
USE MOD_Particle_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: n_in(1:3)
INTEGER,INTENT(IN)                :: locBCID
REAL,INTENT(IN)                   :: PartTrajectory(1:3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL                              :: n_out(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: RandNum(4),random_angle(2),angle,cross_vector(3)
!===================================================================================================================================

cross_vector = CROSSNORM(n_in,PartTrajectory)

! Compute angles
angle = 0.5*PI-ACOS(DOT_PRODUCT(PartTrajectory,n_in))

CALL RANDOM_NUMBER(RandNum)
! Scale to maximal range
RandNum(1) = RandNum(1)*(angle*0.5+0.2*PI) - angle*0.5
RandNum(3) = RandNum(3)*(angle*0.5+0.2*PI) - angle*0.5

CALL RoughnessPDF(RandNum(1:2),PartBound%RoughMeanIC(locBCID),PartBound%RoughVarianceIC(locBCID),angle)
CALL RoughnessPDF(RandNum(3:4),PartBound%RoughMeanIC(locBCID),PartBound%RoughVarianceIC(locBCID),angle)

random_angle(1) = RandNum(1)
random_angle(2) = RandNum(3)

! Modifiy the angle between particle trajectory and normal vector
n_out = n_in*COS(random_angle(1))+(/cross_vector(2)*n_in(3)   - cross_vector(3)*n_in(2),&
                                    cross_vector(3)*n_in(1)   - cross_vector(1)*n_in(3),&
                                    cross_vector(1)*n_in(2)   - cross_vector(2)*n_in(1)/)*SIN(random_angle(1))+&
                                    cross_vector*DOT_PRODUCT(cross_vector,n_in)*(1.-COS(random_angle(1)))
! Rotation of the new normal vector around the old normal vector according to Rodrigues' rotation formula:
n_out = n_out*COS(random_angle(2))+(/n_in(2)*n_out(3)   - n_in(3)*n_out(2),&
                                     n_in(3)*n_out(1)   - n_in(1)*n_out(3),&
                                     n_in(1)*n_out(2)   - n_in(2)*n_out(1)/)*SIN(random_angle(2))+&
                                     n_in(:)*DOT_PRODUCT(n_in(:),n_out)*(1.-COS(random_angle(2)))
END FUNCTION RoughWall


RECURSIVE SUBROUTINE RoughnessPDF(x,MeanIC,VarianceIC,angle)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals_Vars,               ONLY: PI
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: x(2)
REAL,INTENT(IN)                   :: MeanIC
REAL,INTENT(IN)                   :: VarianceIC
REAL,INTENT(IN)                   :: angle
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL                              :: EffectivePDF
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: GaussianNormal
!-----------------------------------------------------------------------------------------------------------------------------------

GaussianNormal = 1./SQRT(2*PI*VarianceIC) * EXP(-0.5*((x(1)-MeanIC)/VarianceIC)**2)
EffectivePDF = SIN(angle+x(1))/SIN(angle)*GaussianNormal
! Draw a sample via the acceptance-rejection method
IF (x(2) .GT. EffectivePDF)THEN
  CALL RANDOM_NUMBER(x)
  CALL RoughnessPDF(x,MeanIC,VarianceIC,angle)
END IF

END SUBROUTINE RoughnessPDF


SUBROUTINE CallANN(nin,nout,xin)
!----------------------------------------------------------------------------------------------------------------------------------!
! Calls the rebound MLP
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars     ,ONLY: PartBoundANN
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: nin
INTEGER,INTENT(IN)               :: nout
REAL(SP),INTENT(IN)              :: xin(nin)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iLayer
REAL                             :: x(PartBoundANN%hgs)
!===================================================================================================================================
x(:)     = 0.
x(1:nin) = xin
! Hidden layers
DO iLayer=1,PartBoundANN%nLayer
  PartBoundANN%output(:) = &
    SILU((MATMUL(x,PartBoundANN%w(:,:,iLayer)) + PartBoundANN%b(:,iLayer)),PartBoundANN%beta(:,iLayer))
  x = PartBoundANN%output(:)
END DO
! Output layer
iLayer = PartBoundANN%nLayer+1
PartBoundANN%output(:)   = SIGMOID((MATMUL(PartBoundANN%output(:),PartBoundANN%w(:,:,iLayer)) + PartBoundANN%b(:,iLayer)))
PartBoundANN%output(1:nout) = LOGNORMINV(PartBoundANN%output(1:nout), PartBoundANN%max_out(1:nout))
END SUBROUTINE CallANN

END MODULE MOD_Particle_Boundary_Condition
