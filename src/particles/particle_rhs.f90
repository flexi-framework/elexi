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
#include "eos.h"
#include "particle.h"

#define PHI(x) 1+0.5*x+1./6*x**2

!===================================================================================================================================
! Subroutine to compute the particle right hand side, therefore the acceleration due to the Lorentz-force with
! respect to the Lorentz factor
!===================================================================================================================================
MODULE MOD_part_RHS
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CalcPartRHS
  MODULE PROCEDURE CalcPartRHS
END INTERFACE

#if USE_EXTEND_RHS || USE_FAXEN_CORR
INTERFACE extRHS
  MODULE PROCEDURE extRHS
END INTERFACE
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */

INTERFACE CalcSourcePart
  MODULE PROCEDURE CalcSourcePart
END INTERFACE

PUBLIC :: CalcPartRHS, InitRHS
#if USE_EXTEND_RHS || USE_FAXEN_CORR
PUBLIC :: extRHS
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */
PUBLIC :: CalcSourcePart
!==================================================================================================================================

CONTAINS

SUBROUTINE InitRHS(drag_factor, FD_Pointer)
!===================================================================================================================================
! Init RHS
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars, ONLY:type_F
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: drag_factor
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(type_F),INTENT(INOUT)  :: FD_Pointer
!-----------------------------------------------------------------------------------------------------------------------------------

SELECT CASE(drag_factor)
  CASE(DF_PART_SCHILLER)
    FD_Pointer%op => DF_SchillerAndNaumann
  CASE(DF_PART_PUTNAM)
    FD_Pointer%op => DF_Putnam
  CASE(DF_PART_HAIDER)
    FD_Pointer%op => DF_Haider
  CASE(DF_PART_HOELZER)
    FD_Pointer%op => DF_Hoelzer
  CASE(DF_PART_LOTH)
    FD_Pointer%op => DF_Loth
  CASE(DF_PART_GANSER)
    FD_Pointer%op => DF_Ganser
END SELECT

END SUBROUTINE InitRHS

SUBROUTINE CalcPartRHS(&
#if USE_BASSETFORCE || ANALYZE_RHS
  t,dt,iStage)
#else
  )
#endif /* USE_BASSETFORCE || ANALYZE_RHS */
!===================================================================================================================================
! Computes the acceleration from the drag force with respect to the species data and velocity
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Interpolation_Vars,  ONLY: FieldAtParticle
#if USE_EXTEND_RHS || USE_FAXEN_CORR
USE MOD_Particle_Interpolation_Vars,  ONLY: GradAtParticle
#endif
USE MOD_Particle_Vars,                ONLY: PDM, Pt
!#if ANALYZE_RHS
!USE MOD_Particle_Vars,                ONLY: tWriteRHS,dtWriteRHS
!#endif /* ANALYZE_RHS */
! #if USE_RW
! USE MOD_Particle_RandomWalk_Vars,     ONLY: RWTime
#if USE_BASSETFORCE
USE MOD_Particle_Vars,                ONLY: bIter,N_Basset
#endif /* USE_BASSETFORCE */
! #endif
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if USE_BASSETFORCE || ANALYZE_RHS
REAL,INTENT(IN)                  :: t
REAL,INTENT(IN)                  :: dt
INTEGER,INTENT(IN),OPTIONAL      :: iStage
#endif /* USE_BASSETFORCE || ANALYZE_RHS */
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE
INTEGER                          :: iPart
!===================================================================================================================================

! Drag force
Pt(:,1:PDM%ParticleVecLength)=0.

DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
!#if USE_RW
!    ! Do not change the particle velocity if RW is working in full Euler mode
!    !> Ideally, this should use tStage. But one cannot start a RK without the first stage and it does not make a difference for Euler
!    IF (RWTime.EQ.'RW') .AND. (t.LT.TurbPartState(4,iPart))) CYCLE
!#endif
    ! Calculate the drag (and gravity) force
#if !USE_FAXEN_CORR
    Pt(1:3,iPart) = ParticlePush(iPart,FieldAtParticle(PRIM,iPart))
#else
    Pt(1:3,iPart) = ParticlePush(iPart,FieldAtParticle(PRIM,iPart),GradAtParticle(1:RHS_GRAD,1:3,iPart))
#endif /* USE_FAXEN_CORR */
#if USE_EXTEND_RHS
#if USE_BASSETFORCE
    bIter(iPart) = MIN(bIter(iPart) + 1,N_Basset + 2)
#endif /* USE_BASSETFORCE */
    ! Calculate other RHS forces and add all forces to compute the particle push
    CALL ParticlePushExtend(iPart,FieldAtParticle(PRIM,iPart)                                     ,&
                                  GradAtParticle (1:RHS_GRAD,1:3,iPart),Pt(1:PP_nVarPartRHS,iPart) &
#if USE_BASSETFORCE || ANALYZE_RHS
                                  ,t,dt,iStage)
#else
                                  )
#endif /* USE_BASSETFORCE || ANALYZE_RHS */
#endif
  END IF
END DO

!#if ANALYZE_RHS
!IF((dtWriteRHS .GT. 0.0) .AND. (tWriteRHS-t .LE. dt*(1.+1.E-4))) tWriteRHS = tWriteRHS + dtWriteRHS
!#endif /* ANALYZE_RHS */

END SUBROUTINE CalcPartRHS


FUNCTION ParticlePush(PartID,FieldAtParticle&
#if USE_FAXEN_CORR
  ,GradAtParticle)
#else
  )
#endif
!===================================================================================================================================
! Push due to Stoke's drag and source terms (gravity)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals
USE MOD_Particle_Vars,     ONLY: Species, PartSpecies, PartGravity
USE MOD_Particle_Vars,     ONLY: PartState
USE MOD_Particle_Vars,     ONLY: TurbPartState
USE MOD_Viscosity
USE MOD_EoS_Vars,          ONLY: kappa
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(PRIM)
#if USE_FAXEN_CORR
REAL,INTENT(IN)     :: GradAtParticle(1:RHS_GRAD,1:3)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: ParticlePush(1:3)           ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Fdm(1:3)
REAL                :: Rep                         ! Particle Reynolds number
REAL                :: Mp                          ! Particle Mach number
!REAL                :: velosqp                    ! v^2 particle
!REAL                :: velosqf                    ! v^2 fluid
REAL                :: udiff(3),urel
REAL                :: f                           ! Drag factor
REAL                :: staup                       ! Inverse of the particle relaxation time
REAL                :: mu                          ! viscosity
!===================================================================================================================================

! Calculate the dyn. viscosity
#if PARABOLIC
mu = VISCOSITY_PRIM(FieldAtParticle)
#endif

SELECT CASE(Species(PartSpecies(PartID))%RHSMethod)

CASE(RHS_NONE)
!===================================================================================================================================
! Debug RHS method for purely inertial particle movement
!===================================================================================================================================
Fdm = 0.

CASE(RHS_TRACER)
!===================================================================================================================================
! Passive tracer moving with fluid velocity
!===================================================================================================================================
IF (ALLOCATED(TurbPartState)) THEN
  Fdm = FieldAtParticle(VELV) + TurbPartState(1:3,PartID)
ELSE
  Fdm = FieldAtParticle(VELV)
END IF

CASE(RHS_TCONVERGENCE)
!===================================================================================================================================
! Special case, drag force only active in x-direction, fixed differential. Gravity in y-direction. Used for convergence tests
!===================================================================================================================================
! Particle relaxation time
Fdm(1) = (FieldAtParticle(VEL1) - PartState(PART_VEL1,PartID)) * 0.5 * 1./Species(PartSpecies(PartID))%StokesIC
Fdm(2) = PartGravity(2)
Fdm(3) = 0.

CASE(RHS_HPCONVERGENCE)
!===================================================================================================================================
! Special case, part push depends only on fluid velocity, which is only active in x-direction, fixed differential.
! Used for h-/p-convergence tests
!===================================================================================================================================
Fdm(1) = FieldAtParticle(VEL1); Fdm(2:3) = 0.

CASE(RHS_SGS1)
!===================================================================================================================================
! Calculation according to Stokes; SGS model provides full fluid velocity seen, see Jean-Pierre Minier & Eric Peirano [2001]
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL Abort(__STAMP__,'Tracking of inertial particles requires mu to be set!')

! In Minier case, TurbPartState it's the FULL  fluid velocity
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)     - PartState(PART_VELV,PartID)
END IF

urel = VECNORM(udiff(1:3))
Rep  = urel*PartState(PART_DIAM,PartID)*FieldAtParticle(DENS)/mu
Mp   = urel/SPEEDOFSOUND_H(FieldAtParticle(PRES),(1./FieldAtParticle(DENS)))

! Empirical relation of nonlinear drag from Clift et al. (1978)
#if USE_SPHERICITY
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, PartState(PART_SPHE,PartID), Mp)
#else
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, Species(PartSpecies(PartID))%SphericityIC, Mp)
#endif

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./PartState(PART_DIAM,PartID)**2

Fdm      = udiff * staup * f

! Add gravity if required
IF(ANY(PartGravity.NE.0)) Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)

CASE(RHS_SGS2)
!===================================================================================================================================
! Calculation according to Stokes; SGS force is added to RHS
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL Abort(__STAMP__,'Tracking of inertial particles requires mu to be set!')

udiff(1:3) = FieldAtParticle(VELV) - PartState(PART_VELV,PartID)
#if USE_FAXEN_CORR
udiff = udiff + (PartState(PART_DIAM,PartID)**2)/6 * GradAtParticle(RHS_LAPLACEVEL,:)
#endif /* USE_FAXEN_CORR */

urel = VECNORM(udiff(1:3))
Rep  = urel*PartState(PART_DIAM,PartID)*FieldAtParticle(DENS)/mu
Mp   = urel/SPEEDOFSOUND_H(FieldAtParticle(PRES),(1./FieldAtParticle(DENS)))

#if USE_SPHERICITY
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, PartState(PART_SPHE,PartID), Mp)
#else
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, Species(PartSpecies(PartID))%SphericityIC, Mp)
#endif

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./PartState(PART_DIAM,PartID)**2

Fdm      = udiff * staup * f

! Add gravity and bouyancy if required
IF(ANY(PartGravity.NE.0)) Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)
IF(ALLOCATED(TurbPartState)) Fdm = Fdm + TurbPartState(:,PartID)

CASE(RHS_INERTIA)
!===================================================================================================================================
! Calculation according to Stokes
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL Abort(__STAMP__,'Tracking of inertial particles requires mu to be set!')

! Assume spherical particles for now
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = FieldAtParticle(VELV) + TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)                             - PartState(PART_VELV,PartID)
END IF
#if USE_FAXEN_CORR
udiff = udiff + (PartState(PART_DIAM,PartID)**2)/6 * GradAtParticle(RHS_LAPLACEVEL,:)
#endif /* USE_FAXEN_CORR */

urel = VECNORM(udiff(1:3))
Rep  = urel*PartState(PART_DIAM,PartID)*FieldAtParticle(DENS)/mu
Mp   = urel/SPEEDOFSOUND_H(FieldAtParticle(PRES),(1./FieldAtParticle(DENS)))

#if USE_SPHERICITY
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, PartState(PART_SPHE,PartID), Mp)
#else
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, Species(PartSpecies(PartID))%SphericityIC, Mp)
#endif

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./PartState(PART_DIAM,PartID)**2

Fdm      = udiff * staup * f

! Add gravity and bouyancy if required
IF(ANY(PartGravity.NE.0)) Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)

CASE DEFAULT
  CALL Abort(__STAMP__, 'No valid RHS method given. Species',IntInfo=PartSpecies(PartID))

END SELECT

ParticlePush(1:3) = Fdm

END FUNCTION ParticlePush


#if PP_nVarPartRHS == 6
FUNCTION ParticlePushRot(PartID,FieldAtParticle,Omega,Rew)
!===================================================================================================================================
! Push due to Stoke's drag and source terms (gravity)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals
USE MOD_Particle_Vars,     ONLY : Species, PartSpecies
USE MOD_Viscosity
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(PRIM)
REAL,INTENT(IN)     :: Omega(1:3),Rew
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: ParticlePushRot(1:3)        ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Cw        ! angular relative fluid-particle velocity
!===================================================================================================================================
Cw = 6.45/SQRT(Rew) + 32.1/SQRT(Rew)
ParticlePushRot(1:3) = FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC * Cw * 15/16 * 1./PI * Omega * VECNORM(Omega)
END FUNCTION ParticlePushRot
#endif


#if USE_EXTEND_RHS
SUBROUTINE ParticlePushExtend(PartID,FieldAtParticle,GradAtParticle,Pt_in&
#if USE_BASSETFORCE || ANALYZE_RHS
  ,t,dt,iStage)
#else
  )
#endif
!===================================================================================================================================
! Push due to additional forces, e.g. Basset force, lift force, added mass effect, viscous and pressure forces
! REMARK:
! The lift and rotational forces (Saffman and Magnus force) are generally small
! The added mass effect and the viscous and pressure forces are similarly calculated and decrease the particle velocity.
! The Basset force has the highest influence off all forces. It decreases the particle velocity as well as it adresses the temporal
! delay due to visous effects.
! All three forces are relevant/significant for unsteady flow independent of rho/rho_p!
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals
USE MOD_Mathtools,              ONLY: CROSS
USE MOD_Particle_Vars,          ONLY: Species,PartSpecies
USE MOD_Particle_Vars,          ONLY: PartState,TurbPartState
#if ANALYZE_RHS
!USE MOD_Particle_Vars,          ONLY: tWriteRHS,FileName_RHS,dtWriteRHS,Pt_ext
USE MOD_Particle_Vars,          ONLY: Pt_ext
!USE MOD_Output,                 ONLY: OutputToFile
#endif /* ANALYZE_RHS */
USE MOD_PreProc,                ONLY: PP_pi
USE MOD_Viscosity
#if USE_BASSETFORCE
USE MOD_Equation_Vars,          ONLY: s43,s23
USE MOD_Particle_Vars,          ONLY: durdt,N_Basset,bIter,FbCoeff,FbCoeffa,Fbdt,FbCoefft!,Fbi,FbCoeffm
USE MOD_TimeDisc_Vars,          ONLY: nRKStages, RKC
#endif /* USE_BASSETFORCE */
#if USE_UNDISTFLOW || USE_VIRTUALMASS
USE MOD_Particle_TimeDisc_Vars, ONLY: useManualTimeStep
#endif /* USE_UNDISTFLOW || USE_VIRTUALMASS */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: PartID
REAL,INTENT(IN)             :: FieldAtParticle(PRIM)
REAL,INTENT(IN)             :: GradAtParticle(1:RHS_GRAD,1:3)
REAL,INTENT(INOUT)          :: Pt_in(1:PP_nVarPartRHS)
#if USE_BASSETFORCE || ANALYZE_RHS
REAL,INTENT(IN)             :: t
REAL,INTENT(IN)             :: dt
INTEGER,INTENT(IN),OPTIONAL :: iStage
#endif /* USE_BASSETFORCE */
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: Pt(1:PP_nVarPartRHS),Fdi(1:3,10),Fre(1:3,10)
REAL                     :: udiff(3)                    ! velocity difference
REAL                     :: mu                          ! viscosity
REAL                     :: globalfactor                ! prefactor of LHS divided by the particle mass
REAL                     :: Flm(1:3)                    ! Saffman force divided by the particle mass
REAL                     :: Fmm(1:3)                    ! Magnus force divided by the particle mass
REAL                     :: Fum(1:3)                    ! undisturbed flow force divided by the particle mass
REAL                     :: Fvm(1:3)                    ! virtual mass force divided by the particle mass
REAL                     :: Fbm(1:3)                    ! Basset force divided by the particle mass
#if (USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE || PP_nVarPartRHS == 6)
REAL                     :: prefactor                   ! factor divided by the particle mass
#endif
#if (USE_UNDISTFLOW || USE_VIRTUALMASS)
REAL                     :: DuDt(1:3)                   ! viscous and pressure forces divided by the particle mass
#endif
#if USE_BASSETFORCE
REAL,PARAMETER           :: s32=3./2.
INTEGER                  :: k,kIndex,nIndex
REAL                     :: dufdt(1:3)                  ! partial derivative of the fluid velocity
REAL                     :: dtk(3),dtn(2)
REAL                     :: tmp(1:3)
! Convergence1
! REAL                     :: Sb,Cb
#endif /* USE_BASSETFORCE */
#if PP_nVarPartRHS == 6
REAL                     :: Rep                         ! particle Reynolds number
REAL                     :: Omega(3),Rew                ! relative fluid-particle Angular velocity, rotational Reynolds number
REAL                     :: rotu(3),rotudiff(3)         ! curl product of the velocity and the velocity difference
REAL                     :: dotp,beta                   ! dot_product, beta=dp*|\omega|/(2*udiff)
#endif /* PP_nVarPartRHS == 6 */
!===================================================================================================================================

SELECT CASE(Species(PartSpecies(PartID))%RHSMethod)
  CASE(RHS_NONE,RHS_TCONVERGENCE,RHS_HPCONVERGENCE,RHS_TRACER)
    RETURN
END SELECT

#if PARABOLIC
! Calculate the dyn. viscosity
mu = VISCOSITY_PRIM(FieldAtParticle)
#endif

! Calcuate velocity difference
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = FieldAtParticle(VELV) + TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)                             - PartState(PART_VELV,PartID)
END IF

! Nullify arrays
Pt(:) = 0.
Flm = 0.; Fbm = 0.; Fvm=0.; Fum=0.; Fmm=0.
#if PP_nVarPartRHS == 6
Rew = 0.; Rep = 0.
#endif
! factor before left hand side to add all dv_p/dt terms of the RHS
globalfactor = 1.

!===================================================================================================================================
! Calculate the Saffman lift force according to:
! Saffman, P.G.: The lift on a small sphere in a slow shear flow. Journal of Fluid Mechanics,
! pp. 385–400, 1965. 10.1017/S0022112065000824.
!===================================================================================================================================
#if PP_nVarPartRHS == 6
! Calculate the RHS of the rotation
IF (Species(PartSpecies(PartID))%CalcMagnusForce .OR. Species(PartSpecies(PartID))%CalcSaffmanForce) THEN
  ! Calculate the rotation: \nabla x u
  rotu = (/GradAtParticle(RHS_GRADVEL3,2)-GradAtParticle(RHS_GRADVEL2,3),&
           GradAtParticle(RHS_GRADVEL1,3)-GradAtParticle(RHS_GRADVEL3,1),&
           GradAtParticle(RHS_GRADVEL2,1)-GradAtParticle(RHS_GRADVEL1,2)/)
  ! Calculate the Re number
  Rep = VECNORM(udiff(1:3))*PartState(PART_DIAM,PartID)*FieldAtParticle(DENS)/mu
END IF

IF (Species(PartSpecies(PartID))%CalcSaffmanForce .AND. .NOT. ALMOSTZERO(Rep)) THEN
  ! Calculate the factor
  prefactor = 9.69/(Species(PartSpecies(PartID))%DensityIC*PartState(PART_DIAM,PartID)*PP_PI)
  beta = PartState(PART_DIAM,PartID) * VECNORM(rotu) * 0.5 / VECNORM(udiff)
  IF (Rep .LE. 40) THEN
    prefactor = prefactor * (1-0.3314*SQRT(beta)*EXP(-0.1*Rep)+0.3314*SQRT(beta))
  ELSE
    prefactor = prefactor * 0.0524*SQRT(beta*Rep)
  END IF
  dotp   = NORM2(rotu(:))
  IF (dotp .GT. 1.E-3) Flm(:) = SQRT(FieldAtParticle(DENS)*mu * 1./dotp) * CROSS(rotu, udiff) * prefactor
END IF

!===================================================================================================================================
! Calculate the Magnus force according to:
! Rubinow, S.I., Keller, J.B.: The transverse force on a spinning sphere moving in a viscous
! fluid. Journal of Fluid Mechanics, pp. 447–459, 1961. 10.1017/S0022112061000640.
!===================================================================================================================================
IF (Species(PartSpecies(PartID))%CalcMagnusForce) THEN
  ! Relative fluid-particle Angular velocity
  Omega = 0.5 * rotu - PartState(PART_AMOMV,PartID)
  ! Prefactor according to Oesterle and Bui Dinh
  IF (ANY(Omega.NE.0)) THEN
    Rew = FieldAtParticle(DENS) * VECNORM(Omega) * PartState(PART_DIAM,PartID)**2 / (4*mu)
    Pt(4:6) = ParticlePushRot(PartID,FieldAtParticle(PRIM),Omega,Rew)
    ! Calculate the rotation: (\nabla x u_p) x udiff
    rotudiff = CROSS(Omega, udiff) * VECNORM(udiff) / VECNORM(Omega)
    IF (Rep.GT.0.) THEN
      prefactor = 0.45 + (4*Rew/Rep-0.45)*EXP(-0.09896*Rew**0.4*Rep**(-0.3))
      Fmm = 3./4 * prefactor * FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC * rotudiff
    END IF
  END IF
END IF
#endif /* PP_nVarPartRHS */

#if USE_UNDISTFLOW || USE_VIRTUALMASS
IF (Species(PartSpecies(PartID))%CalcUndisturbedFlow.OR.Species(PartSpecies(PartID))%CalcVirtualMass) THEN
  ! Material derivative D(u_i)/Dt = \partial (u_i)/\partial t + u_j \partial (u_i)/\partial (x_j) (inkomp.)
  IF (useManualTimeStep) THEN
    DuDt(1)  = DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL1,:))
    DuDt(2)  = DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL2,:))
    DuDt(3)  = DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL3,:))
  ELSE
    DuDt(1)  = GradAtParticle(RHS_dVELdt,1) + DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL1,:))
    DuDt(2)  = GradAtParticle(RHS_dVELdt,2) + DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL2,:))
    DuDt(3)  = GradAtParticle(RHS_dVELdt,3) + DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL3,:))
  END IF
END IF
#endif /* USE_UNDISTFLOW || USE_VIRTUALMASS */

!===================================================================================================================================
! Calculate the viscous and pressure forces (non-conservative compressible form of the NSE):
! 1/\rho_p Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
!===================================================================================================================================
#if USE_UNDISTFLOW
IF (Species(PartSpecies(PartID))%CalcUndisturbedFlow) THEN
  prefactor = FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC

  Fum(1:3) = prefactor * DuDt(1:3)
END IF
#endif /* USE_UNDISTFLOW */

!===================================================================================================================================
! Calculate the added mass force according to:
! Auton, T.R., Hunt, J.C., Prud’Homme, M.: The force exerted on a body in invis-
! cid unsteady non-uniform rotational flow. Journal of Fluid Mechanics, pp. 241–257, 1988.
! DOI: 10.1017/S0022112088003246.
! Non-conservative compressible form of the NSE:
! 1/\rho_p Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
!===================================================================================================================================
#if USE_VIRTUALMASS
IF (Species(PartSpecies(PartID))%CalcVirtualMass) THEN
  prefactor = 0.5*FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC

  Fvm(1:3) = prefactor * DuDt(1:3)

  ! Add to global scaling factor as 0.5*\rho/\rho_p*dv_p/dt is on RHS
  globalfactor = globalfactor + prefactor
END IF
#endif /* USE_VIRTUALMASS */

!===================================================================================================================================
! Calculate the Basset force / history terms following:
! van Hinsberg, M.A., ten Thije Boonkkamp, J.H., Clercx, H.J.: An efficient, second or-
! der method for the approximation of the Basset history force. Journal of Computational Physics,
! pp. 1465–1478, 2011. 10.1016/j.jcp.2010.11.014.
!===================================================================================================================================
#if USE_BASSETFORCE
IF (Species(PartSpecies(PartID))%CalcBassetForce) THEN
  ! Index for previous data
  nIndex = MIN(N_Basset+1, bIter(PartID))
  kIndex = INT((nIndex)*3)

  ! copy previous data
  IF (bIter(PartID) .GT. N_Basset+1) THEN
    tmp(1:3)                 = durdt(1:3,PartID)
    durdt(1:kIndex-3,PartID) = durdt(4:kIndex,PartID)
    Fbdt(1:nIndex,PartID)    = Fbdt(2:nIndex+1,PartID)
  END IF

  IF (PRESENT(iStage)) THEN
    IF (iStage.EQ.1) THEN
      Fbdt(nIndex+1,PartID) = t+RKC(2)*dt
    ELSE
      IF (iStage.NE.nRKStages) THEN
        Fbdt(nIndex+1,PartID) = t+(RKC(iStage+1)-RKC(iStage))*dt
      ELSE
        Fbdt(nIndex+1,PartID) = t+(1.-RKC(nRKStages))*dt
      END IF
    END IF
  ELSE
   Fbdt(nIndex+1,PartID) = t+dt
  END IF

  ! Scaling factor
  prefactor = 9./(PartState(PART_DIAM,PartID)*Species(PartSpecies(PartID))%DensityIC)&
            * SQRT(FieldAtParticle(DENS)*mu/(PP_pi))

  ! d(u_i)/dt = \partial (u_i)/\partial t + v_j \partial (u_i)/\partial (x_j) (inkomp.)
  dufdt(1) = GradAtParticle(RHS_dVELdt,1) + DOT_PRODUCT(PartState(PART_VELV,PartID), GradAtParticle(RHS_GRADVEL1,:))
  dufdt(2) = GradAtParticle(RHS_dVELdt,2) + DOT_PRODUCT(PartState(PART_VELV,PartID), GradAtParticle(RHS_GRADVEL2,:))
  dufdt(3) = GradAtParticle(RHS_dVELdt,3) + DOT_PRODUCT(PartState(PART_VELV,PartID), GradAtParticle(RHS_GRADVEL3,:))

  ! current derivative
  durdt(kIndex-2:kIndex,PartID) = dufdt(:)
  ! NOTE: convergence1
  ! durdt(kIndex-2:kIndex,PartID) = COS(t)
  ! NOTE: convergence2
  ! durdt(kIndex-2:kIndex,PartID) = t**2

  IF (PRESENT(iStage)) THEN
    dtk(1) = Fbdt(nIndex+1,PartID) - Fbdt(1,PartID)
    dtk(2) = Fbdt(nIndex+1,PartID) - Fbdt(2,PartID)
    dtn(1) = Fbdt(nIndex+1,PartID) - Fbdt(nIndex,PartID)
    dtn(2) = Fbdt(2,PartID) - Fbdt(1,PartID)
    Fbm = s43 * durdt(kIndex-2:kIndex,PartID) * SQRT(dtn(1)) + durdt(kIndex-2-(nIndex-1)*3:kIndex-(nIndex-1)*3,PartID) * &
          (2*(SQRT(dtk(1))-SQRT(dtk(2))) + s23/dtn(2)*(dtk(1)**1.5-dtk(2)**1.5) - &
          2*dtk(1)/dtn(2)*(SQRT(dtk(1))-SQRT(dtk(2))))
    DO k=1,nIndex-1
      dtk(1) = Fbdt(nIndex+1,PartID)   - Fbdt(nIndex-k+1,PartID)
      dtk(2) = Fbdt(nIndex+1,PartID)   - Fbdt(nIndex-k+2,PartID)
      dtk(3) = Fbdt(nIndex+1,PartID)   - Fbdt(nIndex-k,PartID)
      dtn(1) = Fbdt(nIndex-k+2,PartID) - Fbdt(nIndex-k+1,PartID)
      dtn(2) = Fbdt(nIndex-k+1,PartID) - Fbdt(nIndex-k,PartID)
      Fbm = Fbm + durdt(kIndex-2-k*3:kIndex-k*3,PartID) * (s23*(dtk(1)**1.5-dtk(2)**1.5)/dtn(1) + &
                  2*(SQRT(dtk(1))-SQRT(dtk(2))) - 2*dtk(1)/dtn(1)*(SQRT(dtk(1))-SQRT(dtk(2))) + &
                  2*dtk(3)/dtn(2)*(SQRT(dtk(3))-SQRT(dtk(1))) - s23/dtn(2)*(dtk(3)**1.5-dtk(1)**1.5))
    END DO
  ELSE
    dtn(1) = Fbdt(nIndex+1,PartID) - Fbdt(nIndex,PartID)
    dtn(2) = Fbdt(2,PartID) - Fbdt(1,PartID)
    Fbm = s43 * durdt(kIndex-2:kIndex,PartID) * SQRT(dtn(1)) + durdt(kIndex-2-(nIndex-1)*3:kIndex-(nIndex-1)*3,PartID) * &
    FbCoeff(N_Basset+nIndex-1) * SQRT(dtn(2))
    DO k=1,nIndex-1
      dtn(1) = Fbdt(nIndex-k+2,PartID) - Fbdt(nIndex-k+1,PartID)
      Fbm = Fbm + durdt(kIndex-2-k*3:kIndex-k*3,PartID) * FbCoeff(k) * SQRT(dtn(1))
    END DO
  END IF

  ! IF (biter(PartID) .GT. N_Basset+1) THEN
  !   DO k=1,FbCoeffm
  !     ! F_i-di(t) = 2 c_B sqrt(e t) exp(-t_win/(2 FbCoefft)) (g_n...)
  !     Fdi(1:3,k) = 2*SQRT(EXP(1.)*FbCoefft(k))*EXP(-SUM(Fbdt(1:nIndex,PartID))/(2*FbCoefft(k)))*(durdt(kIndex-2-(nIndex-1)*3:kIndex-(nIndex-1)*3,PartID)*&
  !       (1-(PHI((-Fbdt(nIndex,PartID)/(2*FbCoefft(k))))))+tmp*EXP(-Fbdt(nIndex,PartID)/(2*FbCoefft(k)))*((PHI((Fbdt(nIndex,PartID)/(2*FbCoefft(k)))))-1.))
  !     ! F_i-re(t) = exp(-dt/(2 FbCoefft)) F_i (t-dt)
  !     Fre(1:3,k) = EXP(-Fbdt(nIndex,PartID)/(2*FbCoefft(k)))*Fbi(1:3,k,PartID)
  !   END DO

  !   DO k=1,FbCoeffm
  !     Fbi(1:3,k,PartID) = FbCoeffa(k) * (prefactor*Fdi(1:3,k) + Fre(1:3,k))
  !     Fbm = Fbm + Fbi(1:3,k,PartID)
  !   END DO
  ! END IF
  Fbm(1:3) = prefactor * Fbm(1:3)

  ! Add to global scaling factor as s43*\rho*prefactor*dv_p/dt is on RHS
  globalfactor     = globalfactor + s43 * prefactor * SQRT(Fbdt(nIndex+1,PartID)-Fbdt(nIndex,PartID))

  ! Correction term if initial particle velocity if different from the surrounding fluid velocity
  IF (biter(PartID) .EQ. 1 .AND. t.GT. 0) THEN
    Fbm(1:3) = Fbm(1:3) + prefactor * (FieldAtParticle(VELV) - PartState(PART_VELV,PartID)) / SQRT(t)
  END IF

  ! NOTE: convergence1
  ! Cb = 0.5 + (1+0.926*SQRT(2/PI*t))/(2+1.792*SQRT(2/PI*t)+3.104*2/PI*t) * SIN(t) -&
  !            1./(2+4.142*SQRT(2/PI*t)+3.492*2/PI*t+6.67*SQRT(2/PI*t)**3) * COS(t)
  ! Sb = 0.5 - (1+0.926*SQRT(2/PI*t))/(2+1.792*SQRT(2/PI*t)+3.104*2/PI*t) * COS(t) -&
  !            1./(2+4.142*SQRT(2/PI*t)+3.492*2/PI*t+6.67*SQRT(2/PI*t)**3) * SIN(t)
  ! tmp(3) = prefactor*SQRT(2*PI)*(Cb*COS(t)+Sb*SIN(t))
  ! NOTE: convergence2
  ! tmp(3) = prefactor*2/15*t**0.5*(8*t**2)

  Pt(1:3) = (Flm + Fmm + Fum + Fvm + Fbm + Pt_in(1:3)) * 1./globalfactor

  ! \rho d(udiff)/dt = \rho d(u)/dt - \rho (dv_p/dt)
  durdt(kIndex-2:kIndex,PartID) = durdt(kIndex-2:kIndex,PartID) - Pt(1:3)
ELSE
#endif /* USE_BASSETFORCE */
  Pt(1:3) = (Flm + Fmm + Fum + Fvm + Fbm + Pt_in(1:3)) * 1./globalfactor
#if USE_BASSETFORCE
END IF
#endif /* USE_BASSETFORCE */

! Output RHS to file
#if ANALYZE_RHS
!IF(dtWriteRHS.GT.0.0)THEN
!  IF(tWriteRHS-t.LE.dt*(1.+1.E-4))THEN
Pt_ext(1:3  ,PartID) = Pt_in(1:3)
#if USE_VIRTUALMASS
IF (Species(PartSpecies(PartID))%CalcVirtualMass) &
  Pt_ext(4:6,PartID) = Fvm-Pt(1:3)*0.5*FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC
#endif
Pt_ext(7:9  ,PartID) = Fum(1:3)
Pt_ext(10:12,PartID) = Flm(1:3)
Pt_ext(13:15,PartID) = Fmm(1:3)
#if USE_BASSETFORCE
IF (Species(PartSpecies(PartID))%CalcBassetForce) THEN
  prefactor = 9./(PartState(PART_DIAM,PartID)*Species(PartSpecies(PartID))%DensityIC)&
            * SQRT(FieldAtParticle(DENS)*mu/(PP_pi))
  Pt_ext(16:18,PartID) = Fbm-Pt(1:3)*s43*prefactor*SQRT(Fbdt(nIndex+1,PartID)-Fbdt(nIndex,PartID))
END IF
#endif /* USE_BASSETFORCE */
!    CALL OutputToFile(FileName_RHS,(/t/),(/23,1/),&
!#if USE_BASSETFORCE
!      (/REAL(PartSpecies(PartID)),Pt(1:3),Pt_in(1:3),Flm(1:3),Fmm(1:3),Fum(1:3),Fvm(1:3),Fbm(1:3),REAL(bIter(PartID),8)/))
!#else
!      (/REAL(PartSpecies(PartID)),Pt(1:3),Pt_in(1:3),Flm(1:3),Fmm(1:3),Fum(1:3),Fvm(1:3),Fbm(1:3),REAL(1.,8)/))
!#endif /* USE_BASSETFORCE */
!  END IF
!END IF
#endif

Pt_in(:) = Pt

END SUBROUTINE ParticlePushExtend
#endif /* USE_EXTEND_RHS */

#if USE_EXTEND_RHS || USE_FAXEN_CORR
SUBROUTINE extRHS(UPrim,Ut,U_RHS)
!===================================================================================================================================
! Compute tau
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
#if USE_FAXEN_CORR
USE MOD_Lifting_BR1_gen,    ONLY: Lifting_BR1_gen
USE MOD_Lifting_Vars,       ONLY: gradUx,gradUy,gradUz
USE MOD_Lifting_Vars,       ONLY: gradUx_master,gradUx_slave
USE MOD_Lifting_Vars,       ONLY: gradUy_master,gradUy_slave
USE MOD_Lifting_Vars,       ONLY: gradUz_master,gradUz_slave
USE MOD_Particle_Vars,      ONLY: gradUx2,gradUy2,gradUz2
USE MOD_Particle_Vars,      ONLY: gradUx_master_loc,gradUx_slave_loc
USE MOD_Particle_Vars,      ONLY: gradUy_master_loc,gradUy_slave_loc
USE MOD_Particle_Vars,      ONLY: gradUz_master_loc,gradUz_slave_loc
#endif
!USE MOD_EoS,                ONLY: ConsToPrim
!USE MOD_Equation_Vars,      ONLY: s13
!USE MOD_Particle_Vars,      ONLY: U_local,gradp_local
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: UPrim(PRIM,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(IN)             :: Ut(   CONS,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)            :: U_RHS(1:RHS_NVARS,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE
INTEGER                     :: i,j,k,iElem
#endif
!===================================================================================================================================

U_RHS             = 0.
#if USE_FAXEN_CORR
gradUx2           = 0.
gradUy2           = 0.
gradUz2           = 0.
gradUx_master_loc = gradUx_master(LIFT_VELV,:,:,:)
gradUx_slave_loc  = gradUx_slave( LIFT_VELV,:,:,:)
gradUy_master_loc = gradUy_master(LIFT_VELV,:,:,:)
gradUy_slave_loc  = gradUy_slave( LIFT_VELV,:,:,:)
gradUz_master_loc = gradUz_master(LIFT_VELV,:,:,:)
gradUz_slave_loc  = gradUz_slave( LIFT_VELV,:,:,:)

! Calculate the second gradient of the velocity and output \nabla \cdot tau
CALL Lifting_BR1_gen(3,3,gradUx(LIFT_VELV,:,:,:,:),gradUx_master_loc,gradUx_slave_loc,&
                    gradUx2(1,:,:,:,:,:),gradUx2(2,:,:,:,:,:),gradUx2(3,:,:,:,:,:))
CALL Lifting_BR1_gen(3,3,gradUy(LIFT_VELV,:,:,:,:),gradUy_master_loc,gradUy_slave_loc,&
                    gradUy2(1,:,:,:,:,:),gradUy2(2,:,:,:,:,:),gradUy2(3,:,:,:,:,:))
CALL Lifting_BR1_gen(3,3,gradUz(LIFT_VELV,:,:,:,:),gradUz_master_loc,gradUz_slave_loc,&
                    gradUz2(1,:,:,:,:,:),gradUz2(2,:,:,:,:,:),gradUz2(3,:,:,:,:,:))

! Compute the Laplacian of the fluid velocity
U_RHS(RHS_LAPLACEVEL1,:,:,:,:) = gradUx2(1,1,:,:,:,:) + gradUy2(2,1,:,:,:,:) + gradUz2(3,1,:,:,:,:)
U_RHS(RHS_LAPLACEVEL2,:,:,:,:) = gradUx2(1,2,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(3,2,:,:,:,:)
U_RHS(RHS_LAPLACEVEL3,:,:,:,:) = gradUx2(1,3,:,:,:,:) + gradUy2(2,3,:,:,:,:) + gradUz2(3,3,:,:,:,:)
#endif /* USE_FAXEN_CORR */

#if USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE
! Time derivative of VELV + GradientVELV or RHS of NS equations (same result) for the calculation of the substantial derivative
!! u_xx + u_yy + u_zz + 1/3 * (u_xx+v_yx+w_zx)
!!divtau(1,:,:,:,:) = gradUx2(1,LIFT_VEL1,:,:,:,:) + gradUy2(2,LIFT_VEL1,:,:,:,:) + gradUz2(3,LIFT_VEL1,:,:,:,:) + &
!!                    s13 * (gradUx2(1,LIFT_VEL1,:,:,:,:) + gradUy2(1,LIFT_VEL2,:,:,:,:) + gradUz2(1,LIFT_VEL3,:,:,:,:))
!U_RHS(RHS_DIVTAU1,:,:,:,:) = gradUx2(1,1,:,:,:,:) + gradUy2(2,1,:,:,:,:) + gradUz2(3,1,:,:,:,:) + &
!                      s13 * (gradUx2(1,1,:,:,:,:) + gradUy2(1,2,:,:,:,:) + gradUz2(1,3,:,:,:,:))
!! v_xx + v_yy + v_zz + 1/3 * (u_xy+v_yy+w_zy)
!!divtau(2,:,:,:,:) = gradUx2(1,LIFT_VEL2,:,:,:,:) + gradUy2(2,LIFT_VEL2,:,:,:,:) + gradUz2(3,LIFT_VEL2,:,:,:,:) + &
!!                    s13 * (gradUx2(2,LIFT_VEL1,:,:,:,:) + gradUy2(2,LIFT_VEL2,:,:,:,:) + gradUz2(2,LIFT_VEL3,:,:,:,:))
!U_RHS(RHS_DIVTAU2,:,:,:,:) = gradUx2(1,2,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(3,2,:,:,:,:) + &
!                      s13 * (gradUx2(2,1,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(2,3,:,:,:,:))
!! w_xx + w_yy + w_zz + 1/3 * (u_xy+v_yy+w_zy)
!!divtau(3,:,:,:,:) = gradUx2(1,LIFT_VEL3,:,:,:,:) + gradUy2(2,LIFT_VEL3,:,:,:,:) + gradUz2(3,LIFT_VEL3,:,:,:,:) + &
!!                    s13 * (gradUx2(3,LIFT_VEL1,:,:,:,:) + gradUy2(3,LIFT_VEL2,:,:,:,:) + gradUz2(3,LIFT_VEL3,:,:,:,:))
!U_RHS(RHS_DIVTAU3,:,:,:,:) = gradUx2(1,3,:,:,:,:) + gradUy2(2,3,:,:,:,:) + gradUz2(3,3,:,:,:,:) + &
!                      s13 * (gradUx2(3,1,:,:,:,:) + gradUy2(3,2,:,:,:,:) + gradUz2(3,3,:,:,:,:))
!
!! Calculate pressure gradient
!DO iElem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
!  CALL ConsToPrim(U_local(:,i,j,k,iElem),U(:,i,j,k,iElem))
!END DO; END DO; END DO; END DO
!CALL Lifting_BR1_gen(1,1,U_local(PRES:PRES,:,:,:,:),gradp_local(:,1,:,:,:,:),gradp_local(:,2,:,:,:,:),gradp_local(:,3,:,:,:,:))
!U_RHS(RHS_GRADP1:RHS_GRADP3,:,:,:,:) = gradp_local(1,:,:,:,:,:)

DO iElem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  U_RHS(RHS_dVELVdt,i,j,k,iElem) = (Ut(MOMV,i,j,k,iElem) - Ut(DENS,i,j,k,iElem)*UPrim(VELV,i,j,k,iElem))*1./UPrim(DENS,i,j,k,iElem)
END DO; END DO; END DO; END DO
#endif /* USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE */

END SUBROUTINE extRHS
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */


FUNCTION DF_SchillerAndNaumann(Rep, SphericityIC, Mp) RESULT(f)
!===================================================================================================================================
! Compute the drag factor according to Schiller and Naumann
!===================================================================================================================================
! MODULES
USE MOD_Globals,           ONLY : UNIT_stdOut
USE MOD_Particle_Vars,     ONLY : RepWarn
#if USE_MPI
USE MOD_Globals,           ONLY : MPIRoot
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: Rep, SphericityIC, Mp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                        :: f
!-----------------------------------------------------------------------------------------------------------------------------------
! Warn when outside valid range of Naumann model
IF(Rep.GT.800) THEN
  IF (RepWarn.EQV..FALSE.) THEN
    SWRITE(UNIT_stdOut,*) 'WARNING: Rep',Rep,'> 800, drag coefficient may not be accurate.'
    RepWarn=.TRUE.
  ENDIF
ENDIF
f  = 1. + 0.15*Rep**0.687

! Suppress compiler warning
NO_OP(SphericityIC)
NO_OP(MP)
END FUNCTION DF_SchillerAndNaumann

FUNCTION DF_Putnam(Rep, SphericityIC, Mp) RESULT(f)
!===================================================================================================================================
! Compute the drag factor according to Putnam et al. (1961)
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: Rep, SphericityIC, Mp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                        :: f
!-----------------------------------------------------------------------------------------------------------------------------------
f = 1. + (Rep**2./3.)/6.
! High Re correction according to Putnam et al. (1961)
IF(Rep .GT. 1000) f = 0.0183*Rep

! Suppress compiler warning
NO_OP(SphericityIC)
NO_OP(MP)
END FUNCTION DF_Putnam

FUNCTION DF_Haider(Rep, SphericityIC, Mp) RESULT(f)
!===================================================================================================================================
! Compute the drag factor according to Haider and Levenspiel (1989) valid for Rep<2.6e5
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: Rep, SphericityIC, Mp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                        :: f
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: k1, k2, k3
!-----------------------------------------------------------------------------------------------------------------------------------
k1 = EXP(2.3288-6.4581*SphericityIC+2.4486*SphericityIC**2)
k2 = EXP(4.905-13.8944*SphericityIC+18.4222*SphericityIC**2-10.2599*SphericityIC**3)
k3 = EXP(1.4681+12.2584*SphericityIC-20.7322*SphericityIC**2+15.8855*SphericityIC**3)
f = (1+k1*Rep**(0.0964+0.5565*SphericityIC))+Rep**2*1./24*k2/(Rep+k3)

! Suppress compiler warning
NO_OP(MP)
END FUNCTION DF_Haider

FUNCTION DF_Hoelzer(Rep, SphericityIC, Mp) RESULT(f)
!===================================================================================================================================
! Compute the drag factor according to Hoelzer et al. (2008)
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars,      ONLY: s13,s23
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: Rep, SphericityIC, Mp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                        :: f
!-----------------------------------------------------------------------------------------------------------------------------------
f =  s13 * 1./SQRT(SphericityIC) + s23 * 1./SQRT(SphericityIC)+&
     SQRT(Rep)/8. * 1./(SphericityIC**(3./4.)) +&
     Rep/24. * 0.4210**(0.4*(-LOG(SphericityIC))**0.2) *1./SphericityIC

! Suppress compiler warning
NO_OP(MP)
END FUNCTION DF_Hoelzer

!FUNCTION DF_Loth(Rep, SphericityIC, Mp) RESULT(f)
!!===================================================================================================================================
!! Compute the drag factor according to Loth (2008)
!! > Loth, E., Compressibility and Rarefaction Effects on Drag of a Spherical Particle,AIAA Journal, 2008, 46, 2219-2228
!!===================================================================================================================================
!! MODULES
!USE MOD_Equation_Vars,      ONLY: s13,s23
!!-----------------------------------------------------------------------------------------------------------------------------------
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)             :: Rep, SphericityIC, Mp
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL                        :: f
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL                        :: Hm, Cm, Gm
!!-----------------------------------------------------------------------------------------------------------------------------------
!IF (Mp .LT. 0.89) THEN
!  Gm = 1. - 1.525*Mp**4                                                      ! (eq. 16a)
!ELSE
!  Gm = 0.0002 + 0.0008*TANH(12.77*(Mp-2.02))                                 ! (eq. 16b)
!END IF
!IF (Mp .LE. 1.45) THEN
!  Cm = 5.*s13 + s23*TANH(3*LOG(Mp+0.1))                                      ! (eq. 14a)
!ELSE
!  Cm = 2.044 + 0.2*EXP(-1.8*(LOG(Mp/1.5))**2)                                ! (eq. 14b)
!END IF
!Hm = 1 - 0.258*Cm/(1+514*Gm)                                                 ! (eq. 16c)
!! Valid up to Rep < 3e5
!f = (1. + 0.15*Rep**0.687) * Hm + Rep/24*0.42*Cm/(1+42500*Gm*Rep**(-1.16))   ! (eq. 15, divided by Rep/24)
!NO_OP(SphericityIC)
!END FUNCTION DF_Loth

FUNCTION DF_Loth(Rep, SphericityIC, Mp) RESULT(f)
!===================================================================================================================================
! Compute the drag factor according to Loth (2021)
! > Loth, E., Supersonic and hypersonic drag coefficients for a sphere, AIAA Journal, 2021, 59, 3261-3274
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: Rep, SphericityIC, Mp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                        :: f
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: Hm, Cm, Gm
!-----------------------------------------------------------------------------------------------------------------------------------
IF (Mp .LT. 0.8) THEN
  Gm = 166*Mp**3+3.29*Mp**2-10.9*Mp+20                                       ! (eq. 8b)
ELSE
  Gm = 5+40*Mp**(-3)                                                         ! (eq. 8c)
END IF
IF (Mp .LE. 1.5) THEN
  Cm = 1.65 + 0.65*TANH(4*Mp-3.4)                                            ! (eq. 7a)
ELSE
  Cm = 2.18 - 0.13*TANH(0.9*Mp-2.7)                                          ! (eq. 7b)
END IF
IF (Mp .LE. 1.0) THEN
  Hm = 0.0239*Mp**3+0.212*Mp**2-0.074*Mp+1                                   ! (eq. 8d)
ELSE
  Hm = 0.93+1/(3.5+Mp**5)                                                    ! (eq. 8e)
END IF
! Valid up to Rep < 3e5
f = (1. + 0.15*Rep**0.687) * Hm + Rep/24*0.42*Cm/(1+42500*Rep**(-1.16*Cm)+Gm*Rep**(-0.5))   ! (eq. 8a, divided by Rep/24)
NO_OP(SphericityIC)
END FUNCTION DF_Loth

FUNCTION DF_Ganser(Rep, SphericityIC, Mp) RESULT(f)
!===================================================================================================================================
! Compute the drag factor according to Loth (2008)
! > Ganser, G., A rational approach to drag prediction of spherical and nonspherical particles, Powder Technology, 1993, 77, 143-152
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars,      ONLY: s13,s23
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: Rep, SphericityIC, Mp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                        :: f
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: K1,K2
!-----------------------------------------------------------------------------------------------------------------------------------
K1 = (s13 + s23*SphericityIC**(-0.5))**(-1)                                          ! (Table 7)
K2 = 1.8148*(-LOG(SphericityIC)**(-0.5743))                                          ! (Table 7)
! Valid up to Rep < 3e5
f = 1./K1*(1. + 0.118*(K1*K2*Rep)**0.6567) + Rep/24*0.4305*K2/(1+3305/(Rep*K1*K2))   ! (eq. 18)
NO_OP(SphericityIC)
END FUNCTION DF_Ganser

!==================================================================================================================================
!> Compute source terms for particles and add them to the nearest DOF
!==================================================================================================================================
SUBROUTINE CalcSourcePart(Ut)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars      ,ONLY: wGPVol
USE MOD_Particle_Globals  ,ONLY: VECNORM
USE MOD_Mesh_Vars         ,ONLY: Elem_xGP,sJ,nElems,offsetElem
USE MOD_Particle_Globals  ,ONLY: PP_nElems
USE MOD_Particle_Mesh_Vars,ONLY: GEO,FIBGM_nElems,FIBGM_offsetElem,FIBGM_Element
USE MOD_Particle_Mesh_Vars,ONLY: Elem_xGP_Shared
USE MOD_Particle_Mesh_Tools,ONLY: GetCNElemID
#if FV_ENABLED
USE MOD_FV_Vars           ,ONLY: FV_Elems
USE MOD_ChangeBasisByDim  ,ONLY: ChangeBasisVolume
USE MOD_FV_Vars           ,ONLY: FV_Vdm
#endif
USE MOD_Particle_Vars     ,ONLY: Species,PartSpecies,PartState,Pt,PEM,PDM
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem,iPart,ijk(3),iBGM,jBGM,kBGM,imin,imax,jmin,jmax,kmin,kmax,ElemID,CNElemID
REAL                :: Ut_src(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)
REAL                :: PartSource(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)
REAL                :: Fp(3),Wp,r,Vol
REAL                :: min_distance_glob,min_distance_loc
#if FV_ENABLED
REAL                :: Ut_src2(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================

Ut_src = 0.
DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    IF (Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER        .OR. &
        Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TCONVERGENCE  .OR. &
        Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_HPCONVERGENCE .OR. &
        Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER) CYCLE

    Vol = 0.
    PartSource = 0.
    ! Radius of particle
    r = PartState(PART_DIAM,iPart)*0.5
    ! Calculate particle force
    Fp(1:3) = -Pt(1:3,iPart)*Species(PartSpecies(iPart))%MassIC
    ! Calculate the work
    Wp = DOT_PRODUCT(Fp,PartState(4:6,iPart))

    ! --- get background mesh cell of point
    imin = FLOOR((PartState(1,iPart)-r-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
    imin = MAX(GEO%FIBGMimin,imin)
    imax = CEILING((PartState(1,iPart)+r-GEO%xminglob)/GEO%FIBGMdeltas(1))
    imax = MIN(GEO%FIBGMimax,imax)
    jmin = FLOOR((PartState(2,iPart)-r-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
    jmin = MAX(GEO%FIBGMjmin,jmin)
    jmax = CEILING((PartState(2,iPart)+r-GEO%yminglob)/GEO%FIBGMdeltas(2))
    jmax = MIN(GEO%FIBGMjmax,jmax)
    kmin = FLOOR((PartState(3,iPart)-r-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
    kmin = MAX(GEO%FIBGMkmin,kmin)
    kmax = CEILING((PartState(3,iPart)+r-GEO%zminglob)/GEO%FIBGMdeltas(3))
    kmax = MIN(GEO%FIBGMkmax,kmax)

    DO kBGM=kmin,kmax; DO jBGM=jmin,jmax; DO iBGM=imin,imax
      !--- check all cells associated with this background mesh cell
      DO iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)
        ElemID = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iElem)
        CNElemID = GetCNElemID(ElemID)
        ElemID = ElemID-offsetElem
        IF (ElemID.LT.1 .OR. ElemID.GT.PP_nElems) CYCLE
        DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
          IF (NORM2(Elem_xGP_Shared(1:2,i,j,k,CNElemID)-PartState(PART_POS1:PART_POS2,iPart)) .LE. r) THEN
            ! TODO: sJ_shared
            Vol = Vol + wGPVol(i,j,k)/sJ(i,j,k,ElemID,0)
            PartSource(MOMV,i,j,k,ElemID) = Fp
            PartSource(ENER,i,j,k,ElemID) = Wp
          END IF
        END DO; END DO; END DO
      END DO
    END DO; END DO; END DO
    IF (Vol .GT. EPSILON(0.)) THEN
      Ut_src = Ut_src + PartSource / Vol
    ELSE
      iElem = PEM%Element(iPart)-offsetElem
      min_distance_glob = VECNORM(Elem_xGP(:,0,0,0,iElem)-PartState(1:3,iPart))
      ijk(:) = 0
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        min_distance_loc = VECNORM(Elem_xGP(:,i,j,k,iElem)-PartState(1:3,iPart))
        IF (min_distance_loc .LT. min_distance_glob) THEN; ijk(:) = (/i,j,k/); min_distance_glob = min_distance_loc; END IF
      END DO; END DO; END DO
      ! Add source term
      Ut_src(MOMV,ijk(1),ijk(2),ijk(3),iElem) = Ut_src(MOMV,ijk(1),ijk(2),ijk(3),iElem)+&
        Fp*sJ(ijk(1),ijk(2),ijk(3),iElem,0)/wGPVol(ijk(1),ijk(2),ijk(3))
      Ut_src(ENER,ijk(1),ijk(2),ijk(3),iElem) = Ut_src(ENER,ijk(1),ijk(2),ijk(3),iElem)+&
        Wp*sJ(ijk(1),ijk(2),ijk(3),iElem,0)/wGPVol(ijk(1),ijk(2),ijk(3))
      ! Add source term
    END IF
  END IF
END DO ! iPart

DO iElem = 1, nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) THEN ! FV elem
    CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:,iElem),Ut_src2(:,:,:,:))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
    END DO; END DO; END DO ! i,j,k
  ELSE
#endif
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k,iElem)/sJ(i,j,k,iElem,0)
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
  END IF
#endif
END DO

END SUBROUTINE CalcSourcePart


END MODULE MOD_part_RHS
