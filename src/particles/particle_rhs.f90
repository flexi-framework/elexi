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
END SELECT

END SUBROUTINE InitRHS

SUBROUTINE CalcPartRHS(&
#if USE_BASSETFORCE || ANALYZE_RHS
  dt,iStage)
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
#if ANALYZE_RHS
USE MOD_Particle_Vars,                ONLY: tWriteRHS,dtWriteRHS
USE MOD_TimeDisc_Vars,                ONLY: t
#endif /* ANALYZE_RHS */
! #if USE_RW
! USE MOD_Particle_RandomWalk_Vars,     ONLY: RWTime
#if USE_BASSETFORCE
USE MOD_Particle_Vars,                ONLY: Species,PartSpecies,bIter,N_Basset
#endif /* USE_BASSETFORCE */
! USE_MOD_Timedisc_Vars,                ONLY: t
! #endif
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if USE_BASSETFORCE
REAL,INTENT(IN)                  :: dt
INTEGER,INTENT(IN),OPTIONAL      :: iStage
#endif /* USE_BASSETFORCE */
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE
INTEGER                          :: iPart
#if USE_MPI && USE_BASSETFORCE
INTEGER                          :: MPIRequest_FB
#endif
!===================================================================================================================================

! Drag force
Pt(:,1:PDM%ParticleVecLength)=0.

#if USE_BASSETFORCE
bIter = bIter + 1 !MIN(bIter + 1,N_Basset + 1)
! Communicate bIter to all other processors (same effect as blocking comm...)
#if USE_MPI
MPIRequest_FB = MPI_REQUEST_NULL
CALL MPI_IBCAST(bIter,1,MPI_INTEGER,0,MPI_COMM_FLEXI,MPIRequest_FB,IERROR)
#endif /* USE_MPI */
#endif /* USE_BASSETFORCE */
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
    ! Calculate other RHS forces and add all forces to compute the particle push
    CALL ParticlePushExtend(iPart,FieldAtParticle(PRIM,iPart)                                     ,&
                                  GradAtParticle (1:RHS_GRAD,1:3,iPart),Pt(1:PP_nVarPartRHS,iPart) &
#if USE_BASSETFORCE || ANALYZE_RHS
                                  ,dt,iStage)
#else
                                  )
#endif /* USE_BASSETFORCE || ANALYZE_RHS */
#endif
  END IF
END DO

#if USE_BASSETFORCE && USE_MPI
CALL MPI_WAIT(MPIRequest_FB,MPI_STATUS_IGNORE,IERROR)
#endif

#if ANALYZE_RHS
IF((dtWriteRHS .GT. 0.0) .AND. (tWriteRHS-t .LE. dt*(1.+1.E-4))) tWriteRHS = tWriteRHS + dtWriteRHS
#endif /* ANALYZE_RHS */

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
!REAL                :: velosqp                    ! v^2 particle
!REAL                :: velosqf                    ! v^2 fluid
REAL                :: udiff(3)
REAL                :: f                           ! Drag factor
REAL                :: staup                       ! Inverse of the particle relaxation time
REAL                :: mu                          ! viscosity
!===================================================================================================================================

! Calculate the dyn. viscosity
mu=VISCOSITY_PRIM(FieldAtParticle)

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
  Fdm         = FieldAtParticle(VELV) + TurbPartState(1:3,PartID)
ELSE
  Fdm         = FieldAtParticle(VELV)
END IF

CASE(RHS_CONVERGENCE)
!===================================================================================================================================
! Special case, drag force only active in x-direction, fixed differential. Gravity in y-direction. Used for convergence tests
!===================================================================================================================================
Fdm(1)      = FieldAtParticle(VEL1) - PartState(PART_VEL1,PartID)

! Gravity fixed to -3
Fdm(2)      = -3.
Fdm(3)      = 0.

CASE(RHS_LI)
!===================================================================================================================================
! Calculation according to AMY LI and GOODARZ AHMADI [1993]
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL ABORT(__STAMP__,'Tracking of inertial particles requires mu to be set!')

udiff(1:3) = FieldAtParticle(VELV) - PartState(PART_VELV,PartID)

Rep     = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu

! Empirical relation of nonlinear drag from Clift et al. (1978)
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, Species(PartSpecies(PartID))%SphericityIC, 0.)

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./Species(PartSpecies(PartID))%DiameterIC**2

! Adding sgs term
IF(ALLOCATED(TurbPartState)) THEN
  Fdm = udiff * staup * f + TurbPartState(1:3,PartID)
ELSE
  Fdm = udiff * staup * f
END IF

! Add gravity if required
IF(ANY(PartGravity.NE.0)) Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)

CASE(RHS_MINIER)
!===================================================================================================================================
! Calculation according to Jean-Pierre Minier & Eric Peirano [2001]
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL ABORT(__STAMP__,'Tracking of inertial particles requires mu to be set!')

! In Minier case, TurbPartState it's the FULL  fluid velocity
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)     - PartState(PART_VELV,PartID)
END IF

Rep     = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu

! Empirical relation of nonlinear drag from Clift et al. (1978)
f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, Species(PartSpecies(PartID))%SphericityIC, 0.)

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./Species(PartSpecies(PartID))%DiameterIC**2

Fdm      = udiff * staup * f

! Add gravity if required
IF(ANY(PartGravity.NE.0)) Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)

CASE(RHS_INERTIA)
!===================================================================================================================================
! Calculation according to Maxey and Riley (1983)
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL ABORT(__STAMP__,'Tracking of inertial particles requires mu to be set!')

! Assume spherical particles for now
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = FieldAtParticle(VELV) + TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)                             - PartState(PART_VELV,PartID)
END IF
#if USE_FAXEN_CORR
udiff = udiff + (Species(PartSpecies(PartID))%DiameterIC**2)/6 * GradAtParticle(RHS_LAPLACEVEL,:)
#endif /* USE_FAXEN_CORR */

Rep     = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu

f = Species(PartSpecies(PartID))%DragFactor_pointer%op(Rep, Species(PartSpecies(PartID))%SphericityIC, 0.)

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./Species(PartSpecies(PartID))%DiameterIC**2

Fdm      = udiff * staup * f

! Add gravity and bouyancy if required
IF(ANY(PartGravity.NE.0)) Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)

CASE DEFAULT
  CALL ABORT(__STAMP__, 'No valid RHS method given. Species',IntInfo=PartSpecies(PartID))

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
  ,dt,iStage)
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
USE MOD_Particle_Vars,          ONLY: Species, PartSpecies
USE MOD_Particle_Vars,          ONLY: PartState, TurbPartState
#if ANALYZE_RHS
USE MOD_Particle_Vars,          ONLY: tWriteRHS,FileName_RHS,dtWriteRHS
USE MOD_Output,                 ONLY: OutputToFile
#endif /* ANALYZE_RHS */
USE MOD_PreProc,                ONLY: PP_pi
USE MOD_Viscosity
USE MOD_TimeDisc_Vars,          ONLY: t
#if USE_BASSETFORCE
USE MOD_Equation_Vars,          ONLY: s43
USE MOD_Particle_Vars,          ONLY: durdt,N_Basset,bIter
USE MOD_TimeDisc_Vars,          ONLY: nRKStages, RKC
#endif /* USE_BASSETFORCE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: PartID
REAL,INTENT(IN)             :: FieldAtParticle(PRIM)
REAL,INTENT(IN)             :: GradAtParticle(1:RHS_GRAD,1:3)
REAL,INTENT(INOUT)          :: Pt_in(1:PP_nVarPartRHS)
#if USE_BASSETFORCE
REAL,INTENT(IN)             :: dt
INTEGER,INTENT(IN),OPTIONAL :: iStage
#endif /* USE_BASSETFORCE */
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: Pt(1:PP_nVarPartRHS)
REAL                     :: udiff(3)                    ! velocity difference
REAL                     :: mu                          ! viscosity
REAL                     :: globalfactor                ! prefactor of LHS divided by the particle mass
REAL                     :: prefactor                   ! factor divided by the particle mass
REAL                     :: Flm(1:3)                    ! Saffman force divided by the particle mass
REAL                     :: Fmm(1:3)                    ! Magnus force divided by the particle mass
REAL                     :: Fum(1:3)                    ! undisturbed flow force divided by the particle mass
REAL                     :: Fvm(1:3)                    ! virtual mass force divided by the particle mass
REAL                     :: Fbm(1:3)                    ! Basset force divided by the particle mass
#if (USE_UNDISTFLOW || USE_VIRTUALMASS)
REAL                     :: DuDt(1:3)                   ! viscous and pressure forces divided by the particle mass
#endif
#if USE_BASSETFORCE
REAL,PARAMETER           :: s32=3./2.
REAL                     :: RKdtFrac                    ! Runge-Kutta time step
INTEGER                  :: k,kIndex,nIndex
REAL                     :: dufdt(1:3)                  ! partial derivative of the fluid velocity
#endif /* USE_BASSETFORCE */
#if PP_nVarPartRHS == 6
REAL                     :: Rep                         ! particle Reynolds number
REAL                     :: Omega(3),Rew                ! relative fluid-particle Angular velocity, rotational Reynolds number
REAL                     :: rotu(3),rotudiff(3)         ! curl product of the velocity and the velocity difference
REAL                     :: dotp,beta                   ! dot_product, beta=dp*|\omega|/(2*udiff)
#endif /* PP_nVarPartRHS == 6 */
!===================================================================================================================================

SELECT CASE(Species(PartSpecies(PartID))%RHSMethod)
  CASE(RHS_NONE,RHS_CONVERGENCE,RHS_TRACER)
    RETURN
END SELECT

! Calculate the dyn. viscosity
mu=VISCOSITY_PRIM(FieldAtParticle)

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
! Calculate the rotation: \nabla x u
rotu = (/GradAtParticle(RHS_GRADVEL3,2)-GradAtParticle(RHS_GRADVEL2,3),&
         GradAtParticle(RHS_GRADVEL1,3)-GradAtParticle(RHS_GRADVEL3,1),&
         GradAtParticle(RHS_GRADVEL2,1)-GradAtParticle(RHS_GRADVEL1,2)/)
! Relative fluid-particle Angular velocity
Omega = 0.5 * rotu - PartState(PART_AMOMV,PartID)
! Prefactor according to Oesterle and Bui Dinh
IF (ANY(Omega.NE.0)) THEN
  Rew = FieldAtParticle(DENS) * VECNORM(Omega) * Species(PartSpecies(PartID))%DiameterIC**2 / (4*mu)
  Pt(4:6) = ParticlePushRot(PartID,FieldAtParticle(PRIM),Omega,Rew)
END IF

! Calculate the Re number
Rep = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu
IF (Species(PartSpecies(PartID))%CalcSaffmanForce .AND. Rep.GT.0.) THEN
  ! Calculate the factor
  prefactor = 9.69/(Species(PartSpecies(PartID))%DensityIC*Species(PartSpecies(PartID))%DiameterIC*PP_PI)
  beta = Species(PartSpecies(PartID))%DiameterIC * VECNORM(PartState(PART_AMOMV,PartID)) * 0.5 / VECNORM(udiff)
  IF (Rep .LE. 40) THEN
    prefactor = prefactor * (1-0.3314*SQRT(beta)*EXP(-0.1*Rep)+0.3314*SQRT(beta))
  ELSE
    prefactor = prefactor * 0.0524*SQRT(beta*Rep)
  END IF
  ! Calculate the rotation: (\nabla x u) x udiff
  rotudiff = CROSS(rotu, udiff)
  dotp   = MAX(SQRT(DOT_PRODUCT(rotu(:),rotu(:))),0.001)
  IF (dotp .NE. 0) Flm(:) = SQRT(2*FieldAtParticle(DENS)*mu * 1./dotp)*rotudiff(:) * prefactor
END IF

!===================================================================================================================================
! Calculate the Magnus force according to:
! Rubinow, S.I., Keller, J.B.: The transverse force on a spinning sphere moving in a viscous
! fluid. Journal of Fluid Mechanics, pp. 447–459, 1961. 10.1017/S0022112061000640.
!===================================================================================================================================
IF (Species(PartSpecies(PartID))%CalcMagnusForce .AND. Rep.GT.0. .AND. ANY(Omega.NE.0)) THEN
  ! Calculate the rotation: (\nabla x u_p) x udiff
  rotudiff = CROSS(Omega, udiff) * VECNORM(udiff) / VECNORM(Omega)
  prefactor = 0.45 + (4*Rew/Rep-0.45)*EXP(-0.05684*Rew**0.4*Rep**0.7)
  Fmm = PP_PI/8 * prefactor * Species(PartSpecies(PartID))%DiameterIC**3 * FieldAtParticle(DENS) * rotudiff
END IF
#endif /* PP_nVarPartRHS */

#if USE_UNDISTFLOW || USE_VIRTUALMASS
IF (Species(PartSpecies(PartID))%CalcUndisturbedFlow.OR.Species(PartSpecies(PartID))%CalcVirtualMass) THEN
  ! Material derivative \rho D(u_i)/Dt = - \nabla p + \nabla \cdot \tau
!  DuDt(:)  = - GradAtParticle(RHS_GRADPRES,:) + mu*GradAtParticle(RHS_GRADTAU,:)
  DuDt(1)  = GradAtParticle(RHS_dVELdt,1) + DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL1,:))
  DuDt(2)  = GradAtParticle(RHS_dVELdt,2) + DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL2,:))
  DuDt(3)  = GradAtParticle(RHS_dVELdt,3) + DOT_PRODUCT(FieldAtParticle(VELV), GradAtParticle(RHS_GRADVEL3,:))
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
  ! Time integration in first RK stage (p. 26)
  IF(PRESENT(iStage))THEN
    IF (iStage.EQ.1) THEN
      RKdtFrac      = RKC(2)*dt
    ELSE
      IF (iStage.NE.nRKStages) THEN
        RKdtFrac      = (RKC(iStage+1)-RKC(iStage))*dt
      ELSE
        RKdtFrac      = (1.-RKC(nRKStages))*dt
      END IF
    END IF
  ELSE
    RKdtFrac = dt
  END IF
  ! Correction term if initial particle velocity if different from the surrounding fluid velocity
  IF(bIter .EQ. 1) durdt(1:3,PartID) = (FieldAtParticle(VELV) - PartState(PART_VELV,PartID)) / SQRT(t)
  IF(MOD(bIter,N_Basset) .EQ. 0) durdt(1:3,PartID) = (FieldAtParticle(VELV) - PartState(PART_VELV,PartID)) / SQRT(t-N_Basset*RKdtFrac)

  ! Scaling factor
  prefactor = 9./(Species(PartSpecies(PartID))%DiameterIC*Species(PartSpecies(PartID))%DensityIC)&
            * SQRT(FieldAtParticle(DENS)*mu/(PP_pi))

  ! Index for previous data
  nIndex = MIN(N_Basset, bIter)
  kIndex = INT((nIndex+1)*3)

  ! copy previous data
  IF(bIter .GT. N_Basset) durdt(4:kIndex-3,PartID) = durdt(7:kIndex,PartID)
  ! \rho d(u)/dt = \rho D(u)/Dt - udiff * (\rho \nabla(u))
!  dufdt(1) = DuDt(1) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(RHS_GRADVEL1,:))
!  dufdt(2) = DuDt(2) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(RHS_GRADVEL2,:))
!  dufdt(3) = DuDt(3) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(RHS_GRADVEL3,:))
  dufdt(1) = GradAtParticle(RHS_dVELdt,1) + DOT_PRODUCT(PartState(PART_VELV,PartID), GradAtParticle(RHS_GRADVEL1,:))
  dufdt(2) = GradAtParticle(RHS_dVELdt,2) + DOT_PRODUCT(PartState(PART_VELV,PartID), GradAtParticle(RHS_GRADVEL2,:))
  dufdt(3) = GradAtParticle(RHS_dVELdt,3) + DOT_PRODUCT(PartState(PART_VELV,PartID), GradAtParticle(RHS_GRADVEL3,:))

  ! \rho d(udiff)/dt = \rho d(u)/dt - \rho (dv_p/dt)
  durdt(kIndex-2:kIndex,PartID) = dufdt(:)

  Fbm = s43 * durdt(kIndex-2:kIndex,PartID) + durdt(4:6,PartID) * (N_Basset-s43)/((N_Basset-1)*SQRT(REAL(N_Basset-1))+(N_Basset-s32)*SQRT(REAL(N_Basset)))
  DO k=1,nIndex-1
    Fbm = Fbm + durdt(kIndex-2-k*3:kIndex-k*3,PartID) * ((k+s43)/((k+1)*SQRT(REAL(k+1))+(k+s32)*SQRT(REAL(k)))+(k-s43)/((k-1)*SQRT(REAL(k-1))+(k-s32)*SQRT(REAL(k))))
  END DO

  ! Add to global scaling factor as s43*\rho*prefactor*dv_p/dt is on RHS
  globalfactor     = globalfactor + s43 * prefactor * SQRT(RKdtFrac)

  Fbm(1:3) = prefactor * (Fbm(1:3) * SQRT(RKdtFrac) + durdt(1:3,PartID))

  Pt(1:3) = (Flm + Fmm + Fum + Fvm + Fbm + Pt_in(1:3)) * 1./globalfactor

  ! Correct durdt with particle push
  durdt(kIndex-2:kIndex,PartID) = durdt(kIndex-2:kIndex,PartID) - Pt(1:3)
ELSE
#endif /* USE_BASSETFORCE */
  Pt(1:3) = (Flm + Fmm + Fum + Fvm + Fbm + Pt_in(1:3)) * 1./globalfactor
#if USE_BASSETFORCE
END IF
#endif /* USE_BASSETFORCE */

! Output RHS to file
#if ANALYZE_RHS
IF(dtWriteRHS.GT.0.0)THEN
  IF(tWriteRHS-t.LE.dt*(1.+1.E-4))THEN
    CALL OutputToFile(FileName_RHS,(/t/),(/19,1/),(/REAL(PartSpecies(PartID)),Pt_in(1:3),Flm(1:3),Fmm(1:3),Fum(1:3),Fvm(1:3),Fbm(1:3)/))
  END IF
END IF
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
  U_RHS(RHS_dVELVdt,i,j,k,iElem) = Ut(MOMV,i,j,k,iElem) - Ut(DENS,i,j,k,iElem)*UPrim(VELV,i,j,k,iElem)
END DO; END DO; END DO; END DO
#endif /* USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE */

END SUBROUTINE extRHS
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */


FUNCTION DF_SchillerAndNaumann(Rep, SphericityIC, Mp) RESULT(f)
!===================================================================================================================================
! Compute the drag factor according to Schiller and Naumann
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,     ONLY : RepWarn
USE MOD_Globals,           ONLY : MPIRoot, UNIT_stdOut
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
END FUNCTION DF_Hoelzer

!==================================================================================================================================
!> Compute source terms for particles and add them to the nearest DOF
!==================================================================================================================================
SUBROUTINE CalcSourcePart(Ut)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Globals ,ONLY: VECNORM
USE MOD_Mesh_Vars        ,ONLY: Elem_xGP,sJ,nElems,offsetElem
#if FV_ENABLED
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisVolume
USE MOD_FV_Vars          ,ONLY: FV_Vdm,FV_Elems
#endif
USE MOD_Particle_Vars    ,ONLY: Species,PartSpecies,PartState,Pt,PEM,PDM
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem,iPart,ijk(3)
REAL                :: Ut_src(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: Fp(3)
#if FV_ENABLED
REAL                :: Ut_src2(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
REAL                :: min_distance_glob,min_distance_loc
!==================================================================================================================================

DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! Calculate particle force
    Fp(1:3) = Pt(1:3,iPart)*Species(PartSpecies(iPart))%MassIC
    ! Determine nearest DOF
    iElem = PEM%Element(iPart)-offsetElem
    min_distance_glob = VECNORM(Elem_xGP(:,0,0,0,iElem)-PartState(1:3,iPart))
    ijk(:) = 0
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      min_distance_loc = VECNORM(Elem_xGP(:,i,j,k,iElem)-PartState(1:3,iPart))
      IF (min_distance_loc .LT. min_distance_glob) THEN; ijk(:) = (/i,j,k/); min_distance_glob = min_distance_loc; END IF
    END DO; END DO; END DO
    ! Add source term
    Ut_src(DENS     ,ijk(1),ijk(2),ijk(3)) = 0.
    Ut_src(MOM1:MOM3,ijk(1),ijk(2),ijk(3)) = Fp
    Ut_src(ENER     ,ijk(1),ijk(2),ijk(3)) = DOT_PRODUCT(Fp,PartState(4:6,iPart))
#if FV_ENABLED
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      Ut(:,ijk(1),ijk(2),ijk(3),iElem) = Ut(:,ijk(1),ijk(2),ijk(3),iElem)+&
                                         Ut_src2(:,ijk(1),ijk(2),ijk(3))/sJ(ijk(1),ijk(2),ijk(3),iElem,1)
    ELSE
#endif
      Ut(:,ijk(1),ijk(2),ijk(3),iElem) = Ut(:,ijk(1),ijk(2),ijk(3),iElem)+&
                                         Ut_src(:,ijk(1),ijk(2),ijk(3))/sJ(ijk(1),ijk(2),ijk(3),iElem,0)
#if FV_ENABLED
    END IF
#endif
  END IF
END DO ! iPart

END SUBROUTINE CalcSourcePart


END MODULE MOD_part_RHS
