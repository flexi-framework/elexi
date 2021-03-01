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
#include "eos.h"

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

PUBLIC :: CalcPartRHS
!==================================================================================================================================

CONTAINS

SUBROUTINE CalcPartRHS(&
#if USE_BASSETFORCE
  dt,iStage)
#else
  )
#endif /* USE_BASSETFORCE */
!===================================================================================================================================
! Computes the acceleration from the drag force with respect to the species data and velocity
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Interpolation_Vars,  ONLY: FieldAtParticle
#if USE_EXTEND_RHS
USE MOD_Particle_Interpolation_Vars,  ONLY: GradAtParticle,TimeDerivAtParticle
USE MOD_Particle_TimeDisc_Vars,       ONLY: UseManualTimestep
#endif
USE MOD_Particle_Vars,                ONLY: PDM, Pt
!#if USE_RW
!USE MOD_Particle_RandomWalk_Vars,     ONLY: RWTime
!USE MOD_Particle_Vars,                ONLY: Species,PartSpecies,TurbPartState
!USE_MOD_Timedisc_Vars,                ONLY: t
!#endif
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
#if USE_EXTEND_RHS
REAL                             :: Fd(1:3)
#endif
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
#if USE_EXTEND_RHS
    ! Calculate the drag (and gravity) force
    Fd(1:3)       = ParticlePush(iPart,FieldAtParticle(PRIM,iPart))
    ! nullify time derivative if steady state computation
    IF(UseManualTimestep) TimeDerivAtParticle(1:RHS_DERIVATIVE,iPart) = 0.
    ! Calculate other RHS forces and add all forces to compute the particle push
    Pt(1:3,iPart) = ParticlePushExtend(iPart,FieldAtParticle(PRIM,iPart),&
                                       GradAtParticle(1:RHS_LIFT,3,iPart),TimeDerivAtParticle(1:RHS_DERIVATIVE,iPart),Fd&
#if USE_BASSETFORCE
                                      ,dt,iStage)
#else
                                      )
#endif /* USE_BASSETFORCE */
#else
    ! Calculate the drag (and gravity) force
    Pt(1:3,iPart) = ParticlePush(iPart,FieldAtParticle(PRIM,iPart))
#endif
  END IF
END DO

END SUBROUTINE CalcPartRHS


FUNCTION ParticlePush(PartID,FieldAtParticle)
!===================================================================================================================================
! Push due to Stoke's drag and source terms (gravity)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals
USE MOD_EOS_Vars,          ONLY : mu0
USE MOD_Particle_Vars,     ONLY : Species, PartSpecies, PartGravity
USE MOD_Particle_Vars,     ONLY : PartState, RepWarn
USE MOD_Particle_Vars,     ONLY : TurbPartState
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(PRIM)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: ParticlePush(1:3)           ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Fdm(1:3)
REAL                :: Rep                         ! Reynolds and Mach number of particle
!REAL                :: velosqp                    ! v^2 particle
!REAL                :: velosqf                    ! v^2 fluid
REAL                :: udiff(3)
REAL                :: Cd                          ! Drag coefficient
REAL                :: staup                       ! Inverse of the particle relaxation time
REAL                :: k1,k2,k3
REAL                :: mu                          ! viscosity
REAL,PARAMETER      :: epsilonRHS=1.0
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

CASE(RHS_VINKOVIC,RHS_WANG)
!===================================================================================================================================
! Calculation according to Vinkovic [2006] or Wang [1996]
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL ABORT(__STAMP__,'Particle tracking with Wang [1996] or Vinkovic [2006] requires mu to be set!')

! Assume spherical particles for now
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = FieldAtParticle(VELV) + TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)                             - PartState(PART_VELV,PartID)
END IF

Rep     = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu

! Empirical relation of nonlinear drag from Clift et al. (1978)
Cd  = 1. + 0.15*Rep**0.687
IF((Species(PartSpecies(PartID))%RHSMethod.EQ.RHS_VINKOVIC).AND.(Rep .LT. epsilonRHS)) Cd  = 1.

! Warn when outside valid range of Wang model
IF(Rep.GT.800) THEN
  IF (RepWarn.EQV..FALSE.) THEN
    SWRITE(UNIT_StdOut,*) 'WARNING: Rep',Rep,'> 800, drag coefficient may not be accurate.'
    RepWarn=.TRUE.
  ENDIF
ENDIF

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./Species(PartSpecies(PartID))%DiameterIC**2

Fdm      = udiff * staup * Cd

! Add gravity if required
IF(ANY(PartGravity.NE.0)) THEN
  Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)
ENDIF

CASE(RHS_JACOBS)
!===================================================================================================================================
! Calculation according to Jacobs [2003]
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL ABORT(__STAMP__,'Particle tracking with Jacobs [2003] requires mu to be set!')

! Assume spherical particles for now
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = FieldAtParticle(VELV) + TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)                             - PartState(PART_VELV,PartID)
END IF

Rep     = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu

! Drag coefficient according to Putnam (1961)
Cd = 1. + (Rep**2./3.)/6.
! High Re correction according to Putnam et al. (1961)
IF(Rep .GT. 1000) Cd = 0.0183*Rep

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./Species(PartSpecies(PartID))%DiameterIC**2

Fdm       = udiff*staup * Cd

! Add gravity if required
IF(ANY(PartGravity.NE.0)) THEN
  Fdm  = Fdm + PartGravity * (1.-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)
END IF

CASE(RHS_HAIDER)
!===================================================================================================================================
! Calculation according to Haider and Levenspiel [1989]
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL ABORT(__STAMP__,'Particle tracking with Haider [1989] requires mu to be set!')

! Assume spherical particles for now
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = FieldAtParticle(VELV) + TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)                             - PartState(PART_VELV,PartID)
END IF

Rep     = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu

! Empirical relation from Haider and Levenspiel (1989) valid for Rep<2.6e5
k1 = EXP(2.3288-6.4581*Species(PartSpecies(PartID))%SphericityIC+2.4486*Species(PartSpecies(PartID))%SphericityIC**2)
k2 = EXP(4.905-13.8944*Species(PartSpecies(PartID))%SphericityIC+18.4222*Species(PartSpecies(PartID))%SphericityIC**2-&
  10.2599*Species(PartSpecies(PartID))%SphericityIC**3)
k3 = EXP(1.4681+12.2584*Species(PartSpecies(PartID))%SphericityIC-20.7322*Species(PartSpecies(PartID))%SphericityIC**2+&
  15.8855*Species(PartSpecies(PartID))%SphericityIC**3)
Cd = (1+k1*Rep**(0.0964+0.5565*Species(PartSpecies(PartID))%SphericityIC))+Rep**2*1./24*k2/(Rep+k3)

! Warn when outside valid range of Haider model
IF(Species(PartSpecies(PartID))%SphericityIC.LT.0.670) THEN
  SWRITE(UNIT_StdOut,*) 'WARNING: SphericityIC',Species(PartSpecies(PartID))%SphericityIC,'< 0.670, drag coefficient may not be accurate.'
ENDIF

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./Species(PartSpecies(PartID))%DiameterIC**2

Fdm       = udiff*staup * Cd

! Add gravity if required
IF(ANY(PartGravity.NE.0)) THEN
    Fdm  = Fdm + PartGravity * (1-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)
ENDIF

CASE(RHS_HOELZER)
!===================================================================================================================================
! Calculation according to Hoelzer [2008]
!===================================================================================================================================
IF(ISNAN(mu) .OR. (mu.EQ.0)) CALL ABORT(__STAMP__,'Particle tracking with Hoelzer [2008] requires mu to be set!')

! Assume spherical particles for now
IF(ALLOCATED(TurbPartState)) THEN
  udiff(1:3) = FieldAtParticle(VELV) + TurbPartState(1:3,PartID) - PartState(PART_VELV,PartID)
ELSE
  udiff(1:3) = FieldAtParticle(VELV)                             - PartState(PART_VELV,PartID)
END IF

Rep     = VECNORM(udiff(1:3))*Species(PartSpecies(PartID))%DiameterIC*FieldAtParticle(DENS)/mu

Cd    =  1./3. * 1./SQRT(Species(PartSpecies(PartID))%SphericityIC) + 2./3. * 1./SQRT(Species(PartSpecies(PartID))%SphericityIC)+&
         1./8. * 1./(Species(PartSpecies(PartID))%SphericityIC**(3./4.)) +&
         Rep/24. * 0.4210**(0.4*(-LOG(Species(PartSpecies(PartID))%SphericityIC))**0.2) *1./Species(PartSpecies(PartID))%SphericityIC

! Particle relaxation time
staup    = (18.*mu) * 1./Species(PartSpecies(PartID))%DensityIC * 1./Species(PartSpecies(PartID))%DiameterIC**2

Fdm       = udiff*staup * Cd

! Add gravity if required
IF(ANY(PartGravity.NE.0)) THEN
    Fdm  = Fdm + PartGravity * (1-FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC)
ENDIF

CASE DEFAULT
  CALL ABORT(__STAMP__, 'No valid RHS method given. Species',IntInfo=PartSpecies(PartID))

END SELECT

ParticlePush = Fdm
!WRITE (*, *) 'Fdm:', Fdm

END FUNCTION ParticlePush

#if USE_EXTEND_RHS
FUNCTION ParticlePushExtend(PartID,FieldAtParticle,GradAtParticle,TimeDerivAtParticle,Fd&
#if USE_BASSETFORCE
  ,dt,iStage)
#else
  )
#endif
!===================================================================================================================================
! Push due to additional forces, e.g. Basset force, lift force, added mass effect, viscous and pressure forces
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals
USE MOD_EOS_Vars,               ONLY: mu0
USE MOD_Particle_Vars,          ONLY: Species, PartSpecies
USE MOD_Particle_Vars,          ONLY: PartState
USE MOD_Particle_Vars,          ONLY: TurbPartState
USE MOD_PreProc,                ONLY: PP_pi
USE MOD_Equation_Vars,          ONLY: s43
#if USE_BASSETFORCE
USE MOD_Particle_Vars,          ONLY: durdt, N_Basset, bIter
USE MOD_TimeDisc_Vars,          ONLY: nRKStages, RKC
#endif /* USE_BASSETFORCE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: PartID
REAL,INTENT(IN)             :: FieldAtParticle(PRIM)
REAL,INTENT(IN)             :: GradAtParticle(1:RHS_LIFT,1:3)
REAL,INTENT(IN)             :: TimeDerivAtParticle(1:RHS_DERIVATIVE)
REAL,INTENT(IN)             :: Fd(1:3)
#if USE_BASSETFORCE
REAL,INTENT(IN)             :: dt
INTEGER,INTENT(IN),OPTIONAL :: iStage
#endif /* USE_BASSETFORCE */
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: ParticlePushExtend(1:3)     ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: Pt(1:3)                     ! Particle push divided by the particle mass
REAL                     :: udiff(3)                    ! velocity difference
REAL                     :: mu                          ! viscosity
REAL                     :: facm,facp                   ! factor divided by the particle mass
#if USE_LIFTFORCE
REAL                     :: Flm(1:3)                    ! lift force divided by the particle mass
REAL                     :: rotu(3), rotudiff(3)        ! curl product of velocity and velocity difference
REAL                     :: dotp                        ! dot_product
#endif /* USE_LIFTFORCE */
#if (USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE)
REAL                     :: DuDt(1:3)                   ! viscous and pressure forces divided by the particle mass
#endif
#if USE_BASSETFORCE
REAL                     :: Fbm(1:3)                    ! Basset force
REAL,PARAMETER           :: s32=3./2.
REAL                     :: RKdtFrac
INTEGER                  :: k,kIndex
REAL                     :: dufdt(1:3)
#endif /* USE_BASSETFORCE */
!===================================================================================================================================

SELECT CASE(Species(PartSpecies(PartID))%RHSMethod)
CASE(RHS_NONE,RHS_CONVERGENCE,RHS_TRACER)
  ParticlePushExtend = Fd
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

Pt(1:3) = 0.
! factor before left hand side
facp = 1.

!===================================================================================================================================
! Calculate the Saffman lift force:
! Saffman, P.G.: The lift on a small sphere in a slow shear flow. Journal of Fluid Mechanics,
! pp. 385–400, 1965. 10.1017/S0022112065000824.
!===================================================================================================================================
#if USE_LIFTFORCE
! Calculate the factor
facm = 9.69./(Species(PartSpecies(PartID))%DensityIC*Species(PartSpecies(PartID))%DiameterIC*PP_PI)
! Calculate the rotation: \nabla x u
rotu     = (/GradAtParticle(3,2)-GradAtParticle(2,3),&
             GradAtParticle(1,3)-GradAtParticle(3,1),&
             GradAtParticle(2,1)-GradAtParticle(1,2)/)
! Calculate the rotation: (\nabla x u) x udiff
rotudiff = (/rotu(2)*udiff(3)-rotu(3)*udiff(2),&
             rotu(3)*udiff(1)-rotu(1)*udiff(3),&
             rotu(1)*udiff(2)-rotu(2)*udiff(1)/)

dotp    = MAX(DOT_PRODUCT(rotu(:),rotu(:)),0.001)
Flm(:)  = SQRT(2*FieldAtParticle(DENS)*mu) * 1./dotp*rotudiff(:)

Pt(1:3) = Pt(1:3) + facm*Flm(1:3)

!WRITE(*,*) 'LIFTFORCE', (facm*Flm(1:3))/Fd*100

#endif /* USE_LIFTFORCE */

!===================================================================================================================================
! Calculate the viscous and pressure forces:
! 1/\rho_p Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
!===================================================================================================================================
#if USE_UNDISTFLOW

facm = 1./Species(PartSpecies(PartID))%DensityIC
! Material derivative Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
DuDt(1) = TimeDerivAtParticle(MOM1) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL1,:))&
                             + FieldAtParticle(VEL1) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
DuDt(2) = TimeDerivAtParticle(MOM2) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL2,:))&
                             + FieldAtParticle(VEL2) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
DuDt(3) = TimeDerivAtParticle(MOM3) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL3,:))&
                                + FieldAtParticle(VEL3) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))

Pt(1:3)  = Pt(1:3) + facm * DuDt(1:3)
!WRITE (*, *) 'UNDISTFLOW:', (facm * DuDt(1:3))/Fd*100
#endif /* USE_UNDISTFLOW */

!===================================================================================================================================
! Calculate the viscous and pressure forces:
! 1/\rho_p Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
!===================================================================================================================================
#if USE_VIRTUALMASS

#if !USE_UNDISTFLOW
! Material derivative Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
DuDt(1) = TimeDerivAtParticle(MOM1) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL1,:))&
                                    + FieldAtParticle(VEL1) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
DuDt(2) = TimeDerivAtParticle(MOM2) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL2,:))&
                                    + FieldAtParticle(VEL2) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
DuDt(3) = TimeDerivAtParticle(MOM3) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL3,:))&
                                    + FieldAtParticle(VEL3) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
#endif

Pt(1:3)  = Pt(1:3) + 0.5 * facm * DuDt(1:3) - 0.5 * facm * PartState(PART_VELV,PartID) * TimeDerivAtParticle(DENS)

facp     = facp + 0.5*FieldAtParticle(DENS)/Species(PartSpecies(PartID))%DensityIC
!WRITE (*, *) 'VIRTUALMASS:', (0.5 * facm * DuDt(1:3)- 0.5 * facm * PartState(PART_VELV,PartID) * TimeDerivAtParticle(DENS))/Fd*100
#endif /* USE_VIRTUALMASS */

!===================================================================================================================================
! Calculate the Basset force / history terms:
! van Hinsberg, M.A., ten Thije Boonkkamp, J.H., Clercx, H.J.: An efficient, second or-
! der method for the approximation of the Basset history force. Journal of Computational Physics,
! pp. 1465–1478, 2011. 10.1016/j.jcp.2010.11.014.
!===================================================================================================================================
#if USE_BASSETFORCE
bIter = bIter + 1
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

facm = 9. * 1./(Species(PartSpecies(PartID))%DiameterIC*Species(PartSpecies(PartID))%DensityIC)&
          * SQRT(mu/(FieldAtParticle(DENS)*PP_pi)) * SQRT(RKdtFrac*dt)

#if !(USE_UNDISTFLOW || USE_VIRTUALMASS)
! Material derivative Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
! D(\rho u)/Dt = \partial(\rho u) / \partial(t) + (\rho u) \cdot \nabla(u) + u \cdot (\nabla(rho) u)
DuDt(1) = TimeDerivAtParticle(MOM1) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL1,:))&
                                    + FieldAtParticle(VEL1) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
DuDt(2) = TimeDerivAtParticle(MOM2) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL2,:))&
                                    + FieldAtParticle(VEL2) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
DuDt(3) = TimeDerivAtParticle(MOM3) + FieldAtParticle(DENS) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(VEL3,:))&
                                    + FieldAtParticle(VEL3) * DOT_PRODUCT(FieldAtParticle(VELV),GradAtParticle(DENS,:))
#endif

!WRITE (*, *) 'TimeDerivAtParticle(DENS:MOM3):', TimeDerivAtParticle(DENS:MOM3)

! Index for previous data
kIndex = INT(MIN(N_Basset, bIter)*3)
! copy previous data
IF(bIter.GT.N_Basset) durdt(1:kIndex-2,PartID) = durdt(4:kIndex,PartID)
! d(\rho u)/dt = D(\rho u)/Dt - udiff * (\rho \nabla(u) + u \nabla(\rho))
dufdt(1) = DuDt(1) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(VEL1,:))&
                   - DOT_PRODUCT(udiff(:),FieldAtParticle(VEL1)*GradAtParticle(DENS,:))
dufdt(2) = DuDt(2) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(VEL2,:))&
                   - DOT_PRODUCT(udiff(:),FieldAtParticle(VEL2)*GradAtParticle(DENS,:))
dufdt(3) = DuDt(3) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(VEL3,:))&
                   - DOT_PRODUCT(udiff(:),FieldAtParticle(VEL3)*GradAtParticle(DENS,:))

! d(\rho udiff)/dt = d(\rho u)/dt - \rho (dv_p/dt) - v_p (d\rho/dt)
durdt(kIndex-2:kIndex,PartID) = dufdt(:) - PartState(PART_VELV,PartID) * TimeDerivAtParticle(DENS)

Fbm = s43 * durdt(kIndex-2:kIndex,PartID) + durdt(1:3,PartID) * (N_Basset-s43)/((N_Basset-1)*SQRT(REAL(N_Basset-1))+(N_Basset-s32)*SQRT(REAL(N_Basset)))
DO k=3,kIndex-3,3
  Fbm = Fbm + durdt(kIndex-2-k:kIndex-k,PartID) * ((k+s43)/((k+1)*SQRT(REAL(k+1))+(k+s32)*SQRT(REAL(k)))+(k-s43)/((k-1)*SQRT(REAL(k-1))+(k-s32)*SQRT(REAL(k))))
END DO

facp     = facp + s43 * facm * FieldAtParticle(DENS)

Pt(1:3)  = Pt(1:3) + facm * Fbm
!WRITE (*, *) 'BASSET_FORCE:', (facm * Fbm)/Fd*100

! Correct durdt with particle push
durdt(kIndex-2:kIndex,PartID) = durdt(kIndex-2:kIndex,PartID) - FieldAtParticle(DENS) * (Pt(1:3) + Fd(1:3)) * 1./facp
#endif /* USE_BASSETFORCE */

ParticlePushExtend(1:3) = (Pt(1:3) + Fd(1:3)) * 1./facp

END FUNCTION ParticlePushExtend
#endif

END MODULE MOD_part_RHS
