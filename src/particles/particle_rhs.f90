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

#if USE_EXTEND_RHS
INTERFACE tauRHS
  MODULE PROCEDURE tauRHS
END INTERFACE
#endif /* USE_EXTEND_RHS */

PUBLIC :: CalcPartRHS
#if USE_EXTEND_RHS
PUBLIC :: tauRHS
#endif /* USE_EXTEND_RHS */
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
USE MOD_Particle_Interpolation_Vars,  ONLY: GradAtParticle
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
    ! Calculate other RHS forces and add all forces to compute the particle push
    Pt(1:3,iPart) = ParticlePushExtend(iPart,FieldAtParticle(PRIM,iPart),&
                                       GradAtParticle(1:RHS_GRAD,1:3,iPart),Fd&
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
USE MOD_Particle_Vars,     ONLY : Species, PartSpecies, PartGravity
USE MOD_Particle_Vars,     ONLY : PartState, RepWarn
USE MOD_Particle_Vars,     ONLY : TurbPartState
USE MOD_Viscosity
USE MOD_Equation_Vars,     ONLY : s13,s23
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

Cd    =  s13 * 1./SQRT(Species(PartSpecies(PartID))%SphericityIC) + s23 * 1./SQRT(Species(PartSpecies(PartID))%SphericityIC)+&
         SQRT(Rep)/8. * 1./(Species(PartSpecies(PartID))%SphericityIC**(3./4.)) +&
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
FUNCTION ParticlePushExtend(PartID,FieldAtParticle,GradAtParticle,Fdm&
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
USE MOD_Particle_Vars,          ONLY: Species, PartSpecies
USE MOD_Particle_Vars,          ONLY: PartState, TurbPartState
#if ANALYZE_RHS
USE MOD_Particle_Vars,          ONLY: tWriteRHS,FileName_RHS,dtWriteRHS
USE MOD_TimeDisc_Vars,          ONLY: t
USE MOD_Output,                 ONLY: OutputToFile
#endif /* ANALYZE_RHS */
USE MOD_PreProc,                ONLY: PP_pi
USE MOD_Equation_Vars,          ONLY: s43
USE MOD_Viscosity
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
REAL,INTENT(IN)             :: GradAtParticle(1:RHS_GRAD,1:3)
REAL,INTENT(IN)             :: Fdm(1:3)
#if USE_BASSETFORCE
REAL,INTENT(IN)             :: dt
INTEGER,INTENT(IN),OPTIONAL :: iStage
#endif /* USE_BASSETFORCE */
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: ParticlePushExtend(1:3)     ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: Pt(1:3)
REAL                     :: udiff(3)                    ! velocity difference
REAL                     :: mu                          ! viscosity
REAL                     :: prefactor,globalfactor                   ! factor divided by the particle mass
REAL                     :: Flm(1:3)                    ! lift force divided by the particle mass
#if USE_LIFTFORCE
REAL                     :: rotu(3), rotudiff(3)        ! curl product of velocity and velocity difference
REAL                     :: dotp                        ! dot_product
#endif /* USE_LIFTFORCE */
REAL                     :: Fum(1:3)
REAL                     :: Fvm(1:3)
REAL                     :: Fbm(1:3)                    ! Basset force
#if (USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE)
REAL                     :: DuDt(1:3)                   ! viscous and pressure forces divided by the particle mass
#endif
#if USE_BASSETFORCE
REAL,PARAMETER           :: s32=3./2.
REAL                     :: RKdtFrac
INTEGER                  :: k,kIndex
REAL                     :: dufdt(1:3)
#endif /* USE_BASSETFORCE */
!===================================================================================================================================

SELECT CASE(Species(PartSpecies(PartID))%RHSMethod)
CASE(RHS_NONE,RHS_CONVERGENCE,RHS_TRACER)
  ParticlePushExtend = Fdm
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
Flm = 0.; Fbm = 0.; Fvm=0.; Fum=0.
! factor before left hand side
globalfactor = 1.

!===================================================================================================================================
! Calculate the Saffman lift force:
! Saffman, P.G.: The lift on a small sphere in a slow shear flow. Journal of Fluid Mechanics,
! pp. 385–400, 1965. 10.1017/S0022112065000824.
!===================================================================================================================================
#if USE_LIFTFORCE
IF (Species(PartSpecies(PartID))%CalcLiftForce) THEN
  ! Calculate the factor
  prefactor = 9.69/(Species(PartSpecies(PartID))%DensityIC*Species(PartSpecies(PartID))%DiameterIC*PP_PI)
  ! Calculate the rotation: \nabla x u
  rotu     = (/GradAtParticle(RHS_GRADVEL3,2)-GradAtParticle(RHS_GRADVEL2,3),&
               GradAtParticle(RHS_GRADVEL1,3)-GradAtParticle(RHS_GRADVEL3,1),&
               GradAtParticle(RHS_GRADVEL2,1)-GradAtParticle(RHS_GRADVEL1,2)/)
  ! Calculate the rotation: (\nabla x u) x udiff
  rotudiff = (/rotu(2)*udiff(3)-rotu(3)*udiff(2),&
               rotu(3)*udiff(1)-rotu(1)*udiff(3),&
               rotu(1)*udiff(2)-rotu(2)*udiff(1)/)

  dotp    = MAX(SQRT(DOT_PRODUCT(rotu(:),rotu(:))),0.001)
  Flm(:)  = SQRT(2*FieldAtParticle(DENS)*mu * 1./dotp)*rotudiff(:)

  Flm     = Flm * prefactor
END IF
#endif /* USE_LIFTFORCE */

#if USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE
IF (Species(PartSpecies(PartID))%CalcUndisturbedFlow.OR.Species(PartSpecies(PartID))%CalcVirtualMass.OR.&
    Species(PartSpecies(PartID))%CalcBassetForce) THEN
  ! Material derivative \rho D(u_i)/Dt = - \nabla p + \nabla \cdot \tau
  DuDt(1)  = - GradAtParticle(RHS_GRADPRES,1) + mu*GradAtParticle(RHS_GRADTAU,1)
  DuDt(2)  = - GradAtParticle(RHS_GRADPRES,2) + mu*GradAtParticle(RHS_GRADTAU,2)
  DuDt(3)  = - GradAtParticle(RHS_GRADPRES,3) + mu*GradAtParticle(RHS_GRADTAU,3)
END IF
#endif /* USE_UNDISTFLOW || USE_VIRTUALMASS || USE_BASSETFORCE */

!===================================================================================================================================
! Calculate the viscous and pressure forces:
! 1/\rho_p Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
!===================================================================================================================================
#if USE_UNDISTFLOW
IF (Species(PartSpecies(PartID))%CalcUndisturbedFlow) THEN
  prefactor = 1./Species(PartSpecies(PartID))%DensityIC

  Fum(1:3) = prefactor * DuDt(1:3)
END IF
#endif /* USE_UNDISTFLOW */

!===================================================================================================================================
! Calculate the viscous and pressure forces:
! 1/\rho_p Du_i/Dt = \partial \rho u_i / \partial t + u_j \partial \rho u_i / \partial x_j
!===================================================================================================================================
#if USE_VIRTUALMASS
IF (Species(PartSpecies(PartID))%CalcVirtualMass) THEN
  prefactor = 0.5/Species(PartSpecies(PartID))%DensityIC

  Fvm(1:3) = prefactor * DuDt(1:3)

  ! Add to global scaling factor
  globalfactor     = globalfactor + prefactor*FieldAtParticle(DENS)
END IF
#endif /* USE_VIRTUALMASS */

!===================================================================================================================================
! Calculate the Basset force / history terms:
! van Hinsberg, M.A., ten Thije Boonkkamp, J.H., Clercx, H.J.: An efficient, second or-
! der method for the approximation of the Basset history force. Journal of Computational Physics,
! pp. 1465–1478, 2011. 10.1016/j.jcp.2010.11.014.
!===================================================================================================================================
#if USE_BASSETFORCE
IF (Species(PartSpecies(PartID))%CalcBassetForce) THEN
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

  ! Scaling factor
  prefactor = 9./(Species(PartSpecies(PartID))%DiameterIC*Species(PartSpecies(PartID))%DensityIC)&
            * SQRT(mu/(FieldAtParticle(DENS)*PP_pi)) * SQRT(RKdtFrac)

  ! Index for previous data
  kIndex = INT(MIN(N_Basset, bIter)*3)
  ! copy previous data
  IF(bIter.GT.N_Basset) durdt(1:kIndex-3,PartID) = durdt(4:kIndex,PartID)
  ! \rho d(u)/dt = \rho D(u)/Dt - udiff * (\rho \nabla(u))
  dufdt(1) = DuDt(1) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(RHS_GRADVEL1,:))
  dufdt(2) = DuDt(2) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(RHS_GRADVEL2,:))
  dufdt(3) = DuDt(3) - DOT_PRODUCT(udiff(:),FieldAtParticle(DENS)*GradAtParticle(RHS_GRADVEL3,:))

  ! \rho d(udiff)/dt = \rho d(u)/dt - \rho (dv_p/dt)
  durdt(kIndex-2:kIndex,PartID) = dufdt(:)

  Fbm = s43 * durdt(kIndex-2:kIndex,PartID) + durdt(1:3,PartID) * (N_Basset-s43)/((N_Basset-1)*SQRT(REAL(N_Basset-1))+(N_Basset-s32)*SQRT(REAL(N_Basset)))
  DO k=3,kIndex-3,3
    Fbm = Fbm + durdt(kIndex-2-k:kIndex-k,PartID) * ((k+s43)/((k+1)*SQRT(REAL(k+1))+(k+s32)*SQRT(REAL(k)))+(k-s43)/((k-1)*SQRT(REAL(k-1))+(k-s32)*SQRT(REAL(k))))
  END DO

  globalfactor     = globalfactor + s43 * prefactor * FieldAtParticle(DENS)

  Fbm(1:3) = prefactor * Fbm(1:3)

  ! Correct durdt with particle push
  durdt(kIndex-2:kIndex,PartID) = durdt(kIndex-2:kIndex,PartID) - FieldAtParticle(DENS) * (Flm + Fum + Fvm + Fbm + Fdm) * 1./globalfactor

  Pt = (Flm + Fum + Fvm + Fbm + Fdm) * 1./globalfactor

  ! Correct durdt with particle push
  durdt(kIndex-2:kIndex,PartID) = durdt(kIndex-2:kIndex,PartID) - FieldAtParticle(DENS) * Pt
ELSE
#endif /* USE_BASSETFORCE */
  Pt = (Flm + Fum + Fvm + Fbm + Fdm) * 1./globalfactor
#if USE_BASSETFORCE
END IF
#endif /* USE_BASSETFORCE */

! Output RHS to file
#if ANALYZE_RHS
IF(dtWriteRHS.GT.0.0)THEN
  IF(tWriteRHS-t.LE.dt*(1.+1.E-4))THEN
    CALL OutputToFile(FileName_RHS,(/t/),(/16,1/),(/REAL(PartSpecies(PartID)),Fdm(1:3),Flm(1:3),Fum(1:3),Fvm(1:3),Fbm(1:3)/))
    tWriteRHS = tWriteRHS + dtWriteRHS
  END IF
END IF
#endif

ParticlePushExtend(1:3) = Pt

END FUNCTION ParticlePushExtend

SUBROUTINE tauRHS(U,divtau,gradp)
!===================================================================================================================================
! Compute tau
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,      ONLY: s13
USE MOD_Lifting_Vars,       ONLY: gradUx,gradUy,gradUz
USE MOD_Lifting_Vars,       ONLY: gradUx_master,gradUx_slave
USE MOD_Lifting_Vars,       ONLY: gradUy_master,gradUy_slave
USE MOD_Lifting_Vars,       ONLY: gradUz_master,gradUz_slave
USE MOD_Lifting_BR1_gen,    ONLY: Lifting_BR1_gen
USE MOD_Mesh_Vars,          ONLY: nElems,nSides
USE MOD_EoS,                ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)             :: U(CONS,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(OUT)            :: divtau(1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(OUT)            :: gradp(1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE            :: gradUx2(:,:,:,:,:,:),gradUy2(:,:,:,:,:,:),gradUz2(:,:,:,:,:,:)
REAL,ALLOCATABLE            :: gradUx_master_loc(:,:,:,:), gradUx_slave_loc(:,:,:,:)
REAL,ALLOCATABLE            :: gradUy_master_loc(:,:,:,:), gradUy_slave_loc(:,:,:,:)
REAL,ALLOCATABLE            :: gradUz_master_loc(:,:,:,:), gradUz_slave_loc(:,:,:,:)
REAL,ALLOCATABLE            :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE            :: gradp_local(:,:,:,:,:,:)
INTEGER                     :: i,j,k,iElem
!===================================================================================================================================
ALLOCATE(gradUx2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
ALLOCATE(gradUy2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
ALLOCATE(gradUz2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
gradUx2 = 0.
gradUy2 = 0.
gradUz2 = 0.
ALLOCATE(gradUx_master_loc(1:3,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUx_slave_loc( 1:3,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_master_loc(1:3,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUy_slave_loc( 1:3,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_master_loc(1:3,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(gradUz_slave_loc( 1:3,0:PP_N,0:PP_NZ,1:nSides))

gradUx_master_loc = gradUx_master(LIFT_VELV,:,:,:)
gradUx_slave_loc  = gradUx_slave(LIFT_VELV,:,:,:)
gradUy_master_loc = gradUy_master(LIFT_VELV,:,:,:)
gradUy_slave_loc  = gradUy_slave(LIFT_VELV,:,:,:)
gradUz_master_loc = gradUz_master(LIFT_VELV,:,:,:)
gradUz_slave_loc  = gradUz_slave(LIFT_VELV,:,:,:)


! Calculate the second gradient of the velocity and output \nabla \cdot tau
CALL Lifting_BR1_gen(3,3,gradUx(LIFT_VELV,:,:,:,:),gradUx_master_loc,gradUx_slave_loc,&
                    gradUx2(1,:,:,:,:,:),gradUx2(2,:,:,:,:,:),gradUx2(3,:,:,:,:,:))
CALL Lifting_BR1_gen(3,3,gradUy(LIFT_VELV,:,:,:,:),gradUy_master_loc,gradUy_slave_loc,&
                    gradUy2(1,:,:,:,:,:),gradUy2(2,:,:,:,:,:),gradUy2(3,:,:,:,:,:))
CALL Lifting_BR1_gen(3,3,gradUz(LIFT_VELV,:,:,:,:),gradUz_master_loc,gradUz_slave_loc,&
                    gradUz2(1,:,:,:,:,:),gradUz2(2,:,:,:,:,:),gradUz2(3,:,:,:,:,:))

! u_xx + u_yy + u_zz + 1/3 * (u_xx+v_yx+w_zx)
!divtau(1,:,:,:,:) = gradUx2(1,LIFT_VEL1,:,:,:,:) + gradUy2(2,LIFT_VEL1,:,:,:,:) + gradUz2(3,LIFT_VEL1,:,:,:,:) + &
!                    s13 * (gradUx2(1,LIFT_VEL1,:,:,:,:) + gradUy2(1,LIFT_VEL2,:,:,:,:) + gradUz2(1,LIFT_VEL3,:,:,:,:))
divtau(1,:,:,:,:) = gradUx2(1,1,:,:,:,:) + gradUy2(2,1,:,:,:,:) + gradUz2(3,1,:,:,:,:) + &
                    s13 * (gradUx2(1,1,:,:,:,:) + gradUy2(1,2,:,:,:,:) + gradUz2(1,3,:,:,:,:))
! v_xx + v_yy + v_zz + 1/3 * (u_xy+v_yy+w_zy)
!divtau(2,:,:,:,:) = gradUx2(1,LIFT_VEL2,:,:,:,:) + gradUy2(2,LIFT_VEL2,:,:,:,:) + gradUz2(3,LIFT_VEL2,:,:,:,:) + &
!                    s13 * (gradUx2(2,LIFT_VEL1,:,:,:,:) + gradUy2(2,LIFT_VEL2,:,:,:,:) + gradUz2(2,LIFT_VEL3,:,:,:,:))
divtau(2,:,:,:,:) = gradUx2(1,2,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(3,2,:,:,:,:) + &
                    s13 * (gradUx2(2,1,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(2,3,:,:,:,:))
! w_xx + w_yy + w_zz + 1/3 * (u_xy+v_yy+w_zy)
!divtau(3,:,:,:,:) = gradUx2(1,LIFT_VEL3,:,:,:,:) + gradUy2(2,LIFT_VEL3,:,:,:,:) + gradUz2(3,LIFT_VEL3,:,:,:,:) + &
!                    s13 * (gradUx2(3,LIFT_VEL1,:,:,:,:) + gradUy2(3,LIFT_VEL2,:,:,:,:) + gradUz2(3,LIFT_VEL3,:,:,:,:))
divtau(3,:,:,:,:) = gradUx2(1,3,:,:,:,:) + gradUy2(2,3,:,:,:,:) + gradUz2(3,3,:,:,:,:) + &
                    s13 * (gradUx2(3,1,:,:,:,:) + gradUy2(3,2,:,:,:,:) + gradUz2(3,3,:,:,:,:))

DEALLOCATE(gradUx2,gradUy2,gradUz2)
DEALLOCATE(gradUx_master_loc,gradUx_slave_loc)
DEALLOCATE(gradUy_master_loc,gradUy_slave_loc)
DEALLOCATE(gradUz_master_loc,gradUz_slave_loc)

! Calculate pressure gradient
ALLOCATE(U_local(PRIM,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
ALLOCATE(gradp_local(1,3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))

DO iElem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  CALL ConsToPrim(U_local(:,i,j,k,iElem),U(:,i,j,k,iElem))
END DO; END DO; END DO; END DO

CALL Lifting_BR1_gen(1,1,U_local(PRES:PRES,:,:,:,:),gradp_local(:,1,:,:,:,:),gradp_local(:,2,:,:,:,:),gradp_local(:,3,:,:,:,:))
gradp = gradp_local(1,:,:,:,:,:)

DEALLOCATE(U_local,gradp_local)

END SUBROUTINE tauRHS
#endif /* USE_EXTEND_RHS */



END MODULE MOD_part_RHS
