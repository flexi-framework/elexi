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

!==================================================================================================================================
!> Module for different Random Walk models of the particle discretization
!>
!> * Random Walk models require knowledge of turbulent kinetic energy and a turbulent length/time scale. Hide everything concerning
!>   RW models if these quantities are not available in FLEXI.
!>
!==================================================================================================================================
MODULE MOD_Particle_RandomWalk
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_RW
INTERFACE ParticleInitRandomWalk
  MODULE PROCEDURE ParticleInitRandomWalk
END INTERFACE

INTERFACE ParticleRandomWalk
  MODULE PROCEDURE ParticleRandomWalk
END INTERFACE

INTERFACE ParticleFinalizeRandomWalk
  MODULE PROCEDURE ParticleFinalizeRandomWalk
END INTERFACE

PUBLIC::ParticleInitRandomWalk
PUBLIC::ParticleRandomWalk
PUBLIC::ParticleFinalizeRandomWalk

!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleInitRandomWalk()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE ParticleInitRandomWalk

!===================================================================================================================================
! Wrapper for Random Walk Models
!===================================================================================================================================
SUBROUTINE ParticleRandomWalk(t)
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Vars,          ONLY : PDM
USE MOD_PICInterpolation_Vars,  ONLY : FieldAtParticle
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                  :: t               !> Current simulation time)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iPart
!===================================================================================================================================

! Iterate over all particles and add random walk to mean push
DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    CALL ParticleRandomWalkPush(iPart,t,FieldAtParticle(1:PP_nVar,iPart))
  END IF
END DO

END SUBROUTINE ParticleRandomWalk

!===================================================================================================================================
! Random Walk Push
!===================================================================================================================================
SUBROUTINE ParticleRandomWalkPush(PartID,t,FieldAtParticle)
! MODULES
USE MOD_Particle_Globals
USE MOD_EOS_Vars,          ONLY:mu0
USE MOD_Equation_Vars,     ONLY: betaStar
USE MOD_Particle_Vars,     ONLY: Species, PartSpecies, PartState, TurbPartState
USE MOD_PICInterpolation_Vars,ONLY: TurbFieldAtParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: t
REAL,INTENT(IN)     :: FieldAtParticle(1:PP_nVar)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: tke,epsturb,C_mu
REAL                :: l_e,t_e,t_r,t_int,tau_e,tau
REAL                :: udiff(1:3),udash(1:3)
REAL,PARAMETER      :: C_L = 0.3                                                                                     ! Debhi (2008)
REAL                :: ReP
REAL                :: Cd
REAL                :: lambda(3)
REAL                :: rho_p,Vol,r
REAL                :: nu
INTEGER             :: i
!===================================================================================================================================

! Check time stepping scheme requested for RW model
SELECT CASE(TRIM(Species(PartSpecies(PartID))%RWTime))

CASE('RK')
  ! Time stepping with every RK step. Do nothing here

CASE('RW')
  ! Time stepping with min(eddy time scale, transit time scale). Return if called too early
  IF (t.LT.TurbPartState(4,PartID)) RETURN

CASE DEFAULT
    CALL abort(__STAMP__, ' No particle random walk time step given. This should not happen.')

END SELECT


SELECT CASE(TRIM(Species(PartSpecies(PartID))%RWModel))

!===================================================================================================================================
! LES/DNS mode. Assume fully resolved turbulence
!===================================================================================================================================
CASE('none')
    ! Do nothing

CASE('Gosman')
    tke     = TurbFieldAtParticle(1,PartID)
    epsturb = TurbFieldAtParticle(2,PartID)

    ! SST renamed C_mu, change back for clarity
    C_mu    = betaStar

    ! Assume spherical particles for now
    Vol     = Species(PartSpecies(PartID))%MassIC/Species(PartSpecies(PartID))%DensityIC
    r       = (3.*Vol/4./pi)**(1./3.)
    rho_p   = Species(PartSpecies(PartID))%DensityIC

    ! Droplet relaxation time
    udiff(1:3) = PartState(4:6,PartID) - (FieldAtParticle(2:4)/FieldAtParticle(1))

    ! Get nu to stay in same equation format
    nu      = mu0/FieldAtParticle(1)

    Rep     = 2.*r*SQRT(udiff(1)**2. + udiff(2)**2. + udiff(3)**2.)/nu
    ! Empirical relation of nonlinear drag from Clift et al. (1978)
!    IF (Rep .LT. 1) THEN
!      Cd  = 1.
!    ELSE
      Cd  = 1. + 0.15*Rep**0.687
!    ENDIF

    ! Division by zero if particle velocity equal fluid velocity. Avoid this!
    IF (ANY(udiff.NE.0)) THEN
      tau     = (4./3.)*FieldAtParticle(1)*(2.*r)/(rho_p*Cd*sqrt(udiff(1)**2. + udiff(2)**2. + udiff(3)**2.))
    ELSE
      tau     = HUGE(1.)
    END IF

    ! Turbulent velocity fluctuation
    udash(1:3)   = (2.*tke/3.)**0.5

    ! Get random number with Gaussian distribution
    DO i = 1,3
      lambda(i) = RandNormal()
    END DO

    ! Catch negative or zero dissipation rate
    ! Eddy length and lifetime
    IF (TKE.GT.0 .AND. epsturb.GT.0) THEN
      l_e     = C_mu**(1./2.)*tke**(3./2.) / epsturb
!      tau_e   = C_L * tke / epsturb

      ! Calculate random walk push. We are isentropic, so use vectors
      TurbPartState(1:3,PartID) = lambda(1:3)*udash(1:3)

      ! Estimate interaction time. The values here are already depend on the random draw.
      ! If TKE!=0, then udash can only be zero if ALL Gaussian random numbers are zero.
      ! What are the chances? Ommit the check for now.
!      IF(ANY(udash.NE.0)) THEN
      t_e     = l_e/     SQRT(udash(1)**2. + udash(2)**2. + udash(3)**2.)
!       END IF
      tau_e   = l_e/(tau*SQRT(udiff(1)**2. + udiff(2)**2. + udiff(3)**2.))

      ! Particles either leaves the eddy or the eddy expires
      IF (tau_e.LT.1.) THEN
        t_r   = -tau*LOG(1.-tau_e)
        t_int = MIN(t_e,t_r)
      ! Particle captured by eddy, so interaction time is eddy life time
      ELSE
        t_int = t_e
      END IF

      ! Save time when eddy expires or particle finished crossing it
      TurbPartState(4,PartID)   = t + t_int

    ! TKE or epsilon is zero. No random velocity but try again in the next RK stage
    ELSE
      TurbPartState(1:3,PartID) = 0.
      TurbPartState(4,  PartID) = t
    ENDIF

CASE DEFAULT
    CALL abort(__STAMP__, ' No particle random walk model given. This should not happen.')

END SELECT


END SUBROUTINE ParticleRandomWalkPush

!===================================================================================================================================
!> Euler particle time integration:
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
SUBROUTINE ParticleFinalizeRandomWalk()
! MODULES

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE ParticleFinalizeRandomWalk

#endif /*USE_RW*/
END MODULE MOD_Particle_RandomWalk
