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
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools,                ONLY: GETSTR
USE MOD_Particle_Vars,              ONLY: PDM,TurbPartState
USE MOD_Particle_RandomWalk_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ALLOCSTAT
INTEGER               :: nRWVars
!----------------------------------------------------------------------------------------------------------------------------------
IF(ParticleRWInitIsDone) RETURN

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE RANDOM WALK ...!'

!--> Random Walk model
RWModel = TRIM(GETSTR('Part-RWModel','none'))
RWTime  = TRIM(GETSTR('Part-RWTime' ,'RW'))

SELECT CASE(RWModel)

  CASE('Mofakham')
    IF (RWTime.EQ.'RK'.AND.MPIROOT) &
      WRITE(*,*) 'WARNING: Modified Mofakham (2020) model is exactly Gosman (1983) with RW time step equal to RK time step.'
  nRWVars = 7

  CASE('Gosman')
    nRWVars = 4

  CASE DEFAULT
    CALL abort(__STAMP__, ' No particle random walk model given.')
END SELECT

! Allocate array to hold the RW properties for every particle
ALLOCATE(TurbPartState(1:nRWVars,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) &
  CALL abort(__STAMP__,'ERROR in particle_randomwalk.f90: Cannot allocate particle random walk arrays!')
TurbPartState = 0.

ParticleRWInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE RANDOM WALK DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE ParticleInitRandomWalk

!===================================================================================================================================
! Wrapper for Random Walk Models
!===================================================================================================================================
SUBROUTINE ParticleRandomWalk(t)
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Vars,                ONLY : PDM
USE MOD_Particle_Interpolation_Vars,  ONLY : FieldAtParticle
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
USE MOD_EOS_Vars,          ONLY: mu0
USE MOD_Equation_Vars,     ONLY: betaStar
USE MOD_Particle_Vars,     ONLY: Species,PartSpecies,PartState,TurbPartState
USE MOD_Particle_Interpolation_Vars, ONLY: TurbFieldAtParticle
USE MOD_Particle_RandomWalk_Vars,    ONLY: RWModel,RWTime
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

SELECT CASE(RWModel)

!===================================================================================================================================
! LES/DNS mode. Assume fully resolved turbulence
!===================================================================================================================================
CASE('none')
    ! Do nothing

CASE('Gosman')
!===================================================================================================================================
! Gosman (1983) model including the correction for exponent of C_mu equal 1/2.
!> Gosman, A. D., and Loannides, E., "Aspects of computer simulation of liquid-fueled combustors." Journal of energy 7.6 (1983):
!> 482-490.
!===================================================================================================================================
  ! Check time stepping scheme requested for RW model
  SELECT CASE(RWTime)

    CASE('RK')
      ! Time stepping with every RK step. Do nothing here

    CASE('RW')
      ! Time stepping with min(eddy time scale, transit time scale). Return if called too early
      !> Ideally, this should use tStage. But one cannot start a RK without the first stage and it does not make a difference for Euler
      IF (t.LT.TurbPartState(4,PartID)) RETURN

    CASE DEFAULT
      CALL abort(__STAMP__, ' No particle random walk time step given. This should not happen.')

  END SELECT

  tke     = TurbFieldAtParticle(1,PartID)
  epsturb = TurbFieldAtParticle(2,PartID)

  ! SST renamed C_mu, change back for clarity
  C_mu    = betaStar

  ! Assume spherical particles for now
  Vol     = Species(PartSpecies(PartID))%MassIC/Species(PartSpecies(PartID))%DensityIC
  r       = (3.*Vol/4./pi)**(1./3.)
  rho_p   = Species(PartSpecies(PartID))%DensityIC

  ! Turbulent velocity fluctuation
  udash(1:3)   = (2.*tke/3.)**0.5

  ! Get random number with Gaussian distribution
  DO i = 1,3
    lambda(i) = RandNormal()
  END DO

  ! Catch negative or zero dissipation rate
  ! Eddy length and lifetime
  IF (TKE.GT.0 .AND. epsturb.GT.0) THEN
    ! Calculate random walk push. We are isentropic, so use vectors
    TurbPartState(1:3,PartID) = lambda(1:3)*udash(1:3)

    ! Droplet relaxation time
    udiff(1:3) = PartState(4:6,PartID) - ((FieldAtParticle(2:4)/FieldAtParticle(1) + TurbPartState(1:3,PartID))

    ! Get nu to stay in same equation format
    nu      = mu0/FieldAtParticle(1)

    Rep     = 2.*r*SQRT(SUM(udiff(1:3)**2.))/nu
    ! Empirical relation of nonlinear drag from Clift et al. (1978)
  !    IF (Rep .LT. 1) THEN
  !      Cd  = 1.
  !    ELSE
      Cd  = 1. + 0.15*Rep**0.687
  !    ENDIF

    ! Division by zero if particle velocity equal fluid velocity. Avoid this!
    IF (ANY(udiff.NE.0)) THEN
      tau     = (4./3.)*rho_p*(2.*r)/(FieldAtParticle(1)*Cd*SQRT(SUM(udiff(1:3)**2.)))
    ELSE
      tau     = HUGE(1.)
    END IF

    l_e     = C_mu**(1./2.)*tke**(3./2.) / epsturb
!    tau_e   = C_L * tke / epsturb

    ! Estimate interaction time. The values here are already depend on the random draw.
    ! If TKE!=0, then udash can only be zero if ALL Gaussian random numbers are zero.
    ! What are the chances? Omit the check for now.
!    IF(ANY(udash.NE.0)) THEN
    t_e     = l_e/     SQRT(SUM(udash(1:3)**2.))
!    END IF
    tau_e   = l_e/(tau*SQRT(SUM(udiff(1:3)**2.)))

    ! Particles either leaves the eddy or the eddy expires
    IF (tau_e.LT.1.) THEN
      t_r   = -tau*LOG(1.-tau_e)
      t_int = MIN(t_e,t_r)
    ! Particle captured by eddy, so interaction time is eddy life time
    ELSE
      t_int = t_e
    END IF

    ! Save time when eddy expires or particle finished crossing it
    !> Ideally, this should use tStage. But one cannot start a RK without the first stage and it does not make a difference for
    !> Euler. Plus, that way we make the same error twice (but in opposite direction in time), so we're regaining some accuracy?
    TurbPartState(4,PartID)   = t + t_int

  ! TKE or epsilon is zero. No random velocity but try again in the next RK stage
  ELSE
    TurbPartState(1:3,PartID) = 0.
    TurbPartState(4,  PartID) = t
  ENDIF

CASE('Mofakham')
!===================================================================================================================================
! Modified Mofakham (2020) model. Calculation of udash and t_int similar to Gosman (1983) but random variables lambda_i fixed for
! t_int. Basis of this reasoning is following passage:
!<<< "It should be emphasized that the cross-correlation between the components of velocity fluctuations in different directions in
!<<<  inhomogeneous flows cannot be included by employing the standard DRW model. In this approach, it is assumed that a particle
!<<<  interacts with a turbulent eddy for a time interval t_int during which [lambda_i] is fixed."
!> Mofakham, Amir A., and Ahmadi, G., "On random walk models for simulation of particle-laden turbulent flows.", International
!> Journal of Multiphase Flow (2019): 103157.
!===================================================================================================================================
  ! Particle still inside of eddy, only update the RMS turbulence
  !> Ideally, this should use tStage. But one cannot start a RK without the first stage and it does not make a difference for Euler
  IF (t.LT.TurbPartState(4,PartID)) THEN
    tke         = TurbFieldAtParticle(1,PartID)
    lambda(1:3) = TurbPartState    (5:7,PartID)

    ! Turbulent velocity fluctuation
    udash(1:3)  = (2.*tke/3.)**0.5

    ! Calculate random walk push. We are isentropic, so use vectors
    TurbPartState(1:3,PartID) = lambda(1:3)*udash(1:3)
  ELSE
    tke     = TurbFieldAtParticle(1,PartID)
    epsturb = TurbFieldAtParticle(2,PartID)

    ! SST renamed C_mu, change back for clarity
    C_mu    = betaStar

    ! Assume spherical particles for now
    Vol     = Species(PartSpecies(PartID))%MassIC/Species(PartSpecies(PartID))%DensityIC
    r       = (3.*Vol/4./pi)**(1./3.)
    rho_p   = Species(PartSpecies(PartID))%DensityIC

    ! Turbulent velocity fluctuation
    udash(1:3)   = (2.*tke/3.)**0.5

    ! Get random number with Gaussian distribution
    DO i = 1,3
      lambda(i) = RandNormal()
    END DO

    ! Catch negative or zero dissipation rate
    ! Eddy length and lifetime
    IF (TKE.GT.0 .AND. epsturb.GT.0) THEN
      ! Calculate random walk push. We are isentropic, so use vectors
      TurbPartState(1:3,PartID) = lambda(1:3)*udash(1:3)

      ! Droplet relaxation time
      udiff(1:3) = PartState(4:6,PartID) - (FieldAtParticle(2:4)/FieldAtParticle(1) + TurbPartState(1:3,PartID))

      ! Get nu to stay in same equation format
      nu      = mu0/FieldAtParticle(1)

      Rep     = 2.*r*SQRT(SUM(udiff(1:3)**2.))/nu
      ! Empirical relation of nonlinear drag from Clift et al. (1978)
!    IF (Rep .LT. 1) THEN
!      Cd  = 1.
!    ELSE
      Cd  = 1. + 0.15*Rep**0.687
!    ENDIF

      ! Division by zero if particle velocity equal fluid velocity. Avoid this!
      IF (ANY(udiff.NE.0)) THEN
        tau     = (4./3.)*rho_p*(2.*r)/(FieldAtParticle(1)*Cd*SQRT(SUM(udiff(1:3)**2.)))
      ELSE
        tau     = HUGE(1.)
      END IF

      l_e     = C_mu**(1./2.)*tke**(3./2.) / epsturb
!      tau_e   = C_L * tke / epsturb

      ! Estimate interaction time. The values here are already depend on the random draw.
      ! If TKE!=0, then udash can only be zero if ALL Gaussian random numbers are zero.
      ! What are the chances? Omit the check for now.
!      IF(ANY(udash.NE.0)) THEN
      t_e     = l_e/     SQRT(SUM(udash(1:3)**2.))
!       END IF
      tau_e   = l_e/(tau*SQRT(SUM(udiff(1:3)**2.)))

      ! Particles either leaves the eddy or the eddy expires
      IF (tau_e.LT.1.) THEN
        t_r   = -tau*LOG(1.-tau_e)
        t_int = MIN(t_e,t_r)
      ! Particle captured by eddy, so interaction time is eddy life time
      ELSE
        t_int = t_e
      END IF

      ! Save time when eddy expires or particle finished crossing it and the random numbers
      TurbPartState(4  ,PartID) = t + t_int
      TurbPartState(5:7,PartID) = lambda(1:3)

    ! TKE or epsilon is zero. No random velocity but try again in the next RK stage
    ELSE
      TurbPartState(1:3,PartID) = 0.
      TurbPartState(4  ,PartID) = t
      ! Nullify random numbers just to make sure
      TurbPartState(5:7,PartID) = 0.
    ENDIF
  END IF ! t.LT.t_int

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
USE MOD_Particle_Vars,              ONLY: TurbPartState
USE MOD_Particle_RandomWalk_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(TurbPartState)

ParticleRWInitIsDone=.FALSE.

END SUBROUTINE ParticleFinalizeRandomWalk

#endif /*USE_RW*/
END MODULE MOD_Particle_RandomWalk
