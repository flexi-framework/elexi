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
#if EQNSYSNR == 4
INTERFACE Particle_InitRandomWalk
  MODULE PROCEDURE Particle_InitRandomWalk
END INTERFACE

INTERFACE Particle_RandomWalk
  MODULE PROCEDURE Particle_RandomWalk
END INTERFACE

INTERFACE Particle_FinalizeRandomWalk
  MODULE PROCEDURE Particle_FinalizeRandomWalk
END INTERFACE

PUBLIC::Particle_InitRandomWalk
PUBLIC::Particle_RandomWalk
PUBLIC::Particle_FinalizeRandomWalk

!===================================================================================================================================

CONTAINS

SUBROUTINE Particle_InitRandomWalk()
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

END SUBROUTINE Particle_InitRandomWalk

!===================================================================================================================================
! Wrapper for Random Walk Models
!===================================================================================================================================
SUBROUTINE Particle_RandomWalk()
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Vars,          ONLY : PDM, Pt
USE MOD_PICInterpolation_Vars,  ONLY : FieldAtParticle
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iPart
!===================================================================================================================================

! Iterate over all particles and add random walk to mean push
DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    Pt(iPart,1:3) = Pt(iPart,1:3) + Particle_RandomWalkPush(iPart,FieldAtParticle(iPart,1:PP_nVar))
  END IF
END DO

END SUBROUTINE Particle_RandomWalk

!===================================================================================================================================
! Random Walk Push
!===================================================================================================================================
FUNCTION Particle_RandomWalkPush(PartID,FieldAtParticle)
! MODULES
USE MOD_Particle_Globals
USE MOD_Particle_Vars,     ONLY : Species, PartSpecies!, PartGravity
!USE MOD_Particle_Vars,     ONLY : PartState, RepWarn
!USE MOD_EOS_Vars,          ONLY : mu0
USE MOD_Equation_Vars,     ONLY : betaStar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(IN)     :: FieldAtParticle(1:PP_nVar)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: Particle_RandomWalkPush(1:3)             ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: tke,omega,epsturb,C_mu
REAL                :: L_e,tau_e
REAL                :: udash,vdash,wdash
REAL,PARAMETER      :: C_L = 0.3                                ! Debhi (2008)
REAL                :: random
REAL                :: lambda(3)
!===================================================================================================================================

SELECT CASE(TRIM(Species(PartSpecies(PartID))%RWModel))

!===================================================================================================================================
! LES/DNS mode. Assume fully resolved turbulence
!===================================================================================================================================
CASE('none')
    ! Do nothing
    
CASE('Gosman')
    tke     = FieldAtParticle(6)
    omega   = FieldAtParticle(7)
    epsturb = omega*tke*betaStar
    
    ! SST renamed C_mu, change back for clarity
    C_mu    = betaStar
    
    ! Eddy length and lifetime
    L_e     = C_mu**(3./4.)*tke**(3./2.) / epsturb
    tau_e   = C_L * tke / epsturb
    
    ! Turbulent velocity fluctuation
    udash   = (2.*tke/3.)**0.5
    vdash   = udash
    wdash   = udash
    
    ! Get random number and rescale to [-1,1]
    CALL RANDOM_NUMBER(random)
    lambda(1) = 2.*random - 1.
    CALL RANDOM_NUMBER(random)
    lambda(2) = 2.*random - 1.
    CALL RANDOM_NUMBER(random)
    lambda(3) = 2.*random - 1.
    
    ! Calculate random walk push. We are isentropic, so use vectors
    Particle_RandomWalkPush(:) = lambda(:)*udash
    
CASE DEFAULT
    CALL abort(__STAMP__, ' No particle random walk model given. This should not happen.')
    
END SELECT


END FUNCTION Particle_RandomWalkPush

!===================================================================================================================================
!> Euler particle time integration:
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
SUBROUTINE Particle_FinalizeRandomWalk()
! MODULES

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE Particle_FinalizeRandomWalk
#endif

END MODULE MOD_Particle_RandomWalk
