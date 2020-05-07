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

MODULE MOD_part_operations
!===================================================================================================================================
! Contains tools for particle related operations. This routine is used in MOD_Particle_Boundary_Tools, but not vice versa!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE CreateParticle
  MODULE PROCEDURE CreateParticle
END INTERFACE

INTERFACE RemoveParticle
  MODULE PROCEDURE RemoveParticle
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: CreateParticle, RemoveParticle
!===================================================================================================================================

CONTAINS

SUBROUTINE CreateParticle(Species,Pos,ElemID,Velocity,NewPartID)
!===================================================================================================================================
!> creates a single particle at correct array position and assign properties
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: PDM,PEM,PartState,LastPartPos,PartSpecies
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: Species
REAL, INTENT(IN)              :: Pos(1:3)
INTEGER, INTENT(IN)           :: ElemID
REAL, INTENT(IN)              :: Velocity(1:3)
INTEGER, INTENT(OUT),OPTIONAL :: NewPartID
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER :: newParticleID
!===================================================================================================================================
PDM%ParticleVecLength = PDM%ParticleVecLength + 1 ! Increase particle vector length
newParticleID         = PDM%ParticleVecLength
IF (newParticleID.GT.PDM%MaxParticleNumber) &
  CALL ABORT(__STAMP__,'CreateParticle: newParticleID.GT.PDM%MaxParticleNumber. newParticleID=',newParticleID)

PartSpecies(newParticleID)     = Species
LastPartPos(1:3,newParticleID) = Pos(1:3)
PartState(1:3,newParticleID)   = Pos(1:3)
PartState(4:6,newParticleID)   = Velocity(1:3)

PDM%ParticleInside(newParticleID) = .TRUE.
PDM%IsNewPart(newParticleID)      = .FALSE.   ! ??????? correct ????
PEM%Element(newParticleID)        = ElemID
PEM%lastElement(newParticleID)    = ElemID

IF (PRESENT(NewPartID)) NewPartID = newParticleID

END SUBROUTINE CreateParticle


SUBROUTINE RemoveParticle(PartID,alpha,crossedBC)
!===================================================================================================================================
!> Removes a single particle "PartID" by setting the required variables.
!> If CalcPartBalance = T: adds/substracts the particle to/from the respective counter
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars           ,ONLY: PDM,PartSpecies
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance,nPartOut,PartEkinOut
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: PartID
REAL, INTENT(OUT),OPTIONAL    :: alpha                   !< if removed during tracking optional alpha can be set to -1
LOGICAL, INTENT(OUT),OPTIONAL :: crossedBC               !< optional flag is needed if particle removed on BC interaction
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iSpec
!===================================================================================================================================

PDM%ParticleInside(PartID) = .FALSE.

iSpec = PartSpecies(PartID)
! Count the number of particles per species and the kinetic energy per species
IF (CalcPartBalance) THEN
  nPartOut(iSpec)    = nPartOut(iSpec) + 1
  PartEkinOut(iSpec) = PartEkinOut(iSpec) + CalcEkinPart(PartID)
END IF ! CalcPartBalance

! Tracking-relevant variables (not required if a particle is removed within the domain, e.g. removal due to radial weighting)
IF (PRESENT(alpha))     alpha     = -1.
IF (PRESENT(crossedBC)) crossedBC = .TRUE.

END SUBROUTINE RemoveParticle

END MODULE MOD_part_operations
