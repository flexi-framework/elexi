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

!===================================================================================================================================
! Contains tools for particle related operations. This routine is used in MOD_Particle_Boundary_Tools, but not vice versa!
!===================================================================================================================================
MODULE MOD_Particle_Operations
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------!

PUBLIC:: CreateParticle
PUBLIC:: RemoveParticle
!===================================================================================================================================

CONTAINS

SUBROUTINE CreateParticle(SpeciesIn,PartStateIn,ElemID,PartID,LastPartPosIn,LastElemID,NewPartID)
!===================================================================================================================================
!> creates a single particle at correct array position and assign properties
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef,LastPartPos,PartSpecies
USE MOD_Particle_Vars          ,ONLY: doPartIndex,PartIndex,PartReflCount!,Species
USE MOD_Particle_Tools         ,ONLY: GetNextFreePosition
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: SpeciesIn
REAL, INTENT(IN)              :: PartStateIn(1:PP_nVarPart)
INTEGER, INTENT(IN)           :: ElemID
INTEGER, INTENT(IN)           :: PartID
REAL, INTENT(IN)              :: LastPartPosIn(1:3)
INTEGER, INTENT(IN)           :: LastElemID
INTEGER, INTENT(OUT),OPTIONAL :: NewPartID
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER :: newParticleID
!===================================================================================================================================
newParticleID = GetNextFreePosition()

PartSpecies(newParticleID)             = SpeciesIn
LastPartPos(1:3,newParticleID)         = LastPartPosIn(1:3)
PartState(1:PP_nVarPart,newParticleID) = PartStateIn(1:PP_nVarPart)
PartReflCount(newParticleID)           = 0

! Set the new reference position here
IF(TrackingMethod.EQ.REFMAPPING)THEN
  CALL GetPositionInRefElem(PartState(PART_POSV,newParticleID),PartPosRef(1:3,newParticleID),ElemID)
END IF ! TrackingMethod.EQ.REFMAPPING

PDM%ParticleInside(newParticleID) = .TRUE.
PDM%IsNewPart(newParticleID)      = .TRUE.
PEM%Element(newParticleID)        = ElemID
PEM%lastElement(newParticleID)    = LastElemID

IF (PRESENT(NewPartID)) NewPartID = newParticleID

!DO iInit = Species(SpeciesIn)%StartnumberOfInits, Species(SpeciesIn)%NumberOfInits
!  Species(SpeciesIn)%Init(iInit)%mySumOfMatchedParticles = Species(SpeciesIn)%Init(iInit)%mySumOfMatchedParticles + 1
!END DO

IF (doPartIndex) PartIndex(newParticleID) = PartIndex(PartID)

END SUBROUTINE CreateParticle


SUBROUTINE RemoveParticle(PartID,BCID,alpha,crossedBC)
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
INTEGER,INTENT(IN)           :: PartID
INTEGER,INTENT(IN) ,OPTIONAL :: BCID                    !< ID of the boundary the particle crossed
REAL   ,INTENT(OUT),OPTIONAL :: alpha                   !< if removed during tracking optional alpha can be set to -1
LOGICAL,INTENT(OUT),OPTIONAL :: crossedBC               !< optional flag is needed if particle removed on BC interaction
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iSpec
!===================================================================================================================================

PDM%ParticleInside(PartID) = .FALSE.

iSpec = PartSpecies(PartID)
! Count the number of particles per species and the kinetic energy per species
IF (CalcPartBalance .AND. PRESENT(BCID)) THEN
  nPartOut(iSpec)    = nPartOut(iSpec) + 1
  PartEkinOut(iSpec) = PartEkinOut(iSpec) + CalcEkinPart(PartID)
END IF ! CalcPartBalance

! Tracking-relevant variables (not required if a particle is removed within the domain, e.g. removal due to radial weighting)
IF (PRESENT(alpha))     alpha     = -1.
IF (PRESENT(crossedBC)) crossedBC = .TRUE.

END SUBROUTINE RemoveParticle

END MODULE MOD_Particle_Operations
