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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_Tracking
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE PerformTracking
  MODULE PROCEDURE PerformTracking
END INTERFACE

INTERFACE ParticleInsideCheck
  MODULE PROCEDURE ParticleInsideCheck
END INTERFACE

PUBLIC::PerformTracking
PUBLIC::ParticleInsideCheck
!===================================================================================================================================

CONTAINS

SUBROUTINE PerformTracking()
!===================================================================================================================================
!> Routine called from the timedisc to call the selected tracking routine
!===================================================================================================================================
! MODULES
USE MOD_Globals                  ,ONLY: Abort
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Tracing         ,ONLY: ParticleTracing
USE MOD_Particle_RefTracking     ,ONLY: ParticleRefTracking
USE MOD_Particle_TriaTracking    ,ONLY: ParticleTriaTracking
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SELECT CASE(TrackingMethod)
  CASE(REFMAPPING)
    CALL ParticleRefTracking()
  CASE(TRACING)
    CALL ParticleTracing()
  CASE(TRIATRACKING)
    CALL ParticleTriaTracking()
  CASE DEFAULT
    CALL Abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfo=TrackingMethod)
END SELECT

END SUBROUTINE PerformTracking


LOGICAL FUNCTION ParticleInsideCheck(Position,iPart,GlobalElemID)
!===================================================================================================================================
!> Checks if the position is inside the element with the appropriate routine depending on the TrackingMethod
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_Particle_Localization   ,ONLY: PartInElemCheck,ParticleInsideQuad3D
USE MOD_Particle_Mesh_Tools     ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Vars           ,ONLY: PartPosRef
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: Position(3)
INTEGER, INTENT(IN)             :: iPart,GlobalElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

ParticleInsideCheck = .FALSE.

SELECT CASE(TrackingMethod)
  CASE(REFMAPPING)
    CALL GetPositionInRefElem(Position,PartPosRef(1:3,iPart),GlobalElemID)
    IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) ParticleInsideCheck=.TRUE.
  CASE(TRACING)
    CALL PartInElemCheck(     Position,iPart,GlobalElemID,ParticleInsideCheck)
  CASE(TRIATRACKING)
    CALL ParticleInsideQuad3D(Position,GlobalElemID,ParticleInsideCheck)
  CASE DEFAULT
    CALL Abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfo=TrackingMethod)
END SELECT

END FUNCTION ParticleInsideCheck

END MODULE MOD_Particle_Tracking
