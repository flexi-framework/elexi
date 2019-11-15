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

!===================================================================================================================================
! Subroutines to interpolate the DG solution to the particle position
!===================================================================================================================================
MODULE  MOD_PICInterpolation
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: InterpolateFieldToParticle,InitializeInterpolation,InterpolateFieldToSingleParticle
!===================================================================================================================================
INTERFACE InitializeInterpolation
  MODULE PROCEDURE InitializeInterpolation
END INTERFACE

INTERFACE InterpolateFieldToParticle
  MODULE PROCEDURE InterpolateFieldToParticle
END INTERFACE

INTERFACE InterpolateFieldToSingleParticle
  MODULE PROCEDURE InterpolateFieldToSingleParticle
END INTERFACE
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeInterpolation
!===================================================================================================================================
! Initialize the interpolation variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals,       ONLY:PP_nElems
USE MOD_ReadInTools
USE MOD_Particle_Vars,          ONLY:PDM
USE MOD_PICInterpolation_Vars
#if USE_RW
USE MOD_Equation_Vars,          ONLY:nVarTurb
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ALLOCSTAT
REAL                      :: scaleExternalField
!===================================================================================================================================

DoInterpolation       = GETLOGICAL('PIC-DoInterpolation','.TRUE.')

! For low number of particles, the loop over all elements becomes quite inefficient. User can opt out with setting
! InterpolationElemLoop = F.
InterpolationElemLoop = GETLOGICAL('PIC-InterpolationElemLoop','.TRUE.')

! Even if the user did not opt out, switch InterpolationElemLoop off for procs with high number of elems.
!>> so far arbitrary treshold of 10 elems per proc
IF (InterpolationElemLoop.AND.(PP_nElems.GT.10)) InterpolationElemLoop=.FALSE.

! Size of external field depends on the equations system. Distinguish here between LES and RANS-SA
#if EQNSYSNR == 3   /*Spalart-Allmaras*/
externalField(1:PP_nVar)= GETREALARRAY('PIC-externalField',PP_nVar,'0.,0.,0.,0.,0.,0.')
#else               /*Navier-Stokes*/
externalField(1:PP_nVar)= GETREALARRAY('PIC-externalField',PP_nVar,'0.,0.,0.,0.,0.')
#endif

! Only consider the external field if it is not equal to zero
IF (ANY(externalField.NE.0)) THEN
  useExternalField   = .TRUE.
  scaleexternalField = GETREAL('PIC-scaleexternalField','1.0')
  externalField      = externalField*ScaleExternalField
ELSE
  useExternalField   = .FALSE.
END IF

!--- Allocate arrays for interpolation of fields to particles
SDEALLOCATE(FieldAtParticle)
! Allocate array for rho,u_x,u_y,u_z,e
ALLOCATE(FieldAtParticle    (1:PP_nVar, 1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
#if USE_RW
SDEALLOCATE(TurbFieldAtParticle)
! Allocate array for TKE,epsilson
ALLOCATE(TurbFieldAtParticle(1:nVarTurb,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
#endif
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in pic_interpolation.f90: Cannot allocate FieldAtParticle array!',ALLOCSTAT)

END SUBROUTINE InitializeInterpolation


SUBROUTINE InterpolateFieldToParticle()
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Eval_xyz,                ONLY: EvaluateFieldAtPhysPos,EvaluateFieldAtRefPos
USE MOD_Mesh_Vars,               ONLY: nElems
USE MOD_Particle_Vars,           ONLY: PartPosRef,PartState,PDM,PEM
USE MOD_Particle_Tracking_Vars,  ONLY: DoRefMapping
USE MOD_PICInterpolation_Vars,   ONLY: FieldAtParticle,useExternalField,externalField
USE MOD_PICInterpolation_Vars,   ONLY: DoInterpolation,InterpolationElemLoop
#if USE_RW
USE MOD_DG_Vars,                 ONLY: UTurb
USE MOD_Equation_Vars,           ONLY: nVarTurb
USE MOD_PICInterpolation_Vars,   ONLY: TurbFieldAtParticle
USE MOD_Restart_Vars,            ONLY: RestartTurb
#endif
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: firstPart,lastPart
REAL                             :: field(PP_nVar)
INTEGER                          :: iPart,iElem,iVar
#if USE_RW
REAL                             :: turbField(nVarTurb)
#endif
!===================================================================================================================================

! Return if no interpolation is wanted
IF (.NOT.DoInterpolation) RETURN

! null field vector
field     = 0.
#if USE_RW
! null turbulent field vector
turbField = 0.
#endif

firstPart = 1
lastPart  = PDM%ParticleVecLength

! Return if there are no particles
IF(firstPart.GT.lastPart) RETURN

! For low number of particles, the loop over all elements becomes quite inefficient. User can opt out with setting
! InterpolationElemLoop = F.
IF (.NOT.InterpolationElemLoop) THEN
  DO iPart = firstPart, LastPart
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:PP_nVar,iPart))
  END DO
  RETURN
END IF

! InterpolationElemLoop is true, so initialize everything for it
FieldAtParticle(:,firstPart:lastPart)   = 0.
IF (useExternalField) THEN
  DO iVar = 1,PP_nVar
    FieldAtParticle(iVar,firstPart:lastPart) = externalField(iVar)
  END DO
END IF

! Loop first over all elements, then over all particles within the element. Ideally, this reduces cache misses as the interpolation
! happens with the same element metrics
DO iElem=1,nElems
  DO iPart=firstPart,LastPart
    ! Particle already left the domain, ignore it
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

    ! Particle is inside and in current element
    IF (PEM%Element(iPart).EQ.iElem) THEN
      ! Not RefMapping, evaluate at physical position
      IF (.NOT.DoRefMapping) THEN
#if USE_RW
        IF (RestartTurb) THEN
          CALL EvaluateFieldAtPhysPos(PartState(1:3,iPart),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,iElem),field    (1:PP_nVar),iElem,iPart &
                                                                       ,UTurb(1:nVarTurb,:,:,:,iElem),turbField(1:nVarTurb))
        ELSE
#endif
          CALL EvaluateFieldAtPhysPos(PartState(1:3,iPart),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,iElem),field    (1:PP_nVar),iElem,iPart)
#if USE_RW
       END IF
#endif
      ! RefMapping, evaluate in reference space
      ELSE
#if USE_RW
        IF (RestartTurb) THEN
          CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,iElem),field    (1:PP_nVar),iElem &
                                                                       ,UTurb(1:nVarTurb,:,:,:,iElem),turbField(1:nVarTurb))
        ELSE
#endif
          CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,iElem),field    (1:PP_nVar),iElem)
#if USE_RW
        END IF
#endif
      END IF ! RefMapping

      ! Add the interpolated field to the background field
      FieldAtParticle(    1:PP_nVar, iPart) = FieldAtParticle(1:PP_nVar,iPart) + field(1:PP_nVar)
#if USE_RW
      TurbFieldAtParticle(1:nVarTurb,iPart) = turbfield(      1:nVarTurb)
#endif
    END IF ! Element(iPart).EQ.iElem
  END DO ! iPart
END DO ! iElem=1,PP_N

END SUBROUTINE InterpolateFieldToParticle


SUBROUTINE InterpolateFieldToSingleParticle(PartID,FieldAtParticle)
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Eval_xyz,                ONLY: EvaluateFieldAtPhysPos,EvaluateFieldAtRefPos
USE MOD_PICInterpolation_Vars,   ONLY: useExternalField,externalField
USE MOD_Particle_Tracking_Vars,  ONLY: DoRefMapping
USE MOD_Particle_Vars,           ONLY: PartPosRef,PartState,PEM
#if USE_MPI
USE MOD_Mesh_Vars,               ONLY: nElems
#endif
#if USE_RW
USE MOD_DG_Vars,                 ONLY: UTurb
USE MOD_Restart_Vars,            ONLY: RestartTurb
USE MOD_Equation_Vars,           ONLY: nVarTurb
USE MOD_PICInterpolation_Vars,   ONLY: TurbFieldAtParticle
#endif
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: FieldAtParticle(1:PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: field(1:PP_nVar)
INTEGER                          :: ElemID
#if USE_RW
REAL                             :: turbField(nVarTurb)
#endif
!===================================================================================================================================

FieldAtParticle(:)   = 0.
IF (useExternalField) THEN
  FieldAtParticle(:) = externalField(:)
END IF

ElemID=PEM%Element(PartID)
! No solution available in the halo region (yet), so return for particles there
#if USE_MPI
IF(ElemID.GT.nElems) RETURN
#endif

IF (.NOT.DoRefMapping) THEN
#if USE_RW
  IF (RestartTurb) THEN
    CALL EvaluateFieldAtPhysPos(PartState(1:3,PartID),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,ElemID),field    (1:PP_nVar),ElemID,PartID &
                                                                 ,UTurb(1:nVarTurb,:,:,:,ElemID),turbField(1:nVarTurb))
  ELSE
#endif
    CALL EvaluateFieldAtPhysPos(PartState(1:3,PartID),PP_nVar,PP_N,U    (1:PP_nVar,:,:,:,ElemID),field     (1:PP_nVar),ElemID,PartID)
#if USE_RW
  END IF
#endif
! RefMapping, evaluate in reference space
ELSE
#if USE_RW
  IF (RestartTurb) THEN
   CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,ElemID),field    (1:PP_nVar),ElemID        &
                                                                ,UTurb(1:nVarTurb,:,:,:,ElemID),turbField(1:nVarTurb))
  ELSE
#endif
   CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),PP_nVar,PP_N,U    (1:PP_nVar,:,:,:,ElemID) ,field    (1:PP_nVar),ElemID)
#if USE_RW
  END IF
#endif
END IF ! RefMapping

! Add the interpolated field to the background field
FieldAtParticle(    1:PP_nVar )        = FieldAtParticle(1:PP_nVar)  + field(1:PP_nVar)
#if USE_RW
TurbFieldAtParticle(1:nVarTurb,PartID) = turbfield(      1:nVarTurb)
#endif

END SUBROUTINE InterpolateFieldToSingleParticle


END MODULE MOD_PICInterpolation
