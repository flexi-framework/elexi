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
! Contains routines for interpolation of one-way coupled particles
!===================================================================================================================================
MODULE MOD_Particle_Interpolation
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleInterpolation
  MODULE PROCEDURE InitParticleInterpolation
END INTERFACE

INTERFACE DefineParametersParticleInterpolation
  MODULE PROCEDURE DefineParametersParticleInterpolation
END INTERFACE

INTERFACE InterpolateFieldToParticle
  MODULE PROCEDURE InterpolateFieldToParticle
END INTERFACE

INTERFACE InterpolateFieldToSingleParticle
  MODULE PROCEDURE InterpolateFieldToSingleParticle
END INTERFACE

PUBLIC :: InitParticleInterpolation
PUBLIC :: DefineParametersParticleInterpolation
PUBLIC :: InterpolateFieldToParticle
PUBLIC :: InterpolateFieldToSingleParticle
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle interpolation
!==================================================================================================================================
SUBROUTINE DefineParametersParticleInterpolation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Interpolation")

CALL prms%CreateLogicalOption(  'Part-DoInterpolation'         , "Compute the DG solution field's influence on the Particle", '.TRUE.')
CALL prms%CreateLogicalOption(  'Part-InterpolationElemLoop'   , 'Interpolate with outer iElem-loop (not'                         //&
                                                                'for many Elems per proc!)', '.TRUE.')
CALL prms%CreateRealArrayOption('Part-externalField'           , 'External field is added to the'                                 //&
                                                                'maxwell-solver-field', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'Part-scaleexternalField'      , 'Scale the provided external field', '1.0')

END SUBROUTINE DefineParametersParticleInterpolation

SUBROUTINE InitParticleInterpolation()
!===================================================================================================================================
! Initialize the interpolation variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals,       ONLY:PP_nElems
USE MOD_ReadInTools
USE MOD_Particle_Vars,          ONLY:PDM
USE MOD_Particle_Interpolation_Vars
#if USE_RW
USE MOD_Equation_Vars,          ONLY:nVarTurb
#endif
#if USE_MPI
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ALLOCSTAT,LoopDisabled,LoopDisabledGlob
REAL                      :: scaleExternalField
!===================================================================================================================================
IF(PartInterpolationInitIsDone)THEN
   SWRITE(*,*) "InitParticleInterpolation already called."
   RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE INTERPOLATION ...'

DoInterpolation       = GETLOGICAL('Part-DoInterpolation','.TRUE.')

! For low number of particles, the loop over all elements becomes quite inefficient. User can opt out with setting
! InterpolationElemLoop = F.
InterpolationElemLoop = GETLOGICAL('Part-InterpolationElemLoop','.TRUE.')

! Even if the user did not opt out, switch InterpolationElemLoop off for procs with high number of elems.
!>> so far arbitrary treshold of 10 elems per proc
IF (InterpolationElemLoop.AND.(PP_nElems.GT.10)) THEN
  InterpolationElemLoop = .FALSE.
  LoopDisabled          = 1

#if USE_MPI
  CALL MPI_ALLREDUCE(LoopDisabled,LoopDisabledGlob,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError)
  SWRITE(*,*) ' InterpolationElemLoop disabled due to high number of elements on ',LoopDisabledGlob,'of',PartMPI%nProcs,'procs'
#else
  LoopDisabledGlob = LoopDisabled
  WRITE(*,*)  ' InterpolationElemLoop disabled due to high number of elements'
#endif
END IF

! Size of external field depends on the equations system. Distinguish here between LES and RANS-SA
#if EQNSYSNR == 3   /*Spalart-Allmaras*/
externalField(1:PP_nVar)= GETREALARRAY('Part-externalField',PP_nVar,'0.,0.,0.,0.,0.,0.')
#else               /*Navier-Stokes*/
externalField(1:PP_nVar)= GETREALARRAY('Part-externalField',PP_nVar,'0.,0.,0.,0.,0.')
#endif

! Only consider the external field if it is not equal to zero
IF (ANY(externalField.NE.0)) THEN
  useExternalField   = .TRUE.
  scaleexternalField = GETREAL('Part-scaleexternalField','1.0')
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
IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in Part_interpolation.f90: Cannot allocate FieldAtParticle array!',ALLOCSTAT)

PartInterpolationInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE INTERPOLATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleInterpolation


SUBROUTINE InterpolateFieldToParticle(nVar,U,FieldAtParticle)
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Eval_xyz,                    ONLY: EvaluateFieldAtPhysPos,EvaluateFieldAtRefPos
USE MOD_Mesh_Vars,                   ONLY: nElems
USE MOD_Particle_Interpolation_Vars, ONLY: useExternalField,externalField
USE MOD_Particle_Interpolation_Vars, ONLY: DoInterpolation,InterpolationElemLoop
USE MOD_Particle_Tracking_Vars,      ONLY: DoRefMapping
USE MOD_Particle_Vars,               ONLY: PartPosRef,PartState,PDM,PEM
USE MOD_Particle_Vars,               ONLY: TurbPartState
#if USE_RW
USE MOD_Particle_RandomWalk_Vars,    ONLY: RWModel,RWTime
USE MOD_TimeDisc_Vars,               ONLY: t
#endif
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: nVar
REAL,INTENT(IN)                     :: U(1:nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                    :: FieldAtParticle(1:nVar,1:PDM%maxParticleNumber)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: firstPart,lastPart
REAL                             :: field(nVar)
INTEGER                          :: iPart,iElem,iVar
!===================================================================================================================================

! Return if no interpolation is wanted
IF (.NOT.DoInterpolation) RETURN

! Null field vector
field     = 0.

! Loop variables
firstPart = 1
lastPart  = PDM%ParticleVecLength

! Return if there are no particles
IF(firstPart.GT.lastPart) RETURN

! For low number of particles, the loop over all elements becomes quite inefficient. User can opt out with setting
! InterpolationElemLoop = F.
IF (.NOT.InterpolationElemLoop) THEN
  DO iPart = firstPart, LastPart
    ! Particle already left the domain, ignore it
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
#if USE_RW
    ! Do not change the particle velocity if RW is working in full Euler mode
    !> Ideally, this should use tStage. But one cannot start a RK without the first stage and it does not make a difference for Euler
    IF ((RWModel.EQ.'Gosman') .AND. (RWTime.EQ.'RW') .AND. (t.LT.TurbPartState(4,iPart))) CYCLE
#endif
    CALL InterpolateFieldToSingleParticle(iPart,nVar,U(:,:,:,:,PEM%Element(iPart)),FieldAtParticle(:,iPart))
  END DO
  RETURN
END IF

! InterpolationElemLoop is true, so initialize everything for it
FieldAtParticle(:,firstPart:lastPart)   = 0.
!IF (useExternalField) THEN
!  DO iVar = 1,PP_nVar
!    FieldAtParticle(iVar,firstPart:lastPart) = externalField(iVar)
!  END DO
!END IF

! Loop first over all elements, then over all particles within the element. Ideally, this reduces cache misses as the interpolation
! happens with the same element metrics
DO iElem=1,nElems
  DO iPart=firstPart,LastPart
    ! Particle already left the domain, ignore it
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
#if USE_RW
    ! Do not change the particle velocity if RW is working in full Euler mode
     !> Ideally, this should use tStage. But one cannot start a RK without the first stage and it does not make a difference for Euler
    IF ((RWTime.EQ.'RW') .AND. (t.LT.TurbPartState(4,iPart))) CYCLE
#endif

    ! Particle is inside and in current element
    IF (PEM%Element(iPart).EQ.iElem) THEN
      ! Not RefMapping, evaluate at physical position
      IF (.NOT.DoRefMapping) THEN
        CALL EvaluateFieldAtPhysPos(PartState(1:3,iPart),nVar,PP_N,U    (:,:,:,:,iElem),field,iElem,iPart)
      ! RefMapping, evaluate in reference space
      ELSE
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),nVar,PP_N,U    (:,:,:,:,iElem),field,iElem)
      END IF ! RefMapping

      ! Add the interpolated field to the background field
      FieldAtParticle(:,iPart) = FieldAtParticle(:,iPart) + field(:)
    END IF ! Element(iPart).EQ.iElem
  END DO ! iPart
END DO ! iElem=1,PP_N

END SUBROUTINE InterpolateFieldToParticle


SUBROUTINE InterpolateFieldToSingleParticle(PartID,nVar,U,FieldAtParticle)
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
!USE MOD_DG_Vars,                 ONLY: U
USE MOD_Eval_xyz,                ONLY: EvaluateFieldAtPhysPos,EvaluateFieldAtRefPos
!USE MOD_Particle_Interpolation_Vars,   ONLY: useExternalField,externalField
USE MOD_Particle_Tracking_Vars,  ONLY: DoRefMapping
USE MOD_Particle_Vars,           ONLY: PartPosRef,PartState,PEM
#if USE_MPI
USE MOD_Mesh_Vars,               ONLY: nElems
#endif
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
INTEGER,INTENT(IN)               :: nVar
REAL,INTENT(IN)                  :: U(1:nVar,0:PP_N,0:PP_N,0:PP_NZ)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: FieldAtParticle(1:nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: field(1:nVar)
INTEGER                          :: ElemID
!===================================================================================================================================

FieldAtParticle(:)   = 0.
!IF (useExternalField) THEN
!  FieldAtParticle(:) = externalField(:)
!END IF

ElemID=PEM%Element(PartID)
! No solution available in the halo region (yet), so return for particles there
#if USE_MPI
IF(ElemID.GT.nElems) RETURN
#endif

IF (.NOT.DoRefMapping) THEN
  CALL EvaluateFieldAtPhysPos(PartState(1:3,PartID),nVar,PP_N,U    (:,:,:,:),field,ElemID,PartID)
! RefMapping, evaluate in reference space
ELSE
  CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),nVar,PP_N,U    (:,:,:,:),field,ElemID)
END IF ! RefMapping

! Add the interpolated field to the background field
FieldAtParticle(:)        = FieldAtParticle(:)  + field(:)

END SUBROUTINE InterpolateFieldToSingleParticle


END MODULE MOD_Particle_Interpolation
