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

!===================================================================================================================================
! Contains utils required by xy- modules
!===================================================================================================================================
MODULE  MOD_PICInterpolation
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: InterpolateFieldToParticle,InitializeInterpolation,InterpolateFieldToSingleParticle,InterpolateVariableExternalField
!===================================================================================================================================
INTERFACE InitializeInterpolation
  MODULE PROCEDURE InitializeInterpolation
END INTERFACE

INTERFACE InterpolateFieldToParticle
  MODULE PROCEDURE InterpolateFieldToParticle
END INTERFACE

INTERFACE InterpolateVariableExternalField
  MODULE PROCEDURE InterpolateVariableExternalField
END INTERFACE

INTERFACE InterpolateFieldToSingleParticle
  MODULE PROCEDURE InterpolateFieldToSingleParticle
END INTERFACE
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeInterpolation
!===================================================================================================================================
! Initialize the interpolation variables first
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
InterpolationType = GETSTR('PIC-Interpolation-Type','particle_position')
InterpolationElemLoop = GETLOGICAL('PIC-InterpolationElemLoop','.TRUE.')
IF (InterpolationElemLoop) THEN !If user-defined F: F for all procs
  IF (PP_nElems.GT.10) THEN !so far arbitrary threshold...
    InterpolationElemLoop=.FALSE. !switch off for procs with high number of Elems
  END IF
END IF
#if EQNSYSNR == 3   /*Spalart-Allmaras*/
externalField(1:PP_nVar)= GETREALARRAY('PIC-externalField',PP_nVar,'0.,0.,0.,0.,0.,0.')
#else               /*Navier-Stokes*/
externalField(1:PP_nVar)= GETREALARRAY('PIC-externalField',PP_nVar,'0.,0.,0.,0.,0.')
#endif
scaleexternalField= GETREAL('PIC-scaleexternalField','1.0')
externalField=externalField*ScaleExternalField
!SWRITE(*,*) " External fied", externalfield
DoInterpolation   = GETLOGICAL('PIC-DoInterpolation','.TRUE.')
useBGField        = GETLOGICAL('PIC-BG-Field','.FALSE.')

! Variable external field
useVariableExternalField = .FALSE.
FileNameVariableExternalField=GETSTR('PIC-curvedexternalField','none')
IF (FileNameVariableExternalField.EQ.'none') THEN
  FileNameVariableExternalField=GETSTR('PIC-variableexternalField','none')
END IF
IF (FileNameVariableExternalField.NE.'none') THEN
  useVariableExternalField = .TRUE.
  CALL ReadVariableExternalField()
END IF

!--- Allocate arrays for interpolation of fields to particles
SDEALLOCATE(FieldAtParticle)
! Allocate array for rho,momentum,energy,(tke,omega)
ALLOCATE(FieldAtParticle(1:PP_nVar,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
#if USE_RW
SDEALLOCATE(TurbFieldAtParticle)
ALLOCATE(TurbFieldAtParticle(1:nVarTurb1,:PDM%maxParticleNumber), STAT=ALLOCSTAT)
#endif
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
  __STAMP__ &
  ,'ERROR in pic_interpolation.f90: Cannot allocate FieldAtParticle array!',ALLOCSTAT)
END IF

SELECT CASE(TRIM(InterpolationType))
CASE('nearest_blurycenter')
   InterpolationType='nearest_blurrycenter'
CASE('nearest_blurrycenter')
CASE('particle_position_slow')
CASE('particle_position')
CASE('nearest_gausspoint')
CASE DEFAULT
  CALL abort(&
  __STAMP__ &
  ,'Unknown InterpolationType in pic_init.f90')
END SELECT
END SUBROUTINE InitializeInterpolation


SUBROUTINE InterpolateFieldToParticle(doInnerParts)
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                 ONLY:U
USE MOD_Particle_Vars,           ONLY:PartPosRef,PDM,PartState,PEM
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
USE MOD_Mesh_Vars,               ONLY:nElems
USE MOD_PIC_Vars
USE MOD_PICInterpolation_Vars,   ONLY:FieldAtParticle,externalField,DoInterpolation,InterpolationType
USE MOD_PICInterpolation_Vars,   ONLY:InterpolationElemLoop
USE MOD_Eval_xyz,                ONLY:TensorProductInterpolation,GetPositionInRefElem,EvaluateFieldAtPhysPos
#if USE_RW
USE MOD_DG_Vars,                 ONLY:UTurb
USE MOD_Restart_Vars,            ONLY:RestartTurb
USE MOD_Equation_Vars,           ONLY:nVarTurb
#endif
#if USE_MPI
! only required for shape function??
USE MOD_Particle_MPI_Vars,       ONLY:PartMPIExchange
#endif
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL                          :: doInnerParts
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: firstPart,lastPart
REAL                             :: field(PP_nVar)
INTEGER                          :: iPart,iElem
#if USE_RW
REAL                             :: turbField(nVarTurb)
#endif
!===================================================================================================================================
! null field vector
field=0.
#if USE_RW
turbField=0.
#endif

IF(doInnerParts)THEN
  firstPart=1
  lastPart =PDM%ParticleVecLength
ELSE
#if USE_MPI
  firstPart=PDM%ParticleVecLength-PartMPIExchange%nMPIParticles+1
  lastPart =PDM%ParticleVecLength
#else
  firstPart=2
  LastPart =1
#endif /*MPI*/
END IF
! that's wrong
IF(firstPart.GT.lastPart) RETURN

IF (.NOT.InterpolationElemLoop) THEN
  DO iPart = firstPart, LastPart
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:PP_nVar,iPart))
  END DO
  RETURN
END IF

FieldAtParticle(:,firstPart:lastPart) = 0.
FieldAtParticle(1,firstPart:lastPart) = externalField(1)
FieldAtParticle(2,firstPart:lastPart) = externalField(2)
FieldAtParticle(3,firstPart:lastPart) = externalField(3)
FieldAtParticle(4,firstPart:lastPart) = externalField(4)
FieldAtParticle(5,firstPart:lastPart) = externalField(5)

IF (DoInterpolation) THEN                 ! skip if no self fields are calculated
  SELECT CASE(TRIM(InterpolationType))
  CASE('particle_position')
!    IF (DoRefMapping) THEN
    ! particles have already been mapped in deposition, other eval routine used
    DO iElem=1,nElems
      DO iPart=firstPart,LastPart
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
        IF (PEM%Element(iPart).EQ.iElem) THEN
          IF (.NOT.DoRefMapping) THEN
            CALL TensorProductInterpolation(PartState(1:3,iPart),PartPosRef(1:3,iPart),iElem,ForceMode=.TRUE.,PartID=iPart)
          END IF
          !--- evaluate at Particle position
#if USE_RW
          IF (RestartTurb) THEN
            CALL EvaluateFieldAtPhysPos(PartPosRef(1:3,iPart),PP_nVar,PP_N,U(1:PP_nVar,:,:,:,iElem),field(1:PP_nVar),iElem,UTurb(1:nVarTurb,:,:,:,iElem),turbField(1:nVarTurb))
          ELSE
#endif
            CALL EvaluateFieldAtPhysPos(PartPosRef(1:3,iPart),PP_nVar,PP_N,U(1:PP_nVar,:,:,:,iElem),field(1:PP_nVar),iElem)
#if USE_RW
          END IF
#endif
          FieldAtParticle(1:PP_nVar,iPart) = FieldAtParticle(1:PP_nVar,iPart) + field(1:PP_nVar)
        END IF ! Element(iPart).EQ.iElem
      END DO ! iPart
    END DO ! iElem=1,PP_N
!    ELSE ! particles are not yet mapped
!      DO iElem=1,nElems
!        DO iPart=firstPart,LastPart
!          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
!          IF(PEM%Element(iPart).EQ.iElem)THEN
!            !--- evaluate at Particle position
!            CALL GetPositionInRefElem(PartState(1:3,iPart),5,PP_N,U(1:5,:,:,:,iElem),field(1:5),iElem,iPart)
!            FieldAtParticle(:,iPart) = FieldAtParticle(:,iPart) + field(1:6)
!          END IF ! Element(iPart).EQ.iElem
!        END DO ! iPart
!      END DO ! iElem=1,PP_nElems
!    END IF ! DoRefMapping .or. Depositiontype=nearest_gausspoint
  CASE DEFAULT
    CALL abort(&
__STAMP__&
       , 'ERROR: Unknown InterpolationType!')
  END SELECT
END IF

RETURN
END SUBROUTINE InterpolateFieldToParticle


SUBROUTINE InterpolateFieldToSingleParticle(PartID,FieldAtParticle)
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Vars,           ONLY:PartPosRef,PartState,PEM
USE MOD_DG_Vars,                 ONLY:U
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
USE MOD_PIC_Vars!,      ONLY:
USE MOD_PICInterpolation_Vars,   ONLY:externalField,DoInterpolation,InterpolationType
USE MOD_Eval_xyz,                ONLY:TensorProductInterpolation,GetPositionInRefElem,EvaluateFieldAtPhysPos
#if USE_MPI
USE MOD_Mesh_Vars,               ONLY:nElems
#endif
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: FieldAtParticle(1:PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: Field(1:PP_nVar) !,Pos(3)
INTEGER                      :: ElemID
!===================================================================================================================================
  FieldAtParticle(:) = 0.
  FieldAtParticle(1) = externalField(1)
  FieldAtParticle(2) = externalField(2)
  FieldAtParticle(3) = externalField(3)
  FieldAtParticle(4) = externalField(4)
  FieldAtParticle(5) = externalField(5)

IF (DoInterpolation) THEN                 ! skip if no self fields are calculated
  field(1:PP_nVar)=0.
  ElemID=PEM%Element(PartID)
#if USE_MPI
  IF(ElemID.GT.nElems) RETURN
#endif
  SELECT CASE(TRIM(InterpolationType))
  CASE('particle_position')
      ! particles have already been mapped in deposition, other eval routine used
      IF(.NOT.DoRefMapping)THEN
        CALL TensorProductInterpolation(PartState(1:3,PartID),PartPosRef(PartID,1:3),ElemID,ForceMode=.TRUE.,PartID=PartID)
      END IF
      !--- evaluate at Particle position
      CALL EvaluateFieldAtPhysPos(PartPosRef(1:3,PartID),PP_nVar,PP_N,U(1:PP_nVar,:,:,:,ElemID),field(1:PP_nVar),ElemID)
      FieldAtParticle(1:PP_nVar) = FieldAtParticle(1:PP_nVar) + field(1:PP_nVar)
  CASE DEFAULT
    CALL abort(&
__STAMP__&
    , 'ERROR: Unknown InterpolationType!')
  END SELECT
END IF

RETURN
END SUBROUTINE InterpolateFieldToSingleParticle


SUBROUTINE ReadVariableExternalField()
!===================================================================================================================================
! ATTENTION: The extrenal field needs to be defined on equidistant data-points
! Usage Information
! The file for the variable Bfield contains only the z coordinates and the static Bz-field
! Use the following format F8.5,1x,F8.5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars, ONLY:VariableExternalField,DeltaExternalField,nIntPoints,FileNameVariableExternalField
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ioUnit, ii, err, ncounts
REAL                  :: dummy, diff_comp, diff_check
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A,3X,A,65X,A)') ' INITIALIZATION OF VARIABLE EXTERNAL FIELD FOR PARTICLES '
!OPEN(NEWUNIT=ioUnit,FILE=VariableExternalField,STATUS='OLD',FORM='FORMATTED')
OPEN(NEWUNIT=ioUnit,FILE=FileNameVariableExternalField,STATUS='OLD')
err = 0
ncounts = 0
DO WHILE (err.EQ.0)
  READ(ioUnit,*,IOSTAT = err) dummy
  IF (err.EQ.-1) THEN
    EXIT
  END IF
  ERR = 0
  ncounts = ncounts + 1
END DO
REWIND(ioUnit)
nIntPoints = ncounts
! allocate needed space
ALLOCATE(VariableExternalField(1:2,1:nIntPoints))
DO ii = 1, ncounts
  read(ioUnit,*) VariableExternalField(1,ii) , VariableExternalField(2,ii)
  IF (ii.GE.2) THEN
    diff_comp  = VariableExternalField(1,2)  - VariableExternalField(1,1)
    diff_check = VariableExternalField(1,ii) - VariableExternalField(1,ii-1)
    IF( (.NOT.ALMOSTEQUALRELATIVE(diff_comp,diff_check,1E-5)) .AND. ((diff_comp.GT.0.0).AND.(diff_check.GT.0.0)) )THEN
      SWRITE(UNIT_stdOut,'(A)') "ReadVariableExternalField: Non-equidistant OR non-increasing points for variable external field."
      ! partns
!      SWRITE(UNIT_stdOut,OUTPUTFORMAT) diff_comp
!      SWRITE(UNIT_stdOut,OUTPUTFORMAT) diff_check
      CALL abort(&
__STAMP__&
        ,' Error in dataset!')
    END IF
  END IF
END DO
CLOSE (ioUnit)

!IF (VariableExternalField(1,1) .NE.0) THEN
  !CALL abort(&
!__STAMP__&
!,  &
      !"ERROR: Points have to start at 0.")
!END IF
IF(ncounts.GT.1) THEN
  DeltaExternalField = VariableExternalField(1,2)  - VariableExternalField(1,1)
  SWRITE(UNIT_stdOut,'(A,1X,E25.14E3)') ' Delta external field: ',DeltaExternalField
  IF(DeltaExternalField.LE.0) THEN
    SWRITE(*,'(A)') ' ERROR: wrong sign in external field delta-x'
  END IF
ELSE
  CALL abort(&
__STAMP__&
, &
      " ERROR: not enough data points in variable external field file!")
END IF
SWRITE(UNIT_stdOut,'(A,I4.0,A)')' Found ', ncounts,' data points.'
SWRITE(UNIT_stdOut,'(A)')' ...VARIABLE EXTERNAL FIELD INITIALIZATION DONE'
END SUBROUTINE ReadVariableExternalField


PURE FUNCTION InterpolateVariableExternalField(Pos)
!===================================================================================================================================
! interpolates Variable external field to z position
! NO z-values smaller than VariableExternalField(1,1) are allowed!
!===================================================================================================================================
! MODULES
!USE MOD_Globals
USE MOD_PICInterpolation_Vars   ,ONLY:DeltaExternalField,nIntPoints,VariableExternalField
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: Pos ! partilce z position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateVariableExternalField  ! Bz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iPos
!===================================================================================================================================
iPos = INT((Pos-VariableExternalField(1,1))/DeltaExternalField) + 1
IF(iPos.GE.nIntPoints)THEN ! particle outside of range (greater -> use constant value)
  InterpolateVariableExternalField = VariableExternalField(2,nIntPoints)
ELSEIF(iPos.LT.1)THEN ! particle outside of range (lower -> use constant value)
  InterpolateVariableExternalField = VariableExternalField(2,1)
ELSE ! Linear Interpolation between iPos and iPos+1 B point
  InterpolateVariableExternalField = (VariableExternalField(2,iPos+1) - VariableExternalField(2,iPos)) & !  dy
                                   / (VariableExternalField(1,iPos+1) - VariableExternalField(1,iPos)) & ! /dx
                             * (Pos - VariableExternalField(1,iPos) ) + VariableExternalField(2,iPos)    ! *(z - z_i) + z_i
END IF
END FUNCTION InterpolateVariableExternalField

END MODULE MOD_PICInterpolation
