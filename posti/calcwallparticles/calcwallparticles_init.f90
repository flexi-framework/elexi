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
! Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_CalcWallParticles_Init
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

! Private Part --------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersParticles
  MODULE PROCEDURE DefineParametersParticles
END INTERFACE

INTERFACE InitParticles
  MODULE PROCEDURE InitParticles
END INTERFACE

INTERFACE FinalizeParticles
  MODULE PROCEDURE FinalizeParticles
END INTERFACE

PUBLIC::DefineParametersParticles
PUBLIC::InitParticles,FinalizeParticles
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particles
!==================================================================================================================================
SUBROUTINE DefineParametersParticles()
! MODULES
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE Mod_Particle_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Particle")

CALL prms%CreateIntOption('BezierClipLineVectorMethod', "")
CALL prms%CreateIntOption('BezierClipMaxIntersec', "")
CALL prms%CreateIntOption('BezierClipMaxIter', "")
CALL prms%CreateIntOption('BezierElevation', "")
CALL prms%CreateIntOption('BezierNewtonGuess', "")
CALL prms%CreateIntOption('BezierNewtonMaxIter', "")
CALL prms%CreateIntOption('BezierSampleN', "")
CALL prms%CreateIntOption('Part-nBounds', "")
CALL prms%CreateIntOption('Part-nSpecies',                  'Number of species in part',                    '1')
CALL prms%CreateIntOption('Part-maxParticleNumber',         'Max number of particles in part',              '1')

CALL prms%CreateLogicalOption('Part-vMPF', "")
CALL prms%CreateLogicalOption('Part-WriteMacroValues',"")
CALL prms%CreateLogicalOption('Part-WriteMacroVolumeValues',"")
CALL prms%CreateLogicalOption('Part-WriteMacroSurfaceValues',"")

CALL prms%SetSection("Particle Boundaries")

CALL prms%CreateStringOption(   'Part-BoundaryType'  &
                                , 'Used boundary condition for boundary[$].\n'//&
                                  '- open\n'//&
                                  '- reflective\n'//&
                                  '- periodic\n'//&
                                 'If condition=open, the following parameters are'//&
                                  ' used: (Part-Boundary[$]-=PB) PB-Ambient ,PB-AmbientTemp,PB-AmbientMeanPartMass,'//&
                                  'PB-AmbientVelo,PB-AmbientDens,PB-AmbientDynamicVisc,PB-AmbientThermalCond,PB-Voltage\n'//&
                                 'If condition=reflective: PB-MomentumACC,PB-WallTemp,PB-TransACC,PB-VibACC,PB-RotACC,'//&
                                  'PB-WallVelo,Voltage,SpeciesSwaps.If condition=periodic:Part-nPeriodicVectors,'//&
                                  'Part-PeriodicVector[$]', multiple=.TRUE.)
CALL prms%CreateStringOption(   'Part-BoundaryName'  &
                                , 'No Default. Source Name of Boundary[i]. Has to be selected for all'//&
                                  'nBounds. Has to be same name as defined in preproc tool', multiple=.TRUE.)

CALL prms%CreateRealArrayOption('Part-PeriodicVector[$]'      , 'TODO-DEFINE-PARAMETER\nVector for periodic boundaries.'//&
                                                                   'Has to be the same as defined in preproc.ini in their'//&
                                                                   ' respective order. ', '1. , 0. , 0.', numberedmulti=.TRUE.)

! Ambient condition=================================================================================================================
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-AmbientCondition'  &
                                , 'Use ambient condition (condition "behind" boundary).', '.FALSE.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-AmbientConditionFix'  &
                                , 'TODO-DEFINE-PARAMETER', '.TRUE.', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption('Part-Boundary[$]-AmbientVelo'  &
                                , 'Ambient velocity', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-WallVelo'  &
                                , 'Velocity (global x,y,z in [m/s]) of reflective particle boundary [$].' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Boundary[$]-AmbientDens'  &
                                , 'Ambient density', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-AmbientDynamicVisc'  &
                                , 'Ambient dynamic viscosity', '1.72326582572253E-5', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Particles-RandomSeed[$]'     , 'Seed [$] for Random Number Generator', '1', numberedmulti=.TRUE.)

END SUBROUTINE DefineParametersParticles

!===================================================================================================================================
! Glue Subroutine for particle initialization
!===================================================================================================================================
SUBROUTINE InitParticles()
! MODULES
USE MOD_Globals
USE Mod_Particle_Globals
USE MOD_ReadInTools
USE MOD_IO_HDF5,                    ONLY: AddToElemData
USE MOD_DSMC_Vars,                  ONLY: DSMC
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Vars,              ONLY: ParticlesInitIsDone,WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Sampling, ONLY: InitParticleBoundarySampling
#if USE_MPI
USE MOD_Particle_MPI,               ONLY: InitParticleCommSize
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ParticlesInitIsDone)THEN
   SWRITE(*,*) "InitParticles already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLES ...'

Call InitParticleGlobals
CALL InitializeVariables()

! Initialize surface sampling
IF (WriteMacroSurfaceValues.OR.DSMC%CalcSurfaceVal) THEN
  CALL InitParticleBoundarySampling()
END IF

#if USE_MPI
! has to be called AFTER InitializeVariables and InitDSMC
CALL InitParticleCommSize()
#endif

ParticlesInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticles

!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
SUBROUTINE InitializeVariables()
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars,ONLY:PartBound,nPartBound
USE MOD_Mesh_Vars,             ONLY:BoundaryName,BoundaryType, nBCs
USE MOD_Particle_Surfaces_Vars,ONLY:BCdata_auxSF
USE MOD_Particle_Tracking_Vars,ONLY:DoRefMapping
USE MOD_Particle_Mesh,         ONLY:InitFIBGM
#if USE_MPI
USE MOD_Particle_MPI,          ONLY:InitEmissionComm
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPartBound
INTEGER               :: iPBC, iBC
INTEGER               :: ALLOCSTAT
CHARACTER(32)         :: hilf
CHARACTER(200)        :: tmpString
CHARACTER(200), ALLOCATABLE        :: tmpStringBC(:)
!===================================================================================================================================

! Read basic particle parameter
PDM%maxParticleNumber = GETINT('Part-maxParticleNumber','1')

IF(DoRefMapping)THEN
  ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,' Cannot allocate partposref!')
  PartPosRef=-888.
END IF

ALLOCATE(PartState(1:PDM%maxParticleNumber,1:6)       , &
         PartReflCount(1:PDM%maxParticleNumber)      , &
         LastPartPos(1:PDM%maxParticleNumber,1:3)     , &
         Pt(1:PDM%maxParticleNumber,1:3)              , &
         PartSpecies(1:PDM%maxParticleNumber)         , &
         PDM%ParticleInside(1:PDM%maxParticleNumber)  , &
         PDM%nextFreePosition(1:PDM%maxParticleNumber), &
         PDM%dtFracPush(1:PDM%maxParticleNumber)      , &
         PDM%IsNewPart(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF

!> output of macroscopic values
WriteMacroValues        = GETLOGICAL('Part-WriteMacroValues',       '.FALSE.')
WriteMacroVolumeValues  = GETLOGICAL('Part-WriteMacroVolumeValues', '.FALSE.')
WriteMacroSurfaceValues = GETLOGICAL('Part-WriteMacroSurfaceValues','.FALSE.')
IF(WriteMacroValues)THEN
  WriteMacroVolumeValues  = .TRUE.
  WriteMacroSurfaceValues = .TRUE.
ELSE IF((WriteMacroVolumeValues.AND.WriteMacroSurfaceValues).AND.(.NOT.WriteMacroValues))THEN
  WriteMacroValues        = .TRUE.
END IF

nSpecies = GETINT('Part-nSpecies','1')
!nSpecies = CNTSTR('Part-Species-SpaceIC')

IF (nSpecies.LE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nSpecies .LE. 0:', nSpecies)
END IF

! Read in boundary parameters
nPartBound = CountOption('Part-BoundaryName')

ALLOCATE(PartBound%SourceBoundName(1:nPartBound))
ALLOCATE(PartBound%TargetBoundCond(1:nPartBound))
ALLOCATE(PartBound%MomentumACC(1:nPartBound))
ALLOCATE(PartBound%WallTemp(1:nPartBound))
ALLOCATE(PartBound%WallVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientCondition(1:nPartBound))
ALLOCATE(PartBound%AmbientConditionFix(1:nPartBound))
ALLOCATE(PartBound%AmbientTemp(1:nPartBound))
ALLOCATE(PartBound%AmbientVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientDens(1:nPartBound))
ALLOCATE(PartBound%SolidState(1:nPartBound))
PartBound%SolidState(1:nPartBound)=.FALSE.

ALLOCATE(tmpStringBC                  (1:nBCs))

! Loop over all particle boundaries and get information
DO iPartBound=1,nPartBound
  tmpStringBC(iPartBound) = TRIM(GETSTR('Part-BoundaryName'))
END DO

!--
DO iBC=1,nBCs
  IF (BoundaryType(iBC,1).EQ.0) THEN
    PartBound%TargetBoundCond(iBC) = -1
    SWRITE(*,*)"... PartBound",iBC,"is internal bound, no mapping needed"
    CYCLE
  END IF

  PartBound%SourceBoundName(iBC) = BoundaryName(iBC)

  SELECT CASE (BoundaryType(iBC, BC_TYPE))
  CASE(1)
    PartBound%SourceBoundType(iBC)='periodic'
  CASE(3,4)
    PartBound%SourceBoundType(iBC)='reflective'
  CASE(9)
    PartBound%SourceBoundType(iBC)='symmetry'
  CASE DEFAULT
    PartBound%SourceBoundType(iBC)='open'
  END SELECT
  DO iPartBound=1,nPartBound
    IF (TRIM(BoundaryName(iBC)).EQ.tmpStringBC(iPartBound)) THEN
      PartBound%SourceBoundType(iBC) = TRIM(GETSTR('Part-BoundaryType'))
      PartBound%SourceBoundName(iBC) = tmpStringBC(iPartBound)
    END IF
  END DO
END DO

DO iBC=1,nBCs
  SELECT CASE (TRIM(tmpString))
  CASE('open')
     PartBound%TargetBoundCond(iPartBound) = PartBound%OpenBC          ! definitions see typesdef_pic
!     PartBound%AmbientCondition(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientCondition','.FALSE.')
!     IF(PartBound%AmbientCondition(iPartBound)) THEN
!       PartBound%AmbientConditionFix(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientConditionFix','.TRUE.')
!       PartBound%AmbientVelo(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientVelo',3,'0. , 0. , 0.')
!       PartBound%AmbientDens(iPartBound) = GETREAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientDens','0')
!       PartBound%AmbientDynamicVisc(iPartBound)=&
!           GETREAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientDynamicVisc','1.72326582572253E-5') ! N2:T=288K
!     END IF
  CASE('reflective')
     PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
!     PartBound%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-WallVelo',3,'0. , 0. , 0.')
  CASE('periodic')
     PartBound%TargetBoundCond(iPartBound) = PartBound%PeriodicBC
  CASE DEFAULT
     SWRITE(*,*) ' Boundary does not exists: ', TRIM(tmpString)
     CALL abort(&
__STAMP__&
         ,'Particle Boundary Condition does not exist')
  END SELECT
END DO

ALLOCATE(PEM%Element(1:PDM%maxParticleNumber), PEM%lastElement(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(&
__STAMP__&
  ,' Cannot allocate PEM arrays!')
END IF

!-- Finalizing InitializeVariables
CALL InitFIBGM()

#if USE_MPI
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/

!!-- Floating Potential
ALLOCATE(BCdata_auxSF(1:nPartBound))
DO iPartBound=1,nPartBound
  BCdata_auxSF(iPartBound)%SideNumber=-1 !init value when not used
  BCdata_auxSF(iPartBound)%GlobalArea=0.
  BCdata_auxSF(iPartBound)%LocalArea=0.
END DO
!nDataBC_CollectCharges=0

SDEALLOCATE(tmpStringBC)

END SUBROUTINE InitializeVariables

!----------------------------------------------------------------------------------------------------------------------------------!
! finalize particle variables
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE FinalizeParticles()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!SDEALLOCATE(Pt_temp)
SDEALLOCATE(PartPosRef)
SDEALLOCATE(RandomVec)
SDEALLOCATE(PartState)
SDEALLOCATE(LastPartPos)
SDEALLOCATE(PartSpecies)
SDEALLOCATE(Pt)
SDEALLOCATE(PDM%ParticleInside)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%dtFracPush)
SDEALLOCATE(PDM%IsNewPart)
SDEALLOCATE(Species)
SDEALLOCATE(PartBound%SourceBoundName)
SDEALLOCATE(PartBound%TargetBoundCond)
SDEALLOCATE(PartBound%MomentumACC)
SDEALLOCATE(PartBound%WallTemp)
SDEALLOCATE(PartBound%WallVelo)
SDEALLOCATE(PartBound%AmbientCondition)
SDEALLOCATE(PartBound%AmbientConditionFix)
SDEALLOCATE(PartBound%AmbientTemp)
SDEALLOCATE(PartBound%AmbientVelo)
SDEALLOCATE(PartBound%AmbientDens)
SDEALLOCATE(PEM%Element)
SDEALLOCATE(PEM%lastElement)
SDEALLOCATE(PEM%pStart)
SDEALLOCATE(PEM%pNumber)
SDEALLOCATE(PEM%pEnd)
SDEALLOCATE(PEM%pNext)
END SUBROUTINE FinalizeParticles


END MODULE MOD_CalcWallParticles_Init
