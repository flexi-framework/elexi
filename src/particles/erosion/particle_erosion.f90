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
! Contains some wrapper functions for analysis of particle erosion
!===================================================================================================================================
MODULE MOD_Particle_Erosion
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersParticleErosion
  MODULE PROCEDURE DefineParametersParticleErosion
END INTERFACE

INTERFACE InitParticleErosion
  MODULE PROCEDURE InitParticleErosion
END INTERFACE

INTERFACE FinalizeParticleErosion
  MODULE PROCEDURE FinalizeParticleErosion
END INTERFACE

PUBLIC :: DefineParametersParticleErosion
PUBLIC :: InitParticleErosion
PUBLIC :: FinalizeParticleErosion
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for erosion if particles
!==================================================================================================================================
SUBROUTINE DefineParametersParticleErosion()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Erosion")

CALL prms%CreateLogicalOption(  'DoErosion',                      'Flag, if the current calculation requires erosion tracking.')
CALL prms%CreateLogicalOption(  'Particles-AnalyzeSurfCollis',    'Output of collided/swaped particles during Sampling period?')
CALL prms%CreateLogicalOption(  'Particles-ErosionOutlet',        'Flag, if erosion should be tracked on outlet.\n'//&
                                                                  'WARNING: Only takes BC named <outlet> into account.')
CALL prms%CreateLogicalOption(  'Particles-CalcSurfaceVal',       'Set [T] to activate sampling, analyze and h5 output for'//&
                                                                  ' surfaces. Therefore either time fraction or iteration'//&
                                                                  'sampling have to be enabled as well.', '.FALSE.')
CALL prms%CreateIntOption(      'Particles-nSurfSample',          'Define polynomial degree of particle BC sampling. Default:'//&
                                                                  ' NGeo', '1')
CALL prms%CreateRealOption(     'Part-TimeFracForSampling',       'Set value greater 0.0 to enable TIME DEPENDANT sampling. The'//&
                                                                  ' given simulation time fraction will be sampled. Sampling'//&
                                                                  ' starts after TEnd*(1-Part-TimefracForSampling).\n'//&
                                                                  'Can not be enabled together with Part-WriteMacroValues.' , '0.0')

END SUBROUTINE DefineParametersParticleErosion


SUBROUTINE InitParticleErosion()
!===================================================================================================================================
! Initializes variables necessary for analyze subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools,            ONLY: GETLOGICAL, GETINT
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Vars,          ONLY: WriteMacroSurfaceValues
!  IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (ParticleErosionInitIsDone) THEN
    CALL abort(__STAMP__,&
    'InitParticleErosion already called.',999,999.)
    RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE EROSION...'

doParticleErosionTrack = GETLOGICAL('DoErosion','F')

! Switch surface macro values flag to .TRUE. for erosion tracking
IF (doParticleErosionTrack) WriteMacroSurfaceValues = .TRUE.

ParticleErosionInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE EROSION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleErosion


SUBROUTINE FinalizeParticleErosion()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Erosion_Vars,ONLY:ParticleErosionInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ParticleErosionInitIsDone = .FALSE.

END SUBROUTINE FinalizeParticleErosion


END MODULE MOD_Particle_Erosion
