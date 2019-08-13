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

!===================================================================================================================================
! Contains global variables used by the Analyze modules.
!===================================================================================================================================
MODULE MOD_Particle_Erosion_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: ParticleErosionInitIsDone = .FALSE.
LOGICAL                       :: PartTrackErosion
LOGICAL                       :: PartTrackReflection = .TRUE.            ! Flag if reflections should be counted
LOGICAL                       :: ErosionOutlet                           ! Flag if erosion should be tracked on outlet
LOGICAL                       :: ErosionRestart                          ! Flag if we are restarting erosion tracking

INTEGER                       :: nErosionVars                            ! Number of Vars = nErosionVars * (nSPecies + 1)
!===================================================================================================================================
END MODULE MOD_Particle_Erosion_Vars
