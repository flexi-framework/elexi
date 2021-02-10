!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
! Module for calculation of particle forces on surfaces for convergence analysis
!===================================================================================================================================
MODULE MOD_CalcWallParticles_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Variables for the specific analyze routines
LOGICAL              :: doCalcWallParticles    =.FALSE.      !< marks if particle properties at walls shall be computed
LOGICAL              :: doWriteWallParticles   =.FALSE.      !< marks if particle properties shall be written to a file

CHARACTER(LEN=255),ALLOCATABLE :: FileName_WallPart(:)       !< output files for wall velocities per BC
!==================================================================================================================================
END MODULE MOD_CalcWallParticles_Vars