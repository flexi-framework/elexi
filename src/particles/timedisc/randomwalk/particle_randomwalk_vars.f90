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

!===================================================================================================================================
! Contains the variables for particle random walk models
!===================================================================================================================================
MODULE MOD_Particle_RandomWalk_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
#if USE_RW
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                                :: ParticleRWInitIsDone=.FALSE.
LOGICAL                                :: RWinUse                          ! Flag for active RW
!Random Walk method
INTEGER                                :: nRWVars                          ! number of variables in TurbPartState
CHARACTER(40)                          :: RWModel                          ! specifying Keyword for RW model
CHARACTER(40)                          :: RWTime                           ! time stepping mode for RW model
!===================================================================================================================================
#endif /*USE_RW*/
END MODULE MOD_Particle_RandomWalk_Vars
