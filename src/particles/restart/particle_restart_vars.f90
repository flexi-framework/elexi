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

!===================================================================================================================================
! Contains the particle restart variables
!===================================================================================================================================
MODULE MOD_Particle_Restart_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                               :: EmissionTime                      ! time at which particle emission started
INTEGER                            :: PP_nVarPartState                  ! Number of particle vars to be copied into PartState
LOGICAL                            :: PartIntExists      = .FALSE.      ! Flag if restart file has saved part ints
LOGICAL                            :: PartDataExists     = .FALSE.      ! Flag if restart file has saved part data
LOGICAL                            :: TurbPartDataExists = .FALSE.      ! Flag if restart file has saved turbulent part data
!===================================================================================================================================
END MODULE MOD_Particle_Restart_Vars
