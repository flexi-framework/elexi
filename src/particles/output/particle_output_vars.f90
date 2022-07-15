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
!> Contains global variables provided by the output routines
!===================================================================================================================================
MODULE MOD_Particle_Output_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                      :: ParticleOutputInitIsDone=.FALSE.
LOGICAL                      :: WriteFieldsToVTK
LOGICAL                      :: OutputMesh

INTEGER                      :: locnPart                                 ! number of particles on processors
INTEGER                      :: offsetnPart                              ! offset of particles on processors

INTEGER                      :: PartDataVarSpecies                       ! Position of the species info in the PartData array
INTEGER                      :: PartDataVarStart                         ! Start of the variable part   in the PartData array
INTEGER                      :: PartDataVarShift                         ! Shift of the variable part   in the PartData array
!===================================================================================================================================
END MODULE MOD_Particle_Output_Vars

