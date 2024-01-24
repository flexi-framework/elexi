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
! Contains the constant Advection Velocity Vector used for the linear scalar advection equation
!===================================================================================================================================
MODULE MOD_Particle_Interpolation_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                               :: PartInterpolationInitIsDone=.FALSE.
LOGICAL                               :: DoInterpolation          ! Flag for interpolation
LOGICAL                               :: InterpolationElemLoop    ! Interpolate with outer iElem-loop (not for many Elems per proc!)
REAL,ALLOCATABLE                      :: FieldAtParticle(:,:)     ! (6,PDM%maxParticleNumber) 1nd index: rho,u_x,u_y,u_z,p,T
#if USE_EXTEND_RHS || USE_FAXEN_CORR
REAL,ALLOCATABLE                      :: GradAtParticle(:,:,:)    ! (3,3,PDM%maxParticleNumber) 1nd index: rho*(u_x,u_y,u_z)
#endif
#if USE_RW
REAL,ALLOCATABLE                      :: turbFieldAtParticle(:,:) ! (PDM%maxParticleNumber,2) 2nd index: k,epsilon
#endif
!===================================================================================================================================
END MODULE MOD_Particle_Interpolation_Vars
