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
!> Contains variables to exchange data using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_Particle_MPI_Shared_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
#if USE_MPI_SHARED
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! Mesh
REAL,POINTER    :: BezierControlPoints3D_Shared(:,:,:,:)      !> BezierControlPoints3D mapping on current node
INTEGER         :: BezierControlPoints3D_Shared_Win           !> Pointer to shared memory window
REAL,POINTER    :: XiEtaZetaBasis_Shared(:,:,:)               !> BezierControlPoints3D mapping on current node
INTEGER         :: XiEtaZetaBasis_Shared_Win                  !> Pointer to shared memory window
REAL,POINTER    :: slenXiEtaZetaBasis_Shared(:,:)             !> BezierControlPoints3D mapping on current node
INTEGER         :: slenXiEtaZetaBasis_Shared_Win              !> Pointer to shared memory window


#endif /* MPI_SHARED */
END MODULE
