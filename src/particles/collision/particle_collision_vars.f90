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
#include "flexi.h"
#include "particle.h"

!===================================================================================================================================
! Module containing general particles collisions routines
!===================================================================================================================================
MODULE MOD_Particle_Collision_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

#if PARTICLES_COUPLING == 4
!----------------------------------------------------------------------------------------------------------------------------------
TYPE tParticleCollisionModel
  REAL                                   :: e        = 1.0              ! restitution coefficient [0;1]
  LOGICAL                                :: Friction = .FALSE.          ! friction during particle collision
  REAL                                   :: f        = 0.0              ! friction coefficient
END TYPE tParticleCollisionModel

TYPE(tParticleCollisionModel)            :: PartCollisionModel          ! parameters for the collision model

! #if USE_MPI
INTEGER                                  :: nComputeNodeNeighElems      !> size of the CN elem to CN neighbor elem mapping
INTEGER,ALLOCPOINT,DIMENSION(:)          :: Neigh_nElems                !>
INTEGER,ALLOCPOINT,DIMENSION(:)          :: Neigh_offsetElem            !>
#if USE_MPI
INTEGER,ALLOCPOINT                       :: Neigh_nElems_Shared(:)      !>
INTEGER,ALLOCPOINT                       :: Neigh_offsetElem_Shared(:)
#endif /*USE_MPI*/
INTEGER,ALLOCPOINT,DIMENSION(:)          :: CNElem2CNNeighElem          !> Compute Nodes mapping 1:nTotal -> 1:nGlobal
INTEGER,ALLOCPOINT,DIMENSION(:)          :: CNNeighElem2CNElem          !> Reverse Mapping
! #endif /*USE_MPI*/

! Shared particle arrays
REAL   ,ALLOCPOINT,DIMENSION(:,:)        :: PartData_Shared             !>
INTEGER,ALLOCPOINT,DIMENSION(:,:)        :: PartInt_Shared              !>
! REAL   ,ALLOCPOINT,DIMENSION(:)          :: PartData_Shared_Flat        !>
! INTEGER,ALLOCPOINT,DIMENSION(:)          :: PartInt_Shared_Flat         !>

#if USE_MPI
INTEGER           :: Neigh_nElems_Shared_Win
INTEGER           :: Neigh_offsetElem_Shared_Win
INTEGER           :: CNElem2CNNeighElem_Win
INTEGER           :: CNNeighElem2CNElem_Win

! Shared particle arrays
INTEGER           :: PartData_Shared_Win
INTEGER           :: PartInt_Shared_Win
#endif /*USE_MPI*/

#endif /*PARTICLES_COUPLING == 4*/

END MODULE MOD_Particle_Collision_Vars
