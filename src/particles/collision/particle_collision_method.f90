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
MODULE MOD_Particle_Collision_Method
! MODULES
IMPLICIT NONE
PRIVATE

#if PARTICLES_COUPLING == 4

INTERFACE ComputeHardSphereCollision
  MODULE PROCEDURE ComputeHardSphereCollision
END INTERFACE

INTERFACE ComputeParticleCollisions
  MODULE PROCEDURE ComputeParticleCollisions
END INTERFACE

PUBLIC :: ComputeParticleCollisions

CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE ComputeParticleCollisions()
!===================================================================================================================================
!> Compute particle collisions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Vars,            ONLY: PDM
#if USE_MPI
USE MOD_MPI_Shared_Vars,          ONLY: MPI_COMM_SHARED
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! ElemID = PDM%Element(PartID)
! DO iElem = Neigh_offsetElem(ElemID)+1:Neigh_offsetElem(ElemID)+Neigh_nElems(ElemID)
!   CNNeighElemID = CNElem2CNNeighElem(iElem)
! END DO
! CALL ComputeHardSphereCollision(iPart1,iPart2,dtColl)

! Done with particle collisions. Deallocate SHM windows
! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

CALL MPI_WIN_UNLOCK_ALL(PartInt_Shared_Win ,iError)
CALL MPI_WIN_UNLOCK_ALL(PartData_Shared_Win,iError)
CALL MPI_WIN_FREE(PartInt_Shared_Win ,iError)
CALL MPI_WIN_FREE(PartData_Shared_Win,iError)

CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

! Then, free the pointers or arrays
MDEALLOCATE(PartInt_Shared)
MDEALLOCATE(PartData_Shared)

END SUBROUTINE ComputeParticleCollisions


SUBROUTINE ComputeHardSphereCollision(iPart1,iPart2,dtColl)
!===================================================================================================================================
!> Apply hard sphere collision model
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Utils                   , ONLY: ALMOSTZERO
USE MOD_Particle_Collision_Vars
USE MOD_Particle_Globals        , ONLY: UNITVECTOR
USE MOD_Particle_Vars           , ONLY: PartState, Species, PartSpecies, LastPartPos!, PDM
!USE MOD_TimeDisc_Vars           , only: CurrentStage, t, iter, dt, RKc
! USE MOD_Particle_Mesh_Vars      , ONLY: GEO
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: iPart1, iPart2
REAL   , INTENT(IN)                :: dtColl
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: n_loc(3),t_loc(3)
REAL                               :: v_rel_contactpoint(3)
REAL                               :: vt_rel_contactpoint(3),vt_rel_contactpoint_norm
REAL                               :: m_red,m1,m2
REAL                               :: Jn, Jt,J(3)
#if USE_PARTROT
REAL                               :: Jxn_loc(3)
#endif /*USE_PARTROT*/
! INTEGER                            :: iPart1Track, iPart2Track
! INTEGER                            :: iPartPeriodicShiftLoc, jPeriodic
!===================================================================================================================================

! ! apply periodic vectors (if there are some)
! IF (GEO%nPeriodicVectors.GT.0) THEN
!   iPartPeriodicShiftLoc = PartCollisionProc%iPartPeriodicShift(iCollision)
!   DO jPeriodic=1,GEO%nPeriodicVectors
!     PartState(PART_POSV,iPart2) = PartState(PART_POSV,iPart2) - REAL(MOD(iPartPeriodicShiftLoc, 3) - 1) * GEO%PeriodicVectors(1:3,jPeriodic)
!     iPartPeriodicShiftLoc = iPartPeriodicShiftLoc / 3
!   END DO
! END IF

! Move the particles back to their positions at the collision time
LastPartPos(PART_POSV,iPart1) = LastPartPos(PART_POSV,iPart1) + dtColl * PartState(PART_VELV,iPart1)
LastPartPos(PART_POSV,iPart2) = LastPartPos(PART_POSV,iPart2) + dtColl * PartState(PART_VELV,iPart2)

! Compute normal vector from iPart1 to iPart2
n_loc = UNITVECTOR(PartState(PART_POSV,iPart2) - PartState(PART_POSV,iPart1))

! Relative velocity of the contact point just before the collision in iPart2's reference frame
v_rel_contactpoint = PartState(PART_VELV,iPart1) - PartState(PART_VELV,iPart2)
#if USE_PARTROT
! Take particle rotation into account
v_rel_contactpoint = v_rel_contactpoint + &
                     PartState(PART_DIAM,iPart1) * 0.5 * CROSSPRODUCT(PartState(PART_AMOMV,iPart1),n_loc) + &
                     PartState(PART_DIAM,iPart2) * 0.5 * CROSSPRODUCT(PartState(PART_AMOMV,iPart2),n_loc)
#endif /*USE_PARTROT*/

! Compute reduced mass of particle pair
m1 = MASS_SPHERE(Species(PartSpecies(iPart1))%DensityIC,PartState(PART_DIAM,iPart1))
m2 = MASS_SPHERE(Species(PartSpecies(iPart2))%DensityIC,PartState(PART_DIAM,iPart2))
m_red = m1*m2/(m1+m2)

! Compute normal component Jn of impulsive contact force J
Jn = - m_red * (1+PartCollisionModel%e) * DOT_PRODUCT(PartState(PART_VELV,iPart1) - PartState(PART_VELV,iPart2),n_loc)

! Compute tangent unit vector t_loc that is in the same plan as n_loc and v_rel_contactpoint
vt_rel_contactpoint = v_rel_contactpoint - DOT_PRODUCT(v_rel_contactpoint,n_loc) * n_loc
vt_rel_contactpoint_norm = NORM2(vt_rel_contactpoint)
IF (ALMOSTZERO(vt_rel_contactpoint_norm)) THEN; t_loc = 0.
ELSE;                                           t_loc = vt_rel_contactpoint / vt_rel_contactpoint_norm
END IF

! Compute friction (Jt) if enabled
IF (PartCollisionModel%Friction .AND. vt_rel_contactpoint_norm.NE.0.) THEN
  ! as particles can slide on each other, the fact that they can stop sliding during the collision due to friction is taken into account
  Jt = MAX(PartCollisionModel%f * Jn, - 2. / 7. * m_red * vt_rel_contactpoint_norm)
ELSE
  Jt = 0.
END IF

! Compute impulsive contact force J
J = Jn * n_loc + Jt * t_loc

! Apply jump relation to momentum equation, i.e., update velocity of particles
PartState(PART_VELV,iPart1) = PartState(PART_VELV,iPart1) + J / m1
PartState(PART_VELV,iPart2) = PartState(PART_VELV,iPart2) - J / m2

! Apply jump relation also to angular momentum if particles can rotate
#if USE_PARTROT
Jxn_loc = CROSSPRODUCT(n_loc, J)
PartState(PART_AMOMV,iPart1) = PartState(PART_AMOMV,iPart1) + Jxn_loc * 0.5 / (0.1 * m1 * PartState(PART_DIAM,iPart1))
PartState(PART_AMOMV,iPart2) = PartState(PART_AMOMV,iPart2) + Jxn_loc * 0.5 / (0.1 * m2 * PartState(PART_DIAM,iPart2))
#endif /*USE_PARTROT*/

! for LSERK algorithm
! PDM%IsNewPart(iPart1) = .TRUE.
! PDM%IsNewPart(iPart2) = .TRUE.

! Update the particle's position at time tStage+dtloc
PartState(PART_POSV,iPart1) = LastPartPos(PART_POSV,iPart1) + dtColl * PartState(PART_VELV,iPart1)
PartState(PART_POSV,iPart2) = LastPartPos(PART_POSV,iPart2) + dtColl * PartState(PART_VELV,iPart2)

! PartState(PART_POSV,iPart1) = LastPartPos(PART_POSV,iPart1) - dtColl * PartState(PART_VELV,iPart1)
! PartState(PART_POSV,iPart2) = LastPartPos(PART_POSV,iPart2) - dtColl * PartState(PART_VELV,iPart2)

! both particles can no longer be considered for collision
! PartCollisionProc%canPartCollide(iPart1Track) = .FALSE.
! PartCollisionProc%canPartCollide(iPart2Track) = .FALSE.

! de-apply periodic vectors (if there are some)
! IF (GEO%nPeriodicVectors.GT.0) THEN
!   iPartPeriodicShiftLoc = PartCollisionProc%iPartPeriodicShift(iCollision)
!   DO jPeriodic=1,GEO%nPeriodicVectors
!     PartState(PART_POSV,iPart2) = PartState(PART_POSV,iPart2) + REAL(MOD(iPartPeriodicShiftLoc, 3) - 1) * GEO%PeriodicVectors(1:3,jPeriodic)
!     iPartPeriodicShiftLoc = iPartPeriodicShiftLoc / 3
!   END DO
! END IF

END SUBROUTINE ComputeHardSphereCollision


#endif /*PARTICLES_COUPLING == 4*/

END MODULE MOD_Particle_Collision_Method
