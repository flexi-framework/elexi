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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_Tracking_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: TrackingMethod                ! Method that is used for tracking of particles
!                                          - REFMAPPING                  - tracking by mapping particle into reference element
!                                          - TRACING                     - tracking by tracing the particle path on curved sides
!                                          - TRIATRACKING                - tracking by using surface representation by triangles
REAL                                    :: tTracking                     ! Tracking time
REAL                                    :: tLocalization                 ! localization time
INTEGER                                 :: nTracks                       ! number of tracked particles
INTEGER                                 :: nCurrentParts                 ! current number of particles
LOGICAL                                 :: MeasureTrackTime              ! switch, if tracking time is measure
LOGICAL                                 :: FastPeriodic                  ! moves the particle along whole periodic vector,
                                                                         ! neglecting possible reflexions
LOGICAL                                 :: CartesianPeriodic             ! old periodic for refmapping and ALL bcs periocic
! LOGICAL                                 :: DoPeriodicCheck               ! flag for reference elem checking after periodic push
! LOGICAL                                 :: DoPeriodicFix                 ! flag for reference elem fixing after periodic push
LOGICAL                                 :: CountNbOfLostParts            ! logical, to count the lost particles
INTEGER                                 :: NbrOfLostParticles            ! Counter for lost particle
INTEGER                                 :: NbrOfLostParticlesTotal       ! Counter for lost particles across all processes
REAL,ALLOCATABLE                        :: Distance(:)                   ! list of distance between particle and element-origin
                                                                         ! to all elements in the same background element
INTEGER,ALLOCATABLE                     :: ListDistance(:)               ! the corresponding element id

TYPE tTrackingInfo
  INTEGER                               :: CurrElem
  INTEGER                               :: LocSide
  INTEGER                               :: GlobSide
  INTEGER                               :: flip
  INTEGER                               :: TriNum
  REAL                                  :: xi
  REAL                                  :: eta
  REAL                                  :: alpha
  REAL                                  :: PartTrajectory(1:3)
  REAL                                  :: LengthPartTrajectory
END TYPE

TYPE(tTrackingInfo)                     :: TrackInfo

#if CODE_ANALYZE
INTEGER                                 :: PartOut
INTEGER                                 :: MPIRankOut
#endif /*CODE_ANALYZE*/
!===================================================================================================================================


END MODULE MOD_Particle_Tracking_Vars
