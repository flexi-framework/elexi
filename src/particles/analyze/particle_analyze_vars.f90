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

!===================================================================================================================================
! Contains global variables used by the Analyze modules.
!===================================================================================================================================
MODULE MOD_Particle_Analyze_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: ParticleAnalyzeInitIsDone = .FALSE.
LOGICAL                       :: DoAnalyze                             ! perform analyze
LOGICAL                       :: CalcPartBalance                       ! Particle Power Balance - input and outflow energy of all
                                                                       ! particles
LOGICAL                       :: CalcEkin                              ! Compute the kinetic energy of each species
LOGICAL                       :: TrackParticlePosition                 ! track the particle movement
                                                                       ! stored in .csv format, debug only, no MPI
LOGICAL                       :: TrackParticleConvergence              ! track the final particle position, stored in .csv format
INTEGER                       :: nSpecAnalyze                          ! number of analyzed species 1 or nSpecies+1
INTEGER,ALLOCATABLE           :: nPartIn(:)                            ! Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartOut(:)                           ! Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartInTmp(:)                         ! Number of entry and leaving particles
REAL,ALLOCATABLE              :: PartEkinIn(:)                         ! energy and temperatur of input particle
REAL,ALLOCATABLE              :: PartEkinOut(:)                        ! energy and temperatur of input particle
REAL,ALLOCATABLE              :: PartEKinInTmp(:)                      ! energy and temperatur of input particle
LOGICAL                       :: printDiff
REAL                          :: printDiffTime
REAL                          :: printDiffVec(6)

REAL                          :: TimeSample
!===================================================================================================================================
END MODULE MOD_Particle_Analyze_Vars
