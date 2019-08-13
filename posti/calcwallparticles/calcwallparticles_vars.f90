!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
!> Global variables for the calcwallparticles tool
!===================================================================================================================================
MODULE MOD_Posti_CalcWallParticles_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

! Variables for spatial averaging
INTEGER, ALLOCATABLE    :: SurfSideIDtoSurfAvgID(:)
!INTEGER, ALLOCATABLE    :: SurfAvgIDtoSurfSideID(:)
LOGICAL                 :: postiAvg                         !< Flag if an averaged statistics should be written
LOGICAL                 :: surfAvg                          !< Flag if an averaged surfState should be written
REAL,DIMENSION(3)       :: surfAvgDir                       !< Vector for averaging direction
LOGICAL                 :: surfRemap                        !< Flag if an surfState should be recalculated from impact files
REAL                    :: dt_remap                         !< Time step used for erosion remapping
INTEGER                 :: surfReflCount                    !< If passed, only particles with matching reflections will be counted

! User-defined parameters
CHARACTER(LEN=255)      :: AnalyzeString                    !< Method for integral calculation

END MODULE MOD_Posti_CalcWallParticles_Vars
