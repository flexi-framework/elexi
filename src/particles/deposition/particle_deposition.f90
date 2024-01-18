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
#include "flexi.h"
#include "particle.h"


!===================================================================================================================================
! Module containing the different deposition methods (NGP, linear (inter-cell) weighting, shape function
!===================================================================================================================================
MODULE MOD_Particle_Deposition
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitializeDeposition
  MODULE PROCEDURE InitializeDeposition
END INTERFACE

PUBLIC :: InitializeDeposition
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the deposition variables first
!===================================================================================================================================
SUBROUTINE InitializeDeposition()
! MODULES
USE MOD_Globals
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE DEPOSITION...'

LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DEPOSITION DONE!'

END SUBROUTINE InitializeDeposition

END MODULE MOD_Particle_Deposition
