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
MODULE MOD_Particle_Deposition_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: DepositionMethod
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  SUBROUTINE DepositionMethodInterface()
  END SUBROUTINE
END INTERFACE

PROCEDURE(DepositionMethodInterface),POINTER :: DepositionMethod             !< pointer defining the deposition method

INTEGER                           :: DepositionType                          ! Type of Deposition-Method
REAL,ALLOCATABLE                  :: PartSource(:,:,:,:,:)
REAL,ALLOCATABLE                  :: PartSource_tmp(:,:,:,:)
REAL,ALLOCATABLE                  :: Ut_src(:,:,:,:,:)

INTEGER                           :: nUniqueFEMNodes
REAL,ALLOCATABLE                  :: FEMNodeSource(:,:)
REAL,ALLOCATABLE                  :: CellVolWeightFac(:)
REAL,ALLOCATABLE                  :: NodeSource_tmp(:,:)

! MPI-SHM arrays
REAL,ALLOCPOINT,DIMENSION(:,:)    :: FEMNodeSource_Shared
#if USE_MPI
! integers to hold shared memory windows
INTEGER                           :: FEMNodeSource_Shared_Win
#endif /*USE_MPI*/

END MODULE MOD_Particle_Deposition_Vars
