!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PUEPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
!==================================================================================================================================
!> Variables needed for the evaluation of the erosion points
!==================================================================================================================================
MODULE MOD_ErosionPoints_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL            :: ErosionPointsInitIsDone = .FALSE. !< mark wheter erosionpoints init routine is finished
LOGICAL            :: doParticleImpactTrack  = .FALSE.     !< mark whether erosionpoints should be evaluated during computation
LOGICAL            :: EP_onProc = .FALSE.     !< marks wheter current proc has EPs
LOGICAL            :: EP_fileExists = .FALSE. !< flag if EP file for analyze level has been created
INTEGER            :: EP_Buffersize           !< no. of time samples (size of EP_Data)
INTEGER            :: EP_MaxBuffersize        !< max. allowed no. of time samples
INTEGER            :: EP_Impacts              !< Impacts on current proc
INTEGER            :: offsetEP                !< offset for each proc in global EP list
INTEGER,ALLOCATABLE:: EP_ElemID(:)            !< mapping from EP->Elem (nEP)
REAL,ALLOCATABLE   :: EP_Data(:,:)            !< solution evaluated at EPs (nvar,nEP,nSamples)
INTEGER            :: EPDataSize = 15
!----------------------------------------------------------------------------------------------------------------------------------
! MPI Communicator for EPs
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER            :: myEPrank                !< rank within EP communicator
INTEGER            :: EP_COMM                 !< MPI EP communicator
INTEGER            :: nEP_Procs               !< number of procs with EPs

END MODULE MOD_erosionPoints_Vars
