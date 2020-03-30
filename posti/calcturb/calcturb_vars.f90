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
MODULE MOD_CalcTurb_Vars
! MODULES
USE ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
LOGICAL                 :: doConservativeDissipation

! Parameters for turbulence calculation
REAL,PARAMETER          :: betaStar = 0.09
REAL                    :: RestartEpsilon                   !< Constant epsilon in the domain for Turb=1,11
REAL                    :: RestartTKE                       !< Constant TKE in the domain for Turb=1,11

!CHARACTER(LEN=255)      :: MeshFile
!CHARACTER(LEN=255)      :: ProjectName
REAL                    :: OutputTime
INTEGER                 :: TurbMode                         !< Flag for mode of calcturb tool

INTEGER                 :: NCalc                            !< polynomial degree to perform the calculations
CHARACTER(LEN=255),ALLOCATABLE :: varnames_loc(:)

! Output format for state visualization
INTEGER                 :: OutputFormat
INTEGER,PARAMETER       :: OUTPUTFORMAT_NONE         = 0
INTEGER,PARAMETER       :: OUTPUTFORMAT_TECPLOT      = 1
INTEGER,PARAMETER       :: OUTPUTFORMAT_TECPLOTASCII = 2
INTEGER,PARAMETER       :: OUTPUTFORMAT_PARAVIEW     = 3

! Variables from time averaged file.
INTEGER                 :: nVarTurb
INTEGER                 :: nVar_HDF5,N_HDF5,nElems_HDF5
CHARACTER(LEN=255)      :: NodeType_HDF5

! Variables for fluid quantities
REAL,ALLOCATABLE        :: TKE        (  :,:,:,:)
REAL,ALLOCATABLE        :: UMean      (:,:,:,:,:)
REAL,ALLOCATABLE        :: gradUMeanx (:,:,:,:,:)           !< gradients in x-dir at degree N
REAL,ALLOCATABLE        :: gradUMeany (:,:,:,:,:)           !< gradients in y-dir at degree N
REAL,ALLOCATABLE        :: gradUMeanz (:,:,:,:,:)           !< gradients in z-dir at degree N
REAL,ALLOCATABLE        :: USolution  (:,:,:,:,:)

REAL,ALLOCATABLE        :: SMean      (:,:,:,:,:,:)         !< Shear tensor for mean flow
REAL,ALLOCATABLE        :: SFluc      (:,:,:,:,:,:)         !< Shear tensor for time-accurate flow
REAL,ALLOCATABLE        :: dFluc      (:,:,:,:)             !< Dilatation for time-accurate flow
REAL,ALLOCATABLE        :: Vorticity  (:,:,:,:,:)           !< Vorticity for time-accurate flow
REAL,ALLOCATABLE        :: EpsilonFin (:,:,:,:)             !< Final epsilon
REAL,ALLOCATABLE        :: EpsilonFluc(:,:,:,:)             !< Epsilon for time-accurate flow
REAL,ALLOCATABLE        :: EpsilonTmp (:,:,:,:)             !< Epsilon temporary storage
REAL,ALLOCATABLE        :: EpsilonSum (:,:,:,:)             !< Final epsilon
REAL,ALLOCATABLE        :: EpsilonMean(:,:,:,:)             !< Epsilon for mean flow

END MODULE MOD_CalcTurb_Vars
