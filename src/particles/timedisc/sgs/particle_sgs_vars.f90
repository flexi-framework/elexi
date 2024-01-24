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

!===================================================================================================================================
! Contains the variables for particle subgrid scale (SGS) models
!===================================================================================================================================
MODULE MOD_Particle_SGS_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                                :: ParticleSGSInitIsDone=.FALSE.
LOGICAL                                :: SGSinUse                            ! Flag for active SGS
!SGS method
!! BREUER SGS MODEL
INTEGER                                :: nSGSVars                            ! number of variables in TurbPartState
INTEGER                                :: nSGSFilter                          ! number of cut-off modes in high-pass filter
CHARACTER(40)                          :: SGSModel                            ! specifying keyword for SGS model
REAL,ALLOCATABLE                       :: USGS(:,:,:,:,:)                     ! Unfiltered state
REAL,ALLOCATABLE                       :: USGSPart(:,:)                       ! Filtered state
REAL,ALLOCATABLE                       :: ElemVolN(:)                         ! ElemVol**(1./3.)/(PP_N+1)
REAL                                   :: randomVar(3)
!! FUKAGATA SGS model
REAL,ALLOCATABLE                       :: gradp(:,:,:,:,:,:)                  ! Pressure gradient
REAL,ALLOCATABLE                       :: divTau(:,:,:,:,:)                   !
!===================================================================================================================================
END MODULE MOD_Particle_SGS_Vars
