!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
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
! Contains the constant Advection Velocity Vector used for the linear scalar advection equation
!===================================================================================================================================
MODULE MOD_PICInterpolation_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                      :: FieldAtParticle(:,:) ! (PIC%maxParticleNumber,5) 2nd index: rho,u_x,u_y,u_z,e
#if USE_RW
REAL,ALLOCATABLE                      :: turbFieldAtParticle(:,:) ! (PIC%maxParticleNumber,2) 2nd index: k,epsilon
#endif
CHARACTER(LEN=256)                    :: InterpolationType    ! Type of Interpolation-Method
LOGICAL                               :: InterpolationElemLoop! Interpolate with outer iElem-loop (not for many Elems per proc!)
REAL                                  :: externalField(PP_nVar)! ext field is added to the maxwell-solver-field
LOGICAL                               :: DoInterpolation      ! Flag for interpolation
LOGICAL                               :: useBGField           ! Flag for BGField via h5-File
INTEGER                               :: NBG                  ! Polynomial degree of BG-Field
INTEGER                               :: BGType               ! Type of BG-Field (Electric,Magnetic,Both)
INTEGER                               :: BGDataSize           ! Type of BG-Field (Electric,Magnetic,Both)
REAL, ALLOCATABLE                     :: BGField(:,:,:,:,:)   ! BGField data
                                                              ! (1:x,0:NBG,0:NBG,0:NBG,1:PP_nElems)
REAL,ALLOCATABLE                      :: BGField_xGP(:)       ! Gauss point coordinates
REAL,ALLOCATABLE                      :: BGField_wGP(:)       ! GP integration weights
REAL,ALLOCATABLE                      :: BGField_wBary(:)     ! barycentric weights


CHARACTER(LEN=256)                    :: FileNameVariableExternalField ! filename containing the externanl field csv table
LOGICAL                               :: useVariableExternalField      ! use given external field. only for Bz variation in z
REAL,ALLOCATABLE                      :: VariableExternalField(:,:)    ! z - Pos , Bz
REAL                                  :: DeltaExternalField
INTEGER                               :: nIntPoints                    ! number of all interpolation points external field
!===================================================================================================================================
END MODULE MOD_PICInterpolation_Vars
