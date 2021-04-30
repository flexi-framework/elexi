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
INTEGER                                :: nSGSVars                            ! number of variables in TurbPartState
INTEGER                                :: nSGSFilter                          ! number of cut-off modes in high-pass filter
CHARACTER(40)                          :: SGSModel                            ! specifying keyword for SGS model
REAL,ALLOCATABLE                       :: USGS(:,:,:,:,:)                     ! Unfiltered state
REAL,ALLOCATABLE                       :: USGSPart(:,:)                       ! Filtered state
!REAL,ALLOCATABLE                       :: kSGS(:,:,:,:,:)                     ! SGS kinetic energy
REAL,ALLOCATABLE                       :: kSGSPart(:)
REAL,ALLOCATABLE                       :: sigmaSGS(:)                         ! SGS kinetic energy standard deviation
REAL,ALLOCATABLE                       :: tauSGS(:)                           ! SGS time scale
REAL,ALLOCATABLE                       :: tauL (:,:)                          ! Parallel and perpendicular SGS time scale
REAL,ALLOCATABLE                       :: B_SGS(:,:,:)                        ! Diffusion matrix
REAL,ALLOCATABLE                       :: E_SGS(:,:,:)                        ! Exponential of the drift matrix
REAL,ALLOCATABLE                       :: G_SGS(:,:,:)                        ! Drift matrix
REAL,ALLOCATABLE                       :: W_SGS(:,:,:)                        ! Covariance matrix
REAL,ALLOCATABLE                       :: ElemVolN(:)                         ! ElemVol**(1./3.)/(PP_N+1)
REAL                                   :: randomVar
!! JIN SGS model
REAL,ALLOCATABLE                       :: deltaT_L(:)                         ! lagrangian time scale
REAL,ALLOCATABLE                       :: k_cf(:)                             ! cutoff wave number
REAL,ALLOCATABLE                       :: deltaT_E(:)                         ! Eulerian integral timescale (Fede & Simonin,2006)
REAL,ALLOCATABLE                       :: beta_Jin(:)                         ! ratio between lagrangian and eulerian time scale
REAL,ALLOCATABLE                       :: taup_Jin(:)                         ! Particle relaxation time
REAL,ALLOCATABLE                       :: St_Jin(:)                           ! Stokes number
REAL,ALLOCATABLE                       :: deltaT_L_p(:)                       ! SGS autocorrelation time seen by a particle
!! FUKAGATA SGS model
REAL,ALLOCATABLE                       :: taup_Fuka(:)                        ! Particle relaxation time
REAL,ALLOCATABLE                       :: Sc_Fuka(:)                          ! Schmidt number of the particles in fluid
!REAL,ALLOCATABLE                       :: ni_B(:)                             ! Brownian force
REAL,ALLOCATABLE                       :: gradUxPart(:,:)                     ! gradient in x-dir interpolate to the particle position
REAL,ALLOCATABLE                       :: gradUyPart(:,:)                     ! gradient in y-dir interpolate to the particle position
REAL,ALLOCATABLE                       :: gradUzPart(:,:)                     ! gradient in z-dir interpolate to the particle position
REAL,ALLOCATABLE                       :: T_f_p(:)                            ! Integral time scale of fluid velocity along a particle trajectory
REAL,ALLOCATABLE                       :: U_M(:)                              ! representative SGS velocity
REAL,ALLOCATABLE                       :: alpha_Fuka(:)                       ! ratio between dt and taup_Fuka 
REAL,ALLOCATABLE                       :: beta_Fuka(:)                        ! ratio between taup_Fuka and dt
REAL,ALLOCATABLE                       :: sigmaS_Fuka(:)                      ! standard deviation of particle velocity due to SGS velocity
REAL,ALLOCATABLE                       :: ni_SGS(:,:)                         ! SGS fluid velocity by Fukagata
REAL,ALLOCATABLE                       :: ni_Fuka(:,:)                        ! stochastic term sum of ni_B and ni_SGS
!===================================================================================================================================
END MODULE MOD_Particle_SGS_Vars
