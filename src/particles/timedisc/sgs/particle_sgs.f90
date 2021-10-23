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
#include "eos.h"
#include "particle.h"

!==================================================================================================================================
!> Module for different SGS models of the particle discretization
!==================================================================================================================================
MODULE MOD_Particle_SGS
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ParticleInitSGS
  MODULE PROCEDURE ParticleInitSGS
END INTERFACE

INTERFACE ParticleSGS
  MODULE PROCEDURE ParticleSGS
END INTERFACE

INTERFACE ParticleFinalizeSGS
  MODULE PROCEDURE ParticleFinalizeSGS
END INTERFACE

PUBLIC::ParticleInitSGS
PUBLIC::ParticleSGS
PUBLIC::ParticleFinalizeSGS
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
! Initialize the particle subgrid-scale model
!===================================================================================================================================
SUBROUTINE ParticleInitSGS()
! MODULES
USE MOD_Globals               ,ONLY: ABORT,UNIT_stdOut
USE MOD_Preproc               ,ONLY: PP_N,PP_NZ
USE MOD_Analyze_Vars          ,ONLY: ElemVol
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_SGS_Vars
USE MOD_Particle_Vars         ,ONLY: PDM,TurbPartState,TurbPt_temp
USE MOD_ReadInTools           ,ONLY: GETINT,GETSTR
#if USE_MPI
USE MOD_Globals               ,ONLY: MPIRoot
#endif /*USE_MPI*/
#if USE_RW
USE MOD_Particle_Randomwalk_Vars,ONLY: RWModel
#endif /*USE_RW*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ALLOCSTAT
!===================================================================================================================================

IF(ParticleSGSInitIsDone) RETURN

!SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SGS MODEL...'

! SGS model
SGSModel = TRIM(GETSTR('Part-SGSModel','none'))

! SGS and RW are not compatible. Go back and think about what you were trying to simulate.
#if USE_RW
IF (SGSModel.NE.'none'.AND.RWModel.NE.'none') &
  CALL abort(__STAMP__,'SGS and RW not compatible!')
#endif

SELECT CASE(SGSModel)
  CASE('Minier')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = GETINT('Part-SGSNFilter','2')
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:4       ,1:PDM%maxParticleNumber),      &
             gradp        (1,1:3     ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             gradpPart    (1:3       ,1:PDM%maxParticleNumber),      &
             gradUxPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart   (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
!    TurbPt_temp     = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)


  CASE('Sommerfeld')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = GETINT('Part-SGSNFilter','2')
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:4       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
!    TurbPt_temp     = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

CASE('Amiri')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = GETINT('Part-SGSNFilter','2')
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:4       ,1:PDM%maxParticleNumber),      &
             gradUxPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart   (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

CASE('Fukagata')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = GETINT('Part-SGSNFilter','2')
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:4       ,1:PDM%maxParticleNumber),      &
             gradUxPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart   (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

  CASE('Jin')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = GETINT('Part-SGSNFilter','2')
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:4       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

  CASE('Breuer')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = GETINT('Part-SGSNFilter','2')
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:4       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    TurbPt_temp     = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

  CASE('Breuer-Analytic')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = GETINT('Part-SGSNFilter','2')
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:4       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState = 0.
    TurbPt_temp   = 0.
    ElemVolN      = ElemVol**(1./3.)/(PP_N+1)

  CASE('none')
    ! Do nothing
    nSGSVars = 0
    SGSinUse = .FALSE.

  CASE DEFAULT
    CALL abort(__STAMP__, ' No valid particle subgrid scale (SGS) model given.')
END SELECT

ParticleSGSInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SGS MODEL DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE ParticleInitSGS


!==================================================================================================================================
!> Initialize all necessary information to perform SGS filtering
!==================================================================================================================================
SUBROUTINE InitSGSFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars
USE MOD_Interpolation_Vars    ,ONLY: Vdm_Leg,sVdm_Leg
USE MOD_Particle_SGS_Vars     ,ONLY: nSGSFilter
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iDeg
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' Init SGS filter...'

! Abort if Navier-Stokes filter is requested in addition to the SST filter
IF(FilterType.GT.0) CALL CollectiveStop(__STAMP__, 'SGS incompatible with Navier-Stokes filter!')

! Prepare modal cut-off filter (low pass)
ALLOCATE(FilterMat(0:PP_N,0:PP_N))
FilterMat = 0.

! Modal Filter, default to cut-off at PP_N-2
NFilter = PP_N - nSGSFilter
DO iDeg = 0,NFilter
  FilterMat(iDeg,iDeg) = 1.
END DO

! Assemble filter matrix in nodal space
FilterMat = MATMUL(MATMUL(Vdm_Leg,FilterMat),sVdm_Leg)

FilterInitIsDone = .TRUE.

END SUBROUTINE InitSGSFilter


!===================================================================================================================================
! SGS deconvolution
!===================================================================================================================================
SUBROUTINE ParticleSGS(dt,iStage) !,b_dt)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars                     ,ONLY: U,UPrim
USE MOD_Filter                      ,ONLY: Filter_Pointer
USE MOD_Filter_Vars                 ,ONLY: FilterMat
USE MOD_Mesh_Vars                   ,ONLY: offsetElem
USE MOD_Particle_Globals
USE MOD_Particle_SGS_Vars
USE MOD_Particle_Interpolation      ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars ,ONLY: FieldAtParticle
USE MOD_Particle_Vars               ,ONLY: TurbPartState,PDM,PEM,PartState
USE MOD_Particle_Vars               ,ONLY: Species,PartSpecies
USE MOD_Lifting_Vars                ,ONLY: gradUx, gradUy, gradUz
USE MOD_Lifting_BR1_gen             ,ONLY: Lifting_BR1_gen
!USE MOD_Particle_Vars               ,ONLY: TurbPt_Temp
USE MOD_TimeDisc_Vars               ,ONLY: nRKStages, RKC
USE MOD_Viscosity
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: dt
INTEGER,INTENT(IN),OPTIONAL   :: iStage
! REAL,INTENT(IN)               :: b_dt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: udiff(3)
REAL                          :: urel(3)    ! Normalized relative velocity vector
INTEGER                       :: ElemID, i, j, k, iPart
REAL                          :: Pt(1:3)
REAL                          :: RKdtFrac
REAL                          :: kSGSPart           ! SGS kinetic energy
REAL                          :: sigmaSGS           ! SGS kinetic energy standard deviation
REAL                          :: epsilonSGS         ! SGS dissipation energy by Shotorban
REAL                          :: taup               ! relative particle relaxation time
REAL                          :: alpha_SGS,beta_SGS
REAL                          :: rot_U(3)           ! rotor of velovity
REAL                          :: mu                 ! viscosity
! BREUER SGS MODEL
REAL,PARAMETER                :: C=1.
REAL,PARAMETER                :: betaSGS=1.
REAL                          :: tauSGS             ! SGS time scale
REAL                          :: tauL(2)            ! Parallel and perpendicular SGS time scale
REAL                          :: G_SGS(3,3)         ! Drift (analytic) / diffusion (num) matrix
REAL                          :: B_SGS(3,3)         ! Covariance (analytic) / exp. of the drift (num) matrix
! LOCAL VARIABLES - SOMMERFELD SGS MODEL
REAL                          :: C_eps_Shotorban=19./12  ! model constant of Shotorban (according to HEINZ)
REAL                          :: C_T=0.24           ! model constant for lagrangian time scale
REAL                          :: C_L=3              ! model constant for eulerian time scale
! LOCAL VARIABLES - MINIER & PEIRANO SGS MODEL
REAL                          :: C_0=19./12         ! model constant for lagrangian time scale MINIER (according to HEINZ)
REAL                          :: Lambda_factor      ! factor in order to calculate diffusion matrix
REAL                          :: trace_HR           ! trace matrix H*R
REAL                          :: trace_H            ! trace matrix H
REAL                          :: aux                ! auxiliary variable
REAL                          :: H_matrix(3,3)      !
REAL                          :: ReynoldsTensor(3,3)! Reynolds Stress Tensor
REAL                          :: H_R_matrix(3,3)    ! H_ij * R_ij
REAL                          :: pg(3)              ! Pressure gradient times dt
REAL                          :: second_drift_term(3) ! term relative at the relative displacement times gradientU
REAL                          :: fluctuationState(3) ! fluctuation velocity: U - <U>
REAL                          :: b_L(2)             ! Parallel and perpendicular Csanady factor
REAL                          :: G_matrix(3,3)      ! Drift matrix
REAL                          :: B_matrix(3,3)      ! Diffusion matrix
REAL                          :: D_matrix(3,3)      ! Symmetric matrix for Cholewski decomposition
! JIN MODEL
REAL                          :: deltaT_L           ! Lagrangian time scale
REAL                          :: deltaT_E           ! Eulerian integral timescale (Fede & Simonin,2006)
REAL                          :: St                 ! Stokes number
! FUKAGATA/AMIRI SGS MODEL
REAL                          :: ni_SGS(3)          ! stochastic term sum of ni_B and ni_SGS
REAL                          :: U_M                ! representative SGS velocity
REAL,PARAMETER                :: C_c= 1.            ! factor of Stokes-Cunningham
REAL,PARAMETER                :: C_Fuka= 0.93       ! constant model of Fukagata
REAL,PARAMETER                :: C_eps=0.08         ! = (2/(3*Ko))**(3/2), Ko: Kolmogorov constant = 1.5
REAL                          :: d_ij_tensor(3,3)   ! rate of deformation tensor
REAL                          :: beta1,beta2
REAL                          :: factor(3)          ! factor to take into account anisotropy (AMIRI MODEL)
REAL                          :: sigmaF(3)
! SOMMERFELD SGS MODEL
REAL                          :: deltaR(3),deltaR_mag ! relative displacement between particle and fluid element (vector/magnitude)
REAL                          :: T_L(3), T_L_f      ! Lagrangian time scale in the three directions (vector/scalar)
REAL                          :: R_L(3)             ! Lagrangian correlation in the three directions
REAL                          :: L_E                ! Integral length scales
REAL                          :: F_dr, G_dr, R_P(3) ! Correlaction vector for Sommerfeld
REAL                          :: diag_R_E(3)        ! Eulerian correlation (only main diagonal, vector form)
REAL                          :: sigmaF_square(3)    ! SGS kinetic energy standard deviation in the three directions
!===================================================================================================================================

SELECT CASE(SGSModel)

CASE('none')
!===================================================================================================================================
! DNS mode. Assume fully resolved turbulence
!===================================================================================================================================
! Do nothing


CASE('Minier')
!===================================================================================================================================
! Breuer (2017) model, first option
!> Breuer, M. and Hoppe, F., "Influence of a cost–efficient Langevin subgrid-scale model on the dispersed phase of large–eddy
!simulations of turbulent bubble–laden and particle–laden flows." International Journal of Multiphase Flow, 89 (2017): 23-44.
!===================================================================================================================================

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat)

! Obtain the high-pass filtered velocity field
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

! Interpolate gradient to particle position
CALL InterpolateFieldToParticle(3,gradUx(1:3,:,:,:,:),3,gradUxPart)
CALL InterpolateFieldToParticle(3,gradUy(1:3,:,:,:,:),3,gradUyPart)
CALL InterpolateFieldToParticle(3,gradUz(1:3,:,:,:,:),3,gradUzPart)

! Calculate pressure gradient
CALL Lifting_BR1_gen(1,1,UPrim(PRES:PRES,:,:,:,:),gradp(:,1,:,:,:,:),gradp(:,2,:,:,:,:),gradp(:,3,:,:,:,:))

! Interpolate pressur gradient to particel position
CALL InterpolateFieldToParticle(3,gradp(1,1:3,:,:,:,:),3,gradpPart)

DO iPart = 1,PDM%ParticleVecLength

  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM((USGSPart(MOMV,iPart)/FieldAtParticle(DENS,iPart))**2.)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Calculate dissipation energy according to  Babak Shotorban & Farzad Mashayek (2006) A stochastic model for particle
  ! motion in large-eddy simulation, Journal of Turbulence, 7, N18, DOI: 10.1080/14685240600595685
  epsilonSGS = C_eps_Shotorban/ElemVolN(ElemID) * (kSGSPart**(3./2))

  ! Relative velocity
  udiff(1:3) = PartState(PART_VELV,iPart) - FieldAtParticle(VELV,iPart) !+ TurbPartState(1:3,iPart) ! only mean velocity

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/SQRT(SUM(udiff**2))
  ELSE
    urel = 0.
  END IF

  ! Calculate fluctuation velocity
  fluctuationState(1:3) = TurbPartState(1:3,iPart) - FieldAtParticle(VELV,iPart)

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(kSGSPart)) THEN ! 1st IF

    ! We ASSUME that these are the CORRECT matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN   ! 2ND IF
      Pt(1:3) = 0.
      second_drift_term(1:3) = 0.
    ELSE     ! 2ND ELSE

      ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j)
      second_drift_term(1) = udiff(1) * gradUxPart(1,iPart) + udiff(2) * gradUyPart(1,iPart) + udiff(3) * gradUzPart(1,iPart)
      second_drift_term(2) = udiff(1) * gradUxPart(2,iPart) + udiff(2) * gradUyPart(2,iPart) + udiff(3) * gradUzPart(2,iPart)
      second_drift_term(3) = udiff(1) * gradUxPart(3,iPart) + udiff(2) * gradUyPart(3,iPart) + udiff(3) * gradUzPart(3,iPart)

      ! parallel Csanady factor
      b_L(1) = betaSGS*SQRT(SUM(udiff**2))
      ! perpendicular Csanady factor
      b_L(2) = betaSGS*2*SQRT(SUM(udiff**2))

      ! Calculate auxliary matrix H_ij, drift matrix G_ij
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            H_matrix(i,j) = b_L(2)  + (b_L(1)  - b_L(2)) *urel(i)*urel(j)
          ELSE
            H_matrix(i,j) =           (b_L(1)  - b_L(2)) *urel(i)*urel(j)
          END IF
          G_matrix(i,j) = -(0.5 + 3*C_0/4) * (C_eps_Shotorban/(SQRT(2./3))*1./ElemVolN(ElemID)) * H_matrix(i,j)
        END DO
      END DO

      ! Calculate the Reynolds stress tensor
      DO i = 1,3
        DO j = 1,3
          ReynoldsTensor(i,j) = (USGSPart(i+1,iPart)/FieldAtParticle(DENS,iPart))*(USGSPart(j+1,iPart)/FieldAtParticle(DENS,iPart))
        END DO
      END DO

      !Calculate factor \lambda
      H_R_matrix = MATMUL(H_matrix,ReynoldsTensor)

      trace_HR = 0.
      trace_H  = 0.

      DO i = 1,3
        trace_HR = trace_HR + H_R_matrix(i,i)
        trace_H  = trace_H  + H_matrix(i,i)
      END DO

      Lambda_factor = trace_HR/trace_H

      ! Calculate simmetrix matrix D_ij for Cholewski decomposition
      DO i = 1,3
        DO j = 1,3
          D_matrix(i,j) = ( C_0/SQRT(8./27) + SQRT(3./2) ) * Lambda_factor * H_matrix(i,j)
        END DO
      END DO
      D_matrix(:,:) = C_eps_Shotorban/ElemVolN(ElemID) * D_matrix(:,:)

      ! Cholewski decomposition to get diffusion matrix B_matrix
      B_matrix(:,:) = 0.

      B_matrix(1,1) = sqrt( D_matrix(1,1) )       ! element (1,1)
      DO i = 2 , 3
        B_matrix(i,1) = D_matrix(1,i) / B_matrix(1,1)  ! first column
      END DO

      DO j = 2 , 3
        DO i = 1 , 3
          IF( i .LT. j ) THEN
            B_matrix(i,j) = 0.0           ! upper triangolar elements
          ELSE
            IF( i .EQ. j) THEN
              aux = 0.0
              DO k = 1 , i -1
                aux = aux + B_matrix(i,k)**2
              END DO
              B_matrix(i,j) = sqrt( D_matrix(i,i) - aux ) ! diagonal elements
            ELSE ! (i .GT. j)
              aux = 0.0
              DO k = 1 , i - 1
                aux = aux + B_matrix(j,k) * B_matrix(i,k)
              END DO
              B_matrix(i,j) = ( D_matrix(j,i) - aux ) / B_matrix(j,j)    ! in our case only element B(3,2)
            END IF
          END IF
        END DO
      END DO

      ! EULER
      ! Sum up turbulent contributions
      Pt(1:3) = 0.
      DO  j=1,3
        Pt(1:3) = Pt(1:3) + G_matrix(1:3,j)*fluctuationState(j)*dt + B_matrix(1:3,j)*SQRT(dt)*RandNormal()
      END DO
    END IF !2ND END IF

  ! Valid SGS turbulent kinetic energy
  ELSE  ! 1ST ELSE

    IF(ALMOSTZERO(MAXVAL(udiff)))THEN  ! 3RD IF

      ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j)
      second_drift_term(1:3) = 0.

      ! parallel Csanady factor
      b_L(1) = 1
      ! perpendicular Csanady factor
      b_L(2) = 1

      ! Calculate auxliary matrix H_ij, drift matrix G_ij
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            H_matrix(i,j) = b_L(2)
          ELSE
            H_matrix(i,j) = 0
          END IF
          G_matrix(i,j) = -(0.5 + 3*C_0/4) * epsilonSGS/kSGSPart * H_matrix(i,j)
        END DO
      END DO

      ! Calculate the Reynolds stress tensor
      DO i = 1,3
        DO j = 1,3
          ReynoldsTensor(i,j) = (USGSPart(i+1,iPart)/FieldAtParticle(DENS,iPart))*(USGSPart(j+1,iPart)/FieldAtParticle(DENS,iPart))
        END DO
      END DO

      !Calculate factor \lambda
      H_R_matrix = MATMUL(H_matrix,ReynoldsTensor)

      trace_HR = 0.
      trace_H  = 0.

      DO i = 1,3
        trace_HR = trace_HR + H_R_matrix(i,i)
        trace_H  = trace_H  + H_matrix(i,i)
      END DO

      Lambda_factor = (3 * trace_HR) / (2 * kSGSPart * trace_H)

      ! Calculate simmetrix matrix D_ij for Cholewski decomposition
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            D_matrix(i,j) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j) -1)/3
          ELSE
            D_matrix(i,j) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j)   )/3
          END IF
        END DO
      END DO
      D_matrix(:,:) = epsilonSGS * D_matrix(:,:)

      ! Cholewski decomposition to get diffusion matrix B_matrix
      B_matrix(:,:) = 0.

      B_matrix(1,1) = sqrt( D_matrix(1,1) )    ! element (1,1)

      DO i = 2 , 3
        B_matrix(i,1) = D_matrix(1,i) / B_matrix(1,1) ! first column
      END DO

      DO j = 2 , 3
        DO i = 1 , 3
          IF( i .LT. j ) THEN
            B_matrix(i,j) = 0.0    ! upper triangolar elements
          ELSE
            IF( i .EQ. j) THEN
              aux = 0.0
              DO k = 1 , i -1
                aux = aux + B_matrix(i,k)**2
              END DO
              B_matrix(i,j) = sqrt( D_matrix(i,i) - aux ) ! diagonal elements
            ELSE ! (i .GT. j)
              aux = 0.0
              DO k = 1 , i - 1
                aux = aux + B_matrix(j,k) * B_matrix(i,k)
              END DO
              B_matrix(i,j) = ( D_matrix(j,i) - aux ) / B_matrix(j,j)    ! in our case only element B(3,2)
            END IF
          END IF
        END DO
      END DO

    ELSE ! with valid KSGS and valid UDIFF      ! 3RD ELSE

      ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j)
      second_drift_term(1) = udiff(1) * gradUxPart(1,iPart) + udiff(2) * gradUyPart(1,iPart) + udiff(3) * gradUzPart(1,iPart)
      second_drift_term(2) = udiff(1) * gradUxPart(2,iPart) + udiff(2) * gradUyPart(2,iPart) + udiff(3) * gradUzPart(2,iPart)
      second_drift_term(3) = udiff(1) * gradUxPart(3,iPart) + udiff(2) * gradUyPart(3,iPart) + udiff(3) * gradUzPart(3,iPart)

      ! parallel Csanady factor
      b_L(1) = (SQRT(1+  betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))
      ! perpendicular Csanady factor
      b_L(2) = (SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))

      ! Calculate auxliary matrix H_ij, drift matrix G_ij
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            H_matrix(i,j) = b_L(2)  + (b_L(1)  - b_L(2)) *urel(i)*urel(j)
          ELSE
            H_matrix(i,j) =           (b_L(1)  - b_L(2)) *urel(i)*urel(j)
          END IF
          G_matrix(i,j) = -(0.5 + 3*C_0/4) * epsilonSGS/kSGSPart * H_matrix(i,j)
        END DO
      END DO

      ! Calculate the Reynolds stress tensor
      DO i = 1,3
        DO j = 1,3
          ReynoldsTensor(i,j) = (USGSPart(i+1,iPart)/FieldAtParticle(DENS,iPart))*(USGSPart(j+1,iPart)/FieldAtParticle(DENS,iPart))
        END DO
      END DO

      !Calculate factor \lambda
      H_R_matrix = MATMUL(H_matrix,ReynoldsTensor)

      trace_HR = 0.
      trace_H  = 0.

      DO i = 1,3
        trace_HR = trace_HR + H_R_matrix(i,i)
        trace_H  = trace_H  + H_matrix(i,i)
      END DO

      Lambda_factor = (3 * trace_HR) / (2 * kSGSPart * trace_H)

      ! Calculate simmetrix matrix D_ij for Cholewski decomposition
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            D_matrix(i,j) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j) -1)/3
          ELSE
            D_matrix(i,j) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j)   )/3
          END IF
        END DO
      END DO
      D_matrix(:,:) = epsilonSGS * D_matrix(:,:)

      ! Cholewski decomposition to get diffusion matrix B_matrix
      B_matrix(:,:) = 0.

      B_matrix(1,1) = sqrt( D_matrix(1,1) )     ! element (1,1)
      DO i = 2 , 3
        B_matrix(i,1) = D_matrix(1,i) / B_matrix(1,1)  ! first column
      END DO

      DO j = 2 , 3
        DO i = 1 , 3
          IF( i .LT. j ) THEN
            B_matrix(i,j) = 0.0    ! upper triangolar elements
          ELSE
            IF( i .EQ. j) THEN
              aux = 0.0
              DO k = 1 , i -1
                aux = aux + B_matrix(i,k)**2
              END DO
              B_matrix(i,j) = SQRT( D_matrix(i,i) - aux ) ! diagonal elements
            ELSE ! (i .GT. j)
              aux = 0.0
              DO k = 1 , i - 1
                aux = aux + B_matrix(j,k) * B_matrix(i,k)
              END DO
              B_matrix(i,j) = ( D_matrix(j,i) - aux ) / B_matrix(j,j)    ! in our case only element B(3,2)
            END IF
          END IF
        END DO
      END DO

    END IF ! 3RD END IF

    ! EULER
    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3)  + G_matrix(1:3,j)*fluctuationState(j)*dt + B_matrix(1:3,j)*SQRT(dt)*RandNormal()
    END DO
  END IF ! 1ST END IF

  ! FULL FLUID VELOCITY SEEN BY THE PARTICLES
  ! 1/\rho_f * ( \partial <P> / \partial x_i ) * dt
  pg(1:3)  = gradpPart(1:3,iPart)/FieldAtParticle(DENS,iPart)*dt

  !  (U_p-U)*(\partial U_i / \partial x_j) * dt
  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) - pg(1:3) + second_drift_term(1:3)*dt + Pt(1:3)

END DO  ! Minier & Peirano


CASE('Sommerfeld')
!===================================================================================================================================
! Martin Sommerfeld 2001
!> Validation of a stochastic Lagrangian modelling approach for inter-particle collisions in homogeneous isotropic turbulence.
!> Institut fur Verfahrenstechnik, Fachbereich Ingenieurwissenschaften, Martin-Luther-Universitat at Halle-Wittenberg,
!> D-06099 Halle (Saale), Germany. Received 30 March 2000; received in revised form 3 May 2001

!> Model implemented from: EVALUATION OF LAGRANGIAN PARTICLE DISPERSION MODELS IN TURBULENT FLOWS (S. Laı́n∗ and C.A. Grillo)
!===================================================================================================================================

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat)

! Obtain the high-pass filtered velocity field
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(DENS,iPart))**2.)

  ! Deviation standard of fluctuation velocity
  sigmaSGS = SQRT(2./3.*kSGSPart)

  ! Deviation standard per direction
  ! considering sigmaF_vector already to the power of 2
  sigmaF_square(1:3) = (USGSPart(2:4,iPart)/FieldAtParticle(DENS,iPart))**2.

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Calculate dissipation energy according to  Babak Shotorban & Farzad Mashayek (2006) A stochastic model for particle
  ! motion in large-eddy simulation, Journal of Turbulence, 7, N18, DOI: 10.1080/14685240600595685
  epsilonSGS = C_eps_Shotorban/ElemVolN(ElemID) * (kSGSPart**(3./2))

  ! Relative velocity
  udiff(1:3) =  (FieldAtParticle(2:4,iPart) + TurbPartState(1:3,iPart)) - PartState(4:6,iPart)

  ! Relative displcement between particle and fluid element
  deltaR(1:3) = udiff(1:3)*dt   ! relative displacement in the three directions
  deltaR_mag = VECNORM(deltaR(1:3)) ! |U-U_p|*dt (scalar value)

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS)) THEN  !1st IF
    ! We ASSUME that these are the CORRECT matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN  !2nd IF
      Pt(1:3) = 0 ! in this case we have 0/0
    ELSE ! 2nd  ELSE IF

      ! Lagrangian correlation
      R_L(1:3) = 1.

      ! Integral length scales L_E when sigmaSGS goes to zero
      L_E = C_L**2./C_eps_Shotorban * SQRT(8./27) * ElemVolN(ElemID)

      ! Function g(dr) and f(dr)
      F_dr = EXP(-deltaR_mag / L_E)
      G_dr = (1 - deltaR_mag / (2*L_E) ) *  EXP(-deltaR_mag/L_E)

      ! I'm considering only the main diagonal of the tensor (treat it as a  vector)
      diag_R_E(1:3) = G_dr + (F_dr - G_dr)*(deltaR(1:3)*deltaR(1:3))/(deltaR_mag**2)

      !Correlaction function: vector R_P
      R_P(1:3) = R_L(1:3) * diag_R_E(1:3)

      ! Calculate the fluctuation in the three directions
      DO i = 1,3
        Pt(i) = R_P(i)*TurbPartState(i,iPart) + SQRT(sigmaF_square(i))*SQRT(1-R_P(i)**2) * RandNormal()  !!!* SQRT(dt)
      END DO

    END IF  !2nd END IF

  ! Valid SGS turbulent kinetic energy
  ELSE    !1st  ELSE IF

    IF(ALMOSTZERO(MAXVAL(udiff)))THEN   ! 3rd IF
      Pt(1:3) = 0 ! in this case we have 0/0
    ELSE ! 3rd  ELSE IF

      ! Lagrangian time scale in the appropriate directions
      T_L(1:3) = C_T * sigmaF_square(1:3)/epsilonSGS

      ! Lagrangian correlation
      R_L(1:3) = EXP(-dt/T_L(1:3))

      ! Integral length scales L_E
      T_L_f = C_T *  (sigmaSGS**2) / epsilonSGS
      L_E = C_L * T_L_f * sigmaSGS

      ! Function g(dr) and f(dr)
      F_dr = EXP(-deltaR_mag / L_E)
      G_dr = (1 - deltaR_mag / (2*L_E) ) *  EXP(-deltaR_mag/L_E)

      diag_R_E(1:3) = G_dr + (F_dr - G_dr)*(deltaR(1:3)*deltaR(1:3))/(deltaR_mag**2)

      !Correlaction function: vector R_P
      R_P(1:3) = R_L(1:3) * diag_R_E(1:3)

      ! Calculate the fluctuation in the three directions
      DO i = 1,3
        Pt(i) = R_P(i)*TurbPartState(i,iPart) + SQRT(sigmaF_square(i))*SQRT(1-R_P(i)**2)* RandNormal()
      END DO
    END IF ! 3rd END IF
  END IF ! 1ST END IF

  TurbPartState(1:3,iPart) =  Pt(1:3)
END DO  ! For Sommerfeld


CASE('Amiri')
!===================================================================================================================================
! A. Elhami Amiri , S. Kazemzadeh Hannani & F. Mashayek (2006)
!>  Large-Eddy Simulation of Heavy-Particle Transport in Turbulent Channel Flow, Numerical Heat Transfer,
!> Part B: Fundamentals, 50:4, 285-313, DOI: 10.1080/10407790600859577

!===================================================================================================================================

! Calculate the dyn. viscosity
mu=VISCOSITY_PRIM(FieldAtParticle)

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat) ! DOUBLE FILTERING

! Obtain the high-pass filtered velocity field
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

! Interpolate gradient to the particle position
CALL InterpolateFieldToParticle(3,gradUx(1:3,:,:,:,:),3,gradUxPart)
CALL InterpolateFieldToParticle(3,gradUy(1:3,:,:,:,:),3,gradUyPart)
CALL InterpolateFieldToParticle(3,gradUz(1:3,:,:,:,:),3,gradUzPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(kSGSPart)) THEN

    !in the limit of sigmaF goes to zero the time sigmaS goes to zero
    ni_SGS(1:3) = 0.

  ! Valid SGS turbulent kinetic energy
  ELSE

    ! Particle relaxation time
    taup = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
    & * 1./(18*mu)

    ! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
    ! Local vorticity
    rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
    rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
    rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy

    ! Calculate the factor f1, f2, f3 for the anisotropy
    factor(1) = SQRT((USGSPart(2,iPart)**2)/(FieldAtParticle(1,iPart))) / SQRT( 2./3 * kSGSPart)   ! f1
    factor(2) = SQRT((USGSPart(3,iPart)**2)/(FieldAtParticle(1,iPart))) / SQRT( 2./3 * kSGSPart)   ! f2
    factor(3) = SQRT(3 - factor(1)**2 - factor(2)**2)            ! f3

    !alpha and beta Amiri
    alpha_SGS =  dt/taup
    beta_SGS  = taup*VECNORM(rot_U)
    beta1 = 1+beta_SGS
    beta2 = 1-beta_SGS

    ! sigmaS_Amiri is the increase of standard deviation of particle velocity due to SGS velocity during time dt
    ! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity
    ! due to SGS velocity in time dt
    ! U_square = (u'')^2
    sigmaF(1:3) = SQRT(factor(1:3)**2 * 2./3 * kSGSPart) * SQRT( 1/beta1 * (1-EXP(-alpha_SGS*beta1))  &
    & -1/beta2 * EXP(-2*alpha_SGS) * (1-EXP(alpha_SGS*beta2))  )

    ! SGS fluid velocity by Amiri
    DO j=1,3
      ni_SGS(j) =  sigmaF(j)/dt * RandNormal()
    END DO

  END IF

  TurbPartState(1:3,iPart) =  ni_SGS(1:3)

END DO ! AMIRI SGS MODEL


CASE('Fukagata')
!===================================================================================================================================
! K. Fukagata (2004)
!> Dynamics of Brownian particles in a turbulent channel flow, K. Fukagata, S. Zahrai, F. H. Bark
!> Heat and Mass Transfer 40 (2004) 715–726. DOI 10.1007/s00231-003-0462-8
!===================================================================================================================================

! Calculate the dyn. viscosity
mu=VISCOSITY_PRIM(FieldAtParticle)

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat) ! DOUBLE FILTERING

! Obtain the high-pass filtered velocity field
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

! Interpolate gradient to the particle position
CALL InterpolateFieldToParticle(3,gradUx(1:3,:,:,:,:),3,gradUxPart)
CALL InterpolateFieldToParticle(3,gradUy(1:3,:,:,:,:),3,gradUyPart)
CALL InterpolateFieldToParticle(3,gradUz(1:3,:,:,:,:),3,gradUzPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Particle relaxation time
  taup = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 * 1./(18*mu)

  ! Schmidt number of the particles in fluid
!  Sc_Fuka = Species(PartSpecies(iPart))%MassIC*nu / ( taup * kSGSPart * FieldAtParticle(TEMP,iPart) )

  ! Brownian force
!  ni_B = RandNormal()*SQRT(2*nu/(Sc_Fuka * taup**2 * dt))

  ! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
  ! Local vorticity
  rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
  rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
  rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy

  ! Calculation of deformation tensor d_ij
  d_ij_tensor(1,:)= 0.5*(/2*gradUxPart(1,iPart), gradUyPart(1,iPart)+gradUxPart(2,iPart), gradUzPart(1,iPart)+gradUxPart(3,iPart)/)
  d_ij_tensor(2,:)= 0.5*(/gradUyPart(1,iPart)+gradUxPart(2,iPart), 2*gradUyPart(2,iPart), gradUzPart(2,iPart)+gradUyPart(3,iPart)/)
  d_ij_tensor(3,:)= 0.5*(/gradUzPart(1,iPart)+gradUxPart(3,iPart), gradUzPart(2,iPart)+gradUyPart(3,iPart), 2*gradUzPart(3,iPart)/)

  ! Calculation of norm of tensor d_ij according to F. Gurniki, K. Fukagata, S. Zahrai and F. H. Bark (2000)
  ! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity due to SGS velocity in time dt
  U_M = 2/(3**0.5) * C_Fuka * C_eps * ElemVolN(ElemID) * SQRT(SUM(d_ij_tensor(:,:)**2))

  !alpha and beta Fukagata
  alpha_SGS =  dt/taup
  ! taup  / T_f_p
  beta_SGS = taup*VECNORM(rot_U)
  beta1 = 1+beta_SGS
  beta2 = 1-beta_SGS

  !sigmaSGS is the increase of standard deviation of particle velocity due to SGS velocity during time dt
  sigmaSGS = U_M * SQRT( 1/beta1 * (1-EXP(-alpha_SGS*beta1)) - 1./beta2 * EXP(-2*alpha_SGS) * (1-EXP(alpha_SGS*beta2)))

  ! SGS fluid velocity by Fukagata
  ! now I'm using only the SGS fluid velocity (NO brownian force)
  ! TODO: check if that is correct
  DO j=1,3
    ni_SGS(j) =  sigmaSGS/dt * RandNormal() !+ ni_B
  END DO

  ! TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)
  TurbPartState(1:3,iPart) =  ni_SGS(1:3)

END DO ! Fukagata


CASE('Jin')
!===================================================================================================================================
! Jin (2010)
!> Guodong Jin a , Guo-Wei He a,*, Lian-Ping Wang b, Jian Zhang a Subgrid scale fluid velocity timescales seen by inertial
!  particles in large-eddy simulation of particle-laden turbulence. International Journal of Multiphase Flow 36 (2010) 432–437
!===================================================================================================================================

! Calculate the dyn. viscosity
mu=VISCOSITY_PRIM(FieldAtParticle)

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat) ! DOUBLE FILTERING

! Obtain the high-pass filtered velocity field
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Time scale of SGS scales = U_rms in case of isotropic turbolence
  sigmaSGS = SQRT(2./3*kSGSPart)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS)) THEN

    ! In the limit of sigmaSGS goes to zero the time scale goes to infinity
    Pt(1:3) = 0.

  ! Valid SGS turbulent kinetic energy
  ELSE

    ! Lagrangian time scale
    deltaT_L   = C*ElemVolN(ElemID)/sigmaSGS

    !  SGS Eulerian integral timescale according to Fede and Simonin, 2006
    deltaT_E = 3*PI/10*1./(PI/ElemVolN(ElemID)*sigmaSGS)

    ! ratio between lagrangian and eulerian time scale
    beta_SGS = deltaT_L/deltaT_E

    ! Particle relaxation time
    taup = Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 * 1./(18*mu)

    ! Stokes number
    St = taup/deltaT_E

    ! SGS autocorrelation time seen by a particle =  deltaT_L
    deltaT_L =  deltaT_L/beta_SGS * (1 - ((1-beta_SGS) / (1+St)**(0.4*(1+0.01*St))))

    ! EULER
    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3) - 1./deltaT_L*TurbPartState(j,iPart)*dt + SQRT( (4*kSGSPart*dt)/(deltaT_L*3) ) * RandNormal()
    END DO

  END IF

  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)

END DO ! Jin


CASE('Breuer')
!===================================================================================================================================
! Breuer (2017) model, first option
!> Breuer, M. and Hoppe, F., "Influence of a cost–efficient Langevin subgrid-scale model on the dispersed phase of large–eddy
!simulations of turbulent bubble–laden and particle–laden flows." International Journal of Multiphase Flow, 89 (2017): 23-44.
!===================================================================================================================================

! Time integration using Euler-Maruyama scheme
IF (iStage.NE.1) RETURN

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat)

! Obtain the high-pass filtered velocity field
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Time scale of SGS scales
  sigmaSGS = SQRT(2./3.*kSGSPart)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Relative velocity
  udiff(1:3) = PartState(PART_VELV,iPart) - (FieldAtParticle(VELV,iPart) + TurbPartState(1:3,iPart))

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/SQRT(SUM(udiff**2))
  ELSE
    urel = 0.
  END IF

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS)) THEN
    ! We ASSUME that these are the CORRECT matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      RETURN
    ELSE
      ! parallel
      tauL(1) = C*ElemVolN(ElemID)/(betaSGS*SQRT(SUM(udiff**2)))
      ! perpendicular
      tauL(2) = C*ElemVolN(ElemID)/(betaSGS*2*SQRT(SUM(udiff**2)))

      ! Calculate drift and diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            G_SGS(i,j) = 1/tauL(2)  + (1/tauL(1) - 1/tauL(2)) *urel(i)*urel(j)
          ELSE
            G_SGS(i,j) =              (1/tauL(1) - 1/tauL(2)) *urel(i)*urel(j)
          END IF
        END DO
      END DO
      ! as sigma is ALMOSTZERO
      B_SGS(:,:) = 0.

      ! EULER
      ! Sum up turbulent contributions
      Pt(1:3) = 0.
      DO j = 1,3
        Pt(1:3) = Pt(1:3) - G_SGS(1:3,j)*TurbPartState(j,iPart)*dt
      END DO
    END IF

  ! Valid SGS turbulent kinetic energy
  ELSE

    tauSGS   = C*ElemVolN(ElemID)/sigmaSGS

    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      ! parallel
      tauL(1) = tauSGS
      ! perpendicular
      tauL(2) = tauSGS

      ! Calculate drift and diffusion matrix
      G_SGS(:,:)=0.
      B_SGS(:,:)=0.
      DO i = 1,3
        G_SGS(i,i) = 1/     tauL(2)
        B_SGS(i,i) = 1/SQRT(tauL(2))
      END DO
      B_SGS(:,:) = SQRT(2*sigmaSGS**2)*B_SGS(:,:)

    ELSE
      ! parallel
      tauL(1) = tauSGS/(SQRT(1+  betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))
      ! perpendicular
      tauL(2) = tauSGS/(SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))

      ! Calculate drift and diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            G_SGS(i,j) = 1/     tauL(2)  + (1/     tauL(1)  - 1/     tauL(2)) *urel(i)*urel(j)
            B_SGS(i,j) = 1/SQRT(tauL(2)) + (1/SQRT(tauL(1)) - 1/SQRT(tauL(2)))*urel(i)*urel(j)
          ELSE
            G_SGS(i,j) =                   (1/     tauL(1)  - 1/     tauL(2)) *urel(i)*urel(j)
            B_SGS(i,j) =                   (1/SQRT(tauL(1)) - 1/SQRT(tauL(2)))*urel(i)*urel(j)
          END IF
        END DO
      END DO

      B_SGS(:,:) = SQRT(2*sigmaSGS**2)*B_SGS(:,:)
    END IF

    ! EULER
    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3) - G_SGS(1:3,j)*TurbPartState(j,iPart)*dt + B_SGS(1:3,j)*RandNormal()*SQRT(dt)
    END DO
  END IF

  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)

  ! RUNGE-KUTTA
  ! Sum up turbulent contributions
!  Pt(1:3)=0.
!  DO j=1,3
!    Pt(1:3) = Pt(1:3) - G_SGS(1:3,j,iPart)*TurbPartState(j,iPart) + B_SGS(1:3,j,iPart)*TurbPartState(j+3,iPart)/SQRT(dt)
!  END DO
!  !--> First RK stage
!  IF (iStage.EQ.1) THEN
!    TurbPt_temp  (1:3,iPart) = Pt
!    TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + TurbPt_temp(1:3,iPart)*b_dt
!  !--> Later RK stage
!  ELSE
!    TurbPt_temp  (1:3,iPart) = Pt(1:3) - RKA(iStage)    * TurbPt_temp(1:3,iPart)
!    TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + TurbPt_temp(1,iPart)*b_dt
!  END IF


!>> In near-wall turbulence, the SGS fluctuations can readily exceed 10*U_resolved
!  ! Sanity check. Use 10*U as arbitrary threshold
!  IF (ANY(ABS(TurbPartState(1:3,iPart)).GT.10.*MAXVAL(ABS(PartState(4:6,iPart))))) THEN
!    IPWRITE(*,*) 'Obtained SGS velocity of',SQRT(SUM(TurbPartState(1:3,iPart)**2)), &
!                 'while resolved velocity is',SQRT(SUM(PartState(4:6,iPart)**2)),   &
!                 'for',iPart
!    TurbPartState(:,iPart) = 0.
!  END IF

  END DO  ! For Breuer

CASE('Breuer-Analytic')
!===================================================================================================================================
! Breuer (2017) model, second option
!> Breuer, M. and Hoppe, F., "Influence of a cost–efficient Langevin subgrid-scale model on the dispersed phase of large–eddy
!simulations of turbulent bubble–laden and particle–laden flows." International Journal of Multiphase Flow, 89 (2017): 23-44.
!===================================================================================================================================

! Time integration in first RK stage (p. 26)
IF(PRESENT(iStage))THEN
  IF (iStage.EQ.1) THEN
    RKdtFrac      = RKC(2)*dt
    DO j=1,3
      randomVar(j)=RandNormal()
    END DO
  ELSE
    IF (iStage.NE.nRKStages) THEN
      RKdtFrac      = (RKC(iStage+1)-RKC(iStage))*dt
    ELSE
      RKdtFrac      = (1.-RKC(nRKStages))*dt
    END IF
  END IF
ELSE
  RKdtFrac = dt
  DO j=1,3
    randomVar(j)=RandNormal()
  END DO
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat)

! Obtain the high-pass filtered velocity field
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Time scale of SGS scales
  sigmaSGS = SQRT(2./3*kSGSPart)

  ! Relative velocity
  udiff(1:3) = PartState(4:6,iPart) - (FieldAtParticle(2:4,iPart) + TurbPartState(1:3,iPart))

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/SQRT(SUM(udiff**2))
  ELSE
    urel = 0.
  END IF

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS)) THEN
    ! We ASSUME that these are the correct matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      RETURN
    ELSE
      ! parallel
      tauL(1) = C*ElemVolN(ElemID)/(betaSGS*SQRT(SUM(udiff**2)))
      ! perpendicular
      tauL(2) = C*ElemVolN(ElemID)/(betaSGS*2*SQRT(SUM(udiff**2)))

      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            G_SGS(i,j) = EXP(-RKdtFrac/tauL(2)) + (EXP(-RKdtFrac/tauL(1))&
              - EXP(-RKdtFrac/tauL(2)))*urel(i)*urel(j)
          ELSE
            G_SGS(i,j) = (EXP(-RKdtFrac/tauL(1)) - EXP(-RKdtFrac/tauL(2)))*urel(i)*urel(j)
          END IF
        END DO
      END DO

      B_SGS(:,:) = 0.

      ! Sum up turbulent contributions
      Pt(1:3) = 0.
      DO j = 1,3
        Pt(1:3) = Pt(1:3) + G_SGS(1:3,j)*TurbPartState(j,iPart)
      END DO
    END IF

  ! Valid SGS turbulent kinetic energy
  ELSE

    tauSGS   = C*ElemVolN(ElemID)/sigmaSGS

    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      ! parallel
      tauL(1) = tauSGS
      ! perpendicular
      tauL(2) = tauSGS

      ! Calculate drift and diffusion matrix
      G_SGS(:,:) = 0.
      DO i = 1,3
        G_SGS(i,i) = EXP(-RKdtFrac/tauL(2))
        B_SGS(i,i) = sigmaSGS*SQRT(1-EXP(-2*RKdtFrac/tauL(2)))
      END DO

    ELSE
      ! parallel
      tauL(1) = tauSGS/(SQRT(1+  betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))
      ! perpendicular
      tauL(2) = tauSGS/(SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))

      ! Calculate drift and diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            G_SGS(i,j) = EXP(-RKdtFrac/tauL(2)) + (EXP(-RKdtFrac/tauL(1)) - EXP(-RKdtFrac/tauL(2)))*urel(i)*urel(j)
            B_SGS(i,j) =  sigmaSGS*SQRT(1-EXP(-2*RKdtFrac/tauL(2)))                                                   &
                             + (sigmaSGS*SQRT(1-EXP(-2*RKdtFrac/tauL(1)))                                                   &
                             -  sigmaSGS*SQRT(1-EXP(-2*RKdtFrac/tauL(2))))                        *urel(i)*urel(j)
          ELSE
            G_SGS(i,j) = (EXP(-RKdtFrac/tauL(1)) - EXP(-RKdtFrac/tauL(2)))*urel(i)*urel(j)
            B_SGS(i,j) = (sigmaSGS*SQRT(1-EXP(-2*RKdtFrac/tauL(1)))                                                   &
                             -  sigmaSGS*SQRT(1-EXP(-2*RKdtFrac/tauL(2))))                        *urel(i)*urel(j)
          END IF
        END DO
      END DO
    END IF

    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3) + G_SGS(1:3,j)*TurbPartState(j,iPart) + B_SGS(1:3,j)*randomVar(j)
    END DO
  END IF

  TurbPartState(1:3,iPart) = Pt(1:3)
END DO

CASE DEFAULT
  CALL ABORT(__STAMP__, ' No particle SGS model given. This should not happen.')

END SELECT

END SUBROUTINE ParticleSGS


!===================================================================================================================================
!> Finalize the SGS model
!===================================================================================================================================
SUBROUTINE ParticleFinalizeSGS()
! MODULES
USE MOD_Particle_Vars,              ONLY: TurbPartState,TurbPt_temp
USE MOD_Particle_SGS_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(USGS)
SDEALLOCATE(USGSPart)
SDEALLOCATE(TurbPartState)
SDEALLOCATE(TurbPt_temp)
SDEALLOCATE(ElemVolN)
SDEALLOCATE(gradUxPart)
SDEALLOCATE(gradUyPart)
SDEALLOCATE(gradUzPart)
SDEALLOCATE(gradp)
SDEALLOCATE(gradpPart)

ParticleSGSInitIsDone=.FALSE.

END SUBROUTINE ParticleFinalizeSGS

END MODULE MOD_Particle_SGS
