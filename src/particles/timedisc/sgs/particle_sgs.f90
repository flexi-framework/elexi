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
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars          ,ONLY: ElemVol
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_SGS_Vars
USE MOD_Particle_Vars         ,ONLY: PDM,TurbPartState,nSpecies,Species!TurbPt_temp
USE MOD_Particle_Vars         ,ONLY: gradUx2,gradUy2,gradUz2
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
INTEGER               :: ALLOCSTAT,iSpec
!===================================================================================================================================

IF(ParticleSGSInitIsDone) RETURN

!LBWRITE(UNIT_stdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SGS MODEL...'

! SGS model
SGSModel = TRIM(GETSTR('Part-SGSModel','none'))

! SGS and RW are not compatible. Go back and think about what you were trying to simulate.
#if USE_RW
IF (SGSModel.NE.'none'.AND.RWModel.NE.'none') &
  CALL Abort(__STAMP__,'SGS and RW not compatible!')
#endif

SELECT CASE(SGSModel)
  CASE('Minier')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = MIN(PP_N-1,GETINT('Part-SGSNFilter'))
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:3       ,1:PDM%maxParticleNumber),      &
             gradp        (1,1:3     ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL Abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

    ! Overwrite RHS if necessary
    DO iSpec = 1, nSpecies
      IF (Species(iSpec)%RHSMethod .EQ. RHS_INERTIA) Species(iSpec)%RHSMethod = RHS_SGS1
    END DO

    ALLOCATE(divTau(1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)       &
#if !FAXEN_CORR
            ,gradUx2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  &
            ,gradUy2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  &
            ,gradUz2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  &
#endif
            )

  CASE('Sommerfeld')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = MIN(PP_N-1,GETINT('Part-SGSNFilter'))
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL Abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

CASE('Amiri')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = MIN(PP_N-1,GETINT('Part-SGSNFilter'))
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL Abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

    ! Overwrite RHS if necessary
    DO iSpec = 1, nSpecies
      Species(iSpec)%RHSMethod = RHS_SGS2
    END DO

CASE('Fukagata')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = MIN(PP_N-1,GETINT('Part-SGSNFilter'))
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL Abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

    ! Overwrite RHS if necessary
    DO iSpec = 1, nSpecies
      Species(iSpec)%RHSMethod = RHS_SGS2
    END DO

  CASE('Jin')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = MIN(PP_N-1,GETINT('Part-SGSNFilter'))
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL Abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

  CASE('Breuer')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = MIN(PP_N-1,GETINT('Part-SGSNFilter'))
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL Abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

  CASE('Breuer-Analytic')
    ! Set number of variables and SGS flag active
    nSGSVars = 3
    SGSinUse =  .TRUE.

    ! Init double-filtering for SGS turbulent kinetic energy
    !>> Modal Filter, default to cut-off at PP_N-2
    nSGSFilter = MIN(PP_N-1,GETINT('Part-SGSNFilter'))
    CALL InitSGSFilter()

    ! Allocate array to hold the SGS properties for every particle
    ALLOCATE(USGS         (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart     (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
!             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (           1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL Abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState = 0.
!    TurbPt_temp   = 0.
    ElemVolN      = ElemVol**(1./3.)/(PP_N+1)

  CASE('none')
    ! Do nothing
    nSGSVars = 0
    SGSinUse = .FALSE.

  CASE DEFAULT
    CALL Abort(__STAMP__, ' No valid particle subgrid scale (SGS) model given.')
END SELECT

ParticleSGSInitIsDone=.TRUE.

LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SGS MODEL DONE!'
LBWRITE(UNIT_stdOut,'(132("-"))')

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
LBWRITE(UNIT_stdOut,'(A)') ' Init SGS filter...'

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
USE MOD_Particle_Interpolation      ,ONLY: InterpolateFieldToParticle,InterpolateFieldToSingleParticle
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
REAL                          :: alpha_SGS
REAL                          :: rotu(3)            ! rotation of velovity
REAL                          :: mu                 ! viscosity
! BREUER SGS MODEL
REAL,PARAMETER                :: C=1.
REAL                          :: betaSGS
REAL                          :: tauSGS             ! SGS time scale
REAL                          :: tauL(2)            ! Parallel and perpendicular SGS time scale
REAL                          :: G_SGS(3,3)         ! Drift (analytic) / diffusion (num) matrix
REAL                          :: B_SGS(3,3)         ! Covariance (analytic) / exp. of the drift (num) matrix
REAL                          :: gradUxPart(1:3), gradUyPart(1:3),gradUzPart(1:3),gradpPart(1:3),divTauPart(1:3)
! SOMMERFELD SGS MODEL
REAL                          :: C_eps_Shotorban=19./12  ! model constant of Shotorban (according to HEINZ)
REAL                          :: C_T=0.24           ! model constant for lagrangian time scale
REAL                          :: C_L=3              ! model constant for eulerian time scale
! MINIER & PEIRANO SGS MODEL
REAL                          :: C_0=2.1            ! model constant for lagrangian time scale MINIER (according to HEINZ)
REAL                          :: trace_HR           ! trace matrix H*R
REAL                          :: trace_H            ! trace matrix H
REAL                          :: aux                ! auxiliary variable
REAL                          :: kSGS_tilde         ! new SGS kinetic energy
REAL                          :: ReynoldsTensor(3,3)! Reynolds Stress Tensor
REAL                          :: H_R_matrix(3,3)    ! H_ij * R_ij
REAL                          :: second_drift_term(3) ! term relative at the relative displacement times gradientU
REAL                          :: fluctuationState(3) ! fluctuation velocity: U - <U>
REAL                          :: D_matrix(3,3)      ! Symmetric matrix for Cholewski decomposition
REAL                          :: D_matrix_L(2)
! JIN MODEL
REAL                          :: deltaT_L           ! Lagrangian time scale
REAL                          :: deltaT_E           ! Eulerian integral timescale (Fede & Simonin,2006)
REAL                          :: St                 ! Stokes number
! FUKAGATA/AMIRI SGS MODEL
REAL                          :: U_M                ! representative SGS velocity
REAL,PARAMETER                :: C_c= 1.            ! factor of Stokes-Cunningham
REAL,PARAMETER                :: C_Fuka= 0.08       ! constant model of Fukagata
REAL,PARAMETER                :: C_eps=0.93         ! = (2/(3*Ko))**(3/2), Ko: Kolmogorov constant = 1.5
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
REAL                          :: R_E_diag(3)        ! Eulerian correlation (only main diagonal, vector form)
!===================================================================================================================================

SELECT CASE(SGSModel)

CASE('none')
!===================================================================================================================================
! DNS mode. Assume fully resolved turbulence
!===================================================================================================================================
! Do nothing


CASE('Minier')
!===================================================================================================================================
!> Minier, J.P., Peirano, E. et al.: "PDF model based on Langevin equation for polydispersed two-phase
!> flows applied to a bluff-body gas-solid flow". Physics of Fluids, 2004.
!
! FOR ANISOTROPIC TURBULENCE
!
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
USGS(2,:,:,:,:) = U(MOM1,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM1,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(3,:,:,:,:) = U(MOM2,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM2,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(4,:,:,:,:) = U(MOM3,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM3,:,:,:,:)/USGS(DENS,:,:,:,:)

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(3,USGS(VELV,:,:,:,:),3,USGSPart)

! Calculate pressure gradient
CALL Lifting_BR1_gen(1,1,UPrim(PRES:PRES,:,:,:,:),gradp(:,1,:,:,:,:),gradp(:,2,:,:,:,:),gradp(:,3,:,:,:,:))

CALL CalcDivTau()

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  kSGSPart = 0.5*SUM(USGSPart(1:3,iPart)**2.)

  ! Lagrangian time scale
  ElemID   = PEM%Element(iPart) - offsetElem

  ! ratio between lagrangian and eulerian time scale
  betaSGS  = 0.356 ! deltaT_L/deltaT_E

  ! Interpolate velocity and pressure gradients to the particle position
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUx(1:3,:,:,:,ElemID) ,3,gradUxPart)
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUy(1:3,:,:,:,ElemID) ,3,gradUyPart)
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUz(1:3,:,:,:,ElemID) ,3,gradUzPart)
  CALL InterpolateFieldToSingleParticle(iPart,3,gradp(1,1:3,:,:,:,ElemID),3,gradpPart )
  CALL InterpolateFieldToSingleParticle(iPart,3,divTau(1:3,:,:,:,ElemID) ,3,divTauPart)

  ! Relative velocity
  udiff(1:3) = PartState(PART_VELV,iPart) - FieldAtParticle(VELV,iPart)

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/VECNORM(udiff)
  ELSE
    urel = 0.
  END IF

  ! Calculate the dyn. viscosity
  mu=VISCOSITY_TEMPERATURE(FieldAtParticle(TEMP,iPart))

  second_drift_term(:) = - gradpPart(:)/FieldAtParticle(DENS,iPart)! + mu * divTauPart(:)
  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(kSGSPart)) THEN
    ! We ASSUME that these are the CORRECT matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      Pt = 0.
    ELSE
      ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j) - nabla p
      second_drift_term(:) = second_drift_term(:) + udiff(1) * gradUxPart(:) + udiff(2) * gradUyPart(:) + udiff(3) * gradUzPart(:)

      ! parallel
      tauL(1) = SQRT(  betaSGS**2*SUM(udiff**2)*3./2)/(ElemVolN(ElemID))*(0.5 + 0.75*C_0)
      ! perpendicular
      tauL(2) = SQRT(4*betaSGS**2*SUM(udiff**2)*3./2)/(ElemVolN(ElemID))*(0.5 + 0.75*C_0)

      ! Calculate drift and diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            G_SGS(i,j) = tauL(2)  + (tauL(1) - tauL(2)) *urel(i)*urel(j)
          ELSE
            G_SGS(i,j) =            (tauL(1) - tauL(2)) *urel(i)*urel(j)
          END IF
        END DO
      END DO
      ! as kSGS is ALMOSTZERO
      B_SGS(:,:) = 0.
    END IF

  ELSE  ! Valid SGS turbulent kinetic energy

    ! Calculate dissipation energy according to  Babak Shotorban & Farzad Mashayek (2006) A stochastic model for particle
    ! motion in large-eddy simulation, Journal of Turbulence, 7, N18, DOI: 10.1080/14685240600595685
    epsilonSGS = 1./ElemVolN(ElemID) * (kSGSPart**(3./2))
    tauSGS = (0.5 + 0.75*C_0)**(-1) * kSGSPart/epsilonSGS

    IF (ALMOSTZERO(MAXVAL(udiff))) THEN
      ! parallel and perpendicular
      tauL(:) = tauSGS

      ! Calculate auxiliary matrix H_ij and drift G_ij
      G_SGS    = 0.
      trace_HR = 0.; trace_H  = 0.
      DO i = 1,3
        DO j = 1,3
          ! Calculate the Reynolds stress tensor
          ReynoldsTensor(i,j) = USGSPart(i,iPart)*USGSPart(j,iPart)
          IF (i.EQ.j) THEN
            G_SGS(i,j) = tauL(2)
            trace_H    = trace_H  + G_SGS(i,j)
          END IF
        END DO
      END DO

      H_R_matrix = MATMUL(G_SGS,ReynoldsTensor)

      DO i = 1,3
        trace_HR = trace_HR + H_R_matrix(i,i)
      END DO

      kSGS_tilde = (3 * trace_HR) / (2 * kSGSPart * trace_H)

      ! parallel (1) and perpendicular (2)
      D_matrix_L(:) = tauL(:)*C_0*kSGS_tilde + 2./3 * (tauL(:)*kSGS_tilde-1)

      ! Calculate diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            B_SGS(i,j) = D_matrix_L(2)
          ELSE
            B_SGS(i,j) = 0.0
          END IF
        END DO
      END DO
      B_SGS(:,:) = epsilonSGS * B_SGS(:,:)

    ELSE ! with valid KSGS and valid UDIFF

      ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j)
      second_drift_term(:) = second_drift_term(:) + udiff(1) * gradUxPart(:) + udiff(2) * gradUyPart(:) + udiff(3) * gradUzPart(:)

      ! parallel
      tauL(1) = (SQRT(1+  betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))/tauSGS
      ! perpendicular
      tauL(2) = (SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart*3/2))/tauSGS

      ! Calculate auxiliary matrix H_ij and drift G_ij
      trace_HR = 0.; trace_H  = 0.
      DO i = 1,3
        DO j = 1,3
          ! Calculate the Reynolds stress tensor
          ReynoldsTensor(i,j) = USGSPart(i,iPart)*USGSPart(j,iPart)
          IF (i.EQ.j) THEN
            G_SGS(i,j) = tauL(2)  + (tauL(1) - tauL(2)) *urel(i)*urel(j)
            trace_H    = trace_H  + G_SGS(i,j)*tauSGS
          ELSE
            G_SGS(i,j) =            (tauL(1) - tauL(2)) *urel(i)*urel(j)
          END IF
        END DO
      END DO

      H_R_matrix = MATMUL(G_SGS*tauSGS,ReynoldsTensor)

      DO i = 1,3
        trace_HR = trace_HR + H_R_matrix(i,i)
      END DO

      kSGS_tilde = (3 * trace_HR) / (2 * kSGSPart * trace_H)

      ! parallel (1) and perpendicular (2)
      D_matrix_L(:) = tauSGS*tauL(:)*C_0*kSGS_tilde + 2./3 * (tauSGS*tauL(:)*kSGS_tilde-1)

      D_matrix = 0.
      ! Calculate symmetry matrix D_ij for Cholewski decomposition
      DO i = 1,3
        DO j = i,3
          IF (i.EQ.j) THEN
            D_matrix(i,j) = D_matrix_L(2) + (D_matrix_L(1) - D_matrix_L(2)) * urel(i)*urel(j)
          ELSE
            D_matrix(i,j) =                 (D_matrix_L(1) - D_matrix_L(2)) * urel(i)*urel(j)
            D_matrix(j,i) = D_matrix(i,j)
          END IF
        END DO
      END DO
      D_matrix(:,:) = epsilonSGS * D_matrix(:,:)

      ! Cholewski decomposition to get diffusion matrix B_SGS
      B_SGS(:,:) = 0.
      DO i = 1 , 3
        DO j = 1 , i
          aux = 0.
          DO k = 1 , j -1
            aux = aux + B_SGS(i,k)*B_SGS(j,k)
          END DO
          IF (i .EQ. j) THEN
            B_SGS(i,j) = SQRT(D_matrix(i,i) - aux) ! diagonal elements
          ELSE ! (i .GT. j)
            B_SGS(i,j) = (D_matrix(i,j) - aux) / B_SGS(j,j)    ! in our case only element B(3,2)
          END IF
        END DO
      END DO
    END IF
  END IF

  ! Calculate fluctuating velocity
  fluctuationState(1:3) = TurbPartState(1:3,iPart) - FieldAtParticle(VELV,iPart)

  Pt(1:3) = 0.
  DO j = 1,3
    Pt(1:3) = Pt(1:3)  - G_SGS(1:3,j)*fluctuationState(j)*dt + B_SGS(j,1:3)*SQRT(dt)*RandNormal()
  END DO

  ! FULL FLUID VELOCITY SEEN BY THE PARTICLES
  !  (U_p-U)*(\partial U_i / \partial x_j) * dt - ( \partial <P> / \partial x_i ) * dt
  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3) + second_drift_term(1:3)*dt

END DO


CASE('Sommerfeld')
!===================================================================================================================================
! Martin Sommerfeld 2001
!> Validation of a stochastic Lagrangian modelling approach for inter-particle collisions in homogeneous isotropic turbulence.
!> Institut fur Verfahrenstechnik, Fachbereich Ingenieurwissenschaften, Martin-Luther-Universitat at Halle-Wittenberg,
!> D-06099 Halle (Saale), Germany. Received 30 March 2000; received in revised form 3 May 2001
!
! FOR ISOTROPIC, HOMOGENEOUS TURBULENCE due to the coefficients f and g (apart from that it is for anisotropic turbulence)

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
USGS(2,:,:,:,:) = U(MOM1,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM1,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(3,:,:,:,:) = U(MOM2,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM2,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(4,:,:,:,:) = U(MOM3,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM3,:,:,:,:)/USGS(DENS,:,:,:,:)

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(3,USGS(MOMV,:,:,:,:),3,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart = 0.5*SUM(USGSPart(1:3,iPart)**2.)

  ! Standard deviation of fluctuation velocity
  sigmaSGS = SQRT(2./3.*kSGSPart)

  ! Standard deviation per direction, sigmaF = sigmaF**2
  sigmaF(1:3) = USGSPart(1:3,iPart)**2.

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Calculate dissipation energy according to  Babak Shotorban & Farzad Mashayek (2006) A stochastic model for particle
  ! motion in large-eddy simulation, Journal of Turbulence, 7, N18, DOI: 10.1080/14685240600595685
  epsilonSGS = C_eps_Shotorban/ElemVolN(ElemID) * (kSGSPart**(3./2))

  ! Relative velocity
  udiff(1:3) =  (FieldAtParticle(2:4,iPart) + TurbPartState(1:3,iPart)) - PartState(4:6,iPart)

  ! Relative displcement between particle and fluid element
  deltaR(1:3) = udiff(1:3)*dt       ! relative displacement in the three directions
  deltaR_mag = VECNORM(deltaR(1:3)) ! |U-U_p|*dt (scalar value)

  ! No SGS turbulent kinetic energy, avoid floating point exception
  IF (ALMOSTZERO(sigmaSGS)) THEN
    ! Lagrangian correlation R_L = 0
    Pt(1:3) = 0

  ! Valid SGS turbulent kinetic energy
  ELSE ! ALMOSTZERO(sigmaSGS)

    IF (ALMOSTZERO(deltaR_mag)) THEN

      ! Lagrangian time scale in the appropriate directions
      T_L(1:3) = C_T * sigmaF(1:3)/epsilonSGS

      ! Lagrangian correlation
      R_L(1:3) = EXP(-dt/T_L(1:3))

      R_E_diag(1:3) = 1.

      !Correlaction function: vector R_P
      R_P(1:3) = R_L(1:3) * R_E_diag(1:3)

      ! Calculate the fluctuation in the three directions
      DO i = 1,3
        Pt(i) = R_P(i)*TurbPartState(i,iPart) + SQRT(sigmaF(i))*SQRT(1-R_P(i)**2)* RandNormal()
      END DO

    ELSE ! ALMOSTZERO(deltaR_mag)

      ! Lagrangian time scale in the appropriate directions
      T_L(1:3) = C_T * sigmaF(1:3)/epsilonSGS

      ! Lagrangian correlation
      R_L(1:3) = EXP(-dt/T_L(1:3))

      ! Integral length scales L_E
      T_L_f = C_T *  (sigmaSGS**2) / epsilonSGS
      L_E = C_L * T_L_f * sigmaSGS

      ! Function g(dr) and f(dr)
      f_dr = EXP( -deltaR_mag / L_E)
      g_dr = (1 - deltaR_mag / (2*L_E) ) *  EXP(-deltaR_mag/L_E)

      R_E_diag(1:3) = g_dr + (f_dr - g_dr)*(deltaR(1:3)*deltaR(1:3))/(deltaR_mag**2)

      !Correlaction function: vector R_P
      R_P(1:3) = R_L(1:3) * R_E_diag(1:3)

      ! Calculate the fluctuation in the three directions
      DO i = 1,3
        Pt(i) = R_P(i)*TurbPartState(i,iPart) + SQRT(sigmaF(i))*SQRT(1.-R_P(i)**2)* RandNormal()
      END DO
    END IF
  END IF

  TurbPartState(1:3,iPart) =  Pt(1:3)
END DO


CASE('Amiri')
!===================================================================================================================================
! A. Elhami Amiri , S. Kazemzadeh Hannani & F. Mashayek (2006)
!>  Large-Eddy Simulation of Heavy-Particle Transport in Turbulent Channel Flow, Numerical Heat Transfer,
!> Part B: Fundamentals, 50:4, 285-313, DOI: 10.1080/10407790600859577
!
! FOR ANISOTROPIC TURBULENCE
!
!===================================================================================================================================

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat) ! DOUBLE FILTERING

! Obtain the high-pass filtered velocity field
USGS(2,:,:,:,:) = U(MOM1,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM1,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(3,:,:,:,:) = U(MOM2,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM2,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(4,:,:,:,:) = U(MOM3,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM3,:,:,:,:)/USGS(DENS,:,:,:,:)

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(3,USGS(MOMV,:,:,:,:),3,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Calculate the dyn. viscosity
  mu=VISCOSITY_TEMPERATURE(FieldAtParticle(TEMP,iPart))

  ElemID   = PEM%Element(iPart) - offsetElem

  ! Interpolate gradient to the particle position
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUx(1:3,:,:,:,ElemID),3,gradUxPart)
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUy(1:3,:,:,:,ElemID),3,gradUyPart)
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUz(1:3,:,:,:,ElemID),3,gradUzPart)

  kSGSPart = 0.5*SUM(USGSPart(1:3,iPart)**2.)

  ! Time scale of SGS scales = U_rms in case of isotropic turbolence
  sigmaSGS = SQRT(2./3*kSGSPart)

  ! No SGS turbulent kinetic energy, avoid floating point exception
  IF (ALMOSTZERO(kSGSPart)) THEN

    ! sigmaF --> 0
    Pt(1:3) = 0.

  ! Valid SGS turbulent kinetic energy
  ELSE

    ! Particle relaxation time
    taup = C_c * Species(PartSpecies(iPart))%DensityIC * PartState(PART_DIAM,iPart)**2 * 1./(18*mu)

    ! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
    rotu(1)=gradUyPart(3) - gradUzPart(2)    ! dw/dy-dv/dz
    rotu(2)=gradUzPart(1) - gradUxPart(3)    ! du/dz-dw/dx
    rotu(3)=gradUxPart(2) - gradUyPart(1)    ! dv/dx-du/dy

    ! Calculate the factor f1, f2, f3 for the anisotropy
    ! TODO: third direction?
    factor(1) = USGSPart(1,iPart) / sigmaSGS
    factor(2) = USGSPart(2,iPart) / sigmaSGS
    factor(3) = SQRT(3. - factor(1)**2 - factor(2)**2)

    alpha_SGS =  dt/taup
    betaSGS  = taup*VECNORM(rotu)
    beta1 = 1+betaSGS
    beta2 = 1-betaSGS

    ! U_square = (u'')^2
    sigmaF(1:3) = factor(1:3) * sigmaSGS * SQRT( 1./beta1 * (1-EXP(-alpha_SGS*beta1))  &
                                                -1./beta2 * EXP(-2*alpha_SGS) * (1-EXP(alpha_SGS*beta2)))

    DO i=1,3
      Pt(i) =  sigmaF(i)/dt * RandNormal()
    END DO

  END IF

  ! The SGS drag force is added to the particle RHS
  TurbPartState(1:3,iPart) =  Pt(1:3)

END DO


CASE('Fukagata')
!===================================================================================================================================
!> K. Fukagata (2004)
!> Dynamics of Brownian particles in a turbulent channel flow, K. Fukagata, S. Zahrai, F. H. Bark
!> Heat and Mass Transfer 40 (2004) 715–726. DOI 10.1007/s00231-003-0462-8
!
! FOR ISOTROPIC TURBULENCE
!
!===================================================================================================================================

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat) ! DOUBLE FILTERING

! Obtain the high-pass filtered velocity field
USGS(2,:,:,:,:) = U(MOM1,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM1,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(3,:,:,:,:) = U(MOM2,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM2,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(4,:,:,:,:) = U(MOM3,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM3,:,:,:,:)/USGS(DENS,:,:,:,:)

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(3,USGS(MOMV,:,:,:,:),3,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Calculate the dyn. viscosity
  mu = VISCOSITY_TEMPERATURE(FieldAtParticle(TEMP,iPart))

  ElemID = PEM%Element(iPart) - offsetElem

  ! Interpolate gradient to the particle position
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUx(1:3,:,:,:,ElemID),3,gradUxPart)
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUy(1:3,:,:,:,ElemID),3,gradUyPart)
  CALL InterpolateFieldToSingleParticle(iPart,3,gradUz(1:3,:,:,:,ElemID),3,gradUzPart)

  kSGSPart = 0.5*SUM(USGSPart(1:3,iPart)**2.)

  ! Particle relaxation time
  taup = C_c * Species(PartSpecies(iPart))%DensityIC * PartState(PART_DIAM,iPart)**2 * 1./(18*mu)

  ! Particle Schmidt number
  ! Sc = Species(MASS_SPHERE(Species(iSpec)%DensityIC, PartState(PART_DIAM,iPart)) * mu / FieldAtParticle(DENS,iPart) &
  !   / (taup * kSGSPart * FieldAtParticle(TEMP,iPart))

  ! Brownian force
  ! ni_B = SQRT(2*mu/FieldAtParticle(DENS,iPart)/(Sc * taup**2 * dt))

  ! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
  rotu(1)=gradUyPart(3) - gradUzPart(2)    ! dw/dy-dv/dz
  rotu(2)=gradUzPart(1) - gradUxPart(3)    ! du/dz-dw/dx
  rotu(3)=gradUxPart(2) - gradUyPart(1)    ! dv/dx-du/dy

  ! Calculation of deformation tensor d_ij
  d_ij_tensor(1,:)= 0.5*(/             2*gradUxPart(1), gradUyPart(1)+gradUxPart(2), gradUzPart(1)+gradUxPart(3)/)
  d_ij_tensor(2,:)= 0.5*(/ gradUyPart(1)+gradUxPart(2),             2*gradUyPart(2), gradUzPart(2)+gradUyPart(3)/)
  d_ij_tensor(3,:)= 0.5*(/ gradUzPart(1)+gradUxPart(3), gradUzPart(2)+gradUyPart(3), 2*gradUzPart(3)/)

  ! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity due to SGS velocity in time dt
  U_M = 2./SQRT(3.) * C_Fuka * C_eps * ElemVolN(ElemID) * VECNORM(d_ij_tensor(:,:))

  alpha_SGS =  dt/taup
  ! beta = taup  / T_f_p, T_f_p = 1./VECNORM(rotu)
  betaSGS = taup*VECNORM(rotu)
  beta1 = 1+betaSGS
  beta2 = 1-betaSGS

  sigmaSGS = U_M * SQRT( 1/beta1 * (1-EXP(-alpha_SGS*beta1)) - 1./beta2 * EXP(-2*alpha_SGS) * (1-EXP(alpha_SGS*beta2)))

  ! SGS fluid velocity by Fukagata, without Brownian force as the particle size is larger then the molecular length scale
  ! n_i^SGS = sigmaSGS/dt * G_i + n_i^B * G_i
  DO i=1,3
    Pt(i) =  sigmaSGS/dt * RandNormal() !+ ni_B * RandNormal()
  END DO

  ! The SGS drag force is added to the particle RHS
  TurbPartState(1:3,iPart) =  Pt(1:3)

END DO


CASE('Jin')
!===================================================================================================================================
! Jin (2010)
!> Guodong Jin a , Guo-Wei He a,*, Lian-Ping Wang b, Jian Zhang a Subgrid scale fluid velocity timescales seen by inertial
!  particles in large-eddy simulation of particle-laden turbulence. International Journal of Multiphase Flow 36 (2010) 432–437
!
! FOR ISOTROPIC TURBULENCE
!
!===================================================================================================================================

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
  IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat) ! DOUBLE FILTERING

! Obtain the high-pass filtered velocity field
USGS(2,:,:,:,:) = U(MOM1,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM1,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(3,:,:,:,:) = U(MOM2,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM2,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(4,:,:,:,:) = U(MOM3,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM3,:,:,:,:)/USGS(DENS,:,:,:,:)

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(3,USGS(MOMV,:,:,:,:),3,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Calculate the dyn. viscosity
  mu=VISCOSITY_TEMPERATURE(FieldAtParticle(TEMP,iPart))

  kSGSPart = 0.5*SUM(USGSPart(1:3,iPart)**2.)

  ! Time scale of SGS scales = U_rms in case of isotropic turbolence
  sigmaSGS = SQRT(2./3*kSGSPart)

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
    betaSGS = deltaT_L/deltaT_E

    ! Particle relaxation time
    taup = Species(PartSpecies(iPart))%DensityIC * PartState(PART_DIAM,iPart)**2 * 1./(18*mu)

    ! Stokes number
    St = taup/deltaT_E

    ! SGS autocorrelation time seen by a particle =  deltaT_L (from Berrouk et al. (2007))
    deltaT_L =  deltaT_E * (1 - (1-betaSGS) / (1+St)**(0.4*(1+0.01*St)))

    ! kSGSPart  = C_0 * kSGS, C_0=f(St), here C_0=1 is assumed
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3) - 1./deltaT_L*TurbPartState(j,iPart)*dt + SQRT(4./3.*kSGSPart*dt/deltaT_L) * RandNormal()
    END DO

  END IF

  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)

END DO


CASE('Breuer')
!===================================================================================================================================
! Breuer (2017) model, first option
!> Breuer, M. and Hoppe, F., "Influence of a cost–efficient Langevin subgrid-scale model on the dispersed phase of large–eddy
!simulations of turbulent bubble–laden and particle–laden flows." International Journal of Multiphase Flow, 89 (2017): 23-44.
!
! FOR ANISOTROPIC TURBULENCE
!
!===================================================================================================================================
betaSGS = 1.

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
  IF (iStage.NE.1) RETURN
END IF

! Filter the velocity field (low-pass)
USGS = U

! Filter overwrites the array in place. FilterMat already filled in InitSGSFilter
CALL Filter_Pointer(USGS,FilterMat)

! Obtain the high-pass filtered velocity field
USGS(2,:,:,:,:) = U(MOM1,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM1,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(3,:,:,:,:) = U(MOM2,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM2,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(4,:,:,:,:) = U(MOM3,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM3,:,:,:,:)/USGS(DENS,:,:,:,:)

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(3,USGS(MOMV,:,:,:,:),3,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  kSGSPart = 0.5*SUM(USGSPart(1:3,iPart)**2.)

  ! Time scale of SGS scales
  sigmaSGS = SQRT(2./3.*kSGSPart)

  ElemID   = PEM%Element(iPart) - offsetElem

  ! Relative velocity
  udiff(1:3) = PartState(PART_VELV,iPart) - (FieldAtParticle(VELV,iPart) + TurbPartState(1:3,iPart))

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/VECNORM(udiff)
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
      tauL(1) = C*ElemVolN(ElemID)/(betaSGS  *VECNORM(udiff))
      ! perpendicular
      tauL(2) = C*ElemVolN(ElemID)/(betaSGS*2*VECNORM(udiff))

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
!    IPWRITE(*,*) 'Obtained SGS velocity of',VECNORM(TurbPartState(1:3,iPart)), &
!                 'while resolved velocity is',VECNORM(PartState(4:6,iPart)),   &
!                 'for',iPart
!    TurbPartState(:,iPart) = 0.
!  END IF

  END DO  ! For Breuer

CASE('Breuer-Analytic')
!===================================================================================================================================
! Breuer (2017) model, second option
!> Breuer, M. and Hoppe, F., "Influence of a cost–efficient Langevin subgrid-scale model on the dispersed phase of large–eddy
!simulations of turbulent bubble–laden and particle–laden flows." International Journal of Multiphase Flow, 89 (2017): 23-44.
!
! FOR ANISOTROPIC TURBULENCE
!
!===================================================================================================================================
betaSGS = 1.

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
USGS(2,:,:,:,:) = U(MOM1,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM1,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(3,:,:,:,:) = U(MOM2,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM2,:,:,:,:)/USGS(DENS,:,:,:,:)
USGS(4,:,:,:,:) = U(MOM3,:,:,:,:)/U(DENS,:,:,:,:) - USGS(MOM3,:,:,:,:)/USGS(DENS,:,:,:,:)

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(3,USGS(MOMV,:,:,:,:),3,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  kSGSPart = 0.5*SUM(USGSPart(1:3,iPart)**2.)

  ! Time scale of SGS scales
  sigmaSGS = SQRT(2./3*kSGSPart)

  ! Relative velocity
  udiff(1:3) = PartState(4:6,iPart) - (FieldAtParticle(2:4,iPart) + TurbPartState(1:3,iPart))

  ElemID   = PEM%Element(iPart) - offsetElem

  IF (ALMOSTZERO(MAXVAL(udiff))) THEN
    urel = 0.
  ELSE
    urel = udiff/VECNORM(udiff)
  END IF

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS)) THEN
    ! We ASSUME that these are the correct matrix indices
    IF (ALL(urel .EQ. 0)) RETURN

    ! parallel
    tauL(1) = C*ElemVolN(ElemID)/(betaSGS  *VECNORM(udiff))
    ! perpendicular
    tauL(2) = C*ElemVolN(ElemID)/(betaSGS*2*VECNORM(udiff))

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

  ! Valid SGS turbulent kinetic energy
  ELSE

    tauSGS   = C*ElemVolN(ElemID)/sigmaSGS

    IF (ALL(urel .EQ. 0)) THEN
      ! parallel
      tauL(1) = tauSGS
      ! perpendicular
      tauL(2) = tauSGS

      ! Calculate drift and diffusion matrix
      G_SGS(:,:) = 0.; B_SGS(:,:) = 0.
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
  CALL Abort(__STAMP__, ' No particle SGS model given. This should not happen.')

END SELECT

END SUBROUTINE ParticleSGS


SUBROUTINE CalcDivTau()
!===================================================================================================================================
! Compute tau
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,      ONLY: s23
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Lifting_BR1_gen,    ONLY: Lifting_BR1_gen
#if !FAXEN_CORR
USE MOD_Lifting_Vars,       ONLY: gradUx,gradUy,gradUz
USE MOD_Lifting_Vars,       ONLY: gradUx_master,gradUx_slave
USE MOD_Lifting_Vars,       ONLY: gradUy_master,gradUy_slave
USE MOD_Lifting_Vars,       ONLY: gradUz_master,gradUz_slave
#endif
USE MOD_Particle_Vars,      ONLY: gradUx2,gradUy2,gradUz2
USE MOD_Particle_SGS_Vars,  ONLY: divTau
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if !FAXEN_CORR
! Calculate the second gradient of the velocity and output \nabla \cdot tau
CALL Lifting_BR1_gen(3,3,gradUx(LIFT_VELV,:,:,:,:),gradUx_master(LIFT_VELV,:,:,:),gradUx_slave(LIFT_VELV,:,:,:),&
                    gradUx2(1,:,:,:,:,:),gradUx2(2,:,:,:,:,:),gradUx2(3,:,:,:,:,:))
CALL Lifting_BR1_gen(3,3,gradUy(LIFT_VELV,:,:,:,:),gradUy_master(LIFT_VELV,:,:,:),gradUy_slave( LIFT_VELV,:,:,:),&
                    gradUy2(1,:,:,:,:,:),gradUy2(2,:,:,:,:,:),gradUy2(3,:,:,:,:,:))
CALL Lifting_BR1_gen(3,3,gradUz(LIFT_VELV,:,:,:,:),gradUz_master(LIFT_VELV,:,:,:),gradUz_slave( LIFT_VELV,:,:,:),&
                    gradUz2(1,:,:,:,:,:),gradUz2(2,:,:,:,:,:),gradUz2(3,:,:,:,:,:))
#endif

! Time derivative of VELV + GradientVELV or RHS of NS equations (same result) for the calculation of the substantial derivative
!divtau(1,:,:,:,:) = gradUx2(1,LIFT_VEL1,:,:,:,:) + gradUy2(2,LIFT_VEL1,:,:,:,:) + gradUz2(3,LIFT_VEL1,:,:,:,:) + &
!                    s13 * (gradUx2(1,LIFT_VEL1,:,:,:,:) + gradUy2(1,LIFT_VEL2,:,:,:,:) + gradUz2(1,LIFT_VEL3,:,:,:,:))
divTau(1,:,:,:,:) = 2*gradUx2(1,1,:,:,:,:) + gradUy2(2,1,:,:,:,:) + gradUz2(3,1,:,:,:,:) + gradUy2(1,2,:,:,:,:) + &
                    + gradUz2(1,3,:,:,:,:) - s23 * (gradUx2(1,1,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(3,3,:,:,:,:))
!divtau(2,:,:,:,:) = gradUx2(1,LIFT_VEL2,:,:,:,:) + gradUy2(2,LIFT_VEL2,:,:,:,:) + gradUz2(3,LIFT_VEL2,:,:,:,:) + &
!                    s13 * (gradUx2(2,LIFT_VEL1,:,:,:,:) + gradUy2(2,LIFT_VEL2,:,:,:,:) + gradUz2(2,LIFT_VEL3,:,:,:,:))
divTau(2,:,:,:,:) = gradUx2(2,1,:,:,:,:) + 2*gradUy2(2,2,:,:,:,:) + gradUz2(2,3,:,:,:,:) + gradUx2(1,2,:,:,:,:) + &
                    + gradUz2(3,2,:,:,:,:) - s23 * (gradUx2(1,1,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(3,3,:,:,:,:))
!divtau(3,:,:,:,:) = gradUx2(1,LIFT_VEL3,:,:,:,:) + gradUy2(2,LIFT_VEL3,:,:,:,:) + gradUz2(3,LIFT_VEL3,:,:,:,:) + &
!                    s13 * (gradUx2(3,LIFT_VEL1,:,:,:,:) + gradUy2(3,LIFT_VEL2,:,:,:,:) + gradUz2(3,LIFT_VEL3,:,:,:,:))
divTau(2,:,:,:,:) = gradUx2(3,1,:,:,:,:) + gradUy2(3,2,:,:,:,:) + 2*gradUz2(3,3,:,:,:,:) + gradUx2(1,3,:,:,:,:) + &
                    + gradUz2(2,3,:,:,:,:) - s23 * (gradUx2(1,1,:,:,:,:) + gradUy2(2,2,:,:,:,:) + gradUz2(3,3,:,:,:,:))

END SUBROUTINE CalcDivTau


!===================================================================================================================================
!> Finalize the SGS model
!===================================================================================================================================
SUBROUTINE ParticleFinalizeSGS()
! MODULES
USE MOD_Particle_Vars,              ONLY: TurbPartState,TurbPt_temp
#if !FAXEN_CORR
USE MOD_Particle_Vars,              ONLY: gradUx2,gradUy2,gradUz2
#endif
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
SDEALLOCATE(gradp)
SDEALLOCATE(divTau)
#if !FAXEN_CORR
SDEALLOCATE(gradUx2)
SDEALLOCATE(gradUy2)
SDEALLOCATE(gradUz2)
#endif

ParticleSGSInitIsDone=.FALSE.

END SUBROUTINE ParticleFinalizeSGS

END MODULE MOD_Particle_SGS
