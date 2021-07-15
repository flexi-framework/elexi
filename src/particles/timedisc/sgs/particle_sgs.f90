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
USE MOD_Globals               ,ONLY: ABORT,Unit_STDOUT
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

!SWRITE(UNIT_StdOut,'(132("-"))')
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
    ALLOCATE(USGS             (1:PP_nVar ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             USGSPart         (1:4       ,1:PDM%maxParticleNumber),      &
!             kSGS             (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart 	      (           1:PDM%maxParticleNumber),      &
	     gradp            (1:3       ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             gradpPart        (1:3       ,1:PDM%maxParticleNumber),      &
             gradUxPart       (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart       (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart       (1:3       ,1:PDM%maxParticleNumber),      &
             epsilonSGS       (           1:PDM%maxParticleNumber),      &
	     pressur_gradient (1:3       ,1:PDM%maxParticleNumber),      &
	     second_drift_term(1:3       ,1:PDM%maxParticleNumber),      &
	     fluctuationState (1:3       ,1:PDM%maxParticleNumber),      &
             b_L              (1:2       ,1:PDM%maxParticleNumber),      &
             G_matrix         (1:3,1:3   ,1:PDM%maxParticleNumber),      &
	     D_matrix         (1:3,1:3   ,1:PDM%maxParticleNumber),      &
	     B_matrix         (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             TurbPartState    (1:nSGSVars,1:PDM%maxParticleNumber),      &
     !        TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
 !   TurbPt_temp     = 0.
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             sigmaF       (           1:PDM%maxParticleNumber),      &
             epsilonSGS   (           1:PDM%maxParticleNumber),      &
             deltaR       (1:3       ,1:PDM%maxParticleNumber),      &
             DELTA_R      (           1:PDM%maxParticleNumber),      &
      sigmaF_vector_square(1:3       ,1:PDM%maxParticleNumber),      &
             T_L          (1:3       ,1:PDM%maxParticleNumber),      &             
 !            T_L          (           1:PDM%maxParticleNumber),      &
             T_L_f        (           1:PDM%maxParticleNumber),      &
             R_L          (1:3       ,1:PDM%maxParticleNumber),      &             
 !            R_L          (           1:PDM%maxParticleNumber),      &
             L_E          (           1:PDM%maxParticleNumber),      &             
             F_dr         (           1:PDM%maxParticleNumber),      &
             G_dr         (           1:PDM%maxParticleNumber),      &                    
  !           R_E          (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             diag_R_E     (1:3       ,1:PDM%maxParticleNumber),      &
  !           R_P          (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             R_P          (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
!    TurbPt_temp     = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)


CASE('Amiri-Analytic')
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             taup_Amiri   (           1:PDM%maxParticleNumber),      &
             gradUxPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart   (1:3       ,1:PDM%maxParticleNumber),      &
             T_f_p        (           1:PDM%maxParticleNumber),      &
             U_square     (1:3       ,1:PDM%maxParticleNumber),      &
             alpha_Amiri  (           1:PDM%maxParticleNumber),      &
             beta_Amiri   (           1:PDM%maxParticleNumber),      &
             sigmaS_Amiri (1:3       ,1:PDM%maxParticleNumber),      &
             ni_SGS_Amiri (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             taup_Amiri   (           1:PDM%maxParticleNumber),      &
             gradUxPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart   (1:3       ,1:PDM%maxParticleNumber),      &
             T_f_p        (           1:PDM%maxParticleNumber),      &
             U_square     (1:3       ,1:PDM%maxParticleNumber),      &
             alpha_Amiri  (           1:PDM%maxParticleNumber),      &
             beta_Amiri   (           1:PDM%maxParticleNumber),      &
             sigmaS_Amiri (1:3       ,1:PDM%maxParticleNumber),      &
             ni_SGS_Amiri (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

CASE('Fukagata-Analytic')
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             sigmaSGS     (           1:PDM%maxParticleNumber),      &
             taup_Fuka    (           1:PDM%maxParticleNumber),      &
             Sc_Fuka      (           1:PDM%maxParticleNumber),      &
             gradUxPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart   (1:3       ,1:PDM%maxParticleNumber),      &
             T_f_p        (           1:PDM%maxParticleNumber),      &
             U_M          (           1:PDM%maxParticleNumber),      & 
             alpha_Fuka   (           1:PDM%maxParticleNumber),      &
             beta_Fuka    (           1:PDM%maxParticleNumber),      &
             sigmaS_Fuka  (	      1:PDM%maxParticleNumber),      &
             ni_SGS       (1:3       ,1:PDM%maxParticleNumber),      &
             ni_Fuka      (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             sigmaSGS     (           1:PDM%maxParticleNumber),      &
             taup_Fuka    (           1:PDM%maxParticleNumber),      &
             Sc_Fuka      (           1:PDM%maxParticleNumber),      &
  !           ni_B         (           1:PDM%maxParticleNumber),      &
             gradUxPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUyPart   (1:3       ,1:PDM%maxParticleNumber),      &
             gradUzPart   (1:3       ,1:PDM%maxParticleNumber),      &
             T_f_p        (           1:PDM%maxParticleNumber),      &
             U_M          (           1:PDM%maxParticleNumber),      & 
             alpha_Fuka   (           1:PDM%maxParticleNumber),      &
             beta_Fuka    (           1:PDM%maxParticleNumber),      &
             sigmaS_Fuka  (	      1:PDM%maxParticleNumber),      &
             ni_SGS       (1:3       ,1:PDM%maxParticleNumber),      &
             ni_Fuka      (1:3       ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')

    TurbPartState   = 0.
    ElemVolN        = ElemVol**(1./3.)/(PP_N+1)

  CASE('Jin-Analytic')
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             sigmaSGS     (           1:PDM%maxParticleNumber),      &
             deltaT_L     (           1:PDM%maxParticleNumber),      &
             k_cf         (           1:PDM%maxParticleNumber),      &
             deltaT_E     (           1:PDM%maxParticleNumber),      &
             beta_Jin     (           1:PDM%maxParticleNumber),      &
             taup_Jin     (           1:PDM%maxParticleNumber),      &
             St_Jin       (           1:PDM%maxParticleNumber),      &
             deltaT_L_p   (           1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             sigmaSGS     (           1:PDM%maxParticleNumber),      &
             deltaT_L     (           1:PDM%maxParticleNumber),      &
             k_cf         (           1:PDM%maxParticleNumber),      &
             deltaT_E     (           1:PDM%maxParticleNumber),      &
             beta_Jin     (           1:PDM%maxParticleNumber),      &
             taup_Jin     (           1:PDM%maxParticleNumber),      &
             St_Jin       (           1:PDM%maxParticleNumber),      &
             deltaT_L_p   (           1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             sigmaSGS     (           1:PDM%maxParticleNumber),      &
             tauSGS       (           1:PDM%maxParticleNumber),      &
             tauL         (1:2       ,1:PDM%maxParticleNumber),      &
             G_SGS        (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             B_SGS        (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
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
!             kSGS         (           0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (           1:PDM%maxParticleNumber),      &
             sigmaSGS     (           1:PDM%maxParticleNumber),      &
             tauSGS       (           1:PDM%maxParticleNumber),      &
             tauL         (1:2       ,1:PDM%maxParticleNumber),      &
             E_SGS        (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             W_SGS        (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),      &
             ElemVolN     (         1:nElems),STAT=ALLOCSTAT)
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
SWRITE(UNIT_StdOut,'(132("-"))')

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
USE MOD_DG_Vars                     ,ONLY: U
USE MOD_Filter                      ,ONLY: Filter_Pointer
USE MOD_Filter_Vars                 ,ONLY: FilterMat
USE MOD_Mesh_Vars                   ,ONLY: offsetElem, nElems
USE MOD_Particle_Globals
USE MOD_Particle_SGS_Vars
USE MOD_Particle_Interpolation      ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars ,ONLY: FieldAtParticle
USE MOD_Particle_Vars               ,ONLY: TurbPartState,PDM,PEM,PartState
!!! other SGS MODELS !!!
USE MOD_EOS_Vars,          ONLY : mu0
USE MOD_Particle_Vars,     ONLY : Species, PartSpecies
USE MOD_Lifting_Vars, 	   ONLY : gradUx, gradUy, gradUz
! M&P SGS MODEL - Modules used in order to calculate pressur gradient
USE MOD_Lifting_BR1_gen,    ONLY: Lifting_BR1_gen
USE MOD_EoS,                ONLY: ConsToPrim
!USE MOD_Particle_Vars               ,ONLY: TurbPt_Temp
USE MOD_TimeDisc_Vars               ,ONLY: nRKStages, RKC
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
REAL,PARAMETER                :: C=1.
REAL,PARAMETER                :: betaSGS=1.
REAL                          :: udiff(3)
REAL                          :: urel(3)    ! Normalized relative velocity vector
INTEGER                       :: ElemID, i, j, k, iPart, iElem
REAL                          :: Pt(1:3)
REAL                          :: dxminus(1:3)
REAL                          :: dxplus(1:3)
REAL                          :: RKdtFrac
REAL 			      :: filterwidth ! IN ORDER TO VALIDATE FLEXI IMPLEMENTATION
! LOCAL VARIABLES - FUKAGATA/AMIRI SGS MODEL
REAL                          :: nu  		    ! kynematic viscosity
REAL,PARAMETER                :: C_c= 1. 	    ! factor of Stokes-Cunningham
REAL,PARAMETER                :: C_Fuka= 0.93 	    ! constant model of Fukagata
REAL,PARAMETER                :: C_eps=0.08 	    ! = (2/(3*Ko))**(3/2), Ko: Kolmogorov constant = 1.5
REAL                          :: beta1
REAL                          :: beta2
REAL    		      :: rot_U(3) 	    ! rotor of velovity
REAL                          :: d_ij_tensor(3,3)   ! rate of deformation tensor
REAL                          :: d_ij  		    ! used for the calculation of representative u_sgs
REAL    		      :: factor(3)          ! factor to take into account anisotropy (AMIRI MODEL)
! LOCAL VARIABLES - SOMMERFELD SGS MODEL
REAL                          :: C_eps_Shotorban=19./12  ! model constant of Shotorban (according to HEINZ)
REAL                          :: C_T=0.24           ! model constant for lagrangian time scale
REAL                          :: C_L=3              ! model constant for eulerian time scale
!REAL                          :: M_ones(3,3)        ! identity matrix
! LOCAL VARIABLES - MINIER & PEIRANO SGS MODEL		
REAL                          :: C_0=19./12          ! model constant for lagrangian time scale MINIER (according to HEINZ)
REAL			      :: Lambda_factor      ! factor in order to calculate diffusion matrix
REAL                          :: trace_HR           ! trace matrix H*R
REAL                          :: trace_H            ! trace matrix H
REAL                          :: aux                ! auxiliary variable 
REAL                          :: gradient_U(3,3)    ! Velocity gradient
REAL                          :: H_matrix(3,3)      ! 
REAL                          :: ReynoldsTensor(3,3)! Reynolds Stress Tensor
REAL                          :: H_R_matrix(3,3)    ! H_ij * R_ij
REAL    		      :: pg(3)		! Pressure gradient times dt
REAL    		      :: sdt(3)		! Second drift term times dt
!!! only for validation MINIER
!REAL    		      :: ddrift(3)		! 
!REAL    		      :: diffusion(3)		! 
! LOCAL ALLOCATABLE VARIABLES FOR M&P SGS model
REAL,ALLOCATABLE              :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE              :: gradp_local(:,:,:,:,:,:)
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
ALLOCATE(U_local(PRIM,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
ALLOCATE(gradp_local(1,3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
DO iElem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  CALL ConsToPrim(U_local(:,i,j,k,iElem),U(:,i,j,k,iElem))
END DO; END DO; END DO; END DO
CALL Lifting_BR1_gen(1,1,U_local(PRES:PRES,:,:,:,:),gradp_local(:,1,:,:,:,:),gradp_local(:,2,:,:,:,:),gradp_local(:,3,:,:,:,:))
gradp = gradp_local(1,:,:,:,:,:)
DEALLOCATE(U_local,gradp_local)

! Interpolate pressur gradient to particel position
CALL InterpolateFieldToParticle(3,gradp(1:3,:,:,:,:),3,gradpPart) 

DO iPart = 1,PDM%ParticleVecLength

  !filterwidth = 2*PI/32    	! to validate the implementation
  !USGSPart(2:4,iPart) = 0.1	! to validate the implementation
!PRINT *, 'filterwidth:  ', filterwidth
!PRINT *, 'USGSPart:  ', USGSPart(2:4,iPart)

  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.) 
!PRINT *, 'kSGSPart:  ', kSGSPart(iPart)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Calculate dissipation energy according to  Babak Shotorban & Farzad Mashayek (2006) A stochastic model for particle
  ! motion in large-eddy simulation, Journal of Turbulence, 7, N18, DOI: 10.1080/14685240600595685
  epsilonSGS(iPart) = C_eps_Shotorban/ElemVolN(ElemID) * (kSGSPart(iPart)**(3./2)) 
  !epsilonSGS(iPart) = C_eps_Shotorban/filterwidth * (kSGSPart(iPart)**(3./2))  ! to validate the implementation
!PRINT *, 'epsilonSGS:  ', epsilonSGS(iPart)
  
  ! Relative velocity
  udiff(1:3) = PartState(4:6,iPart) - FieldAtParticle(2:4,iPart) !+ TurbPartState(1:3,iPart) ! only mean velocity
   !udiff(1:3) = 1. ! to validate the implementation

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/SQRT(SUM(udiff**2))
  ELSE
    urel = 0.
  END IF

  ! Calculate fluctuation velocity
  !FieldAtParticle(2:4,iPart) = 0.5 ! to validate the implementation
  fluctuationState(1:3,iPart) = TurbPartState(1:3,iPart) - FieldAtParticle(2:4,iPart)
!PRINT *, 'fluctuationState:  ', fluctuationState(1:3,iPart)

  ! Calculate gradient_U for the drift vector: \partial U_i / \partial x_j
  gradient_U(1,:) = (/	gradUxPart(1,iPart), 	gradUyPart(1,iPart),	 gradUzPart(1,iPart)/) 	
  gradient_U(2,:) = (/	gradUxPart(2,iPart), 	gradUyPart(2,iPart),	 gradUzPart(2,iPart)/) 	
  gradient_U(3,:) = (/	gradUxPart(3,iPart), 	gradUyPart(3,iPart),	 gradUzPart(3,iPart)/) 	       
  !gradient_U(1:3,1:3) = 1. ! to validate the implementation

  !  Pressur gradient: 1/\rho_f * \partial <P> / \partial x_i
  pressur_gradient(1:3,iPart) = gradpPart(1:3,iPart)/FieldAtParticle(1,iPart)
  !pressur_gradient(1:3,iPart) = 15  ! to validate the implementation
  !pressur_gradient(:,iPart) = (/1, 1, 15/)  ! to validate the implementation
!PRINT *, 'pressur_gradient:  ', pressur_gradient(1:3,iPart)

  ! No SGS turbulent kinetic energy, avoid float error

  IF (ALMOSTZERO(kSGSPart(iPart))) THEN ! 1st IF
    ! We ASSUME that these are the CORRECT matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN   ! 2ND IF
      Pt(1:3) = 0.
      second_drift_term(1:3,iPart) = 0.
    ELSE     ! 2ND ELSE

      ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j)
      second_drift_term(1:3,iPart) = MATMUL(udiff,gradient_U)
	  
      ! parallel Csanady factor
      b_L(1,iPart) = betaSGS*SQRT(SUM(udiff**2))
      ! perpendicular Csanady factor
      b_L(2,iPart) = betaSGS*2*SQRT(SUM(udiff**2))

      ! Calculate auxliary matrix H_ij, drift matrix G_ij  
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            H_matrix(i,j) = b_L(2,iPart)  + (b_L(1,iPart)  - b_L(2,iPart)) *urel(i)*urel(j)
            ELSE
            H_matrix(i,j) =                 (b_L(1,iPart)  - b_L(2,iPart)) *urel(i)*urel(j)
          END IF
	G_matrix(i,j,iPart) = -(0.5 + 3*C_0/4) * (C_eps_Shotorban/(SQRT(2./3))*1./ElemVolN(ElemID)) * H_matrix(i,j)
	END DO
      END DO
    
	! Calculate the Reynolds stress tensor
      DO i = 1,3 
        DO j = 1,3
          ReynoldsTensor(i,j) = (USGSPart(i+1,iPart)/FieldAtParticle(1,iPart))*(USGSPart(j+1,iPart)/FieldAtParticle(1,iPart))
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
            D_matrix(i,j,iPart) = ( C_0/SQRT(8./27) + SQRT(3./2) ) * Lambda_factor * H_matrix(i,j) 
          END DO
        END DO
	D_matrix(:,:,iPart) = C_eps_Shotorban/ElemVolN(ElemID) * D_matrix(:,:,iPart)

	! Cholewski decomposition to get diffusion matrix B_matrix
	B_matrix(:,:,iPart) = 0.	
	
	B_matrix(1,1,iPart) = sqrt( D_matrix(1,1,iPart) ) 						! element (1,1)
	DO i = 2 , 3									
	  B_matrix(i,1,iPart) = D_matrix(1,i,iPart) / B_matrix(1,1,iPart)		! first column
	END DO

	DO j = 2 , 3
	  DO i = 1 , 3
	    IF( i .LT. j ) THEN
	      B_matrix(i,j,iPart) = 0.0   								! upper triangolar elements
  	    ELSE
	      IF( i .EQ. j) THEN
		aux = 0.0
		DO k = 1 , i -1
		  aux = aux + B_matrix(i,k,iPart)**2
		END DO
		B_matrix(i,j,iPart) = sqrt( D_matrix(i,i,iPart) - aux ) ! diagonal elements
	      ELSE ! (i .GT. j)
		aux = 0.0
		DO k = 1 , i - 1
		  aux = aux + B_matrix(j,k,iPart) * B_matrix(i,k,iPart)
		END DO
	        B_matrix(i,j,iPart) = ( D_matrix(j,i,iPart) - aux ) / B_matrix(j,j,iPart)    ! in our case only element B(3,2)
	      END IF
            END IF
	  END DO
	END DO

      ! EULER
      ! Sum up turbulent contributions
      Pt(1:3) = 0.
      DO  j=1,3					
      Pt(1:3) = Pt(1:3) + G_matrix(1:3,j,iPart)*fluctuationState(j,iPart)*dt + B_matrix(1:3,j,iPart)*SQRT(dt)*RandNormal()
      END DO
    END IF !2ND END IF 

  ! Valid SGS turbulent kinetic energy
  ELSE  ! 1ST ELSE

    IF(ALMOSTZERO(MAXVAL(udiff)))THEN  ! 3RD IF
   
     ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j)
     second_drift_term(1:3,iPart) = 0
	  
     ! parallel Csanady factor
     b_L(1,iPart) = 1
     ! perpendicular Csanady factor
     b_L(2,iPart) = 1

      ! Calculate auxliary matrix H_ij, drift matrix G_ij  
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            H_matrix(i,j) = b_L(2,iPart) 
            ELSE
            H_matrix(i,j) = 0
          END IF
	G_matrix(i,j,iPart) = -(0.5 + 3*C_0/4) * epsilonSGS(iPart)/kSGSPart(iPart) * H_matrix(i,j)
	END DO
      END DO
    
	! Calculate the Reynolds stress tensor
      DO i = 1,3 
        DO j = 1,3
          ReynoldsTensor(i,j) = (USGSPart(i+1,iPart)/FieldAtParticle(1,iPart))*(USGSPart(j+1,iPart)/FieldAtParticle(1,iPart))
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
	
	Lambda_factor = (3 * trace_HR) / (2 * kSGSPart(iPart) * trace_H)

	! Calculate simmetrix matrix D_ij for Cholewski decomposition
	DO i = 1,3
          DO j = 1,3
            IF (i.EQ.j) THEN
            D_matrix(i,j,iPart) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j) -1)/3
            ELSE
            D_matrix(i,j,iPart) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j)   )/3
            END IF
          END DO
        END DO
	D_matrix(:,:,iPart) = epsilonSGS(iPart) * D_matrix(:,:,iPart)

	! Cholewski decomposition to get diffusion matrix B_matrix
	B_matrix(:,:,iPart) = 0.	
	
	B_matrix(1,1,iPart) = sqrt( D_matrix(1,1,iPart) ) 			! element (1,1)

	DO i = 2 , 3									
	  B_matrix(i,1,iPart) = D_matrix(1,i,iPart) / B_matrix(1,1,iPart)	! first column
	END DO

	
	DO j = 2 , 3
	  DO i = 1 , 3
	    IF( i .LT. j ) THEN
	      B_matrix(i,j,iPart) = 0.0   	! upper triangolar elements
  	    ELSE
	      IF( i .EQ. j) THEN
		aux = 0.0
		DO k = 1 , i -1
		  aux = aux + B_matrix(i,k,iPart)**2
		END DO
		B_matrix(i,j,iPart) = sqrt( D_matrix(i,i,iPart) - aux ) ! diagonal elements
	      ELSE ! (i .GT. j)
		aux = 0.0
		DO k = 1 , i - 1
		  aux = aux + B_matrix(j,k,iPart) * B_matrix(i,k,iPart)
		END DO
	        B_matrix(i,j,iPart) = ( D_matrix(j,i,iPart) - aux ) / B_matrix(j,j,iPart)    ! in our case only element B(3,2)
	      END IF
            END IF
	  END DO
	END DO

    ELSE ! with valid KSGS and valid UDIFF      ! 3RD ELSE
	
	  ! Calculate second term of drift vector: (U_p-U)*(\partial U_i / \partial x_j)
	  second_drift_term(1:3,iPart) = MATMUL(udiff,gradient_U)
!PRINT *, 'second_drift_term:  ', second_drift_term(1:3,iPart)
	  
      ! parallel Csanady factor
      b_L(1,iPart) = (SQRT(1+  betaSGS**2*SUM(udiff**2)/kSGSPart(iPart)*3/2))
      ! perpendicular Csanady factor
      b_L(2,iPart) = (SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart(iPart)*3/2))
!PRINT *, 'b_L:  ', b_L(1:2,iPart)
	
      ! Calculate auxliary matrix H_ij, drift matrix G_ij  
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            H_matrix(i,j) = b_L(2,iPart)  + (b_L(1,iPart)  - b_L(2,iPart)) *urel(i)*urel(j)
            ELSE
            H_matrix(i,j) =                 (b_L(1,iPart)  - b_L(2,iPart)) *urel(i)*urel(j)
          END IF
	G_matrix(i,j,iPart) = -(0.5 + 3*C_0/4) * epsilonSGS(iPart)/kSGSPart(iPart) * H_matrix(i,j)
	END DO
      END DO
!PRINT *, 'H_matrix:  ', H_matrix(1:3,1:3)	
!PRINT *, 'G_matrix:  ', G_matrix(1:3,1:3,iPart)	
    
	! Calculate the Reynolds stress tensor
    DO i = 1,3 
      DO j = 1,3
        ReynoldsTensor(i,j) = (USGSPart(i+1,iPart)/FieldAtParticle(1,iPart))*(USGSPart(j+1,iPart)/FieldAtParticle(1,iPart))
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
	
	Lambda_factor = (3 * trace_HR) / (2 * kSGSPart(iPart) * trace_H)
!PRINT *, 'Lambda_factor:  ', Lambda_factor	

	! Calculate simmetrix matrix D_ij for Cholewski decomposition
	DO i = 1,3
          DO j = 1,3
            IF (i.EQ.j) THEN
              D_matrix(i,j,iPart) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j) -1)/3
            ELSE
              D_matrix(i,j,iPart) = C_0*Lambda_factor*H_matrix(i,j) + 2 * ( Lambda_factor*H_matrix(i,j)   )/3
            END IF
          END DO
        END DO
	D_matrix(:,:,iPart) = epsilonSGS(iPart) * D_matrix(:,:,iPart)
!PRINT *, 'D_matrix:  ', D_matrix(1:3,1:3,iPart)	
	
	! Cholewski decomposition to get diffusion matrix B_matrix
	B_matrix(:,:,iPart) = 0.	
	
	B_matrix(1,1,iPart) = sqrt( D_matrix(1,1,iPart) ) 				! element (1,1)
	DO i = 2 , 3									
	  B_matrix(i,1,iPart) = D_matrix(1,i,iPart) / B_matrix(1,1,iPart)		! first column
	END DO

	DO j = 2 , 3
	  DO i = 1 , 3
	    IF( i .LT. j ) THEN
	      B_matrix(i,j,iPart) = 0.0   	! upper triangolar elements
  	    ELSE
	      IF( i .EQ. j) THEN
		aux = 0.0
		DO k = 1 , i -1
		  aux = aux + B_matrix(i,k,iPart)**2
		END DO
		B_matrix(i,j,iPart) = SQRT( D_matrix(i,i,iPart) - aux ) ! diagonal elements
	      ELSE ! (i .GT. j)
		aux = 0.0
		DO k = 1 , i - 1
		  aux = aux + B_matrix(j,k,iPart) * B_matrix(i,k,iPart)
		END DO
	        B_matrix(i,j,iPart) = ( D_matrix(j,i,iPart) - aux ) / B_matrix(j,j,iPart)    ! in our case only element B(3,2)
	      END IF
            END IF
	  END DO
	END DO

!PRINT *, 'B_matrix:  ', B_matrix(1:3,1:3,iPart)
	
	! Check FOR THE CHOLEWSKI DECOMPOSITION
	!B_transpose = transpose(B_matrix)
	!D_check = matmul(B_matrix,B_transpose)
	!tolerance = 10e-4
	!DO i = 1,3
	!  DO j = 1,3
	!    aux = D_check(i,j,iPart) - D_matrix(i,j,iPart)
	!  	 IF (aux .GT. tolerance) THEN
	!    	PRINT *, 'ERROR IN CHOLEWSKI DECOMPOSITION '
	!    END IF 
  	!  END DO
  	!END DO
  	
    END IF ! 3RD END IF

    ! EULER
    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3	
      !ddrift(1:3)    = G_matrix(1:3,j,iPart)*fluctuationState(j,iPart)*dt
      !diffusion(1:3) = B_matrix(1:3,j,iPart)*SQRT(dt)!*RandNormal()
      !Pt(1:3) = Pt(1:3)  + ddrift(1:3) + diffusion(1:3)
      Pt(1:3) = Pt(1:3)  + G_matrix(1:3,j,iPart)*fluctuationState(j,iPart)*dt + B_matrix(1:3,j,iPart)*SQRT(dt)!*RandNormal()
    END DO					!CHECK FLUCTUATION PART
  END IF	! 1ST END IF
  
  ! FULL FLUID VELOCITY SEEN BY THE PARTICLES

  pg(1:3)  = pressur_gradient(1:3,iPart)*dt ! 1/\rho_f * ( \partial <P> / \partial x_i ) * dt 
  sdt(1:3) = second_drift_term(1:3,iPart)*dt !  (U_p-U)*(\partial U_i / \partial x_j) * dt

  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) - pg(1:3) + sdt(1:3) + Pt(1:3)

!PRINT *, 'pg:  ', pg(1:3)
!PRINT *, 'sdt: ', sdt(1:3)
!PRINT *, 'ddrift:  ', ddrift(1:3)
!PRINT *, 'diffusion:  ', diffusion(1:3)
!PRINT *, 'Pt:  ', Pt(1:3)
!PRINT *, 'TurbPartState:  ', TurbPartState(1:3,iPart)
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


  !filterwidth = 2*PI/32    	! to validate the implementation
  !USGSPart(2:4,iPart) = 0.1	! to validate the implementation
  !FieldAtParticle(1,iPart) = 1.	! to validate the implementation
!PRINT *, 'filterwidth:  ', filterwidth
!PRINT *, 'USGSPart:  ', USGSPart(2:4,iPart)

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)
!PRINT *, 'kSGSPart:  ', kSGSPart(iPart)

  ! Deviation standard of fluctuation velocity
  sigmaF(iPart) = SQRT(2./3.*kSGSPart(iPart))
!PRINT *, 'sigmaF:  ', sigmaF(iPart)
 
  ! Deviation standard per direction
  sigmaF_vector_square(1:3,iPart) = (USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2. !considering sigmaF_vector already to the power of 2
!PRINT *, 'sigmaF_vector_square:  ', sigmaF_vector_square(1:3,iPart)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Calculate dissipation energy according to  Babak Shotorban & Farzad Mashayek (2006) A stochastic model for particle
  ! motion in large-eddy simulation, Journal of Turbulence, 7, N18, DOI: 10.1080/14685240600595685
  epsilonSGS(iPart) = C_eps_Shotorban/ElemVolN(ElemID) * (kSGSPart(iPart)**(3./2)) 
  !epsilonSGS(iPart) = C_eps_Shotorban/filterwidth * (kSGSPart(iPart)**(3./2))  ! to validate the implementation
!PRINT *, 'epsilonSGS:  ', epsilonSGS(iPart)
 
  ! Relative velocity
  udiff(1:3) =  (FieldAtParticle(2:4,iPart) + TurbPartState(1:3,iPart)) - PartState(4:6,iPart) 
	!udiff(1:3) = 1.	! to validate the implementation
	
  ! Relative displcement between particle and fluid element
  deltaR(1:3,iPart) = udiff(1:3)*dt   ! relative displacement in the three directions
  DELTA_R(iPart) = VECNORM(deltaR(1:3,iPart)) ! |U-U_p|*dt (scalar value)
	
  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaF(iPart))) THEN 	!1st IF
    ! We ASSUME that these are the CORRECT matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN 	!2nd IF
      Pt(1:3) = 0 ! in this case we have 0/0 
    ELSE	! 2nd 	ELSE IF
      
	! Lagrangian correlation 
	R_L(1:3,iPart) = 1.
			
	! Integral length scales L_E when sigmaF goes to zero 
	L_E(iPart) = C_L**2./C_eps_Shotorban * SQRT(8./27) * ElemVolN(ElemID)

	! Function g(dr) and f(dr)
	F_dr(iPart) = EXP(-DELTA_R(iPart)/L_E(iPart))
	G_dr(iPart) = (1 - DELTA_R(iPart) / (2*L_E(iPart)) ) *  EXP(-DELTA_R(iPart)/L_E(iPart))

		! Eulerian  correlation 
!	DO i = 1,3
!          DO j = 1,3
!            IF (i.EQ.j) THEN
!              R_E(i,j,iPart) = G_dr(iPart) + (F_dr(iPart) - G_dr(iPart))*(deltaR(i,iPart)*deltaR(j,iPart))/(DELTA_R(iPart)**2)
!            ELSE
!              R_E(i,j,iPart) = (F_dr(iPart) - G_dr(iPart))*(deltaR(i,iPart)*deltaR(j,iPart))/(DELTA_R(iPart)**2)
!            END IF
!          END DO
!        END DO

	DO i = 1,3 ! I'm considering only the main diagonal of the tensor (treat it as a  vector)
	   diag_R_E(i,iPart) = G_dr(iPart) + (F_dr(iPart) - G_dr(iPart))*(deltaR(i,iPart)*deltaR(i,iPart))/(DELTA_R(iPart)**2)
	END DO
			
        ! Identity matrix
!        DO i = 1,3
!          DO j = 1,3
!            M_ones(i,j) = 1.
!          END DO
!        END DO
     	 
       	!Correlaction function: vector R_P
	DO i = 1,3
	  R_P(i,iPart) = R_L(i,iPart) * diag_R_E(i,iPart)
	END DO 
      
    ! Calculate the fluctuation in the three directions
    DO i = 1,3
      Pt(1:3) = R_P(i,iPart)*TurbPartState(i,iPart) + SQRT(sigmaF_vector_square(i,iPart))*SQRT(1-R_P(i,iPart)**2) !* RandNormal()  !!!* SQRT(dt)
    END DO

    END IF  !2nd END IF 

  ! Valid SGS turbulent kinetic energy
  ELSE   	!1st 	ELSE IF
   
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN  	! 3rd IF
     Pt(1:3) = 0 ! in this case we have 0/0 
    ELSE	! 3rd 	ELSE IF
      	! Lagrangian time scale in the appropriate directions
	T_L(1:3,iPart) = C_T * sigmaF_vector_square(1:3,iPart)/epsilonSGS(iPart)
!PRINT *, 'T_L:  ', T_L(1:3,iPart)

	! Lagrangian correlation 
	R_L(1:3,iPart) = EXP(-dt/T_L(1:3,iPart))
!PRINT *, 'R_L:  ', R_L(1:3,iPart)
	
	! Integral length scales L_E
	T_L_f(iPart) = C_T *  (sigmaF(iPart)**2) / epsilonSGS(iPart)
	L_E(iPart) = C_L * T_L_f(iPart) * sigmaF(iPart)
!PRINT *, 'T_L_f:  ', T_L_f(iPart)
!PRINT *, 'L_E:  ', L_E(iPart)

	! Function g(dr) and f(dr)
	F_dr(iPart) = EXP(-DELTA_R(iPart)/L_E(iPart))
	G_dr(iPart) = (1 - DELTA_R(iPart) / (2*L_E(iPart)) ) *  EXP(-DELTA_R(iPart)/L_E(iPart))
!PRINT *, 'F_dr:  ', F_dr(iPart)	
!PRINT *, 'G_dr:  ', G_dr(iPart)
	
	! Eulerian  correlation 
!	DO i = 1,3
!          DO j = 1,3
!            IF (i.EQ.j) THEN
!              R_E(i,j,iPart) = G_dr(iPart) + (F_dr(iPart) - G_dr(iPart))*(deltaR(i,iPart)*deltaR(j,iPart))/(DELTA_R(iPart)**2)
!            ELSE
!              R_E(i,j,iPart) = 	      (F_dr(iPart) - G_dr(iPart))*(deltaR(i,iPart)*deltaR(j,iPart))/(DELTA_R(iPart)**2)
!            END IF
!          END DO
!        END DO

	DO i = 1,3 ! I'm considering only the main diagonal of the tensor (treating it as a  vector)
	   diag_R_E(i,iPart) = G_dr(iPart) + (F_dr(iPart) - G_dr(iPart))*(deltaR(i,iPart)*deltaR(i,iPart))/(DELTA_R(iPart)**2)
	END DO
!PRINT *, 'diag_R_E:  ', diag_R_E(1:3,iPart)				
      
       ! Identity matrix
!       DO i = 1,3
!          DO j = 1,3
!            M_ones(i,j) = 1.
!          END DO
!       END DO
     	 
       	!Correlaction function: vector R_P
	DO i = 1,3
	  R_P(i,iPart) = R_L(i,iPart) * diag_R_E(i,iPart)
	END DO 
!PRINT *, 'R_P:  ', R_P(1:3,iPart)	      

    
    ! Calculate the fluctuation in the three directions
    DO i = 1,3
      Pt(1:3) = R_P(i,iPart)*TurbPartState(i,iPart) + SQRT(sigmaF_vector_square(i,iPart))*SQRT(1-R_P(i,iPart)**2)* RandNormal()  !!!* SQRT(dt)
    END DO
      
    END IF ! 3rd END IF

    
  END IF ! 1ST END IF 
  
  TurbPartState(1:3,iPart) =  Pt(1:3)

!PRINT *, 'TurbPartState:  ', TurbPartState(1:3,iPart)

  END DO  ! For Sommerfeld


CASE('Amiri-Analytic')
!===================================================================================================================================
! A. Elhami Amiri , S. Kazemzadeh Hannani & F. Mashayek (2006)
!>  Large-Eddy Simulation of Heavy-Particle Transport in Turbulent Channel Flow, Numerical Heat Transfer, 
!> Part B: Fundamentals, 50:4, 285-313, DOI: 10.1080/10407790600859577

!===================================================================================================================================

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF 

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

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

   ! Setting fixed values of USGSpart in order to validate the code 
  USGSPart(2:4,iPart)  = 0.01 ! Fluctuation in order to validate the code
  PRINT *, 'USGSPart:  ', USGSPart(2:4,ipart)

  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)
PRINT *, 'kSGSPart:  ', kSGSPart(ipart)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem
  
  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(kSGSPart(iPart))) THEN
    
	!in the limit of sigmaSGS goes to zero the time sigmaS goes to zero 
	ni_SGS(1:3,iPart) = 0.

     ! Valid SGS turbulent kinetic energy
     ELSE
	
	! Particle relaxation time
	taup_Amiri(iPart)    = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0)		
PRINT *, 'taup_Amiri:  ', taup_Amiri(ipart)

	! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
		! Local vorticity 
		rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
		rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
		rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy 		

	!T_f_p(iPart) =1./VECNORM(rot_U)
	T_f_p(iPart) = 1.
PRINT *, 'T_f_p:  ', T_f_p(ipart)

	! Calculate the factor f1, f2, f3 for the anisotropy
	factor(1) = SQRT(USGSPart(2,iPart)**2) / SQRT( 2./3 * kSGSPart(iPart) )   ! f1 
	factor(2) = SQRT(USGSPart(3,iPart)**2) / SQRT( 2./3 * kSGSPart(iPart) )   ! f2
	factor(3) = SQRT(3 - factor(1)**2 - factor(2)**2)			  ! f3
PRINT *, 'factor:  ', factor(1:3)

	! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity 
	! due to SGS velocity in time dt
 	U_square(1:3,iPart) = factor(1:3)**2 * 2./3 * kSGSPart(iPart)  ! U_square = (u'')^2
PRINT *, 'U_square:  ', U_square(1:3,iPart)

	!alpha and beta Amiri
	alpha_Amiri(ipart) =  dt/taup_Amiri(iPart)
	beta_Amiri(iPart)  = taup_Amiri(iPart)/T_f_p(iPart) 
	beta1 = 1+beta_Amiri(iPart)
	beta2 = 1-beta_Amiri(iPart)
PRINT *, 'alpha_Amiri:  ', alpha_Amiri(ipart)
PRINT *, 'beta_Amiri:  ', beta_Amiri(ipart)

	! sigmaS_Amiri is the increase of standard deviation of particle velocity due to SGS velocity during time dt
	
	sigmaS_Amiri(1:3,iPart) = SQRT(U_square(1:3,iPart)) * SQRT(	1/beta1 * (1-EXP(-alpha_Amiri(ipart)*beta1)) 	&
				&	-1/beta2 * EXP(-2*alpha_Amiri(ipart)) * (1-EXP(alpha_Amiri(ipart)*beta2))		)
PRINT *, 'sigmaS_Amiri:  ', sigmaS_Amiri(1:3,iPart)

	! SGS fluid velocity by Amiri
	ni_SGS_Amiri(1:3,iPart) =  sigmaS_Amiri(1:3,iPart)/dt !* RandNormal()
  
    END IF
	TurbPartState(1:3,iPart) =  ni_SGS_Amiri(1:3,iPart)
PRINT *, 'ni_SGS_Amiri:  ', ni_SGS_Amiri(1:3,iPart)

END DO ! AMIRI-ANALYTIC SGS MODEL


CASE('Amiri')
!===================================================================================================================================
! A. Elhami Amiri , S. Kazemzadeh Hannani & F. Mashayek (2006)
!>  Large-Eddy Simulation of Heavy-Particle Transport in Turbulent Channel Flow, Numerical Heat Transfer, 
!> Part B: Fundamentals, 50:4, 285-313, DOI: 10.1080/10407790600859577

!===================================================================================================================================

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF 

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

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
  kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem
  
  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(kSGSPart(iPart))) THEN
    
	!in the limit of sigmaSGS goes to zero the time sigmaS goes to zero 
	ni_SGS(1:3,iPart) = 0.

     ! Valid SGS turbulent kinetic energy
     ELSE
	
	! Particle relaxation time
	taup_Amiri(iPart)    = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0)		

	! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
		! Local vorticity 
		rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
		rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
		rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy 		

	T_f_p(iPart) =1./VECNORM(rot_U)

	! Calculate the factor f1, f2, f3 for the anisotropy
	factor(1) = SQRT((USGSPart(2,iPart)**2)/(FieldAtParticle(1,iPart))) / SQRT( 2./3 * kSGSPart(iPart) )   ! f1 
	factor(2) = SQRT((USGSPart(3,iPart)**2)/(FieldAtParticle(1,iPart))) / SQRT( 2./3 * kSGSPart(iPart) )   ! f2
	factor(3) = SQRT(3 - factor(1)**2 - factor(2)**2)		          ! f3
 

	! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity 
	! due to SGS velocity in time dt
 	U_square(1:3,iPart) = factor(1:3)**2 * 2./3 * kSGSPart(iPart)  ! U_square = (u'')^2

	!alpha and beta Amiri
	alpha_Amiri(ipart) =  dt/taup_Amiri(iPart)
	beta_Amiri(iPart)  = taup_Amiri(iPart)/T_f_p(iPart) 
	beta1 = 1+beta_Amiri(iPart)
	beta2 = 1-beta_Amiri(iPart)

	! sigmaS_Amiri is the increase of standard deviation of particle velocity due to SGS velocity during time dt
	
	sigmaS_Amiri(1:3,iPart) = SQRT(U_square(1:3,iPart)) * SQRT(	1/beta1 * (1-EXP(-alpha_Amiri(ipart)*beta1)) 	&
				&	-1/beta2 * EXP(-2*alpha_Amiri(ipart)) * (1-EXP(alpha_Amiri(ipart)*beta2))		)

	! SGS fluid velocity by Amiri
	ni_SGS_Amiri(1:3,iPart) =  sigmaS_Amiri(1:3,iPart)/dt * RandNormal()
  
    END IF
	TurbPartState(1:3,iPart) =  ni_SGS_Amiri(1:3,iPart)

END DO ! AMIRI SGS MODEL


CASE('Fukagata-Analytic') ! CASE USED IN ORDER TO VALIDATE THE CODE
!===================================================================================================================================
! K. Fukagata (2004)
!> Dynamics of Brownian particles in a turbulent channel flow, K. Fukagata, S. Zahrai, F. H. Bark
!> Heat and Mass Transfer 40 (2004) 715–726. DOI 10.1007/s00231-003-0462-8
!===================================================================================================================================

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

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

CALL InterpolateFieldToParticle(3,gradUx(1:3,:,:,:,:),3,gradUxPart) 
CALL InterpolateFieldToParticle(3,gradUy(1:3,:,:,:,:),3,gradUyPart)
CALL InterpolateFieldToParticle(3,gradUz(1:3,:,:,:,:),3,gradUzPart)

filterwidth = 2*PI/32  ! validation of the implementation


DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

PRINT *, 'filterwidth:  ', filterwidth

	! Particle relaxation time
	taup_Fuka(iPart)    = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0)			
PRINT *, 'taup_Fuka:  ', taup_Fuka(iPart)

	! Integral time scale of fluid velocity along a particle trajectory
		rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
		rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
		rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy 		

	!T_f_p(iPart) = 1./VECNORM(rot_U)
	T_f_p(iPart) =1
PRINT *, 'T_f_p:  ', T_f_p(iPart)

	! Calculation of deformation tensor d_ij

	d_ij_tensor(1,:)= 0.5*(/2*gradUxPart(1,iPart), gradUyPart(1,iPart)+gradUxPart(2,iPart),	gradUzPart(1,iPart)+gradUxPart(3,iPart)/) 	
	
	d_ij_tensor(2,:)= 0.5*(/gradUyPart(1,iPart)+gradUxPart(2,iPart), 2*gradUyPart(2,iPart),	gradUzPart(2,iPart)+gradUyPart(3,iPart)/) 
	
	d_ij_tensor(3,:)= 0.5*(/gradUzPart(1,iPart)+gradUxPart(3,iPart), gradUzPart(2,iPart)+gradUyPart(3,iPart), 2*gradUzPart(3,iPart)/) 

	! Calculation of norm of tensor d_ij according to F. Gurniki, K. Fukagata, S. Zahrai and F. H. Bark (2000)
      d_ij = 0 ! inizialization
      DO i = 1,3
        DO j = 1,3
          	d_ij= d_ij + d_ij_tensor(i,j)**2
        END DO
      END DO
	d_ij = SQRT(d_ij)
	d_ij = 1
PRINT *, 'd_ij:  ', d_ij

	! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity due to SGS velocity in time dt
 	U_M(iPart) = 2/(3**0.5) * C_Fuka * C_eps * filterwidth * d_ij	
PRINT *, 'U_M:  ', U_M(iPart)
	
!alpha and beta Fukagata
	alpha_Fuka(ipart) =  dt/taup_Fuka(iPart)
	beta_Fuka(iPart)  = taup_Fuka(iPart)/T_f_p(iPart) 
	beta1 = 1+beta_Fuka(iPart)
	beta2 = 1-beta_Fuka(iPart)
PRINT *, 'alpha_Fuka:  ', alpha_Fuka(iPart)
PRINT *, 'beta_Fuka:  ', beta_Fuka(iPart)

	! sigmaS_Fuka is the increase of standard deviation of particle velocity due to SGS velocity during time dt
	sigmaS_Fuka(iPart) = U_M(iPart) * SQRT(	1/beta1 * (1-EXP(-alpha_Fuka(ipart)*beta1)) 	&
				        &      -1/beta2 * EXP(-2*alpha_Fuka(ipart)) * (1-EXP(alpha_Fuka(ipart)*beta2))		)
PRINT *, 'sigmaS_Fuka:  ', sigmaS_Fuka(iPart)

	! SGS fluid velocity by Fukagata
	DO j = 1,3
	ni_SGS(1:3,iPart) = sigmaS_Fuka(iPart)/dt !* RandNormal()
	END DO 
PRINT *, 'ni_SGS:  ', ni_SGS(1:3,iPart)

    	! sum up the two terms 
	ni_Fuka(1:3,iPart) = ni_SGS(1:3,iPart)

	TurbPartState(1:3,iPart) =  ni_Fuka(1:3,iPart)
PRINT *, 'TurbPartState:  ', TurbPartState(1:3,iPart)
END DO ! Fukagata-Analytic


CASE('Fukagata')
!===================================================================================================================================
! K. Fukagata (2004)
!> Dynamics of Brownian particles in a turbulent channel flow, K. Fukagata, S. Zahrai, F. H. Bark
!> Heat and Mass Transfer 40 (2004) 715–726. DOI 10.1007/s00231-003-0462-8
!===================================================================================================================================

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

! Time integration using Euler-Maruyama scheme
IF(PRESENT(iStage))THEN
IF (iStage.NE.1) RETURN
END IF 

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

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
 kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)
	
  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  
  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(kSGSPart(iPart))) THEN
    
	! Particle relaxation time
	taup_Fuka(iPart)    = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0)	

	! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
		! Local vorticity 
		rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
		rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
		rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy 		

	T_f_p(iPart) = 1./VECNORM(rot_U)

	! Calculation of deformation tensor d_ij
	d_ij_tensor(1,:)= 0.5*(/2*gradUxPart(1,iPart), gradUyPart(1,iPart)+gradUxPart(2,iPart),	gradUzPart(1,iPart)+gradUxPart(3,iPart)/) 	
	
	d_ij_tensor(2,:)= 0.5*(/gradUyPart(1,iPart)+gradUxPart(2,iPart), 2*gradUyPart(2,iPart),	gradUzPart(2,iPart)+gradUyPart(3,iPart)/) 
	
	d_ij_tensor(3,:)= 0.5*(/gradUzPart(1,iPart)+gradUxPart(3,iPart), gradUzPart(2,iPart)+gradUyPart(3,iPart), 2*gradUzPart(3,iPart)/) 

	! Calculation of norm of tensor d_ij according to F. Gurniki, K. Fukagata, S. Zahrai and F. H. Bark (2000)
      d_ij = 0 ! inizialization
      DO i = 1,3
        DO j = 1,3
          	d_ij= d_ij + d_ij_tensor(i,j)**2
        END DO
      END DO
	d_ij = SQRT(d_ij)
	
	! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity due to SGS velocity in time dt
 	U_M(iPart) = 2/(3**0.5) * C_Fuka * C_eps * ElemVolN(ElemID) * d_ij	

	!alpha and beta Fukagata
	alpha_Fuka(ipart) =  dt/taup_Fuka(iPart)
	beta_Fuka(iPart)  = taup_Fuka(iPart)/T_f_p(iPart) 
	beta1 = 1+beta_Fuka(iPart)
	beta2 = 1-beta_Fuka(iPart)

	! sigmaS_Fuka is the increase of standard deviation of particle velocity due to SGS velocity during time dt
	sigmaS_Fuka(iPart) = U_M(iPart) * SQRT(	1/beta1 * (1-EXP(-alpha_Fuka(ipart)*beta1)) 	&
				 &	-1/beta2 * EXP(-2*alpha_Fuka(ipart)) * (1-EXP(alpha_Fuka(ipart)*beta2))		)

	! SGS fluid velocity by Fukagata
	DO j = 1,3
	ni_SGS(1:3,iPart) = sigmaS_Fuka(iPart)/dt * RandNormal()
	END DO 

    	! sum up the two terms 
	ni_Fuka(1:3,iPart) = ni_SGS(1:3,iPart)


     ! Valid SGS turbulent kinetic energy
     ELSE
	
	! Particle relaxation time
	taup_Fuka(iPart)    = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0)		

	! Schmidt number of the particles in fluid
!	Sc_Fuka(iPart) = Species(PartSpecies(iPart))%MassIC*nu / ( taup_Fuka(iPart) * kSGSPart(iPart) * FieldAtParticle(TEMP,iPart) ) 

	! Brownian force 
!	ni_B(iPart) = RandNormal()*SQRT(2*nu/(Sc_Fuka(iPart) * taup_Fuka(iPart)**2 * dt))

	! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
		! Local vorticity 
		rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
		rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
		rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy 		

	T_f_p(iPart) =1./VECNORM(rot_U)

	! Calculation of deformation tensor d_ij
	d_ij_tensor(1,:)= 0.5*(/2*gradUxPart(1,iPart), gradUyPart(1,iPart)+gradUxPart(2,iPart),	gradUzPart(1,iPart)+gradUxPart(3,iPart)/) 	
	
	d_ij_tensor(2,:)= 0.5*(/gradUyPart(1,iPart)+gradUxPart(2,iPart), 2*gradUyPart(2,iPart),	gradUzPart(2,iPart)+gradUyPart(3,iPart)/) 
	
	d_ij_tensor(3,:)= 0.5*(/gradUzPart(1,iPart)+gradUxPart(3,iPart), gradUzPart(2,iPart)+gradUyPart(3,iPart), 2*gradUzPart(3,iPart)/) 

	! Calculation of norm of tensor d_ij according to F. Gurniki, K. Fukagata, S. Zahrai and F. H. Bark (2000)
     	 d_ij = 0 
     	 DO i = 1,3
      	   DO j = 1,3
          	d_ij= d_ij + d_ij_tensor(i,j)**2
      	   END DO
    	 END DO
	d_ij = SQRT(d_ij)

	
	! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity due to SGS velocity in time dt
 	U_M(iPart) = 2/(3**0.5) * C_Fuka * C_eps * ElemVolN(ElemID) * d_ij	
	

	!alpha and beta Fukagata
	alpha_Fuka(ipart) =  dt/taup_Fuka(iPart)
	beta_Fuka(iPart)  = taup_Fuka(iPart)/T_f_p(iPart) 
	beta1 = 1+beta_Fuka(iPart)
	beta2 = 1-beta_Fuka(iPart)

	!sigmaS_Fuka is the increase of standard deviation of particle velocity due to SGS velocity during time dt
	
	sigmaS_Fuka(iPart) = U_M(iPart) * SQRT(	1/beta1 * (1-EXP(-alpha_Fuka(ipart)*beta1)) 	&
				&	-1/beta2 * EXP(-2*alpha_Fuka(ipart)) * (1-EXP(alpha_Fuka(ipart)*beta2))		)
!PRINT *, 'sigmaS_Fuka:  ', sigmaS_Fuka

	! SGS fluid velocity by Fukagata
	DO j = 1,3
	ni_SGS(1:3,iPart) =  sigmaS_Fuka(iPart)/dt * RandNormal()
	END DO  

    	! sum up the two terms
	!!!!!! now I'm using only te SGS fluid velocity and NO the brownian force
	ni_Fuka(1:3,iPart) = ni_SGS(1:3,iPart) !+ ni_B(iPart)

  END IF
	! TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)
	TurbPartState(1:3,iPart) =  ni_Fuka(1:3,iPart)
PRINT *, 'TurbPartState:  ', TurbPartState(:,iPart)
END DO ! Fukagata


CASE('Jin-Analytic') ! CASE USED IN ORDER TO VALIDATE THE CODE
!===================================================================================================================================
! Jin (2010)
!> Guodong Jin a , Guo-Wei He a,*, Lian-Ping Wang b, Jian Zhang a Subgrid scale fluid velocity timescales seen by inertial 
!  particles in large-eddy simulation of particle-laden turbulence. International Journal of Multiphase Flow 36 (2010) 432–437
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
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

filterwidth = 2*PI/32  ! validation of the implementation

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

PRINT *, 'filterwidth:  ', filterwidth
    ! Use unfiltered density to obtain velocity in primitive variables
    ! kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)
    kSGSPart(iPart)= 3./2! validation of the implementation
PRINT *, 'kSGSPart:  ', kSGSPart(iPart)

  ! Time scale of SGS scales = U_rms in case of isotropic turbolence
   sigmaSGS(iPart) = SQRT(2./3*kSGSPart(iPart))
PRINT *, 'sigmaSGS:  ', sigmaSGS(iPart)

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS(iPart))) THEN
    
    	!in the limit of sigmaSGS goes to zero the time scale goes to infinity
	Pt(1:3) = 0.

     ! Valid SGS turbulent kinetic energy
     ELSE
		
	! Lagrangian time scale 
	!deltaT_L(iPart)   = C*ElemVolN(ElemID)/sigmaSGS(iPart)
	 deltaT_L(iPart) = C*filterwidth/sigmaSGS(iPart) ! validation of the implementation
PRINT *, 'deltaT_L:  ', deltaT_L(iPart)

	! cutoff wave number 
	!K_cf(iPart) = PI/ElemVolN(ElemID)
	 K_cf(iPart) = PI/filterwidth  ! validation of the implementation
PRINT *, 'K_cf:  ', K_cf(iPart)

	!  SGS Eulerian integral timescale according to Fede and Simonin, 2006
	 deltaT_E(iPart) = 3*PI/10*1./(K_cf(iPart)*sigmaSGS(iPart))
PRINT *, 'deltaT_E:  ', deltaT_E(iPart)	

	! ratio between lagrangian and eulerian time scale
	beta_Jin(iPart) = deltaT_L(iPart)/deltaT_E(iPart)
PRINT *, 'beta_Jin:  ', beta_Jin(iPart)

    	! Particle relaxation time
	taup_Jin(iPart)    = Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0) 
PRINT *, 'taup_Jin:  ', taup_Jin(iPart)
   
    	! Stokes number 
	St_Jin(iPart) = taup_Jin(iPart)/deltaT_E(iPart)
PRINT *, 'St_Jin:  ', St_Jin(iPart)

	! SGS autocorrelation time seen by a particle =  deltaT_L_p
	deltaT_L_p(iPart) =  deltaT_L(iPart)/beta_Jin(iPart) * &
				& ( 1 - ( (1-beta_Jin(iPart)) / ( 1+St_Jin(iPart))**(0.4*(1+0.01*St_Jin(iPart)))  )	)
PRINT *, 'deltaT_L_p:  ', deltaT_L_p(iPart)

    ! EULER
    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
	!dxminus(1:3) = - 1./deltaT_L_p(iPart)*TurbPartState(j,iPart)*dt 
	!dxplus(1:3)  = SQRT( (4*kSGSPart(iPart)*dt)/(deltaT_L_p(iPart)*3) )
      Pt(1:3) = Pt(1:3) + - 1./deltaT_L_p(iPart)*TurbPartState(j,iPart)*dt + SQRT( (4*kSGSPart(iPart)*dt)/(deltaT_L_p(iPart)*3) )!*RandNormal()
    END DO

  END IF

  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)

!PRINT *, 'dxminus:  ', dxminus(1:3)
!PRINT *, 'dxplus:  ', dxplus(1:3)
PRINT *, 'Pt:  ', Pt(1:3)
PRINT *, 'TurbPartState:  ', TurbPartState(1:3,iPart)

END DO ! for Jin-Analytic


CASE('Jin')
!===================================================================================================================================
! Jin (2010)
!> Guodong Jin a , Guo-Wei He a,*, Lian-Ping Wang b, Jian Zhang a Subgrid scale fluid velocity timescales seen by inertial 
!  particles in large-eddy simulation of particle-laden turbulence. International Journal of Multiphase Flow 36 (2010) 432–437
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
USGS = U - USGS

! Interpolate SGS kinetic energy to particle position
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),4,USGSPart)

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  
  ! Use unfiltered density to obtain velocity in primitive variables
   kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)
   
  ! Time scale of SGS scales = U_rms in case of isotropic turbolence
   sigmaSGS(iPart) = SQRT(2./3*kSGSPart(iPart))

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS(iPart))) THEN
    
    	!in the limit of sigmaSGS goes to zero the time scale goes to infinity
	Pt(1:3) = 0.

     ! Valid SGS turbulent kinetic energy
     ELSE
		
	! Lagrangian time scale 
	deltaT_L(iPart)   = C*ElemVolN(ElemID)/sigmaSGS(iPart)
	 
	! cutoff wave number 
	K_cf(iPart) = PI/ElemVolN(ElemID)
	 
	!  SGS Eulerian integral timescale according to Fede and Simonin, 2006
	deltaT_E(iPart) = 3*PI/10*1./(K_cf(iPart)*sigmaSGS(iPart))

	! ratio between lagrangian and eulerian time scale
	beta_Jin(iPart) = deltaT_L(iPart)/deltaT_E(iPart)

    	! Particle relaxation time
	taup_Jin(iPart)    = Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0) 
   
    	! Stokes number 
	St_Jin(iPart) = taup_Jin(iPart)/deltaT_E(iPart)

	! SGS autocorrelation time seen by a particle =  deltaT_L_p
	deltaT_L_p(iPart) =  deltaT_L(iPart)/beta_Jin(iPart) * &
				& ( 1 - ( (1-beta_Jin(iPart)) / ( 1+St_Jin(iPart))**(0.4*(1+0.01*St_Jin(iPart)))  )	)

	! JIN simulation results can be fitted by the next empirical curve:
	!deltaT_L_p(iPart) =  deltaT_L(iPart)/beta_Jin(iPart) * ((0.444-0.7*ElemVolN(ElemID))*EXP(-(LOG(St_Jin(iPart)/0.5))**2) &
	!& + 1 - (1-beta_Jin(iPart))*EXP(-St_Jin(iPart)/5.15))
        !deltaT_L_p(iPart) =  deltaT_L(iPart)/beta_Jin(iPart) * ((0.444-0.7*filterwidth)*EXP(-(LOG(St_Jin(iPart)/0.5))**2) &
	!& + 1 - (1-beta_Jin(iPart))*EXP(-St_Jin(iPart)/5.15)) ! validation of the implementation

    ! EULER
    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3) - 1./deltaT_L_p(iPart)*TurbPartState(j,iPart)*dt + SQRT( (4*kSGSPart(iPart)*dt)/(deltaT_L_p(iPart)*3) ) * RandNormal()
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
  kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Time scale of SGS scales
  sigmaSGS(iPart) = SQRT(2./3.*kSGSPart(iPart))

  ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  ! Relative velocity
  udiff(1:3) = PartState(4:6,iPart) - (FieldAtParticle(2:4,iPart) + TurbPartState(1:3,iPart))

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/SQRT(SUM(udiff**2))
  ELSE
    urel = 0.
  END IF

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS(iPart))) THEN
    ! We ASSUME that these are the CORRECT matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      RETURN
    ELSE
      ! parallel
      tauL(1,iPart) = C*ElemVolN(ElemID)/(betaSGS*SQRT(SUM(udiff**2)))
      ! perpendicular
      tauL(2,iPart) = C*ElemVolN(ElemID)/(betaSGS*2*SQRT(SUM(udiff**2)))

      ! Calculate drift and diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            G_SGS(i,j,iPart) = 1/     tauL(2,iPart)  + (1/     tauL(1,iPart)  - 1/     tauL(2,iPart)) *urel(i)*urel(j)
          ELSE
            G_SGS(i,j,iPart) =                         (1/     tauL(1,iPart)  - 1/     tauL(2,iPart)) *urel(i)*urel(j)
          END IF
        END DO
      END DO
      ! as sigma is ALMOSTZERO
      B_SGS(:,:,iPart) = 0.

      ! EULER
      ! Sum up turbulent contributions
      Pt(1:3) = 0.
      DO j = 1,3
        Pt(1:3) = Pt(1:3) - G_SGS(1:3,j,iPart)*TurbPartState(j,iPart)*dt
      END DO
    END IF

  ! Valid SGS turbulent kinetic energy
  ELSE

    tauSGS   = C*ElemVolN(ElemID)/sigmaSGS(iPart)

    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      ! parallel
      tauL(1,iPart) = tauSGS(iPart)
      ! perpendicular
      tauL(2,iPart) = tauSGS(iPart)

      ! Calculate drift and diffusion matrix
      G_SGS(:,:,iPart)=0.
      B_SGS(:,:,iPart)=0.
      DO i = 1,3
        G_SGS(i,i,iPart) = 1/     tauL(2,iPart)
        B_SGS(i,i,iPart) = 1/SQRT(tauL(2,iPart))
      END DO
      B_SGS(:,:,iPart) = SQRT(2*sigmaSGS(iPart)**2)*B_SGS(:,:,iPart)

    ELSE
      ! parallel
      tauL(1,iPart) = tauSGS(iPart)/(SQRT(1+  betaSGS**2*SUM(udiff**2)/kSGSPart(iPart)*3/2))
      ! perpendicular
      tauL(2,iPart) = tauSGS(iPart)/(SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart(iPart)*3/2))

      ! Calculate drift and diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            G_SGS(i,j,iPart) = 1/     tauL(2,iPart)  + (1/     tauL(1,iPart)  - 1/     tauL(2,iPart)) *urel(i)*urel(j)
            B_SGS(i,j,iPart) = 1/SQRT(tauL(2,iPart)) + (1/SQRT(tauL(1,iPart)) - 1/SQRT(tauL(2,iPart)))*urel(i)*urel(j)
          ELSE
            G_SGS(i,j,iPart) =                         (1/     tauL(1,iPart)  - 1/     tauL(2,iPart)) *urel(i)*urel(j)
            B_SGS(i,j,iPart) =                         (1/SQRT(tauL(1,iPart)) - 1/SQRT(tauL(2,iPart)))*urel(i)*urel(j)
          END IF
        END DO
      END DO

      B_SGS(:,:,iPart) = SQRT(2*sigmaSGS(iPart)**2)*B_SGS(:,:,iPart)
    END IF

    ! EULER
    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3) - G_SGS(1:3,j,iPart)*TurbPartState(j,iPart)*dt + B_SGS(1:3,j,iPart)*RandNormal()*SQRT(dt)
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
    randomVar=RandNormal()
  ELSE
    IF (iStage.NE.nRKStages) THEN
      RKdtFrac      = (RKC(iStage+1)-RKC(iStage))*dt
    ELSE
      RKdtFrac      = (1.-RKC(nRKStages))*dt
    END IF
  END IF
ELSE
  RKdtFrac = dt
  randomVar=RandNormal()
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
  	 kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)
!kSGSPart(iPart) = 3./2
!PRINT *, 'kSGSPart:  ', kSGSPart
  
	! Time scale of SGS scales
 	sigmaSGS(iPart) = SQRT(2./3*kSGSPart(iPart))
!PRINT *, 'sigmaSGS:  ', sigmaSGS
 
  	! Relative velocity
  	udiff(1:3) = PartState(4:6,iPart) - (FieldAtParticle(2:4,iPart) + TurbPartState(1:3,iPart))
!udiff(1:3) = 1
!PRINT *, 'udiff:  ', udiff
 
 ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
  ElemID   = PEM%Element(iPart) - offsetElem

  IF (ANY(udiff.NE.0)) THEN
    urel = udiff/SQRT(SUM(udiff**2))
PRINT *, 'urel:  ', urel
  ELSE
    urel = 0.
  END IF

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS(iPart))) THEN
    ! We ASSUME that these are the correct matrix indices
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      RETURN
    ELSE
      ! parallel
      tauL(1,iPart) = C*ElemVolN(ElemID)/(betaSGS*SQRT(SUM(udiff**2)))
      ! perpendicular
      tauL(2,iPart) = C*ElemVolN(ElemID)/(betaSGS*2*SQRT(SUM(udiff**2)))

      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            E_SGS(i,j,iPart) = EXP(-RKdtFrac/tauL(2,iPart)) + (EXP(-RKdtFrac/tauL(1,iPart))&
              - EXP(-RKdtFrac/tauL(2,iPart)))*urel(i)*urel(j)
          ELSE
            E_SGS(i,j,iPart) =                                   (EXP(-RKdtFrac/tauL(1,iPart))&
              - EXP(-RKdtFrac/tauL(2,iPart)))*urel(i)*urel(j)
          END IF
        END DO
      END DO

      W_SGS(:,:,iPart) = 0.

      ! Sum up turbulent contributions
      Pt(1:3) = 0.
      DO j = 1,3
        Pt(1:3) = Pt(1:3) + E_SGS(1:3,j,iPart)*TurbPartState(j,iPart)
      END DO
    END IF

  ! Valid SGS turbulent kinetic energy
  ELSE

    tauSGS   = C*ElemVolN(ElemID)/sigmaSGS(iPart)
!tauSGS = 2*PI/32
!PRINT *, 'tauSGS:  ', tauSGS
    IF(ALMOSTZERO(MAXVAL(udiff)))THEN
      ! parallel
      tauL(1,iPart) = tauSGS(iPart)
      ! perpendicular
      tauL(2,iPart) = tauSGS(iPart)

      ! Calculate drift and diffusion matrix
      E_SGS(:,:,iPart) = 0.
      DO i = 1,3
        E_SGS(i,i,iPart) = EXP(-RKdtFrac/tauL(2,iPart))
        W_SGS(i,i,iPart) = sigmaSGS(iPart)*SQRT(1-EXP(-2*RKdtFrac/tauL(2,iPart)))
      END DO

    ELSE
      ! parallel
      tauL(1,iPart) = tauSGS(iPart)/(SQRT(1+  betaSGS**2*SUM(udiff**2)/kSGSPart(iPart)*3/2))
      ! perpendicular
      tauL(2,iPart) = tauSGS(iPart)/(SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart(iPart)*3/2))
!PRINT *, 'tauL:  ', tauL

      ! Calculate drift and diffusion matrix
      DO i = 1,3
        DO j = 1,3
          IF (i.EQ.j) THEN
            E_SGS(i,j,iPart) = EXP(-RKdtFrac/tauL(2,iPart)) + (EXP(-RKdtFrac/tauL(1,iPart))&
              - EXP(-RKdtFrac/tauL(2,iPart)))*urel(i)*urel(j)
            W_SGS(i,j,iPart) =  sigmaSGS(iPart)*SQRT(1-EXP(-2*RKdtFrac/tauL(2,iPart)))                                                   &
                             + (sigmaSGS(iPart)*SQRT(1-EXP(-2*RKdtFrac/tauL(1,iPart)))                                                   &
                             -  sigmaSGS(iPart)*SQRT(1-EXP(-2*RKdtFrac/tauL(2,iPart))))                        *urel(i)*urel(j)
          ELSE
            E_SGS(i,j,iPart) =                                   (EXP(-RKdtFrac/tauL(1,iPart))&
              - EXP(-RKdtFrac/tauL(2,iPart)))*urel(i)*urel(j)
            W_SGS(i,j,iPart) = (sigmaSGS(iPart)*SQRT(1-EXP(-2*RKdtFrac/tauL(1,iPart)))                                                   &
                             -  sigmaSGS(iPart)*SQRT(1-EXP(-2*RKdtFrac/tauL(2,iPart))))                        *urel(i)*urel(j)
          END IF
        END DO
      END DO
    END IF

    ! Sum up turbulent contributions
    Pt(1:3) = 0.
    DO j = 1,3
      Pt(1:3) = Pt(1:3) + E_SGS(1:3,j,iPart)*TurbPartState(j,iPart) + W_SGS(1:3,j,iPart)*randomVar
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
!SDEALLOCATE(kSGS)
SDEALLOCATE(kSGSPart)
SDEALLOCATE(sigmaSGS)
SDEALLOCATE(tauSGS)
SDEALLOCATE(tauL)
SDEALLOCATE(B_SGS)
SDEALLOCATE(E_SGS)
SDEALLOCATE(G_SGS)
SDEALLOCATE(W_SGS)
SDEALLOCATE(TurbPartState)
SDEALLOCATE(TurbPt_temp)
SDEALLOCATE(ElemVolN)
!!! JIN SGS model
SDEALLOCATE(deltaT_L)
SDEALLOCATE(k_cf)
SDEALLOCATE(deltaT_E)
SDEALLOCATE(beta_Jin)
SDEALLOCATE(taup_Jin)
SDEALLOCATE(St_Jin)
SDEALLOCATE(deltaT_L_p)
!!! FUKAGATA SGS model
SDEALLOCATE(taup_Fuka)
SDEALLOCATE(Sc_Fuka)
!SDEALLOCATE(ni_B)
SDEALLOCATE(gradUxPart)
SDEALLOCATE(gradUyPart)
SDEALLOCATE(gradUzPart)
SDEALLOCATE(T_f_p)
SDEALLOCATE(U_M)
SDEALLOCATE(alpha_Fuka)
SDEALLOCATE(beta_Fuka)
SDEALLOCATE(sigmaS_Fuka)
SDEALLOCATE(ni_SGS)
SDEALLOCATE(ni_Fuka)
!!! AMIRI MODEL 
SDEALLOCATE(taup_Amiri)
SDEALLOCATE(U_square)
SDEALLOCATE(alpha_Amiri)
SDEALLOCATE(beta_Amiri)
SDEALLOCATE(sigmaS_Amiri)
SDEALLOCATE(ni_SGS_Amiri)
!!! SOMMERFELD MODEL
SDEALLOCATE(sigmaF)
SDEALLOCATE(epsilonSGS)
SDEALLOCATE(deltaR)
SDEALLOCATE(DELTA_R)
SDEALLOCATE(sigmaF_vector_square)
SDEALLOCATE(T_L)
SDEALLOCATE(T_L_f)
SDEALLOCATE(R_L)
SDEALLOCATE(L_E)
SDEALLOCATE(F_dr)
SDEALLOCATE(G_dr)
!SDEALLOCATE(R_E)
SDEALLOCATE(diag_R_E)
SDEALLOCATE(R_P)
! MINIER&PEIRANO MODEL
SDEALLOCATE(pressur_gradient)
SDEALLOCATE(second_drift_term)
SDEALLOCATE(b_L)
SDEALLOCATE(G_matrix)
SDEALLOCATE(D_matrix)
SDEALLOCATE(B_matrix)
SDEALLOCATE(gradp)
SDEALLOCATE(gradpPart)

ParticleSGSInitIsDone=.FALSE.

END SUBROUTINE ParticleFinalizeSGS

END MODULE MOD_Particle_SGS
