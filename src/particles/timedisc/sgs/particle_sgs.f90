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
             taup_Fuka    (           1:PDM%maxParticleNumber),      &
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
USE MOD_Mesh_Vars                   ,ONLY: offsetElem
USE MOD_Particle_Globals
USE MOD_Particle_SGS_Vars
USE MOD_Particle_Interpolation      ,ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars ,ONLY: FieldAtParticle
USE MOD_Particle_Vars               ,ONLY: TurbPartState,PDM,PEM,PartState
!!! JIN and FUKAGATA MODEL !!!
USE MOD_EOS_Vars,          ONLY : mu0
USE MOD_Particle_Vars,     ONLY : Species, PartSpecies
USE MOD_Lifting_Vars, 	   ONLY : gradUx, gradUy, gradUz !! in order to calculate rot_U (only Fukagata model)
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
INTEGER                       :: ElemID, i, j, iPart
REAL                          :: Pt(1:3)
REAL                          :: RKdtFrac
REAL 			      :: filterwidth ! IN ORDER TO VALIDATE FLEXI IMPLEMENTATION
! LOCAL VARIABLES - FUKAGATA SGS MODEL
REAL                          :: nu  			! kynematic viscosity
REAL,PARAMETER                :: C_c= 1. 		! factor of Stokes-Cunningham
REAL,PARAMETER                :: C_Fuka= 0.93 		! constant model of Fukagata
REAL,PARAMETER                :: C_eps=0.08 		! = (2/(3*Ko))**(3/2), Ko: Kolmogorov constant = 1.5
REAL                          :: beta1
REAL                          :: beta2
REAL    		      :: rot_U(3)
REAL                          :: d_ij_tensor(3,3) 	! rate of deformation tensor
REAL                          :: d_ij  			! used for the calculation of representative u_sgs
!===================================================================================================================================

SELECT CASE(SGSModel)

CASE('none')
!===================================================================================================================================
! DNS mode. Assume fully resolved turbulence
!===================================================================================================================================
! Do nothing


CASE('Fukagata-Analytic')
!===================================================================================================================================
! K. Fukagata (2004)
!> Dynamics of Brownian particles in a turbulent channel flow, K. Fukagata, S. Zahrai, F. H. Bark
!> Heat and Mass Transfer 40 (2004) 715–726. DOI 10.1007/s00231-003-0462-8
!===================================================================================================================================

! nu = mu0/FieldAtParticle(1,iPart)  !! for the Schmidt number 

! Time integration using Euler-Maruyama scheme
IF (iStage.NE.1) RETURN

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
PRINT *, 'filterwidth:  ', filterwidth

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  
	! Particle relaxation time
	taup_Fuka(iPart)    = C_c * Species(PartSpecies(iPart))%DensityIC * Species(PartSpecies(iPart))%DiameterIC**2 &
	& * 1./(18*mu0)			
PRINT *, 'taup_Fuka:  ', taup_Fuka

	! Integral time scale of fluid velocity along a particle trajectory
	T_f_p(iPart) =1
PRINT *, 'T_f_p:  ', T_f_p

	! Calculation of deformation tensor d_ij
	d_ij = 1
PRINT *, 'd_ij:  ', d_ij

	! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity due to SGS velocity in time dt
 	U_M(iPart) = 2/(3**0.5) * C_Fuka * C_eps * filterwidth * d_ij	
PRINT *, 'U_M:  ', U_M
	
!alpha and beta Fukagata
	alpha_Fuka(ipart) =  dt/taup_Fuka(iPart)
	beta_Fuka(iPart)  = taup_Fuka(iPart)/T_f_p(iPart) 
	beta1 = 1+beta_Fuka(iPart)
	beta2 = 1-beta_Fuka(iPart)
PRINT *, 'alpha_Fuka:  ', alpha_Fuka
PRINT *, 'beta_Fuka:  ', beta_Fuka

	! sigmaS_Fuka is the increase of standard deviation of particle velocity due to SGS velocity during time dt
	sigmaS_Fuka(iPart) = U_M(iPart) * SQRT(	1/beta1 * (1-EXP(-alpha_Fuka(ipart)*beta1)) 	&
				        &      -1/beta2 * EXP(-2*alpha_Fuka(ipart)) * (1-EXP(alpha_Fuka(ipart)*beta2))		)
PRINT *, 'sigmaS_Fuka:  ', sigmaS_Fuka

	! SGS fluid velocity by Fukagata
	DO j = 1,3
	ni_SGS(1:3,iPart) = sigmaS_Fuka(iPart)/dt !* RandNormal()
	END DO 
PRINT *, 'ni_SGS:  ', ni_SGS(1,iPart)

    	! sum up the two terms 
	ni_Fuka(1:3,iPart) = ni_SGS(1:3,iPart)


	! TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)
	TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + ni_Fuka(1:3,iPart)

END DO


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

filterwidth = 2*PI/32  ! validation of the implementation
!PRINT *, 'filterwidth:  ', filterwidth

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Use unfiltered density to obtain velocity in primitive variables
  ! kSGSPart(iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)
	kSGSPart(iPart) = 1;

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
!PRINT *, 'taup_Fuka:  ', taup_Fuka

	! Schmidt number of the particles in fluid
!	Sc_Fuka(iPart) = Species(PartSpecies(iPart))%MassIC*nu / ( taup_Fuka(iPart) * kSGSPart(iPart) * FieldAtParticle(TEMP,iPart) ) 

	! Brownian force 
!	ni_B(iPart) = RandNormal()*SQRT(2*nu/(Sc_Fuka(iPart) * taup_Fuka(iPart)**2 * dt))

	! Integral time scale of fluid velocity along a particle trajectory by estimated by an inverse of resolved local vorticity
		! Local vorticity 
		rot_U(1)=gradUyPart(3,iPart) - gradUzPart(2,iPart)    ! dw/dy-dv/dz
		rot_U(2)=gradUzPart(1,iPart) - gradUxPart(3,iPart)    ! du/dz-dw/dx
		rot_U(3)=gradUxPart(2,iPart) - gradUyPart(1,iPart)    ! dv/dx-du/dy 		

	!T_f_p(iPart) =1./VECNORM(rot_U)
	T_f_p(iPart) = 1.
!PRINT *, 'T_f_p:  ', T_f_p

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
	!d_ij = SQRT(d_ij)
	d_ij = 1.
!PRINT *, 'd_ij:  ', d_ij
	
	! Calculate the representative SGS velocity in order to model the increase of variance of particle velocity due to SGS velocity in time dt
 	!U_M(iPart) = 2/(3**0.5) * C_Fuka * C_eps * ElemVolN(ElemID) * d_ij	
	U_M(iPart) = 2/(3**0.5) * C_Fuka * C_eps * filterwidth * d_ij
!PRINT *, 'U_M:  ', U_M

	!alpha and beta Fukagata
	alpha_Fuka(ipart) =  dt/taup_Fuka(iPart)
	beta_Fuka(iPart)  = taup_Fuka(iPart)/T_f_p(iPart) 
	beta1 = 1+beta_Fuka(iPart)
	beta2 = 1-beta_Fuka(iPart)
!PRINT *, 'alpha_Fuka:  ', alpha_Fuka
!PRINT *, 'beta_Fuka:  ', beta_Fuka

	! sigmaS_Fuka is the increase of standard deviation of particle velocity due to SGS velocity during time dt
	
	sigmaS_Fuka(iPart) = U_M(iPart) * SQRT(	1/beta1 * (1-EXP(-alpha_Fuka(ipart)*beta1)) 	&
				&	-1/beta2 * EXP(-2*alpha_Fuka(ipart)) * (1-EXP(alpha_Fuka(ipart)*beta2))		)
!PRINT *, 'sigmaS_Fuka:  ', sigmaS_Fuka

	! SGS fluid velocity by Fukagata
	DO j = 1,3
	ni_SGS(1:3,iPart) =  sigmaS_Fuka(iPart)/dt !* RandNormal()
	END DO  

    	! sum up the two terms
	!!!!!! now I'm using only te SGS fluid velocity and NO the brownian force
	ni_Fuka(1:3,iPart) = ni_SGS(1:3,iPart) !+ ni_B(iPart)

  END IF
	! TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)
	TurbPartState(1:3,iPart) =  ni_Fuka(1:3,iPart)

END DO


CASE('Jin-Analytic')
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
PRINT *, 'filterwidth:  ', filterwidth

DO iPart = 1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  
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
      Pt(1:3) = Pt(1:3) - 1./deltaT_L_p(iPart)*TurbPartState(j,iPart)*dt + (4/3*kSGSPart(iPart)*dt/deltaT_L_p(iPart))**0.5!*RandNormal()
    END DO

  END IF

  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)

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
      Pt(1:3) = Pt(1:3) - 1./deltaT_L_p(iPart)*TurbPartState(j,iPart)*dt + ((4/3*kSGSPart(iPart)*dt/deltaT_L_p(iPart))**0.5)*RandNormal()
    END DO

  END IF

  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)

END DO


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

  END DO

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
        W_SGS(i,j,iPart) = sigmaSGS(iPart)*SQRT(1-EXP(-2*RKdtFrac/tauL(2,iPart)))
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
! JIN SGS model
SDEALLOCATE(deltaT_L)
SDEALLOCATE(k_cf)
SDEALLOCATE(deltaT_E)
SDEALLOCATE(beta_Jin)
SDEALLOCATE(taup_Jin)
SDEALLOCATE(St_Jin)
SDEALLOCATE(deltaT_L_p)
! FUKAGATA SGS model
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

ParticleSGSInitIsDone=.FALSE.

END SUBROUTINE ParticleFinalizeSGS

END MODULE MOD_Particle_SGS
