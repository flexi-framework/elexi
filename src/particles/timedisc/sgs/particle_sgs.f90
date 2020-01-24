!=================================================================================================================================
! Copyright (c) 2010-2020  Prof. Claus-Dieter Munz
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
!
!===================================================================================================================================
SUBROUTINE ParticleInitSGS()
! MODULES
USE MOD_Globals,                    ONLY: ABORT, MPIRoot, Unit_STDOUT
USE MOD_Preproc,                    ONLY: PP_N, PP_NZ
USE MOD_Mesh_Vars,                  ONLY: nElems
USE MOD_ReadInTools,                ONLY: GETINT,GETSTR
USE MOD_Particle_SGS_Vars
USE MOD_Particle_Vars,              ONLY: PDM,TurbPartState,TurbPt_temp
#if USE_RW
USE MOD_Particle_Randomwalk_Vars,   ONLY: RWModel
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ALLOCSTAT
!----------------------------------------------------------------------------------------------------------------------------------
IF(ParticleSGSInitIsDone) RETURN

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SGS MODEL ...'

! SGS model
SGSModel = TRIM(GETSTR('Part-SGSModel','none'))

! SGS and RW are not compatible. Go back and think about what you were trying to simulate.
#if USE_RW
IF (SGSModel.NE.'none'.AND.RWModel.NE.'none') &
  CALL abort(__STAMP__,'SGS and RW not compatible!')
#endif

SELECT CASE(SGSModel)
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
             kSGS         (1         ,0:PP_N,0:PP_N,0:PP_NZ,nElems), &
             kSGSPart     (1         ,1:PDM%maxParticleNumber),      &
             sigmaSGS     (1         ,1:PDM%maxParticleNumber),      &
             tauSGS       (1         ,1:PDM%maxParticleNumber),      &
             tauL         (1:2       ,1:PDM%maxParticleNumber),      &
             G_SGS        (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             B_SGS        (1:3,1:3   ,1:PDM%maxParticleNumber),      &
             TurbPartState(1:nSGSVars,1:PDM%maxParticleNumber),      &
             TurbPt_temp  (1:3       ,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) &
      CALL abort(__STAMP__,'ERROR in particle_sgs.f90: Cannot allocate particle SGS arrays!')
    TurbPartState = 0.
    TurbPt_temp   = 0.

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
IF(FilterType.GT.0) CALL CollectiveStop(__STAMP__,"SGS incompatible with Navier-Stokes filter!")

! Prepare modal cut-off filter (low pass)
ALLOCATE(FilterMat(0:PP_N,0:PP_N))
FilterMat = 0.

! Modal Filter, default to cut-off at PP_N-2
NFilter = PP_N - nSGSFilter
DO iDeg=0,NFilter
  FilterMat(iDeg,iDeg) = 1.
END DO

! Assemble filter matrix in nodal space
FilterMat=MATMUL(MATMUL(Vdm_Leg,FilterMat),sVdm_Leg)

FilterInitIsDone = .TRUE.
END SUBROUTINE InitSGSFilter

!===================================================================================================================================
! SGS deconvolution
!===================================================================================================================================
SUBROUTINE ParticleSGS(iStage,dt,b_dt)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Analyze_Vars          ,ONLY: ElemVol
USE MOD_DG_Vars               ,ONLY: U
USE MOD_Filter                ,ONLY: Filter_Pointer
USE MOD_Filter_Vars           ,ONLY: FilterMat
USE MOD_Particle_Globals
USE MOD_Particle_SGS_Vars
USE MOD_Particle_Interpolation, ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars, ONLY: FieldAtParticle
USE MOD_Particle_Vars         ,ONLY: TurbPartState, PDM, PEM, PartState, TurbPt_Temp
USE MOD_TimeDisc_Vars         ,ONLY: RKA
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iStage
REAL,INTENT(IN)               :: dt
REAL,INTENT(IN)               :: b_dt
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
!===================================================================================================================================

SELECT CASE(SGSModel)

CASE('none')
!===================================================================================================================================
! DNS mode. Assume fully resolved turbulence
!===================================================================================================================================
  ! Do nothing

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
CALL InterpolateFieldToParticle(4,USGS(1:4,:,:,:,:),USGSPart)

DO iPart=1,PDM%ParticleVecLength
  ! Only consider particles
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  ! Use unfiltered density to obtain velocity in primitive variables
  kSGSPart(1,iPart) = 0.5*SUM((USGSPart(2:4,iPart)/FieldAtParticle(1,iPart))**2.)

  ! Time scale of SGS scales
  sigmaSGS(1,iPart) = SQRT(2./3.*kSGSPart(1,iPart))

  ! No SGS turbulent kinetic energy, avoid float error
  IF (ALMOSTZERO(sigmaSGS(1,iPart))) THEN
    ! We ASSUME that these are the correct matrix indices
    G_SGS(:,:,iPart) = 0.
    DO i=1,3
      G_SGS(i,i,iPart) = 1.
    END DO
    B_SGS(:,:,iPart) = 0.

  ! Valid SGS turbulent kinetic energy
  ELSE
    ! Estimate the filter width with the equivalent cell length and polynominal degree, see Flad (2017)
    ElemID   = PEM%Element(iPart)
    tauSGS   = C*(ElemVol(ElemID)**(1./3.)/(PP_N+1))/sigmaSGS(1,iPart)

    ! Relative velocity
    udiff(1:3) = PartState(4:6,iPart) - (FieldAtParticle(2:4,iPart)/FieldAtParticle(1,iPart) + TurbPartState(1:3,iPart))
    urel       = udiff/SQRT(SUM(udiff**2))
    ! parallel
    tauL(1,iPart) = tauSGS(1,iPart)/(SQRT(1+betaSGS**2*SUM(udiff**2)/kSGSPart(1,iPart)*3/2))
    ! perpendicular
    tauL(2,iPart) = tauSGS(1,iPart)/(SQRT(1+4*betaSGS**2*SUM(udiff**2)/kSGSPart(1,iPart)*3/2))

    ! Calculate drift and diffusion matrix
    DO i=1,3
      DO j=1,3
        IF (i.EQ.j) THEN
          G_SGS(i,j,iPart) = 1/tauL(2,iPart) + (1/tauL(1,iPart) - 1/tauL(2,iPart))*urel(i)*urel(j)
          B_SGS(i,j,iPart) = 1/SQRT(tauL(2,iPart)) + (1/SQRT(tauL(1,iPart)) - 1/SQRT(tauL(2,iPart)))*urel(i)*urel(j)
        ELSE
          G_SGS(i,j,iPart) =                   (1/tauL(1,iPart) - 1/tauL(2,iPart))*urel(i)*urel(j)
          B_SGS(i,j,iPart) =                   (1/SQRT(tauL(1,iPart)) - 1/SQRT(tauL(2,iPart)))*urel(i)*urel(j)
        END IF
      END DO
    END DO

    B_SGS(:,:,iPart)=SQRT(2*sigmaSGS(1,iPart)**2)*B_SGS(:,:,iPart)
  END IF

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

  ! EULER
  ! Sum up turbulent contributions
  Pt(1:3)=0.
  DO j=1,3
    Pt(1:3) = Pt(1:3) - G_SGS(1:3,j,iPart)*TurbPartState(j,iPart)*dt + B_SGS(1:3,j,iPart)*RandNormal()*SQRT(dt)
  END DO
  TurbPartState(1:3,iPart) = TurbPartState(1:3,iPart) + Pt(1:3)

  ! Sanity check. Use 10*U as arbitrary threshold
  IF (ANY(ABS(TurbPartState(1:3,iPart)).GT.10.*MAXVAL(ABS(PartState(4:6,iPart))))) THEN
    IPWRITE(*,*) 'Obtained SGS velocity of',SQRT(SUM(TurbPartState(1:3,iPart)**2)), &
                 'while resolved velocity is',SQRT(SUM(PartState(4:6,iPart)**2)),   &
                 'for',iPart
!    TurbPartState(:,iPart) = 0.
  END IF
END DO

CASE DEFAULT
    CALL abort(__STAMP__, ' No particle SGS model given. This should not happen.')

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
SDEALLOCATE(kSGS)
SDEALLOCATE(kSGSPart)
SDEALLOCATE(sigmaSGS)
SDEALLOCATE(tauSGS)
SDEALLOCATE(tauL)
SDEALLOCATE(G_SGS)
SDEALLOCATE(B_SGS)
SDEALLOCATE(TurbPartState)
SDEALLOCATE(TurbPt_temp)

ParticleSGSInitIsDone=.FALSE.

END SUBROUTINE ParticleFinalizeSGS

END MODULE MOD_Particle_SGS
