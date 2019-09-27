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
#include "flexi.h"

!===================================================================================================================================
! Module for DSMC Sampling and Output
!===================================================================================================================================
MODULE MOD_CalcWallParticles_Analyze
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE CalcWallSurfaceValues
  MODULE PROCEDURE CalcWallSurfaceValues
END INTERFACE

INTERFACE InitAnalyze
  MODULE PROCEDURE InitAnalyze
END INTERFACE

INTERFACE AnalyzeEquation
  MODULE PROCEDURE AnalyzeEquation
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: CalcWallSurfaceValues,InitAnalyze,AnalyzeEquation
!===================================================================================================================================

CHARACTER(LEN=255),ALLOCATABLE :: FileName_Wall(:,:)       !< output files for wall velocities per BC

CONTAINS


SUBROUTINE CalcWallSurfaceValues(during_dt_opt,restart_opt)
!===================================================================================================================================
!> Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Restart_Vars               ,ONLY: DoRestart,RestartTime
USE MOD_Analyze_Vars               ,ONLY: Analyze_dt
USE MOD_Timedisc_Vars              ,ONLY: t,dt,Tend
USE MOD_Restart_Vars               ,ONLY: RestartTime
USE MOD_DSMC_Vars                  ,ONLY: MacroSurfaceVal ,MacroSurfaceSpecVal, DSMC
USE MOD_Particle_Vars              ,ONLY: WriteMacroSurfaceValues, nSpecies, MacroValSampTime
USE MOD_Particle_Analyze_Vars      ,ONLY: TimeSample
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfMesh,nSurfSample,SampWall,CalcSurfCollis
USE MOD_Particle_Erosion_Vars
USE MOD_CalcWallParticles_Vars
USE MOD_Posti_CalcWallParticles_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL      :: during_dt_opt !routine was called during tstep (i.e. before iter=iter+1, t=t+dt...)
LOGICAL, INTENT(IN), OPTIONAL      :: restart_opt   !routine was called during tstep (i.e. before iter=iter+1, t=t+dt...)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iSpec,iSurfSide,p,q,nShift,nShiftRHS
REAL                               :: ActualTime
INTEGER, ALLOCATABLE               :: CounterTotal(:), SumCounterTotal(:)              ! Total Wall-Collision counter
LOGICAL                            :: during_dt
!===================================================================================================================================

IF (PRESENT(during_dt_opt)) THEN
  during_dt=during_dt_opt
ELSE
  during_dt=.FALSE.
END IF
IF (during_dt) THEN
  ActualTime=t+dt
ELSE
  ActualTime=t
END IF

IF (WriteMacroSurfaceValues) THEN
  TimeSample = t - MacroValSampTime !elapsed t since last sampling (variable dt's possible!)
  MacroValSampTime = t
ELSE IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
  TimeSample = t - RestartTime
ELSE
  TimeSample = (t-(1-DSMC%TimeFracSamp)*TEnd)
END IF
IF (PRESENT(restart_opt)) THEN
    IF (restart_opt) THEN
        TimeSample = Analyze_dt
        t          = MERGE(RestartTime,0.,DoRestart)
        ActualTime = t
    END IF
ELSE
  IF(ALMOSTZERO(TimeSample)) RETURN
END IF

IF(.NOT.SurfMesh%SurfOnProc) RETURN

! Allocate N+1 Species to have space for average
IF (nSpecies.EQ.1) THEN
    ALLOCATE(MacroSurfaceVal(nErosionVars-1,1:nSurfSample,1:nSurfSample,SurfMesh%nSides))
ELSE
    ALLOCATE(MacroSurfaceVal((nErosionVars-1)*(nSpecies+1),1:nSurfSample,1:nSurfSample,SurfMesh%nSides))
END IF
ALLOCATE(MacroSurfaceSpecVal(1,1:nSurfSample,1:nSurfSample,SurfMesh%nSides,nSpecies))

MacroSurfaceVal    = 0.
MacroSurfaceSpecVal= 0.

IF (CalcSurfCollis%Output) THEN
  ALLOCATE(CounterTotal(1:nSpecies))
  ALLOCATE(SumCounterTotal(1:nSpecies+1))
  CounterTotal(1:nSpecies)=0
  SumCounterTotal(1:nSpecies+1)=0
END IF

!> Erosion tracking
iSpec = 1
!---- Only one species. Only total values necessary
!===================================================================================================================================
DO iSurfSide=1,SurfMesh%nSides
  DO q=1,nSurfSample
    DO p=1,nSurfSample
      !---- 1. - .. / Impact Counter
        MacroSurfaceVal(1,p,q,iSurfSide) = SampWall(iSurfSide)%State(1,p,q)
        MacroSurfaceSpecVal(1,p,q,iSurfSide,iSpec) = SampWall(iSurfSide)%State(1,p,q) / TimeSample
      !---- 2. - .. / Impact Counter per AREA
        MacroSurfaceVal(2,p,q,iSurfSide) = SampWall(iSurfSide)%State(1,p,q) / SurfMesh%SurfaceArea(p,q,iSurfSide)
      !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
        MacroSurfaceVal(3,p,q,iSurfSide) = SampWall(iSurfSide)%State(2,p,q)
        MacroSurfaceVal(4,p,q,iSurfSide) = SampWall(iSurfSide)%State(3,p,q)
        MacroSurfaceVal(5,p,q,iSurfSide) = SampWall(iSurfSide)%State(4,p,q)
        MacroSurfaceVal(6,p,q,iSurfSide) = SampWall(iSurfSide)%State(6,p,q)
      !---- 7. - 10 / Impact angle (mean, min, max, variance)
        MacroSurfaceVal(7,p,q,iSurfSide) = SampWall(iSurfSide)%State(7,p,q)
        MacroSurfaceVal(8,p,q,iSurfSide) = SampWall(iSurfSide)%State(8,p,q)
        MacroSurfaceVal(9,p,q,iSurfSide) = SampWall(iSurfSide)%State(9,p,q)
        MacroSurfaceVal(10,p,q,iSurfSide)= SampWall(iSurfSide)%State(11,p,q)
      !---- 11 - 13 / Sampling Current Forces at walls
        MacroSurfaceVal(11,p,q,iSurfSide) = SampWall(iSurfSide)%State(12,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(12,p,q,iSurfSide) = SampWall(iSurfSide)%State(13,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(13,p,q,iSurfSide) = SampWall(iSurfSide)%State(14,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
      !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        SampWall(iSurfSide)%State(12,p,q) = 0.
!        SampWall(iSurfSide)%State(13,p,q) = 0.
!        SampWall(iSurfSide)%State(14,p,q) = 0.
      !---- 14 - 16 / Sampling Average Forces at walls
        MacroSurfaceVal(14,p,q,iSurfSide) = SampWall(iSurfSide)%State(15,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(15,p,q,iSurfSide) = SampWall(iSurfSide)%State(16,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(16,p,q,iSurfSide) = SampWall(iSurfSide)%State(17,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
    END DO
  END DO
END DO

!---- Multiple species. All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
IF (nSpecies.GT.1) THEN
  DO iSurfSide=1,SurfMesh%nSides
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        DO iSpec=1,nSpecies
          nShift    = iSpec * (nErosionVars-1)
          nShiftRHS = iSpec * nErosionVars
        !---- 1. - .. / Impact Counter
          MacroSurfaceVal(1+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(1+nShiftRHS,p,q)
          MacroSurfaceSpecVal(1,p,q,iSurfSide,iSpec)= SampWall(iSurfSide)%State(1+nShiftRHS,p,q) / TimeSample
        !---- 2. - .. / Impact Counter per AREA
          MacroSurfaceVal(2+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(1+nShiftRHS,p,q) / SurfMesh%SurfaceArea(p,q,iSurfSide)
        !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
          MacroSurfaceVal(3+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(2+nShiftRHS,p,q)
          MacroSurfaceVal(4+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(3+nShiftRHS,p,q)
          MacroSurfaceVal(5+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(4+nShiftRHS,p,q)
          MacroSurfaceVal(6+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(6+nShiftRHS,p,q)
        !---- 7. - 10 / Impact angle (mean, min, max, variance)
          MacroSurfaceVal(7+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(7+nShiftRHS,p,q)
          MacroSurfaceVal(8+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(8+nShiftRHS,p,q)
          MacroSurfaceVal(9+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(9+nShiftRHS,p,q)
          MacroSurfaceVal(10+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(11+nShiftRHS,p,q)
        !---- 11 - 13 / Sampling Current Forces at walls
          MacroSurfaceVal(11+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(12+nShiftRHS,p,q)                                    &
                                                    / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
          MacroSurfaceVal(12+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(13+nShiftRHS,p,q)                                    &
                                                    / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
          MacroSurfaceVal(13+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(14+nShiftRHS,p,q)                                    &
                                                    / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        SampWall(iSurfSide)%State(12+nShiftRHS,p,q) = 0.
!        SampWall(iSurfSide)%State(13+nShiftRHS,p,q) = 0.
!        SampWall(iSurfSide)%State(14+nShiftRHS,p,q) = 0.
        !---- 14 - 16 / Sampling Average Forces at walls
          MacroSurfaceVal(14+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(15+nShiftRHS,p,q)                                    &
                                                    / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
          MacroSurfaceVal(15+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(16+nShiftRHS,p,q)                                    &
                                                    / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
          MacroSurfaceVal(16+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(17+nShiftRHS,p,q)                                    &
                                                    / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        END DO
      END DO
    END DO
  END DO
END IF

! Only deallocate if we don't need in for wall calculations
IF (.NOT.doCalcWallParticles) THEN
    DEALLOCATE(MacroSurfaceVal,MacroSurfaceSpecVal)
ELSE
    IF(.NOT.surfAvg) THEN
        DEALLOCATE(MacroSurfaceSpecVal)
    END IF
END IF

END SUBROUTINE CalcWallSurfaceValues


!==================================================================================================================================
!> Initializes variables necessary for analyze subroutines
!> - provides basic quantities like global domain volume, surface area of boundary conditions
!>   or precomputed surface and volume integration weights
!> - initializes other specific analysis and benchmarking routines
!==================================================================================================================================
SUBROUTINE InitAnalyze()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_ReadInTools,        ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_StringTools,        ONLY: INTTOSTR
USE MOD_Interpolation_Vars, ONLY: xGP,wGP,wBary,InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: nBCs,SurfElem,nSides,AnalyzeSide,sJ,nElems
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Benchmarking,       ONLY: InitBenchmarking
USE MOD_Timedisc_Vars,      ONLY: TEnd
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,j,k,iSurf,iElem,iSide
!==================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.AnalyzeInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,'InitAnalyse not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ANALYZE...'

! Get the various analysis/output variables
Analyze_dt        =GETREAL('Analyze_dt','0.0')
nWriteData        =GETINT('nWriteData' ,'1')
NAnalyze          =GETINT('NAnalyze'   ,INTTOSTR(2*(PP_N+1)))
#if PP_dim == 3
NAnalyzeZ         =NAnalyze
#else
NAnalyzeZ         =0
#endif
! If Analyze_dt is set to 0 (default) or to a negative value, no analyze calls should be performed at all.
! To achieve this, Analyze_dt is set to the final simulation time. This will prevent any calls of the analyze routine
! except at the beginning and the end of the simulation.
IF (Analyze_dt.LE.0.) THEN
  Analyze_dt = TEnd
  nWriteData = 1
END IF

WriteData_dt = Analyze_dt*nWriteData

! precompute integration weights
ALLOCATE(wGPSurf(0:PP_N,0:PP_NZ),wGPVol(0:PP_N,0:PP_N,0:PP_NZ))
#if PP_dim == 3
DO j=0,PP_N; DO i=0,PP_N
  wGPSurf(i,j)  = wGP(i)*wGP(j)
END DO; END DO
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  wGPVol(i,j,k) = wGP(i)*wGP(j)*wGP(k)
END DO; END DO; END DO
#else
DO i=0,PP_N
  wGPSurf(i,0)  = wGP(i)
END DO
DO j=0,PP_N; DO i=0,PP_N
  wGPVol(i,j,0) = wGP(i)*wGP(j)
END DO; END DO
#endif

! precompute volume of the domain
ALLOCATE(ElemVol(nElems))
ElemVol=0.
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    ElemVol(iElem)=ElemVol(iElem)+wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO !i,j,k
END DO ! iElem
Vol=SUM(ElemVol)


! compute surface of each boundary
ALLOCATE(Surf(nBCs))
Surf=0.
DO iSide=1,nSides
  iSurf=AnalyzeSide(iSide)
  IF(iSurf.EQ.0) CYCLE
  DO j=0,PP_NZ; DO i=0,PP_N
    Surf(iSurf)=Surf(iSurf)+wGPSurf(i,j)*SurfElem(i,j,0,iSide)
  END DO; END DO
END DO
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Vol ,1   ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Surf,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/

! Initialize eval routines
CALL InitAnalyzeBasis(PP_N,NAnalyze,xGP,wBary)

CALL InitAnalyzeEquation()
CALL InitBenchmarking()

AnalyzeInitIsDone=.TRUE.
SWRITE(UNIT_StdOut,'(A,ES18.9)')' Volume of computational domain : ',Vol
SWRITE(UNIT_stdOut,'(A)')' INIT ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitAnalyze


!==================================================================================================================================
!> Initializes variables necessary for NavierStokes specific analyze subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_Equation_Vars,      ONLY: StrVarNamesPrim,StrVarNames
USE MOD_ReadInTools,        ONLY: GETLOGICAL
USE MOD_Mesh_Vars,          ONLY: nBCs,BoundaryType,BoundaryName
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_TimeAverage,        ONLY: InitCalcTimeAverage
#if USE_PARTICLES
USE MOD_Particle_Vars,      ONLY: nSpecies
USE MOD_Particle_Vars,      ONLY: WriteMacroSurfaceValues
USE MOD_CalcWallParticles
USE MOD_CalcWallParticles_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)               :: formatStr
INTEGER                         :: i,iSpec
!==================================================================================================================================
#if USE_PARTICLES
IF (doCalcWallParticles) WriteMacroSurfaceValues = .TRUE.
#endif

! Generate wallmap
ALLOCATE(isWall(nBCs))
DO i=1,nBCs
  SELECT CASE(BoundaryType(i,BC_TYPE))
  CASE(3,4,9)
    isWall(i)=.TRUE.
  CASE DEFAULT
    isWall(i)=.FALSE.
  END SELECT
END DO
maxlen=MAX(MAXVAL(LEN_TRIM(BoundaryName))+1,14)

IF(.NOT.ANY(isWall))THEN
  doCalcBodyForces=.FALSE.
  doCalcWallVelocity=.FALSE.
END IF

! Initialize eval routines
IF(MPIRoot)THEN
#if USE_PARTICLES
    IF (nSpecies.GT.1) THEN
        ALLOCATE(Filename_Wall(nBCs,0:nSpecies))
            DO i=1,nBCs
                IF(.NOT.isWall(i)) CYCLE
                FileName_Wall(i,0) = TRIM(ProjectName)//'_WallPart_'//TRIM(BoundaryName(i))
                CALL InitOutputToFile(FileName_Wall(i,0),TRIM(BoundaryName(i)),6,&
                [CHARACTER(9) :: "AlphaMean","AlphaVar","EkinMean","EkinVar","PartForce","MaxForce"])
            END DO

        DO iSpec=1,nSpecies
            DO i=1,nBCs
                IF(.NOT.isWall(i)) CYCLE
                IF(iSpec.LT.10) THEN
                    WRITE(formatStr,'(I1)') iSpec
                ELSE
                    WRITE(formatStr,'(I2)') iSpec
                END IF
                FileName_Wall(i,iSpec) = TRIM(ProjectName)//'_WallPart_Spec'//TRIM(formatStr)//'_'//TRIM(BoundaryName(i))
                CALL InitOutputToFile(FileName_Wall(i,iSpec),TRIM(BoundaryName(i)),6,&
                [CHARACTER(9) :: "AlphaMean","AlphaVar","EkinMean","EkinVar","PartForce","MaxForce"])
            END DO
        END DO
    ELSE
        ALLOCATE(Filename_Wall(nBCs,1))
        DO i=1,nBCs
                IF(.NOT.isWall(i)) CYCLE
                FileName_Wall(i,1) = TRIM(ProjectName)//'_WallPart_'//TRIM(BoundaryName(i))
                CALL InitOutputToFile(FileName_Wall(i,1),TRIM(BoundaryName(i)),6,&
                [CHARACTER(9) :: "AlphaMean","AlphaVar","EkinMean","EkinVar","PartForce","MaxForce"])
        END DO
    END IF
#endif
END IF

IF(doCalcTimeAverage)  CALL InitCalcTimeAverage()

END SUBROUTINE InitAnalyzeEquation

!==================================================================================================================================
!> Wrapper routine for the equation system specific analyze routines. Will call the specific subroutines to calculate the quantities
!> set in the parameter file and the respective output routines.
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation_Vars
USE MOD_DSMC_Vars,          ONLY: MacroSurfaceVal
USE MOD_Mesh_Vars,          ONLY: nBCs
USE MOD_Output,             ONLY: OutputToFile
#if USE_PARTICLES
USE MOD_Particle_Vars,      ONLY: nSpecies
USE MOD_Posti_CalcWallParticles
USE MOD_Posti_CalcWallParticles_Vars
USE MOD_CalcWallParticles_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                              !< Current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(LEN=40)               :: formatStr
INTEGER                         :: i,iSpec
#if USE_PARTICLES
REAL,DIMENSION(nBCs)            :: AlphaMean,AlphaVar,EkinMean,EkinVar,PartForce,MaxForce
#endif
!==================================================================================================================================
! Calculate derived quantities
#if USE_PARTICLES
IF(MPIRoot.AND.doCalcWallParticles)THEN
    IF (nSpecies.GT.1) THEN
        CALL CalcWallParticles(AlphaMean,AlphaVar,EkinMean,EkinVar,partForce,maxForce,DoSpecies=.FALSE.)
        DO i=1,nBCs
            IF(.NOT.isWall(i)) CYCLE
            CALL OutputToFile(FileName_Wall(i,0),(/Time/),(/6,1/),(/AlphaMean(i),AlphaVar(i),EkinMean(i),EkinVar(i),&
                                                               PartForce(i),MaxForce(i)/))
        END DO

        DO iSpec=1,nSpecies
            CALL CalcWallParticles(AlphaMean,AlphaVar,EkinMean,EkinVar,partForce,maxForce,DoSpecies=.TRUE.,Species_opt=iSpec)
            DO i=1,nBCs
                IF(.NOT.isWall(i)) CYCLE
                CALL OutputToFile(FileName_Wall(i,iSpec),(/Time/),(/6,1/),(/AlphaMean(i),AlphaVar(i),EkinMean(i),EkinVar(i),&
                                                                PartForce(i),MaxForce(i)/))
            END DO
        END DO
    ELSE
        CALL CalcWallParticles(AlphaMean,AlphaVar,EkinMean,EkinVar,partForce,maxForce,DoSpecies=.FALSE.)
        DO i=1,nBCs
            IF(.NOT.isWall(i)) CYCLE
            CALL OutputToFile(FileName_Wall(i,1),(/Time/),(/6,1/),(/AlphaMean(i),AlphaVar(i),EkinMean(i),EkinVar(i),&
                                                               PartForce(i),MaxForce(i)/))
        END DO
    END IF
END IF
#endif

! Deallocate Macrosurfaces we couldn't deallocate earlier
IF(.NOT.surfAvg) THEN
    SDEALLOCATE(MacroSurfaceVal)
END IF

END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!> - Builds Vandermonde to interpolate the solution onto a Gauss-Lobatto mesh at a higher polynomial degree
!> - Precomputes volume interpolation weights
!==================================================================================================================================
SUBROUTINE InitAnalyzeBasis(N_in,Nloc,xGP,wBary)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars, ONLY: wGPVolAnalyze,Vdm_GaussN_NAnalyze
USE MOD_Basis,        ONLY: InitializeVandermonde
USE MOD_Interpolation,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars,ONLY: NodeTypeGL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: N_in                  !< input polynomial degree
INTEGER,INTENT(IN)               :: Nloc                  !< polynomial degree of analysis polynomial
REAL,INTENT(IN)                  :: xGP(0:N_in)           !< interpolation points
REAL,INTENT(IN)                  :: wBary(0:N_in)         !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: XiAnalyze(0:Nloc)
REAL                             :: wAnalyze( 0:Nloc)  ! GL integration weights used for the analyze
INTEGER                          :: i,j
#if PP_dim == 3
INTEGER                          :: k
#endif
!==================================================================================================================================
ALLOCATE(wGPVolAnalyze(0:Nloc,0:Nloc,0:ZDIM(Nloc)),Vdm_GaussN_NAnalyze(0:Nloc,0:N_in))
CALL GetNodesAndWeights(Nloc,NodeTypeGL,XiAnalyze,wAnalyze)
CALL InitializeVandermonde(N_in,Nloc,wBary,xGP,XiAnalyze,Vdm_GaussN_NAnalyze)

#if PP_dim == 3
DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
  wGPVolAnalyze(i,j,k) = wAnalyze(i)*wAnalyze(j)*wAnalyze(k)
END DO; END DO; END DO
#else
DO j=0,Nloc; DO i=0,Nloc
  wGPVolAnalyze(i,j,0) = wAnalyze(i)*wAnalyze(j)
END DO; END DO
#endif

END SUBROUTINE InitAnalyzeBasis


END MODULE MOD_CalcWallParticles_Analyze
