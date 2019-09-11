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
MODULE MOD_Particle_Erosion_Analyze
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE CalcSurfaceValues
  MODULE PROCEDURE CalcSurfaceValues
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: CalcSurfaceValues
!===================================================================================================================================

CONTAINS


SUBROUTINE CalcSurfaceValues(during_dt_opt,restart_opt,remap_opt)
!===================================================================================================================================
!> Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Restart_Vars               ,ONLY: DoRestart,RestartTime
USE MOD_Analyze_Vars               ,ONLY: Analyze_dt
USE MOD_Mesh_Vars                  ,ONLY: MeshFile
USE MOD_Timedisc_Vars              ,ONLY: t,dt,Tend
USE MOD_Restart_Vars               ,ONLY: RestartTime
USE MOD_DSMC_Vars                  ,ONLY: MacroSurfaceVal ,MacroSurfaceSpecVal, DSMC
USE MOD_Particle_Vars              ,ONLY: WriteMacroSurfaceValues, nSpecies, MacroValSampTime
USE MOD_Particle_Analyze_Vars      ,ONLY: TimeSample
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfMesh,nSurfSample,SampWall,CalcSurfCollis
USE MOD_Particle_Boundary_Sampling ,ONLY: WriteSurfSampleToHDF5
USE MOD_Particle_Erosion_Vars
USE MOD_CalcWallParticles_Vars
#if USE_MPI
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfCOMM
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
LOGICAL, INTENT(IN), OPTIONAL      :: during_dt_opt !routine was called during tstep (i.e. before iter=iter+1, t=t+dt...)
LOGICAL, INTENT(IN), OPTIONAL      :: restart_opt   !routine was called during tstep (i.e. before iter=iter+1, t=t+dt...)
CHARACTER(LEN=*),INTENT(IN),OPTIONAL::remap_opt     !routine was called from posti. Change output file name
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
    TimeSample = t - MacroValSampTime               ! Elapsed time since last sampling (variable dt's possible!)
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
ELSE ! Avoid division by zero if no simulation time has passed
    IF(ALMOSTZERO(TimeSample)) RETURN
END IF

!IF (CalcSurfCollis%AnalyzeSurfCollis.AND.(t.NE.0.)) THEN
!    CALL WriteAnalyzeSurfCollisToHDF5(ActualTime,TimeSample)
!END IF

! Do not try to record impacts if there are no walls on the current proc
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
        !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        MacroSurfaceVal(11+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(12+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(12+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(13+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(13+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(14+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        SampWall(iSurfSide)%State(12+nShiftRHS,p,q) = 0.
!        SampWall(iSurfSide)%State(13+nShiftRHS,p,q) = 0.
!        SampWall(iSurfSide)%State(14+nShiftRHS,p,q) = 0.
        !---- 14 - 16 / Sampling Average Forces at walls
        MacroSurfaceVal(14+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(15+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(15+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(16+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(16+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(17+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        END DO
      END DO
    END DO
  END DO
END IF

!===================================================================================================================================
! ANALYZE SURF COLLIS
! LEGACY CODE: IGNORE FOR NOW
!===================================================================================================================================
IF (CalcSurfCollis%Output) THEN
#if USE_MPI
    CALL MPI_REDUCE(CounterTotal,SumCounterTotal(1:nSpecies),nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%COMM,iError)
#else
    SumCounterTotal(1:nSpecies)=CounterTotal
#endif
    DO iSpec=1,nSpecies
        IF (CalcSurfCollis%SpeciesFlags(iSpec)) THEN ! Sum up all Collisions with SpeciesFlags for output
            SumCounterTotal(nSpecies+1) = SumCounterTotal(nSpecies+1) + SumCounterTotal(iSpec)
        END IF
    END DO
    SWRITE(UNIT_stdOut,'(A)') ' The following species swaps at walls have been sampled:'
    DO iSpec=1,nSpecies
        SWRITE(*,'(A9,I2,A2,E16.9,A6)') ' Species ',iSpec,': ',REAL(SumCounterTotal(iSpec)) / TimeSample,' MP/s;'
    END DO
    SWRITE(*,'(A23,E16.9,A6)') ' All with SpeciesFlag: ', REAL(SumCounterTotal(nSpecies+1)) / TimeSample,' MP/s.'
    DEALLOCATE(CounterTotal)
    DEALLOCATE(SumCounterTotal)
END IF

CALL WriteSurfSampleToHDF5(TRIM(MeshFile),ActualTime,remap_opt)

! Only deallocate if we don't need the values for wall calculations
IF (.NOT.doCalcWallParticles) THEN
    DEALLOCATE(MacroSurfaceVal,MacroSurfaceSpecVal)
ELSE
    DEALLOCATE(MacroSurfaceSpecVal)
END IF

END SUBROUTINE CalcSurfaceValues


!SUBROUTINE WriteAnalyzeSurfCollisToHDF5(OutputTime,TimeSample)
!!===================================================================================================================================
!!> Wrinting AnalyzeSurfCollis-Data to hdf5 file (based on WriteParticleToHDF5 and WriteDSMCHOToHDF5)
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Particle_Vars          ,ONLY: nSpecies
!USE MOD_Output_Vars            ,ONLY: ProjectName
!USE MOD_io_HDF5
!USE MOD_PICDepo_Vars           ,ONLY: SFResampleAnalyzeSurfCollis, LastAnalyzeSurfCollis, r_SF
!USE MOD_Particle_Boundary_Vars ,ONLY: nPartBound, AnalyzeSurfCollis
!USE MOD_Particle_HDF5_Output   ,ONLY: WriteAttributeToHDF5, WriteHDF5Header, WriteArrayToHDF5
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)                :: OutputTime, TimeSample
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!CHARACTER(LEN=255)             :: Filename, TypeString, H5_Name
!INTEGER,ALLOCATABLE            :: SpeciesPositions(:,:)
!CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
!#if USE_MPI
!INTEGER,ALLOCATABLE            :: sendbuf(:),recvbuf(:)
!REAL,ALLOCATABLE               :: sendbuf2(:),recvbuf2(:)
!INTEGER                        :: iProc
!INTEGER                        :: globalNum(0:nProcessors-1), Displace(0:nProcessors-1), RecCount(0:nProcessors-1)
!#endif
!INTEGER                        :: TotalNumberMPF, counter2, BCTotalNumberMPF
!INTEGER,ALLOCATABLE            :: locnPart(:),offsetnPart(:),nPart_glob(:),minnParts(:), iPartCount(:)
!INTEGER                        :: iPart, iSpec, counter
!REAL,ALLOCATABLE               :: PartData(:,:)
!INTEGER                        :: PartDataSize       !number of entries in each line of PartData
!REAL                           :: TotalFlowrateMPF, RandVal, BCTotalFlowrateMPF
!LOGICAL,ALLOCATABLE            :: PartDone(:)
!!===================================================================================================================================
!SWRITE(*,*) 'WRITE EROSION SURFACE COLLISIONS TO FILE...'
!
!TypeString='DSMCSurfCollis'
!FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
!PartDataSize=10
!ALLOCATE(StrVarNames(PartDataSize))
!StrVarNames(1)='ParticlePositionX'
!StrVarNames(2)='ParticlePositionY'
!StrVarNames(3)='ParticlePositionZ'
!StrVarNames(4)='VelocityX'
!StrVarNames(5)='VelocityY'
!StrVarNames(6)='VelocityZ'
!StrVarNames(7)='OldParticlePositionX'
!StrVarNames(8)='OldParticlePositionY'
!StrVarNames(9)='OldParticlePositionZ'
!StrVarNames(10)='BCid'
!ALLOCATE(locnPart(1:nSpecies) &
!        ,offsetnPart(1:nSpecies) &
!        ,nPart_glob(1:nSpecies) &
!        ,minnParts(1:nSpecies) &
!        ,iPartCount(1:nSpecies) )
!#if USE_MPI
!ALLOCATE(sendbuf(1:nSpecies) &
!        ,recvbuf(1:nSpecies) )
!#endif
!ALLOCATE(SpeciesPositions( 1:nSpecies,1:MAXVAL(AnalyzeSurfCollis%Number(1:nSpecies)) ))
!
!iPartCount(:)=0
!DO iPart=1,AnalyzeSurfCollis%Number(nSpecies+1)
!  IF (AnalyzeSurfCollis%Spec(iPart).LT.1 .OR. AnalyzeSurfCollis%Spec(iPart).GT.nSpecies) THEN
!    CALL Abort(&
!      __STAMP__,&
!      'Error 1 in AnalyzeSurfCollis!')
!  ELSE
!    iPartCount(AnalyzeSurfCollis%Spec(iPart))=iPartCount(AnalyzeSurfCollis%Spec(iPart))+1
!    SpeciesPositions(AnalyzeSurfCollis%Spec(iPart),iPartCount(AnalyzeSurfCollis%Spec(iPart)))=iPart
!  END IF
!END DO
!DO iSpec=1,nSpecies
!  locnPart(iSpec) = AnalyzeSurfCollis%Number(iSpec)
!  IF (iPartCount(iSpec).NE.locnPart(iSpec)) CALL Abort(&
!    __STAMP__,&
!    'Error 2 in AnalyzeSurfCollis!')
!END DO
!
!#if USE_MPI
!sendbuf(:)=locnPart(:)
!recvbuf(:)=0
!CALL MPI_EXSCAN(sendbuf,recvbuf,nSpecies,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
!offsetnPart(:)=recvbuf(:)
!sendbuf(:)=recvbuf(:)+locnPart(:)
!CALL MPI_BCAST(sendbuf(:),nSpecies,MPI_INTEGER,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!!global numbers
!nPart_glob(:)=sendbuf(:)
!DEALLOCATE(sendbuf &
!          ,recvbuf )
!!LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
!CALL MPI_ALLREDUCE(locnPart(:),minnParts(:),nSpecies,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,IERROR)
!IF (SFResampleAnalyzeSurfCollis) THEN
!  CALL MPI_ALLGATHER(AnalyzeSurfCollis%Number(nSpecies+1), 1, MPI_INTEGER, globalNum, 1, MPI_INTEGER, MPI_COMM_WORLD, IERROR)
!  TotalNumberMPF = SUM(globalNum)
!ELSE
!  CALL MPI_ALLREDUCE(AnalyzeSurfCollis%Number(nSpecies+1),TotalNumberMPF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
!END IF
!#else
!offsetnPart(:)=0
!nPart_glob(:)=locnPart(:)
!minnParts(:)=locnPart(:)
!TotalNumberMPF=AnalyzeSurfCollis%Number(nSpecies+1)
!#endif
!! determine number of parts at BC of interest
!BCTotalNumberMPF=0
!IF (SFResampleAnalyzeSurfCollis) THEN
!  DO iPart=1,AnalyzeSurfCollis%Number(nSpecies+1)
!    IF (AnalyzeSurfCollis%BCid(iPart).LT.1 .OR. AnalyzeSurfCollis%BCid(iPart).GT.nPartBound) THEN
!      CALL Abort(&
!        __STAMP__,&
!        'Error 3 in AnalyzeSurfCollis!')
!    ELSE IF ( ANY(LastAnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(LastAnalyzeSurfCollis%BCs.EQ.AnalyzeSurfCollis%BCid(iPart)) ) THEN
!      BCTotalNumberMPF = BCTotalNumberMPF + 1
!    END IF
!  END DO
!#if USE_MPI
!  CALL MPI_ALLREDUCE(MPI_IN_PLACE,BCTotalNumberMPF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
!#endif
!  BCTotalFlowrateMPF=REAL(BCTotalNumberMPF)/TimeSample
!END IF
!TotalFlowrateMPF=REAL(TotalNumberMPF)/TimeSample
!
!IF(MPIRoot) THEN !create File-Skeleton
!  ! Create file
!  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
!
!  ! Write file header
!  CALL WriteHDF5Header(TRIM(TypeString),File_ID)
!
!  ! Write dataset properties "Time","VarNames","nSpecies","TotalFlowrateMPF"
!  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
!  CALL WriteAttributeToHDF5(File_ID,'VarNames',PartDataSize,StrArray=StrVarNames)
!  CALL WriteAttributeToHDF5(File_ID,'NSpecies',1,IntegerScalar=nSpecies)
!  CALL WriteAttributeToHDF5(File_ID,'TotalFlowrateMPF',1,RealScalar=TotalFlowrateMPF)
!
!  CALL CloseDataFile()
!END IF
!
!#if USE_MPI
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
!#else
!CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
!#endif
!
!
!IF (SFResampleAnalyzeSurfCollis) THEN
!  IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts to MaxPartNumber
!    LastAnalyzeSurfCollis%PartNumberSamp=MIN(BCTotalNumberMPF,LastAnalyzeSurfCollis%PartNumberReduced)
!    ALLOCATE(PartDone(1:TotalNumberMPF))
!    PartDone(:)=.FALSE.
!  ELSE
!    LastAnalyzeSurfCollis%PartNumberSamp=BCTotalNumberMPF
!  END IF
!  SWRITE(*,*) 'Number of saved particles for SFResampleAnalyzeSurfCollis: ',LastAnalyzeSurfCollis%PartNumberSamp
!  SDEALLOCATE(LastAnalyzeSurfCollis%WallState)
!  SDEALLOCATE(LastAnalyzeSurfCollis%Species)
!  ALLOCATE(LastAnalyzeSurfCollis%WallState(6,LastAnalyzeSurfCollis%PartNumberSamp))
!  ALLOCATE(LastAnalyzeSurfCollis%Species(LastAnalyzeSurfCollis%PartNumberSamp))
!  LastAnalyzeSurfCollis%pushTimeStep = HUGE(LastAnalyzeSurfCollis%pushTimeStep)
!#if USE_MPI
!  IF (BCTotalNumberMPF.GT.0) THEN
!    ALLOCATE(sendbuf2(1:AnalyzeSurfCollis%Number(nSpecies+1)*8))
!    ALLOCATE(recvbuf2(1:TotalNumberMPF*8))
!    ! Fill sendbufer
!    counter2 = 0
!    DO iPart=1,AnalyzeSurfCollis%Number(nSpecies+1)
!      sendbuf2(counter2+1:counter2+6) = AnalyzeSurfCollis%Data(iPart,1:6)
!      sendbuf2(counter2+7)           = REAL(AnalyzeSurfCollis%Spec(iPart))
!      sendbuf2(counter2+8)           = REAL(AnalyzeSurfCollis%BCid(iPart))
!      counter2 = counter2 + 8
!    END DO
!    ! Distribute particles to all procs
!    counter2 = 0
!    DO iProc = 0, nProcessors-1
!      RecCount(iProc) = globalNum(iProc) * 8
!      Displace(iProc) = counter2
!      counter2 = counter2 + globalNum(iProc)*8
!    END DO
!    CALL MPI_ALLGATHERV(sendbuf2, 8*globalNum(myRank), MPI_DOUBLE_PRECISION, &
!      recvbuf2, RecCount, Displace, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERROR)
!    ! Add them to particle list
!    counter2 = -8 !moved increment before usage, thus: -8 instead of 0
!    DO counter = 1, LastAnalyzeSurfCollis%PartNumberSamp
!      IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts (differently in each proc. Could be changed)
!        DO !get random (equal!) position between 8*[0,TotalNumberMPF-1] and accept if .NOT.PartDone and with right BC
!          CALL RANDOM_NUMBER(RandVal)
!          counter2 = MIN(1+INT(RandVal*REAL(TotalNumberMPF)),TotalNumberMPF) !( MIN(1+INT(RandVal*REAL(TotalNumberMPF)),TotalNumberMPF) - 1) *8
!          IF (.NOT.PartDone(counter2) .AND. &
!            ( ANY(LastAnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(LastAnalyzeSurfCollis%BCs.EQ.INT(recvbuf2(8*counter2))) )) THEN
!            PartDone(counter2)=.TRUE.
!            counter2 = 8*(counter2-1)
!            EXIT
!          END IF
!        END DO
!      ELSE
!        counter2 = counter2 + 8
!      END IF
!      LastAnalyzeSurfCollis%WallState(:,counter) = recvbuf2(counter2+1:counter2+6)
!      LastAnalyzeSurfCollis%Species(counter) = INT(recvbuf2(counter2+7))
!      LastAnalyzeSurfCollis%pushTimeStep = MIN( LastAnalyzeSurfCollis%pushTimeStep &
!        , DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%WallState(4:6,counter)) )
!    END DO
!    DEALLOCATE(sendbuf2 &
!              ,recvbuf2 )
!  END IF
!#else
!  ! Add particle to list
!  counter2 = 0
!  DO counter = 1, LastAnalyzeSurfCollis%PartNumberSamp
!    IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts (differently for each proc. Could be changed)
!      DO !get random (equal!) position between [1,TotalNumberMPF] and accept if .NOT.PartDone and with right BC
!        CALL RANDOM_NUMBER(RandVal)
!        counter2 = MIN(1+INT(RandVal*REAL(TotalNumberMPF)),TotalNumberMPF)
!        IF (.NOT.PartDone(counter2) .AND. &
!          ( ANY(LastAnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(LastAnalyzeSurfCollis%BCs.EQ.AnalyzeSurfCollis%BCid(counter2)) )) THEN
!          PartDone(counter2)=.TRUE.
!          EXIT
!        END IF
!      END DO
!    ELSE
!      counter2 = counter2 + 1
!    END IF
!    LastAnalyzeSurfCollis%WallState(:,counter) = AnalyzeSurfCollis%Data(counter2,1:6)
!    LastAnalyzeSurfCollis%Species(counter) = AnalyzeSurfCollis%Spec(counter2)
!    LastAnalyzeSurfCollis%pushTimeStep = MIN( LastAnalyzeSurfCollis%pushTimeStep &
!      , DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%WallState(4:6,counter)) )
!  END DO
!#endif
!  IF (LastAnalyzeSurfCollis%pushTimeStep .LE. 0.) THEN
!    CALL Abort(&
!      __STAMP__,&
!      'Error with SFResampleAnalyzeSurfCollis. Something is wrong with velocities or NormVecOfWall!')
!  ELSE
!    LastAnalyzeSurfCollis%pushTimeStep = r_SF / LastAnalyzeSurfCollis%pushTimeStep !dt required for smallest projected velo to cross r_SF
!    LastAnalyzeSurfCollis%PartNumberDepo = NINT(BCTotalFlowrateMPF * LastAnalyzeSurfCollis%pushTimeStep)
!    SWRITE(*,'(A,E12.5,x,I0)') 'Total Flowrate and to be inserted number of MP for SFResampleAnalyzeSurfCollis: ' &
!      ,BCTotalFlowrateMPF, LastAnalyzeSurfCollis%PartNumberDepo
!    IF (LastAnalyzeSurfCollis%PartNumberDepo .GT. LastAnalyzeSurfCollis%PartNumberSamp) THEN
!      SWRITE(*,*) 'WARNING: PartNumberDepo .GT. PartNumberSamp!'
!    END IF
!    IF (LastAnalyzeSurfCollis%PartNumberDepo .GT. LastAnalyzeSurfCollis%PartNumThreshold) THEN
!      CALL Abort(&
!        __STAMP__,&
!        'Error with SFResampleAnalyzeSurfCollis: PartNumberDepo .gt. PartNumThreshold',LastAnalyzeSurfCollis%PartNumberDepo)
!    END IF
!  END IF
!END IF !SFResampleAnalyzeSurfCollis
!
!DO iSpec=1,nSpecies
!  ALLOCATE(PartData(offsetnPart(iSpec)+1:offsetnPart(iSpec)+locnPart(iSpec),PartDataSize))
!  DO iPart=1,locnPart(iSpec)
!    PartData(offsetnPart(iSpec)+iPart,1)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),1)
!    PartData(offsetnPart(iSpec)+iPart,2)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),2)
!    PartData(offsetnPart(iSpec)+iPart,3)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),3)
!    PartData(offsetnPart(iSpec)+iPart,4)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),4)
!    PartData(offsetnPart(iSpec)+iPart,5)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),5)
!    PartData(offsetnPart(iSpec)+iPart,6)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),6)
!    PartData(offsetnPart(iSpec)+iPart,7)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),7)
!    PartData(offsetnPart(iSpec)+iPart,8)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),8)
!    PartData(offsetnPart(iSpec)+iPart,9)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),9)
!    PartData(offsetnPart(iSpec)+iPart,10)=REAL(AnalyzeSurfCollis%BCid(SpeciesPositions(iSpec,iPart)))
!  END DO
!  WRITE(H5_Name,'(A,I3.3)') 'SurfCollisData_Spec',iSpec
!  IF(minnParts(iSpec).EQ.0)THEN
!    CALL WriteArrayToHDF5(DataSetName=TRIM(H5_Name), rank=2,&
!                          nValGlobal=(/nPart_glob(iSpec),PartDataSize/),&
!                          nVal=      (/locnPart(iSpec),PartDataSize  /),&
!                          offset=    (/offsetnPart(iSpec) , 0  /),&
!                          collective=.FALSE., RealArray=PartData)
!  ELSE
!    CALL WriteArrayToHDF5(DataSetName=TRIM(H5_Name), rank=2,&
!                          nValGlobal=(/nPart_glob(iSpec),PartDataSize/),&
!                          nVal=      (/locnPart(iSpec),PartDataSize  /),&
!                          offset=    (/offsetnPart(iSpec) , 0  /),&
!                          collective=.TRUE., RealArray=PartData)
!  END IF
!  DEALLOCATE(PartData)
!END DO !iSpec
!
!CALL CloseDataFile()
!DEALLOCATE(locnPart &
!          ,offsetnPart &
!          ,nPart_glob &
!          ,minnParts &
!          ,iPartCount )
!DEALLOCATE(SpeciesPositions)
!DEALLOCATE(StrVarNames)
!
!END SUBROUTINE WriteAnalyzeSurfCollisToHDF5


END MODULE MOD_Particle_Erosion_Analyze
