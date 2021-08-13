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

!===================================================================================================================================
! Module for particle analysis and output
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CalcSurfaceValues
  MODULE PROCEDURE CalcSurfaceValues
END INTERFACE

INTERFACE WriteBoundaryParticleToHDF5
  MODULE PROCEDURE WriteBoundaryParticleToHDF5
END INTERFACE

PUBLIC :: CalcSurfaceValues
PUBLIC :: WriteBoundaryParticleToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcSurfaceValues(restart_opt,remap_opt)
!===================================================================================================================================
!> Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Restart_Vars               ,ONLY: DoRestart,RestartTime
USE MOD_Analyze_Vars               ,ONLY: Analyze_dt
USE MOD_Mesh_Vars                  ,ONLY: MeshFile
USE MOD_Timedisc_Vars              ,ONLY: t
USE MOD_Restart_Vars               ,ONLY: RestartTime
USE MOD_Particle_Analyze_Vars      ,ONLY: TimeSample
USE MOD_Particle_Boundary_Vars     ,ONLY: WriteMacroSurfaceValues,MacroValSampTime
USE MOD_Particle_Boundary_Vars     ,ONLY: nSurfSample
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfOnNode,SurfSideArea
USE MOD_Particle_Boundary_Vars     ,ONLY: nComputeNodeSurfSides,SampWallState_Shared
USE MOD_Particle_Boundary_Sampling ,ONLY: WriteSurfSample
USE MOD_Particle_Boundary_Vars     ,ONLY: MacroSurfaceVal,MacroSurfaceSpecVal,nImpactVars
USE MOD_Particle_Vars              ,ONLY: nSpecies
USE MOD_CalcWallParticles_Vars
#if USE_MPI
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfOnNode
USE MOD_Particle_MPI_Boundary_Sampling,ONLY: ExchangeSurfData
USE MOD_Particle_MPI_Shared_Vars   ,ONLY: MPI_COMM_LEADERS_SURF
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL      :: restart_opt   !routine was called during tstep (i.e. before iter=iter+1, t=t+dt...)
CHARACTER(LEN=*),INTENT(IN),OPTIONAL::remap_opt     !routine was called from posti. Change output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iSpec,iSurfSide
INTEGER                            :: p,q
INTEGER                            :: nShift,nShiftRHS
REAL                               :: ActualTime
!===================================================================================================================================

! Set ActualTime to current time
ActualTime = t

IF (WriteMacroSurfaceValues) THEN
  TimeSample       = t - MacroValSampTime
  MacroValSampTime = t
END IF

! Update values if we are called from a restart
IF (PRESENT(restart_opt)) THEN
  IF (restart_opt) THEN
    TimeSample = Analyze_dt
    t          = MERGE(RestartTime,0.,DoRestart)
    ActualTime = t
  END IF
! Avoid division by zero if no simulation time has passed
ELSE
  IF(ALMOSTZERO(TimeSample)) RETURN
END IF

! Do not try to record impacts if there are no walls on the current proc
IF(.NOT.SurfOnNode) RETURN

#if USE_MPI
CALL ExchangeSurfData()

! Only surface sampling leaders take part in the remainder of this routine
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
#endif /*USE_MPI*/

! Allocate N+1 Species to have space for average
IF (nSpecies.EQ.1) THEN
  ALLOCATE(MacroSurfaceVal(nImpactVars-1,1:nSurfSample,1:nSurfSample,nComputeNodeSurfSides))
ELSE
  ALLOCATE(MacroSurfaceVal((nImpactVars-1)*(nSpecies+1),1:nSurfSample,1:nSurfSample,nComputeNodeSurfSides))
END IF
ALLOCATE(MacroSurfaceSpecVal(1,1:nSurfSample,1:nSurfSample,nComputeNodeSurfSides,nSpecies))

MacroSurfaceVal    = 0.
MacroSurfaceSpecVal= 0.

!> Impact sampling
iSpec = 1

ASSOCIATE(SampWallState => SampWallState_Shared)
!===================================================================================================================================
DO iSurfSide = 1,nComputeNodeSurfSides; DO q = 1,nSurfSample; DO p = 1,nSurfSample
  !---- 1. - .. / Impact Counter
  MacroSurfaceVal(1    ,p,q,iSurfSide) =       SampWallState(1 ,p,q,iSurfSide)
  !---- 2. - .. / Impact Counter per AREA
  MacroSurfaceVal(2    ,p,q,iSurfSide) =       SampWallState(1 ,p,q,iSurfSide) / SurfSideArea(p,q,iSurfSide)
  !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
  MacroSurfaceVal(3    ,p,q,iSurfSide) =       SampWallState(2 ,p,q,iSurfSide)
  MacroSurfaceVal(4    ,p,q,iSurfSide) = MERGE(SampWallState(3 ,p,q,iSurfSide),0.,SampWallState(3,p,q,iSurfSide).LT. HUGE(1.))
  MacroSurfaceVal(5    ,p,q,iSurfSide) = MERGE(SampWallState(4 ,p,q,iSurfSide),0.,SampWallState(4,p,q,iSurfSide).GT.-HUGE(1.))
  MacroSurfaceVal(6    ,p,q,iSurfSide) =       SampWallState(6 ,p,q,iSurfSide)
  !---- 7. - 10 / Impact angle (mean, min, max, variance)
  MacroSurfaceVal(7    ,p,q,iSurfSide) =       SampWallState(7 ,p,q,iSurfSide)
  MacroSurfaceVal(8    ,p,q,iSurfSide) = MERGE(SampWallState(8 ,p,q,iSurfSide),0.,SampWallState(3,p,q,iSurfSide).LT. HUGE(1.))
  MacroSurfaceVal(9    ,p,q,iSurfSide) = MERGE(SampWallState(9 ,p,q,iSurfSide),0.,SampWallState(4,p,q,iSurfSide).GT.-HUGE(1.))
  MacroSurfaceVal(10   ,p,q,iSurfSide) =       SampWallState(11,p,q,iSurfSide)
  !---- 11 - 13 / Sampling Current Forces at walls
  MacroSurfaceVal(11:13,p,q,iSurfSide) =       SampWallState(12:14,p,q,iSurfSide) / (SurfSideArea(p,q,iSurfSide) * TimeSample)
  !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  SampWall(iSurfSide)%State(12:14,p,q) = 0.
  !---- 14 - 16 / Sampling Average Forces at walls
  MacroSurfaceVal(14:16,p,q,iSurfSide) =       SampWallState(15:17,p,q,iSurfSide) / (SurfSideArea(p,q,iSurfSide) * t)
END DO; END DO; END DO

!---- Multiple species. All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
IF (nSpecies.GT.1) THEN
  DO iSurfSide = 1,nComputeNodeSurfSides; DO q = 1,nSurfSample; DO p = 1,nSurfSample; DO iSpec=1,nSpecies
    nShift    = iSpec * (nImpactVars-1)
    nShiftRHS = iSpec *  nImpactVars
    !---- 1. - .. / Impact Counter
    MacroSurfaceVal(1+nShift ,p,q,iSurfSide) =       SampWallState(1 +nShiftRHS,p,q,iSurfSide)
    MacroSurfaceSpecVal(1    ,p,q,iSurfSide,iSpec) = SampWallState(1 +nShiftRHS,p,q,iSurfSide) / TimeSample
    !---- 2. - .. / Impact Counter per AREA
    MacroSurfaceVal(2+nShift ,p,q,iSurfSide) =       SampWallState(1 +nShiftRHS,p,q,iSurfSide) / SurfSideArea(p,q,iSurfSide)
    !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
    MacroSurfaceVal(3+nShift ,p,q,iSurfSide) =       SampWallState(2 +nShiftRHS,p,q,iSurfSide)
    MacroSurfaceVal(4+nShift ,p,q,iSurfSide) = MERGE(SampWallState(3 +nShiftRHS,p,q,iSurfSide),0.,SampWallState(3+nShiftRHS,p,q,iSurfSide).LT. HUGE(1.))
    MacroSurfaceVal(5+nShift ,p,q,iSurfSide) = MERGE(SampWallState(4 +nShiftRHS,p,q,iSurfSide),0.,SampWallState(4+nShiftRHS,p,q,iSurfSide).GT.-HUGE(1.))
    MacroSurfaceVal(6+nShift ,p,q,iSurfSide) =       SampWallState(6 +nShiftRHS,p,q,iSurfSide)
    !---- 7. - 10 / Impact angle (mean, min, max, variance)
    MacroSurfaceVal(7+nShift ,p,q,iSurfSide) =       SampWallState(7 +nShiftRHS,p,q,iSurfSide)
    MacroSurfaceVal(8+nShift ,p,q,iSurfSide) = MERGE(SampWallState(8 +nShiftRHS,p,q,iSurfSide),0.,SampWallState(8+nShiftRHS,p,q,iSurfSide).LT. HUGE(1.))
    MacroSurfaceVal(9+nShift ,p,q,iSurfSide) = MERGE(SampWallState(9 +nShiftRHS,p,q,iSurfSide),0.,SampWallState(9+nShiftRHS,p,q,iSurfSide).GT.-HUGE(1.))
    MacroSurfaceVal(10+nShift,p,q,iSurfSide) =       SampWallState(11+nShiftRHS,p,q,iSurfSide)
    !---- 11 - 13 / Sampling Current Forces at walls
    MacroSurfaceVal(11+nShift:13+nShift,p,q,iSurfSide) = SampWallState(12+nShiftRHS:14+nShiftRHS,p,q,iSurfSide)             &
                                                       / (SurfSideArea(p,q,iSurfSide) * TimeSample)
    !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!    SampWall(iSurfSide)%State(12+nShiftRHS,p,q) = 0.
!     SampWall(iSurfSide)%State(13+nShiftRHS,p,q) = 0.
!     SampWall(iSurfSide)%State(14+nShiftRHS,p,q) = 0.
    !---- 14 - 16 / Sampling Average Forces at walls
    MacroSurfaceVal(14+nShift:16+nShift,p,q,iSurfSide) = SampWallState(15+nShiftRHS:17+nShiftRHS,p,q,iSurfSide)             &
                                                       / (SurfSideArea(p,q,iSurfSide) * TimeSample)
  END DO; END DO; END DO; END DO
ELSE
  DO iSurfSide = 1,nComputeNodeSurfSides; DO q = 1,nSurfSample; DO p = 1,nSurfSample
    !---- 1. - .. / Impact Counter
    MacroSurfaceSpecVal(1,p,q,iSurfSide,iSpec) = SampWallState(1,p,q,iSurfSide) / TimeSample
  END DO; END DO; END DO
END IF
END ASSOCIATE

CALL WriteSurfSample(TRIM(MeshFile),ActualTime,remap_opt)

! Only deallocate if we don't need the values for wall calculations
IF (.NOT.doCalcWallParticles) THEN
    DEALLOCATE(MacroSurfaceVal,MacroSurfaceSpecVal)
ELSE
    DEALLOCATE(MacroSurfaceSpecVal)
END IF

END SUBROUTINE CalcSurfaceValues


SUBROUTINE WriteBoundaryParticleToHDF5(OutputTime)
!===================================================================================================================================
! Write data of impacting particles on specific boundary conditions of .h5 file (position, velocity, species ID, kinetic energy, time of impact, impact obliqueness angle)
!===================================================================================================================================
! MODULES
USE MOD_Globals
! USE MOD_PreProc
USE MOD_HDF5_WriteArray        ,ONLY: WriteArray,GatheredWriteArray
USE MOD_HDF5_Output            ,ONLY: WriteAttribute
USE MOD_IO_HDF5                ,ONLY: GatheredWrite
USE MOD_IO_HDF5                ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Output_Vars            ,ONLY: ProjectName,WriteStateFiles
USE MOD_Particle_Analyze_Vars  ,ONLY: doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength
USE MOD_Particle_Boundary_Vars ,ONLY: ImpactDataSize,ImpactnGlob
USE MOD_Particle_Vars          ,ONLY: doPartIndex
#if USE_MPI
! USE MOD_Particle_Boundary_Vars ,ONLY: MPI_COMM_IMPACT
USE MOD_Particle_HDF5_Output   ,ONLY: DistributedWriteArray
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,   INTENT(IN)             :: OutputTime            !< time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName,FileString
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
LOGICAL                        :: reSwitch
REAL                           :: startT,endT
INTEGER                        :: ImpactnLoc,ImpactOffset,dims(2)
#if USE_MPI
INTEGER                        :: sendbuf(2),recvbuf(2)
! INTEGER                        :: nImpacts(0:nProcessors-1)
#endif
!===================================================================================================================================
IF (.NOT.WriteStateFiles) RETURN

! Find amount of recorded impacts on current proc
ImpactnLoc  = PartStateBoundaryVecLength

!>> Sum up particles from the other procs
#if USE_MPI
sendbuf(1) = ImpactnLoc
recvbuf    = 0
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
!>> Offset of each proc is the sum of the particles on the previous procs
ImpactOffset = recvbuf(1)
sendbuf(1)   = recvbuf(1) + ImpactnLoc
!>> Last proc knows the global number
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER,nProcessors-1,MPI_COMM_FLEXI,iError)
!>> Gather the global number and communicate to root (MPIRank.EQ.0)
ImpactnGlob  = sendbuf(1)
#else
ImpactOffset   = 0
ImpactnGlob    = PartStateBoundaryVecLength
#endif

IF (MPIRoot) THEN
  IF (ImpactnGlob.NE.0) THEN
    WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='NO')' WRITE PARTICLE IMPACTS STATE TO HDF5 FILE...'
    GETTIME(startT)
  ELSE
    WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='NO')' SKIP  PARTICLE IMPACTS STATE TO HDF5 FILE...'
    GETTIME(startT)
  END IF
END IF

! Regenerate state file skeleton
FileName   = TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime))
FileString = TRIM(FileName)//'.h5'

reSwitch = .FALSE.
IF (gatheredWrite) THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch      = .TRUE.
  gatheredWrite = .FALSE.
END IF

! Array for impact tracking variable
ALLOCATE(StrVarNames(ImpactDataSize))
StrVarNames(1)  = 'ParticlePositionX'
StrVarNames(2)  = 'ParticlePositionY'
StrVarNames(3)  = 'ParticlePositionZ'
StrVarNames(4)  = 'VelocityX'
StrVarNames(5)  = 'VelocityY'
StrVarNames(6)  = 'VelocityZ'
StrVarNames(7)  = 'Species'
StrVarNames(8)  = 'BoundaryNumber'
StrVarNames(9)  = 'ImpactTime'
StrVarNames(10) = 'ReflectionCount'
StrVarNames(11) = 'E_kin_impact'
StrVarNames(12) = 'E_kin_reflected'
StrVarNames(13) = 'Alpha_impact'
StrVarNames(14) = 'Alpha_reflected'
IF (doPartIndex) StrVarNames(15)= 'Index'
IF (doParticleDispersionTrack) THEN
  StrVarNames(ImpactDataSize-2) = 'PartPathAbsX'
  StrVarNames(ImpactDataSize-1) = 'PartPathAbsY'
  StrVarNames(ImpactDataSize  ) = 'PartPathAbsZ'
END IF
IF (doParticlePathTrack) THEN
  StrVarNames(ImpactDataSize-2) = 'PartPathX'
  StrVarNames(ImpactDataSize-1) = 'PartPathY'
  StrVarNames(ImpactDataSize  ) = 'PartPathZ'
END IF

IF (MPIRoot) THEN
  CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'VarNamesImpactTracking',ImpactDataSize,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

!> Writing empty arrays can cause problems with HDF5
! IF (ImpactnGlob.EQ.0) THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
!   IF (MPIRoot) THEN ! only root writes the container
!     CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
!     CALL WriteArray(         DataSetName = 'PartData'                      ,&
!                              rank        = 2                               ,&
!                              nValGlobal  = (/ImpactDataSize,ImpactnGlob/)  ,&
!                              nVal        = (/ImpactDataSize,ImpactnLoc /)  ,&
!                              offset      = (/ 0             , 0        /)  ,&
!                              collective  = .FALSE.                         ,&
!                              RealArray   = PartStateBoundary(1:ImpactDataSize,1:ImpactnLoc))
!     CALL CloseDataFile()
!     GETTIME(EndT)
!     WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES') 'DONE  [',EndT-StartT,'s] [NO IMPACTS WRITTEN]'
!   END IF ! MPIRoot
! ELSE
#if USE_MPI
  CALL DistributedWriteArray(FileString                                    ,&
                             DataSetName  = 'ImpactData'                   ,&
                             rank         = 2                              ,&
                             nValGlobal   = (/ImpactDataSize,ImpactnGlob /),&
                             nVal         = (/ImpactDataSize,ImpactnLoc  /),&
                             offset       = (/0             ,ImpactOffset/),&
                             collective   = .FALSE.                        ,&
                             offSetDim    = 2                              ,&
                             communicator = PartMPI%COMM                   ,&
                             RealArray    = PartStateBoundary(1:ImpactDataSize,1:ImpactnLoc))
#else
  CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArray(           DataSetName  = 'ImpactData'                   ,&
                             rank         = 2                              ,&
                             nValGlobal   = (/ImpactDataSize,ImpactnLoc/)  ,&
                             nVal         = (/ImpactDataSize,ImpactnLoc/)  ,&
                             offset       = (/0             ,0         /)  ,&
                             collective   = .TRUE.                         ,&
                             RealArray    = PartStateBoundary(1:ImpactDataSize,1:ImpactnLoc))
  CALL CloseDataFile()
#endif /*USE_MPI*/
! END IF ! ImpactnGlob.EQ.0

! reswitch
IF (reSwitch) gatheredWrite = .TRUE.

DEALLOCATE(StrVarNames)

! Check if PartStateBoundary was grown beyond initial size
dims = SHAPE(PartStateBoundary)
! Re-allocate PartStateBoundary for a small number of particles and double the array size each time the
! maximum is reached
IF (dims(2).GT.10) THEN
  DEALLOCATE(PartStateBoundary)
  ALLOCATE(PartStateBoundary(1:ImpactDataSize,1:10))
END IF

! Nullify and reset boundary parts container after write out
PartStateBoundaryVecLength = 0
PartStateBoundary          = 0.

IF (MPIROOT) THEN
!  CALL MarkWriteSuccessfull(FileString)
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES') 'DONE  [',EndT-StartT,'s]'
END IF

END SUBROUTINE WriteBoundaryParticleToHDF5

END MODULE MOD_Particle_Boundary_Analyze
