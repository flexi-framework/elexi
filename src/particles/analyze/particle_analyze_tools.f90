!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
#include "particle.h"

!===================================================================================================================================
!> Contains helper routines for statistical analysis of particle behavior
!===================================================================================================================================
MODULE MOD_Particle_Analyze_Tools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: CalcEkinPart
PUBLIC:: ParticleRecordPath
PUBLIC:: ParticleRecord
PUBLIC:: TrackingParticlePath
!===================================================================================================================================

CONTAINS

PURE FUNCTION CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,           ONLY: PI
USE MOD_Particle_Vars,          ONLY : PartState, PartSpecies, Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: iPart
REAL                               :: CalcEkinPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

CalcEkinPart = 0.5*PI/6*Species(PartSpecies(iPart))%DensityIC*PartState(PART_DIAM,iPart)**3*&
               DOT_PRODUCT(PartState(PART_VELV,iPart),PartState(PART_VELV,iPart))

END FUNCTION CalcEkinPart


SUBROUTINE TrackingParticlePath()
!===================================================================================================================================
! Outputs the particle position and velocity at every time step to determine the absolute or relative (since emission) path
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars,  ONLY: PartPath,doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Vars,          ONLY: PartState,PDM,LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iPart
!===================================================================================================================================

! No BC interaction expected, so path can be calculated here. Periodic BCs are ignored purposefully
IF (doParticleDispersionTrack) THEN
  DO iPart = 1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) PartPath(1:3,iPart) = PartPath(1:3,iPart) + ABS(PartState(1:3,iPart) - LastPartPos(1:3,iPart))
  END DO
ELSEIF (doParticlePathTrack) THEN
  DO iPart = 1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) PartPath(1:3,iPart) = PartPath(1:3,iPart) +    (PartState(1:3,iPart) - LastPartPos(1:3,iPart))
  END DO
END IF

END SUBROUTINE TrackingParticlePath


SUBROUTINE ParticleRecordPath()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals        ,ONLY: VECNORM
USE MOD_HDF5_Output             ,ONLY: WriteArray
USE MOD_Output_Vars             ,ONLY: WriteStateFiles
USE MOD_Particle_Analyze_Vars   ,ONLY: RPP_Plane,RecordPart,RPP_Records,RPP_Records_Glob
USE MOD_Particle_Boundary_Vars  ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Vars           ,ONLY: PartState,PDM,LastPartPos,PartSpecies
USE MOD_Particle_Vars           ,ONLY: doPartIndex,PartIndex
USE MOD_Particle_Vars           ,ONLY: PartReflCount
USE MOD_TimeDisc_Vars           ,ONLY: t,dt,currentStage,RKC,nRKStages
USE MOD_Utils                   ,ONLY: ALMOSTZERO
#if USE_MPI
USE MOD_Particle_Analyze_Vars   ,ONLY: RPP_MPI_Request
#endif /*USE_MPI*/
#if ANALYZE_RHS
USE MOD_Particle_Vars           ,ONLY: Pt_ext
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPart,iRecord,m
REAL                           :: t_loc
! RecordPlane
REAL                           :: PartTrajectory(3),lengthPartTrajectory,Inter1(3)
REAL                           :: locOrigin(1:3),locNormVec(1:3),locDistance
REAL                           :: coeffA,locPlaneDistance,alpha,alphaNorm
!===================================================================================================================================
IF (.NOT.WriteStateFiles) RETURN

!IF(RPP_Type.EQ.'plane')THEN
DO iRecord = 1,RecordPart
  DO iPart=1,PDM%ParticleVecLength
    IF (.NOT. PDM%ParticleInside(iPart)) CYCLE

    ! Compute particle trajectory
    PartTrajectory       = PartState(1:3,iPart) - LastPartPos(1:3,iPart)
    lengthPartTrajectory = VECNORM(PartTrajectory(1:3))
    PartTrajectory       = PartTrajectory/lengthPartTrajectory

    ! Compute planar rect intersection
    locOrigin   = RPP_Plane(iRecord)%pos
    locNormVec  = RPP_Plane(iRecord)%dir
    locDistance = RPP_Plane(iRecord)%dist
    coeffA      = DOT_PRODUCT(locNormVec,PartTrajectory)

    ! Particle moving parallel to plane
    IF (ALMOSTZERO(coeffA)) CYCLE

    ! Difference between SideDistance (distance from origin to side) and the dot product is the distance of the particle to the side
    locPlaneDistance = locDistance-DOT_PRODUCT(LastPartPos(1:3,iPart),locNormVec)

    ! Length of particle vector until side intersection in physical space
    alpha     = locPlaneDistance/coeffA

    ! Calculate normalized alpha, i.e. length of particle vector until intersection in reference element
    alphaNorm = alpha/lengthPartTrajectory

    IF (alphaNorm.GE.0 .AND. alphaNorm.LE.1.) THEN
      ! Calculate intersection point
      Inter1 = LastPartPos(1:3,iPart) + alphaNorm*PartTrajectory*lengthPartTrajectory

      ! Calculate exact impact time
      IF (CurrentStage.EQ.1) THEN
        t_loc = t                                                       & ! current physical time
              +  RKc(2)                                *dt*alphaNorm      ! relative time till intersection
      ELSE IF (CurrentStage.GT.1 .AND. currentStage.LT.nRKStages) THEN
        t_loc = t                                                       & ! current physical time
              +  RKc(CurrentStage)                     *dt              & ! current stage time
              + (RKc(CurrentStage+1)-RKc(currentStage))*dt*alphaNorm      ! relative time till intersection
      ELSE ! nRKStages
        t_loc = t                                                       & ! current physical time
              +  RKc(CurrentStage)                     *dt              & ! current stage time
              + (1.                 -RKc(currentStage))*dt*alphaNorm      ! relative time till intersection
      END IF

      RPP_Plane(iRecord)%RPP_Records = RPP_Plane(iRecord)%RPP_Records+1
      ! Part intersection point
      RPP_Plane(iRecord)%RPP_Data(1:3,RPP_Plane(iRecord)%RPP_Records) = Inter1
      ! Part velocity
      RPP_Plane(iRecord)%RPP_Data(4:6,RPP_Plane(iRecord)%RPP_Records) = PartState(4:6,iPart)
      ! dp
      RPP_Plane(iRecord)%RPP_Data(7  ,RPP_Plane(iRecord)%RPP_Records) = PartState(PART_DIAM,iPart)
      ! Species
      RPP_Plane(iRecord)%RPP_Data(8  ,RPP_Plane(iRecord)%RPP_Records) = PartSpecies(iPart)
      m = 8
      ! Reflection
      IF (doParticleReflectionTrack) THEN; m = m + 1
        RPP_Plane(iRecord)%RPP_Data(m,RPP_Plane(iRecord)%RPP_Records) = PartReflCount(iPart)
      END IF
      ! Time
      m = m + 1
      RPP_Plane(iRecord)%RPP_Data(m  ,RPP_Plane(iRecord)%RPP_Records) = t_loc
#if USE_SPHERICITY
      ! Sphericity
      m = m + 1
      RPP_Plane(iRecord)%RPP_Data(m  ,RPP_Plane(iRecord)%RPP_Records) = PartState(PART_SPHE,iPart)
#endif
      ! Index
      IF(doPartIndex) THEN; m = m + 1
        RPP_Plane(iRecord)%RPP_Data(m,RPP_Plane(iRecord)%RPP_Records) = PartIndex(iPart)
      END IF
#if ANALYZE_RHS
      m = m + 1
      RPP_Plane(iRecord)%RPP_Data(m:m+18,RPP_Plane(iRecord)%RPP_Records) = Pt_ext(:,iPart)
#endif /*ANALYZE_RHS*/
    END IF
  END DO
  RPP_Records(iRecord) = RPP_Plane(iRecord)%RPP_Records
END DO

#if USE_MPI
RPP_Records_Glob = 0
CALL MPI_IALLREDUCE(RPP_Records,RPP_Records_Glob,RecordPart,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,RPP_MPI_Request,iError)
#else
RPP_Records_Glob = RPP_Records
#endif /*USE_MPI*/

END SUBROUTINE ParticleRecordPath


SUBROUTINE ParticleRecord(OutputTime,writeToBinary)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Output             ,ONLY: WriteAttribute,WriteArray
USE MOD_IO_HDF5                 ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Output_Vars             ,ONLY: WriteStateFiles
USE MOD_Particle_Analyze_Vars   ,ONLY: RPP_MaxBufferSize,RPP_Plane,RecordPart,RPP_nVarNames,RPP_Records_Glob
USE MOD_Particle_Boundary_Vars  ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Vars           ,ONLY: Species,nSpecies,doPartIndex
#if USE_MPI
USE MOD_Particle_Analyze_Vars   ,ONLY: RPP_MPI_Request
USE MOD_Particle_HDF5_Output    ,ONLY: DistributedWriteArray
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: OutputTime
LOGICAL,OPTIONAL,INTENT(IN)    :: writeToBinary
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: RPP_Output
INTEGER                        :: iRecord,iSpecies,m
CHARACTER(LEN=200)             :: FileName_loc
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: locRPP,offsetRPP,RPP_glob
REAL                           :: SphericityIC(nSpecies)
#if USE_EXTEND_RHS
INTEGER                        :: ForceIC(6,nSpecies)
#else
INTEGER                        :: ForceIC(1,nSpecies)
#endif
#if USE_MPI
INTEGER                        :: sendbuf(2),recvbuf(2)
! INTEGER                        :: nRecords(0:nProcessors-1)
#endif
CHARACTER(LEN=32)              :: tmpStr
! Timers
REAL                           :: StartT,EndT
!===================================================================================================================================
IF (.NOT.WriteStateFiles) RETURN

#if USE_MPI
IF (RPP_MPI_Request.NE.MPI_REQUEST_NULL) THEN
  CALL MPI_WAIT(RPP_MPI_Request,MPI_STATUS_IGNORE,iError)
END IF
#endif /*USE_MPI*/

RPP_Output = .FALSE.

DO iRecord = 1,RecordPart
  IF((RPP_Records_Glob(iRecord) .GE. RPP_MaxBufferSize .OR. PRESENT(writeToBinary)))THEN

    ! Cycle if we reached analyze_dt but no record parts
    IF(RPP_Records_Glob(iRecord).EQ.0) CYCLE

    RPP_Output = .TRUE.
    locRPP     = RPP_Plane(iRecord)%RPP_Records
    RPP_glob   = RPP_Records_Glob(iRecord)

    !>> Sum up particles from the other procs
#if USE_MPI
    sendbuf(1) = locRPP
    recvbuf    = 0
    CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
    !>> Offset of each proc is the sum of the particles on the previous procs
    offsetRPP  = recvbuf(1)
    sendbuf(1) = recvbuf(1)+locRPP
    !>> Last proc knows the global number
    CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER,nProcessors-1,MPI_COMM_FLEXI,iError)
    !>> Gather the global number and communicate to root (MPIRank.EQ.0)
    ! RPP_glob    = sendbuf(1)
    ! CALL MPI_GATHER(locRPP,1,MPI_INTEGER,nRecords,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
#else
    offsetRPP  = 0
    ! RPP_glob   = locRPP
#endif

    ALLOCATE(StrVarNames(RPP_nVarNames))
    StrVarNames(1)   ='PartPosX'
    StrVarNames(2)   ='PartPosY'
    StrVarNames(3)   ='PartPosZ'
    StrVarNames(4)   ='VelocityX'
    StrVarNames(5)   ='VelocityY'
    StrVarNames(6)   ='VelocityZ'
    StrVarNames(7)   ='PartDiam'
    StrVarNames(8)   ='Species'
    m = 8
    IF (doParticleReflectionTrack) THEN; m = m + 1
      StrVarNames(m) ='ReflectionCount'
    END IF
    m = m + 1
    StrVarNames(m)   ='Time'
#if USE_SPHERICITY
    m = m + 1
    StrVarNames(m)   ='Sphericity'
#endif
    IF(doPartIndex) THEN; m = m + 1
      StrVarNames(m) ='Index'
    END IF
#if ANALYZE_RHS
    m = m + 1; StrVarNames(m)   ='FdX'
    m = m + 1; StrVarNames(m)   ='FdY'
    m = m + 1; StrVarNames(m)   ='FdZ'
    m = m + 1; StrVarNames(m)   ='FvX'
    m = m + 1; StrVarNames(m)   ='FvY'
    m = m + 1; StrVarNames(m)   ='FvZ'
    m = m + 1; StrVarNames(m)   ='FuX'
    m = m + 1; StrVarNames(m)   ='FuY'
    m = m + 1; StrVarNames(m)   ='FuZ'
    m = m + 1; StrVarNames(m)   ='FlX'
    m = m + 1; StrVarNames(m)   ='FlY'
    m = m + 1; StrVarNames(m)   ='FlZ'
    m = m + 1; StrVarNames(m)   ='FmX'
    m = m + 1; StrVarNames(m)   ='FmY'
    m = m + 1; StrVarNames(m)   ='FmZ'
    m = m + 1; StrVarNames(m)   ='FbX'
    m = m + 1; StrVarNames(m)   ='FbY'
    m = m + 1; StrVarNames(m)   ='FbZ'
#endif

    ForceIC = 0

    WRITE(UNIT=tmpStr,FMT='(I0)') iRecord
    FileName_loc = TRIM(TIMESTAMP('recordpoints/recordpoints_part'//TRIM(ADJUSTL(tmpStr)),OutputTime))//'.h5'
    ! SWRITE(UNIT_stdOut,*)' Opening file '//TRIM(FileName_loc)

    IF(MPIRoot)THEN
      WRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' WRITE PARTICLE RECORD  PLANE TO HDF5 FILE...'
      GETTIME(startT)

      CALL OpenDataFile(FileName_loc,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttribute(File_ID,'File_Type'   ,1            ,StrScalar=(/CHARACTER(LEN=255)::'RecordPlane'/))
      CALL WriteAttribute(File_ID,'VarNamesPart',RPP_nVarNames,StrArray =StrVarNames)
      CALL WriteAttribute(File_ID,'nSpecies'    ,1            ,IntScalar=nSpecies)
      DO iSpecies=1,nSpecies
        SphericityIC(iSpecies) = Species(iSpecies)%SphericityIC
        IF(Species(iSpecies)%RHSMethod .NE. RHS_TRACER) ForceIC(1,iSpecies)    = 1
#if USE_EXTEND_RHS
        IF(Species(iSpecies)%CalcSaffmanForce)    ForceIC(2,iSpecies) = 1
        IF(Species(iSpecies)%CalcBassetForce)     ForceIC(3,iSpecies) = 1
        IF(Species(iSpecies)%CalcVirtualMass)     ForceIC(4,iSpecies) = 1
        IF(Species(iSpecies)%CalcUndisturbedFlow) ForceIC(5,iSpecies) = 1
        IF(Species(iSpecies)%CalcMagnusForce)     ForceIC(6,iSpecies) = 1
#endif
      END DO
#if !USE_SPHERICITY
      CALL WriteAttribute(File_ID,'SphericityIC',nSpecies,RealArray=SphericityIC)
#endif
      CALL WriteAttribute(File_ID,'DragForce'   ,nSpecies,IntArray=ForceIC(1,:))
#if USE_EXTEND_RHS
      CALL WriteAttribute(File_ID,'SaffmanForce',nSpecies,IntArray=ForceIC(2,:))
      CALL WriteAttribute(File_ID,'BassForce'   ,nSpecies,IntArray=ForceIC(3,:))
      CALL WriteAttribute(File_ID,'VirtForce'   ,nSpecies,IntArray=ForceIC(4,:))
      CALL WriteAttribute(File_ID,'UndiForce'   ,nSpecies,IntArray=ForceIC(5,:))
      CALL WriteAttribute(File_ID,'MagnusForce' ,nSpecies,IntArray=ForceIC(6,:))
#endif
      CALL CloseDataFile()
    END IF ! MPIRoot

#if USE_MPI
    CALL DistributedWriteArray(FileName_loc                                       ,&
                               DataSetName  = 'RecordData'                        ,&
                               rank         = 2                                   ,&
                               nValGlobal   = (/RPP_nVarNames ,RPP_Glob        /) ,&
                               nVal         = (/RPP_nVarNames ,locRPP          /) ,&
                               offset       = (/0             ,offsetRPP       /) ,&
                               collective   = .FALSE.                             ,&
                               offSetDim    = 2                                   ,&
                               communicator = MPI_COMM_FLEXI                      ,&
                               RealArray    = RPP_Plane(iRecord)%RPP_Data(1:RPP_nVarNames,1:locRPP))
#else
    CALL OpenDataFile(FileName_loc,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArray(           DataSetName  = 'RecordData'                        ,&
                               rank         = 2                                   ,&
                               nValGlobal   = (/RPP_nVarNames ,RPP_Glob        /) ,&
                               nVal         = (/RPP_nVarNames ,locRPP          /) ,&
                               offset       = (/0             ,offsetRPP       /) ,&
                               collective   = .TRUE.                              ,&
                               RealArray    = RPP_Plane(iRecord)%RPP_Data(1:RPP_nVarNames,1:locRPP))
    CALL CloseDataFile()
#endif /*MPI*/

    RPP_Plane(iRecord)%RPP_Data=0.0

    IF (MPIRoot) THEN
      GETTIME(EndT)
      CALL DisplayMessageAndTime(EndT-StartT,'DONE!',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
    END IF
    DEALLOCATE(StrVarNames)

    RPP_Plane(iRecord)%RPP_Records=0
  END IF
END DO

IF (RPP_Output .AND. MPIRoot) WRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE ParticleRecord

END MODULE MOD_Particle_Analyze_Tools
