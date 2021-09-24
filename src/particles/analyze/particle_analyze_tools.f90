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
#include "particle.h"


!===================================================================================================================================
!> Contains helper routines for statistical analysis of particle behavior
!===================================================================================================================================
MODULE MOD_Particle_Analyze_Tools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CalcEkinPart
  MODULE PROCEDURE CalcEkinPart
END INTERFACE

INTERFACE TrackingParticlePath
  MODULE PROCEDURE TrackingParticlePath
END INTERFACE

INTERFACE ParticleRecord
  MODULE PROCEDURE ParticleRecord
END INTERFACE

PUBLIC :: CalcEkinPart
PUBLIC :: ParticleRecord
PUBLIC :: TrackingParticlePath
!===================================================================================================================================

CONTAINS

PURE FUNCTION CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Particle_Globals,       ONLY : pi
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
USE MOD_Particle_Vars,          ONLY: PartState,PDM,LastPartPos
USE MOD_Particle_Analyze_Vars,  ONLY: PartPath,doParticleDispersionTrack,doParticlePathTrack
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


SUBROUTINE ParticleRecord(OutputTime,writeToBinary)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_WriteArray         ,ONLY: WriteArray
USE MOD_IO_HDF5                 ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Output_Vars             ,ONLY: WriteStateFiles
USE MOD_Particle_Vars           ,ONLY: PartState,PDM,LastPartPos,PartSpecies,Species,nSpecies,PartIndex,doPartIndex
USE MOD_Particle_Analyze_Vars   ,ONLY: RPP_MaxBufferSize,RPP_Plane,RecordPart,RPP_nVarNames
USE MOD_HDF5_Output             ,ONLY: WriteAttribute
#if USE_MPI
USE MOD_Particle_HDF5_output    ,ONLY: DistributedWriteArray
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: OutputTime
LOGICAL,OPTIONAL,INTENT(IN)    :: writeToBinary
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPart,iRecord,iSpecies
CHARACTER(LEN=200)             :: FileName_loc
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: locRPP,RPP_glob,offsetRPP
#if USE_MPI
INTEGER                        :: sendbuf(2),recvbuf(2)
INTEGER                        :: nRecords(0:nProcessors-1)
#endif
CHARACTER(LEN=32)              :: tmpStr
REAL                           :: SphericityIC(nSpecies)
#if USE_EXTEND_RHS
INTEGER                        :: ForceIC(6,nSpecies)
#else
INTEGER                        :: ForceIC(1,nSpecies)
#endif
!===================================================================================================================================
IF (.NOT.WriteStateFiles) RETURN

!IF(RPP_Type.EQ.'plane')THEN
DO iRecord = 1,RecordPart
  DO iPart=1,PDM%ParticleVecLength
    IF ((PartState(  RPP_Plane(iRecord)%dir,iPart).GE.RPP_Plane(iRecord)%pos) .AND. &
        (LastPartPos(RPP_Plane(iRecord)%dir,iPart).LT.RPP_Plane(iRecord)%pos)) THEN
      RPP_Plane(iRecord)%RPP_Records = RPP_Plane(iRecord)%RPP_Records+1
      ! Part pos and vel
      RPP_Plane(iRecord)%RPP_Data(1:6,RPP_Plane(iRecord)%RPP_Records) = PartState(1:6,iPart)
      ! dp
      RPP_Plane(iRecord)%RPP_Data(7,RPP_Plane(iRecord)%RPP_Records)   = PartState(PART_DIAM,iPart)
      ! Species
      RPP_Plane(iRecord)%RPP_Data(8,RPP_Plane(iRecord)%RPP_Records)   = PartSpecies(iPart)
      ! Index
      IF(doPartIndex) RPP_Plane(iRecord)%RPP_Data(9,RPP_Plane(iRecord)%RPP_Records)   = PartIndex(iPart)
    END IF
  END DO
!END IF

  IF((RPP_Plane(iRecord)%RPP_Records .GE. RPP_MaxBufferSize .OR. PRESENT(writeToBinary)))THEN

    locRPP    = RPP_Plane(iRecord)%RPP_Records
    RPP_glob  = 0

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
    RPP_glob    = sendbuf(1)
    CALL MPI_GATHER(locRPP,1,MPI_INTEGER,nRecords,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
#else
    offsetRPP  = 0
    RPP_glob   = locRPP
#endif

    IF(RPP_glob.EQ.0) RETURN

    ALLOCATE(StrVarNames(RPP_nVarNames))
    StrVarNames(1) ='PartPosX'
    StrVarNames(2) ='PartPosY'
    StrVarNames(3) ='PartPosZ'
    StrVarNames(4) ='VelocityX'
    StrVarNames(5) ='VelocityY'
    StrVarNames(6) ='VelocityZ'
    StrVarNames(7) ='PartDiam'
    StrVarNames(8) ='Species'
    IF(doPartIndex) StrVarNames(9) ='Index'

    ForceIC = 0

    WRITE(UNIT=tmpStr,FMT='(I0)') iRecord
    FileName_loc = TRIM(TIMESTAMP('recordpoints/recordpoints_part'//TRIM(ADJUSTL(tmpStr)),OutputTime))//'.h5'
    SWRITE(UNIT_stdOut,*)' Opening file '//TRIM(FileName_loc)
    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName_loc,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttribute(File_ID,'VarNamesPart',RPP_nVarNames,StrArray=StrVarNames)
      CALL WriteAttribute(File_ID,'nSpecies',1,IntScalar=nSpecies)
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
      CALL WriteAttribute(File_ID,'SphericityIC',nSpecies,RealArray=SphericityIC)
      CALL WriteAttribute(File_ID,'DragForce'   ,nSpecies,IntArray=ForceIC(1,:))
#if USE_EXTEND_RHS
      CALL WriteAttribute(File_ID,'SaffmanForce',nSpecies,IntArray=ForceIC(2,:))
      CALL WriteAttribute(File_ID,'BassForce'   ,nSpecies,IntArray=ForceIC(3,:))
      CALL WriteAttribute(File_ID,'VirtForce'   ,nSpecies,IntArray=ForceIC(4,:))
      CALL WriteAttribute(File_ID,'UndiForce'   ,nSpecies,IntArray=ForceIC(5,:))
      CALL WriteAttribute(File_ID,'MagnusForce' ,nSpecies,IntArray=ForceIC(6,:))
#endif
      CALL CloseDataFile()
    END IF

#if USE_MPI
    CALL DistributedWriteArray(FileName_loc                                  ,&
                               DataSetName  = 'RecordData'                   ,&
                               rank         = 2                              ,&
                               nValGlobal   = (/RPP_nVarNames ,RPP_glob  /)    ,&
                               nVal         = (/RPP_nVarNames ,locRPP    /)    ,&
                               offset       = (/0           ,offsetRPP/)     ,&
                               collective   = .FALSE.                        ,&
                               offSetDim=2                                   ,&
                               communicator = MPI_COMM_FLEXI                 ,&
                               RealArray    = RPP_Plane(iRecord)%RPP_Data(1:RPP_nVarNames,1:locRPP))
  !CALL MPI_BARRIER(PartMPI%COMM,iERROR)
#else
    CALL OpenDataFile(FileName_loc,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArray(           DataSetName  = 'RecordData'                   ,&
                               rank         = 2                              ,&
                               nValGlobal   = (/RPP_nVarNames ,RPP_glob  /)    ,&
                               nVal         = (/RPP_nVarNames ,locRPP    /)    ,&
                               offset       = (/0           ,offsetRPP/)     ,&
                               collective   = .TRUE.                         ,&
                               RealArray    = RPP_Plane(iRecord)%RPP_Data(1:RPP_nVarNames,1:locRPP))
    CALL CloseDataFile()
#endif /*MPI*/

    RPP_Plane(iRecord)%RPP_Data=0.0

    SWRITE(UNIT_stdOut,'(A)') ' WRITE PARTICLE RECORD PLANE TO HDF5 ... DONE'
    SWRITE(UNIT_StdOut,'(132("-"))')
    DEALLOCATE(StrVarNames)

    RPP_Plane(iRecord)%RPP_Records=0
  END IF
END DO

END SUBROUTINE ParticleRecord

END MODULE MOD_Particle_Analyze_Tools
