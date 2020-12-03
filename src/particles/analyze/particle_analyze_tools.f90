!=================================================================================================================================
! Copyright (c) 2010-2019  Prof. Claus-Dieter Munz
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

INTERFACE ParticleRecord
  MODULE PROCEDURE ParticleRecord
END INTERFACE

PUBLIC :: CalcEkinPart
PUBLIC :: ParticleRecord
!===================================================================================================================================

CONTAINS

PURE FUNCTION CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
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
REAL                               :: partV2
!===================================================================================================================================

partV2 = PartState(4,iPart) * PartState(4,iPart) &
       + PartState(5,iPart) * PartState(5,iPart) &
       + PartState(6,iPart) * PartState(6,iPart)

CalcEkinPart = 0.5*Species(PartSpecies(iPart))%MassIC*partV2

END FUNCTION CalcEkinPart


SUBROUTINE ParticleRecord(OutputTime,writeToBinary)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,           ONLY: PartState,PDM,LastPartPos,PartSpecies,nSpecies,PartIndex
USE MOD_Particle_Analyze_Vars,   ONLY: RPP_MaxBufferSize,RPP_Plane,RPP_Type
USE MOD_io_bin,                  ONLY: prt_bin!,load_bin
USE MOD_HDF5_Output             ,ONLY: WriteAttribute
USE MOD_IO_HDF5                 ,ONLY: File_ID,OpenDataFile,CloseDataFile
#if USE_MPI
USE MOD_Particle_HDF5_output    ,ONLY: DistributedWriteArray
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: OutputTime
LOGICAL,OPTIONAL,INTENT(IN) :: writeToBinary
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart,RP_glob
!REAL,ALLOCATABLE            :: buffOuttmp(:,:)
CHARACTER(LEN=200)          :: FileName_loc     ! FileName with data type extension
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER, DIMENSION(2)       :: M_shape
INTEGER, PARAMETER          :: type_label = 8
INTEGER, PARAMETER          :: RPPDataSize = 8
INTEGER                     :: locRPP,RPP_glob,offsetRPP
#if USE_MPI
INTEGER                     :: sendbuf(2),recvbuf(2)
INTEGER                     :: nRecords(0:nProcessors-1)
#endif
!===================================================================================================================================

IF(RPP_Type.EQ.'plane')THEN
  DO iPart=1,PDM%ParticleVecLength
    IF ((PartState(1,iPart).GE.RPP_Plane%x(1,1)) .AND. (LastPartPos(1,iPart).LT.RPP_Plane%x(1,1)))THEN
      RPP_Plane%RPP_Records=RPP_Plane%RPP_Records+1
      RPP_Plane%RPP_Data(1:6,RPP_Plane%RPP_Records) = PartState(1:6,iPart)
      RPP_Plane%RPP_Data(7,RPP_Plane%RPP_Records)   = PartSpecies(iPart)
      RPP_Plane%RPP_Data(8,RPP_Plane%RPP_Records)   = PartIndex(iPart)
    END IF
  END DO
END IF

IF((RPP_Plane%RPP_Records .GE. RPP_MaxBufferSize .OR. PRESENT(writeToBinary)))THEN

  locRPP    = RPP_Plane%RPP_Records
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

  ALLOCATE(StrVarNames(RPPDataSize))
  StrVarNames(1) ='ParticlePositionX'
  StrVarNames(2) ='ParticlePositionY'
  StrVarNames(3) ='ParticlePositionZ'
  StrVarNames(4) ='VelocityX'
  StrVarNames(5) ='VelocityY'
  StrVarNames(6) ='VelocityZ'
  StrVarNames(7) ='Species'
  StrVarNames(8) ='Index'

  FileName_loc = TRIM(TIMESTAMP('recordpoints_part',OutputTime))//'.h5'
  SWRITE(UNIT_stdOut,*)' Opening file '//TRIM(FileName_loc)
  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName_loc,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'VarNamesPart',RPPDataSize,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

#if USE_MPI
  CALL DistributedWriteArray(FileName_loc                                  ,&
                             DataSetName  = 'RecordData'                   ,&
                             rank         = 2                              ,&
                             nValGlobal   = (/RPPDataSize ,RPP_glob  /)    ,&
                             nVal         = (/RPPDataSize ,locRPP    /)    ,&
                             offset       = (/0           ,offsetRPP/)     ,&
                             collective   = .FALSE.                        ,&
                             offSetDim=2                                   ,&
                             communicator = MPI_COMM_FLEXI                 ,&
                             RealArray    = RPP_Plane%RPP_Data(1:RPPDataSize,1:locRPP))
!CALL MPI_BARRIER(PartMPI%COMM,iERROR)
#else
  CALL OpenDataFile(FileName_loc,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArray(           DataSetName  = 'RecordData'                   ,&
                             rank         = 2                              ,&
                             nValGlobal   = (/RPPDataSize ,RPP_glob  /)    ,&
                             nVal         = (/RPPDataSize ,locRPP    /)    ,&
                             offset       = (/0           ,offsetRPP/)     ,&
                             collective   = .TRUE.                         ,&
                             RealArray    = RPP_Plane%RPP_Data(1:RPPDataSize,1:locRPP))
  CALL CloseDataFile()
#endif /*MPI*/

  RPP_Plane%RPP_Data=0.0

  SWRITE(UNIT_stdOut,'(A)') ' WRITE STATE TO BINARY... DONE'
  SWRITE(UNIT_StdOut,'(132("-"))')
  DEALLOCATE(StrVarNames)

  RPP_Plane%RPP_Records=0
END IF

END SUBROUTINE ParticleRecord

END MODULE MOD_Particle_Analyze_Tools
