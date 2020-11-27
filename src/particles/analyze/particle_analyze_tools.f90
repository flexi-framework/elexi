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
USE MOD_Particle_Vars,           ONLY: PartState,PDM,LastPartPos,PartSpecies,nSpecies
USE MOD_Particle_Analyze_Vars,   ONLY: RPP_MaxBufferSize,RPP_Plane,RPP_Type
USE MOD_io_bin,                  ONLY: prt_bin!,load_bin
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: OutputTime
LOGICAL,OPTIONAL,INTENT(IN) :: writeToBinary
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart,iElem,ElemID,RP_glob
!REAL,ALLOCATABLE            :: buffOuttmp(:,:)
CHARACTER(LEN=200)          :: FileName_loc, FileString_loc        ! FileName with data type extension
#if USE_MPI
INTEGER,ALLOCATABLE         :: displacement(:), nRPRecords_per_proc(:)
INTEGER                     :: i
REAL,ALLOCATABLE            :: RPP_Data_glob(:,:)
#endif
!===================================================================================================================================

IF(RPP_Type.EQ.'plane')THEN
  DO iPart=1,PDM%ParticleVecLength
    IF ((PartState(1,iPart).GE.RPP_Plane%x(1,1)) .AND. (LastPartPos(1,iPart).LT.RPP_Plane%x(1,1)))THEN
      RPP_Plane%RPP_Records=RPP_Plane%RPP_Records+1
      RPP_Plane%RPP_Data(1:6,RPP_Plane%RPP_Records) = PartState(1:6,iPart)
      RPP_Plane%RPP_Data(7,RPP_Plane%RPP_Records)   = PartSpecies(iPart)
    END IF
  END DO
END IF

IF((RPP_Plane%RPP_Records .GE. RPP_MaxBufferSize .OR. PRESENT(writeToBinary)))THEN
  !>> Sum up particles from the other procs
#if USE_MPI
  IF(MPIRoot)THEN
    ALLOCATE(nRPRecords_per_proc(0:nProcessors-1))
  ELSE
    ALLOCATE(nRPRecords_per_proc(1))
  END IF
  CALL MPI_GATHER(RPP_Plane%RPP_Records, 1, MPI_INTEGER, nRPRecords_per_proc,&
    1, MPI_INTEGER, 0, MPI_COMM_FLEXI, iERROR)

  !>> Gather data to root
  IF (MPIRoot) THEN
    ALLOCATE(displacement(1:nProcessors))
    displacement(1)=0
    DO i=2,nProcessors
      displacement(i)=SUM(7*nRPRecords_per_proc(0:i-2))
    END DO
    ALLOCATE(RPP_Data_glob(1:7,1:SUM(nRPRecords_per_proc)))
  ELSE
    ALLOCATE(displacement(1))
    ALLOCATE(RPP_Data_glob(1,1))
  END IF

  CALL MPI_GATHERV(RPP_Plane%RPP_Data(:,1:RPP_Plane%RPP_Records),RPP_Plane%RPP_Records*7,MPI_DOUBLE_PRECISION,RPP_Data_glob,nRPRecords_per_proc*7,displacement,&
    MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
#endif

#if USE_MPI
  IF(MPIRoot .AND. SUM(nRPRecords_per_proc).GT.0)THEN
#endif
    ! Write to binary file
    FileName_loc = TRIM(TIMESTAMP('recordpoints_part_',OutputTime))
    FileString_loc=TRIM(FileName_loc)//'.dat'
    SWRITE(UNIT_stdOut,*)' Opening file '//TRIM(FileString_loc)
    SWRITE(UNIT_stdOut,*)' Number of particles ',SUM(nRPRecords_per_proc)

#if USE_MPI
    CALL prt_bin(RPP_Data_glob, FileString_loc)
#else
    CALL prt_bin(RPP_Plane%RPP_Data, FileString_loc)
#endif
!    CALL load_bin(FileName_loc,buffOuttmp)
!    print *, buffOuttmp

    RPP_Plane%RPP_Data=0.0

    SWRITE(UNIT_stdOut,'(A)') ' WRITE STATE TO BINARY... DONE'
    SWRITE(UNIT_StdOut,'(132("-"))')
    DEALLOCATE(RPP_Data_glob)
    DEALLOCATE(displacement,nRPRecords_per_proc)
#if USE_MPI
  END IF
#endif
  RPP_Plane%RPP_Records=0
END IF

END SUBROUTINE ParticleRecord

END MODULE MOD_Particle_Analyze_Tools
