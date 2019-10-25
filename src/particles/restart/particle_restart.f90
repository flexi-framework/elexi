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

!==================================================================================================================================
!> \brief Routines that handle restart capabilities.
!>
!> With this feature a simulation can be resumed from a state file that has been created during a previous
!> simulation (restart file). The restart file is passed to FLEXI as a second command line argument.
!> The restart can also be performed from a file with a different polynomial degree or node type than the current simulation.
!==================================================================================================================================
MODULE MOD_Particle_Restart
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ParticleRestart
  MODULE PROCEDURE ParticleRestart
END INTERFACE

PUBLIC :: ParticleRestart
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE ParticleRestart(doFlushFiles)
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Restart_Vars
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_HDF5_Output,             ONLY:FlushFiles
USE MOD_Mesh_Vars,               ONLY:offsetElem
USE MOD_Restart_Vars,            ONLY:RestartTime,RestartFile
USE MOD_Particle_Vars,           ONLY:PartState, PartSpecies, PEM, PDM, Species, nSpecies, PartPosRef,PartReflCount
USE MOD_Part_Tools,              ONLY:UpdateNextFreePosition
USE MOD_Eval_XYZ,                ONLY:TensorProductInterpolation
USE MOD_Particle_Localization,   ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Mesh_Vars,      ONLY:epsOneCell
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
USE MOD_Particle_Erosion_Analyze,ONLY:CalcSurfaceValues
USE MOD_Particle_Erosion_Vars,   ONLY:ErosionRestart,PartTrackReflection
USE MOD_ErosionPoints,           ONLY:RestartErosionPoint
USE MOD_ErosionPoints_Vars,      ONLY:EP_inUse
#if USE_MPI
USE MOD_Particle_MPI_Vars,       ONLY:PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: doFlushFiles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: FirstElemInd,LastelemInd,iInit
INTEGER,ALLOCATABLE      :: PartInt(:,:)
INTEGER,PARAMETER        :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                  :: PartDim              !dummy for rank of partData
INTEGER                  :: PartDataSize         !number of entries in each line of PartData
INTEGER                  :: locnPart,offsetnPart
INTEGER,PARAMETER        :: ELEM_FirstPartInd=1
INTEGER,PARAMETER        :: ELEM_LastPartInd=2
REAL,ALLOCATABLE         :: PartData(:,:)
REAL                     :: xi(3)
LOGICAL                  :: InElementCheck
INTEGER                  :: COUNTER, COUNTER2
#if USE_MPI
REAL, ALLOCATABLE        :: SendBuff(:), RecBuff(:)
INTEGER                  :: LostParts(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                  :: NbrOfFoundParts, CompleteNbrOfFound, RecCount(0:PartMPI%nProcs-1)
#endif /*MPI*/
!REAL                     :: VFR_total
INTEGER                  :: i
INTEGER                  :: iElem
LOGICAL                  :: doFlushFiles_loc
!===================================================================================================================================
doFlushFiles_loc = MERGE(doFlushFiles, .TRUE., PRESENT(doFlushFiles))

IF (LEN_TRIM(RestartFile).GT.0) THEN
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Reading Particles from Restartfile...'

  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  !read local ElemInfo from HDF5
  FirstElemInd=offsetElem+1
  LastElemInd=offsetElem+PP_nElems
  ! read local ParticleInfo from HDF5
  CALL DatasetExists(File_ID,'PartData',PartDataExists)
  IF(PartDataExists)THEN
    ! get PartInt
    ALLOCATE(PartInt(PartIntSize,FirstElemInd:LastElemInd))
    CALL ReadArray('PartInt',2,(/PartIntSize,PP_nElems/),offsetElem,2,IntArray=PartInt)
    ! read local Particle Data from HDF5
    locnPart=PartInt(ELEM_LastPartInd,LastElemInd)-PartInt(ELEM_FirstPartInd,FirstElemInd)
    offsetnPart=PartInt(ELEM_FirstPartInd,FirstElemInd)
    ! get size of PartData
    CALL GetDataSize(File_ID,'PartData',PartDim,HSize)
    CHECKSAFEINT(HSize(2),4)
    PartDataSize = INT(HSize(1))
    SWRITE(UNIT_stdOut,'(A3,A30,A3,I33)')' | ','Number of particle variables',' | ',PartDataSize
    ! Do not start counting reflections mid-simulation
    IF (PartDataSize.EQ.7) THEN
        PartTrackReflection = .FALSE.
        SWRITE(UNIT_stdOut,'(A3,A30,A3,I33)')' | ','Reflections not tracked previously. Disabling reflection counter',' | ',PartDataSize
    END IF
    ! get PartData
    ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
    CALL ReadArray('PartData',2,(/PartDataSize,locnPart/),offsetnPart,2,RealArray=PartData)!,&
                           !xfer_mode_independent=.TRUE.)
    IF (locnPart.GT.0) THEN
      PartState(1,1:locnPart)   = PartData(1,offsetnPart+1:offsetnPart+locnPart)
      PartState(2,1:locnPart)   = PartData(2,offsetnPart+1:offsetnPart+locnPart)
      PartState(3,1:locnPart)   = PartData(3,offsetnPart+1:offsetnPart+locnPart)
      PartState(4,1:locnPart)   = PartData(4,offsetnPart+1:offsetnPart+locnPart)
      PartState(5,1:locnPart)   = PartData(5,offsetnPart+1:offsetnPart+locnPart)
      PartState(6,1:locnPart)   = PartData(6,offsetnPart+1:offsetnPart+locnPart)
      PartSpecies(1:locnPart)   = INT(PartData(7,offsetnPart+1:offsetnPart+locnPart))
      IF (PartDataSize.EQ.8) THEN
          PartReflCount(1:locnPart) = INT(PartData(8,offsetnPart+1:offsetnPart+locnPart))
      END IF
      DO iElem=FirstElemInd,LastElemInd
        IF (PartInt(ELEM_LastPartInd,iElem).GT.PartInt(ELEM_FirstPartInd,iElem)) THEN
          PEM%Element(PartInt(ELEM_FirstPartInd,iElem)-offsetnPart+1 : &
                                 PartInt(ELEM_LastPartInd,iElem) -offsetnPart)  = iElem-offsetElem
          PEM%LastElement(PartInt(ELEM_FirstPartInd,iElem)-offsetnPart+1 : &
                                 PartInt(ELEM_LastPartInd,iElem) -offsetnPart)  = iElem-offsetElem
        END IF
      END DO
      PDM%ParticleInside(1:locnPart) = .TRUE.
    END IF
    DEALLOCATE(PartInt,PartData)
    PDM%ParticleVecLength = PDM%ParticleVecLength + locnPart
    CALL UpdateNextFreePosition()

    SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Reading Particles from Restartfile... DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')

    DO i=1,nSpecies
      DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
        Species(i)%Init(iInit)%InsertedParticle = INT(Species(i)%Init(iInit)%ParticleEmission * RestartTime,8)
      END DO
    END DO

    ! if ParticleVecLength GT maxParticleNumber: Stop
    IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
      CALL abort(&
  __STAMP__&
  ,' Number of Particles in Restart higher than MaxParticleNumber!')
    END IF

    ! Since the elementside-local node number are NOT persistant and dependent on the location
    ! of the MPI borders, all particle-element mappings need to be checked after a restart
    ! Step 1: Identify particles that are not in the element in which they were before the restart
    COUNTER = 0
    COUNTER2 = 0
    IF(DoRefMapping) THEN
      DO i = 1,PDM%ParticleVecLength
        CALL TensorProductInterpolation(PartState(1:3,i),Xi,PEM%Element(i))
        IF(ALL(ABS(Xi).LE.EpsOneCell(PEM%Element(i)))) THEN ! particle inside
          InElementCheck=.TRUE.
          PartPosRef(1:3,i)=Xi
        ELSE
          InElementCheck=.FALSE.
        END IF
        IF (.NOT.InElementCheck) THEN  ! try to find them within MyProc
          COUNTER = COUNTER + 1
          !CALL SingleParticleToExactElement(i)
          CALL SingleParticleToExactElement(i,doHALO=.FALSE.,initFix=.FALSE.,doRelocate=.FALSE.)
          IF (.NOT.PDM%ParticleInside(i)) THEN
            COUNTER2 = COUNTER2 + 1
            PartPosRef(1:3,i) = -888.
          ELSE
            PEM%LastElement(i) = PEM%Element(i)
          END IF
        END IF
      END DO
    ELSE ! no Ref Mapping
      DO i = 1,PDM%ParticleVecLength
        CALL TensorProductInterpolation(PartState(1:3,i),Xi,PEM%Element(i))
        IF(ALL(ABS(Xi).LE.1.0)) THEN ! particle inside
          InElementCheck=.TRUE.
          IF(ALLOCATED(PartPosRef)) PartPosRef(1:3,i)=Xi
        ELSE
          InElementCheck=.FALSE.
        END IF
        IF (.NOT.InElementCheck) THEN  ! try to find them within MyProc
          COUNTER = COUNTER + 1
          !CALL SingleParticleToExactElement(i)
          CALL SingleParticleToExactElementNoMap(i,doHALO=.FALSE.,doRelocate=.FALSE.)
          IF (.NOT.PDM%ParticleInside(i)) THEN
            COUNTER2 = COUNTER2 + 1
          ELSE
            PEM%LastElement(i) = PEM%Element(i)
          END IF
        END IF
      END DO
    END IF

#if USE_MPI
    ! Step 2: All particles that are not found withing MyProc need to be communicated to the others and located there
    ! Combine number of lost particles of all processes and allocate variables
    CALL MPI_ALLGATHER(COUNTER2, 1, MPI_INTEGER, LostParts, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
    IF (SUM(LostParts).GT.0) THEN
      ALLOCATE(SendBuff(1:COUNTER2*PartDataSize))
      ALLOCATE(RecBuff(1:SUM(LostParts)*PartDataSize))
      ! Fill SendBuffer
      COUNTER = 0
      DO i = 1, PDM%ParticleVecLength
        IF (.NOT.PDM%ParticleInside(i)) THEN
          SendBuff(COUNTER+1:COUNTER+6) = PartState(1:6,i)
          SendBuff(COUNTER+7)           = REAL(PartSpecies(i))
          COUNTER = COUNTER + PartDataSize
        END IF
      END DO
      ! Distribute lost particles to all procs
      COUNTER = 0
      DO i = 0, PartMPI%nProcs-1
        RecCount(i) = LostParts(i) * PartDataSize
        Displace(i) = COUNTER
        COUNTER = COUNTER + LostParts(i)*PartDataSize
      END DO
      CALL MPI_ALLGATHERV(SendBuff, PartDataSize*LostParts(PartMPI%MyRank), MPI_DOUBLE_PRECISION, &
           RecBuff, RecCount, Displace, MPI_DOUBLE_PRECISION, PartMPI%COMM, IERROR)
      ! Add them to particle list and check if they are in MyProcs domain
      NbrOfFoundParts = 0
      CurrentPartNum = PDM%ParticleVecLength+1
      COUNTER = 0
      DO i = 1, SUM(LostParts)
        PartState(1:6,CurrentPartNum) = RecBuff(COUNTER+1:COUNTER+6)
        PDM%ParticleInside(CurrentPartNum) = .true.
        IF(DoRefMapping)THEN
          CALL SingleParticleToExactElement(CurrentPartNum,doHALO=.FALSE.,initFix=.FALSE.,doRelocate=.FALSE.)
        ELSE
          CALL SingleParticleToExactElementNoMap(CurrentPartNum,doHALO=.FALSE.,doRelocate=.FALSE.)
        END IF
        !CALL SingleParticleToExactElement(CurrentPartNum)
        IF (PDM%ParticleInside(CurrentPartNum)) THEN
          PEM%LastElement(CurrentPartNum) = PEM%Element(CurrentPartNum)
          NbrOfFoundParts = NbrOfFoundParts + 1
          PartSpecies(CurrentPartNum) = INT(RecBuff(COUNTER+7))
          CurrentPartNum = CurrentPartNum + 1
        END IF
        COUNTER = COUNTER + PartDataSize
      END DO
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfFoundParts
      ! Combine number of found particles to make sure none are lost completely
      CALL MPI_ALLREDUCE(NbrOfFoundParts, CompleteNbrOfFound, 1, MPI_INTEGER, MPI_SUM, PartMPI%COMM, IERROR)
      SWRITE(UNIT_stdOut,*) SUM(LostParts),'were not in the correct proc after restart.'
      SWRITE(UNIT_stdOut,*) CompleteNbrOfFound,'of these were found in other procs.'
      SWRITE(UNIT_stdOut,*) SUM(LostParts)-CompleteNbrOfFound,'were not found and have been removed.'
    END IF
#else
    IF (COUNTER.NE.0) WRITE(*,*) COUNTER,'Particles are in different element after restart!'
    IF (COUNTER2.NE.0) WRITE(*,*) COUNTER2,'of which could not be found and are removed!'
#endif
    CALL UpdateNextFreePosition()
  ELSE
      SWRITE(UNIT_stdOut,*)'PartData does not exists in restart file'
      SWRITE(UNIT_stdOut,'(132("-"))')
  END IF ! PartDataExists
CALL CloseDataFile()

  ! Get individual impact data
IF (EP_inUse)           CALL RestartErosionPoint

  ! Delete all files that will be rewritten --> moved from restart.f90 since we need it here
  IF (doFlushFiles_loc) CALL FlushFiles(RestartTime)
ELSE
  ! Delete all files since we are doing a fresh start --> moved from restart.f90 since we need it here
  IF (doFlushFiles_loc) CALL FlushFiles()
END IF

! Make sure we update our surface data after reading the sampling data
#if USE_MPI
! Only procs with SurfMesh know about the restart so far. Communicate to all here
CALL MPI_ALLREDUCE(MPI_IN_PLACE,ErosionRestart,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError)
#endif

IF (ErosionRestart) CALL CalcSurfaceValues(restart_opt=.TRUE.)

!#if USE_MPI
!CALL MPI_BARRIER(PartMPI%COMM,iError)
!#endif /*MPI*/

END SUBROUTINE ParticleRestart

END MODULE MOD_Particle_Restart
