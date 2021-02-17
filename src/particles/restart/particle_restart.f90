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
USE MOD_PreProc
USE MOD_ErosionPoints,           ONLY: RestartErosionPoint
USE MOD_ErosionPoints_Vars,      ONLY: doParticleImpactTrack
USE MOD_Eval_XYZ,                ONLY: GetPositionInRefElem
USE MOD_HDF5_Input
USE MOD_HDF5_Output,             ONLY: FlushFiles
USE MOD_Mesh_Vars,               ONLY: offsetElem
USE MOD_Part_Tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Boundary_Analyze,ONLY:CalcSurfaceValues
USE MOD_Particle_Boundary_Vars,  ONLY: ImpactRestart,doParticleReflectionTrack
USE MOD_Particle_Globals
USE MOD_Particle_Localization,   ONLY: LocateParticleInElement
USE MOD_Particle_Mesh_Vars,      ONLY: ElemEpsOneCell
USE MOD_Particle_Mesh_Tools,     ONLY: GetCNElemID
USE MOD_Particle_Restart_Vars
USE MOD_Particle_Tracking_Vars,  ONLY: TrackingMethod
USE MOD_Particle_Vars,           ONLY: PartState,PartSpecies,PEM,PDM,Species,nSpecies,PartPosRef,PartReflCount
USE MOD_Restart_Vars,            ONLY: RestartTime,RestartFile
#if USE_MPI
USE MOD_Mesh_Vars,               ONLY: nGlobalElems
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Vars,       ONLY: PartMPI
#endif /*MPI*/
! Particle turbulence models
USE MOD_Particle_Vars,           ONLY: TurbPartState
USE MOD_Particle_SGS_Vars,       ONLY: nSGSVars
#if USE_RW
USE MOD_Particle_RandomWalk_Vars,ONLY: nRWVars
#endif
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
REAL, ALLOCATABLE        :: SendBuff(:),RecvBuff(:)
INTEGER                  :: LostParts(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                  :: NbrOfFoundParts, CompleteNbrOfFound, RecCount(0:PartMPI%nProcs-1)
#endif /*MPI*/
!REAL                     :: VFR_total
INTEGER                  :: i
INTEGER                  :: iElem,CNElemID
LOGICAL                  :: doFlushFiles_loc
! Particle turbulence models
INTEGER                  :: TurbPartSize         !number of turbulent properties with curent setup
INTEGER                  :: TurbPartDataSize     !number of turbulent properties in HDF5 file
REAL,ALLOCATABLE         :: TurbPartData(:,:)    !number of entries in each line of TurbPartData
!===================================================================================================================================
doFlushFiles_loc = MERGE(doFlushFiles, .TRUE., PRESENT(doFlushFiles))

IF (LEN_TRIM(RestartFile).GT.0) THEN
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' READING PARTICLES FROM RESTARTFILE...'
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

  ! Read the first/last ElemID on the local proc
  FirstElemInd = offsetElem+1
  LastElemInd  = offsetElem+PP_nElems

  ! Check if file contains particle information
  CALL DatasetExists(File_ID,'PartData',PartDataExists)
  IF(PartDataExists)THEN
    ! PartInt contains the indices of the particles in each element
    ALLOCATE(PartInt(PartIntSize,FirstElemInd:LastElemInd))

    ! Read number of local particles and their offset from HDF5
    CALL ReadArray('PartInt',2,(/PartIntSize,PP_nElems/),offsetElem,2,IntArray=PartInt)
    locnPart    = PartInt(ELEM_LastPartInd,LastElemInd)-PartInt(ELEM_FirstPartInd,FirstElemInd)
    offsetnPart = PartInt(ELEM_FirstPartInd,FirstElemInd)

    ! Get size of PartData
    CALL GetDataSize(File_ID,'PartData',PartDim,HSize)
    CHECKSAFEINT(HSize(2),4)
    PartDataSize = INT(HSize(1))
    SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of particle variables:           ', PartDataSize

    ! Reflections are stored in the 8th data column. Do not start counting reflections mid-simulation
    IF (PartDataSize.EQ.7) THEN
        doParticleReflectionTrack = .FALSE.
        SWRITE(UNIT_stdOut,'(A3,A30,A3,I33)')' | ','Reflections not tracked previously. Disabling reflection counter',' | ',PartDataSize
    END IF

    ! Get PartData
    ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
    CALL ReadArray('PartData',2,(/PartDataSize,locnPart/),offsetnPart,2,RealArray=PartData)

    ! Loop over all particles on local proc and fill the PartState
    IF (locnPart.GT.0) THEN
      PartState(1:6,1:locnPart)     =     PartData(1:6,offsetnPart+1:offsetnPart+locnPart)
      PartSpecies(  1:locnPart)     = INT(PartData(7  ,offsetnPart+1:offsetnPart+locnPart))
      ! Reflections were tracked previously and are therefore still enabled
      IF (PartDataSize.EQ.8) THEN
          PartReflCount(1:locnPart) = INT(PartData(8,offsetnPart+1:offsetnPart+locnPart))
      END IF

      ! Fill the particle-to-element-mapping (PEM) with the information from HDF5
      DO iElem = FirstElemInd,LastElemInd
        ! If the element contains a particle, add its ID to the element. Particles are still ordered along the SFC at this point
        IF (PartInt(ELEM_LastPartInd,iElem).GT.PartInt(ELEM_FirstPartInd,iElem)) THEN
          PEM%Element    (PartInt(ELEM_FirstPartInd,iElem)-offsetnPart+1 : &
                                 PartInt(ELEM_LastPartInd,iElem) -offsetnPart)  = iElem
          PEM%LastElement(PartInt(ELEM_FirstPartInd,iElem)-offsetnPart+1 : &
                                 PartInt(ELEM_LastPartInd,iElem) -offsetnPart)  = iElem
        END IF
      END DO

      ! All particle properties restored. Particle ready to be tracked
      PDM%ParticleInside(1:locnPart) = .TRUE.
    END IF

    ! Total size of turbulent properties
    TurbPartDataSize = nSGSVars
#if USE_RW
    TurbPartDataSize = TurbPartDataSize + nRWVars
#endif

    ! Check if TurbPartData exists in HDF5 file
    CALL DatasetExists(File_ID,'TurbPartData',TurbPartDataExists)
    IF(TurbPartDataExists)THEN

      ! Get size of TurbPartData, reuse PartDim dummy
      CALL GetDataSize(File_ID,'TurbPartData',PartDim,HSize)
      CHECKSAFEINT(HSize(2),4)
      TurbPartSize = INT(HSize(1))
      SWRITE(UNIT_stdOut,'(A3,A38,A3,I25)')' | ','Number of turbulent particle variables',' | ',TurbPartSize

      ! Compare number of turbulent properties of current run against HDF5
      IF ((TurbPartDataSize.EQ.0).AND.(TurbPartSize.NE.0)) THEN
        SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES') ' | HDF5 state file containing SGS/RW data but current run in DNS mode. Ignoring...'
      ELSEIF ((TurbPartDataSize.NE.0).AND.(TurbPartDataSize.NE.TurbPartSize)) THEN
        CALL abort(__STAMP__,' Number of turbulent variables in HDF5 does not match requested SGS/RW model!')
      ELSEIF ((TurbPartDataSize.NE.0).AND.(TurbPartDataSize.EQ.TurbPartSize)) THEN
        ! Maybe add check that compares the models here. Should resort to a warning in case we want to change the model midrun, so
        ! nothing urgent
        SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='YES') ' | HDF5 state file containing SGS/RW data with ',TurbPartDataSize, &
                                                    ' variables and matching current setup. Continuing run...'
        ! Get TurbPartData
        ALLOCATE(TurbPartData(TurbPartDataSize,offsetnPart+1:offsetnPart+locnPart))
        CALL ReadArray('TurbPartData',2,(/TurbPartDataSize,locnPart/),offsetnPart,2,RealArray=TurbPartData)

        ! Loop over all particles on local proc and fill the TurbPartState
        IF (locnPart.GT.0) THEN
          TurbPartState(1:TurbPartDataSize,1:locnPart) = TurbPartData(1:TurbPartDataSize,offsetnPart+1:offsetnPart+locnPart)
        END IF

        ! De-allocate array used for readin
        DEALLOCATE(TurbPartData)

      ! Last case is no turbulent properties in current run or HDF5. Do nothing ...
      END IF
    ! No turbulent particle properties in HDF5
    ELSE
      IF (TurbPartDataSize.NE.0) THEN
        SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES') ' | HDF5 state file not containing SGS/RW data but model active in current run.'
        SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES') ' | SGS/RW model will start without history...'
      END IF
    END IF

    ! De-allocate arrays used for read-in
    DEALLOCATE(PartInt,PartData)
    PDM%ParticleVecLength = PDM%ParticleVecLength + locnPart
    CALL UpdateNextFreePosition()

    SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' READING PARTICLES FROM RESTARTFILE DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')

    ! Reconstruct the number of particles inserted before restart from the emission rate
    DO i=1,nSpecies
      DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
        Species(i)%Init(iInit)%InsertedParticle = INT(Species(i)%Init(iInit)%ParticleEmission * RestartTime,8)
      END DO
    END DO

    ! Abort if HDF5 file contains more particles than memory allocated for in current run
    IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
      CALL abort(__STAMP__,' Number of Particles in Restart higher than MaxParticleNumber!')
    END IF

    ! Since the elementside-local node number are NOT persistant and dependent on the location
    ! of the MPI borders, all particle-element mappings need to be checked after a restart
    ! Step 1: Identify particles that are not in the element in which they were before the restart
    COUNTER  = 0
    COUNTER2 = 0

    DO i = 1,PDM%ParticleVecLength
      CALL GetPositionInRefElem(PartState(1:3,i),Xi,PEM%Element(i))
      ! Particle already inside the correct element
      CNElemID = GetCNElemID(PEM%Element(i))

      SELECT CASE(TrackingMethod)
        ! CASE(REFMAPPING)
        CASE(TRACING,REFMAPPING)
          IF(ALL(ABS(Xi).LE.ElemEpsOneCell(CNElemID))) THEN
            InElementCheck=.TRUE.
            PartPosRef(1:3,i)=Xi
          ELSE
            InElementCheck=.FALSE.
          END IF

        ! FLEXI has ElemEpsOneCell also with TRACING
        ! CASE(TRACING,TRIATRACKING)
        CASE(TRIATRACKING)
          IF(ALL(ABS(Xi).LE.1.0)) THEN ! particle inside
            InElementCheck=.TRUE.
            IF(ALLOCATED(PartPosRef)) PartPosRef(1:3,i)=Xi
          ELSE
            InElementCheck=.FALSE.
          END IF
      END SELECT

      ! Particle not inside the correct element, try to find them within the current proc
      IF (.NOT.InElementCheck) THEN
        COUNTER = COUNTER + 1
        CALL LocateParticleInElement(i,doHALO=.TRUE.)
        ! Particle located on the current node
        IF (PDM%ParticleInside(i)) THEN
          ! particle on the current node but on the wrong proc
          IF ((PEM%Element(i).GT.offsetElem).AND.(PEM%Element(i).LE.offsetElem+PP_nElems)) THEN
            COUNTER2 = COUNTER2 + 1
            PartPosRef(1:3,i) = -888.
          ! Particle on the correct proc but inside the wrong element
          ELSE
            PEM%LastElement(i) = PEM%Element(i)
          END IF
        ! Particle located on the current node
        ELSE
          ! Every element was checked, the particle is truly lost
#if USE_MPI
          IF (nComputeNodeTotalElems.EQ.nGlobalElems) THEN
#endif /*USE_MPI*/
            CALL ABORT(__STAMP__,'Particle not located on the compute-node during restart. PartID:',i)
          ! Still elements unchecked, let the procs on the other nodes figure it out
#if USE_MPI
          ELSE
            COUNTER2 = COUNTER2 + 1
            PartPosRef(1:3,i) = -888.
          END IF
#endif /*USE_MPI*/
        END IF
      END IF
    END DO

#if USE_MPI
    ! Step 2: All particles that are not found withing MyProc need to be communicated to the others and located there
    ! Combine number of lost particles of all processes and allocate variables
    CALL MPI_ALLGATHER(COUNTER2, 1, MPI_INTEGER, LostParts, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
    IF (SUM(LostParts).GT.0) THEN
      ALLOCATE(SendBuff(1:COUNTER2*PartDataSize))
      ALLOCATE(RecvBuff(1:SUM(LostParts)*PartDataSize))

      ! Fill SendBuffer
      COUNTER = 0
      DO i = 1, PDM%ParticleVecLength
        IF (.NOT.PDM%ParticleInside(i)) THEN
          SendBuff(COUNTER+1:COUNTER+6) = PartState(1:6,i)
          SendBuff(COUNTER+7)           = REAL(PartSpecies(i))
          ! Reflections were tracked previously and are therefore still enabled
          IF (PartDataSize.EQ.8) THEN
          SendBuff(COUNTER+8)           = REAL(PartReflCount(i))
          END IF
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
           RecvBuff, RecCount, Displace, MPI_DOUBLE_PRECISION, PartMPI%COMM, IERROR)

      ! Add them to particle list and check if they are in MyProcs domain
      NbrOfFoundParts = 0
      CurrentPartNum  = PDM%ParticleVecLength+1
      COUNTER         = 0
      DO i = 1, SUM(LostParts)
        PartState(1:6,CurrentPartNum) = RecvBuff(COUNTER+1:COUNTER+6)
        PDM%ParticleInside(CurrentPartNum) = .TRUE.
        CALL LocateParticleInElement(i,doHALO=.FALSE.)
        ! Particle is located on current node. Only keep it if it belongs on this proc
        IF (PDM%ParticleInside(CurrentPartNum).AND.(PEM%Element(i).GT.offsetElem).AND.(PEM%Element(i).LE.offsetElem+PP_nElems)) THEN
          PEM%LastElement(CurrentPartNum) = PEM%Element(CurrentPartNum)
          PartSpecies(    CurrentPartNum) = INT(RecvBuff(COUNTER+7))
          ! Reflections were tracked previously and are therefore still enabled
          IF (PartDataSize.EQ.8) THEN
            PartReflCount(CurrentPartNum) = INT(RecvBuff(COUNTER+8))
          END IF
          NbrOfFoundParts = NbrOfFoundParts + 1
          CurrentPartNum  = CurrentPartNum + 1
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
    IF (COUNTER .NE.0) WRITE(*,*) COUNTER ,'Particles are in different element after restart!'
    IF (COUNTER2.NE.0) WRITE(*,*) COUNTER2,'of which could not be found and are removed!'
#endif
    CALL UpdateNextFreePosition()
  ELSE
      SWRITE(UNIT_stdOut,*)'PartData does not exists in restart file'
      SWRITE(UNIT_stdOut,'(132("-"))')
  END IF ! PartDataExists
CALL CloseDataFile()

  ! Get individual impact data
  IF (doParticleImpactTrack) CALL RestartErosionPoint

  ! Delete all files that will be rewritten --> moved from restart.f90 since we need it here
  IF (doFlushFiles_loc) CALL FlushFiles(RestartTime)
! No restart
ELSE
  ! Delete all files since we are doing a fresh start --> moved from restart.f90 since we need it here
  IF (doFlushFiles_loc) CALL FlushFiles()
END IF

! Make sure we update our surface data after reading the sampling data
#if USE_MPI
! Only procs with SurfMesh know about the restart so far. Communicate to all here
CALL MPI_ALLREDUCE(MPI_IN_PLACE,ImpactRestart,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError)
#endif

IF (ImpactRestart) CALL CalcSurfaceValues(restart_opt=.TRUE.)

!#if USE_MPI
!CALL MPI_BARRIER(PartMPI%COMM,iError)
!#endif /*MPI*/

END SUBROUTINE ParticleRestart

END MODULE MOD_Particle_Restart
