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
USE MOD_Eval_XYZ,                   ONLY: GetPositionInRefElem
USE MOD_HDF5_Input
USE MOD_HDF5_Output,                ONLY: FlushFiles
USE MOD_Mesh_Vars,                  ONLY: offsetElem
USE MOD_Part_Tools,                 ONLY: UpdateNextFreePosition
USE MOD_Particle_Boundary_Analyze,  ONLY: CalcSurfaceValues
USE MOD_Particle_Boundary_Tracking, ONLY: RestartParticleBoundaryTracking
USE MOD_Particle_Boundary_Vars,     ONLY: ImpactRestart,doParticleReflectionTrack
USE MOD_Particle_Boundary_Vars,     ONLY: doParticleImpactTrack
USE MOD_Particle_Globals
USE MOD_Particle_Interpolation_Vars,ONLY: DoInterpolation
USE MOD_Particle_Localization,      ONLY: LocateParticleInElement,ParticleInsideQuad3D
USE MOD_Particle_Mesh_Vars,         ONLY: ElemEpsOneCell
USE MOD_Particle_Mesh_Tools,        ONLY: GetCNElemID
USE MOD_Particle_Restart_Vars
USE MOD_Particle_Tracking_Vars,     ONLY: TrackingMethod,NbrOfLostParticles,CountNbrOfLostParts
USE MOD_Particle_Tracking_Vars,     ONLY: NbrOfLostParticlesTotal,TotalNbrOfMissingParticlesSum,NbrOfLostParticlesTotal_old
USE MOD_Particle_Vars,              ONLY: PartState,PartSpecies,PEM,PDM,Species,nSpecies,PartPosRef,PartReflCount
USE MOD_Restart_Vars,               ONLY: RestartTime,RestartFile
#if USE_MPI
USE MOD_Particle_MPI_Vars,          ONLY: PartMPI
#endif /*MPI*/
! Particle turbulence models
USE MOD_Particle_Vars,              ONLY: TurbPartState
USE MOD_Particle_SGS_Vars,          ONLY: nSGSVars
#if USE_RW
USE MOD_Particle_RandomWalk_Vars,   ONLY: nRWVars
#endif
USE MOD_StringTools,                ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL     :: doFlushFiles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: FirstElemInd,LastelemInd,iInit
INTEGER,ALLOCATABLE             :: PartInt(:,:)
INTEGER,PARAMETER               :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                         :: PartDim              !dummy for rank of partData
INTEGER                         :: PartDataSize         !number of entries in each line of PartData
INTEGER                         :: locnPart,offsetnPart
INTEGER,PARAMETER               :: ELEM_FirstPartInd=1
INTEGER,PARAMETER               :: ELEM_LastPartInd=2
REAL,ALLOCATABLE                :: PartData(:,:)
REAL                            :: xi(3),det(6,2)
LOGICAL                         :: InElementCheck,EmissionTimeExists
INTEGER                         :: NbrOfMissingParticles
#if USE_MPI
INTEGER                         :: iProc
INTEGER,ALLOCATABLE             :: IndexOfFoundParticles(:),CompleteIndexOfFoundParticles(:)
INTEGER                         :: CompleteNbrOfLost,CompleteNbrOfFound,CompleteNbrOfDuplicate
REAL,ALLOCATABLE                :: RecBuff(:,:)
INTEGER                         :: TotalNbrOfMissingParticles(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                         :: OffsetTotalNbrOfMissingParticles(0:PartMPI%nProcs-1)
INTEGER                         :: NbrOfFoundParts, RecCount(0:PartMPI%nProcs-1)
#endif /*MPI*/
!REAL                            :: VFR_total
INTEGER                         :: iSpec,iPart,iStr
INTEGER                         :: iElem,CNElemID
LOGICAL                         :: doFlushFiles_loc
! Particle turbulence models
INTEGER                         :: TurbPartSize         !number of turbulent properties with curent setup
INTEGER                         :: TurbPartDataSize     !number of turbulent properties in HDF5 file
REAL,ALLOCATABLE                :: TurbPartData(:,:)    !number of entries in each line of TurbPartData
INTEGER                         :: PP_nVarPart_loc
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
!===================================================================================================================================
doFlushFiles_loc = MERGE(doFlushFiles, .TRUE., PRESENT(doFlushFiles))

IF (LEN_TRIM(RestartFile).GT.0) THEN
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' READING PARTICLES FROM RESTARTFILE...'
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

  ! Read the emission time
  CALL DatasetExists(File_ID,'EmissionTime',EmissionTimeExists,attrib=.TRUE.)
  ! Reset EmissionTime for ResetTime
  IF (EmissionTimeExists .AND. RestartTime.GT.0) THEN
    CALL ReadAttribute(File_ID,'EmissionTime',1,RealScalar=EmissionTime)
    SWRITE(UNIT_stdOut,'(A,F16.6)')' | Resuming particle emission from time    ',EmissionTime
  ELSE
    EmissionTime = 0.
  END IF

  ! Read the first/last ElemID on the local proc
  FirstElemInd = offsetElem+1
  LastElemInd  = offsetElem+PP_nElems

  ! Check if file contains particle information
  CALL DatasetExists(File_ID,'PartData',PartDataExists)
  IF (PartDataExists) THEN
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

    ! For files, where no particle diameter was saved
    ALLOCATE(StrVarNames(PartDataSize))
    CALL ReadAttribute(File_ID,'VarNamesParticles',PartDataSize,StrArray=StrVarNames)
    ! Species and PartDiam
    PP_nVarPart_loc = PartDataSize-1
    DO iStr=1,PartDataSize
      ! Reflection
      IF (STRICMP('ReflectionCount',TRIM(StrVarNames(iStr)))) PP_nVarPart_loc = PP_nVarPart_loc - 1
    END DO
    DEALLOCATE(StrVarNames)

    ! Reflections are stored in the 8th data column. Do not start counting reflections mid-simulation
    IF (PartDataSize.EQ.PP_nVarPart_loc+1) THEN
        doParticleReflectionTrack = .FALSE.
        SWRITE(UNIT_stdOut,'(A3,A30,A3,I33)')' | ','Reflections not tracked previously. Disabling reflection counter',' | ',PartDataSize
    END IF

    ! Get PartData
    ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
    CALL ReadArray('PartData',2,(/PartDataSize,locnPart/),offsetnPart,2,RealArray=PartData)

    ! Loop over all particles on local proc and fill the PartState
    IF (locnPart.GT.0) THEN
      PartState(1:PP_nVarPart_loc,1:locnPart) = PartData(1:PP_nVarPart_loc,offsetnPart+1:offsetnPart+locnPart)
      PartSpecies(                1:locnPart) = INT(PartData(PP_nVarPart_loc+1,offsetnPart+1:offsetnPart+locnPart))
      ! For old files, where no particle diameter was saved
      IF (ANY(PartState(PART_DIAM,1:locnPart).EQ.0.)) THEN
        DO iPart = 1,locnPart
          PartState(PART_DIAM,iPart) = Species(PartSpecies(iPart))%DiameterIC
        END DO
      END IF
      ! Reflections were tracked previously and are therefore still enabled
      IF (PartDataSize.EQ.PP_nVarPart_loc+2) THEN
          PartReflCount(1:locnPart) = INT(PartData(PP_nVarPart_loc+2,offsetnPart+1:offsetnPart+locnPart))
      END IF

      ! Fill the particle-to-element-mapping (PEM) with the information from HDF5
      DO iElem = FirstElemInd,LastElemInd
        ! If the element contains a particle, add its ID to the element. Particles are still ordered along the SFC at this point
        IF (PartInt(ELEM_LastPartInd,iElem).GT.PartInt(ELEM_FirstPartInd,iElem)) THEN
          PEM%Element    (PartInt(ELEM_FirstPartInd,iElem)-offsetnPart+1 : &
                          PartInt(ELEM_LastPartInd ,iElem)-offsetnPart)  = iElem
          PEM%LastElement(PartInt(ELEM_FirstPartInd,iElem)-offsetnPart+1 : &
                          PartInt(ELEM_LastPartInd ,iElem)-offsetnPart)  = iElem
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
        SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' | HDF5 state file containing SGS/RW data but current run in DNS mode. Ignoring...'
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
        SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' | HDF5 state file not containing SGS/RW data but model active in current run.'
        SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' | SGS/RW model will start without history...'
      END IF
    END IF

    ! De-allocate arrays used for read-in
    DEALLOCATE(PartInt,PartData)
    PDM%ParticleVecLength = PDM%ParticleVecLength + locnPart
    CALL UpdateNextFreePosition()

    SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' READING PARTICLES FROM RESTARTFILE DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')

    ! Reconstruct the number of particles inserted before restart from the emission rate
    DO iSpec=1,nSpecies
      DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
        Species(iSpec)%Init(iInit)%InsertedParticle = INT(Species(iSpec)%Init(iInit)%ParticleEmission * RestartTime / Species(iSpec)%Init(iInit)%ParticleEmissionTime,8)
      END DO
    END DO

    ! Abort if HDF5 file contains more particles than memory allocated for in current run
    IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
      SWRITE (UNIT_stdOut,'(A,I0)') "PDM%ParticleVecLength =", PDM%ParticleVecLength
      SWRITE (UNIT_stdOut,'(A,I0)') "PDM%maxParticleNumber =", PDM%maxParticleNumber
      CALL abort(__STAMP__&
          ,' Number of Particles in Restart file is higher than MaxParticleNumber! Increase MaxParticleNumber!')
    END IF ! PDM%ParticleVecLength.GT.PDM%maxParticleNumber

    ! Since the elementside-local node number are NOT persistant and dependent on the location
    ! of the MPI borders, all particle-element mappings need to be checked after a restart
    ! Step 1: Identify particles that are not in the element in which they were before the restart
    NbrOfMissingParticles = 0
    NbrOfLostParticles    = 0

    SELECT CASE(TrackingMethod)
      CASE(TRIATRACKING)
        DO iPart = 1,PDM%ParticleVecLength
          IF (DoInterpolation) THEN
            CALL GetPositionInRefElem(PartState(1:3,iPart),Xi,PEM%Element(iPart))
            CNElemID = GetCNElemID(PEM%Element(iPart))

            ! Particle already inside the correct element
            IF(ALL(ABS(Xi).LE.1.0)) THEN
              InElementCheck=.TRUE.
              IF(ALLOCATED(PartPosRef)) PartPosRef(1:3,iPart)=Xi
            ELSE
              InElementCheck=.FALSE.
            END IF
          ELSE
            ! Check if particle is inside the correct element
            CALL ParticleInsideQuad3D(PartState(1:3,iPart),PEM%Element(iPart),InElementCheck,Det)
          END IF ! DoInterpolation

          ! Particle not in correct element, try to find them within MyProc
          IF (.NOT.InElementCheck) THEN
            NbrOfMissingParticles = NbrOfMissingParticles + 1
            CALL LocateParticleInElement(iPart,doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (.NOT.PDM%ParticleInside(iPart)) THEN
              NbrOfLostParticles = NbrOfLostParticles + 1
! #if !(USE_MPI)
!               IF (CountNbrOfLostParts) CALL StoreLostParticleProperties(iPart, PEM%Element(iPart), UsePartState_opt=.TRUE.)
! #endif /*!(USE_MPI)*/
            ELSE
              PEM%LastElement(iPart) = PEM%Element(iPart)
            END IF
          END IF
        END DO ! iPart = 1,PDM%ParticleVecLength

      CASE(TRACING)
        DO iPart = 1,PDM%ParticleVecLength
          ! Check if particle is inside the correct element
          CALL GetPositionInRefElem(PartState(1:3,iPart),Xi,PEM%Element(iPart))

          ! Particle already inside the correct element
          ! FLEXI has ElemEpsOneCell also with TRACING
          ! IF(ALL(ABS(Xi).LE.ElemEpsOneCell(CNElemID))) THEN
          IF (ALL(ABS(Xi).LE.1.0)) THEN
            InElementCheck = .TRUE.
            IF(ALLOCATED(PartPosRef)) PartPosRef(1:3,iPart)=Xi
          ELSE
            InElementCheck = .FALSE.
          END IF

          ! Particle not in correct element, try to find them within MyProc
          IF (.NOT.InElementCheck) THEN
            NbrOfMissingParticles = NbrOfMissingParticles + 1
            CALL LocateParticleInElement(iPart,doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (.NOT.PDM%ParticleInside(iPart)) THEN
              NbrOfLostParticles = NbrOfLostParticles + 1
! #if !(USE_MPI)
!               IF (CountNbrOfLostParts) CALL StoreLostParticleProperties(iPart, PEM%Element(iPart), UsePartState_opt=.TRUE.)
! #endif /*!(USE_MPI)*/
            ELSE
              PEM%LastElement(iPart) = PEM%Element(iPart)
            END IF ! .NOT.PDM%ParticleInside(iPart)
          END IF ! .NOT.InElementCheck
        END DO ! iPart = 1,PDM%ParticleVecLength

      CASE(REFMAPPING)
        DO iPart = 1,PDM%ParticleVecLength
          ! Check if particle is inside the correct element
          CALL GetPositionInRefElem(PartState(1:3,iPart),Xi,PEM%Element(iPart))
          CNElemID = GetCNElemID(PEM%Element(iPart))

          ! Particle already inside the correct element
          IF (ALL(ABS(Xi).LE.ElemEpsOneCell(CNElemID))) THEN ! particle inside
            InElementCheck    = .TRUE.
            PartPosRef(1:3,iPart) = Xi
          ELSE
            InElementCheck    = .FALSE.
          END IF

          ! Particle not in correct element, try to find them within MyProc
          IF (.NOT.InElementCheck) THEN
            NbrOfMissingParticles = NbrOfMissingParticles + 1
            CALL LocateParticleInElement(iPart,doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (.NOT.PDM%ParticleInside(iPart)) THEN
              NbrOfLostParticles = NbrOfLostParticles + 1
! #if !(USE_MPI)
!               IF (CountNbrOfLostParts) CALL StoreLostParticleProperties(iPart, PEM%Element(iPart), UsePartState_opt=.TRUE.)
! #endif /*!(USE_MPI)*/
              PartPosRef(1:3,iPart) = -888.
            ELSE
              PEM%LastElement(iPart) = PEM%Element(iPart)
            END IF
          END IF
        END DO ! iPart = 1,PDM%ParticleVecLength
    END SELECT

#if USE_MPI
    ! Step 2: All particles that are not found within MyProc need to be communicated to the others and located there
    ! Combine number of lost particles of all processes and allocate variables
    ! Note: Particles that are lost on MyProc are also searched for here again
    CALL MPI_ALLGATHER(NbrOfLostParticles, 1, MPI_INTEGER, TotalNbrOfMissingParticles, 1, MPI_INTEGER, PartMPI%COMM, IERROR)

    ! Check total number of missing particles and start re-locating them on other procs
    IF (TotalNbrOfMissingParticlesSum.GT.0) THEN
      ! Set offsets
      OffsetTotalNbrOfMissingParticles(0) = 0
      DO iProc = 1, PartMPI%nProcs-1
        OffsetTotalNbrOfMissingParticles(iProc) = OffsetTotalNbrOfMissingParticles(iProc-1) + TotalNbrOfMissingParticles(iProc-1)
      END DO ! iProc = 0, PartMPI%nProcs-1

      ALLOCATE(RecBuff(PartDataSize,1:TotalNbrOfMissingParticlesSum))

      ! Fill SendBuffer
      NbrOfMissingParticles = OffsetTotalNbrOfMissingParticles(PartMPI%MyRank) + 1
      DO iPart = 1, PDM%ParticleVecLength
        IF (.NOT.PDM%ParticleInside(iPart)) THEN
          RecBuff(1:PP_nVarPart,NbrOfMissingParticles) = PartState(1:PP_nVarPart,iPart)
          RecBuff(PP_nVarPart+1,NbrOfMissingParticles) = REAL(PartSpecies(iPart))
          ! Reflections were tracked previously and are therefore still enabled
          IF (PartDataSize.EQ.8) THEN
            RecBuff(PP_nVarPart+2,NbrOfMissingParticles) = REAL(PartReflCount(iPart))
          END IF
          NbrOfMissingParticles = NbrOfMissingParticles + 1
        END IF ! .NOT.PDM%ParticleInside(iPart)
      END DO ! iPart = 1, PDM%ParticleVecLength

      ! Distribute lost particles to all procs
      NbrOfMissingParticles = 0

      DO iProc = 0, PartMPI%nProcs-1
        RecCount(iProc) = TotalNbrOfMissingParticles(iProc)
        Displace(iProc) = NbrOfMissingParticles
        NbrOfMissingParticles = NbrOfMissingParticles + TotalNbrOfMissingParticles(iProc)
      END DO ! iProc = 0, PartMPI%nProcs-1

      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                     &
                         , 0                                                &
                         , MPI_DATATYPE_NULL                                &
                         , RecBuff                                          &
                         , PartDataSize*TotalNbrOfMissingParticles(:)       &
                         , PartDataSize*OffsetTotalNbrOfMissingParticles(:) &
                         , MPI_DOUBLE_PRECISION                             &
                         , PartMPI%COMM                                     &
                         , IERROR)

      ! Keep track which particles are found on the current proc
      ALLOCATE(IndexOfFoundParticles          (TotalNbrOfMissingParticlesSum))
      IF (MPIRoot) THEN
        ALLOCATE(CompleteIndexOfFoundParticles(TotalNbrOfMissingParticlesSum))
      END IF ! MPIRoot
      IndexOfFoundParticles = -1

      ! Free lost particle positions in local array to make room for missing particles that are tested
      CALL UpdateNextFreePosition()

      ! Add them to particle list and check if they are in MyProcs domain
      NbrOfFoundParts = 0
      CurrentPartNum  = PDM%ParticleVecLength+1

      DO iPart = 1, TotalNbrOfMissingParticlesSum
        ! Sanity check
        IF(CurrentPartNum.GT.PDM%maxParticleNumber)THEn
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " CurrentPartNum        = ",  CurrentPartNum
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " PDM%maxParticleNumber = ",  PDM%maxParticleNumber
          CALL abort(__STAMP__,'Missing particle ID > PDM%maxParticleNumber. Increase Part-MaxParticleNumber!')
        END IF ! CurrentPartNum.GT.PDM%maxParticleNumber

        ! Do not search particles twice: Skip my own particles, because these have already been searched for before they are
        ! sent to all other procs
        ASSOCIATE( myFirst => OffsetTotalNbrOfMissingParticles(PartMPI%MyRank) + 1 ,&
                   myLast  => OffsetTotalNbrOfMissingParticles(PartMPI%MyRank) + TotalNbrOfMissingParticles(PartMPI%MyRank))
          IF((iPart.GE.myFirst).AND.(iPart.LE.myLast))THEN
            IndexOfFoundParticles(iPart) = 0
            CYCLE
          END IF
        END ASSOCIATE

        PartState(1:PP_nVarPart,CurrentPartNum) = RecBuff(1:PP_nVarPart,iPart)
        PDM%ParticleInside(CurrentPartNum) = .true.

        CALL LocateParticleInElement(CurrentPartNum,doHALO=.FALSE.)
        IF (PDM%ParticleInside(CurrentPartNum)) THEN
          IndexOfFoundParticles(iPart) = 1
          PEM%LastElement(CurrentPartNum) = PEM%Element(CurrentPartNum)
          ! Set particle properties (if the particle is lost, it's properties are written to a .h5 file)
          PartSpecies(CurrentPartNum)     = INT(RecBuff(PP_nVarPart+1,iPart))
          ! Reflections were tracked previously and are therefore still enabled
          IF (PartDataSize.EQ.8) THEN
            PartReflCount(CurrentPartNum) = INT(RecBuff(PP_nVarPart+2,iPart))
          END IF
          NbrOfFoundParts = NbrOfFoundParts + 1
          CurrentPartNum  = CurrentPartNum  + 1
        ELSE ! Lost
          IndexOfFoundParticles(iPart) = 0
        END IF

        ! Sanity Check
        IF(IndexOfFoundParticles(iPart).EQ.-1)THEN
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " iPart                        : ",  iPart
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " IndexOfFoundParticles(iPart) : ",  IndexOfFoundParticles(iPart)
          CALL abort(__STAMP__,'IndexOfFoundParticles(iPart) was not set correctly)')
        END IF ! IndexOfFoundParticles(iPart)
      END DO ! iPart = 1, TotalNbrOfMissingParticlesSum

      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfFoundParts

      ! Combine number of found particles to make sure none are lost completely or found twice
      IF(MPIroot)THEN
        CALL MPI_REDUCE(IndexOfFoundParticles,CompleteIndexOfFoundParticles,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
      ELSE
        CALL MPI_REDUCE(IndexOfFoundParticles,0                            ,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
      END IF

      CompleteNbrOfFound      = 0
      CompleteNbrOfLost       = 0
      CompleteNbrOfDuplicate  = 0
      ! Only mpi root: Check if particle are not found or found twice
      IF (MPIRoot) THEN
        DO iPart = 1,TotalNbrOfMissingParticlesSum
          IF    (CompleteIndexOfFoundParticles(iPart).EQ.0) THEN ! permanently lost
            CompleteNbrOfLost       = CompleteNbrOfLost      + 1
          ELSEIF(CompleteIndexOfFoundParticles(iPart).EQ.1) THEN ! lost and found on one other processor
            CompleteNbrOfFound      = CompleteNbrOfFound     + 1
          ELSE ! lost and found multiple times on other processors
            CompleteNbrOfDuplicate  = CompleteNbrOfDuplicate + 1
          END IF

          ! Store the particle info
          IF(CountNbrOfLostParts)THEN
            CurrentPartNum = PDM%ParticleVecLength+1

            ! Set properties of the "virtual" particle (only for using the routine StoreLostParticleProperties to store this info
            ! in the .h5 container)
            PartState(1:PP_nVarPart,CurrentPartNum) = RecBuff(1:PP_nVarPart,iPart)
            PartSpecies(CurrentPartNum)             = INT(RecBuff(PP_nVarPart+1,iPart))
            ! Reflections were tracked previously and are therefore still enabled
            IF (PartDataSize.EQ.8) THEN
              PartReflCount(CurrentPartNum)      = INT(RecBuff(PP_nVarPart+2,iPart))
            END IF
            PEM%LastElement(CurrentPartNum)      = -1
            PDM%ParticleInside(CurrentPartNum)   = .FALSE.

            ! CALL StoreLostParticleProperties(CurrentPartNum, PEM%Element(CurrentPartNum), &
            !                                  UsePartState_opt=.TRUE., PartMissingType_opt=CompleteIndexOfFoundParticles(iPart))
          END IF ! CountNbrOfLostParts
        END DO

        WRITE(UNIT_stdOut,'(A,I0)') ' Particles initially lost during restart    : ',TotalNbrOfMissingParticlesSum
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles found on (other) procs : ',CompleteNbrOfFound
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles permanently lost       : ',CompleteNbrOfLost
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles found multiple times   : ',CompleteNbrOfDuplicate
        NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + CompleteNbrOfLost

        DEALLOCATE(CompleteIndexOfFoundParticles)

      END IF ! MPIRoot
      CALL MPI_BCAST(NbrOfLostParticlesTotal,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
      NbrOfLostParticlesTotal_old = NbrOfLostParticlesTotal

    END IF ! TotalNbrOfMissingParticlesSum.GT.0
#else /*not USE_MPI*/
    IF(NbrOfMissingParticles.GT.0)THEN
      NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + NbrOfLostParticles
      TotalNbrOfMissingParticlesSum = NbrOfMissingParticles
      NbrOfLostParticlesTotal_old = NbrOfLostParticlesTotal
      WRITE(UNIT_stdOut,'(A,I0)') ' Particles initially lost during restart : ',NbrOfMissingParticles
      WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles permanently lost    : ',NbrOfLostParticles
      WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles relocated           : ',NbrOfMissingParticles-NbrOfLostParticles
    END IF ! NbrOfMissingParticles.GT.0
#endif /*USE_MPI*/

    CALL UpdateNextFreePosition()
  ELSE
      SWRITE(UNIT_stdOut,'(A)') ' PartData does not exists in restart file'
      SWRITE(UNIT_stdOut,'(132("-"))')
  END IF ! PartDataExists
  CALL CloseDataFile()

  ! Get individual impact data
  IF (doParticleImpactTrack) CALL RestartParticleBoundaryTracking

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
