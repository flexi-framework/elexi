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
!> \brief Routines that handle particle restart capabilities.
!>
!> With this feature a simulation can be resumed from a state file that has been created during a previous
!> simulation (restart file). The restart file is passed to FLEXI as a second command line argument.
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
USE MOD_Particle_Globals
USE MOD_Particle_Readin
USE MOD_Particle_Restart_Vars
USE MOD_ReadInTools,                ONLY: PrintOption
! Boundary
USE MOD_Particle_Boundary_Analyze,  ONLY: CalcSurfaceValues
USE MOD_Particle_Boundary_Vars,     ONLY: ImpactRestart
! Localization
USE MOD_Particle_Localization,      ONLY: LocateParticleInElement,ParticleInsideQuad3D
USE MOD_Particle_Tracking_Vars,     ONLY: TrackingMethod,NbrOfLostParticles
USE MOD_Particle_Tracking_Vars,     ONLY: NbrOfLostParticlesTotal,TotalNbrOfMissingParticlesSum,NbrOfLostParticlesTotal_old
! Interpolation
USE MOD_Eval_XYZ,                   ONLY: GetPositionInRefElem
USE MOD_Particle_Interpolation_Vars,ONLY: DoInterpolation
USE MOD_Particle_Mesh_Vars,         ONLY: ElemEpsOneCell
! Mesh
USE MOD_Mesh_Vars,                  ONLY: offsetElem
USE MOD_Particle_Mesh_Tools,        ONLY: GetCNElemID
! Particles
USE MOD_Particle_Output,            ONLY: GetOffsetAndGlobalNumberOfParts
USE MOD_Part_Tools,                 ONLY: UpdateNextFreePosition
USE MOD_Particle_Vars,              ONLY: PartState,PartSpecies,PEM,PDM,Species,nSpecies
USE MOD_Particle_Vars,              ONLY: PartInt,PartData,TurbPartData
USE MOD_Particle_Vars,              ONLY: PartDataSize,TurbPartDataSize
USE MOD_Particle_Vars,              ONLY: PartPosRef,PartReflCount,doPartIndex,PartIndex
#if USE_MPI
USE MOD_Particle_Tracking_Vars,     ONLY: CountNbrOfLostParts
#endif /*USE_MPI*/
! Particle turbulence models
USE MOD_Particle_Vars,              ONLY: TurbPartState
#if PARABOLIC
USE MOD_Particle_SGS_Vars,          ONLY: nSGSVars
#if USE_RW
USE MOD_Particle_RandomWalk_Vars,   ONLY: nRWVars
#endif /*USE_RW*/
#endif
! Restart
USE MOD_Restart_Vars,               ONLY: RestartTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL     :: doFlushFiles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Parameters
INTEGER,PARAMETER                  :: ELEM_FirstPartInd = 1
INTEGER,PARAMETER                  :: ELEM_LastPartInd  = 2
! Counters
INTEGER                            :: iElem,CNElemID
INTEGER                            :: FirstElemInd,LastelemInd
INTEGER                            :: locnPart,offsetnPart
INTEGER                            :: iSpec,iPart,iInit
! Localization
LOGICAL                            :: InElementCheck
REAL                               :: xi(3)
REAL                               :: det(6,2)
INTEGER                            :: NbrOfMissingParticles
! Timers
REAL                               :: StartT,EndT
! MPI
#if USE_MPI
INTEGER,ALLOCATABLE                :: IndexOfFoundParticles(:),CompleteIndexOfFoundParticles(:)
INTEGER                            :: CompleteNbrOfLost,CompleteNbrOfFound,CompleteNbrOfDuplicate
REAL,ALLOCATABLE                   :: RecBuff(:,:)
INTEGER                            :: CurrentPartNum
INTEGER                            :: TotalNbrOfMissingParticles(      0:nProcessors-1)
INTEGER                            :: OffsetTotalNbrOfMissingParticles(0:nProcessors-1)
INTEGER                            :: Displace(                        0:nProcessors-1)
INTEGER                            :: RecCount(                        0:nProcessors-1)
INTEGER                            :: NbrOfFoundParts
INTEGER                            :: iProc
#endif /*USE_MPI*/
!===================================================================================================================================

! ===========================================================================
! Distribute or read the particle solution
! ===========================================================================
CALL ParticleReadin(doFlushFiles)

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' RESTARTING PARTICLES...'
GETTIME(StartT)

IF(PartIntExists)THEN
  IF(PartDataExists)THEN
    FirstElemInd = offsetElem+1
    LastElemInd  = offsetElem+PP_nElems
    locnPart     = PartInt(ELEM_LastPartInd,LastElemInd)-PartInt(ELEM_FirstPartInd,FirstElemInd)
    offsetnPart  = PartInt(ELEM_FirstPartInd,FirstElemInd)

    ! Loop over all particles on local proc and fill the PartState
    IF (locnPart.GT.0) THEN
      PartState(1:PP_nVarPartState,1:locnPart) = PartData(  1:PP_nVarPartState  ,offsetnPart+1:offsetnPart+locnPart)
      PartSpecies(                 1:locnPart) = INT(PartData(PP_nVarPartState+1,offsetnPart+1:offsetnPart+locnPart))
      ! Sanity check: SpecID > 0
      IF (ANY(PartSpecies(1:locnPart).LE.0)) THEN
        IPWRITE(UNIT_StdOut,'(I0,A,I0,A,I0,A,3(ES25.14E3),A,I0)') 'Warning: Found particle in restart file with SpecID.LE.0!'
        CALL Abort(__STAMP__,'Found particle in restart file with species ID zero, which indicates a corrupted restart file.')
      END IF
      ! For old files, where no particle diameter was saved
      IF (ANY(PartState(PART_DIAM,1:locnPart).EQ.0.)) THEN
        DO iPart = 1,locnPart
          PartState(PART_DIAM,iPart) = Species(PartSpecies(iPart))%DiameterIC
        END DO
      END IF
#if USE_SPHERICITY
      IF (ANY(PartState(PART_SPHE,1:locnPart).EQ.0.)) THEN
        DO iPart = 1,locnPart
          PartState(PART_SPHE,iPart) = Species(PartSpecies(iPart))%SphericityIC
        END DO
      END IF
#endif
      ! Particle index
      IF (doPartIndex) PartIndex(1:locnPart) = INT(PartData(PP_nVarPartState+2,offsetnPart+1:offsetnPart+locnPart))
      ! Reflections were tracked previously and are therefore still enabled
      IF (PartDataSize.EQ.PP_nVarPartState+2) &
        PartReflCount(1:locnPart) = INT(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))

      FirstElemInd = offsetElem+1
      LastElemInd  = offsetElem+PP_nElems

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
    END IF ! PartDataExists

#if PARABOLIC
    ! Total size of turbulent properties
    TurbPartDataSize = nSGSVars
#if USE_RW
    TurbPartDataSize = TurbPartDataSize + nRWVars
#endif
#endif

    ! Check if TurbPartData exists in HDF5 file
    IF(TurbPartDataExists)THEN
      ! Loop over all particles on local proc and fill the TurbPartState
      IF (locnPart.GT.0) THEN
        TurbPartState(1:TurbPartDataSize,1:locnPart) = TurbPartData(1:TurbPartDataSize,offsetnPart+1:offsetnPart+locnPart)
      END IF

      ! De-allocate array used for readin
      DEALLOCATE(TurbPartData)
    ! No turbulent particle properties in HDF5
    END IF ! TurbPartDataExists

    ! De-allocate arrays used for read-in
    DEALLOCATE(PartInt)
    DEALLOCATE(PartData)
    PDM%ParticleVecLength = PDM%ParticleVecLength + locnPart
    CALL UpdateNextFreePosition()

    ! Reconstruct the number of particles inserted before restart from the emission rate
    DO iSpec=1,nSpecies
      DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
        Species(iSpec)%Init(iInit)%InsertedParticle = INT(Species(iSpec)%Init(iInit)%ParticleEmission * (RestartTime - EmissionTime) / Species(iSpec)%Init(iInit)%ParticleEmissionTime,8)
      END DO
    END DO

    ! Abort if HDF5 file contains more particles than memory allocated for in current run
    IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
      SWRITE (UNIT_stdOut,'(A,I0)') "PDM%ParticleVecLength =", PDM%ParticleVecLength
      SWRITE (UNIT_stdOut,'(A,I0)') "PDM%maxParticleNumber =", PDM%maxParticleNumber
      CALL Abort(__STAMP__&
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
    CALL MPI_ALLGATHER(NbrOfLostParticles, 1, MPI_INTEGER, TotalNbrOfMissingParticles, 1, MPI_INTEGER, MPI_COMM_FLEXI, IERROR)

    ! Check total number of missing particles and start re-locating them on other procs
    IF (TotalNbrOfMissingParticlesSum.GT.0) THEN
      ! Set offsets
      OffsetTotalNbrOfMissingParticles(0) = 0
      DO iProc = 1,nProcessors-1
        OffsetTotalNbrOfMissingParticles(iProc) = OffsetTotalNbrOfMissingParticles(iProc-1) + TotalNbrOfMissingParticles(iProc-1)
      END DO ! iProc = 0,nProcessors-1

      ALLOCATE(RecBuff(PartDataSize,1:TotalNbrOfMissingParticlesSum))

      ! Fill SendBuffer
      NbrOfMissingParticles = OffsetTotalNbrOfMissingParticles(myRank) + 1
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

      DO iProc = 0,nProcessors-1
        RecCount(iProc) = TotalNbrOfMissingParticles(iProc)
        Displace(iProc) = NbrOfMissingParticles
        NbrOfMissingParticles = NbrOfMissingParticles + TotalNbrOfMissingParticles(iProc)
      END DO ! iProc = 0,nProcessors-1

      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                     &
                         , 0                                                &
                         , MPI_DATATYPE_NULL                                &
                         , RecBuff                                          &
                         , PartDataSize*TotalNbrOfMissingParticles(:)       &
                         , PartDataSize*OffsetTotalNbrOfMissingParticles(:) &
                         , MPI_DOUBLE_PRECISION                             &
                         , MPI_COMM_FLEXI                                   &
                         , IERROR)

      ! Keep track which particles are found on the current proc
      ALLOCATE(IndexOfFoundParticles          (TotalNbrOfMissingParticlesSum))
      IF (MPIRoot) &
        ALLOCATE(CompleteIndexOfFoundParticles(TotalNbrOfMissingParticlesSum))
      IndexOfFoundParticles = -1

      ! Free lost particle positions in local array to make room for missing particles that are tested
      CALL UpdateNextFreePosition()

      ! Add them to particle list and check if they are in MyProcs domain
      NbrOfFoundParts = 0
      CurrentPartNum  = PDM%ParticleVecLength+1

      DO iPart = 1, TotalNbrOfMissingParticlesSum
        ! Sanity check
        IF(CurrentPartNum.GT.PDM%maxParticleNumber)THEn
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') " CurrentPartNum        = ",  CurrentPartNum
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') " PDM%maxParticleNumber = ",  PDM%maxParticleNumber
          CALL Abort(__STAMP__,'Missing particle ID > PDM%maxParticleNumber. Increase Part-MaxParticleNumber!')
        END IF ! CurrentPartNum.GT.PDM%maxParticleNumber

        ! Do not search particles twice: Skip my own particles, because these have already been searched for before they are
        ! sent to all other procs
        ASSOCIATE( myFirst => OffsetTotalNbrOfMissingParticles(myRank) + 1 ,&
                   myLast  => OffsetTotalNbrOfMissingParticles(myRank) + TotalNbrOfMissingParticles(myRank))
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
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') " iPart                        : ",  iPart
          IPWRITE(UNIT_stdOut,'(I0,A,I0)') " IndexOfFoundParticles(iPart) : ",  IndexOfFoundParticles(iPart)
          CALL Abort(__STAMP__,'IndexOfFoundParticles(iPart) was not set correctly)')
        END IF ! IndexOfFoundParticles(iPart)
      END DO ! iPart = 1, TotalNbrOfMissingParticlesSum

      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfFoundParts

      ! Combine number of found particles to make sure none are lost completely or found twice
      IF(MPIRoot)THEN
        CALL MPI_REDUCE(IndexOfFoundParticles,CompleteIndexOfFoundParticles,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
      ELSE
        CALL MPI_REDUCE(IndexOfFoundParticles,0                            ,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
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

      CALL MPI_BCAST(NbrOfLostParticlesTotal,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
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
      SWRITE(UNIT_stdOut,'(A)') ' | PartData does not exists in restart file. No particles restarted.'
  END IF ! PartDataExists
END IF ! PartIntExists


! Make sure we update our surface data after reading the sampling data
#if USE_MPI
! Only procs with SurfMesh know about the restart so far. Communicate to all here
CALL MPI_ALLREDUCE(MPI_IN_PLACE,ImpactRestart,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_FLEXI,iError)
#endif

IF (ImpactRestart) CALL CalcSurfaceValues(restart_opt=.TRUE.)

! Communicate the total number and offset
CALL GetOffsetAndGlobalNumberOfParts('ParticleRestart',offsetnPart,nGlobalNbrOfParticles,locnPart,.TRUE.)

GETTIME(EndT)
SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' RESTARTING PARTICLES DONE! [',EndT-StartT,'s]'
SWRITE(UNIT_stdOut,'(132("-"))')

!#if USE_MPI
!CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
!#endif /*MPI*/

END SUBROUTINE ParticleRestart

END MODULE MOD_Particle_Restart
