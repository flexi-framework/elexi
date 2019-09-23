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
! Module for particle emission
!===================================================================================================================================
MODULE MOD_part_emission
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitializeParticleEmission
  MODULE PROCEDURE InitializeParticleEmission
END INTERFACE

INTERFACE ParticleInserting
  MODULE PROCEDURE ParticleInserting
END INTERFACE

INTERFACE SetParticleMass
  MODULE PROCEDURE SetParticleMass
END INTERFACE

INTERFACE SetParticleVelocity
  MODULE PROCEDURE SetParticleVelocity
END INTERFACE


!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC         :: InitializeParticleEmission, ParticleInserting
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeParticleEmission()
!===================================================================================================================================
! Initialize particles / Insert initial particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,      ONLY : Species,nSpecies,PDM,PEM
USE MOD_Part_Tools,         ONLY : UpdateNextFreePosition
USE MOD_Restart_Vars,       ONLY : DoRestart
#if USE_MPI
USE MOD_Particle_MPI_Vars,  ONLY : PartMPI
#endif /* MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i,NbrOfParticle,iInit,insertParticles
LOGICAL               :: EmType6
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' Initial particle inserting... '

CALL UpdateNextFreePosition()
EmType6=.false.

DO i=1, nSpecies
    DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
        IF ((Species(i)%Init(iInit)%ParticleEmissionType.EQ.6)) THEN
            EmType6=.true.
            EXIT
        END IF
    END DO
    IF (EmType6) EXIT
END DO

! FLEXI gained restart capability from pure fluid solution
IF (.NOT.DoRestart) THEN
  ! for the case of particle insertion per time, the inserted particle number for the current time must
  ! be updated. Otherwise, at the first timestep after restart, these particles will be inserted again
  ! Do insanity check of max. particle number compared to the number that is to be inserted for certain insertion types
  insertParticles = 0

  ! Get number of particles to be inserted. Divide it between MPI procs if required.
  DO i=1,nSpecies
      DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
          IF (TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
              IF (Species(i)%Init(iInit)%PartDensity.EQ.0) THEN
#if USE_MPI
                insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs)
#else
                insertParticles = insertParticles + Species(i)%Init(iInit)%initialParticleNumber
#endif
              ELSE
                  insertParticles = insertParticles + Species(i)%Init(iInit)%initialParticleNumber
              END IF
          ELSE IF ((TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cuboid') .OR.(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
#if USE_MPI
            insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs)
#else
            insertParticles = insertParticles + Species(i)%Init(iInit)%initialParticleNumber
#endif
          END IF
      END DO
  END DO

  ! Check if requested to insert more particles than the particle array can hold
  IF (insertParticles.GT.PDM%maxParticleNumber) THEN
#if USE_MPI
    WRITE(UNIT_stdOut,'(I0,A40,I0)')PartMPI%MyRank,' Maximum particle number : ',PDM%maxParticleNumber
    WRITE(UNIT_stdOut,'(I0,A40,I0)')PartMPI%MyRank,' To be inserted particles: ',insertParticles
#else
    WRITE(UNIT_stdOut,'(A40,I0)')                  ' Maximum particle number : ',PDM%maxParticleNumber
    WRITE(UNIT_stdOut,'(A40,I0)')                  ' To be inserted particles: ',insertParticles
#endif
    CALL abort(__STAMP__,'Number of to be inserted particles per init-proc exceeds max. particle number! ')
  END IF

  DO i = 1,nSpecies
    DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
      ! special emission type: constant density in cell, + to be used for init
      IF (((Species(i)%Init(iInit)%ParticleEmissionType .EQ. 4).OR.(Species(i)%Init(iInit)%ParticleEmissionType .EQ. 6)) .AND. &
           (Species(i)%Init(iInit)%UseForInit)) THEN
        CALL abort(__STAMP__,' particle pressure not moved to picasso!')
      ! no special emissiontype to be used
      ELSE IF (Species(i)%Init(iInit)%UseForInit) THEN
        ! initial particle number too large to handle
        IF(Species(i)%Init(iInit)%initialParticleNumber.GT.HUGE(1)) &
          CALL abort(__STAMP__,' Integer for initialParticleNumber too large!')

        NbrOfParticle = INT(Species(i)%Init(iInit)%initialParticleNumber,4)
#if USE_MPI
        ! set the particles at the correct position
        CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
        CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
        CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif /*MPI*/
        ! give the particles their correct velocity
        SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle velocities for species ',i,' ... '
        CALL SetParticleVelocity(i,iInit,NbrOfParticle)
        ! give the particles their correct (species) mass
        SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle and mass for species ',i,' ... '
        CALL SetParticleMass(i,NbrOfParticle)
        ! update number of particles on proc and find next free position in particle array
        PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
        CALL UpdateNextFreePosition()
      END IF
    END DO ! inits
  END DO ! species
END IF ! doRestart

!--- set last element to current element (needed when ParticlePush is not executed, e.g. "delay")
DO i = 1,PDM%ParticleVecLength
    PEM%lastElement(i) = PEM%Element(i)
END DO

SWRITE(UNIT_stdOut,'(A)') ' ...DONE '

END SUBROUTINE InitializeParticleEmission


#if USE_MPI
SUBROUTINE ParticleInserting(mode_opt)
#else
SUBROUTINE ParticleInserting()
#endif
!===================================================================================================================================
! Particle Inserting
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars     , ONLY: PartMPI
#endif /* MPI*/
USE MOD_Globals
USE MOD_Restart_Vars          , ONLY:RestartTime
USE MOD_Particle_Restart_Vars
USE MOD_Particle_Globals
USE MOD_Timedisc_Vars         , ONLY: dt,t,RKc,nRKStages,currentStage
USE MOD_Particle_Timedisc_Vars, ONLY: RKdtFrac,RKdtFracTotal
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_Part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze_Vars
USE MOD_Particle_Analyze       ,ONLY: CalcEkinPart
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if USE_MPI
INTEGER,INTENT(IN),OPTIONAL      :: mode_opt
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                          :: i,iInit,IntSample
INTEGER                , SAVE    :: NbrOfParticle=0
INTEGER(KIND=8)                  :: inserted_Particle_iter,inserted_Particle_time
INTEGER(KIND=8)                  :: inserted_Particle_diff
REAL                             :: PartIns, RandVal1
REAL                             :: RiseFactor, RiseTime
#if USE_MPI
INTEGER                          :: mode
INTEGER                          :: InitGroup
#endif
!===================================================================================================================================

!----------------------------------------------------------------------------------------------------------------------------------
!!! VORSICHT: FUNKTIONIERT SO MOMENTAN NUR MIT 1 SPEZIES!!!!
! --- fuer mehr als eine Spezies gibt es bei der Benutzung des mode_opt Flags Probleme mit den non-blocking communications. Es
!     koennte dann passieren, dass Nachrichten falsch zugeordnet werden. Sicherheitshalber sollte man kein mode_opt Argument bei
!     mehreren Spezies uebergeben.
#if USE_MPI
IF (PRESENT(mode_opt)) THEN
  mode=mode_opt
ELSE
  mode=0
END IF
#endif

! We are calling particleInserting differently in FLEXI. Calculate RKdtFrac here using Timedisc_Vars.
IF (currentStage.EQ.1) THEN
  RKdtFrac      = RKC(2)
  RKdtFracTotal = RKdtFrac
ELSE
  IF (currentStage.NE.nRKStages) THEN
    RKdtFrac      = RKC(currentStage+1)-RKC(currentStage)
    RKdtFracTotal = RKdtFracTotal+RKdtFrac
  ELSE
    RKdtFrac      = 1.-RKC(nRKStages)
    RKdtFracTotal = 1.
  END IF
END IF

!---  Emission at time step (initial emission see particle_init.f90: InitializeParticleEmission)
DO i=1,nSpecies
  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
    ! species to be used for init
    IF (Species(i)%Init(iInit)%UseForEmission) THEN
#if USE_MPI
      IF (mode.NE.2) THEN
#endif
        SELECT CASE(Species(i)%Init(iInit)%ParticleEmissionType)

          ! Emission Type: Particles per !!!!!SECOND!!!!!!!! (not per ns)
          CASE(1)
            IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
              ! requested particles during time-slab
              PartIns=Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac
              ! integer number of particles to be inserted
              inserted_Particle_iter = INT(PartIns,8)

              ! Attention: FLEXI uses a different definition for elapsed time due to restart capability!
              IF ((RestartTime.GT.0).AND.(.NOT.PartDataExists)) THEN
                ! total number of emitted particles over simulation
                PartIns=Species(i)%Init(iInit)%ParticleEmission * (t-RestartTime + dt*RKdtFracTotal)
              ELSE
                PartIns=Species(i)%Init(iInit)%ParticleEmission * (t + dt*RKdtFracTotal)
              END IF

              !-- random-round the inserted_Particle_time for preventing periodicity
              ! PO & SC: why, sometimes we do not want this add, TB is bad!
              CALL RANDOM_NUMBER(RandVal1)
              !-- repeat if inserting multiple particles
              IF (inserted_Particle_iter.GE.1) THEN
                CALL RANDOM_NUMBER(RandVal1)
                inserted_Particle_time = INT(PartIns + RandVal1,8) ! adds up to ONE
              ! needed, since InsertedParticleSurplus can increase
              ! and _iter > 1 needs to be possible for preventing periodicity
              ELSE IF((inserted_Particle_iter.GE.0).AND.(inserted_Particle_iter.LT.1)) THEN
                ! needed, since InsertedParticleSurplus can increase
                ! and _iter > 1 needs to be possible for preventing periodicity
                IF (ALMOSTEQUAL(PartIns,0.)) THEN !dummy
                  inserted_Particle_time = INT(PartIns,8)
                ELSE
                  !poisson-distri of PartIns-INT(PartIns)
                  CALL SamplePoissonDistri( PartIns-INT(PartIns) , IntSample )
                inserted_Particle_time = INT(INT(PartIns)+IntSample,8)
                END IF
              ELSE !dummy
                inserted_Particle_time = INT(PartIns,8)
              END IF

              !-- evaluate inserted_Particle_time and inserted_Particle_iter
              inserted_Particle_diff = inserted_Particle_time - Species(i)%Init(iInit)%InsertedParticle                  &
                                     - inserted_Particle_iter - Species(i)%Init(iInit)%InsertedParticleSurplus           &
                                     + Species(i)%Init(iInit)%InsertedParticleMisMatch

              Species(i)%Init(iInit)%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0))
              NbrOfParticle = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)

              !-- if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
              IF (TRIM(Species(i)%Init(iInit)%velocityDistribution).EQ.'maxwell') THEN
                IF (NbrOfParticle.LT.5) NbrOfParticle=0
              END IF

            ELSE IF (DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
              ! linear rise of inflow
              RiseTime=Species(i)%Init(iInit)%InflowRiseTime

              ! ramp up the particle insertion by increasing the RiseFactor depending on time
              IF(RiseTime.GT.0.)THEN
                IF(t-DelayTime.LT.RiseTime)THEN
                    RiseFactor=(t-DelayTime)/RiseTime
                ELSE
                    RiseFactor=1.
                END IF
              ELSE
                RiseFactor=1.
              END IF

              ! emitted particles during time-slab
              PartIns=Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac * RiseFactor
              CALL RANDOM_NUMBER(RandVal1)

              IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
                IPWRITE(*,*)'WARNING: target is too large for poisson sampling: switching now to Random rounding...'
                NbrOfParticle = INT(PartIns + RandVal1)
                DoPoissonRounding = .FALSE.
              ELSE
                !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects)
                !> [Tysanner and Garcia 2004]
                CALL SamplePoissonDistri( PartIns , NbrOfParticle , DoPoissonRounding)
              END IF

            ! DoTimeDepInflow
            ELSE
              ! linear rise of inflow
              RiseTime=Species(i)%Init(iInit)%InflowRiseTime
              IF (RiseTime.GT.0.) THEN
                IF (t-DelayTime.LT.RiseTime) THEN
                  RiseFactor=(t-DelayTime)/RiseTime
                ELSE
                  RiseFactor=1.
                END IF
              ELSE
                 RiseFactor=1.
              END IF

              ! emitted particles during time-slab
              PartIns = Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac * RiseFactor &
                      + Species(i)%Init(iInit)%InsertedParticleMisMatch
              CALL RANDOM_NUMBER(RandVal1)
              NbrOfParticle = INT(PartIns + RandVal1)
            END IF
#if USE_MPI
            ! communicate number of particles with all procs in the same init group
            InitGroup=Species(i)%Init(iInit)%InitCOMM
            IF(PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL) THEN
              ! only procs which are part of group take part in the communication
              ! NbrOfParticle based on RandVals!
              CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
            ELSE
              NbrOfParticle=0
            END IF
#endif
            Species(i)%Init(iInit)%InsertedParticle = Species(i)%Init(iInit)%InsertedParticle + INT(NbrOfParticle,8)

          ! Emission Type: Particles per Iteration
          CASE(2)
            ! insert in last stage only, so that no reconstruction is necessary and number/iter matches
            IF (RKdtFracTotal .EQ. 1.) THEN
              NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission)
            ELSE
              NbrOfParticle = 0
            END IF

          CASE DEFAULT
            CALL abort(__STAMP__,'Unknown particle emission type')
            ! Line below not required anymore, but might want to change back to silent fail later
            ! NbrOfParticle = 0
          END SELECT

#if USE_MPI
          CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
        END IF

        IF (mode.NE.1) THEN
          CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
          CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif
          CALL SetParticleVelocity(i,iInit,NbrOfParticle)
          CALL SetParticleMass(i,NbrOfParticle)

         ! instead of UpdateNextfreePosition we update the particleVecLength only and doing it later, after CalcPartBalance
         PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
         PDM%ParticleVecLength       = PDM%ParticleVecLength       + NbrOfParticle
         !CALL UpdateNextFreePosition()
#if USE_MPI
      END IF
#endif
    END IF ! UseForEmission
  END DO
END DO

END SUBROUTINE ParticleInserting


#if USE_MPI
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle,mode)
#else
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle)
#endif /* MPI*/
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
#if USE_MPI
USE MOD_Particle_MPI_Vars,     ONLY: PartMPI,PartMPIInsert
#endif /* MPI*/
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_PIC_Vars
USE MOD_Particle_Vars,         ONLY: Species,PDM,BoltzmannConst,PartState
USE MOD_Particle_Mesh_Vars,    ONLY: GEO
USE MOD_Timedisc_Vars,         ONLY: dt
USE MOD_Particle_TimeDisc_Vars,ONLY: RKdtFrac
USE MOD_Particle_Mesh,         ONLY: SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Tracking_Vars,ONLY: DoRefMapping, TriaTracking
USE MOD_PICInterpolation,      ONLY: InterpolateVariableExternalField
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr, iInit
#if USE_MPI
INTEGER,INTENT(IN)                       :: mode
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                                  :: iProc,tProc, CellX, CellY, CellZ
INTEGER                                  :: msg_status(1:MPI_STATUS_SIZE)
INTEGER                                  :: MessageSize
LOGICAL                                  :: InsideMyBGM
#endif
REAL,ALLOCATABLE                         :: particle_positions(:)
INTEGER                                  :: allocStat
INTEGER                                  :: i,j,k,ParticleIndexNbr
INTEGER                                  :: mySumOfMatchedParticles, sumOfMatchedParticles
INTEGER                                  :: nChunks, chunkSize, chunkSize2
REAL                                     :: lineVector(3),VectorGap(3)
REAL                                     :: RandVal(3), Particle_pos(3),lineVector2(3)
REAL                                     :: RandVal1
REAL                                     :: radius, argumentTheta
REAL                                     :: pilen
REAL                                     :: x_step, y_step, z_step,  x_pos , y_pos
REAL                                     :: xlen, ylen, zlen
INTEGER                                  :: iPart
REAL                                     :: PartIns
INTEGER                                  :: DimSend
LOGICAL                                  :: insideExcludeRegion
#if USE_MPI
INTEGER                                  :: InitGroup,nChunksTemp,mySumOfRemovedParticles
INTEGER,ALLOCATABLE                      :: PartFoundInProc(:,:) ! 1 proc id, 2 local part id
#endif
!===================================================================================================================================

! emission group communicator
#if USE_MPI
InitGroup=Species(FractNbr)%Init(iInit)%InitCOMM
IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
  NbrofParticle=0
  RETURN
END IF
#endif /*MPI*/

PartIns      = 0.
lineVector   = 0.
Particle_pos = 0.
DimSend      = 3            !save (and send) only positions

IF ( (NbrOfParticle .LE. 0).AND.(PartIns .LE. 0.).AND. (ABS(Species(FractNbr)%Init(iInit)%PartDensity).LE.0.) ) &
  RETURN !0<Partins<1: statistical handling of exact REAL-INT-conv. below!

! Standard: Non-MPI
nChunks                 = 1
sumOfMatchedParticles   = 0
mySumOfMatchedParticles = 0
chunkSize               = nbrOfParticle

! for equidistant distribution, process myRank=0 generates the complete list of random positions for all emitted particles
#if USE_MPI
IF(( (nbrOfParticle.GT.PartMPI%InitGroup(InitGroup)%nProcs*10                             ) .AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'circle_equidistant'                 ) .AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'sin_deviation'                      ) .AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'cuboid_with_equidistant_distribution').AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'line_with_equidistant_distribution' )))THEN
   nChunks = PartMPI%InitGroup(InitGroup)%nProcs
ELSE
   nChunks = 1
END IF

! communication
IF(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle') nChunks=1

IF (mode.EQ.1) THEN
  chunkSize = INT(nbrOfParticle/nChunks)
  IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
    chunkSize = chunkSize + ( nbrOfParticle - (nChunks*chunkSize) )
  END IF
  IF (PartMPI%InitGroup(InitGroup)%MPIROOT .OR. nChunks.GT.1) THEN
#endif
    ALLOCATE( particle_positions(1:chunkSize*DimSend), STAT=allocStat )
    IF (allocStat .NE. 0) THEN
      CALL abort(&
__STAMP__&
,'ERROR in SetParticlePosition: cannot allocate particle_positions!')
    END IF

    chunkSize2=chunkSize !will be changed during insertion for:
                         !  1.: excludeRegions (orig. chunksize is for SpaceIC without taking excludeRegions into account)
    !------------------SpaceIC-cases: start-----------------------------------------------------------!
    SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
    !------------------SpaceIC-case: point------------------------------------------------------------------------------------------
    CASE ('point')
       Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC
       DO i=1,chunkSize
          particle_positions(i*3-2) = Particle_pos(1)
          particle_positions(i*3-1) = Particle_pos(2)
          particle_positions(i*3  ) = Particle_pos(3)
       END DO
    !------------------SpaceIC-case: line_with_equidistant_distribution-------------------------------------------------------------
    CASE ('line_with_equidistant_distribution')
      IF(NbrOfParticle.EQ.1)THEN
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + 0.5 * Species(FractNbr)%Init(iInit)%BaseVector1IC
      ELSE
        VectorGap = Species(FractNbr)%Init(iInit)%BaseVector1IC/(NbrOfParticle-1)
        DO i=1,chunkSize
          Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + (i-1)*VectorGap
          particle_positions(i*3-2) = Particle_pos(1)
          particle_positions(i*3-1) = Particle_pos(2)
          particle_positions(i*3  ) = Particle_pos(3)
        END DO
      END IF
    !------------------SpaceIC-case: line-------------------------------------------------------------------------------------------
    CASE ('line')
      DO i=1,chunkSize
        CALL RANDOM_NUMBER(RandVal1)
        Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC*RandVal1
        particle_positions(i*3-2) = Particle_pos(1)
        particle_positions(i*3-1) = Particle_pos(2)
        particle_positions(i*3  ) = Particle_pos(3)
      END DO
    !------------------SpaceIC-case: disc-------------------------------------------------------------------------------------------
    CASE('disc')
      IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
        lineVector(1) = 1.0
        lineVector(2) = 1.0
        lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                         Species(FractNbr)%Init(iInit)%NormalIC(3)
      ELSE
        IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
          lineVector(1) = 1.0
          lineVector(3) = 1.0
          lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                            Species(FractNbr)%Init(iInit)%NormalIC(2)
        ELSE
          IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
            lineVector(2) = 1.0
            lineVector(3) = 1.0
            lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                 Species(FractNbr)%Init(iInit)%NormalIC(1)
          ELSE
            CALL abort(&
__STAMP__&
,'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein')
          END IF
        END IF
      END IF

      lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
           lineVector(2) + lineVector(3) * lineVector(3))

      lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
           Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
      lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
           Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
      lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
           Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

      lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
           lineVector2(2) + lineVector2(3) * lineVector2(3))

      DO i=1,chunkSize
         radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1
         DO WHILE(radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC)
            CALL RANDOM_NUMBER(RandVal)
            RandVal = RandVal * 2. - 1.
            Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * &
                     (RandVal(1) * lineVector + RandVal(2) *lineVector2)

            radius = SQRT( (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) * &
                           (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) + &
                           (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) * &
                           (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) + &
                           (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) * &
                           (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) )
         END DO
         particle_positions(i*3-2) = Particle_pos(1)
         particle_positions(i*3-1) = Particle_pos(2)
         particle_positions(i*3  ) = Particle_pos(3)
      END DO
    CASE('circle')
      IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
         lineVector(1) = 1.0
         lineVector(2) = 1.0
         lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                           Species(FractNbr)%Init(iInit)%NormalIC(3)
      ELSE
         IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
            lineVector(1) = 1.0
            lineVector(3) = 1.0
            lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                              Species(FractNbr)%Init(iInit)%NormalIC(2)
         ELSE
            IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
               lineVector(2) = 1.0
               lineVector(3) = 1.0
               lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                 Species(FractNbr)%Init(iInit)%NormalIC(1)
            ELSE
              CALL abort(&
__STAMP__&
,'Error in SetParticlePosition, NormalIC should not be zero')
            END IF
         END IF
      END IF

      lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
           lineVector(2) + lineVector(3) * lineVector(3))

      lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
           Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
      lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
           Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
      lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
           Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

      lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
           lineVector2(2) + lineVector2(3) * lineVector2(3))

      radius = Species(FractNbr)%Init(iInit)%RadiusIC
      DO i=1,chunkSize
         CALL RANDOM_NUMBER(RandVal1)
         argumentTheta = 2.*pi*RandVal1
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                        linevector * cos(argumentTheta) * radius +  &
                        linevector2 * sin(argumentTheta) * radius
         particle_positions(i*3-2) = Particle_pos(1)
         particle_positions(i*3-1) = Particle_pos(2)
         particle_positions(i*3  ) = Particle_pos(3)
      END DO
    !------------------SpaceIC-case: circle_equidistant-----------------------------------------------------------------------------
    CASE('circle_equidistant')
      IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
         lineVector(1) = 1.0
         lineVector(2) = 1.0
         lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                           Species(FractNbr)%Init(iInit)%NormalIC(3)
      ELSE
         IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
            lineVector(1) = 1.0
            lineVector(3) = 1.0
            lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                              Species(FractNbr)%Init(iInit)%NormalIC(2)
         ELSE
            IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
               lineVector(2) = 1.0
               lineVector(3) = 1.0
               lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                 Species(FractNbr)%Init(iInit)%NormalIC(1)
            ELSE
              CALL abort(&
__STAMP__&
,'Error in SetParticlePosition, NormalIC should not be zero')
            END IF
         END IF
      END IF

      lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
           Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
      lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
           Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
      lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
           Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

      lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
           lineVector(2) + lineVector(3) * lineVector(3))

      lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
           lineVector2(2) + lineVector2(3) * lineVector2(3))

      radius = Species(FractNbr)%Init(iInit)%RadiusIC
      DO i=1,chunkSize
         argumentTheta = 2.*pi*i/chunkSize
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                        linevector * cos(argumentTheta) * radius +  &
                        linevector2 * sin(argumentTheta) * radius
         particle_positions(i*3-2) = Particle_pos(1)
         particle_positions(i*3-1) = Particle_pos(2)
         particle_positions(i*3  ) = Particle_pos(3)
      END DO
    !------------------SpaceIC-case: cuboid-----------------------------------------------------------------------------------------
    CASE('cuboid')
      lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
        CALL abort(&
__STAMP__&
,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
         CALL RANDOM_NUMBER(RandVal)
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1)
         Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2)
         IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
         ELSE
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * RandVal(3)
         END IF
         IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
           CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
           IF (insideExcludeRegion) THEN
             i=i+1
             CYCLE !particle is in excluded region
           END IF
         END IF
         particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
         particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
         particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
         i=i+1
         chunkSize2=chunkSize2+1
      END DO
    !------------------SpaceIC-case: cylinder---------------------------------------------------------------------------------------
    CASE('cylinder')
      lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
        CALL abort(&
__STAMP__&
,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
         radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
         DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
            CALL RANDOM_NUMBER(RandVal)
            Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC * (RandVal(1)*2-1) &
                         + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2-1)
            radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                           Particle_pos(2) * Particle_pos(2) + &
                           Particle_pos(3) * Particle_pos(3) )
         END DO
         Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
         IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
         ELSE
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal(3)
         END IF
         IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
           CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
           IF (insideExcludeRegion) THEN
             i=i+1
             CYCLE !particle is in excluded region
           END IF
         END IF
         particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
         particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
         particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
         i=i+1
         chunkSize2=chunkSize2+1
      END DO
    !------------------SpaceIC-case: cuboid_equal-----------------------------------------------------------------------------------
    CASE('cuboid_equal')
#if USE_MPI
      IF (PartMPI%InitGroup(InitGroup)%nProcs.GT. 1) THEN
        SWRITE(UNIT_stdOut,*)'WARNING in SetParticlePosition:'
        SWRITE(UNIT_stdOut,*)'cannot fully handle Particle Initial Condition \"cuboid equal\"'
        SWRITE(UNIT_stdOut,*)'in parallel mode (with more than one CPU)!'
        SWRITE(UNIT_stdOut,*)'USE WITH CARE!!!'
      END IF
      j=0
      mySumOfMatchedParticles = 0
      DO i=1,PDM%ParticleVecLength
         j=j+1
         ParticleIndexNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (ParticleIndexNbr .ne. 0) THEN
            PartState(ParticleIndexNbr,1:3) = PartState(j,1:3)
            PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
            IF(DoRefMapping.OR.TriaTracking)THEN
              CALL SingleParticleToExactElement(ParticleIndexNbr,doHALO=.FALSE.,initFix=.TRUE.,doRelocate=.FALSE.)
            ELSE
              CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.,doRelocate=.FALSE.)
            END IF
            IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
               mySumOfMatchedParticles = mySumOfMatchedParticles + 1
            ELSE
               PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
            END IF
            IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
              PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
              PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
            END IF
         ELSE
           CALL abort(&
__STAMP__&
,'ERROR in SetParticlePosition: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
         END IF
      END DO
      CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER &
          , MPI_SUM, PartMPI%InitGroup(InitGroup)%COMM, IERROR)
      nbrOfParticle = NbrOfParticle - sumOfMatchedParticles
      IF (nbrOfParticle .NE. 0) THEN
        IPWRITE(UNIT_stdOut,*)'ERROR in ParticleEmission_parallel:'
        IPWRITE(UNIT_stdOut,'(I4,A,I8,A)')'matched ', sumOfMatchedParticles, ' particles'
        IPWRITE(UNIT_stdOut,'(I4,A,I8,A)')'when ', NbrOfParticle+sumOfMatchedParticles, ' particles were required!'
        CALL abort(&
__STAMP__&
,'ERROR in ParticleEmission_parallel')
      END IF
      NbrOfParticle = mySumOfMatchedParticles
      DEALLOCATE( particle_positions, STAT=allocStat )
      IF (allocStat .NE. 0) THEN
        CALL abort(&
__STAMP__&
,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
      END IF
      RETURN
#else
      DO i=1,chunkSize
        ParticleIndexNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        particle_positions(i*3-2 : i*3) = PartState(ParticleIndexNbr-Species(FractNbr)%Init(iInit)%initialParticleNumber,1:3)
      END DO
#endif
    !------------------SpaceIC-case: cuboid_with_equidistant_distribution-----------------------------------------------------------
    CASE ('cuboid_with_equidistant_distribution')
       IF(Species(FractNbr)%Init(iInit)%initialParticleNumber.NE. &
            (Species(FractNbr)%Init(iInit)%maxParticleNumberX * Species(FractNbr)%Init(iInit)%maxParticleNumberY &
            * Species(FractNbr)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',FractNbr,' does not match number of particles in each direction!'
         CALL abort(&
__STAMP__&
,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = SQRT(Species(FractNbr)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(FractNbr)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(FractNbr)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(FractNbr)%Init(iInit)%BaseVector1IC(1)).OR. &
          (ylen.NE.Species(FractNbr)%Init(iInit)%BaseVector2IC(2)).OR. &
          (zlen.NE.Species(FractNbr)%Init(iInit)%CuboidHeightIC)) THEN
         CALL abort(&
__STAMP__&
,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
        END IF
       x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
       y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
       z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
       iPart = 1
       DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
         x_pos = (i-0.5) * x_step + Species(FractNbr)%Init(iInit)%BasePointIC(1)
         DO j=1,Species(FractNbr)%Init(iInit)%maxParticleNumberY
           y_pos =  Species(FractNbr)%Init(iInit)%BasePointIC(2) + (j-0.5) * y_step
           DO k=1,Species(FractNbr)%Init(iInit)%maxParticleNumberZ
             particle_positions(iPart*3-2) = x_pos
             particle_positions(iPart*3-1) = y_pos
             particle_positions(iPart*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                  + (k-0.5) * z_step
             iPart = iPart + 1
           END DO
         END DO
       END DO
    !------------------SpaceIC-case: sin_deviation----------------------------------------------------------------------------------
    CASE('sin_deviation')
       IF(Species(FractNbr)%Init(iInit)%initialParticleNumber.NE. &
            (Species(FractNbr)%Init(iInit)%maxParticleNumberX * Species(FractNbr)%Init(iInit)%maxParticleNumberY &
            * Species(FractNbr)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',FractNbr,' does not match number of particles in each direction!'
         CALL abort(&
         __STAMP__&
         ,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = abs(GEO%xmaxglob  - GEO%xminglob)
       ylen = abs(GEO%ymaxglob  - GEO%yminglob)
       zlen = abs(GEO%zmaxglob  - GEO%zminglob)
       pilen=2.0*PI/xlen
       x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
       y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
       z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
       iPart = 1
       DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
          x_pos = (i * x_step - x_step*0.5)
          x_pos = GEO%xminglob + x_pos + Species(FractNbr)%Init(iInit)%Amplitude &
                  * sin(Species(FractNbr)%Init(iInit)%WaveNumber * pilen * x_pos)
          DO j=1,Species(FractNbr)%Init(iInit)%maxParticleNumberY
            y_pos =  GEO%yminglob + j * y_step - y_step * 0.5
            DO k=1,Species(FractNbr)%Init(iInit)%maxParticleNumberZ
              particle_positions(iPart*3-2) = x_pos
              particle_positions(iPart*3-1) = y_pos
              particle_positions(iPart*3  ) = GEO%zminglob &
                                        + k * z_step - z_step * 0.5
              iPart = iPart + 1
            END DO
          END DO
       END DO
    END SELECT
    !------------------SpaceIC-cases: end-----------------------------------------------------------!
    chunkSize=chunkSize2

#if USE_MPI
 ELSE !no mpi root, nchunks=1
   chunkSize=0
 END IF
 IF(nChunks.GT.1) THEN
   ALLOCATE( PartMPIInsert%nPartsSend  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   ALLOCATE( PartMPIInsert%nPartsRecv  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   ALLOCATE( PartMPIInsert%SendRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
   ALLOCATE( PartMPIInsert%RecvRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
   ALLOCATE( PartMPIInsert%send_message(0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   PartMPIInsert%nPartsSend(:)=0
   DO i=1,chunkSize
     CellX = INT((particle_positions(DimSend*(i-1)+1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
     CellY = INT((particle_positions(DimSend*(i-1)+2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
     CellZ = INT((particle_positions(DimSend*(i-1)+3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
     InsideMyBGM=.TRUE.
     IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
         (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
         (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
       InsideMyBGM=.FALSE.
     END If
     IF (InsideMyBGM) THEN
       IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
     END IF
     IF (InsideMyBGM) THEN
       DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
         iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
         tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(iProc)
         IF(tProc.EQ.-1)CYCLE
         !IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
         PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
       END DO
       PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)=&
              PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)+1
     ELSE
       DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
!         IF (iProc.EQ.PartMPI%iProc) CYCLE
         PartMPIInsert%nPartsSend(iProc)=PartMPIInsert%nPartsSend(iProc)+1
       END DO
     END IF
   END DO
 ELSE
    IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
      ALLOCATE( PartMPIInsert%send_message(0:0), STAT=allocStat )
      MessageSize=DimSend*chunkSize
      ALLOCATE( PartMPIInsert%send_message(0)%content(1:MessageSize), STAT=allocStat )
      PartMPIInsert%send_message(0)%content(:)=particle_positions(1:DimSend*chunkSize)
      DEALLOCATE(particle_positions, STAT=allocStat)
    END IF
 END IF
 IF (nChunks.GT.1) THEN
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      ! sent particles
      !--- MPI_ISEND lengths of lists of particles leaving local mesh
      CALL MPI_ISEND(PartMPIInsert%nPartsSend(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, &
                     PartMPIInsert%SendRequest(iProc,1), IERROR)
      !--- MPI_IRECV lengths of lists of particles entering local mesh
      CALL MPI_IRECV(PartMPIInsert%nPartsRecv(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, &
                     PartMPIInsert%RecvRequest(iProc,1), IERROR)
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        ALLOCATE( PartMPIInsert%send_message(iProc)%content(1:DimSend*PartMPIInsert%nPartsSend(iProc)), STAT=allocStat )
      END IF
    END DO
    PartMPIInsert%nPartsSend(:)=0
    DO i=1,chunkSize
      CellX = INT((particle_positions(DimSend*(i-1)+1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
      CellY = INT((particle_positions(DimSend*(i-1)+2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
      CellZ = INT((particle_positions(DimSend*(i-1)+3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
      InsideMyBGM=.TRUE.
      IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
          (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
          (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
        InsideMyBGM=.FALSE.
      END If
      IF (InsideMyBGM) THEN
        IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
      END IF
      IF (InsideMyBGM) THEN
        DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
          iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
          tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(iProc)
          IF(tProc.EQ.-1)CYCLE
          PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
          k=PartMPIInsert%nPartsSend(tProc)
          PartMPIInsert%send_message(tProc)%content(DimSend*(k-1)+1:DimSend*k)=particle_positions(DimSend*(i-1)+1:DimSend*i)
        END DO
        PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)= &
            PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)+1
        k=PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)
        PartMPIInsert%send_message(PartMPI%InitGroup(InitGroup)%MyRank)%content(DimSend*(k-1)+1:DimSend*k)=&
                                                                          particle_positions(DimSend*(i-1)+1:DimSend*i)
      ELSE
        DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
 !         IF (iProc.EQ.PartMPI%iProc) CYCLE
          PartMPIInsert%nPartsSend(iProc)=PartMPIInsert%nPartsSend(iProc)+1
          k=PartMPIInsert%nPartsSend(iProc)
          PartMPIInsert%send_message(iProc)%content(DimSend*(k-1)+1:DimSend*k)=particle_positions(DimSend*(i-1)+1:DimSend*i)
        END DO
      END IF
    END DO
    DEALLOCATE(particle_positions, STAT=allocStat)
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      !--- (non-blocking:) send messages to all procs receiving particles from myself
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        CALL MPI_ISEND(PartMPIInsert%send_message(iProc)%content, DimSend*PartMPIInsert%nPartsSend(iProc),&
         MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%SendRequest(iProc,2), IERROR)
      END IF
    END DO
  END IF
ELSE ! mode.NE.1:
!--- RECEIVE:
  nChunksTemp=0
  IF(nChunks.EQ.1) THEN
    IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN !chunkSize can be 1 higher than NbrOfParticle for VPI+PartDens
       chunkSize=INT( REAL(SIZE(PartMPIInsert%send_message(0)%content)) / REAL(DimSend) )
       ALLOCATE(particle_positions(1:chunkSize*DimSend), STAT=allocStat)
       particle_positions(:)=PartMPIInsert%send_message(0)%content(:)
       DEALLOCATE( PartMPIInsert%send_message(0)%content )
       DEALLOCATE( PartMPIInsert%send_message )
    END IF
    chunkSize=NbrOfParticle
    IF(.NOT.PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
      ALLOCATE(particle_positions(1:chunkSize*DimSend), STAT=allocStat)
    END IF
    CALL MPI_BCAST(particle_positions, chunkSize*DimSend, MPI_DOUBLE_PRECISION,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
    nChunksTemp=1
  ELSE
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,1),msg_status(:),IERROR)
    END DO
    k=SUM(PartMPIInsert%nPartsRecv)
    ALLOCATE(particle_positions(1:k*DimSend), STAT=allocStat)
    k=0
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      IF (PartMPIInsert%nPartsRecv(iProc).GT.0) THEN
      !--- MPI_IRECV lengths of lists of particles entering local mesh
        CALL MPI_IRECV(particle_positions(k*DimSend+1), DimSend*PartMPIInsert%nPartsRecv(iProc),&
                                                  MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr,   &
                                                  PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%RecvRequest(iProc,2), IERROR)
        CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,2),msg_status(:),IERROR)
        k=k+PartMPIInsert%nPartsRecv(iProc)
      END IF
    END DO
    DEALLOCATE( PartMPIInsert%nPartsRecv )
    DEALLOCATE( PartMPIInsert%RecvRequest )
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      CALL MPI_WAIT(PartMPIInsert%SendRequest(iProc,1),msg_status(:),IERROR)
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        CALL MPI_WAIT(PartMPIInsert%SendRequest(iProc,2),msg_status(:),IERROR)
        DEALLOCATE( PartMPIInsert%send_message(iProc)%content )
      END IF
    END DO
    DEALLOCATE( PartMPIInsert%nPartsSend )
    DEALLOCATE( PartMPIInsert%send_message )
    DEALLOCATE( PartMPIInsert%SendRequest )
    chunkSize=k
    nChunks=1
  END IF
#endif
   ! each process checks which particle can be matched to its elements, counting the elements inside (local particles)
!   WRITE(*,*)'locating',chunkSize,'*',nChunks,' particles...'
!   WRITE(UNIT=debugFileName,FMT='(A,I2.2)')'prtcls_',PartMPI%iProc
!   OPEN(UNIT=130+PartMPI%iProc,FILE=debugFileName)
!   DO i=1,chunkSize*nChunks
!      WRITE(130+PartMPI%iProc,'(3(ES15.8))')particle_positions(i*3-2:i*3)
!   END DO
!   CLOSE(130+PartMPI%iProc)

#if USE_MPI
  ! in order to remove duplicated particles
  IF(nChunksTemp.EQ.1) THEN
    ALLOCATE(PartFoundInProc(1:2,1:ChunkSize),STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) THEN
        CALL abort(&
__STAMP__,&
"abort: Error during emission in PartFoundInProc allocation")
      END IF
    PartFoundInProc=-1
  END IF
#endif /*MPI*/

  mySumOfMatchedParticles=0
  ParticleIndexNbr = 1
  DO i=1,chunkSize*nChunks
    IF ((i.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
       ParticleIndexNbr = PDM%nextFreePosition(mySumOfMatchedParticles + 1 &
                                             + PDM%CurrentNextFreePosition)
    END IF
    IF (ParticleIndexNbr .ne. 0) THEN
       PartState(ParticleIndexNbr,1:DimSend) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+DimSend)
       PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
       IF(DoRefMapping.OR.TriaTracking)THEN
         CALL SingleParticleToExactElement(ParticleIndexNbr,doHALO=.FALSE.,InitFix=.TRUE.,doRelocate=.FALSE.)
       ELSE
         CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.,doRelocate=.FALSE.)
       END IF
       IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
          mySumOfMatchedParticles = mySumOfMatchedParticles + 1
#if USE_MPI
          IF(nChunksTemp.EQ.1) THEN
            ! mark elements with Rank and local found particle index
            PartFoundInProc(1,i)=MyRank
            PartFoundInProc(2,i)=mySumOfMatchedParticles
          END IF ! nChunks.EQ.1
#endif /*MPI*/
       ELSE
          PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
       END IF
       IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
         PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
         PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
       END IF
    ELSE
      CALL abort(&
__STAMP__&
,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
    END IF
  END DO

! we want always warnings to know if the emission has failed. if a timedisc does not require this, this
! timedisc has to be handled separately
#if USE_MPI
  mySumOfRemovedParticles=0
  IF(nChunksTemp.EQ.1) THEN
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,PartfoundInProc(1,:), ChunkSize, MPI_INTEGER, MPI_MAX &
                                                        , PartMPI%InitGroup(InitGroup)%COMM, IERROR)
    ! loop over all particles and check, if particle is found in my proc
    ! proc with LARGES id gets the particle, all other procs remove the duplicated
    ! particle from their list
    DO i=1,chunkSize
      IF(PartFoundInProc(2,i).GT.-1)THEN ! particle has been previously found by MyRank
        IF(PartFoundInProc(1,i).NE.MyRank)THEN ! particle should not be found by MyRank
          !ParticleIndexNbr = PartFoundInProc(2,i)
          ParticleIndexNbr = PDM%nextFreePosition(PartFoundInProc(2,i) + PDM%CurrentNextFreePosition)
          IF(.NOT.PDM%ParticleInside(ParticleIndexNbr)) WRITE(UNIT_stdOut,*) ' Error in emission in parallel!!'
          PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
          PDM%IsNewPart(ParticleIndexNbr)=.FALSE.
          ! correct number of found particles
          mySumOfRemovedParticles = mySumOfRemovedParticles +1
          ! set update next free position to zero for removed particle
          PDM%nextFreePosition(PartFoundInProc(2,i) + PDM%CurrentNextFreePosition) = 0
          !mySumOfMatchedParticles = mySumOfMatchedParticles -1
        END IF
      END IF
    END DO ! i=1,chunkSize
    DEALLOCATE(PartFoundInProc)
    mySumOfMatchedParticles = mySumOfMatchedParticles - mySumOfRemovedParticles
  END IF

  ! check the sum of the matched particles: did each particle find its "home"-CPU?
  CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER, MPI_SUM &
                                           , PartMPI%InitGroup(InitGroup)%COMM, IERROR)
#else
  ! im seriellen Fall kommen alle Partikel auf einen CPU,
  ! daher ist PIC%maxParticleNumber die harte Grenze
  sumOfMatchedParticles = mySumOfMatchedParticles
#endif

#if USE_MPI
  IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
#endif
      ! add number of matching error to particle emission to fit
      ! number of added particles
      Species(FractNbr)%Init(iInit)%InsertedParticleMisMatch = nbrOfParticle  - sumOfMatchedParticles
      IF (nbrOfParticle .GT. sumOfMatchedParticles) THEN
          SWRITE(UNIT_StdOut,'(A)')'WARNING in ParticleEmission_parallel:'
          SWRITE(UNIT_StdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
        SWRITE(UNIT_StdOut,'(A,I0,A)')'matched only ', sumOfMatchedParticles, ' particles'
        SWRITE(UNIT_StdOut,'(A,I0,A)')'when ', NbrOfParticle, ' particles were required!'
      ELSE IF (nbrOfParticle .LT. sumOfMatchedParticles) THEN
            SWRITE(UNIT_StdOut,'(A)')'ERROR in ParticleEmission_parallel:'
            SWRITE(UNIT_StdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
            SWRITE(UNIT_StdOut,'(A,I8,A)')'matched ', sumOfMatchedParticles, ' particles'
            SWRITE(UNIT_StdOut,'(A,I8,A)')'when ', NbrOfParticle, ' particles were required!'
      ELSE IF (nbrOfParticle .EQ. sumOfMatchedParticles) THEN
        !  WRITE(UNIT_stdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
        !  WRITE(UNIT_stdOut,'(A,I0,A)')'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
      END IF
#if USE_MPI
  END IF ! PartMPI%iProc.EQ.0
#endif

  ! Return the *local* NbrOfParticle so that the following Routines only fill in
  ! the values for the local particles
#if USE_MPI
  NbrOfParticle = mySumOfMatchedParticles + mySumOfRemovedParticles
#else
  NbrOfParticle = mySumOfMatchedParticles
#endif

  DEALLOCATE( particle_positions, STAT=allocStat )
  IF (allocStat .NE. 0) THEN
    CALL abort(&
__STAMP__&
,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
  END IF
#if USE_MPI
END IF ! mode 1/2
#endif

END SUBROUTINE SetParticlePosition


SUBROUTINE SetParticleVelocity(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Particle_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Tracking_Vars,  ONLY: DoRefMapping
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_PICInterpolation_Vars,   ONLY: externalField, FieldAtParticle, InterpolationType
USE MOD_Eval_xyz,                ONLY: TensorProductInterpolation, GetPositionInRefElem, EvaluateFieldAtPhysPos
USE MOD_PIC_Vars
!IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr
INTEGER,INTENT(IN)               :: iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: PositionNbr,i,iElem,j
REAL                             :: field(PP_nVar)
REAL                             :: Radius(3), n_vec(3), tan_vec(3), Velo1, Velo2
REAL                             :: Vec3D(3), RandVal(3), Vec1D
INTEGER                          :: Rotation
REAL                             :: Velosq
REAL                             :: RandN_save
LOGICAL                          :: RandN_in_Mem
CHARACTER(30)                    :: velocityDistribution             ! specifying keyword for velocity distribution
REAL                             :: RadiusIC                         ! Radius for IC circle
REAL                             :: NormalIC(3)                      ! Normal / Orientation of circle
REAL                             :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
REAL                             :: VeloIC                           ! velocity for inital Data
REAL                             :: VeloTurbIC                       ! turbulent fluctuation for velocity for initial Data
REAL                             :: VeloIC2                          ! square of velocity for inital Data
REAL                             :: VeloVecIC(3)                     ! normalized velocity vector
REAL                             :: Alpha                            ! WaveNumber for sin-deviation initiation.
REAL                             :: VelocitySpread                   ! widening of init velocity
REAL                             :: vMag2                            ! magnitude of velocity
!===================================================================================================================================
! Abort if we don't have any/too many particles
IF(NbrOfParticle.lt.1) RETURN
IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
    CALL abort(&
    __STAMP__&
    ,'NbrOfParticle > PIC%maxParticleNumber!')
END IF
RandN_in_Mem=.FALSE.

velocityDistribution  = Species(FractNbr)%Init(iInit)%velocityDistribution
VeloVecIC             = Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
VeloIC                = Species(FractNbr)%Init(iInit)%VeloIC
VeloTurbIC            = Species(FractNbr)%Init(iInit)%VeloTurbIC
BasePointIC           = Species(FractNbr)%Init(iInit)%BasePointIC(1:3)
NormalIC              = Species(FractNbr)%Init(iInit)%NormalIC(1:3)
RadiusIC              = Species(FractNbr)%Init(iInit)%RadiusIC
Alpha                 = Species(FractNbr)%Init(iInit)%alpha

SELECT CASE(TRIM(velocityDistribution))

CASE('tangential_constant')
    Rotation       = Species(FractNbr)%Init(iInit)%Rotation
    VelocitySpread = Species(FractNbr)%Init(iInit)%VelocitySpread
    IF(VelocitySpread.GT.0)THEN
        IF(Species(FractNbr)%Init(iInit)%VelocitySpreadMethod.EQ.0)THEN
            ! sigma of normal Distribution, Kostas proposal
            VelocitySpread = VelocitySpread * VeloIC   !/(2.*SQRT(2.*LOG(10.)))
        ELSE IF(Species(FractNbr)%Init(iInit)%VelocitySpreadMethod.EQ.1)THEN
            ! sigma is defined by changing the width of the distribution function at 10% of its maxima
            ! the input value is the spread in percent, hence, 5% => v = v +- 0.05*v at 10% of maximum value
            ! width of the velocity spread, deltaV:
            VelocitySpread = 2.0*VelocitySpread * VeloIC
            ! computing the corresponding sigma
            VelocitySpread = VelocitySpread / (2.*SQRT(2.*LOG(10.)))
        ELSE
            CALL abort(&
            __STAMP__&
            ,' This method for the velocity spread is not implemented.')
        END IF
        IF(alpha.GT.0) THEN
            vMag2 = (1.0+1./(alpha*alpha)) * VeloIC*VeloIC
        ELSE
            vMag2 = VeloIC*VeloIC
        END IF
    END IF
    VeloIC2        = VeloIC*VeloIC
END SELECT

SELECT CASE(TRIM(velocityDistribution))
!===================================================================================================================================
! Emission with random distribution (zero mean, VeloIC standard deviation)
!===================================================================================================================================
CASE('random')
    i = 1
    DO WHILE (i .le. NbrOfParticle)
        PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        IF (PositionNbr .ne. 0) THEN
            CALL RANDOM_NUMBER(RandVal)
            RandVal(:) = RandVal(:) - 0.5
            RandVal(:) = RandVal(:)/SQRT(RandVal(1)**2+RandVal(2)**2+RandVal(3)**2)
            PartState(PositionNbr,4:6) = RandVal(1:3) * VeloIC
            ! New particles have not been reflected
            PartReflCount(PositionNbr) = 0
        END IF
        i = i + 1
  END DO

!===================================================================================================================================
! Emission with constant velocity
!===================================================================================================================================
CASE('constant')
    i = 1
    DO WHILE (i .le. NbrOfParticle)
        PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        IF (PositionNbr .ne. 0) THEN
            PartState(PositionNbr,4:6) = VeloVecIC(1:3) * VeloIC
            ! New particles have not been reflected
            PartReflCount(PositionNbr) = 0
        END IF
        i = i + 1
    END DO

!===================================================================================================================================
! Emission with constant velocity plus Gaussian random fluctuations
!>  Dehbi, A., "A CFD model for particle dispersion in turbulent boundary layer flow", 2006
!===================================================================================================================================
CASE('constant_turbulent')
    i = 1
    DO WHILE (i .le. NbrOfParticle)
        PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        IF (PositionNbr .ne. 0) THEN
            ! Get Gaussian distribution for 2D/3D
            RandVal(:) = 0.
            DO j = 1,PP_dim
                RandVal(j) = RandNormal()
            END DO
            ! Superposition mean and turbulent velocity components
            PartState(PositionNbr,4:6) = VeloVecIC(1:3) * VeloIC + RandVal(1:3) * VeloTurbIC * VeloIC
            ! New particles have not been reflected
            PartReflCount(PositionNbr) = 0
        END IF
        i = i + 1
    END DO

!===================================================================================================================================
!
!===================================================================================================================================
CASE('radial_constant')
    i = 1
    DO WHILE (i .le. NbrOfParticle)
        PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        IF (PositionNbr .ne. 0) THEN
            Radius(1:3) = PartState(PositionNbr,1:3) - BasePointIC(1:3)
            ! Unity radius
            ! Radius(1:3) = Radius(1:3) / RadiusIC
            Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
            PartState(PositionNbr,4:6) = Radius(1:3) * VeloIC
            ! New particles have not been reflected
            PartReflCount(PositionNbr) = 0
        END IF
        i = i + 1
    END DO

!===================================================================================================================================
!
!===================================================================================================================================
CASE('tangential_constant')
    i = 1
    DO WHILE (i .le. NbrOfParticle)
        PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        IF (PositionNbr .ne. 0) THEN
            Radius(1:3) = PartState(PositionNbr,1:3) - BasePointIC(1:3)
            ! Normal Vector of circle
            n_vec(1:3) = NormalIC(1:3)

            ! Unity radius
            Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
            !  Vector Product rxn
            tan_vec(1) = Radius(2)*n_vec(3) - Radius(3)*n_vec(2)
            tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3)
            tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)

            IF(VelocitySpread.GT.0.)THEN
                IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
                    Vec1D = RandN_save
                    RandN_in_Mem=.FALSE.
                ELSE
                    CALL RANDOM_NUMBER(RandVal)
                    Velo1 = 2.0*RandVal(1)-1.0
                    Velo2 = 2.0*RandVal(2)-1.0
                    Velosq= Velo1**2+Velo2**2
                    DO WHILE ((Velosq.LE.0).OR.(Velosq.GE.1))
                        CALL RANDOM_NUMBER(RandVal)
                        Velo1 = 2.0*RandVal(1)-1.0
                        Velo2 = 2.0*RandVal(2)-1.0
                        Velosq= Velo1**2+Velo2**2
                    END DO
                Vec1D = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
                RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
                RandN_in_Mem=.TRUE.
                END IF
                ! velocity spread of tangential velocity
                IF(Rotation.EQ.1)THEN
                    Vec3D  = tan_vec(1:3) * (VeloIC+Vec1D*VelocitySpread)
                ELSE
                    Vec3D = -tan_vec(1:3) * (VeloIC+Vec1D*VelocitySpread)
                END IF
                ! compute axial velocity
                Vec1D = vMag2  - DOT_PRODUCT(Vec3D,Vec3D)
                IF(Vec1D.LT.0) CALL abort(&
                               __STAMP__&
                              ,' Error in set velocity!',PositionNbr)
                Vec1D=SQRT(Vec1D)
                PartState(PositionNbr,4:6) = Vec3D+n_vec(1:3) * Vec1D
            ELSE ! no velocity spread
                ! If Gyrotron resonator: Add velocity in normal direction!
                IF (Alpha .gt. 0.) THEN
                    n_vec = n_vec * ( 1 / Alpha )
                ELSE
                    n_vec = 0
                END IF
                !  And finally the velocities
                IF(Rotation.EQ.1)THEN
                    PartState(PositionNbr,4:6) = tan_vec(1:3) * VeloIC + n_vec(1:3) * VeloIC
                ELSE
                PartState(PositionNbr,4:6) = -tan_vec(1:3) * VeloIC + n_vec(1:3) * VeloIC
                END IF
            END IF
        END IF
        i = i + 1
    END DO

!===================================================================================================================================
! Emission with local fluid velocity. Currently not usable for MPI
!===================================================================================================================================
CASE('fluid')
    ! null field vector
#if USE_MPI
    CALL abort(&
         __STAMP__&
         , 'ERROR: Fluid velocity distribution unstable for MPI case!')
#endif

    field=0.

    FieldAtParticle(:,:) = 0.
    FieldAtParticle(:,1) = externalField(1)
    FieldAtParticle(:,2) = externalField(2)
    FieldAtParticle(:,3) = externalField(3)
    FieldAtParticle(:,4) = externalField(4)
    FieldAtParticle(:,5) = externalField(5)

    DO i = 1,NbrOfParticle
        SELECT CASE(TRIM(InterpolationType))
            CASE('particle_position')
                ! particles have already been mapped in deposition, other eval routine used
                IF (.NOT.DoRefMapping) THEN
                    CALL TensorProductInterpolation(PartState(i,1:3),PartPosRef(1:3,i),iElem)
                END IF
                IPWRITE(*,*) 'min U', minval(U(1:5,:,:,:,iElem)), 'max U', maxval(U(1:5,:,:,:,iElem))
                !--- evaluate at Particle position
                CALL EvaluateFieldAtPhysPos(PartPosRef(1:3,i),PP_nVar,PP_N,U(1:PP_nVar,:,:,:,iElem),field(1:PP_nVar),iElem)
                FieldAtParticle(i,1:PP_nVar) = FieldAtParticle(i,1:PP_nVar) + field(1:PP_nVar)
            CASE DEFAULT
                CALL abort(&
                __STAMP__&
                , 'ERROR: Unknown InterpolationType!')
         END SELECT

        PartState(i,4) = FieldAtParticle(i,2)/FieldAtParticle(i,1)
        PartState(i,5) = FieldAtParticle(i,3)/FieldAtParticle(i,1)
        PartState(i,6) = FieldAtParticle(i,4)/FieldAtParticle(i,1)

        ! New particles have not been reflected
        PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        PartReflCount(PositionNbr) = 0
   END DO

CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'Wrong particle velocity distribution!')

END SELECT

END SUBROUTINE SetParticleVelocity


SUBROUTINE SetParticleMass(FractNbr,NbrOfParticle)
!===================================================================================================================================
! And partilces mass and charge
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,    ONLY : PDM, PartSpecies
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: i,PositionNbr
!===================================================================================================================================

IF(NbrOfParticle.gt.PDM%maxParticleNumber) THEN
    NbrOfParticle = PDM%maxParticleNumber
END IF
i = 1
DO WHILE (i .le. NbrOfParticle)
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .ne. 0) THEN
        PartSpecies(PositionNbr) = FractNbr
    END IF
    i = i + 1
END DO

END SUBROUTINE SetParticleMass


SUBROUTINE InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
!===================================================================================================================================
! Subroutine for checking if calculated particle position would be inside user-defined ExcludeRegion (cuboid or cylinder)
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY : abort
USE MOD_Particle_Vars,          ONLY : Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr, iInit
REAL,INTENT(IN)                  :: Particle_pos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: insideExcludeRegion
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: VecExclude(3), DistExclude
INTEGER                          :: iExclude
!===================================================================================================================================
insideExcludeRegion=.FALSE.

DO iExclude=1,Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions
    VecExclude = Particle_pos - Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC

    SELECT CASE (TRIM(Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC))

        CASE ('cuboid')
            !--check normal direction
            DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
                        + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
                        + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
            IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC) &
            .AND.(DistExclude .GE. 0.) ) THEN
                insideExcludeRegion = .TRUE.
            ELSE
                insideExcludeRegion = .FALSE.
                CYCLE
            END IF

            !--check BV1 direction
            DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
                        + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
                        + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)
            IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1)**2) &
            .AND.(DistExclude .GE. 0.) ) THEN
                insideExcludeRegion = .TRUE.
            ELSE
                insideExcludeRegion = .FALSE.
                CYCLE
            END IF

            !--check BV2 direction
            DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
                        + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
                        + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
            IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2)**2) &
            .AND.(DistExclude .GE. 0.) ) THEN
                insideExcludeRegion = .TRUE.
                RETURN !particle is inside current ExcludeRegion based an all dimensions
            ELSE
                insideExcludeRegion = .FALSE.
                CYCLE
            END IF

        CASE ('cylinder')
            !--check normal direction
            DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
                        + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
                        + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
            IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC) &
            .AND.(DistExclude .GE. 0.) ) THEN
                insideExcludeRegion = .TRUE.
            ELSE
                insideExcludeRegion = .FALSE.
                CYCLE
            END IF

            !--check radial direction
            DistExclude = SQRT( VecExclude(1)**2 + VecExclude(2)**2 + VecExclude(3)**2 - DistExclude**2 )
            IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC) &
            .AND.(DistExclude .GE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC) ) THEN
                insideExcludeRegion = .TRUE.
                RETURN !particle is inside current ExcludeRegion based an all dimensions
            ELSE
                insideExcludeRegion = .FALSE.
                CYCLE
            END IF

        CASE DEFAULT
            CALL abort(&
            __STAMP__&
            ,'wrong SpaceIC for ExcludeRegion!')

    END SELECT
END DO

END SUBROUTINE InsideExcludeRegionCheck


SUBROUTINE SamplePoissonDistri(RealTarget,IntSample,Flag_opt)
!===================================================================================================================================
! Sample IntSample from Poisson-Distri around RealTarget (if Flag present it will be turned off at sample limit, otherwise abort)
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: RealTarget
LOGICAL,INTENT(INOUT),OPTIONAL :: Flag_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: IntSample
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL         :: Flag
INTEGER         :: Npois
REAL            :: Tpois, RandVal1
!===================================================================================================================================

IF (PRESENT(Flag_opt)) THEN
  Flag=Flag_opt
ELSE
  Flag=.FALSE.
END IF

Npois=0
Tpois=1.0
CALL RANDOM_NUMBER(RandVal1)
DO
  Tpois=RandVal1*Tpois
  IF (Tpois.LT.TINY(Tpois)) THEN
    IF (Flag) THEN !Turn off Poisson Sampling and "sample" by random-rounding
      IPWRITE(*,*)'WARNING: target is too large for poisson sampling: switching now to Random rounding...'
      IntSample = INT(RealTarget + RandVal1)
      Flag = .FALSE.
      EXIT
    ELSE !Turning off not allowed: abort (RealTarget must be decreased ot PoissonSampling turned off manually)
      CALL abort(&
__STAMP__&
,'ERROR in SamplePoissonDistri: RealTarget (e.g. flux) is too large for poisson sampling!')
    END IF
  END IF
  IF (Tpois.GT.EXP(-RealTarget)) THEN
    Npois=Npois+1
    CALL RANDOM_NUMBER(RandVal1)
  ELSE
    IntSample = Npois
    EXIT
  END IF
END DO

END SUBROUTINE SamplePoissonDistri


END MODULE MOD_part_emission
