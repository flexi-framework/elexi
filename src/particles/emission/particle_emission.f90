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
#include "particle.h"

!===================================================================================================================================
! Module for particle emission
!===================================================================================================================================
MODULE MOD_Part_Emission
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitializeParticleEmission
  MODULE PROCEDURE InitializeParticleEmission
END INTERFACE

INTERFACE ParticleInserting
  MODULE PROCEDURE ParticleInserting
END INTERFACE

PUBLIC :: InitializeParticleEmission
PUBLIC :: ParticleInserting
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeParticleEmission()
!===================================================================================================================================
! Initialize particles / Insert initial particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,      ONLY: Species,nSpecies,PDM,PEM
USE MOD_Part_Emission_Tools,ONLY: SetParticleMass
USE MOD_Part_Pos_and_Velo,  ONLY: SetParticlePosition,SetParticleVelocity
USE MOD_Part_Tools,         ONLY: UpdateNextFreePosition
USE MOD_Restart_Vars,       ONLY: DoRestart
#if USE_MPI
USE MOD_Particle_MPI_Vars,  ONLY: PartMPI
#endif /* MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i,NbrOfParticle,iInit
INTEGER               :: insertParticles
INTEGER               :: ParticleVecLengthGlob
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INITIAL PARTICLE INSERTING... '

! Update next free position in particle array to get insertion position
CALL UpdateNextFreePosition()

! FLEXI gained restart capability from pure fluid solution. Do initial particle inserting only if we are not restarting or if there
! are no particles present in the current file. Check globally, otherwise some procs might block
#if USE_MPI
CALL MPI_ALLREDUCE(PDM%ParticleVecLength,ParticleVecLengthGlob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
#else
ParticleVecLengthGlob = PDM%ParticleVecLength
#endif /* MPI */

IF (.NOT.DoRestart .OR. (ParticleVecLengthGlob.EQ.0)) THEN
  ! for the case of particle insertion per time, the inserted particle number for the current time must
  ! be updated. Otherwise, at the first timestep after restart, these particles will be inserted again
  ! Do insanity check of max. particle number compared to the number that is to be inserted for certain insertion types
  insertParticles = 0

  ! Get number of particles to be inserted. Divide it between MPI procs if required.
  DO i=1,nSpecies
      DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
          IF ((TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cuboid') .OR.(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
#if USE_MPI
            insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs)
#else
            insertParticles = insertParticles + INT(Species(i)%Init(iInit)%initialParticleNumber)
#endif
          END IF
      END DO
  END DO

  ! Check if requested to insert more particles than the particle array can hold
  IF (insertParticles.GT.PDM%maxParticleNumber) THEN
    IPWRITE(UNIT_stdOut,'(A40,I0)') ' Maximum particle number : ',PDM%maxParticleNumber
    IPWRITE(UNIT_stdOut,'(A40,I0)') ' To be inserted particles: ',insertParticles
    CALL abort(__STAMP__,'Number of to be inserted particles per init-proc exceeds max. particle number! ')
  END IF

  DO i = 1,nSpecies
    DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
      ! no special emissiontype to be used
      IF (Species(i)%Init(iInit)%UseForInit) THEN
        ! initial particle number too large to handle
        IF(Species(i)%Init(iInit)%initialParticleNumber.GT.HUGE(1)) &
          CALL abort(__STAMP__,' Integer for initialParticleNumber too large!')

        NbrOfParticle = INT(Species(i)%Init(iInit)%initialParticleNumber,4)
        SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle position for species ',i,' ... '
        CALL SetParticlePosition(i,iInit,NbrOfParticle)
        ! give the particles their correct velocity
        SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle velocities for species ',i,' ... '
        CALL SetParticleVelocity(i,iInit,NbrOfParticle,1)
        ! give the particles their correct (species) mass
        SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle mass for species ',i,' ... '
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

SWRITE(UNIT_stdOut,'(A)') ' INITIAL PARTICLE INSERTING DONE '
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitializeParticleEmission


SUBROUTINE ParticleInserting()
!===================================================================================================================================
! Particle Inserting
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /* MPI*/
USE MOD_Globals
USE MOD_Restart_Vars           ,ONLY: RestartTime
USE MOD_Part_Emission_Tools    ,ONLY: SamplePoissonDistri
USE MOD_Part_Emission_Tools    ,ONLY: SetParticleMass
USE MOD_Part_Pos_and_Velo      ,ONLY: SetParticlePosition,SetParticleVelocity
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPartBalance,nPartIn,PartEkinIn
USE MOD_Particle_Globals       ,ONLY: ALMOSTEQUAL
USE MOD_Particle_Restart_Vars  ,ONLY: PartDataExists
USE MOD_Particle_Timedisc_Vars ,ONLY: RKdtFrac,RKdtFracTotal
USE MOD_Particle_Vars          ,ONLY: Species,nSpecies,PartSpecies,PDM,DelayTime,DoPoissonRounding,DoTimeDepInflow
USE MOD_Timedisc_Vars          ,ONLY: dt,t,RKc,nRKStages,currentStage
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                          :: i,iPart,PositionNbr,iInit,IntSample
INTEGER                , SAVE    :: NbrOfParticle=0
INTEGER(KIND=8)                  :: inserted_Particle_iter,inserted_Particle_time
INTEGER(KIND=8)                  :: inserted_Particle_diff
REAL                             :: PartIns, RandVal1
REAL                             :: RiseFactor, RiseTime
#if USE_MPI
INTEGER                          :: InitGroup
#endif
!===================================================================================================================================

! We are calling particleInserting differently in FLEXI. Calculate RKdtFrac here using Timedisc_Vars
!> Basically a fancy way to reconstruct tStage for the next Runge-Kutta stage
IF (currentStage.EQ.1) THEN
  RKdtFrac      = RKC(2)
  RKdtFracTotal = RKdtFrac
ELSE
  IF (currentStage.NE.nRKStages) THEN
    RKdtFrac      = RKC(currentStage+1) - RKC(currentStage)
    RKdtFracTotal = RKdtFracTotal + RKdtFrac
  ELSE
    RKdtFrac      = 1. - RKC(nRKStages)
    RKdtFracTotal = 1.
  END IF
END IF

!---  Emission at time step (initial emission see particle_init.f90: InitializeParticleEmission)
DO i = 1,nSpecies
  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
    ! Reset the number of particles per species AND init region
    NbrOfParticle = 0

    ! species to be used for init
    IF (Species(i)%Init(iInit)%UseForEmission) THEN

      SELECT CASE(Species(i)%Init(iInit)%ParticleEmissionType)
        ! Emission Type: Particles per second
        CASE(1)
          IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
            ! requested particles during time-slab
            PartIns = Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac * 1/Species(i)%Init(iInit)%ParticleEmissionTime
            ! integer number of particles to be inserted
            inserted_Particle_iter = INT(PartIns,8)

            ! Attention: FLEXI uses a different definition for elapsed time due to restart capability!
            IF ((RestartTime.GT.0) .AND. (.NOT.PartDataExists)) THEN
              ! total number of emitted particles over simulation
              PartIns = Species(i)%Init(iInit)%ParticleEmission * (t-RestartTime + dt*RKdtFracTotal) * 1/Species(i)%Init(iInit)%ParticleEmissionTime
            ELSE
              PartIns = Species(i)%Init(iInit)%ParticleEmission * (t             + dt*RKdtFracTotal) * 1/Species(i)%Init(iInit)%ParticleEmissionTime
            END IF

            !-- random-round the inserted_Particle_time for preventing periodicity
            CALL RANDOM_NUMBER(RandVal1)
            !-- repeat if inserting multiple particles
            IF (inserted_Particle_iter.GE.1) THEN
              CALL RANDOM_NUMBER(RandVal1)
              inserted_Particle_time = INT(PartIns + RandVal1,8) ! adds up to ONE
            ! needed, since InsertedParticleSurplus can increase
            ! and _iter > 1 needs to be possible for preventing periodicity
            ELSE IF ((inserted_Particle_iter.GE.0).AND.(inserted_Particle_iter.LT.1)) THEN
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
            inserted_Particle_diff = inserted_Particle_time - Species(i)%Init(iInit)%InsertedParticle        &
                                   - inserted_Particle_iter - Species(i)%Init(iInit)%InsertedParticleSurplus &
                                   + Species(i)%Init(iInit)%InsertedParticleMisMatch

            Species(i)%Init(iInit)%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0))
            NbrOfParticle = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)

          ELSE IF (DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
            ! linear rise of inflow
            RiseTime = Species(i)%Init(iInit)%InflowRiseTime

            ! ramp up the particle insertion by increasing the RiseFactor depending on time
            IF (RiseTime.GT.0.) THEN
              IF (t-DelayTime.LT.RiseTime) THEN
                  RiseFactor = (t-DelayTime)/RiseTime
              ELSE
                  RiseFactor = 1.
              END IF
            ELSE
              RiseFactor = 1.
            END IF

            ! emitted particles during time-slab
            PartIns = Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac * RiseFactor
            CALL RANDOM_NUMBER(RandVal1)

            IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
              IPWRITE(UNIT_StdOut,'(A)') ' WARNING: target is too large for poisson sampling: switching now to Random rounding...'
              NbrOfParticle     = INT(PartIns + RandVal1)
              DoPoissonRounding = .FALSE.
            ELSE
              !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects)
              !> [Tysanner and Garcia 2004]
              CALL SamplePoissonDistri( PartIns , NbrOfParticle , DoPoissonRounding)
            END IF

          ! DoTimeDepInflow
          ELSE
            ! linear rise of inflow
            RiseTime = Species(i)%Init(iInit)%InflowRiseTime
            IF (RiseTime.GT.0.) THEN
              IF (t-DelayTime .LT. RiseTime) THEN
                RiseFactor = (t-DelayTime)/RiseTime
              ELSE
                RiseFactor = 1.
              END IF
            ELSE
               RiseFactor = 1.
            END IF

            ! emitted particles during time-slab
            PartIns = Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac * RiseFactor &
                    + Species(i)%Init(iInit)%InsertedParticleMisMatch
            CALL RANDOM_NUMBER(RandVal1)
            NbrOfParticle = INT(PartIns + RandVal1)
          END IF
#if USE_MPI
          ! communicate number of particles with all procs in the same init group
          InitGroup = Species(i)%Init(iInit)%InitCOMM
          IF (PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL) THEN
            ! only procs which are part of group take part in the communication
            ! NbrOfParticle based on RandVals!
            CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
          ELSE
            NbrOfParticle = 0
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

        ! Emission Type: Particles at time x
        CASE(3)
          ! insert in last stage only, so that no reconstruction is necessary and number/iter matches
          IF (RKdtFracTotal .EQ. 1..AND.(t-RestartTime).GT.0) THEN
!            print *, 't', 'MOD', t-RestartTime,dt, MOD((t-RestartTime),Species(i)%Init(iInit)%ParticleEmissionTime)
            IF (MOD((t-RestartTime),Species(i)%Init(iInit)%ParticleEmissionTime).LT.dt) THEN
              NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission)
            ELSE
              NbrOfParticle = 0
            END IF
          ELSE
            NbrOfParticle = 0
          END IF

        CASE DEFAULT
          CALL ABORT(__STAMP__,'Unknown particle emission type')
          ! Line below not required anymore, but might want to change back to silent fail later
          ! NbrOfParticle = 0
      END SELECT

      CALL SetParticlePosition(i,iInit,NbrOfParticle)
      CALL SetParticleVelocity(i,iInit,NbrOfParticle,1)
      CALL SetParticleMass(    i      ,NbrOfParticle)

      ! instead of UpdateNextfreePosition we update the particleVecLength only and doing it later, after CalcPartBalance
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
      PDM%ParticleVecLength       = PDM%ParticleVecLength       + NbrOfParticle
      !CALL UpdateNextFreePosition()

      ! Compute number of input particles and energy
      IF (CalcPartBalance) THEN
        ! Alter history, dirty hack for balance calculation
        PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition - NbrOfParticle
        IF(NbrOfParticle.GT.0) THEN
          nPartIn(i) = nPartIn(i) + NBrofParticle
          DO iPart = 1,NbrOfparticle
              PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
              IF (PositionNbr .NE. 0) PartEkinIn(PartSpecies(PositionNbr)) = PartEkinIn(PartSpecies(PositionNbr)) + &
                                                                             CalcEkinPart(PositionNbr)
          END DO ! iPart
        END IF
        ! alter history, dirty hack for balance calculation
        PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
      END IF ! CalcPartBalance

      ! Complete check if all particles were emitted successfully
#if USE_MPI
      InitGroup = Species(i)%Init(iInit)%InitCOMM
      IF (PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL .AND. Species(i)%Init(iInit)%sumOfRequestedParticles.GT.0) THEN
        CALL MPI_WAIT(PartMPI%InitGroup(InitGroup)%Request, MPI_STATUS_IGNORE, iError)

        IF (PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
#endif
          ! add number of matching error to particle emission to fit
          ! number of added particles
          Species(i)%Init(iInit)%InsertedParticleMisMatch = Species(i)%Init(iInit)%sumOfRequestedParticles - Species(i)%Init(iInit)%sumOfMatchedParticles
          IF (Species(i)%Init(iInit)%sumOfRequestedParticles.GT. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
            SWRITE(UNIT_StdOut,'(A)')      'WARNING in ParticleEmission_parallel:'
            SWRITE(UNIT_StdOut,'(A,I0)')   'Fraction Nbr: '  , i
            SWRITE(UNIT_StdOut,'(A,I0,A)') 'matched only '   , Species(i)%Init(iInit)%sumOfMatchedParticles  , ' particles'
            SWRITE(UNIT_StdOut,'(A,I0,A)') 'when '           , Species(i)%Init(iInit)%sumOfRequestedParticles, ' particles were required!'
          ELSE IF (Species(i)%Init(iInit)%sumOfRequestedParticles.LT. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
            SWRITE(UNIT_StdOut,'(A)')      'ERROR in ParticleEmission_parallel:'
            SWRITE(UNIT_StdOut,'(A,I0)')   'Fraction Nbr: '  , i
            SWRITE(UNIT_StdOut,'(A,I0,A)') 'matched '        , Species(i)%Init(iInit)%sumOfMatchedParticles  , ' particles'
            SWRITE(UNIT_StdOut,'(A,I0,A)') 'when '           , Species(i)%Init(iInit)%sumOfRequestedParticles, ' particles were required!'
          ! ELSE IF (nbrOfParticle .EQ. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
          !  WRITE(UNIT_stdOut,'(A,I0)')   'Fraction Nbr: '  , FractNbr
          !  WRITE(UNIT_stdOut,'(A,I0,A)') 'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
          END IF
#if USE_MPI
        END IF ! PartMPI%InitGroup(InitGroup)%MPIRoot
      END IF ! PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL
#endif
    END IF ! UseForEmission
  END DO ! iInit
END DO ! i=1,nSpecies

END SUBROUTINE ParticleInserting



END MODULE MOD_Part_Emission
