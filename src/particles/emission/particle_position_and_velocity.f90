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

MODULE MOD_Part_Pos_and_Velo
!===================================================================================================================================
! module for particle emission
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE SetParticlePosition
  MODULE PROCEDURE SetParticlePosition
END INTERFACE

INTERFACE SetParticleVelocity
  MODULE PROCEDURE SetParticleVelocity
END INTERFACE

!===================================================================================================================================
PUBLIC         :: SetParticleVelocity,SetParticlePosition
!===================================================================================================================================
CONTAINS


SUBROUTINE SetParticlePositionCellLocal(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
USE MOD_Part_Emission_Tools    ,ONLY: IntegerDivide,SetCellLocalParticlePosition
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: chunkSize
LOGICAL                                  :: DoExactPartNumInsert
#if USE_MPI
INTEGER                                  :: InitGroup
REAL,ALLOCATABLE                         :: ProcMeshVol(:)
INTEGER,ALLOCATABLE                      :: ProcNbrOfParticle(:)
#endif
!===================================================================================================================================/

  DoExactPartNumInsert =  .FALSE.
  ! check if particle inserting during simulation or initial inserting and also if via partdensity or exact particle number
  ! nbrOfParticles is set for initial inserting if initialPartNum or partdensity is set in ini
  ! ParticleEmission and Partdensity not working together
  IF (NbrofParticle.EQ.0.AND.(Species(FractNbr)%Init(iInit)%ParticleEmission.EQ.0)) RETURN
  IF ((NbrofParticle.GT.0).AND.(Species(FractNbr)%Init(iInit)%PartDensity.LE.0.)) THEN
    DoExactPartNumInsert =  .TRUE.
  END IF
  chunksize = 0
#if USE_MPI
! emission group communicator
  InitGroup=Species(FractNbr)%Init(iInit)%InitCOMM
  IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
    NbrofParticle=0
    RETURN
  END IF
  IF (PartMPI%InitGroup(InitGroup)%nProcs.GT.1) THEN
    IF (DoExactPartNumInsert) THEN !###$ ToDo
      IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
        ALLOCATE(ProcMeshVol(0:PartMPI%InitGroup(InitGroup)%nProcs-1))
        ALLOCATE(ProcNbrOfParticle(0:PartMPI%InitGroup(InitGroup)%nProcs-1))
        ProcMeshVol=0.
        ProcNbrOfParticle=0
      ELSE ! to reduce global memory allocation if a lot of procs are used
        ALLOCATE(ProcMeshVol(1))
        ALLOCATE(ProcNbrOfParticle(1))
        ProcMeshVol=0.
        ProcNbrOfParticle=0
      END IF !InitGroup%MPIroot
      CALL MPI_GATHER(LocalVolume,1,MPI_DOUBLE_PRECISION &
                     ,ProcMeshVol,1,MPI_DOUBLE_PRECISION,0,PartMPI%InitGroup(InitGroup)%COMM,iError)
      IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
        CALL IntegerDivide(NbrOfParticle,PartMPI%InitGroup(InitGroup)%nProcs,ProcMeshVol,ProcNbrOfParticle)
      END IF
      CALL MPI_SCATTER(ProcNbrOfParticle, 1, MPI_INTEGER, chunksize, 1, MPI_INTEGER, 0, PartMPI%InitGroup(InitGroup)%COMM, IERROR)
      SDEALLOCATE(ProcMeshVol)
      SDEALLOCATE(ProcNbrOfParticle)
    END IF
  ELSE
    chunksize = NbrOfParticle
  END IF
#else
  IF (DoExactPartNumInsert) chunksize = NbrOfParticle
#endif /*USE_MPI*/
  IF ((chunksize.GT.0).OR.(Species(FractNbr)%Init(iInit)%PartDensity.GT.0.)) THEN
    CALL SetCellLocalParticlePosition(chunkSize,FractNbr,iInit,DoExactPartNumInsert)
  END IF
  NbrOfParticle = chunksize

END SUBROUTINE SetParticlePositionCellLocal


SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Globals,      ONLY: PI
USE MOD_Particle_Vars,         ONLY: Species,PDM,PartState
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
USE MOD_Timedisc_Vars,         ONLY: dt
USE MOD_Particle_TimeDisc_Vars,ONLY: RKdtFrac
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_Part_Emission_Tools    ,ONLY: IntegerDivide,SetCellLocalParticlePosition,SetParticlePositionPoint
USE MOD_Part_Emission_Tools    ,ONLY: SetParticlePositionEquidistLine, SetParticlePositionLine, SetParticlePositionDisk
USE MOD_Part_Emission_Tools    ,ONLY: SetParticlePositionCuboidCylinder,SetParticlePositionCircle
USE MOD_Part_Emission_Tools    ,ONLY: SetParticlePositionSphere,SetParticlePositionSinDeviation
#if USE_MPI
USE MOD_Particle_MPI_Emission  ,ONLY: SendEmissionParticlesToProcs
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                         :: particle_positions(:)
INTEGER                                  :: i,ParticleIndexNbr,allocStat,nChunks, chunkSize
INTEGER                                  :: mySumOfMatchedParticles, sumOfMatchedParticles, DimSend
LOGICAL                                  :: DoExactPartNumInsert
#if USE_MPI
INTEGER                                  :: InitGroup
REAL,ALLOCATABLE                         :: ProcMeshVol(:)
INTEGER,ALLOCATABLE                      :: ProcNbrOfParticle(:)
#endif
!===================================================================================================================================
IF (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
  CALL SetParticlePositionCellLocal(FractNbr,iInit,NbrOfParticle)
  RETURN
END IF
IF ( (NbrOfParticle .LE. 0).AND. (ABS(Species(FractNbr)%Init(iInit)%PartDensity).LE.0.) ) RETURN

! emission group communicator for the current iInit
#if USE_MPI
InitGroup = Species(FractNbr)%Init(iInit)%InitCOMM
IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
  NbrofParticle=0
  RETURN
END IF
#endif /*USE_MPI*/

DimSend      = 3            !save (and send) only positions
nChunks      = 1  
sumOfMatchedParticles   = 0
mySumOfMatchedParticles = 0
chunkSize    = nbrOfParticle
! process myRank=0 generates the complete list of random positions for all emitted particles
#if USE_MPI
IF ( (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle_equidistant'                 ) .OR.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'sin_deviation'                      ) .OR.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle'                             ) .OR.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'line_with_equidistant_distribution' )) THEN
  nChunks = 1
ELSE IF (nbrOfParticle.GT.(PartMPI%InitGroup(InitGroup)%nProcs*10)) THEN
  nChunks = PartMPI%InitGroup(InitGroup)%nProcs
ELSE
  nChunks = 1
END IF
chunkSize = INT(nbrOfParticle/nChunks)
IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
  chunkSize = chunkSize + ( nbrOfParticle - (nChunks*chunkSize) )
END IF
! all proc taking part in particle inserting
IF (PartMPI%InitGroup(InitGroup)%MPIROOT.OR.nChunks.GT.1) THEN
#endif
  ALLOCATE( particle_positions(1:chunkSize*DimSend), STAT=allocStat )
  IF (allocStat .NE. 0) &
    CALL abort(__STAMP__,'ERROR in SetParticlePosition: cannot allocate particle_positions!')

  !------------------SpaceIC-cases: start-----------------------------------------------------------!
  SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
  CASE ('point')
    CALL SetParticlePositionPoint(FractNbr,iInit,chunkSize,particle_positions)
  CASE ('line_with_equidistant_distribution')
    CALL SetParticlePositionEquidistLine(FractNbr,iInit,chunkSize,particle_positions)
  CASE ('line')
    CALL SetParticlePositionLine(FractNbr,iInit,chunkSize,particle_positions)
  CASE ('Gaussian')
    ! TODO
  CASE('disc')
    CALL SetParticlePositionDisk(FractNbr,iInit,chunkSize,particle_positions)
  CASE('circle', 'circle_equidistant')
    CALL SetParticlePositionCircle(FractNbr,iInit,chunkSize,particle_positions)
  CASE('cuboid','cylinder')
    CALL SetParticlePositionCuboidCylinder(FractNbr,iInit,chunkSize,particle_positions)
  CASE('sphere')
    CALL SetParticlePositionSphere(FractNbr,iInit,chunkSize,particle_positions)
  CASE('sin_deviation')
    CALL SetParticlePositionSinDeviation(FractNbr,iInit,chunkSize,particle_positions)
  END SELECT
  !------------------SpaceIC-cases: end-------------------------------------------------------------------------------------------
#if USE_MPI
ELSE !no mpi root, nchunks=1
  chunkSize=0
END IF
! Need to open MPI communication regardless of the chunk number. Make it only dependent on the number of procs
IF (PartMPI%InitGroup(InitGroup)%nProcs.GT.1) THEN
  CALL SendEmissionParticlesToProcs(chunkSize,DimSend,FractNbr,iInit,mySumOfMatchedParticles,particle_positions)
! Finish emission on local proc
ELSE
  mySumOfMatchedParticles = 0
  ParticleIndexNbr        = 1
  DO i=1,chunkSize
    ! Find a free position in the PDM array
    IF ((i.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
      ParticleIndexNbr = PDM%nextFreePosition(mySumOfMatchedParticles + 1 + PDM%CurrentNextFreePosition)
    END IF
    IF (ParticleIndexNbr.NE.0) THEN
      PartState(1:DimSend,ParticleIndexNbr) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+DimSend)
      PDM%ParticleInside( ParticleIndexNbr) = .TRUE.
      CALL LocateParticleInElement(ParticleIndexNbr,doHALO=.FALSE.)
      IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
        mySumOfMatchedParticles = mySumOfMatchedParticles + 1
        PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
      END IF
    ELSE
          CALL abort(&
    __STAMP__&
    ,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
    END IF
  END DO
END IF
! we want always warnings to know if the emission has failed. if a timedisc does not require this, this
! timedisc has to be handled separately
! check the sum of the matched particles: did each particle find its "home"-CPU?
CALL MPI_ALLREDUCE( mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER, MPI_SUM &
                  , PartMPI%InitGroup(InitGroup)%COMM, IERROR)
#else
! in the seriell case, particles are only emitted on the current proc
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
        SWRITE(UNIT_StdOut,'(A,I0,A)')'matched ', sumOfMatchedParticles, ' particles'
        SWRITE(UNIT_StdOut,'(A,I0,A)')'when ', NbrOfParticle, ' particles were required!'
  ELSE IF (nbrOfParticle .EQ. sumOfMatchedParticles) THEN
  !  WRITE(UNIT_stdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
  !  WRITE(UNIT_stdOut,'(A,I0,A)')'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
  END IF
#if USE_MPI
END IF ! PartMPI%iProc.EQ.0
#endif
! Return the *local* NbrOfParticle so that the following Routines only fill in
! the values for the local particles
NbrOfParticle = mySumOfMatchedParticles

IF (chunkSize.GT.0) THEN
  DEALLOCATE(particle_positions, STAT=allocStat)
  IF (allocStat .NE. 0) &
    CALL ABORT(__STAMP__,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
END IF

END SUBROUTINE SetParticlePosition

SUBROUTINE SetParticleVelocity(FractNbr,iInit,NbrOfParticle,init_or_sf)
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Globals,        ONLY: RandNormal
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Eval_xyz,                ONLY: EvaluateFieldAtPhysPos,EvaluateFieldAtRefPos
USE MOD_Particle_Interpolation,      ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars, ONLY: DoInterpolation,FieldAtParticle,externalField
USE MOD_Particle_Tracking_Vars,  ONLY: DoRefMapping
USE MOD_Particle_Vars,           ONLY: PartState,PDM,PEM,Species,PartPosRef,PartReflCount
#if USE_RW
USE MOD_DG_Vars,                 ONLY: UTurb
USE MOD_Restart_Vars,            ONLY: RestartTurb
USE MOD_Equation_Vars,           ONLY: nVarTurb
USE MOD_Particle_Interpolation_Vars, ONLY: TurbFieldAtParticle
#endif /* USE_RW */
!IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr
INTEGER,INTENT(IN)               :: iInit
INTEGER,INTENT(IN)               :: init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: PositionNbr,i,iElem,j
REAL                             :: field(PP_nVar)
REAL                             :: Radius(3)
REAL                             :: RandVal(3)
CHARACTER(30)                    :: velocityDistribution             ! specifying keyword for velocity distribution
REAL                             :: RadiusIC                         ! Radius for IC circle
REAL                             :: NormalIC(3)                      ! Normal / Orientation of circle
REAL                             :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
REAL                             :: VeloIC                           ! velocity for inital Data
REAL                             :: VeloTurbIC                       ! turbulent fluctuation for velocity for initial Data
REAL                             :: VeloVecIC(3)                     ! normalized velocity vector
REAL                             :: Alpha                            ! WaveNumber for sin-deviation initiation.
#if USE_RW
REAL                             :: turbField(nVarTurb)
#endif
!===================================================================================================================================
! Abort if we don't have any/too many particles
IF(NbrOfParticle.LT.1) RETURN
IF(NbrOfParticle.GT.PDM%maxParticleNumber) &
    CALL abort(__STAMP__,'NbrOfParticle > PDM%maxParticleNumber!')

SELECT CASE (init_or_sf)
  CASE(1) !iInit
  ! Collect velocity parameters for current Init
  velocityDistribution  = Species(FractNbr)%Init(iInit)%velocityDistribution
  VeloVecIC             = Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
  VeloIC                = Species(FractNbr)%Init(iInit)%VeloIC
  VeloTurbIC            = Species(FractNbr)%Init(iInit)%VeloTurbIC
  BasePointIC           = Species(FractNbr)%Init(iInit)%BasePointIC(1:3)
  NormalIC              = Species(FractNbr)%Init(iInit)%NormalIC(1:3)
  RadiusIC              = Species(FractNbr)%Init(iInit)%RadiusIC
  Alpha                 = Species(FractNbr)%Init(iInit)%alpha
CASE(2) !SurfaceFlux
  ! TODO
  CALL ABORT(__STAMP__,'surface flux not yet implemented in FLEXI')
  
  IF (TRIM(Species(FractNbr)%Surfaceflux(iInit)%velocityDistribution).EQ.'constant') THEN
    velocityDistribution=Species(FractNbr)%Surfaceflux(iInit)%velocityDistribution
  ELSE
    CALL ABORT(__STAMP__,'only constant velo-distri implemented in SetParticleVelocity for surfaceflux!')
  END IF
  VeloVecIC             = Species(FractNbr)%Surfaceflux(iInit)%VeloVecIC(1:3)
  VeloIC                = Species(FractNbr)%Surfaceflux(iInit)%VeloIC
CASE DEFAULT
  CALL ABORT(__STAMP__,'neither iInit nor Surfaceflux defined as reference!')
END SELECT

! Set velocity according to velocityDistribution
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
      PartState(4:6,PositionNbr) = RandVal(1:3) * VeloIC
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
      PartState(4:6,PositionNbr) = VeloVecIC(1:3) * VeloIC
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
      PartState(4:6,PositionNbr) = VeloVecIC(1:3) * VeloIC + RandVal(1:3) * VeloTurbIC * VeloIC
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
      Radius(1:3) = PartState(1:3,PositionNbr) - BasePointIC(1:3)
      ! Unity radius
      Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
      PartState(4:6,PositionNbr) = Radius(1:3) * VeloIC
      ! New particles have not been reflected
      PartReflCount(PositionNbr) = 0
    END IF
    i = i + 1
  END DO

!===================================================================================================================================
! Emission with local fluid velocity. Currently untested for MPI
!===================================================================================================================================
CASE('fluid')
#if USE_MPI
!  CALL abort(&
!       __STAMP__&
!       , 'ERROR: Fluid velocity distribution unstable for MPI case!')
  SWRITE(*,*) 'Emission with local fluid velocity currently untested for MPI case.'
#endif

  ! null field vector
  field = 0.
  i     = 1

  ! Get background field
  FieldAtParticle(:,:) = 0.
  FieldAtParticle(1,:) = externalField(1)
  FieldAtParticle(2,:) = externalField(2)
  FieldAtParticle(3,:) = externalField(3)
  FieldAtParticle(4,:) = externalField(4)
  FieldAtParticle(5,:) = externalField(5)

  ! Interpolate fluid and field to particle
  DO WHILE (i .LE. NbrOfParticle)
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    ! Valid particle which got recently position
    IF (PositionNbr .NE. 0) THEN

    ! Return if no interpolation is wanted
    IF (.NOT.DoInterpolation) RETURN

      iElem = PEM%Element(PositionNbr)
      IF (.NOT.DoRefMapping) THEN
#if USE_RW
        IF (RestartTurb) THEN
          CALL EvaluateFieldAtPhysPos(PartState(1:3,PositionNbr),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,iElem),field    (1:PP_nVar) ,iElem,PositionNbr &
                                                                             ,UTurb(1:nVarTurb,:,:,:,iElem),turbField(1:nVarTurb))
        ELSE
#endif
          CALL EvaluateFieldAtPhysPos(PartState(1:3,PositionNbr),PP_nVar,PP_N,U(1:PP_nVar,:,:,:,iElem),field    (1:PP_nVar) ,iElem,PositionNbr)
#if USE_RW
        END IF ! RestartTurb
#endif
      ! RefMapping, evaluate in reference space
      ELSE
#if USE_RW
        IF (RestartTurb) THEN
          CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PositionNbr),PP_nVar,PP_N,U    (1:PP_nVar ,:,:,:,iElem),field    (1:PP_nVar)  &
                                                                             ,UTurb(1:nVarTurb,:,:,:,iElem),turbField(1:nVarTurb))
        ELSE
#endif
          CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PositionNbr),PP_nVar,PP_N,U    (1:PP_nVar,:,:,:,iElem),field    (1:PP_nVar))
#if USE_RW
        END IF ! RestartTurb
#endif
      END IF ! RefMapping

      ! Add the interpolated field to the background field
      FieldAtParticle(    1:PP_nVar, PositionNbr) = FieldAtParticle(1:PP_nVar,PositionNbr) + field(1:PP_nVar)
#if USE_RW
      TurbFieldAtParticle(1:nVarTurb,PositionNbr) = turbfield(      1:nVarTurb)
#endif

      ! Calculate velocity from momentum and density
      PartState(4,PositionNbr) = FieldAtParticle(2,PositionNbr)/FieldAtParticle(1,PositionNbr)
      PartState(5,PositionNbr) = FieldAtParticle(3,PositionNbr)/FieldAtParticle(1,PositionNbr)
      PartState(6,PositionNbr) = FieldAtParticle(4,PositionNbr)/FieldAtParticle(1,PositionNbr)

      ! New particles have not been reflected
      PartReflCount(PositionNbr) = 0
    END IF
  END DO

CASE DEFAULT
    CALL abort(__STAMP__,'Wrong particle velocity distribution!')

END SELECT

END SUBROUTINE SetParticleVelocity


END  MODULE MOD_Part_Pos_and_Velo
