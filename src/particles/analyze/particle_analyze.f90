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
!> Contains routines for statistical analysis of particle behavior
!===================================================================================================================================
MODULE MOD_Particle_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleAnalyze
  MODULE PROCEDURE InitParticleAnalyze
END INTERFACE

INTERFACE ParticleAnalyze
  MODULE PROCEDURE ParticleAnalyze
END INTERFACE

INTERFACE ParticleInformation
  MODULE PROCEDURE ParticleInformation
END INTERFACE

INTERFACE FinalizeParticleAnalyze
  MODULE PROCEDURE FinalizeParticleAnalyze
END INTERFACE

INTERFACE CalcKineticEnergy
  MODULE PROCEDURE CalcKineticEnergy
END INTERFACE

INTERFACE TrackingParticlePosition
  MODULE PROCEDURE TrackingParticlePosition
END INTERFACE

INTERFACE TrackingParticlePath
  MODULE PROCEDURE TrackingParticlePath
END INTERFACE

PUBLIC :: InitParticleAnalyze
PUBLIC :: ParticleAnalyze
PUBLIC :: ParticleInformation
PUBLIC :: FinalizeParticleAnalyze
PUBLIC :: DefineParametersParticleAnalyze
PUBLIC :: CalcKineticEnergy
PUBLIC :: TrackingParticlePosition
PUBLIC :: TrackingParticlePath
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for analyze if particles (.csv output)
!==================================================================================================================================
SUBROUTINE DefineParametersParticleAnalyze()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
!USE MOD_AnalyzeEquation ,ONLY: DefineParametersAnalyzeEquation
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Analyze")

CALL prms%CreateIntOption(      'Part-AnalyzeStep'        , 'Analyze is performed each Nth time step','1')
CALL prms%CreateLogicalOption(  'CalcKineticEnergy'       , 'Calculate Kinetic Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPartBalance'         , 'Calculate the Particle Kinetic Energy Balance'//&
                                                            '- input and outflow kinetic energy of all particles','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackPosition'      , 'Track particle position','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackConvergence'   , 'Track particle convergence (i.e. final position)','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackDispersion'    , 'Track particle convergence radius (i.e. absolute path)','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackPath'          , 'Track particle path (i.e. relative path)','.FALSE.')
#if CODE_ANALYZE
CALL prms%CreateIntOption(      'PartOut'                 , 'If compiled with CODE_ANALYZE flag: For This particle number'   //&
                                                            ' every tracking information is written as STDOUT.'                &
                                                          , '0')
CALL prms%CreateIntOption(      'MPIRankOut'              , 'If compiled with CODE_ANALYZE flag: This MPI-Proc writes the'   //&
                                                            ' tracking information for the defined PartOut.'
                                                          , '0')
#endif /*CODE_ANALYZE*/
CALL prms%CreateIntOption(      'Part-nRPs'               , 'Number of record planes'                                          &
                                                          , '0')
CALL prms%CreateIntOption(      'Part-RPMemory'           , 'Record particles memory'                                          &
                                                          , '100')
!CALL prms%CreateStringOption(  'Part-RecordType[$]'        , 'Type of record plane.\n'                                        //&
!                                                             ' - plane\n'                                                       &
!                                                           , 'none', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-RPDirection[$]'     , 'Direction of the normal vector of the record plane'               &
                                                          , '1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-RPPosition[$]'      , 'Position of the record plane in RPDirection.'                     &
                                                          , '0.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-FilenameRecordPart' , 'Specifying filename for load_from_file init.'                     &
                                                          , 'data/recordplane_')
CALL prms%CreateLogicalOption(  'doPartIndex'             , 'Flag to write out unique part index'                              &
                                                          , '.FALSE.')

END SUBROUTINE DefineParametersParticleAnalyze


SUBROUTINE InitParticleAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars
USE MOD_Particle_Boundary_Vars  ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Vars  ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Tracking_Vars  ,ONLY: CountNbOfLostParts
USE MOD_Particle_Vars           ,ONLY: nSpecies,PDM
USE MOD_ReadInTools             ,ONLY: GETLOGICAL,GETINT,GETSTR,GETINTARRAY,GETREALARRAY,GETREAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF (ParticleAnalyzeInitIsDone) THEN
  CALL abort(__STAMP__,'InitParticleAnalyse already called.',999,999.)
  RETURN
END IF

!SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' INIT PARTICLE ANALYZE...'

IF (WriteMacroSurfaceValues .OR. doParticleImpactTrack .OR. RecordPart.GT.0 .OR. CountNbOfLostParts .OR. CalcEkin .OR. DoAnalyze) THEN
  doParticleAnalyze = .TRUE.
ELSE
  doParticleAnalyze = .FALSE.
END IF

DoAnalyze    = .FALSE.
nSpecAnalyze = MERGE(nSpecies + 1,1,nSpecies.GT.1)

CalcEkin = GETLOGICAL('CalcKineticEnergy','.FALSE.')
IF (CalcEkin) THEN
  DoAnalyze = .TRUE.
  SDEALLOCATE(PartEkin)
  ALLOCATE( PartEkin (nSpecAnalyze))
  PartEkin = 0.
END IF

! Calculate number and kinetic energy of particles entering / leaving the domain
CalcPartBalance = GETLOGICAL('CalcPartBalance','.FALSE.')
IF (CalcPartBalance) THEN
  DoAnalyze = .TRUE.
  SDEALLOCATE(nPartIn)
  SDEALLOCATE(nPartOut)
  SDEALLOCATE(PartEkinIn)
  SDEALLOCATE(PartEkinOut)
  ALLOCATE( nPartIn    (nSpecAnalyze)    &
          , nPartOut   (nSpecAnalyze)    &
          , PartEkinOut(nSpecAnalyze)    &
          , PartEkinIn (nSpecAnalyze))
  nPartIn     = 0
  nPartOut    = 0
  PartEkinOut = 0.
  PartEkinIn  = 0.

  SDEALLOCATE( nPartInTmp)
  SDEALLOCATE( PartEkinInTmp)
  ALLOCATE( nPartInTmp(nSpecies)     &
          , PartEkinInTmp(nSpecies))
  PartEkinInTmp = 0.
  nPartInTmp    = 0
END IF

! Track and write position of each particle. Primarily intended for debug purposes
doParticlePositionTrack    = GETLOGICAL('Part-TrackPosition',   '.FALSE.')
doParticleConvergenceTrack = GETLOGICAL('Part-TrackConvergence','.FALSE.')

doParticleDispersionTrack    = GETLOGICAL('Part-TrackDispersion','.FALSE.')
doParticlePathTrack          = GETLOGICAL('Part-TrackPath'      ,'.FALSE.')
IF (doParticleDispersionTrack .AND. doParticlePathTrack) &
  CALL CollectiveStop(__STAMP__,'doParticleDispersionTrack and doParticlePathTrack cannot be enabled at the same time!')

IF(doParticleDispersionTrack .OR. doParticlePathTrack) THEN
  ALLOCATE(PartPath(1:3,PDM%maxParticleNumber))
  PartPath = 0.
END IF

ParticleAnalyzeInitIsDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTCILE ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleAnalyze


SUBROUTINE ParticleAnalyze(t    &
#if USE_LOADBALANCE
                          ,iter &
#endif /* USE_LOADBALANCE */
  )
!==================================================================================================================================
!> Controls particle analysis routines and is called at analyze time levels
!> - writes load balancing statistics
!> - calls impact sampling output
!> - calls impact tracking output
!> - informs about lost particles
!> - calls kinetic energy balance
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Analyze_Tools    ,ONLY: ParticleRecord
USE MOD_Particle_Analyze_Vars     ,ONLY: DoAnalyze,DoParticleAnalyze,CalcEkin
USE MOD_Particle_Analyze_Vars     ,ONLY: RecordPart
USE MOD_Particle_Boundary_Vars    ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars    ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Analyze ,ONLY: CalcSurfaceValues,WriteBoundaryParticleToHDF5
USE MOD_Particle_Output           ,ONLY: WriteParticleAnalyze,WriteInfoStdOut
USE MOD_Particle_Tracking_Vars    ,ONLY: CountNbOfLostParts
#if USE_LOADBALANCE
USE MOD_LoadDistribution          ,ONLY: WriteElemTimeStatistics
#endif /*LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< current simulation time
#if USE_LOADBALANCE
INTEGER(KIND=8),INTENT(IN)      :: iter                   !< current iteration
#endif /* USE_LOADBALANCE */
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

#if USE_LOADBALANCE
! Create .csv file for performance analysis and load balance: write header line
CALL WriteElemTimeStatistics(WriteHeader=.TRUE.,iter=iter)
#endif /*LOADBALANCE*/

! Prettify output
IF (doParticleAnalyze) THEN
  SWRITE(UNIT_StdOut,'(132("-"))')
END IF

! Calculate cumulative particle surface impact sampling data
IF (WriteMacroSurfaceValues) THEN
  CALL CalcSurfaceValues
END IF

! Write individual particle surface impact tracking data
IF (doParticleImpactTrack) THEN
  CALL WriteBoundaryParticleToHDF5(OutputTime=t)
END IF

! Write individual particle record plane data
IF (RecordPart.GT.0) THEN
  CALL ParticleRecord(t,writeToBinary=.TRUE.)
END IF

! Write information to console output
IF (CountNbOfLostParts) THEN
  CALL WriteInfoStdOut()
END IF

IF (CalcEkin) THEN
  CALL CalcKineticEnergy()
END IF

IF (DoAnalyze) THEN
  CALL WriteParticleAnalyze()
END IF

END SUBROUTINE ParticleAnalyze


SUBROUTINE ParticleInformation()
!==================================================================================================================================
!> Displays the actual status of the particle phase in the simulation and counts the number of impacts
!==================================================================================================================================
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Boundary_Vars    ,ONLY: doParticleImpactTrack,ImpactnGlob
USE MOD_Particle_Vars             ,ONLY: PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: nParticleOnProc,iPart
INTEGER :: nParticleInDomain
!==================================================================================================================================

! Count number of particles on local proc
nParticleOnProc = 0
DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) nParticleOnProc = nParticleOnProc + 1
END DO
#if USE_MPI
! Gather number of particles on all procs
CALL MPI_REDUCE(nParticleOnProc,nParticleInDomain,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
#else
nParticleInDomain = nParticleOnProc
#endif

! Output particle and impact information
!SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(UNIT_StdOut,'(A14,I16)')' # Particle : ', nParticleInDomain
IF (doParticleImpactTrack) THEN
SWRITE(UNIT_StdOut,'(A14,I16)')' # Impacts  : ', ImpactnGlob
END IF

END SUBROUTINE ParticleInformation


SUBROUTINE CalcKineticEnergy()
!===================================================================================================================================
! compute the kinetic energy of particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Globals      ,ONLY: DOTPRODUCT
USE MOD_Particle_Analyze_Vars ,ONLY: PartEkin,nSpecAnalyze
USE MOD_Particle_Vars         ,ONLY: PartState,PartSpecies,Species,PDM
#if USE_MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
REAL(KIND=8)      :: partV2
!===================================================================================================================================

PartEkin = 0.

! nSpecAnalyze > 1, calculate species and sum
IF (nSpecAnalyze.GT.1) THEN
  DO i = 1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      partV2   = DOTPRODUCT(PartState(4:6,i))
      PartEkin(nSpecAnalyze)   = PartEkin(nSpecAnalyze)   + 0.5 *  Species(PartSpecies(i))%MassIC * partV2
      PartEkin(PartSpecies(i)) = PartEkin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
! nSpecAnalyze = 1 : only 1 species
ELSE
  DO i = 1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      partV2   = DOTPRODUCT(PartState(4:6,i))
      PartEkin(PartSpecies(i)) = PartEkin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2
    END IF ! particle inside
  END DO ! 1,PDM%ParticleVecLength
END IF

#if USE_MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PartEkin,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(PartEkin    ,0       ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcKineticEnergy


SUBROUTINE TrackingParticlePosition(time)
!===================================================================================================================================
! Outputs the particle position and velocity at every time step. Use only for debugging purposes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,          ONLY:PartState, PDM, PEM
#if USE_MPI
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
#endif
USE MOD_Particle_Globals,       ONLY:GETFREEUNIT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iunit,iPartState
CHARACTER(LEN=60)  :: TrackingFilename
LOGICAL            :: fexist
REAL               :: diffPos,diffVelo
!===================================================================================================================================

iunit=GETFREEUNIT()
TrackingFilename = ('ParticlePosition.csv')

INQUIRE(FILE = TrackingFilename, EXIST=fexist)

IF(.NOT.fexist) THEN
#if USE_MPI
 IF(PartMPI%MPIRoot)THEN
#endif
   iunit=GETFREEUNIT()
   OPEN(UNIT=iunit,FILE=TrackingFilename,FORM='FORMATTED',STATUS='UNKNOWN')
   !CALL FLUSH (iunit)
    ! writing header
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'TIME', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartNum', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosX', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosY', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosZ', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelX', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelY', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelZ', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartElem', ' '
    CLOSE(iunit)
#if USE_MPI
  END IF
  ! Wait with all procs until the file is available.
  ! WARNING: Global sync point, but this routine is only supposed to work for debug anyways
  CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif
END IF

! Check again, we might be on the other proc and only see the file now
INQUIRE(FILE = TrackingFilename, EXIST=fexist)
IF(fexist) THEN
  iunit=GETFREEUNIT()
  OPEN(unit=iunit,FILE=TrackingFileName,FORM='Formatted',POSITION='APPEND',STATUS='old')
  !CALL FLUSH (iunit)
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      WRITE(iunit,104,ADVANCE='NO') TIME
      WRITE(iunit,'(A1)',ADVANCE='NO') ','
      WRITE(iunit,'(I12)',ADVANCE='NO') i
      DO iPartState=1,6
        WRITE(iunit,'(A1)',ADVANCE='NO') ','
        WRITE(iunit,104,ADVANCE='NO') PartState(iPartState,i)
      END DO
      WRITE(iunit,'(A1)',ADVANCE='NO') ','
      WRITE(iunit,'(I12)',ADVANCE='NO') PEM%Element(i)
      WRITE(iunit,'(A)') ' '
     END IF
  END DO
  CLOSE(iunit)
END IF

104    FORMAT (e25.14)

END SUBROUTINE TrackingParticlePosition


SUBROUTINE TrackingParticlePath()
!===================================================================================================================================
! Outputs the particle position and velocity at every time step. Use only for debugging purposes
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



SUBROUTINE FinalizeParticleAnalyze()
!===================================================================================================================================
! Finalizes variables necessary for analyze subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(PartEkin)

SDEALLOCATE(nPartIn)
SDEALLOCATE(nPartOut)
SDEALLOCATE(PartEkinIn)
SDEALLOCATE(PartEkinOut)

SDEALLOCATE(nPartInTmp)
SDEALLOCATE(PartEkinInTmp)

SDEALLOCATE(RPP_Plane)

ParticleAnalyzeInitIsDone = .FALSE.

END SUBROUTINE FinalizeParticleAnalyze


END MODULE MOD_Particle_Analyze
