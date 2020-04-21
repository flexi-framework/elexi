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
!> Contains routines for statistical analysis of particle behavior
!===================================================================================================================================
MODULE MOD_Particle_Analyze
! IMPLICIT VARIABLE HANDLING
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

INTERFACE CalcEkinPart
  MODULE PROCEDURE CalcEkinPart
END INTERFACE

INTERFACE TrackingParticlePosition
  MODULE PROCEDURE TrackingParticlePosition
END INTERFACE

PUBLIC :: InitParticleAnalyze
PUBLIC :: ParticleAnalyze
PUBLIC :: ParticleInformation
PUBLIC :: FinalizeParticleAnalyze
PUBLIC :: DefineParametersParticleAnalyze
PUBLIC :: CalcKineticEnergy
PUBLIC :: CalcEkinPart
PUBLIC :: TrackingParticlePosition
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

CALL prms%CreateIntOption(      'Part-AnalyzeStep'      , 'Analyze is performed each Nth time step','1')
CALL prms%CreateLogicalOption(  'CalcKineticEnergy'     , 'Calculate Kinetic Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPartBalance'       , 'Calculate the Particle Power Balance'//&
                                                          '- input and outflow energy of all particles','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackPosition'    , 'Track particle position','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackConvergence' , 'Track particle convergence (i.e. final position)','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackDispersion'  , 'Track particle convergence radius (i.e. absolute path)','.FALSE.')
CALL prms%CreateLogicalOption(  'printDiff'             , '.FALSE.')
CALL prms%CreateRealOption(     'PrintDiffTime'         , '12')
CALL prms%CreateRealArrayOption('printDiffVec'          , '0. , 0. , 0. , 0. , 0. , 0.')

END SUBROUTINE DefineParametersParticleAnalyze


SUBROUTINE InitParticleAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars
USE MOD_ReadInTools             ,ONLY: GETLOGICAL, GETINT, GETSTR, GETINTARRAY, GETREALARRAY, GETREAL
USE MOD_Particle_Vars           ,ONLY: nSpecies,PDM
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

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ANALYZE...'

DoAnalyze = .FALSE.
  CalcEkin = GETLOGICAL('CalcKineticEnergy','.FALSE.')
  IF(CalcEkin) DoAnalyze = .TRUE.

IF(nSpecies.GT.1) THEN
  nSpecAnalyze = nSpecies + 1
ELSE
  nSpecAnalyze = 1
END IF

! Calculate number and kinetic energy of particles entering / leaving the domain
CalcPartBalance = GETLOGICAL('CalcPartBalance','.FALSE.')
IF (CalcPartBalance) THEN
  DoAnalyze = .TRUE.
  SDEALLOCATE(nPartIn)
  SDEALLOCATE(nPartOut)
  SDEALLOCATE(PartEkinIn)
  SDEALLOCATE(PartEkinOut)
  ALLOCATE( nPartIn(nSpecies)        &
          , nPartOut(nSpecies)       &
          , PartEkinOut(nSpecies)    &
          , PartEkinIn(nSpecies))
  nPartIn=0
  nPartOut=0
  PartEkinOut=0.
  PartEkinIn=0.

  SDEALLOCATE( nPartInTmp)
  SDEALLOCATE( PartEkinInTmp)
  ALLOCATE( nPartInTmp(nSpecies)     &
          , PartEkinInTmp(nSpecies))
  PartEkinInTmp=0.
  nPartInTmp=0
END IF

! Track and write position of each particle. Primarily intended for debug purposes
doParticlePositionTrack    = GETLOGICAL('Part-TrackPosition',   '.FALSE.')
doParticleConvergenceTrack = GETLOGICAL('Part-TrackConvergence','.FALSE.')

IF(doParticlePositionTrack) THEN
  printDiff=GETLOGICAL('printDiff','.FALSE.')
  IF(printDiff) THEN
    printDiffTime = GETREAL(     'printDiffTime','12.')
    printDiffVec  = GETREALARRAY('printDiffVec',6,'0.,0.,0.,0.,0.,0.')
  END IF
END IF

doParticleDispersionTrack = GETLOGICAL('Part-TrackDispersion','.FALSE.')
IF(doParticleDispersionTrack) THEN
  ALLOCATE(PartPath(1:3,PDM%maxParticleNumber))
  PartPath = 0.
END IF

ParticleAnalyzeInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTCILE ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleAnalyze


SUBROUTINE ParticleAnalyze(iter)
!==================================================================================================================================
!> Controls particle analysis routines and is called at analyze time levels
!> - calls erosion impcat output
!> - informs about lost particles
!> - writes load balancing statistics
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars    ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Analyze ,ONLY: CalcSurfaceValues
USE MOD_Particle_Tracking_Vars    ,ONLY: CountNbOfLostParts
USE MOD_Particle_Output           ,ONLY: WriteInfoStdOut
#if USE_LOADBALANCE
USE MOD_LoadDistribution          ,ONLY: WriteElemTimeStatistics
#endif /*LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)      :: iter                   !< current iteration
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Calculate particle surface erosion data
IF (WriteMacroSurfaceValues) THEN
  CALL CalcSurfaceValues
END IF

! Write information to console output
IF (CountNbOfLostParts) THEN
  CALL WriteInfoStdOut()
END IF

#if USE_LOADBALANCE
! Create .csv file for performance analysis and load balance: write header line
CALL WriteElemTimeStatistics(WriteHeader=.TRUE.,iter=iter)
#endif /*LOADBALANCE*/

END SUBROUTINE ParticleAnalyze


SUBROUTINE ParticleInformation()
!==================================================================================================================================
!> Displays the actual status of the particle phase in the simulation and counts the number of impacts
!==================================================================================================================================
USE MOD_Globals
USE MOD_PreProc
USE MOD_Erosionpoints_Vars        ,ONLY: EP_Buffersize
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
DO iPart=1,PDM%ParticleVecLength
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
SWRITE(UNIT_StdOut,'(A14,I16)')' # Impacts  : ', EP_Buffersize

END SUBROUTINE ParticleInformation



SUBROUTINE CalcKineticEnergy(Ekin)
!===================================================================================================================================
! compute the kinetic energy of particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,          ONLY : PartState, PartSpecies, Species, PDM
USE MOD_Particle_Analyze_Vars,  ONLY : nSpecAnalyze
#if USE_MPI
USE MOD_Particle_MPI_Vars,      ONLY : PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Ekin(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
REAL(KIND=8)              :: partV2
#if USE_MPI
REAL                      :: RD(nSpecAnalyze)
#endif /*MPI*/
!===================================================================================================================================

Ekin = 0.!d0
IF (nSpecAnalyze.GT.1) THEN
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      partV2 = PartState(4,i) * PartState(4,i) &
             + PartState(5,i) * PartState(5,i) &
             + PartState(6,i) * PartState(6,i)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + 0.5 *  Species(PartSpecies(i))%MassIC * partV2
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength

ELSE ! nSpecAnalyze = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      partV2 = PartState(4,i) * PartState(4,i) &
             + PartState(5,i) * PartState(5,i) &
             + PartState(6,i) * PartState(6,i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2
    END IF ! particle inside
  END DO ! particleveclength
END IF

#if USE_MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Ekin,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(Ekin  ,RD        ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcKineticEnergy


Function CalcEkinPart(iPart)
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
REAL                               :: partV2, Ekin !,gamma1
!===================================================================================================================================

  partV2 = PartState(4,iPart) * PartState(4,iPart) &
         + PartState(5,iPart) * PartState(5,iPart) &
         + PartState(6,iPart) * PartState(6,iPart)

  Ekin= 0.5*Species(PartSpecies(iPart))%MassIC*partV2
  CalcEkinPart=Ekin

END FUNCTION CalcEkinPart


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
USE MOD_Particle_Analyze_Vars,  ONLY:printDiff,printDiffVec,printDiffTime
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

IF (printDiff) THEN
  diffPos=0.
  diffVelo=0.
  IF (time.GE.printDiffTime) THEN
    printDiff=.FALSE.
    DO iPartState=1,3
      diffPos=diffPos+(printDiffVec(iPartState)-PartState(1,iPartState))**2
      diffVelo=diffVelo+(printDiffVec(iPartState+3)-PartState(1,iPartState+3))**2
    END DO
    WRITE(*,'(A,e24.14,x,e24.14)') 'L2-norm from printDiffVec: ',SQRT(diffPos),SQRT(diffVelo)
  END IF
END IF
104    FORMAT (e25.14)

END SUBROUTINE TrackingParticlePosition

SUBROUTINE FinalizeParticleAnalyze()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars,ONLY:ParticleAnalyzeInitIsDone
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ParticleAnalyzeInitIsDone = .FALSE.
END SUBROUTINE FinalizeParticleAnalyze


END MODULE MOD_Particle_Analyze
