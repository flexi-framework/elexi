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
!> Contains routines for statistical analysis of particle behavior
!===================================================================================================================================
MODULE MOD_Particle_Analyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitParticleAnalyze
  MODULE PROCEDURE InitParticleAnalyze
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

PUBLIC:: InitParticleAnalyze
PUBLIC:: FinalizeParticleAnalyze
PUBLIC:: DefineParametersParticleAnalyze
PUBLIC:: CalcKineticEnergy
PUBLIC:: CalcEkinPart
PUBLIC:: TrackingParticlePosition
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

CALL prms%CreateIntOption(      'Part-AnalyzeStep'   , 'Analyze is performed each Nth time step','1') 
CALL prms%CreateLogicalOption(  'CalcKineticEnergy'  , 'Calculate Kinetic Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPartBalance'    , 'Calculate the Particle Power Balance'//&
                                                       '- input and outflow energy of all particles','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcVelos'          , 'Calculate thermal and flow velocities.'//&
                                                       'if CalcVelos = T VelocityDirections = (/[int],[int],[int],[int]/)  '//&
                                                       'Switching dimensions for CalcVelos on (1) or off (0)\n'//&
                                                       '(/v_x,v_y,v_z,|v|/) ','.FALSE.')
CALL prms%CreateIntArrayOption( 'VelocityDirections' , 'x,y,z,abs -> 0/1 = T/F. (please note: CalcVelos)'&
                                                     ,'1 , 1 , 1 , 1')
CALL prms%CreateLogicalOption(  'Part-TrackPosition' , 'Track particle position','.FALSE.')
CALL prms%CreateLogicalOption(  'printDiff'          , '.FALSE.')
CALL prms%CreateRealOption(     'PrintDiffTime'      , '12')
CALL prms%CreateRealArrayOption('printDiffVec'       , '0. , 0. , 0. , 0. , 0. , 0.')
CALL prms%CreateLogicalOption(  'IsRestart'          , 'Flag, if the current calculation is a restart. ','.FALSE.')

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
  USE MOD_Particle_Vars           ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: dir, VeloDirs_hilf(4)
!===================================================================================================================================
  IF (ParticleAnalyzeInitIsDone) THEN
CALL abort(__STAMP__,&
'InitParticleAnalyse already called.',999,999.)
    RETURN
  END IF
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ANALYZE...'

  PartAnalyzeStep = GETINT('Part-AnalyzeStep','1')
IF (PartAnalyzeStep.EQ.0) PartAnalyzeStep = 123456789

  DoAnalyze = .FALSE.
  CalcEkin = GETLOGICAL('CalcKineticEnergy','.FALSE.')
  IF(CalcEkin) DoAnalyze = .TRUE.
  IF(nSpecies.GT.1) THEN
    nSpecAnalyze = nSpecies + 1
  ELSE
    nSpecAnalyze = 1
  END IF
  CalcPartBalance = GETLOGICAL('CalcPartBalance','.FALSE.')
  IF (CalcPartBalance) THEN
    DoAnalyze = .TRUE.
    SDEALLOCATE(nPartIn)
    SDEALLOCATE(nPartOut)
    SDEALLOCATE(PartEkinIn)
    SDEALLOCATE(PartEkinOut)
    ALLOCATE( nPartIn(nSpecies)     &
            , nPartOut(nSpecies)     &
            , PartEkinOut(nSpecies) &
            , PartEkinIn(nSpecies)  )
    nPartIn=0
    nPartOut=0
    PartEkinOut=0.
    PartEkinIn=0.

    SDEALLOCATE( nPartInTmp)
    SDEALLOCATE( PartEkinInTmp)
    ALLOCATE( nPartInTmp(nSpecies)     &
            , PartEkinInTmp(nSpecies)  )
    PartEkinInTmp=0.
    nPartInTmp=0

  END IF
TrackParticlePosition = GETLOGICAL('Part-TrackPosition','.FALSE.')
IF(TrackParticlePosition)THEN
  printDiff=GETLOGICAL('printDiff','.FALSE.')
  IF(printDiff)THEN
    printDiffTime=GETREAL('printDiffTime','12.')
    printDiffVec=GETREALARRAY('printDiffVec',6,'0.,0.,0.,0.,0.,0.')
  END IF
END IF

  CalcVelos = GETLOGICAL('CalcVelos','.FALSE')
  IF (CalcVelos) THEN
    DoAnalyze=.TRUE.
    VeloDirs_hilf = GetIntArray('VelocityDirections',4,'1,1,1,1') ! x,y,z,abs -> 0/1 = T/F
    VeloDirs(:) = .FALSE.
  IF(.NOT.CalcNumSpec)THEN
    SWRITE(UNIT_stdOut,'(A)') ' Velocity computation requires NumSpec and SimNumSpec. Setting CalcNumSpec=.TRUE.'
    CalcNumSpec = .TRUE.
  END IF
    DO dir = 1,4
      IF (VeloDirs_hilf(dir) .EQ. 1) THEN
        VeloDirs(dir) = .TRUE.
      END IF
    END DO
    IF ((.NOT. VeloDirs(1)) .AND. (.NOT. VeloDirs(2)) .AND. &
        (.NOT. VeloDirs(3)) .AND. (.NOT. VeloDirs(4))) THEN
      CALL abort(&
        __STAMP__&
        ,'No VelocityDirections set in CalcVelos!')
    END IF
  END IF

  IsRestart = GETLOGICAL('IsRestart','.FALSE.')

  ParticleAnalyzeInitIsDone=.TRUE.

  SWRITE(UNIT_stdOut,'(A)')' INIT PARTCILE ANALYZE DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleAnalyze


SUBROUTINE CalcKineticEnergy(Ekin) 
!===================================================================================================================================
! compute the kinetic energy of particles
! for velocity <1e3 non-relativistic formula is used, for larger velocities the relativistic kinetic energy is computed
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
      partV2 = PartState(i,4) * PartState(i,4) &
              + PartState(i,5) * PartState(i,5) &
              + PartState(i,6) * PartState(i,6)
          Ekin(nSpecAnalyze) = Ekin(nSpecAnalyze) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
ELSE ! nSpecAnalyze = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      partV2 = PartState(i,4) * PartState(i,4) &
             + PartState(i,5) * PartState(i,5) &
             + PartState(i,6) * PartState(i,6)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
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


!SUBROUTINE CalcNumPartsOfSpec(NumSpec,SimNumSpec) 
!!===================================================================================================================================
!! computes the number of simulated particles AND number of real particles within the domain
!! CAUTION: SimNumSpec equals NumSpec only for constant MPF, not vMPF
!! Last section of the routine contains the MPI-communication
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!USE MOD_Globals
!USE MOD_Particle_Vars,          ONLY: PDM,PartSpecies
!USE MOD_Particle_Analyze_Vars,  ONLY: nSpecAnalyze
!USE MOD_Particle_Vars,          ONLY: nSpecies
!#if USE_MPI
!USE MOD_Particle_MPI_Vars,      ONLY: PartMPI
!USE MOD_Particle_Analyze_Vars,  ONLY: CalcNumSpec
!#endif /*MPI*/
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES 
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)                   :: NumSpec(nSpecAnalyze)
!INTEGER(KIND=8),INTENT(OUT)        :: SimNumSpec(nSpecAnalyze)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                            :: iPart
!#if USE_MPI
!REAL                               :: RD(nSpecAnalyze)
!INTEGER(KIND=8)                    :: ID(nSpecAnalyze)
!#endif /*MPI*/
!!===================================================================================================================================
!
!  SimNumSpec = 0
!  DO iPart=1,PDM%ParticleVecLength
!    IF (PDM%ParticleInside(iPart)) THEN
!      SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1                      ! NumSpec =  particle number
!    END IF
!  END DO
!  NumSpec = REAL(SimNumSpec)
!  IF(nSpecAnalyze.GT.1)THEN
!    SimNumSpec(nSpecAnalyze) = SUM(SimNumSpec(1:nSpecies))
!    NumSpec(nSpecAnalyze) = SUM(NumSpec(1:nSpecies))
!  END IF
!!END IF
!
!#if USE_MPI
!IF (PartMPI%MPIRoot) THEN
!  CALL MPI_REDUCE(MPI_IN_PLACE,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
!  IF(CalcNumSpec) & 
!  CALL MPI_REDUCE(MPI_IN_PLACE,SimNumSpec ,nSpecAnalyze,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
!ELSE
!  CALL MPI_REDUCE(NumSpec     ,RD         ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
!  IF(CalcNumSpec) & 
!  CALL MPI_REDUCE(SimNumSpec  ,ID         ,nSpecAnalyze,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
!END IF
!#endif /*MPI*/
!
!END SUBROUTINE CalcNumPartsOfSpec

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

  partV2 = PartState(iPart,4) * PartState(iPart,4) &
         + PartState(iPart,5) * PartState(iPart,5) &
         + PartState(iPart,6) * PartState(iPart,6)
         
  Ekin= 0.5*Species(PartSpecies(iPart))%MassIC*partV2* Species(PartSpecies(iPart))%MacroParticleFactor
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
USE MOD_Particle_Globals,       ONLY: GETFREEUNIT
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
CHARACTER(LEN=60) :: TrackingFilename
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
    CLOSE(iunit)
#if USE_MPI
  END IF
#endif
ELSE
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
        WRITE(iunit,104,ADVANCE='NO') PartState(i,iPartState)
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
