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
!> Write particle data to Tecplot ASCII file
!===================================================================================================================================
MODULE MOD_Particle_Output
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

!INTERFACE InitParticleOutput
!  MODULE PROCEDURE InitParticleOutput
!END INTERFACE

INTERFACE WriteInfoStdOut
  MODULE PROCEDURE WriteInfoStdOut
END INTERFACE

INTERFACE WriteParticleAnalyze
  MODULE PROCEDURE WriteParticleAnalyze
END INTERFACE

INTERFACE Visualize_Particles
  MODULE PROCEDURE Visualize_Particles
END INTERFACE

!PUBLIC :: InitParticleOutput
PUBLIC :: WriteInfoStdOut
PUBLIC :: WriteParticleAnalyze
PUBLIC :: Visualize_Particles
!===================================================================================================================================

CONTAINS

!SUBROUTINE InitParticleOutput()
!!===================================================================================================================================
!! Initialize all output (and analyze) variables.
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Particle_Output_Vars,ONLY:ParticleOutputInitIsDone
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!===================================================================================================================================
!IF (ParticleOutputInitIsDone) &
!  CALL ABORT(__STAMP__,'InitOutput not ready to be called or already called.',999,999.)
!
!SWRITE(UNIT_StdOut,'(132("-"))')
!SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE OUTPUT...'
!
!ParticleOutputInitIsDone = .TRUE.
!SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE OUTPUT DONE!'
!SWRITE(UNIT_StdOut,'(132("-"))')
!END SUBROUTINE InitParticleOutput


SUBROUTINE WriteInfoStdOut()
!===================================================================================================================================
! Writes particle information to standard output
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Tracking_Vars,ONLY: NbrOfLostParticles,countNbOfLostParts
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: nLostPartsTot
!===================================================================================================================================

IF (CountNbOfLostParts) THEN
#if USE_MPI
  CALL MPI_REDUCE(NbrOfLostParticles,nLostPartsTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
#else
  nLostPartsTot = NbrOfLostParticles
#endif /*MPI*/
  WRITE(UNIT_stdOut,'(A,I12)')' NbOfLostParticle : ',nLostPartsTot
END IF

END SUBROUTINE WriteInfoStdOut


SUBROUTINE WriteParticleAnalyze()
!===================================================================================================================================
! Writes particle information to standard PartAnalyze file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Analyze_Vars
USE MOD_Particle_Vars             ,ONLY: nSpecies
USE MOD_Restart_Vars              ,ONLY: DoRestart
USE MOD_TimeDisc_Vars             ,ONLY: t
#if USE_MPI
USE MOD_Particle_MPI_Vars         ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: isOpen
CHARACTER(LEN=350)  :: outfile
INTEGER             :: unit_index,iSpec,OutputCounter
!===================================================================================================================================

OutputCounter = 2
unit_index = 535

IF (MPIRoot) THEN
  INQUIRE(UNIT=unit_index, OPENED=isOpen)
  ! Continue analyze file if opened, otherwise append or create
  IF (.NOT.isOpen) THEN
    outfile = 'PartAnalyze.csv'
    IF (DoRestart .and. FILEEXISTS(outfile)) THEN
      OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
    ELSE
      OPEN(unit_index,file=TRIM(outfile))
      !--- insert header
      WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
      IF (CalcPartBalance) THEN
        DO iSpec = 1,nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A14,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartIn-Spec-',iSpec,' '
          OutputCounter = OutputCounter + 1
        END DO
        DO iSpec=1, nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A15,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartOut-Spec-',iSpec,' '
          OutputCounter = OutputCounter + 1
        END DO
!     END IF
!     IF (CalcPartBalance) THEN
        DO iSpec = 1,nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A8,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinIn-',iSpec,' '
          OutputCounter = OutputCounter + 1
        END DO
        DO iSpec=1, nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A9,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinOut-',iSpec,' '
          OutputCounter = OutputCounter + 1
        END DO
      END IF

      WRITE(unit_index,'(A1)') ' '
    END IF
  END IF
END IF

! Communicate analysis variables
#if USE_MPI
IF (MPIRoot) THEN
  IF (CalcPartBalance) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,nPartIn    (1:nSpecAnalyze),nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,nPartOut   (1:nSpecAnalyze),nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinIn (1:nSpecAnalyze),nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinOut(1:nSpecAnalyze),nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
ELSE
  IF (CalcPartBalance) THEN
    CALL MPI_REDUCE(nPartIn    (1:nSpecAnalyze),0           ,nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(nPartOut   (1:nSpecAnalyze),0           ,nSpecAnalyze,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(PartEkinIn (1:nSpecAnalyze),0           ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(PartEkinOut(1:nSpecAnalyze),0           ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
END IF
#endif /*USE_MPI*/

! Perform averaging/summation of the MPI communicated variables on the root only
IF (MPIRoot) THEN
  IF (CalcPartBalance) THEN
    IF(nSpecies.GT.1) THEN
      nPartIn    (nSpecAnalyze) = SUM(nPartIn    (1:nSpecies))
      nPartOut   (nSpecAnalyze) = SUM(nPartOut   (1:nSpecies))
      PartEkinIn (nSpecAnalyze) = SUM(PartEkinIn (1:nSpecies))
      PartEkinOut(nSpecAnalyze) = SUM(PartEkinOut(1:nSpecies))
    END IF
  END IF
END IF

! Output
IF (MPIRoot) THEN
  ! Time
  WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') t
  IF (CalcPartBalance) THEN
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartIn(iSpec))
    END DO
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartOut(iSpec))
    END DO
!  END IF
!  IF (CalcPartBalance) THEN
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkinIn(iSpec)
    END DO
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkinOut(iSpec)
    END DO
  END IF

  WRITE(unit_index,'(A1)') ' '
END IF

! Reset the particle counter
IF (CalcPartBalance) THEN
  nPartIn     = 0
  nPartOut    = 0
  PartEkinIn  = 0.
  PartEkinOut = 0.
END IF


END SUBROUTINE WriteParticleAnalyze


SUBROUTINE Visualize_Particles(OutputTime)
!===================================================================================================================================
! Simple visualization of conservative variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars     ,ONLY: ProjectName
USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nParts, i, index_unit
CHARACTER(LEN=255)            :: FileString
!===================================================================================================================================

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Particles',OutputTime))//'.dat'
! Visualize data
nParts = 0
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    nParts = nParts + 1
  END IF
END DO
!CALL WriteDataToTecplotBinary(NVisu,PP_nElems,PP_nVar,0,VarNames,Coords_NVisu(1:3,:,:,:,:),U_NVisu,TRIM(FileString))
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE PARTICLE DATA TO TECPLOT ASCII FILE..."

index_unit = 45
OPEN(index_unit,FILE=TRIM(FileString),Status="REPLACE")
WRITE(index_unit,*)'TITLE = "',TRIM(ProjectName),'"'
WRITE(index_unit,'(A)')'VARIABLES = "x[m]" ,"y[m]" ,"z[m]" ,"v_x[m/s]" ,"v_y[m/s]" ,"v_z[m/s]" ,"m[kg]" ,"Particle_Number"'
WRITE(index_unit,*)'ZONE T= "',TRIM(TIMESTAMP('Particles',OutputTime)),'"'
WRITE(index_unit,*)'I=',nParts,' J=1, K=1, F=POINT'
DO i=1,nParts
  WRITE(index_unit,'(8(1X,e19.12),1X,i0)')PartState(1,i),&
                                          PartState(2,i),&
                                          PartState(3,i),&
                                          PartState(4,i),&
                                          PartState(5,i),&
                                          PartState(6,i),&
                                          Species(PartSpecies(i))%MassIC,&
                                          i
END DO
CLOSE(index_unit)
SWRITE(UNIT_stdOut,'(A)')"DONE!"

END SUBROUTINE Visualize_Particles


END MODULE MOD_Particle_Output
