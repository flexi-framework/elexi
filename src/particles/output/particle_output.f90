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

INTERFACE FillParticleData
  MODULE PROCEDURE FillParticleData
END INTERFACE

INTERFACE WriteInfoStdOut
  MODULE PROCEDURE WriteInfoStdOut
END INTERFACE

INTERFACE WriteParticleAnalyze
  MODULE PROCEDURE WriteParticleAnalyze
END INTERFACE

INTERFACE Visualize_Particles
  MODULE PROCEDURE Visualize_Particles
END INTERFACE

INTERFACE GetOffsetAndGlobalNumberOfParts
  MODULE PROCEDURE GetOffsetAndGlobalNumberOfParts
END INTERFACE

!PUBLIC :: InitParticleOutput
PUBLIC :: FillParticleData
PUBLIC :: WriteInfoStdOut
PUBLIC :: WriteParticleAnalyze
PUBLIC :: Visualize_Particles
PUBLIC :: GetOffsetAndGlobalNumberOfParts
!===================================================================================================================================

CONTAINS

SUBROUTINE FillParticleData()
!===================================================================================================================================
! Fills the particle data arrays required for loadbalance/output
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Mesh_Vars               ,ONLY: offsetElem
USE MOD_Part_Tools              ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze_Vars   ,ONLY: PartPath,doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars  ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Output_Vars
USE MOD_Particle_Vars           ,ONLY: PartDataSize,TurbPartDataSize
USE MOD_Particle_Vars           ,ONLY: PDM,PEM,PartState,PartSpecies,PartReflCount,PartIndex,nSpecies
USE MOD_Particle_Vars           ,ONLY: PartInt,PartData,TurbPartData
USE MOD_Particle_Vars           ,ONLY: useLinkedList,doPartIndex,doWritePartDiam
! Particle turbulence models
USE MOD_Particle_Vars           ,ONLY: TurbPartState
USE MOD_Particle_SGS_Vars       ,ONLY: nSGSVars
#if USE_RW
USE MOD_Particle_RandomWalk_Vars,ONLY: nRWVars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: PartIntSize=2
INTEGER                        :: pcount,nInvalidPart
INTEGER                        :: iPart,iElem
INTEGER                        :: PP_nVarPart_loc
!===================================================================================================================================

! Return if running particle code without any species
IF (nSpecies.LE.0) RETURN

! Size and location of particle data
PP_nVarPart_loc = PP_nVarPart-1
PartDataSize    = PP_nVarPart_loc + 1
PartDataVarSpecies       = PartDataSize
IF (doWritePartDiam) THEN
  PP_nVarPart_loc = PP_nVarPart
  PartDataSize    = PartDataSize + 1
  PartDataVarSpecies       = PartDataSize
END IF
PartDataVarStart     = PartDataSize + 1
! Increase size if index is tracked
IF (doPartIndex) THEN
  PartDataSize = PartDataSize + 1
  PartDataVarStart     = PartDataVarStart     + 1
END IF
PartDataVarShift     = 0
! Increase size if reflections are tracked
IF (doParticleReflectionTrack) THEN
  PartDataSize = PartDataSize + 1
  PartDataVarShift     = 1
END IF
! Increase size if the absolute particle path is tracked
IF (doParticleDispersionTrack.OR.doParticlePathTrack) &
  PartDataSize = PartDataSize + 3

! Add turbulent dispersion data to output
IF (ALLOCATED(TurbPartState)) THEN
  TurbPartDataSize = nSGSVars
#if USE_RW
  TurbPartDataSize = TurbPartDataSize + nRWVars
#endif
END IF

nInvalidPart = 0
! Make sure to eliminate invalid particles as we cannot restart from NaN
DO pcount = 1,PDM%ParticleVecLength
  IF (.NOT. PDM%ParticleInside(pcount)) CYCLE
  IF (ANY(IEEE_IS_NAN(PartState(:,pcount)))) THEN
    nInvalidPart = nInvalidPart + 1
    PDM%ParticleInside(pcount) = .FALSE.
  END IF
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nInvalidPart,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(nInvalidPart,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif /*USE_MPI*/
IF (nInvalidPart.GT.0 .AND. MPIRoot) THEN
  WRITE(UNIT_stdOut,'(A,I0,A)') ' Detected ',nInvalidPart,' invalid particles during output. Removing ...'
END IF

! Determine number of particles in the complete domain
locnPart =   0
!>> Count number of particle on local proc
DO pcount = 1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(pcount)) THEN
    locnPart = locnPart + 1
  END IF
END DO

! Communicate the total number and offset
CALL GetOffsetAndGlobalNumberOfParts('WriteParticleToHDF5',offsetnPart,nGlobalNbrOfParticles,locnPart,.TRUE.)

! Allocate data arrays for mean particle quantities
#if USE_LOADBALANCE
! Arrays might still be allocated from previous loadbalance step
SDEALLOCATE(PartInt)
SDEALLOCATE(PartData)
#endif /*USE_LOADBALANCE*/
ALLOCATE(PartInt( PartIntSize ,offsetElem +1:offsetElem +PP_nElems))
ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
! Allocate data arrays for turbulent particle quantities
IF (ALLOCATED(TurbPartState)) ALLOCATE(TurbPartData(TurbPartDataSize,offsetnPart+1:offsetnPart+locnPart))

! Order the particles along the SFC using a linked list
ALLOCATE(PEM%pStart (offsetElem+1:offsetElem+PP_nElems) , &
         PEM%pNumber(offsetElem+1:offsetElem+PP_nElems) , &
         PEM%pNext  (1           :PDM%maxParticleNumber), &
         PEM%pEnd   (offsetElem+1:offsetElem+PP_nElems))
useLinkedList = .TRUE.
CALL UpdateNextFreePosition()

! Walk along the linked list and fill the data arrays
iPart = offsetnPart
! Walk over all elements on local proc
DO iElem = offsetElem+1,offsetElem+PP_nElems
  ! Set start of particle numbers in current element
  PartInt(1,iElem) = iPart
  ! Find all particles in current element
  IF (ALLOCATED(PEM%pNumber)) THEN
    PartInt(2,iElem) = PartInt(1,iElem) + PEM%pNumber(iElem)
    ! Sum up particles and add properties to output array
    pcount = PEM%pStart(iElem)
    DO iPart = PartInt(1,iElem)+1,PartInt(2,iElem)
      PartData(1:PartDataVarSpecies-1,iPart) = PartState(1:PartDataVarSpecies-1,pcount)
      PartData(PartDataVarSpecies,iPart)     = REAL(PartSpecies(pcount))
      IF (doPartIndex)                                      PartData(PP_nVarPart_loc+2                    ,iPart) = REAL(PartIndex(pcount))
      IF (doParticleReflectionTrack)                        PartData(PartDataVarStart                             ,iPart) = REAL(PartReflCount(pcount))
      IF (doParticleDispersionTrack.OR.doParticlePathTrack) PartData(PartDataVarStart+PartDataVarShift:PartDataVarStart+2+PartDataVarShift,iPart) = PartPath(1:3,pcount)

      ! Turbulent particle properties
      IF (ALLOCATED(TurbPartState))  TurbPartData(:,iPart) = TurbPartState(:,pcount)

      ! Set the index to the next particle
      pcount = PEM%pNext(pcount)
    END DO
    ! Set counter to the end of particle number in the current element
    iPart = PartInt(2,iElem)
  ELSE
    CALL Abort(__STAMP__, " Particle HDF5-Output method not supported! PEM%pNumber not associated")
  END IF
  PartInt(2,iElem)=iPart
END DO ! iElem = offsetElem+1,offsetElem+PP_nElems

! De-allocate linked list and return to normal particle array mode
useLinkedList=.FALSE.
DEALLOCATE( PEM%pStart   &
          , PEM%pNumber  &
          , PEM%pNext    &
          , PEM%pEnd)

END SUBROUTINE FillParticleData


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
unit_index    = 535

IF (MPIRoot) THEN
  INQUIRE(UNIT=unit_index, OPENED=isOpen)
  ! Continue analyze file if opened, otherwise append or create
  IF (.NOT.isOpen) THEN
    outfile = 'PartAnalyze.csv'
    IF (DoRestart .and. FILEEXISTS(outfile)) THEN
      OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
    ELSE
      !--- insert header
      OPEN(unit_index,file=TRIM(outfile))
      ! Time
      WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
      ! Particle number
      IF (CalcPartNumber) THEN
        DO iSpec = 1,nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A14,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPart-Spec-',iSpec,' '
          OutputCounter = OutputCounter + 1
        END DO
     END IF
      ! Particle number balance
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
     END IF
     ! Particle kinetic energy balance
     IF (CalcPartBalance) THEN
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
      ! Instantaneous particle kinetic energy
      IF (CalcEkin) THEN
        DO iSpec=1, nSpecAnalyze
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Ekin-',iSpec,' '
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
  ! Particle number
  IF (CalcPartNumber) THEN
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPart(iSpec))
    END DO
  END IF
  ! Particle number balance
  IF (CalcPartBalance) THEN
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartIn(iSpec))
    END DO
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') REAL(nPartOut(iSpec))
    END DO
  END IF
  ! Particle kinetic energy balance
  IF (CalcPartBalance) THEN
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkinIn(iSpec)
    END DO
    DO iSpec = 1,nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkinOut(iSpec)
    END DO
  END IF
  ! Instantaneous particle kinetic energy
  IF (CalcEkin) THEN
    DO iSpec=1, nSpecAnalyze
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PartEkin(iSpec)
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
USE MOD_Particle_Globals,ONLY: pi
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
                                          MASS_SPHERE(Species(PartSpecies(i))%DensityIC,PartState(PART_DIAM,i)),&
                                          i
END DO
CLOSE(index_unit)
SWRITE(UNIT_stdOut,'(A)')"DONE!"

END SUBROUTINE Visualize_Particles


!===================================================================================================================================
!> Calculate the particle offset and global number of particles across all processors
!> In this routine the number are calculated using integer KIND=8, but are returned with integer KIND=ICC in order to test if using
!> integer KIND=8 is required for total number of particles, particle boundary state, lost particles or clones
!===================================================================================================================================
SUBROUTINE GetOffsetAndGlobalNumberOfParts(CallingRoutine,offsetnPart,globnPart,locnPart,GetMinMaxNbrOfParticles)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Globals
#if USE_MPI
USE MOD_Particle_MPI_Vars ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: CallingRoutine
INTEGER(KIND=IK),INTENT(IN)  :: locnPart
LOGICAL,INTENT(IN)           :: GetMinMaxNbrOfParticles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=IK),INTENT(OUT) :: offsetnPart
INTEGER(KIND=IK),INTENT(OUT) :: globnPart(6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER(KIND=8)              :: globnPart8                         ! always integer KIND=8
INTEGER(KIND=8)              :: locnPart8,locnPart8Recv            ! always integer KIND=8
INTEGER(KIND=IK)             :: SimNumSpecMin,SimNumSpecMax
#endif
!===================================================================================================================================
#if USE_MPI
locnPart8     = INT(locnPart,8)
locnPart8Recv = 0_IK
CALL MPI_EXSCAN(locnPart8,locnPart8Recv,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_FLEXI,iError)
offsetnPart   = INT(locnPart8Recv,KIND=IK)

! Last proc calculates the global number and broadcasts it
IF(myRank.EQ.nProcessors-1) locnPart8=locnPart8Recv+locnPart8
CALL MPI_BCAST(locnPart8,1,MPI_INTEGER8,nProcessors-1,MPI_COMM_FLEXI,iError)

! Global numbers
globnPart8   = locnPart8
LOGWRITE(*,*) TRIM(CallingRoutine)//'offsetnPart,locnPart,globnPart8',offsetnPart,locnPart,globnPart8

! Sanity check: Add up all particles with integer KIND=8 and compare
IF (MPIRoot) THEN
  ! Check if offsetnPart is kind=8 is the number of particles is larger than integer KIND=4
  IF (globnPart8.GT.INT(HUGE(offsetnPart),8)) THEN
    WRITE(UNIT_stdOut,'(A,I0)') '\n\n\nTotal number of particles  : ',globnPart8
    WRITE(UNIT_stdOut,'(A,I0)')       'Maximum number of particles: ',HUGE(offsetnPart)
    CALL Abort(__STAMP__,TRIM(CallingRoutine)//' has encountered more than integer KIND=4 particles!')
  END IF
END IF ! MPIRoot

! Get min/max number of particles
SimNumSpecMin = 0
SimNumSpecMax = 0
IF(GetMinMaxNbrOfParticles)THEN
  IF (PartMPI%MPIRoot) THEN
    CALL MPI_REDUCE(locnPart,SimNumSpecMin,1,MPI_INTEGER_INT_KIND,MPI_MIN,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(locnPart,SimNumSpecMax,1,MPI_INTEGER_INT_KIND,MPI_MAX,0,PartMPI%COMM,IERROR)
  ELSE
    CALL MPI_REDUCE(locnPart,0            ,1,MPI_INTEGER_INT_KIND,MPI_MIN,0,PartMPI%COMM,IERROR)
    CALL MPI_REDUCE(locnPart,0            ,1,MPI_INTEGER_INT_KIND,MPI_MAX,0,PartMPI%COMM,IERROR)
  END IF
END IF ! GetMinMaxNbrOfParticles

! Cast to Kind=IK before returning the number
globnPart(1) = INT(SimNumSpecMin , KIND = IK)
globnPart(2) = INT(SimNumSpecMax , KIND = IK)
globnPart(3) = INT(globnPart8    , KIND = IK)
#else
offsetnPart=0_IK
globnPart(1:3)=INT(locnPart,KIND=IK)

! Suppress compiler warning
NO_OP(CallingRoutine)
#endif

! Get extrema over the complete simulation only during WriteParticleToHDF5
IF(GetMinMaxNbrOfParticles)THEN
  globnPart(4) = MIN(globnPart(1),globnPart(4))
  globnPart(5) = MAX(globnPart(2),globnPart(5))
  globnPart(6) = MAX(globnPart(3),globnPart(6))
END IF ! GetMinMaxNbrOfParticles

END SUBROUTINE GetOffsetAndGlobalNumberOfParts

END MODULE MOD_Particle_Output
