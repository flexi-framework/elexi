!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PUEPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
!> Module that provides functions for computing the solutions time history at a defined set of points ("erosionpoints")
!==================================================================================================================================
MODULE MOD_ErosionPoints
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersErosionPoints
  MODULE PROCEDURE DefineParametersErosionPoints
END INTERFACE

INTERFACE InitErosionPoints
  MODULE PROCEDURE InitErosionPoints
END INTERFACE

INTERFACE WriteEP
  MODULE PROCEDURE WriteEP
END INTERFACE

INTERFACE RecordErosionPoint
  MODULE PROCEDURE RecordErosionPoint
END INTERFACE

!INTERFACE RestartErosionPoint
!  MODULE PROCEDURE RestartErosionPoint
!END INTERFACE

INTERFACE FinalizeErosionPoints
  MODULE PROCEDURE FinalizeErosionPoints
END INTERFACE

PUBLIC :: DefineParametersErosionPoints
PUBLIC :: InitErosionPoints
PUBLIC :: FinalizeErosionPoints
PUBLIC :: RecordErosionPoint
!PUBLIC :: RestartErosionPoint
PUBLIC :: WriteEP
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersErosionPoints()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("ErosionPoints")
CALL prms%CreateLogicalOption('Part-TrackImpacts',     "Set true to record individual particle impact data.",                         &
                                                   '.FALSE.')
CALL prms%CreateIntOption(    'Part-TrackImpactsMemory', "Maximum memory in MiB to be used for storing erosionpoint state history. ",   &!//&
!                                                   "If memory is exceeded before regular IO level states are written to file.",&
                                                   '100')
END SUBROUTINE DefineParametersErosionPoints


!==================================================================================================================================
!> Init EP tracking
!==================================================================================================================================
SUBROUTINE InitErosionPoints()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools            ,ONLY: GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Interpolation_Vars     ,ONLY: InterpolationInitIsDone
USE MOD_ErosionPoints_Vars     ,ONLY: doParticleImpactTrack
USE MOD_ErosionPoints_Vars     ,ONLY: EP_Data,EPDataSize,EP_Impacts
USE MOD_ErosionPoints_Vars     ,ONLY: EP_onProc,nEP_Procs,EP_MaxBufferSize
USE MOD_ErosionPoints_Vars,     ONLY: ErosionPointsInitIsDone
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfTotalSides,doParticleImpactSample
USE MOD_Particle_Vars          ,ONLY: doPartIndex
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: EP_maxMemory
INTEGER(KIND=8)       :: EP_maxBufferSize_glob
!==================================================================================================================================

IF((.NOT.InterpolationInitIsDone) .OR. ErosionPointsInitIsDone)THEN
   CALL Abort(__STAMP__,"InitErosionPoints not ready to be called or already called.")
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT IMPACT TRACKING...'

! check if erosionpoints are activated
doParticleImpactTrack = GETLOGICAL('Part-TrackImpacts','.FALSE.')
IF (.NOT.doParticleImpactTrack) THEN
  SWRITE(UNIT_stdOut,'(A)')' INIT EROSIONPOINTS DONE!'
  RETURN
END IF

! surfaces sides are determined in particle_boundary_sampling.f90!
IF (.NOT.doParticleImpactSample) &
 CALL COLLECTIVESTOP(__STAMP__,'Impact tracking only available with Part-SurfaceSampling=T!')

IF (doPartIndex) THEN
  EPDataSize = 15
ELSE
  EPDataSize = 14
END IF

EP_maxMemory     = GETINT('Part-TrackImpactsMemory','100')           ! Max buffer (100MB)
EP_MaxBufferSize = EP_MaxMemory*131072/EPDataSize    != size in bytes/(real*EPDataSize)

IF(nSurfTotalSides.NE.0) THEN
  EP_onProc        = .TRUE.
END IF

! Initialize impact tracking communicator (needed for output)
#if USE_MPI
CALL InitEPCommunicator()
#endif /*USE_MPI*/

! Allocate and nullify impact tracking arrays
ALLOCATE(EP_Data(EPDataSize,EP_MaxBufferSize))
EP_Data    = 0.
EP_Impacts = 0

! This might overflow a kind=4 integer, so use a larger on to be sure
EP_maxBufferSize_glob   = EP_MaxBufferSize*nEP_Procs

ErosionPointsInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A,I0,A,I0,A)') ' | Buffer allocatated for max. ',EP_maxBufferSize_glob, ' impacts (',     &
                                                                      EP_MaxBufferSize,      ' impacts/proc)'
SWRITE(UNIT_stdOut,'(A)')' INIT IMPACT TRACKING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitErosionPoints


#if USE_MPI
!==================================================================================================================================
!> Read EP parameters from ini file and EP definitions from HDF5
!==================================================================================================================================
SUBROUTINE InitEPCommunicator()
! MODULES
USE MOD_Globals
USE MOD_ErosionPoints_Vars   ,ONLY: EP_onProc,myEPrank,EP_COMM,nEP_Procs
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: color
!==================================================================================================================================
color = MPI_UNDEFINED
IF(EP_onProc) color=2

! create new EP communicator for EP communication. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,color,MPI_INFO_NULL,EP_COMM,iError)

! Find my rank on the shared communicator, comm size and proc name
IF (EP_onProc) THEN
  CALL MPI_COMM_RANK(EP_COMM, myEPrank ,iError)
  CALL MPI_COMM_SIZE(EP_COMM, nEP_Procs,iError)

  IF (myEPrank.EQ.0) &
    WRITE(UNIT_StdOut,'(A,I0,A)') ' Starting impact tracking communications between ', nEP_Procs, ' procs'
END IF

CALL MPI_BARRIER(MPI_COMM_FLEXI,iERROR)

END SUBROUTINE InitEPCommunicator
#endif /*USE_MPI*/


SUBROUTINE RecordErosionPoint(BCSideID,PartID,PartFaceAngle,v_old,PartFaceAngle_old,PartReflCount,alpha,n_loc)
!----------------------------------------------------------------------------------------------------------------------------------!
! Combined routine to add erosion impacts to tracking array
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_TimeDisc_Vars,           ONLY: t,CurrentStage,dt,RKc
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Boundary_Vars
USE MOD_ErosionPoints_Vars
USE MOD_Particle_Vars,           ONLY: Species,PartState,PartSpecies,LastPartPos,PartIndex,doPartIndex
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: PartFaceAngle, v_old(1:3)
REAL,INTENT(IN)                   :: PartFaceAngle_old
REAL,INTENT(IN)                   :: alpha
INTEGER,INTENT(IN)                :: BCSideID,PartID,PartReflCount
REAL,INTENT(IN)                   :: n_loc(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_magnitude_old,v_magnitude_new
REAL                              :: e_kin_old,e_kin_new,v_norm(3)!,v_tang(3)
REAL                              :: t_loc
!===================================================================================================================================

!----  Calculating values before and after reflection
v_magnitude_old   = SQRT(DOT_PRODUCT(v_old(1:3),v_old(1:3)))
v_magnitude_new   = SQRT(DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID)))
e_kin_old         = .5*Species(PartSpecies(PartID))%MassIC*v_magnitude_old**2.
e_kin_new         = .5*Species(PartSpecies(PartID))%MassIC*v_magnitude_new**2.
!e_kin_loss        = e_kin_old-e_kin_new
v_norm            = DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc
!v_tang            = PartState(4:6,PartID) - v_norm

!IF ((e_kin_new.GT.e_kin_old).AND.(.NOT.ALMOSTEQUAL(e_kin_new,e_kin_old))) THEN
!    CALL abort(&
!      __STAMP__&
!      ,'Increase in kinetic energy upon reflection! Aborting ...')
!END IF

! LastParPos is set to impact location!

! Calculate exact impact time
IF (CurrentStage.EQ.1) THEN
    t_loc = t                                                         ! current physical time
ELSEIF (CurrentStage.EQ.2) THEN
    t_loc = t                                                       & ! current physical time
          +  RKc(CurrentStage)                     *dt*alpha          ! current stage time
ELSE
    t_loc = t                                                       & ! current physical time
          +  RKc(CurrentStage-1)                   *dt              & ! current stage time
          + (RKc(CurrentStage)-RKc(currentStage-1))*dt*alpha
END IF

! Record individual impact
EP_Impacts = EP_Impacts + 1

EP_Data(1:3,EP_Impacts) = LastPartPos(1:3,PartID)
EP_Data(4:6,EP_Impacts) = v_old(1:3)
EP_Data(7,EP_Impacts)   = REAL(PartSpecies(PartID))
EP_Data(8,EP_Impacts)   = REAL(BCSideID)
EP_Data(9,EP_Impacts)   = t_loc
EP_Data(10,EP_Impacts)  = REAL(PartReflCount)
EP_Data(11,EP_Impacts)  = e_kin_old
EP_Data(12,EP_Impacts)  = e_kin_new
EP_Data(13,EP_Impacts)  = PartFaceAngle_old
EP_Data(14,EP_Impacts)  = PartFaceAngle
IF(doPartIndex) EP_Data(15,EP_Impacts)  = PartIndex(PartID)

END SUBROUTINE RecordErosionPoint


! TODO: A state file only contains the impacts between the outputtime and the last write out. No need to restart this!
!SUBROUTINE RestartErosionPoint
!!==================================================================================================================================
!!> Restarts the impact tracking. Needed before files are flushed
!!==================================================================================================================================
!USE MOD_Globals
!USE MOD_IO_HDF5
!USE MOD_HDF5_Input
!USE MOD_HDF5_Output
!USE MOD_Restart_Vars,               ONLY: RestartFile
!USE MOD_Particle_Boundary_Vars
!USE MOD_Erosionpoints_Vars
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!LOGICAL                        :: ImpactDataExists
!INTEGER                        :: EP_glob
!INTEGER                        :: ImpactDim              !dummy for rank of ImpactData
!!==================================================================================================================================
!!#if USE_MPI
!! Ignore procs without erosion surfaces on them
!!IF(SurfMesh%nSides.EQ.0) RETURN
!!#endif
!SWRITE(UNIT_stdOut,'(A,F0.3,A)')' RESTARTING PARTICLE IMPACT RECORDING FROM HDF5 FILE...'
!
!EP_Impacts = 0
!
!! Open the restart file and search for ImpactData
!IF (MPIRoot) THEN
!  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
!  CALL DatasetExists(File_ID,'ImpactData',ImpactDataExists)
!
!  IF (ImpactDataExists) THEN
!    CALL GetDataSize(File_ID,'ImpactData',ImpactDim,HSize)
!    CHECKSAFEINT(HSize(2),4)
!    EP_glob    = INT(HSize(2))
!    WRITE(UNIT_StdOut,'(A,I8)') ' | Number of impacts:                      ', EP_glob
!    ! We lost the impact <-> proc association, so fill the entire array
!    CALL ReadArray(ArrayName  = 'ImpactData'          , &
!                   rank       = 2                      , &
!                   nVal       = (/EPDataSize,EP_glob/) , &
!                   offset_in  = 0                      , &
!                   offset_dim = 2                      , &
!                   RealArray  = EP_Data(1:EPDataSize,1:EP_glob))
!  END IF
!END IF
!
!! Root has significantly more load here, keep everything in sync
!CALL MPI_BARRIER(MPI_COMM_FLEXI,iERROR)
!
!! Pretend all impacts happened on MPI_ROOT, so we can write out in the next state file write-out
!IF(MPIRoot) THEN
!  EP_Impacts = EP_glob
!  WRITE(UNIT_stdOut,'(A,F0.3,A)')' RESTARTING PARTICLE IMPACT RECORDING FROM HDF5 FILE DONE'
!  WRITE(UNIT_StdOut,'(132("-"))')
!END IF
!
!END SUBROUTINE RestartErosionPoint


!==================================================================================================================================
!> Writes the time history of the solution at the erosionpoints to an HDF5 file
!==================================================================================================================================
SUBROUTINE WriteEP(OutputTime,resetCounters)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Erosionpoints_Vars    ,ONLY: EP_Buffersize
USE MOD_Erosionpoints_Vars    ,ONLY: EPDataSize,offsetEP
USE MOD_Erosionpoints_Vars    ,ONLY: EP_Data,EP_Impacts
USE MOD_HDF5_Output           ,ONLY: WriteAttribute,MarkWriteSuccessfull
USE MOD_HDF5_WriteArray       ,ONLY: WriteArray
USE MOD_IO_HDF5               ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Output_Vars           ,ONLY: ProjectName
USE MOD_Particle_Vars         ,ONLY: doPartIndex
#if USE_MPI
USE MOD_Erosionpoints_Vars    ,ONLY: EP_COMM
USE MOD_Particle_HDF5_output  ,ONLY: DistributedWriteArray
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,   INTENT(IN)             :: OutputTime            !< time
LOGICAL,INTENT(IN)             :: resetCounters         !< flag to reset sample counters and reallocate buffers, once file is done
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: locEP,EP_glob
CHARACTER(LEN=255)             :: FileName,FileString
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL                           :: startT,endT
#if USE_MPI
INTEGER                        :: sendbuf(2),recvbuf(2)
INTEGER                        :: nImpacts(0:nProcessors-1)
#endif
!==================================================================================================================================

! Find amount of recorded impacts on current proc
locEP   = EP_Impacts
EP_glob = 0

! IF(MPIroot)THEN
!  WRITE(UNIT_StdOut,'(132("-"))')
!  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='NO')' WRITE EROSION IMPACTS STATE TO HDF5 FILE...'
!  WRITE(UNIT_stdOut,'(a,I4,a,I4,a)')' EP Buffer  : ',locEP,'/',EP_Buffersize,' impacts.'
  ! GETTIME(startT)
! END IF

!>> Sum up particles from the other procs
#if USE_MPI
sendbuf(1) = locEP
recvbuf    = 0
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
!>> Offset of each proc is the sum of the particles on the previous procs
offsetEP   = recvbuf(1)
sendbuf(1) = recvbuf(1)+locEP
!>> Last proc knows the global number
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER,nProcessors-1,MPI_COMM_FLEXI,iError)
!>> Gather the global number and communicate to root (MPIRank.EQ.0)
EP_glob    = sendbuf(1)
CALL MPI_GATHER(locEP,1,MPI_INTEGER,nImpacts,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
#else
offsetEP   = 0
EP_glob    = locEP
#endif

!>> No impacts recorded
IF (EP_glob.LE.0) THEN
  ! IF(MPIroot)THEN
  !  ! CALL MarkWriteSuccessfull(FileName)
  !  GETTIME(EndT)
  !  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES') 'DONE  [',EndT-StartT,'s]'
  !  WRITE(UNIT_stdOut,'(a,I4,a,I4,a)')' EP Buffer  : ',locEP,' impacts local / ',EP_Buffersize,' impacts global.'
  !  WRITE(UNIT_StdOut,'(132("-"))')
  ! END IF
  RETURN
END IF

EP_Buffersize = EP_glob

! Array for erosion point vars
ALLOCATE(StrVarNames(EPDataSize))
StrVarNames(1) ='ParticlePositionX'
StrVarNames(2) ='ParticlePositionY'
StrVarNames(3) ='ParticlePositionZ'
StrVarNames(4) ='VelocityX'
StrVarNames(5) ='VelocityY'
StrVarNames(6) ='VelocityZ'
StrVarNames(7) ='Species'
StrVarNames(8) ='BoundaryNumber'
StrVarNames(9) ='ImpactTime'
StrVarNames(10)='ReflectionCount'
StrVarNames(11)='E_kin_impact'
StrVarNames(12)='E_kin_reflected'
StrVarNames(13)='Alpha_impact'
StrVarNames(14)='Alpha_reflected'
IF (doPartIndex) StrVarNames(15)='Index'

! Regenerate state file skeleton
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime))
FileString=TRIM(FileName)//'.h5'

IF(MPIRoot)THEN
  CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'VarNamesErosion',EPDataSize,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

#if USE_MPI
CALL DistributedWriteArray(FileString                                    ,&
                           DataSetName  = 'ImpactData'                   ,&
                           rank         = 2                              ,&
                           nValGlobal   = (/EPDataSize  ,EP_glob  /)     ,&
                           nVal         = (/EPDataSize  ,locEP    /)     ,&
                           offset       = (/0           ,offsetEP /)     ,&
                           collective   = .FALSE.                        ,&
                           offSetDim=2                                   ,&
                           communicator = EP_COMM                        ,&
                           RealArray    = EP_Data(1:EPDataSize,1:locEP))
!CALL MPI_BARRIER(PartMPI%COMM,iERROR)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
CALL WriteArray(           DataSetName  = 'ImpactData'                   ,&
                           rank         = 2                              ,&
                           nValGlobal   = (/EPDataSize     ,EP_glob  /)  ,&
                           nVal         = (/EPDataSize     ,locEP    /)  ,&
                           offset       = (/0              ,offsetEP /)  ,&
                           collective   = .TRUE.                         ,&
                           RealArray    = EP_Data(1:EPDataSize,1:locEP))
CALL CloseDataFile()
#endif /*MPI*/

! Deallocate everything
DEALLOCATE(StrVarNames)

! Erase record variables
IF (resetCounters) THEN
  EP_Impacts = 0
  locEP      = 0
  EP_glob    = 0
  EP_Data    = 0.
END IF

! IF(MPIroot)THEN
!  CALL MarkWriteSuccessfull(FileName)
!  GETTIME(EndT)
!  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES') 'DONE  [',EndT-StartT,'s]'
!  WRITE(UNIT_stdOut,'(a,I4,a,I4,a)')' EP Buffer  : ',locEP,' impacts local / ',EP_Buffersize,' impacts global.'
!  WRITE(UNIT_StdOut,'(132("-"))')
! END IF
END SUBROUTINE WriteEP


!==================================================================================================================================
!> Deallocate erosionpoint arrays
!==================================================================================================================================
SUBROUTINE FinalizeErosionPoints()
! MODULES
USE MOD_ErosionPoints_Vars
IMPLICIT NONE
!==================================================================================================================================

SDEALLOCATE(EP_Data)

ErosionPointsInitIsDone = .FALSE.

END SUBROUTINE FinalizeErosionPoints


END MODULE MOD_ErosionPoints
