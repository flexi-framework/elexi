!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
!> Define and init parameters for particle boundaries
!==================================================================================================================================
MODULE MOD_Particle_Boundary_Tracking
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersParticleBoundaryTracking
  MODULE PROCEDURE DefineParametersParticleBoundaryTracking
END INTERFACE

INTERFACE InitParticleBoundaryTracking
  MODULE PROCEDURE InitParticleBoundaryTracking
END INTERFACE

INTERFACE RestartParticleBoundaryTracking
  MODULE PROCEDURE RestartParticleBoundaryTracking
END INTERFACE

INTERFACE StoreBoundaryParticleProperties
  MODULE PROCEDURE StoreBoundaryParticleProperties
END INTERFACE

INTERFACE FinalizeParticleBoundaryTracking
  MODULE PROCEDURE FinalizeParticleBoundaryTracking
END INTERFACE

PUBLIC :: DefineParametersParticleBoundaryTracking
PUBLIC :: InitParticleBoundaryTracking
PUBLIC :: RestartParticleBoundaryTracking
PUBLIC :: StoreBoundaryParticleProperties
PUBLIC :: FinalizeParticleBoundaryTracking
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersParticleBoundaryTracking()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================

CALL prms%SetSection("Particle Impact Tracking")
CALL prms%CreateLogicalOption('Part-TrackImpacts'      , "Set true to record individual particle impact data.",                    &
                                                         '.FALSE.')
CALL prms%CreateIntOption(    'Part-TrackImpactsMemory', "Maximum memory in MiB to be used for storing particle impact history. ", &!//&
!                                                        "If memory is exceeded before regular IO level states are written to file.",&
                                                         '100')

END SUBROUTINE DefineParametersParticleBoundaryTracking


!==================================================================================================================================
!> Init particle impact tracking
!==================================================================================================================================
SUBROUTINE InitParticleBoundaryTracking()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools            ,ONLY: GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Interpolation_Vars     ,ONLY: InterpolationInitIsDone
USE MOD_Particle_Analyze_Vars  ,ONLY: doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,ImpactDataSize
USE MOD_Particle_Boundary_Vars ,ONLY: ImpactTrackInitIsDone,ImpactSideOnProc
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfTotalSides,doParticleImpactSample
USE MOD_Particle_Vars          ,ONLY: doPartIndex
! #if USE_MPI
! USE MOD_Particle_Boundary_Vars ,ONLY: nImpactProcs
! #endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if !USE_MPI
INTEGER,PARAMETER     :: nImpactProcs = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

IF((.NOT.InterpolationInitIsDone) .OR. ImpactTrackInitIsDone)THEN
   CALL Abort(__STAMP__,"InitializeParticleBoundary not ready to be called or already called.")
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT IMPACT TRACKING...'

! Check if impact tracking is activated
doParticleImpactTrack = GETLOGICAL('Part-TrackImpacts','.FALSE.')
IF (.NOT.doParticleImpactTrack) THEN
  SWRITE(UNIT_stdOut,'(A)')' INIT IMPACT TRACKING DONE!'
  RETURN
END IF

! surfaces sides are determined in particle_boundary_sampling.f90!
IF (.NOT.doParticleImpactSample) &
 CALL COLLECTIVESTOP(__STAMP__,'Impact tracking only available with Part-SurfaceSampling=T!')

ImpactDataSize = 14
IF (doPartIndex)                                      ImpactDataSize = ImpactDataSize + 1
IF (doParticleDispersionTrack.OR.doParticlePathTrack) ImpactDataSize = ImpactDataSize + 3

ImpactSideOnProc = MERGE(.TRUE.,.FALSE.,nSurfTotalSides.NE.0)
#if USE_MPI
CALL MPI_ALLREDUCE(ImpactSideOnProc,doParticleImpactTrack,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_FLEXI,iError)
! Initialize impact tracking communicator (needed for output)
! CALL InitImpactCommunication()
#else
doParticleImpactTrack = ImpactSideOnProc
#endif /*USE_MPI*/

! Allocate and nullify impact tracking arrays
ALLOCATE(PartStateBoundary(ImpactDataSize,1:10))
PartStateBoundary          = 0.
PartStateBoundaryVecLength = 0
SWRITE(UNIT_StdOut,'(A)')' | Starting impact tracking ...'

SWRITE(UNIT_stdOut,'(A)')' INIT IMPACT TRACKING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleBoundaryTracking


! #if USE_MPI
! !==================================================================================================================================
! !> Initialize communication for impact tracking
! !==================================================================================================================================
! SUBROUTINE InitImpactCommunication()
! ! MODULES
! USE MOD_Globals
! USE MOD_Particle_Boundary_Vars ,ONLY: ImpactSideOnProc,myImpactRank,nImpactProcs
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT/OUTPUT VARIABLES
! !----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER                   :: color
! !==================================================================================================================================

! color = MERGE(1,MPI_UNDEFINED,ImpactSideOnProc)

! ! create new EP communicator for EP communication. Pass MPI_INFO_NULL as rank to follow the original ordering
! CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,color,MPI_INFO_NULL,MPI_COMM_IMPACT,iError)

! ! Find my rank on the shared communicator, comm size and proc name
! IF (ImpactSideOnProc) THEN
!   CALL MPI_COMM_RANK(EP_COMM, myImpactRank ,iError)
!   CALL MPI_COMM_SIZE(EP_COMM, nImpactProcs, iError)

!   IF (myImpactRank.EQ.0) &
!     WRITE(UNIT_StdOut,'(A,I0,A)') ' Starting impact tracking communications between ', nImpactProcs, ' procs'
! END IF

! CALL MPI_BARRIER(MPI_COMM_FLEXI,iERROR)

! END SUBROUTINE InitImpactCommunication
! #endif /*USE_MPI*/


SUBROUTINE RestartParticleBoundaryTracking
!==================================================================================================================================
!> Restarts the impact tracking. Needed before files are flushed
!==================================================================================================================================
USE MOD_Globals
USE MOD_Particle_Boundary_Vars     ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,ImpactDataSize
USE MOD_HDF5_Input                 ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,GetDataSize
USE MOD_HDF5_Input                 ,ONLY: ReadArray,File_ID,HSize
USE MOD_Particle_Memory            ,ONLY: Allocate_Safe
USE MOD_Restart_Vars               ,ONLY: RestartFile
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: ImpactDataExists
INTEGER                        :: j,offsetImpact
INTEGER                        :: ImpactDim              ! dummy for rank of ImpactData
INTEGER                        :: PartStateBoundaryVecLength_glob
INTEGER                        :: ALLOCSTAT
#if USE_MPI
INTEGER                        :: iProc
INTEGER                        :: offsetImpactsProcCount,offsetImpacts(0:nProcessors)
#endif /* USE_MPI */
!==================================================================================================================================

SWRITE(UNIT_stdOut,'(A,F0.3,A)')' READING PARTICLE IMPACT DATA FROM HDF5 FILE...'

PartStateBoundaryVecLength = 0

! Open the restart file and search for ImpactData
CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,'ImpactData',ImpactDataExists)

IF (ImpactDataExists) THEN
  CALL GetDataSize(File_ID,'ImpactData',ImpactDim,HSize)
  CHECKSAFEINT(HSize(2),4)
  PartStateBoundaryVecLength_glob    = INT(HSize(2))

#if USE_MPI
  ! Distribute impacts between procs
  offsetImpacts              = 0
  PartStateBoundaryVecLength = PartStateBoundaryVecLength_glob/nProcessors
  offsetImpactsProcCount     = PartStateBoundaryVecLength_glob-PartStateBoundaryVecLength*nProcessors
  DO iProc = 0,nProcessors-1
    offsetImpacts(iProc) = PartStateBoundaryVecLength*iProc+MIN(iProc,offsetImpactsProcCount)
  END DO
  offsetImpacts(nProcessors) = PartStateBoundaryVecLength_glob

  ! local impacts and offset
  PartStateBoundaryVecLength = offsetImpacts(myRank+1)-offsetImpacts(myRank)
  offsetImpact               = offsetImpacts(myRank)
#else
  PartStateBoundaryVecLength = PartStateBoundaryVecLength
  offsetImpact               = 0
#endif /* USE_MPI */

  ! Check if PartStateBoundary is sufficiently large
  ASSOCIATE( iMax => PartStateBoundaryVecLength )

    ! Check if array maximum is reached.
    ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
    IF (iMax.GT.10) THEN
      j = 1
      DO WHILE (iMax.GT.j*10)
        j = j*2
      END DO

      ! Check if the new PartStateBoundary size can be safely allocated
      CALL Allocate_Safe(PartStateBoundary,(/ImpactDataSize,j*10/), STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'ERROR in particle_boundary_init.f90: Cannot allocate PartStateBoundary array!')
      PartStateBoundary(:,:) = 0.
    END IF

    ! We lost the impact <-> proc association, so read in according to calculated distribution
    CALL ReadArray(ArrayName  = 'ImpactData'                                  ,&
                   rank       = 2                                             ,&
                   nVal       = (/ImpactDataSize,PartStateBoundaryVecLength/) ,&
                   offset_in  = offsetImpact                                  ,&
                   offset_dim = 2                                             ,&
                   RealArray  = PartStateBoundary(1:ImpactDataSize,1:PartStateBoundaryVecLength))
  END ASSOCIATE
END IF

CALL CloseDataFile()

! Keep everything in sync
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iERROR)
#endif /* USE_MPI */

! Write out in the next state file write-out
IF(MPIRoot) THEN
  WRITE(UNIT_StdOut,'(A,I0,A)')  ' | ',PartStateBoundaryVecLength,' impact read from restart file'
  WRITE(UNIT_StdOut,'(A,F0.3,A)')' READING PARTICLE IMPACT DATA FROM HDF5 FILE DONE'
  WRITE(UNIT_StdOut,'(132("-"))')
END IF

END SUBROUTINE RestartParticleBoundaryTracking


SUBROUTINE StoreBoundaryParticleProperties(BCSideID,PartID,PartFaceAngle,v_old,PartFaceAngle_old,PartReflCount,alpha,n_loc)
!----------------------------------------------------------------------------------------------------------------------------------!
! Save particle position, velocity, and species to PartStateBoundary container for writing to .h5 later
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Analyze_Vars,   ONLY: PartPath,doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars,  ONLY: PartStateBoundary,PartStateBoundaryVecLength,ImpactDataSize
USE MOD_Particle_Memory,         ONLY: Allocate_Safe
USE MOD_Particle_Vars,           ONLY: Species,PartState,PartSpecies,LastPartPos,PartIndex,doPartIndex
USE MOD_TimeDisc_Vars,           ONLY: t,CurrentStage,dt,RKc
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
INTEGER              :: dims(2)
REAL,ALLOCATABLE     :: PartStateBoundary_tmp(:,:) !                 2nd index: 1 to number of boundary-crossed particles
INTEGER              :: ALLOCSTAT
!
REAL                 :: v_magnitude_old,v_magnitude_new
REAL                 :: e_kin_old,e_kin_new,v_norm(3)!,v_tang(3)
REAL                 :: t_loc
!===================================================================================================================================

!----  Calculating values before and after reflection
v_magnitude_old   = SQRT(DOT_PRODUCT(v_old(1:3),v_old(1:3)))
v_magnitude_new   = SQRT(DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID)))
e_kin_old         = .5*Species(PartSpecies(PartID))%MassIC*v_magnitude_old**2.
e_kin_new         = .5*Species(PartSpecies(PartID))%MassIC*v_magnitude_new**2.
!e_kin_loss        = e_kin_old-e_kin_new
v_norm            = DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc
!v_tang            = PartState(4:6,PartID) - v_norm

! Check if PartStateBoundary is sufficiently large
dims = SHAPE(PartStateBoundary)

ASSOCIATE( iMax => PartStateBoundaryVecLength )
  ! Increase maximum number of boundary-impact particles
  iMax = iMax + 1

  ! Check if array maximum is reached.
  ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
  IF (iMax.GT.dims(2)) THEN
    ! Check if the new PartStateBoundary size can be safely allocated
    CALL Allocate_Safe(PartStateBoundary_tmp,(/ImpactDataSize,dims(2)/), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary_tmp temporary array!')
    ! Save old data
    PartStateBoundary_tmp(1:ImpactDataSize,1:dims(2)) = PartStateBoundary(1:ImpactDataSize,1:dims(2))

    ! Re-allocate PartStateBoundary to twice the size
    DEALLOCATE(PartStateBoundary)
    CALL Allocate_Safe(PartStateBoundary,(/ImpactDataSize,2*dims(2)/), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary array!')
    PartStateBoundary(1:ImpactDataSize,        1:  dims(2)) = PartStateBoundary_tmp(1:ImpactDataSize,1:dims(2))
    PartStateBoundary(1:ImpactDataSize,dims(2)+1:2*dims(2)) = 0.
  END IF

  ! LastPartPos is set to impact location!

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
  PartStateBoundary(1:3,iMax) = LastPartPos(1:3,PartID)
  PartStateBoundary(4:6,iMax) = v_old(1:3)
  PartStateBoundary(7  ,iMax) = REAL(PartSpecies(PartID))
  PartStateBoundary(8  ,iMax) = REAL(BCSideID)
  PartStateBoundary(9  ,iMax) = t_loc
  PartStateBoundary(10 ,iMax) = REAL(PartReflCount)
  PartStateBoundary(11 ,iMax) = e_kin_old
  PartStateBoundary(12 ,iMax) = e_kin_new
  PartStateBoundary(13 ,iMax) = PartFaceAngle_old
  PartStateBoundary(14 ,iMax) = PartFaceAngle
  IF(doPartIndex)                                      PartStateBoundary(15                             ,iMax) = PartIndex(    PartID)
  IF(doParticleDispersionTrack.OR.doParticlePathTrack) PartStateBoundary(ImpactDataSize-2:ImpactDataSize,iMax) = PartPath (1:3,PartID)
END ASSOCIATE

END SUBROUTINE StoreBoundaryParticleProperties


!==================================================================================================================================
!> Deallocate impact tracking arrays
!==================================================================================================================================
SUBROUTINE FinalizeParticleBoundaryTracking()
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars
IMPLICIT NONE
!==================================================================================================================================

SDEALLOCATE(PartStateBoundary)

! #if USE_MPI
! ! Free MPI communicator
! IF (MPI_COMM_IMPACT.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_IMPACT, IERROR)
! #endif /* USE_MPI */

ImpactTrackInitIsDone = .FALSE.

END SUBROUTINE FinalizeParticleBoundaryTracking


END MODULE MOD_Particle_Boundary_Tracking