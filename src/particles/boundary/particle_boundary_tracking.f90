!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
#include "particle.h"

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

INTERFACE StoreBoundaryParticleProperties
  MODULE PROCEDURE StoreBoundaryParticleProperties
END INTERFACE

INTERFACE FinalizeParticleBoundaryTracking
  MODULE PROCEDURE FinalizeParticleBoundaryTracking
END INTERFACE

PUBLIC :: DefineParametersParticleBoundaryTracking
PUBLIC :: InitParticleBoundaryTracking
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
USE MOD_ReadInTools              ,ONLY: GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Interpolation_Vars       ,ONLY: InterpolationInitIsDone
USE MOD_Particle_Analyze_Vars    ,ONLY: doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars   ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Vars   ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,ImpactDataSize
USE MOD_Particle_Boundary_Vars   ,ONLY: ImpactTrackInitIsDone,ImpactSideOnProc
USE MOD_Particle_Boundary_Vars   ,ONLY: nSurfTotalSides
USE MOD_Particle_Boundary_Vars   ,ONLY: SurfSide2GlobalSide,SurfSide2GlobalSide_Shared
USE MOD_Particle_Boundary_Vars   ,ONLY: GlobalSide2SurfSide,GlobalSide2SurfSide_Shared
USE MOD_Particle_Boundary_Vars   ,ONLY: SurfBCName,SurfSampleBCs
USE MOD_Particle_Vars            ,ONLY: doPartIndex,doWritePartDiam
#if USE_MPI
USE MOD_MPI_Shared_Vars          ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Boundary_Vars   ,ONLY: SurfSide2GlobalSide_Shared_Win,GlobalSide2SurfSide_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL               :: ImpactSideGlobal
#if !USE_MPI
INTEGER,PARAMETER     :: nImpactProcs = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

IF((.NOT.InterpolationInitIsDone) .OR. ImpactTrackInitIsDone)THEN
   CALL Abort(__STAMP__,"InitializeParticleBoundary not ready to be called or already called.")
END IF

LBWRITE(UNIT_stdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT IMPACT TRACKING...'

! Check if impact tracking is activated
IF (.NOT.doParticleImpactTrack) THEN
  LBWRITE(UNIT_stdOut,'(A)')' INIT IMPACT TRACKING DONE!'
  RETURN
END IF

#if USE_PARTROT
ImpactDataSize = 19
#else
ImpactDataSize = 17
#endif
IF (doWritePartDiam)                                  ImpactDataSize = ImpactDataSize + 2
#if USE_SPHERICITY
ImpactDataSize = ImpactDataSize + 1
#endif
IF (doPartIndex)                                      ImpactDataSize = ImpactDataSize + 1
IF (doParticleDispersionTrack.OR.doParticlePathTrack) ImpactDataSize = ImpactDataSize + 3

! surfaces sides are determined in particle_boundary_sampling.f90!
ImpactSideOnProc = MERGE(.TRUE.,.FALSE.,nSurfTotalSides.NE.0)
#if USE_MPI
CALL MPI_ALLREDUCE(ImpactSideOnProc,ImpactSideGlobal,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_FLEXI,iError)
! Initialize impact tracking communicator (needed for output)
! CALL InitImpactCommunication()
#else
ImpactSideGlobal = ImpactSideOnProc
#endif /*USE_MPI*/

! Free the arrays if they were previously associated
IF (doParticleImpactTrack .AND. .NOT.ImpactSideGlobal) THEN
#if USE_MPI
  ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

  CALL MPI_WIN_UNLOCK_ALL(SurfSide2GlobalSide_Shared_Win,iError)
  CALL MPI_WIN_FREE(      SurfSide2GlobalSide_Shared_Win,iError)
  CALL MPI_WIN_UNLOCK_ALL(GlobalSide2SurfSide_Shared_Win,iError)
  CALL MPI_WIN_FREE(      GlobalSide2SurfSide_Shared_Win,iError)

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

  ! Then, free the pointers or arrays
  MDEALLOCATE(GlobalSide2SurfSide)
  MDEALLOCATE(GlobalSide2SurfSide_Shared)
  MDEALLOCATE(SurfSide2GlobalSide)
  MDEALLOCATE(SurfSide2GlobalSide_Shared)

  ! Finally, free the local arrays
  SDEALLOCATE(SurfBCName)
  SDEALLOCATE(SurfSampleBCs)

  doParticleImpactTrack = .FALSE.
  LBWRITE(UNIT_stdOut,'(A)')' | Impact tracking  requested but no impact faces  found! Disabling...'
  LBWRITE(UNIT_stdOut,'(A)')' INIT IMPACT TRACKING DONE!'
  LBWRITE(UNIT_stdOut,'(132("-"))')

  RETURN
END IF

! Allocate and nullify impact tracking arrays
IF (.NOT.ALLOCATED(PartStateBoundary)) THEN
  ALLOCATE(PartStateBoundary(ImpactDataSize,1:10))
  PartStateBoundary          = 0.
  PartStateBoundaryVecLength = 0
  LBWRITE(UNIT_stdOut,'(A)')' | Starting impact tracking ...'
ELSE
  SWRITE(UNIT_stdOut,'(A)')' | Restarting impact tracking ...'
END IF

LBWRITE(UNIT_stdOut,'(A)')' INIT IMPACT TRACKING DONE!'
LBWRITE(UNIT_stdOut,'(132("-"))')

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
!     WRITE(UNIT_stdOut,'(A,I0,A)') ' Starting impact tracking communications between ', nImpactProcs, ' procs'
! END IF

! CALL MPI_BARRIER(MPI_COMM_FLEXI,iERROR)

! END SUBROUTINE InitImpactCommunication
! #endif /*USE_MPI*/


SUBROUTINE StoreBoundaryParticleProperties(BCSideID,PartID,PartFaceAngle,v_old,PartFaceAngle_old,PartReflCount,alpha,dp_old&
#if USE_PARTROT
    ,rot_old&
#endif
    )
!----------------------------------------------------------------------------------------------------------------------------------!
! Save particle position, velocity, and species to PartStateBoundary container for writing to .h5 later
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Memory,                  ONLY: Allocate_Safe
USE MOD_Particle_Globals,        ONLY: PI
USE MOD_Particle_Analyze_Vars,   ONLY: PartPath,doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars,  ONLY: PartStateBoundary,PartStateBoundaryVecLength,ImpactDataSize
USE MOD_Particle_Vars,           ONLY: Species,PartState,PartSpecies,LastPartPos,PartIndex,doPartIndex,doWritePartDiam
USE MOD_TimeDisc_Vars,           ONLY: t,CurrentStage,dt,RKc,nRKStages
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: PartFaceAngle, v_old(1:3)
REAL,INTENT(IN)                   :: PartFaceAngle_old
REAL,INTENT(IN)                   :: alpha
INTEGER,INTENT(IN)                :: BCSideID,PartID,PartReflCount
REAL,INTENT(IN)                   :: dp_old
#if USE_PARTROT
REAL,INTENT(IN)                   :: rot_old(1:3)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: dims(2),tmp
REAL,ALLOCATABLE     :: PartStateBoundary_tmp(:,:) !                 2nd index: 1 to number of boundary-crossed particles
INTEGER              :: ALLOCSTAT
REAL                 :: e_kin_old,e_kin_new
#if USE_PARTROT
REAL                 :: e_rot_old, e_rot_new
#endif
REAL                 :: t_loc
!===================================================================================================================================

!----  Calculating values before and after reflection
e_kin_old = ENERGY_KINETIC(Species(PartSpecies(PartID))%DensityIC,dp_old,v_old(1:3))
e_kin_new = ENERGY_KINETIC(Species(PartSpecies(PartID))%DensityIC,PartState(PART_DIAM,PartID),PartState(PART_VELV,PartID))
#if USE_PARTROT
e_rot_old = ENERGY_ROTATION(Species(PartSpecies(PartID))%DensityIC,dp_old,rot_old(1:3))
e_rot_new = ENERGY_ROTATION(Species(PartSpecies(PartID))%DensityIC,PartState(PART_DIAM,PartID),PartState(PART_AMOMV,PartID))
#endif

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
    IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary_tmp temporary array!')
    ! Save old data
    PartStateBoundary_tmp(1:ImpactDataSize,1:dims(2)) = PartStateBoundary(1:ImpactDataSize,1:dims(2))

    ! Re-allocate PartStateBoundary to twice the size
    DEALLOCATE(PartStateBoundary)
    CALL Allocate_Safe(PartStateBoundary,(/ImpactDataSize,2*dims(2)/), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary array!')
    PartStateBoundary(1:ImpactDataSize,        1:  dims(2)) = PartStateBoundary_tmp(1:ImpactDataSize,1:dims(2))
    PartStateBoundary(1:ImpactDataSize,dims(2)+1:2*dims(2)) = 0.
  END IF

  ! LastPartPos is set to impact location!

  ! Calculate exact impact time
  IF (CurrentStage.EQ.1) THEN
    t_loc = t                                                       & ! current physical time
          +  RKc(2)                                *dt*alpha          ! relative time till intersection
  ELSE IF (CurrentStage.GT.1 .AND. currentStage.LT.nRKStages) THEN
    t_loc = t                                                       & ! current physical time
          +  RKc(CurrentStage)                     *dt              & ! current stage time
          + (RKc(CurrentStage+1)-RKc(currentStage))*dt*alpha          ! relative time till intersection
  ELSE ! nRKStages
    t_loc = t                                                       & ! current physical time
          +  RKc(CurrentStage)                     *dt              & ! current stage time
          + (1.                 -RKc(currentStage))*dt*alpha          ! relative time till intersection
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
  PartStateBoundary(15:17,iMax) = PartState(PART_VELV,PartID)
  IF (doWritePartDiam) THEN
    PartStateBoundary(18 ,iMax) = dp_old
    PartStateBoundary(19 ,iMax) = PartState(PART_DIAM,PartID)
    tmp = 20
  ELSE
    tmp = 18
  END IF
#if USE_PARTROT
  PartStateBoundary(tmp ,iMax) = e_rot_old
  PartStateBoundary(tmp+1 ,iMax) = e_rot_new
  tmp = tmp+2
#endif
#if USE_SPHERICITY
  PartStateBoundary(tmp ,iMax) = REAL(PartState(PART_SPHE,PartID))
  tmp = tmp+1
#endif
  IF(doPartIndex)                                      PartStateBoundary(tmp                            ,iMax) = PartIndex(    PartID)
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
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,    ONLY :PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!==================================================================================================================================

#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  SDEALLOCATE(PartStateBoundary)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

! #if USE_MPI
! ! Free MPI communicator
! IF (MPI_COMM_IMPACT.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_IMPACT, IERROR)
! #endif /* USE_MPI */

ImpactTrackInitIsDone = .FALSE.

END SUBROUTINE FinalizeParticleBoundaryTracking


END MODULE MOD_Particle_Boundary_Tracking
