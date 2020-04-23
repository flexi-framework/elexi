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
!> Samples the particle conditions upon encountering a boundary
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Sampling
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersParticleBoundarySampling
  MODULE PROCEDURE DefineParametersParticleBoundarySampling
END INTERFACE

INTERFACE InitParticleBoundarySampling
  MODULE PROCEDURE InitParticleBoundarySampling
END INTERFACE

INTERFACE RestartParticleBoundarySampling
  MODULE PROCEDURE RestartParticleBoundarySampling
END INTERFACE

INTERFACE RecordParticleBoundarySampling
  MODULE PROCEDURE RecordParticleBoundarySampling
END INTERFACE

INTERFACE SideErosion
  MODULE PROCEDURE SideErosion
END INTERFACE

INTERFACE FinalizeParticleBoundarySampling
  MODULE PROCEDURE FinalizeParticleBoundarySampling
END INTERFACE

INTERFACE WriteSurfSample
  MODULE PROCEDURE WriteSurfSample
END INTERFACE

PUBLIC :: DefineParametersParticleBoundarySampling
PUBLIC :: InitParticleBoundarySampling
PUBLIC :: RestartParticleBoundarySampling
PUBLIC :: RecordParticleBoundarySampling
PUBLIC :: SideErosion
PUBLIC :: FinalizeParticleBoundarySampling
PUBLIC :: WriteSurfSample
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle surface sampling
!==================================================================================================================================
SUBROUTINE DefineParametersParticleBoundarySampling()
! MODULES
USE MOD_ReadInTools             ,ONLY: PRMS
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Boundary Sampling")
CALL prms%CreateLogicalOption(  'DoErosion',                      'Flag, if the current calculation requires erosion tracking.')
CALL prms%CreateLogicalOption(  'Part-WriteMacroSurfaceValues',   'Set [T] to activate iteration dependant sampling and h5 output'//&
                                                                  ' surfaces.',                                                    &
                                                                  '.FALSE.')
CALL prms%CreateIntOption(      'Particles-nSurfSample',          'Define polynomial degree of particle BC sampling. Default:'   //&
                                                                  ' NGeo', '1')

END SUBROUTINE DefineParametersParticleBoundarySampling


SUBROUTINE InitParticleBoundarySampling()
!===================================================================================================================================
! Initialization of particle boundary sampling
! default: use for sampling same polynomial degree as NGeo
! 1) all procs identify surfaces for sampling on the node (plus halo region)
! 2) the compute-node leaders communicate the number of sampling surfaces
! 3) every proc holds a copy of the complete sampling region on the node
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Basis                   ,ONLY:LegendreGaussNodesAndWeights
USE MOD_Mesh_Vars               ,ONLY: NGeo,nBCs,BoundaryName
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides,nComputeNodeSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide,SurfSide2GlobalSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSideArea,SurfSampSize
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSampleBCs
USE MOD_Particle_Boundary_Vars  ,ONLY: nErosionVars,doParticleErosionTrack,WriteMacroSurfaceValues
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeElem,nComputeNodeElems
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemSideNodeID_Shared
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TriaTracking
USE MOD_Particle_Vars           ,ONLY: nSpecies
USE MOD_ReadInTools             ,ONLY: GETINT,GETLOGICAL,GETINTARRAY,GETSTR,COUNTOPTION
USE MOD_StringTools             ,ONLY: LowCase
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: nNonUniqueGlobalSides
USE MOD_Particle_MPI_Shared     ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPI_COMM_SHARED,MPIRankLeader,nLeaderGroupProcs
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPI_COMM_LEADERS_SURF,mySurfRank
USE MOD_Particle_MPI_Shared_Vars,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeTotalSides
USE MOD_Particle_MPI_Shared_Vars,ONLY: myLeaderGroupRank,nLeaderGroupProcs
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide_Shared,GlobalSide2SurfSide_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2GlobalSide_Shared,SurfSide2GlobalSide_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSideArea_Shared,SurfSideArea_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState_Shared_Win
USE MOD_Particle_MPI_Boundary_Sampling,ONLY: InitSurfCommunication
#else
USE MOD_Particle_Boundary_Vars  ,ONLY: mySurfRank
USE MOD_Particle_Mesh_Vars      ,ONLY: nTotalSides
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: iBC
INTEGER                                :: iSide,firstSide,lastSide
INTEGER                                :: nSurfSidesProc
INTEGER                                :: offsetSurfTotalSidesProc
INTEGER,ALLOCATABLE                    :: GlobalSide2SurfSideProc(:,:)
INTEGER,ALLOCATABLE                    :: SurfSide2GlobalSideProc(:,:)
! user defined sampling surfaces
LOGICAL                                :: DoSide
INTEGER                                :: iSurfBC,nSurfSampleBC
CHARACTER(20)                          :: tmpStr,tmpStrBC
CHARACTER(LEN=255),ALLOCATABLE         :: BCName(:)
! surface area
INTEGER                                :: SideID,ElemID,LocSideID
INTEGER                                :: p,q,iSample,jSample
INTEGER                                :: TriNum, Node1, Node2
REAL                                   :: area,nVal
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,DIMENSION(:),ALLOCATABLE          :: Xi_NGeo,wGP_NGeo
REAL                                   :: XiOut(1:2),E,F,G,D,tmp1,tmpI2,tmpJ2
REAL                                   :: xNod, zNod, yNod, Vector1(3), Vector2(3), nx, ny, nz
#if USE_MPI
INTEGER                                :: iLeader
INTEGER                                :: offsetSurfSidesProc
INTEGER                                :: GlobalElemID,GlobalElemRank
INTEGER(KIND=MPI_ADDRESS_KIND)         :: MPISharedSize
INTEGER                                :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

! Get input parameters
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING ...'

! Output of macroscopic surface values
!> Double Variable because we follow both FLEXI and PICLas style
doParticleErosionTrack  = GETLOGICAL('DoErosion','F')
WriteMacroSurfaceValues = GETLOGICAL('Part-WriteMacroSurfaceValues','.FALSE.')

! Switch surface macro values flag to .TRUE. for erosion tracking
IF (doParticleErosionTrack) WriteMacroSurfaceValues = .TRUE.

! standard is sampling on NGeo
WRITE(UNIT=tmpStr,FMT='(I0)') NGeo
nSurfSample             = GETINT    ('Particles-nSurfSample',TRIM(tmpStr))

IF((nSurfSample.GT.1).AND.(TriaTracking)) &
  CALL abort(__STAMP__,'nSurfSample cannot be >1 if TriaTracking=T')

! Calculate equidistant surface points
ALLOCATE(XiEQ_SurfSample(0:nSurfSample))
dXiEQ_SurfSample =2./REAL(nSurfSample)
DO q=0,nSurfSample
  XiEQ_SurfSample(q) = dXiEQ_SurfSample * REAL(q) - 1.
END DO

! Allocate shared array for surf sides
#if USE_MPI
MPISharedSize = INT((3*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalSides/),GlobalSide2SurfSide_Shared_Win,GlobalSide2SurfSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,GlobalSide2SurfSide_Shared_Win,IERROR)
GlobalSide2SurfSide => GlobalSide2SurfSide_Shared
#else
ALLOCATE(GlobalSide2SurfSide(1:3,1:nTotalSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  GlobalSide2SurfSide = -1.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(GlobalSide2SurfSide_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

! create boundary name mapping for surfaces SurfaceBC number mapping
nSurfSampleBC  = CountOption('Particles-SurfSampleBC')
ALLOCATE(SurfSampleBCs(nSurfSampleBC))
DO iSurfBC=1,nSurfSampleBC
  SurfSampleBCs(iSurfBC) = TRIM(GETSTR('Particles-SurfSampleBC'))
END DO

! compare BCs against requested BCs
ALLOCATE(BCName(1:nBCs))
nSurfBC = 0
BCName  = ''

DO iBC=1,nBCs
  ! inner side
  IF (PartBound%TargetBoundCond(iBC).EQ.-1) CYCLE

  ! count number of reflective BCs
  IF (PartBound%TargetBoundCond(iBC).EQ.PartBound%ReflectiveBC) THEN
    nSurfBC = nSurfBC + 1
    BCName(nSurfBC) = BoundaryName(iBC)
    CYCLE
  END IF

  ! check if BC is explicitly requested
  CALL lowcase(TRIM(BoundaryName(iBC)),tmpStrBC)
  DO iSurfBC=1,nSurfSampleBC
    IF (tmpStrBC.EQ.SurfSampleBCs(iSurfBC)) THEN
      nSurfBC = nSurfBC + 1
      BCName(nSurfBC) = BoundaryName(iBC)
      EXIT
    END IF
  END DO
END DO

! write back found BC names
IF (nSurfBC.GE.1) THEN
  ALLOCATE(SurfBCName(1:nSurfBC))
  DO iBC=1,nSurfBC
    SurfBCName(iBC) = BCName(iBC)
  END DO
END IF
DEALLOCATE(BCName)

! get number of BC-Sides
#if USE_MPI
! NO HALO REGION REDUCTION
firstSide = INT(REAL( myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
ALLOCATE(GlobalSide2SurfSideProc(1:3,firstSide:lastSide) &
        ,SurfSide2GlobalSideProc(1:3,1         :INT(nNonUniqueGlobalSides/REAL(nComputeNodeProcessors))))
#else
firstSide = 1
lastSide  = nTotalSides
ALLOCATE(GlobalSide2SurfSideProc(1:3,1:nTotalSides) &
        ,SurfSide2GlobalSideProc(1:3,1:nTotalSides))
#endif /*USE_MPI*/

GlobalSide2SurfSideProc    = -1
SurfSide2GlobalSideProc    = -1
nComputeNodeSurfSides      = 0
nSurfSidesProc             = 0

! check every BC side
DO iSide = firstSide,lastSide
  ! ignore non-BC sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  ! ignore sides outside of halo region
  IF (ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.0) CYCLE

  ! check if BC is explicitly requested
  DoSide = .FALSE.
  CALL lowcase(TRIM(BoundaryName(SideInfo_Shared(SIDE_BCID,iSide))),tmpStrBC)
  DO iSurfBC=1,nSurfSampleBC
    IF (tmpStrBC.EQ.SurfSampleBCs(iSurfBC)) THEN
      DoSide = .TRUE.
      EXIT
    END IF
  END DO

  ! count number of reflective and inner BC sides
  IF ((PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,iSide)).EQ.PartBound%ReflectiveBC).OR.DoSide) THEN
    nSurfSidesProc = nSurfSidesProc + 1
    ! check if element for this side is on the current compute-node
    IF ((SideInfo_Shared(SIDE_ELEMID,iSide).GE.offsetComputeNodeElem+1).AND.SideInfo_Shared(SIDE_ELEMID,iSide).LE.offsetComputeNodeElem+nComputeNodeElems) THEN
      nComputeNodeSurfSides  = nComputeNodeSurfSides + 1
    END IF

    ! TODO: Add another check to determine the surface side in halo_eps from current proc. Node-wide halo can become quite large with
    !       with 128 procs!

    ! Write local mapping from Side to Surf side. The rank is already correct, the offset must be corrected by the proc offset later
    GlobalSide2SurfSideProc(SURF_SIDEID,iSide) = nSurfSidesProc
    GlobalSide2SurfSideProc(SURF_RANK  ,iSide) = ElemInfo_Shared(ELEM_RANK,SideInfo_Shared(SIDE_ELEMID,iSide))
#if USE_MPI
    ! get global Elem ID
    GlobalElemID   = SideInfo_Shared(SIDE_ELEMID,iSide)
    GlobalElemRank = ElemInfo_Shared(ELEM_RANK,GlobalElemID)
    ! running on one node, everything belongs to us
    IF (nLeaderGroupProcs.EQ.1) THEN
      GlobalSide2SurfSideProc(SURF_LEADER,iSide) = myLeaderGroupRank
    ELSE
      ! find the compute node
      DO iLeader = 0,nLeaderGroupProcs-2
        ! The last proc is not a leader proc, so catch it separately
        IF ((GlobalElemRank.GE.MPIRankLeader(iLeader)).AND.(GlobalElemRank.LT.MPIRankLeader(iLeader+1))) THEN
          GlobalSide2SurfSideProc(SURF_LEADER,iSide) = iLeader
          EXIT
        ELSE
          GlobalSide2SurfSideProc(SURF_LEADER,iSide) = nLeaderGroupProcs-1
          EXIT
        END IF
      END DO
    END IF
#else
    GlobalSide2SurfSideProc(SURF_LEADER,iSide) = GlobalSide2SurfSideProc(SURF_RANK,iSide)
#endif /*USE_MPI*/

    ! check if element for this side is on the current compute-node. Alternative version to the check above
!    IF (GlobalSide2SurfSideProc(SURF_LEADER,iSide).EQ.myLeaderGroupRank) THEN
!      nComputeNodeSurfSides  = nComputeNodeSurfSides + 1
!    END IF

    SurfSide2GlobalSideProc(SURF_SIDEID,nSurfSidesProc) = iSide
    SurfSide2GlobalSideProc(2:3        ,nSurfSidesProc) = GlobalSide2SurfSideProc(2:3,iSide)
  END IF ! reflective side
END DO

! Find CN global number of total surf sides and write Side to Surf Side mapping into shared array
#if USE_MPI
sendbuf = nSurfSidesProc
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetSurfTotalSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetSurfTotalSidesProc + nSurfSidesProc
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeSurfTotalSides = sendbuf

! Find CN global number of local surf sides and write Side to Surf Side mapping into shared array
sendbuf = nComputeNodeSurfSides
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetSurfSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetSurfSidesProc + nComputeNodeSurfSides
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeSurfSides = sendbuf

! increment SURF_SIDEID by offset
DO iSide = firstSide,lastSide
  IF (GlobalSide2SurfSideProc(SURF_SIDEID,iSide).EQ.-1) CYCLE

  ! sort compute-node local sides first
  IF (GlobalSide2SurfSideProc(SURF_LEADER,iSide).EQ.myLeaderGroupRank) THEN
    GlobalSide2SurfSideProc(SURF_SIDEID,iSide) = GlobalSide2SurfSideProc(SURF_SIDEID,iSide) + offsetSurfSidesProc
    GlobalSide2SurfSide    (:          ,iSide) = GlobalSide2SurfSideProc(:          ,iSide)
  ! sampling sides in halo region follow at the end
  ELSE
    GlobalSide2SurfSideProc(SURF_SIDEID,iSide) = GlobalSide2SurfSideProc(SURF_SIDEID,iSide) + nComputeNodeSurfSides     &
                                                                                            + offsetSurfTotalSidesProc  &
                                                                                            - offsetSurfSidesProc
    GlobalSide2SurfSide    (:          ,iSide) = GlobalSide2SurfSideProc(:          ,iSide)
  END IF
END DO

#else
offsetSurfTotalSidesProc   = 0
nComputeNodeSurfTotalSides = nSurfSidesProc
GlobalSide2SurfSide(:,firstSide:lastSide) = GlobalSide2SurfSide(:,firstSide:lastSide)
#endif /*USE_MPI*/

! Build inverse mapping
#if USE_MPI
MPISharedSize = INT((3*nComputeNodeSurfTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeSurfTotalSides/),SurfSide2GlobalSide_Shared_Win,SurfSide2GlobalSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SurfSide2GlobalSide_Shared_Win,IERROR)
SurfSide2GlobalSide => SurfSide2GlobalSide_Shared

! This barrier MIGHT not be required
CALL MPI_WIN_SYNC(GlobalSide2SurfSide_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SurfSide2GlobalSide_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

SurfSide2GlobalSide(        :,offsetSurfTotalSidesProc+1:offsetSurfTotalSidesProc+nSurfSidesProc) &
  = SurfSide2GlobalSideProc(:,1                         :                         nSurfSidesProc)

CALL MPI_WIN_SYNC(GlobalSide2SurfSide_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SurfSide2GlobalSide_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ALLOCATE(SurfSide2GlobalSide(1:3,1:nComputeNodeSurfTotalSides))
SurfSide2GlobalSide = SurfSide2GlobalSideProc(:,1:nComputeNodeSurfTotalSides)
#endif /*USE_MPI*/

! free temporary arrays
DEALLOCATE(GlobalSide2SurfSideProc)
DEALLOCATE(SurfSide2GlobalSideProc)

! flag if there is at least one surf side on the node (sides in halo region do also count)
SurfOnNode = .FALSE.
IF (nComputeNodeSurfTotalSides.GT.0) SurfOnNode = .TRUE.

! Set size of erosion sampling array
nErosionVars        = 17
IF (nSpecies.EQ.1) THEN
    SurfSampSize   = nErosionVars
ELSE
    SurfSampSize   = nErosionVars*(nSpecies+1)
END IF

!> Leader communication
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  CALL InitSurfCommunication()
END IF
! The leaders are synchronized at this point, but behind the other procs. Bring them back into sync
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
CALL MPI_BCAST(nSurfTotalSides,1,MPI_INTEGER,0,MPI_COMM_SHARED,iError)
#else
mySurfRank      = 0
nSurfTotalSides = nComputeNodeSurfTotalSides
#endif /* USE_MPI */

!> User sanity check
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of sampling sides:               ', nSurfTotalSides

! surface sampling array do not need to be allocated if there are no sides within halo_eps range
IF(.NOT.SurfOnNode) RETURN

! allocate arrays to hold data. If nSpecies > 1, one more species is added to hold the global average
! Caution: SideID is surfSideID, not global ID
!> First local arrays for boundary sampling
!>>> Pointer structure type ruins possibility of ALLREDUCE, need separate arrays for everything
ALLOCATE(SampWallState(1:SurfSampSize,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
SampWallState = 0.

#if USE_MPI
!> Then shared arrays for boundary sampling
MPISharedSize = INT((SurfSampSize*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/SurfSampSize,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SampWallState_Shared_Win,SampWallState_Shared)
CALL MPI_WIN_LOCK_ALL(0,SampWallState_Shared_Win,IERROR)
#else
ALLOCATE(SampWallState_Shared(1:SurfSampSize,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
SampWallState_Shared = 0.
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  SampWallState_Shared = 0.
#if USE_MPI
END IF
CALL MPI_WIN_SYNC(SampWallState_Shared_Win,IERROR)
#endif /*USE_MPI*/

! Surf sides are shared, array calculation can be distributed
#if USE_MPI
MPISharedSize = INT((nSurfSample*nSurfSample*nComputeNodeSurfTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SurfSideArea_Shared_Win,SurfSideArea_Shared)
CALL MPI_WIN_LOCK_ALL(0,SurfSideArea_Shared_Win,IERROR)
SurfSideArea => SurfSideArea_Shared

firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(SurfSideArea(1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))

firstSide = 1
lastSide  = nSurfTotalSides
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  SurfSideArea=0.
#if USE_MPI
END IF
CALL MPI_WIN_SYNC(SurfSideArea_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/

! get interpolation points and weights
ALLOCATE( Xi_NGeo( 0:NGeo)  &
        , wGP_NGeo(0:NGeo) )
CALL LegendreGaussNodesAndWeights(NGeo,Xi_NGeo,wGP_NGeo)

! compute area of sub-faces
tmp1=dXiEQ_SurfSample/2.0 !(b-a)/2

DO iSide = firstSide,LastSide
  ! get global SideID. This contains only nonUniqueSide, no special mortar treatment required
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)

  IF (TriaTracking) THEN
    ElemID    = SideInfo_Shared(SIDE_ELEMID ,SideID)
    LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
    area = 0.
    xNod = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,LocSideID,ElemID)+1)
    yNod = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,LocSideID,ElemID)+1)
    zNod = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,LocSideID,ElemID)+1)

    DO TriNum = 1,2
      Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
      Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle
        Vector1(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node1,LocSideID,ElemID)+1) - xNod
        Vector1(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node1,LocSideID,ElemID)+1) - yNod
        Vector1(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node1,LocSideID,ElemID)+1) - zNod
        Vector2(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node2,LocSideID,ElemID)+1) - xNod
        Vector2(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node2,LocSideID,ElemID)+1) - yNod
        Vector2(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node2,LocSideID,ElemID)+1) - zNod
      nx = - Vector1(2) * Vector2(3) + Vector1(3) * Vector2(2) !NV (inwards)
      ny = - Vector1(3) * Vector2(1) + Vector1(1) * Vector2(3)
      nz = - Vector1(1) * Vector2(2) + Vector1(2) * Vector2(1)
      nVal = SQRT(nx*nx + ny*ny + nz*nz)
        area = area + nVal/2.
    END DO
    SurfSideArea(1,1,iSide) = area
  ! .NOT. TriaTracking
  ELSE
    ! call here stephens algorithm to compute area
    DO jSample=1,nSurfSample
      DO iSample=1,nSurfSample
        area=0.
        tmpI2=(XiEQ_SurfSample(iSample-1)+XiEQ_SurfSample(iSample))/2. ! (a+b)/2
        tmpJ2=(XiEQ_SurfSample(jSample-1)+XiEQ_SurfSample(jSample))/2. ! (a+b)/2
        DO q=0,NGeo
          DO p=0,NGeo
            XiOut(1)=tmp1*Xi_NGeo(p)+tmpI2
            XiOut(2)=tmp1*Xi_NGeo(q)+tmpJ2
            CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
                                                    ,Gradient=gradXiEta3D)
            ! calculate first fundamental form
            E=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
            F=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
            G=DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
            D=SQRT(E*G-F*F)
            area = area+tmp1*tmp1*D*wGP_NGeo(p)*wGP_NGeo(q)
          END DO
        END DO
        SurfSideArea(iSample,jSample,iSide) = area
      END DO ! iSample=1,nSurfSample
    END DO ! jSample=1,nSurfSample
  END IF
END DO ! iSide = firstSide,lastSide

#if USE_MPI
CALL MPI_WIN_SYNC(SurfSideArea_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/

! get the full area of all surface sides
area=0.

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  DO iSide = 1,nComputeNodeSurfTotalSides
#if USE_MPI
    ! ignore sides in halo region
    IF (SurfSide2GlobalSide(SURF_LEADER,iSide).NE.myLeaderGroupRank) CYCLE
#endif /*USE_MPI*/

    area = area + SUM(SurfSideArea(:,:,iSide))
  END DO ! iSide = 1,nComputeNodeSurfTotalSides
#if USE_MPI
  ! surf leaders communicate total surface area
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,area,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SURF,IERROR)
END IF

! surf leaders inform the other procs on their node
CALL MPI_BCAST(area,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! de-allocate temporary interpolation points and weights
DEALLOCATE(Xi_NGeo,wGP_NGeo)

#if USE_MPI
IF (mySurfRank.EQ.0) THEN
#endif
  WRITE(UNIT_StdOut,'(A,ES10.4E2)') ' | Surface-Area:                         ', Area
  WRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING DONE'
#if USE_MPI
END IF
#endif

END SUBROUTINE InitParticleBoundarySampling


SUBROUTINE SideErosion(PartTrajectory,n_loc,xi,eta,PartID,SideID,alpha)
!----------------------------------------------------------------------------------------------------------------------------------!
! Tracks erosion on designated sides other than reflective wall
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
!USE MOD_ErosionPoints,          ONLY:RecordErosionPoint
!USE MOD_ErosionPoints_Vars,     ONLY:EP_inUse
!USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Tracking_Vars, ONLY:TrackingMethod
!USE MOD_Particle_Vars,          ONLY:PartReflCount
USE MOD_Particle_Vars,          ONLY:PartState
!USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: PartTrajectory(1:3)
REAL,INTENT(IN)                   :: n_loc(1:3)
REAL,INTENT(IN)                   :: xi,eta,alpha
INTEGER,INTENT(IN)                :: PartID,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3)
REAL                              :: Xitild,EtaTild
INTEGER                           :: p,q,SurfSideID
REAL                              :: PartFaceAngle
!===================================================================================================================================

IF (WriteMacroSurfaceValues) THEN
  ! Find correct boundary on SurfMesh
  SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

  ! Not a sampling surface (e.g. call for outlet without sampling)
  IF (SurfSideID.EQ.-1) RETURN

  ! Make sure we have to old velocity safe
  v_old = PartState(4:6,PartID)

  ! compute p and q for supersampling
  SELECT CASE(TrackingMethod)
    CASE(REFMAPPING,TRACING)
      Xitild  = MIN(MAX(-1.,xi ),0.99)
      Etatild = MIN(MAX(-1.,eta),0.99)
      p = INT((Xitild +1.0)/dXiEQ_SurfSample)+1
      q = INT((Etatild+1.0)/dXiEQ_SurfSample)+1
    CASE(TRIATRACKING)
  END SELECT

  PartFaceAngle=ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

  CALL RecordParticleBoundarySampling(PartID                                              &
                                     ,SurfSideID                                          &
                                     ,p                                                   &
                                     ,q                                                   &
                                     ,v_old                                               &
                                     ,PartFaceAngle)
END IF

!IF (EP_inUse) CALL RecordErosionPoint(BCSideID        = BC(SideID),                       &
!                                      PartID          = PartID,                           &
!                                      PartFaceAngle   = PartFaceAngle,                    &
!                                      v_old           = v_old,                            &
!                                      PartFaceAngle_old = PartFaceAngle,                  &
!                                      PartReflCount   = PartReflCount(PartID),            &
!                                      alpha           = alpha)

END SUBROUTINE SideErosion


SUBROUTINE RecordParticleBoundarySampling(PartID,SurfSideID,p,q,v_old,PartFaceAngle)
!----------------------------------------------------------------------------------------------------------------------------------!
! Combined routine to add calculated erosion variables to tracking array
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Vars,          ONLY:Species,PartSpecies,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: PartFaceAngle, v_old(1:3)
INTEGER,INTENT(IN)                :: PartID, SurfSideID, p, q
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: delta                          ! Reusable variable for variance calculation
REAL                              :: delta2                         ! Reusable variable for variance calculation
INTEGER                           :: nShift                         ! Shift amount for species tracking
REAL                              :: v_magnitude
REAL                              :: e_kin
!===================================================================================================================================
nShift        = PartSpecies(PartID) * nErosionVars

!===================================================================================================================================
! SAMP WALL
!===================================================================================================================================

!  Sampling kinetic energy at walls
v_magnitude   = SQRT(DOT_PRODUCT(v_old(1:3),v_old(1:3)))
e_kin         = .5*Species(PartSpecies(PartID))%MassIC*v_magnitude**2.

! All Variables are saved DOUBLE. First Total, then per SPECIES
!-- 1. - .. / Impact Counter
    SampWallState(1,p,q,surfSideID)                       = SampWallState(1,p,q,surfSideID) + 1
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN
    SampWallState(1+nShift,p,q,surfSideID)                = SampWallState(1+nShift,p,q,surfSideID) + 1
END IF

!===================================================================================================================================
!---- 2. - 6. / Kinetic energy on impact (mean, min, max, M2, variance)
!<<< Welford's algorithm for variance
!    (count, mean, M2) = existingAggregate
!    count = count + 1
!    delta = newValue - mean
!    mean = mean + delta / count
!    delta2 = newValue - mean
!    M2 = M2 + delta * delta2

!   Record first impact, otherwise min will be frozen at zero
    IF (SampWallState(1,p,q,surfSideID).EQ.1) THEN
!        SampWallState(2,p,q,surfSideID)                   = e_kin
        SampWallState(3,p,q,surfSideID)                   = e_kin
        SampWallState(4,p,q,surfSideID)                   = e_kin
    END IF

!   All subsequent impacts
    delta                                                   = e_kin - SampWallState(2,p,q,surfSideID)
!   Update mean
    SampWallState(2,p,q,surfSideID)                      = SampWallState(2,p,q,surfSideID) + delta /                          &
                                                              SampWallState(1,p,q,surfSideID)
!    Find min/max of distribution
    IF (e_kin.LT.SampWallState(3,p,q,surfSideID)) THEN
        SampWallState(3,p,q,surfSideID)               = e_kin
    END IF
    IF (e_kin.GT.SampWallState(4,p,q,surfSideID)) THEN
        SampWallState(4,p,q,surfSideID)               = e_kin
    END IF
!   delta2 = newValue - mean
    delta2                                                  = e_kin - SampWallState(2,p,q,surfSideID)
!   M2 = M2 + delta * delta2
    SampWallState(5,p,q,surfSideID)                       = SampWallState(5,p,q,surfSideID) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWallState(6,p,q,surfSideID)                       = SampWallState(5,p,q,surfSideID)/SampWallState(1,p,q,surfSideID)
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN

!   Record first impact, otherwise min will be frozen at zero
    IF (SampWallState(1+nShift,p,q,surfSideID).EQ.1) THEN
!        SampWallState(2+nShift,p,q,surfSideID)            = e_kin
        SampWallState(3+nShift,p,q,surfSideID)            = e_kin
        SampWallState(4+nShift,p,q,surfSideID)            = e_kin
    END IF

!   All subsequent impacts
    delta                                                   = e_kin - SampWallState(2+nShift,p,q,surfSideID)
!   Update mean
    SampWallState(2+nShift,p,q,surfSideID)                = SampWallState(2+nShift,p,q,surfSideID) + delta /                   &
                                                              SampWallState(1+nShift,p,q,surfSideID)
!   Find min/max of distribution
    IF (e_kin.LT.SampWallState(3+nShift,p,q,surfSideID)) THEN
        SampWallState(3+nShift,p,q,surfSideID)        = e_kin
    END IF
    IF (e_kin.GT.SampWallState(4+nShift,p,q,surfSideID)) THEN
        SampWallState(4+nShift,p,q,surfSideID)        = e_kin
    END IF
!   delta2 = newValue - mean
    delta2                                                  = e_kin - SampWallState(2+nShift,p,q,surfSideID)
!   M2 = M2 + delta * delta2
    SampWallState(5+nShift,p,q,surfSideID)                = SampWallState(5+nShift,p,q,surfSideID) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWallState(6+nShift,p,q,surfSideID)                = SampWallState(5+nShift,p,q,surfSideID) /                           &
                                                              SampWallState(1+nShift,p,q,surfSideID)
END IF

!-- 7. - 11 / Impact angle (mean, min, max, M2, variance)
!   Record first impact, otherwise min will be frozen at zero
    IF (SampWallState(1,p,q,surfSideID).EQ.1) THEN
!        SampWallState(7,p,q,surfSideID)                   = PartFaceAngle
        SampWallState(8,p,q,surfSideID)                   = PartFaceAngle
        SampWallState(9,p,q,surfSideID)                   = PartFaceAngle
    END IF

!   All subsequent impacts
    delta                                                   = PartFaceAngle - SampWallState(7,p,q,surfSideID)
!   Update mean
    SampWallState(7,p,q,surfSideID)                       = SampWallState(7,p,q,surfSideID) + delta /                          &
                                                              SampWallState(1,p,q,surfSideID)
!    Find min/max of distribution
    IF (PartFaceAngle.LT.SampWallState(8,p,q,surfSideID)) THEN
        SampWallState(8,p,q,surfSideID)                = PartFaceAngle
    END IF
    IF (PartFaceAngle.GT.SampWallState(9,p,q,surfSideID)) THEN
        SampWallState(9,p,q,surfSideID)                = PartFaceAngle
    END IF
!   delta2 = newValue - mean
    delta2                                                  = PartFaceAngle - SampWallState(7,p,q,surfSideID)
!   M2 = M2 + delta * delta2
    SampWallState(10,p,q,surfSideID)                      = SampWallState(10,p,q,surfSideID) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWallState(11,p,q,surfSideID)                      = SampWallState(10,p,q,surfSideID)/SampWallState(1,p,q,surfSideID)
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN

!   Record first impact, otherwise min will be frozen at zero
    IF (SampWallState(1+nShift,p,q,surfSideID).EQ.1) THEN
!        SampWallState(7+nShift,p,q,surfSideID)            = PartFaceAngle
        SampWallState(8+nShift,p,q,surfSideID)            = PartFaceAngle
        SampWallState(9+nShift,p,q,surfSideID)            = PartFaceAngle
    END IF

!   All subsequent impacts
    delta                                                   = PartFaceAngle - SampWallState(7+nShift,p,q,surfSideID)
!   Update mean
    SampWallState(7+nShift,p,q,surfSideID)                = SampWallState(7+nShift,p,q,surfSideID) + delta /                   &
                                                              SampWallState(1+nShift,p,q,surfSideID)
!   Find min/max of distribution
    IF (PartFaceAngle.LT.SampWallState(8+nShift,p,q,surfSideID)) THEN
        SampWallState(8+nShift,p,q,surfSideID)        = PartFaceAngle
    END IF
    IF (PartFaceAngle.GT.SampWallState(9+nShift,p,q,surfSideID)) THEN
        SampWallState(9+nShift,p,q,surfSideID)        = PartFaceAngle
    END IF
!   delta2 = newValue - mean
    delta2                                                  = PartFaceAngle - SampWallState(7+nShift,p,q,surfSideID)
!   M2 = M2 + delta * delta2
    SampWallState(10+nShift,p,q,surfSideID)               = SampWallState(10+nShift,p,q,surfSideID) + delta * delta2
!   Update sample variance here so we can check values at runtime, variance       = M2/count
!                                                                  sampleVariance = M2/(count - 1)
    SampWallState(11+nShift,p,q,surfSideID)               = SampWallState(10+nShift,p,q,surfSideID) /                          &
                                                              SampWallState(1+nShift,p,q,surfSideID)
END IF

!-- 12 - 14 / Sampling Current Forces at walls
    SampWallState(12,p,q,surfSideID)= SampWallState(12,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(1))
    SampWallState(13,p,q,surfSideID)= SampWallState(13,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(2))
    SampWallState(14,p,q,surfSideID)= SampWallState(14,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(3))
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN
    SampWallState(12+nShift,p,q,surfSideID)= SampWallState(12+nShift,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC     &
                                               * (v_old(1))
    SampWallState(13+nShift,p,q,surfSideID)= SampWallState(13+nShift,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC     &
                                               * (v_old(2))
    SampWallState(14+nShift,p,q,surfSideID)= SampWallState(14+nShift,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC     &
                                               * (v_old(3))
END IF

!-- 15 - 17 / Sampling Average Forces at walls
    SampWallState(15,p,q,surfSideID)= SampWallState(15,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(1))
    SampWallState(16,p,q,surfSideID)= SampWallState(16,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(2))
    SampWallState(17,p,q,surfSideID)= SampWallState(17,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC                   &
                                        * (v_old(3))
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
IF (nSpecies.GT.1) THEN
    SampWallState(15+nShift,p,q,surfSideID)= SampWallState(15+nShift,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC     &
                                        * (v_old(1))
    SampWallState(16+nShift,p,q,surfSideID)= SampWallState(16+nShift,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC     &
                                        * (v_old(2))
    SampWallState(17+nShift,p,q,surfSideID)= SampWallState(17+nShift,p,q,surfSideID) + Species(PartSpecies(PartID))%MassIC     &
                                        * (v_old(3))
END IF

END SUBROUTINE RecordParticleBoundarySampling


SUBROUTINE RestartParticleBoundarySampling(remap_opt)
!===================================================================================================================================
!> Restart the particle boundary sampling (assuming RestartData array exists)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Input
USE MOD_HDF5_Output
USE MOD_Output_Vars                ,ONLY: ProjectName
USE MOD_Restart_Vars               ,ONLY: RestartTime
USE MOD_Particle_Vars              ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfOnNode
USE MOD_Particle_Boundary_Vars     ,ONLY: nSurfSample,offsetComputeNodeSurfSide,nComputeNodeSurfSides
USE MOD_Particle_Boundary_Vars,     ONLY: SampWallState_Shared
USE MOD_Particle_Boundary_Vars
#if USE_MPI
USE MOD_Particle_Boundary_Vars,     ONLY: SampWallState_Shared_Win
USE MOD_Particle_MPI_Shared_Vars   ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_Particle_MPI_Shared_Vars   ,ONLY: myComputeNodeRank,mySurfRank
#else
USE MOD_Particle_Boundary_vars     ,ONLY: mySurfRank
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: remap_opt       ! routine was called from posti. Change input file name
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString
CHARACTER(LEN=255)                  :: H5_Name
LOGICAL                             :: ErosionDataExists
INTEGER                             :: RestartVarCount,NRestartFile
REAL,ALLOCATABLE                    :: RestartArray(:,:,:,:)
!===================================================================================================================================

#if USE_MPI
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfOnNode) RETURN

! surface sampling is not restarted if not required in current run
IF(.NOT.WriteMacroSurfaceValues) THEN
    ErosionRestart = .FALSE.
    RETURN
END IF

IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  IF (mySurfRank.EQ.0) WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' RESTARTING EROSION TRACKING FROM HDF5 FILE...'

  ! Get number of required erosion vars
  IF (nSpecies.EQ.1) THEN
      RestartVarCount =  nErosionVars
  ELSE
      RestartVarCount = (nErosionVars*(nSpecies+1))
  END IF

  ! Allocate array for restart on ALL procs
  IF (nSpecies.EQ.1) THEN
    ALLOCATE(RestartArray (nErosionVars              ,1:nSurfSample,1:nSurfSample,nComputeNodeSurfSides))
  ELSE
    ALLOCATE(RestartArray((nErosionVars)*(nSpecies+1),1:nSurfSample,1:nSurfSample,nComputeNodeSurfSides))
  END IF
  RestartArray = 0.

  ! Open restart array in general erosion file
  IF (PRESENT(remap_opt)) THEN
      FileString=remap_opt
  ELSE
      FileName=TIMESTAMP(TRIM(ProjectName)//'_ErosionSurfState',RestartTime)
      FileString=TRIM(FileName)//'.h5'
  END IF
#if USE_MPI
END IF
#endif /*USE_MPI*/

! Return if the file was not found. This can happen quite often, so it's generally no reason to abort
IF ((mySurfRank.EQ.0).AND..NOT.FILEEXISTS(FileString)) THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' | Erosion File does not exist for current time. Aborting erosion restart ...'
  WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' RESTARTING EROSION TRACKING FROM HDF5 FILE DONE'
  WRITE(UNIT_StdOut,'(132("-"))')
ELSE
  ! Set flag to take restart time into account
  ErosionRestart = .TRUE.
END IF

! Only surf root knows about this so far, first inform the other surf leaders. They inform every proc
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_BCAST(ErosionRestart,1,MPI_LOGICAL,0,MPI_COMM_LEADERS_SURF,iError)
END IF
CALL MPI_BCAST(ErosionRestart,1,MPI_LOGICAL,0,MPI_COMM_SHARED      ,iError)
#endif
IF (.NOT.ErosionRestart) RETURN

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*USE_MPI*/

  ! Check if erosion file has restart data
  CALL DatasetExists(File_ID,'RestartData',ErosionDataExists)

  IF (ErosionDataExists) THEN
    ! Check if size of restart array is correct
    CALL ReadAttribute(File_ID,'NRestart',1,IntScalar=NRestartFile)
    IF (NRestartFile.NE.RestartVarCount) THEN
      IF (mySurfRank.EQ.0) THEN
        WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' | Number of variables in erosion file does not match. Aborting erosion restart ...'
        WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' RESTARTING EROSION TRACKING FROM HDF5 FILE DONE'
        WRITE(UNIT_StdOut,'(132("-"))')
      END IF
      ErosionRestart = .FALSE.
      CALL CloseDataFile()
    END IF
  ELSE
    IF (mySurfRank.EQ.0) THEN
      WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' | RestartData does not exists in erosion tracking file.'
      WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' RESTARTING EROSION TRACKING FROM HDF5 FILE DONE'
      WRITE(UNIT_StdOut,'(132("-"))')
    END IF
    ErosionRestart = .FALSE.
    CALL CloseDataFile()
  END IF
#if USE_MPI
END IF
#endif /*USE_MPI*/

! Only surf leader knows about this so far, inform other procs
#if USE_MPI
CALL MPI_BCAST(ErosionRestart,1,MPI_LOGICAL,0,MPI_COMM_SHARED      ,iError)
#endif
IF (.NOT.ErosionRestart) RETURN

! Dataset is valid and still open on node roots
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  ASSOCIATE (&
    nLocalSides          => nComputeNodeSurfSides     , &
    offsetSurfSide       => offsetComputeNodeSurfSide)

  ! Read array from restart file
  WRITE(H5_Name,'(A)') 'RestartData'
  CALL ReadArray(ArrayName  = H5_Name                                                           , &
                 rank       = 4                                                                 , &
                 nVal       = (/RestartVarCount,nSurfSample,nSurfSample,nLocalSides/)           , &
                 offset_in  = offsetSurfSide                                                    , &
                 offset_dim = 4                                                                 , &
                 RealArray  = RestartArray)

  END ASSOCIATE

  ! We got the array, close the file
  CALL CloseDataFile()

  ! Only loop over sides on current proc (meaningless in single core case)
  SampWallState_Shared(:,:,:,1:nComputeNodeSurfSides) = RestartArray(:,:,:,:)

  IF (mySurfRank.EQ.0) THEN
    WRITE(UNIT_stdOut,'(a,E11.3)',ADVANCE='YES')' | Erosion tracking restart successful at t =',RestartTime
    WRITE(UNIT_stdOut,'(a)'      ,ADVANCE='YES')' RESTARTING EROSION TRACKING FROM HDF5 FILE DONE'
    WRITE(UNIT_StdOut,'(132("-"))')
  END IF

  DEALLOCATE(RestartArray)
#if USE_MPI
END IF ! myComputeNodeRank.EQ.0
#endif /*USE_MPI*/

! Synchronize data on the compute node
#if USE_MPI
CALL MPI_WIN_SYNC(SampWallState_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

END SUBROUTINE RestartParticleBoundarySampling


SUBROUTINE WriteSurfSample(MeshFileName,OutputTime,remap_opt)
!===================================================================================================================================
!> write the final values of the surface sampling to a HDF5 state file
!> additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Output             ,ONLY: WriteHeader,WriteAttribute
USE MOD_HDF5_WriteArray         ,ONLY: WriteArray
USE MOD_Output_Vars             ,ONLY: ProjectName
USE MOD_Particle_Vars           ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample,offSetSurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState_Shared
USE MOD_Particle_boundary_Vars  ,ONLY: nComputeNodeSurfSides,offsetComputeNodeSurfSide
USE MOD_Particle_Boundary_Vars
#if USE_MPI
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPI_COMM_LEADERS_SURF,mySurfRank
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: remap_opt       !routine was called from posti. Change output file name
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString,Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255)                  :: SpecID
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: nVar2D, nVar2D_Spec, nVar2D_Total, nVarCount, iSpec, nShift
INTEGER                             :: nShiftRHS
!INTEGER                             :: iSurfSide
INTEGER                             :: OutputVarCount
REAL,ALLOCATABLE                    :: OutputArray(:,:,:,:)
REAL                                :: startT,endT
!===================================================================================================================================

#if USE_MPI
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nSurfTotalSides      .EQ.0) RETURN
#endif /*USE_MPI*/

IF (mySurfRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE EROSION SURFACE STATE TO HDF5 FILE...'
  GETTIME(startT)
END IF

! Allocate and fill output array
IF (nSpecies.EQ.1) THEN
    OutputVarCount = nErosionVars
ELSE
    OutputVarCount = (nErosionVars*(nSpecies+1))
END IF
ALLOCATE(OutputArray (OutputVarCount,1:nSurfSample,1:nSurfSample,nComputeNodeSurfSides))
OutputArray(:,:,:,1:nComputeNodeSurfSides) = SampWallState_Shared(:,:,:,1:nComputeNodeSurfSides)

IF (PRESENT(remap_opt)) THEN
  FileName=TIMESTAMP(TRIM(ProjectName)//remap_opt,OutputTime)
ELSE
  FileName=TIMESTAMP(TRIM(ProjectName)//'_ErosionSurfState',OutputTime)
END IF
FileString=TRIM(FileName)//'.h5'

! Create dataset attribute "SurfVarNames"
IF (nSpecies.EQ.1) THEN
    nVar2D = nErosionVars - 1
ELSE
    nVar2D = (nErosionVars - 1) * (nSpecies+1)
END IF
nVar2D_Spec=1
nVar2D_Total = nVar2D + nVar2D_Spec*nSpecies

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF (mySurfRank.EQ.0) THEN
#endif
  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'DSMCSurfState'

  ! Write file header
  CALL WriteHeader(Statedummy,File_ID)
  CALL WriteAttribute(File_ID,'DSMC_nSurfSample',1       ,IntScalar =nSurfSample)
  CALL WriteAttribute(File_ID,'DSMC_nSpecies'   ,1       ,IntScalar =nSpecies)
  CALL WriteAttribute(File_ID,'MeshFile'        ,1       ,StrScalar =(/MeshFileName/))
  CALL WriteAttribute(File_ID,'Time'            ,1       ,RealScalar=OutputTime)
  CALL WriteAttribute(File_ID,'BC_Surf'         ,nSurfBC ,StrArray  =SurfBCName)
  CALL WriteAttribute(File_ID,'N'               ,1       ,IntScalar =nSurfSample)
  NodeTypeTemp='VISU'
  CALL WriteAttribute(File_ID,'NodeType'        ,1       ,StrScalar =(/NodeTypeTemp/))

  ALLOCATE(Str2DVarNames(1:nVar2D_Total))
  nVarCount=0

  DO iSpec=1,nSpecies
    WRITE(SpecID,'(I3.3)') iSpec
    Str2DVarNames(nVarCount+1) ='Spec'//TRIM(SpecID)//'_Counter'
    nVarCount=nVarCount+nVar2D_Spec
  END DO ! iSpec=1,nSpecies

  ! fill varnames for total values
  Str2DVarNames(nVarCount+1) ='Impacts'
  Str2DVarNames(nVarCount+2) ='ImpactsPerArea'
  Str2DVarNames(nVarCount+3) ='E_kin(mean)'
  Str2DVarNames(nVarCount+4) ='E_kin(min)'
  Str2DVarNames(nVarCount+5) ='E_kin(max)'
  Str2DVarNames(nVarCount+6) ='E_kin(variance)'
  Str2DVarNames(nVarCount+7) ='ImpactAngle(mean)'
  Str2DVarNames(nVarCount+8) ='ImpactAngle(min)'
  Str2DVarNames(nVarCount+9) ='ImpactAngle(max)'
  Str2DVarNames(nVarCount+10)='ImpactAngle(variance)'
  Str2DVarNames(nVarCount+11)='CurrentForceX'
  Str2DVarNames(nVarCount+12)='CurrentForceY'
  Str2DVarNames(nVarCount+13)='CurrentForceZ'
  Str2DVarNames(nVarCount+14)='AverageForceX'
  Str2DVarNames(nVarCount+15)='AverageForceY'
  Str2DVarNames(nVarCount+16)='AverageForceZ'

  IF (nSpecies.GT.1) THEN
      DO iSpec=1,nSpecies
          WRITE(SpecID,'(I3.3)') iSpec
          nShift = iSpec * (nErosionVars - 1)

          Str2DVarNames(nVarCount+nShift+1) ='Species'//TRIM(SpecID)//'_Impacts'
          Str2DVarNames(nVarCount+nShift+2) ='Species'//TRIM(SpecID)//'_ImpactsPerArea'
          Str2DVarNames(nVarCount+nShift+3) ='Species'//TRIM(SpecID)//'_E_kin(mean)'
          Str2DVarNames(nVarCount+nShift+4) ='Species'//TRIM(SpecID)//'_E_kin(min)'
          Str2DVarNames(nVarCount+nShift+5) ='Species'//TRIM(SpecID)//'_E_kin(max)'
          Str2DVarNames(nVarCount+nShift+6) ='Species'//TRIM(SpecID)//'_E_kin(variance)'
          Str2DVarNames(nVarCount+nShift+7) ='Species'//TRIM(SpecID)//'_ImpactAngle(mean)'
          Str2DVarNames(nVarCount+nShift+8) ='Species'//TRIM(SpecID)//'_ImpactAngle(min)'
          Str2DVarNames(nVarCount+nShift+9) ='Species'//TRIM(SpecID)//'_ImpactAngle(max)'
          Str2DVarNames(nVarCount+nShift+10)='Species'//TRIM(SpecID)//'_ImpactAngle(variance)'
          Str2DVarNames(nVarCount+nShift+11)='Species'//TRIM(SpecID)//'_CurrentForceX'
          Str2DVarNames(nVarCount+nShift+12)='Species'//TRIM(SpecID)//'_CurrentForceY'
          Str2DVarNames(nVarCount+nShift+13)='Species'//TRIM(SpecID)//'_CurrentForceZ'
          Str2DVarNames(nVarCount+nShift+14)='Species'//TRIM(SpecID)//'_AverageForceX'
          Str2DVarNames(nVarCount+nShift+15)='Species'//TRIM(SpecID)//'_AverageForceY'
          Str2DVarNames(nVarCount+nShift+16)='Species'//TRIM(SpecID)//'_AverageForceZ'
      END DO
  END IF

  CALL WriteAttribute(File_ID,'VarNamesSurface',nVar2D_Total,StrArray=Str2DVarNames)

  CALL CloseDataFile()
  DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF

CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

nVarCount=0
WRITE(H5_Name,'(A)') 'SurfaceData'
ASSOCIATE (&
      nGlobalSides         => nSurfTotalSides           , &
      nLocalSides          => nComputeNodeSurfSides     , &
      offsetSurfSide       => offsetComputeNodeSurfSide)

DO iSpec = 1,nSpecies
    CALL WriteArray(DataSetName = H5_Name                                                 , &
                    rank        = 4                                                       , &
                    nValGlobal  = (/nVar2D_Total,nSurfSample,nSurfSample,nGlobalSides  /) , &
                    nVal        = (/nVar2D_Spec ,nSurfSample,nSurfSample,nLocalSides   /) , &
                    offset      = (/nVarCount   ,          0,          0,offsetSurfSide/) , &
                    collective  = .TRUE.                                                  , &
                    RealArray   = MacroSurfaceSpecVal(:,:,:,:,iSpec))
    nVarCount = nVarCount + nVar2D_Spec
END DO

CALL WriteArray(    DataSetName = H5_Name                                                 , &
                    rank        = 4                                                       , &
                    nValGlobal  = (/NVar2D_Total,nSurfSample,nSurfSample,nGlobalSides  /) , &
                    nVal        = (/nVar2D      ,nSurfSample,nSurfSample,nLocalSides   /) , &
                    offset      = (/nVarCount   ,          0,          0,offsetSurfSide/) , &
                    collective  = .TRUE.                                                  , &
                    RealArray   = MacroSurfaceVal)
CALL CloseDataFile()

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)
IF (mySurfRank.EQ.0) THEN
#endif
  CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'RestartState'
  CALL WriteAttribute(File_ID,'NRestart',1,IntScalar=OutputVarCount)
  CALL CloseDataFile()
#if USE_MPI
END IF

CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*MPI*/

! This array is purely for restart
WRITE(H5_Name,'(A)') 'RestartData'
CALL WriteArray    (DataSetName = H5_Name                                                   , &
                    rank        = 4                                                         , &
                    nValGlobal  = (/OutputVarCount,nSurfSample,nSurfSample,nGlobalSides  /) , &
                    nVal        = (/OutputVarCount,nSurfSample,nSurfSample,nLocalSides   /) , &
                    offset      = (/0             ,          0,          0,offsetSurfSide/) , &
                    collective  = .TRUE.                                                    , &
                    RealArray   = OutputArray)
END ASSOCIATE
CALL CloseDataFile()

!> Only reset current forces here so we have them in case of restart
iSpec = 1
!-- Only one species. Only total values necessary
SampWallState_Shared(12:14,:,:,:) = 0.
!-- Multiple species. All Variables are saved DOUBLE. First Total, then per SPECIES
IF (nSpecies.GT.1) THEN
  DO iSpec=1,nSpecies
    nShift    = iSpec * (nErosionVars-1)
    nShiftRHS = iSpec * nErosionVars
    SampWallState_Shared(12+nShiftRHS:14+nShiftRHS,:,:,:) = 0.
  END DO
END IF

IF (mySurfRank.EQ.0) THEN
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES') 'DONE  [',EndT-StartT,'s]'
END IF

END SUBROUTINE WriteSurfSample


SUBROUTINE FinalizeParticleBoundarySampling()
!===================================================================================================================================
! deallocate everything
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Boundary_Vars
!USE MOD_Particle_Vars               ,ONLY: WriteMacroSurfaceValues
#if USE_MPI
USE MOD_Particle_MPI_Boundary_Sampling,ONLY: FinalizeSurfCommunication
USE MOD_Particle_MPI_Vars           ,ONLY: SurfSendBuf,SurfRecvBuf
USE MOD_Particle_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF,nSurfLeaders
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: iSurfSide
#if USE_MPI
INTEGER :: iProc,iError
#endif /*USE_MPI*/
!===================================================================================================================================

! Return if nothing was allocated
IF (.NOT.WriteMacroSurfaceValues) RETURN

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

CALL MPI_WIN_UNLOCK_ALL(SampWallState_Shared_Win      ,iError)
CALL MPI_WIN_FREE(      SampWallState_Shared_Win      ,iError)
CALL MPI_WIN_UNLOCK_ALL(GlobalSide2SurfSide_Shared_Win,iError)
CALL MPI_WIN_FREE(      GlobalSide2SurfSide_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(SurfSide2GlobalSide_Shared_Win,iError)
CALL MPI_WIN_FREE(      SurfSide2GlobalSide_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(SurfSideArea_Shared_Win       ,iError)
CALL MPI_WIN_FREE(      SurfSideArea_Shared_Win       ,iError)

CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! Deallocate surf leader mapping and free surf leader communicator
IF (ALLOCATED(SurfMapping)) THEN
  DO iProc = 0,nSurfLeaders-1
    SDEALLOCATE(SurfMapping(iProc)%RecvSurfGlobalID)
    SDEALLOCATE(SurfMapping(iProc)%SendSurfGlobalID)
  END DO
  SDEALLOCATE(SurfMapping)
END IF

SDEALLOCATE(SurfSendBuf)
SDEALLOCATE(SurfRecvBuf)

CALL FinalizeSurfCommunication()

IF(MPI_COMM_LEADERS_SURF.NE.MPI_COMM_NULL) THEN
  CALL MPI_COMM_FREE(MPI_COMM_LEADERS_SURF,iERROR)
END IF

MDEALLOCATE(SampWallState_Shared)
MDEALLOCATE(SurfSideArea_Shared)
MDEALLOCATE(GlobalSide2SurfSide_Shared)
MDEALLOCATE(SurfSide2GlobalSide_Shared)
#endif /*USE_MPI*/

! Then, free the pointers or arrays
SDEALLOCATE(XiEQ_SurfSample)
SDEALLOCATE(SurfBCName)
SDEALLOCATE(SurfSampleBCs)
SDEALLOCATE(SampWallState)
MDEALLOCATE(SurfSideArea)
MDEALLOCATE(GlobalSide2SurfSide)
MDEALLOCATE(SurfSide2GlobalSide)

END SUBROUTINE FinalizeParticleBoundarySampling

END MODULE MOD_Particle_Boundary_Sampling
