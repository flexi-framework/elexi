!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
! Contains routines to build the halo exchange
!===================================================================================================================================
MODULE MOD_Particle_MPI_Halo
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
PUBLIC:: IdentifyPartExchangeProcs
PUBLIC:: FinalizePartExchangeProcs
!===================================================================================================================================

CONTAINS

SUBROUTINE IdentifyPartExchangeProcs
!===================================================================================================================================
! Identifies processors in physical range for particle exchange communication. This communication has to occur at every RK step and
! would be too costly if done as an all-to-all communication
!> Procs on the same compute node are always assumed to be in range since communication is handled on proc
!> Procs in the compute node halo region are only considered in range if they lie within halo_eps of the mesh on the current proc
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_CalcTimeStep            ,ONLY: CalcTimeStep
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
USE MOD_Mesh_Vars               ,ONLY: ElemInfo_Shared,SideInfo_Shared!,NodeCoords_Shared
USE MOD_MPI_Vars                ,ONLY: offsetElemMPI
USE MOD_Particle_Globals        ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: BoundsOfElem_Shared
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID,GetGlobalNonUniqueSideID
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeTotalElems
USE MOD_MPI_Shared_Vars         ,ONLY: nComputenodeProcessors,nProcessors_Global,myLeaderGroupRank
USE MOD_Particle_MPI_Vars       ,ONLY: SafetyFactor,halo_eps,halo_eps_velo
USE MOD_Particle_MPI_Vars       ,ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc,CheckExchangeProcs
! USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Timedisc_Vars  ,ONLY: ManualTimeStep
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_ReadInTools             ,ONLY: PrintOption
USE MOD_TimeDisc_Vars           ,ONLY: nRKStages,RKc
USE MOD_Utils                   ,ONLY: ALMOSTZERO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Partner identification
LOGICAL                        :: ProcInRange
INTEGER                        :: iPeriodicVector,jPeriodicVector
INTEGER                        :: iDir,jDir,kDir
REAL,DIMENSION(6)              :: xCoordsProc,xCoordsOrigin
REAL                           :: origin(3),radius
INTEGER                        :: iElem,ElemID,firstElem,lastElem,NbElemID,HaloElem
INTEGER                        :: iSide,SideID,firstSide,lastSide,iLocSide
INTEGER                        :: iMortar,nMortarElems,NbSideID
INTEGER                        :: iProc,HaloProc
INTEGER                        :: nExchangeSides
INTEGER,ALLOCATABLE            :: ExchangeSides(:)
REAL                           :: BoundsOfElemCenter(1:4)
REAL,ALLOCATABLE               :: MPISideBoundsOfElemCenter(:,:)
INTEGER                        :: ExchangeProcLeader
LOGICAL,ALLOCATABLE            :: MPISideElem(:)
LOGICAL                        :: ProcHasExchangeElem
! halo_eps reconstruction
REAL                           :: MPI_halo_eps,MPI_halo_eps_velo,MPI_halo_diag,vec(1:3),deltaT
INTEGER                        :: iStage,errType
! Non-symmetric particle exchange
INTEGER,ALLOCATABLE            :: SendRequest(:),RecvRequest(:)
LOGICAL,ALLOCATABLE            :: GlobalProcToRecvProc(:)
LOGICAL,ALLOCATABLE            :: CommFlag(:)
INTEGER                        :: nNonSymmetricExchangeProcs,nNonSymmetricExchangeProcsGlob
INTEGER                        :: nExchangeProcessorsGlobal
! Timers
REAL                           :: StartT,EndT
!=================================================================================================================================

! Keep everything in sync here
! CALL MPI_BARRIER(MPI_COMM_FLEXI,IERROR)

IF (MPIRoot) THEN
  ! WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors ...'
  StartT = MPI_WTIME()
END IF ! MPIRoot

! Allocate arrays
ALLOCATE(GlobalProcToExchangeProc(EXCHANGE_PROC_SIZE,0:nProcessors_Global-1))
GlobalProcToExchangeProc(:,:) = -1

! Identify all procs on same node
nExchangeProcessors = 0

! This is generally not required, keep communication to a minimum
!DO iProc = ComputeNodeRootRank,ComputeNodeRootRank+nComputeNodeProcessors-1
!  ! Do not attempt to communicate with myself
!  IF (iProc.EQ.myRank) CYCLE
!
!  ! Build mapping global to compute-node
!  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
!  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
!  nExchangeProcessors = nExchangeProcessors + 1
!END DO

! Communicate non-symmetric part exchange partners to catch non-symmetric proc identification due to inverse distance calculation
ALLOCATE(GlobalProcToRecvProc(0:nProcessors_Global-1), &
         SendRequest         (0:nProcessors_Global-1), &
         RecvRequest         (0:nProcessors_Global-1), &
         CommFlag            (0:nProcessors_Global-1))

GlobalProcToRecvProc = .FALSE.

! Identify all procs with elements in range. This includes checking the procs on the compute-node as they might lie far apart
IF (nProcessors.EQ.1) THEN
  SWRITE(UNIT_stdOut,'(A)') ' | Running on one processor. Particle exchange communication disabled.'
  SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors DONE!'
  SWRITE(UNIT_stdOut,'(132("-"))')
  RETURN
END IF

! Open receive buffer for non-symmetric exchange identification
IF (CheckExchangeProcs) THEN
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    CALL MPI_IRECV( GlobalProcToRecvProc(iProc)  &
                  , 1                            &
                  , MPI_LOGICAL                  &
                  , iProc                        &
                  , 1999                         &
                  , MPI_COMM_FLEXI               &
                  , RecvRequest(iProc)           &
                  , IERROR)
  END DO
END IF

!> Count all MPI sides on current proc.
firstElem = offsetElem+1
lastElem  = offsetElem+nElems

! This approach does not work, we get only the MPI sides pointing into the compute-node halo region
!DO iSide = firstSide,lastSide
!  IF (SideInfo_Shared(SIDE_NBELEMTYPE,iSide).EQ.2) THEN
!    nExchangeSides = nExchangeSides + 1
!  END IF
!END DO

!>>> For all elements, loop over the six sides and check if the neighbor element is on the current proc
!>>> Special care for big mortar sides, here the SIDE_ELEMID must be used
ALLOCATE(MPISideElem(offsetElem+1:offsetElem+nElems))
nExchangeSides = 0
MPISideElem    = .FALSE.

DO iElem = firstElem,lastElem
  ! Element already flagged
  IF (MPISideElem(iElem)) CYCLE

  DO iLocSide = 1,6
    SideID   = GetGlobalNonUniqueSideID(iElem,iLocSide)
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

    ! Mortar side
    IF (NbElemID.LT.0) THEN
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

      DO iMortar = 1,nMortarElems
        NbSideID = -SideInfo_Shared(SIDE_LOCALID,SideID + iMortar)
        ! If small mortar side not defined, skip it for now, likely not inside the halo region
        IF (NbSideID.LT.1) CYCLE

        NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID + iMortar)
        ! If small mortar element not defined, skip it for now, likely not inside the halo region
        IF (NbElemID.LT.1) CYCLE

        ! If any of the small mortar sides is not on the local proc, the side is a MPI side
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides     = nExchangeSides + 1
          MPISideElem(iElem) = .TRUE.
          EXIT
        END IF
      END DO

    ! regular side or small mortar side
    ELSE
      ! Only check inner (MPI interfaces) and boundary sides (NbElemID.EQ.0)
      ! Boundary sides cannot be discarded because one proc might have MPI interfaces and the other might not
      ! NbElemID.LT.firstElem is always true for NbElemID.EQ.0 because firstElem.GE.1
      IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
        nExchangeSides     = nExchangeSides + 1
        MPISideElem(iElem) = .TRUE.
      END IF
    END IF
  END DO
END DO

IF (nComputeNodeProcessors.GT.1.AND.nExchangeSides.EQ.0) &
  CALL Abort(__STAMP__,'Found no side connectivity between processor domains')

!> Build mapping for all MPI sides on current proc
ALLOCATE(ExchangeSides(1:nExchangeSides))

nExchangeSides = 0
MPISideElem    = .FALSE.

DO iElem = firstElem,lastElem
  ! Element already flagged
  IF (MPISideElem(iElem)) CYCLE

  DO iLocSide = 1,6
    SideID   = GetGlobalNonUniqueSideID(iElem,iLocSide)
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

    ! Mortar side
    IF (NbElemID.LT.0) THEN
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

      DO iMortar = 1,nMortarElems
        NbSideID = -SideInfo_Shared(SIDE_LOCALID,SideID + iMortar)
        ! If small mortar side not defined, skip it for now, likely not inside the halo region
        IF (NbSideID.LT.1) CYCLE

        NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID + iMortar)
        ! If small mortar element not defined, skip it for now, likely not inside the halo region
        IF (NbElemID.LT.1) CYCLE

        ! If any of the small mortar sides is not on the local proc, the side is a MPI side
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides = nExchangeSides + 1
          ExchangeSides(nExchangeSides) = SideID
          MPISideElem(iElem) = .TRUE.
          EXIT
        END IF
      END DO

    ! regular side or small mortar side
    ELSE
      ! Only check inner (MPI interfaces) and boundary sides (NbElemID.EQ.0)
      ! Boundary sides cannot be discarded because one proc might have MPI interfaces and the other might not
      ! NbElemID.LT.firstElem is always true for NbElemID.EQ.0 because firstElem.GE.1
      IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
        nExchangeSides = nExchangeSides + 1
        ExchangeSides(nExchangeSides) = SideID
        MPISideElem(iElem) = .TRUE.
      END IF
    END IF
  END DO
END DO

DEALLOCATE(MPISideElem)

!> Build metrics for all MPI sides on current proc
ALLOCATE(MPISideBoundsOfElemCenter(1:4,1:nExchangeSides))

DO iSide = 1, nExchangeSides
  SideID = ExchangeSides(iSide)
  ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
  MPISideBoundsOfElemCenter(1:3,iSide) = (/   SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
                                              SUM(BoundsOfElem_Shared(1:2,2,ElemID)), &
                                              SUM(BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
  MPISideBoundsOfElemCenter(4,iSide) = VECNORM ((/BoundsOfElem_Shared(2  ,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                                  BoundsOfElem_Shared(2  ,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                                  BoundsOfElem_Shared(2  ,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
END DO

!> Check all elements in the CN halo region against local MPI sides. Check is identical to particle_bgm.f90
!>>> Check the bounding box of each element in compute-nodes' halo domain against the bounding boxes of the elements of the
!>>> MPI-surface (local proc MPI sides)

! if running on one node, halo_eps is meaningless. Get a representative MPI_halo_eps for MPI proc identification
IF (halo_eps.LE.0.) THEN
  ! reconstruct halo_eps_velo
  IF (halo_eps_velo.EQ.0) THEN
    ! Set to speed of sound
    MPI_halo_eps_velo = 343
  ELSE
    MPI_halo_eps_velo = halo_eps_velo
  END IF

  ! reconstruct deltaT
  deltaT = 0.
  IF (ManualTimeStep.GT.0.) THEN
    deltaT    = ManualTimeStep
  ELSE
    deltaT    = CalcTimeStep(errType)
  END IF

  ! calculate halo_eps
  MPI_halo_eps = RKc(2)
  DO iStage=2,nRKStages-1
    MPI_halo_eps = MAX(MPI_halo_eps,RKc(iStage+1)-RKc(iStage))
  END DO
  MPI_halo_eps = MAX(MPI_halo_eps,1.-RKc(nRKStages))
  MPI_halo_eps = MPI_halo_eps*MPI_halo_eps_velo*deltaT*SafetyFactor

  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  MPI_halo_diag = VECNORM(vec)

  ! compare halo_eps against global diagonal and reduce if necessary
  IF (.NOT.ALMOSTZERO(MPI_halo_eps).AND.(MPI_halo_diag.GE.MPI_halo_eps)) THEN
    CALL PrintOption('No halo_eps given. Reconstructed','CALC',RealOpt=MPI_halo_eps)
  ELSEIF (.NOT.ALMOSTZERO(MPI_halo_eps).AND.(MPI_halo_diag.LT.MPI_halo_eps)) THEN
    MPI_halo_eps = MPI_halo_diag
    CALL PrintOption('No halo_eps given. Reconstructed to global diag','CALC',RealOpt=MPI_halo_eps)
  ! halo_eps still at zero. Set it to global diagonal
  ELSE
    MPI_halo_eps = MPI_halo_diag
    CALL PrintOption('No halo_eps given and could not be reconstructed. Using global diag','CALC',RealOpt=MPI_halo_eps)
  END IF
ELSE
  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  MPI_halo_diag = VECNORM(vec)

  MPI_halo_eps = halo_eps
END IF

! Identify all procs with elements in range
! xCoordsProc(1) = GEO%xmin
! xCoordsProc(2) = GEO%xmax
! xCoordsProc(3) = GEO%ymin
! xCoordsProc(4) = GEO%ymax
! xCoordsProc(5) = GEO%zmin
! xCoordsProc(6) = GEO%zmax
xCoordsProc(1) = -HUGE(1.)
xCoordsProc(2) =  HUGE(1.)
xCoordsProc(3) = -HUGE(1.)
xCoordsProc(4) =  HUGE(1.)
xCoordsProc(5) = -HUGE(1.)
xCoordsProc(6) =  HUGE(1.)

DO iElem = 1, nElems
  ElemID    = iElem + offsetElem
  ! Flag elements depending on radius
  origin(1:3) = (/ SUM(   BoundsOfElem_Shared(1:2,1,ElemID)), &
                   SUM(   BoundsOfElem_Shared(1:2,2,ElemID)), &
                   SUM(   BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
  ! Calculate halo element outer radius
  radius    = VECNORM ((/ BoundsOfElem_Shared(2  ,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                          BoundsOfElem_Shared(2  ,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                          BoundsOfElem_Shared(2  ,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)

  xCoordsProc(1) = MIN(xCoordsProc(1), origin(1) - radius)
  xCoordsProc(2) = MAX(xCoordsProc(2), origin(1) + radius)
  xCoordsProc(3) = MIN(xCoordsProc(3), origin(2) - radius)
  xCoordsProc(4) = MAX(xCoordsProc(4), origin(2) + radius)
  xCoordsProc(5) = MIN(xCoordsProc(5), origin(3) - radius)
  xCoordsProc(6) = MAX(xCoordsProc(6), origin(3) + radius)
END DO ! iElem = 1, nElems

! Use a named loop so the entire element can be cycled
ElemLoop:  DO iElem = 1,nComputeNodeTotalElems
  ElemID   = GetGlobalElemID(iElem)
  HaloProc = ElemInfo_Shared(ELEM_RANK,ElemID)

  ! Skip elements on same rank
  IF (HaloProc.EQ.myRank) CYCLE ElemLoop

  ! Skip tracing mortar elements
  IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).EQ.4) CYCLE ElemLoop

  BoundsOfElemCenter(1:3) = (/SUM(      BoundsOfElem_Shared(1:2,1,ElemID)), &
                              SUM(      BoundsOfElem_Shared(1:2,2,ElemID)), &
                              SUM(      BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
  BoundsOfElemCenter(4)   = VECNORM ((/ BoundsOfElem_Shared(2  ,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                        BoundsOfElem_Shared(2  ,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                        BoundsOfElem_Shared(2  ,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)

!#if CODE_ANALYZE
!    ! Sanity checks. Elems in halo region must have ELEM_HALOFLAG=2 and the proc must not be flagged yet
!    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.2) THEN
!      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Flag: ',ElemInfo_Shared(ELEM_HALOFLAG,ElemID)
!      CALL Abort(__STAMP__,  'Element found in range of halo elements while not flagged as such!')
!    END IF
!
!    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.1) THEN
!      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Proc: ',HaloProc
!      CALL Abort(__STAMP__, 'Proc claimed to have elements both on compute node and in halo region!')
!    END IF
!#endif

  SELECT CASE(GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc))
    ! Proc not previously encountered, check if possibly in range
    CASE(-1)
      firstElem = offsetElemMPI(HaloProc)+1
      lastElem  = offsetElemMPI(HaloProc +1)

      SELECT CASE(TrackingMethod)
        ! TriaTracking does not have curved elements, nodeCoords are sufficient
        CASE(TRIATRACKING)
          ! xCoordsOrigin(1) = MINVAL(NodeCoords_Shared(1,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
          !                                              :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
          ! xCoordsOrigin(2) = MAXVAL(NodeCoords_Shared(1,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
          !                                              :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
          ! xCoordsOrigin(3) = MINVAL(NodeCoords_Shared(2,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
          !                                              :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
          ! xCoordsOrigin(4) = MAXVAL(NodeCoords_Shared(2,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
          !                                              :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
          ! xCoordsOrigin(5) = MINVAL(NodeCoords_Shared(3,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
          !                                              :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
          ! xCoordsOrigin(6) = MAXVAL(NodeCoords_Shared(3,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
          !                                              :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))

          xCoordsOrigin(1) = -HUGE(1.)
          xCoordsOrigin(2) =  HUGE(1.)
          xCoordsOrigin(3) = -HUGE(1.)
          xCoordsOrigin(4) =  HUGE(1.)
          xCoordsOrigin(5) = -HUGE(1.)
          xCoordsOrigin(6) =  HUGE(1.)

          DO HaloElem = firstElem, lastElem
            ! Flag elements depending on radius
            origin(1:3) = (/ SUM(   BoundsOfElem_Shared(1:2,1,HaloElem)), &
                             SUM(   BoundsOfElem_Shared(1:2,2,HaloElem)), &
                             SUM(   BoundsOfElem_Shared(1:2,3,HaloElem)) /) / 2.
            ! Calculate halo element outer radius
            radius    = VECNORM ((/ BoundsOfElem_Shared(2  ,1,HaloElem)-BoundsOfElem_Shared(1,1,HaloElem), &
                                    BoundsOfElem_Shared(2  ,2,HaloElem)-BoundsOfElem_Shared(1,2,HaloElem), &
                                    BoundsOfElem_Shared(2  ,3,HaloElem)-BoundsOfElem_Shared(1,3,HaloElem) /) / 2.)

            xCoordsOrigin(1) = MIN(xCoordsOrigin(1), origin(1) - radius)
            xCoordsOrigin(2) = MAX(xCoordsOrigin(2), origin(1) + radius)
            xCoordsOrigin(3) = MIN(xCoordsOrigin(3), origin(2) - radius)
            xCoordsOrigin(4) = MAX(xCoordsOrigin(4), origin(2) + radius)
            xCoordsOrigin(5) = MIN(xCoordsOrigin(5), origin(3) - radius)
            xCoordsOrigin(6) = MAX(xCoordsOrigin(6), origin(3) + radius)
        END DO ! HaloElem = firstElem, lastElem

      ! Build mesh min/max on BezierControlPoints for possibly curved elements
      CASE(REFMAPPING,TRACING)
        firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,firstElem)+1
        lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND ,lastElem)

        ! xCoordsOrigin(1) = MINVAL(BezierControlPoints3D(1,:,:,firstSide:lastSide))
        ! xCoordsOrigin(2) = MAXVAL(BezierControlPoints3D(1,:,:,firstSide:lastSide))
        ! xCoordsOrigin(3) = MINVAL(BezierControlPoints3D(2,:,:,firstSide:lastSide))
        ! xCoordsOrigin(4) = MAXVAL(BezierControlPoints3D(2,:,:,firstSide:lastSide))
        ! xCoordsOrigin(5) = MINVAL(BezierControlPoints3D(3,:,:,firstSide:lastSide))
        ! xCoordsOrigin(6) = MAXVAL(BezierControlPoints3D(3,:,:,firstSide:lastSide))

        xCoordsOrigin(1) = -HUGE(1.)
        xCoordsOrigin(2) =  HUGE(1.)
        xCoordsOrigin(3) = -HUGE(1.)
        xCoordsOrigin(4) =  HUGE(1.)
        xCoordsOrigin(5) = -HUGE(1.)
        xCoordsOrigin(6) =  HUGE(1.)

        DO HaloElem = firstElem, lastElem
          ! Flag elements depending on radius
          origin(1:3) = (/ SUM(   BoundsOfElem_Shared(1:2,1,HaloElem)), &
                           SUM(   BoundsOfElem_Shared(1:2,2,HaloElem)), &
                           SUM(   BoundsOfElem_Shared(1:2,3,HaloElem)) /) / 2.
          ! Calculate halo element outer radius
          radius    = VECNORM ((/ BoundsOfElem_Shared(2  ,1,HaloElem)-BoundsOfElem_Shared(1,1,HaloElem), &
                                  BoundsOfElem_Shared(2  ,2,HaloElem)-BoundsOfElem_Shared(1,2,HaloElem), &
                                  BoundsOfElem_Shared(2  ,3,HaloElem)-BoundsOfElem_Shared(1,3,HaloElem) /) / 2.)

          xCoordsOrigin(1) = MIN(xCoordsOrigin(1), origin(1) - radius)
          xCoordsOrigin(2) = MAX(xCoordsOrigin(2), origin(1) + radius)
          xCoordsOrigin(3) = MIN(xCoordsOrigin(3), origin(2) - radius)
          xCoordsOrigin(4) = MAX(xCoordsOrigin(4), origin(2) + radius)
          xCoordsOrigin(5) = MIN(xCoordsOrigin(5), origin(3) - radius)
          xCoordsOrigin(6) = MAX(xCoordsOrigin(6), origin(3) + radius)
        END DO ! HaloElem = firstElem, lastElem
      END SELECT

      ! Keep direction to account for accuracy issues
      IF (myRank.LT.HaloProc) THEN
        ProcInRange = HaloBoxInProc(xCoordsOrigin,xCoordsProc,MPI_halo_eps,GEO%nPeriodicVectors,GEO%PeriodicVectors)
      ELSE
        ProcInRange = HaloBoxInProc(xCoordsProc,xCoordsOrigin,MPI_halo_eps,GEO%nPeriodicVectors,GEO%PeriodicVectors)
      END IF

      ! Check if proc is in range
      IF (.NOT.ProcInRange) THEN
        ! Proc definitely not in range
        GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = -2
        CYCLE ElemLoop
      ELSE
        ! Proc possible in range
        GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 0
      END IF

    ! Proc definitely not in range or already flagged
    CASE(-2,1,2)
      CYCLE ElemLoop
  END SELECT

  DO iSide = 1, nExchangeSides
    ! compare distance of centers with sum of element outer radii+halo_eps
    IF (VECNORM(BoundsOfElemCenter(1:3)-MPISideBoundsOfElemCenter(1:3,iSide)) &
      .GT. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide)) THEN

      ! Also check periodic directions. Only MPI sides of the local proc are
      ! taken into account, so do not perform additional case distinction
      SELECT CASE(GEO%nPeriodicVectors)
        ! One periodic vector
        CASE(1)
          DO iDir = -1,1,2
            IF (VECNORM( BoundsOfElemCenter(1:3)                                                       &
                       + GEO%PeriodicVectors(1:3,1) * REAL(iDir)                                       &
                       - MPISideBoundsOfElemCenter(1:3,iSide))                                         &
                    .GT. MPI_halo_eps+BoundsOfElemCenter(4)                                            & !-BoundsOfElemCenter(5) &
                       + MPISideBoundsOfElemCenter(4,iSide) ) CYCLE

              ! flag the proc as exchange proc (in halo region)
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ElemLoop
            END DO

        ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
        CASE(2)
          DO iPeriodicVector = 1,2
            DO iDir = -1,1,2
              ! check if element is within halo_eps of periodically displaced element
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                         + GEO%PeriodicVectors(1:3,iPeriodicVector) * REAL(iDir)                      &
                         - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                      .GT. MPI_halo_eps+BoundsOfElemCenter(4)                                         & !-BoundsOfElemCenter(5) &
                         + MPISideBoundsOfElemCenter(4,iSide)) CYCLE

              ! flag the proc as exchange proc (in halo region)
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ElemLoop
            END DO
          END DO

          DO iDir = -1,1,2
            DO jDir = -1,1,2
              ! check if element is within halo_eps of periodically displaced element
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                         + GEO%PeriodicVectors(1:3,1) * REAL(iDir)                                    &
                         + GEO%PeriodicVectors(1:3,2) * REAL(jDir)                                    &
                         - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                        .GT. MPI_halo_eps+BoundsOfElemCenter(4)                                       & !-BoundsOfElemCenter(5) &
                           + MPISideBoundsOfElemCenter(4,iSide)) CYCLE

              ! flag the proc as exchange proc (in halo region)
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ElemLoop
            END DO
          END DO

        ! Three periodic vectors. Also check linear combination, see particle_bgm.f90
        CASE(3)
          ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
          ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
          DO iPeriodicVector = 1,3
            DO iDir = -1,1,2
              ! element might be already added back
              ! check if element is within halo_eps of periodically displaced element
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                      &
                         + GEO%PeriodicVectors(1:3,iPeriodicVector) * REAL(iDir)                        &
                         - MPISideBoundsOfElemCenter(1:3,iSide))                                        &
                      .GT. MPI_halo_eps+BoundsOfElemCenter(4)                                           & !-BoundsOfElemCenter(5) &
                         + MPISideBoundsOfElemCenter(4,iSide)) CYCLE

              ! flag the proc as exchange proc (in halo region)
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ElemLoop
            END DO
          END DO

          DO iPeriodicVector = 1,3
            DO iDir = -1,1,2
              DO jPeriodicVector = 1,3
                DO jDir = -1,1,2
                  IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

                  ! check if element is within halo_eps of periodically displaced element
                  IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                             + GEO%PeriodicVectors(1:3,iPeriodicVector) * REAL(iDir)                        &
                             + GEO%PeriodicVectors(1:3,jPeriodicVector) * REAL(jDir)                        &
                             - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                          .GT. MPI_halo_eps+BoundsOfElemCenter(4)                                           & !-BoundsOfElemCenter(5) &
                             + MPISideBoundsOfElemCenter(4,iSide)) CYCLE

                  ! flag the proc as exchange proc (in halo region)
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  CYCLE ElemLoop
                END DO
              END DO
            END DO
          END DO

          ! check if element is within halo_eps of periodically displaced element
          DO iDir = -1,1,2
            DO jDir = -1,1,2
              DO kDir = -1,1,2
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                        &
                           + GEO%PeriodicVectors(1:3,1) * REAL(iDir)                                        &
                           + GEO%PeriodicVectors(1:3,2) * REAL(jDir)                                        &
                           + GEO%PeriodicVectors(1:3,3) * REAL(kDir)                                        &
                           - MPISideBoundsOfElemCenter(1:3,iSide))                                          &
                        .GT. MPI_halo_eps+BoundsOfElemCenter(4)                                             & !-BoundsOfElemCenter(5) &
                           + MPISideBoundsOfElemCenter(4,iSide) ) CYCLE

                ! flag the proc as exchange proc (in halo region)
                GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                nExchangeProcessors = nExchangeProcessors + 1
                CYCLE ElemLoop
              END DO
            END DO
          END DO

        ! No periodic vectors, element out of range
        CASE(0)
          ! Do nothing

        CASE DEFAULT
          CALL Abort(__STAMP__,'Invalid number of periodic vectors in particle_mpi_halo.f90')

      END SELECT

    ! Element is in range of not-periodically displaced MPI side
    ELSE
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
      GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
      nExchangeProcessors = nExchangeProcessors + 1
      CYCLE ElemLoop
    END IF
  END DO ! iSide = 1, nExchangeSides
END DO ElemLoop

! Remove elements if the halo proc contains only internal elements, i.e. we cannot possibly reach the halo element
!
!   CN1     CN2    > If a processor contains large changes in element size, internal elements might intersect with
!  _ _ _    _ _ _  > the MPI sides. Since a processor checks all potential halo elements against only its own exchange
! |_|_|_|  |_|_|_| > sides, the large elements will only take effect for the proc not containing it. However, if the
! |_|_|_|  |_| | | > proc flags only large internal elements without flagging a single exchange element, there is no
! |_|_|_|  |_|_|_| > way for a particle to actually reach.
! |_|_|_|  |_| | |
! |_|_|_|  |_|_|_| > This routine therefore checks for the presence of exchange sides on the procs and unflags the
!                  > proc if none is found.
!
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).LE.0) CYCLE

  ProcHasExchangeElem = .FALSE.
  ! Use a named loop so the entire element can be cycled
ExchangeLoop: DO iElem = offsetElemMPI(iProc)+1,offsetElemMPI(iProc+1)
    ! Ignore elements outside nComputeNodeTotalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).LT.0) CYCLE ExchangeLoop

    DO iLocSide = 1,6
      SideID   = GetGlobalNonUniqueSideID(iElem,iLocSide)
      NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

      ! Mortar side
      IF (NbElemID.LT.0) THEN
        nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

        DO iMortar = 1,nMortarElems
          NbSideID = -SideInfo_Shared(SIDE_LOCALID,SideID + iMortar)
          ! If small mortar side not defined, skip it for now, likely not inside the halo region
          IF (NbSideID.LT.1) CYCLE

          NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID + iMortar)
          ! If small mortar element not defined, skip it for now, likely not inside the halo region
          IF (NbElemID.LT.1) CYCLE

          ! If any of the small mortar sides is not on the halo proc, the side is a MPI side
          IF (NbElemID.LT.offsetElemMPI(iProc)+1 .OR. NbElemID.GT.offsetElemMPI(iProc+1)) THEN
            ProcHasExchangeElem = .TRUE.
            EXIT ExchangeLoop
          END IF
        END DO

      ! regular side or small mortar side
      ELSE
        ! Only check inner (MPI interfaces) and boundary sides (NbElemID.EQ.0)
        ! Boundary sides cannot be discarded because one proc might have MPI interfaces and the other might not
        IF (NbElemID.LT.offsetElemMPI(iProc)+1 .OR. NbElemID.GT.offsetElemMPI(iProc+1)) THEN
          ProcHasExchangeElem = .TRUE.
          EXIT ExchangeLoop
        END IF
      END IF
    END DO
  END DO ExchangeLoop ! iElem = offsetElemMPI(iProc)+1,offsetElemMPI(iProc+1)

  ! Processor has halo elements but no MPI sides, remove the exchange processor
  IF (.NOT.ProcHasExchangeElem) THEN
    GlobalProcToExchangeProc(:,iProc) = 0
    nExchangeProcessors = nExchangeProcessors - 1
  END IF
END DO ! iProc = 1,nExchangeProcessors

! Notify every proc if it was identified by the local proc
IF(CheckExchangeProcs)THEN
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    ! CommFlag holds the information if the local proc wants to communicate with iProc. Cannot be a logical because ISEND might not
    ! return before the next value is written
    CommFlag(iProc) = MERGE(.TRUE.,.FALSE.,GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).GT.0)
    CALL MPI_ISEND( CommFlag(iProc)              &
                  , 1                            &
                  , MPI_LOGICAL                  &
                  , iProc                        &
                  , 1999                         &
                  , MPI_COMM_FLEXI               &
                  , SendRequest(iProc)           &
                  , IERROR)
  END DO

  ! Finish communication
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)
    CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL Abort(__STAMP__,' MPI Communication error', IERROR)
  END DO

  ! Append previously not found procs to list of exchange processors
  nNonSymmetricExchangeProcs = 0
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    ! Ignore procs that are already flagged or not requesting communication
    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).GT.0) CYCLE
    ! GlobalProcToRecvProc(iProc) is true if iProc has flagged myRank
    IF (.NOT.GlobalProcToRecvProc(iProc)) CYCLE

    ! Found a previously missing proc
    GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
    GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors

    nNonSymmetricExchangeProcs = nNonSymmetricExchangeProcs + 1
    nExchangeProcessors        = nExchangeProcessors + 1
  END DO

  ! On smooth grids, nNonSymmetricExchangeProcs should be zero. Only output if previously missing particle exchange procs are found
  CALL MPI_ALLREDUCE(nNonSymmetricExchangeProcs,nNonSymmetricExchangeProcsGlob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
  IF (nNonSymmetricExchangeProcsGlob.GT.0) THEN
    SWRITE(UNIT_stdOut,'(132("~"))')
    SWRITE(UNIT_stdOut,'(A,I0,A)') ' | Found ',nNonSymmetricExchangeProcsGlob, &
                                   ' previously missing non-symmetric particle exchange procs'
    IF(CheckExchangeProcs) CALL CollectiveStop(__STAMP__,&
      ' Non-symmetric particle exchange procs > 0. This check is optional. You can disable it via CheckExchangeProcs = F')
  END IF ! nNonSymmetricExchangeProcsGlob.GT.0
END IF

SDEALLOCATE(GlobalProcToRecvProc)
SDEALLOCATE(RecvRequest)
SDEALLOCATE(SendRequest)
SDEALLOCATE(CommFlag)

! Build reverse mapping
!-- EXCHANGE_PROC_TYPE information is currently unused and either -1 (no communication) or 2 (communication). Can be used to
!-- implement check if exchange partner is on the same compute node, so build it here
ALLOCATE(ExchangeProcToGlobalProc(2,0:nExchangeProcessors-1))

DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).GT.0) THEN
    ! Find it the other proc is on the same compute node
    ExchangeProcLeader = INT(GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)/nComputeNodeProcessors)
    IF (ExchangeProcLeader.EQ.myLeaderGroupRank) THEN
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = iProc
    ELSE
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = 2
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = iProc
    END IF
  END IF
END DO

! -- Average number of exchange processors
CALL MPI_REDUCE(nExchangeProcessors,nExchangeProcessorsGlobal,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
SWRITE(UNIT_stdOut,'(A,I0,A)') ' | Started particle exchange communication with average ', &
                                 nExchangeProcessorsGlobal/nProcessors_Global            , &
                                 ' partners per proc'

IF (MPIRoot) THEN
  EndT = MPI_WTIME()
  CALL DisplayMessageAndTime(EndT-StartT,'IDENTIFYING Particle Exchange Processors DONE!',DisplayDespiteLB=.TRUE.,DisplayLine=.TRUE.)
END IF ! MPIRoot

END SUBROUTINE IdentifyPartExchangeProcs


!===================================================================================================================================
! Deallocates arrays for halo exchange
!===================================================================================================================================
SUBROUTINE FinalizePartExchangeProcs()
! MODULES
USE MOD_Particle_MPI_Vars       ,ONLY: ExchangeProcToGlobalProc,GlobalProcToExchangeProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(ExchangeProcToGlobalProc)
SDEALLOCATE(GlobalProcToExchangeProc)

END SUBROUTINE FinalizePartExchangeProcs


PURE FUNCTION HaloBoxInProc(CartNodes,CartProc,halo_eps,nPeriodicVectors,PeriodicVectors)
!===================================================================================================================================
! Check if bounding box is on proc by comparing against the other bounding box extended by halo_eps
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: CartNodes(6)
REAL,INTENT(IN)             :: CartProc( 6)
REAL,INTENT(IN)             :: halo_eps
INTEGER,INTENT(IN)          :: nPeriodicVectors
REAL,INTENT(IN),OPTIONAL    :: PeriodicVectors(3,nPeriodicVectors)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                     :: HaloBoxInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPeriodicVector,jPeriodicVector
INTEGER                        :: iDir,jDir,kDir
REAL,DIMENSION(1:6)            :: xCordsPeri
!===================================================================================================================================

HaloBoxInProc = .FALSE.

! Check whether the bounding boxes intersect
IF (   ((CartNodes(1).LE.CartProc(2)+halo_eps).AND.(CartNodes(2).GE.CartProc(1)-halo_eps))  &
  .AND.((CartNodes(3).LE.CartProc(4)+halo_eps).AND.(CartNodes(4).GE.CartProc(3)-halo_eps))  &
  .AND.((CartNodes(5).LE.CartProc(6)+halo_eps).AND.(CartNodes(6).GE.CartProc(5)-halo_eps))) HaloBoxInProc = .TRUE.

! Also check periodic directions. Only MPI sides of the local proc are taken into account, so do not perform additional case distinction
SELECT CASE(nPeriodicVectors)
  ! One periodic vector
  CASE(1)
    DO iDir = -1,1,2
      ! check if element is within halo_eps of periodically displaced element
      xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,1) * REAL(iDir)
      xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,1) * REAL(iDir)
      xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,1) * REAL(iDir)

      ! Check whether the bounding boxes intersect
      IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
        .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
        .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
        HaloBoxInProc = .TRUE.
        RETURN
      END IF
    END DO

  ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
  CASE(2)
    DO iPeriodicVector = 1,2
      DO iDir = -1,1,2
        ! check if element is within halo_eps of periodically displaced element
        xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,iPeriodicVector) * REAL(iDir)
        xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,iPeriodicVector) * REAL(iDir)
        xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,iPeriodicVector) * REAL(iDir)

        ! Check whether the bounding boxes intersect
        IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
          .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
          .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
          HaloBoxInProc = .TRUE.
          RETURN
        END IF
      END DO
    END DO

    DO iDir = -1,1,2
      DO jDir = -1,1,2
        ! check if element is within halo_eps of periodically displaced element
        xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,1) * REAL(iDir) &
                                         + PeriodicVectors(1,2) * REAL(jDir)
        xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,1) * REAL(iDir) &
                                         + PeriodicVectors(2,2) * REAL(jDir)
        xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,1) * REAL(iDir) &
                                         + PeriodicVectors(3,2) * REAL(jDir)

        ! Check whether the bounding boxes intersect
        IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
          .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
          .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
          HaloBoxInProc = .TRUE.
          RETURN
        END IF
      END DO
    END DO

  ! Three periodic vectors. Also check linear combination, see particle_bgm.f90
  CASE(3)
    ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
    ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
    DO iPeriodicVector = 1,3
      DO iDir = -1,1,2
        ! check if element is within halo_eps of periodically displaced element
        xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,iPeriodicVector) * REAL(iDir)
        xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,iPeriodicVector) * REAL(iDir)
        xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,iPeriodicVector) * REAL(iDir)

        ! Check whether the bounding boxes intersect
        IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
          .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
          .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
          HaloBoxInProc = .TRUE.
          RETURN
        END IF

        DO jPeriodicVector = 1,3
          DO jDir = -1,1,2
            IF (iPeriodicVector.GE.jPeriodicVector) CYCLE
            ! check if element is within halo_eps of periodically displaced element
            xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,iPeriodicVector) * REAL(iDir) &
                                             + PeriodicVectors(1,jPeriodicVector) * REAL(jDir)
            xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,iPeriodicVector) * REAL(iDir) &
                                             + PeriodicVectors(2,jPeriodicVector) * REAL(jDir)
            xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,iPeriodicVector) * REAL(iDir) &
                                             + PeriodicVectors(3,jPeriodicVector) * REAL(jDir)

            ! Check whether the bounding boxes intersect
            IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
              .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
              .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
              HaloBoxInProc = .TRUE.
              RETURN
            END IF
          END DO
        END DO
      END DO
    END DO

    ! check if element is within halo_eps of periodically displaced element
    DO iDir = -1,1,2
      DO jDir = -1,1,2
        DO kDir = -1,1,2
          ! check if element is within halo_eps of periodically displaced element
          xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,1) * REAL(iDir) &
                                           + PeriodicVectors(1,2) * REAL(jDir) &
                                           + PeriodicVectors(1,3) * REAL(kDir)
          xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,1) * REAL(iDir) &
                                           + PeriodicVectors(2,2) * REAL(jDir) &
                                           + PeriodicVectors(2,3) * REAL(kDir)
          xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,1) * REAL(iDir) &
                                           + PeriodicVectors(3,2) * REAL(jDir) &
                                           + PeriodicVectors(3,3) * REAL(kDir)

          ! Check whether the bounding boxes intersect
          IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
            .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
            .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
            HaloBoxInProc = .TRUE.
            RETURN
          END IF
        END DO
      END DO
    END DO

END SELECT

END FUNCTION HaloBoxInProc
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI_Halo
