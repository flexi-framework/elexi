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
#include "particle.h"

!===================================================================================================================================
! Contains routines to build the halo exchange
!===================================================================================================================================
MODULE MOD_Particle_MPI_Halo
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE IdentifyPartExchangeProcs
  MODULE PROCEDURE IdentifyPartExchangeProcs
END INTERFACE

INTERFACE FinalizePartExchangeProcs
  MODULE PROCEDURE FinalizePartExchangeProcs
END INTERFACE

PUBLIC :: IdentifyPartExchangeProcs
PUBLIC :: FinalizePartExchangeProcs
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
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
!USE MOD_MPI_Vars                ,ONLY: offsetElemMPI
USE MOD_Particle_Globals        ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
!USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeElems
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,SideInfo_Shared!,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: BoundsOfElem_Shared
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID,GetGlobalNonUniqueSideID
USE MOD_Particle_MPI_Vars       ,ONLY: halo_eps
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputenodeProcessors,nProcessors_Global,myLeaderGroupRank!,ComputeNodeRootRank
USE MOD_Particle_MPI_Vars       ,ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Partner identification
INTEGER                        :: iPeriodicVector,jPeriodicVector,iPeriodicDir
INTEGER,DIMENSION(2)           :: DirPeriodicVector = [-1,1]
INTEGER                        :: iElem,ElemID,firstElem,lastElem,NbElemID
INTEGER                        :: iSide,SideID,iLocSide
!INTEGER                        :: firstSide,lastSide
INTEGER                        :: iMortar,nMortarElems,NbSideID
INTEGER                        :: iProc,HaloProc
!INTEGER                        :: GlobalProcID
INTEGER                        :: nExchangeSides
!INTEGER                        :: nExchangeProcs
INTEGER,ALLOCATABLE            :: ExchangeSides(:)
REAL,ALLOCATABLE               :: BoundsOfElemCenter(:),MPISideBoundsOfElemCenter(:,:)
INTEGER                        :: ExchangeProcLeader
! Non-symmetric particle exchange
INTEGER,ALLOCATABLE            :: SendRequest(:),RecvRequest(:)
LOGICAL,ALLOCATABLE            :: GlobalProcToRecvProc(:)
LOGICAL                        :: CommFlag
INTEGER                        :: nNonSymmetricExchangeProcs,nNonSymmetricExchangeProcsGlob
INTEGER                        :: nExchangeProcessorsGlobal
!=================================================================================================================================

! Keep everything in sync here
! CALL MPI_BARRIER(MPI_COMM_FLEXI,IERROR)

!SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors...'

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

! Identify all procs with elements in range. This includes checking the procs on the compute-node as they might lie far apart
IF (nProcessors.GT.1) THEN
  !> Count all MPI sides on current proc.
  firstElem = offsetElem+1
  lastElem  = offsetElem+nElems

  ! This approach does not work, we get only the MPI sides pointing into the compute-node halo region
  !DO iSide = firstSide,lastSide
  !  IF (SideInfo_Shared(SIDE_NBELEMTYPE,iSide).EQ.2) THEN
  !    nExchangeSides = nExchangeSides + 1
  !  END IF
  !END DO

  !>>> For all element, loop over the six sides and check if the neighbor element is on the current proc
  !>>> Special care for big mortar sides, here the SIDE_ELEMID must be used
  nExchangeSides = 0

  DO iElem = firstElem,lastElem
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

          NbElemID = SideInfo_Shared(SIDE_ELEMID,NbSideID)
          ! If small mortar element not defined, skip it for now, likely not inside the halo region
          IF (NbElemID.LT.1) CYCLE

          ! If any of the small mortar sides is not on the local proc, the side is a MPI side
          IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
            nExchangeSides = nExchangeSides + 1
            EXIT
          END IF
        END DO

      ! regular side or small mortar side
      ELSE
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides = nExchangeSides + 1
        END IF
      END IF
    END DO
  END DO

  IF (nComputeNodeProcessors.GT.1.AND.nExchangeSides.EQ.0) &
    CALL ABORT(__STAMP__,'Found no side connectivity between processor domains')

  !> Build mapping for all MPI sides on current proc
  ALLOCATE(ExchangeSides(1:nExchangeSides))

  nExchangeSides = 0

  ! This approach does not work, we get only the MPI sides pointing into the compute-node halo region
  !DO iSide = firstSide,lastSide
  !  IF (SideInfo_Shared(SIDE_NBELEMTYPE,iSide).EQ.2) THEN
  !    nExchangeSides = nExchangeSides + 1
  !    ExchangeSides(nExchangeSides) = SideInfo_Shared(SIDE_ID,iSide)
  !  END IF
  !END DO

  DO iElem = firstElem,lastElem
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

          NbElemID = SideInfo_Shared(SIDE_ELEMID,NbSideID)
          ! If small mortar element not defined, skip it for now, likely not inside the halo region
          IF (NbElemID.LT.1) CYCLE

          ! If any of the small mortar sides is not on the local proc, the side is a MPI side
          IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
            nExchangeSides = nExchangeSides + 1
            ExchangeSides(nExchangeSides) = SideID
            EXIT
          END IF
        END DO

      ! regular side or small mortar side
      ELSE
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides = nExchangeSides + 1
          ExchangeSides(nExchangeSides) = SideID
        END IF
      END IF
    END DO
  END DO

  !> Build metrics for all MPI sides on current proc
  ALLOCATE(BoundsOfElemCenter(1:4))
  ALLOCATE(MPISideBoundsOfElemCenter(1:4,1:nExchangeSides))

  DO iSide = 1, nExchangeSides
    SideID = ExchangeSides(iSide)
    ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
    MPISideBoundsOfElemCenter(1:3,iSide) = (/ SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
                                              SUM(BoundsOfElem_Shared(1:2,2,ElemID)), &
                                              SUM(BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
    MPISideBoundsOfElemCenter(4,iSide) = VECNORM ((/BoundsOfElem_Shared(2,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                                    BoundsOfElem_Shared(2,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                                    BoundsOfElem_Shared(2,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
  END DO

  !> Check all elements in the CN halo region against local MPI sides. Check is identical to particle_bgm.f90
  !>>> Check the bounding box of each element in compute-nodes' halo domain against the bounding boxes of the elements of the
  !>>> MPI-surface (local proc MPI sides)

  ! Use a named loop so the entire element can be cycled
ElemLoop:  DO iElem = 1,nComputeNodeTotalElems
    ElemID   = GetGlobalElemID(iElem)
    HaloProc = ElemInfo_Shared(ELEM_RANK,ElemID)

    IF (HaloProc.EQ.myRank) CYCLE

!#if CODE_ANALYZE
!    ! Sanity checks. Elems in halo region must have ELEM_HALOFLAG=2 and the proc must not be flagged yet
!    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.2) THEN
!      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Flag: ',ElemInfo_Shared(ELEM_HALOFLAG,ElemID)
!      CALL ABORT(__STAMP__,  'Element found in range of halo elements while not flagged as such!')
!    END IF
!
!    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.1) THEN
!      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Proc: ',HaloProc
!      CALL ABORT(__STAMP__, 'Proc claimed to have elements both on compute node and in halo region!')
!    END IF
!#endif

    ! Skip if the proc is already flagged
    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).NE.-1) CYCLE

    BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
                                SUM(BoundsOfElem_Shared(1:2,2,ElemID)), &
                                SUM(BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
    BoundsOfElemCenter(4)   = VECNORM ((/ BoundsOfElem_Shared(2,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                          BoundsOfElem_Shared(2,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                          BoundsOfElem_Shared(2,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
    DO iSide = 1, nExchangeSides
      ! compare distance of centers with sum of element outer radii+halo_eps
      IF (VECNORM(BoundsOfElemCenter(1:3)-MPISideBoundsOfElemCenter(1:3,iSide)) &
        .GT. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide)) THEN

        ! Also check periodic directions. Only MPI sides of the local proc are
        ! taken into account, so do not perform additional case distinction
        SELECT CASE(GEO%nPeriodicVectors)
          ! One periodic vector
          CASE(1)
            DO iPeriodicDir = 1,2
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                       &
                         + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)                  &
                         - MPISideBoundsOfElemCenter(1:3,iSide))                                         &
                .LE. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                  ! flag the proc as exchange proc (in halo region)
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  CYCLE ElemLoop
                END IF
              END DO

          ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
          CASE(2)
            DO iPeriodicVector = 1,2
              DO iPeriodicDir = 1,2
                ! check if element is within halo_eps of periodically displaced element
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                           + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                           - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                          .LE. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                  ! flag the proc as exchange proc (in halo region)
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  CYCLE ElemLoop
                END IF

                DO jPeriodicVector = 1,2

                  ! check if element is within halo_eps of periodically displaced element
                  IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                             + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                             + GEO%PeriodicVectors(1:3,jPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                             - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                          .LE. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                    ! flag the proc as exchange proc (in halo region)
                    GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                    GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                    nExchangeProcessors = nExchangeProcessors + 1
                    CYCLE ElemLoop
                  END IF
                END DO
              END DO
            END DO

          ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
          CASE(3)
            ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
            ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
            DO iPeriodicVector = 1,3
              DO iPeriodicDir = 1,2
                ! element might be already added back
                ! check if element is within halo_eps of periodically displaced element
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                      &
                           + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)   &
                           - MPISideBoundsOfElemCenter(1:3,iSide))                                        &
                          .LE. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                  ! flag the proc as exchange proc (in halo region)
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  CYCLE ElemLoop
                END IF

                DO jPeriodicVector = 1,3
                  IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

                  ! check if element is within halo_eps of periodically displaced element
                  IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                             + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                             + GEO%PeriodicVectors(1:3,jPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                             - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                          .LE. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                    ! flag the proc as exchange proc (in halo region)
                    GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                    GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                    nExchangeProcessors = nExchangeProcessors + 1
                    CYCLE ElemLoop
                  END IF

                END DO
              END DO
            END DO

            ! check if element is within halo_eps of periodically displaced element
            DO iPeriodicDir = 1,2
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                        &
                         + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)                   &
                         + GEO%PeriodicVectors(1:3,2) * DirPeriodicVector(iPeriodicDir)                   &
                         + GEO%PeriodicVectors(1:3,3) * DirPeriodicVector(iPeriodicDir)                   &
                         - MPISideBoundsOfElemCenter(1:3,iSide))                                          &
                      .LE. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                ! flag the proc as exchange proc (in halo region)
                GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                nExchangeProcessors = nExchangeProcessors + 1
                CYCLE ElemLoop
              END IF
            END DO

          ! No periodic vectors, element out of range
          CASE(0)
            ! Do nothing

          CASE DEFAULT
            CALL ABORT(__STAMP__,'Invalid number of periodic vectors in particle_mpi_halo.f90')

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
END IF

! Communicate non-symmetric part exchange partners to catch non-symmetric proc identification due to inverse distance calculation
ALLOCATE(GlobalProcToRecvProc(0:nProcessors_Global-1), &
         SendRequest         (0:nProcessors_Global-1), &
         RecvRequest         (0:nProcessors_Global-1))

GlobalProcToRecvProc = .FALSE.

! Notify every proc which was identified by the local proc
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  CALL MPI_IRECV( GlobalProcToRecvProc(iProc)  &
                , 1                            &
                , MPI_INTEGER                  &
                , iProc                        &
                , 1999                         &
                , MPI_COMM_WORLD               &
                , RecvRequest(iProc)           &
                , IERROR)

  ! CommFlag holds the information if the local proc wants to communicate with iProc
  CommFlag = MERGE(.TRUE.,.FALSE.,GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).NE.-1)
  CALL MPI_ISEND( CommFlag                     &
                , 1                            &
                , MPI_INTEGER                  &
                , iProc                        &
                , 1999                         &
                , MPI_COMM_WORLD               &
                , SendRequest(iProc)           &
                , IERROR)
END DO

! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

! Append previously not found procs to list of exchange processors
nNonSymmetricExchangeProcs = 0
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  ! Ignore procs that are already flagged or not requesting communication
  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) .NE.-1) CYCLE
  IF (.NOT.GlobalProcToRecvProc(iProc)) CYCLE

  ! Found a previously missing proc
  nNonSymmetricExchangeProcs = nNonSymmetricExchangeProcs + 1
  nExchangeProcessors        = nExchangeProcessors + 1

  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
END DO

DEALLOCATE(GlobalProcToRecvProc,RecvRequest,SendRequest)

! On smooth grids, nNonSymmetricExchangeProcs should be zero. Only output if previously missing particle exchange procs are found
CALL MPI_REDUCE(nNonSymmetricExchangeProcs,nNonSymmetricExchangeProcsGlob,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
IF (nNonSymmetricExchangeProcsGlob.GT.1) THEN
  SWRITE(Unit_StdOut,'(A,I0,A)') ' | Found ',nNonSymmetricExchangeProcsGlob, &
                                 ' previously missing non-symmetric particle exchange procs'
END IF

! Build reverse mapping
!-- EXCHANGE_PROC_TYPE information is currently unused and either -1 (no communication) or 2 (communication). Can be used to
!-- implement check if exchange partner is on the same compute node, so build it here
ALLOCATE(ExchangeProcToGlobalProc(2,0:nExchangeProcessors-1))

! Loop through all procs and build reverse mapping
!DO iProc = 0,nComputeNodeProcessors-1
!  ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,iProc) = 1
!  ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc) = iProc + ComputeNodeRootRank
!END DO
!nExchangeProcessors = nComputeNodeProcessors

nExchangeProcessors = 0
DO iProc = 0,nProcessors_Global-1
  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).NE.-1) THEN
    ! Find it the other proc is on the same compute node
    ExchangeProcLeader = INT(GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)/nComputeNodeProcessors)
    IF (ExchangeProcLeader.EQ.myLeaderGroupRank) THEN
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,nExchangeProcessors) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,nExchangeProcessors) = iProc
    ELSE
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,nExchangeProcessors) = 2
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,nExchangeProcessors) = iProc
    END IF
    nExchangeProcessors = nExchangeProcessors + 1
  END IF
END DO

! -- Average number of exchange processors
CALL MPI_REDUCE(nExchangeProcessors,nExchangeProcessorsGlobal,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
SWRITE(UNIT_stdOut,'(A,I0,A)') ' | Started particle exchange communication with average ', &
                                 nExchangeProcessorsGlobal/nProcessors_Global            , &
                                 ' partners per proc'

SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

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


!PURE FUNCTION HaloBoxInProc(CartNodes,CartProc,halo_eps)
!!===================================================================================================================================
!! check if bounding box is on proc
!!===================================================================================================================================
!! MODULES
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)       :: CartNodes(6)
!REAL,INTENT(IN)       :: CartProc( 6)
!REAL,INTENT(IN)       :: halo_eps
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!LOGICAL               :: HaloBoxInProc
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER               :: iNode
!REAL,DIMENSION(1:3,8) :: xCoordsProc,xCordsTest
!!===================================================================================================================================
!
!HaloBoxInProc = .FALSE.
!
!! Reconstruct the eight corner nodes of the cuboid
!xCoordsProc(1:3,1) = (/CartProc(1), CartProc(3), CartProc(5)/)
!xCoordsProc(1:3,2) = (/CartProc(2), CartProc(3), CartProc(5)/)
!xCoordsProc(1:3,3) = (/CartProc(2), CartProc(4), CartProc(5)/)
!xCoordsProc(1:3,4) = (/CartProc(1), CartProc(4), CartProc(5)/)
!xCoordsProc(1:3,5) = (/CartProc(1), CartProc(3), CartProc(6)/)
!xCoordsProc(1:3,6) = (/CartProc(2), CartProc(3), CartProc(6)/)
!xCoordsProc(1:3,7) = (/CartProc(2), CartProc(4), CartProc(6)/)
!xCoordsProc(1:3,8) = (/CartProc(1), CartProc(4), CartProc(6)/)
!
!xCordsTest(1:3,1) = (/CartNodes(1), CartNodes(3), CartNodes(5)/)
!xCordsTest(1:3,2) = (/CartNodes(2), CartNodes(3), CartNodes(5)/)
!xCordsTest(1:3,3) = (/CartNodes(2), CartNodes(4), CartNodes(5)/)
!xCordsTest(1:3,4) = (/CartNodes(1), CartNodes(4), CartNodes(5)/)
!xCordsTest(1:3,5) = (/CartNodes(1), CartNodes(3), CartNodes(6)/)
!xCordsTest(1:3,6) = (/CartNodes(2), CartNodes(3), CartNodes(6)/)
!xCordsTest(1:3,7) = (/CartNodes(2), CartNodes(4), CartNodes(6)/)
!xCordsTest(1:3,8) = (/CartNodes(1), CartNodes(4), CartNodes(6)/)
!
!! Check if any of the eight test corner nodes is within the current proc
!DO iNode = 1,8
!  IF (   ((xCordsTest (1,iNode).LE.CartProc (2)+halo_eps).AND.(xCordsTest (1,iNode).GE.CartProc (1)-halo_eps))  &
!    .AND.((xCordsTest (2,iNode).LE.CartProc (4)+halo_eps).AND.(xCordsTest (2,iNode).GE.CartProc (3)-halo_eps))  &
!    .AND.((xCordsTest (3,iNode).LE.CartProc (6)+halo_eps).AND.(xCordsTest (3,iNode).GE.CartProc (5)-halo_eps))) THEN
!    HaloBoxInProc = .TRUE.
!  END IF
!END DO
!
!! Reverse check if any of the proc corner nodes is within the test proc
!DO iNode = 1,8
!  IF (   ((xCoordsProc(1,iNode).LE.CartNodes(2)+halo_eps).AND.(xCoordsProc(1,iNode).GE.CartNodes(1)-halo_eps))  &
!    .AND.((xCoordsProc(2,iNode).LE.CartNodes(4)+halo_eps).AND.(xCoordsProc(2,iNode).GE.CartNodes(3)-halo_eps))  &
!    .AND.((xCoordsProc(3,iNode).LE.CartNodes(6)+halo_eps).AND.(xCoordsProc(3,iNode).GE.CartNodes(5)-halo_eps))) THEN
!    HaloBoxInProc = .TRUE.
!  END IF
!END DO
!
!END FUNCTION HaloBoxInProc
#endif /*MPI*/

END MODULE MOD_Particle_MPI_Halo
