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
USE MOD_Particle_Globals        ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeElems
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,SideInfo_Shared,BoundsOfElem_Shared
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID,GetGlobalNonUniqueSideID
USE MOD_Particle_MPI_Vars       ,ONLY: halo_eps
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputenodeProcessors,nProcessors_Global,ComputeNodeRootRank
USE MOD_Particle_MPI_Vars       ,ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPeriodicVector,jPeriodicVector,iPeriodicDir
INTEGER,DIMENSION(2)           :: DirPeriodicVector = [-1,1]
INTEGER                        :: iElem,ElemID,firstElem,lastElem,NbElemID
INTEGER                        :: iSide,SideID,iLocSide !,firstSide,lastSide
INTEGER                        :: iMortar,nMortarElems,NbSideID
INTEGER                        :: iProc,HaloProc !,GlobalProcID
INTEGER                        :: nExchangeSides !,nExchangeProcs
INTEGER,ALLOCATABLE            :: ExchangeSides(:)
REAL,ALLOCATABLE               :: BoundsOfElemCenter(:),MPISideBoundsOfElemCenter(:,:)
!=================================================================================================================================

!SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors...'

! Allocate arrays
ALLOCATE(GlobalProcToExchangeProc(EXCHANGE_PROC_SIZE,0:nProcessors_Global-1))
GlobalProcToExchangeProc(:,:) = -1

! Identify all procs on same node
nExchangeProcessors = 0
DO iProc = ComputeNodeRootRank,ComputeNodeRootRank+nComputeNodeProcessors-1
  ! Build mapping global to compute-node
  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
  nExchangeProcessors = nExchangeProcessors + 1
END DO

! Identify all procs with elements in range. If all elements are on the current proc, they are already added
! and we are done here
IF (nComputeNodeElems.NE.nComputeNodeTotalElems) THEN
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
    MPISideBoundsOfElemCenter(4,iSide) = VECNORM ((/ BoundsOfElem_Shared(2,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                                     BoundsOfElem_Shared(2,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                                     BoundsOfElem_Shared(2,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
  END DO

  !> Check all elements in the CN halo region against local MPI sides. Check is identical to particle_bgm.f90
  !>>> Check the bounding box of each element in compute-nodes' halo domain against the bounding boxes of the
  !>>> of the elements of the MPI-surface (local proc MPI sides)

  ! Use a named loop so the entire element can be cycled
ElemLoop:  DO iElem = nComputeNodeElems+1,nComputeNodeTotalElems
    ElemID   = GetGlobalElemID(iElem)
    HaloProc = ElemInfo_Shared(ELEM_RANK,ElemID)

!#if CODE_ANALYZE
    ! Sanity checks. Elems in halo region must have ELEM_HALOFLAG=2 and the proc must not be flagged yet
    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.2) THEN
      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Flag: ',ElemInfo_Shared(ELEM_HALOFLAG,ElemID)
      CALL ABORT(__STAMP__,  'Element found in range of halo elements while not flagged as such!')
    END IF

    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.1) THEN
      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Proc: ',HaloProc
      CALL ABORT(__STAMP__, 'Proc claimed to have elements both on compute node and in halo region!')
    END IF
!#endif

    ! Skip if the proc is already flagged
    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE

    BoundsOfElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
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

!
ALLOCATE(ExchangeProcToGlobalProc(2,0:nExchangeProcessors-1))

! Loop through all procs and build reverse mapping
DO iProc = 0,nComputeNodeProcessors-1
  ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,iProc) = 1
  ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc) = iProc + ComputeNodeRootRank
END DO

nExchangeProcessors = nComputeNodeProcessors
DO iProc = 0,nProcessors_Global-1
  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).EQ.2) THEN
    ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,nExchangeProcessors) = 2
    ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,nExchangeProcessors) = iProc
    nExchangeProcessors = nExchangeProcessors +1
  END IF
END DO

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


#endif /*MPI*/

END MODULE MOD_Particle_MPI_Halo
