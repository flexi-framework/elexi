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
! Module contains the routines to partition the mesh according to load balancing
!===================================================================================================================================
MODULE MOD_LoadMesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE LoadPartition
  MODULE PROCEDURE LoadPartition
END INTERFACE

INTERFACE LoadElemTime
  MODULE PROCEDURE LoadElemTime
END INTERFACE

PUBLIC  :: LoadPartition
PUBLIC  :: LoadElemTime
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

#if USE_MPI
!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE LoadPartition(FileString)
! MODULES                                                                                                                         !
USE MOD_Globals
USE MOD_LoadDistribution
USE MOD_LoadBalance_Vars
USE MOD_HDF5_Input
USE MOD_Mesh_Vars,      ONLY: nElems,nGlobalElems,offsetElem
USE MOD_Restart_Vars,   ONLY: DoRestart,RestartFile
USE MOD_StringTools,    ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileString                          !> mesh filename
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nVal(15),iVar
REAL,ALLOCATABLE               :: ElemTime_local(:)
REAL,ALLOCATABLE               :: ElemData_loc(:,:),tmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesElemData_loc(:)
!===================================================================================================================================
! Fresh computation cannot have load inbalance
IF (.NOT.DoRestart) RETURN

!----------------------------------------------------------------------------------------------------------------------------------!
IF(MPIRoot) THEN
    WRITE(UNIT_stdOut,'(A,I0,A)',ADVANCE='YES') ' Starting load partitioning on ',nProcessors,' processors.'
END IF

! Clean up variables usually done in mesh_reading
SDEALLOCATE(LoadDistri)
ALLOCATE(   LoadDistri(0:nProcessors-1))
LoadDistri(:) = 0.

SDEALLOCATE(PartDistri)
ALLOCATE(   PartDistri(0:nProcessors-1))
PartDistri(:) = 0

ElemTimeExists=.FALSE.

!----------------------------------------------------------------------------------------------------------------------------------!
! Readin of ElemTime: Read in only by MPIRoot in single mode, only communicate logical ElemTimeExists
SDEALLOCATE(ElemGlobalTime)
ALLOCATE(   ElemGlobalTime(1:nGlobalElems))
ElemGlobalTime = 0.

! Close the still open mesh file
CALL CloseDataFile()

! 1) Only MPIRoot does readin of ElemTime
IF(MPIRoot)THEN
    ALLOCATE(ElemTime_local(1:nGlobalElems))
    ElemTime_local = 0.0
    nElems         = nGlobalElems    ! Temporary set nElems as nGlobalElems for GetArrayAndName
    offsetElem     = 0               ! Offset is the index of first entry, hdf5 array starts at 0-.GT. -1

    ! Ready to get the restart data
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
    CALL GetArrayAndName('ElemData','VarNamesAdd',nVal,tmp,VarNamesElemData_loc)
    CALL CloseDataFile()

    ! Check if restart file has load balancing data
    IF (ALLOCATED(VarNamesElemData_loc)) THEN
        ALLOCATE(ElemData_loc(nVal(1),nVal(2)))
        ElemData_loc = RESHAPE(tmp,(/nVal(1),nVal(2)/))

        ! Search for ElemTime
        DO iVar=1,nVal(1)
            IF (STRICMP(VarNamesElemData_loc(iVar),"ElemTime")) THEN
                ElemTime_local = REAL(ElemData_loc(iVar,:))
                ElemTimeExists = .TRUE.
            END IF
        END DO

        ! De-allocate temporary variables
        DEALLOCATE(ElemData_loc,VarNamesElemData_loc,tmp)
    END IF

    ! Set the global elem time
    ElemGlobalTime = ElemTime_local
    DEALLOCATE(ElemTime_local)

    ! if the elemtime is 0.0, the value must be changed in order to prevent a division by zero
    IF(MAXVAL(ElemGlobalTime).LE.0.0) THEN
        ElemGlobalTime = 1.0
        ElemTimeExists = .FALSE.
    END IF
END IF

! 2) Distribute logical information ElemTimeExists
CALL MPI_BCAST (ElemTimeExists,1,MPI_LOGICAL,0,MPI_COMM_FLEXI,iError)

! Distribute the elements according to the selected distribution method
CALL ApplyWeightDistributionMethod(ElemTimeExists)

! Restore the mesh file pointer
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
CHECKSAFEINT(HSize(2),4)
DEALLOCATE(HSize)

! Wait for all procs before communicating the final load partitioning
IF(MPIRoot) THEN
    WRITE(UNIT_StdOut,'(132("-"))')
END IF

END SUBROUTINE LoadPartition


SUBROUTINE LoadElemTime()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_LoadDistribution
USE MOD_LoadBalance_Vars
USE MOD_IO_HDF5,        ONLY:ElementOut,AddToElemData
USE MOD_Mesh_Vars,      ONLY:nElems
USE MOD_MPI_Vars,       ONLY:offsetElemMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iProc
REAL,ALLOCATABLE  :: WeightSum_proc(:)
!===================================================================================================================================

! Set new ElemTime depending on new load distribution
SDEALLOCATE(ElemTime)
ALLOCATE(ElemTime(nElems))
ElemTime = 0.

! Write it in every state file
CALL AddToElemData(ElementOut,'ElemTime',RealArray=ElemTime)

! Calculate new (theoretical) imbalance with offsetElemMPI information
IF(ElemTimeExists.AND.MPIRoot)THEN
    ALLOCATE(WeightSum_proc(0:nProcessors-1))
    DO iProc=0,nProcessors-1
        WeightSum_proc(iProc) = SUM(ElemGlobalTime(1+offsetElemMPI(iProc):offsetElemMPI(iProc+1)))
    END DO
    MaxWeight = MAXVAL(WeightSum_proc)
    MinWeight = MINVAL(WeightSum_proc)
    ! WeightSum (Mesh global value) is already set in BalanceMethod scheme

    ! new computation of current imbalance
    TargetWeight=SUM(WeightSum_proc)/nProcessors
    NewImbalance =  (MaxWeight-TargetWeight ) / TargetWeight
    IF(TargetWeight.LE.0.0) CALL abort(__STAMP__,' LoadBalance: TargetWeight = ',RealInfo=TargetWeight)

    SWRITE(UNIT_stdOut,'(A)') ' Calculated new (theoretical) imbalance with offsetElemMPI information'
    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MaxWeight:        ', MaxWeight
    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MinWeight:        ', MinWeight
    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' TargetWeight:     ', TargetWeight
    SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' NewImbalance:     ', NewImbalance
ELSE
    SWRITE(UNIT_stdOut,'(A)') ' No ElemTime found in restart file'
    NewImbalance = -1.
    MaxWeight = -1.
    MinWeight = -1.
END IF

SDEALLOCATE(ElemGlobalTime)

! Re-allocate nPartsPerElem depending on new number of elements
IF(.NOT.ALLOCATED(nPartsPerElem))THEN
    ALLOCATE(nPartsPerElem(1:nElems))
ELSE
    SDEALLOCATE(nPartsPerElem)
    ALLOCATE(nPartsPerElem(1:nElems))
END IF
nPartsPerElem=0
CALL AddToElemData(ElementOut,'nPartsPerElem',IntArray=nPartsPerElem(:))

SDEALLOCATE(nTracksPerElem)
ALLOCATE(nTracksPerElem(1:nElems))
nTracksPerElem=0

! Clean up local variables
SDEALLOCATE(WeightSum_proc)

IF(MPIRoot) THEN
    WRITE(UNIT_StdOut,'(132("-"))')
END IF

END SUBROUTINE LoadElemTime
#endif /*MPI*/

END MODULE MOD_LoadMesh
