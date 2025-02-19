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

!===================================================================================================================================
!> Contains the routines to exchange data using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_MPI_Shared
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
#if DEBUG_MEMORY
INTERFACE Allocate_Shared_DEBUG
#else
INTERFACE Allocate_Shared
#endif /*DEBUG_MEMORY*/
  MODULE PROCEDURE Allocate_Shared_Logical_1
  MODULE PROCEDURE Allocate_Shared_Logical_2
  MODULE PROCEDURE Allocate_Shared_Int_1
  MODULE PROCEDURE Allocate_Shared_Int_2
  MODULE PROCEDURE Allocate_Shared_Int_3
  MODULE PROCEDURE Allocate_Shared_Int_4
  MODULE PROCEDURE Allocate_Shared_Real_1
  MODULE PROCEDURE Allocate_Shared_Real_2
  MODULE PROCEDURE Allocate_Shared_Real_3
  MODULE PROCEDURE Allocate_Shared_Real_4
  MODULE PROCEDURE Allocate_Shared_Real_5
  MODULE PROCEDURE Allocate_Shared_Real_6
END INTERFACE

PUBLIC:: DefineParametersMPIShared
PUBLIC:: InitMPIShared
PUBLIC:: FinalizeMPIShared
#if DEBUG_MEMORY
PUBLIC:: Allocate_Shared_DEBUG
#else
PUBLIC:: Allocate_Shared
#endif /*DEBUG_MEMORY*/
! PUBLIC:: UpdateDGShared
PUBLIC:: BARRIER_AND_SYNC
PUBLIC:: MPI_SIZE
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE DefineParametersMPIShared()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools               ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection     ('MPI Shared')
#if ! (CORE_SPLIT==0)
CALL prms%CreateIntOption('NbrOfPhysicalNodes' , 'Number of physical nodes (as opposed to virtual nodes). Required for RAM'    //&
                                                 ' monitoring when using more than one virtual compute node per physical node.','-1')
#endif /*! (CORE_SPLIT==0)*/

END SUBROUTINE DefineParametersMPIShared


!==================================================================================================================================
!> Initialize for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE InitMPIShared()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: MemoryMonitor
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
USE MOD_ReadInTools               ,ONLY: GETINT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: sharedGroup
INTEGER                         :: worldGroup
INTEGER                         :: i
INTEGER                         :: color
!==================================================================================================================================

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MPI SHARED COMMUNICATION...'

! Save the global number of procs
nProcessors_Global = nProcessors

MemoryMonitor = .TRUE.
#if ! (CORE_SPLIT==0)
! When core-level splitting is used, it is not clear how many cores are on the same physical compute node.
#if USE_MPI
NbrOfPhysicalNodes =  GETINT('NbrOfPhysicalNodes')
IF(NbrOfPhysicalNodes.LE.0) MemoryMonitor = .FALSE.
#endif /*USE_MPI*/
#endif /*! (CORE_SPLIT==0)*/

! Split the node communicator (shared memory) from the global communicator on physical processor or node level
#if (CORE_SPLIT==1)
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,myRank,0,MPI_COMM_SHARED,iError)
#elif (CORE_SPLIT==0)
  ! Note that using SharedMemoryMethod=OMPI_COMM_TYPE_CORE somehow does not work in every case (intel/amd processors)
  ! Also note that OMPI_COMM_TYPE_CORE is undefined when not using OpenMPI
  CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_FLEXI,SHARED_MEMORY_METHOD,0,MPI_INFO_NULL,MPI_COMM_SHARED,IERROR)
#else
  ! Check if more nodes than procs are required or
  ! if the resulting split would create unequal procs per node
  IF((CORE_SPLIT.GT.nProcessors_Global).OR.(MOD(nProcessors_Global,CORE_SPLIT).GT.0))THEN
    SWRITE (UNIT_stdOut,'(A,I0,A,I0,A,F0.2,A)') ' WARNING: Either more nodes than cores selected (nodes: ',CORE_SPLIT,', cores: ',&
        nProcessors_Global,') OR unequal number of cores per node (=',REAL(nProcessors_Global)/REAL(CORE_SPLIT),&
        '). Setting 1 core per node for MPI_COMM_SHARED!'
    color = myRank
  ELSE
    ! Group procs so that every CORE_SPLIT procs are in the same group
    color = INT(REAL(myRank)/REAL(CORE_SPLIT))
  END IF ! (CORE_SPLIT.GE.nProcessors_Global).OR.(MOD().GT.0)
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,color,0,MPI_COMM_SHARED,iError)
#endif

! Find my rank on the shared communicator, comm size and proc name
CALL MPI_COMM_RANK(MPI_COMM_SHARED, myComputeNodeRank,IERROR)
CALL MPI_COMM_SIZE(MPI_COMM_SHARED, nComputeNodeProcessors,IERROR)

! MPI3 shared implementation currently only works with equal procs per node
IF (MOD(nProcessors_Global,nComputeNodeProcessors).NE.0) &
  CALL Abort(__STAMP__,'MPI shared communication currently only supported with equal procs per node!')

IF (nProcessors_Global.EQ.nComputeNodeProcessors) THEN
  nComputeNodes = 1
  SWRITE(UNIT_stdOut,'(A,I0,A,I0,A)') ' | Starting shared communication with ',nComputeNodeProcessors,' procs on ',1,' node'
ELSE
  nComputeNodes = INT(nProcessors_Global/nComputeNodeProcessors)
  SWRITE(UNIT_stdOut,'(A,I0,A,I0,A,I0,A)') ' | Starting shared communication with ',nComputeNodeProcessors,' procs on ',         &
                                                         nComputeNodes     ,' nodes for a total number of ',&
                                                         nProcessors_Global,' procs'
END IF

! Map global rank number into shared rank number. Returns MPI_UNDEFINED if not on the same node
ALLOCATE(MPIRankGlobal(0:nProcessors-1))
ALLOCATE(MPIRankShared(0:nProcessors-1))
DO i=0,nProcessors-1
  MPIRankGlobal(i) = i
END DO

! Get handles for each group
sharedGroup = MPI_GROUP_NULL
worldGroup  = MPI_GROUP_NULL
CALL MPI_COMM_GROUP(MPI_COMM_FLEXI , worldGroup,IERROR)
CALL MPI_COMM_GROUP(MPI_COMM_SHARED,sharedGroup,IERROR)

! Finally translate global rank to local rank
CALL MPI_GROUP_TRANSLATE_RANKS(worldGroup,nProcessors,MPIRankGlobal,sharedGroup,MPIRankShared,IERROR)

! Free group handles
CALL MPI_GROUP_FREE( worldGroup,IERROR)
CALL MPI_GROUP_FREE(sharedGroup,IERROR)

! Send rank of compute node root to all procs on shared comm
IF (myComputeNodeRank.EQ.0) ComputeNodeRootRank = myRank
CALL MPI_BCAST(ComputeNodeRootRank,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS_SHARED=MPI_COMM_NULL
myLeaderGroupRank=-1
color = MERGE(101,MPI_UNDEFINED,myComputeNodeRank.EQ.0)
CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,color,0,MPI_COMM_LEADERS_SHARED,IERROR)
IF(myComputeNodeRank.EQ.0)THEN
  CALL MPI_COMM_RANK(MPI_COMM_LEADERS_SHARED,myLeaderGroupRank,IERROR)
  CALL MPI_COMM_SIZE(MPI_COMM_LEADERS_SHARED,nLeaderGroupProcs,IERROR)
END IF

! Leaders inform every processor on their node about the leader group rank and group size
CALL MPI_BCAST(myLeaderGroupRank,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)
CALL MPI_BCAST(nLeaderGroupProcs,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

! communicate ranks for each leader rank
ALLOCATE(MPIRankLeader(0:nLeaderGroupProcs-1))
IF(myComputeNodeRank.EQ.0)THEN
  CALL MPI_ALLGATHER(ComputeNodeRootRank,1,MPI_INTEGER,MPIRankLeader(0:nLeaderGroupProcs-1),1,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
END IF
CALL MPI_BCAST(MPIRankLeader,nLeaderGroupProcs,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

! Create MPI_Info for shared memory windows
CALL MPI_INFO_CREATE(MPI_INFO_SHARED_LOOSE,IERROR)
CALL MPI_INFO_SET(   MPI_INFO_SHARED_LOOSE,'accumulate_ordering','none',IERROR)
! Only root allocates, size differs between ranks
!CALL MPI_INFO_SET(   MPI_INFO_SHARED_LOOSE,'same_size'          ,'true',IERROR)
CALL MPI_INFO_SET(   MPI_INFO_SHARED_LOOSE,'same_disp_unit'     ,'true',IERROR)

! Create a group from the communicator
CALL MPI_COMM_GROUP(MPI_COMM_SHARED, MPI_GROUP_SHARED,IERROR)
IF(myComputeNodeRank.EQ.0)THEN
  CALL MPI_COMM_GROUP(MPI_COMM_LEADERS_SHARED, MPI_GROUP_LEADERS_SHARED,IERROR)
END IF

! synchronize everything or bad things will happen
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

MPISharedInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')      ' INIT MPI SHARED COMMUNICATION DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitMPIShared


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Logical_1(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
LOGICAL,INTENT(OUT),POINTER               :: DataPointer(:)           !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_LOGICAL_1'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Logical1) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')


END SUBROUTINE ALLOCATE_SHARED_LOGICAL_1

!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Logical_2(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
LOGICAL,INTENT(OUT),POINTER               :: DataPointer(:,:)         !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_LOGICAL_2'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Logical1) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Logical1) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_LOGICAL_2

!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Int_1(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:)           !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_INT_1'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Int1) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int1) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_INT_1


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================

SUBROUTINE Allocate_Shared_Int_2(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:,:)         !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_INT_2'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Int2) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int2) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_INT_2


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Int_3(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(3)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:,:,:)       !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_INT_3'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Int3) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int3) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_INT_3


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Int_4(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(4)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:)     !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_INT_4'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Int3) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int3) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_INT_4


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_1(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:)         !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_REAL_1'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Real1) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real1) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_REAL_1


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_2(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:)         !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_REAL_2'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Real2) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real2) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_REAL_2


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_3(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(3)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:)       !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_REAL_3'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Real3) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real3) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_REAL_3


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_4(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(4)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:)     !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_REAL_4'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Real4) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real4) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_REAL_4


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_5(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(5)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:,:)   !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_REAL_5'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Real5) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real5) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_REAL_5


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_6(nVal,SM_WIN,DataPointer&
#if DEBUG_MEMORY
        ,SM_WIN_NAME,SM_CALL_FILE,SM_CALL_FILE_LINE&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(6)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:,:,:) !> Pointer to the RMA window
#if DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
CHARACTER(LEN=*),INTENT(IN)               :: SM_CALL_FILE             !> File from which the call comes
INTEGER,INTENT(IN)                        :: SM_CALL_FILE_LINE        !> Line in file from which the call comes
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
CHARACTER(32),PARAMETER                   :: FuncName='ALLOCATE_SHARED_REAL_6'
#if DEBUG_MEMORY
CHARACTER(32)                             :: hilf                     !> Convert line number from int to string
#endif /*DEBUG_MEMORY*/
!==================================================================================================================================
IF (ASSOCIATED(DataPointer)) CALL Abort(__STAMP__,'ERROR: Datapointer (Real6) already associated')

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#if DEBUG_MEMORY
WRITE(UNIT=hilf,FMT='(I0)') SM_CALL_FILE_LINE
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real6) already associated'&
#if DEBUG_MEMORY
//' for '//TRIM(SM_WIN_NAME)//' in file '//TRIM(SM_CALL_FILE)//' in line '//TRIM(hilf)&
#endif /*DEBUG_MEMORY*/
)

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN, IERROR)
IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_ALLOCATE_SHARED returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR, IERROR)
  IF (IERROR.NE.MPI_SUCCESS)  CALL abort(__STAMP__,'ERROR: MPI_WIN_SHARED_QUERY returned IERROR.NE.MPI_SUCCESS ('//TRIM(FuncName)//')')
END IF
IF (WIN_SIZE.LE.0) CALL abort(__STAMP__,'ERROR: WIN_SIZE ('//TRIM(FuncName)//') is zero')

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer, nVal)
IF (.NOT.ASSOCIATED(DataPointer)) CALL abort(__STAMP__,'ERROR: Datapointer ('//TRIM(FuncName)//') could not be associated')

END SUBROUTINE ALLOCATE_SHARED_REAL_6


!==================================================================================================================================
!> Unlock and free shared memory array
!==================================================================================================================================
SUBROUTINE BARRIER_AND_SYNC(SharedWindow,Communicator) !,Barrier_Opt)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(INOUT)       :: SharedWindow !> Shared memory window
INTEGER,INTENT(INOUT)       :: Communicator !> Shared memory communicator
! LOGICAL,INTENT(IN)          :: Barrier_Opt  !
! LOCAL VARIABLES
! LOGICAL                     :: Barrier
!==================================================================================================================================
! Barrier = MERGE(Barrier_Opt,.TRUE.,PRESENT(Barrier_Opt)

CALL MPI_WIN_SYNC(SharedWindow,iError)
! IF (Barrier) CALL MPI_BARRIER (Communicator,iError)
CALL MPI_BARRIER (Communicator,iError)
CALL MPI_WIN_SYNC(SharedWindow,iError)

! IF(iError.NE.MPI_SUCCESS)THEN
!   CALL Abort(__STAMP__,'ERROR in MPI_WIN_SYNC() for '//TRIM(SM_WIN_NAME)//': iError returned non-zero value =',IntInfoOpt=iError)
! END IF ! iError.NE.MPI_SUCCESS

END SUBROUTINE BARRIER_AND_SYNC


!==================================================================================================================================
!> Finalize for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE FinalizeMPIShared()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Free arrays
SDEALLOCATE(MPIRankGlobal)
SDEALLOCATE(MPIRankShared)
SDEALLOCATE(MPIRankLeader)

! Free MPI_INFO objects
IF(MPI_INFO_SHARED_LOOSE   .NE.MPI_INFO_NULL ) CALL MPI_INFO_FREE( MPI_INFO_SHARED_LOOSE   ,IERROR)

! Free the shared communicator
IF(MPI_COMM_SHARED         .NE.MPI_COMM_NULL ) CALL MPI_COMM_FREE( MPI_COMM_SHARED         ,IERROR)
IF(MPI_COMM_LEADERS_SHARED .NE.MPI_COMM_NULL ) CALL MPI_COMM_FREE( MPI_COMM_LEADERS_SHARED ,IERROR)

! Free the shared groups
IF(MPI_GROUP_SHARED        .NE.MPI_GROUP_NULL) CALL MPI_GROUP_FREE(MPI_GROUP_SHARED        ,IERROR)
IF(MPI_GROUP_LEADERS_SHARED.NE.MPI_GROUP_NULL) CALL MPI_GROUP_FREE(MPI_GROUP_LEADERS_SHARED,IERROR)

MPISharedInitIsDone=.FALSE.

END SUBROUTINE FinalizeMPIShared


! SUBROUTINE UpdateDGShared(U)
! !===================================================================================================================================
! ! Updates the DG solution in the MPI-3 shared memory window in the current node
! !===================================================================================================================================
! ! MODULES
! USE MOD_Globals,            ONLY: iError
! USE MOD_Preproc,            ONLY: N
! USE MOD_Mesh_Vars,          ONLY: nElems,offsetElem
! USE MOD_MPI_Shared_Vars
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT/OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! REAL,INTENT(IN)                 :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER                         :: FirstElemShared,LastElemShared
! !=================================================================================================================================
!
! !> Calculate the local offset relative to the node MPI root
! FirstElemShared = offsetElem-offsetComputeNodeElem+1
! LastElemShared  = offsetElem-offsetComputeNodeElem+nElems
!
! !> Update the DG solution in the RMA window
! U_Shared(:,:,:,:,FirstElemShared:LastElemShared) = U(:,:,:,:,:)
!
! ! Synchronize all RMA communication
! CALL BARRIER_AND_SYNC(U_Shared_Win,MPI_COMM_SHARED)
! CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
! ! Intel documentation claims this is required on "certain architectures". Whatever this means...
! ! https://software.intel.com/en-us/articles/an-introduction-to-mpi-3-shared-memory-programming
! CALL BARRIER_AND_SYNC(U_Shared_Win,MPI_COMM_SHARED)
!
! END SUBROUTINE UpdateDGShared


!===================================================================================================================================
!
!===================================================================================================================================
FUNCTION MPI_SIZE(nVal,VarSize)
! MODULES
USE MOD_Globals
USE MOD_MPI_Shared_Vars              ,ONLY: myComputeNodeRank
USE MOD_Memory                       ,ONLY: ProcessMemUsage
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=DP),INTENT(IN)    :: nVal
INTEGER,INTENT(IN)             :: VarSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPI_SIZE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                      :: memory(3)
INTEGER(KIND=DP),PARAMETER                :: kByte = 1024
!===================================================================================================================================

IF (INT(nVal*INT(VarSize,KIND=8),KIND=8).LT.INT(HUGE(INT(1,KIND=MPI_ADDRESS_KIND)),KIND=8)) THEN
  ! Compute node root also checks available memory
  IF (myComputeNodeRank.EQ.0) THEN
    ! Find memory usage and requirements
    CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

    ! Compare requested size against available memory
    IF (INT(nVal*INT(VarSize,KIND=8),KIND=8).GT.memory(2)*kByte) &
      CALL Abort(__STAMP__,'Trying to allocate shared array larger than available memory!')
  END IF

  MPI_SIZE = INT(nVal,KIND=MPI_ADDRESS_KIND) * INT(VarSize,KIND=MPI_ADDRESS_KIND)
ELSE
  CALL Abort(__STAMP__,'MPI_SIZE for shared array too large!')
  ! Suppress compiler warning
  MPI_SIZE = -1
END IF

END FUNCTION MPI_SIZE
#endif /* USE_MPI */

END MODULE MOD_MPI_Shared
