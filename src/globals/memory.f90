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

!==================================================================================================================================
!> Provides memory routines to particle code
!============================================================================================================   ===================
MODULE MOD_Memory
!===================================================================================================================================
! MODULES
USE MOD_Globals
#if USE_MPI
USE MPI
#endif /*MPI*/
IMPLICIT NONE
PRIVATE
SAVE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL :: MemoryStatFileExists = .FALSE.

!=================================================================================================================================

INTERFACE ProcessMemUsage
  ! SUBROUTINE processmemusage(memUsed,memAvail,memTotal) BIND(C, name='processmemusage')
  !   USE ISO_C_BINDING,   ONLY : c_double
  !   real(c_double) :: memUsed
  !   real(c_double) :: memAvail
  !   real(c_double) :: memTotal
  ! END SUBROUTINE processmemusage
  MODULE PROCEDURE ProcessMemUsage
END INTERFACE

INTERFACE GetMemUsage
  MODULE PROCEDURE GetMemUsage
END INTERFACE

INTERFACE VerifyMemUsage
  MODULE PROCEDURE VerifyMemUsage
END INTERFACE

INTERFACE Allocate_Safe
  MODULE PROCEDURE Allocate_Safe_Logical_1
  MODULE PROCEDURE Allocate_Safe_Int_1
  MODULE PROCEDURE Allocate_Safe_Int_2
  MODULE PROCEDURE Allocate_Safe_Real_1
  MODULE PROCEDURE Allocate_Safe_Real_2
  MODULE PROCEDURE Allocate_Safe_Real_3
END INTERFACE Allocate_Safe

PUBLIC :: ProcessMemUsage
PUBLIC :: GetMemUsage
PUBLIC :: VerifyMemUsage
PUBLIC :: Allocate_Safe
!===================================================================================================================================

CONTAINS

SUBROUTINE ProcessMemUsage(memory)
!==================================================================================================================================
! processmemusage(REAL,REAL,REAL) - takes three REALS, attemps to read the system-dependent data for available and total memory as
! well as the system-dependent data for a process' unique set size (USS), and return the results in KB.
!
! On failure, returns -1 on the failed value
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_StringTools               ,ONLY: split_string,STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)                          :: memory(3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
REAL                                      :: memSize
INTEGER                                   :: stat,memstat,ioUnit
INTEGER                                   :: memCount
CHARACTER(LEN=255)                        :: memName
CHARACTER(LEN=255)                        :: memStrings(2),memArray(2)
!==================================================================================================================================

stat   = 0
memory = 0

! Open the magic symbolic link /proc/self, which resolves to the process's own /proc/[pid] directory
! /proc/pid/smaps_rollup provides pre-summed memory information for a process.
! > https://www.kernel.org/doc/Documentation/ABI/testing/procfs-smaps_rollup
OPEN(NEWUNIT= ioUnit                   , &
     FILE   = '/proc/self/smaps_rollup', &
     STATUS = 'OLD'                    , &
     ACTION = 'READ'                   , &
     ACCESS = 'SEQUENTIAL'             , &
     IOSTAT = stat)
IF (stat.NE.0) CALL Abort(__STAMP__,"Could not open /proc/self/smaps_rollup.")

! Skip the first line, containing the hash
READ (ioUnit,"(A)",IOSTAT=stat)

DO
  READ (ioUnit,"(A)",IOSTAT=stat) memName
  ! Split once to get the memory array name
  CALL split_string(memName      ,':',memStrings,memCount)
  ! Split again to get the memory array size
  CALL split_string(TRIM(ADJUSTL(memStrings(2))),' ',memArray,memCount)
  ! Read the memSize into an integer
  READ(memArray(1),*,IOSTAT=memstat) memSize
  IF (memstat.NE.0) EXIT
  IF (STRICMP(memStrings(1),'Private_Clean')) memory(1) = memory(1) + memSize
  IF (STRICMP(memStrings(1),'Private_Dirty')) memory(1) = memory(1) + memSize
  IF(IS_IOSTAT_END(stat)) EXIT
END DO
CLOSE(ioUnit)

! /proc/meminfo provides information about distribution and utilization of memory.
! > https://www.kernel.org/doc/Documentation/filesystems/proc.txt
OPEN(NEWUNIT= ioUnit                   , &
     FILE   = '/proc/meminfo'          , &
     STATUS = 'OLD'                    , &
     ACTION = 'READ'                   , &
     ACCESS = 'SEQUENTIAL'             , &
     IOSTAT = stat)
IF (stat.NE.0) CALL Abort(__STAMP__,"Could not open /proc/meminfo.")

DO
  READ (ioUnit,"(A)",IOSTAT=stat) memName
  ! Split once to get the memory array name
  CALL split_string(memName      ,':',memStrings,memCount)
  ! Split again to get the memory array size
  CALL split_string(TRIM(ADJUSTL(memStrings(2))),' ',memArray,memCount)
  ! Read the memSize into an integer
  READ(memArray(1),*,IOSTAT=memstat) memSize
  IF (memstat.NE.0) EXIT
  IF (STRICMP(memStrings(1),'MemTotal'    )) memory(3) = memSize
  IF (STRICMP(memStrings(1),'MemAvailable')) memory(2) = memSize
  IF(IS_IOSTAT_END(stat)) EXIT
END DO
CLOSE(ioUnit)

END SUBROUTINE ProcessMemUsage


SUBROUTINE GetMemUsage(memory,caller)
!==================================================================================================================================
!> Verifies sufficient memory is available to allocate
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: StartTime
USE MOD_StringTools               ,ONLY: split_string,STRICMP
#if USE_MPI
USE MOD_MPI_Shared_Vars           ,ONLY: myComputeNodeRank,myLeaderGroupRank
! USE MOD_MPI_Shared_Vars           ,ONLY: nComputeNodeProcessors
USE MOD_MPI_Shared_Vars           ,ONLY: MPI_COMM_LEADERS_SHARED,MPI_COMM_SHARED
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: nLoadBalanceSteps
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)                          :: memory(3)
CHARACTER(LEN=*),INTENT(IN)               :: caller
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
REAL                                      :: buffers,cached,memSize
INTEGER                                   :: stat,memstat,ioUnit
INTEGER                                   :: memCount
CHARACTER(LEN=255)                        :: memName
CHARACTER(LEN=255)                        :: memStrings(2),memArray(2)
REAL                                      :: currentT
#if USE_MPI
REAL                                      :: ProcMemoryUsed    ! Used memory on a single proc
REAL                                      :: NodeMemoryUsed    ! Sum of used memory across one compute node
#endif /*USE_MPI*/
!==================================================================================================================================

GETTIME(CurrentT)
! IF (nProcessors.NE.nComputeNodeProcessors) &
!   CALL Abort(__STAMP__,'GetMemUsage routine only suitable for single node executation!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

! Only CN roots communicate available and total memory info (count once per node)
#if USE_MPI
IF (nProcessors.GT.1) THEN
  ! Collect data on node roots
  ProcMemoryUsed = memory(1)
  IF (myComputeNodeRank.EQ.0) THEN
    CALL MPI_REDUCE(ProcMemoryUsed , NodeMemoryUsed , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_SHARED , IERROR)
    memory(1) = NodeMemoryUsed
  ELSE
    CALL MPI_REDUCE(ProcMemoryUsed , 0              , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_SHARED , IERROR)
  END IF

  ! collect data from node roots on first root node
  IF (myComputeNodeRank.EQ.0) THEN ! only leaders
    IF (myLeaderGroupRank.EQ.0) THEN ! first node leader MUST be MPIRoot
      CALL MPI_REDUCE(MPI_IN_PLACE , memory , 3 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_LEADERS_SHARED , IERROR)
    ELSE
      CALL MPI_REDUCE(memory       , 0      , 3 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_LEADERS_SHARED , IERROR)
    END IF ! myLeaderGroupRank.EQ.0
  END IF ! myComputeNodeRank.EQ.0
END IF ! nProcessors.GT.1
#endif /*USE_MPI*/

! Find cache usage
IF (MPIRoot) THEN
  OPEN(NEWUNIT= ioUnit         , &
       FILE   = '/proc/meminfo', &
       STATUS = 'OLD'          , &
       ACTION = 'READ'         , &
       ACCESS = 'SEQUENTIAL'   , &
       IOSTAT = stat)
  IF (stat.NE.0) CALL Abort(__STAMP__,"Could not open /proc/meminfo.")

  ! Parallel IO: CN ROOT reads file and sends it to all other procs
  stat   = 0
  DO
    READ (ioUnit,"(A)",IOSTAT=stat) memName
    ! Split once to get the memory array name
    CALL split_string(memName      ,':',memStrings,memCount)
    ! Split again to get the memory array size
    CALL split_string(TRIM(ADJUSTL(memStrings(2))),' ',memArray,memCount)
    ! Read the memSize into an integer
    READ(memArray(1),*,IOSTAT=memstat) memSize
    IF (memstat.NE.0) EXIT
    IF (STRICMP(memStrings(1),'Buffers')) buffers = memSize
    IF (STRICMP(memStrings(1),'Cached' )) cached  = memSize
    IF(IS_IOSTAT_END(stat)) EXIT
  END DO
  CLOSE(ioUnit)
END IF ! MPIRoot

! Block, so no new caching/buffering occurs
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

IF (MPIRoot) THEN
  ! Check if the file exists
  IF (.NOT.MemoryStatFileExists) THEN
    INQUIRE(FILE='MemoryStatistics.csv',EXIST=MemoryStatFileExists)
    IF (.NOT.MemoryStatFileExists) THEN
      ! File really does not exist
      OPEN(NEWUNIT  = ioUnit                 , &
           FILE     = 'MemoryStatistics.csv' , &
           FORM     = 'FORMATTED'            , &
           STATUS   = 'NEW'                  , &
           RECL     = 50000                  , &
           IOSTAT   = stat)
      ! Write header
#if USE_LOADBALANCE
      WRITE(ioUnit,'(A,7(A1,A))') 'Walltime',',','Function',',','memUsed',',','memAvail',',','memTotal',',','Cached',',','Buffers',',','nLoadBalanceSteps'
#else
      WRITE(ioUnit,'(A,6(A1,A))') 'Walltime',',','Function',',','memUsed',',','memAvail',',','memTotal',',','Cached',',','Buffers'
#endif /*USE_LOADBALANCE*/
      CLOSE(ioUnit)
    END IF
  END IF

  OPEN(NEWUNIT  = ioUnit                 , &
       FILE     = 'MemoryStatistics.csv' , &
       FORM     = 'FORMATTED'            , &
       STATUS   = 'OLD'                  , &
       POSITION = 'APPEND'               , &
       RECL     = 50000                  , &
       IOSTAT   = stat)
#if USE_LOADBALANCE
  WRITE(ioUnit,'(E21.14E3,A1,A,6(A1,E21.14E3))') currentT-StartTime,',',TRIM(caller),',',memory(1),',',memory(2),',',memory(3),',',cached,',',buffers,',',REAL(nLoadBalanceSteps)
#else
  WRITE(ioUnit,'(E21.14E3,A1,A,5(A1,E21.14E3))') currentT-StartTime,',',TRIM(caller),',',memory(1),',',memory(2),',',memory(3),',',cached,',',buffers
#endif /*USE_LOADBALANCE*/
  CLOSE(ioUnit)
END IF

! Block, so writing is finished
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

END SUBROUTINE GetMemUsage


FUNCTION VerifyMemUsage(ArraySize)
!==================================================================================================================================
!> Verifies sufficient memory is available to allocate
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=8)                           :: ArraySize
LOGICAL                                   :: VerifyMemUsage
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                      :: memory(3)
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
#if !USE_MPI
INTEGER(KIND=8),PARAMETER                 :: nComputeNodeProcessors = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nProc => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (ArraySize*nProc .LT. memory(2)*kByte) THEN
  VerifyMemUsage = .TRUE.
ELSE
  VerifyMemUsage = .FALSE.
END IF

END ASSOCIATE

END FUNCTION VerifyMemUsage


!==================================================================================================================================
!> Allocate data after checking for available memory
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
SUBROUTINE Allocate_Safe_Logical_1(Array,nVal,STAT)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(OUT)           :: Array(:)                 !> Array to be allocated
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Number of variables in each rank
INTEGER,OPTIONAL,INTENT(OUT)              :: STAT                     !> Allocation status
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ALLOCSTAT                !> Allocation status
REAL                                      :: memory(3)
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
#if !USE_MPI
INTEGER(KIND=8),PARAMETER                 :: nComputeNodeProcessors = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

! Check if array is already allocated
IF (ALLOCATED(ARRAY)) CALL Abort(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL Abort(__STAMP__,'Trying to allocate array larger than available memory!')
END IF

END ASSOCIATE

END SUBROUTINE Allocate_Safe_Logical_1


!==================================================================================================================================
!> Allocate data after checking for available memory
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
SUBROUTINE Allocate_Safe_Int_1(Array,nVal,STAT)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,ALLOCATABLE,INTENT(OUT)           :: Array(:)                 !> Array to be allocated
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Number of variables in each rank
INTEGER,OPTIONAL,INTENT(OUT)              :: STAT                     !> Allocation status
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ALLOCSTAT                !> Allocation status
REAL                                      :: memory(3)
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
#if !USE_MPI
INTEGER(KIND=8),PARAMETER                 :: nComputeNodeProcessors = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

! Check if array is already allocated
IF (ALLOCATED(ARRAY)) CALL Abort(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL Abort(__STAMP__,'Trying to allocate array larger than available memory!')
END IF

END ASSOCIATE

END SUBROUTINE Allocate_Safe_Int_1


!==================================================================================================================================
!> Allocate data after checking for available memory
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
SUBROUTINE Allocate_Safe_Int_2(Array,nVal,STAT)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,ALLOCATABLE,INTENT(OUT)           :: Array(:,:)               !> Array to be allocated
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Number of variables in each rank
INTEGER,OPTIONAL,INTENT(OUT)              :: STAT                     !> Allocation status
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ALLOCSTAT                !> Allocation status
REAL                                      :: memory(3)
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
#if !USE_MPI
INTEGER(KIND=8),PARAMETER                 :: nComputeNodeProcessors = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

! Check if array is already allocated
IF (ALLOCATED(ARRAY)) CALL Abort(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1),nVal(2)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL Abort(__STAMP__,'Trying to allocate array larger than available memory!')
END IF

END ASSOCIATE

END SUBROUTINE Allocate_Safe_Int_2


!==================================================================================================================================
!> Allocate data after checking for available memory
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
SUBROUTINE Allocate_Safe_Real_1(Array,nVal,STAT)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,ALLOCATABLE,INTENT(OUT)              :: Array(:)                 !> Array to be allocated
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Number of variables in each rank
INTEGER,OPTIONAL,INTENT(OUT)              :: STAT                     !> Allocation status
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ALLOCSTAT                !> Allocation status
REAL                                      :: memory(3)
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
#if !USE_MPI
INTEGER(KIND=8),PARAMETER                 :: nComputeNodeProcessors = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

! Check if array is already allocated
IF (ALLOCATED(ARRAY)) CALL Abort(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL Abort(__STAMP__,'Trying to allocate array larger than available memory!')
END IF

END ASSOCIATE

END SUBROUTINE Allocate_Safe_Real_1


!==================================================================================================================================
!> Allocate data after checking for available memory
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
SUBROUTINE Allocate_Safe_Real_2(Array,nVal,STAT)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,ALLOCATABLE,INTENT(OUT)              :: Array(:,:)               !> Array to be allocated
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Number of variables in each rank
INTEGER,OPTIONAL,INTENT(OUT)              :: STAT                     !> Allocation status
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ALLOCSTAT                !> Allocation status
REAL                                      :: memory(3)
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
#if !USE_MPI
INTEGER(KIND=8),PARAMETER                 :: nComputeNodeProcessors = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

! Check if array is already allocated
IF (ALLOCATED(ARRAY)) CALL Abort(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1),nVal(2)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL Abort(__STAMP__,'Trying to allocate array larger than available memory!')
END IF

END ASSOCIATE

END SUBROUTINE Allocate_Safe_Real_2


!==================================================================================================================================
!> Allocate data after checking for available memory
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
SUBROUTINE Allocate_Safe_Real_3(Array,nVal,STAT)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,ALLOCATABLE,INTENT(OUT)              :: Array(:,:,:)             !> Array to be allocated
INTEGER,INTENT(IN)                        :: nVal(3)                  !> Number of variables in each rank
INTEGER,OPTIONAL,INTENT(OUT)              :: STAT                     !> Allocation status
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: ALLOCSTAT                !> Allocation status
REAL                                      :: memory(3)
INTEGER(KIND=8),PARAMETER                 :: kByte = 1024
#if !USE_MPI
INTEGER(KIND=8),PARAMETER                 :: nComputeNodeProcessors = 1
#endif /*!USE_MPI*/
!==================================================================================================================================

! Check if array is already allocated
IF (ALLOCATED(ARRAY)) CALL Abort(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1),nVal(2),nVal(3)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL Abort(__STAMP__,'Trying to allocate array larger than available memory!')
END IF

END ASSOCIATE

END SUBROUTINE Allocate_Safe_Real_3

END MODULE MOD_Memory
