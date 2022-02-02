!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
MODULE MOD_Particle_Memory
!===================================================================================================================================
! MODULES
USE MOD_Globals
#if USE_MPI
USE MPI
#endif /*MPI*/
IMPLICIT NONE
PRIVATE

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

INTERFACE
  SUBROUTINE processmemusage(memUsed,memAvail,memTotal) BIND(C, name='processmemusage')
    USE ISO_C_BINDING,   ONLY : c_double
    real(c_double) :: memUsed
    real(c_double) :: memAvail
    real(c_double) :: memTotal
  END SUBROUTINE processmemusage
END INTERFACE

PUBLIC :: VerifyMemUsage
PUBLIC :: Allocate_Safe
PUBLIC :: processmemusage
!===================================================================================================================================

CONTAINS

FUNCTION VerifyMemUsage(ArraySize)
!==================================================================================================================================
!> Verifies sufficient memory is available to allocate
!> CAVE: Currently assumes each rank is calling with the same size
!==================================================================================================================================
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
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
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal in kB

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
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
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
IF (ALLOCATED(ARRAY)) CALL ABORT(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL ABORT(__STAMP__,'Trying to allocate array larger than available memory!')
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
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
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
IF (ALLOCATED(ARRAY)) CALL ABORT(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL ABORT(__STAMP__,'Trying to allocate array larger than available memory!')
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
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
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
IF (ALLOCATED(ARRAY)) CALL ABORT(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1),nVal(2)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL ABORT(__STAMP__,'Trying to allocate array larger than available memory!')
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
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
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
IF (ALLOCATED(ARRAY)) CALL ABORT(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL ABORT(__STAMP__,'Trying to allocate array larger than available memory!')
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
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
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
IF (ALLOCATED(ARRAY)) CALL ABORT(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1),nVal(2)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL ABORT(__STAMP__,'Trying to allocate array larger than available memory!')
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
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors
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
IF (ALLOCATED(ARRAY)) CALL ABORT(__STAMP__,'Trying to allocate already allocated array!')

! Find memory usage and requirements
CALL ProcessMemUsage(memory(1),memory(2),memory(3)) ! memUsed,memAvail,memTotal in kB

ASSOCIATE(nVal     => INT(nVal       ,KIND=8),          &
          VarSize  => INT(KIND(Array),KIND=8),          &
          nProc    => INT(nComputeNodeProcessors,KIND=8))

! Compare requested size against available memory
IF (PRODUCT(nVal)*VarSize*nProc .LT. memory(2)*kByte) THEN
  ALLOCATE(Array(nVal(1),nVal(2),nVal(3)),STAT=ALLOCSTAT)
  IF (PRESENT(STAT)) STAT=ALLOCSTAT
ELSE
  CALL ABORT(__STAMP__,'Trying to allocate array larger than available memory!')
END IF

END ASSOCIATE

END SUBROUTINE Allocate_Safe_Real_3

END MODULE MOD_Particle_Memory
