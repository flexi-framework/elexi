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
!> Contains the routines to exchange data using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_MPI_Shared
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE

#if USE_MPI_SHARED
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersMPIShared
  MODULE PROCEDURE DefineParametersMPIShared
END INTERFACE

INTERFACE InitMPIShared
  MODULE PROCEDURE InitMPIShared
END INTERFACE

INTERFACE FinalizeMPIShared
  MODULE PROCEDURE FinalizeMPIShared
END INTERFACE

INTERFACE Allocate_Shared
  MODULE PROCEDURE Allocate_Shared_Real_5
END INTERFACE

PUBLIC::DefineParametersMPIShared
PUBLIC::InitMPIShared
PUBLIC::FinalizeMPIShared
PUBLIC::Allocate_Shared
!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE DefineParametersMPIShared()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection         ("MPI Shared")

END SUBROUTINE DefineParametersMPIShared


!==================================================================================================================================
!> Initialize for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE InitMPIShared()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
integer :: ierr
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MPI Shared ...'

! Split the node communicator (shared memory) from the global communicator
CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, myRank, MPI_INFO_NULL, MPI_COMM_SHARED, ierr)

! Find my rank on the shared communicator, comm size and proc name
CALL MPI_COMM_RANK(MPI_COMM_SHARED, myRank_Shared, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_SHARED, nProcessors_Shared, ierr)
!CALL MPI_GET_PROCESSOR_NAME(myNode_Shared, myNode_Shared_Len, ierr)

MPISharedInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MPI Shared DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitMPIShared


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_5(Datasize_Byte,nVal,SM_WIN,DataPointer)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER    ,INTENT(IN)                    :: Datasize_Byte            !> Length of the data size in bytes
INTEGER    ,INTENT(IN)                    :: nVal(5)                  !> Local number of variables in each rank
INTEGER    ,INTENT(OUT)                   :: SM_WIN                   !> Shared memory window
REAL       ,INTENT(OUT),POINTER           :: DataPointer(:,:,:,:,:)   !> Pointer to the RMA window
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================

! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
IF (myRank_Shared.EQ.0) THEN
  WIN_SIZE  = datasize_byte
ELSE
  WIN_SIZE  = 0
END IF
DISP_UNIT = 1

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_NULL, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myRank_Shared.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_REAL_5

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
CALL MPI_COMM_FREE(MPI_COMM_SHARED, IERROR)
MPISharedInitIsDone=.FALSE.

END SUBROUTINE FinalizeMPIShared

#endif /* MPI_SHARED */
END MODULE
