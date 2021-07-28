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
!> Module for the Particle HDF5 Output
!==================================================================================================================================
MODULE MOD_Particle_HDF5_output
! MODULES
USE MOD_IO_HDF5
USE MOD_HDF5_WriteArray
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! Do not define an interface for this, cast information directly into subroutine
!INTERFACE DistributedWriteArray
!  MODULE PROCEDURE DistributedWriteArray
!END INTERFACE

#if USE_LOADBALANCE
!INTERFACE WriteElemDataToSeparateContainer
!  MODULE PROCEDURE WriteElemDataToSeparateContainer
!END INTERFACE

INTERFACE WriteElemTime
  MODULE PROCEDURE WriteElemTime
END INTERFACE
#endif /*USE_LOADBALANCE*/

#if USE_MPI
PUBLIC :: DistributedWriteArray
#endif
#if USE_LOADBALANCE
!PUBLIC :: WriteElemDataToSeparateContainer
PUBLIC :: WriteElemTime
#endif
!==================================================================================================================================

CONTAINS

#if USE_MPI
SUBROUTINE DistributedWriteArray(FileName,DataSetName,rank,nValGlobal,nVal,offset,collective,&
                                 offSetDim,communicator,RealArray,IntArray,StrArray)
!===================================================================================================================================
! Write distributed data to proc, e.g. particles which are not hosted by each proc
! a new output-communicator is build and afterwards killed
! offset is in the last dimension
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Output_Vars,       ONLY: WriteStateFiles
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName,DataSetName
LOGICAL,INTENT(IN)             :: collective
INTEGER,INTENT(IN)             :: offSetDim,communicator
INTEGER,INTENT(IN)             :: rank,nVal(rank),nValGlobal(rank),offset(rank)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray( PRODUCT(nVal))
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: color,OutPutCOMM,nOutPutProcs
LOGICAL                        :: DataOnProc,DoNotSplit,OutputCollective,OutputSingle
!===================================================================================================================================
IF (.NOT.WriteStateFiles) RETURN

! 1: check if every proc of given communicator has data
DataOnProc = MERGE(.TRUE.,.FALSE.,nVal(offSetDim).GT.0)
CALL MPI_ALLREDUCE(DataOnProc,DoNotSplit,1,MPI_LOGICAL,MPI_LAND,COMMUNICATOR,IERROR)

! 2: if any proc has no data, split the communicator and write only with the new communicator
IF (.NOT.DoNotSplit) THEN
  color = MERGE(87,MPI_UNDEFINED,DataOnProc)
  CALL MPI_COMM_SPLIT(COMMUNICATOR,color,MPI_INFO_NULL,OutputCOMM,iError)

  IF(DataOnProc) THEN
    CALL MPI_COMM_SIZE(OutputCOMM,nOutPutProcs,iError)

    OutputCollective = MERGE(.FALSE.,collective,nOutPutProcs.EQ.1)
    OutputSingle     = MERGE(.TRUE.,.FALSE.    ,nOutPutProcs.EQ.1)

    CALL OpenDataFile(FileName,create=.FALSE.,single=OutputSingle,readOnly=.FALSE.,communicatorOpt=OutputCOMM)
    IF(PRESENT(RealArray)) CALL WriteArray( DataSetName,rank,nValGlobal,nVal,                        &
                                            offset,collective=OutputCollective,RealArray=RealArray)
    IF(PRESENT(IntArray))  CALL WriteArray( DataSetName,rank,nValGlobal,nVal,                        &
                                            offset,collective=OutputCollective,IntArray =IntArray)
    IF(PRESENT(StrArray))  CALL WriteArray( DataSetName,rank,nValGlobal,nVal,                        &
                                            offset,collective=OutputCollective,StrArray =StrArray)
    CALL CloseDataFile()
    CALL MPI_BARRIER  (OutputCOMM,IERROR)
    CALL MPI_COMM_FREE(OutputCOMM,iERROR)
  END IF
  ! MPI Barrier is required so procs don't open the datafile while the output procs are still writing
  CALL MPI_BARRIER(COMMUNICATOR,IERROR)
  ! Communicator was realized, also nullify variable
  OutputCOMM = MPI_UNDEFINED

! 3: else write with all procs of the given communicator
! communicator_opt has to be the given communicator or else procs that are not in the given communicator might block the write out
! e.g. surface communicator contains only procs with physical surface and MPI_COMM_WORLD contains every proc
!      Consequently, MPI_COMM_WORLD would block communication
ELSE
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
  IF(PRESENT(RealArray)) CALL WriteArray( DataSetName,rank,nValGlobal,nVal,                          &
                                          offset,collective,RealArray=RealArray)
  IF(PRESENT(IntArray))  CALL WriteArray( DataSetName,rank,nValGlobal,nVal,                          &
                                          offset,collective,IntArray =IntArray)
  IF(PRESENT(StrArray))  CALL WriteArray( DataSetName,rank,nValGlobal,nVal,                          &
                                          offset,collective,StrArray =StrArray)
  CALL CloseDataFile()
END IF

END SUBROUTINE DistributedWriteArray
#endif /*USE_MPI*/


#if USE_LOADBALANCE
!SUBROUTINE WriteElemDataToSeparateContainer(FileName,ElemList,ElemDataName)
!!===================================================================================================================================
!!> Similar to WriteAdditionalElemData() but only writes one of the fields to a separate container
!!> ----------------
!!> Write additional data for analyze purpose to HDF5.
!!> The data is taken from a lists, containing either pointers to data arrays or pointers
!!> to functions to generate the data, along with the respective varnames.
!!>
!!> Two options are available:
!!>    1. WriteAdditionalElemData:
!!>       Element-wise scalar data, e.g. the timestep or indicators.
!!>       The data is collected in a single array and written out in one step.
!!>       DO NOT MISUSE NODAL DATA FOR THIS! IT WILL DRASTICALLY INCREASE FILE SIZE AND SLOW DOWN IO!
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_PreProc
!USE MOD_Mesh_Vars        ,ONLY: nElems
!USE MOD_HDF5_Input       ,ONLY: ReadArray
!USE MOD_LoadBalance_Vars ,ONLY: ElemTime,ElemTime_tmp
!USE MOD_Restart_Vars     ,ONLY: DoRestart
!USE MOD_Mesh_Vars        ,ONLY: nGlobalElems,offsetelem
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!CHARACTER(LEN=255),INTENT(IN)        :: FileName
!TYPE(tElementOut),POINTER,INTENT(IN) :: ElemList            !< Linked list of arrays to write to file
!CHARACTER(LEN=*),INTENT(IN)          :: ElemDataName
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!CHARACTER(LEN=255)                   :: StrVarNames
!REAL,ALLOCATABLE                     :: ElemData(:,:)
!INTEGER                              :: nVar
!TYPE(tElementOut),POINTER            :: e
!!===================================================================================================================================
!
!IF(.NOT. ASSOCIATED(ElemList)) RETURN
!
!! Allocate variable names and data array
!ALLOCATE(ElemData(1,nElems))
!
!! Fill the arrays
!nVar = 0
!e=>ElemList
!DO WHILE(ASSOCIATED(e))
!  StrVarNames=e%VarName
!  IF(StrVarNames.EQ.TRIM(ElemDataName))THEN
!    nVar = nVar+1
!    IF(ASSOCIATED(e%RealArray))    ElemData(nVar,:)=e%RealArray(1:nElems)
!    IF(ASSOCIATED(e%RealScalar))   ElemData(nVar,:)=e%RealScalar
!    IF(ASSOCIATED(e%IntArray))     ElemData(nVar,:)=REAL(e%IntArray(1:nElems))
!    IF(ASSOCIATED(e%IntScalar))    ElemData(nVar,:)=REAL(e%IntScalar)
!    IF(ASSOCIATED(e%eval))         CALL e%eval(ElemData(nVar,:)) ! function fills elemdata
!    EXIT
!  END IF ! StrVarNames.EQ.TRIM(ElemDataName)
!  e=>e%next
!END DO
!
!IF(nVar.NE.1) CALL ABORT(__STAMP__,'WriteElemDataToSeparateContainer: Array not found in ElemData = '//TRIM(ElemDataName))
!
!! Check if ElemTime is all zeros and if this is a restart (save the old values)
!IF((MAXVAL(ElemData).LE.0.0)          .AND.& ! Restart
!    DoRestart                         .AND.& ! Restart
!    (TRIM(ElemDataName).EQ.'ElemTime').AND.& ! only for ElemTime array
!    ALLOCATED(ElemTime_tmp)) THEN            ! only allocated when not starting simulation from zero
!  ! Additionally, store old values in ElemData container
!  ElemTime = ElemTime_tmp
!
!    ! Write 'ElemTime' container
!    CALL GatheredWriteArray(FileName                               ,&
!                            create          = .FALSE.              ,&
!                            DataSetName     = TRIM(ElemDataName)   ,&
!                            rank            = 2                    ,&
!                            nValGlobal      = (/nVar,nGlobalElems/),&
!                            nVal            = (/nVar,nElems      /),&
!                            offset          = (/0   ,offsetElem  /),&
!                            collective      = .TRUE.               ,&
!                            RealArray       = ElemTime_tmp)
!ELSE
!    CALL GatheredWriteArray(FileName                               ,&
!                            create          = .FALSE.              ,&
!                            DataSetName     = TRIM(ElemDataName)   ,&
!                            rank            = 2                    ,&
!                            nValGlobal      = (/nVar,nGlobalElems/),&
!                            nVal            = (/nVar,nElems      /),&
!                            offset          = (/0   ,offsetElem  /),&
!                            collective      = .TRUE.               ,&
!                            RealArray       = ElemData)
!END IF ! (MAXVAL(ElemData).LE.0.0).AND.DoRestart.AND.(TRIM(ElemDataName).EQ.'ElemTime')
!
!DEALLOCATE(ElemData)
!
!END SUBROUTINE WriteElemDataToSeparateContainer


SUBROUTINE WriteElemTime(FileName)
!===================================================================================================================================
!> Similar to WriteAdditionalElemData() but only writes one of the fields to a separate container
!> ----------------
!> Write additional data for analyze purpose to HDF5.
!> The data is taken from a lists, containing either pointers to data arrays or pointers
!> to functions to generate the data, along with the respective varnames.
!>
!> Two options are available:
!>    1. WriteAdditionalElemData:
!>       Element-wise scalar data, e.g. the timestep or indicators.
!>       The data is collected in a single array and written out in one step.
!>       DO NOT MISUSE NODAL DATA FOR THIS! IT WILL DRASTICALLY INCREASE FILE SIZE AND SLOW DOWN IO!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: nElems
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,ElemTime_tmp
USE MOD_Restart_Vars     ,ONLY: DoRestart
USE MOD_Mesh_Vars        ,ONLY: nGlobalElems,offsetelem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)        :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Check if ElemTime is all zeros and if this is a restart (save the old values)
IF((MAXVAL(ElemTime).LE.0.0)          .AND.& ! Restart
    DoRestart                         .AND.& ! Restart
    ALLOCATED(ElemTime_tmp)) THEN            ! only allocated when not starting simulation from zero
    ! Additionally, store old values in ElemData container
    ElemTime = ElemTime_tmp

    ! Write 'ElemTime' container
    CALL GatheredWriteArray(FileName                               ,&
                            create          = .FALSE.              ,&
                            DataSetName     = 'ElemTime'           ,&
                            rank            = 2                    ,&
                            nValGlobal      = (/1   ,nGlobalElems/),&
                            nVal            = (/1   ,nElems      /),&
                            offset          = (/0   ,offsetElem  /),&
                            collective      = .TRUE.               ,&
                            RealArray       = ElemTime_tmp)
ELSE
    CALL GatheredWriteArray(FileName                               ,&
                            create          = .FALSE.              ,&
                            DataSetName     = 'ElemTime'           ,&
                            rank            = 2                    ,&
                            nValGlobal      = (/1   ,nGlobalElems/),&
                            nVal            = (/1   ,nElems      /),&
                            offset          = (/0   ,offsetElem  /),&
                            collective      = .TRUE.               ,&
                            RealArray       = ElemTime)
END IF ! (MAXVAL(ElemData).LE.0.0).AND.DoRestart)

END SUBROUTINE WriteElemTime
#endif /*USE_LOADBALANCE*/

END MODULE MOD_Particle_HDF5_output
