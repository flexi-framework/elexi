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
USE MOD_HDF5_Output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! Do not define an interface for this, cast information directly into subroutine
!INTERFACE DistributedWriteArray
!  MODULE PROCEDURE DistributedWriteArray
!END INTERFACE

INTERFACE WriteParticle
  MODULE PROCEDURE WriteParticle
END INTERFACE

#if USE_LOADBALANCE
!INTERFACE WriteElemDataToSeparateContainer
!  MODULE PROCEDURE WriteElemDataToSeparateContainer
!END INTERFACE

INTERFACE WriteElemTime
  MODULE PROCEDURE WriteElemTime
END INTERFACE
#endif /*USE_LOADBALANCE*/

PUBLIC :: WriteParticle
#if USE_MPI
PUBLIC :: DistributedWriteArray
#endif
#if USE_LOADBALANCE
!PUBLIC :: WriteElemDataToSeparateContainer
PUBLIC :: WriteElemTime
#endif
!==================================================================================================================================

CONTAINS


!===================================================================================================================================
! Subroutine that write the particle information to the state file
!> PartInt  contains the index of particles in each global element
!> PartData contains the indidividual properties of each particle
!===================================================================================================================================
SUBROUTINE WriteParticle(FileName)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,             ONLY: nGlobalElems, offsetElem
USE MOD_Part_Tools,            ONLY: UpdateNextFreePosition
USE MOD_Particle_Globals
USE MOD_Particle_Analyze_Vars, ONLY: PartPath,doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars,ONLY: doParticleReflectionTrack
! USE MOD_Particle_HDF5_Output
USE MOD_Particle_Restart_Vars, ONLY: EmissionTime
USE MOD_Particle_Vars,         ONLY: PDM,PEM,PartState,PartSpecies,PartReflCount,PartIndex
USE MOD_Particle_Vars,         ONLY: useLinkedList,doPartIndex,doWritePartDiam
#if USE_MPI
USE MOD_Particle_MPI_Vars,     ONLY: PartMPI
#endif /*MPI*/
! Particle turbulence models
USE MOD_Particle_Vars,         ONLY: TurbPartState
USE MOD_Particle_SGS_Vars,     ONLY: nSGSVars
#if USE_RW
USE MOD_Particle_RandomWalk_Vars,ONLY: nRWVars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVar,VarShift
LOGICAL                        :: reSwitch
INTEGER                        :: pcount,nInvalidPart
INTEGER                        :: locnPart,offsetnPart
INTEGER                        :: iPart,nGlobalNbrOfParticles,iElem
INTEGER,ALLOCATABLE            :: PartInt(:,:)
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER,PARAMETER              :: PartIntSize=2      !number of entries in each line of PartInt
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER                        :: tmpIndex,tmpIndex2,PP_nVarPart_loc
! Particle turbulence models
INTEGER                        :: TurbPartDataSize
REAL,ALLOCATABLE               :: TurbPartData(:,:)
!===================================================================================================================================

! Size and location of particle data
PP_nVarPart_loc = PP_nVarPart-1
PartDataSize    = PP_nVarPart_loc + 1
tmpIndex2       = PartDataSize
IF (doWritePartDiam) THEN
  PP_nVarPart_loc = PP_nVarPart
  PartDataSize    = PartDataSize + 1
  tmpIndex2       = PartDataSize
END IF
tmpIndex     = PartDataSize + 1
! Increase size if index is tracked
IF (doPartIndex) THEN
  PartDataSize = PartDataSize + 1
  tmpIndex     = tmpIndex     + 1
END IF
varShift     = 0
! Increase size if reflections are tracked
IF (doParticleReflectionTrack) THEN
  PartDataSize = PartDataSize + 1
  varShift     = 1
END IF
! Increase size if the absolute particle path is tracked
IF (doParticleDispersionTrack.OR.doParticlePathTrack) &
  PartDataSize = PartDataSize + 3

! Add turbulent dispersion data to output
IF (ALLOCATED(TurbPartState)) THEN
  TurbPartDataSize = nSGSVars
#if USE_RW
  TurbPartDataSize = TurbPartDataSize + nRWVars
#endif
END IF

nInvalidPart = 0
! Make sure to eliminate invalid particles as we cannot restart from NaN
DO pcount = 1,PDM%ParticleVecLength
  IF (.NOT. PDM%ParticleInside(pcount)) CYCLE
  IF (ANY(IEEE_IS_NAN(PartState(:,pcount)))) THEN
    nInvalidPart = nInvalidPart + 1
    PDM%ParticleInside(pcount) = .FALSE.
  END IF
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nInvalidPart,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(nInvalidPart,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif /*USE_MPI*/
IF (nInvalidPart.GT.0 .AND. MPIRoot) THEN
  WRITE(UNIT_stdOut,'(A,I0,A)') ' Detected ',nInvalidPart,' invalid particles during output. Removing ...'
END IF

! Determine number of particles in the complete domain
locnPart =   0
!>> Count number of particle on local proc
DO pcount = 1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(pcount)) THEN
    locnPart = locnPart + 1
  END IF
END DO

! Communicate the total number and offset
CALL GetOffsetAndGlobalNumberOfParts('WriteParticleToHDF5',offsetnPart,nGlobalNbrOfParticles,locnPart)

! Allocate data arrays for mean particle quantities
ALLOCATE(PartInt( PartIntSize ,offsetElem +1:offsetElem +PP_nElems))
ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
! Allocate data arrays for turbulent particle quantities
IF (ALLOCATED(TurbPartState)) ALLOCATE(TurbPartData(TurbPartDataSize,offsetnPart+1:offsetnPart+locnPart))

! Update next free position using a linked list
ALLOCATE(PEM%pStart (offsetElem+1:offsetElem+PP_nElems) , &
         PEM%pNumber(offsetElem+1:offsetElem+PP_nElems) , &
         PEM%pNext  (1           :PDM%maxParticleNumber), &
         PEM%pEnd   (offsetElem+1:offsetElem+PP_nElems))
useLinkedList = .TRUE.
CALL UpdateNextFreePosition()

! Walk along the linked list and fill the data arrays
iPart = offsetnPart
! Walk over all elements on local proc
DO iElem = offsetElem+1,offsetElem+PP_nElems
  ! Set start of particle numbers in current element
  PartInt(1,iElem) = iPart
  ! Find all particles in current element
  IF (ALLOCATED(PEM%pNumber)) THEN
    PartInt(2,iElem) = PartInt(1,iElem) + PEM%pNumber(iElem)
    ! Sum up particles and add properties to output array
    pcount = PEM%pStart(iElem)
    DO iPart = PartInt(1,iElem)+1,PartInt(2,iElem)
      PartData(1:tmpIndex2-1,iPart) = PartState(1:tmpIndex2-1,pcount)
      PartData(tmpIndex2,iPart)     = REAL(PartSpecies(pcount))
      IF (doPartIndex)                                      PartData(PP_nVarPart_loc+2                    ,iPart) = REAL(PartIndex(pcount))
      IF (doParticleReflectionTrack)                        PartData(tmpIndex                             ,iPart) = REAL(PartReflCount(pcount))
      IF (doParticleDispersionTrack.OR.doParticlePathTrack) PartData(tmpIndex+varShift:tmpIndex+2+varShift,iPart) = PartPath(1:3,pcount)

      ! Turbulent particle properties
      IF (ALLOCATED(TurbPartState))  TurbPartData(:,iPart) = TurbPartState(:,pcount)

      ! Set the index to the next particle
      pcount = PEM%pNext(pcount)
    END DO
    ! Set counter to the end of particle number in the current element
    iPart = PartInt(2,iElem)
  ELSE
    CALL Abort(__STAMP__, " Particle HDF5-Output method not supported! PEM%pNumber not associated")
  END IF
  PartInt(2,iElem)=iPart
END DO ! iElem = offsetElem+1,offsetElem+PP_nElems

! Allocate PartInt varnames array and fill it
nVar=2
ALLOCATE(StrVarNames(nVar))
StrVarNames(1)='FirstPartID'
StrVarNames(2)='LastPartID'

IF (MPIRoot) THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'VarNamesPartInt',nVar,StrArray=StrVarNames)
  ! Write the beginning of the emission for restart from fluid solution
  IF (EmissionTime .NE. 0) &
    CALL WriteAttribute(File_ID,'EmissionTime',1,RealScalar=EmissionTime)
  CALL CloseDataFile()
END IF

reSwitch=.FALSE.
IF(gatheredWrite)THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch=.TRUE.
  gatheredWrite=.FALSE.
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems)                              ,&
      nVar            => INT(nVar)                                      ,&
      PP_nElems       => INT(PP_nElems)                                 ,&
      offsetElem      => INT(offsetElem)                                ,&
      PartDataSize    => INT(PartDataSize))

  CALL GatheredWriteArray(FileName                                      ,&
                          create      = .FALSE.                         ,&
                          DataSetName = 'PartInt'                       ,&
                          rank        = 2                               ,&
                          nValGlobal  = (/nVar,nGlobalElems/)           ,&
                          nVal        = (/nVar,PP_nElems   /)           ,&
                          offset      = (/0   ,offsetElem  /)           ,&
                          collective  = .TRUE.                          ,&
                          IntArray    = PartInt)
  DEALLOCATE(StrVarNames)

  ! Allocate PartData varnames array and fill it
  ALLOCATE(StrVarNames(PartDataSize))
  StrVarNames(1:3) = (/'ParticlePositionX','ParticlePositionY','ParticlePositionZ'/)
  StrVarNames(4:6) = (/'VelocityX'        ,'VelocityY'        ,'VelocityZ'        /)
#if PP_nVarPartRHS == 6
  StrVarNames(7:9) = (/'AngularVelX'      ,'AngularVelY'      ,'AngularVelZ'      /)
#endif
  IF (doWritePartDiam) StrVarNames(PP_nVarPart) = 'PartDiam'
  StrVarNames(tmpIndex2) = 'Species'
  IF (doPartIndex) StrVarNames(tmpIndex2+1) = 'Index'
  IF (doParticleReflectionTrack) &
    StrVarNames(tmpIndex) = 'ReflectionCount'
  IF (doParticleDispersionTrack.OR.doParticlePathTrack) &
    StrVarNames(tmpIndex+varShift:tmpIndex+2+varShift)=(/'PartPathX','PartPathY','PartPathZ'/)

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'VarNamesParticles',PartDataSize,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

  ! Zero particles present in the complete domain.
  ! > Root writes empty dummy container to .h5 file (required for subsequent file access in ParaView)
  IF (nGlobalNbrOfParticles.EQ.0 .AND. MPIRoot) THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArray(      DataSetName = 'PartData'                                      ,&
                            rank        = 2                                               ,&
                            nValGlobal  = (/PartDataSize,nGlobalNbrOfParticles /)         ,&
                            nVal        = (/PartDataSize,locnPart   /)                    ,&
                            offset      = (/0           ,offsetnPart/)                    ,&
                            collective  = .FALSE.                                         ,&
                            RealArray   = PartData)
      CALL CloseDataFile()
  END IF ! nGlobalNbrOfParticles.EQ.0 .AND. MPIRoot
#if USE_MPI
 CALL DistributedWriteArray(FileName                                                      ,&
                            DataSetName  = 'PartData'                                     ,&
                            rank         = 2                                              ,&
                            nValGlobal   = (/PartDataSize,nGlobalNbrOfParticles /)        ,&
                            nVal         = (/PartDataSize,locnPart   /)                   ,&
                            offset       = (/0           ,offsetnPart/)                   ,&
                            collective   =.FALSE.                                         ,&
                            offSetDim    = 2                                              ,&
                            communicator = PartMPI%COMM                                   ,&
                            RealArray    = PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArray(          DataSetName  = 'PartData'                                     ,&
                            rank         = 2                                              ,&
                            nValGlobal   = (/PartDataSize,nGlobalNbrOfParticles /)        ,&
                            nVal         = (/PartDataSize,locnPart   /)                   ,&
                            offset       = (/0           ,offsetnPart/)                   ,&
                            collective   = .TRUE.                                         ,&
                            RealArray    = PartData)
  CALL CloseDataFile()
#endif /*MPI*/

  ! Turbulent particle properties currently not supported to be read directly. Do not associate varnames
#if USE_MPI
  IF (ALLOCATED(TurbPartState)) &
    CALL DistributedWriteArray(FileName                                                   ,&
                               DataSetName  = 'TurbPartData'                              ,&
                               rank         = 2                                           ,&
                               nValGlobal   = (/TurbPartDataSize,nGlobalNbrOfParticles /) ,&
                               nVal         = (/TurbPartDataSize,locnPart   /)            ,&
                               offset       = (/0               ,offsetnPart/)            ,&
                               collective   = .FALSE.                                     ,&
                               offSetDim    = 2                                           ,&
                               communicator = PartMPI%COMM                                ,&
                               RealArray    = TurbPartData)
#else
  IF (ALLOCATED(TurbPartState)) THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArray(           DataSetName  = 'TurbPartData'                              ,&
                               rank         = 2                                           ,&
                               nValGlobal   = (/TurbPartDataSize,nGlobalNbrOfParticles/)  ,&
                               nVal         = (/TurbPartDataSize,locnPart/)               ,&
                               offset       = (/0               ,offsetnPart/)            ,&
                               collective   = .TRUE.                                      ,&
                               RealArray    = TurbPartData)
    CALL CloseDataFile()
  END IF
#endif /*MPI*/

END ASSOCIATE
  ! reswitch
IF(reSwitch) gatheredWrite=.TRUE.

! De-allocate linked list and return to normal particle array mode
useLinkedList=.FALSE.
DEALLOCATE( StrVarNames  &
          , PartInt      &
          , PartData     &
          , PEM%pStart   &
          , PEM%pNumber  &
          , PEM%pNext    &
          , PEM%pEnd)

END SUBROUTINE WriteParticle


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


!===================================================================================================================================
!> Calculate the particle offset and global number of particles across all processors
!> In this routine the number are calculated using integer KIND=8, but are returned with integer KIND=ICC in order to test if using
!> integer KIND=8 is required for total number of particles, particle boundary state, lost particles or clones
!===================================================================================================================================
SUBROUTINE GetOffsetAndGlobalNumberOfParts(CallingRoutine,offsetnPart,globnPart,locnPart)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: CallingRoutine
INTEGER(KIND=IK),INTENT(IN)  :: locnPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=IK),INTENT(OUT) :: offsetnPart
INTEGER(KIND=IK),INTENT(OUT) :: globnPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)              :: globnPart8                         ! always integer KIND=8
#if USE_MPI
INTEGER(KIND=8)              :: locnPart8,locnPart8Recv            ! always integer KIND=8
#endif
!===================================================================================================================================
#if USE_MPI
locnPart8     = INT(locnPart,8)
locnPart8Recv = 0_IK
CALL MPI_EXSCAN(locnPart8,locnPart8Recv,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart   = INT(locnPart8Recv,KIND=IK)

! Last proc calculates the global number and broadcasts it
IF(myRank.EQ.nProcessors-1) locnPart8=locnPart8Recv+locnPart8
CALL MPI_BCAST(locnPart8,1,MPI_INTEGER8,nProcessors-1,MPI_COMM_WORLD,iError)

! Global numbers
globnPart8   = locnPart8
LOGWRITE(*,*) TRIM(CallingRoutine)//'offsetnPart,locnPart,globnPart8',offsetnPart,locnPart,globnPart8
#else
offsetnPart  = 0_IK
globnPart8   = INT(locnPart,8)
#endif

! Sanity check: Add up all particles with integer KIND=8 and compare
IF (MPIRoot) THEN
  ! Check if offsetnPart is kind=8 is the number of particles is larger than integer KIND=4
  IF (globnPart8.GT.INT(HUGE(offsetnPart),8)) THEN
    WRITE(UNIT_stdOut,'(A,I0)') '\n\n\nTotal number of particles  : ',globnPart8
    WRITE(UNIT_stdOut,'(A,I0)')       'Maximum number of particles: ',HUGE(offsetnPart)
    CALL Abort(__STAMP__,TRIM(CallingRoutine)//' has encountered more than integer KIND=4 particles!')
  END IF
END IF ! MPIRoot

! Cast to Kind=IK before returning the number
globnPart = INT(globnPart8,KIND=IK)

END SUBROUTINE GetOffsetAndGlobalNumberOfParts


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
END IF ! (MAXVAL(ElemData).LE.0.0).AND.DoRestart)

! Write 'ElemTime' container
CALL GatheredWriteArray(FileName                               ,&
                        create          = .FALSE.              ,&
                        DataSetName     = 'ElemTime'           ,&
                        rank            = 2                    ,&
                        nValGlobal      = (/1   ,nGlobalElems/),&
                        nVal            = (/1   ,nElems      /),&
                        offset          = (/0   ,offsetElem  /),&
                        collective      = .TRUE.               ,&
                        RealArray       = ElemTime)

END SUBROUTINE WriteElemTime
#endif /*USE_LOADBALANCE*/

END MODULE MOD_Particle_HDF5_output
