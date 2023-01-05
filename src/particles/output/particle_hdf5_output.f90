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
#include "particle.h"

!==================================================================================================================================
!> Module for the Particle HDF5 Output
!==================================================================================================================================
MODULE MOD_Particle_HDF5_Output
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
USE MOD_Particle_Globals
USE MOD_Mesh_Vars               ,ONLY: nGlobalElems,offsetElem
USE MOD_Part_Tools              ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze_Vars   ,ONLY: doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Boundary_Vars  ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Output_Vars
USE MOD_Particle_Restart_Vars   ,ONLY: EmissionTime
USE MOD_Particle_Vars           ,ONLY: nSpecies
USE MOD_Particle_Vars           ,ONLY: PartInt,PartData,TurbPartData
USE MOD_Particle_Vars           ,ONLY: PartDataSize,TurbPartDataSize
USE MOD_Particle_Vars           ,ONLY: doPartIndex,doWritePartDiam
#if USE_MPI
USE MOD_Particle_MPI_Vars       ,ONLY: PartMPI
#endif /*MPI*/
! Load balance
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: PartIntSize=2
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
LOGICAL                        :: reSwitch
! Particle turbulence models
!===================================================================================================================================

! Return if running particle code without any species
IF (nSpecies.LE.0) RETURN

! Allocate PartInt varnames array and fill it
IF (MPIRoot) THEN
  ALLOCATE(StrVarNames(PartIntSize))
  StrVarNames(1) = 'FirstPartID'
  StrVarNames(2) = 'LastPartID'

  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'VarNamesPartInt',PartIntSize,StrArray=StrVarNames)
  ! Write the beginning of the emission for restart from fluid solution
  IF (EmissionTime .NE. 0) &
    CALL WriteAttribute(File_ID,'EmissionTime',1,RealScalar=EmissionTime)
  CALL CloseDataFile()

  DEALLOCATE(StrVarNames)
END IF

reSwitch=.FALSE.
IF (gatheredWrite) THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch      = .TRUE.
  gatheredWrite = .FALSE.
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems)                              ,&
      nVar            => INT(PartIntSize)                               ,&
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

  IF (MPIRoot) THEN
    ! Allocate PartData varnames array and fill it
    ALLOCATE(StrVarNames(PartDataSize))
    StrVarNames(1:3) = (/'ParticlePositionX','ParticlePositionY','ParticlePositionZ'/)
    StrVarNames(4:6) = (/'VelocityX'        ,'VelocityY'        ,'VelocityZ'        /)
#if PP_nVarPartRHS == 6
    StrVarNames(7:9) = (/'AngularVelX'      ,'AngularVelY'      ,'AngularVelZ'      /)
#endif
    IF (doWritePartDiam) StrVarNames(PP_nVarPart) = 'PartDiam'
    StrVarNames(PartDataVarSpecies) = 'Species'
    IF (doPartIndex) StrVarNames(PartDataVarSpecies+1) = 'Index'
    IF (doParticleReflectionTrack) &
      StrVarNames(PartDataVarStart) = 'ReflectionCount'
    IF (doParticleDispersionTrack.OR.doParticlePathTrack) &
      StrVarNames(PartDataVarStart+PartDataVarShift:PartDataVarStart+2+PartDataVarShift)=(/'PartPathX','PartPathY','PartPathZ'/)

    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'VarNamesParticles',PartDataSize,StrArray=StrVarNames)
    CALL CloseDataFile()

    DEALLOCATE(StrVarNames)
  END IF

  ! Zero particles present in the complete domain.
  ! > Root writes empty dummy container to .h5 file (required for subsequent file access in ParaView)
  IF (nGlobalNbrOfParticles(3).EQ.0 .AND. MPIRoot) THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArray(      DataSetName  = 'PartData'     , rank = 2                          ,&
                            nValGlobal   = (/PartDataSize , nGlobalNbrOfParticles(3)/)        ,&
                            nVal         = (/PartDataSize , locnPart                /)        ,&
                            offset       = (/0            , offsetnPart             /)        ,&
                            collective   = .FALSE.        , RealArray = PartData)
      CALL CloseDataFile()
  END IF ! nGlobalNbrOfParticles(3).EQ.0 .AND. MPIRoot
#if USE_MPI
 CALL DistributedWriteArray(FileName                                                         ,&
                            DataSetName  = 'PartData'     , rank = 2                         ,&
                            nValGlobal   = (/PartDataSize , nGlobalNbrOfParticles(3)/)       ,&
                            nVal         = (/PartDataSize , locnPart                /)       ,&
                            offset       = (/0            , offsetnPart             /)       ,&
                            collective   = UseCollectiveIO, offSetDim = 2                    ,&
                            communicator = PartMPI%COMM   , RealArray = PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArray(          DataSetName  = 'PartData'     , rank = 2                         ,&
                            nValGlobal   = (/PartDataSize , nGlobalNbrOfParticles(3)/)       ,&
                            nVal         = (/PartDataSize , locnPart                /)       ,&
                            offset       = (/0            , offsetnPart             /)       ,&
                            collective   = .FALSE.        , RealArray = PartData)
  CALL CloseDataFile()
#endif /*MPI*/

  ! Turbulent particle properties currently not supported to be read directly. Do not associate varnames
#if USE_MPI
  IF (ALLOCATED(TurbPartData)) &
    CALL DistributedWriteArray(FileName                                                      ,&
                               DataSetName  = 'TurbPartData'     , rank = 2                  ,&
                               nValGlobal   = (/TurbPartDataSize , nGlobalNbrOfParticles(3)/),&
                               nVal         = (/TurbPartDataSize , locnPart                /),&
                               offset       = (/0                , offsetnPart             /),&
                               collective   = UseCollectiveIO    , offSetDim = 2             ,&
                               communicator = PartMPI%COMM       , RealArray = TurbPartData)
#else
  IF (ALLOCATED(TurbPartData)) THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArray(           DataSetName  = 'TurbPartData'     , rank = 2                  ,&
                               nValGlobal   = (/TurbPartDataSize , nGlobalNbrOfParticles(3)/),&
                               nVal         = (/TurbPartDataSize , locnPart                /),&
                               offset       = (/0                , offsetnPart             /),&
                               collective   = .FALSE.            , RealArray = TurbPartData)
    CALL CloseDataFile()
  END IF
#endif /*MPI*/

END ASSOCIATE
  ! reswitch
IF(reSwitch) gatheredWrite=.TRUE.

#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
DEALLOCATE( PartInt      &
          , PartData)
SDEALLOCATE(TurbPartData)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

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


#if USE_LOADBALANCE
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

! ElemTime might not be allocated, e.g, when running in posti mode
IF (.NOT.ALLOCATED(ElemTime)) RETURN

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

END MODULE MOD_Particle_HDF5_Output
