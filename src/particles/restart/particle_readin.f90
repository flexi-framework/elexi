!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "flexi.h"
#include "particle.h"

MODULE MOD_Particle_Readin
!===================================================================================================================================
! Module to handle particle arrays upon loadbalance / restart
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ParticleReadin
  MODULE PROCEDURE ParticleReadin
END INTERFACE

PUBLIC :: ParticleReadin
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleReadin(doFlushFiles)
!===================================================================================================================================
! Distribute or readin particle data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Globals
USE MOD_Particle_Restart_Vars
USE MOD_ReadInTools            ,ONLY: PrintOption
USE MOD_StringTools            ,ONLY: STRICMP
! Boundary
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,ImpactDataSize
! HDF5
USE MOD_IO_HDF5
USE MOD_HDF5_Input             ,ONLY: OpenDataFile,CloseDataFile,ReadArray,ReadAttribute,GetDataSize
USE MOD_HDF5_Input             ,ONLY: File_ID,DatasetExists,nDims,HSize
USE MOD_HDF5_Output            ,ONLY: FlushFiles
! Memory
USE MOD_Particle_Memory        ,ONLY: Allocate_Safe
! Mesh
USE MOD_Mesh_Vars              ,ONLY: OffsetElem,nGlobalElems
! Particles
USE MOD_Particle_Boundary_Vars ,ONLY: doParticleReflectionTrack,doParticleImpactTrack
USE MOD_Particle_Vars          ,ONLY: PartInt,PartData,TurbPartData
USE MOD_Particle_Vars          ,ONLY: PartDataSize,TurbPartDataSize
! Restart
USE MOD_Restart_Vars           ,ONLY: RestartFile,RestartTime
! LoadBalance
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: nElemsOld,offsetElemOld,ElemInfoRank_Shared
USE MOD_LoadBalance_Vars       ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_LoadBalance_Vars       ,ONLY: MPInPartSend,MPInPartRecv,MPIoffsetPartSend,MPIoffsetPartRecv
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
USE MOD_Particle_Output_Vars   ,ONLY: PartDataVarSpecies
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL        :: doFlushFiles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Parameters
INTEGER,PARAMETER                  :: ELEM_FirstPartInd = 1
INTEGER,PARAMETER                  :: ELEM_LastPartInd  = 2
INTEGER,PARAMETER                  :: PartIntSize       = 2           ! number of entries in each line of PartInt
! Counters
INTEGER                            :: locnPart,offsetnPart
INTEGER                            :: iStr
INTEGER                            :: FirstElemInd,LastelemInd
! HDF5
INTEGER                            :: PartDim                         ! dummy for rank of partData
LOGICAL                            :: doFlushFiles_loc
LOGICAL                            :: EmissionTimeExists
! Impacts
LOGICAL                            :: ImpactDataExists
INTEGER                            :: ALLOCSTAT,j
INTEGER(KIND=IK)                   :: offsetImpact
INTEGER(KIND=IK)                   :: PartStateBoundaryVecLengthGlob
#if USE_MPI
INTEGER(KIND=IK)                   :: offsetImpactsProcCount,offsetImpacts(0:nProcessors)
#endif /* USE_MPI */
! LoadBalance
#if USE_LOADBALANCE
INTEGER                            :: iElem,PartRank
INTEGER                            :: offsetPartSend,offsetPartRecv
! Temporary arrays
INTEGER,ALLOCATABLE                :: PartIntTmp(:,:)
REAL,ALLOCATABLE                   :: PartDataTmp(:,:)
REAL,ALLOCATABLE                   :: TurbPartDataTmp(:,:)
#endif /*USE_LOADBALANCE*/
! MPI
#if USE_MPI
INTEGER                            :: iProc
#endif /*USE_MPI*/
! Particle turbulence models
INTEGER                            :: TurbPartHDF5Size                ! number of turbulent properties in HDF5 file
CHARACTER(LEN=255),ALLOCATABLE     :: StrVarNames(:)
! Timers
REAL                               :: StartT,EndT
!===================================================================================================================================

FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+PP_nElems

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  ! PartInt and PartData are still allocated from last WriteState
  StartT = MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' REDISTRIBUTING PARTICLES DURING LOADBALANCE...'

  ! ------------------------------------------------
  ! Check and set sizes
  ! ------------------------------------------------
  IF (PartDataSize.EQ.0) &
    CALL Abort(__STAMP__,'PartDataSize.EQ.0 but should have been set before loadbalance!')

  ! ------------------------------------------------
  ! PartInt
  ! ------------------------------------------------
  ALLOCATE(PartIntTmp(PartIntSize,FirstElemInd:LastElemInd))
  ASSOCIATE (&
          counts_send  => (PartIntSize*MPInElemSend     ) ,&
          disp_send    => (PartIntSize*MPIoffsetElemSend) ,&
          counts_recv  => (PartIntSize*MPInElemRecv     ) ,&
          disp_recv    => (PartIntSize*MPIoffsetElemRecv))
    ! Communicate PartInt over MPI
    CALL MPI_ALLTOALLV(PartInt,counts_send,disp_send,MPI_INTEGER_INT_KIND,PartIntTmp,counts_recv,disp_recv,MPI_INTEGER_INT_KIND,MPI_COMM_FLEXI,iError)
  END ASSOCIATE

  ! Calculate the PartInt deltas
  MPInPartSend      = 0
  MPIoffsetPartSend = 0
  ! Calculate the particles to send
  ! Loop with the old element over the new particle distribution
  DO iElem = 1,nElemsOld
    PartRank               = ElemInfo_Shared(ELEM_RANK,offsetElemOld+iElem)+1
    MPInPartSend(PartRank) = MPInPartSend(PartRank) + PartInt(2,offsetElemOld+iElem) - PartInt(1,offsetElemOld+iElem)
  END DO

  offsetPartSend = 0
  DO iProc = 2,nProcessors
    MPIoffsetPartSend(iProc) = SUM(MPInPartSend(1:iProc-1))
  END DO

  ! Calculate the elements to send
  MPInPartRecv      = 0
  MPIoffsetPartRecv = 0
  ! Loop with the new element over the old particle distribution
  DO iElem = 1,nElems
    PartRank               = ElemInfoRank_Shared(offsetElem+iElem)+1
    MPInPartRecv(PartRank) = MPInPartRecv(PartRank) + PartIntTmp(2,offsetElem+iElem) - PartIntTmp(1,offsetElem+iElem)
  END DO

  offsetPartRecv = 0
  DO iProc = 2,nProcessors
    MPIoffsetPartRecv(iProc) = SUM(MPInPartRecv(1:iProc-1))
  END DO
  CALL MOVE_ALLOC(PartIntTmp,PartInt)
  PartIntExists = .TRUE.

  ! ------------------------------------------------
  ! PartData
  ! ------------------------------------------------
  locnPart    = PartInt(ELEM_LastPartInd,LastElemInd)-PartInt(ELEM_FirstPartInd,FirstElemInd)
  offsetnPart = PartInt(ELEM_FirstPartInd,FirstElemInd)
  ALLOCATE(PartDataTmp(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
  ASSOCIATE (&
          counts_send  => (PartDataSize*MPInPartSend     ) ,&
          disp_send    => (PartDataSize*MPIoffsetPartSend) ,&
          counts_recv  => (PartDataSize*MPInPartRecv     ) ,&
          disp_recv    => (PartDataSize*MPIoffsetPartRecv))
    ! Communicate PartInt over MPI
    CALL MPI_ALLTOALLV(PartData,counts_send,disp_send,MPI_DOUBLE_PRECISION,PartDataTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
  END ASSOCIATE
  CALL MOVE_ALLOC(PartDataTmp,PartData)
  PartDataExists   = .TRUE.
  PP_nVarPartState = PartDataVarSpecies-1 !PP_nVarPart-1

  ! ------------------------------------------------
  ! TurbPartData
  ! ------------------------------------------------
  IF (ALLOCATED(TurbPartData)) THEN
    ALLOCATE(TurbPartDataTmp(TurbPartDataSize,offsetnPart+1:offsetnPart+locnPart))
    ASSOCIATE (&
            counts_send  => (TurbPartDataSize*MPInPartSend     ) ,&
            disp_send    => (TurbPartDataSize*MPIoffsetPartSend) ,&
            counts_recv  => (TurbPartDataSize*MPInPartRecv     ) ,&
            disp_recv    => (TurbPartDataSize*MPIoffsetPartRecv))
      ! Communicate PartInt over MPI
      CALL MPI_ALLTOALLV(TurbPartData,counts_send,disp_send,MPI_DOUBLE_PRECISION,TurbPartDataTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
    END ASSOCIATE
    CALL MOVE_ALLOC(TurbPartDataTmp,TurbPartData)
    TurbPartDataExists = .TRUE.
  END IF ! ALLOCATED(TurbPartData)

  ! ------------------------------------------------
  ! ImpactData
  ! ------------------------------------------------
  ! Impact data has no processor association, can remain on the previous processor
  GETTIME(EndT)
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' DONE! [',EndT-StartT,'s]'
  SWRITE(UNIT_stdOut,'(132("-"))')

! NOT. PerformLoadBalance
ELSE
#endif /*USE_LOADBALANCE*/
  doFlushFiles_loc = MERGE(doFlushFiles, .TRUE., PRESENT(doFlushFiles))

  IF (LEN_TRIM(RestartFile).EQ.0) THEN
    ! Delete all files since we are doing a fresh start --> moved from restart.f90 since we need it here
    IF (doFlushFiles_loc) CALL FlushFiles()
    RETURN
  END IF

  ! FIXME: Deallocate PartInt/PartData until loadbalance is always handled with MPI
  ! SDEALLOCATE(PartInt)
  ! SDEALLOCATE(PartData)

  GETTIME(StartT)
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' READING PARTICLES FROM RESTARTFILE...'
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

  ! Read the emission time
  CALL DatasetExists(File_ID,'EmissionTime',EmissionTimeExists,attrib=.TRUE.)
  ! Reset EmissionTime for ResetTime
  IF (EmissionTimeExists .AND. RestartTime.GT.0) THEN
    CALL ReadAttribute(File_ID,'EmissionTime',1,RealScalar=EmissionTime)
    SWRITE(UNIT_stdOut,'(A,F16.6)')' | Resuming particle emission from time    ',EmissionTime
  ELSE
    EmissionTime = 0.
  END IF

  ! ------------------------------------------------
  ! PartInt
  ! ------------------------------------------------
  CALL DatasetExists(File_ID,'PartInt',PartIntExists)
  IF(PartIntExists)THEN
    ALLOCATE(PartInt(PartIntSize,FirstElemInd:LastElemInd))

    ! Read number of local particles and their offset from HDF5
    CALL ReadArray('PartInt',2,(/PartIntSize,PP_nElems/),offsetElem,2,IntArray=PartInt)

    ! ------------------------------------------------
    ! PartData
    ! ------------------------------------------------
    locnPart    = PartInt(ELEM_LastPartInd ,LastElemInd )-PartInt(ELEM_FirstPartInd,FirstElemInd)
    offsetnPart = PartInt(ELEM_FirstPartInd,FirstElemInd)

    CALL DatasetExists(File_ID,'PartData',PartDataExists)
    IF(PartDataExists)THEN
      ! Get size of PartData
      CALL GetDataSize(File_ID,'PartData',PartDim,HSize)
      CHECKSAFEINT(HSize(2),4)
      PartDataSize = INT(HSize(1))
      CALL PrintOption('Number of particle variables','INFO',IntOpt=PartDataSize)

      ! For files, where no particle diameter was saved
      ALLOCATE(StrVarNames(PartDataSize))
      CALL ReadAttribute(File_ID,'VarNamesParticles',PartDataSize,StrArray=StrVarNames)
      ! Subtract variables that go into other arrays: PartSpecies
      PP_nVarPartState = PartDataSize-1
      DO iStr=1,PartDataSize
        ! Reflection
        IF (STRICMP('ReflectionCount',TRIM(StrVarNames(iStr)))) PP_nVarPartState = PP_nVarPartState - 1
        ! PartIndex
        IF (STRICMP('Index'          ,TRIM(StrVarNames(iStr)))) PP_nVarPartState = PP_nVarPartState - 1
      END DO
      DEALLOCATE(StrVarNames)

      ! Reflections are stored in the 8th data column. Do not start counting reflections mid-simulation
      IF (PartDataSize.EQ.PP_nVarPartState+1) THEN
        doParticleReflectionTrack = .FALSE.
        SWRITE(UNIT_stdOut,'(A3,A30,A3,I33)')' | ','Reflections not tracked previously. Disabling reflection counter',' | ',PartDataSize
      END IF

      ! Get PartData
      ALLOCATE(PartData(PartDataSize,offsetnPart+1:offsetnPart+locnPart))
      CALL ReadArray('PartData',2,(/PartDataSize,locnPart/),offsetnPart,2,RealArray=PartData)

      ! ------------------------------------------------
      ! TurbPartData
      ! ------------------------------------------------
      CALL DatasetExists(File_ID,'TurbPartData',TurbPartDataExists)
      IF(TurbPartDataExists)THEN
      ! Get size of TurbPartData, reuse PartDim dummy
        CALL GetDataSize(File_ID,'TurbPartData',PartDim,HSize)
        CHECKSAFEINT(HSize(2),4)
        TurbPartHDF5Size = INT(HSize(1))
        SWRITE(UNIT_stdOut,'(A3,A38,A3,I25)')' | ','Number of turbulent particle variables',' | ',TurbPartHDF5Size

        ! Compare number of turbulent properties of current run against HDF5
        IF ((TurbPartDataSize.EQ.0).AND.(TurbPartHDF5Size.NE.0)) THEN
          SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' | HDF5 state file containing SGS/RW data but current run in DNS mode. Ignoring...'
          TurbPartDataExists = .FALSE.
        ELSEIF ((TurbPartDataSize.NE.0).AND.(TurbPartDataSize.NE.TurbPartHDF5Size)) THEN
          CALL Abort(__STAMP__,' Number of turbulent variables in HDF5 does not match requested SGS/RW model!')
        ELSEIF ((TurbPartDataSize.NE.0).AND.(TurbPartDataSize.EQ.TurbPartHDF5Size)) THEN
          ! Maybe add check that compares the models here. Should resort to a warning in case we want to change the model midrun, so
          ! nothing urgent
          SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='YES') ' | HDF5 state file containing SGS/RW data with ',TurbPartDataSize, &
                                                      ' variables and matching current setup. Continuing run...'
          ! Get TurbPartData
          ALLOCATE(TurbPartData(TurbPartDataSize,offsetnPart+1:offsetnPart+locnPart))
          CALL ReadArray('TurbPartData',2,(/TurbPartDataSize,locnPart/),offsetnPart,2,RealArray=TurbPartData)
        END IF
      ELSE
        IF (TurbPartDataSize.NE.0) THEN
          SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' | HDF5 state file not containing SGS/RW data but model active in current run.'
          SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' | SGS/RW model will start without history...'
        END IF
      END IF ! TurbPartDataExists

      IF (myRank.EQ.nProcessors-1) WRITE(UNIT_stdOut,'(A,I0)')  ' | Particle(s) read from restart file: ',PartInt(ELEM_LastPartInd,nGlobalElems)
    END IF ! PartDataExists
  END IF ! PartIntExists

  ! ------------------------------------------------
  ! ImpactData
  ! ------------------------------------------------
  IF (doParticleImpactTrack) THEN
    PartStateBoundaryVecLength = 0
    CALL DatasetExists(File_ID,'ImpactData',ImpactDataExists)

    IF (ImpactDataExists) THEN
      CALL GetDataSize(File_ID,'ImpactData',PartDim,HSize)
      CHECKSAFEINT(HSize(2),4)
      PartStateBoundaryVecLengthGlob    = INT(HSize(2))

#if USE_MPI
      ! Distribute impacts between procs
      offsetImpacts              = 0
      PartStateBoundaryVecLength = PartStateBoundaryVecLengthGlob/nProcessors
      offsetImpactsProcCount     = PartStateBoundaryVecLengthGlob-PartStateBoundaryVecLength*nProcessors
      DO iProc = 0,nProcessors-1
        offsetImpacts(iProc) = PartStateBoundaryVecLength*iProc+MIN(iProc,offsetImpactsProcCount)
      END DO
      offsetImpacts(nProcessors) = PartStateBoundaryVecLengthGlob

      ! local impacts and offset
      PartStateBoundaryVecLength = offsetImpacts(myRank+1)-offsetImpacts(myRank)
      offsetImpact               = offsetImpacts(myRank)
#else
      PartStateBoundaryVecLength = PartStateBoundaryVecLength
      offsetImpact               = 0
#endif /* USE_MPI */

      ! Check if PartStateBoundary is sufficiently large
      ASSOCIATE( iMax => PartStateBoundaryVecLength )

        ! Check if array maximum is reached.
        ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
        IF (iMax.GT.10) THEN
          j = 1
          DO WHILE (iMax.GT.j*10)
            j = j*2
          END DO

          ! Check if the new PartStateBoundary size can be safely allocated
          CALL Allocate_Safe(PartStateBoundary,(/ImpactDataSize,j*10/), STAT=ALLOCSTAT)
          IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'ERROR in particle_boundary_init.f90: Cannot allocate PartStateBoundary array!')
          PartStateBoundary(:,:) = 0.
        END IF

        ! We lost the impact <-> proc association, so read in according to calculated distribution
        CALL ReadArray(ArrayName  = 'ImpactData'                                  ,&
                       rank       = 2                                             ,&
                       nVal       = (/ImpactDataSize,PartStateBoundaryVecLength/) ,&
                       offset_in  = offsetImpact                                  ,&
                       offset_dim = 2                                             ,&
                       RealArray  = PartStateBoundary(1:ImpactDataSize,1:PartStateBoundaryVecLength))
      END ASSOCIATE

      SWRITE(UNIT_stdOut,'(A,I0)')  ' | Impact(s)   read from restart file: ',PartStateBoundaryVecLengthGlob
    END IF ! ImpactDataExists
  END IF ! doParticleImpactTrack

  ! Keep everything in sync
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_FLEXI,iERROR)
#endif /* USE_MPI */

  CALL CloseDataFile()
  GETTIME(EndT)
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' READING PARTICLES FROM RESTARTFILE DONE! [',EndT-StartT,'s]'
  SWRITE(UNIT_stdOut,'(132("-"))')

  ! Delete all files that will be rewritten --> moved from restart.f90 since we need it here
  IF (doFlushFiles_loc) CALL FlushFiles(RestartTime)
#if USE_LOADBALANCE
END IF ! PerformLoadBalance
#endif /*USE_LOADBALANCE*/

END SUBROUTINE ParticleReadin

END MODULE MOD_Particle_Readin
