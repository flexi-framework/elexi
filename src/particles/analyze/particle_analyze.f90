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

!===================================================================================================================================
!> Contains routines for statistical analysis of particle behavior
!===================================================================================================================================
MODULE MOD_Particle_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleAnalyze
  MODULE PROCEDURE InitParticleAnalyze
END INTERFACE

INTERFACE ParticleAnalyze
  MODULE PROCEDURE ParticleAnalyze
END INTERFACE

INTERFACE ParticleInformation
  MODULE PROCEDURE ParticleInformation
END INTERFACE

INTERFACE FinalizeParticleAnalyze
  MODULE PROCEDURE FinalizeParticleAnalyze
END INTERFACE

INTERFACE CalcKineticEnergy
  MODULE PROCEDURE CalcKineticEnergy
END INTERFACE

PUBLIC :: InitParticleAnalyze
PUBLIC :: ParticleAnalyze
PUBLIC :: ParticleInformation
PUBLIC :: FinalizeParticleAnalyze
PUBLIC :: DefineParametersParticleAnalyze
PUBLIC :: CalcKineticEnergy
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for analyze if particles (.csv output)
!==================================================================================================================================
SUBROUTINE DefineParametersParticleAnalyze()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
!USE MOD_AnalyzeEquation ,ONLY: DefineParametersAnalyzeEquation
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Analyze")

CALL prms%CreateIntOption(      'Part-AnalyzeStep'        , 'Analyze is performed each Nth time step','1')
CALL prms%CreateLogicalOption(  'CalcKineticEnergy'       , 'Calculate Kinetic Energy. ','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPartBalance'         , 'Calculate the Particle Kinetic Energy Balance'//&
                                                            '- input and outflow kinetic energy of all particles','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackDispersion'    , 'Track particle convergence radius (i.e. absolute path)','.FALSE.')
CALL prms%CreateLogicalOption(  'Part-TrackPath'          , 'Track particle path (i.e. relative path)','.FALSE.')
#if CODE_ANALYZE
CALL prms%CreateIntOption(      'PartOut'                 , 'If compiled with CODE_ANALYZE flag: For This particle number'   //&
                                                            ' every tracking information is written as STDOUT.'                &
                                                          , '0')
CALL prms%CreateIntOption(      'MPIRankOut'              , 'If compiled with CODE_ANALYZE flag: This MPI-Proc writes the'   //&
                                                            ' tracking information for the defined PartOut.'                   &
                                                          , '0')
#endif /*CODE_ANALYZE*/
CALL prms%CreateIntOption(      'Part-nRPs'               , 'Number of record planes'                                          &
                                                          , '0')
CALL prms%CreateIntOption(      'Part-RPMemory'           , 'Record particles memory'                                          &
                                                          , '100')
!CALL prms%CreateStringOption(  'Part-RecordType[$]'        , 'Type of record plane.\n'                                        //&
!                                                             ' - plane\n'                                                       &
!                                                           , 'none', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-RPDirection[$]'     , 'Direction of the normal vector of the record plane'               &
                                                          , '1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-RPPosition[$]'      , 'Position of the record plane in RPDirection.'                     &
                                                          , '0.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-FilenameRecordPart' , 'Specifying filename for load_from_file init.'                     &
                                                          , 'data/recordplane_')
CALL prms%CreateLogicalOption(  'Part-PartIndex'          , 'Flag to give each particle an unique index'                      &
                                                               , '.FALSE.')

END SUBROUTINE DefineParametersParticleAnalyze


SUBROUTINE InitParticleAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars   ,ONLY: doParticleAnalyze
USE MOD_Particle_Analyze_Vars   ,ONLY: doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcEkin,CalcPartBalance
USE MOD_Particle_Analyze_Vars   ,ONLY: ParticleAnalyzeInitIsDone,nSpecAnalyze
USE MOD_Particle_Analyze_Vars   ,ONLY: nPartIn,nPartOut,PartPath,PartEkin,PartEkinOut,PartEkinIn,nPartInTmp,PartEkinInTmp
USE MOD_Particle_Vars           ,ONLY: nSpecies,PDM
USE MOD_ReadInTools             ,ONLY: GETLOGICAL,GETINT,GETSTR,GETINTARRAY,GETREALARRAY,GETREAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF (ParticleAnalyzeInitIsDone) THEN
  CALL abort(__STAMP__,'InitParticleAnalyse already called.',999,999.)
  RETURN
END IF

!SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ANALYZE...'

doParticleAnalyze    = .FALSE.
nSpecAnalyze = MERGE(nSpecies + 1,1,nSpecies.GT.1)

CalcEkin = GETLOGICAL('CalcKineticEnergy')
IF (CalcEkin) THEN
  doParticleAnalyze = .TRUE.
  SDEALLOCATE(PartEkin)
  ALLOCATE( PartEkin (nSpecAnalyze))
  PartEkin = 0.
END IF

! Calculate number and kinetic energy of particles entering / leaving the domain
CalcPartBalance = GETLOGICAL('CalcPartBalance')
IF (CalcPartBalance) THEN
  doParticleAnalyze = .TRUE.
  SDEALLOCATE(nPartIn)
  SDEALLOCATE(nPartOut)
  SDEALLOCATE(PartEkinIn)
  SDEALLOCATE(PartEkinOut)
  ALLOCATE( nPartIn    (nSpecAnalyze)    &
          , nPartOut   (nSpecAnalyze)    &
          , PartEkinOut(nSpecAnalyze)    &
          , PartEkinIn (nSpecAnalyze))
  nPartIn     = 0
  nPartOut    = 0
  PartEkinOut = 0.
  PartEkinIn  = 0.

  SDEALLOCATE( nPartInTmp)
  SDEALLOCATE( PartEkinInTmp)
  ALLOCATE( nPartInTmp(nSpecies)     &
          , PartEkinInTmp(nSpecies))
  PartEkinInTmp = 0.
  nPartInTmp    = 0
END IF

doParticleDispersionTrack    = GETLOGICAL('Part-TrackDispersion')
doParticlePathTrack          = GETLOGICAL('Part-TrackPath'      )
IF (doParticleDispersionTrack .AND. doParticlePathTrack) &
  CALL CollectiveStop(__STAMP__,'doParticleDispersionTrack and doParticlePathTrack cannot be enabled at the same time!')

IF(doParticleDispersionTrack .OR. doParticlePathTrack) THEN
  ALLOCATE(PartPath(1:3,PDM%maxParticleNumber))
  PartPath = 0.
END IF

ParticleAnalyzeInitIsDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE ANALYZE DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitParticleAnalyze


SUBROUTINE ParticleAnalyze(t,iter)
!==================================================================================================================================
!> Controls particle analysis routines and is called at analyze time levels
!> - writes load balancing statistics
!> - calls impact sampling output
!> - calls impact tracking output
!> - informs about lost particles
!> - calls kinetic energy balance
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars              ,ONLY: nWriteData
USE MOD_Particle_Analyze_Tools    ,ONLY: ParticleRecord
USE MOD_Particle_Analyze_Vars     ,ONLY: doParticleAnalyze,CalcEkin
USE MOD_Particle_Analyze_Vars     ,ONLY: RecordPart
USE MOD_Particle_Boundary_Vars    ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars    ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Analyze ,ONLY: CalcSurfaceValues,WriteBoundaryParticleToHDF5
USE MOD_Particle_Output           ,ONLY: WriteParticleAnalyze,WriteInfoStdOut
USE MOD_Restart_Vars              ,ONLY: RestartTime
USE MOD_TimeDisc_Vars             ,ONLY: doFinalize,writeCounter
USE MOD_Particle_Tracking_Vars    ,ONLY: CountNbOfLostParts
#if USE_LOADBALANCE
USE MOD_LoadDistribution          ,ONLY: WriteElemTimeStatistics
#endif /*LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< current simulation time
INTEGER(KIND=8),INTENT(IN)      :: iter                   !< current iteration
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

#if USE_LOADBALANCE
! Create .csv file for performance analysis and load balance: write header line
CALL WriteElemTimeStatistics(WriteHeader=.TRUE.,iter=iter)
#endif /*LOADBALANCE*/

! Write information to console output
IF (CountNbOfLostParts) THEN
  CALL WriteInfoStdOut()
END IF

IF (CalcEkin) THEN
  CALL CalcKineticEnergy()
END IF

IF (doParticleAnalyze) THEN
  CALL WriteParticleAnalyze()
END IF

! Keep DG and particle output synchronous for nWriteData > 1
IF ((writeCounter.NE.nWriteData) .AND. .NOT.doFinalize .AND. t.NE.RestartTime) RETURN

! Calculate cumulative particle surface impact sampling data
IF (WriteMacroSurfaceValues) THEN
  CALL CalcSurfaceValues
END IF

! Write individual particle surface impact tracking data
IF (doParticleImpactTrack) THEN
  CALL WriteBoundaryParticleToHDF5(OutputTime=t)
END IF

! Write individual particle record plane data
IF (RecordPart.GT.0) THEN
  CALL ParticleRecord(t,writeToBinary=.TRUE.)
END IF

END SUBROUTINE ParticleAnalyze


SUBROUTINE ParticleInformation()
!==================================================================================================================================
!> Displays the actual status of the particle phase in the simulation and counts the number of impacts
!==================================================================================================================================
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Boundary_Vars    ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Vars    ,ONLY: PartStateBoundaryVecLength
USE MOD_Particle_Boundary_Vars    ,ONLY: ImpactnGlob,ImpactnLoc,ImpactOffset
USE MOD_Particle_Vars             ,ONLY: PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nParticleOnProc,iPart
INTEGER                        :: nParticleInDomain
#if USE_MPI
INTEGER                        :: sendbuf(2),recvbuf(2)
#endif
!==================================================================================================================================

! Count number of particles on local proc
nParticleOnProc = 0
DO iPart = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) nParticleOnProc = nParticleOnProc + 1
END DO
#if USE_MPI
! Gather number of particles on all procs
CALL MPI_REDUCE(nParticleOnProc,nParticleInDomain,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
#else
nParticleInDomain = nParticleOnProc
#endif

! Find amount of recorded impacts on current proc
ImpactnLoc  = PartStateBoundaryVecLength

!>> Sum up particles from the other procs
#if USE_MPI
sendbuf(1) = ImpactnLoc
recvbuf    = 0
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
!>> Offset of each proc is the sum of the particles on the previous procs
ImpactOffset = recvbuf(1)
sendbuf(1)   = recvbuf(1) + ImpactnLoc
!>> Last proc knows the global number
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER,nProcessors-1,MPI_COMM_FLEXI,iError)
!>> Gather the global number and communicate to root (MPIRank.EQ.0)
ImpactnGlob  = sendbuf(1)
#else
ImpactOffset   = 0
ImpactnGlob    = PartStateBoundaryVecLength
#endif

! Output particle and impact information
SWRITE(UNIT_stdOut,'(A14,I16)')' # Particle : ', nParticleInDomain
IF (doParticleImpactTrack) THEN
SWRITE(UNIT_stdOut,'(A14,I16)')' # Impacts  : ', ImpactnGlob
END IF

END SUBROUTINE ParticleInformation


SUBROUTINE CalcKineticEnergy()
!===================================================================================================================================
! compute the kinetic energy of particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Tools,ONLY: CalcEkinPart
USE MOD_Particle_Analyze_Vars ,ONLY: PartEkin,nSpecAnalyze
USE MOD_Particle_Vars         ,ONLY: PartSpecies,PDM
#if USE_MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
REAL              :: EkinPart
!===================================================================================================================================

PartEkin = 0.

! nSpecAnalyze > 1, calculate species and sum
IF (nSpecAnalyze.GT.1) THEN
  DO i = 1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      EkinPart = CalcEkinPart(i)
      PartEkin(nSpecAnalyze)   = PartEkin(nSpecAnalyze)   + EkinPart
      PartEkin(PartSpecies(i)) = PartEkin(PartSpecies(i)) + EkinPart
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
! nSpecAnalyze = 1 : only 1 species
ELSE
  DO i = 1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      PartEkin(PartSpecies(i)) = PartEkin(PartSpecies(i)) + CalcEkinPart(i)
    END IF ! particle inside
  END DO ! 1,PDM%ParticleVecLength
END IF

#if USE_MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PartEkin,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(PartEkin    ,0       ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcKineticEnergy


SUBROUTINE FinalizeParticleAnalyze()
!===================================================================================================================================
! Finalizes variables necessary for analyze subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(PartEkin)

SDEALLOCATE(nPartIn)
SDEALLOCATE(nPartOut)
SDEALLOCATE(PartEkinIn)
SDEALLOCATE(PartEkinOut)

SDEALLOCATE(nPartInTmp)
SDEALLOCATE(PartEkinInTmp)

SDEALLOCATE(RPP_Plane)

ParticleAnalyzeInitIsDone = .FALSE.

END SUBROUTINE FinalizeParticleAnalyze


END MODULE MOD_Particle_Analyze
