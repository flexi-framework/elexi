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

MODULE MOD_Flexi

IMPLICIT NONE
PRIVATE
SAVE

INTERFACE InitFlexi
   MODULE PROCEDURE InitFlexi
END INTERFACE

INTERFACE FinalizeFlexi
   MODULE PROCEDURE FinalizeFlexi
END INTERFACE

PUBLIC::InitFlexi,FinalizeFlexi

CONTAINS

!==================================================================================================================================
!> Initialization of the computation
!==================================================================================================================================
SUBROUTINE InitFlexi(nArgs_In,Args_In &
#if USE_MPI
                    ,mpi_comm_loc     &
#endif /*USE_MPI*/
                    )
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,      ONLY:InitializationWallTime,StartTime
USE MOD_Commandline_Arguments
USE MOD_PreProc
USE MOD_Analyze,           ONLY:DefineParametersAnalyze,InitAnalyze
USE MOD_DG,                ONLY:InitDG
USE MOD_Eos,               ONLY:DefineParametersEos
USE MOD_Equation,          ONLY:DefineParametersEquation,InitEquation
USE MOD_Exactfunc,         ONLY:DefineParametersExactFunc
USE MOD_Filter,            ONLY:DefineParametersFilter,InitFilter
USE MOD_Implicit,          ONLY:DefineParametersImplicit,InitImplicit
USE MOD_Interpolation,     ONLY:DefineParametersInterpolation,InitInterpolation
USE MOD_IO_HDF5,           ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Mesh,              ONLY:DefineParametersMesh,InitMesh
USE MOD_Mesh_Check,        ONLY:CheckMesh
USE MOD_Mortar,            ONLY:InitMortar
USE MOD_MPI,               ONLY:DefineParametersMPI,InitMPI
USE MOD_Output,            ONLY:DefineParametersOutput,InitOutput
USE MOD_Overintegration,   ONLY:DefineParametersOverintegration,InitOverintegration
USE MOD_Precond,           ONLY:DefineParametersPrecond
USE MOD_ReadInTools,       ONLY:prms,PrintDefaultParameterFile,ExtractParameterFile
USE MOD_RecordPoints,      ONLY:DefineParametersRecordPoints,InitRecordPoints
USE MOD_Restart,           ONLY:DefineParametersRestart,InitRestart,Restart
USE MOD_Sponge,            ONLY:DefineParametersSponge,InitSponge
USE MOD_StringTools,       ONLY:STRICMP, GetFileExtension
USE MOD_TestCase,          ONLY:DefineParametersTestcase
USE MOD_TimeDisc,          ONLY:TimeDisc
USE MOD_TimeDisc_Functions,ONLY:DefineParametersTimedisc,InitTimeDisc
USE MOD_Unittest,          ONLY:GenerateUnittestReferenceData
#if FV_ENABLED
USE MOD_FV,                ONLY:DefineParametersFV,InitFV
USE MOD_FV_Basis,          ONLY:InitFV_Basis
USE MOD_Indicator,         ONLY:DefineParametersIndicator,InitIndicator
#endif /*FV_ENABLED*/
#if PARABOLIC
USE MOD_Lifting,           ONLY:DefineParametersLifting,InitLifting
#endif /*PARABOLIC*/
#if USE_MPI
USE MOD_MPI,               ONLY:InitMPIvars
USE MOD_MPI_Shared,        ONLY:DefineParametersMPIShared,InitMPIShared
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance,       ONLY:DefineParametersLoadBalance,InitLoadBalance
USE MOD_Restart_Vars,      ONLY:DoRestart,RestartFile
#endif /*USE_LOADBALANCE*/
#if USE_PARTICLES
USE MOD_Particle_Init,     ONLY:DefineParametersParticles,InitParticleGlobals,InitParticles
#if USE_MPI
USE MOD_Particle_MPI,      ONLY:DefineParticleMPI,InitParticleMPI
#endif /*USE_MPI*/
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                     :: nArgs_In
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: Args_In(*)
#if USE_MPI
MPI_TYPE_COMM,INTENT(IN),OPTIONAL      :: mpi_comm_loc
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                 :: userblockFound
CHARACTER(LEN=255)      :: RestartFile_loc = ''
!==================================================================================================================================
CALL SetStackSizeUnlimited()

#if USE_MPI
IF(PRESENT(mpi_comm_loc))THEN
  CALL InitMPI(mpi_comm_loc)
ELSE
#endif /*USE_MPI*/
  CALL InitMPI()
#if USE_MPI
END IF
#endif /*USE_MPI*/

IF(nArgs_In.EQ.0)THEN
  CALL ParseCommandlineArguments()
ELSE
  CALL ParseCommandlineArguments(Args_In(1:nArgs_In))
END IF

! Check if the number of arguments is correct
IF (nArgs.GT.2) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: flexi parameter.ini [restart.h5] or flexi --help'// &
  '[option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF

ParameterFile = Args(1)
IF (nArgs.GT.1) THEN
  RestartFile_loc = Args(2)
#if USE_LOADBALANCE
  ! LoadBalance needs the information that we are performing a restart
  DoRestart = .TRUE.
  ! LoadBalance needs the information which is the RestartFile
  RestartFile = RestartFile_loc
#endif
ELSE IF (STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
  ParameterFile = ".flexi.ini"
  CALL ExtractParameterFile(Args(1), ParameterFile, userblockFound)
  IF (.NOT.userblockFound) THEN
    CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(Args(1))//"'")
  END IF
  RestartFile_loc = Args(1)
#if USE_LOADBALANCE
  ! LoadBalance needs the information that we are performing a restart
  DoRestart = .TRUE.
  ! LoadBalance needs the information which is the RestartFile
  RestartFile = RestartFile_loc
#endif
END IF

CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersEos()
CALL DefineParametersEquation()
CALL DefineParametersExactFunc()
CALL DefineParametersTestcase()
CALL DefineParametersFilter()
CALL DefineParametersOverintegration()
#if FV_ENABLED
CALL DefineParametersIndicator()
CALL DefineParametersFV()
#endif
#if PARABOLIC
CALL DefineParametersLifting ()
#endif /*PARABOLIC*/
CALL DefineParametersSponge()
CALL DefineParametersTimedisc()
CALL DefineParametersImplicit()
CALL DefineParametersPrecond()
CALL DefineParametersAnalyze()
CALL DefineParametersRecordPoints()
#if USE_MPI
CALL DefineParametersMPIShared()
#if USE_LOADBALANCE
CALL DefineParametersLoadBalance()
#endif /*USE_LOADBALANCE*/
#if USE_PARTICLES
CALL DefineParticleMPI
#endif /*USE_PARTICLES*/
#endif /*USE_MPI*/
#if USE_PARTICLES
CALL DefineParametersParticles()
#endif /*USE_PARTICLES*/

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
#if USE_MPI
  ! free the communicator
  CALL MPI_BARRIER  (MPI_COMM_FLEXI,iError)
  CALL MPI_COMM_FREE(MPI_COMM_FLEXI,iError)
  ! we also have to finalize MPI itself here
  CALL MPI_FINALIZE(iError)
  IF(iError.NE.MPI_SUCCESS) STOP 'MPI finalize error'
#endif
  STOP
END IF

CALL prms%read_options(ParameterFile)

CALL InitIOHDF5()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
"           __________________   _______              __________________   ______      ______   __________________ "
SWRITE(UNIT_stdOut,'(A)') &
"          /                 /) /      /)            /                 /) /      |   _/     /) /                 /)"
SWRITE(UNIT_stdOut,'(A)') &
"         /       __________// /      //            /      ___________// /__     |_/´     _// /____       ______// "
SWRITE(UNIT_stdOut,'(A)') &
"        /      /)__________) /      //            /      /)__________)  (__|          _/´_)  )___/      /)_____)  "
SWRITE(UNIT_stdOut,'(A)') &
"       /      //___         /      //            /      //___              |       _/´_/´       /      //         "
SWRITE(UNIT_stdOut,'(A)') &
"      /           /)       /      //            /           /)             |     /´ /´         /      //          "

SWRITE(UNIT_stdOut,'(A)') &
"     /      _____//       /      //            /      _____//            _/´     |/´          /      //           "
SWRITE(UNIT_stdOut,'(A)') &
"    /      /)____)       /      //            /      /)____)          _/´        |           /      //            "

SWRITE(UNIT_stdOut,'(A)') &
"   /      //            /      //_________   /      //_________   __/´     _     |__   _____/      //____         "
SWRITE(UNIT_stdOut,'(A)') &
"  /      //            /                 /) /                 /) /      _/´ |      /) /                 /)        "
SWRITE(UNIT_stdOut,'(A)') &
" /______//            /_________________// /_________________// /_____/` _/´|_____// /_________________//         "
SWRITE(UNIT_stdOut,'(A)') &
" )______)             )_________________)  )_________________)  )_____)/´   )_____)  )_________________)          "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')
! Measure init duration
StartTime=FLEXITIME()

! Initialization
CALL InitInterpolation()
#if FV_ENABLED
CALL InitFV_Basis()
#endif
CALL InitMortar()
CALL InitOutput()
#if USE_PARTICLES
CALL InitParticleGlobals
#endif /*USE_PARTICLES*/
#if USE_MPI
CALL InitMPIShared()
#endif /*USE_MPÌ*/
CALL InitMesh(meshMode=2)
CALL InitRestart(RestartFile_loc)
CALL InitFilter()
CALL InitOverintegration()
#if USE_MPI
CALL InitMPIvars()
#endif
CALL CheckMesh()
CALL InitEquation()
CALL InitDG()
#if FV_ENABLED
CALL InitIndicator()
CALL InitFV()
#endif
#if PARABOLIC
CALL InitLifting()
#endif /*PARABOLIC*/
CALL InitSponge()
CALL InitTimeDisc()
CALL InitAnalyze()
CALL InitImplicit()
CALL InitRecordpoints()
CALL Restart()
#if USE_PARTICLES && USE_MPI
CALL InitParticleMPI()
#endif /*USE_PARTICLES && USE_MPI*/
#if USE_LOADBALANCE
CALL InitLoadBalance()
#endif /*USE_LOADBALANCE*/
#if USE_PARTICLES
CALL InitParticles()
#endif /*USE_PARTICLES*/
CALL prms%WriteUnused()

! Measure init duration
InitializationWallTime = FLEXITIME()-StartTime
SWRITE(UNIT_stdOut,'(132("="))')
CALL DisplayMessageAndTime(InitializationWallTime,'INITIALIZATION DONE!',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("="))')

! Generate Unittest Data
IF (doGenerateUnittestReferenceData) THEN
#if USE_MPI
  CALL CollectiveStop(__STAMP__, "Can't generate Unittest Reference Data with MPI enabled.")
#endif
  CALL GenerateUnittestReferenceData()
  CALL FinalizeFlexi()
  CALL EXIT()
END IF

END SUBROUTINE InitFlexi


!==================================================================================================================================
!> Finalize the computation.
!==================================================================================================================================
SUBROUTINE FinalizeFlexi()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,      ONLY:StartTime
USE MOD_Analyze,           ONLY:FinalizeAnalyze
USE MOD_Commandline_Arguments,ONLY:FinalizeCommandlineArguments
USE MOD_DG,                ONLY:FinalizeDG
USE MOD_Equation,          ONLY:FinalizeEquation
USE MOD_Filter,            ONLY:FinalizeFilter
USE MOD_Implicit,          ONLY:FinalizeImplicit
USE MOD_Interpolation,     ONLY:FinalizeInterpolation
USE MOD_IO_HDF5,           ONLY:FinalizeIOHDF5
USE MOD_Mesh,              ONLY:FinalizeMesh
USE MOD_Mortar,            ONLY:FinalizeMortar
USE MOD_Output,            ONLY:FinalizeOutput
USE MOD_Overintegration,   ONLY:FinalizeOverintegration
USE MOD_ReadInTools,       ONLY:FinalizeParameters
USE MOD_RecordPoints,      ONLY:FinalizeRecordPoints
USE MOD_Restart,           ONLY:FinalizeRestart
USE MOD_Sponge,            ONLY:FinalizeSponge
USE MOD_TimeDisc_Functions,ONLY:FinalizeTimeDisc
#if PARABOLIC
USE MOD_Lifting,           ONLY:FinalizeLifting
#endif /*PARABOLIC*/
#if USE_MPI
USE MOD_MPI,               ONLY:FinalizeMPI
USE MOD_MPI_Shared        ,ONLY:FinalizeMPIShared
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance,       ONLY:FinalizeLoadBalance
#endif /*USE_LOADBALANCE*/
#if FV_ENABLED
USE MOD_FV,                ONLY:FinalizeFV
USE MOD_FV_Basis,          ONLY:FinalizeFV_Basis
USE MOD_Indicator,         ONLY:FinalizeIndicator
#endif /*FV_ENABLED*/
#if USE_PARTICLES
USE MOD_Particle_Init,     ONLY:FinalizeParticles
#if USE_MPI
USE MOD_Particle_MPI,      ONLY:FinalizeParticleMPI
#endif /*USE_MPI*/
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Time                              !< Used to measure simulation time
!==================================================================================================================================
!Finalize
CALL FinalizeIOHDF5()
CALL FinalizeOutput()
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
#if PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/
#if FV_ENABLED
CALL FinalizeFV()
#endif
CALL FinalizeDG()
CALL FinalizeEquation()
CALL FinalizeInterpolation()
CALL FinalizeImplicit
CALL FinalizeTimeDisc()
CALL FinalizeRestart()
CALL FinalizeMesh()
CALL FinalizeMortar()
CALL FinalizeSponge()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
#if FV_ENABLED
CALL FinalizeFV_Basis()
CALL FinalizeIndicator()
#endif
#if USE_PARTICLES && USE_MPI
CALL FinalizeParticleMPI()
#endif /*USE_PARTICLES && USE_MPI*/
#if USE_LOADBALANCE
CALL FinalizeLoadBalance()
#endif /*USE_LOADBALANCE*/
#if USE_PARTICLES
CALL FinalizeParticles()
#endif /*USE_PARTICLES*/
#if USE_MPI
! Must be called last because MPI3 shared windows are deallocated on MPI_COMM_SHARED
CALL FinalizeMPIShared()
#endif /*MPI*/
! Measure simulation duration
Time=FLEXITIME()
CALL FinalizeParameters()
CALL FinalizeCommandlineArguments()
#if USE_MPI
! For flexilib MPI init/finalize is controlled by main program
CALL FinalizeMPI()
#endif

CALL DisplaySimulationTime(Time, StartTime, 'FINISHED!')
END SUBROUTINE FinalizeFlexi

END MODULE MOD_Flexi
