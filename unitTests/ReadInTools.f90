#include "flexi.h"

!==================================================================================================================================
!> Unit test 'ReadInToolsUnitTest'
!> Test the module: MOD_ReadInTools
!==================================================================================================================================
PROGRAM ReadInToolsUnitTest
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_MPI,         ONLY: InitMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nArgs,i
CHARACTER(LEN=255),PARAMETER   :: FileName='ReadInTools.ini'
INTEGER :: intOpt
INTEGER :: intOpt_def
INTEGER :: intOpt_mult
INTEGER :: intOpt_mult_A
REAL    :: realOpt
REAL    :: realOptsci
REAL    :: realOpt_def
REAL    :: realOpt_defsci
REAL    :: realOpt_mult
REAL    :: realOpt_mult_A
REAL    :: realOpt_multsci
REAL    :: realOpt_multsci_A
LOGICAL :: logOpt
LOGICAL :: logOpt_def
LOGICAL :: logOpt_mult
LOGICAL :: logOpt_mult_A
character(len=255) :: strOpt
character(len=255) :: strOpt_def
character(len=255) :: strOpt_mult
character(len=255) :: strOpt_mult_A

INTEGER :: intArrayOpt(3)
INTEGER :: intArrayOpt_def(3)
INTEGER :: intArrayOpt_mult(2)
INTEGER :: intArrayOpt_mult_A(2)
REAL    :: realArrayOpt(3)
REAL    :: realArrayOptsci(2)
REAL    :: realArrayOpt_def(3)
REAL    :: realArrayOpt_defsci(2)
REAL    :: realArrayOpt_mult(2)
REAL    :: realArrayOpt_mult_A(2)
REAL    :: realArrayOpt_multsci(2)
REAL    :: realArrayOpt_multsci_A(2)
LOGICAL :: logArrayOpt(2)
LOGICAL :: logArrayOpt_def(2)
LOGICAL :: logArrayOpt_mult(2)
LOGICAL :: logArrayOpt_mult_A(2)
!character(len=255) :: strArrayOpt(2)
!character(len=255) :: strArrayOpt_def(2)
!character(len=255) :: strArrayOpt_mult(2)
!character(len=255) :: strArrayOpt_mult_A(2)
#if USE_PARTICLES
INTEGER,ALLOCATABLE  :: Dollar(:)
INTEGER              :: nDollar,iDollar,nDollar2,iDollar2
CHARACTER(32)        :: hilf,hilf2
INTEGER, ALLOCATABLE :: IntOption(:)
REAL, ALLOCATABLE    :: RealOption(:,:)
#endif /*USE_PARTICLES*/
!==================================================================================================================================
CALL InitMPI()
! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) CALL Abort(__STAMP__,'ERROR - Unknown command line argument.')

CALL prms%SetSection("UnitTest")
CALL prms%CreateIntOption('intOpt'                     , "Description IntOpt")
CALL prms%CreateIntOption('intOpt_def'                 , "Description IntOpt with default value"             , '4')
CALL prms%CreateIntOption('intOpt_mult'                , "Description IntOpt multiple"                       , multiple=.TRUE.)
CALL prms%CreateRealOption('realOpt'                   , "Description RealOpt")
CALL prms%CreateRealOption('realOptsci'                , "Description RealOpt")
CALL prms%CreateRealOption('realOpt_def'               , "Description RealOpt with default value"            , '-1.00')
CALL prms%CreateRealOption('realOpt_defsci'            , "Description RealOpt with default value scientific" , '0.3e-7')
CALL prms%CreateRealOption('realOpt_mult'              , "Description RealOpt multiple"                      , multiple=.TRUE.)
CALL prms%CreateRealOption('realOpt_multsci'           , "Description RealOpt multiple"                      , multiple=.TRUE.)
CALL prms%CreateLogicalOption('logOpt'                 , "Description LogOpt")
CALL prms%CreateLogicalOption('logOpt_def'             , "Description LogOpt with default value"             , 'T')
CALL prms%CreateLogicalOption('logOpt_mult'            , "Description LogOpt multiple"                       , multiple=.TRUE.)
CALL prms%CreateStringOption('strOpt'                  , "Description StrOpt")
CALL prms%CreateStringOption('strOpt_def'              , "Description StrOpt with default value"             , 'dummyValue')
CALL prms%CreateStringOption('strOpt_mult'             , "Description StrOpt multiple"                       , multiple=.TRUE.)

CALL prms%CreateIntArrayOption('intArrayOpt'           , "Description IntOpt")
CALL prms%CreateIntArrayOption('intArrayOpt_def'       , "Description IntOpt with default value"             , '-1,0,-3')
CALL prms%CreateIntArrayOption('intArrayOpt_mult'      , "Description IntOpt multiple"                       , multiple=.TRUE.)
CALL prms%CreateRealArrayOption('realArrayOpt'         , "Description RealOpt")
CALL prms%CreateRealArrayOption('realArrayOptsci'      , "Description RealOpt")
CALL prms%CreateRealArrayOption('realArrayOpt_def'     , "Description RealOpt with default value"            , '-1.00,5.,22')
CALL prms%CreateRealArrayOption('realArrayOpt_defsci'  , "Description RealOpt with default value scientific" , '0.3e-7,-5e2')
CALL prms%CreateRealArrayOption('realArrayOpt_mult'    , "Description RealOpt multiple"                      , multiple=.TRUE.)
CALL prms%CreateRealArrayOption('realArrayOpt_multsci' , "Description RealOpt multiple"                      , multiple=.TRUE.)
CALL prms%CreateLogicalArrayOption('logArrayOpt'       , "Description LogOpt")
CALL prms%CreateLogicalArrayOption('logArrayOpt_def'   , "Description LogOpt with default value"             , '(/T,F/)')
CALL prms%CreateLogicalArrayOption('logArrayOpt_mult'  , "Description LogOpt multiple"                       , multiple=.TRUE.)
!CALL prms%CreateStringArrayOption('strArrayOpt'        , "Description StrOpt")
!CALL prms%CreateStringArrayOption('strArrayOpt_def'    , "Description StrOpt with default value", 'dum1,dum2')
!CALL prms%CreateStringArrayOption('strArrayOpt_mult'   , "Description StrOpt multiple", multiple=.TRUE.)

#if USE_PARTICLES
CALL prms%CreateIntOption('Dollar[$]-IntOptionUndef'             , 'Single dollar variable undefined.'                    , numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Dollar[$]-Dollar[$]-RealOptionUndef' , 'Double dollar variable undefined.'                    , numberedmulti=.TRUE.)

CALL prms%CreateIntOption('Dollar[$]-IntOption'                  , 'Single dollar variable with regular default.', '-1'   , numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Dollar[$]-Dollar[$]-RealOption'      , 'Double dollar variable with regular default.', '-1.0' , numberedmulti=.TRUE.)

CALL prms%CreateIntOption('Dollar[$]-IntOptionDef'               , 'Single dollar variable with both defaults.'  , '-1'   , numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Dollar[$]-Dollar[$]-RealOptionDef'   , 'Double dollar variable with both defaults.'  , '-1.0' , numberedmulti=.TRUE.)
#endif /*USE_PARTICLES*/

!CALL PrintDefaultParameterFile(.FALSE.)
CALL prms%read_options(FileName)

intOpt           = GETINT('intOpt'        )
intOpt_def       = GETINT('intOpt_def'    )
IF (intOpt.NE.intOpt_def) CALL Abort(__STAMP__,"intOpt failed")
intOpt_mult      = GETINT('intOpt_mult'   )
intOpt_mult_A    = GETINT('intOpt_mult'   )
IF (intOpt_mult.NE.intOpt_mult_A) CALL Abort(__STAMP__,"intOpt_mult failed")
realOpt          = GETREAL('realOpt'        )
realOpt_def      = GETREAL('realOpt_def'    )
IF (realOpt.NE.realOpt_def) CALL Abort(__STAMP__,"realOpt failed")
realOptsci         = GETREAL('realOptsci'       )
realOpt_defsci     = GETREAL('realOpt_defsci'    )
IF (realOptsci.NE.realOpt_defsci) CALL Abort(__STAMP__,"realOptsci failed")
realOpt_mult     = GETREAL('realOpt_mult'   )
realOpt_mult_A   = GETREAL('realOpt_mult'   )
IF (realOpt_mult.NE.realOpt_mult_A) CALL Abort(__STAMP__,"realOpt_mult failed")
realOpt_multsci    = GETREAL('realOpt_multsci'  )
realOpt_multsci_A  = GETREAL('realOpt_multsci'  )
IF (realOpt_multsci.NE.realOpt_multsci_A) CALL Abort(__STAMP__,"realOpt_multsci failed")
logOpt           = GETLOGICAL('logOpt'        )
logOpt_def       = GETLOGICAL('logOpt_def'    )
IF (logOpt.NEQV.logOpt_def) CALL Abort(__STAMP__,"logOpt failed")
logOpt_mult      = GETLOGICAL('logOpt_mult'   )
logOpt_mult_A    = GETLOGICAL('logOpt_mult'   )
IF (logOpt_mult.NEQV.logOpt_mult_A) CALL Abort(__STAMP__,"logOpt_mult failed")
strOpt           = GETSTR('strOpt'        )
strOpt_def       = GETSTR('strOpt_def'    )
IF (strOpt.NE.strOpt_def) CALL Abort(__STAMP__,"strOpt failed")
strOpt_mult      = GETSTR('strOpt_mult'   )
strOpt_mult_A    = GETSTR('strOpt_mult'   )
IF (strOpt_mult.NE.strOpt_mult_A) CALL Abort(__STAMP__,"strOpt_mult failed")

intArrayOpt             = GETINTARRAY('intArrayOpt',3    )
intArrayOpt_def         = GETINTARRAY('intArrayOpt_def',3)
DO i=1,3
  IF (intArrayOpt(i).NE.intArrayOpt_def(i)) CALL Abort(__STAMP__,"intArrayOpt failed")
END DO
intArrayOpt_mult        = GETINTARRAY('intArrayOpt_mult', 2)
intArrayOpt_mult_A      = GETINTARRAY('intArrayOpt_mult', 2)
DO i=1,2
  IF (intArrayOpt_mult(i).NE.intArrayOpt_mult_A(i)) CALL Abort(__STAMP__,"intArrayOpt_mult failed")
END DO
realArrayOpt            = GETREALARRAY('realArrayOpt', 3)
realArrayOpt_def        = GETREALARRAY('realArrayOpt_def', 3)
DO i=1,3
  IF (realArrayOpt(i).NE.realArrayOpt_def(i)) CALL Abort(__STAMP__,"realArrayOpt failed")
END DO
realArrayOptsci         = GETREALARRAY('realArrayOptsci', 2)
realArrayOpt_defsci     = GETREALARRAY('realArrayOpt_defsci', 2)
DO i=1,2
  IF (realArrayOptsci(i).NE.realArrayOpt_defsci(i)) CALL Abort(__STAMP__,"realArrayOptsci failed")
END DO
realArrayOpt_mult       = GETREALARRAY('realArrayOpt_mult', 2)
realArrayOpt_mult_A     = GETREALARRAY('realArrayOpt_mult', 2)
DO i=1,2
  IF (realArrayOpt_mult(i).NE.realArrayOpt_mult_A(i)) CALL Abort(__STAMP__,"realArrayOpt_mult failed")
END DO
realArrayOpt_multsci    = GETREALARRAY('realArrayOpt_multsci', 2)
realArrayOpt_multsci_A  = GETREALARRAY('realArrayOpt_multsci', 2)
DO i=1,2
  IF (realArrayOpt_multsci(i).NE.realArrayOpt_multsci_A(i)) CALL Abort(__STAMP__,"realArrayOpt_multsci failed")
END DO
logArrayOpt             = GETLOGICALARRAY('logArrayOpt', 2)
logArrayOpt_def         = GETLOGICALARRAY('logArrayOpt_def', 2)
DO i=1,2
  IF (logArrayOpt(i).NEQV.logArrayOpt_def(i)) CALL Abort(__STAMP__,"logArrayOpt failed")
END DO
logArrayOpt_mult        = GETLOGICALARRAY('logArrayOpt_mult', 2)
logArrayOpt_mult_A      = GETLOGICALARRAY('logArrayOpt_mult', 2)
DO i=1,2
  IF (logArrayOpt_mult(i).NEQV.logArrayOpt_mult_A(i)) CALL Abort(__STAMP__,"logArrayOpt_mult failed")
END DO
!strArrayOpt             = GETSTRARRAY('strArrayOpt', 2)
!strArrayOpt_def         = GETSTRARRAY('strArrayOpt_def', 2)
!DO i=1,2
  !IF (strArrayOpt(i).NE.strArrayOpt_def(i)) CALL Abort(__STAMP__,"strArrayOpt failed")
!END DO
!strArrayOpt_mult        = GETSTRARRAY('strArrayOpt_mult', 2)
!strArrayOpt_mult_A      = GETSTRARRAY('strArrayOpt_mult', 2)
!DO i=1,2
  !IF (strArrayOpt_mult(i).NE.strArrayOpt_mult_A(i)) CALL Abort(__STAMP__,"strArrayOpt_mult failed")
!END DO

#if USE_PARTICLES
! Test single and double dollar variables
nDollar  = ABS(intOpt)
nDollar2 = ABS(intOpt)
ALLOCATE(Dollar(1:nDollar))
Dollar   = 0
ALLOCATE(RealOption(nDollar,nDollar2))
RealOption = 0
ALLOCATE(IntOption(nDollar))
IntOption  = 0
! Loop twice. The default value <= zero is not allowed to appear
DO i = 1, 2
  WRITE(UNIT_StdOut,'(A,I0,A)') "--------------",i,"--------------"
  ! Test parameters read-in without a default in CreateOption nor GETVAR (Option to force a user input)
  DO iDollar=1, nDollar
    WRITE(UNIT=hilf,FMT='(I0)') iDollar
    IntOption(iDollar) = GETINT('Dollar'//TRIM(hilf)//'-IntOptionUndef')
    DO iDollar2 = 1, nDollar2
      WRITE(UNIT=hilf2,FMT='(I0)') iDollar2
      RealOption(iDollar,iDollar2) = GETREAL('Dollar'//TRIM(hilf)//'-Dollar'//TRIM(hilf2)//'-RealOptionUndef')
    END DO ! iDollar2 = 1, nDollar2
  END DO ! iDollar=1, nDollar
  IF(ANY(IntOption .LE.0)) CALL Abort(__STAMP__,'ERROR: IntOption cannot be greater than zero!')
  IF(ANY(RealOption.LE.0)) CALL Abort(__STAMP__,'ERROR: RealOption cannot be greater than zero!')
  ! Test parameters read-in with a default in CreateOption and not in GETVAR (Preferred option)
  RealOption = 0
  IntOption  = 0
  DO iDollar=1, nDollar
    WRITE(UNIT=hilf,FMT='(I0)') iDollar
    IntOption(iDollar) = GETINT('Dollar'//TRIM(hilf)//'-IntOption')
    DO iDollar2 = 1, nDollar2
      WRITE(UNIT=hilf2,FMT='(I0)') iDollar2
      RealOption(iDollar,iDollar2) = GETREAL('Dollar'//TRIM(hilf)//'-Dollar'//TRIM(hilf2)//'-RealOption')
    END DO ! iDollar2 = 1, nDollar2
  END DO ! iDollar=1, nDollar
  IF(ANY(IntOption .LE.0)) CALL Abort(__STAMP__,'ERROR: IntOption cannot be greater than zero!')
  IF(ANY(RealOption.LE.0)) CALL Abort(__STAMP__,'ERROR: RealOption cannot be greater than zero!')
  ! Test parameters read-in with a default in CreateOption and in GETVAR (Rare cases, where e.g. HUGE is required)
  RealOption = 0
  IntOption = 0
  DO iDollar=1, nDollar
    WRITE(UNIT=hilf,FMT='(I0)') iDollar
    IntOption(iDollar) = GETINT('Dollar'//TRIM(hilf)//'-IntOptionDef','-1')
    DO iDollar2 = 1, nDollar2
      WRITE(UNIT=hilf2,FMT='(I0)') iDollar2
      RealOption(iDollar,iDollar2) = GETREAL('Dollar'//TRIM(hilf)//'-Dollar'//TRIM(hilf2)//'-RealOptionDef','-1.')
    END DO ! iDollar2 = 1, nDollar2
  END DO ! iDollar=1, nDollar
  IF(ANY(IntOption .LE.0)) CALL Abort(__STAMP__,'ERROR: IntOption cannot be greater than zero!')
  IF(ANY(RealOption.LE.0)) CALL Abort(__STAMP__,'ERROR: RealOption cannot be greater than zero!')
  WRITE (*,*) "IntOption =", IntOption
  DO iDollar=1, nDollar
    WRITE (*,*) "RealOption(iDollar,:) =", RealOption(iDollar,:)
  END DO
  IF(i.EQ.1) THEN
    CALL prms%WriteUnused()
    IF(prms%count_unread().GT.0) CALL Abort(__STAMP__,'ERROR: Unused parameters were found but are not allowed!')
  END IF
  CALL prms%finalize(.TRUE.) ! is the same as CALL FinalizeParameters(), but considers load balancing
END DO ! i = 1, 2

DEALLOCATE(Dollar)
DEALLOCATE(RealOption)
DEALLOCATE(IntOption)

CALL prms%finalize(.FALSE.) ! is the same as CALL FinalizeParameters(), but considers load balancing
#endif /*USE_PARTICLES*/

#if USE_MPI
! free the communicator
CALL MPI_BARRIER  (MPI_COMM_FLEXI,IERROR)
CALL MPI_COMM_FREE(MPI_COMM_FLEXI,IERROR)
! we also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'MPI finalize error')
#endif

END PROGRAM ReadInToolsUnitTest
