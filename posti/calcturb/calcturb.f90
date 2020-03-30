!===================================================================================================================================
! Copyright (c) 2010-2019  Prof. Claus-Dieter Munz
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
!===================================================================================================================================
#include "flexi.h"

!===================================================================================================================================
!> This tool will take a TimeAvg (and multiple state) file and calculate turbulence quantities from Mean and MeanSquare values
!===================================================================================================================================
!
! There are quantities that require the computation of gradients. In this case the DG operator 'DGTimeDerivative_weakForm' is
! called once to fill the gradients and the reconstruction of the FV subcell method. This requires the initialization of several
! modules of FLEXI. U is read via a call of 'Restart'. In the DGTimeDerivative_weakForm the primitive quantities U_Prim and
! gradUx/y/z as well as gradUxi/eta/zeta are filled. These are used to calculate the remaining quantities.
!
! * The calculation of derived quantities is performed on a arbitrary polynomial degree NCalc and afterwards interpolated to NVisu.
!   Default is PP_N. These require to reconstruct the solution first to the visu grid and afterwards can calculate the derived
!   quantities on the NVisu_FV grid.
!
!===================================================================================================================================
MODULE MOD_CalcTurb
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineCalcTurb
  MODULE PROCEDURE DefineCalcTurb
END INTERFACE

INTERFACE InitCalcTurb
  MODULE PROCEDURE InitCalcTurb
END INTERFACE

INTERFACE FinalizeCalcTurb
  MODULE PROCEDURE FinalizeCalcTurb
END INTERFACE

INTERFACE ALMOSTEQUAL
  MODULE PROCEDURE ALMOSTEQUAL
END INTERFACE

INTERFACE ALMOSTZERO
  MODULE PROCEDURE ALMOSTZERO
END INTERFACE ALMOSTZERO

PUBLIC :: DefineCalcTurb
PUBLIC :: InitCalcTurb
PUBLIC :: FinalizeCalcTurb
PUBLIC :: ALMOSTEQUAL
PUBLIC :: ALMOSTZERO
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> This routine is used to prepare everything we need for the calcturb tool.
!===================================================================================================================================
SUBROUTINE DefineCalcTurb
! MODULES
USE MOD_ReadInTools
USE MOD_CalcTurb_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Define parameters needed
CALL FinalizeParameters()
! Read Varnames to visualize and build calc and visu dependencies
CALL prms%SetSection("posti")
CALL prms%CreateIntOption(    "NVisu"           ,"Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateIntOption(    "NCalc"           ,"Polynomial degree at which calculations are done.")
CALL prms%CreateIntOption(    "TurbMode"        ,"Mode for turbulence quantities calculation:\n"//&
                                                 " 1 - RANS mode. Using lifted tau and S for epsilon.\n"//&
                                                 " 2 - LES mode.  Using lifted gradients for direct Reynolds stress calculation."//&
                                                 " 3 - Statefile mode (incomp.). Using one timeAvg file and several state files."//&
                                                 " 4 - Statefile mode (comp.). Using one timeAvg file and several state files.")
CALL prms%CreateLogicalOption("doConservativeDissipation", "Calculate turbulent dissipation by eps = eps_mean - eps_turb instead"//&
                                                 " of using the filtered velocity gradients.")
CALL prms%CreateIntFromStringOption('OutputFormat',"File format for visualization: None, Tecplot, TecplotASCII, ParaView. "//&
                                                 " Note: Tecplot output is currently unavailable due to licensing issues.", 'None')
CALL prms%CreateRealOption(   'RestartTKE'      ,"Constant TKE to be applied throughout the domain.", '0.')
CALL prms%CreateRealOption(   'Restartepsilon'  ,"Constant epsilon to be applied throughout the domain.", '0.')
CALL prms%CreateStringOption( 'NodeTypeVisu'    ,"NodeType for visualization. Visu, Gauss,Gauss-Lobatto,Visu_inner"    ,"VISU")
CALL addStrListEntry(         'OutputFormat'    ,'none',        OUTPUTFORMAT_NONE)
CALL addStrListEntry(         'OutputFormat'    ,'tecplot',     OUTPUTFORMAT_TECPLOT)
CALL addStrListEntry(         'OutputFormat'    ,'tecplotascii',OUTPUTFORMAT_TECPLOTASCII)
CALL addStrListEntry(         'OutputFormat'    ,'paraview',    OUTPUTFORMAT_PARAVIEW)
END SUBROUTINE DefineCalcTurb

!===================================================================================================================================
!> This routine is used to get everything ready for the first read of the data file.
!> - Find the mesh file
!> - Read the desired polynomial degree
!===================================================================================================================================
SUBROUTINE InitCalcTurb(ParameterFile,StateFile,ArrayName)
! MODULES
USE MOD_Globals,            ONLY: UNIT_stdOut,CollectiveStop,MPIROOT
USE MOD_PreProc,            ONLY: PP_N
USE MOD_CalcTurb_Vars
USE MOD_CalcTurb_Restart,   ONLY: GetDataProps
USE HDF5
USE MOD_HDF5_Input,         ONLY: ReadAttribute,File_ID,OpenDataFile,CloseDataFile,ReadArray,DatasetExists
USE MOD_HDF5_Input,         ONLY: GetArrayAndName
USE MOD_Interpolation_Vars, ONLY: NodeType
USE MOD_Mesh_Vars,          ONLY: MeshFile
USE MOD_Output_Vars,        ONLY: ProjectName,NOut
USE MOD_ReadInTools,        ONLY: GETINT,GETINTFROMSTR,GETLOGICAL,PRMS
USE MOD_StringTools,        ONLY: STRICMP,INTTOSTR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: ParameterFile
CHARACTER(LEN=*),INTENT(IN)   :: StateFile
CHARACTER(LEN=*),INTENT(IN)   :: ArrayName
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
!===================================================================================================================================
! open state file to be able to read attributes
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! read the meshfile attribute from statefile
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar = MeshFile)

! get properties from state file
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,ArrayName)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar = ProjectName)
CALL ReadAttribute(File_ID,'Time',        1,RealScalar= OutputTime)

! check if we are running on the current node type
IF (.NOT.STRICMP(NodeType_HDF5, NodeType)) THEN
    CALL CollectiveStop(__STAMP__, "NodeType of state does not match with NodeType the turbulence tool is compiled with!")
END IF

! read options from parameter file
CALL prms%read_options(ParameterFile)
! We might only need this on the first run
TurbMode                  = GETINT    ('TurbMode'                 ,'-1')
doConservativeDissipation = GETLOGICAL('doConservativeDissipation','.FALSE.')

! Polynomial degree for calculations. Set to N_HDF5 if not given otherwise
NCalc                 = GETINT("NCalc",INTTOSTR(N_HDF5))
IF (NCalc.LE.0) NCalc = N_HDF5
IF (NOut .LE.0) NOut  = PP_N

CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(132("="))')

END SUBROUTINE InitCalcTurb


!===================================================================================================================================
!> This routine deallocates everything used the calcturb tool.
!===================================================================================================================================
SUBROUTINE FinalizeCalcTurb
! MODULES
USE MOD_ReadInTools
USE MOD_CalcTurb_Vars
! External calls
USE MOD_DG,                     ONLY: FinalizeDG
#if USE_MPI
USE MOD_MPI,                    ONLY: FinalizeMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! External arrays
CALL FinalizeDG()
#if USE_MPI
CALL FinalizeMPI()
#endif

END SUBROUTINE FinalizeCalcTurb


!==================================================================================================================================
!> Determines if two real numbers are equal up to a specified tolerance (=PP_RealTolerance, normaly set to machine precision)
!> Takes into account that x,y are located in-between [-1;1] for additional accuracy
!> Based on Algorithm 139, Kopriva
!==================================================================================================================================
ELEMENTAL FUNCTION ALMOSTEQUAL(x,y)
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: x                !< (IN)  first scalar to be compared
REAL,INTENT(IN) :: y                !< (IN)  second scalar to be compared
LOGICAL         :: AlmostEqual      !< (OUT) TRUE if |x-y| < 2*PP_RealTolerance
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
AlmostEqual=.FALSE.
IF((x.EQ.0.).OR.(y.EQ.0.)) THEN
  IF(ABS(x-y).LE.2.*PP_RealTolerance) AlmostEqual=.TRUE.
ELSE ! x, y not zero
  IF((ABS(x-y).LE.PP_RealTolerance*ABS(x)).AND.((ABS(x-y).LE.PP_RealTolerance*ABS(y)))) AlmostEqual=.TRUE.
END IF ! x,y zero
END FUNCTION ALMOSTEQUAL

!===================================================================================================================================
! Performe an almost zero check. But ...
! Bruce Dawson quote:
! "There is no silver bullet. You have to choose wisely."
!    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
!      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
!      to your calculation. Maybe."
!    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
!      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
!      An absolute epsilon could be used if you knew exactly what number you were comparing against."
!    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
!      Good luck and God speed."
!===================================================================================================================================
FUNCTION ALMOSTZERO(Num)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL            :: Num ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL         :: ALMOSTZERO
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

ALMOSTZERO=.FALSE.
IF(ABS(Num).LE.epsilon(0.)) ALMOSTZERO=.TRUE.

END FUNCTION ALMOSTZERO


END MODULE MOD_CalcTurb
