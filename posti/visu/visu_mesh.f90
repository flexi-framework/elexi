!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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

!=================================================================================================================================
!> Routines to build the mesh for visualization.
!=================================================================================================================================
MODULE MOD_Posti_ReadMesh
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE VisualizeMesh
  MODULE PROCEDURE VisualizeMesh
END INTERFACE

PUBLIC:: VisualizeMesh
!=================================================================================================================================

CONTAINS

!=================================================================================================================================
!> Visualize mesh only
!> 1. read mesh
!> 2. BuildVisuCoords
!> 3. Convert scaled jacobian
!> 4. write mesh to VTK array
!> 5. set length of all other output arrays to zero
!=================================================================================================================================
SUBROUTINE VisualizeMesh(postifile,meshfile_in,useCurveds)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_HDF5_Input    ,ONLY: ReadAttribute,ReadArray,File_ID,OpenDataFile,CloseDataFile
USE MOD_Mesh          ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars     ,ONLY: nElems,offsetElem,Ngeo,scaledJac
USE MOD_Interpolation ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG
USE MOD_Posti_VisuMesh      ,ONLY: BuildVisuCoords
USE MOD_ReadInTools   ,ONLY: prms,GETINT,GETSTR,GETLOGICAL,CountOption
USE MOD_ReadInTools   ,ONLY: FinalizeParameters
USE MOD_StringTools   ,ONLY: STRICMP
USE MOD_VTK           ,ONLY: WriteCoordsToVTK_array
#if USE_MPI
USE MOD_MPI           ,ONLY: FinalizeMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN):: postifile
CHARACTER(LEN=255),INTENT(IN):: meshfile_in
LOGICAL,INTENT(IN),OPTIONAL  :: UseCurveds
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iVar,jVar,iVarVisu,meshModeLoc
CHARACTER(LEN=255)  :: VarName
!===================================================================================================================================
#if USE_MPI
CALL FinalizeMPI()
#endif
CALL FinalizeMesh()
CALL FinalizeInterpolation()

CALL OpenDataFile(meshfile_in,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=Ngeo)
CALL CloseDataFile()

IF (LEN_TRIM(postifile).GT.0) THEN
  ! read options from parameter file
  CALL DefineParametersInterpolation()
  CALL DefineParametersMesh()
  CALL prms%SetSection("posti")
  CALL prms%CreateIntOption('NVisu', "Number of points at which solution is sampled for visualization.")
  CALL prms%read_options(postifile)
  NVisu     = GETINT('NVisu','1') ! Degree of visualization basis
  HighOrder = GETLOGICAL('HighOrder')
  ! Get number of variables to be visualized
  nVarIni = CountOption("VarName")
ELSE
  NVisu = 2*NGeo ! TODO: correct?
END IF
NVisu_FV = 1

! read mesh, depending if we should visualize the Jacobian or not different mesh modes are needed (calculate metrics or not)
meshModeLoc = 0
IF (nVarIni.GT.0) meshModeLoc=2
CALL InitInterpolation(NVisu)
IF (PRESENT(UseCurveds)) THEN; CALL InitMesh(meshMode=meshModeLoc,MeshFile_IN=MeshFile_in,UseCurveds_IN=UseCurveds)
ELSE;                          CALL InitMesh(meshMode=meshModeLoc,MeshFile_IN=MeshFile_in)
END IF

! convert to visu grid
nElems_DG = nElems
nElems_FV = 0
SDEALLOCATE(mapDGElemsToAllElems)
ALLOCATE(mapDGElemsToAllElems(nElems))
DO iElem=1,nElems
  mapDGElemsToAllElems(iElem) = iElem
END DO
CALL BuildVisuCoords()
DEALLOCATE(mapDGElemsToAllElems)

! Do we need to visualize the scaled Jacobian, or the max scaled Jacobian?
IF (nVarIni.GT.0) THEN
  ! A very simple mapping is build: There are two depending variables, either one or both of them can be visualized
  NCalc = PP_N
  nVarVisu = nVarIni
  IF (IJK_exists) THEN
    nVarDep = 5
    nVarAll = 5
  ELSE
    nVarDep = 2
    nVarAll = 2
  END IF
  SDEALLOCATE(mapDepToCalc)
  SDEALLOCATE(mapAllVarsToVisuVars)
  SDEALLOCATE(mapAllVarsToSurfVisuVars)
  ALLOCATE(mapDepToCalc(nVarDep))
  DO iVar = 1, nVarDep
    mapDepToCalc(iVar) = iVar
  END DO ! iVar
  ALLOCATE(mapAllVarsToVisuVars(nVarAll))
  mapAllVarsToVisuVars = 0
  ALLOCATE(mapAllVarsToSurfVisuVars(1:nVarAll))
  mapAllVarsToSurfVisuVars = 0
  iVarVisu = 1
  DO iVar = 1, nVarIni
    VarName = GETSTR("VarName")
    DO jVar = 1, nVarAll
      IF (STRICMP(VarNamesAll(jVar),VarName)) THEN
        mapAllVarsToVisuVars(jVar) = iVarVisu
        iVarVisu = iVarVisu + 1
      END IF
    END DO ! jVar = 1, nVarAll
  END DO ! iVar = 1, nVarIni
  SDEALLOCATE(UCalc_DG)
  ALLOCATE(UCalc_DG(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems_DG,nVarDep))

  ! Fill the array
  UCalc_DG(:,:,:,:,1) = scaledJac
  DO iElem=1,nElems
    UCalc_DG(:,:,:,iElem,2) = MINVAL(UCalc_DG(:,:,:,iElem,1))
  END DO ! iElem

  IF (IJK_exists) THEN
    ALLOCATE(Elem_IJK(3,nElems))
    CALL OpenDataFile(meshfile_in,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL ReadArray('Elem_IJK',2,(/3,nElems/),offsetElem,2,IntArray=Elem_IJK)
    DO iElem=1,nElems
      UCalc_DG(:,:,:,iElem,3) = Elem_IJK(1,iElem)
      UCalc_DG(:,:,:,iElem,4) = Elem_IJK(2,iElem)
      UCalc_DG(:,:,:,iElem,5) = Elem_IJK(3,iElem)
    END DO
    DEALLOCATE(Elem_IJK)
    CALL CloseDataFile()
  END IF

  CALL ConvertToVisu_DG()
ELSE
  nVarVisu = 0
END IF

CALL FinalizeInterpolation()
CALL FinalizeParameters()
END SUBROUTINE VisualizeMesh

END MODULE MOD_Posti_ReadMesh
