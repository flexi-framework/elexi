!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
!> This module contains all the routines that provide the interfaces between the FORTRAN visu tool and the C ParaView plugin.
!> These routines include:
!> * RequestInformation: Called by ParaView when a new state is loaded into the pipeline. Will return the available variable and
!>   boundary names to display them in ParaView for the user to choose from.
!> * visuCwrapper: Called by ParaView when data is requested after he apply button has been pressed. In this routine, the actual
!>   visu main routine is called with the parameter file created by the ParaView reader (based on the settings in ParaView choosen
!>   by the user). After the visu routine, the data and coordinate arrays are converted to a C pointer that is passed to ParaView.
!===================================================================================================================================
MODULE MOD_Visu_Cwrapper
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: visu_requestInformation
PUBLIC:: visu_CWrapper
PUBLIC:: visu_dealloc_nodeids
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Function to convert a C string with length strlen to a FORTRAN character array with length 255.
!===================================================================================================================================
FUNCTION cstrToChar255(cstr, strlen)
! MODULES
USE ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
TYPE(C_PTR),TARGET,INTENT(IN)  :: cstr
INTEGER,INTENT(IN)             :: strlen
CHARACTER(LEN=255)             :: cstrToChar255
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(KIND=C_CHAR),POINTER :: tmp(:)
!===================================================================================================================================
CALL C_F_POINTER(C_LOC(cstr), tmp, [strlen])
cstrToChar255 = TRANSFER(tmp(1:strlen), cstrToChar255)
cstrToChar255(strlen+1:255) = ' '
END FUNCTION cstrToChar255


!===================================================================================================================================
!> Wrapper to visu_InitFile for Paraview plugin, returns the available variable names and boundary names.
!===================================================================================================================================
SUBROUTINE visu_requestInformation(mpi_comm_IN, strlen_state, statefile_IN, strlen_mesh, meshfile_IN, varnames, bcnames &
#if USE_PARTICLES
                                  ,partnames &
#endif /*USE_PARTICLES*/
                                  )
USE ISO_C_BINDING
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5    ,ONLY: InitMPIInfo
USE MOD_MPI        ,ONLY: InitMPI
USE MOD_Visu_Init  ,ONLY: visu_getVarNamesAndFileType
USE MOD_Visu_Vars  ,ONLY: VarnamesAll,BCNamesAll,nVarIni
USE MOD_VTK        ,ONLY: CARRAY
#if USE_PARTICLES
USE MOD_Visu_Vars  ,ONLY: PartNamesAll
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                    :: mpi_comm_IN
INTEGER,INTENT(IN)                    :: strlen_state
TYPE(C_PTR),TARGET,INTENT(IN)         :: statefile_IN
INTEGER,INTENT(IN)                    :: strlen_mesh
TYPE(C_PTR),TARGET,INTENT(IN)         :: meshfile_IN
TYPE (CARRAY), INTENT(INOUT)          :: varnames
TYPE (CARRAY), INTENT(INOUT)          :: bcnames
#if USE_PARTICLES
TYPE (CARRAY), INTENT(INOUT)          :: partnames
#endif /*USE_PARTICLES*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
CHARACTER(LEN=255)                    :: statefile
CHARACTER(LEN=255)                    :: meshfile
CHARACTER(LEN=255),POINTER            :: varnames_pointer(:)
CHARACTER(LEN=255),POINTER            :: bcnames_pointer(:)
#if USE_PARTICLES
CHARACTER(LEN=255),POINTER            :: partnames_pointer(:)
#endif
!===================================================================================================================================
statefile = cstrToChar255(statefile_IN, strlen_state)
meshfile  = cstrToChar255(meshfile_IN , strlen_mesh)

! Set dummy value to force visu to read the actual values
nVarIni = -1

CALL InitMPIInfo()
CALL InitMPI( &
#if USE_MPI
             mpi_comm_IN &
#endif /*USE_MPI*/
            )
CALL visu_getVarNamesAndFileType(statefile,meshfile,VarnamesAll,BCNamesAll)
IF (ALLOCATED(VarnamesAll)) THEN
  varnames_pointer => VarnamesAll
  varnames%len  = SIZE(varnames_pointer)*255
  varnames%data = C_LOC(varnames_pointer(1))
ELSE
  varnames%len  = 0
  varnames%data = C_NULL_PTR
END IF
IF (ALLOCATED(BCNamesAll)) THEN
  bcnames_pointer => BCNamesAll
  bcnames%len  = SIZE(bcnames_pointer)*255
  bcnames%data = C_LOC(bcnames_pointer(1))
ELSE
  bcnames%len  = 0
  bcnames%data = C_NULL_PTR
END IF
#if USE_PARTICLES
IF (ALLOCATED(PartNamesAll)) THEN
  partnames_pointer => PartNamesAll
  partnames%len  = SIZE(partnames_pointer)*255
  partnames%data = C_LOC(partnames_pointer(1))
ELSE
  partnames%len  = 0
  partnames%data = C_NULL_PTR
END IF
#endif
END SUBROUTINE visu_requestInformation


!===================================================================================================================================
!> C wrapper routine for the visu call from ParaView. The main visu routine is called with the parameter file created by the
!> ParaView reader, and afterwards the data and coordinate arrays as well as the variable names are converted to C arrays since
!> ParaView needs the data in this format.
!===================================================================================================================================
SUBROUTINE visu_CWrapper(mpi_comm_IN,  &
    UseHighOrder,UseCurveds_IN,                                                                                     &
    strlen_prm      ,prmfile_IN      ,strlen_posti     ,postifile_IN       ,strlen_state ,statefile_IN,             &
    coordsDG_out    ,valuesDG_out    ,nodeidsDG_out    ,&
    coordsFV_out    ,valuesFV_out    ,nodeidsFV_out    ,&
    varnames_out,                                                                                                   &
    coordsSurfDG_out,valuesSurfDG_out,nodeidsSurfDG_out,&
    coordsSurfFV_out,valuesSurfFV_out,nodeidsSurfFV_out,&
    varnamesSurf_out                                    &
#if USE_PARTICLES
   ,coordsPart_out  ,valuesPart_out  ,nodeidsPart_out  ,varnamesPart_out  ,componentsPart_out,                      &
    coordsImpact_out,valuesImpact_out,nodeidsImpact_out,varnamesImpact_out,componentsImpact_out                     &
#endif /*USE_PARTICLES*/
    )
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_Output_Vars ,ONLY: doPrintStatusLine
USE MOD_Visu_Vars
USE MOD_Visu        ,ONLY: visu
USE MOD_VTK         ,ONLY: WriteCoordsToVTK_array,WriteDataToVTK_array,WriteVarnamesToVTK_array,CARRAY
#if USE_PARTICLES
USE MOD_VTK         ,ONLY: WritePartDataToVTK_array
#endif /*USE_PARTICLES*/
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: mpi_comm_IN
INTEGER,INTENT(IN)            :: UseHighOrder
INTEGER,INTENT(IN)            :: UseCurveds_IN
INTEGER,INTENT(IN)            :: strlen_prm
INTEGER,INTENT(IN)            :: strlen_posti
INTEGER,INTENT(IN)            :: strlen_state
TYPE(C_PTR),TARGET,INTENT(IN) :: prmfile_IN
TYPE(C_PTR),TARGET,INTENT(IN) :: postifile_IN
TYPE(C_PTR),TARGET,INTENT(IN) :: statefile_IN
TYPE (CARRAY), INTENT(INOUT)  :: coordsDG_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesDG_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsDG_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsFV_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesFV_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsFV_out
TYPE (CARRAY), INTENT(INOUT)  :: varnames_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsSurfDG_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesSurfDG_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsSurfDG_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsSurfFV_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesSurfFV_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsSurfFV_out
TYPE (CARRAY), INTENT(INOUT)  :: varnamesSurf_out
#if USE_PARTICLES
TYPE (CARRAY), INTENT(INOUT)  :: coordsPart_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesPart_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsPart_out
TYPE (CARRAY), INTENT(INOUT)  :: varnamesPart_out
TYPE (CARRAY), INTENT(INOUT)  :: componentsPart_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsImpact_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesImpact_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsImpact_out
TYPE (CARRAY), INTENT(INOUT)  :: varnamesImpact_out
TYPE (CARRAY), INTENT(INOUT)  :: componentsImpact_out
#endif /*USE_PARTICLES*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                      :: prmfile
CHARACTER(LEN=255)                      :: postifile
CHARACTER(LEN=255)                      :: statefile
LOGICAl                                 :: UseCurveds
!===================================================================================================================================
prmfile    = cstrToChar255(prmfile_IN,   strlen_prm)
postifile  = cstrToChar255(postifile_IN, strlen_posti)
statefile  = cstrToChar255(statefile_IN, strlen_state)
UseCurveds = MERGE(.TRUE.,.FALSE.,UseCurveds_IN.GT.0)

! Enable progress indicator
doPrintStatusLine = .TRUE.

CALL visu(mpi_comm_IN, prmfile, postifile, statefile, UseCurveds)

! Map Fortran arrays to C pointer
IF (MeshFileMode) THEN
  ! Write only the DG coordinates to the VTK file
  CALL WriteCoordsToVTK_array(NVisu,nElems_DG,coordsDG_out,nodeidsDG_out,CoordsVisu_DG,nodeids_DG,dim=PP_dim,DGFV=0,HighOrder=UseHighOrder)
  ! We may visualize the scaled Jacobian for debug purposes
  IF (nVarVisu.GT.0) THEN
    CALL WriteDataToVTK_array(nVarVisu,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG,PP_dim)
    CALL WriteVarnamesToVTK_array(nVarAll,mapAllVarsToVisuVars,varnames_out,VarnamesAll,nVarVisu)
  ELSE
    ! Otherwise, no output of data
    valuesDG_out%len      = 0
    varnames_out%len      = 0
  END IF

  ! set length of all other output arrays to zero so they are not used in the reader
  coordsFV_out%dim      = 3
  coordsFV_out%len      = 0
  valuesFV_out%len      = 0
  nodeidsFV_out%len     = 0
  coordsSurfDG_out%dim  = 2
  coordsSurfDG_out%len  = 0
  valuesSurfDG_out%len  = 0
  nodeidsSurfDG_out%len = 0
  coordsSurfFV_out%dim  = 2
  coordsSurfFV_out%len  = 0
  valuesSurfFV_out%len  = 0
  nodeidsSurfFV_out%len = 0
  varnamesSurf_out%len  = 0
#if USE_PARTICLES
  coordsPart_out%len    = 0
  coordsPart_out%dim    = 0
  valuesPart_out%len    = 0
  nodeidsPart_out%len   = 0
  varnamesPart_out%len  = 0
  componentsPart_out%len  = 0
  coordsImpact_out%len    = 0
  coordsImpact_out%dim    = 0
  valuesImpact_out%len    = 0
  nodeidsImpact_out%len   = 0
  varnamesImpact_out%len  = 0
  componentsImpact_out%len  = 0
#endif /*USE_PARTICLES*/

  RETURN
END IF

! write UVisu to VTK 2D / 3D arrays (must be done always!)
! write coords, UVisu to VTK  2D / 3D arrays (must be done always!)
IF (Avg2D) THEN
  CALL WriteDataToVTK_array(nVarVisu,NVisu   ,nElemsAvg2D_DG,valuesDG_out,UVisu_DG,2)
  CALL WriteDataToVTK_array(nVarVisu,NVisu_FV,nElemsAvg2D_FV,valuesFV_out,UVisu_FV,2)
  CALL WriteCoordsToVTK_array(NVisu,nElemsAvg2D_DG,coordsDG_out,nodeidsDG_out,CoordsVisu_DG,nodeids_DG,dim=2,DGFV=0,HighOrder=UseHighOrder)
  CALL WriteCoordsToVTK_array(NVisu_FV,nElemsAvg2D_FV,coordsFV_out,nodeidsFV_out,CoordsVisu_FV,nodeids_FV,dim=2,DGFV=1,HighOrder=UseHighOrder)
ELSE
  CALL WriteDataToVTK_array(nVarVisu,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG,PP_dim)
  CALL WriteDataToVTK_array(nVarVisu,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV,PP_dim)
  CALL WriteCoordsToVTK_array(NVisu,nElems_DG,coordsDG_out,nodeidsDG_out,CoordsVisu_DG,nodeids_DG,dim=PP_dim,DGFV=0,HighOrder=UseHighOrder)
  CALL WriteCoordsToVTK_array(NVisu_FV,nElems_FV,coordsFV_out,nodeidsFV_out,CoordsVisu_FV,nodeids_FV,dim=PP_dim,DGFV=1,HighOrder=UseHighOrder)
END IF

#if USE_PARTICLES
ALLOCATE(PD%PartIds_Visu (                   1:PD%nPart_visu) &
        ,PD%Part_Pos_visu(1:3               ,1:PD%nPart_Visu) &
        ,PD%Part_visu    (1:PD%nPartVar_visu,1:PD%nPart_visu))
IF (ALLOCATED(PD%PartData_HDF5)) THEN
  PD%Part_Pos_visu = PD%PartData_HDF5(1:3,:)
  PD%Part_visu     = PD%PartData_HDF5(4:,:)
END IF

CALL WritePartDataToVTK_array(PD%nPart_visu,PD%nPartVar_visu,coordsPart_out,valuesPart_out,nodeidsPart_out,varnamesPart_out,&
                              componentsPart_out,PD%Part_Pos_visu,PD%Part_visu,PD%PartIds_Visu,PD%VarNamePartCombine,&
                              PD%VarNamePartCombineLen,PD%VarNamePartVisu,PD%PartCPointers_allocated)

ALLOCATE(PDE%PartIds_Visu (                    1:PDE%nPart_visu) &
        ,PDE%Part_Pos_visu(1:3                ,1:PDE%nPart_Visu) &
        ,PDE%Part_visu    (1:PDE%nPartVar_visu,1:PDE%nPart_visu))
IF(ALLOCATED(PDE%PartData_HDF5))THEN
  PDE%Part_Pos_visu = PDE%PartData_HDF5(1:3,:)
  PDE%Part_visu     = PDE%PartData_HDF5(4:,:)
END IF

CALL WritePartDataToVTK_array(PDE%nPart_visu,PDE%nPartVar_visu,coordsImpact_out,valuesImpact_out,nodeidsImpact_out,varnamesImpact_out,&
                              componentsImpact_out,PDE%Part_Pos_visu,PDE%Part_visu,PDE%PartIds_Visu,PDE%VarNamePartCombine,&
                              PDE%VarNamePartCombineLen,PDE%VarNamePartVisu,PDE%PartCPointers_allocated)
#endif

CALL WriteVarnamesToVTK_array(nVarAll,mapAllVarsToVisuVars,varnames_out,VarnamesAll,nVarVisu)

! Surface
CALL WriteDataToVTK_array(nVarSurfVisuAll,NVisu   ,nBCSidesVisu_DG,valuesSurfDG_out,USurfVisu_DG,PP_dim-1)
CALL WriteDataToVTK_array(nVarSurfVisuAll,NVisu_FV,nBCSidesVisu_FV,valuesSurfFV_out,USurfVisu_FV,PP_dim-1)

CALL WriteCoordsToVTK_array(NVisu   ,nBCSidesVisu_DG,coordsSurfDG_out,nodeidsSurfDG_out,&
    CoordsSurfVisu_DG,nodeidsSurf_DG,dim=PP_dim-1,DGFV=0,HighOrder=UseHighOrder)
CALL WriteCoordsToVTK_array(NVisu_FV,nBCSidesVisu_FV,coordsSurfFV_out,nodeidsSurfFV_out,&
    CoordsSurfVisu_FV,nodeidsSurf_FV,dim=PP_dim-1,DGFV=1,HighOrder=UseHighOrder)

CALL WriteVarnamesToVTK_array(nVarAll,mapAllVarsToSurfVisuVars,varnamesSurf_out,VarnamesAll,nVarSurfVisuAll)

END SUBROUTINE visu_CWrapper


!===================================================================================================================================
!> Deallocate the different NodeID arrays.
!===================================================================================================================================
SUBROUTINE visu_dealloc_nodeids()
! MODULES
USE MOD_Visu_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! volume
SDEALLOCATE(nodeids_DG)
SDEALLOCATE(nodeids_FV)
! surface
SDEALLOCATE(nodeidsSurf_DG)
SDEALLOCATE(nodeidsSurf_FV)
! particles
#if USE_PARTICLES
IF (PD%PartCPointers_allocated) THEN
  DEALLOCATE(PD%Part_visu)
  DEALLOCATE(PD%Part_Pos_visu)
  DEALLOCATE(PD%PartIDs_visu)
  PD%PartCPointers_allocated=.FALSE.
END IF
IF (PDE%PartCPointers_allocated) THEN
  DEALLOCATE(PDE%Part_visu)
  DEALLOCATE(PDE%Part_Pos_visu)
  DEALLOCATE(PDE%PartIDs_visu)
  PDE%PartCPointers_allocated=.FALSE.
END IF
#endif

END SUBROUTINE visu_dealloc_nodeids

END MODULE MOD_Visu_Cwrapper
