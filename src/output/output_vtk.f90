!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
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

!===================================================================================================================================
!> Module for generic data output in vtk xml fromat
!> WARNING: WriteDataToVTK works only for POSTPROCESSING or for debug output during runtime
!===================================================================================================================================
MODULE MOD_VTK
USE ISO_C_BINDING
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE, BIND(C) :: CARRAY
  INTEGER (C_INT) :: len
  INTEGER (C_INT) :: dim
  TYPE (C_PTR)    :: data
END TYPE CARRAY

INTERFACE WriteDataToVTK
  MODULE PROCEDURE WriteDataToVTK
END INTERFACE

INTERFACE WriteVTKMultiBlockDataSet
  MODULE PROCEDURE WriteVTKMultiBlockDataSet
END INTERFACE

INTERFACE WriteCoordsToVTK_array
  MODULE PROCEDURE WriteCoordsToVTK_array
END INTERFACE

INTERFACE WriteDataToVTK_array
  MODULE PROCEDURE WriteDataToVTK_array
END INTERFACE

INTERFACE WriteVarnamesToVTK_array
  MODULE PROCEDURE WriteVarnamesToVTK_array
END INTERFACE

#if USE_PARTICLES
INTERFACE WritePartDataToVTK_array
  MODULE PROCEDURE WritePartDataToVTK_array
END INTERFACE

INTERFACE WriteDataToVTKPart
  MODULE PROCEDURE WriteDataToVTKPart
END INTERFACE
#endif

PUBLIC::WriteDataToVTK
PUBLIC::WriteVTKMultiBlockDataSet
PUBLIC::WriteCoordsToVTK_array
PUBLIC::WriteDataToVTK_array
PUBLIC::WriteVarnamesToVTK_array
#if USE_PARTICLES
PUBLIC::WriteDataToVTKPart
PUBLIC::WritePartDataToVTK_array
#endif
PUBLIC::CARRAY
!===================================================================================================================================

CONTAINS

SUBROUTINE CreateConnectivity(NVisu,nElems,nodeids,dim,DGFV)
USE ISO_C_BINDING
USE MOD_Globals
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                       :: NVisu
INTEGER,INTENT(IN)                       :: nElems
INTEGER,ALLOCATABLE,TARGET,INTENT(INOUT) :: nodeids(:)        !< stores the connectivity
INTEGER,INTENT(IN)                       :: dim               !< 3 = 3d connectivity, 2 = 2d connectivity
INTEGER,INTENT(IN)                       :: DGFV              !< flag indicating DG = 0 or FV = 1 data
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,iElem
INTEGER           :: NodeID,NodeIDElem
INTEGER           :: NVisu_k, NVisu_j, NVisu_elem, NVisu_p1_2
INTEGER           :: nVTKCells
!===================================================================================================================================
IF (dim.EQ.3) THEN
  NVisu_k = NVisu
  NVisu_j = NVisu
ELSE IF (dim.EQ.2) THEN
  NVisu_k = 1
  NVisu_j = NVisu
ELSE IF (dim.EQ.1) THEN
  NVisu_k = 1
  NVisu_j = 1
ELSE
  CALL Abort(__STAMP__, &
      "Only 2D and 3D connectivity can be created. dim must be 2 or 3.")
END IF

NVisu_elem = (NVisu+1)**dim
NVisu_p1_2 = (NVisu+1)**2

nVTKCells  = ((NVisu+DGFV)/(1+DGFV))**dim*nElems
SDEALLOCATE(nodeids)
ALLOCATE(nodeids((2**dim)*nVTKCells))


! create connectivity
NodeID = 0
NodeIDElem = 0
DO iElem=1,nElems
  DO k=1,NVisu_k,(DGFV+1)
    DO j=1,NVisu_j,(DGFV+1)
      DO i=1,NVisu,(DGFV+1)
        IF (dim.GE.2) THEN
          NodeID=NodeID+1
          nodeids(NodeID) = NodeIDElem+i+   j   *(NVisu+1)+(k-1)*NVisu_p1_2-1 !P4(CGNS=tecVisu standard)
        END IF
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+  (j-1)*(NVisu+1)+(k-1)*NVisu_p1_2-1 !P1
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+1+(j-1)*(NVisu+1)+(k-1)*NVisu_p1_2-1 !P2
        IF (dim.GE.2) THEN
          NodeID=NodeID+1
          nodeids(NodeID) = NodeIDElem+i+1+ j   *(NVisu+1)+(k-1)*NVisu_p1_2-1 !P3
        END IF
        IF (dim.EQ.3) THEN
          NodeID=NodeID+1
          nodeids(NodeID)=NodeIDElem+i+   j   *(NVisu+1)+ k   *NVisu_p1_2-1 !P8
          NodeID=NodeID+1
          nodeids(NodeID)=NodeIDElem+i+  (j-1)*(NVisu+1)+ k   *NVisu_p1_2-1 !P5
          NodeID=NodeID+1
          nodeids(NodeID)=NodeIDElem+i+1+(j-1)*(NVisu+1)+ k   *NVisu_p1_2-1 !P6
          NodeID=NodeID+1
          nodeids(NodeID)=NodeIDElem+i+1+ j   *(NVisu+1)+ k   *NVisu_p1_2-1 !P7
        END IF
      END DO
    END DO
  END DO
  NodeIDElem=NodeIDElem+NVisu_elem
END DO
END SUBROUTINE CreateConnectivity

!===================================================================================================================================
!> Subroutine to write 2D or 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteDataToVTK(nVal,NVisu,nElems,VarNames,Coord,Value,FileString,dim,DGFV,nValAtLastDimension)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nVal                 !< Number of nodal output variables
INTEGER,INTENT(IN)          :: NVisu                !< Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)          :: nElems               !< Number of output elements
INTEGER,INTENT(IN)          :: dim                  !< dimension: 2 or 3
REAL,INTENT(IN)             :: Coord(1:3,0:NVisu,0:NVisu*(MERGE(1,0,dim.GT.1)),0:NVisu*(MERGE(1,0,dim.GT.2)),nElems)     !< CoordsVector
CHARACTER(LEN=*),INTENT(IN) :: VarNames(nVal)       !< Names of all variables that will be written out
REAL,INTENT(IN)             :: Value(:,:,:,:,:)     !< Statevector
CHARACTER(LEN=*),INTENT(IN) :: FileString           !< Output file name
INTEGER,OPTIONAL,INTENT(IN) :: DGFV                 !< flag indicating DG = 0 or FV =1 data
LOGICAL,OPTIONAL,INTENT(IN) :: nValAtLastDimension  !< if TRUE, nVal is stored in the last index of value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iVal,ivtk
INTEGER                     :: nElems_glob(0:nProcessors-1)
INTEGER                     :: NVisu_elem,nVTKPoints,nVTKCells
INTEGER                     :: nTotalElems
INTEGER                     :: nBytes,Offset
INTEGER                     :: INTdummy
REAL(KIND=4)                :: FLOATdummy
CHARACTER(LEN=35)           :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
INTEGER                     :: ElemType,iElem
INTEGER,ALLOCATABLE,TARGET  :: nodeids(:)
INTEGER                     :: NVisu_k,NVisu_j,PointsPerVTKCell
#if USE_MPI
INTEGER                     :: iProc,nElems_proc,nElemsMax
REAL,ALLOCATABLE            :: buf(:,:,:,:), buf2(:,:,:,:,:)
#endif /*USE_MPI*/
INTEGER                     :: DGFV_loc
LOGICAL                     :: nValAtLastDimension_loc
!===================================================================================================================================
IF (PRESENT(DGFV)) THEN
  DGFV_loc = DGFV
ELSE
  DGFV_loc = 0
END IF
IF (PRESENT(nValAtLastDimension)) THEN
  nValAtLastDimension_loc = nValAtLastDimension
ELSE
  nValAtLastDimension_loc = .FALSE.
END IF

IF (dim.EQ.3) THEN
  NVisu_k = NVisu
  NVisu_j = NVisu
  PointsPerVTKCell = 8
ELSE IF (dim.EQ.2) THEN
  NVisu_k = 0
  NVisu_j = NVisu
  PointsPerVTKCell = 4
ELSE IF (dim.EQ.1) THEN
  NVisu_k = 0
  NVisu_j = 0
  PointsPerVTKCell = 2
ELSE
  CALL Abort(__STAMP__, &
      "Only 2D and 3D connectivity can be created. dim must be 1, 2 or 3.")
END IF

SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO')"   WRITE ",dim,"D DATA TO VTX XML BINARY (VTU) FILE..."

! get total number of elements on all processors
#if USE_MPI
CALL MPI_GATHER(nElems,1,MPI_INTEGER,nElems_glob,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
#else
nElems_glob(0) = nElems
#endif
nTotalElems = SUM(nElems_glob)

NVisu_elem = (NVisu+1)**dim
nVTKPoints = NVisu_elem * nTotalElems
nVTKCells  = ((NVisu+DGFV_loc)/(1+DGFV_loc))**dim*nTotalElems

! write header of VTK file
IF(MPIROOT)THEN
  ! Line feed character
  lf = char(10)

  ! Write file
  OPEN(NEWUNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  WRITE(TempStr1,'(I16)')nVTKPoints
  WRITE(TempStr2,'(I16)')nVTKCells
  Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
         &NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  DO iVal=1,nVal
    Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'" '// &
                     'format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    Offset=Offset+SIZEOF_F(INTdummy)+nVTKPoints*SIZEOF_F(FLOATdummy)
    WRITE(StrOffset,'(I16)')Offset
  END DO
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+3*nVTKPoints*SIZEOF_F(FLOATdummy)
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify necessary cell data
  Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Connectivity
  Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+PointsPerVTKCell*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Offsets
  Buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Elem types
  Buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)
END IF



#if USE_MPI
IF(MPIroot)THEN
  !ALLOCATE buffer for Root
  nElemsMax=MAXVAL(nElems_glob)
  ALLOCATE(buf(   0:NVisu,0:NVisu_j,0:NVisu_k,nElemsMax))
END IF
#endif

! Write binary raw data into append section
! Solution data
DO iVal=1,nVal
  IF(MPIroot)THEN
    nBytes = nVTKPoints*SIZEOF_F(FLOATdummy)
    IF (nValAtLastDimension_loc) THEN
      WRITE(ivtk) nBytes,REAL(Value(:,:,:,:,iVal),4)
    ELSE
      WRITE(ivtk) nBytes,REAL(Value(iVal,:,:,:,:),4)
    END IF
#if USE_MPI
    DO iProc=1,nProcessors-1
      nElems_proc=nElems_glob(iProc)
      IF (nElems_proc.GT.0) THEN
        CALL MPI_RECV(buf(:,:,:,1:nElems_proc),nElems_proc*NVisu_elem,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_FLEXI,MPIstatus,iError)
        WRITE(ivtk) REAL(buf(:,:,:,1:nElems_proc),4)
      END IF
    END DO !iProc
  ELSE
    IF (nElems.GT.0) THEN
      IF (nValAtLastDimension_loc) THEN
        CALL MPI_SEND(Value(:,:,:,:,iVal),nElems*NVisu_elem,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
      ELSE
        CALL MPI_SEND(Value(iVal,:,:,:,:),nElems*NVisu_elem,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
      END IF
    END IF
#endif /*USE_MPI*/
  END IF !MPIroot
END DO       ! iVar

#if USE_MPI
IF(MPIroot)THEN
  SDEALLOCATE(buf)
  ALLOCATE(buf2(3,0:NVisu,0:NVisu_j,0:NVisu_k,nElemsMax))
END IF
#endif

! Coordinates
IF(MPIRoot)THEN
  nBytes = nVTKPoints*SIZEOF_F(FLOATdummy) * 3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(Coord(:,:,:,:,:),4)
#if USE_MPI
  DO iProc=1,nProcessors-1
    nElems_proc=nElems_glob(iProc)
    IF (nElems_proc.GT.0) THEN
      CALL MPI_RECV(buf2(:,:,:,:,1:nElems_proc),nElems_proc*NVisu_elem*3,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_FLEXI,MPIstatus,iError)
      WRITE(ivtk) REAL(buf2(:,:,:,:,1:nElems_proc),4)
    END IF
  END DO !iProc
ELSE
  IF (nElems.GT.0) THEN
    CALL MPI_SEND(Coord(:,:,:,:,:),nElems*NVisu_elem*3,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
  END IF
#endif /*USE_MPI*/
END IF !MPIroot

#if USE_MPI
IF(MPIroot)THEN
  SDEALLOCATE(buf2)
END IF
#endif

! Connectivity and footer
IF(MPIROOT)THEN
  CALL CreateConnectivity(NVisu,nTotalElems,nodeids,dim,DGFV_loc)

  nBytes = PointsPerVTKCell*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) nodeids
  ! Offset
  nBytes = nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) (Offset,Offset=PointsPerVTKCell,PointsPerVTKCell*nVTKCells,PointsPerVTKCell)
  ! Elem type
  IF (dim.EQ.3) THEN
    ElemType = 12 ! VTK_HEXAHEDRON
  ELSE IF (dim.EQ.2) THEN
    ElemType = 9  ! VTK_QUAD
  ELSE IF (dim.EQ.1) THEN
    ElemType = 3  ! VTK_LINE
  END IF
  WRITE(ivtk) nBytes
  WRITE(ivtk) (ElemType,iElem=1,nVTKCells)

  DEALLOCATE(nodeids)

  ! Footer
  lf = char(10)
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
ENDIF
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK

!===================================================================================================================================
!> Links DG and FV VTK files together
!===================================================================================================================================
SUBROUTINE WriteVTKMultiBlockDataSet(FileString,FileString_DG,FileString_FV)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileString     !< Output file name
CHARACTER(LEN=*),INTENT(IN) :: FileString_DG  !< Filename of DG VTU file
CHARACTER(LEN=*),INTENT(IN) :: FileString_FV  !< Filename of FV VTU file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ivtk
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
!===================================================================================================================================
IF (MPIRoot) THEN
  ! write multiblock file
  OPEN(NEWUNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
  ! Line feed character
  lf = char(10)
  Buffer='<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf
  WRITE(ivtk) TRIM(BUFFER)
  Buffer='  <vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    <DataSet index="0" name="DG" file="'//TRIM(FileString_DG)//'">'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    </DataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    <DataSet index="1" name="FV" file="'//TRIM(FileString_FV)//'">'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    </DataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='  </vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(BUFFER)
  CLOSE(ivtk)
ENDIF
END SUBROUTINE WriteVTKMultiBlockDataSet

!===================================================================================================================================
!> Subroutine to write 2D or 3D coordinates to VTK format
!===================================================================================================================================
SUBROUTINE WriteCoordsToVTK_array(NVisu,nElems,coords_out,nodeids_out,coords,nodeids,dim,DGFV)
USE ISO_C_BINDING
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: NVisu                        !< Polynomial degree for visualization
INTEGER,INTENT(IN)                   :: nElems                       !< Number of elements
INTEGER,INTENT(IN)                   :: dim                          !< Spacial dimension (2D or 3D)
INTEGER,INTENT(IN)                   :: DGFV                         !< flag indicating DG = 0 or FV =1 data
REAL,ALLOCATABLE,TARGET,INTENT(IN)    :: coords(:,:,:,:,:) !< Array containing coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,ALLOCATABLE,TARGET,INTENT(INOUT) :: nodeids(:)
TYPE (CARRAY), INTENT(INOUT)         :: coords_out
TYPE (CARRAY), INTENT(INOUT)         :: nodeids_out

!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
coords_out%dim  = dim
IF (nElems.EQ.0) THEN
  coords_out%len  = 0
  nodeids_out%len = 0
  RETURN
END IF
SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO')"   WRITE ",dim,"D COORDS TO VTX XML BINARY (VTU) ARRAY..."
! values and coords are already in the correct structure of VTK/Paraview

! create connectivity
CALL CreateConnectivity(NVisu,nElems,nodeids,dim,DGFV)

! set the sizes of the arrays
coords_out%len = 3*(NVisu+1)**dim*nElems
nodeids_out%len = (2**dim)*((NVisu+DGFV)/(1+DGFV))**dim*nElems

! assign data to the arrays (no copy!!!)
coords_out%data = C_LOC(Coords(1,0,0,0,1))
nodeids_out%data = C_LOC(nodeids(1))

SWRITE(UNIT_stdOut,'(A)')" Done!"
END SUBROUTINE WriteCoordsToVTK_array

!===================================================================================================================================
!> Subroutine to write actual 2D or 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteDataToVTK_array(nVal,NVisu,nElems,Values_out,values,dim)
USE ISO_C_BINDING
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: nVal                         !> Number of nodal output variables
INTEGER,INTENT(IN)                :: NVisu                        !> Polynomial degree for visualization
INTEGER,INTENT(IN)                :: nElems                       !> Number of elements
INTEGER,INTENT(IN)                :: dim                          !> Spacial dimension (2D or 3D)
REAL(C_DOUBLE),ALLOCATABLE,TARGET,INTENT(IN) :: values(:,:,:,:,:) !> Array containing the points values
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE (CARRAY), INTENT(INOUT)      :: values_out

!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
values_out%dim  = dim
IF (nElems.EQ.0) THEN
  values_out%len  = 0
  RETURN
END IF
SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO')"   WRITE ",dim,"D DATA TO VTX XML BINARY (VTU) ARRAY..."

! values and coords are already in the correct structure of VTK/Paraview
! set the sizes of the arrays
values_out%len = nVal*(NVisu+1)**dim*nElems

IF (nVal*(NVisu+1)**dim*nElems.GT.0) THEN
  ! assign data to the arrays (no copy!!!)
  values_out%data = C_LOC(values(0,0,0,1,1))
END IF

SWRITE(UNIT_stdOut,'(A)')" Done!"
END SUBROUTINE WriteDataToVTK_array

!===================================================================================================================================
!> Subroutine to write variable names to VTK format
!===================================================================================================================================
SUBROUTINE WriteVarnamesToVTK_array(nVarTotal,mapVisu,varnames_out,VarNamesTotal,nVarVisu)
USE ISO_C_BINDING
! MODULES
USE MOD_Globals
USE MOD_StringTools    ,ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVarTotal
INTEGER,INTENT(IN)             :: mapVisu(nVarTotal)
TYPE (CARRAY), INTENT(INOUT)   :: varnames_out
CHARACTER(LEN=255),INTENT(IN)  :: VarNamesTotal(nVarTotal)
INTEGER,INTENT(IN)             :: nVarVisu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(C_CHAR),POINTER    :: VarNames_loc(:,:)
INTEGER                      :: i,iVar
!===================================================================================================================================
! copy varnames
ALLOCATE(VarNames_loc(255,nVarVisu))
varnames_out%len  = nVarVisu*255
IF (nVarVisu.GT.0) THEN
  varnames_out%data = C_LOC(VarNames_loc(1,1))

  DO iVar=1,nVarTotal
    IF (mapVisu(iVar).GT.0) THEN
      DO i=1,255
        VarNames_loc(i,mapVisu(iVar)) = VarNamesTotal(iVar)(i:i)
      END DO
    END IF
  END DO
END IF

END SUBROUTINE WriteVarnamesToVTK_array

#if USE_PARTICLES
SUBROUTINE WritePartDataToVTK_array(nParts_out,nVar_out,coords_out,values_out,nodeids_out,&
                                      varnamespart_out, componentspart_out, coords, values,nodeids,&
                                      VarNamePartCombine,VarNamePartCombineLen,VarNamePartVisu,PartCPointers_allocated) 
USE ISO_C_BINDING
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: nParts_out                    ! Number of nodal output variables
INTEGER,INTENT(IN)                :: nVar_out
REAL(C_DOUBLE),POINTER,INTENT(IN) :: coords(:,:)
REAL(C_DOUBLE),POINTER,INTENT(IN) :: values(:,:)
INTEGER(C_INT),POINTER,INTENT(IN) :: nodeids(:)
INTEGER,INTENT(IN)                :: VarNamePartCombine(:)
CHARACTER(LEN=255),INTENT(IN)     :: VarNamePartVisu(:)
INTEGER,INTENT(IN)                :: VarNamePartCombineLen(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE (CARRAY), INTENT(INOUT)      :: coords_out
TYPE (CARRAY), INTENT(INOUT)      :: values_out
TYPE (CARRAY), INTENT(INOUT)      :: nodeids_out
TYPE (CARRAY), INTENT(INOUT)      :: varnamespart_out
TYPE (CARRAY), INTENT(INOUT)      :: componentspart_out
LOGICAL,INTENT(INOUT)             :: PartCPointers_allocated 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iPart, nPartVarCombine, iVar, iVar2
INTEGER                           :: NodeID, str_len, i
CHARACTER(C_CHAR),POINTER         :: VarNamesPart_loc(:,:)
INTEGER(C_INT),POINTER            :: componentspart_loc(:)
!===================================================================================================================================
IF(nVar_out.EQ.0)THEN
  coords_out%len=0
  RETURN
END IF
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE PARTICEL DATA TO VTX XML BINARY (VTU) ARRAY..."

PartCPointers_allocated=.TRUE.
! values and coords are already in the correct structure of VTK/Paraview 
! set the sizes of the arrays
coords_out%len = 3*nParts_Out
values_out%len = nVar_out*nParts_out
nodeids_out%len = nParts_Out

! assign data to the arrays (no copy!!!)
coords_out%data = C_LOC(Coords(1,1))
values_out%data = C_LOC(values(1,1))
nodeids_out%data = C_LOC(nodeids(1))

! copy varnames
nPartVarCombine = 0
DO iVar=1,nVar_out
  IF (VarNamePartCombine(iVar).LE.1) THEN
    nPartVarCombine = nPartVarCombine + 1
  END IF
END DO

ALLOCATE(VarNamesPart_loc(255,nPartVarCombine))
varnamespart_out%len = nPartVarCombine*255
varnamespart_out%data = C_LOC(VarNamesPart_loc(1,1))

ALLOCATE(componentspart_loc(nPartVarCombine))
componentspart_out%len = nPartVarCombine
componentspart_out%data = C_LOC(componentspart_loc(1))

iVar2 = 1
DO iVar=1,nVar_out
  IF (VarNamePartCombine(iVar).EQ.0) THEN
    DO i=1,255
      VarNamesPart_loc(i,iVar2) = VarNamePartVisu(iVar)(i:i)
    END DO
    componentspart_loc(iVar2) = 1
    iVar2 = iVar2 + 1
  ELSE IF (VarNamePartCombine(iVar).EQ.1) THEN
    str_len = LEN_TRIM(VarNamePartVisu(iVar))
    DO i=1,255
      VarNamesPart_loc(i,iVar2) = VarNamePartVisu(iVar)(i:i)
    END DO
    VarNamesPart_loc(str_len,iVar2) =  ' '
    componentspart_loc(iVar2) = VarNamePartCombineLen(iVar)
    iVar2 = iVar2 + 1
  END IF
END DO

! create connectivity
NodeID = 0
DO iPart=1,nParts_out
  NodeID=NodeID+1
  nodeids(NodeID)=NodeID
END DO ! iPart=1,nParts_out

SWRITE(UNIT_stdOut,'(A)')" Done!"
END SUBROUTINE WritePartDataToVTK_array

SUBROUTINE WriteDataToVTKPart(nParts,nVal,Coord,Value,FileString,VarNamePartVisu,VarNamePartCombine,VarNamePartCombineLen,&
    nGlobalParts)
!===================================================================================================================================
! Subroutine to write unstructured 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
INTEGER,INTENT(IN)            :: nParts                  ! Number of output Particle
REAL,INTENT(IN)               :: Coord(3,nParts)         ! CoordsVector 
REAL,INTENT(IN)               :: Value(nVal,nParts)      ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
CHARACTER(LEN=255),INTENT(IN) :: VarNamePartVisu(:)
INTEGER,INTENT(IN)            :: VarNamePartCombine(:)
INTEGER,INTENT(IN)            :: VarNamePartCombineLen(:)
INTEGER,INTENT(IN)            :: nGlobalParts 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iPart,Offset,nBytes,nVTKElems,nVTKCells,ivtk=44,iVar,str_len
INTEGER            :: INTdummy
INTEGER            :: Vertex(nGlobalParts)
INTEGER            :: ElemType
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf,components_string
CHARACTER(LEN=255) :: VarNameString
REAL(KIND=4)       :: FLOATdummy
REAL,ALLOCATABLE   :: buf(:), buf2(:,:)
INTEGER            :: nParts_glob(1:nProcessors), nParts_proc
INTEGER            :: iProc, nMaxParts
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE PART/EROSION DATA TO VTX XML BINARY (VTU) FILE..."

! collect number of particles on each proc
CALL MPI_ALLGATHER(nParts,1,MPI_INTEGER,nParts_glob(:),1,MPI_INTEGER,MPI_COMM_FLEXI,iError)

IF(nParts.LT.1)THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
  RETURN
END IF

nVTKElems=nGlobalParts
nVTKCells=nGlobalParts

! Write file
IF(MPIRoot)THEN
  ! Line feed character
  lf = char(10)
 
  OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
   
  Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  WRITE(TempStr1,'(I16)')nVTKElems
  WRITE(TempStr2,'(I16)')nVTKCells
  Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//&
  '" NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  IF (nVal .GT.0)THEN
    DO iVar=1,nVal
      IF (VarNamePartCombine(iVar).EQ.0) THEN
        Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNamePartVisu(iVar))//&
        '" NumberOfComponents="1" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
        Offset=Offset+SIZEOF_F(INTdummy)+nVTKElems*SIZEOF_F(FLOATdummy)
        WRITE(StrOffset,'(I16)')Offset
      ELSE IF (VarNamePartCombine(iVar).EQ.1) THEN
        str_len = LEN_TRIM(VarNamePartVisu(iVar))
        WRITE(components_string,'(I1)') VarNamePartCombineLen(iVar)
        VarNameString = VarNamePartVisu(iVar)(1:str_len-1)
        Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNameString)//&
        '" NumberOfComponents="'//components_string//'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf
        WRITE(ivtk) TRIM(Buffer)
        Offset=Offset+SIZEOF_F(INTdummy)+nVTKElems*SIZEOF_F(FLOATdummy)*VarNamePartCombineLen(iVar)
        WRITE(StrOffset,'(I16)')Offset
      END IF
    END DO
  END IF
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended"'// &
         ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+3*nVTKElems*SIZEOF_F(FLOATdummy)
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify necessary cell data
  Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Connectivity
  Buffer='        <DataArray type="Int32" Name="connectivity" format="appended"'// &
         ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Offsets
  Buffer='        <DataArray type="Int32" Name="offsets" format="appended"'// &
         ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Elem types
  Buffer='        <DataArray type="Int32" Name="types" format="appended"'// &
         ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)
END IF

#if USE_MPI
IF(MPIroot)THEN
  !ALLOCATE buffer for Root
  nMaxParts=MAXVAL(nParts_glob)
  ALLOCATE(buf(nMaxParts))
END IF
#endif

! Write binary raw data into append section
! Point data
DO iVar=1,nVal
  IF(MPIRoot)THEN
    nBytes = nVTKElems*SIZEOF_F(FLOATdummy)
    IF (VarNamePartCombine(iVar).EQ.0) THEN
      WRITE(ivtk) nBytes,REAL(Value(iVar,:),4)
    ELSEIF(VarNamePartCombine(iVar).EQ.1) THEN
      WRITE(ivtk) nBytes*VarNamePartCombineLen(iVar),REAL(Value(iVar:iVar+VarNamePartCombineLen(iVar)-1,:),4)
    END IF
#if USE_MPI
    DO iProc=1,nProcessors-1
      nParts_proc=nParts_glob(iProc+1)
      IF (nParts_proc.GT.0) THEN
        CALL MPI_RECV(buf(1:nParts_proc),nParts_proc,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_FLEXI,MPIstatus,iError)
        WRITE(ivtk) REAL(buf(1:nParts_proc),4)
      END IF
    END DO !iProc
  ELSE
    IF (nParts.GT.0) THEN
      CALL MPI_SEND(Value(iVar,:),nParts,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
    END IF
#endif /*USE_MPI*/
  END IF
END DO

#if USE_MPI
IF(MPIroot)THEN
  SDEALLOCATE(buf)
  nMaxParts=MAXVAL(nParts_glob)
  ALLOCATE(buf2(3,nMaxParts))
END IF
#endif

! Coordinates
IF(MPIRoot)THEN
  nBytes = nVTKElems*SIZEOF_F(FLOATdummy)*3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(Coord(:,:),4)
#if USE_MPI
  DO iProc=1,nProcessors-1
    nParts_proc=nParts_glob(iProc+1)
    IF (nParts_proc.GT.0) THEN
      CALL MPI_RECV(buf2(:,1:nParts_proc),nParts_proc*3,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_FLEXI,MPIstatus,iError)
      WRITE(ivtk) REAL(buf2(:,1:nParts_proc),4)
    END IF
  END DO !iProc
ELSE
  IF (nParts.GT.0) THEN
    CALL MPI_SEND(Coord(:,:),nParts*3,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
  END IF
#endif /*USE_MPI*/
END IF !MPIroot

#if USE_MPI
IF(MPIroot)THEN
  SDEALLOCATE(buf2)
END IF
#endif

IF(MPIRoot)THEN
  ! Connectivity
  DO iPart=1,nGlobalParts
    Vertex(iPart)=iPart-1
  END DO
  nBytes = nVTKElems*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) Vertex(:)
  ! Offset
  nBytes = nVTKElems*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) (Offset,Offset=1,nVTKElems,1)
  ! Elem type
  ElemType = 2  ! VTK_VERTEX (POINT)
  !ElemType = 12 ! VTK_HEXAHEDRON
  WRITE(ivtk) nBytes
  WRITE(ivtk) (ElemType,iPart=1,nVTKElems)
  ! Write footer
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END IF
END SUBROUTINE WriteDataToVTKPart
#endif

END MODULE MOD_VTK
