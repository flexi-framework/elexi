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
! Calculate GCC version
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

!=================================================================================================================================
!> Routines to build the mesh for visualization.
!=================================================================================================================================
MODULE MOD_Posti_VisuMesh
#if !FV_ENABLED
USE ISO_C_BINDING
#endif
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if !FV_ENABLED
TYPE, BIND(C) :: CARRAY
  INTEGER (C_INT) :: len
  INTEGER (C_INT) :: dim
  TYPE (C_PTR)    :: data
END TYPE CARRAY

INTEGER,PARAMETER :: SHA256_LENGTH=64
#endif

INTERFACE BuildVisuCoords
  MODULE PROCEDURE BuildVisuCoords
END INTERFACE

INTERFACE BuildSurfVisuCoords
  MODULE PROCEDURE BuildSurfVisuCoords
END INTERFACE

INTERFACE VisualizeMesh
  MODULE PROCEDURE VisualizeMesh
END INTERFACE

#if !FV_ENABLED
INTERFACE WriteGlobalNodeIDsToVTK_array
  MODULE PROCEDURE WriteGlobalNodeIDsToVTK_array
END INTERFACE
#endif

PUBLIC:: BuildVisuCoords
PUBLIC:: BuildSurfVisuCoords
PUBLIC:: VisualizeMesh
#if !FV_ENABLED
PUBLIC:: WriteGlobalNodeIDsToVTK_array
#endif

CONTAINS

!=================================================================================================================================
!> Converts the coordinates of the mesh to the visu-mesh.
!=================================================================================================================================
SUBROUTINE BuildVisuCoords()
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeVISUFVEqui,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP,nGlobalElems,NodeCoords,NGeo
USE MOD_Visu_Vars          ,ONLY: CoordsVisu_DG
USE MOD_Visu_Vars          ,ONLY: NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu,nElems_DG,mapDGElemsToAllElems
USE MOD_Visu_Vars          ,ONLY: nElemsAvg2D_DG,Avg2D
USE MOD_Visu_Vars          ,ONLY: Elem_IJK_glob,mapElemIJToDGElemAvg2D
#if FV_ENABLED
USE MOD_Visu_Vars          ,ONLY: FVAmountAvg2D,mapElemIJToFVElemAvg2D,nElemsAvg2D_FV
USE MOD_Visu_Vars          ,ONLY: NVisu_FV,nElems_FV,mapFVElemsToAllElems,hasFV_Elems
USE MOD_Visu_Vars          ,ONLY: CoordsVisu_FV,changedMeshFile,changedFV_Elems,changedAvg2D
#endif
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY: offsetElemMPI
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem, iElem_DG
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)
#if FV_ENABLED
INTEGER            :: iElem_FV
REAL,ALLOCATABLE   :: Vdm_N_NVisu_FV(:,:)
#endif
INTEGER            :: iElemAvg,ii,jj,kk
REAL,ALLOCATABLE   :: Elem_xGP_glob(:,:,:,:,:)
#if USE_MPI
INTEGER            :: nDOFPerProc(0:nProcessors-1),offsetDOF(0:nProcessors-1)
INTEGER            :: iProc
#endif
!===================================================================================================================================

! Convert coordinates to visu grid
SWRITE (*,*) "[MESH] Convert coordinates to visu grid (DG)"

ALLOCATE(Vdm_N_NVisu(0:NVisu,0:NGeo))
CALL GetVandermonde(NGeo,NodeTypeVISU,NVisu   ,NodeTypeVisuPosti  ,Vdm_N_NVisu   ,modal=.FALSE.)
#if FV_ENABLED
ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu_FV,NodeTypeVISUFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
#endif

! convert coords of DG elements
IF (Avg2D) THEN
#if USE_MPI
  ! For parallel averaging, the root will gather the whole mesh and convert the first layer to the
  ! visu grid.
  ! For the gather operation, we need to know the number of DOFs per processor
  DO iProc = 0, nProcessors-1
    nDOFPerProc(iProc) = (offsetElemMPI(iProc+1) - offsetElemMPI(iProc)) * 3*(PP_N+1)**(PP_dim)
  END DO ! iProc = 1, nProcessors
  IF (MPIRoot) THEN
    ! On the root, we need the recieve array and the offset for each proc
    ALLOCATE(Elem_xGP_glob(1:3,0:PP_N,0:PP_N,0:PP_NZ,nGlobalElems))
    DO iProc = 0, nProcessors-1
      offsetDOF(iProc) = offsetElemMPI(iProc) * 3*(PP_N+1)**(PP_dim)
    END DO ! iProc = 1, nProcessors
  END IF
  CALL MPI_GATHERV(Elem_xGP,nDOFPerProc(myRank),MPI_DOUBLE_PRECISION,&
                   Elem_xGP_glob,nDOFPerProc,offsetDOF,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
#else
  ALLOCATE(Elem_xGP_glob(1:3,0:PP_N,0:PP_N,0:PP_NZ,nGlobalElems))
  Elem_xGP_glob = Elem_xGP
#endif /* USE_MPI */

  IF (MPIRoot) THEN
    SDEALLOCATE(CoordsVisu_DG)
    ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:0,nElemsAvg2D_DG))
#if FV_ENABLED
    SDEALLOCATE(CoordsVisu_FV)
    ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV))
#endif
    DO iElem = 1,nGlobalElems
      ii = Elem_IJK_glob(1,iElem)
      jj = Elem_IJK_glob(2,iElem)
      kk = Elem_IJK_glob(3,iElem)
      IF (kk.EQ.1) THEN
#if FV_ENABLED
        IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN ! DG
#endif
          iElemAvg = mapElemIJToDGElemAvg2D(ii,jj)
          CALL ChangeBasis2D(3,PP_N,NVisu,Vdm_N_NVisu,Elem_xGP_glob(:,:,:,0,iElem),CoordsVisu_DG(:,:,:,0,iElemAvg))
#if FV_ENABLED
        ELSE ! FV
          iElemAvg = mapElemIJToFVElemAvg2D(ii,jj)
          CALL ChangeBasis2D(3,PP_N,NVisu_FV,Vdm_N_NVisu_FV,Elem_xGP_glob(:,:,:,0,iElem),CoordsVisu_FV(:,:,:,0,iElemAvg))
        END IF
#endif
      END IF
    END DO
    DEALLOCATE(Elem_xGP_glob)
  ELSE ! MPIRoot
    ! All other procs will have 0 average elements
    SDEALLOCATE(CoordsVisu_DG)
    ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:0,nElemsAvg2D_DG))
#if FV_ENABLED
    SDEALLOCATE(CoordsVisu_FV)
    ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV))
#endif
  END IF ! MPIRoot
ELSE
  SDEALLOCATE(CoordsVisu_DG)
  ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems_DG))
  DO iElem_DG = 1,nElems_DG
    iElem = mapDGElemsToAllElems(iElem_DG)
    ! CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_N_NVisu,Elem_xGP(:,:,:,:,iElem),CoordsVisu_DG(:,:,:,:,iElem_DG))
    CALL ChangeBasisVolume(3,NGeo,NVisu,Vdm_N_NVisu,NodeCoords(:,:,:,:,iElem),CoordsVisu_DG(:,:,:,:,iElem_DG))
  END DO

#if FV_ENABLED
  IF (hasFV_Elems) THEN
    SWRITE (*,*) "[MESH] Convert coordinates to visu grid (FV)"
    ! only NVisu changed, but NVisu_FV is independent of NVisu
    IF ((.NOT.changedMeshFile).AND.(.NOT.changedFV_Elems).AND.(.NOT.changedAvg2D)) RETURN
    !ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:PP_N))
    !CALL GetVandermonde(PP_N,NodeType,NVisu_FV,NodeTypeVISUFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
    ! convert coords of FV elements
    SDEALLOCATE(CoordsVisu_FV)
    ALLOCATE(CoordsVisu_FV(3,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),nElems_FV))
    DO iElem_FV = 1,nElems_FV
      iElem = mapFVElemsToAllElems(iElem_FV)
      CALL ChangeBasisVolume(3,PP_N,NVisu_FV,Vdm_N_NVisu_FV,Elem_xGP(:,:,:,:,iElem),CoordsVisu_FV(:,:,:,:,iElem_FV))
    END DO
    SDEALLOCATE(Vdm_N_NVisu_FV)
  END IF
#endif

END IF
SDEALLOCATE(Vdm_N_NVisu)

END SUBROUTINE BuildVisuCoords


!=================================================================================================================================
!> Converts the coordinates of the surface mesh to the surface visu-mesh.
!=================================================================================================================================
SUBROUTINE BuildSurfVisuCoords()
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars          ,ONLY: CoordsSurfVisu_DG,nBCSidesVisu_DG,mapAllBCSidesToDGVisuBCSides
USE MOD_Visu_Vars          ,ONLY: NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu
#if FV_ENABLED
USE MOD_Visu_Vars          ,ONLY: CoordsSurfVisu_FV,nBCSidesVisu_FV,mapAllBCSidesToFVVisuBCSides
USE MOD_Visu_Vars          ,ONLY: NVisu_FV,hasFV_Elems
USE MOD_Visu_Vars          ,ONLY: changedMeshFile,changedFV_Elems,changedBCnames
#endif
USE MOD_Interpolation_Vars ,ONLY: NodeTypeVisu,NodeTypeVISUFVEqui,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
USE MOD_Mesh_Vars          ,ONLY: Face_xGP,nBCSides
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSide,iSideVisu
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)
#if FV_ENABLED
REAL,ALLOCATABLE   :: Vdm_N_NVisu_FV(:,:)
#endif
CHARACTER(LEN=255) :: NodeType_loc
INTEGER            :: Nloc
!===================================================================================================================================
Nloc = PP_N
NodeType_loc = NodeType

! Convert coordinates to visu grid
SWRITE (*,*) "[MESH] Convert coordinates to surface visu grid (DG)"
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:Nloc))
CALL GetVandermonde(Nloc,NodeType_loc,NVisu   ,NodeTypeVisuPosti  ,Vdm_N_NVisu   ,modal=.FALSE.)
! convert coords of DG elements
SDEALLOCATE(CoordsSurfVisu_DG)
ALLOCATE(CoordsSurfVisu_DG(3,0:NVisu,0:ZDIM(NVisu),0:0,nBCSidesVisu_DG))
DO iSide=1,nBCSides
  iSideVisu = mapAllBCSidesToDGVisuBCSides(iSide)
  IF (iSideVisu.GT.0)THEN
    CALL ChangeBasisSurf(3,Nloc,NVisu,   Vdm_N_NVisu, Face_xGP(:,:,:,0,iSide),CoordsSurfVisu_DG(:,:,:,0,iSideVisu))
  END IF
END DO
SDEALLOCATE(Vdm_N_NVisu)

#if FV_ENABLED
IF (hasFV_Elems) THEN
  SWRITE (*,*) "[MESH] Convert coordinates to surface visu grid (FV)"
  IF ((.NOT.changedMeshFile).AND.(.NOT.changedFV_Elems).AND.(.NOT.changedBCnames)) RETURN
  ALLOCATE(Vdm_N_NVisu_FV(0:NVisu_FV,0:Nloc))
  CALL GetVandermonde(Nloc,NodeType_loc,NVisu_FV,NodeTypeVISUFVEqui,Vdm_N_NVisu_FV,modal=.FALSE.)
  ! convert coords of FV elements
  SDEALLOCATE(CoordsSurfVisu_FV)
  ALLOCATE(CoordsSurfVisu_FV(3,0:NVisu_FV,0:ZDIM(NVisu_FV),0:0,nBCSidesVisu_FV))
  DO iSide=1,nBCSides
    iSideVisu = mapAllBCSidesToFVVisuBCSides(iSide)
    IF (iSideVisu.GT.0)THEN
      CALL ChangeBasisSurf(3,Nloc,NVisu_FV,Vdm_N_NVisu_FV,Face_xGP(:,:,:,0,iSide),CoordsSurfVisu_FV(:,:,:,0,iSideVisu))
    END IF
  END DO
  SDEALLOCATE(Vdm_N_NVisu_FV)
END IF
#endif

END SUBROUTINE BuildSurfVisuCoords


#if !FV_ENABLED
!===================================================================================================================================
!> Subroutine to write 2D or 3D coordinates to VTK format
!===================================================================================================================================
SUBROUTINE WriteGlobalNodeIDsToVTK_array(NVisu,nElems,globalnodeids_out,globalnodeids,dim,DGFV)
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,ALLOCATABLE,TARGET,INTENT(INOUT) :: globalnodeids(:)
TYPE (CARRAY), INTENT(INOUT)             :: globalnodeids_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (nElems.EQ.0) THEN
  globalnodeids_out%len = 0
  RETURN
END IF

! create global node ID
CALL BuildGlobalNodeIDs(NVisu=NVisu,nElems=nElems,globalnodeids=globalnodeids,dim=dim,DGFV=DGFV)

! set the sizes of the arrays
! globalnodeids_out%len = (2**dim)*((NVisu+DGFV)/(1+DGFV))**dim*nElems
globalnodeids_out%len = (nVisu+1)**dim*nElems

! assign data to the arrays (no copy!!!)
globalnodeids_out%data = C_LOC(globalnodeids(1))

END SUBROUTINE WriteGlobalNodeIDsToVTK_array


!=================================================================================================================================
!> Calculates the global node ID for VTK D3 filter
!=================================================================================================================================
SUBROUTINE BuildGlobalNodeIDs(NVisu,nElems,globalnodeids,dim,DGFV)
! MODULES
USE MOD_Globals
USE MOD_Visu_Vars    ,ONLY: CoordsVisu_DG
USE MOD_SHA256       ,ONLY: SHA256
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                       :: NVisu
INTEGER,INTENT(IN)                       :: nElems
INTEGER,ALLOCATABLE,TARGET,INTENT(INOUT) :: globalnodeids(:)  !< stores the unique global node ID
INTEGER,INTENT(IN)                       :: dim               !< 3 = 3d connectivity, 2 = 2d connectivity
INTEGER,INTENT(IN)                       :: DGFV              !< flag indicating DG = 0 or FV = 1 data
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,iElem
INTEGER           :: NodeID
INTEGER           :: NVisu_k,NVisu_j
INTEGER           :: fVTKCells!,nVTKCells
INTEGER           :: iProc
INTEGER           :: nUniqueNodeHashes
INTEGER           :: nUniqueNodeHashesProc(0:nProcessors-1),offsetUniqueHashesProc(0:nProcessors-1)
CHARACTER(LEN=48) :: NodePosChar48
CHARACTER(LEN=SHA256_LENGTH),ALLOCATABLE :: NodeHash(:,:,:,:),NodeHashGlob(:),unique(:),sorted(:)
CHARACTER(LEN=SHA256_LENGTH)             :: min_val,max_val

!=================================================================================================================================

SELECT CASE(dim)
  CASE(3)
    NVisu_k = NVisu
    NVisu_j = NVisu
  CASE(2)
    NVisu_k = 1
    NVisu_j = NVisu
  CASE(1)
    NVisu_k = 1
    NVisu_j = 1
  CASE DEFAULT
    ! Dummy variables to make GCC happy
    NVisu_k = -1
    NVisu_j = -1
    CALL Abort(__STAMP__, "Only 2D and 3D connectivity can be created. dim must be 2 or 3.")
END SELECT

fVTKCells  = ((NVisu+DGFV)/(1+DGFV))**dim
! nVTKCells  = fVTKCells*nElems
SDEALLOCATE(globalnodeids)

! Reuse the existing elem distribution
ALLOCATE(NodeHash(0:NVisu,0:NVisu_j,0:NVisu_k,nElems))

! Transfer node coordinates to continous data string and compute hash
DO iElem = 1,nElems
  DO k = 0,NVisu_k!,(DGFV+1)
    DO j = 0,NVisu_j!,(DGFV+1)
      DO i = 0,NVisu!,(DGFV+1)
        ! Just write the coordinates into the string without space
        WRITE(NodePosChar48,"(3E16.10)") CoordsVisu_DG(1:3,i,j,k,iElem)

        ! Hash the string (strictly we don't even need to hash it?)
        !> SHA256 is definitely overkill but Mikael Leetmaa provided a library
        !> Keep it allocated because we need it again when assigning unique global node IDs
        NodeHash(i,j,k,iElem) = SHA256(NodePosChar48)
      END DO
    END DO
  END DO
END DO

! Sort the NodeHash array
ALLOCATE(unique((nElems)*(2**dim)*fVTKCells))
min_val = MINVAL(NodeHash)
max_val = MAXVAL(NodeHash)

! Initialize
nUniqueNodeHashes         = 1
unique(nUniqueNodeHashes) = min_val

DO WHILE(LLT(min_val,max_val))
  nUniqueNodeHashes = nUniqueNodeHashes+1
  min_val = minval(NodeHash, mask=LGT(NodeHash,min_val))
  unique(nUniqueNodeHashes) = min_val
END DO

! Write the sorted values back
ALLOCATE(sorted(nUniqueNodeHashes))
sorted = unique(1:nUniqueNodeHashes)
DEALLOCATE(unique)

#if USE_MPI
! Communicate the size of the unique nodes to the root
CALL MPI_GATHER(nUniqueNodeHashes,1,MPI_INTEGER,nUniqueNodeHashesProc,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
! CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
! Allocate recv buffer on the root
IF (MPIRoot) THEN
  ALLOCATE(NodeHashGlob(SUM(nUniqueNodeHashesProc)))
  ! Calculate offsets for procs
  offsetUniqueHashesProc(0) = 0
  DO iProc = 1,nProcessors-1
    offsetUniqueHashesProc(iProc) = nUniqueNodeHashesProc(iProc-1) + offsetUniqueHashesProc(iProc-1)
  END DO
END IF

CALL MPI_GATHERV(sorted                ,nUniqueNodeHashes    *SHA256_LENGTH                                     ,MPI_CHARACTER ,&
                 NodeHashGlob          ,nUniqueNodeHashesProc*SHA256_LENGTH,offsetUniqueHashesProc*SHA256_LENGTH,MPI_CHARACTER ,&
                 0,MPI_COMM_FLEXI,iError)

! Root sort the unique hashes from the nodes
IF (MPIRoot) THEN
  ALLOCATE(unique(SUM(nUniqueNodeHashesProc)))
  min_val = MINVAL(NodeHashGlob)
  max_val = MAXVAL(NodeHashGlob)

  nUniqueNodeHashes         = 1
  unique(nUniqueNodeHashes) = min_val
  DO WHILE(LLT(min_val,max_val))
    nUniqueNodeHashes = nUniqueNodeHashes+1
    min_val = minval(NodeHashGlob, mask=LGT(NodeHashGlob,min_val))
    unique(nUniqueNodeHashes) = min_val
  END DO
  DEALLOCATE(NodeHashGlob)
END IF

! Send the information back to the procs
SWRITE(UNIT_stdOut,'(A,I0)') ' Number of unique visu nodes: ',nUniqueNodeHashes
CALL MPI_BCAST(nUniqueNodeHashes,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
ALLOCATE(NodeHashGlob(nUniqueNodeHashes))

! Write the sorted values back
IF (MPIRoot) THEN
  NodeHashGlob = unique(1:nUniqueNodeHashes)
  DEALLOCATE(unique)
END IF

CALL MPI_BCAST(NodeHashGlob,nUniqueNodeHashes*SHA256_LENGTH,MPI_CHARACTER,0,MPI_COMM_FLEXI,iError)
#else
ALLOCATE(NodeHashGlob(nUniqueNodeHashes))
NodeHashGlob = sorted
#endif
DEALLOCATE(sorted)

ALLOCATE(GlobalNodeIDs((nVisu+1)**dim*nElems))
GlobalNodeIDs = -1
NodeID = 0

! Finally, the unique node ID is just the position in the sorted NodeHashGlob array
DO iElem = 1,nElems
  DO k = 0,NVisu_k!,(DGFV+1)
    DO j = 0,NVisu_j!,(DGFV+1)
      DO i = 0,NVisu!,(DGFV+1)
        NodeID = NodeID + 1
        GlobalNodeIDs(NodeID) = FINDLOC(NodeHashGlob,NodeHash(i,j,k,iElem),1)
      END DO
    END DO
  END DO
END DO

DEALLOCATE(NodeHash)
DEALLOCATE(NodeHashGlob)

END SUBROUTINE BuildGlobalNodeIDs
#endif


!=================================================================================================================================
!> Visualize mesh only
!> 1. read mesh
!> 2. BuildVisuCoords
!> 3. Convert scaled jacobian
!> 4. write mesh to VTK array
!> 5. set length of all other output arrays to zero
!=================================================================================================================================
SUBROUTINE VisualizeMesh(postifile,meshfile_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_ReadInTools   ,ONLY: prms,GETINT,GETSTR,CountOption
USE MOD_ReadInTools   ,ONLY: FinalizeParameters
USE MOD_StringTools   ,ONLY: STRICMP
#if USE_MPI
USE MOD_MPI           ,ONLY: FinalizeMPI
#endif
USE MOD_Interpolation ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh_Vars     ,ONLY: nElems,Ngeo,scaledJac
USE MOD_Mesh          ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_VTK           ,ONLY: WriteCoordsToVTK_array
USE MOD_HDF5_Input    ,ONLY: ReadAttribute,File_ID,OpenDataFile,CloseDataFile
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN):: postifile
CHARACTER(LEN=255),INTENT(IN):: meshfile_in
! LOCAL VARIABLES
INTEGER             :: iElem,nVarIni,iVar,jVar,iVarVisu,meshModeLoc
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
  NVisu = GETINT('NVisu','1') ! Degree of visualization basis
ELSE
  NVisu = 2*NGeo ! TODO: correct?
END IF
NVisu_FV = 1

! read mesh, depending if we should visualize the Jacobian or not different mesh modes are needed (calculate metrics or not)
nVarIni=CountOption("VarName")
meshModeLoc = 0
IF (nVarIni.GT.0) meshModeLoc=2
CALL InitInterpolation(Ngeo)
CALL InitMesh(meshMode=meshModeLoc, MeshFile_IN=meshfile_in)

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
  nVarDep = 2
  nVarAll = 2
  SDEALLOCATE(mapDepToCalc)
  SDEALLOCATE(mapAllVarsToVisuVars)
  SDEALLOCATE(mapAllVarsToSurfVisuVars)
  ALLOCATE(mapDepToCalc(nVarDep))
  mapDepToCalc(1) = 1
  mapDepToCalc(2) = 2
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
  UCalc_DG(:,:,:,:,1) = scaledJac
  DO iElem=1,nElems
    UCalc_DG(:,:,:,iElem,2) = MINVAL(UCalc_DG(:,:,:,iElem,1))
  END DO ! iElem

  CALL ConvertToVisu_DG()
ELSE
  nVarVisu = 0
END IF

CALL FinalizeInterpolation()
CALL FinalizeParameters()
END SUBROUTINE VisualizeMesh


#if GCC_VERSION < 100000
PURE FUNCTION FINDLOC(Array,Value,Dim)
!===================================================================================================================================
!> Implements a subset of the intrinsic FINDLOC function for Fortran < 2008
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
CHARACTER(LEN=SHA256_LENGTH),INTENT(IN) :: Array(:)
CHARACTER(LEN=SHA256_LENGTH),INTENT(IN) :: Value
INTEGER,INTENT(IN)                      :: Dim
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER                        :: FINDLOC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVar
!===================================================================================================================================
DO iVar = 1,SIZE(ARRAY,1)
  IF (Array(iVar).EQ.Value) THEN
    FINDLOC = iVar
    RETURN
  END IF
END DO

! Return error code -1 if the value was not found
FINDLOC = -1

END FUNCTION FINDLOC
#endif /*GCC_VERSION < 100000*/

END MODULE MOD_Posti_VisuMesh
