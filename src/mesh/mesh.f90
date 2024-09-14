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

!==================================================================================================================================
!> Contains control routines to read high-order meshes, provide mesh data to the solver, build the metrics, partition the domain.
!==================================================================================================================================
MODULE MOD_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::FinalizeMesh
!==================================================================================================================================

PUBLIC::DefineParametersMesh
CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateStringOption(  'MeshFile',            "(relative) path to meshfile (mandatory).")
CALL prms%CreateLogicalOption( 'useCurveds',          "Controls usage of high-order information in mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption( 'interpolateFromTree', "For non-conforming meshes, built by refinement from a tree structure, "//&
                                                      "the metrics can be built from the tree geometry if it is contained "//&
                                                      "in the mesh. Can improve free-stream preservation.",&
                                                      '.TRUE.')
CALL prms%CreateRealOption(    'meshScale',           "Scale the mesh by this factor (shrink/enlarge).",&
                                                      '1.0')
CALL prms%CreateLogicalOption( 'meshdeform',          "Apply simple sine-shaped deformation on cartesion mesh (for testing).",&
                                                      '.FALSE.')
CALL prms%CreateLogicalOption( 'meshCheckRef',        "Flag if the mesh Jacobians should be checked in the reference system in "//&
                                                      "addition to the computational system.",'.TRUE.')
CALL prms%CreateLogicalOption( 'meshCheckConnectivity',"Flag if the mesh connectivity should be checked by comparing the face "//&
                                                      "Gauss points.",'.FALSE.')
#if (PP_dim == 3)
CALL prms%CreateLogicalOption( 'crossProductMetrics', "Compute mesh metrics using cross product form. Caution: in this case "//&
                                                      "free-stream preservation is only guaranteed for N=3*NGeo.",&
                                                      '.FALSE.')
#endif
CALL prms%CreateIntOption(     'debugmesh',           "Output file with visualization and debug information for the mesh. "//&
                                                      "0: no visualization, 3: Paraview binary",'0')
CALL prms%CreateStringOption(  'BoundaryName',        "Names of boundary conditions to be set (must be present in the mesh!)."//&
                                                      "For each BoundaryName a BoundaryType needs to be specified.",&
                                                      multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType',        "Type of boundary conditions to be set. Format: (BC_TYPE,BC_STATE)",&
                                                      multiple=.TRUE.)
CALL prms%CreateLogicalOption( 'writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
CALL prms%CreateIntOption(     'NGeoOverride',        "Override switch for NGeo. Interpolate mesh to different NGeo." //&
                                                      "<1: off, >0: Interpolate",'-1')
END SUBROUTINE DefineParametersMesh


!==================================================================================================================================
!> Routine controlling the initialization of the mesh.
!> - parameter and mesh reading
!> - domain partitioning
!> - allocation of mesh arrays
!> - build mesh mappings to handle volume/surface operations
!> - compute the mesh metrics
!> - provide mesh metrics for overintegration
!==================================================================================================================================
SUBROUTINE InitMesh(meshMode,MeshFile_IN,UseCurveds_IN)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_HDF5_Input
USE MOD_ChangeBasisByDim   ,ONLY:ChangeBasisVolume
USE MOD_Interpolation      ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone,NodeType,NodeTypeVISU
USE MOD_IO_HDF5,            ONLY:AddToElemData,ElementOut
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR,GETREAL,GETINT
USE MOD_Metrics,            ONLY:BuildCoords,CalcMetrics,CalcSurfMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Mappings,           ONLY:buildMappings
#if USE_MPI
USE MOD_Prepare_Mesh,       ONLY:exchangeFlip
#endif /*USE_MPI*/
#if FV_ENABLED
USE MOD_FV_Metrics,         ONLY:InitFV_Metrics,FV_CalcMetrics
#endif /*FV_ENABLED*/
#if (PP_dim == 2)
USE MOD_2D
#endif /*(PP_dim == 2)*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Metrics,ONLY:MoveCoords,MoveMetrics
USE MOD_LoadBalance_Vars   ,ONLY:PerformLoadBalance,UseH5IOLoadBalance
! USE MOD_Interpolation_Vars, ONLY:xGP
#endif /*USE_LOADBALANCE*/
#if USE_PARTICLES
USE MOD_Particle_Mesh      ,ONLY:InitParticleMeshBasis
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode !< 0: only read and build Elem_xGP,
                               !< 1: as 0 + build connectivity, 2: as 1 + calc metrics
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: MeshFile_IN !< file name of mesh to be read
LOGICAL,INTENT(IN),OPTIONAL            :: UseCurveds_IN
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: x(3)
REAL,POINTER      :: coords(:,:,:,:,:),coordsTmp(:,:,:,:,:),Vdm_EQNGeo_EQNGeoOverride(:,:)
INTEGER           :: iElem,i,j,k,nElemsLoc
LOGICAL           :: validMesh
INTEGER           :: firstMasterSide     ! lower side ID of array U_master/gradUx_master...
INTEGER           :: lastMasterSide      ! upper side ID of array U_master/gradUx_master...
INTEGER           :: firstSlaveSide      ! lower side ID of array U_slave/gradUx_slave...
INTEGER           :: lastSlaveSide       ! upper side ID of array U_slave/gradUx_slave...
! Timer
REAL              :: StartT,EndT,WallTime
!==================================================================================================================================
IF (.NOT.InterpolationInitIsDone .OR. MeshInitIsDone) &
  CALL CollectiveStop(__STAMP__,'InitMesh not ready to be called or already called.')

LBWRITE(UNIT_stdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A,I1,A)') ' INIT MESH IN MODE ',meshMode,'...'
GETTIME(StartT)

! prepare pointer structure (get nElems, etc.)
IF (PRESENT(MeshFile_IN)) THEN
  MeshFile = MeshFile_IN
ELSE
  MeshFile = GETSTR('MeshFile')
END IF

#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  validMesh = ISVALIDMESHFILE(MeshFile)
  IF(.NOT.validMesh) CALL CollectiveStop(__STAMP__,'ERROR - Mesh file not a valid HDF5 mesh.')
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

IF (PRESENT(UseCurveds_IN)) THEN
  useCurveds = UseCurveds_IN
ELSE
  useCurveds = GETLOGICAL('useCurveds')
END IF

#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo)
  CALL CloseDataFile()
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

IF(useCurveds.AND.(PP_N.LT.NGeo))THEN
#if USE_LOADBALANCE
  IF (.NOT.PerformLoadBalance) &
#endif /*USE_LOADBALANCE*/
  CALL PrintWarning('N<NGeo, for curved hexa normals are only approximated, can cause problems on periodic boundaries! Set N>=NGeo')
ENDIF

#if USE_PARTICLES
NGeoOverride = GETINT ('NGeoOverride','-1' )
IF (.NOT.UseCurveds) THEN
  LBWRITE(UNIT_stdOut,'(A)') ' NGeoOverride set but UseCurveds = F, setting NGeoOverride = -1'
  NGeoOverride = -1
END IF
meshScale    = GETREAL('meshScale'   ,'1.0')
#endif /*USE_PARTICLES*/
CALL ReadMesh(MeshFile) !set nElems

#if (PP_dim == 2)
! If this is a two dimensional calculation, all subsequent operations are performed on the reduced mesh.
LBWRITE(UNIT_stdOut,'(A)') " RUNNING A 2D SIMULATION! "
! The mesh coordinates read in by the readMesh routine are therefore reduced by one dimension.
CALL to2D_rank5((/1,0,0,0,1/),(/3,NGeo,NGeo,NGeo,nElems/),4,NodeCoords)
NodeCoords(3,:,:,:,:) = 0.
#endif

! if trees are available: compute metrics on tree level and interpolate to elements
interpolateFromTree=.FALSE.
IF(isMortarMesh) interpolateFromTree=GETLOGICAL('interpolateFromTree')
IF(interpolateFromTree)THEN
#if (PP_dim == 2)
  CALL CollectiveStop(__STAMP__,&
      "interpolateFromTree not supported in 2D.")
#endif
  coords    => TreeCoords
  NGeo      =  NGeoTree
  nElemsLoc =  nTrees
ELSE
  coords    => NodeCoords
  nElemsLoc =  nElems
ENDIF

#if !USE_PARTICLES
NGeoOverride=GETINT('NGeoOverride')
#endif /*!USE_PARTICLES*/
IF(NGeoOverride.GT.0)THEN
  ALLOCATE(CoordsTmp(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
  ALLOCATE(Vdm_EQNGeo_EQNGeoOverride(0:NGeoOverride,0:NGeo))
  CALL GetVandermonde(Ngeo, NodeTypeVISU, NgeoOverride, NodeTypeVISU, Vdm_EQNgeo_EQNgeoOverride, modal=.FALSE.)
  DO iElem=1,nElemsLoc
    CALL ChangeBasisVolume(3,Ngeo,NgeoOverride,Vdm_EQNGeo_EQNGeoOverride,coords(:,:,:,:,iElem),coordsTmp(:,:,:,:,iElem))
  END DO
  ! cleanup
  IF(interpolateFromTree)THEN
    DEALLOCATE(TreeCoords)
    ALLOCATE(TreeCoords(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
    coords=>TreeCoords
  ELSE
    DEALLOCATE(NodeCoords)
    ALLOCATE(NodeCoords(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
    coords=>NodeCoords
  END IF
  Coords = CoordsTmp
  NGeo = NGeoOverride
  DEALLOCATE(CoordsTmp, Vdm_EQNGeo_EQNGeoOverride)
END IF

#if USE_PARTICLES
IF(NGeoOverride.GT.0)THEN
#endif /*!USE_PARTICLES*/
LBWRITE(UNIT_stdOut,'(A3,A30,A3,I0)')' | ','Ngeo',' | ', Ngeo
#if USE_PARTICLES
End IF
#endif /*!USE_PARTICLES*/

! scale and deform mesh if desired (warning: no mesh output!)
#if !USE_PARTICLES
meshScale=GETREAL('meshScale')
#endif /*!USE_PARTICLES*/
IF(ABS(meshScale-1.).GT.1e-14)&
  Coords = Coords*meshScale

IF(GETLOGICAL('meshdeform'))THEN
#if USE_PARTICLES
  CALL COLLECTIVESTOP(__STAMP__,'Mesh deformation currently not supported in conjunction with particles!')
#endif
  DO iElem=1,nElemsLoc
    DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
      x(:)=Coords(:,i,j,k,iElem)
#if PP_dim==3
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))*SIN(PP_Pi*x(3))
#else
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))
#endif
    END DO; END DO; END DO;
  END DO
END IF

#if USE_LOADBALANCE
IF (PerformLoadBalance .AND. .NOT.UseH5IOLoadBalance) THEN
  CALL MoveCoords()
ELSE
! Build the coordinates of the solution gauss points in the volume
#endif /*USE_LOADBALANCE*/
  SDEALLOCATE(Elem_xGP)
  ALLOCATE(Elem_xGP      (3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  IF(interpolateFromTree)THEN
    CALL BuildCoords(NodeCoords,NodeType,PP_N,Elem_xGP,TreeCoords)
  ELSE
    CALL BuildCoords(NodeCoords,NodeType,PP_N,Elem_xGP)
  END IF
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

! Return if no connectivity and metrics are required (e.g. for visualization mode)
IF (meshMode.GT.0) THEN

  LBWRITE(UNIT_stdOut,'(A)') " NOW CALLING setLocalSideIDs..."
  CALL setLocalSideIDs()

#if USE_MPI
  ! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
  LBWRITE(UNIT_stdOut,'(A)') " NOW CALLING exchangeFlip..."
  CALL exchangeFlip()
#endif

  !RANGES
  !-----------------|-----------------|-------------------|
  !    U_master     | U_slave         |    FLUX           |
  !-----------------|-----------------|-------------------|
  !  BCsides        |                 |    BCSides        |
  !  InnerMortars   |                 |    InnerMortars   |
  !  InnerSides     | InnerSides      |    InnerSides     |
  !  MPI_MINE sides | MPI_MINE sides  |    MPI_MINE sides |
  !                 | MPI_YOUR sides  |    MPI_YOUR sides |
  !  MPIMortars     |                 |    MPIMortars     |
  !-----------------|-----------------|-------------------|

  firstBCSide          = 1
  firstMortarInnerSide = firstBCSide         +nBCSides
  firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
  firstMPISide_MINE    = firstInnerSide      +nInnerSides
  firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
  firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR

  lastBCSide           = firstMortarInnerSide-1
  lastMortarInnerSide  = firstInnerSide    -1
  lastInnerSide        = firstMPISide_MINE -1
  lastMPISide_MINE     = firstMPISide_YOUR -1
  lastMPISide_YOUR     = firstMortarMPISide-1
  lastMortarMPISide    = nSides

  firstMasterSide = 1
  lastMasterSide  = nSides
  firstSlaveSide  = firstInnerSide
  lastSlaveSide   = lastMPISide_YOUR
  nSidesMaster    = lastMasterSide-firstMasterSide+1
  nSidesSlave     = lastSlaveSide -firstSlaveSide+1

  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8)')   'first/lastMasterSide     ', firstMasterSide,lastMasterSide
  LOGWRITE(*,'(A25,I8)')   'first/lastSlaveSide      ', firstSlaveSide, lastSlaveSide
  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8,I8)')'first/lastBCSide         ', firstBCSide         ,lastBCSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMortarInnerSide', firstMortarInnerSide,lastMortarInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastInnerSide      ', firstInnerSide      ,lastInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_MINE   ', firstMPISide_MINE   ,lastMPISide_MINE
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_YOUR   ', firstMPISide_YOUR   ,lastMPISide_YOUR
  LOGWRITE(*,'(A30,I8,I8)')'first/lastMortarMPISide  ', firstMortarMPISide  ,lastMortarMPISide
  LOGWRITE(*,*)'-------------------------------------------------------'

  ! fill ElemToSide, SideToElem,BC
  ALLOCATE(ElemToSide(2,6,nElems))
  ALLOCATE(SideToElem(5,nSides))
  ALLOCATE(BC(1:nBCSides))
  ALLOCATE(AnalyzeSide(1:nSides))
  ! ALLOCATE(SideToGlobalSide(1:nSides))
  ElemToSide  = 0
  SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
  BC          = 0
  AnalyzeSide = 0

  !NOTE: nMortarSides=nMortarInnerSides+nMortarMPISides
  ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo
  ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides)) ! [1]: 1: Neighbour sides, 2: Flip, [2]: small sides
  MortarType=0
  MortarInfo=-1

  LBWRITE(UNIT_stdOut,'(A)') " NOW CALLING fillMeshInfo..."
  CALL fillMeshInfo()

#if (PP_dim ==2)
  ! In 2D, there is only one flip for the slave sides (1)
  SideToElem(S2E_FLIP,:)   = MIN(1,SideToElem(S2E_FLIP,:))
  ElemToSide(:,1,:)        = -999
  ElemToSide(:,6,:)        = -999
  ElemToSide(E2S_FLIP,:,:) = MIN(1,ElemToSide(E2S_FLIP,:,:))
  MortarInfo(MI_FLIP,:,:)  = MIN(1,MortarInfo(MI_FLIP,:,:))
#endif

  ! Build necessary mappings
  CALL buildMappings(PP_N,V2S=V2S,S2V=S2V,S2V2=S2V2,FS2M=FS2M,dim=PP_dim)
END IF
! deallocate pointers
LBWRITE(UNIT_stdOut,'(A)') " NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----

IF (meshMode.GT.1) THEN
#if USE_LOADBALANCE
  IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
    ! Shift metric arrays during load balance
    CALL MoveMetrics()
  ELSE
    ! deallocate existing arrays
    ! mesh basis
    ! SDEALLOCATE(Xi_NGeo)
    ! SDEALLOCATE(XiCL_NGeo)
    ! SDEALLOCATE(DCL_N)
    ! SDEALLOCATE(DCL_NGeo)
    ! SDEALLOCATE(Vdm_CLN_GaussN)
    ! SDEALLOCATE(Vdm_CLNGeo_GaussN)
    ! SDEALLOCATE(Vdm_CLNGeo_CLN)
    ! SDEALLOCATE(Vdm_NGeo_CLNGeo)
    ! SDEALLOCATE(wBaryCL_NGeo)

    ! mesh metrics
    SDEALLOCATE(      dXCL_N)
    SDEALLOCATE(      JaCL_N)
    SDEALLOCATE(Metrics_fTilde)
    SDEALLOCATE(Metrics_gTilde)
    SDEALLOCATE(Metrics_hTilde)
    SDEALLOCATE(sJ)
    SDEALLOCATE(scaledJac)
    SDEALLOCATE(DetJac_Ref)

    ! Vandermonde
    SDEALLOCATE(Vdm_CLN_N)
    SDEALLOCATE(XCL_N)

!     SDEALLOCATE(XCL_NGeo)
! #if USE_PARTICLES
!     SDEALLOCATE(dXCL_NGeo)
! #endif /*USE_PARTICLES*/

    ! Calculate metric arrays
#endif /*USE_LOADBALANCE*/
    ! volume data
    ALLOCATE(       XCL_N(  3,0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems))
    ALLOCATE(      dXCL_N(3,3,0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems))
    ALLOCATE(      JaCL_N(3,3,0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems))
    ALLOCATE(Metrics_fTilde(3,0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems,0:FV_SIZE))
    ALLOCATE(Metrics_gTilde(3,0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems,0:FV_SIZE))
    ALLOCATE(Metrics_hTilde(3,0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems,0:FV_SIZE))
    ALLOCATE(sJ            (  0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems,0:FV_SIZE))
    ALLOCATE(     scaledJac(  0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems))
    ALLOCATE(    DetJac_N(  1,0:PP_N,   0:PP_N,   0:PP_NZ        ,nElems))
    NGeoRef = 3*NGeo ! build jacobian at higher degree
    JaCL_N  = 0.
    ALLOCATE(    DetJac_Ref(1,0:NGeoRef,0:NGeoRef,0:ZDIM(NGeoRef),nElems))

    ! Vandermonde
    ALLOCATE(Vdm_CLN_N(       0:PP_N,0:PP_N))

#if USE_PARTICLES
    CALL InitParticleMeshBasis()
#endif

    ! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
    ! assign 1/detJ (sJ)
    ! assign normal and tangential vectors and surfElems on faces
#if (PP_dim == 3)
    ! compute metrics using cross product instead of curl form (warning: no free stream preservation!)
    crossProductMetrics = GETLOGICAL('crossProductMetrics')
#endif
    LBWRITE(UNIT_stdOut,'(A)') ' NOW CALLING calcMetrics...'
    ! CALL InitMeshBasis(NGeo,PP_N,xGP)

    CALL CalcMetrics()
#if USE_LOADBALANCE
  END IF
#endif /*USE_LOADBALANCE*/

  ! surface data
  ALLOCATE(      Face_xGP(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
  ALLOCATE(       NormVec(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
  ALLOCATE(      TangVec1(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
  ALLOCATE(      TangVec2(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
  ALLOCATE(      SurfElem(  0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
  ALLOCATE(     Ja_Face(3,3,0:PP_N,0:PP_NZ,          1:nSides)) ! temp
  ALLOCATE(    Ja_slave(3,3,0:PP_N,0:PP_NZ,          1:nSides)) ! temp
  Face_xGP = 0.
  NormVec  = 0.
  TangVec1 = 0.
  TangVec2 = 0.
  SurfElem = 0.

#if FV_ENABLED
  ! NOTE: initialize with 1 and not with 0
  ALLOCATE(sJ_master(1,0:PP_N,0:PP_NZ,1:nSides,0:FV_SIZE))
  ALLOCATE(sJ_slave( 1,0:PP_N,0:PP_NZ,1:nSides,0:FV_SIZE))
  sJ_master = 1.
  sJ_slave  = 1.
#endif

  ! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
  ! assign 1/detJ (sJ)
  ! assign normal and tangential vectors and surfElems on faces
  DO iElem = 1,nElems
    CALL CalcSurfMetrics(JaCL_N(:,:,:,:,:,iElem),iElem)
  END DO

#if FV_ENABLED
  CALL InitFV_Metrics()  ! Init FV metrics
  CALL FV_CalcMetrics()  ! FV metrics
#endif

! debugmesh: param specifies format to output, 0: no output, 1: tecplot ascii, 2: tecplot binary, 3: paraview binary
  CALL WriteDebugMesh(GETINT('debugmesh'))
END IF ! meshMode.GT.1

! IF (meshMode.GT.0) THEN
!   ALLOCATE(SideToGlobalSide(nSides))
!   DO iElem=1,nElems
! #if PP_dim == 3
!     DO LocSideID=1,6
! #else
!     DO LocSideID=2,5
! #endif
!       SideID = ElemToSide(E2S_SIDE_ID,LocSideID,iElem)
!       iSide  = ElemInfo(3,iElem+offsetElem) + LocSideID
!       SideToGlobalSide(SideID) = ABS(SideInfo(2,iSide))
!     END DO
!   END DO ! iElem
! END IF

! SDEALLOCATE(dXCL_N)
! SDEALLOCATE(Ja_Face)
! SDEALLOCATE(Ja_Slave)
SDEALLOCATE(TreeCoords)
SDEALLOCATE(xiMinMax)
SDEALLOCATE(ElemToTree)
IF ((.NOT.postiMode).AND.(ALLOCATED(scaledJac))) DEALLOCATE(scaledJac)

! Write debug information
CALL AddToElemData(ElementOut,'myRank',IntScalar=myRank)

! ALLOCATE(ElemGlobalID(1:nElems))
! DO iElem = 1,nElems
!   ElemGlobalID(iElem) = OffsetElem + iElem
! END DO ! iElem=1,nElems
! CALL AddToElemData(ElementOut,'ElemID',LongIntArray=ElemGlobalID)

MeshInitIsDone=.TRUE.

GETTIME(EndT)
WallTime = EndT-StartT
CALL DisplayMessageAndTime(WallTime,'INIT MESH DONE!',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitMesh

!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE FinalizeMesh()
! MODULES
USE MOD_Mesh_Vars
USE MOD_Mappings             ,ONLY: FinalizeMappings
#if FV_ENABLED
USE MOD_FV_Vars              ,ONLY: FV_Elems_master
USE MOD_FV_Metrics           ,ONLY: FinalizeFV_Metrics
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars     ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_PARTICLES
USE MOD_Particle_Mesh        ,ONLY: FinalizeParticleMeshBasis
#endif /*USE_PARTICLES*/
#if USE_LOADBALANCE || USE_PARTICLES
USE MOD_Mesh_Shared          ,ONLY: FinalizeMeshShared
#endif /*USE_LOADBALANCE || USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(ElemInfo)
! mapping from elems to sides and vice-versa
SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
SDEALLOCATE(AnalyzeSide)

! mortars
SDEALLOCATE(MortarType)
SDEALLOCATE(MortarInfo)

! allocated during ReadMesh
SDEALLOCATE(NodeCoords)

! Volume
! SDEALLOCATE(Elem_xGP)
! SDEALLOCATE(Metrics_fTilde)
! SDEALLOCATE(Metrics_gTilde)
! SDEALLOCATE(Metrics_hTilde)
! SDEALLOCATE(sJ)
SDEALLOCATE(scaledJac)
#if FV_ENABLED
SDEALLOCATE(sJ_master)
SDEALLOCATE(sJ_slave)
#endif

! surface
SDEALLOCATE(Face_xGP)
SDEALLOCATE(NormVec)
SDEALLOCATE(TangVec1)
SDEALLOCATE(TangVec2)
SDEALLOCATE(SurfElem)

SDEALLOCATE(Ja_Face)
SDEALLOCATE(Ja_slave)

! ijk sorted mesh
SDEALLOCATE(Elem_IJK)
SDEALLOCATE(ElemInfo)
SDEALLOCATE(SideInfo)
! SDEALLOCATE(SideToGlobalSide)

! mappings
CALL FinalizeMappings()

#if FV_ENABLED
SDEALLOCATE(FV_Elems_master) ! moved here from fv.f90
CALL FinalizeFV_Metrics()
#endif

#if USE_PARTICLES
CALL FinalizeParticleMeshBasis()
#endif /*USE_PARTICLES*/

MeshInitIsDone = .FALSE.

! Shared mesh readin happens during mesh readin, finalize with gathered routine here
#if USE_PARTICLES || USE_LOADBALANCE
CALL FinalizeMeshShared()
#endif /*USE_PARTICLES || USE_LOADBALANCE*/

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

! geometry information and VDMS
! SDEALLOCATE(DCL_N)
! SDEALLOCATE(DCL_NGeo)
! SDEALLOCATE(Vdm_CLN_GaussN)
! SDEALLOCATE(Vdm_CLNGeo_GaussN)
! SDEALLOCATE(Vdm_CLNGeo_CLN)
! SDEALLOCATE(Vdm_NGeo_CLNGeo)
! SDEALLOCATE(Xi_NGeo)
! SDEALLOCATE(XiCL_NGeo)
! SDEALLOCATE(wBaryCL_NGeo)
! particle input/output
SDEALLOCATE(Vdm_N_EQ)
SDEALLOCATE(Vdm_EQ_N)
! superb
SDEALLOCATE(Vdm_GL_N)
SDEALLOCATE(Vdm_N_GL)
! BCS
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
! elem-xgp and metrics
! SDEALLOCATE(XCL_NGeo)
! SDEALLOCATE(dXCL_NGeo)
SDEALLOCATE(Metrics_fTilde)
SDEALLOCATE(Metrics_gTilde)
SDEALLOCATE(Metrics_hTilde)
SDEALLOCATE(Elem_xGP)
SDEALLOCATE(sJ)
SDEALLOCATE(DetJac_N)
SDEALLOCATE(DetJac_Ref)
SDEALLOCATE(JaCL_N)
SDEALLOCATE(Vdm_CLN_N)
SDEALLOCATE(XCL_N)
! #if USE_LOADBALANCE
SDEALLOCATE(dXCL_N)
! SDEALLOCATE(Ja_Face)
! #endif /*USE_LOADBALANCE*/

END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
