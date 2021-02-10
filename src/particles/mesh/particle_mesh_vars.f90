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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_Mesh_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

LOGICAL            :: ParticleMeshInitIsDone   = .FALSE.
INTEGER            :: NGeoOverride
REAL               :: meshScale
!-----------------------------------------------------------------------------------------------------------------------------------
! Mesh info
REAL               :: MeshVolume                            ! total Volume of mesh
REAL               :: LocalVolume                           ! volume of proc
INTEGER            :: nNonUniqueGlobalSides                 !> total nb. of non-unique sides of mesh (hexahedral: 6*nElems)
INTEGER            :: nNonUniqueGlobalNodes                 !> total nb. of non-unique nodes of mesh (hexahedral: 8**NGeo * nElems)
INTEGER            :: nNonUniqueGlobalTrees                 !> total nb. of trees
INTEGER            :: nUniqueMasterMortarSides              !> total nb. of master mortar sides in the mesh
INTEGER            :: nUniqueBCSides                        !> total nb. of BC sides in the mesh
INTEGER            :: nComputeNodeElems                     !> Number of elems on current compute-node
INTEGER            :: nComputeNodeSides                     !> Number of sides on current compute-node
INTEGER            :: nComputeNodeNodes                     !> Number of nodes on current compute-node
INTEGER            :: nComputeNodeTrees                     !> Number of trees on current compute-node
INTEGER            :: offsetComputeNodeElem                 !> elem offset of compute-node root
INTEGER            :: offsetComputeNodeSide                 !> side offset of compute-node root
INTEGER            :: offsetComputeNodeNode                 !> node offset of compute-node root
INTEGER            :: offsetComputeNodeTree                 !> tree offset of compute-node root

! ====================================================================
! MPI3 shared variables
REAL,ALLOCPOINT,DIMENSION(:,:)           :: ElemBaryNGeo       ! element local basis: origin
REAL,ALLOCPOINT,DIMENSION(:)             :: ElemRadiusNGeo     ! radius of element
REAL,ALLOCPOINT,DIMENSION(:)             :: ElemRadius2NGeo    ! radius of element + 2% tolerance
REAL,ALLOCPOINT,DIMENSION(:,:)           :: slenXiEtaZetaBasis ! inverse of length of basis vector
REAL,ALLOCPOINT,DIMENSION(:,:,:)         :: XiEtaZetaBasis     ! element local basis vector (linear elem)

! XCL_NGeo and dXCL_NGeo always exist for DG mesh
REAL,POINTER,DIMENSION(:,:,:,:,:)        :: XCL_NGeo_Shared
REAL,POINTER,DIMENSION(:,:,:,:,:)        :: Elem_xGP_Shared
REAL,POINTER,DIMENSION(:,:,:,:,:,:)      :: dXCL_NGeo_Shared   ! Jacobi matrix of the mapping P\in NGeo

! FIBGM
INTEGER,ALLOCPOINT,DIMENSION(:,:,:)      :: FIBGM_nTotalElems  !> FastInitBackgroundMesh of global domain
INTEGER,ALLOCPOINT,DIMENSION(:,:,:)      :: FIBGM_nElems       !> FastInitBackgroundMesh of compute node
INTEGER,ALLOCPOINT,DIMENSION(:,:,:)      :: FIBGM_offsetElem   !> element offsets in 1D FIBGM_Element_Shared array
INTEGER,ALLOCPOINT,DIMENSION(:)          :: FIBGM_Element      !> element offsets in 1D FIBGM_Element_Shared array

LOGICAL,ALLOCPOINT,DIMENSION(:)          :: ElemCurved         !> flag if an element is curved

INTEGER,ALLOCPOINT,DIMENSION(:)          :: ElemToBCSides(:,:) !> Mapping from elem to BC sides within halo eps
REAL,ALLOCPOINT,DIMENSION(:,:)           :: SideBCMetrics(:,:) !> Metrics for BC sides, see piclas.h

REAL,ALLOCPOINT,DIMENSION(:,:,:,:)       :: ElemsJ             !> 1/DetJac for each Gauss Point
REAL,ALLOCPOINT,DIMENSION(:)             :: ElemEpsOneCell     !> tolerance for particle in inside ref element 1+epsinCell

! Boundary sides
INTEGER,ALLOCPOINT,DIMENSION(:)          :: BCSide2SideID      !> Mapping from compute-node BC side ID to global Side ID
INTEGER,ALLOCPOINT,DIMENSION(:)          :: SideID2BCSide      !> Inverse mapping
REAL,ALLOCPOINT,DIMENSION(:,:)           :: BCSideMetrics      !> Side origin and radius for each compute-node BC side

! FIBGM to proc mapping
INTEGER,ALLOCPOINT,DIMENSION(:,:,:,:)    :: FIBGMToProc
INTEGER,ALLOCPOINT,DIMENSION(:)          :: FIBGMProcs

! Shared arrays containing information for complete mesh
INTEGER,ALLOCPOINT,DIMENSION(:)          :: ElemToTree_Shared
INTEGER,ALLOCPOINT,DIMENSION(:,:)        :: ElemInfo_Shared
INTEGER,ALLOCPOINT,DIMENSION(:,:)        :: SideInfo_Shared
INTEGER,ALLOCPOINT,DIMENSION(:)          :: NodeInfo_Shared
REAL,ALLOCPOINT,DIMENSION(:,:)           :: NodeCoords_Shared
REAL,ALLOCPOINT,DIMENSION(:,:,:,:,:)     :: TreeCoords_Shared

! Shared arrays for halo debug information
INTEGER,ALLOCPOINT,DIMENSION(:,:)        :: ElemHaloInfo_Shared
INTEGER,ALLOCPOINT,DIMENSION(:)          :: ElemHaloInfo_Array

REAL,ALLOCPOINT    :: xiMinMax_Shared(:,:,:)

INTEGER,ALLOCPOINT :: ElemToBCSides_Shared(:,:)                !> Mapping from elem to BC sides within halo eps
REAL,ALLOCPOINT    :: SideBCMetrics_Shared(:,:)                !> Metrics for BC sides, see piclas.h
                                                               !> 1 - Global SideID
                                                               !> 2 - ElemID for BC side (non-unique)
                                                               !> 3 - Distance from BC side to element origin
                                                               !> 4 - Radius of BC Side
                                                               !> 5 - Origin of BC Side, x-coordinate
                                                               !> 6 - Origin of BC Side, y-coordinate
                                                               !> 7 - Origin of BC Side, z-coordinate

INTEGER,ALLOCPOINT :: ElemToBGM_Shared(:,:)                    !> BGM Bounding box around element (respective BGM indices) of compute node
INTEGER,ALLOCPOINT :: FIBGM_nTotalElems_Shared(:)              !> FastInitBackgroundMesh of global domain
INTEGER,ALLOCPOINT :: FIBGM_nElems_Shared(:)                   !> FastInitBackgroundMesh of compute node
INTEGER,ALLOCPOINT :: FIBGM_Element_Shared(:)                  !> FastInitBackgroundMesh of compute node
INTEGER,ALLOCPOINT :: FIBGM_offsetElem_Shared(:)

INTEGER,ALLOCPOINT :: FIBGMToProc_Shared(:,:,:,:)
INTEGER,ALLOCPOINT :: FIBGMProcs_Shared(:)

REAL,ALLOCPOINT    :: BoundsOfElem_Shared(:,:,:)               !> Cartesian bounding box around element

REAL,ALLOCPOINT    :: XCL_NGeo_Array(:)                        !> 1D array, pointer changes to proper array bounds
REAL,ALLOCPOINT    :: Elem_xGP_Array(:)                        !> 1D array, pointer changes to proper array bounds
REAL,ALLOCPOINT    :: dXCL_NGeo_Array(:)                       !> 1D array, pointer changes to proper array bounds
REAL,ALLOCPOINT    :: BezierControlPoints3D_Shared(:)          !> BezierControlPoints in 1D array. Pointer changes to proper array bounds
REAL,ALLOCPOINT    :: BezierControlPoints3DElevated_Shared(:)  !> BezierControlPoints in 1D array. Pointer changes to proper array bounds
REAL,ALLOCPOINT    :: ElemsJ_Shared(:)                         !> 1/DetJac for each Gauss Point. 1D array, pointer changes to proper array bounds
REAL,ALLOCPOINT    :: ElemEpsOneCell_Shared(:)                 !> tolerance for particle in inside ref element 1+epsinCell

REAL,ALLOCPOINT    :: ElemBaryNGeo_Shared(:,:)
REAL,ALLOCPOINT    :: ElemRadiusNGeo_Shared(:)
REAL,ALLOCPOINT    :: ElemRadius2NGeo_Shared(:)
REAL,ALLOCPOINT    :: XiEtaZetaBasis_Shared(:,:,:)
REAL,ALLOCPOINT    :: slenXiEtaZetaBasis_Shared(:,:)

LOGICAL,ALLOCPOINT :: ElemCurved_Shared(:)                     !> Flag if an element is curved
LOGICAL,ALLOCPOINT :: ConcaveElemSide_Shared(:,:)
INTEGER,ALLOCPOINT :: ElemNodeID_Shared(:,:)                   !> Contains the 8 corner nodes of an element, important for NGeo > 1
INTEGER,ALLOCPOINT :: ElemSideNodeID_Shared(:,:,:)             !> Contains the 4 corner nodes of the local sides in an element
REAL,ALLOCPOINT    :: ElemMidPoint_Shared(:,:)

REAL,ALLOCPOINT    :: SideSlabNormals_Shared(:,:,:)
REAL,ALLOCPOINT    :: SideSlabIntervals_Shared(:,:)
LOGICAL,ALLOCPOINT :: BoundingBoxIsEmpty_Shared(:)

INTEGER,ALLOCPOINT :: SideType_Shared(:)
REAL,ALLOCPOINT    :: SideDistance_Shared(:)
REAL,ALLOCPOINT    :: SideNormVec_Shared(:,:)

REAL,ALLOCPOINT    :: BaseVectors0_Shared(:,:)
REAL,ALLOCPOINT    :: BaseVectors1_Shared(:,:)
REAL,ALLOCPOINT    :: BaseVectors2_Shared(:,:)
REAL,ALLOCPOINT    :: BaseVectors3_Shared(:,:)
!REAL,ALLOCPOINT    :: BaseVectorsScale_Shared(:)

! Boundary sides
INTEGER,ALLOCPOINT :: BCSide2SideID_Shared(:)
INTEGER,ALLOCPOINT :: SideID2BCSide_Shared(:)
REAL,ALLOCPOINT    :: BCSideMetrics_Shared(:,:)

! Shared arrays containing information for mesh on compute node
REAL,ALLOCPOINT    :: ElemVolume_Shared(:)
!REAL,ALLOCPOINT    :: ElemCharLength_Shared(:)
!REAL,ALLOCPOINT    :: ElemCharLengthX_Shared(:)
!REAL,ALLOCPOINT    :: ElemCharLengthY_Shared(:)
!REAL,ALLOCPOINT    :: ElemCharLengthZ_Shared(:)

#if USE_MPI
! integers to hold shared memory windows
INTEGER           :: ElemToTree_Shared_Win
INTEGER           :: ElemInfo_Shared_Win
INTEGER           :: SideInfo_Shared_Win
INTEGER           :: NodeInfo_Shared_Win
INTEGER           :: NodeCoords_Shared_Win
INTEGER           :: TreeCoords_Shared_Win

LOGICAL           :: CalcHaloInfo                          !> Output halo element information to ElemData
INTEGER           :: ElemHaloInfo_Shared_Win

INTEGER           :: xiMinMax_Shared_Win

INTEGER           :: ElemToBCSides_Shared_Win
INTEGER           :: SideBCMetrics_Shared_Win

INTEGER           :: ElemToBGM_Shared_Win
INTEGER           :: FIBGM_nTotalElems_Shared_Win
INTEGER           :: FIBGM_nElems_Shared_Win
INTEGER           :: FIBGM_Element_Shared_Win
INTEGER           :: FIBGM_offsetElem_Shared_Win

INTEGER           :: FIBGMToProc_Shared_Win
INTEGER           :: FIBGMProcs_Shared_Win

INTEGER           :: BoundsOfElem_Shared_Win

INTEGER           :: XCL_NGeo_Shared_Win
INTEGER           :: Elem_xGP_Shared_Win
INTEGER           :: dXCL_NGeo_Shared_Win
INTEGER           :: BezierControlPoints3D_Shared_Win
INTEGER           :: BezierControlPoints3DElevated_Shared_Win
INTEGER           :: ElemsJ_Shared_Win
INTEGER           :: ElemEpsOneCell_Shared_Win

INTEGER           :: ElemBaryNGeo_Shared_Win
INTEGER           :: ElemRadiusNGeo_Shared_Win
INTEGER           :: ElemRadius2NGeo_Shared_Win
INTEGER           :: XiEtaZetaBasis_Shared_Win
INTEGER           :: slenXiEtaZetaBasis_Shared_Win

INTEGER           :: ElemCurved_Shared_Win
INTEGER           :: ConcaveElemSide_Shared_Win
INTEGER           :: ElemNodeID_Shared_Win
INTEGER           :: ElemSideNodeID_Shared_Win
INTEGER           :: ElemMidPoint_Shared_Win

INTEGER           :: SideSlabNormals_Shared_Win
INTEGER           :: SideSlabIntervals_Shared_Win
INTEGER           :: BoundingBoxIsEmpty_Shared_Win

INTEGER           :: SideType_Shared_Win
INTEGER           :: SideDistance_Shared_Win
INTEGER           :: SideNormVec_Shared_Win

INTEGER           :: BaseVectors0_Shared_Win
INTEGER           :: BaseVectors1_Shared_Win
INTEGER           :: BaseVectors2_Shared_Win
INTEGER           :: BaseVectors3_Shared_Win
!INTEGER           :: BaseVectorsScale_Shared_Win

! Boundary sides
INTEGER           :: BCSide2SideID_Shared_Win
INTEGER           :: SideID2BCSide_Shared_Win
INTEGER           :: BCSideMetrics_Shared_Win

! Shared arrays containing information for mesh on compute node
INTEGER           :: ElemVolume_Shared_Win
INTEGER           :: ElemMPVolumePortion_Shared_Win
INTEGER           :: ElemCharLength_Shared_Win
INTEGER           :: ElemCharLengthX_Shared_Win
INTEGER           :: ElemCharLengthY_Shared_Win
INTEGER           :: ElemCharLengthZ_Shared_Win
#endif

! ElemID for WriteHaloInfo
INTEGER,ALLOCATABLE                      :: ElemHaloID(:)

! ====================================================================
LOGICAL,ALLOCATABLE :: ElemHasAuxBCs(:,:)

INTEGER                                  :: RefMappingGuess                   ! select guess for mapping into reference element
                                                                              ! 1 - Linear, cubical element
                                                                              ! 2 - closest Gauss-Point
                                                                              ! 3 - closest XCL-point
                                                                              ! 4 - trivial guess - element origin
REAL                                     :: RefMappingEps                     ! tolerance for Netwton to get xi from X
REAL                                     :: epsInCell                         ! tolerance for eps for particle
                                                                              ! inside of ref element
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tFastInitBGM
  INTEGER                                :: nElem                             ! Number of elements in background mesh cell
  INTEGER, ALLOCATABLE                   :: Element(:)                        ! List of elements/physical cells in BGM cell
END TYPE

TYPE tGeometry
  REAL                                   :: CNxmin                            ! minimum x coord of all compute-node nodes
  REAL                                   :: CNxmax                            ! minimum y coord of all compute-node nodes
  REAL                                   :: CNymin                            ! minimum z coord of all compute-node nodes
  REAL                                   :: CNymax                            ! max x coord of all compute-node nodes
  REAL                                   :: CNzmin                            ! max y coord of all compute-node nodes
  REAL                                   :: CNzmax                            ! max z coord of all compute-node nodes
  REAL                                   :: xminglob                          ! global minimum x coord of all nodes
  REAL                                   :: yminglob                          ! global minimum y coord of all nodes
  REAL                                   :: zminglob                          ! global minimum z coord of all nodes
  REAL                                   :: xmaxglob                          ! global max x coord of all nodes
  REAL                                   :: ymaxglob                          ! global max y coord of all nodes
  REAL                                   :: zmaxglob                          ! global max z coord of all nodes
  REAL                                   :: xmin                              ! minimum x coord of all nodes
  REAL                                   :: xmax                              ! maximum x coord of all nodes
  REAL                                   :: ymin                              ! minimum y coord of all nodes
  REAL                                   :: ymax                              ! maximum y coord of all nodes
  REAL                                   :: zmin                              ! minimum z coord of all nodes
  REAL                                   :: zmax                              ! maximum z coord of all nodes
  ! periodic
  INTEGER                                :: nPeriodicVectors                  ! Number of periodic Vectors
  REAL, ALLOCATABLE                      :: PeriodicVectors(:,:)              ! PeriodicVectors(1:3,1:nPeriodicVectors), 1:3=x,y,z
  INTEGER,ALLOCATABLE                    :: DirPeriodicVectors(:)             ! direction of periodic vectors
  LOGICAL                                :: directions(3)                     ! flag for direction
  ! required for cartesian BGM for desposition
  INTEGER, ALLOCATABLE                   :: PeriodicBGMVectors(:,:)           ! = periodic vectors in backgroundmesh coords
  ! FIBGM
  REAL                                   :: FIBGMdeltas(3)                    ! size of background mesh cell for particle init
  REAL                                   :: FactorFIBGM(3)                    ! scaling factor for FIBGM

  ! caution, possible pointer
  TYPE (tFastInitBGM),ALLOCATABLE        :: FIBGM(:,:,:)  !        =>NULL()   ! FastInitBackgroundMesh
  INTEGER                                :: FIBGMimin                         ! smallest index of FastInitBGM (x)
  INTEGER                                :: FIBGMimax                         ! biggest index of FastInitBGM (x)
  INTEGER                                :: FIBGMjmin                         ! smallest index of FastInitBGM (y)
  INTEGER                                :: FIBGMjmax                         ! biggest index of FastInitBGM (y)
  INTEGER                                :: FIBGMkmin                         ! smallest index of FastInitBGM (z)
  INTEGER                                :: FIBGMkmax                         ! biggest index of FastInitBGM (z)
  INTEGER                                :: FIBGMiminglob
  INTEGER                                :: FIBGMimaxglob
  INTEGER                                :: FIBGMjminglob
  INTEGER                                :: FIBGMjmaxglob
  INTEGER                                :: FIBGMkminglob
  INTEGER                                :: FIBGMkmaxglob

  REAL, ALLOCATABLE                      :: CharLength(:)                     ! Characteristic length for each cell: L=V^(1/3)
  INTEGER, ALLOCATABLE                   :: ElemToRegion(:)                   ! ElemToRegion(1:nElems)
END TYPE

TYPE (tGeometry)                         :: GEO

INTEGER                                  :: WeirdElems                        ! Number of Weird Elements (=Elements which are folded
                                                                              ! into themselves)

INTEGER                                  :: NbrOfRegions                      ! Nbr of regions to be mapped to Elems
REAL, ALLOCATABLE                        :: RegionBounds(:,:)                 ! RegionBounds ((xmin,xmax,ymin,...)|1:NbrOfRegions)

!===================================================================================================================================
INTEGER                                  :: NGeoElevated                      !< polynomial degree of elevated geometric transformation
!-----------------------------------------------------------------------------------------------------------------------------------
! Interpolation - for Newton localisation of particles in curved Elements
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                         :: XiCL_NGeo(:)
REAL,ALLOCATABLE                         :: XiCL_NGeo1(:)
REAL,ALLOCATABLE,TARGET                  :: XCL_NGeo(:,:,:,:,:)               !> mapping X(xi) P\in Ngeo
REAL,ALLOCATABLE,TARGET                  :: dXCL_NGeo(:,:,:,:,:,:)            !> Jacobi matrix of the mapping P\in NGeo
REAL,ALLOCATABLE                         :: wBaryCL_NGeo(:)
REAL,ALLOCATABLE                         :: wBaryCL_NGeo1(:)
REAL,ALLOCATABLE                         :: DCL_NGeo(:,:)
REAL,ALLOCATABLE                         :: DCL_N(:,:)
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                         :: Xi_NGeo(:)                        !> 1D equidistant point positions for curved elements (during readin)
REAL                                     :: DeltaXi_NGeo
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                         :: Vdm_CLNGeo_CLN(:,:)
REAL,ALLOCATABLE                         :: Vdm_CLNGeo_GaussN(:,:)
REAL,ALLOCATABLE                         :: Vdm_CLN_GaussN(:,:)
REAL,ALLOCATABLE                         :: Vdm_CLNGeo1_CLNGeo(:,:)
REAL,ALLOCATABLE                         :: Vdm_NGeo_CLNGeo(:,:)

END MODULE MOD_Particle_Mesh_Vars
