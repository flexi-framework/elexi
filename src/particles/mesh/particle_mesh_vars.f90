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

LOGICAL             :: ParticleMeshInitIsDone
!-----------------------------------------------------------------------------------------------------------------------------------
! Mesh info
REAL                                     :: MeshVolume         ! total Volume of mesh
REAL                                     :: LocalVolume        ! volume of proc
! ====================================================================
! MPI3 shared variables
REAL,ALLOCPOINT,DIMENSION(:,:)           :: ElemBaryNGeo       ! element local basis: origin
REAL,ALLOCPOINT,DIMENSION(:)             :: ElemRadiusNGeo     ! radius of element
REAL,ALLOCPOINT,DIMENSION(:)             :: ElemRadius2NGeo    ! radius of element + 2% tolerance
REAL,ALLOCPOINT,DIMENSION(:,:)           :: slenXiEtaZetaBasis ! inverse of length of basis vector
REAL,ALLOCPOINT,DIMENSION(:,:,:,:,:)     :: XCL_NGeo_Shared
REAL,ALLOCPOINT,DIMENSION(:,:,:,:,:,:)   :: dXCL_NGeo_Shared   ! Jacobi matrix of the mapping P\in NGeo
REAL,ALLOCPOINT,DIMENSION(:,:,:)         :: XiEtaZetaBasis     ! element local basis vector (linear elem)

! FIBGM
INTEGER,ALLOCPOINT,DIMENSION(:,:,:)      :: FIBGM_nElems       !> FastInitBackgroundMesh of compute node
INTEGER,ALLOCPOINT,DIMENSION(:,:,:)      :: FIBGM_offsetElem   !> element offsets in 1D FIBGM_Element_Shared array
INTEGER,ALLOCPOINT,DIMENSION(:)          :: FIBGM_Element      !> element offsets in 1D FIBGM_Element_Shared array

LOGICAL,ALLOCPOINT,DIMENSION(:)          :: ElemCurved         ! flag if an element is curved

INTEGER,ALLOCPOINT,DIMENSION(:)          :: ElemToBCSides(:,:) ! Mapping from elem to BC sides within halo eps
REAL,ALLOCPOINT,DIMENSION(:,:)           :: SideBCMetrics(:,:) ! Metrics for BC sides, see piclas.h

REAL,ALLOCPOINT,DIMENSION(:,:,:,:)       :: ElemsJ             !< 1/DetJac for each Gauss Point
REAL,ALLOCPOINT,DIMENSION(:)             :: ElemEpsOneCell     ! tolerance for particle in inside ref element 1+epsinCell

! periodic case
INTEGER, ALLOCATABLE                     :: casematrix(:,:)                              ! matrix to compute periodic cases
INTEGER                                  :: NbrOfCases                                   ! Number of periodic cases
LOGICAL                                  :: PeriodicCheck                                ! Flag if periodic vectors are checked


! general: periodic sides have to be Cartesian
INTEGER,ALLOCATABLE :: SidePeriodicType(:)                                                ! 1:nTotalSides, periodic type of side
                                                                                          ! 0 - normal or BC side
                                                                                          ! >0 type of periodic displacement
REAL,ALLOCATABLE    :: SidePeriodicDisplacement(:,:)                                      ! displacement vector

! ====================================================================
INTEGER,ALLOCATABLE :: PartElemToSide(:,:,:)                                              ! extended list: 1:2,1:6,1:nTotalElems
                                                                                          ! ElemToSide: my geometry + halo
                                                                                          ! geometry + halo information
! -> this information is now in ElemInfo_Shared
! ====================================================================


! ====================================================================
INTEGER,ALLOCATABLE :: PartSideToElem(:,:)                                                ! extended list: 1:5,1:6,1:nTotalSides
                                                                                          ! SideToElem: my geometry + halo
                                                                                          ! geometry + halo information
! -> this information is now in SideInfo_Shared
! ====================================================================


! ====================================================================
INTEGER(KIND=8),ALLOCATABLE :: PartElemToElemGlob(:,:,:)                                  ! Mapping from ElemToElem
                                                                                          ! 1:4,1:6,1:nTotalElems
                                                                                          ! now in global-elem-ids !!!
INTEGER(KIND=4),ALLOCATABLE :: PartElemToElemAndSide(:,:,:)                               ! Mapping from ElemToElem
                                                                                          ! 1:8,1:6,1:nTotalElems
                                                                                          ! [1]1:4 - MortarNeighborElemID
                                                                                          ! [1]5:8 -       Neighbor locSideID
                                                                                          ! [2]1:6 - locSideID
                                                                                          ! [3]    - nTotalElems
                                                                                          ! now in global-elem-ids !!!
! -> this information is now deprecated
! ====================================================================


! ====================================================================
INTEGER             :: nPartSides                                                         ! nPartSides - nSides+nPartPeriodicSides
INTEGER             :: nTotalSides                                                        ! total nb. of sides (my+halo)
INTEGER             :: nPartPeriodicSides                                                 ! total nb. of sides (my+halo)
INTEGER             :: nTotalElems                                                        ! total nb. of elems (my+halo)
INTEGER             :: nTotalNodes                                                        ! total nb. of nodes (my+halo)
! -> this information is now shared (nComputeNodeSides)
! ====================================================================

INTEGER,ALLOCATABLE :: TracingBCInnerSides(:)                                             ! number of local element boundary faces
                                                                                          ! used for tracing (connected to element)
INTEGER,ALLOCATABLE :: TracingBCTotalSides(:)                                             ! total number of element boundary faces
                                                                                          ! used for tracing (loc faces + other
                                                                                          ! element faces that are possibly reached)
INTEGER,ALLOCATABLE :: ElemType(:)              !< Type of Element 1: only planar side, 2: one bilinear side 3. one curved side
LOGICAL,ALLOCATABLE :: ElemHasAuxBCs(:,:)
INTEGER             :: nTotalBCSides                                                      ! total number of BC sides (my+halo)
INTEGER             :: nTotalBCElems                                                      ! total number of bc elems (my+halo)
INTEGER,ALLOCATABLE :: PartBCSideList(:)                                                  ! mapping from SideID to BCSideID

INTEGER                                 :: RefMappingGuess                                ! select guess for mapping into reference
                                                                                          ! element
                                                                                          ! 1 - Linear, cubical element
                                                                                          ! 2 - closest Gauss-Point
                                                                                          ! 3 - closest XCL-point
                                                                                          ! 4 - trivial guess - element origin
REAL                                    :: RefMappingEps                                  ! tolerance for Netwton to get xi from X
REAL                                    :: epsInCell                                      ! tolerance for eps for particle
                                                                                          ! inside of ref element

!LOGICAL                                 :: DoRefMapping                  ! tracking by mapping particle into reference element
!-----------------------------------------------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tFastInitBGM
  INTEGER                                :: nElem                             ! Number of elements in background mesh cell
  INTEGER, ALLOCATABLE                   :: Element(:)                        ! List of elements/physical cells in BGM cell
#if USE_MPI
  INTEGER, ALLOCATABLE                   :: SharedProcs(:)                    ! first Entry: Number of Sharedprocs,
                                                                              ! following: SharedProcs
  !INTEGER                                :: nBCSides                          ! number BC sides in BGM cell
#endif
END TYPE

INTEGER                                  :: FIBGMCellPadding(1:3)
! -> this is now shared and correct (large enough halo region)
! ====================================================================

! ====================================================================
TYPE tNodeToElem
  INTEGER, ALLOCATABLE                   :: ElemID(:)
END TYPE
! -> this should be replaced with NodeInfo_Shared
! ====================================================================
TYPE tGeometry
  REAL                                   :: CNxmin                   ! minimum x coord of all compute-node nodes
  REAL                                   :: CNxmax                   ! minimum y coord of all compute-node nodes
  REAL                                   :: CNymin                   ! minimum z coord of all compute-node nodes
  REAL                                   :: CNymax                   ! max x coord of all compute-node nodes
  REAL                                   :: CNzmin                   ! max y coord of all compute-node nodes
  REAL                                   :: CNzmax                   ! max z coord of all compute-node nodes
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

  TYPE (tFastInitBGM),ALLOCATABLE        :: TFIBGM(:,:,:)  !       =>NULL()   ! FastInitBackgroundMesh
  INTEGER                                :: TFIBGMimin                        ! smallest index of FastInitBGM (x)
  INTEGER                                :: TFIBGMimax                        ! biggest index of FastInitBGM (x)
  INTEGER                                :: TFIBGMjmin                        ! smallest index of FastInitBGM (y)
  INTEGER                                :: TFIBGMjmax                        ! biggest index of FastInitBGM (y)
  INTEGER                                :: TFIBGMkmin                        ! smallest index of FastInitBGM (z)
  INTEGER                                :: TFIBGMkmax                        ! biggest index of FastInitBGM (z)

  INTEGER,ALLOCATABLE                    :: ElemToFIBGM(:,:)                  ! range of FIGMB cells per element
                                                                              ! 1:6,1:nTotalElems, xmin,max,yminmax,...
!  REAL, ALLOCATABLE                      :: Volume(:)                         ! Volume(nElems) for nearest_blurrycenter
  REAL, ALLOCATABLE                      :: CharLength(:)                     ! Characteristic length for each cell: L=V^(1/3)
!  REAL                                   :: MeshVolume                        ! Total Volume of mesh
!  REAL                                   :: LocalVolume                       ! Volume of proc
  INTEGER, ALLOCATABLE                   :: ElemToRegion(:)                   ! ElemToRegion(1:nElems)

  LOGICAL                                :: SelfPeriodic                      ! does process have periodic bounds with itself?
!  INTEGER, ALLOCATABLE                   :: ElemToNodeID(:,:)                 ! ElemToNodeID(1:nElemNodes,1:nElems)
!  INTEGER, ALLOCATABLE                   :: ElemSideNodeID(:,:,:)             ! ElemSideNodeID(1:nSideNodes,1:nLocSides,1:nElems)
                                                                              ! From element sides to node IDs
!  INTEGER, ALLOCATABLE                   :: ElemsOnNode(:)                    ! ElemSideNodeID(1:nSideNodes,1:nLocSides,1:nElems)
!  TYPE(tNodeToElem), ALLOCATABLE         :: NodeToElem(:)
!  INTEGER, ALLOCATABLE                   :: NumNeighborElems(:)
!  TYPE(tNodeToElem), ALLOCATABLE         :: ElemToNeighElems(:)
                                                                              ! From element sides to node IDs
  INTEGER, ALLOCATABLE                   :: PeriodicElemSide(:,:)             ! 0=not periodic side, others=PeriodicVectorsNum
!  LOGICAL, ALLOCATABLE                   :: ConcaveElemSide(:,:)              ! Whether LocalSide of Element is concave side
!  REAL, ALLOCATABLE                      :: NodeCoords(:,:)                   ! Node Coordinates (1:nDim,1:nNodes)
END TYPE

TYPE (tGeometry)                         :: GEO

INTEGER                                  :: WeirdElems                        ! Number of Weird Elements (=Elements which are folded
                                                                              ! into themselves)
LOGICAL                                  :: FindNeighbourElems=.FALSE.        ! Flag defining if mapping for neighbour elements
                                                                              ! is built via nodes

TYPE tBCElem
  INTEGER                                :: nInnerSides                       ! Number of BC-Sides of Element
  INTEGER                                :: lastSide                          ! total number of BC-Sides in eps-vicinity of element
  INTEGER, ALLOCATABLE                   :: BCSideID(:)                       ! List of elements in BGM cell
  REAL,ALLOCATABLE                       :: ElemToSideDistance(:)             ! stores the distance between each element and the
                                                                              ! sides associated with this element
END TYPE

TYPE (tBCElem),ALLOCATABLE               :: BCElem(:)

INTEGER                                  :: NbrOfRegions      ! Nbr of regions to be mapped to Elems
REAL, ALLOCATABLE                        :: RegionBounds(:,:) ! RegionBounds ((xmin,xmax,ymin,...)|1:NbrOfRegions)
LOGICAL,ALLOCATABLE                      :: isTracingTrouble(:)
REAL,ALLOCATABLE                         :: ElemTolerance(:)
INTEGER, ALLOCATABLE                     :: ElemToGlobalElemID(:)  ! mapping form local-elemid to global-id

!===================================================================================================================================
INTEGER          :: NGeoElevated                !< polynomial degree of elevated geometric transformation
!-----------------------------------------------------------------------------------------------------------------------------------
! Interpolation - for Newton localisation of particles in curved Elements
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: XiCL_NGeo(:)
REAL,ALLOCATABLE :: XiCL_NGeo1(:)
REAL,ALLOCATABLE :: XCL_NGeo(:,:,:,:,:)         !> mapping X(xi) P\in Ngeo
REAL,ALLOCATABLE :: dXCL_NGeo(:,:,:,:,:,:)      !> Jacobi matrix of the mapping P\in NGeo
REAL,ALLOCATABLE :: wBaryCL_NGeo(:)
REAL,ALLOCATABLE :: wBaryCL_NGeo1(:)
REAL,ALLOCATABLE :: DCL_NGeo(:,:)
REAL,ALLOCATABLE :: DCL_N(:,:)
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: Xi_NGeo(:)                  !> 1D equidistant point positions for curved elements (during readin)
REAL             :: DeltaXi_NGeo
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE :: Vdm_CLNGeo_CLN(:,:)
REAL,ALLOCATABLE :: Vdm_CLNGeo_GaussN(:,:)
REAL,ALLOCATABLE :: Vdm_CLN_GaussN(:,:)
REAL,ALLOCATABLE :: Vdm_CLNGeo1_CLNGeo(:,:)
REAL,ALLOCATABLE :: Vdm_NGeo_CLNGeo(:,:)

!-----------------------------------------------------------------------------------------------------------------------------------
! mapping from GaussPoints to Side or Neighbor Volume
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: VolToSideA(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: VolToSideIJKA(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: VolToSide2A(:,:,:,:,:)
INTEGER,ALLOCATABLE :: CGNS_VolToSideA(:,:,:,:,:)
INTEGER,ALLOCATABLE :: CGNS_SideToVol2A(:,:,:,:)
INTEGER,ALLOCATABLE :: SideToVolA(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: SideToVol2A(:,:,:,:,:)

END MODULE MOD_Particle_Mesh_Vars
