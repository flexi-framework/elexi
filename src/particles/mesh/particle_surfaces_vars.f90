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
#include "particle.h"

!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_Surfaces_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE,DIMENSION(:,:)         :: XiBuf
REAL,ALLOCATABLE,DIMENSION(:,:)         :: minmax
REAL,ALLOCATABLE,DIMENSION(:)           :: xiup
REAL,ALLOCATABLE,DIMENSION(:)           :: xidown
REAL,ALLOCATABLE,DIMENSION(:,:)         :: BezierControlPoints1D
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: BezierControlPoints2D        ! Bezier basis control points of degree equal to NGeo
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: BezierControlPoints2D_temp
REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: BezierControlPoints2D_temp2
REAL,ALLOCPOINT,DIMENSION(:,:,:,:)      :: BezierControlPoints3D        ! Bezier basis control points of degree equal to NGeo
REAL,ALLOCPOINT,DIMENSION(:,:,:,:)      :: BezierControlPoints3DElevated! Bezier basis control points of degree equal to NGeoElevated

REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: BiLinearCoeff                ! contains the bi-linear coefficients for each side
REAL,ALLOCPOINT,DIMENSION(:,:)          :: BaseVectors0                 ! vectors for building intersectionsurfaces for particle
                                                                        ! from Bezierpoints (1:3,1:nBCSurfaces)
REAL,ALLOCPOINT,DIMENSION(:,:)          :: BaseVectors1                 ! vectors for building intersectionsurfaces for particle
                                                                        ! from Bezierpoints (1:3,1:nBCSurfaces)
REAL,ALLOCPOINT,DIMENSION(:,:)          :: BaseVectors2                 ! vectors for building intersectionsurfaces for particle
                                                                        ! from Bezierpoints (1:3,1:nBCSurfaces)
REAL,ALLOCPOINT,DIMENSION(:,:)          :: BaseVectors3                 ! additional vector for bilinear intersection
                                                                        ! from Bezierpoints (1:3,1:nBCSurfaces)
!REAL,ALLOCPOINT,DIMENSION(:)            :: BaseVectorsScale             ! approx. size of face for bilinear intersection
!                                                                        ! from Bezierpoints (1:nBCSurfaces)

REAL,ALLOCATABLE,DIMENSION(:,:)         :: ElevationMatrix              ! array for binomial coefficients used for Bezier Elevation

REAL,ALLOCPOINT,DIMENSION(:,:,:)        :: SideSlabNormals              ! normal vectors of bounding slab box (Sides)
REAL,ALLOCPOINT,DIMENSION(:,:)          :: SideSlabIntervals            ! intervals beta1, beta2, beta3 (Sides)
LOGICAL,ALLOCPOINT,DIMENSION(:)         :: BoundingBoxIsEmpty           ! logical if Side bounding box is empty

REAL,ALLOCATABLE,DIMENSION(:,:)         :: Vdm_Bezier,sVdm_Bezier       ! Vdm from/to Bezier Polynomial from BC representation
REAL,ALLOCATABLE,DIMENSION(:,:)         :: D_Bezier                     ! D-Matrix of Bezier Polynomial from BC representation
REAL,ALLOCATABLE,DIMENSION(:,:)         :: arrayNchooseK                ! array for binomial coefficients
REAL,ALLOCATABLE,DIMENSION(:,:)         :: FacNchooseK                  ! array for binomial coefficients times prefactor

INTEGER,ALLOCPOINT,DIMENSION(:)         :: SideType                     ! integer array with side type - planar - bilinear - curved
REAL,ALLOCPOINT,DIMENSION(:,:)          :: SideNormVec                  ! normal Vector of planar sides
REAL,ALLOCPOINT,DIMENSION(:)            :: SideDistance                 ! distance of planar base from origin

INTEGER,ALLOCATABLE,DIMENSION(:)        :: gElemBCSides                 ! number of BC-Sides of element
REAL                                    :: epsilontol                   ! epsilon for setting the tolerance
REAL                                    :: OneMinusEps                  ! 1 - eps: epsilontol
REAL                                    :: OnePlusEps                   ! 1 + eps: epsilontol for setting the boundary tolerance
REAL                                    :: MinusEps                     ! - eps: epsilontol
LOGICAL                                 :: ParticleSurfaceInitIsDone=.FALSE.
! settings for Bezier-Clipping and definition of maximal number of intersections
REAL                                    :: BezierNewtonAngle            ! switch for intersection with bezier newton algorithm
                                                                        ! smallest angle of impact of particle trajectory on face
REAL                                    :: BezierClipHit                ! value for clip hit
REAL                                    :: BezierClipTolerance          ! tolerance for root of bezier clipping
REAL                                    :: BezierClipLocalTol           ! abort tolerance for bezier clipping
REAL                                    :: BezierNewtonTolerance2       ! tolerance for root of bezier Newton
INTEGER                                 :: BezierNewtonGuess            ! guess of bezier newton
INTEGER                                 :: BezierNewtonMaxIter          ! maximum iterations for bezier newton
REAL                                    :: BezierNewtonHit              ! value for bezier Newton hit
REAL                                    :: BezierSplitLimit             ! clip if remaining area after clip is > clipforce %
INTEGER                                 :: BezierClipMaxIntersec        ! maximal possible intersections for Bezier clipping
INTEGER                                 :: BezierClipMaxIter            ! maximal iterations per intersections
INTEGER                                 :: BezierClipLineVectorMethod   ! recompute method for Lu,Lv
                                                                        ! 0 - once
                                                                        ! 1 - after each clip
                                                                        ! 2 - after each xi,eta pair
INTEGER                                 :: BezierElevation              ! elevate polynomial degree to NGeo+BezierElevation
REAL,ALLOCATABLE,DIMENSION(:)           :: locAlpha,locXi,locEta        ! position of trajectory-patch
REAL,ALLOCATABLE,DIMENSION(:,:)         :: XiArray,EtaArray             ! xi and eta history for computation of intersection
!LOGICAL                                 :: MultipleBCs                  ! allow for multiple BC during one tracking step
                                                                        ! only for do-ref-mapping required
#if CODE_ANALYZE
REAL                                    :: rBoundingBoxChecks           ! number of bounding box checks
REAL(KIND=QP)                           :: rTotalBBChecks               ! total number of bounding box checks
REAL                                    :: rPerformBezierClip           ! number of performed bezier clips
REAL                                    :: rPerformBezierNewton         ! number of performed bezier newton intersections
REAL(KIND=QP)                           :: rTotalBezierClips            ! total number of performed bezier clips
REAL(KIND=QP)                           :: rTotalBezierNewton           ! total number of performed bezier newton intersections
! REAL,ALLOCATABLE,DIMENSION(:)           :: SideBoundingBoxVolume        ! Bounding Box volume
#endif /*CODE_ANALYZE*/

! Surface flux
INTEGER                                 :: BezierSampleN                ! equidistant sampling of bezier surface for emission
INTEGER                                 :: SurfFluxSideSize(2)          ! discretization of sides for Surfaceflux
REAL,ALLOCATABLE,DIMENSION(:)           :: BezierSampleXi               ! ref coordinate for equidistant bezier surface sampling
LOGICAL                                 :: TriaSurfaceFlux, WriteTriaSurfaceFluxDebugMesh

REAL,ALLOCATABLE,DIMENSION(:)           :: SurfMeshSideAreas            ! areas of of sides of surface mesh (1:nBCSides)
TYPE tSurfMeshSubSideData
  REAL                                  :: vec_nIn(3)                  ! inward directed normal of sub-sides of surface mesh
  REAL                                  :: vec_t1(3)                   ! first orth. vector in sub-sides of surface mesh
  REAL                                  :: vec_t2(3)                   ! second orth. vector in sub-sides of surface mesh
  REAL                                  :: area                        ! area of sub-sides of surface mesh
END TYPE tSurfMeshSubSideData
TYPE(tSurfMeshSubSideData),ALLOCATABLE  :: SurfMeshSubSideData(:,:,:)  ! geodata of sub-sides of surface mesh (:,:,1:nBCSides)

TYPE tTriaSwapGeo
  REAL                                  :: midpoint(3)                 ! midpoint for tria swapping in parallelogram
  REAL                                  :: ndist(3)                    ! normal vector for tria swapping in parallelogram
END TYPE tTriaSwapGeo
TYPE tTriaSideGeo
  REAL                                  :: xyzNod(3)                   ! first corner
  REAL                                  :: Vectors(3,3)                ! vectors from xyzNod to the 3 other corners of side
END TYPE tTriaSideGeo

TYPE tBCdata_auxSF
  INTEGER                               :: SideNumber                  ! Number of Particles in Sides in SurfacefluxBC
  REAL                                  :: GlobalArea, LocalArea       ! Sum of global and local tria-areas
  INTEGER                , ALLOCATABLE  :: SideList(:)                 ! List of Sides in BC (1:SideNumber)
  TYPE(tTriaSwapGeo)     , ALLOCATABLE  :: TriaSwapGeo(:,:,:)          ! data for tria-swapping in surfflux (:,:,1:SideNumber)
  TYPE(tTriaSideGeo)     , ALLOCATABLE  :: TriaSideGeo(:)              ! data for trias in surfflux (1:SideNumber)
END TYPE tBCdata_auxSF
TYPE(tBCdata_auxSF),ALLOCATABLE         :: BCdata_auxSF(:)             ! aux. data of BCs for surfacefluxes, (1:nBCs) (!!!)


!===================================================================================================================================

END MODULE MOD_Particle_Surfaces_Vars
