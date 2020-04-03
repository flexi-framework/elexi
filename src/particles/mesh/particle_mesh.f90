!=================================================================================================================================
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
!=================================================================================================================================
#include "flexi.h"
#include "particle.h"

!===================================================================================================================================
! Builds the mesh for particle tracking, separate from the DG mesh
!===================================================================================================================================
MODULE MOD_Particle_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleMeshBasis
    MODULE PROCEDURE InitParticleMeshBasis
END INTERFACE

INTERFACE InitParticleMesh
  MODULE PROCEDURE InitParticleMesh
END INTERFACE

INTERFACE InitParticleGeometry
  MODULE PROCEDURE InitParticleGeometry
END INTERFACE

INTERFACE FinalizeParticleMesh
  MODULE PROCEDURE FinalizeParticleMesh
END INTERFACE

INTERFACE GetMeshMinMax
  MODULE PROCEDURE GetMeshMinMax
END INTERFACE

INTERFACE MapRegionToElem
  MODULE PROCEDURE MapRegionToElem
END INTERFACE

INTERFACE MarkAuxBCElems
  MODULE PROCEDURE MarkAuxBCElems
END INTERFACE

INTERFACE GetGlobalNonUniqueSideID
  MODULE PROCEDURE GetGlobalNonUniqueSideID
END INTERFACE

PUBLIC :: DefineParametersParticleMesh
PUBLIC :: InitParticleMeshBasis
PUBLIC :: InitParticleMesh
PUBLIC :: InitParticleGeometry
PUBLIC :: FinalizeParticleMesh
PUBLIC :: MapRegionToElem
PUBLIC :: MarkAuxBCElems
PUBLIC :: GetMeshMinMax
PUBLIC :: GetGlobalNonUniqueSideID
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle tracking
!==================================================================================================================================
SUBROUTINE DefineParametersParticleMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('Tracking')

CALL prms%CreateLogicalOption( 'DoRefMapping'&
  , 'Refmapping [T] or Tracing [F] algorithms are used for tracking of particles.'&
  , '.TRUE.')

CALL prms%CreateLogicalOption( 'TriaTracking'&
  , 'Using Triangle-aproximation [T] or (bi-)linear and bezier (curved) description [F] of sides for tracing algorithms.'//&
  ' Requires DoRefMapping=F.'&
  ,'.FALSE.')
CALL prms%CreateLogicalOption( 'Write-Tria-DebugMesh'&
  , 'Writes per proc triangulated Surfacemesh used for Triatracking. Requires TriaTracking=T.'&
  ,'.FALSE.')
CALL prms%CreateLogicalOption( 'TriaSurfaceFlux'&
  , 'Using Triangle-aproximation [T] or (bi-)linear and bezier (curved) description [F] of sides for surfaceflux.'//&
  ' Default is set to TriaTracking')
CALL prms%CreateLogicalOption( 'Write-TriaSurfaceFlux-DebugMesh'&
  , 'Writes per proc triangulated Surfacemesh used for TriaSurfaceFlux. Requires TriaSurfaceFlux=T.'&
  ,'.FALSE.')

CALL prms%CreateLogicalOption( 'CountNbOfLostParts'&
  , 'Count number of lost particles during tracking that can not be found with fallbacks.','.FALSE.')
#if CODE_ANALYZE
CALL prms%CreateIntOption(     'PartOut'&
  , 'If compiled with CODE_ANALYZE flag: For This particle number every tracking information is written as STDOUT.','0')
CALL prms%CreateIntOption(     'MPIRankOut'&
  , 'If compiled with CODE_ANALYZE flag: This MPI-Proc writes the tracking information for the defined PartOut.','0')
#endif /*CODE_ANALYZE*/
CALL prms%CreateLogicalOption( 'MeasureTrackTime'&
  , 'If .TRUE. then the time how long the tracking routines are called are sampled and written for each MPI-Proc.','.FALSE.')
CALL prms%CreateLogicalOption( 'CartesianPeriodic'&
    , ' Simplified treatment for periodic box with Refmapping. Not computation of intersection points at periodic BCs.','.FALSE.')
CALL prms%CreateLogicalOption( 'FastPeriodic'&
  , ' Further simplification by directly moving particle into grid. Instead of moving the particle several times the periodic'//&
    ' displacements, the particle is mapped directly back into the domain. ','.FALSE.')
CALL prms%CreateIntOption(     'RefMappingGuess'&
  , ' Initial guess of the Newton for mapping the particle into reference coordinates.\n'//&
    '1 -linear pseudo-Cartesian coordinates\n'//&
    '2 - Xi of closest Gauss point\n'//&
    '3 - Xi of closest XCL_ngeo point\n'//&
    '4 -trival guess (0,0,0)^t')
CALL prms%CreateRealOption(    'RefMappingEps'&
  , ' Tolerance for mapping particle into reference element measured as L2-norm of deltaXi' , '1e-4')
CALL prms%CreateIntOption(     'BezierElevation'&
  , ' Use BezierElevation>0 to tighten the bounding box. Typical values>10','0')
CALL prms%CreateIntOption(     'BezierSampleN'&
  , 'TODO-DEFINE-PARAMETER\n'//&
    'Default value: NGeo equidistant sampling of bezier surface for emission','0')


! Background mesh init variables
CALL prms%CreateRealArrayOption('Part-FIBGMdeltas'&
  , 'Define the deltas for the cartesian Fast-Init-Background-Mesh.'//&
  ' They should be of the similar size as the smallest cells of the used mesh for simulation.'&
  , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption('Part-FactorFIBGM'&
  , 'Factor with which the background mesh will be scaled.'&
  , '1. , 1. , 1.')
CALL prms%CreateLogicalOption( 'printMPINeighborWarnings'&
    ,  ' Print warning if the MPI-Halo-region between to procs are not overlapping. Only one proc find the other in halo ' &
    ,'.FALSE.')
CALL prms%CreateLogicalOption( 'printBezierControlPointsWarnings'&
    ,  ' Print warning if MINVAL(BezierControlPoints3d(iDir,:,:,newSideID)) and global boundaries are too close ' &
    ,'.FALSE.')

! Mesh check of periodic vectors
CALL prms%CreateLogicalOption( 'Part-PeriodicCheck'&
    ,  ' Disable check of periodic vectors in case of borderline conformal mesh ' &
    ,'.TRUE.')
CALL prms%CreateLogicalOption( 'Part-PeriodicReorder'&
    ,  ' Flag if FLEXI should try to reorder particle periodic vectors ' &
    ,'.FALSE.')
CALL prms%CreateLogicalOption( 'Part-PeriodicStrict'&
    ,  ' Flag if FLEXI should abort in case no periodic vectors are found ' &
    ,'.FALSE.')
CALL prms%CreateRealOption('Part-PeriodicVecTol'&
    ,   'Tolerance angle when looking for periodic vectors ' &
    ,'0')

! Wall models
!CALL prms%CreateStringOption('Part-WallModel'&
!  , ' Wall model to be used for particle tracking Available options:.\n'//&
!    'perfRef  - perfect reflection\n'//&
!    'coeffRes - Coefficient of restitution')
!CALL prms%CreateRealOption('Part-WallCoeff-Tang'        , 'Coefficient of restituation in tangential direction', '1.0')
!CALL prms%CreateRealOption('Part-WallCoeff-Norm'        , 'Coefficient of restituation in tangential direction', '1.0')

CALL prms%CreateRealOption(    'BezierNewtonAngle'      , ' BoundingBox intersection angle for switching between '//&
'Bezierclipping and BezierNewton.' , '1.570796326')
CALL prms%CreateRealOption(    'BezierClipTolerance'    , ' Tolerance for BezierClipping' , '1e-8')
CALL prms%CreateRealOption(    'BezierNewtonTolerance'  , ' Tolerance for BezierNewton' , '1e-4')
CALL prms%CreateIntOption(     'BezierNewtonGuess'      , ' Initial guess for BezierNewton\n'// &
    '1 - linear projected face\n'//&
    '2 - closest projected BeziercontrolPoint\n'//&
    '4 - (0,0)^t' , '1')
CALL prms%CreateIntOption(     'BezierNewtonMaxIter'    , ' TODO-DEFINE-PARAMETER' , '100')
CALL prms%CreateRealOption(    'BezierSplitLimit'       , ' Limit for splitting in BezierClipping.'// &
   ' Value allows to detect multiple intersections and speed up computation. Parameter is multiplied by 2' , '0.6')
CALL prms%CreateIntOption(     'BezierClipMaxIter'      , ' Max iteration of BezierClipping' , '100')
CALL prms%CreateRealOption(    'epsilontol'             , ' Tolerance (absolute) for comparison against zero', '0.')
CALL prms%CreateRealOption(    'BezierClipHit'          , ' Tolerance in [-1,1] of BezierFace' , '0.')
CALL prms%CreateRealOption(    'BezierNewtonHit'        , ' Tolerance in [-1,1] of BezierNewton' , '0.')
CALL prms%CreateIntOption(     'BezierClipMaxIntersec'  , ' Max. number of multiple intersections. Default: 2*(NGeo+1)')

END SUBROUTINE DefineParametersParticleMesh


SUBROUTINE InitParticleMeshBasis()
!===================================================================================================================================
! Read Parameter from inputfile
!===================================================================================================================================
! MODULES
USE MOD_PreProc                ,ONLY: N
USE MOD_Basis
USE MOD_Interpolation_Vars     ,ONLY: xGP
USE MOD_Mesh_Vars
USE MOD_Particle_Basis
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N)                     :: XiCL_N,wBaryCL_N
REAL,DIMENSION(0:NGeo)                  :: wBary_NGeo
INTEGER                                 :: i
!===================================================================================================================================
ALLOCATE(DCL_N             (0:N   ,0:N)     &
        ,Vdm_CLN_GaussN    (0:N   ,0:N)     &
        ,Vdm_CLNGeo_GaussN (0:N   ,0:NGeo)  &
        ,Vdm_CLNGeo_CLN    (0:N   ,0:NGeo)  &
        ,Vdm_CLNGeo1_CLNGeo(0:NGeo,0:1)     &
        ,Vdm_NGeo_CLNGeo   (0:NGeo,0:NGeo)  &
        ,Vdm_Bezier        (0:NGeo,0:NGeo)  &
        ,sVdm_Bezier       (0:NGeo,0:NGeo)  &
        ,wBaryCL_NGeo      (0:NGeo)         &
        ,wBaryCL_NGeo1     (0:1)            &
        ,Xi_NGeo           (0:NGeo)         &
        ,XiCL_NGeo         (0:NGeo)         &
        ,XiCL_NGeo1        (0:1)            &
        ,DCL_NGeo          (0:NGeo,0:NGeo)  &
        ,D_Bezier          (0:NGeo,0:NGeo))

! Chebyshev-Lobatto N
CALL ChebyGaussLobNodesAndWeights(N,XiCL_N)
CALL BarycentricWeights          (N,XiCL_N,wBaryCL_N)
CALL PolynomialDerivativeMatrix  (N,XiCL_N,DCL_N)
CALL InitializeVandermonde       (N,N,wBaryCL_N,XiCL_N,xGP,Vdm_CLN_GaussN)

!equidistant-Lobatto NGeo
DO i=0,NGeo
  Xi_NGeo(i) = 2./REAL(NGeo) * REAL(i) - 1.
END DO
DeltaXi_NGeo = 2./NGeo
CALL BarycentricWeights          (NGeo,Xi_NGeo,wBary_NGeo)

! Chebyshev-Lobatto NGeo
CALL ChebyGaussLobNodesAndWeights(NGeo,XiCL_NGeo)
CALL BarycentricWeights          (NGeo,XiCL_NGeo,wBaryCL_NGeo)
CALL PolynomialDerivativeMatrix  (NGeo,XiCL_NGeo,DCL_NGeo)

CALL InitializeVandermonde       (NGeo,N   ,wBaryCL_NGeo,XiCL_NGeo,xGP      ,Vdm_CLNGeo_GaussN)
CALL InitializeVandermonde       (NGeo,N   ,wBaryCL_NGeo,XiCL_NGeo,XiCL_N   ,Vdm_CLNGeo_CLN   )
CALL InitializeVandermonde       (NGeo,NGeo,wBary_NGeo  ,Xi_NGeo  ,XiCL_NGeo,Vdm_NGeo_CLNGeo  )

! small wBaryCL_NGEO
ALLOCATE                         (wBaryCL_NGeo1(0:1),XiCL_NGeo1(0:1))
CALL ChebyGaussLobNodesAndWeights(1                 ,XiCL_NGeo1)
CALL BarycentricWeights          (1                 ,XiCL_NGeo1,wBaryCL_NGeo1)
CALL InitializeVandermonde       (1, NGeo,wBaryCL_NGeo1,XiCL_NGeo1,XiCL_NGeo ,Vdm_CLNGeo1_CLNGeo)

! initialize Vandermonde for Bezier basis surface representation (particle tracking with curved elements)
CALL BuildBezierVdm              (NGeo,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier)
CALL BuildBezierDMat             (NGeo,Xi_NGeo,D_Bezier)

! allocate Chebyshev-Lobatto physical coordinates
ALLOCATE( XCL_NGeo(1:3,    0:Ngeo,0:Ngeo,0:ZDIM(NGeo),1:nElems)   &
        ,dXCL_NGeo(1:3,1:3,0:Ngeo,0:Ngeo,0:ZDIM(NGeo),1:nElems))
XCL_NGeo  = 0.
dXCL_NGeo = 0.

END SUBROUTINE InitParticleMeshBasis


SUBROUTINE InitParticleMesh()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: NGeo,nElems,useCurveds
USE MOD_Particle_BGM           ,ONLY: BuildBGMAndIdentifyHaloRegion
USE MOD_Particle_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: InitGetGlobalElemID,InitGetCNElemID
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation,BezierControlPoints3D,BezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping,FastPeriodic,CountNbOfLostParts,nLostParts,CartesianPeriodic
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GetRealArray
#if CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: SideBoundingBoxVolume
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif /* USE_MPI */
#if USE_LOADBALANCE
USE MOD_Particle_Tracking_Vars ,ONLY: MeasureTrackTime
#else
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: RefMappingGuessProposal
INTEGER          :: iSample
INTEGER          :: firstSide,lastSide,iSide
CHARACTER(LEN=2) :: hilf
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#else
INTEGER          :: ALLOCSTAT
#endif
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH ...'
IF(ParticleMeshInitIsDone) CALL ABORT(__STAMP__, ' Particle-Mesh is already initialized.')

!===================================================================================================================================
! Collect particle variables that are not initialized somewhere else
!===================================================================================================================================
PP_nElems = nElems

#if USE_LOADBALANCE
! ElemTime already set in loadmesh
#else
SDEALLOCATE(ElemTime)
ALLOCATE(ElemTime(1:nElems))
#endif
!===================================================================================================================================

DoRefMapping       = GETLOGICAL('DoRefMapping',".TRUE.")
TriaTracking       = GETLOGICAL('TriaTracking','.FALSE.')

! Build BGM to Element mapping and identify which of the elements, sides and nodes are in the compute-node local and halo region
CALL BuildBGMAndIdentifyHaloRegion()

! Initialize mapping function: GetGlobalElemID()
CALL InitGetGlobalElemID()

! Initialize mapping function: GetCNElemID()
CALL InitGetCNElemID()

IF ((DoRefMapping.OR.UseCurveds.OR.(NGeo.GT.1)).AND.(TriaTracking)) THEN
  CALL abort(&
__STAMP__&
,'DoRefMapping=T .OR. UseCurveds=T .OR. NGEO>1! Not possible with TriaTracking=T at the same time!')
END IF
CountNbOfLostParts = GETLOGICAL('CountNbOfLostParts',".FALSE.")
nLostParts         = 0

PeriodicCheck      = GETLOGICAL('Part-PeriodicCheck',".TRUE.")

#if CODE_ANALYZE
PARTOUT            = GETINT('PartOut','0')
MPIRankOut         = GETINT('MPIRankOut','0')
#endif /*CODE_ANALYZE*/

! Particle wall model
!PartWallModel      = GETSTR('Part-WallModel','perfRef')
!PartWallTang       = GETREAL('Part-WallCoeff-Tang', '1.0')
!PartWallNorm       = GETREAL('Part-WallCoeff-Norm', '1.0')

!IF(.NOT.DoRefMapping) THEN
!  SDEALLOCATE(nTracksPerElem)
!END IF
CartesianPeriodic = GETLOGICAL('CartesianPeriodic','.FALSE.')
IF(CartesianPeriodic) FastPeriodic = GETLOGICAL('FastPeriodic','.FALSE.')

! method from xPhysic to parameter space

IF(UseCurveds)THEN ! don't use RefMappingGuess=1, because RefMappingGuess is only best for linear cubical elements
  ! curved elements can be stronger deformed, hence, a better guess can be used
  ! RefMappingGuess 2,3 searches the closest Gauss/CL points of the considered element. This point is used as the initial value for
  ! the mapping. Note, that the position of the CL points can still be advantageous for the initial guess.
  RefMappingGuessProposal=2
  IF(PP_N.GT.NGeo)THEN ! there are more Gauss points within an element then CL-points
                       ! Gauss points sample the element finer
                       ! Note: the Gauss points does not exist for HALO elements, here, the closest CL point is used
    RefMappingGuessProposal=2
  ELSE ! more CL-points than Gauss points, hence, better sampling of the element
    RefMappingGuessProposal=3
  END IF
ELSE
  RefMappingGuessProposal=1 ! default for linear meshes. Guess is exact for cubical, non-twisted elements
END IF
WRITE(hilf,'(I2.2)') RefMappingGuessProposal
RefMappingGuess = GETINT('RefMappingGuess',hilf)
IF((RefMappingGuess.LT.1).AND.(UseCurveds)) THEN ! this might cause problems
  SWRITE(UNIT_stdOut,'(A)')' WARNING: read-in [RefMappingGuess=1] when using [UseCurveds=T] may create problems!'
END IF
RefMappingEps   = GETREAL('RefMappingEps','1e-4')

epsInCell       = SQRT(3.0*RefMappingEps)
!epsOneCell      = 1.0+epsInCell

IF((RefMappingGuess.LT.1).OR.(RefMappingGuess.GT.4))THEN
   CALL abort(&
__STAMP__ &
,'Wrong guessing method for mapping from physical space in reference space.',RefMappingGuess,999.)
END IF

IF (TriaTracking) THEN
  CALL InitParticleGeometry()
  ! Compute convex element radius^2
  CALL BuildElementRadiusTria()
ELSE
  CALL CalcParticleMeshMetrics()

  BezierElevation = GETINT('BezierElevation','0')
  NGeoElevated    = NGeo + BezierElevation
  CALL CalcBezierControlPoints()

#if USE_MPI
  MPISharedSize = INT((3**2*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,3,nComputeNodeTotalSides/),SideSlabNormals_Shared_Win,SideSlabNormals_Shared)
  CALL MPI_WIN_LOCK_ALL(0,SideSlabNormals_Shared_Win,IERROR)
  MPISharedSize = INT((6*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalSides/),SideSlabIntervals_Shared_Win,SideSlabIntervals_Shared)
  CALL MPI_WIN_LOCK_ALL(0,SideSlabIntervals_Shared_Win,IERROR)
  MPISharedSize = INT((nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalSides/),BoundingBoxIsEmpty_Shared_Win,BoundingBoxIsEmpty_Shared)
  CALL MPI_WIN_LOCK_ALL(0,BoundingBoxIsEmpty_Shared_Win,IERROR)
  firstSide = INT(REAL (myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
  lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
  SideSlabNormals    => SideSlabNormals_Shared
  SideSlabIntervals  => SideSlabIntervals_Shared
  BoundingBoxIsEmpty => BoundingBoxIsEmpty_Shared
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
  ALLOCATE(SideSlabNormals(1:3,1:3,1:nNonUniqueGlobalSides) &
          ,SideSlabIntervals(  1:6,1:nNonUniqueGlobalSides) &
          ,BoundingBoxIsEmpty(     1:nNonUniqueGlobalSides) &
          ,STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate SideMetrics arrays!')
  firstSide = 1
  lastSide  = nNonUniqueGlobalSides
#endif /* USE_MPI */
#if CODE_ANALYZE
  ALLOCATE(SideBoundingBoxVolume(nSides))
#endif /*CODE_ANALYZE*/

  IF (BezierElevation.GT.0) THEN
    DO iSide=firstSide,LastSide
      ! Ignore small mortar sides attached to big mortar sides
      IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE
      CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide)     &
                                         ,SideSlabNormals(   1:3,1:3,iSide)                          &
                                         ,SideSlabInterVals( 1:6    ,iSide)                          &
                                         ,BoundingBoxIsEmpty(iSide))
    END DO
  ELSE
    DO iSide=firstSide,LastSide
      ! Ignore small mortar sides attached to big mortar sides
      IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE
      CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)             &
                                         ,SideSlabNormals(   1:3,1:3,iSide)                          &
                                         ,SideSlabInterVals( 1:6    ,iSide)                          &
                                         ,BoundingBoxIsEmpty(iSide))
    END DO
  END IF
#if USE_MPI
  CALL MPI_WIN_SYNC(SideSlabNormals_Shared_Win,IERROR)
  CALL MPI_WIN_SYNC(SideSlabIntervals_Shared_Win,IERROR)
  CALL MPI_WIN_SYNC(BoundingBoxIsEmpty_Shared_Win,IERROR)
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI */
#if CODE_ANALYZE
  ! TODO: bounding box volumes must be calculated for all unique sides.
  offsetSideID = ElemInfo_Shared(SideIf
  DO iSide=offsetMPISides_YOUR,LastSide
    dx=ABS(SideSlabIntervals(2)-SideSlabIntervals(1))
    dy=ABS(SideSlabIntervals(4)-SideSlabIntervals(3))
    dz=ABS(SideSlabIntervals(6)-SideSlabIntervals(5))
    SideID = SideInfo
    SideBoundingBoxVolume(SideID)=dx*dy*dz
  END DO
#endif /*CODE_ANALYZE*/

  ! Compute element bary and element radius for node elements (with halo region)
  CALL BuildElementOriginShared()

  ! Check the side type (planar, bilinear, curved)
  CALL IdentifyElemAndSideType()

  ! Compute the element XiEtaZetaBasis and the radius of the convex hull
  CALL BuildElementBasisAndRadius()

  ! Get basevectors for (bi-)linear sides
  CALL GetLinearSideBaseVectors()

  IF (DoRefMapping) THEN
    ! Identify BCElems
    CALL BuildBCElemDistance()
  END IF

  CALL BuildEpsOneCell()
END IF


! BezierAreaSample stuff:
WRITE(hilf,'(L1)') TriaTracking
TriaSurfaceFlux = GETLOGICAL('TriaSurfaceFlux',TRIM(hilf))
IF (TriaSurfaceFlux) THEN
  BezierSampleN = 1
  SurfFluxSideSize=(/1,2/)
ELSE
  WRITE(hilf,'(I2.2)') NGeo
  BezierSampleN = GETINT('BezierSampleN',hilf)
  SurfFluxSideSize=BezierSampleN
  ALLOCATE(BezierSampleXi(0:BezierSampleN))!,STAT=ALLOCSTAT)
  DO iSample=0,BezierSampleN
    BezierSampleXi(iSample)=-1.+2.0/BezierSampleN*iSample
  END DO
END IF

ParticleMeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMesh


SUBROUTINE CalcParticleMeshMetrics()
!===================================================================================================================================
!> calculates XCL_Ngeo and dXCL_Ngeo for compute node mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Basis                  ,ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights,InitializeVandermonde
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Interpolation          ,ONLY: GetDerivativeMatrix
USE MOD_Interpolation_Vars     ,ONLY: NodeTypeCL
USE MOD_Mesh_Vars              ,ONLY: NGeo,InterpolateFromTree
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Array,dXCL_NGeo_Array
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo,dXCL_NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared,dXCL_NGeo_Shared
USE MOD_Particle_MPI_Shared_Vars,ONLY:GlobalElem2CNTotalElem
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared_Win,dXCL_NGeo_Shared_Win
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                        :: iElem
REAL                           :: Vdm_NGeo_CLNGeo(0:NGeo,0:NGeo)
REAL                           :: DCL_NGeo(0:Ngeo,0:Ngeo)
INTEGER                        :: firstHaloElem,lastHaloElem,nComputeNodeHaloElems
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: firstNodeID,nodeID,i,j,k,ll
REAL                           :: NodeCoordstmp(1:3,0:NGeo,0:NGeo,0:NGeo)
!===================================================================================================================================
#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
MPISharedSize = INT((3*(NGeo+1)**3*nComputeNodeElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*  (NGeo+1)*(NGeo+1)*(NGeo+1)*nComputeNodeElems/), XCL_NGeo_Shared_Win,XCL_NGeo_Array)
MPISharedSize = INT((3*3*(NGeo+1)**3*nComputeNodeElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*3*(NGeo+1)*(NGeo+1)*(NGeo+1)*nComputeNodeElems/),dXCL_NGeo_Shared_Win,dXCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,dXCL_NGeo_Shared_Win,IERROR)
XCL_NGeo_Shared (1:3    ,0:NGeo,0:NGeo,0:NGeo,1:nComputeNodeElems) => XCL_NGeo_Array
dXCL_NGeo_Shared(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nComputeNodeElems) => dXCL_NGeo_Array
! Copy local XCL and dXCL into shared
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  DO iElem = 1, nElems
    XCL_NGeo_Shared (:  ,:,:,:,offsetElem+iElem) = XCL_NGeo (:  ,:,:,:,iElem)
    dXCL_NGeo_Shared(:,:,:,:,:,offsetElem+iElem) = dXCL_NGeo(:,:,:,:,:,iElem)
  END DO ! iElem = 1, nElems
ELSE
  DO iElem = 1, nElems
    XCL_NGeo_Shared (:,  :,:,:,GlobalElem2CNTotalElem(offsetElem+iElem)) = XCL_NGeo (:,  :,:,:,iElem)
    dXCL_NGeo_Shared(:,:,:,:,:,GlobalElem2CNTotalElem(offsetElem+iElem)) = dXCL_NGeo(:,:,:,:,:,iElem)
  END DO ! iElem = 1, nElems
END IF
nComputeNodeHaloElems = nComputeNodeTotalElems - nComputeNodeElems
IF (nComputeNodeHaloElems.GT.nComputeNodeProcessors) THEN
  firstHaloElem = INT(REAL( myComputeNodeRank   *nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))+1
  lastHaloElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))
ELSE
  firstHaloElem = myComputeNodeRank + 1
  IF (myComputeNodeRank.LT.nComputeNodeHaloElems) THEN
    lastHaloElem = myComputeNodeRank + 1
  ELSE
    lastHaloElem = 0
  END IF
END IF

! Build XCL and dXCL for compute node halo region (each proc of compute-node build only its fair share)
IF(interpolateFromTree) THEN
  CALL abort(__STAMP__,'ERROR: InterpolateFromTree not yet implemented for new halo region!')
ELSE
  CALL GetDerivativeMatrix(NGeo  , NodeTypeCL  , DCL_Ngeo)

  DO iElem = firstHaloElem, lastHaloElem
    firstNodeID=ElemInfo_Shared(ELEM_FIRSTNODEIND,GetGlobalElemID(nComputeNodeElems+iElem))+1
    nodeID = 0
    DO i = 0, NGeo
      DO j = 0, NGeo
        DO k = 0, NGeo
          NodeCoordstmp(:,i,j,k) = NodeCoords_Shared(1,firstNodeID+NodeID)
          nodeID = nodeID + 1
        END DO
      END DO
    END DO ! i = 0, NGeo
    CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_NGeo_CLNGeo,NodeCoordstmp,XCL_NGeo_Shared(:,:,:,:,nComputeNodeElems+iElem))

    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      ! Matrix-vector multiplication
      DO ll=0,Ngeo
        dXCL_NGeo_Shared(1,:,i,j,k,iElem) = dXCL_NGeo_Shared(1,:,i,j,k,iElem) + DCL_NGeo(i,ll)*XCL_NGeo_Shared(:,ll,j,k,iElem)
        dXCL_NGeo_Shared(2,:,i,j,k,iElem) = dXCL_NGeo_Shared(2,:,i,j,k,iElem) + DCL_NGeo(j,ll)*XCL_NGeo_Shared(:,i,ll,k,iElem)
        dXCL_NGeo_Shared(3,:,i,j,k,iElem) = dXCL_NGeo_Shared(3,:,i,j,k,iElem) + DCL_NGeo(k,ll)*XCL_NGeo_Shared(:,i,j,ll,iElem)
      END DO !l=0,N
    END DO; END DO; END DO !i,j,k=0,Ngeo
    END DO ! iElem = firstHaloElem, lastHaloElem
END IF

CALL MPI_WIN_SYNC(XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(dXCL_NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ALLOCATE(XCL_NGeo_Shared (3,  0:NGeo,0:NGeo,0:NGeo,nElems))
ALLOCATE(dXCL_NGeo_Shared(3,3,0:NGeo,0:NGeo,0:NGeo,nElems))
XCL_NGeo_Shared  = XCL_NGeo
dXCL_NGeo_Shared = dXCL_NGeo
#endif /*USE_MPI*/

END SUBROUTINE CalcParticleMeshMetrics


SUBROUTINE CalcBezierControlPoints()
!===================================================================================================================================
!> calculate the bezier control point (+elevated) for shared compute-node mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_Mappings               ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: NGeoElevated,XCL_NGeo_Shared
USE MOD_Particle_Surfaces      ,ONLY: GetBezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,sVdm_Bezier
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3DElevated,BezierElevation
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID
#endif
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iSide,ilocSide
INTEGER                        :: firstElem,lastElem,firstSide,lastSide
INTEGER                        :: p,q,pq(2), SideID
REAL                           :: tmp(3,0:NGeo,0:NGeo)
REAL                           :: tmp2(3,0:NGeo,0:NGeo)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#else
INTEGER                        :: ALLOCSTAT
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' CALCULATING BezierControlPoints ...'

! Build BezierControlPoints3D (compute-node local+halo)
#if USE_MPI
MPISharedSize = INT((3*(NGeo+1)**2*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared(MPISharedSize,(/3*(NGeo+1)*(NGeo+1)*nComputeNodeTotalSides/),BezierControlPoints3D_Shared_Win,BezierControlPoints3D_Shared)
CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3D_Shared_Win,IERROR)
BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nComputeNodeTotalSides) => BezierControlPoints3D_Shared
IF (myComputeNodeRank.EQ.0) THEN
  BezierControlPoints3D         = 0.
END IF
IF (BezierElevation.GT.0) THEN
  MPISharedSize = INT((3*(NGeoElevated+1)**2*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3*(NGeoElevated+1)*(NGeoElevated+1)*nComputeNodeTotalSides/), &
                                      BezierControlPoints3DElevated_Shared_Win,BezierControlPoints3DElevated_Shared)
  CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3DElevated_Shared_Win,IERROR)
  BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nComputeNodeTotalSides) => BezierControlPoints3DElevated_Shared
  IF (myComputeNodeRank.EQ.0) THEN
    BezierControlPoints3DElevated = 0.
  END IF
END IF
#else
ALLOCATE(BezierControlPoints3D(1:3,1:NGeo+1,1:NGeo+1,1:nNonUniqueGlobalSides) &
        ,STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3D!')
BezierControlPoints3D         = 0.

IF (BezierElevation.GT.0) THEN
  ALLOCATE(BezierControlPoints3DElevated(1:3,1:NGeoElevated+1,1:NGeoElevated+1,1:nNonUniqueGlobalSides) &
          ,STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3DElevated!')
  BezierControlPoints3DElevated = 0.
END IF
#endif

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
firstSide = INT(REAL (myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif

DO iElem = firstElem, lastElem
  DO ilocSide=1,6
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , 0    , :    , :   ,iElem )
    CASE(XI_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , NGeo , :    , :   ,iElem )
    CASE(ETA_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , 0    , :   ,iElem )
    CASE(ETA_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , NGeo , :   ,iElem )
    CASE(ZETA_MINUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , :    , 0   ,iElem )
    CASE(ZETA_PLUS)
      tmp=XCL_NGeo_Shared(1:3 , :    , :    , NGeo,iElem )
    END SELECT
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,tmp,tmp2)
    ! get global SideID of local side
    SideID = GetGlobalNonUniqueSideID(GetGlobalElemID(iElem),iLocSide)
    DO q=0,NGeo; DO p=0,NGeo
      ! turn into right hand system of side
      pq=CGNS_SideToVol2(NGeo,p,q,iLocSide,3)
      BezierControlPoints3D(1:3,p,q,SideID)=tmp2(1:3,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! ilocSide=1,6
END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

! Calculate elevated BezierControlPoints
IF (BezierElevation.GT.0) THEN
  DO iSide=firstSide,LastSide
    ! Ignore small mortar sides attached to big mortar sides
    IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).LT.6) CYCLE
    ! Indices in shared arrays are shifted by 1
    CALL GetBezierControlPoints3DElevated( NGeo,NGeoElevated                                                       &
                                         , BezierControlPoints3D        (1:3,0:NGeo        ,0:NGeo        ,iSide)  &
                                         , BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide))
  END DO

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3DElevated_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif
END IF

END SUBROUTINE CalcBezierControlPoints


SUBROUTINE InitParticleGeometry()
!===================================================================================================================================
! Subroutine for particle geometry initialization (GEO container)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
#else
USE MOD_Mesh_Vars              ,ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER            :: ALLOCSTAT
INTEGER            :: iElem,iLocSide,iNode
!INTEGER            :: jNode
INTEGER            :: nStart, NodeNum
INTEGER            :: NodeMap(4,6)
!INTEGER            :: nSides
REAL               :: A(3,3),detcon
!REAL,ALLOCATABLE   :: Coords(:,:,:,:)
!CHARACTER(32)      :: hilf
!CHARACTER(LEN=255) :: FileString
INTEGER            :: FirstElem, LastElem, GlobalSideID, nlocSides, localSideID
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
INTEGER            :: CornerNodeIDswitch(8)
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE GEOMETRY INFORMATION...'


! the cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
CornerNodeIDswitch(1)=1
CornerNodeIDswitch(2)=(Ngeo+1)
CornerNodeIDswitch(3)=(Ngeo+1)**2
CornerNodeIDswitch(4)=(Ngeo+1)*Ngeo+1
CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2
CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1

! New crazy corner node switch (philipesque)
ASSOCIATE(CNS => CornerNodeIDswitch )
  ! CGNS Mapping
  NodeMap(:,1)=(/CNS(1),CNS(4),CNS(3),CNS(2)/)
  NodeMap(:,2)=(/CNS(1),CNS(2),CNS(6),CNS(5)/)
  NodeMap(:,3)=(/CNS(2),CNS(3),CNS(7),CNS(6)/)
  NodeMap(:,4)=(/CNS(3),CNS(4),CNS(8),CNS(7)/)
  NodeMap(:,5)=(/CNS(1),CNS(5),CNS(8),CNS(4)/)
  NodeMap(:,6)=(/CNS(5),CNS(6),CNS(7),CNS(8)/)

!ALLOCATE(GEO%ElemToNodeID(1:8,1:nElems),       &
!         GEO%ElemSideNodeID(1:4,1:6,1:nElems), &
!         GEO%NodeCoords(1:3,1:nNodes),         &
!         GEO%ConcaveElemSide(1:6,1:nElems), &
!         GEO%ElemMidPoint(1:3,nElems), STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) THEN
! CALL abort(__STAMP__&
! ,'ERROR in InitParticleGeometry: Cannot allocate GEO%... stuff!')
!END IF

#if USE_MPI
MPISharedSize = INT(6*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),ConcaveElemSide_Shared_Win,ConcaveElemSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,ConcaveElemSide_Shared_Win,IERROR)
firstElem = INT(REAL(myComputeNodeRank*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))

MPISharedSize = INT(8*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/8,nComputeNodeTotalElems/),ElemNodeID_Shared_Win,ElemNodeID_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemNodeID_Shared_Win,IERROR)

MPISharedSize = INT(4*6*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/4,6,nComputeNodeTotalElems/),ElemSideNodeID_Shared_Win,ElemSideNodeID_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemSideNodeID_Shared_Win,IERROR)

MPISharedSize = INT(3*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemMidPoint_Shared_Win,ElemMidPoint_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemMidPoint_Shared_Win,IERROR)
#else
ALLOCATE(ConcaveElemSide_Shared(   1:6,1:nElems))
ALLOCATE(ElemNodeID_Shared(        1:8,1:nElems))
ALLOCATE(ElemSideNodeID_Shared(1:4,1:6,1:nElems))
ALLOCATE(ElemMidPoint_Shared(      1:3,1:nElems))
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

!ALLOCATE(GEO%ElemToNodeIDGlobal(1:8,1:nElems))

!GEO%ElemToNodeID(:,:)=0
!GEO%ElemSideNodeID(:,:,:)=0
!GEO%NodeCoords(:,:)=0.
!GEO%ConcaveElemSide(:,:)=.FALSE.
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif
  ElemNodeID_Shared      = 0
  ElemSideNodeID_Shared  = 0
  ConcaveElemSide_Shared = .FALSE.
#if USE_MPI
END IF
#endif

!iNode=0
!DO iElem=1,nElems
!  DO jNode=1,8
!    Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=0
!  END DO
!END DO
!DO iElem=1,nElems
!  !--- Save corners of sides
!  DO jNode=1,8
!    IF (Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID.EQ.0) THEN
!      iNode=iNode+1
!      Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=iNode
!      GEO%NodeCoords(1:3,iNode)=Elems(iElem+offsetElem)%ep%node(jNode)%np%x(1:3)
!    END IF
!    GEO%ElemToNodeID(jNode,iElem)=Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID
!    !GEO%ElemToNodeIDGlobal(jNode,iElem) = Elems(iElem+offsetElem)%ep%node(jNode)%np%ind
!  END DO
!END DO
!
!DO iElem=1,nElems
!  DO iLocSide=1,6
!    nStart=MAX(0,ElemToSide(E2S_FLIP,iLocSide,iElem)-1)
!    GEO%ElemSideNodeID(1:4,iLocSide,iElem)=(/Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart  ,4)+1,iLocSide))%np%NodeID,&
!                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+1,4)+1,iLocSide))%np%NodeID,&
!                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+2,4)+1,iLocSide))%np%NodeID,&
!                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+3,4)+1,iLocSide))%np%NodeID/)
!  END DO
!END DO
DO iElem = firstElem,lastElem
  DO iNode = 1,8
    ElemNodeID_Shared(iNode,iElem) = ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) + CNS(iNode)
  END DO

  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  DO iLocSide = 1,nlocSides
    ! Get global SideID
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
    IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    ! Find start of CGNS mapping from flip
    nStart = MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,GlobalSideID),10)-1)
!    IF(iElem.EQ.1) THEN
!      SWRITE(*,*) GlobalSideID
!      SWRITE(*,*) MOD(SideInfo_Shared(SIDE_FLIP,GlobalSideID),10)
!    END IF
    ! Shared memory array starts at 1, but NodeID at 0
    ElemSideNodeID_Shared(1:4,localSideID,iElem) = (/ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart  ,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart+1,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart+2,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart+3,4)+1,localSideID)-1/)
  END DO
END DO
END ASSOCIATE

!--- Save whether Side is concave or convex
DO iElem = firstElem,lastElem
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
  DO iLocSide = 1,nlocSides
    !--- Check whether the bilinear side is concave
    !--- Node Number 4 and triangle 1-2-3
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
    IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    DO NodeNum = 1,3               ! for all 3 nodes of triangle
!      A(:,NodeNum) = GEO%NodeCoords(:,GEO%ElemSideNodeID(NodeNum,iLocSide,iElem)) &
!                   - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,iLocSide,iElem))
       A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,localSideID,iElem)+1) &
                    - NodeCoords_Shared(:,ElemSideNodeID_Shared(4      ,localSideID,iElem)+1)
    END DO
    !--- concave if detcon < 0:
    detcon = ((A(2,1) * A(3,2) - A(3,1) * A(2,2)) * A(1,3) +     &
              (A(3,1) * A(1,2) - A(1,1) * A(3,2)) * A(2,3) +     &
              (A(1,1) * A(2,2) - A(2,1) * A(1,2)) * A(3,3))
!    IF (detcon.LT.0) GEO%ConcaveElemSide(iLocSide,iElem)=.TRUE.
    IF (detcon.LT.0) ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
  END DO
END DO
!SWRITE(*,*) ConcaveElemSide_Shared(:,25)
!SideInfo_Shared(SIDE_NBELEMID,ElemInfo_Shared(ELEM_FIRSTSIDEIND,1)+1)

DO iElem = firstElem,lastElem
  ElemMidPoint_Shared(:,iElem) = 0.
  DO iNode = 1,8
    ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+iNode)
  END DO
  ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) / 8.
END DO

!--- check for elements with intersecting sides (e.g. very flat elements)
CALL WeirdElementCheck()

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GEOMETRY INFORMATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleGeometry


SUBROUTINE WeirdElementCheck()
!===================================================================================================================================
! Calculate whether element edges intersect other sides
! If this is the case it means that part of the element is turned inside-out
! which results in a warning so the user can decide whether it is a problem that
! necessitates a new mesh.
! Fixing the problem would involve defining the bilinear edge between nodes 2 and 4
! (instead of 1 and 3). This information would need to be stored and used throughout
! the particle treatment. Additionally, since the edge would need to be changed
! for both neighboring elements, it is possible that both element might have the problem
! hence no solution exists.
! tl;dr: Hard/maybe impossible to fix, hence only a warning is given so the user can decide
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Vars ,ONLY: WeirdElems
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, iLocSide, kLocSide, iNode
INTEGER,ALLOCATABLE :: WeirdElemNbrs(:)
REAL              :: vec(1:3), Node(1:3,1:4),det(1:3)
LOGICAL           :: WEIRD, TRICHECK, TRIABSCHECK
INTEGER           :: firstElem,lastElem
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' CHECKING FOR WEIRD ELEMENTS...'

#if USE_MPI
firstElem = INT(REAL(myComputeNodeRank*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

ALLOCATE(WeirdElemNbrs(1:lastElem-firstElem+1))

WeirdElems = 0
DO iElem = firstElem,lastElem ! go through all elements
  WEIRD = .FALSE.
  DO iLocSide = 1,5  ! go through local sides
    IF (.not.WEIRD) THEN  ! if one is found there is no need to continue
      IF (ConcaveElemSide_Shared(iLocSide,iElem)) THEN  ! only concave elements need to be checked
        ! build vector from node 1 to node 3
        vec(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3,iLocSide,iElem)+1) &
               - NodeCoords_Shared(:,ElemSideNodeID_Shared(1,iLocSide,iElem)+1)
        ! check all other sides
        DO kLocSide = iLocSide + 1, 6
          IF (ConcaveElemSide_Shared(kLocSide,iElem)) THEN  ! only concave elements need to be checked
            ! build 4 vectors from point 1 of edge to 4 nodes of kLocSide
            DO iNode = 1,4
              Node(:,iNode) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1    ,iLocSide,iElem)+1) &
                            - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
            END DO
            ! Compute whether any of the triangle intersects with the vector vec:
            ! If all three volumes built by the vector vec and the vectors Node
            ! are either positive or negative then there is an intersection

            ! Triangle 1 (Nodes 1,2,3)
            ! Only check this if neither point of vec is part of the triangle.
            ! If points of vec correspont to point 1 or 3 or triangle then both
            ! triangles can be skipped (triabscheck = true), else point 4 needs to be checked
            ! separately for triangle 2 (see below)
            TRICHECK = .FALSE.
            TRIABSCHECK = .FALSE.
            DO iNode = 1,3
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1    ,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) THEN
                TRICHECK = .TRUE.
                IF(iNode.NE.2)TRIABSCHECK = .TRUE.
              END IF
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3    ,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(iNode,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) THEN
                TRICHECK = .TRUE.
                IF(iNode.NE.2)TRIABSCHECK = .TRUE.
              END IF
            END DO
            IF (.not.TRICHECK) THEN
              det(1) = ((Node(2,1) * Node(3,2) - Node(3,1) * Node(2,2)) * vec(1)  + &
                        (Node(3,1) * Node(1,2) - Node(1,1) * Node(3,2)) * vec(2)  + &
                        (Node(1,1) * Node(2,2) - Node(2,1) * Node(1,2)) * vec(3))
              det(2) = ((Node(2,2) * Node(3,3) - Node(3,2) * Node(2,3)) * vec(1)  + &
                        (Node(3,2) * Node(1,3) - Node(1,2) * Node(3,3)) * vec(2)  + &
                        (Node(1,2) * Node(2,3) - Node(2,2) * Node(1,3)) * vec(3))
              det(3) = ((Node(2,3) * Node(3,1) - Node(3,3) * Node(2,1)) * vec(1)  + &
                        (Node(3,3) * Node(1,1) - Node(1,3) * Node(3,1)) * vec(2)  + &
                        (Node(1,3) * Node(2,1) - Node(2,3) * Node(1,1)) * vec(3))
              IF ((det(1).LT.0).AND.(det(2).LT.0).AND.(det(3).LT.0)) WEIRD = .TRUE.
              IF ((det(1).GT.0).AND.(det(2).GT.0).AND.(det(3).GT.0)) WEIRD = .TRUE.
            END IF

            ! Triangle 2 (Nodes 1,3,4)
            TRICHECK = .FALSE.
            IF (.not.TRIABSCHECK) THEN
              ! Node 4 needs to be checked separately (see above)
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(1,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(4,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) TRICHECK = .TRUE.
              det(:) = NodeCoords_Shared(:,ElemSideNodeID_Shared(3,iLocSide,iElem)+1) &
                     - NodeCoords_Shared(:,ElemSideNodeID_Shared(4,kLocSide,iElem)+1)
              IF (SUM(abs(det(:))).EQ.0) TRICHECK = .TRUE.
              IF (.not.TRICHECK) THEN
                det(1) = ((Node(2,1) * Node(3,3) - Node(3,1) * Node(2,3)) * vec(1)  + &
                          (Node(3,1) * Node(1,3) - Node(1,1) * Node(3,3)) * vec(2)  + &
                          (Node(1,1) * Node(2,3) - Node(2,1) * Node(1,3)) * vec(3))
                det(2) = ((Node(2,3) * Node(3,4) - Node(3,3) * Node(2,4)) * vec(1)  + &
                          (Node(3,3) * Node(1,4) - Node(1,3) * Node(3,4)) * vec(2)  + &
                          (Node(1,3) * Node(2,4) - Node(2,3) * Node(1,4)) * vec(3))
                det(3) = ((Node(2,4) * Node(3,1) - Node(3,4) * Node(2,1)) * vec(1)  + &
                          (Node(3,4) * Node(1,1) - Node(1,4) * Node(3,1)) * vec(2)  + &
                          (Node(1,4) * Node(2,1) - Node(2,4) * Node(1,1)) * vec(3))
                IF ((det(1).LT.0).AND.(det(2).LT.0).AND.(det(3).LT.0)) WEIRD = .TRUE.
                IF ((det(1).GT.0).AND.(det(2).GT.0).AND.(det(3).GT.0)) WEIRD = .TRUE.
              END IF
            END IF
          END IF
        END DO
      END IF
    END IF
  END DO
  IF (WEIRD) THEN
    SWRITE(*,*) iElem
    DO iNode=1,8
      SWRITE(*,*) NodeCoords_Shared(:,ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+iNode)
    END DO
    EXIT

    WeirdElems = WeirdElems + 1
    WeirdElemNbrs(WeirdElems) = iElem
  END IF
END DO
!STOP

SWRITE(UNIT_stdOut,'(A)')' CHECKING FOR WEIRD ELEMENTS DONE!'
IF(WeirdElems.GT.0) THEN
  IPWRITE(UNIT_stdOut,*)' FOUND', WeirdElems, 'ELEMENTS!'
  IPWRITE(UNIT_stdOut,*)' WEIRD ELEM NUMBERS:'
  DO iElem = 1,WeirdElems
    IPWRITE(UNIT_stdOut,*) WeirdElemNbrs(iElem)
  END DO
END IF

DEALLOCATE(WeirdElemNbrs)

SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE WeirdElementCheck


SUBROUTINE MapRegionToElem()
!----------------------------------------------------------------------------------------------------------------------------------!
! map a particle region to element
! check only element barycenter, nothing else
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Mesh_Vars ,ONLY: NbrOfRegions,RegionBounds,GEO
USE MOD_Particle_Mesh_Vars ,ONLY: ElemBaryNGeo
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
 INTEGER                :: iElem, iRegions
!===================================================================================================================================
SDEALLOCATE(GEO%ElemToRegion)
ALLOCATE(GEO%ElemToRegion(1:nElems))
GEO%ElemToRegion=0

DO iElem=1,nElems
  DO iRegions=1,NbrOfRegions
    IF ((ElemBaryNGeo(1,iElem).LT.RegionBounds(1,iRegions)).OR.(ElemBaryNGEO(1,iElem).GE.RegionBounds(2,iRegions))) CYCLE
    IF ((ElemBaryNGeo(2,iElem).LT.RegionBounds(3,iRegions)).OR.(ElemBaryNGEO(2,iElem).GE.RegionBounds(4,iRegions))) CYCLE
    IF ((ElemBaryNGeo(3,iElem).LT.RegionBounds(5,iRegions)).OR.(ElemBaryNGEO(3,iElem).GE.RegionBounds(6,iRegions))) CYCLE
    IF (GEO%ElemToRegion(iElem).EQ.0) THEN
      GEO%ElemToRegion(iElem)=iRegions
    ELSE
      CALL abort(&
__STAMP__&
,'Defined regions are overlapping')
    END IF
  END DO ! iRegions=1,NbrOfRegions
END DO ! iElem=1,nElems


END SUBROUTINE MapRegionToElem


SUBROUTINE BuildElementBasisAndRadius()
!================================================================================================================================
! Build the element local basis system, where the origin is located at xi=(0,0,0)^T and each local coordinate system is pointing
! to an element side
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis                  ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Basis         ,ONLY: DeCasteljauInterpolation
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Vars     ,ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo,XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadiusNGEO_Shared,ElemRadiusNGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo_Shared,ElemRadius2NGeo_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
#else
USE MOD_Mesh_Vars              ,ONLY: ElemBaryNGeo,XCL_NGeo,nELems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                        :: ALLOCSTAT
INTEGER                        :: iElem,SideID,i,j,k,ilocSide
REAL                           :: Xi(3,6),xPos(3),Radius
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem, lastElem, iDir
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
!================================================================================================================================
#if USE_MPI
  MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemRadiusNGeo_Shared_Win,ElemRadiusNGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadiusNGeo_Shared_Win,IERROR)
  CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemRadius2NGeo_Shared_Win,ElemRadius2NGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadius2NGeo_Shared_Win,IERROR)
  MPISharedSize = INT((3*6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,6,nComputeNodeTotalElems/),XiEtaZetaBasis_Shared_Win,XiEtaZetaBasis_Shared)
  CALL MPI_WIN_LOCK_ALL(0,XiEtaZetaBasis_Shared_Win,IERROR)
  MPISharedSize = INT((6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),slenXiEtaZetaBasis_Shared_Win,slenXiEtaZetaBasis_Shared)
  CALL MPI_WIN_LOCK_ALL(0,slenXiEtaZetaBasis_Shared_Win,IERROR)
  ElemRadiusNGeo     => ElemRadiusNGeo_Shared
  ElemRadius2NGeo    => ElemRadius2NGeo_Shared
  XiEtaZetaBasis     => XiEtaZetaBasis_Shared
  slenXiEtaZetaBasis => slenXiEtaZetaBasis_Shared

  ASSOCIATE(ElemBaryNGeo => ElemBaryNGeo_Shared, &
            XCL_NGeo     => XCL_NGeo_Shared)
#else
  ALLOCATE(ElemRadiusNGeo(          nElems) &
          ,ElemRadius2NGeo(         nElems) &
          ,XiEtaZetaBasis(1:3,1:6,1:nElems) &
          ,slenXiEtaZetaBasis(1:6,1:nElems))
#endif

ElemRadiusNGeo =0.
ElemRadius2NGeo=0.

#if USE_MPI
firstElem=INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

DO iElem=firstElem,lastElem
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
    xPos=0.
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          xPos=xPos+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
        END DO !i=0,NGeo
      END DO !j=0,NGeo
    END DO !k=0,NGeo
    XiEtaZetaBasis(1:3,iDir,iElem)=xPos
    ! compute vector from each barycenter to sidecenter
    XiEtaZetaBasis(:,iDir,iElem)=XiEtaZetaBasis(:,iDir,iElem)-ElemBaryNGeo(:,iElem)
    ! compute length: The root is omitted here due to optimization
    slenXiEtaZetaBasis(iDir,iElem)=1.0/DOT_PRODUCT(XiEtaZetaBasis(:,iDir,iElem),XiEtaZetaBasis(:,iDir,iElem))
  END DO ! iDir = 1, 6

  Radius=0.
  DO ilocSide=1,6
    SideID = GetGlobalNonUniqueSideID(GetGlobalElemID(iElem),iLocSide)
    IF(SideID.EQ.-1) CYCLE
    DO j=0,NGeo
      DO i=0,NGeo
        xPos=BezierControlPoints3D(:,i,j,SideID)-ElemBaryNGeo(:,iElem)
        Radius=MAX(Radius,SQRT(DOT_PRODUCT(xPos,xPos)))
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO ! ilocSide
  ElemRadiusNGeo (iElem)=Radius
  ElemRadius2NGeo(iElem)=Radius*Radius
END DO ! iElem

#if USE_MPI
END ASSOCIATE
CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(XiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(slenXiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

END SUBROUTINE BuildElementBasisAndRadius


SUBROUTINE BuildElementRadiusTria()
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Globals       ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo,ElemBaryNGeo_Shared,ElemRadius2NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared_Win,ElemRadius2NGeo_Shared_Win
#else
USE MOD_Mesh_Vars              ,ONLY: nELems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                        :: ALLOCSTAT
INTEGER                        :: iElem,iNode
REAL                           :: xPos(3),Radius
INTEGER                        :: firstElem, lastElem
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
!================================================================================================================================
#if USE_MPI
  MPISharedSize = INT((3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemBaryNGeo_Shared_Win,ElemBaryNGeo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemBaryNGeo_Shared_Win,IERROR)
  MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemRadius2NGeo_Shared_Win,ElemRadius2NGEO_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemRadius2NGeo_Shared_Win,IERROR)
  ElemRadius2NGeo    => ElemRadius2NGeo_Shared
  ElemBaryNGeo       => ElemBaryNGeo_Shared
#else
  ALLOCATE(ElemBaryNGeo(1:3,nElems) &
          ,ElemRadius2NGeo( nElems)
#endif

ElemRadius2NGeo=0.

#if USE_MPI
firstElem=INT(REAL(myComputeNodeRank*    nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem =INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem=1
lastElem=nElems
#endif

DO iElem=firstElem,lastElem
  Radius=0.
  xPos  =0.
  DO iNode=1,8
    xPos = xPos + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+iNode)
  END DO
    ElemBaryNGeo(:,iElem) = xPos/8.
  DO iNode=1,8
    xPos   = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+iNode) - ElemBaryNGeo(:,iElem)
    Radius = MAX(Radius,VECNORM(xPos))
  END DO
  ElemRadius2NGeo(iElem) = Radius*Radius
END DO ! iElem

#if USE_MPI
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

END SUBROUTINE BuildElementRadiusTria


SUBROUTINE PointsEqual(N,Points1,Points2,IsNotEqual)
!===================================================================================================================================
! compute the distance between two data sets
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: N
REAL,INTENT(IN)           :: Points1(1:3,1:N)
REAL,INTENT(IN)           :: Points2(1:3,1:N)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL                   :: IsNotEqual
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i
!===================================================================================================================================

IsNotEqual=.FALSE.

DO i=1,N
  IF( ABS(Points1(1,i)-Points2(1,i)).GT.1e-14 .OR. &
      ABS(Points1(2,i)-Points2(2,i)).GT.1e-14 .OR. &
      ABS(Points1(3,i)-Points2(3,i)).GT.1e-14 ) THEN
    IsNotEqual=.TRUE.
    RETURN
  END IF
END DO ! i=0,N

END SUBROUTINE PointsEqual


SUBROUTINE BuildElementOriginShared()
!================================================================================================================================
! compute the element origin at xi=(0,0,0)^T and set it as ElemBaryNGeo
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars          ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars ,ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars,ONLY: XCL_NGeo_Shared,ElemBaryNGeo_Shared
#if USE_MPI
USE MOD_Particle_MPI_Shared,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars,ONLY: ElemBaryNGeo_Shared_Win
#else
USE MOD_Particle_Mesh_Vars ,ONLY: XCL_NGeo
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,i,j,k
REAL                           :: XPos(3),buf
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem,lastElem
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
!================================================================================================================================
#if USE_MPI
MPISharedSize = INT((3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemBaryNGeo_Shared_Win,ElemBaryNGeo_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemBaryNGeo_Shared_Win,IERROR)

ASSOCIATE(XCL_NGeo => XCL_NGeo_Shared)

! Set ranges
firstElem = INT(REAL(myComputeNodeRank*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

ElemBaryNGeo_Shared=0.
DO iElem=firstElem,lastElem
  ! evaluate the polynomial at origin: Xi=(/0.0,0.0,0.0/)
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
  xPos=0.
  DO k=0,NGeo
    DO j=0,NGeo
      buf=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        xPos=xPos+XCL_NGeo(1:3,i,j,k,iElem)*Lag(1,i)*buf
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo
  ElemBaryNGeo_Shared(:,iElem)=xPos
END DO ! iElem

#if USE_MPI
END ASSOCIATE
CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementOriginShared


SUBROUTINE IdentifyElemAndSideType()
!===================================================================================================================================
!> get the element and side type of each element
!> 1) Get Elem Type (curved_elem)
!> 2) Get Side Type (planar_rect, planar_nonrect, bilineard, curved, planar_curved)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_ChangeBasis            ,ONLY: changeBasis3D
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Globals       ,ONLY: ALMOSTZERO,CROSSNORM,UNITVECTOR
USE MOD_Particle_Mesh_Vars     ,ONLY: Vdm_CLNGeo1_CLNGeo,Vdm_CLNGeo1_CLNGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: XCL_NGeo_Shared,ElemBaryNGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared,ElemCurved
USE MOD_Particle_Mesh_Vars     ,ONLY: BoundingBoxIsEmpty_Shared
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalElemID
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,SideType,SideNormVec,SideDistance
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideDistance_Shared,SideDistance_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideType_Shared,SideType_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: SideNormVec_Shared,SideNormVec_Shared_Win
#endif /* USE_MPI */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: ilocSide,SideID,flip
INTEGER                                  :: iElem,firstElem,lastElem
REAL,DIMENSION(1:3)                      :: v1,v2,v3
LOGICAL,ALLOCATABLE                      :: SideIsDone(:)
REAL                                     :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                                     :: XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                                     :: XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                                     :: BezierControlPoints_loc(1:3,0:NGeo,0:NGeo)
INTEGER                                  :: NGeo3,NGeo2
REAL                                     :: XCL_NGeoSideNew(1:3,0:NGeo,0:NGeo)
REAL                                     :: XCL_NGeoSideOld(1:3,0:NGeo,0:NGeo)
LOGICAL                                  :: isCurvedSide,isRectangular
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND)           :: MPISharedSize
#endif /* USE_MPI */
! output and sanity check
INTEGER                                  :: nPlanarRectangular,   nPlanarNonRectangular,   nPlanarCurved,   nBilinear,   nCurved
INTEGER                                  :: nPlanarRectangularTot,nPlanarNonRectangularTot,nPlanarCurvedTot,nBilinearTot,nCurvedTot
INTEGER                                  :: nLinearElems,   nCurvedElems
INTEGER                                  :: nLinearElemsTot,nCurvedElemsTot
#if USE_MPI
INTEGER                                  :: nDummy
#endif /* USE_MPI */
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying side types and whether elements are curved ...'

! elements
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemCurved_Shared_Win,ElemCurved_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemCurved_Shared_Win,IERROR)
ElemCurved => ElemCurved_Shared
#else
ALLOCATE(ElemCurved(1:nComputeNodeTotalElems))
#endif /*USE_MPI*/
!IF (.NOT.DoRefMapping) THEN
!  ALLOCATE(ElemType(1:nTotalElems))
!  ElemType=-1
!END IF

! sides
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalSides/),SideType_Shared_Win,SideType_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideType_Shared_Win,IERROR)
SideType => SideType_Shared
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalSides/),SideDistance_Shared_Win,SideDistance_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideDistance_Shared_Win,IERROR)
SideDistance => SideDistance_Shared
MPISharedSize = INT((3*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalSides/),SideNormVec_Shared_Win,SideNormVec_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideNormVec_Shared_Win,IERROR)
SideNormVec => SideNormVec_Shared
#else
ALLOCATE(SideType(       nComputeNodeTotalSides))
ALLOCATE(SideDistance(   nComputeNodeTotalSides))
ALLOCATE(SideNormVec(1:3,nComputeNodeTotalSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemCurved   = .FALSE.
  SideType     = -1
  SideDistance = -0.
  SideNormVec  = 0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemCurved_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideType_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideDistance_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideNormVec_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

ALLOCATE(SideIsDone(nComputeNodeTotalSides))
SideIsDone = .FALSE.

NGeo2 = (NGeo+1)*(NGeo+1)
NGeo3 = NGeo2   *(NGeo+1)

! decide if element is (bi-)linear or curved
! decide if sides are planar-rect, planar-nonrect, planar-curved, bilinear or curved
#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

DO iElem=firstElem,lastElem
  XCL_NGeoLoc = XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
  ! 1) check if elem is curved
  !   a) get the coordinates of the eight nodes of the hexahedral
  XCL_NGeo1(1:3,0,0,0) = XCL_NGeoLoc(1:3, 0  , 0  , 0  )
  XCL_NGeo1(1:3,1,0,0) = XCL_NGeoLoc(1:3,NGeo, 0  , 0  )
  XCL_NGeo1(1:3,0,1,0) = XCL_NGeoLoc(1:3, 0  ,NGeo, 0  )
  XCL_NGeo1(1:3,1,1,0) = XCL_NGeoLoc(1:3,NGeo,NGeo, 0  )
  XCL_NGeo1(1:3,0,0,1) = XCL_NGeoLoc(1:3, 0  , 0  ,NGeo)
  XCL_NGeo1(1:3,1,0,1) = XCL_NGeoLoc(1:3,NGeo, 0  ,NGeo)
  XCL_NGeo1(1:3,0,1,1) = XCL_NGeoLoc(1:3, 0  ,NGeo,NGeo)
  XCL_NGeo1(1:3,1,1,1) = XCL_NGeoLoc(1:3,NGeo,NGeo,NGeo)

  !  b) interpolate from the nodes to NGeo
  !     Compare the bi-linear mapping with the used mapping
  !     For NGeo=1, this should always be true, because the mappings are identical
  CALL ChangeBasis3D(3,1,NGeo,Vdm_CLNGeo1_CLNGeo,XCL_NGeo1,XCL_NGeoNew)
  ! check the coordinates of all Chebychev-Lobatto geometry points between the bi-linear and used
  ! mapping
  CALL PointsEqual(NGeo3,XCL_NGeoNew,XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0:NGeo),ElemCurved(iElem))

  ! 2) check sides
  ! loop over all 6 sides of element
  ! a) check if the sides are straight
  ! b) use curved information to decide side type
  DO ilocSide=1,6
    SideID = GetGlobalNonUniqueSideID(GetGlobalElemID(iElem),iLocSide)

    ! Why were only flips LT. 0 considered? All flips are equal!
!    IF (SideInfo_Shared(SIDE_ID,SideID).GT.0) THEN
!      flip = 0
!    ELSE
!      flip = MOD(Sideinfo_Shared(SIDE_FLIP,SideID),10)
!    END IF
    flip = Sideinfo_Shared(SIDE_FLIP,SideID)

    IF(.NOT.ElemCurved(iElem))THEN
      BezierControlPoints_loc(1:3,0:NGeo,0:NGeo) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      ! linear element
      IF(BoundingBoxIsEmpty_Shared(SideID))THEN
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,SideID) = CROSSNORM(v1,v2)
        v1=0.25*(BezierControlPoints_loc(:,0,0      )     &
                +BezierControlPoints_loc(:,NGeo,0   )  &
                +BezierControlPoints_loc(:,0,NGeo   )  &
                +BezierControlPoints_loc(:,NGeo,NGeo))
        ! check if normal vector points outwards
        v2=v1-ElemBaryNGeo_Shared(:,iElem)
        IF(flip.EQ.0)THEN
          IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).LT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
        ELSE
          IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).GT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
        END IF
        SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
        ! check if it is rectangular
        isRectangular=.TRUE.
        v1=UNITVECTOR(BezierControlPoints_loc(:,0   ,NGeo)-BezierControlPoints_loc(:,0   ,0   ))
        v2=UNITVECTOR(BezierControlPoints_loc(:,NGeo,0   )-BezierControlPoints_loc(:,0   ,0   ))
        v3=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,0   ,NGeo))
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
        IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        IF(isRectangular)THEN
          v1=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,NGeo,0))
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
        END IF
        IF(isRectangular)THEN
          SideType(SideID)=PLANAR_RECT
        ELSE
          SideType(SideID)=PLANAR_NONRECT
        END IF
      ELSE
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,SideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
        SideType(SideID)=BILINEAR
      END IF
    ELSE
      ! possible curved face
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0,0:NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0,0:NGeo,0:NGeo)
      CASE(XI_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,NGeo,0:NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,NGeo,0:NGeo,0:NGeo)
      CASE(ETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0,0:NGeo)
      CASE(ETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,NGeo,0:NGeo)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,NGeo,0:NGeo)
      CASE(ZETA_MINUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,0)
        XCL_NGeoSideNew=XCL_NGeoNew(1:3,0:NGeo,0:NGeo,0)
      CASE(ZETA_PLUS)
        XCL_NGeoSideOld=XCL_NGeoLoc(1:3,0:NGeo,0:NGeo,NGeo)
        XCL_NGeoSideNew=XCL_NGeoNEw(1:3,0:NGeo,0:NGeo,NGeo)
      END SELECT
      CALL PointsEqual(NGeo2,XCL_NGeoSideNew,XCL_NGeoSideOld,isCurvedSide)
      IF(isCurvedSide)THEn
        IF(BoundingBoxIsEmpty_Shared(SideID))THEN
          SideType(SideID)=PLANAR_CURVED
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,SideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0   )  &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo_Shared(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).LT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).GT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          END IF
          SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
        ELSE
          SideType(SideID)=CURVED
        END IF
      ELSE
        IF (BoundingBoxIsEmpty_Shared(SideID)) THEN
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,SideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0)     &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo_Shared(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).LT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,SideID)).GT.0) SideNormVec(:,SideID)=-SideNormVec(:,SideID)
          END IF
          SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
          ! check if it is rectangular
          isRectangular=.TRUE.
          v1=UNITVECTOR(BezierControlPoints_loc(:,0   ,NGeo)-BezierControlPoints_loc(:,0   ,0   ))
          v2=UNITVECTOR(BezierControlPoints_loc(:,NGeo,0   )-BezierControlPoints_loc(:,0   ,0   ))
          v3=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,0   ,NGeo))
!          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
!          IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
          IF(DOT_PRODUCT(v1,v2).GT.1E-14) isRectangular=.FALSE.
          IF(DOT_PRODUCT(v1,v3).GT.1E-14) isRectangular=.FALSE.
          IF(isRectangular)THEN
            v1=UNITVECTOR(BezierControlPoints_loc(:,NGeo,NGeo)-BezierControlPoints_loc(:,NGeo,0))
!            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v2))) isRectangular=.FALSE.
!            IF(.NOT.ALMOSTZERO(DOT_PRODUCT(v1,v3))) isRectangular=.FALSE.
            IF(DOT_PRODUCT(v1,v2).GT.1E-14) isRectangular=.FALSE.
            IF(DOT_PRODUCT(v1,v3).GT.1E-14) isRectangular=.FALSE.
          END IF
          IF(isRectangular)THEN
            SideType(SideID)=PLANAR_RECT
          ELSE
            SideType(SideID)=PLANAR_NONRECT
          END IF
        ELSE
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          SideNormVec(:,SideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
          SideType(SideID)=BILINEAR
        END IF
      END IF
    END IF
    SideIsDone(SideID)=.TRUE.
  END DO ! ilocSide=1,6
END DO ! iElem=1,nTotalElems

! TODO with PERIODIC SIDES
! sanity check for side periodic type
!DO iSide=1,nPartSides
!  IF(DoRefmapping)THEN
!    BCSideID  =PartBCSideList(iSide)
!    IF(BCSideID.LE.0) CYCLE
!  ELSE
!    BCSideID  =iSide
!  END IF
!  PVID=SidePeriodicType(iSide)
!  IF(PVID.EQ.0) CYCLE
!  IF(.NOT.PartMeshHasPeriodicBCs) CYCLE
!  Vec1=SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
!  ScalarProduct=DOT_PRODUCT(SideNormVec(1:3,BCSideID),Vec1)
!  IF(ALMOSTEQUAL(ScalarProduct,GEO%PeriodicVectorsLength(ABS(PVID))))THEN
!    SidePeriodicType(iSide)=-SidePeriodicType(iSide)
!  ELSEIF(.NOT.ALMOSTEQUAL(ScalarProduct,-GEO%PeriodicVectorsLength(ABS(PVID))))THEN
!    WRITE (*,*) "BCSideID                  : ", BCSideID
!    WRITE (*,*) "SideNormVec(1:3,BCSideID) : ", SideNormVec(1:3,BCSideID)
!    WRITE (*,*) "Vec1                      : ", Vec1
!    CALL abort(&
!__STAMP__&
!        , ' Missalignment between SideNormVec and PeriodicVector!',ABS(PVID),ScalarProduct)
!  END IF
!END DO ! iSide=1,nPartSides

!! fill Element type checking sides
!IF (.NOT.DoRefMapping) THEN
!  DO iElem=1,nTotalElems
!    DO ilocSide=1,6
!      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!      SELECT CASE(SideType(SideID))
!      CASE(PLANAR_RECT,PLANAR_NONRECT)
!        IF (ElemType(iElem).GE.1) THEN
!          CYCLE
!        ELSE
!          ElemType(iElem) = 1
!        END IF
!      CASE(BILINEAR)
!        IF (ElemType(iElem).GE.2) THEN
!          CYCLE
!        ELSE
!          ElemType(iElem) = 2
!        END IF
!      CASE(PLANAR_CURVED,CURVED)
!        ElemType(iElem) = 3
!        EXIT
!      END SELECT
!    END DO ! ilocSide=1,6
!  END DO ! iElem=1,nTotalElems
!END IF

#if USE_MPI
CALL MPI_WIN_SYNC(ElemCurved_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI */

! zero counter for side types
nPlanarRectangular         = 0
nPlanarNonRectangular      = 0
nPlanarCurved              = 0
nBilinear                  = 0
nCurved                    = 0
! zero counter for elem types
nCurvedElems               = 0
nLinearElems               = 0

DO iElem = firstElem,lastElem
  IF (ElemCurved(iElem)) THEN
    nCurvedElems = nCurvedElems+1
  ELSE
    nLinearElems = nLinearElems+1
  END IF

  DO ilocSide = 1,6
    ! ignore small mortar sides attached to big mortar sides
    SideID = GetGlobalNonUniqueSideID(iElem,ilocSide)
    SELECT CASE(SideType(SideID))
      CASE (PLANAR_RECT)
        nPlanarRectangular    = nPlanarRectangular   +1
      CASE (PLANAR_NONRECT)
        nPlanarNonRectangular = nPlanarNonRectangular+1
      CASE (BILINEAR)
        nBilinear             = nBilinear            +1
      CASE (PLANAR_CURVED)
        nPlanarCurved         = nPlanarCurved        +1
      CASE (CURVED)
        nCurved               = nCurved              +1
    END SELECT
  END DO
END DO

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(nPlanarRectangular   ,nPlanarRectangularTot   ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nPlanarNonRectangular,nPlanarNonRectangularTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nBilinear            ,nBilinearTot            ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nPlanarCurved        ,nPlanarCurvedTot        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nCurved              ,nCurvedTot              ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nLinearElems         ,nLinearElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nCurvedElems         ,nCurvedElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(nPlanarRectangular   ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nPlanarNonRectangular,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nBilinear            ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nPlanarCurved        ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nCurved              ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nLinearElems         ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
  CALL MPI_REDUCE(nCurvedElems         ,nDummy                  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
END IF
#else
nPlanarRectangularTot    = nPlanarRectangular
nPlanarNonRectangularTot = nPlanarNonRectangular
nBilinearTot             = nBilinear
nPlanarCurvedTot         = nPlanarCurved
nCurvedTot               = nCurved
nLinearElemsTot          = nLinearElems
nCurvedElemsTot          = nCurvedElems
#endif /* USE_MPI */

! sanity check
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI */
  IF ((nComputeNodeTotalElems-nCurvedElemsTot).NE.nLinearElemsTot) &
    CALL ABORT(__STAMP__, 'Error in particle mesh: lost elements while trying to dermine if elements are curved')
#if USE_MPI
END IF
#endif /* USE_MPI */


SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-rectangular     faces: ', nPlanarRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-non-rectangular faces: ', nPlanarNonRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of bi-linear              faces: ', nBilineartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-curved          faces: ', nPlanarCurvedtot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 faces: ', nCurvedtot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of (bi-)linear            elems: ', nLinearElemsTot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 elems: ', nCurvedElemsTot

END SUBROUTINE IdentifyElemAndSideType


SUBROUTINE BuildBCElemDistance()
!===================================================================================================================================
! get the distance of each BC face
!> 1) identify BC elements to be handled by Tracing
!> 2) build mapping global elem ID to BC elem ID
!> 3) calculate distance from element origin to each side
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Basis         ,ONLY: DeCasteljauInterpolation
USE MOD_Particle_Globals       ,ONLY: ALMOSTZERO
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides,SideBCMetrics
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared,ElemRadiusNGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides_Shared,SideBCMetrics_Shared
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Utils         ,ONLY: InsertionSort
USE MOD_Particle_Vars          ,ONLY: ManualTimeStep
#if USE_MPI
USE MOD_CalcTimeStep           ,ONLY: CalcTimeStep
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBCSides_Shared_Win,SideBCMetrics_Shared_Win
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps,halo_eps_velo
USE MOD_TimeDisc_Vars          ,ONLY: nRKStages,RKc
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: p,q,r,s
INTEGER                        :: firstBezierPoint,lastBezierPoint
INTEGER                        :: ElemID,SideID
INTEGER                        :: iElem,firstElem,lastElem
INTEGER                        :: iSide,firstSide,lastSide,iLocSide
INTEGER                        :: nComputeNodeBCSides
INTEGER                        :: nBCSidesElem,nBCSidesProc,offsetBCSidesProc,offsetBCSides
INTEGER,ALLOCATABLE            :: ElemToBCSidesProc(:,:)
LOGICAL,ALLOCATABLE            :: BCSideElem(:)
LOGICAL                        :: SideInRange
REAL                           :: dX,dY,dZ
REAL                           :: origin(1:3),xi(1:2),radius,radiusMax,vec(1:3)
REAL                           :: NodesLocSide(1:3,0:NGeo,0:NGeo),NodeBCSide(1:3)
REAL                           :: BC_halo_eps,BC_halo_eps_velo,deltaT
REAL,ALLOCATABLE               :: tmpSideBCMetrics(:,:)
REAL,ALLOCATABLE               :: tmpSideBCDistance(:)
INTEGER,ALLOCATABLE            :: intSideBCMetrics(:)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: errType
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
INTEGER                        :: iStage
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying BC sides and calculating side metrics ...'

! elements
#if USE_MPI
MPISharedSize = INT((2*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,nComputeNodeTotalElems/),ElemToBCSides_Shared_Win,ElemToBCSides_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToBCSides_Shared_Win,IERROR)
ElemToBCSides => ElemToBCSides_Shared
#else
ALLOCATE(ElemToBCSides(1:2,1:nComputeNodeTotalElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemToBCSides = -1
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemToBCSides_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
firstSide = 1
lastSide  = nComputeNodeTotalSides
ALLOCATE(ElemToBCSidesProc(1:2,1:nComputeNodeTotalElems))
ALLOCATE(BCSideElem(1:nComputeNodeTotalSides))

! if running on one node, halo_eps is meaningless. Get a representative BC_halo_eps for BC side identification
IF (halo_eps.EQ.0) THEN
  ! reconstruct halo_eps_velo
  IF (halo_eps_velo.EQ.0) THEN
    ! Set to speed of sound
    BC_halo_eps_velo = 343
  ELSE
    BC_halo_eps_velo = halo_eps_velo
  END IF

  ! reconstruct deltaT
  deltaT = 0.
  IF (ManualTimeStep.GT.0.) THEN
    deltaT    = ManualTimeStep
  ELSE
    deltaT    = CalcTimeStep(errType)
  END IF

  ! calculate halo_eps
  BC_halo_eps = RKc(2)
  DO iStage=2,nRKStages-1
    BC_halo_eps = MAX(BC_halo_eps,RKc(iStage+1)-RKc(iStage))
  END DO
  BC_halo_eps = MAX(BC_halo_eps,1.-RKc(nRKStages))
  BC_halo_eps = BC_halo_eps*BC_halo_eps_velo*deltaT

  ! sanity check the calculated halo_eps
  IF (.NOT.ALMOSTZERO(BC_halo_eps)) THEN
    SWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to ',BC_halo_eps
  ! halo_eps still at zero. Set it to global diagonal
  ELSE
    vec(1)   = GEO%xmaxglob-GEO%xminglob
    vec(2)   = GEO%ymaxglob-GEO%yminglob
    vec(3)   = GEO%zmaxglob-GEO%zminglob
    BC_halo_eps = DOT_PRODUCT(vec,vec)
    SWRITE(UNIT_stdOUt,'(A,F11.3)') ' | No halo_eps given and could not be reconstructed. Using global diag with ',BC_halo_eps
  END IF
END IF

#else
! get distance of diagonal of mesh
vec(1)   = GEO%xmaxglob-GEO%xminglob
vec(2)   = GEO%ymaxglob-GEO%yminglob
vec(3)   = GEO%zmaxglob-GEO%zminglob
BC_halo_eps = DOT_PRODUCT(vec,vec)

firstElem = 1
lastElem  = nElems
firstSide = 1
lastSide  = nTotalSides
ALLOCATE(ElemToBCSidesProc(1:2,1:nElems))
ALLOCATE(BCSideElem(1:nTotalSides))
#endif /*USE_MPI*/

ElemToBCSidesProc = -1
nBCSidesProc      = 0
offsetBCSides     = 0

! sum up all BC sides in range of BC_halo_eps
DO iElem = firstElem,lastElem
  nBCSidesElem  = 0
  BCSideElem    = .FALSE.

  ! check local side of an element
  DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)
    IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
    nBCSidesElem = nBCSidesElem + 1
    nBCSidesProc = nBCSidesProc + 1
  END DO

  ! loop over all sides. Check distance from every local side to total sides. Once a side has been flagged,
  ! it must not be counted again
  DO ilocSide = 1,6
    SideID = GetGlobalNonUniqueSideID(iElem,ilocSide)
    ! Get BezierControlPoints for local side. BezierControlPoints3D available for ALL sides in shared memory
    NodesLocSide(:,:,:)= BezierControlPoints3D(:,:,:,SideID)
    SELECT CASE(ilocSide)
      CASE(XI_MINUS,XI_PLUS)
        firstBezierPoint = 0
        lastBezierPoint  = NGeo
      CASE DEFAULT
        firstBezierPoint = 1
        lastBezierPoint  = NGeo-1
    END SELECT

    ! check every BC side
    DO iSide = firstSide,lastSide
      ! ignore non-BC sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
      ! ignore sides of the same element
      IF (SideInfo_Shared(SIDE_ELEMID,iSide).EQ.iElem) CYCLE
      ! ignore already flagged sides
      IF (BCSideElem(iSide).EQV..TRUE.) CYCLE

      ! compare all nodes between local and BC side to check if within BC_halo_eps. Once one node pair is in range, flag the entire
      ! side and stop checking
      SideInRange = .FALSE.
      DO q = firstBezierPoint,lastBezierPoint
        DO p = firstBezierPoint,lastBezierPoint
          ! get all nodes for BC side
          NodeBCSide(:) = BezierControlPoints3D(:,p,q,iSide)
          ! finally compare the node coords
          DO s = firstBezierPoint,lastBezierPoint
            DO r = firstBezierPoint,lastBezierPoint
              dX = ABS(NodesLocSide(1,r,s)-NodeBCSide(1))
              IF (dX.GT.BC_halo_eps) CYCLE
              dY = ABS(NodesLocSide(2,r,s)-NodeBCSide(2))
              IF (dY.GT.BC_halo_eps) CYCLE
              dZ = ABS(NodesLocSide(3,r,s)-NodeBCSide(3))
              IF (dZ.GT.BC_halo_eps) CYCLE

              IF (SQRT(dX*dX+dY*dY+dZ*dZ).LE.BC_halo_eps) THEN
                nBCSidesElem = nBCSidesElem + 1
                nBCSidesProc = nBCSidesProc + 1
                BCSideElem(iSide) = .TRUE.
                SideInRange       = .TRUE.
                EXIT
              END IF
            END DO ! r
            IF (SideInRange) EXIT
          END DO ! s
          IF (SideInRange) EXIT
        END DO ! p
        IF (SideInRange) EXIT
      END DO ! q
    END DO ! iSide
  END DO ! ilocSide

  ! Write local mapping from Elem to BC sides. The number is already correct, the offset must be corrected later
  IF (nBCSidesElem.GT.0) THEN
    ElemToBCSidesProc(ELEM_NBR_BCSIDES ,iElem) = nBCSidesElem
    ElemToBCSidesProc(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides + 1
  END IF

  offsetBCSides = nBCSidesProc
END DO

! Find CN global number of BC sides and write Elem to BC Side mapping into shared array
#if USE_MPI
sendbuf = nBCSidesProc
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetBCSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetBCSidesProc + nBCSidesProc
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeBCSides = sendbuf

ElemToBCSides(1,firstElem:lastElem) = ElemToBCSidesProc(1,firstElem:lastElem)
ElemToBCSides(2,firstElem:lastElem) = ElemToBCSidesProc(2,firstElem:lastElem) + offsetBCSidesProc
#else
offsetBCSidesProc   = 0
nComputeNodeBCSides = nBCSidesProc

ElemToBCSides(:,firstElem:lastElem) = ElemToBCSidesProc(:,firstElem:lastElem)
#endif /*USE_MPI*/
DEALLOCATE(ElemToBCSidesProc)

! Allocate shared array for BC sides
#if USE_MPI
MPISharedSize = INT((7*nComputeNodeBCSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/7,nComputeNodeBCSides/),SideBCMetrics_Shared_Win,SideBCMetrics_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideBCMetrics_Shared_Win,IERROR)
SideBCMetrics => SideBCMetrics_Shared
#else
ALLOCATE(SideBCMetrics(1:7,1:nComputeNodeBCSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  SideBCMetrics = -1.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

nBCSidesProc      = 0

! We did not know the number of BC sides before. Therefore, we need to do the check again and build the final mapping
DO iElem = firstElem,lastElem
  BCSideElem    = .FALSE.

  ! check local side of an element
  DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)
    IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
    nBCSidesProc = nBCSidesProc + 1
    SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
    SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(iElem)
  END DO

  ! loop over all sides. Check distance from every local side to total sides. Once a side has been flagged,
  ! it must not be counted again
  DO ilocSide = 1,6
    SideID = GetGlobalNonUniqueSideID(iElem,ilocSide)
      ! Get BezierControlPoints for local side. BezierControlPoints3D available for ALL sides in shared memory
      NodesLocSide(:,:,:)= BezierControlPoints3D(:,:,:,SideID)
      SELECT CASE(ilocSide)
        CASE(XI_MINUS,XI_PLUS)
          firstBezierPoint = 0
          lastBezierPoint  = NGeo
        CASE DEFAULT
          firstBezierPoint = 1
          lastBezierPoint  = NGeo-1
      END SELECT

      ! check every BC side
      DO iSide = firstSide,lastSide
        ! ignore non-BC sides
        IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
        ! ignore sides of the same element
        IF (SideInfo_Shared(SIDE_ELEMID,iSide).EQ.iElem) CYCLE
        ! ignore already flagged sides
        IF (BCSideElem(iSide).EQV..TRUE.) CYCLE

        ! compare all nodes between local and BC side to check if within halo_eps. Once one node pair is in range, flag the entire
        ! side and stop checking
        SideInRange = .FALSE.
        DO q = firstBezierPoint,lastBezierPoint
          DO p = firstBezierPoint,lastBezierPoint
            ! get all nodes for BC side
            NodeBCSide(:) = BezierControlPoints3D(:,p,q,iSide)
            ! finally compare the node coords
            DO s = firstBezierPoint,lastBezierPoint
              DO r = firstBezierPoint,lastBezierPoint
                dX = ABS(NodesLocSide(1,r,s)-NodeBCSide(1))
                IF (dX.GT.BC_halo_eps) CYCLE
                dY = ABS(NodesLocSide(2,r,s)-NodeBCSide(2))
                IF (dY.GT.BC_halo_eps) CYCLE
                dZ = ABS(NodesLocSide(3,r,s)-NodeBCSide(3))
                IF (dZ.GT.BC_halo_eps) CYCLE

                IF (SQRT(dX*dX+dY*dY+dZ*dZ).LE.BC_halo_eps) THEN
                  nBCSidesProc = nBCSidesProc + 1
                  SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
                  SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(iElem)
                  BCSideElem(iSide) = .TRUE.
                  SideInRange       = .TRUE.
                  EXIT
                END IF
              END DO ! r
              IF (SideInRange) EXIT
            END DO ! s
            IF (SideInRange) EXIT
          END DO ! p
          IF (SideInRange) EXIT
        END DO ! q
      END DO ! iSide
    END DO ! ilocSide
END DO

#if USE_MPI
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeBCSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeBCSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nComputeNodeBCSides
#endif /*USE_MPI*/

! calculate origin, radius and distance to sides
DO iSide = firstSide,lastSide
  SideID = INT(SideBCMetrics(BCSIDE_SIDEID,iSide))
  ElemID = INT(SideBCMetrics(BCSIDE_ELEMID,iSide))

  !> build side origin
  xi     = 0.
  ! TODO: BezierControlPoints are alloced with global side ID, so this SHOULD work. Breaks if we reduce the halo region
  CALL DeCasteljauInterpolation(NGeo,xi,SideID,origin)
  SideBCMetrics(5:7,iSide) = origin(1:3)

  !> build side radius
  radiusMax = 0.
  DO q = 0,NGeo
    DO p = 0,NGeo
      vec(1:3) = BezierControlPoints3D(:,p,q,SideID) - origin
      radius   = DOT_PRODUCT(Vec,Vec)
      radiusMax= MAX(radiusMax,radius)
    END DO
  END DO
  SideBCMetrics(BCSIDE_RADIUS,iSide) = SQRT(RadiusMax)

  !> build side distance
  origin(1:3) = ElemBaryNGeo_Shared(1:3,ElemID)
  vec(1:3)    = origin(1:3) - SideBCMetrics(5:7,iSide)
  SideBCMetrics(BCSIDE_DISTANCE,iSide) = SQRT(DOT_PRODUCT(vec,vec))-ElemRadiusNGeo_Shared(ElemID)-SideBCMetrics(BCSIDE_RADIUS,iSide)
END DO ! iSide

#if USE_MPI
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

! finally, sort by distance to help speed up BC tracking
!> allocate dummy array to hold variables
ALLOCATE(tmpSideBCMetrics(1:7,1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))
ALLOCATE(tmpSideBCDistance(   1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))
ALLOCATE(intSideBCMetrics(    1:MAXVAL(ElemToBCSides(ELEM_NBR_BCSIDES,:))))

DO iElem = firstElem,lastElem
  ! skip elements with no BC sides
  IF (ElemToBCSides(ELEM_NBR_BCSIDES,iElem).LE.0) CYCLE

  ! save values in temporary array
  firstSide    = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem)
  lastSide     = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem)+ElemToBCSides(ELEM_NBR_BCSIDES,iElem)-1
  nBCSidesElem = ElemToBCSides(ELEM_NBR_BCSIDES,iElem)

  tmpSideBCMetrics(:,1:nBCSidesElem) = SideBCMetrics(:,firstSide:lastSide)
  tmpSideBCDistance( 1:nBCSidesElem) = SideBCMetrics(BCSIDE_DISTANCE,firstSide:lastSide)
  intSideBCMetrics(  1:nBCSidesElem) = (/((p),p=1,nBCSidesElem)/)

  ! sort SideID according to distance
  CALL InsertionSort(tmpSideBCDistance(1:nBCSidesElem),intSideBCMetrics(1:nBCSidesElem),nBCSidesElem)

  ! write back dummy array with variables
  DO iSide = 1,nBCSidesElem
    SideID = intSideBCMetrics(iSide)
    SideBCMetrics(:,firstSide+iSide-1) = tmpSideBCMetrics(:,SideID)
  END DO
END DO

DEALLOCATE(tmpSideBCMetrics)
DEALLOCATE(tmpSideBCDistance)
DEALLOCATE(intSideBCMetrics)

#if USE_MPI
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif

END SUBROUTINE BuildBCElemDistance


SUBROUTINE BuildEpsOneCell()
!===================================================================================================================================
! Build epsOneCell for each element
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars     ,ONLY: NodeTypeCL,NodeType
USE MOD_Mesh_Vars              ,ONLY: NGeo,NGeoRef
USE MOD_Mesh_Vars              ,ONLY: sJ
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: dXCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ_Shared,ElemEpsOneCell_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: RefMappingEps
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemsJ_Shared_Win,ElemEpsOneCell_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k
! Vandermonde matrices
REAL                           :: Vdm_CLNGeo_NGeoRef(0:NgeoRef,0:NGeo)
REAL                           :: Vdm_NGeoRef_N(     0:PP_N   ,0:NGeoRef)
! Jacobian on CL N and NGeoRef
REAL                           :: detJac_Ref(1  ,0:NGeoRef,0:NGeoRef,0:NGeoRef)
REAL                           :: DetJac_N(  1  ,0:PP_N,   0:PP_N,   0:PP_N)
! interpolation points and derivatives on CL N
REAL                           :: dX_NGeoRef(3,3,0:NGeoRef,0:NGeoRef,0:NGeoRef)

INTEGER                        :: iElem,firstElem,lastElem,ElemLocID
REAL                           :: scaleJ,maxScaleJ
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Building EpsOneCell for all elements ...'

! build sJ for all elements not on local proc
#if USE_MPI
MPISharedSize = INT(((PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/(PP_N+1)*(PP_N+1)*(PP_N+1)*nComputeNodeTotalElems/),ElemsJ_Shared_Win,ElemsJ_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemsJ_Shared_Win,IERROR)
ElemsJ(0:PP_N,0:PP_N,0:PP_N,1:nComputeNodeTotalElems) => ElemsJ_Shared

IF (myComputeNodeRank.EQ.0) THEN
  ElemsJ_Shared = 0.
END IF

CALL MPI_WIN_SYNC(ElemsJ_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeElems)/REAL(nComputeNodeProcessors))

! Calculate sJ for elements not inside current proc, otherwise copy local values
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , NgeoRef , NodeType  , Vdm_CLNGeo_NGeoRef, modal=.FALSE.)
CALL GetVandermonde(    NgeoRef, NodeType    , PP_N    , NodeType  , Vdm_NGeoRef_N     , modal=.TRUE.)

DO iElem = firstElem,lastElem
  ElemLocID = iElem-offsetElem
  ! element on local proc, sJ already calculated in metrics.f90
  IF ((ElemLocID.GT.0) .AND. (ElemLocID.LE.nElems)) THEN
    ElemsJ(:,:,:,iElem) = sJ(:,:,:,ElemLocID,0)

  ! element not on local proc, calculate sJ frm dXCL_NGeo_Shared
  ELSE
    detJac_Ref = 0.
    ! Compute Jacobian on NGeo and then interpolate:
    ! required to guarantee conservativity when restarting with N<NGeo
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,1,:,:,:,iElem),dX_NGeoRef(:,1,:,:,:))
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,2,:,:,:,iElem),dX_NGeoRef(:,2,:,:,:))
    CALL ChangeBasis3D(3,Ngeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo_Shared(:,3,:,:,:,iElem),dX_NGeoRef(:,3,:,:,:))
    DO k=0,NGeoRef; DO j=0,NGeoRef; DO i=0,NGeoRef
      detJac_Ref(1,i,j,k)=detJac_Ref(1,i,j,k) &
        + dX_NGeoRef(1,1,i,j,k)*(dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(3,3,i,j,k) - dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(2,3,i,j,k))  &
        + dX_NGeoRef(2,1,i,j,k)*(dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(1,3,i,j,k) - dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(3,3,i,j,k))  &
        + dX_NGeoRef(3,1,i,j,k)*(dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(2,3,i,j,k) - dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(1,3,i,j,k))
    END DO; END DO; END DO !i,j,k=0,NgeoRef

    ! interpolate detJac_ref to the solution points
    CALL ChangeBasis3D(1,NgeoRef,PP_N,Vdm_NgeoRef_N,DetJac_Ref(:,:,:,:),DetJac_N)

    ! assign to global Variable sJ
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ElemsJ(i,j,k,iElem)=1./DetJac_N(1,i,j,k)
    END DO; END DO; END DO !i,j,k=0,PP_N
  END IF
END DO

CALL MPI_WIN_SYNC(ElemsJ_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ALLOCATE(ElemsJ(0:PP_N,0:PP_N,0:PP_N,1:nElems))
ElemsJ => sJ
#endif /* USE_MPI*/

! allocate epsOneCell
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemEpsOneCell_Shared_Win,ElemEpsOneCell_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemEpsOneCell_Shared_Win,IERROR)
ElemEpsOneCell => ElemEpsOneCell_Shared
#else
ALLOCATE(ElemEpsOneCell(1:nComputeNodeTotalElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemEpsOneCell = -1.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemEpsOneCell_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

maxScaleJ = 0.
DO iElem = firstElem,lastElem
  scaleJ = MAXVAL(ElemsJ(:,:,:,iElem))/MINVAL(ElemsJ(:,:,:,iElem))
  ElemepsOneCell(iElem) = 1.0 + SQRT(3.0*scaleJ*RefMappingEps)
  maxScaleJ  =MAX(scaleJ,maxScaleJ)
END DO ! iElem = firstElem,lastElem

!IF(CalcMeshInfo)THEN
!  CALL AddToElemData(ElementOut,'epsOneCell',RealArray=epsOneCell(1:nElems))
!END IF

END SUBROUTINE BuildEpsOneCell


SUBROUTINE GetLinearSideBaseVectors()
!===================================================================================================================================
! computes the face base vector for linear (planar or bilinear) face intersection calculation
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mathtools              ,ONLY: CROSS
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: NGeoElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,BaseVectorsScale
#if USE_MPI
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectorsScale_Shared,BaseVectorsScale_Shared_Win
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectors0_Shared,BaseVectors1_Shared,BaseVectors2_Shared,BaseVectors3_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: BaseVectors0_Shared_Win,BaseVectors1_Shared_Win,BaseVectors2_Shared_Win,BaseVectors3_Shared_Win
#endif /* USE_MPI */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSide,firstSide,lastSide
REAL                           :: crossVec(3)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' GET LINEAR SIDE BASEVECTORS...'
#if USE_MPI
MPISharedSize = INT((3*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalSides/),BaseVectors0_Shared_Win,BaseVectors0_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors0_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalSides/),BaseVectors1_Shared_Win,BaseVectors1_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors1_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalSides/),BaseVectors2_Shared_Win,BaseVectors2_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors2_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalSides/),BaseVectors3_Shared_Win,BaseVectors3_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors3_Shared_Win,IERROR)
MPISharedSize = INT((nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalSides/),BaseVectorsScale_Shared_Win,BaseVectorsScale_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectorsScale_Shared_Win,IERROR)
BaseVectors0 => BaseVectors0_Shared
BaseVectors1 => BaseVectors1_Shared
BaseVectors2 => BaseVectors2_Shared
BaseVectors3 => BaseVectors3_Shared
BaseVectorsScale => BaseVectorsScale_Shared

firstSide = INT(REAL (myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE( BaseVectors0(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors1(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors2(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors3(1:3,1:nNonUniqueGlobalSides),&
          BaseVectorsScale(1:nNonUniqueGlobalSides))

firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

IF (BezierElevation.GT.0) THEN
  DO iSide=firstSide,lastSide
    BaseVectors0(:,iSide) = (+BezierControlPoints3DElevated(:,0,0           ,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             +BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors1(:,iSide) = (-BezierControlPoints3DElevated(:,0,0           ,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             -BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors2(:,iSide) = (-BezierControlPoints3DElevated(:,0,0           ,iSide)-BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             +BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    BaseVectors3(:,iSide) = (+BezierControlPoints3DElevated(:,0,0           ,iSide)-BezierControlPoints3DElevated(:,NGeoElevated,0           ,iSide)   &
                             -BezierControlPoints3DElevated(:,0,NGeoElevated,iSide)+BezierControlPoints3DElevated(:,NGeoElevated,NGeoElevated,iSide) )
    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
ELSE
  DO iSide=firstSide,lastSide
    BaseVectors0(:,iSide) = (+BezierControlPoints3D(:,0,0   ,iSide)+BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors1(:,iSide) = (-BezierControlPoints3D(:,0,0   ,iSide)+BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors2(:,iSide) = (-BezierControlPoints3D(:,0,0   ,iSide)-BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    BaseVectors3(:,iSide) = (+BezierControlPoints3D(:,0,0   ,iSide)-BezierControlPoints3D(:,NGeo,0   ,iSide)   &
                             -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
END IF

#if USE_MPI
CALL MPI_WIN_SYNC(BaseVectors0_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors1_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors2_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors3_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectorsScale_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI */

SWRITE(UNIT_stdOut,'(A)')' GET LINEAR SIDE BASEVECTORS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE GetLinearSideBaseVectors


SUBROUTINE GetMeshMinMax()
!===================================================================================================================================
! computes the minimum and maximum value of the mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Mesh_Vars          ,ONLY: NodeCoords
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
#if USE_MPI
USE MOD_Particle_Mesh_Vars ,ONLY: NodeCoords_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,POINTER                   :: NodeCoordsPointer(:,:,:,:,:)
!===================================================================================================================================

NodeCoordsPointer => NodeCoords
GEO%xmin     = MINVAL(NodeCoordsPointer(1,:,:,:,:))
GEO%xmax     = MAXVAL(NodeCoordsPointer(1,:,:,:,:))
GEO%ymin     = MINVAL(NodeCoordsPointer(2,:,:,:,:))
GEO%ymax     = MAXVAL(NodeCoordsPointer(2,:,:,:,:))
GEO%zmin     = MINVAL(NodeCoordsPointer(3,:,:,:,:))
GEO%zmax     = MAXVAL(NodeCoordsPointer(3,:,:,:,:))

#if USE_MPI
GEO%CNxmin   = MINVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNxmax   = MAXVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNymin   = MINVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNymax   = MAXVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNzmin   = MINVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%CNzmax   = MAXVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
GEO%xminglob = MINVAL(NodeCoords_Shared(1,:))
GEO%xmaxglob = MAXVAL(NodeCoords_Shared(1,:))
GEO%yminglob = MINVAL(NodeCoords_Shared(2,:))
GEO%ymaxglob = MAXVAL(NodeCoords_Shared(2,:))
GEO%zminglob = MINVAL(NodeCoords_Shared(3,:))
GEO%zmaxglob = MAXVAL(NodeCoords_Shared(3,:))
#else
GEO%CNxmin   = GEO%xmin
GEO%CNxmax   = GEO%xmax
GEO%CNymin   = GEO%ymin
GEO%CNymax   = GEO%ymax
GEO%CNzmin   = GEO%zmin
GEO%CNzmax   = GEO%zmax
GEO%xminglob = GEO%xmin
GEO%xmaxglob = GEO%xmax
GEO%yminglob = GEO%ymin
GEO%ymaxglob = GEO%ymax
GEO%zminglob = GEO%zmin
GEO%zmaxglob = GEO%zmax
#endif /*USE_MPI*/

END SUBROUTINE GetMeshMinMax


SUBROUTINE FinalizeParticleMesh()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
!USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Tracking_Vars ,ONLY: Distance,ListDistance
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                             :: iELem,iNode
!===================================================================================================================================

!SDEALLOCATE(PartElemToSide)
!SDEALLOCATE(PartSideToElem)
!SDEALLOCATE(PartElemIsMortar)
!SDEALLOCATE(PartElemToElemGlob)
!SDEALLOCATE(PartElemToElemAndSide)
SDEALLOCATE(PartBCSideList)
SDEALLOCATE(SidePeriodicType)
SDEALLOCATE(SidePeriodicDisplacement)
!SDEALLOCATE(IsTracingBCElem)
!SDEALLOCATE(TracingBCInnerSides)
!SDEALLOCATE(TracingBCTotalSides)
SDEALLOCATE(ElemType)
SDEALLOCATE(GEO%PeriodicVectors)
!SDEALLOCATE(GEO%PeriodicVectorsLength)
SDEALLOCATE(GEO%FIBGM)
!SDEALLOCATE(GEO%Volume)
!SDEALLOCATE(GEO%MPVolumePortion)
!SDEALLOCATE(GEO%CharLength)
SDEALLOCATE(GEO%ElemToFIBGM)
SDEALLOCATE(GEO%TFIBGM)

!SDEALLOCATE(GEO%ElemToNodeID)
!SDEALLOCATE(GEO%ElemSideNodeID)
!SDEALLOCATE(GEO%ElemToNodeIDGlobal)
!SDEALLOCATE(GEO%NodeCoords)
!SDEALLOCATE(GEO%ElemsOnNode)
!SDEALLOCATE(GEO%NeighNodesOnNode)
!SDEALLOCATE(GEO%NumNeighborElems)
!IF (ALLOCATED(GEO%ElemToNeighElems)) THEN
!  DO iElem=1,nElems
!    SDEALLOCATE(GEO%ElemToNeighElems(iElem)%ElemID)
!  END DO
!END IF
!SDEALLOCATE(GEO%ElemToNeighElems)
!IF (ALLOCATED(GEO%NodeToElem)) THEN
!  DO iNode=1,nNodes
!    SDEALLOCATE(GEO%NodeToElem(iNode)%ElemID)
!  END DO
!END IF
!SDEALLOCATE(GEO%NodeToElem)
!IF (ALLOCATED(GEO%NodeToNeighNode)) THEN
!  DO iNode=1,nNodes
!    SDEALLOCATE(GEO%NodeToNeighNode(iNode)%ElemID)
!  END DO
!END IF
!SDEALLOCATE(GEO%NodeToNeighNode)
!SDEALLOCATE(GEO%ConcaveElemSide)
!SDEALLOCATE(GEO%ElemMidPoint)
!SDEALLOCATE(GEO%BoundsOfElem)

SDEALLOCATE(BCElem)
MDEALLOCATE(XiEtaZetaBasis)
MDEALLOCATE(slenXiEtaZetaBasis)
MDEALLOCATE(ElemRadiusNGeo)
MDEALLOCATE(ElemRadius2NGeo)
MDEALLOCATE(ElemEpsOneCell)
SDEALLOCATE(Distance)
SDEALLOCATE(ListDistance)
!SDEALLOCATE(isTracingTrouble)
SDEALLOCATE(ElemTolerance)
SDEALLOCATE(ElemToGlobalElemID)
!SDEALLOCATE(ElemHaloInfoProc)
MDEALLOCATE(ElemCurved)

ParticleMeshInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMesh


SUBROUTINE MarkAuxBCElems()
!===================================================================================================================================
! check if auxBCs are inside BoundingBox of Elems
! -- plane: use plane equation f=a1*x+a2*y+a3*z+a4=0 and insert corresponding intervals of box -> fmin and fmax
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars                          ,ONLY:nElems
USE MOD_Particle_Globals
USE MOD_Particle_Mesh_Vars                 ,ONLY:ElemHasAuxBCs
USE MOD_Particle_Mesh_Vars                 ,ONLY: BoundsOfElem_Shared
USE MOD_Particle_Boundary_Vars             ,ONLY:nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone!,AuxBC_parabol
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,iAuxBC,icoord,dir(3),positiontype,positiontype_tmp
REAL                     :: r_vec(3),n_vec(3),fmin,fmax,radius,BoundsBC(1:2,1:3)
REAL                     :: lmin,lmax,deltamin,deltamax,origin(2),halfangle
LOGICAL                  :: cartesian, backwards
!===================================================================================================================================

ALLOCATE(ElemHasAuxBCs(1:nElems , 1:nAuxBCs))
ElemHasAuxBCs=.FALSE.

DO iAuxBC=1,nAuxBCs
  SELECT CASE (TRIM(AuxBCType(iAuxBC)))
  CASE ('plane')
    r_vec=AuxBC_plane(AuxBCMap(iAuxBC))%r_vec
    n_vec=AuxBC_plane(AuxBCMap(iAuxBC))%n_vec
    radius=AuxBC_plane(AuxBCMap(iAuxBC))%radius
    ! loop over all  elements
    DO iElem=1,nElems
      ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
      fmin=-DOT_PRODUCT(r_vec,n_vec)
      fmax=fmin
      DO icoord=1,3
        IF (n_vec(icoord).GE.0) THEN
          fmin = fmin + n_vec(icoord)*Bounds(1,icoord)
          fmax = fmax + n_vec(icoord)*Bounds(2,icoord)
        ELSE
          fmin = fmin + n_vec(icoord)*Bounds(2,icoord)
          fmax = fmax + n_vec(icoord)*Bounds(1,icoord)
        END IF
      END DO
      IF ((fmin.LE.0 .AND. fmax.GT.0).OR.(fmin.LT.0 .AND. fmax.GE.0)) THEN !plane intersects the box!
        !radius check needs to be implemented (compute intersection polygon and minimum radii): would sort out further elements!!!
        !quick, conservative solution: calculate bounding box of disc in space and compare with bb of element
        ElemHasAuxBCs(iElem,iAuxBC)=.TRUE.
        IF (radius .LT. 0.5*HUGE(radius)) THEN !huge was default
          BoundsBC(1,1:3) = r_vec - radius * SQRT(1.-(n_vec*n_vec))
          BoundsBC(2,1:3) = r_vec + radius * SQRT(1.-(n_vec*n_vec))
          DO icoord=1,3
            IF ( BoundsBC(2,icoord).LT.Bounds(1,icoord) .OR. BoundsBC(1,icoord).GT.Bounds(2,icoord) ) THEN
              ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
              EXIT
            END IF
          END DO
        END IF
      ELSE IF ((fmin.LT.0 .AND. fmax.LT.0).OR.(fmin.GT.0 .AND. fmax.GT.0)) THEN !plane does not intersect the box!
        ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
      ELSE !e.g. if elem has zero volume...
        CALL abort(&
          __STAMP__&
          ,'Error in MarkAuxBCElems for AuxBC:',iAuxBC)
      END IF
      END ASSOCIATE
    END DO
  CASE ('cylinder','cone')
    IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
      r_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%r_vec
      n_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%axis
      radius=AuxBC_cylinder(AuxBCMap(iAuxBC))%radius
      lmin=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin
      lmax=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmax
    ELSE !cone
      r_vec=AuxBC_cone(AuxBCMap(iAuxBC))%r_vec
      n_vec=AuxBC_cone(AuxBCMap(iAuxBC))%axis
      halfangle=AuxBC_cone(AuxBCMap(iAuxBC))%halfangle
      lmin=AuxBC_cone(AuxBCMap(iAuxBC))%lmin
      lmax=AuxBC_cone(AuxBCMap(iAuxBC))%lmax
    END IF
    cartesian=.TRUE.
    backwards=.FALSE.
    IF (ABS(n_vec(1)).EQ.1.) THEN
      dir(1)=1
      dir(2)=2
      dir(3)=3
      IF (n_vec(1).LT.0.) backwards=.TRUE.
    ELSE IF (ABS(n_vec(2)).EQ.1.) THEN
      dir(1)=2
      dir(2)=3
      dir(3)=1
      IF (n_vec(2).LT.0.) backwards=.TRUE.
    ELSE IF (ABS(n_vec(3)).EQ.1.) THEN
      dir(1)=3
      dir(2)=1
      dir(3)=2
      IF (n_vec(3).LT.0.) backwards=.TRUE.
    ELSE
      cartesian=.FALSE.
      SWRITE(*,*) 'WARNING in MarkAuxBCElems: all Elems are set to ElemHasAuxBCs=.TRUE. for AuxBC:',iAuxBC
      ElemHasAuxBCs(:,iAuxBC)=.TRUE. !actual intersection with box check to-be implemented!!!
    END IF
    IF (cartesian) THEN
      IF (backwards) THEN
        deltamin = -lmax
        deltamax = -lmin
      ELSE
        deltamin = lmin
        deltamax = lmax
      END IF
      origin(1) = r_vec(dir(2))
      origin(2) = r_vec(dir(3))
      ! loop over all  elements
      DO iElem=1,nElems
        ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
        ! check for lmin and lmax
        IF ( r_vec(dir(1))+deltamax.LT.Bounds(1,dir(1)) .OR. r_vec(dir(1))+deltamin.GT.Bounds(2,dir(1)) ) THEN
          ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
        ELSE !between lmin and lmax
          IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
            CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
          ELSE !cone
            !local minimum radius
            IF (backwards) THEN
              radius = MAX(-Bounds(2,dir(1))+r_vec(dir(1)),lmin)*TAN(halfangle)
            ELSE
              radius = MAX(Bounds(1,dir(1))-r_vec(dir(1)),lmin)*TAN(halfangle)
            END IF
            CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype_tmp)
            !local maximum radius
            IF (backwards) THEN
              radius = MIN(-Bounds(1,dir(1))+r_vec(dir(1)),lmax)*TAN(halfangle)
            ELSE
              radius = MIN(Bounds(2,dir(1))-r_vec(dir(1)),lmax)*TAN(halfangle)
            END IF
            CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
            !if both are type 0 or both are type 1 than the "total" type is not 2:
            IF ( .NOT.(positiontype_tmp.EQ.0 .AND. positiontype.EQ.0) &
              .AND. .NOT.(positiontype_tmp.EQ.1 .AND. positiontype.EQ.1) ) THEN
              positiontype=2
            END IF
          END IF
          IF (positiontype.EQ.2) THEN
            ElemHasAuxBCs(iElem,iAuxBC)=.TRUE.
          ELSE
            ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
          END IF
        END IF !check for lmin and lmax
        END ASSOCIATE
      END DO !iElem
    END IF !cartesian
  CASE('parabol')
    ElemHasAuxBCs(:,iAuxBC)=.TRUE. ! to be implemented!!!
  CASE DEFAULT
    SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
    CALL abort(&
      __STAMP__&
      ,'AuxBC does not exist')
  END SELECT
END DO

END SUBROUTINE MarkAuxBCElems


SUBROUTINE CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
!===================================================================================================================================
! checks how a cartesian bb is located with regard to a radius with cartesian axis (dir is cartesian axis and origin in orth. dirs)
!- positiontype=0 : complete bb is inside of radius
!- positiontype=1 : complete bb is outside of radius
!- positiontype=2 : bb is partly inside of radius
! (based on "check where the sides are located relative to rmax" in particle_emission for SimpleRadialVeloFit)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
!
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)           :: Bounds(1:2,1:3), origin(2), radius
INTEGER,INTENT(IN)        :: dir(3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)       :: positiontype
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iDir1, iDir2, iDir3, iPoint
REAL                      :: BoundingBox(1:3,1:8), point(2), pointRadius
LOGICAL                   :: done, insideBound
!===================================================================================================================================
!-- convert minmax-values to bb-points
DO iDir1=0,1
  DO iDir2=0,1
      DO iDir3=0,1
        BoundingBox(1,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir1+1,1)
        BoundingBox(2,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir2+1,2)
        BoundingBox(3,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir3+1,3)
      END DO
  END DO
END DO

!-- check where the points are located relative to radius
done=.FALSE.
DO iDir1=0,1
  IF(done) EXIT
  DO iDir2=0,1
    IF(done) EXIT
    DO iDir3=0,1
      !-- coords orth. to axis of point:
      iPoint=iDir1*4 + iDir2*2 + iDir3+1
      point(1) = BoundingBox(dir(2),iPoint)-origin(1)
      point(2) = BoundingBox(dir(3),iPoint)-origin(2)
      pointRadius = SQRT( (point(1))**2+(point(2))**2 )
      IF (iPoint.EQ.1) THEN
        IF (pointRadius.LE.radius) THEN
          insideBound=.TRUE.
        ELSE !outside
          insideBound=.FALSE.
        END IF !in-/outside?
      ELSE !iPoint.GT.1: type must be 2 if state of point if different from last point
        IF (pointRadius.LE.radius) THEN
          IF (.NOT.insideBound) THEN !different from last point
            positiontype=2
            done=.TRUE.
            EXIT
        END IF
        ELSE !outside
          IF (insideBound) THEN !different from last point
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        END IF !in-/outside?
      END IF !iPoint.EQ.1
    END DO !iDir3
  END DO !iDir2
END DO !iDir1
IF (.NOT.done) THEN
  IF (insideBound) THEN
    positiontype=0
  ELSE
    ! all points are outside of radius, but when radius is smaller than box, it can intersect it:
    IF ( origin(1) + radius .GE. Bounds(1,dir(2)) .AND. &
         origin(1) - radius .LE. Bounds(2,dir(2)) .AND. &
         origin(2) + radius .GE. Bounds(1,dir(3)) .AND. &
         origin(2) - radius .LE. Bounds(2,dir(3)) ) THEN !circle completely or partly inside box
      positiontype=2
    ELSE !points are really outside
      positiontype=1
    END IF
  END IF
END IF

END SUBROUTINE CheckBoundsWithCartRadius


FUNCTION GetGlobalNonUniqueSideID(ElemID,localSideID)
!===================================================================================================================================
!> Determines the non-unique global side ID of the local side in global element ElemID
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars   ,ONLY: ElemInfo_Shared, SideInfo_Shared
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID                              !< global element ID
INTEGER,INTENT(IN) :: localSideID                         !< local side id of an element (1:6)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER :: GetGlobalNonUniqueSideID
INTEGER :: iSide,firstSide,lastSide
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + 1
lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND, ElemID)

! Small mortar sides are added after
DO iSide = firstSide,lastSide
  IF (SideInfo_Shared(SIDE_LOCALID,iSide).EQ.localSideID) THEN
    GetGlobalNonUniqueSideID = iSide
    RETURN
  END IF
END DO

! We should never arrive here
CALL ABORT(__STAMP__,'GlobalSideID not found for Elem',ElemID)
END FUNCTION GetGlobalNonUniqueSideID

END MODULE MOD_Particle_Mesh
