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

PUBLIC :: DefineParametersParticleMesh
PUBLIC :: InitParticleMeshBasis
PUBLIC :: InitParticleMesh
PUBLIC :: InitParticleGeometry
PUBLIC :: FinalizeParticleMesh
PUBLIC :: MapRegionToElem
PUBLIC :: MarkAuxBCElems
PUBLIC :: GetMeshMinMax
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle tracking
!==================================================================================================================================
SUBROUTINE DefineParametersParticleMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Halo region
CALL prms%CreateRealOption(         'Part-SafetyFactor'         , 'Factor to scale the halo region with MPI'                       &
                                                                , '1.')
CALL prms%CreateRealOption(         'Part-HaloEpsVelo'          , 'Maximum velocity to be considered for halo region'              &
                                                                , '0.')
CALL prms%CreateLogicalOption(      'CalcHaloInfo'              , 'Output halo info to ElemData','.FALSE.')
! Background mesh init variables
CALL prms%CreateRealArrayOption(    'Part-FIBGMdeltas'          , 'Define the deltas for the cartesian Fast-Init-Background-Mesh.'//&
                                                                  ' They should be of the similar size as the smallest cells of' //&
                                                                  ' the used mesh for simulation.'                                 &
                                                                , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption(    'Part-FactorFIBGM'          , 'Factor with which the background mesh will be scaled.'          &
                                                                , '1. , 1. , 1.')

! Periodic vectors
CALL prms%CreateIntOption(          'Part-nPeriodicVectors'     , 'Number of the periodic vectors j=1,...,n.'                    //&
                                                                  ' Value has to be the same as defined in preprog.ini'            &
                                                                , '0')
CALL prms%CreateLogicalOption(      'CartesianPeriodic'         , ' Simplified treatment for periodic box with Refmapping. Not'  //&
                                                                  ' computation of intersection points at periodic BCs.'           &
                                                                , '.FALSE.')
CALL prms%CreateLogicalOption(      'FastPeriodic'              , ' Further simplification by directly moving particle into'     //&
                                                                  ' grid. Instead of moving the particle several times the'      //&
                                                                  ' periodic displacements, the particle is mapped directly back'//&
                                                                  ' into the domain.'                                              &
                                                                , '.FALSE.')

! RefMapping
CALL prms%CreateIntOption(          'RefMappingGuess'           , ' Initial guess of the Newton for mapping the particle into'   //&
                                                                  ' reference coordinates.\n'                                    //&
                                                                  '1 -linear pseudo-Cartesian coordinates\n'                     //&
                                                                  '2 - Xi of closest Gauss point\n'                              //&
                                                                  '3 - Xi of closest XCL_ngeo point\n'                           //&
                                                                  '4 -trival guess (0,0,0)^t')
CALL prms%CreateRealOption(         'RefMappingEps'             , ' Tolerance for mapping particle into reference element'       //&
                                                                  ' measured as L2-norm of deltaXi.'                               &
                                                                , '1e-4')

! Bezier surfaces
CALL prms%CreateIntOption(          'BezierElevation'           , ' Use BezierElevation>0 to tighten the bounding box. Typical'  //&
                                                                  ' values>10'                                                     &
                                                                , '0')
CALL prms%CreateIntOption(          'BezierSampleN'             , 'TODO-DEFINE-PARAMETER\n'                                      //&
                                                                  'Default value: NGeo equidistant sampling of bezier surface'   //&
                                                                  ' for emission'                                                  &
                                                                , '0')

CALL prms%CreateRealOption(         'BezierNewtonAngle'         , ' BoundingBox intersection angle for switching between '       //&
                                                                  'Bezierclipping and BezierNewton.'                               &
                                                                , '1.570796326')
CALL prms%CreateRealOption(         'BezierClipTolerance'       , ' Tolerance for BezierClipping'                                  &
                                                                , '1e-8')
CALL prms%CreateRealOption(         'BezierNewtonTolerance'     , ' Tolerance for BezierNewton'                                    &
                                                                , '1e-4')
CALL prms%CreateIntOption(          'BezierNewtonGuess'         , ' Initial guess for BezierNewton\n'                            //&
                                                                  ' 1 - linear projected face\n'                                 //&
                                                                  ' 2 - closest projected BeziercontrolPoint\n'                  //&
                                                                  ' 4 - (0,0)^t'                                                   &
                                                                , '1')
CALL prms%CreateIntOption(          'BezierNewtonMaxIter'       , ' TODO-DEFINE-PARAMETER'                                         &
                                                                , '100')
CALL prms%CreateRealOption(         'BezierSplitLimit'          , ' Limit for splitting in BezierClipping.'                      //&
                                                                  ' Value allows to detect multiple intersections and speed up ' //&
                                                                  'computation. Parameter is multiplied by 2'                      &
                                                                , '0.6')
CALL prms%CreateIntOption(          'BezierClipMaxIter'         , ' Max iteration of BezierClipping'                               &
                                                                , '100')
CALL prms%CreateRealOption(         'epsilontol'                , ' Tolerance (absolute) for comparison against zero'              &
                                                                , '0.')
CALL prms%CreateRealOption(         'epsilonrel'                , ' Tolerance (relative) for comparison against zero'              &
                                                                , '0.')
CALL prms%CreateRealOption(         'BezierClipHit'             , ' Tolerance in [-1,1] of BezierFace'                             &
                                                                , '0.')
CALL prms%CreateRealOption(         'BezierNewtonHit'           , ' Tolerance in [-1,1] of BezierNewton'                           &
                                                                , '0.')
CALL prms%CreateIntOption(          'BezierClipMaxIntersec'     , ' Max. number of multiple intersections. Default: 2*(NGeo+1)')

END SUBROUTINE DefineParametersParticleMesh


SUBROUTINE InitParticleMeshBasis()
!===================================================================================================================================
! Init basic metrics needed for particle mesh
!===================================================================================================================================
! MODULES
USE MOD_PreProc                ,ONLY: N
USE MOD_Basis
USE MOD_Interpolation_Vars     ,ONLY: xGP
USE MOD_Mesh_Vars
USE MOD_Particle_Basis
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars
USE MOD_ReadInTools            ,ONLY: GETINT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
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
CALL ChebyGaussLobNodesAndWeights(1                 ,XiCL_NGeo1)
CALL BarycentricWeights          (1                 ,XiCL_NGeo1,wBaryCL_NGeo1)
CALL InitializeVandermonde       (1, NGeo,wBaryCL_NGeo1,XiCL_NGeo1,XiCL_NGeo ,Vdm_CLNGeo1_CLNGeo)

! initialize Vandermonde for Bezier basis surface representation (particle tracking with curved elements)
BezierElevation = GETINT('BezierElevation','0')
NGeoElevated    = NGeo + BezierElevation

CALL BuildBezierVdm              (NGeo,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier)
CALL BuildBezierDMat             (NGeo,Xi_NGeo,D_Bezier)

! allocate Chebyshev-Lobatto physical coordinates
ALLOCATE( XCL_NGeo(1:3,    0:NGeo,0:NGeo,0:ZDIM(NGeo),1:nElems)   &
        ,dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:ZDIM(NGeo),1:nElems))
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
USE MOD_Particle_Interpolation_Vars, ONLY: DoInterpolation
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: InitGetGlobalElemID,InitGetCNElemID,GetCNElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: InitGetGlobalSideID,InitGetCNSideID,GetGlobalSideID
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation,BezierControlPoints3D,BezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Tracking_Vars ,ONLY: FastPeriodic,CountNbOfLostParts,NbrOfLostParticles,NbrOfLostParticlesTotal,CartesianPeriodic
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GetRealArray
#if CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: SideBoundingBoxVolume
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_Particle_BGM           ,ONLY: WriteHaloInfo
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars
#endif /* USE_MPI */
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#else
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: RefMappingGuessProposal
INTEGER          :: iSample
INTEGER          :: firstSide,lastSide,iSide,SideID
CHARACTER(LEN=2) :: tmpStr
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#else
INTEGER          :: ALLOCSTAT
#endif
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH...'
IF(ParticleMeshInitIsDone) CALL ABORT(__STAMP__, ' Particle-Mesh is already initialized.')

!===================================================================================================================================
! Collect particle variables that are not initialized somewhere else
!===================================================================================================================================
PP_nElems = nElems

# if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  CALL ComputePeriodicVec()
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

! ElemTime already set in loadbalance.f90
#if !USE_LOADBALANCE
SDEALLOCATE(ElemTime)
ALLOCATE(ElemTime(1:nElems))
ELemTime = 0.
#endif /*!USE_LOADBALANCE*/
!===================================================================================================================================

! Potentially curved elements. FIBGM needs to be built on BezierControlPoints rather than NodeCoords to avoid missing elements
IF (TrackingMethod.EQ.TRACING .OR. TrackingMethod.EQ.REFMAPPING) THEN
  CALL CalcParticleMeshMetrics()

  CALL CalcBezierControlPoints()
END IF

! Mesh min/max must be built on BezierControlPoint for possibly curved elements
CALL GetMeshMinMax()

! Build BGM to Element mapping and identify which of the elements, sides and nodes are in the compute-node local and halo region
CALL BuildBGMAndIdentifyHaloRegion()

#if USE_MPI
CalcHaloInfo = GETLOGICAL('CalcHaloInfo')
IF (CalcHaloInfo) THEN
  CALL WriteHaloInfo
END IF
#endif /*USE_MPI*/

! Initialize mapping function: GetGlobalElemID(1:nComputeNodeTotalElems)
CALL InitGetGlobalElemID()

! Initialize mapping function: GetCNElemID()
CALL InitGetCNElemID()

! Initialize mapping function: GetGlobalSideID(1:nComputeNodeTotalSides)
CALL InitGetGlobalSideID()

! Initialize mapping function: GetCNSideID(1:GlobalSideID)
CALL InitGetCNSideID()

CountNbOfLostParts      = GETLOGICAL('CountNbOfLostParts',".FALSE.")
IF (CountNbOfLostParts) THEN
  NbrOfLostParticles      = 0
  NbrOfLostParticlesTotal = 0
!  DisplayLostParticles    = GETLOGICAL('DisplayLostParticles')
END IF

#if CODE_ANALYZE
PARTOUT            = GETINT('PartOut','0')
MPIRankOut         = GETINT('MPIRankOut','0')
#endif /*CODE_ANALYZE*/

CartesianPeriodic = GETLOGICAL('CartesianPeriodic','.FALSE.')
IF (CartesianPeriodic) FastPeriodic = GETLOGICAL('FastPeriodic','.FALSE.')

IF (UseCurveds) THEN ! don't use RefMappingGuess=1, because RefMappingGuess is only best for linear cubical elements
  ! curved elements can be stronger deformed, hence, a better guess can be used
  ! RefMappingGuess 2,3 searches the closest Gauss/CL points of the considered element. This point is used as the initial value for
  ! the mapping. Note, that the position of the CL points can still be advantageous for the initial guess.
  RefMappingGuessProposal = 2
  ! there are more Gauss points within an element then CL-points. Gauss points sample the element finer
  IF (PP_N.GT.NGeo) THEN
    RefMappingGuessProposal = 2
  ! more CL-points than Gauss points, hence, better sampling of the element
  ELSE
    RefMappingGuessProposal = 3
  END IF
! default for linear meshes. Guess is exact for cubical, non-twisted elements
ELSE
  RefMappingGuessProposal = 1
END IF

WRITE(tmpStr,'(I2.2)') RefMappingGuessProposal
RefMappingGuess = GETINT('RefMappingGuess',tmpStr)
! Linear intial guess on curved meshes might cause problems.
IF((RefMappingGuess.EQ.1).AND.(UseCurveds)) THEN
  SWRITE(UNIT_stdOut,'(A)')' WARNING: read-in [RefMappingGuess=1] when using [UseCurveds=T] may create problems!'
END IF

RefMappingEps   = GETREAL('RefMappingEps','1e-4')
epsInCell       = SQRT(3.0*RefMappingEps)

IF((RefMappingGuess.LT.1).OR.(RefMappingGuess.GT.4))THEN
   CALL ABORT(__STAMP__ ,'Wrong guessing method for mapping from physical space in reference space.',RefMappingGuess,999.)
END IF

CALL InitElemVolumes()

! TriaSurfaceFlux needs to be determined before particle mesh init to build all required information. TriaSurfaceFlux is enabled by
! default for TriaTracking and disabled otherwise
WRITE(tmpStr,'(L1)') (TrackingMethod.EQ.TRIATRACKING)
TriaSurfaceFlux = GETLOGICAL('TriaSurfaceFlux',TRIM(tmpStr))
IF((TrackingMethod.EQ.TRIATRACKING).AND.(.NOT.TriaSurfaceFlux)) THEN
  CALL ABORT(__STAMP__,'TriaSurfaceFlux explicitly turned off for TriaTracking!')
END IF

SELECT CASE(TrackingMethod)

  CASE(TRIATRACKING)
    CALL InitParticleGeometry()
    ! Compute convex element radius^2
    CALL BuildElementRadiusTria()

    ! Interpolation needs coordinates in reference system
    IF (DoInterpolation) THEN
      CALL CalcParticleMeshMetrics()

      CALL BuildElemTypeAndBasisTria()
    END IF

  CASE(TRACING,REFMAPPING)
    ! TriaSurfaceFlux needs ElemNodeID_Shared
    IF (TriaSurfaceFlux) THEN
      CALL InitParticleGeometry()
    END IF

     ! Calculated before BGM to find convex hull
!    CALL CalcParticleMeshMetrics()

!    CALL CalcBezierControlPoints()

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
      DO iSide = firstSide,LastSide
        ! ignore sides that are not on the compute node
        ! IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

        SideID = GetGlobalSideID(iSide)

        ! Ignore small mortar sides attached to big mortar sides
        IF (SideInfo_Shared(SIDE_LOCALID,SideID).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,SideID).GT.6) CYCLE

        ! BezierControlPoints are always on nonUniqueGlobalSide
        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,SideID) &
                                           ,SideSlabNormals(   1:3,1:3,iSide)                                       &
                                           ,SideSlabInterVals( 1:6    ,iSide)                                       &
                                           ,BoundingBoxIsEmpty(iSide))
      END DO
    ELSE
      DO iSide=firstSide,LastSide
        ! ignore sides that are not on the compute node
        ! IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

        SideID = GetGlobalSideID(iSide)

        ! Ignore small mortar sides attached to big mortar sides
        IF (SideInfo_Shared(SIDE_LOCALID,SideID).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,SideID).GT.6) CYCLE

        ! BezierControlPoints are always on nonUniqueGlobalSide
        CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)                         &
                                           ,SideSlabNormals(   1:3,1:3,iSide)                                       &
                                           ,SideSlabInterVals( 1:6    ,iSide)                                       &
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
    DO iSide = offsetMPISides_YOUR,LastSide
      dx = ABS(SideSlabIntervals(2)-SideSlabIntervals(1))
      dy = ABS(SideSlabIntervals(4)-SideSlabIntervals(3))
      dz = ABS(SideSlabIntervals(6)-SideSlabIntervals(5))
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

    IF (TrackingMethod.EQ.REFMAPPING) THEN
      ! Identify BCSides and build side origin and radius
      CALL GetBCSidesAndOrgin()

      ! Identify BCElems
      CALL BuildBCElemDistance()
    END IF

    CALL BuildEpsOneCell()

  CASE DEFAULT
    CALL ABORT(__STAMP__,'Invalid tracking method in particle_mesh.f90!')

END SELECT

! Equidistant BezierControlPoints, required for BezierAreaSample. Dummy in case surfaceFlux is calculated on triangles
IF (TriaSurfaceFlux) THEN
  BezierSampleN    = 1
  SurfFluxSideSize = (/1,2/)
ELSE
  WRITE(tmpStr,'(I2.2)') NGeo
  BezierSampleN    = GETINT('BezierSampleN',tmpStr)
  SurfFluxSideSize = BezierSampleN
  ALLOCATE(BezierSampleXi(0:BezierSampleN))
  DO iSample = 0,BezierSampleN
    BezierSampleXi(iSample) = -1. + 2.0/BezierSampleN*iSample
  END DO
END IF

ParticleMeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMesh


SUBROUTINE ComputePeriodicVec()
!===================================================================================================================================
! Init of Particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars                ,ONLY: NGeo,BoundaryType
USE MOD_Particle_Boundary_Vars   ,ONLY: PartBound
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars       ,ONLY: GEO,ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared,nNonUniqueGlobalSides
#if USE_MPI
USE MOD_Mesh_Vars                ,ONLY: nGlobalElems
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: myComputeNodeRank,nComputeNodeProcessors,MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars                ,ONLY: nElems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: iNode=1
INTEGER                        :: firstElem,lastElem,iNbSide,iVec,BCALPHA
INTEGER                        :: SideID,ElemID,GlobalSideID,NbElemID,localSideID,localSideNbID,nStart
INTEGER                        :: CornerNodeIDswitch(8),NodeMap(4,6)
REAL,DIMENSION(3)              :: MasterCoords,SlaveCoords,PeriodicTmp
LOGICAL,ALLOCPOINT             :: PeriodicFound(:)
#if USE_MPI
REAL,ALLOCATABLE               :: sendbuf(:),recvbuf(:,:)
INTEGER                        :: PeriodicFound_Win
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

! the cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
CornerNodeIDswitch(1)=1
CornerNodeIDswitch(2)=(Ngeo+1)
CornerNodeIDswitch(3)=(Ngeo+1)**2
CornerNodeIDswitch(4)=(Ngeo+1)*Ngeo+1
CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2
CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1

! Corner node switch to order HOPR coordinates in CGNS format
ASSOCIATE(CNS => CornerNodeIDswitch )
! CGNS Mapping
NodeMap(:,1)=(/CNS(1),CNS(4),CNS(3),CNS(2)/)
NodeMap(:,2)=(/CNS(1),CNS(2),CNS(6),CNS(5)/)
NodeMap(:,3)=(/CNS(2),CNS(3),CNS(7),CNS(6)/)
NodeMap(:,4)=(/CNS(3),CNS(4),CNS(8),CNS(7)/)
NodeMap(:,5)=(/CNS(1),CNS(5),CNS(8),CNS(4)/)
NodeMap(:,6)=(/CNS(5),CNS(6),CNS(7),CNS(8)/)

! Find number of periodic vectors
GEO%nPeriodicVectors = MAXVAL(BoundaryType(:,BC_ALPHA))
IF (GEO%nPeriodicVectors.EQ.0) RETURN

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank*   nGlobalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nGlobalElems)/REAL(nComputeNodeProcessors))

MPISharedSize = INT(GEO%nPeriodicVectors,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
! Somehow the pointer is associated at this point, nullify it
NULLIFY(PeriodicFound)
CALL Allocate_Shared(MPISharedSize,(/GEO%nPeriodicVectors/),PeriodicFound_Win,PeriodicFound)
CALL MPI_WIN_LOCK_ALL(0,PeriodicFound_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
firstElem = 1
lastElem  = nElems

ALLOCATE(PeriodicFound(1:GEO%nPeriodicVectors))
#endif /*USE_MPI*/

! Only root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif  /*USE_MPI*/
  PeriodicFound(:) = .FALSE.
#if USE_MPI
END IF
CALL MPI_WIN_SYNC(PeriodicFound_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

ALLOCATE(GEO%PeriodicVectors(1:3,GEO%nPeriodicVectors))
GEO%PeriodicVectors = 0.

DO ElemID = firstElem,lastElem
  ! Every periodic vector already found
  IF (ALL(PeriodicFound(:))) EXIT

  DO SideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)
    IF (SideInfo_Shared(SIDE_BCID,SideID).EQ.0) CYCLE

    ! Boundary is a periodic boundary
    IF (PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,SideID)).NE.3) CYCLE

    ! Check if side is master side
    BCALPHA = BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA)

    IF (BCALPHA.GT.0) THEN
      ! Periodic vector already found
      IF (PeriodicFound(BCALPHA)) CYCLE

      ! Periodic slave side has same ID, but negative sign
      GlobalSideID = SideInfo_Shared(SIDE_ID,SideID)
      DO iNbSide = 1,nNonUniqueGlobalSides
        IF (SideInfo_Shared(SIDE_ID,iNbSide).EQ.-GlobalSideID) THEN
          NbElemID      = SideInfo_Shared(SIDE_ELEMID,iNbSide)
          localSideID   = SideInfo_Shared(SIDE_LOCALID,SideID)
          localSideNbID = SideInfo_Shared(SIDE_LOCALID,iNbSide)
          nStart        = MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,iNbSide),10)-1)

          ! Only take the first node into account, no benefit in accuracy if running over others as well
!          DO iNode = 1,4
            MasterCoords = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)  +NodeMap(iNode                  ,localSideID))
            SlaveCoords  = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,NbElemID)+NodeMap(MOD(nStart+5-iNode,4)+1,localSideNbID))
            PeriodicTmp  = SlaveCoords - MasterCoords
            GEO%PeriodicVectors(:,BCALPHA) = PeriodicTmp
            PeriodicFound(BCALPHA) = .TRUE.
            EXIT
!          END DO
        END IF
      END DO
    END IF
  END DO
END DO

END ASSOCIATE

#if USE_MPI
ALLOCATE(sendbuf(GEO%nPeriodicVectors)                           ,&
         recvbuf(GEO%nPeriodicVectors,0:nComputeNodeProcessors-1))
sendbuf = 0.
recvbuf = 0.

DO iVec = 1,GEO%nPeriodicVectors
  sendbuf(iVec) = MERGE(VECNORM(GEO%PeriodicVectors(:,iVec)),HUGE(1.),VECNORM(GEO%PeriodicVectors(:,iVec)).GT.0)
END DO

! Do it by hand, MPI_ALLREDUCE seems problematic with MPI_2DOUBLE_PRECISION and MPI_MINLOC
! https://stackoverflow.com/questions/56307320/mpi-allreduce-not-synchronizing-properly
!CALL MPI_ALLREDUCE(MPI_IN_PLACE,sendbuf,GEO%nPeriodicVectors,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_SHARED,iERROR)
DO iVec = 1,GEO%nPeriodicVectors
  CALL MPI_ALLGATHER(sendbuf(iVec),1,MPI_DOUBLE_PRECISION,recvbuf(iVec,:),1,MPI_DOUBLE_PRECISION,MPI_COMM_SHARED,iERROR)
  ! MINLOC does not follow array bounds, so root rank = 1
  CALL MPI_BCAST(GEO%PeriodicVectors(:,iVec),3,MPI_DOUBLE_PRECISION,MINLOC(recvbuf(iVec,:),1)-1,MPI_COMM_SHARED,iError)
END DO

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
! Deallocate only after MPI_BCAST safely returned
DEALLOCATE(sendbuf,recvbuf)
CALL MPI_WIN_UNLOCK_ALL(PeriodicFound_Win,iError)
CALL MPI_WIN_FREE(      PeriodicFound_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/
MDEALLOCATE(PeriodicFound)

SWRITE(UNIT_StdOut,'(A,I0,A)') ' | Found ',GEO%nPeriodicVectors,' periodic vectors for particle tracking'

END SUBROUTINE ComputePeriodicVec


SUBROUTINE CalcParticleMeshMetrics()
!===================================================================================================================================
!> calculates XCL_Ngeo and dXCL_Ngeo for compute node mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Basis                    ,ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights,InitializeVandermonde
USE MOD_ChangeBasis              ,ONLY: ChangeBasis3D
USE MOD_Interpolation            ,ONLY: GetDerivativeMatrix
USE MOD_Interpolation            ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars       ,ONLY: NodeType,NodeTypeCL,NodeTypeVISU
USE MOD_Mesh_Vars                ,ONLY: NGeo
USE MOD_Mesh_Vars                ,ONLY: nElems,nGlobalElems
USE MOD_Mesh_Vars                ,ONLY: Elem_xGP
USE MOD_Particle_Mesh_Vars       ,ONLY: Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: XCL_NGeo,dXCL_NGeo
USE MOD_Particle_Mesh_Vars       ,ONLY: XCL_NGeo_Shared,dXCL_NGeo_Shared
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_PreProc                  ,ONLY: N
USE MOD_Mesh_Vars                ,ONLY: offsetElem
!USE MOD_Mesh_Vars                ,ONLY: InterpolateFromTree
!USE MOD_Mesh_Vars                ,ONLY: UseCurveds
USE MOD_Metrics                  ,ONLY: BuildCoords
!USE MOD_Particle_Mesh_Vars       ,ONLY: nComputeNodeElems
!USE MOD_Particle_Mesh_Vars       ,ONLY: ElemInfo_Shared,NodeCoords_Shared,TreeCoords_Shared
!USE MOD_Particle_Mesh_Vars       ,ONLY: Vdm_NGeo_CLNGeo
USE MOD_Particle_Mesh_Vars       ,ONLY: XCL_NGeo_Array,dXCL_NGeo_Array
USE MOD_Particle_Mesh_Vars       ,ONLY: Elem_xGP_Array
USE MOD_Particle_Mesh_Vars       ,ONLY: XCL_NGeo_Shared_Win,dXCL_NGeo_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: Elem_xGP_Shared_Win
USE MOD_Particle_Mesh_Tools      ,ONLY: GetGlobalElemID
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
!USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,nProcessors_Global,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SHARED,displsElem,recvcountelem
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
#if USE_MPI
INTEGER                        :: iElem!,ElemID
!INTEGER                        :: firstHaloElem,lastHaloElem,nComputeNodeHaloElems
!INTEGER                        :: firstNodeID,nodeID,i,j,k,ll
!REAL                           :: NodeCoordstmp(1:3,0:NGeo,0:NGeo,0:NGeo)
!REAL                           :: DCL_NGeo(         0:Ngeo,0:Ngeo)
!REAL                           :: Vdm_EQNGeo_CLN(   0:PP_N,0:NGeo)
!REAL                           :: Vdm_CLNloc_N(     0:PP_N,0:PP_N)
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_LOADBALANCE
! XCL and dXCL are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

#if USE_MPI
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
MPISharedSize = INT((3*(NGeo+1)**3*nGlobalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*  (NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/), XCL_NGeo_Shared_Win,XCL_NGeo_Array)
MPISharedSize = INT((3*(PP_N+1)**3*nGlobalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*  (PP_N+1)*(PP_N+1)*(PP_N+1)*nGlobalElems/), Elem_xGP_Shared_Win,Elem_xGP_Array)
MPISharedSize = INT((3*3*(NGeo+1)**3*nGlobalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3*3*(NGeo+1)*(NGeo+1)*(NGeo+1)*nGlobalElems/),dXCL_NGeo_Shared_Win,dXCL_NGeo_Array)
CALL MPI_WIN_LOCK_ALL(0,XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,dXCL_NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
XCL_NGeo_Shared (1:3    ,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => XCL_NGeo_Array
Elem_xGP_Shared (1:3    ,0:PP_N,0:PP_N,0:PP_N,1:nGlobalElems) => Elem_xGP_Array
dXCL_NGeo_Shared(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nGlobalElems) => dXCL_NGeo_Array

! Copy local XCL and dXCL into shared memory
!IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  DO iElem = 1, nElems
    XCL_NGeo_Shared (:  ,:,:,:,offsetElem+iElem) = XCL_NGeo (:  ,:,:,:,iElem)
    Elem_xGP_Shared (:  ,:,:,:,offsetElem+iElem) = Elem_xGP (:  ,:,:,:,iElem)
    dXCL_NGeo_Shared(:,:,:,:,:,offsetElem+iElem) = dXCL_NGeo(:,:,:,:,:,iElem)
  END DO ! iElem = 1, nElems
!ELSE
!  DO iElem = 1, nElems
!    XCL_NGeo_Shared (:,  :,:,:,offsetElem+iElem) = XCL_NGeo (:,  :,:,:,iElem)
!    Elem_xGP_Shared (:,  :,:,:,offsetElem+iElem) = Elem_xGP (:,  :,:,:,iElem)
!    dXCL_NGeo_Shared(:,:,:,:,:,offsetElem+iElem) = dXCL_NGeo(:,:,:,:,:,iElem)
!  END DO ! iElem = 1, nElems
!END IF

! XCL_NGeo\dXCL_NGeo can be deallocated in MPI case
DEALLOCATE(XCL_Ngeo)
DEALLOCATE(dXCL_Ngeo)

! Communicate XCL and dXCL between compute node rootss instead of calculating globally
CALL MPI_WIN_SYNC(XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(dXCL_NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

IF (nComputeNodeProcessors.NE.nProcessors_Global .AND. myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , XCL_NGeo_Shared               &
                     , 3*(NGeo+1)**3*recvcountElem   &
                     , 3*(NGeo+1)**3*displsElem      &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)
!  CALL MPI_ALLREDUCE(MPI_IN_PLACE,XCL_NGeo_Shared,3*(NGeo+1)**3*4096,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)

  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , Elem_xGP_Shared               &
                     , 3*(PP_N+1)**3*recvcountElem   &
                     , 3*(PP_N+1)**3*displsElem      &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)
!  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Elem_xGP_Array,3*(PP_N+1)**3*4096,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)

  CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , dXCL_NGeo_Shared              &
                     , 3*3*(NGeo+1)**3*recvcountElem &
                     , 3*3*(NGeo+1)**3*displsElem    &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)
!  CALL MPI_ALLREDUCE(MPI_IN_PLACE,dXCL_NGeo_Shared,3*3*(NGeo+1)**3*4096,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)
END IF

!nComputeNodeHaloElems = nComputeNodeTotalElems - nComputeNodeElems
!IF (nComputeNodeHaloElems.GT.nComputeNodeProcessors) THEN
!  firstHaloElem = INT(REAL( myComputeNodeRank   *nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))+1
!  lastHaloElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeHaloElems)/REAL(nComputeNodeProcessors))
!ELSE
!  firstHaloElem = myComputeNodeRank + 1
!  IF (myComputeNodeRank.LT.nComputeNodeHaloElems) THEN
!    lastHaloElem = myComputeNodeRank + 1
!  ELSE
!    lastHaloElem = 0
!  END IF
!END IF
!
!! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
!!       Important for curved meshes if NGeo<N, no effect for N>=NGeo
!CALL GetVandermonde(    NGeo, NodeTypeVISU, PP_N, NodeTypeCL, Vdm_EQNGeo_CLN,  modal=.FALSE.)
!CALL GetVandermonde(    PP_N, NodeTypeCL  , PP_N, NodeType  , Vdm_CLNloc_N,    modal=.FALSE.)
!
!!Transform from EQUI_NGeo to solution points on Nloc
!Vdm_EQNGeo_CLN = MATMUL(Vdm_CLNloc_N,Vdm_EQNGeo_CLN)
!
!! Build XCL and dXCL for compute node halo region (each proc of compute-node build only its fair share)
!IF(interpolateFromTree) THEN
!  CALL abort(__STAMP__,'ERROR: InterpolateFromTree not yet implemented for new halo region!')
!ELSE
!  CALL GetDerivativeMatrix(NGeo  , NodeTypeCL  , DCL_Ngeo)
!
!  DO iElem = firstHaloElem,lastHaloElem
!    ElemID = GetGlobalElemID(nComputeNodeElems+iElem)
!!    firstNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+1
!    firstNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)
!!    nodeID = 0
!    nodeID = 1
!    IF (useCurveds) THEN
!      DO k = 0, NGeo; DO j = 0, NGeo; DO i = 0, NGeo
!        NodeCoordstmp(:,i,j,k) = NodeCoords_Shared(:,firstNodeID+NodeID)
!        nodeID = nodeID + 1
!      END DO; END DO; END DO ! i = 0, NGeo
!    ELSE
!      NodeCoordstmp(:,0,0,0) = NodeCoords_Shared(:,firstNodeID+1)
!      NodeCoordstmp(:,1,0,0) = NodeCoords_Shared(:,firstNodeID+2)
!      NodeCoordstmp(:,0,1,0) = NodeCoords_Shared(:,firstNodeID+3)
!      NodeCoordstmp(:,1,1,0) = NodeCoords_Shared(:,firstNodeID+4)
!      NodeCoordstmp(:,0,0,1) = NodeCoords_Shared(:,firstNodeID+5)
!      NodeCoordstmp(:,1,0,1) = NodeCoords_Shared(:,firstNodeID+6)
!      NodeCoordstmp(:,0,1,1) = NodeCoords_Shared(:,firstNodeID+7)
!      NodeCoordstmp(:,1,1,1) = NodeCoords_Shared(:,firstNodeID+8)
!    END IF
!    CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_NGeo_CLNGeo,NodeCoordstmp,XCL_NGeo_Shared(:,:,:,:,nComputeNodeElems+iElem))
!    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_EQNGeo_CLN ,NodeCoordstmp,Elem_xGP_Shared(:,:,:,:,nComputeNodeElems+iElem))
!
!    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
!      ! Matrix-vector multiplication
!      DO ll=0,Ngeo
!        dXCL_NGeo_Shared(1,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(1,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(i,ll) &
!                                                            *  XCL_NGeo_Shared(: ,ll,j,k,nComputeNodeElems+iElem)
!        dXCL_NGeo_Shared(2,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(2,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(j,ll) &
!                                                            *  XCL_NGeo_Shared(: ,i,ll,k,nComputeNodeElems+iElem)
!        dXCL_NGeo_Shared(3,:,i,j,k,nComputeNodeElems+iElem) = dXCL_NGeo_Shared(3,:,i,j,k,nComputeNodeElems+iElem) + DCL_NGeo(k,ll) &
!                                                            *  XCL_NGeo_Shared(: ,i,j,ll,nComputeNodeElems+iElem)
!      END DO
!    END DO; END DO; END DO !i,j,k=0,Ngeo
!    END DO ! iElem = firstHaloElem, lastHaloElem
!END IF

CALL MPI_WIN_SYNC(XCL_NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(Elem_xGP_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(dXCL_NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
XCL_NGeo_Shared  => XCL_NGeo
Elem_xGP_Shared  => Elem_xGP
dXCL_NGeo_Shared => dXCL_NGeo
#endif /*USE_MPI*/

END SUBROUTINE CalcParticleMeshMetrics


SUBROUTINE CalcBezierControlPoints()
!===================================================================================================================================
!> calculate the Bezier control point (+elevated) for shared compute-node mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis              ,ONLY: ChangeBasis2D
USE MOD_Mappings                 ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars                ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars       ,ONLY: nNonUniqueGlobalSides,SideInfo_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: NGeoElevated,XCL_NGeo_Shared
USE MOD_Particle_Mesh_Tools      ,ONLY: GetGlobalNonUniqueSideID!,GetGlobalElemID
USE MOD_Particle_Surfaces        ,ONLY: GetBezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierControlPoints3D,sVdm_Bezier
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierControlPoints3DElevated,BezierElevation
#if USE_MPI
USE MOD_Mesh_Vars                ,ONLY: nGlobalElems
USE MOD_Particle_Mesh_Vars       ,ONLY: BezierControlPoints3D_Shared,BezierControlPoints3D_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: BezierControlPoints3DElevated_Shared,BezierControlPoints3DElevated_Shared_Win
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
!USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars                ,ONLY: nElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iSide,ilocSide
INTEGER                        :: SideID
INTEGER                        :: firstElem,lastElem,firstSide,lastSide
INTEGER                        :: p,q,pq(2)
REAL                           :: tmp(3,0:NGeo,0:NGeo)
REAL                           :: tmp2(3,0:NGeo,0:NGeo)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#else
INTEGER                        :: ALLOCSTAT
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_LOADBALANCE
! BezierControlPoints are global and do not change during load balance, return
IF (PerformLoadBalance) RETURN
#endif

SWRITE(UNIT_stdOut,'(A)') ' CALCULATING BezierControlPoints...'

! Build BezierControlPoints3D (compute-node local+halo)
#if USE_MPI
MPISharedSize = INT((3*(NGeo+1)**2*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
! This is a trick. Allocate as 1D array and then set a pointer with the proper array bounds
CALL Allocate_Shared(MPISharedSize,(/3*(NGeo+1)*(NGeo+1)*nNonUniqueGlobalSides/),BezierControlPoints3D_Shared_Win,BezierControlPoints3D_Shared)
CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3D_Shared_Win,IERROR)
BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nNonUniqueGlobalSides) => BezierControlPoints3D_Shared
IF (myComputeNodeRank.EQ.0) THEN
  BezierControlPoints3D         = 0.
END IF
IF (BezierElevation.GT.0) THEN
  MPISharedSize = INT((3*(NGeoElevated+1)**2*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
  CALL Allocate_Shared(MPISharedSize,(/3*(NGeoElevated+1)*(NGeoElevated+1)*nNonUniqueGlobalSides/), &
                                      BezierControlPoints3DElevated_Shared_Win,BezierControlPoints3DElevated_Shared)
  CALL MPI_WIN_LOCK_ALL(0,BezierControlPoints3DElevated_Shared_Win,IERROR)
  BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nNonUniqueGlobalSides) => BezierControlPoints3DElevated_Shared
  IF (myComputeNodeRank.EQ.0) THEN
    BezierControlPoints3DElevated = 0.
  END IF
END IF
#else
ALLOCATE(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nNonUniqueGlobalSides) &
        ,STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3D!')
BezierControlPoints3D         = 0.

IF (BezierElevation.GT.0) THEN
  ALLOCATE(BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,1:nNonUniqueGlobalSides) &
          ,STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(__STAMP__,'  Cannot allocate BezierControlPoints3DElevated!')
  BezierControlPoints3DElevated = 0.
END IF
#endif /*USE_MPI*/

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
IF (BezierElevation.GT.0) THEN
  CALL MPI_WIN_SYNC(BezierControlPoints3DElevated_Shared_Win,IERROR)
END IF
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

firstElem = INT(REAL( myComputeNodeRank*   nGlobalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nGlobalElems)/REAL(nComputeNodeProcessors))
!firstSide = INT(REAL (myComputeNodeRank   *nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))+1
!lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalSides)/REAL(nComputeNodeProcessors))
firstSide = INT(REAL (myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

! iElem is CN elem
DO iElem = firstElem, lastElem
  DO ilocSide=1,6
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp = XCL_NGeo_Shared(1:3 , 0    , :    , :   ,iElem )
    CASE(XI_PLUS)
      tmp = XCL_NGeo_Shared(1:3 , NGeo , :    , :   ,iElem )
    CASE(ETA_MINUS)
      tmp = XCL_NGeo_Shared(1:3 , :    , 0    , :   ,iElem )
    CASE(ETA_PLUS)
      tmp = XCL_NGeo_Shared(1:3 , :    , NGeo , :   ,iElem )
    CASE(ZETA_MINUS)
      tmp = XCL_NGeo_Shared(1:3 , :    , :    , 0   ,iElem )
    CASE(ZETA_PLUS)
      tmp = XCL_NGeo_Shared(1:3 , :    , :    , NGeo,iElem )
    END SELECT
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,tmp,tmp2)

    ! get global SideID of local side
    SideID = GetGlobalNonUniqueSideID(iElem,iLocSide)

    DO q = 0,NGeo; DO p = 0,NGeo
      ! turn into right hand system of side
      pq = CGNS_SideToVol2(NGeo,p,q,iLocSide,3)
      BezierControlPoints3D(1:3,p,q,SideID) = tmp2(1:3,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! ilocSide=1,6
END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3D_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! Calculate elevated BezierControlPoints
IF (BezierElevation.GT.0) THEN
  DO iSide=firstSide,LastSide
    ! Ignore small mortar sides attached to big mortar sides
    IF (SideInfo_Shared(SIDE_LOCALID,iSide).LT.1 .OR. SideInfo_Shared(SIDE_LOCALID,iSide).GT.6) CYCLE
    ! Indices in shared arrays are shifted by 1
    CALL GetBezierControlPoints3DElevated( NGeo,NGeoElevated                                                       &
                                         , BezierControlPoints3D        (1:3,0:NGeo        ,0:NGeo        ,iSide)  &
                                         , BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide))
  END DO

#if USE_MPI
CALL MPI_WIN_SYNC(BezierControlPoints3DElevated_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/
END IF

END SUBROUTINE CalcBezierControlPoints


SUBROUTINE InitElemVolumes()
!===================================================================================================================================
! Calculate Element volumes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars       ,ONLY: wGP
USE MOD_Mesh_Vars                ,ONLY: nElems
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Particle_Mesh_Vars       ,ONLY: LocalVolume,MeshVolume
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemVolume_Shared
#if USE_MPI
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Mesh_Vars       ,ONLY: nComputeNodeElems,offsetComputeNodeElem
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemVolume_Shared_Win
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k
INTEGER            :: iElem
REAL               :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
INTEGER            :: offsetElemCNProc
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
!===================================================================================================================================
#if USE_MPI
! J_N is only build for local DG elements. Therefore, array is only filled for elements on the same compute node
offsetElemCNProc = offsetElem - offsetComputeNodeElem
#else
offsetElemCNProc = 0
#endif  /*USE_MPI*/

#if USE_MPI
MPISharedSize = INT(nComputeNodeElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeElems/),ElemVolume_Shared_Win,ElemVolume_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemVolume_Shared_Win,IERROR)
!CALL Allocate_Shared(MPISharedSize,(/nComputeNodeElems/),ElemCharLength_Shared_Win,ElemCharLength_Shared)
!CALL MPI_WIN_LOCK_ALL(0,ElemCharLength_Shared_Win,IERROR)

! Only root nullifies
IF (myComputeNodeRank.EQ.0) THEN
  ElemVolume_Shared(:)          = 0.
!  ElemCharLength_Shared(:)      = 0.
END IF
CALL MPI_WIN_SYNC(ElemVolume_Shared_Win,IERROR)
!CALL MPI_WIN_SYNC(ElemCharLength_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#else
ALLOCATE(ElemVolume_Shared(nElems))
!ALLOCATE(ElemCharLength_Shared(nElems))
#endif  /*USE_MPI*/

! Only root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif  /*USE_MPI*/
  ElemVolume_Shared(:)     = 0.
!  ElemCharLength_Shared(:) = 0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemVolume_Shared_Win,IERROR)
!CALL MPI_WIN_SYNC(ElemCharLength_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif  /*USE_MPI*/

! Calculate element volumes and characteristic lengths
DO iElem = 1,nElems
  !--- Calculate and save volume of element iElem
  J_N(1,0:PP_N,0:PP_N,0:PP_N) = 1./sJ(:,:,:,iElem,0)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    ElemVolume_Shared(iElem+offsetElemCNProc) = ElemVolume_Shared(iElem+offsetElemCNProc) + wGP(i)*wGP(j)*wGP(k)*J_N(1,i,j,k)
  END DO; END DO; END DO
!  !---- Calculate characteristic cell length: V^(1/3)
!  ElemCharLength_Shared(iElem+offsetElemCNProc) = ElemVolume_Shared(iElem+offsetElemCNProc)**(1./3.)
END DO

#if USE_MPI
CALL MPI_WIN_SYNC(ElemVolume_Shared_Win,IERROR)
!CALL MPI_WIN_SYNC(ElemCharLength_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

LocalVolume = SUM(ElemVolume_Shared(offsetElemCNProc+1:offsetElemCNProc+nElems))
MeshVolume  = SUM(ElemVolume_Shared(:))

SWRITE(UNIT_StdOut,'(A,E18.8)') ' | Total MESH Volume: ', MeshVolume
END SUBROUTINE InitElemVolumes


SUBROUTINE InitParticleGeometry()
!===================================================================================================================================
! Subroutine for particle geometry initialization in case of TriaTracking.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                ,ONLY: NGeo
USE MOD_Particle_Mesh_Tools      ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: ConcaveElemSide_Shared,ElemNodeID_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemSideNodeID_Shared,ElemMidPoint_Shared
#if USE_MPI
USE MOD_Particle_Mesh_Vars       ,ONLY: ConcaveElemSide_Shared_Win,ElemNodeID_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemSideNodeID_Shared_Win,ElemMidPoint_Shared_Win
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars                ,ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,FirstElem,LastElem,GlobalElemID
INTEGER            :: GlobalSideID,nlocSides,localSideID,iLocSide
INTEGER            :: iNode
INTEGER            :: nStart, NodeNum
INTEGER            :: NodeMap(4,6)
REAL               :: A(3,3),detcon
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

! Corner node switch to order HOPR coordinates in CGNS format
ASSOCIATE(CNS => CornerNodeIDswitch )
! CGNS Mapping
NodeMap(:,1)=(/CNS(1),CNS(4),CNS(3),CNS(2)/)
NodeMap(:,2)=(/CNS(1),CNS(2),CNS(6),CNS(5)/)
NodeMap(:,3)=(/CNS(2),CNS(3),CNS(7),CNS(6)/)
NodeMap(:,4)=(/CNS(3),CNS(4),CNS(8),CNS(7)/)
NodeMap(:,5)=(/CNS(1),CNS(5),CNS(8),CNS(4)/)
NodeMap(:,6)=(/CNS(5),CNS(6),CNS(7),CNS(8)/)

#if USE_MPI
MPISharedSize = INT(6*nComputeNodeTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),ConcaveElemSide_Shared_Win,ConcaveElemSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,ConcaveElemSide_Shared_Win,IERROR)
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
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

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif
  ElemNodeID_Shared      = 0
  ElemSideNodeID_Shared  = 0
  ConcaveElemSide_Shared = .FALSE.
#if USE_MPI
END IF
CALL MPI_WIN_SYNC(ElemNodeID_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemSideNodeID_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ConcaveElemSide_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

! iElem is CNElemID
DO iElem = firstElem,lastElem
  GlobalElemID = GetGlobalElemID(iElem)

  DO iNode = 1,8
    ElemNodeID_Shared(iNode,iElem) = ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID) + CNS(iNode)
  END DO

  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
  DO iLocSide = 1,nlocSides
    ! Get global SideID
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide

    IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    ! Find start of CGNS mapping from flip
    nStart = MERGE(0,MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,GlobalSideID),10)-1),SideInfo_Shared(SIDE_ID,GlobalSideID).GT.0)
    ! Shared memory array starts at 1, but NodeID at 0
    ElemSideNodeID_Shared(1:4,localSideID,iElem) = (/ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart  ,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+1,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+2,4)+1,localSideID)-1, &
                                                     ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+NodeMap(MOD(nStart+3,4)+1,localSideID)-1/)
  END DO
END DO
END ASSOCIATE
#if USE_MPI
CALL MPI_WIN_SYNC(ElemNodeID_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemSideNodeID_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

!--- Save whether Side is concave or convex
DO iElem = firstElem,lastElem
  ! iElem is CNElemID
  GlobalElemID = GetGlobalElemID(iElem)

  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
  DO iLocSide = 1,nlocSides
    !--- Check whether the bilinear side is concave
    !--- Node Number 4 and triangle 1-2-3
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
    IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
    localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
    DO NodeNum = 1,3               ! for all 3 nodes of triangle
       A(:,NodeNum) = NodeCoords_Shared(:,ElemSideNodeID_Shared(NodeNum,localSideID,iElem)+1) &
                    - NodeCoords_Shared(:,ElemSideNodeID_Shared(4      ,localSideID,iElem)+1)
    END DO
    !--- concave if detcon < 0:
    detcon = ((A(2,1) * A(3,2) - A(3,1) * A(2,2)) * A(1,3) +     &
              (A(3,1) * A(1,2) - A(1,1) * A(3,2)) * A(2,3) +     &
              (A(1,1) * A(2,2) - A(2,1) * A(1,2)) * A(3,3))
    !--- arbitrary choice if detcon exactly zero, define the one with lower SideID as concave
    IF (detcon.LT.0) THEN
      ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
    ELSE IF (detcon.EQ.0.0) THEN
      IF (GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) ConcaveElemSide_Shared(localSideID,iElem) = .TRUE.
    END IF
  END DO
END DO

DO iElem = firstElem,lastElem
  ! iElem is CNElemID
  GlobalElemID = GetGlobalElemID(iElem)

  ElemMidPoint_Shared(:,iElem) = 0.
  DO iNode = 1,8
    ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,GlobalElemID)+iNode)
  END DO
  ElemMidPoint_Shared(1:3,iElem) = ElemMidPoint_Shared(1:3,iElem) / 8.
END DO

#if USE_MPI
CALL MPI_WIN_SYNC(ConcaveElemSide_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemMidPoint_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif

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
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared,ConcaveElemSide_Shared,ElemSideNodeID_Shared
USE MOD_Particle_Mesh_Vars        ,ONLY: WeirdElems
USE MOD_Particle_Mesh_Tools       ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: nComputeNodeTotalElems,nComputeNodeProcessors,myComputeNodeRank
#else
USE MOD_Mesh_Vars                 ,ONLY: nElems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, iLocSide, kLocSide, iNode
INTEGER,ALLOCATABLE :: WeirdElemNbrs(:)
REAL                :: vec(1:3), Node(1:3,1:4),det(1:3)
LOGICAL             :: WEIRD, TRICHECK, TRIABSCHECK
INTEGER             :: firstElem,lastElem
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' CHECKING FOR WEIRD ELEMENTS...'

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif

ALLOCATE(WeirdElemNbrs(1:lastElem-firstElem+1))

WeirdElems = 0

! go through all CN elements
DO iElem = firstElem,lastElem
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
    WeirdElems = WeirdElems + 1
    WeirdElemNbrs(WeirdElems) = GetGlobalElemID(iElem)
  END IF
END DO

IF(WeirdElems.GT.0) THEN
  IPWRITE(UNIT_stdOut,*)' FOUND', WeirdElems, 'ELEMENTS!'
  IPWRITE(UNIT_stdOut,*)' WEIRD ELEM NUMBERS:'
  DO iElem = 1,WeirdElems
    IPWRITE(UNIT_stdOut,*) WeirdElemNbrs(iElem)
  END DO
END IF

SWRITE(UNIT_stdOut,'(A)')' CHECKING FOR WEIRD ELEMENTS DONE!'

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
USE MOD_Mesh_Vars                 ,ONLY: nElems
USE MOD_Particle_Mesh_Vars        ,ONLY: NbrOfRegions,RegionBounds,GEO
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemBaryNGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
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
      GEO%ElemToRegion(iElem) = iRegions
    ELSE
      CALL ABORT(__STAMP__,'Defined regions are overlapping')
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
USE MOD_Basis                     ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars                 ,ONLY: NGeo
USE MOD_Particle_Basis            ,ONLY: DeCasteljauInterpolation
USE MOD_Particle_Surfaces_Vars    ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Vars        ,ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars        ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis,ElemRadiusNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Tools       ,ONLY: GetGlobalElemID,GetGlobalNonUniqueSideID
#if USE_MPI
USE MOD_Particle_MPI_Shared       ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars        ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemRadiusNGEO_Shared,ElemRadiusNGeo_Shared_Win
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemRadius2NGeo_Shared,ElemRadius2NGeo_Shared_Win
USE MOD_Particle_Mesh_Vars        ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars        ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
#else
USE MOD_Particle_Mesh_Vars        ,ONLY: XCL_NGeo
USE MOD_Mesh_Vars                 ,ONLY: nELems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,SideID
INTEGER                        :: i,j,k,ilocSide
INTEGER                        :: iDir
REAL                           :: Xi(3,6),xPos(3),Radius
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem,lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
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

ASSOCIATE(XCL_NGeo     => XCL_NGeo_Shared)

#else
ALLOCATE(ElemRadiusNGeo(          nElems) &
        ,ElemRadius2NGeo(         nElems) &
        ,XiEtaZetaBasis(1:3,1:6,1:nElems) &
        ,slenXiEtaZetaBasis(1:6,1:nElems))
#endif /*USE_MPI*/

#if USE_MPI
IF(myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  ElemRadiusNGeo =0.
  ElemRadius2NGeo=0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemRadiusNGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

! iElem is CN elem
DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

    xPos = 0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos = xPos + XCL_NGeo(:,i,j,k,ElemID)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO

    XiEtaZetaBasis(1:3,iDir,iElem) = xPos
    ! compute vector from each barycenter to sidecenter
    XiEtaZetaBasis(  :,iDir,iElem) = XiEtaZetaBasis(:,iDir,iElem)-ElemBaryNGeo(:,iElem)
    ! compute length: The root is omitted here due to optimization
    slenXiEtaZetaBasis(iDir,iElem) = 1.0/DOT_PRODUCT(XiEtaZetaBasis(:,iDir,iElem),XiEtaZetaBasis(:,iDir,iElem))
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
#endif /*USE_MPI*/

END SUBROUTINE BuildElementBasisAndRadius


SUBROUTINE BuildElementRadiusTria()
!================================================================================================================================
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Globals          ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemBaryNGeo,ElemRadius2NGeo
USE MOD_Particle_Mesh_Tools       ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemBaryNGeo_Shared,ElemRadius2NGeo_Shared
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemBaryNGeo_Shared_Win,ElemRadius2NGeo_Shared_Win
USE MOD_Particle_MPI_Shared       ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_Particle_MPI_Shared_Vars  ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars                 ,ONLY: nELems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,iNode
REAL                           :: xPos(3),Radius
INTEGER                        :: firstElem, lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
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
        ,ElemRadius2NGeo( nElems))
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemRadius2NGeo = 0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  Radius = 0.
  xPos   = 0.

  DO iNode = 1,8
    xPos = xPos + NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode)
  END DO
    ElemBaryNGeo(:,iElem) = xPos/8.
  DO iNode = 1,8
    xPos   = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode) - ElemBaryNGeo(:,iElem)
    Radius = MAX(Radius,VECNORM(xPos))
  END DO
  ElemRadius2NGeo(iElem) = Radius*Radius
END DO ! iElem

#if USE_MPI
CALL MPI_WIN_SYNC(ElemRadius2NGeo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildElementRadiusTria


SUBROUTINE BuildElemTypeAndBasisTria()
!===================================================================================================================================
!> Dummy routine to fill the ElemCurved array with TriaTracking
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis                   ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: XiEtaZetaBasis_Shared,XiEtaZetaBasis_Shared_Win
USE MOD_Particle_Mesh_Vars      ,ONLY: slenXiEtaZetaBasis_Shared,slenXiEtaZetaBasis_Shared_Win
USE MOD_Particle_MPI_Shared     ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPI_COMM_SHARED
#else
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo
USE MOD_Mesh_Vars               ,ONLY: nElems
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,iDir
INTEGER                        :: i,j,k
REAL                           :: xPos(3)
REAL                           :: Xi(3,6),Lag(1:3,0:NGeo)
INTEGER                        :: firstElem, lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /* USE_MPI */
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying side types and whether elements are curved...'

! elements
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemCurved_Shared_Win,ElemCurved_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemCurved_Shared_Win,IERROR)
MPISharedSize = INT((3*6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,6,nComputeNodeTotalElems/),XiEtaZetaBasis_Shared_Win,XiEtaZetaBasis_Shared)
CALL MPI_WIN_LOCK_ALL(0,XiEtaZetaBasis_Shared_Win,IERROR)
MPISharedSize = INT((6*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nComputeNodeTotalElems/),slenXiEtaZetaBasis_Shared_Win,slenXiEtaZetaBasis_Shared)
CALL MPI_WIN_LOCK_ALL(0,slenXiEtaZetaBasis_Shared_Win,IERROR)
ElemCurved         => ElemCurved_Shared
XiEtaZetaBasis     => XiEtaZetaBasis_Shared
slenXiEtaZetaBasis => slenXiEtaZetaBasis_Shared

ASSOCIATE(XCL_NGeo     => XCL_NGeo_Shared)

#else
ALLOCATE(ElemCurved(            1:nElems) &
        ,XiEtaZetaBasis(1:3,1:6,1:nElems) &
        ,slenXiEtaZetaBasis(1:6,1:nElems))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  ElemCurved   = .FALSE.
#if USE_MPI
END IF
#endif /*USE_MPI*/

#if USE_MPI
firstElem = INT(REAL( myComputeNodeRank*   nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

Xi(:,1) = (/ 1.0 , 0.0  ,  0.0/) ! xi plus
Xi(:,2) = (/ 0.0 , 1.0  ,  0.0/) ! eta plus
Xi(:,3) = (/ 0.0 , 0.0  ,  1.0/) ! zeta plus
Xi(:,4) = (/-1.0 , 0.0  ,  0.0/) ! xi minus
Xi(:,5) = (/ 0.0 , -1.0 ,  0.0/) ! eta minus
Xi(:,6) = (/ 0.0 , 0.0  , -1.0/) ! zeta minus

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  ! get point on each side
  DO iDir = 1, 6
    CALL LagrangeInterpolationPolys(Xi(1,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3,iDir),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

    xPos = 0.
    DO k = 0,NGeo; DO j = 0,NGeo; DO i = 0,NGeo
      xPos = xPos+XCL_NGeo(:,i,j,k,ElemID)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO

    XiEtaZetaBasis(1:3,iDir,iElem) = xPos
    ! compute vector from each barycenter to sidecenter
    XiEtaZetaBasis(:,iDir,iElem)   = XiEtaZetaBasis(:,iDir,iElem)-ElemBaryNGeo(:,iElem)
    ! compute length: The root is omitted here due to optimization
    slenXiEtaZetaBasis(iDir,iElem) = 1.0/DOT_PRODUCT(XiEtaZetaBasis(:,iDir,iElem),XiEtaZetaBasis(:,iDir,iElem))
  END DO ! iDir = 1, 6
END DO

#if USE_MPI
END ASSOCIATE
CALL MPI_WIN_SYNC(ElemCurved_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(XiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(slenXiEtaZetaBasis_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE BuildElemTypeAndBasisTria


SUBROUTINE PointsEqual(N,Points1,Points2,IsNotEqual)
!===================================================================================================================================
! compute the distance between two data sets
!===================================================================================================================================
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
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
USE MOD_Basis                   ,ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo_Shared,ElemBaryNGeo_Shared_Win
USE MOD_Particle_MPI_Shared     ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars,ONLY: MPI_COMM_SHARED
#else
USE MOD_Mesh_Vars               ,ONLY: nElems
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeElems
USE MOD_Particle_Mesh_Vars      ,ONLY: XCL_NGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,ElemID,i,j,k
REAL                           :: XPos(3),buf
REAL                           :: Lag(1:3,0:NGeo)
INTEGER                        :: firstElem,lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!================================================================================================================================

#if USE_MPI
MPISharedSize = INT((3*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalElems/),ElemBaryNGeo_Shared_Win,ElemBaryNGeo_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemBaryNGeo_Shared_Win,IERROR)
ElemBaryNGeo => ElemBaryNGeo_Shared

ASSOCIATE(XCL_NGeo => XCL_NGeo_Shared)

! Set ranges
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(ElemBaryNGeo(1:3,nComputeNodeElems))
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  ElemBaryNGeo = 0.
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(ElemBaryNGeo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

! evaluate the polynomial at origin: Xi=(/0.0,0.0,0.0/)
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(0.0,NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  xPos = 0.

  DO k = 0,NGeo
    DO j = 0,NGeo
      buf = Lag(2,j)*Lag(3,k)
      DO i = 0,NGeo
        xPos = xPos + XCL_NGeo(1:3,i,j,k,ElemID)*Lag(1,i)*buf
      END DO ! i = 0,NGeo
    END DO ! j = 0,NGeo
  END DO ! k = 0,NGeo

  ElemBaryNGeo(:,iElem)=xPos
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
USE MOD_ChangeBasis              ,ONLY: changeBasis3D
USE MOD_Mesh_Vars                ,ONLY: NGeo,nElems
USE MOD_Particle_Globals         ,ONLY: ALMOSTZERO,CROSSNORM,UNITVECTOR
USE MOD_Particle_Mesh_Vars       ,ONLY: nNonUniqueGlobalSides
USE MOD_Particle_Mesh_Vars       ,ONLY: Vdm_CLNGeo1_CLNGeo,Vdm_CLNGeo1_CLNGeo
USE MOD_Particle_Mesh_Vars       ,ONLY: XCL_NGeo_Shared,ElemBaryNGeo
USE MOD_Particle_Mesh_Vars       ,ONLY: SideInfo_Shared,ElemCurved
USE MOD_Particle_Mesh_Tools      ,ONLY: GetGlobalElemID,GetCNElemID,GetCNSideID,GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars   ,ONLY: BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierControlPoints3D,SideType,SideNormVec,SideDistance
#if USE_MPI
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemCurved_Shared,ElemCurved_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: SideDistance_Shared,SideDistance_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: SideType_Shared,SideType_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: SideNormVec_Shared,SideNormVec_Shared_Win
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeTotalElems,nComputeNodeTotalSides
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
#else
USE MOD_Particle_Mesh_Vars       ,ONLY: nComputeNodeElems
#endif /* USE_MPI */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: ilocSide,SideID,CNSideID,flip
INTEGER                                  :: iElem,firstElem,lastElem,ElemID
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
!#if USE_MPI
!INTEGER                                  :: nDummy
!#endif /* USE_MPI */
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying side types and whether elements are curved...'

! elements
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemCurved_Shared_Win,ElemCurved_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemCurved_Shared_Win,IERROR)
ElemCurved => ElemCurved_Shared
#else
ALLOCATE(ElemCurved(1:nComputeNodeElems))
#endif /*USE_MPI*/

! sides
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalSides/),  SideType_Shared_Win,    SideType_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideType_Shared_Win,IERROR)
SideType => SideType_Shared
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalSides/),  SideDistance_Shared_Win,SideDistance_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideDistance_Shared_Win,IERROR)
SideDistance => SideDistance_Shared
MPISharedSize = INT((3*nComputeNodeTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nComputeNodeTotalSides/),SideNormVec_Shared_Win, SideNormVec_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideNormVec_Shared_Win,IERROR)
SideNormVec => SideNormVec_Shared
#else
ALLOCATE(SideType(       nNonUniqueGlobalSides))
ALLOCATE(SideDistance(   nNonUniqueGlobalSides))
ALLOCATE(SideNormVec(1:3,nNonUniqueGlobalSides))
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

ALLOCATE(SideIsDone(nNonUniqueGlobalSides))
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

DO iElem = firstElem,lastElem
  ElemID = GetGlobalElemID(iElem)
  XCL_NGeoLoc = XCL_NGeo_Shared(1:3,0:NGeo,0:NGeo,0:NGeo,ElemID)
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
    SideID   = GetGlobalNonUniqueSideID(ElemID,iLocSide)
    CNSideID = GetCNSideID(SideID)
    flip     = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

    IF(.NOT.ElemCurved(iElem))THEN
      BezierControlPoints_loc(1:3,0:NGeo,0:NGeo) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
      ! linear element
      IF(BoundingBoxIsEmpty(CNSideID))THEN
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,CNSideID) = CROSSNORM(v1,v2)
        v1=0.25*(BezierControlPoints_loc(:,0,0      )     &
                +BezierControlPoints_loc(:,NGeo,0   )  &
                +BezierControlPoints_loc(:,0,NGeo   )  &
                +BezierControlPoints_loc(:,NGeo,NGeo))
        ! check if normal vector points outwards
        v2=v1-ElemBaryNGeo(:,iElem)
        IF(flip.EQ.0)THEN
          IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).LT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
        ELSE
          IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).GT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
        END IF
        SideDistance(CNSideID)=DOT_PRODUCT(v1,SideNormVec(:,CNSideID))
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
          SideType(CNSideID)=PLANAR_RECT
        ELSE
          SideType(CNSideID)=PLANAR_NONRECT
        END IF
      ELSE
        v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
            -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
            +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
        SideNormVec(:,CNSideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
        SideType(CNSideID)=BILINEAR
      END IF
    ELSE
      BezierControlPoints_loc(1:3,0:NGeo,0:NGeo) = BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)
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
        IF(BoundingBoxIsEmpty(CNSideID))THEN
          SideType(CNSideID)=PLANAR_CURVED
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,CNSideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0   )  &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).LT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).GT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          END IF
          SideDistance(CNSideID)=DOT_PRODUCT(v1,SideNormVec(:,CNSideID))
        ELSE
          SideType(CNSideID)=CURVED
        END IF
      ELSE
        IF (BoundingBoxIsEmpty(CNSideID)) THEN
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )

          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo) )
          SideNormVec(:,CNSideID) = CROSSNORM(v1,v2)
          v1=0.25*(BezierControlPoints_loc(:,0,0)     &
                  +BezierControlPoints_loc(:,NGeo,0)  &
                  +BezierControlPoints_loc(:,0,NGeo)  &
                  +BezierControlPoints_loc(:,NGeo,NGeo))
          ! check if normal vector points outwards
          v2=v1-ElemBaryNGeo(:,iElem)
          IF(flip.EQ.0)THEN
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).LT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          ELSE
            IF(DOT_PRODUCT(v2,SideNormVec(:,CNSideID)).GT.0) SideNormVec(:,CNSideID)=-SideNormVec(:,CNSideID)
          END IF
          SideDistance(CNSideID)=DOT_PRODUCT(v1,SideNormVec(:,CNSideID))
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
            SideType(CNSideID)=PLANAR_RECT
          ELSE
            SideType(CNSideID)=PLANAR_NONRECT
          END IF
        ELSE
          v1=(-BezierControlPoints_loc(:,0,0   )+BezierControlPoints_loc(:,NGeo,0   )   &
              -BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          v2=(-BezierControlPoints_loc(:,0,0   )-BezierControlPoints_loc(:,NGeo,0   )   &
              +BezierControlPoints_loc(:,0,NGeo)+BezierControlPoints_loc(:,NGeo,NGeo))
          SideNormVec(:,CNSideID) = CROSSNORM(v1,v2) !non-oriented, averaged normal vector based on all four edges
          SideType(CNSideID)=BILINEAR
        END IF
      END IF
    END IF
    SideIsDone(SideID)=.TRUE.
  END DO ! ilocSide=1,6
END DO ! iElem=1,nTotalElems

DEALLOCATE(SideIsDone)

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

#if USE_MPI
firstElem = offsetElem + 1
lastElem  = offsetElem + nElems
#else
firstElem = 1
lastElem  = nElems
#endif


DO iElem = firstElem,lastElem
  ElemID = GetCNElemID(iElem)
  IF (ElemCurved(ElemID)) THEN
    nCurvedElems = nCurvedElems+1
  ELSE
    nLinearElems = nLinearElems+1
  END IF

  DO ilocSide = 1,6
    ! ignore small mortar sides attached to big mortar sides
    SideID   = GetGlobalNonUniqueSideID(iElem,ilocSide)
    CNSideID = GetCNSideID(SideID)
    SELECT CASE(SideType(CNSideID))
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
CALL MPI_REDUCE(nPlanarRectangular   ,nPlanarRectangularTot   ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
CALL MPI_REDUCE(nPlanarNonRectangular,nPlanarNonRectangularTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
CALL MPI_REDUCE(nBilinear            ,nBilinearTot            ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
CALL MPI_REDUCE(nPlanarCurved        ,nPlanarCurvedTot        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
CALL MPI_REDUCE(nCurved              ,nCurvedTot              ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
CALL MPI_REDUCE(nLinearElems         ,nLinearElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
CALL MPI_REDUCE(nCurvedElems         ,nCurvedElemsTot         ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,IERROR)
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
! This only works if the full mesh is built on one node!
!#if USE_MPI
!IF (myComputeNodeRank.EQ.0) THEN
!#endif /* USE_MPI */
!  IF ((nComputeNodeTotalElems-nCurvedElemsTot).NE.nLinearElemsTot) &
!    CALL ABORT(__STAMP__, 'Error in particle mesh: lost elements while trying to dermine if elements are curved')
!#if USE_MPI
!END IF
!#endif /* USE_MPI */

SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-rectangular     faces: ', nPlanarRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-non-rectangular faces: ', nPlanarNonRectangulartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of bi-linear              faces: ', nBilineartot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of planar-curved          faces: ', nPlanarCurvedtot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 faces: ', nCurvedtot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of (bi-)linear            elems: ', nLinearElemsTot
SWRITE(UNIT_StdOut,'(A,I8)') ' | Number of curved                 elems: ', nCurvedElemsTot

END SUBROUTINE IdentifyElemAndSideType


SUBROUTINE GetBCSidesAndOrgin()
!===================================================================================================================================
! Globally identifies all BC sides and build side origin and radius
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars                ,ONLY: NGeo
USE MOD_Particle_Basis           ,ONLY: DeCasteljauInterpolation
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars       ,ONLY: SideInfo_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars       ,ONLY: nNonUniqueGlobalSides,nUniqueBCSides
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierControlPoints3D
#if USE_MPI
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars       ,ONLY: BCSide2SideID_Shared,SideID2BCSide_Shared,BCSideMetrics_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: BCSide2SideID_Shared_Win,SideID2BCSide_Shared_Win,BCSideMetrics_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: p,q
INTEGER                        :: iSide,firstSide,lastSide,BCSideID
INTEGER                        :: nUniqueBCSidesProc,offsetUniqueBCSidesProc
REAL                           :: origin(1:3),xi(1:2),radius,radiusMax,vec(1:3)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nNonUniqueGlobalSides
#endif /*USE_MPI*/

! Count number of BC sides in range
nUniqueBCSidesProc = 0
DO iSide = firstSide,lastSide
  ! ignore elements not on the compute node
  IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

  ! ignore inner and virtual (mortar) sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  nUniqueBCSidesProc = nUniqueBCSidesProc + 1
END DO

! Find global number of BC sides and write side <=> BCSide mapping into shared array
#if USE_MPI
sendbuf = nUniqueBCSidesProc
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetUniqueBCSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetUniqueBCSidesProc + nUniqueBCSidesProc
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nUniqueBCSides = sendbuf

MPISharedSize = INT((nUniqueBCSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nUniqueBCSides/),BCSide2SideID_Shared_Win,BCSide2SideID_Shared)
CALL MPI_WIN_LOCK_ALL(0,BCSide2SideID_Shared_Win,IERROR)
MPISharedSize = INT((nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalSides/),SideID2BCSide_Shared_Win,SideID2BCSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideID2BCSide_Shared_Win,IERROR)
BCSide2SideID => BCSide2SideID_Shared
SideID2BCSide => SideID2BCSide_Shared

! Also allocate array to hold BC Side metrics
MPISharedSize = INT((4*nUniqueBCSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/4,nUniqueBCSides/),BCSideMetrics_Shared_Win,BCSideMetrics_Shared)
CALL MPI_WIN_LOCK_ALL(0,BCSideMetrics_Shared_Win,IERROR)
BCSideMetrics => BCSideMetrics_Shared
#else
offsetUniqueBCSidesProc = 0
nUniqueBCSides = nUniqueBCSidesProc

ALLOCATE(BCSide2SideID(    1:nUniqueBCSides)        &
        ,SideID2BCSide(    1:nNonUniqueGlobalSides))
! Also allocate array to hold BC Side metrics
ALLOCATE(BCSideMetrics(1:4,1:nUniqueBCSides))
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  BCSide2SideID = -1
  SideID2BCSide = -1
  BCSideMetrics = -1
#if USE_MPI
END IF

CALL MPI_WIN_SYNC(BCSide2SideID_Shared_Win,iError)
CALL MPI_WIN_SYNC(SideID2BCSide_Shared_Win,iError)
CALL MPI_WIN_SYNC(BCSideMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

nUniqueBCSidesProc = 0
DO iSide = firstSide,lastSide
  ! ignore elements not on the compute node
  IF (GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.-1) CYCLE

  ! ignore inner and virtual (mortar) sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

  nUniqueBCSidesProc = nUniqueBCSidesProc + 1
  BCSideID           = offsetUniqueBCSidesProc + nUniqueBCSidesProc
  BCSide2SideID(BCSideID) = iSide
  SideID2BCSide(iSide)    = BCSideID

  ! calculate origin, radius for all BC sides
  !> build side origin
  xi     = 0.
  ! TODO: BezierControlPoints are allocated with global side ID, so this SHOULD work. Breaks if we reduce the halo region
  CALL DeCasteljauInterpolation(NGeo,xi,iSide,origin)
  BCSideMetrics(1:3,BCSideID) = origin(1:3)

  !> build side radius
  radiusMax = 0.
  DO q = 0,NGeo
    DO p = 0,NGeo
      vec(1:3) = BezierControlPoints3D(:,p,q,iSide) - origin
      radius   = DOT_PRODUCT(Vec,Vec)
      radiusMax= MAX(radiusMax,radius)
    END DO
  END DO
  BCSideMetrics(4,BCSideID) = SQRT(RadiusMax)
END DO

#if USE_MPI
CALL MPI_WIN_SYNC(BCSide2SideID_Shared_Win,iError)
CALL MPI_WIN_SYNC(SideID2BCSide_Shared_Win,iError)
CALL MPI_WIN_SYNC(BCSideMetrics_Shared_Win,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

END SUBROUTINE GetBCSidesAndOrgin


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
USE MOD_Mesh_Vars                ,ONLY: NGeo
USE MOD_Particle_Globals         ,ONLY: ALMOSTZERO,VECNORM
USE MOD_Particle_Mesh_Vars       ,ONLY: GEO
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemInfo_Shared,SideInfo_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemToBCSides,SideBCMetrics
USE MOD_Particle_Mesh_Vars       ,ONLY: BCSide2SideID,SideID2BCSide,BCSideMetrics
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemBaryNGeo,ElemRadiusNGeo
USE MOD_Particle_Mesh_Vars       ,ONLY: nNonUniqueGlobalSides,nUniqueBCSides
USE MOD_Particle_Mesh_Tools      ,ONLY: GetGlobalElemID,GetCNElemID,GetCNElemID,GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierControlPoints3D
USE MOD_Particle_Utils           ,ONLY: InsertionSort
#if USE_MPI
USE MOD_CalcTimeStep             ,ONLY: CalcTimeStep
USE MOD_Particle_Mesh_Vars       ,ONLY: nNonUniqueGlobalSides
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemToBCSides_Shared,SideBCMetrics_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemToBCSides_Shared_Win,SideBCMetrics_Shared_Win
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_MPI_Vars        ,ONLY: halo_eps,halo_eps_velo
USE MOD_Particle_Timedisc_Vars   ,ONLY: ManualTimeStep
USE MOD_TimeDisc_Vars            ,ONLY: nRKStages,RKc
#else
USE MOD_Mesh_Vars                ,ONLY: nElems
USE MOD_Particle_Mesh_Vars       ,ONLY: nComputeNodeElems
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
INTEGER                        :: iBCSide,BCElemID,BCSideID
INTEGER                        :: CNElemID,BCCNElemID
REAL                           :: dX,dY,dZ
REAL                           :: origin(1:3),vec(1:3)
REAL                           :: BC_halo_eps
LOGICAL                        :: fullMesh
REAL,ALLOCATABLE               :: tmpSideBCMetrics(:,:)
REAL,ALLOCATABLE               :: tmpSideBCDistance(:)
INTEGER,ALLOCATABLE            :: intSideBCMetrics(:)
#if USE_MPI
REAL                           :: BC_halo_eps_velo,BC_halo_diag,deltaT
INTEGER                        :: iStage
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: errType
INTEGER                        :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Identifying BC sides and calculating side metrics...'

! elements
#if USE_MPI
MPISharedSize = INT((2*nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,nComputeNodeTotalElems/),ElemToBCSides_Shared_Win,ElemToBCSides_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToBCSides_Shared_Win,IERROR)
ElemToBCSides => ElemToBCSides_Shared
#else
ALLOCATE(ElemToBCSides(1:2,1:nComputeNodeElems))
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

! if running on one node, halo_eps is meaningless. Get a representative BC_halo_eps for BC side identification
fullMesh = .FALSE.
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

  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  BC_halo_diag = VECNORM(vec)

  ! compare halo_eps against global diagonal and reduce if necessary
  IF (.NOT.ALMOSTZERO(BC_halo_eps).AND.(BC_halo_diag.GE.BC_halo_eps)) THEN
    SWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to ',BC_halo_eps
  ELSEIF (.NOT.ALMOSTZERO(BC_halo_eps).AND.(BC_halo_diag.LT.BC_halo_eps)) THEN
    fullMesh = .TRUE.
    BC_halo_eps = BC_halo_diag
    SWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to global diag with ',BC_halo_eps
  ! halo_eps still at zero. Set it to global diagonal
  ELSE
    fullMesh = .TRUE.
    BC_halo_eps = BC_halo_diag
    SWRITE(UNIT_stdOUt,'(A,F11.3)') ' | No halo_eps given and could not be reconstructed. Using global diag with ',BC_halo_eps
  END IF
ELSE
  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  BC_halo_diag = VECNORM(vec)

  IF (BC_halo_diag.LE.halo_eps) fullMesh = .TRUE.
  BC_halo_eps = halo_eps
END IF

#else
! get distance of diagonal of mesh
vec(1)   = GEO%xmaxglob-GEO%xminglob
vec(2)   = GEO%ymaxglob-GEO%yminglob
vec(3)   = GEO%zmaxglob-GEO%zminglob
BC_halo_eps = DOT_PRODUCT(vec,vec)
fullMesh = .TRUE.

firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

nBCSidesProc      = 0
offsetBCSides     = 0

! for fullMesh, each element requires ALL BC faces
IF (fullMesh) THEN
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO

    DO iBCSide = 1,nUniqueBCSides
      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO ! iBCSide

    ! Write local mapping from Elem to BC sides. The number is already correct, the offset must be corrected later
    IF (nBCSidesElem.GT.0) THEN
      ElemToBCSides(ELEM_NBR_BCSIDES ,iElem) = nBCSidesElem
      ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides
    END IF
  END DO ! iElem

  offsetBCSides = nBCSidesProc

! .NOT. fullMesh
ELSE
  ! sum up all BC sides in range of BC_halo_eps
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)
        ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesElem = nBCSidesElem + 1
      nBCSidesProc = nBCSidesProc + 1
    END DO

      ! loop over all sides. Check distance from every local side to total sides.
      DO iBCSide = 1,nUniqueBCSides

        BCSideID   = BCSide2SideID(iBCSide)
        BCElemID   = SideInfo_Shared(SIDE_ELEMID,BCSideID)
        BCCNElemID = GetCNElemID(BCElemID)

        ! Ignore elements not on the compute node
        IF (BCCNElemID.EQ.-1) CYCLE

        ! Ignore the same element
        IF (BCElemID.EQ.iElem) CYCLE

        ! Check if barycenter of element is in range
        IF (VECNORM(ElemBaryNGeo(:,iElem) - ElemBaryNGeo(:,BCCNElemID)) &
          .GT. (BC_halo_eps + ElemRadiusNGeo(iElem) + ElemRadiusNGeo(BCCNElemID))) CYCLE

      ! loop over all local sides of the element. Use a named loop so the entire element can be cycled
Check1: DO ilocSide = 1,6
          SideID = GetGlobalNonUniqueSideID(ElemID,ilocSide)

          ! compare all nodes between local and BC side to check if within BC_halo_eps. Once one node pair is in range, flag the entire
          ! side and stop checking. First, get BezierControlPoints for local side. BezierControlPoints3D available for ALL sides in shared memory
          SELECT CASE(ilocSide)
            CASE(XI_MINUS,XI_PLUS)
              firstBezierPoint = 0
              lastBezierPoint  = NGeo
            CASE DEFAULT
              firstBezierPoint = 1
              lastBezierPoint  = NGeo-1
          END SELECT

            ! finally compare the node coords
          DO q = firstBezierPoint,lastBezierPoint
            DO p = firstBezierPoint,lastBezierPoint
!           ! get all nodes for BC side
!           NodeBCSide(:) = BezierControlPoints3D(:,p,q,BCSideID)
            ! finally compare the node coords
            DO s = firstBezierPoint,lastBezierPoint
              DO r = firstBezierPoint,lastBezierPoint
                  dX = ABS(BezierControlPoints3D(1,r,s,SideID)-BezierControlPoints3D(1,p,q,BCSideID))
                IF (dX.GT.BC_halo_eps) CYCLE
                  dY = ABS(BezierControlPoints3D(2,r,s,SideID)-BezierControlPoints3D(2,p,q,BCSideID))
                IF (dY.GT.BC_halo_eps) CYCLE
                  dZ = ABS(BezierControlPoints3D(3,r,s,SideID)-BezierControlPoints3D(3,p,q,BCSideID))
                IF (dZ.GT.BC_halo_eps) CYCLE

                IF (SQRT(dX*dX+dY*dY+dZ*dZ).LE.BC_halo_eps) THEN
                  nBCSidesElem = nBCSidesElem + 1
                  nBCSidesProc = nBCSidesProc + 1
                  EXIT Check1
                END IF
              END DO ! r
            END DO ! s
          END DO ! p
        END DO ! q
      END DO Check1 ! ilocSide

!      ! Compare distance of element center to side center while also taking the radii into account
!      !-- This approach gives almost triple the amount of BCSides of the method above!
!      IF (VECNORM(ElemBaryNGeo(:,iElem) - BCSideMetrics(1:3,iBCSide)) &
!        .LE. (BC_halo_eps + ElemRadiusNGeo(iElem) + BCSideMetrics(4,iBCSide))) THEN
!        nBCSidesElem = nBCSidesElem + 1
!        nBCSidesProc = nBCSidesProc + 1
!      END IF
    END DO ! iBCSide

    ! Write local mapping from Elem to BC sides. The number is already correct, the offset must be corrected later
    IF (nBCSidesElem.GT.0) THEN
      ElemToBCSides(ELEM_NBR_BCSIDES ,iElem) = nBCSidesElem
      ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) = offsetBCSides
    END IF

    offsetBCSides = nBCSidesProc
  END DO ! iElem
END IF ! fullMesh

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

ElemToBCSides(ELEM_FIRST_BCSIDE,firstElem:lastElem) = ElemToBCSides(ELEM_FIRST_BCSIDE,firstElem:lastElem) + offsetBCSidesProc
#else
offsetBCSidesProc   = 0
nComputeNodeBCSides = nBCSidesProc
#endif /*USE_MPI*/

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

CALL MPI_WIN_SYNC(ElemToBCSides_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(SideBCMetrics_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI*/

nBCSidesProc      = 0

! We did not know the number of BC sides before. Therefore, we need to do the check again and build the final mapping
! for fullMesh, each element requires ALL BC faces
IF (fullMesh) THEN
  DO iElem = firstElem,lastElem
    ElemID     = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(iElem)
    END DO

    DO iBCSide = 1,nUniqueBCSides
      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(BCSideID)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(BCElemID)
    END DO ! iBCSide
  END DO ! iElem

! .NOT. fullMesh
ELSE
  ! sum up all BC sides in range of BC_halo_eps
  DO iElem = firstElem,lastElem
    ElemID = GetGlobalElemID(iElem)
    nBCSidesElem  = 0

    ! check local side of an element
    DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID)

      ! ignore inner and virtual (mortar) sides
      IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

      nBCSidesProc = nBCSidesProc + 1
      SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(iSide)
      SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(ElemID)
    END DO

    ! loop over all sides. Check distance from every local side to total sides. Once a side has been flagged,
    ! it must not be counted again
    DO iBCSide = 1,nUniqueBCSides

      BCSideID = BCSide2SideID(iBCSide)
      BCElemID = SideInfo_Shared(SIDE_ELEMID,BCSideID)
      BCCNElemID = GetCNElemID(BCElemID)

      ! Ignore elements not on the compute node
      IF (BCCNElemID.EQ.-1) CYCLE

      ! Ignore the same element
      IF (BCElemID.EQ.iElem) CYCLE

      ! Check if barycenter of element is in range
      IF (VECNORM(ElemBaryNGeo(:,iElem) - ElemBaryNGeo(:,BCCNElemID)) &
        .GT. (BC_halo_eps + ElemRadiusNGeo(iElem) + ElemRadiusNGeo(BCCNElemID))) CYCLE

      ! loop over all local sides of the element. Use a named loop so the entire element can be cycled
Check2: DO ilocSide = 1,6
    SideID = GetGlobalNonUniqueSideID(ElemID,ilocSide)

        ! compare all nodes between local and BC side to check if within BC_halo_eps. Once one node pair is in range, flag the entire
        ! side and stop checking. First, get BezierControlPoints for local side. BezierControlPoints3D available for ALL sides in shared memory
      SELECT CASE(ilocSide)
        CASE(XI_MINUS,XI_PLUS)
          firstBezierPoint = 0
          lastBezierPoint  = NGeo
        CASE DEFAULT
          firstBezierPoint = 1
          lastBezierPoint  = NGeo-1
      END SELECT

        ! finally compare the node coords
        DO q = firstBezierPoint,lastBezierPoint
          DO p = firstBezierPoint,lastBezierPoint
!           ! get all nodes for BC side
!           NodeBCSide(:) = BezierControlPoints3D(:,p,q,BCSideID)
            ! finally compare the node coords
            DO s = firstBezierPoint,lastBezierPoint
              DO r = firstBezierPoint,lastBezierPoint
                dX = ABS(BezierControlPoints3D(1,r,s,SideID)-BezierControlPoints3D(1,p,q,BCSideID))
                IF (dX.GT.BC_halo_eps) CYCLE
                dY = ABS(BezierControlPoints3D(2,r,s,SideID)-BezierControlPoints3D(2,p,q,BCSideID))
                IF (dY.GT.BC_halo_eps) CYCLE
                dZ = ABS(BezierControlPoints3D(3,r,s,SideID)-BezierControlPoints3D(3,p,q,BCSideID))
                IF (dZ.GT.BC_halo_eps) CYCLE

                IF (SQRT(dX*dX+dY*dY+dZ*dZ).LE.BC_halo_eps) THEN
                  nBCSidesProc = nBCSidesProc + 1
                  SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(BCSideID)
                  SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(BCElemID)
                  EXIT Check2
                END IF
              END DO ! r
            END DO ! s
          END DO ! p
        END DO ! q
      END DO Check2 ! ilocSide

!      ! Compare distance of element center to side center while also taking the radii into account
!      !-- This approach gives almost triple the amount of BCSides of the method above!
!      IF (VECNORM(ElemBaryNGeo(:,iElem) - BCSideMetrics(1:3,iBCSide)) &
!        .LE. (BC_halo_eps + ElemRadiusNGeo(iElem) + BCSideMetrics(4,iBCSide))) THEN
!        nBCSidesProc = nBCSidesProc + 1
!        SideBCMetrics(BCSIDE_SIDEID,nBCSidesProc+offsetBCSidesProc) = REAL(BCSideID)
!        SideBCMetrics(BCSIDE_ELEMID,nBCSidesProc+offsetBCSidesProc) = REAL(BCElemID)
!      END IF
    END DO ! iBCSide
  END DO ! iElem
END IF ! fullMesh

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
  SideID   = INT(SideBCMetrics(BCSIDE_SIDEID,iSide))
  BCSideID = SideID2BCSide(SideID)
  CNElemID = GetCNElemID(INT(SideBCMetrics(BCSIDE_ELEMID,iSide)))

  !> get origin and radius from BC Side
  SideBCMetrics(5:7          ,iSide) = BCSideMetrics(1:3,BCSideID)
  SideBCMetrics(BCSIDE_RADIUS,iSide) = BCSideMetrics(4  ,BCSideID)

  !> build side distance
  origin(1:3) = ElemBaryNGeo(1:3,CNElemID)
  vec(1:3)    = origin(1:3) - SideBCMetrics(5:7,iSide)
  SideBCMetrics(BCSIDE_DISTANCE,iSide) = SQRT(DOT_PRODUCT(vec,vec))-ElemRadiusNGeo(CNElemID)-SideBCMetrics(BCSIDE_RADIUS,iSide)
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
  firstSide    = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) + 1
  lastSide     = ElemToBCSides(ELEM_FIRST_BCSIDE,iElem) + ElemToBCSides(ELEM_NBR_BCSIDES,iElem)
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
USE MOD_ChangeBasis              ,ONLY: ChangeBasis3D
USE MOD_Interpolation            ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars       ,ONLY: NodeTypeCL,NodeType
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Mesh_Vars                ,ONLY: nElems
USE MOD_Particle_Mesh_Vars       ,ONLY: nComputeNodeElems
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemsJ,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars       ,ONLY: RefMappingEps
USE MOD_Particle_Mesh_Tools      ,ONLY: GetGlobalElemID
#if USE_MPI
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Mesh_Vars                ,ONLY: NGeo,NGeoRef
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars       ,ONLY: dXCL_NGeo_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemsJ_Shared,ElemEpsOneCell_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemsJ_Shared_Win,ElemEpsOneCell_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,firstElem,lastElem
REAL                           :: scaleJ,maxScaleJ
#if USE_MPI
INTEGER                        :: ElemID
INTEGER                        :: i,j,k
! Vandermonde matrices
REAL                           :: Vdm_CLNGeo_NGeoRef(0:NgeoRef,0:NGeo)
REAL                           :: Vdm_NGeoRef_N(     0:PP_N   ,0:NGeoRef)
! Jacobian on CL N and NGeoRef
REAL                           :: detJac_Ref(1  ,0:NGeoRef,0:NGeoRef,0:NGeoRef)
REAL                           :: DetJac_N(  1  ,0:PP_N,   0:PP_N,   0:PP_N)
! interpolation points and derivatives on CL N
REAL                           :: dX_NGeoRef(3,3,0:NGeoRef,0:NGeoRef,0:NGeoRef)

INTEGER                        :: ElemLocID
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Building EpsOneCell for all elements...'

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
#else
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

#if USE_MPI
! Calculate sJ for elements not inside current proc, otherwise copy local values
CALL GetVandermonde(    Ngeo   , NodeTypeCL  , NgeoRef , NodeType  , Vdm_CLNGeo_NGeoRef, modal=.FALSE.)
CALL GetVandermonde(    NgeoRef, NodeType    , PP_N    , NodeType  , Vdm_NGeoRef_N     , modal=.TRUE.)

DO iElem = firstElem,lastElem
  ElemID    = GetGlobalElemID(iElem)
  ElemLocID = ElemID-offsetElem
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
ElemsJ = sJ(:,:,:,:,0)
#endif /* USE_MPI*/

! allocate epsOneCell
#if USE_MPI
MPISharedSize = INT((nComputeNodeTotalElems),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeTotalElems/),ElemEpsOneCell_Shared_Win,ElemEpsOneCell_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemEpsOneCell_Shared_Win,IERROR)
ElemEpsOneCell => ElemEpsOneCell_Shared
#else
ALLOCATE(ElemEpsOneCell(1:nComputeNodeElems))
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
USE MOD_Mathtools                ,ONLY: CROSS
USE MOD_Mesh_Vars                ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars       ,ONLY: NGeoElevated
USE MOD_Particle_Mesh_Vars       ,ONLY: nNonUniqueGlobalSides
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierElevation
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierControlPoints3D,BezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars   ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3!,BaseVectorsScale
#if USE_MPI
USE MOD_Particle_MPI_Shared      ,ONLY: Allocate_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,myComputeNodeRank
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_SHARED
!USE MOD_Particle_Mesh_Vars       ,ONLY: BaseVectorsScale_Shared,BaseVectorsScale_Shared_Win
USE MOD_Particle_Mesh_Vars       ,ONLY: BaseVectors0_Shared,BaseVectors1_Shared,BaseVectors2_Shared,BaseVectors3_Shared
USE MOD_Particle_Mesh_Vars       ,ONLY: BaseVectors0_Shared_Win,BaseVectors1_Shared_Win,BaseVectors2_Shared_Win,BaseVectors3_Shared_Win
#endif /* USE_MPI */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSide,firstSide,lastSide
!REAL                           :: crossVec(3)
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' GET LINEAR SIDE BASEVECTORS...'
#if USE_MPI
MPISharedSize = INT((3*nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors0_Shared_Win,BaseVectors0_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors0_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors1_Shared_Win,BaseVectors1_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors1_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors2_Shared_Win,BaseVectors2_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors2_Shared_Win,IERROR)
CALL Allocate_Shared(MPISharedSize,(/3,nNonUniqueGlobalSides/),BaseVectors3_Shared_Win,BaseVectors3_Shared)
CALL MPI_WIN_LOCK_ALL(0,BaseVectors3_Shared_Win,IERROR)
!MPISharedSize = INT((nNonUniqueGlobalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
!CALL Allocate_Shared(MPISharedSize,(/nNonUniqueGlobalSides/),BaseVectorsScale_Shared_Win,BaseVectorsScale_Shared)
!CALL MPI_WIN_LOCK_ALL(0,BaseVectorsScale_Shared_Win,IERROR)
BaseVectors0 => BaseVectors0_Shared
BaseVectors1 => BaseVectors1_Shared
BaseVectors2 => BaseVectors2_Shared
BaseVectors3 => BaseVectors3_Shared
!BaseVectorsScale => BaseVectorsScale_Shared

firstSide = INT(REAL (myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE( BaseVectors0(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors1(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors2(1:3,1:nNonUniqueGlobalSides),&
          BaseVectors3(1:3,1:nNonUniqueGlobalSides)) !,&
!          BaseVectorsScale(1:nNonUniqueGlobalSides))

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
!    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
!    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
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
!    crossVec = CROSS(BaseVectors1(:,iSide),BaseVectors2(:,iSide)) !vector with length of approx. 4x area (BV12 have double length)
!    BaseVectorsScale(iSide) = 0.25*SQRT(DOT_PRODUCT(crossVec,crossVec))
  END DO ! iSide
END IF

#if USE_MPI
CALL MPI_WIN_SYNC(BaseVectors0_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors1_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors2_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BaseVectors3_Shared_Win,IERROR)
!CALL MPI_WIN_SYNC(BaseVectorsScale_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /* USE_MPI */

SWRITE(UNIT_stdOut,'(A)')' GET LINEAR SIDE BASEVECTORS DONE!'
!SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE GetLinearSideBaseVectors


SUBROUTINE GetMeshMinMax()
!===================================================================================================================================
! computes the minimum and maximum value of the mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Mesh_Vars               ,ONLY: offsetElem,nElems
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeSide,nComputeNodeSides
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeNode,nComputeNodeNodes
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: offsetLocalSide,nLocalSides
INTEGER                        :: offsetLocalNode,nLocalNodes
!===================================================================================================================================

SELECT CASE(TrackingMethod)
  ! Build mesh min/max on BezierControlPoints for possibly curved elements
  CASE(REFMAPPING,TRACING)
    ! calculate all offsets
    offsetLocalSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1)
    nLocalSides     = ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElem+nElems)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1)

    ! proc local
    GEO%xmin     = MINVAL(BezierControlPoints3D(1,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%xmax     = MAXVAL(BezierControlPoints3D(1,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%ymin     = MINVAL(BezierControlPoints3D(2,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%ymax     = MAXVAL(BezierControlPoints3D(2,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%zmin     = MINVAL(BezierControlPoints3D(3,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))
    GEO%zmax     = MAXVAL(BezierControlPoints3D(3,:,:,offsetLocalSide+1:offsetLocalSide+nLocalSides))

#if USE_MPI
    ! compute-node local
    GEO%CNxmin   = MINVAL(BezierControlPoints3D(1,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNxmax   = MAXVAL(BezierControlPoints3D(1,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNymin   = MINVAL(BezierControlPoints3D(2,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNymax   = MAXVAL(BezierControlPoints3D(2,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNzmin   = MINVAL(BezierControlPoints3D(3,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))
    GEO%CNzmax   = MAXVAL(BezierControlPoints3D(3,:,:,offsetComputeNodeSide+1:offsetComputeNodeSide+nComputeNodeSides))

    ! global
    GEO%xminglob = MINVAL(BezierControlPoints3D(1,:,:,:))
    GEO%xmaxglob = MAXVAL(BezierControlPoints3D(1,:,:,:))
    GEO%yminglob = MINVAL(BezierControlPoints3D(2,:,:,:))
    GEO%ymaxglob = MAXVAL(BezierControlPoints3D(2,:,:,:))
    GEO%zminglob = MINVAL(BezierControlPoints3D(3,:,:,:))
    GEO%zmaxglob = MAXVAL(BezierControlPoints3D(3,:,:,:))
#endif /*USE_MPI*/
  ! TriaTracking does not have curved elements, nodeCoords are sufficient
  CASE(TRIATRACKING)
    ! calculate all offsets
    offsetLocalNode = ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)
    nLocalNodes     = ElemInfo_Shared(ELEM_LASTNODEIND ,offsetElem+nElems)-ElemInfo_Shared(ELEM_FIRSTNODEIND,offsetElem+1)

    ! proc local
    GEO%xmin     = MINVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%xmax     = MAXVAL(NodeCoords_Shared(1,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%ymin     = MINVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%ymax     = MAXVAL(NodeCoords_Shared(2,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%zmin     = MINVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))
    GEO%zmax     = MAXVAL(NodeCoords_Shared(3,offsetLocalNode+1:offsetLocalNode+nLocalNodes))

#if USE_MPI
    ! compute-node local
    GEO%CNxmin   = MINVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNxmax   = MAXVAL(NodeCoords_Shared(1,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNymin   = MINVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNymax   = MAXVAL(NodeCoords_Shared(2,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNzmin   = MINVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))
    GEO%CNzmax   = MAXVAL(NodeCoords_Shared(3,offsetComputeNodeNode+1:offsetComputeNodeNode+nComputeNodeNodes))

    ! global
    GEO%xminglob = MINVAL(NodeCoords_Shared(1,:))
    GEO%xmaxglob = MAXVAL(NodeCoords_Shared(1,:))
    GEO%yminglob = MINVAL(NodeCoords_Shared(2,:))
    GEO%ymaxglob = MAXVAL(NodeCoords_Shared(2,:))
    GEO%zminglob = MINVAL(NodeCoords_Shared(3,:))
    GEO%zmaxglob = MAXVAL(NodeCoords_Shared(3,:))
#endif /*USE_MPI*/
END SELECT

#if !USE_MPI
! compute-node local (dummy)
GEO%CNxmin   = GEO%xmin
GEO%CNxmax   = GEO%xmax
GEO%CNymin   = GEO%ymin
GEO%CNymax   = GEO%ymax
GEO%CNzmin   = GEO%zmin
GEO%CNzmax   = GEO%zmax

! global (dummy)
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
! Deallocates variables for the particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_BGM                ,ONLY: FinalizeBGM
USE MOD_Particle_Interpolation_Vars ,ONLY: DoInterpolation
USE MOD_Particle_Mesh_Readin        ,ONLY: FinalizeMeshReadin
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars
USE MOD_Particle_Tracking_Vars      ,ONLY: TrackingMethod,Distance,ListDistance
#if USE_MPI
USE MOD_Particle_MPI_Shared_vars    ,ONLY: MPI_COMM_SHARED
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Particle mesh readin happens during mesh readin, finalize with gathered routine here
CALL FinalizeMeshReadin()

CALL FinalizeBGM()

! Particle mesh metrics
SDEALLOCATE(DCL_N)
SDEALLOCATE(Vdm_CLN_GaussN)
SDEALLOCATE(Vdm_CLNGeo_GaussN)
SDEALLOCATE(Vdm_CLNGeo_CLN)
SDEALLOCATE(Vdm_CLNGeo1_CLNGeo)
SDEALLOCATE(Vdm_NGeo_CLNGeo)
SDEALLOCATE(Vdm_Bezier)
SDEALLOCATE(sVdm_Bezier)
SDEALLOCATE(wBaryCL_NGeo)
SDEALLOCATE(wBaryCL_NGeo1)
SDEALLOCATE(Xi_NGeo)
SDEALLOCATE(XiCL_NGeo)
SDEALLOCATE(XiCL_NGeo1)
SDEALLOCATE(DCL_NGeo)
SDEALLOCATE(D_Bezier)
SDEALLOCATE(XCL_NGeo)
SDEALLOCATE(dXCL_NGeo)

SELECT CASE(TrackingMethod)

  ! RefMapping, Tracing
  CASE(REFMAPPING,TRACING)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! InitElemVolumes
    CALL MPI_WIN_UNLOCK_ALL(ElemVolume_Shared_Win           ,iError)
    CALL MPI_WIN_FREE(      ElemVolume_Shared_Win           ,iError)

    ! InitParticleGeometry
    IF (TriaSurfaceFlux) THEN
      CALL MPI_WIN_UNLOCK_ALL(ConcaveElemSide_Shared_Win      ,iError)
      CALL MPI_WIN_FREE(      ConcaveElemSide_Shared_Win      ,iError)
      CALL MPI_WIN_UNLOCK_ALL(ElemNodeID_Shared_Win           ,iError)
      CALL MPI_WIN_FREE(      ElemNodeID_Shared_Win           ,iError)
      CALL MPI_WIN_UNLOCK_ALL(ElemSideNodeID_Shared_Win       ,iError)
      CALL MPI_WIN_FREE(      ElemSideNodeID_Shared_Win       ,iError)
      CALL MPI_WIN_UNLOCK_ALL(ElemMidPoint_Shared_Win         ,iError)
      CALL MPI_WIN_FREE(      ElemMidPoint_Shared_Win         ,iError)
    END IF

    ! GetBCSidesAndOrgin
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      CALL MPI_WIN_UNLOCK_ALL(BCSide2SideID_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      BCSide2SideID_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(SideID2BCSide_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      SideID2BCSide_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(BCSideMetrics_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      BCSideMetrics_Shared_Win        ,iError)
    END IF

#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ! CalcParticleMeshMetrics
      CALL MPI_WIN_UNLOCK_ALL(XCL_NGeo_Shared_Win             ,iError)
      CALL MPI_WIN_FREE(      XCL_NGeo_Shared_Win             ,iError)
      CALL MPI_WIN_UNLOCK_ALL(Elem_xGP_Shared_Win             ,iError)
      CALL MPI_WIN_FREE(      Elem_xGP_Shared_Win             ,iError)
      CALL MPI_WIN_UNLOCK_ALL(dXCL_NGeo_Shared_Win            ,iError)
      CALL MPI_WIN_FREE(      dXCL_NGeo_Shared_Win            ,iError)

      ! CalcBezierControlPoints
      CALL MPI_WIN_UNLOCK_ALL(BezierControlPoints3D_Shared_Win,iError)
      CALL MPI_WIN_FREE(      BezierControlPoints3D_Shared_Win,iError)
      IF (BezierElevation.GT.0) THEN
        CALL MPI_WIN_UNLOCK_ALL(BezierControlPoints3DElevated_Shared_Win,iError)
        CALL MPI_WIN_FREE(      BezierControlPoints3DElevated_Shared_Win,iError)
      END IF
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    ! GetSideSlabNormalsAndIntervals (allocated in particle_mesh.f90)
    CALL MPI_WIN_UNLOCK_ALL(SideSlabNormals_Shared_Win      ,iError)
    CALL MPI_WIN_FREE(      SideSlabNormals_Shared_Win      ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideSlabIntervals_Shared_Win    ,iError)
    CALL MPI_WIN_FREE(      SideSlabIntervals_Shared_Win    ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BoundingBoxIsEmpty_Shared_Win   ,iError)
    CALL MPI_WIN_FREE(      BoundingBoxIsEmpty_Shared_Win   ,iError)

    ! BuildElementOriginShared
    CALL MPI_WIN_UNLOCK_ALL(ElemBaryNGeo_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      ElemBaryNGeo_Shared_Win         ,iError)

    ! IdentifyElemAndSideType
    CALL MPI_WIN_UNLOCK_ALL(ElemCurved_Shared_Win           ,iError)
    CALL MPI_WIN_FREE(      ElemCurved_Shared_Win           ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideType_Shared_Win             ,iError)
    CALL MPI_WIN_FREE(      SideType_Shared_Win             ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideDistance_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      SideDistance_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(SideNormVec_Shared_Win          ,iError)
    CALL MPI_WIN_FREE(      SideNormVec_Shared_Win          ,iError)

    ! BuildElementBasisAndRadius
    CALL MPI_WIN_UNLOCK_ALL(ElemRadiusNGeo_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      ElemRadiusNGeo_Shared_Win       ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemRadius2NGeo_Shared_Win      ,iError)
    CALL MPI_WIN_FREE(      ElemRadius2NGeo_Shared_Win      ,iError)
    CALL MPI_WIN_UNLOCK_ALL(XiEtaZetaBasis_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      XiEtaZetaBasis_Shared_Win       ,iError)
    CALL MPI_WIN_UNLOCK_ALL(slenXiEtaZetaBasis_Shared_Win   ,iError)
    CALL MPI_WIN_FREE(      slenXiEtaZetaBasis_Shared_Win   ,iError)

    ! GetLinearSideBaseVectors
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors0_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors0_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors1_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors1_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors2_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors2_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(BaseVectors3_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      BaseVectors3_Shared_Win         ,iError)
!    CALL MPI_WIN_UNLOCK_ALL(BaseVectorsScale_Shared_Win     ,iError)
!    CALL MPI_WIN_FREE(      BaseVectorsScale_Shared_Win     ,iError)

    ! BuildBCElemDistance
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      CALL MPI_WIN_UNLOCK_ALL(ElemToBCSides_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      ElemToBCSides_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(SideBCMetrics_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      SideBCMetrics_Shared_Win        ,iError)
    END IF

    ! BuildEpsOneCell
    CALL MPI_WIN_UNLOCK_ALL(ElemsJ_Shared_Win               ,iError)
    CALL MPI_WIN_FREE(      ElemsJ_Shared_Win               ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemEpsOneCell_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      ElemEpsOneCell_Shared_Win       ,iError)

    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

    ! Then, free the pointers or arrays
    ! InitElemVolumes
    MDEALLOCATE(ElemVolume_Shared)

    ! InitParticleGeometry
    IF (TriaSurfaceFlux) THEN
      MDEALLOCATE(ConcaveElemSide_Shared)
      MDEALLOCATE(ElemNodeID_Shared)
      MDEALLOCATE(ElemSideNodeID_Shared)
      MDEALLOCATE(ElemMidPoint_Shared)
    END IF

    ! GetBCSidesAndOrgin
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      MDEALLOCATE(BCSide2SideID_Shared)
      MDEALLOCATE(SideID2BCSide_Shared)
      MDEALLOCATE(BCSideMetrics_Shared)
      MDEALLOCATE(BCSide2SideID)
      MDEALLOCATE(SideID2BCSide)
      MDEALLOCATE(BCSideMetrics)
    END IF

#if USE_LOADBALANCE
    IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
      ! CalcParticleMeshMetrics
      MDEALLOCATE(XCL_NGeo_Array)
      MDEALLOCATE(Elem_xGP_Array)
      MDEALLOCATE(dXCL_NGeo_Array)
      IF(ASSOCIATED(XCL_NGeo_Shared))  NULLIFY(XCL_NGeo_Shared)
      IF(ASSOCIATED(Elem_xGP_Shared))  NULLIFY(Elem_xGP_Shared)
      IF(ASSOCIATED(dXCL_NGeo_Shared)) NULLIFY(dXCL_NGeo_Shared)

      ! CalcBezierControlPoints
      MDEALLOCATE(BezierControlPoints3D)
      MDEALLOCATE(BezierControlPoints3D_Shared)
      MDEALLOCATE(BezierControlPoints3DElevated)
      MDEALLOCATE(BezierControlPoints3DElevated_Shared)
#if USE_LOADBALANCE
    END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

    ! GetSideSlabNormalsAndIntervals (allocated in particle_mesh.f90)
    MDEALLOCATE(SideSlabNormals)
    MDEALLOCATE(SideSlabNormals_Shared)
    MDEALLOCATE(SideSlabIntervals)
    MDEALLOCATE(SideSlabIntervals_Shared)
    MDEALLOCATE(BoundingBoxIsEmpty)
    MDEALLOCATE(BoundingBoxIsEmpty_Shared)

    ! BuildElementOriginShared
    MDEALLOCATE(ElemBaryNGeo)
    MDEALLOCATE(ElemBaryNGeo_Shared)

    ! IdentifyElemAndSideType
    MDEALLOCATE(ElemCurved)
    MDEALLOCATE(ElemCurved_Shared)
    MDEALLOCATE(SideType)
    MDEALLOCATE(SideType_Shared)
    MDEALLOCATE(SideDistance)
    MDEALLOCATE(SideDistance_Shared)
    MDEALLOCATE(SideNormVec)
    MDEALLOCATE(SideNormVec_Shared)

    ! BuildElementBasisAndRadius
    MDEALLOCATE(ElemRadiusNGeo)
    MDEALLOCATE(ElemRadiusNGeo_Shared)
    MDEALLOCATE(ElemRadius2NGeo)
    MDEALLOCATE(ElemRadius2NGeo_Shared)
    MDEALLOCATE(XiEtaZetaBasis)
    MDEALLOCATE(XiEtaZetaBasis_Shared)
    MDEALLOCATE(slenXiEtaZetaBasis)
    MDEALLOCATE(slenXiEtaZetaBasis_Shared)

    ! GetLinearSideBaseVectors
    MDEALLOCATE(BaseVectors0)
    MDEALLOCATE(BaseVectors0_Shared)
    MDEALLOCATE(BaseVectors1)
    MDEALLOCATE(BaseVectors1_Shared)
    MDEALLOCATE(BaseVectors2)
    MDEALLOCATE(BaseVectors2_Shared)
    MDEALLOCATE(BaseVectors3)
    MDEALLOCATE(BaseVectors3_Shared)
!    MDEALLOCATE(BaseVectorsScale)
!    MDEALLOCATE(BaseVectorsScale_Shared)

    ! BuildBCElemDistance
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      MDEALLOCATE(ElemToBCSides)
      MDEALLOCATE(ElemToBCSides_Shared)
      MDEALLOCATE(SideBCMetrics)
      MDEALLOCATE(SideBCMetrics_Shared)
    END IF

    ! BuildEpsOneCell
    MDEALLOCATE(ElemsJ)
    MDEALLOCATE(ElemsJ_Shared)
    MDEALLOCATE(ElemEpsOneCell)
    MDEALLOCATE(ElemEpsOneCell_Shared)

!  ! Tracing
!  CASE(TRACING)

  ! TriaTracking
  CASE(TRIATRACKING)
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! InitElemVolumes
    CALL MPI_WIN_UNLOCK_ALL(ElemVolume_Shared_Win           ,iError)
    CALL MPI_WIN_FREE(      ElemVolume_Shared_Win           ,iError)

    IF (DoInterpolation) THEN
#if USE_LOADBALANCE
      IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
        ! CalcParticleMeshMetrics
        CALL MPI_WIN_UNLOCK_ALL(XCL_NGeo_Shared_Win             ,iError)
        CALL MPI_WIN_FREE(      XCL_NGeo_Shared_Win             ,iError)
        CALL MPI_WIN_UNLOCK_ALL(Elem_xGP_Shared_Win             ,iError)
        CALL MPI_WIN_FREE(      Elem_xGP_Shared_Win             ,iError)
        CALL MPI_WIN_UNLOCK_ALL(dXCL_NGeo_Shared_Win            ,iError)
        CALL MPI_WIN_FREE(      dXCL_NGeo_Shared_Win            ,iError)
#if USE_LOADBALANCE
      END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

      ! BuildElemTypeAndBasisTria
      CALL MPI_WIN_UNLOCK_ALL(ElemCurved_Shared_Win           ,iError)
      CALL MPI_WIN_FREE(      ElemCurved_Shared_Win           ,iError)
      CALL MPI_WIN_UNLOCK_ALL(XiEtaZetaBasis_Shared_Win       ,iError)
      CALL MPI_WIN_FREE(      XiEtaZetaBasis_Shared_Win       ,iError)
      CALL MPI_WIN_UNLOCK_ALL(slenXiEtaZetaBasis_Shared_Win   ,iError)
      CALL MPI_WIN_FREE(      slenXiEtaZetaBasis_Shared_Win   ,iError)
    END IF ! DoInterpolation

    ! InitParticleGeometry
    CALL MPI_WIN_UNLOCK_ALL(ConcaveElemSide_Shared_Win      ,iError)
    CALL MPI_WIN_FREE(      ConcaveElemSide_Shared_Win      ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemNodeID_Shared_Win           ,iError)
    CALL MPI_WIN_FREE(      ElemNodeID_Shared_Win           ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemSideNodeID_Shared_Win       ,iError)
    CALL MPI_WIN_FREE(      ElemSideNodeID_Shared_Win       ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemMidPoint_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      ElemMidPoint_Shared_Win         ,iError)

    ! BuildElementRadiusTria
    CALL MPI_WIN_UNLOCK_ALL(ElemBaryNGeo_Shared_Win         ,iError)
    CALL MPI_WIN_FREE(      ElemBaryNGeo_Shared_Win         ,iError)
    CALL MPI_WIN_UNLOCK_ALL(ElemRadius2NGeo_Shared_Win      ,iError)
    CALL MPI_WIN_FREE(      ElemRadius2NGeo_Shared_Win      ,iError)

    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

    ! Then, free the pointers or arrays
    ! InitElemVolumes
    MDEALLOCATE(ElemVolume_Shared)

    IF (DoInterpolation) THEN
#if USE_LOADBALANCE
      IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
        ! CalcParticleMeshMetrics
        MDEALLOCATE(XCL_NGeo_Array)
        MDEALLOCATE(Elem_xGP_Array)
        MDEALLOCATE(dXCL_NGeo_Array)
        IF(ASSOCIATED(XCL_NGeo_Shared))  NULLIFY(XCL_NGeo_Shared)
        IF(ASSOCIATED(Elem_xGP_Shared))  NULLIFY(Elem_xGP_Shared)
        IF(ASSOCIATED(dXCL_NGeo_Shared)) NULLIFY(dXCL_NGeo_Shared)
#if USE_LOADBALANCE
      END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

      ! BuildElemTypeAndBasisTria
      MDEALLOCATE(ElemCurved)
      MDEALLOCATE(ElemCurved_Shared)
      MDEALLOCATE(XiEtaZetaBasis)
      MDEALLOCATE(XiEtaZetaBasis_Shared)
      MDEALLOCATE(slenXiEtaZetaBasis)
      MDEALLOCATE(slenXiEtaZetaBasis_Shared)
    END IF

    ! InitParticleGeometry
    MDEALLOCATE(ConcaveElemSide_Shared)
    MDEALLOCATE(ElemNodeID_Shared)
    MDEALLOCATE(ElemSideNodeID_Shared)
    MDEALLOCATE(ElemMidPoint_Shared)

    ! BuildElementRadiusTria
    MDEALLOCATE(ElemBaryNGeo)
    MDEALLOCATE(ElemBaryNGeo_Shared)
    MDEALLOCATE(ElemRadius2NGeo)
    MDEALLOCATE(ElemRadius2NGeo_Shared)

END SELECT

! Background mesh
#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
  SDEALLOCATE(GEO%PeriodicVectors)
END IF
#endif /*USE_LOADBALANCE*/
SDEALLOCATE(GEO%FIBGM)

! Tracking Vars
SDEALLOCATE(Distance)
SDEALLOCATE(ListDistance)

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
      r_vec  = AuxBC_plane(AuxBCMap(iAuxBC))%r_vec
      n_vec  = AuxBC_plane(AuxBCMap(iAuxBC))%n_vec
      radius = AuxBC_plane(AuxBCMap(iAuxBC))%radius

      ! loop over all  elements
      DO iElem = 1,nElems
        ! 1-2: Min, Max value; 1-3: x,y,z
        ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) )

        fmin = -DOT_PRODUCT(r_vec,n_vec)
        fmax = fmin
        DO icoord=1,3
          IF (n_vec(icoord).GE.0) THEN
            fmin = fmin + n_vec(icoord)*Bounds(1,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(2,icoord)
          ELSE
            fmin = fmin + n_vec(icoord)*Bounds(2,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(1,icoord)
          END IF
        END DO

        ! plane intersects the box!
        IF ((fmin.LE.0 .AND. fmax.GT.0).OR.(fmin.LT.0 .AND. fmax.GE.0)) THEN
          !radius check needs to be implemented (compute intersection polygon and minimum radii): would sort out further elements!!!
          !quick, conservative solution: calculate bounding box of disc in space and compare with bb of element
          ElemHasAuxBCs(iElem,iAuxBC) = .TRUE.

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
        ! plane does not intersect the box!
        ELSE IF ((fmin.LT.0 .AND. fmax.LT.0).OR.(fmin.GT.0 .AND. fmax.GT.0)) THEN
          ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
        ! e.g. if elem has zero volume...
        ELSE
          CALL abort(__STAMP__,'Error in MarkAuxBCElems for AuxBC:',iAuxBC)
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

      cartesian = .TRUE.
      backwards = .FALSE.
      IF (ABS(n_vec(1)).EQ.1.) THEN
        dir(1)=1
        dir(2)=2
        dir(3)=3
        IF (n_vec(1).LT.0.) backwards = .TRUE.
      ELSE IF (ABS(n_vec(2)).EQ.1.) THEN
        dir(1)=2
        dir(2)=3
        dir(3)=1
        IF (n_vec(2).LT.0.) backwards = .TRUE.
      ELSE IF (ABS(n_vec(3)).EQ.1.) THEN
        dir(1)=3
        dir(2)=1
        dir(3)=2
        IF (n_vec(3).LT.0.) backwards = .TRUE.
      ELSE
        cartesian=.FALSE.
        SWRITE(*,*) 'WARNING in MarkAuxBCElems: all Elems are set to ElemHasAuxBCs=.TRUE. for AuxBC:',iAuxBC
        ElemHasAuxBCs(:,iAuxBC) = .TRUE. !actual intersection with box check to-be implemented!!!
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
          ! 1-2: Min, Max value; 1-3: x,y,z
          ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) )

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
      ! to be implemented!!!
      ElemHasAuxBCs(:,iAuxBC) = .TRUE.

    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(__STAMP__,'AuxBC does not exist')
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
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
DO iDir1 = 0,1
  DO iDir2 = 0,1
      DO iDir3 = 0,1
        BoundingBox(1,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir1+1,1)
        BoundingBox(2,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir2+1,2)
        BoundingBox(3,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir3+1,3)
      END DO
  END DO
END DO

!-- check where the points are located relative to radius
done = .FALSE.
DO iDir1 = 0,1
  IF(done) EXIT

  DO iDir2  =0,1
    IF(done) EXIT

    DO iDir3 = 0,1
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

      ! iPoint.GT.1: type must be 2 if state of point if different from last point
      ELSE
        IF (pointRadius.LE.radius) THEN
          ! different from last point
          IF (.NOT.insideBound) THEN
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        ! outside
        ELSE
          ! different from last point
          IF (insideBound) THEN
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
    ! points are really outside
    ELSE
      positiontype=1
    END IF
  END IF
END IF

END SUBROUTINE CheckBoundsWithCartRadius


END MODULE MOD_Particle_Mesh
