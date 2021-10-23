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
! Builds the mesh for particle tracking, separate from the DG mesh
!===================================================================================================================================
MODULE MOD_Particle_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersParticleMesh
    MODULE PROCEDURE DefineParametersParticleMesh
END INTERFACE

INTERFACE InitParticleMeshBasis
    MODULE PROCEDURE InitParticleMeshBasis
END INTERFACE

INTERFACE InitParticleMesh
  MODULE PROCEDURE InitParticleMesh
END INTERFACE

INTERFACE FinalizeParticleMeshBasis
    MODULE PROCEDURE FinalizeParticleMeshBasis
END INTERFACE

INTERFACE FinalizeParticleMesh
  MODULE PROCEDURE FinalizeParticleMesh
END INTERFACE

PUBLIC :: DefineParametersParticleMesh
PUBLIC :: InitParticleMeshBasis
PUBLIC :: InitParticleMesh
PUBLIC :: FinalizeParticleMeshBasis
PUBLIC :: FinalizeParticleMesh
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
! CALL prms%CreateIntOption(          'Part-nPeriodicVectors'     , 'Number of the periodic vectors j=1,...,n.'                    //&
!                                                                   ' Value has to be the same as defined in preprog.ini'            &
!                                                                 , '0')
CALL prms%CreateLogicalOption(      'Part-CartesianPeriodic'    , ' Simplified treatment for periodic box with Refmapping. Not'  //&
                                                                  ' computation of intersection points at periodic BCs.'           &
                                                                , '.FALSE.')
CALL prms%CreateLogicalOption(      'Part-FastPeriodic'           , ' Further simplification by directly moving particle into'   //&
                                                                  ' grid. Instead of moving the particle several times the'      //&
                                                                  ' periodic displacements, the particle is mapped directly back'//&
                                                                  ' into the domain.'                                              &
                                                                , '.FALSE.')
! CALL prms%CreateLogicalOption(      'Part-DoPeriodicCheck'      , ' Check position in reference element after periodic displacement'&
!                                                                 , '.FALSE.')
! CALL prms%CreateLogicalOption(      'Part-DoPeriodicFix'        , ' Fix position in reference element after periodic displacement' &
!                                                                 , '.FALSE.')

! Mesh checking
CALL prms%CreateLogicalOption(      'meshCheckWeirdElements'    , 'Abort when weird elements are found: it means that part of'   //&
                                                                  ' the element is turned inside-out. '                            &
                                                                , '.TRUE.')

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
USE MOD_Particle_Basis         ,ONLY: BuildBezierVdm,BuildBezierDMat
USE MOD_Particle_BGM           ,ONLY: BuildBGMAndIdentifyHaloRegion
USE MOD_Particle_Globals
USE MOD_Particle_Interpolation_Vars, ONLY: DoInterpolation
USE MOD_Particle_Mesh_Build    ,ONLY: GetMeshMinMax,IdentifyElemAndSideType,ComputePeriodicVec
USE MOD_Particle_Mesh_Build    ,ONLY: BuildElementRadiusTria,BuildElemTypeAndBasisTria,BuildEpsOneCell,BuildBCElemDistance
USE MOD_Particle_Mesh_Build    ,ONLY: BuildElementOriginShared,BuildElementBasisAndRadius
USE MOD_Particle_Mesh_Build    ,ONLY: BuildSideOriginAndRadius,BuildLinearSideBaseVectors
USE MOD_Particle_Mesh_Build    ,ONLY: CalcParticleMeshMetrics,InitParticleGeometry,CalcBezierControlPoints
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: InitGetGlobalElemID,InitGetCNElemID,GetCNElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: InitGetGlobalSideID,InitGetCNSideID,GetGlobalSideID
USE MOD_Particle_Surfaces      ,ONLY: GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierElevation,BezierControlPoints3D,BezierControlPoints3DElevated
USE MOD_Particle_Surfaces_Vars ,ONLY: SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars ,ONLY: D_Bezier,Vdm_Bezier,sVdm_Bezier
USE MOD_Particle_Tracking_Vars ,ONLY: FastPeriodic,CartesianPeriodic,CountNbOfLostParts
USE MOD_Particle_Tracking_Vars ,ONLY: NbrOfLostParticles,NbrOfLostParticlesTotal,NbrOfLostParticlesTotal_old
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierSampleN,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod, DisplayLostParticles
USE MOD_Particle_Tracking_Vars ,ONLY: PartStateLostVecLength,PartStateLost,PartLostDataSize
! USE MOD_Particle_Tracking_Vars ,ONLY: DoPeriodicCheck,DoPeriodicFix
USE MOD_ReadInTools            ,ONLY: GETREAL,GETINT,GETLOGICAL,GetRealArray
#if CODE_ANALYZE
! USE MOD_Particle_Surfaces_Vars ,ONLY: SideBoundingBoxVolume
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_Particle_BGM           ,ONLY: WriteHaloInfo
USE MOD_Particle_MPI_Shared    ,ONLY: Allocate_Shared,MPI_SIZE,BARRIER_AND_SYNC
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
#if !USE_MPI
INTEGER          :: ALLOCSTAT
#endif
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(132("-"))')
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

! Initialize Vandermonde for Bezier basis surface representation (particle tracking with curved elements)
BezierElevation = GETINT('BezierElevation','0')
NGeoElevated    = NGeo + BezierElevation

CALL BuildBezierVdm              (NGeo,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier)
CALL BuildBezierDMat             (NGeo,Xi_NGeo,D_Bezier)

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
  ! Nullify and reset lost parts container after write out
  PartStateLostVecLength = 0

  ! Allocate PartStateLost for a small number of particles and double the array size each time the
  ! maximum is reached
  ALLOCATE(PartStateLost(1:PartLostDataSize,1:10))
  PartStateLost=0.
END IF ! CountNbrOfLostParts
NbrOfLostParticles      = 0
#if USE_LOADBALANCE
! Nullify only once at the beginning of a simulation (not during load balance steps!)
IF(.NOT.PerformLoadBalance)THEN
#endif
  NbrOfLostParticlesTotal = 0
  NbrOfLostParticlesTotal_old = 0
#if USE_LOADBALANCE
END IF
#endif
DisplayLostParticles    = GETLOGICAL('DisplayLostParticles')

#if CODE_ANALYZE
PARTOUT            = GETINT('PartOut','0')
MPIRankOut         = GETINT('MPIRankOut','0')
#endif /*CODE_ANALYZE*/

CartesianPeriodic = GETLOGICAL('Part-CartesianPeriodic','.FALSE.')
IF (CartesianPeriodic) FastPeriodic = GETLOGICAL('Part-FastPeriodic','.FALSE.')

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
    ! TriaSurfaceFlux needs ElemMidPoint_Shared
    IF (TriaSurfaceFlux) THEN
      CALL InitParticleGeometry()
    END IF

     ! Calculated before BGM to find convex hull
!    CALL CalcParticleMeshMetrics()

!    CALL CalcBezierControlPoints()

#if USE_MPI
    CALL Allocate_Shared((/3,3,nComputeNodeTotalSides/),SideSlabNormals_Shared_Win,SideSlabNormals_Shared)
    CALL MPI_WIN_LOCK_ALL(0,SideSlabNormals_Shared_Win,IERROR)
    CALL Allocate_Shared((/6,nComputeNodeTotalSides/),SideSlabIntervals_Shared_Win,SideSlabIntervals_Shared)
    CALL MPI_WIN_LOCK_ALL(0,SideSlabIntervals_Shared_Win,IERROR)
    CALL Allocate_Shared((/nComputeNodeTotalSides/),BoundingBoxIsEmpty_Shared_Win,BoundingBoxIsEmpty_Shared)
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
! #if CODE_ANALYZE
!     ALLOCATE(SideBoundingBoxVolume(nSides))
! #endif /*CODE_ANALYZE*/

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
    CALL BARRIER_AND_SYNC(SideSlabNormals_Shared_Win   ,MPI_COMM_SHARED)
    CALL BARRIER_AND_SYNC(SideSlabIntervals_Shared_Win ,MPI_COMM_SHARED)
    CALL BARRIER_AND_SYNC(BoundingBoxIsEmpty_Shared_Win,MPI_COMM_SHARED)
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
      ! SideBoundingBoxVolume(SideID)=dx*dy*dz
    END DO
#endif /*CODE_ANALYZE*/

    ! Compute element bary and element radius for node elements (with halo region)
    CALL BuildElementOriginShared()

    ! Check the side type (planar, bilinear, curved)
    CALL IdentifyElemAndSideType()

    ! Compute the element XiEtaZetaBasis and the radius of the convex hull
    CALL BuildElementBasisAndRadius()

    ! Get basevectors for (bi-)linear sides
    CALL BuildLinearSideBaseVectors()

    IF (TrackingMethod.EQ.REFMAPPING) THEN
      ! Identify BCSides and build side origin and radius
      CALL BuildSideOriginAndRadius()

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
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitParticleMesh


SUBROUTINE FinalizeParticleMeshBasis()
!===================================================================================================================================
! Finalize basic metrics needed for particle mesh
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

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

END SUBROUTINE FinalizeParticleMeshBasis


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
USE MOD_LoadBalance_Vars            ,ONLY: PerformLoadBalance
#else
USE MOD_LoadBalance_Vars            ,ONLY: ElemTime
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

    ! BuildSideOriginAndRadius
    IF (TrackingMethod.EQ.REFMAPPING) THEN
      CALL MPI_WIN_UNLOCK_ALL(BCSide2SideID_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      BCSide2SideID_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(SideID2BCSide_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      SideID2BCSide_Shared_Win        ,iError)
      CALL MPI_WIN_UNLOCK_ALL(BCSideMetrics_Shared_Win        ,iError)
      CALL MPI_WIN_FREE(      BCSideMetrics_Shared_Win        ,iError)
    END IF

#if USE_LOADBALANCE
    ! BezierControlPoints are global and do not change during load balance
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

    ! BuildLinearSideBaseVectors
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

    ! BuildSideOriginAndRadius
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

    ! BuildLinearSideBaseVectors
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

! Load Balance
#if !USE_LOADBALANCE
SDEALLOCATE(ElemTime)
#endif /* !USE_LOADBALANCE */

ParticleMeshInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleMesh

END MODULE MOD_Particle_Mesh
