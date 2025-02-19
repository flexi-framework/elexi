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
! Contains subroutines to build the requiered data to track particles on (curviilinear) meshes, etc.
!===================================================================================================================================
MODULE MOD_Particle_Surfaces
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: InitParticleSurfaces
PUBLIC:: FinalizeParticleSurfaces
PUBLIC:: GetBezierControlPoints3DElevated
PUBLIC:: GetSideSlabNormalsAndIntervals
PUBLIC:: GetSideBoundingBox
PUBLIC:: GetBezierSampledAreas
PUBLIC:: EvaluateBezierPolynomialAndGradient
PUBLIC:: CalcNormAndTangBilinear
PUBLIC:: CalcNormAndTangBezier
PUBLIC:: CalcNormAndTangTriangle
PUBLIC:: RotateMasterToSlave
#if CODE_ANALYZE
PUBLIC:: OutputBezierControlPoints
#endif /*CODE_ANALYZE*/
!===================================================================================================================================


CONTAINS

SUBROUTINE InitParticleSurfaces()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,               ONLY: epsMach
USE MOD_PreProc
USE MOD_Mesh_Vars,                  ONLY: NGeo
USE MOD_Particle_Surfaces_Vars
USE MOD_ReadInTools,                ONLY: GETREAL,GETINT,GETLOGICAL,PrintOption
USE MOD_Utils,                      ONLY: ALMOSTZERO
#if CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,     ONLY: rBoundingBoxChecks,rPerformBezierClip,rPerformBezierNewton
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: tmp
REAL                            :: epsilonrel
CHARACTER(LEN=2)                :: dummy
!===================================================================================================================================

IF(ParticleSurfaceInitIsDone) RETURN
! LBWRITE(UNIT_stdOut,'(132("."))')
LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SURFACES...'

BezierNewtonAngle          = GETREAL('BezierNewtonAngle'         )  ! 1°=0.01754 (in rad)
BezierClipTolerance        = GETREAL('BezierClipTolerance'       )
BezierNewtonTolerance2     = GETREAL('BezierNewtonTolerance'     )
BezierNewtonGuess          = GETINT( 'BezierNewtonGuess'         )
BezierNewtonMaxIter        = GETINT( 'BezierNewtonMaxIter'       )
BezierSplitLimit           = GETREAL('BezierSplitLimit'          )
BezierSplitLimit           = 2.*BezierSplitLimit
BezierClipMaxIter          = GETINT( 'BezierClipMaxIter'         )
BezierClipLineVectorMethod = GETINT( 'BezierClipLineVectorMethod')
BezierClipHit              = GETREAL('BezierClipHit'             )
IF(ALMOSTZERO(BezierClipHit)) BezierClipHit = 100.*BezierClipTolerance
BezierClipHit   = 1.0 + BezierClipHit
BezierNewtonHit            = GETREAL('BezierNewtonHit'           )
IF(ALMOSTZERO(BezierNewtonHit)) BezierNewtonHit=BezierNewtonTolerance2
BezierNewtonHit = 1.0 + (BezierNewtonHit) ! it is not clear how this value should be determined
BezierNewtonTolerance2 = BezierNewtonTolerance2**2
tmp             = 2*(NGeo+1)
WRITE(dummy,'(I2.2)') tmp
BezierClipMaxIntersec      = GETINT( 'BezierClipMaxIntersec'     ,dummy)

epsilontol            = GETREAL('epsilontol')
! it might be smarter to give the tolerance as a relationship to the machine epsilon (epsMach)
IF(ALMOSTZERO(epsilontol)) THEN
  epsilonrel          = GETREAL('epsilonrel')
  IF(.NOT.ALMOSTZERO(epsilonrel)) THEN
    epsilontol=epsilonrel*EpsMach
  ! if nothing is entered, than a default value is used
  ! for tolerance issuses see, e.g. Haselbacher et al., "An efficient and robust particle-localization algorithm for unstructured grids"
  ! epsilon approx 100*tolerance of the algorithm
  ELSE
    epsilontol=100.*EpsMach
  END IF
END IF
CALL PrintOption('epsilon (abs.,rel. to epsMach)','INFO',RealArrayOpt=(/epsilontol,epsilontol/EpsMach/))
MinusEps              = -epsilontol
OnePlusEps            = 1.0 + 100.*epsilontol
OneMinusEps           = 1.0 - epsilontol

#if CODE_ANALYZE
rBoundingBoxChecks   = 0.
rPerformBezierClip   = 0.
rPerformBezierNewton = 0.
rTotalBBChecks       = 0.
rTotalBezierClips    = 0.
rTotalBezierNewton   = 0.
#endif /*CODE_ANALYZE*/

ALLOCATE( locAlpha(1:BezierClipMaxIntersec) &
        , locXi   (1:BezierClipMaxIntersec) &
        , locEta  (1:BezierClipMaxIntersec) &
        , XiArray (1:2,1:BezierClipMaxIter) &
        , EtaArray(1:2,1:BezierClipMaxIter) )

ParticleSurfaceInitIsDone = .TRUE.

LBWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SURFACES DONE!'
LBWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitParticleSurfaces


SUBROUTINE FinalizeParticleSurfaces()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars            ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! InitParticleMesh
SDEALLOCATE(BezierSampleXi)
! #if CODE_ANALYZE
! SDEALLOCATE(SideBoundingBoxVolume)
! #endif

#if USE_LOADBALANCE
! BezierControlPoints are global and do not change during load balance
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  ! InitParticleMeshBasis
  SDEALLOCATE(Vdm_Bezier)
  SDEALLOCATE(sVdm_Bezier)
  SDEALLOCATE(D_Bezier)

  ! CalcBezierControlPoints (MPI3 shared freed in FinalizeParticleMesh)
  SDEALLOCATE(XiBuf)
  SDEALLOCATE(MinMax)
  SDEALLOCATE(XiUp)
  SDEALLOCATE(XiDown)
  SDEALLOCATE(BezierControlPoints1D)
  SDEALLOCATE(BezierControlPoints2D)
  SDEALLOCATE(BezierControlPoints2D_temp)
  SDEALLOCATE(BezierControlPoints2D_temp2)

  MDEALLOCATE(BezierControlPoints3D)
  MDEALLOCATE(BezierControlPoints3DElevated)
#if USE_LOADBALANCE
END IF !PerformLoadBalance
#endif /*USE_LOADBALANCE*/

! GetSideSlabNormalsAndIntervals (MPI3 shared freed in FinalizeParticleMesh)
MDEALLOCATE(SideSlabNormals)
MDEALLOCATE(SideSlabIntervals)
MDEALLOCATE(BoundingBoxIsEmpty)

! GetLinearSideBaseVectors (MPI3 shared freed in FinalizeParticleMesh)
MDEALLOCATE(BaseVectors0)
MDEALLOCATE(BaseVectors1)
MDEALLOCATE(BaseVectors2)
MDEALLOCATE(BaseVectors3)
!MDEALLOCATE(BaseVectorsScale)

! BuildBezierVdm
SDEALLOCATE(arrayNChooseK)
SDEALLOCATE(FacNchooseK)
SDEALLOCATE(ElevationMatrix)

! InitParticleSurfaces
SDEALLOCATE(locAlpha)
SDEALLOCATE(locXi)
SDEALLOCATE(locEta)
SDEALLOCATE(XiArray)
SDEALLOCATE(EtaArray)

! Surface flux
SDEALLOCATE(SurfMeshSubSideData)
SDEALLOCATE(SurfMeshSideAreas)
SDEALLOCATE(BCdata_auxSF)

ParticleSurfaceInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleSurfaces


SUBROUTINE CalcNormAndTangTriangle(nVec,tang1,tang2,area,midpoint,ndist,xyzNod,Vectors,TriNum,SideID,ElemID_opt,LocSideID_opt)
!================================================================================================================================
! function to compute the geo-data of a triangulated surface
!================================================================================================================================
USE MOD_Globals,                              ONLY: Abort
USE MOD_PreProc
USE MOD_Mesh_Vars,                            ONLY: SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars,                   ONLY: ElemSideNodeID_Shared
USE MOD_Particle_Mesh_Tools,                  ONLY: GetCNElemID,GetCNSideID
USE MOD_Particle_Surfaces_Vars,               ONLY: SideNormVec,SideType
USE MOD_Particle_Tracking_Vars,               ONLY: TrackingMethod
USE MOD_Utils,                                ONLY: ALMOSTZERO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                     :: TriNum
INTEGER,INTENT(IN),OPTIONAL            :: SideID
INTEGER,INTENT(IN),OPTIONAL            :: ElemID_opt,LocSideID_opt
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,INTENT(OUT),OPTIONAL              :: nVec(3)
REAL,INTENT(OUT),OPTIONAL              :: tang1(3), tang2(3), area, midpoint(3), ndist(3)
REAL,INTENT(INOUT),OPTIONAL            :: xyzNod(3) ,Vectors(3,3)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: CNSideID
INTEGER                                :: CNElemID
INTEGER                                :: LocSideID
INTEGER                                :: Node1,Node2
REAL                                   :: xNod,zNod,yNod,Vector1(3),Vector2(3)
REAL                                   :: nVal,ndistVal,nx,ny,nz,dotpr
!================================================================================================================================

IF (PRESENT(ElemID_opt).AND.PRESENT(LocSideID_opt)) THEN
  CNElemID  = GetCNElemID(ElemID_opt)
  LocSideID = LocSideID_opt
ELSE IF (PRESENT(SideID)) THEN
  CNElemID  = GetCNElemID(SideInfo_Shared(SIDE_ELEMID ,SideID))
  LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
ELSE
  CNElemID  = -1 ! Initial value to eliminate compiler warning
  LocSideID = -1 ! Initial value to eliminate compiler warning
  CALL Abort(__STAMP__, 'Either SideID or ElemID+LocSideID have to be given to CalcNormAndTangTriangle!')
END IF

xNod = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)
yNod = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)
zNod = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)

IF(PRESENT(xyzNod) .AND. TriNum.EQ.1) THEN !only write for first Tria
  xyzNod = (/xNod,yNod,zNod/)
END IF

Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

Vector1(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node1,LocSideID,CNElemID)+1) - xNod
Vector1(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node1,LocSideID,CNElemID)+1) - yNod
Vector1(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node1,LocSideID,CNElemID)+1) - zNod
Vector2(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node2,LocSideID,CNElemID)+1) - xNod
Vector2(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node2,LocSideID,CNElemID)+1) - yNod
Vector2(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node2,LocSideID,CNElemID)+1) - zNod

nx = - Vector1(2) * Vector2(3) + Vector1(3) * Vector2(2) !NV (inwards)
ny = - Vector1(3) * Vector2(1) + Vector1(1) * Vector2(3)
nz = - Vector1(1) * Vector2(2) + Vector1(2) * Vector2(1)
nVal = SQRT(nx*nx + ny*ny + nz*nz)
nx = -nx / nVal
ny = -ny / nVal
nz = -nz / nVal

IF (TrackingMethod.NE.TRIATRACKING) THEN
  CNSideID = GetCNSideID(SideID)
  IF ((SideType(CNSideID).EQ.PLANAR_RECT .OR. SideType(CNSideID).EQ.PLANAR_NONRECT)) THEN
    !if surfflux-side are planar, TriaSurfaceFlux can be also used for tracing or Refmapping (for which SideNormVec exists)!
    !warning: these values go into SurfMeshSubSideData and if TriaSurfaceFlux they should be used only for planar_rect/_nonrect sides
    dotpr = DOT_PRODUCT(SideNormVec(1:3,CNSideID),(/nx,ny,nz/))
    IF ( .NOT.ALMOSTEQUALRELATIVE(dotpr,1.,1.0E-2) ) &
      CALL Abort(__STAMP__, 'SideNormVec is not identical with V1xV2!')
  END IF
END IF

IF (PRESENT(Vectors) .AND. TriNum.EQ.1) THEN
  Vectors(:,1) = Vector1
  Vectors(:,2) = Vector2
ELSE IF (PRESENT(Vectors)) THEN !TriNum.EQ.2
  Vectors(:,3) = Vector2
END IF

IF(PRESENT(nVec)) THEN
  nVec = (/nx,ny,nz/)
END IF

IF(PRESENT(area)) THEN
  area = nVal*0.5
END IF

IF(PRESENT(midpoint)) THEN
  midpoint(1) = xNod + (Vector1(1)+Vector2(1))/2.
  midpoint(2) = yNod + (Vector1(2)+Vector2(2))/2.
  midpoint(3) = zNod + (Vector1(3)+Vector2(3))/2.
END IF

IF(PRESENT(ndist)) THEN
  ndist(1) = ny * (Vector1(3)-Vector2(3)) - nz * (Vector1(2)-Vector2(2)) !NV to Vec1-Vec2 in plane (outwards from tria)
  ndist(2) = nz * (Vector1(1)-Vector2(1)) - nx * (Vector1(3)-Vector2(3))
  ndist(3) = nx * (Vector1(2)-Vector2(2)) - ny * (Vector1(1)-Vector2(1))
  ndistVal = SQRT(ndist(1)*ndist(1) + ndist(2)*ndist(2) + ndist(3)*ndist(3))
  IF (ALMOSTZERO(ndistVal)) CALL Abort(&
__STAMP__&
, 'ndistVal=0 in InitSurfaceflux!')
  ndist(1:3) = ndist(1:3) / ndistVal
END IF

IF(PRESENT(tang1) .OR. PRESENT(tang2)) THEN
  IF (.NOT.ALMOSTZERO(nz)) THEN
    tang1(1) = 1.0
    tang1(2) = 1.0
    tang1(3) = -(nx+ny)/nz
    tang2(1) = ny * tang1(3) - nz
    tang2(2) = nz - nx * tang1(3)
    tang2(3) = nx - ny
    tang1 = tang1 / SQRT(2.0 + tang1(3)*tang1(3))
  ELSE
    IF (.NOT.ALMOSTZERO(ny)) THEN
      tang1(1) = 1.0
      tang1(3) = 1.0
      tang1(2) = -(nx+nz)/ny
      tang2(1) = ny - nz * tang1(2)
      tang2(2) = nz - nx
      tang2(3) = nx * tang1(2) - ny
      tang1 = tang1 / SQRT(2.0 + tang1(2)*tang1(2))
    ELSE
      IF (.NOT.ALMOSTZERO(nx)) THEN
        tang1(2) = 1.0
        tang1(3) = 1.0
        tang1(1) = -(ny+nz)/nx
        tang2(1) = ny - nz
        tang2(2) = nz * tang1(1) - nx
        tang2(3) = nx - ny * tang1(1)
        tang1 = tang1 / SQRT(2.0 + tang1(1)*tang1(1))
      ELSE
        CALL Abort(__STAMP__,&
          'Error in InitParticles: normal vector is zero!')
      END IF
    END IF
  END IF
tang2 = tang2 / SQRT(tang2(1)*tang2(1) + tang2(2)*tang2(2) + tang2(3)*tang2(3))
END IF

END SUBROUTINE CalcNormAndTangTriangle


SUBROUTINE CalcNormAndTangBilinear(nVec,tang1,tang2,xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface described by 4 corner nodes
!================================================================================================================================
USE MOD_Particle_Globals,                     ONLY:CROSSNORM,UNITVECTOR
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,INTENT(OUT)                       :: nVec(3)
REAL,INTENT(OUT),OPTIONAL              :: tang1(3), tang2(3)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)                      :: a,b
!================================================================================================================================

b=xi*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)  &
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
   +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)   &
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)+BezierControlPoints3D(:,0   ,NGeo,SideID) )

a=eta*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
    +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) )

nVec=CROSSNORM(a,b)
IF(PRESENT(tang1)) tang1=UNITVECTOR(a)
IF(PRESENT(tang2)) tang2=CROSSNORM(nVec,a)

END SUBROUTINE CalcNormAndTangBilinear


SUBROUTINE CalcNormAndTangBezier(nVec,tang1,tang2,xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface described by a Bezier polynomial
!================================================================================================================================
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Particle_Globals,                     ONLY:CROSSNORM,UNITVECTOR
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,DIMENSION(3),INTENT(OUT)          :: nVec
REAL,DIMENSION(3),INTENT(OUT),OPTIONAL :: tang1,tang2
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(2,3)                    :: gradXiEta
!================================================================================================================================

! caution we require the formula in [0;1]
CALL EvaluateBezierPolynomialAndGradient((/xi,eta/),NGeo,3,BezierControlPoints3D(:,:,:,SideID),Gradient=gradXiEta)
nVec =CROSSNORM(gradXiEta(1,:),gradXiEta(2,:))
IF(PRESENT(tang1)) tang1=UNITVECTOR(gradXiEta(1,:))
IF(PRESENT(tang2)) tang2=CROSSNORM(nVec,gradXiEta(1,:))

END SUBROUTINE CalcNormAndTangBezier


SUBROUTINE EvaluateBezierPolynomialAndGradient(Xi,N_in,iSize,BezierControlPoints,dBezierControlPoints,Point,Gradient)
!===================================================================================================================================
! evaluate the BezierPolynomial and/or gradient by means of direct evaluation and help of a pseudo-horner scheme
! faster for the evaluation of gradients
! results depend on mode
! mode=1 compute/interpolate point in patch at xi
! mode=2 compute gradient in xi at xi
! mode=3 point and gradient || both
! range: [-1,1]
! info for gradient
! grad(1,:) are all gradients in u direction,...
! book:
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
! 1) computation of Bezier Polynomial is from Farin, Chap. 5, p. 57, EQ. 1
!    precompute t^i (1-t)^(n-i) by means of Horner Schema and store values direcly sorted
! 2) point value
! 3) evaluate gradients
! optimized in terms of floating point operations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars,               ONLY:facNchooseK
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)              :: Xi(2)                 ! value of the plane (u,v)
INTEGER,INTENT(IN)           :: N_in
INTEGER,INTENT(IN)           :: iSize                 ! depending if two or three dimensional
REAL,INTENT(IN)              :: BezierControlPoints(1:iSize,0:N_in,0:N_in)
REAL,INTENT(IN),OPTIONAL     :: dBezierControlPoints(1:2,1:iSize,0:N_in,0:N_in)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT),OPTIONAL    :: Point(1:iSize)        ! point value
REAL,INTENT(OUT),OPTIONAL    :: Gradient(1:2,1:iSize) ! gradient in u,v
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: Mode,iMode
INTEGER                      :: p,q,m
REAL                         :: MinusXi,PlusXI,MinusEta,PlusEta
REAL                         :: xiup(0:N_in),etaup(0:N_in),xidown(0:N_in),etadown(0:N_in)
REAL                         :: xitild(0:N_in), etatild(0:N_in)
!===================================================================================================================================

Mode=0
IF(PRESENT(Point)) Mode=1
iMode=0
IF(Present(Gradient)) iMode=2
Mode=Mode+iMode

M=N_in-1
MinusXi =1.0-xi(1)
PlusXi  =1.0+xi(1)
MinusEta=1.0-xi(2)
PlusEta =1.0+xi(2)

! compute the required stuff || pseudo Horner or precomputation
xiup(0)=1.0
etaup(0)=1.0
! caution, here the indicies are switched from n-j to j  for **down
xidown(N_in)=1.0
etadown(N_in)=1.0
DO p=1,N_in
  xiup   (p)     =xiup   (p-1)*PlusXi
  etaup  (p)     =etaup  (p-1)*PlusEta
  xidown (N_in-p)=xidown (N_in-p+1)*MinusXi
  etadown(N_in-p)=etadown(N_in-p+1)*MinusEta
END DO ! p

! only once
DO p=0,N_in
  xitild (p)=xiup (p)*xidown (p)
  etatild(p)=etaup(p)*etadown(p)
END DO ! p=0,N_in

IF ((Mode.EQ.1).OR.(Mode.EQ.3)) THEN
  Point=0.
  DO q=0,N_in
    DO p=0,N_in
      ! derivative in xi: d/dxi = u
      !Point(:)=Point(:)+BezierControlPoints(:,p,q)              &
      !                 *facNchooseK(N_In,p)*xiup(p) *xidown(p)  &
      !                 *facNChooseK(N_In,q)*etaup(q)*etadown(q)
      Point(:)=Point(:)+BezierControlPoints(:,p,q)     &
                       *facNchooseK(N_In,p)*xitild (p) &
                       *facNChooseK(N_In,q)*etatild(q)
    END DO ! p
  END DO ! q
  !Point=Point
END IF

! gradient
IF (Mode.GE.2) THEN
  gradient=0.
  IF(PRESENT(dBezierControlPoints))THEN
    DO q=0,N_in
      DO p=0,N_in
        gradient(:,:)=gradient(:,:)+dBezierControlPoints(:,:,p,q)  &
                                   *facNchooseK(N_In,p)*xitild (p) &
                                   *facNChooseK(N_In,q)*etatild(q)
      END DO ! p
    END DO ! q
  ELSE
    DO q=0,N_in
      DO p=0,M
        ! derivative in xi: d/dxi = u
        gradient(1,:)=gradient(1,:)+(BezierControlPoints(:,p+1,q)-BezierControlPoints(:,p,q)) &
                                   *facNchooseK(M,p)*xiup(p)*xidown(p+1)                      &
                                   *facNChooseK(N_in,q)*etatild(q)
        ! derivative in eta: d/deta = v ! caution - exchange indicies
        gradient(2,:)=gradient(2,:)+(BezierControlPoints(:,q,p+1)-BezierControlPoints(:,q,p)) &
                                   *facNchooseK(N_in,q)*xitild(q)                             &
                                   *facNchooseK(M,p)*etaup(p)*etadown(p+1)
      END DO ! p
    END DO ! q
    ! caution, maybe one 1/(b-a) missing??
    gradient=0.5*REAL(N_in)*gradient
  END IF
END IF

END SUBROUTINE EvaluateBezierPolynomialAndGradient


SUBROUTINE GetSideSlabNormalsAndIntervals(BezierControlPoints3D,SideSlabNormals,SideSlabInterVals,BoundingBoxIsEmpty)
!===================================================================================================================================
! computes the oriented-slab box for each bezier basis surface (i.e. 3 slab normals + 3 intervals)
! see article:
!    author = {Shyue-wu Wang and Zen-chung Shih and Ruei-chuan Chang},
!    title = {An Efficient and Stable Ray Tracing Algorithm for Parametric Surfaces},
!    year = {2001},
! original article: oriented-slab box (cartesian bounding box)
!   author = {Yen, Jonathan and Spach, Susan and Smith, Mark and Pulleyblank, Ron},
!   title = {Parallel Boxing in B-Spline Intersection},
!   issue_date = {January 1991},
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Globals,         ONLY: CROSSNORM
USE MOD_Particle_Mesh_Vars,       ONLY: NGeoElevated
USE MOD_Utils,                    ONLY: ALMOSTZERO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: BezierControlPoints3D(1:3,0:NGeoElevated,0:NGeoElevated)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: SideSlabNormals(1:3,1:3)
REAL,INTENT(OUT)    :: SideSlabInterVals(1:6)
LOGICAL,INTENT(OUT) :: BoundingBoxIsEmpty
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q, i
REAL               :: skalprod(3),dx,dy,dz,dMax,w,h,l
LOGICAL            :: SideIsCritical
!===================================================================================================================================

!-----------------------------------------------------------------------------------------------------------------------------------
! 0.) check if side is planar
!-----------------------------------------------------------------------------------------------------------------------------------
! done in particle_mesh

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) slab normal vectors
!-----------------------------------------------------------------------------------------------------------------------------------
! n_1=V_1+V_2 (V: corner vectors in xi-direction)
SideSlabNormals(:,1) = BezierControlPoints3D(:,NGeoElevated,0)              &
                     - BezierControlPoints3D(:,0,0)                         &
                     + BezierControlPoints3D(:,NGeoElevated,NGeoElevated)   &
                     - BezierControlPoints3D(:,0,NGeoElevated)

IF (ALL(SideSlabNormals(:,1).EQ.0)) &
  CALL Abort(__STAMP__,             &
  'Error while calculating side slab normals and intervals. Normal vector length zero. Possibly wrong BezierControlPoints')

SideSlabNormals(:,1) = SideSlabNormals(:,1)/SQRT(DOT_PRODUCT(SideSlabNormals(:,1),SideSlabNormals(:,1)))
! n_2=n_1 x (U_1+U_2) (U: corner vectors in eta-direction)
SideSlabNormals(:,2) = BezierControlPoints3D(:,0,NGeoElevated)                      &
                     - BezierControlPoints3D(:,0,0)                                         &
                     + BezierControlPoints3D(:,NGeoElevated,NGeoElevated)   &
                     - BezierControlPoints3D(:,NGeoElevated,0)

!fehlt das?
SideSlabNormals(:,2) = CROSSNORM(SideSlabNormals(:,1),SideSlabNormals(:,2))

! n_3=n_1 x n_2
SideSlabNormals(:,3) = CROSSNORM(SideSlabNormals(:,2),SideSlabNormals(:,1))

! check vector length=1
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1),SideSlabNormals(:,1))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 1 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,1),SideSlabNormals(:,1)))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,2),SideSlabNormals(:,2))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 2 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,2),SideSlabNormals(:,2)))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,3),SideSlabNormals(:,3))-1.)).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 3 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,3),SideSlabNormals(:,3)))

! check perpendicularity
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1),SideSlabNormals(:,2)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 1 and 2 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,1),SideSlabNormals(:,2))))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1),SideSlabNormals(:,3)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 1 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,1),SideSlabNormals(:,3))))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,2),SideSlabNormals(:,3)))).GT.1.E-6) CALL Abort(&
__STAMP__&
,'Side slab normal 2 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,2),SideSlabNormals(:,3))))
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) slab box intervals beta_1, beta_2, beta_3
!-----------------------------------------------------------------------------------------------------------------------------------
!SideSlabIntervals(x- x+ y- y+ z- z+, SideID)


! Interval beta_1
!print*,"SideID",SideID
SideSlabIntervals(:)=0.

DO q=0,NGeoElevated
  DO p=0,NGeoElevated
    IF((p.EQ.0).AND.(q.EQ.0))CYCLE
    skalprod(1)=DOT_PRODUCT(BezierControlPoints3D(:,p,q)-&
                            BezierControlPoints3D(:,0,0),SideSlabNormals(:,1))
    skalprod(2)=DOT_PRODUCT(BezierControlPoints3D(:,p,q)-&
                            BezierControlPoints3D(:,0,0),SideSlabNormals(:,2))
    skalprod(3)=DOT_PRODUCT(BezierControlPoints3D(:,p,q)-&
                            BezierControlPoints3D(:,0,0),SideSlabNormals(:,3))
    IF    (skalprod(1).LT.0.)THEN
      SideSlabIntervals(1)=MIN(SideSlabIntervals(1),skalprod(1))
    ELSEIF(skalprod(1).GT.0.)THEN
      SideSlabIntervals(2)=MAX(SideSlabIntervals(2),skalprod(1))
    END IF
    IF    (skalprod(2).LT.0.)THEN
      SideSlabIntervals(3)=MIN(SideSlabIntervals(3),skalprod(2))
    ELSEIF(skalprod(2).GT.0.)THEN
      SideSlabIntervals(4)=MAX(SideSlabIntervals(4),skalprod(2))
    END IF
    IF    (skalprod(3).LT.0.)THEN
      SideSlabIntervals(5)=MIN(SideSlabIntervals(5),skalprod(3))
    ELSEIF(skalprod(3).GT.0.)THEN
      SideSlabIntervals(6)=MAX(SideSlabIntervals(6),skalprod(3))
    END IF
  END DO !p
END DO !q

!-----------------------------------------------------------------------------------------------------------------------------------
! 2-b.) sanity check
!-----------------------------------------------------------------------------------------------------------------------------------

DO i = 1,3
  IF(SideSlabIntervals(2*i).LT.SideSlabIntervals(2*i-1))  CALL Abort(&
__STAMP__&
,' SideSlabIntervals are corrupted! ')
END DO

!-----------------------------------------------------------------------------------------------------------------------------------
! 2-c.) bounding box extension
!-----------------------------------------------------------------------------------------------------------------------------------

dx=ABS(SideSlabIntervals(2)-SideSlabIntervals(1))
dy=ABS(SideSlabIntervals(4)-SideSlabIntervals(3))
dz=ABS(SideSlabIntervals(6)-SideSlabIntervals(5))

!-----------------------------------------------------------------------------------------------------------------------------------
! 3.) Is Side critical? (particle path parallel to the larger surface, therefore numerous intersections are possilbe)
! from Wang et al. 2001 An Efficient and Stable ray tracing algorithm for parametric surfaces
! critical sides: h << w <= l
! h: height
! w: width (short side)
! l: length (long side)
!-----------------------------------------------------------------------------------------------------------------------------------

h=dy
w=MIN(dx,dz)
l=MAX(dx,dz)

IF((h.LT.w*1.E-2))THEN!.AND.(w.LE.l))THEN ! critical case possible: second condition is obsolete
  SideIsCritical=.TRUE.
ELSE
  SideIsCritical=.FALSE.
END IF


!-----------------------------------------------------------------------------------------------------------------------------------
! 4.) determine is side is planar -> "BoundingBoxIsEmpty(SideID)"
!     this results also in the decision whether a side is also considered flat or bilinear!
!-----------------------------------------------------------------------------------------------------------------------------------

! check dimensions of sideslabs in x, y and z direction where the slab area (w)x(l) is defined by x-z and the height of the side-slab is y
dMax=MAX(h,l,w)
IF(l/dMax.LT.1.0e-6)THEN
  ! side is almost a line
  CALL Abort(__STAMP__,'ERROR: found degenerated side. length/dMax of a side slab is ->',0,l/dMax)
END IF
IF(w/dMax.LT.1.0e-6)THEN
  ! side is almost a line
  CALL Abort(__STAMP__,'ERROR: found degenerated side. width/dMax of a side slab is ->',0,w/dMax)
END IF
IF(h/dMax.LT.1.0e-6)THEN
  ! the height of the sideslab is almost zero --> flat in regards to machine precision
  SideSlabIntervals(3:4)=0.
  h=0.
END IF

IF(l*w*h.LT.0) THEN
  CALL Abort(__STAMP__,'A bounding box (for sides) is negative!?. length*width*height.LT.0 ->',0,(l*w*h))
END IF

IF(ALMOSTZERO(h/SQRT(l*l+w*w+h*h)))THEN ! bounding box volume is approx zeros
  BoundingBoxIsEmpty=.TRUE.
ELSE
  BoundingBoxIsEmpty=.FALSE.
END IF

END SUBROUTINE GetSideSlabNormalsAndIntervals



SUBROUTINE GetSideBoundingBox(SideID, BoundingBox)
!===================================================================================================================================
! computes the 8 corners of bounding box of bezier basis surface (based on values from GetSideSlabNormalsAndIntervals)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Tools,        ONLY: GetCNSideID
USE MOD_Particle_Surfaces_Vars,     ONLY: BezierControlPoints3D,SideSlabIntervals,SideSlabNormals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: BoundingBox(1:3,1:8)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: CNSideID
INTEGER            :: iDir1,iDir2,iDir3
!===================================================================================================================================

CNSideID = GetCNSideID(SideID)

! BezierControlPoints are always on nonUniqueGlobalSide
DO iDir1 = 0,1
  DO iDir2 = 0,1
      DO iDir3 = 0,1
        BoundingBox(1:3,iDir1*4 + iDir2*2 + iDir3+1) = BezierControlPoints3D(:,0,0,SideID) &
          + SideSlabNormals(:,1,CNSideID)*SideSlabIntervals(2*1-iDir1,CNSideID) &
          + SideSlabNormals(:,2,CNSideID)*SideSlabIntervals(2*2-iDir2,CNSideID) &
          + SideSlabNormals(:,3,CNSideID)*SideSlabIntervals(2*3-iDir3,CNSideID)
      END DO
  END DO
END DO

END SUBROUTINE GetSideBoundingBox


SUBROUTINE GetBezierControlPoints3DElevated(NGeo,NGeoElevated,BezierControlPoints,BezierControlPointsElevated)
!===================================================================================================================================
! compute the elevated bezier control points (use pre-computed ElevationMatrix with binomial coefficicents)
! uses the tensor-product basis and combines two 1-D elevation steps
!===================================================================================================================================
! MODULES
!USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,BezierControlPoints3DElevated,BezierElevation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: NGeo,NGeoElevated
REAL,INTENT(IN)    :: BezierControlPoints(1:3,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: BezierControlPointsElevated(1:3,0:NGeoElevated,0:NGeoElevated)
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
REAL               :: temp(1:3,0:NGeoElevated,0:NGeo)
!===================================================================================================================================
temp=0.
! p-direction
DO q=0,NGeo
  temp(:,:,q) = ElevateBezierPolynomial(NGeo,BezierControlPoints(:,:,q))
END DO
! q-direction
DO p=0,NGeoElevated
  BezierControlPointsElevated(:,p,:) = ElevateBezierPolynomial(NGeo,temp(:,p,:))
END DO
END SUBROUTINE GetBezierControlPoints3DElevated


SUBROUTINE GetBezierSampledAreas(SideID,BezierSampleN,BezierSurfFluxProjection_opt,SurfMeshSubSideAreas,SurfMeshSideArea_opt &
                                ,SurfMeshSubSideVec_nOut_opt,SurfMeshSubSideVec_t1_opt,SurfMeshSubSideVec_t2_opt &
                                ,DmaxSampleN_opt,Dmax_opt,BezierControlPoints2D_opt)
!===================================================================================================================================
! equidistantly super-sampled Bezier surface area and vector calculation. Required for surface flux
! --------------------------------------
! book: see also for general remarks
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
! --------------------------------------
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Surfaces_Vars, ONLY:epsilontol,BezierControlPoints3D
USE MOD_Basis,                  ONLY:LegendreGaussNodesAndWeights
USE MOD_Mesh_Vars,              ONLY:NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: SideID,BezierSampleN
LOGICAL,INTENT(IN),OPTIONAL         :: BezierSurfFluxProjection_opt
INTEGER,INTENT(IN),OPTIONAL         :: DmaxSampleN_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:BezierSampleN,1:BezierSampleN),INTENT(OUT)              :: SurfMeshSubSideAreas
REAL,INTENT(OUT),OPTIONAL                                                :: SurfMeshSideArea_opt
REAL,DIMENSION(1:3,1:BezierSampleN,1:BezierSampleN),INTENT(OUT),OPTIONAL :: SurfMeshSubSideVec_nOut_opt &
                                                                           ,SurfMeshSubSideVec_t1_opt &
                                                                           ,SurfMeshSubSideVec_t2_opt
REAL,DIMENSION(1:BezierSampleN,1:BezierSampleN),INTENT(OUT),OPTIONAL     :: Dmax_opt
REAL,DIMENSION(1:2,0:NGeo,0:NGeo,1:BezierSampleN,1:BezierSampleN),INTENT(OUT),OPTIONAL       :: BezierControlPoints2D_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: p,q
INTEGER                                :: I,J,iSample,jSample,DmaxSampleN
REAL                                   :: areaTotal,areaTotalAbs,area,deltaXi,tmp1,E,F,G,D
REAL                                   :: tmpI2,tmpJ2
REAL                                   :: BezierControlPoints2D(1:2,0:NGeo,0:NGeo,1:BezierSampleN,1:BezierSampleN)
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,DIMENSION(2,2)                    :: gradXiEta2D,xiab
REAL,DIMENSION(2)                      :: XiOut(2)
REAL,DIMENSION(2)                      :: Xi
REAL,DIMENSION(3)                      :: n1,n2
REAL,ALLOCATABLE,DIMENSION(:)          :: Xi_NGeo,wGP_NGeo
LOGICAL                                :: BezierSurfFluxProjection, CalcDmax
REAL,DIMENSION(3)                      :: ProjectionVector
REAL,DIMENSION(0:BezierSampleN)        :: BezierSampleXi
REAL,DIMENSION(1:3,1:BezierSampleN,1:BezierSampleN) :: SurfMeshSubSideVec_nOut,SurfMeshSubSideVec_t1,SurfMeshSubSideVec_t2
REAL,DIMENSION(1:BezierSampleN,1:BezierSampleN) :: Dmax
!===================================================================================================================================

IF (PRESENT(BezierSurfFluxProjection_opt)) THEN; BezierSurfFluxProjection = BezierSurfFluxProjection_opt
ELSE                                           ; BezierSurfFluxProjection = .FALSE.
END IF
CalcDmax = .FALSE.

IF (PRESENT(Dmax_opt)) THEN
  IF (PRESENT(DmaxSampleN_opt) .AND. DmaxSampleN_opt.GT.0) THEN
    DmaxSampleN = DmaxSampleN_opt
    CalcDmax    = .TRUE.
  ELSE
    CALL Abort(__STAMP__,'DmaxSampleN not defined in GetBezierSampledAreas!')
  END IF
END IF

Xi(1)   =-SQRT(1./3.)
Xi(2)   = SQRT(1./3.)
deltaXi = 2.0/BezierSampleN
tmp1    = deltaXi/2.0 !(b-a)/2
DO jSample = 0,BezierSampleN
  BezierSampleXi(jSample) = -1.+deltaXi*jSample
END DO

!===================================================================================================================================

SurfMeshSubSideAreas    = 0.
SurfMeshSubSideVec_nOut = 0.
SurfMeshSubSideVec_t1   = 0.
SurfMeshSubSideVec_t2   = 0.
Dmax=1. !dummy

ALLOCATE(Xi_NGeo( 0:NGeo)  &
        ,wGP_NGeo(0:NGeo) )
CALL LegendreGaussNodesAndWeights(NGeo,Xi_NGeo,wGP_NGeo)

areaTotal    = 0.
areaTotalAbs = 0.

! loop through Sub-Elements
DO jSample = 1,BezierSampleN; DO iSample = 1,BezierSampleN

  area  = 0.
  tmpI2 = (BezierSampleXi(iSample-1)+BezierSampleXi(iSample))/2. ! (a+b)/2
  tmpJ2 = (BezierSampleXi(jSample-1)+BezierSampleXi(jSample))/2. ! (a+b)/2

  IF (BezierSurfFluxProjection             .OR. &
      PRESENT(SurfMeshSubSideVec_nOut_opt) .OR. &
      PRESENT(SurfMeshSubSideVec_t1_opt)   .OR. &
      PRESENT(SurfMeshSubSideVec_t2_opt))  THEN
    CALL CalcNormAndTangBezier( nVec=SurfMeshSubSideVec_nOut(:,iSample,jSample) &
                              ,tang1=SurfMeshSubSideVec_t1(  :,iSample,jSample) &
                              ,tang2=SurfMeshSubSideVec_t2(  :,iSample,jSample) &
                              ,xi=tmpI2,eta=tmpJ2,SideID=SideID )
  END IF

  IF (BezierSurfFluxProjection) THEN
    ! inwards normal is ProjVec
    ProjectionVector(1:3) = -SurfMeshSubSideVec_nOut(1:3,iSample,jSample)

    ! transformation of bezier patch 3D->2D
    IF(ABS(ProjectionVector(3)).LT.epsilontol)THEN
      n1 = (/ -ProjectionVector(2)-ProjectionVector(3) ,&
               ProjectionVector(1),ProjectionVector(1) /)
    ELSE
      n1 = (/  ProjectionVector(3),ProjectionVector(3) ,&
              -ProjectionVector(1)-ProjectionVector(2) /)
    END IF

    n1    = n1/SQRT(DOT_PRODUCT(n1,n1))
    n2(:) = CROSSNORM(ProjectionVector,n1)

    DO q = 0,NGeo; DO p = 0,NGeo
      BezierControlPoints2D(1,p,q,iSample,jSample) = DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID),n1)
      ! origin is (0,0,0)^T
      BezierControlPoints2D(2,p,q,iSample,jSample) = DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID),n2)
      ! origin is (0,0,0)^T
    END DO; END DO
  END IF ! BezierSurfFluxProjection

  ! calc integral
  DO I = 0,NGeo; DO J = 0,NGeo
    XiOut(1) = tmp1*Xi_NGeo(I)+tmpI2
    XiOut(2) = tmp1*Xi_NGeo(J)+tmpJ2

    IF (BezierSurfFluxProjection) THEN
      ! get gradients
      CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,2,BezierControlPoints2D(1:2,0:NGeo,0:NGeo,iSample,jSample) &
        ,Gradient=gradXiEta2D)
      ! calculate first fundamental form
      E = DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(1,1:2))
      F = DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(2,1:2))
      G = DOT_PRODUCT(gradXiEta2D(2,1:2),gradXiEta2D(2,1:2))
    ELSE
      CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(:,:,:,SideID) &
                                              ,Gradient=gradXiEta3D)
      ! calculate first fundamental form
      E = DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
      F = DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
      G = DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
    END IF

    D    = SQRT(E*G-F*F)
    area = area+tmp1*tmp1*D*wGP_NGeo(i)*wGP_NGeo(j)
  END DO; END DO ! I,J

  ! calc Dmax
  IF (PRESENT(Dmax_opt) .AND. CalcDmax) THEN
    Dmax(iSample,jSample) = -1.

    DO I = 0,DmaxSampleN; DO J = 0,DmaxSampleN
      xiab(1,1:2) = (/BezierSampleXi(ISample-1),BezierSampleXi(ISample)/) !correct order?!?
      xiab(2,1:2) = (/BezierSampleXi(JSample-1),BezierSampleXi(JSample)/) !correct order?!?
      XiOut = (xiab(:,2)-xiab(:,1))*(/ REAL(J) , REAL(I) /)/REAL(DmaxSampleN) + xiab(:,1)

      IF (BezierSurfFluxProjection) THEN
        ! get gradients
        CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,2,BezierControlPoints2D(1:2,0:NGeo,0:NGeo,iSample,jSample) &
          ,Gradient=gradXiEta2D)
        ! calculate first fundamental form
        E = DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(1,1:2))
        F = DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(2,1:2))
        G = DOT_PRODUCT(gradXiEta2D(2,1:2),gradXiEta2D(2,1:2))
      ELSE
        CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(:,:,:,SideID) &
          ,Gradient=gradXiEta3D)
        ! calculate first fundamental form
        E = DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
        F = DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
        G = DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
      END IF

      D = SQRT(E*G-F*F)
      Dmax(iSample,jSample) = MAX(Dmax(iSample,jSample),D)
    END DO; END DO ! I,J

    IF (Dmax(iSample,jSample).LT.0.) CALL Abort(__STAMP__,'ERROR in GetBezierSampledAreas: No Dmax found!?')
  END IF

  IF (BezierSurfFluxProjection) THEN
    ! add facing sides and substract non-facing sides
    SurfMeshSubSideAreas(iSample,jSample) = SIGN(area,-DOT_PRODUCT(ProjectionVector,SurfMeshSubSideVec_nOut(:,iSample,jSample)))
  ELSE
    SurfMeshSubSideAreas(iSample,jSample) = area
  END IF
  areaTotal    = areaTotal+SurfMeshSubSideAreas(iSample,jSample)
  areaTotalAbs = areaTotalAbs+area
END DO; END DO ! jSample=1,BezierSampleN;iSample=1,BezierSampleN: loop through Sub-Elements

DEALLOCATE(Xi_NGeo,wGP_NGeo)

IF (PRESENT(SurfMeshSideArea_opt))        SurfMeshSideArea_opt        =areaTotal
IF (PRESENT(SurfMeshSubSideVec_nOut_opt)) SurfMeshSubSideVec_nOut_opt =SurfMeshSubSideVec_nOut
IF (PRESENT(SurfMeshSubSideVec_t1_opt))   SurfMeshSubSideVec_t1_opt   =SurfMeshSubSideVec_t1
IF (PRESENT(SurfMeshSubSideVec_t2_opt))   SurfMeshSubSideVec_t2_opt   =SurfMeshSubSideVec_t2
IF (PRESENT(Dmax_opt))                    Dmax_opt                    =Dmax
IF (PRESENT(BezierControlPoints2D_opt)) THEN
  BezierControlPoints2D_opt = MERGE(BezierControlPoints2D,1.,BezierSurfFluxProjection)
END IF

#if CODE_ANALYZE
IPWRITE(UNIT_stdOut,'(I0,A)')            ' ===== SINGLE AREA ====='
IPWRITE(UNIT_stdOut,'(I0,A,I0,A,F16.7)') 'SideID:',SideID,'areaTotal   =',areaTotal
IPWRITE(UNIT_stdOut,'(I0,A,I0,A,F16.7)') 'SideID:',SideID,'areaTotalAbs=',areaTotalAbs
IPWRITE(UNIT_stdOut,'(I0,A)')            ' ============================'
#endif /*CODE_ANALYZE*/

END SUBROUTINE GetBezierSampledAreas


FUNCTION ElevateBezierPolynomial(NGeo,BezierPolynomial)
!===================================================================================================================================
! this function creates a new equidistantly distributed set of control points (bézier polynomial basis coefficients) based on the
! control points "BezierPolynomial" with order "N" and elevates them to "ElevateBezierPolynomial" on "Xi_NGeo_elevated" with order
! "p+BezierElevation"
! book: single step procedure (1D)
!   author = {Piegl, Les and Tiller, Wayne},
!   title = {The NURBS Book (2Nd Ed.)},
!   year = {1997},
! book: see also for general remarks
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
!===================================================================================================================================
! MODULES
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierElevation,ElevationMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: NGeo
REAL,INTENT(IN)    :: BezierPolynomial(1:3,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL               :: ElevateBezierPolynomial(1:3,0:NGeo+BezierElevation)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL            :: length
INTEGER            :: i,j,jStart,jEnd
!===================================================================================================================================
! algorithm originally from "The NURBS Book" by Les Piegl, Wayne Tiller (p.205)
ElevateBezierPolynomial = 0.
! the first and last points remain the same!
! edge points remain: P_0 and P_p
ElevateBezierPolynomial(:,0)                    = BezierPolynomial(:,0)
ElevateBezierPolynomial(:,NGeo+BezierElevation) = BezierPolynomial(:,NGeo)
! inner points change: P_1,...,P_p-1
DO i=1,NGeo+BezierElevation-1
  jStart = MAX(0,i-BezierElevation)
  jEnd   = MIN(NGeo,i)
  DO j=jStart,jEnd
    ElevateBezierPolynomial(:,i)=ElevateBezierPolynomial(:,i)+ElevationMatrix(i,j)*BezierPolynomial(:,j)
  END DO
END DO
END FUNCTION ElevateBezierPolynomial


SUBROUTINE RotateMasterToSlave(flip,BezierControlPoints3D)
!===================================================================================================================================
! Rotate BezierControlPoints3D from Master to Slave
! The resulting controlpoints are orientated as a "new master" for the slave element
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Mesh_Vars,      ONLY:nGeo
USE MOD_Mappings,       ONLY:Flip_M2S
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: flip
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)            :: BezierControlPoints3D(1:3,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: BezierControlPoints3D_tmp(1:3,0:NGeo,0:NGeo)
INTEGER                       :: p,q,pq(2)
!===================================================================================================================================

BezierControlPoints3D_tmp=BezierControlPoints3D

DO q=0,NGeo; DO p=0,NGeo
  pq = Flip_M2S(NGeo,p,q,flip,3)
  BezierControlPoints3D(:,pq(1),pq(2))=BezierControlPoints3d_tmp(:,p,q)
END DO; END DO ! p,q

END SUBROUTINE RotateMasterToSlave

#if CODE_ANALYZE
SUBROUTINE OutputBezierControlPoints(BezierControlPoints3D_in,BezierControlPoints2D_in)
!===================================================================================================================================
! dump the BezierControlPoints
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY: NGeo
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL  :: BezierControlPoints3d_in(1:3,0:NGeo,0:NGeo)
REAL,INTENT(IN),OPTIONAL  :: BezierControlPoints2d_in(1:2,0:NGeo,0:NGeo)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: k,i,j
REAL             :: BezierControlPoints3D(1:3,0:NGeo,0:NGeo)
!===================================================================================================================================

BezierControlPoints3D=0.
IF(PRESENT(BezierControlPoints3d_in))THEN
  BezierControlPoints3D=BezierControlPoints3d_in
ELSE IF(PRESENT(BezierControlPoints2d_in))THEN
  BezierControlPoints3D(1:2,:,:)=BezierControlPoints2d_in(:,:,:)
ELSE
  CALL Abort(&
__STAMP__&
,' OutputBezierControlPoints(): Something went wrong.!')
END IF

DO K=1,3
  WRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO')' P(:,:,',K,') = [ '
  DO I=0,NGeo ! output for MATLAB
    DO J=0,NGeo
      WRITE(UNIT_stdOut,'(E24.12)',ADVANCE='NO') BezierControlPoints3D(K,J,I)
      IF(J.EQ.NGeo)THEN
        IF(I.EQ.NGeo)THEN
          WRITE(UNIT_stdOut,'(A)')' ];'
        ELSE
          WRITE(UNIT_stdOut,'(A)')' ;...'
        END IF
      ELSE ! comma
        WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' , '
      END IF
    END DO
  END DO
END DO

END SUBROUTINE OutputBezierControlPoints
#endif /*CODE_ANALYZE*/

END MODULE MOD_Particle_Surfaces
