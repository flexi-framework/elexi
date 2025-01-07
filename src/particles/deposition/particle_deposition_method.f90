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
#include "eos.h"
#include "particle.h"

!===================================================================================================================================
! Module containing the different deposition methods (NGP, linear (inter-cell) weighting, shape function
!===================================================================================================================================
MODULE MOD_Particle_Deposition_Method
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

#if PARTICLES_COUPLING >= 2
PUBLIC:: DefineParametersDepositionMethod
PUBLIC:: InitDepositionMethod
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersDepositionMethod()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools              ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection('Particle Deposition')

CALL prms%CreateLogicalOption(      'Part-CalcSource'          , 'Flag to enable two-way coupling'                                 &
                                                               , '.FALSE.')
CALL prms%CreateIntFromStringOption('Part-DepositionType'      , 'Deposition method used for two-way coupling .\n'               //&
                                                      ' - 1) cell_vol        : imposed to element the particle resides in \n'    //&
                                                      ' - 2) cell_vol_linear : imposed linear continuous distribution \n'        //&
                                                      ' - 3) step            : step function \n'                                 //&
                                                      ' - 4) dirac           : dirac impulse \n'                                 //&
                                                      ' - 5) shapefunc_gauss : Gaussian shape function \n'                       //&
                                                      ' - 6) shapefunc_poly  : Polynomial shape function \n'                     &
                                                    , 'cell_vol')
CALL addStrListEntry('Part-DepositionType' , 'cell_vol'           , DEPO_CV)
CALL addStrListEntry('Part-DepositionType' , 'cell_vol_linear'    , DEPO_CVLM)
CALL addStrListEntry('Part-DepositionType' , 'step'               , DEPO_STEP)
CALL addStrListEntry('Part-DepositionType' , 'dirac'              , DEPO_DIRAC)
CALL addStrListEntry('Part-DepositionType' , 'shapefunc_gauss'    , DEPO_SF_GAUSS)
CALL addStrListEntry('Part-DepositionType' , 'shapefunc_poly'     , DEPO_SF_POLY)

END SUBROUTINE DefineParametersDepositionMethod

!==================================================================================================================================!
!> Initialize deposition method function pointer
!==================================================================================================================================!
SUBROUTINE InitDepositionMethod()
! MODULES
USE MOD_Globals
USE MOD_Particle_Deposition_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Select the deposition method function pointer
SELECT CASE(DepositionType)
  CASE(DEPO_CV)
    DepositionMethod => DepositionMethod_CV
  CASE(DEPO_CVLM)
    DepositionMethod => DepositionMethod_CVLM
  CASE(DEPO_STEP)
    DepositionMethod => DepositionMethod_Step
  CASE(DEPO_DIRAC)
    DepositionMethod => DepositionMethod_Dirac
  CASE(DEPO_SF_GAUSS)
    DepositionMethod => DepositionMethod_SF_Gauss
  CASE(DEPO_SF_POLY)
    DepositionMethod => DepositionMethod_SF_Poly
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'Unknown DepositionMethod!', IntInfo=DepositionType)
END SELECT

END SUBROUTINE InitDepositionMethod


!===================================================================================================================================
! Deposition method 'cell_vol_linear'
! > Linear deposition with deposition radius equal to cell size
!===================================================================================================================================
SUBROUTINE DepositionMethod_CV()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars             ,ONLY: wGPVol
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBPauseTime
USE MOD_LoadBalance_Vars         ,ONLY: nDeposPerElem
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Globals         ,ONLY: ElementOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iPart
INTEGER             :: i,j,k
REAL                :: Vol
REAL                :: Source(PP_nVar)
! Timers
#if USE_LOADBALANCE
REAL                :: tLBStart
#endif /*USE_LOADBALANCE*/
!==================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

DO iPart = 1,PDM%ParticleVecLength
  ! Particle not inside
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Particle has no RHS
  IF (Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER        .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TCONVERGENCE  .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_HPCONVERGENCE .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER) CYCLE

  CALL CalcPartSourceTerm(iPart,Source)
  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING,TRACING)
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(:,iPart),PEM%Element(iPart))
    END SELECT

  iElem = PEM%Element(iPart) - offsetElem
  Vol   = 0.
  PartSource_tmp = 0.

  DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
    Vol = Vol + wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
    PartSource_tmp(:,i,j,k) = Source(:)
  END DO; END DO; END DO

  PartSource_tmp            = PartSource_tmp / Vol
  PartSource(:,:,:,:,iElem) = PartSource(:,:,:,:,iElem) + PartSource_tmp

#if USE_LOADBALANCE
  ! Cell is on current proc, assign load to this cell
  IF (ElementOnProc(iElem)) nDeposPerElem(iElem-offsetElem) = nDeposPerElem(iElem-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
END DO ! iPart = 1,PDM%ParticleVecLength

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DEPOSITION,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositionMethod_CV


!===================================================================================================================================
! Deposition method 'cell_vol_linear'
! > Linear deposition with deposition radius equal to cell size
!===================================================================================================================================
SUBROUTINE DepositionMethod_CVLM()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Mesh_Vars       ,ONLY: VertexInfo_Shared,VertexVol_Shared
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
USE MOD_Particle_Vars            ,ONLY: PartNodeSource
#if USE_MPI
USE MOD_MPI_Shared               ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars          ,ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars          ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SHARED
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBPauseTime
USE MOD_LoadBalance_Vars         ,ONLY: nDeposPerElem
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Globals         ,ONLY: ElementOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iPart,iNode
INTEGER             :: FEMNodeID(   1:8)
REAL                :: PartDistDepo(1:8)
REAL                :: alpha1,alpha2,alpha3
REAL                :: Vol
REAL                :: Source(PP_nVar)
! Timers
#if USE_LOADBALANCE
REAL                :: tLBStart
#endif /*USE_LOADBALANCE*/
!==================================================================================================================================
! Nullify
#if USE_MPI
IF (myComputeNodeRank.EQ.0) FEMNodeSource_Shared(:,:) = 0.
CALL BARRIER_AND_SYNC(FEMNodeSource_Shared_Win,MPI_COMM_SHARED)
CALL MPI_WIN_UNLOCK_ALL(FEMNodeSource_Shared_Win,iError)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

!> Start an RMA exposure epoch
! MPI_WIN_POST must complete first as per https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node281.htm
! > "The call to MPI_WIN_START can block until the matching call to MPI_WIN_POST occurs at all target processes."
! No local operations prior to this epoch, so give an assertion
!> > ACTIVE SYNCHRONIZATION
! CALL MPI_WIN_FENCE(    MPI_MODE_NOPRECEDE                                                &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
!> > PASSIVE SYNCHRONIZATION
CALL MPI_WIN_LOCK_ALL(0,FEMNodeSource_Shared_Win,iError)
#else
FEMNodeSource_Shared(:,:) = 0.
#endif /*USE_MPI*/
PartNodeSource            = 0.

! Loop all particles and deposit their force contribution
DO iPart = 1,PDM%ParticleVecLength
  ! Particle not inside
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Particle has no RHS
  IF (Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER        .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TCONVERGENCE  .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_HPCONVERGENCE .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER) CYCLE

  CALL CalcPartSourceTerm(iPart,Source)
  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING,TRACING)
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(:,iPart),PEM%Element(iPart))
    END SELECT

  ! Distance to the 8 corner nodes
  alpha1          = 0.5*(PartPosRef(1,iPart)+1.0)
  alpha2          = 0.5*(PartPosRef(2,iPart)+1.0)
  alpha3          = 0.5*(PartPosRef(3,iPart)+1.0)
  PartDistDepo(1) = (1-alpha1)*(1-alpha2)*(1-alpha3)
  PartDistDepo(2) = (  alpha1)*(1-alpha2)*(1-alpha3)
  PartDistDepo(3) = (  alpha1)*  (alpha2)*(1-alpha3)
  PartDistDepo(4) = (1-alpha1)*  (alpha2)*(1-alpha3)
  PartDistDepo(5) = (1-alpha1)*(1-alpha2)*  (alpha3)
  PartDistDepo(6) = (  alpha1)*(1-alpha2)*  (alpha3)
  PartDistDepo(7) = (  alpha1)*  (alpha2)*  (alpha3)
  PartDistDepo(8) = (1-alpha1)*  (alpha2)*  (alpha3)

  iElem = PEM%Element(iPart)
  Vol   = 0.

  DO iNode = 1,8
    FEMNodeID(iNode) = VertexInfo_Shared(1,(iElem-1)*8 + iNode)
    ! FEMNodeSource(MOMV,FEMNodeID(iNode)) = FEMNodeSource(MOMV,FEMNodeID(iNode)) + Source(MOMV)*PartDistDepo(iNode)
    ! FEMNodeSource(ENER,FEMNodeID(iNode)) = FEMNodeSource(ENER,FEMNodeID(iNode)) + Source(ENER)*PartDistDepo(iNode)
    Vol = Vol + VertexVol_Shared(FEMNodeID(iNode))
    PartNodeSource(MOMV,iNode,iPart) = PartNodeSource(MOMV,iNode,iPart) + Source(MOMV)*PartDistDepo(iNode)
    PartNodeSource(ENER,iNode,iPart) = PartNodeSource(ENER,iNode,iPart) + Source(ENER)*PartDistDepo(iNode)
  END DO ! iNode = 1,8

  ! Normalize deposition with volume
  PartNodeSource(:,:,iPart) = PartNodeSource(:,:,iPart) / Vol

  ! Loop over all nodes and add normalized deposition to shared array
  DO iNode = 1,8
#if USE_MPI
    CALL MPI_ACCUMULATE(PartNodeSource(:,iNode,iPart)                                                                        ,&
                        PP_nVar, MPI_DOUBLE_PRECISION, 0,       INT(PP_nVar*(FEMNodeID(iNode)-1)*SIZE_REAL,MPI_ADDRESS_KIND) ,&
                        PP_nVar, MPI_DOUBLE_PRECISION, MPI_SUM, FEMNodeSource_Shared_Win, iError)
#else
    FEMNodeSource_Shared(:,FEMNodeID(iNode)) = FEMNodeSource_Shared(:,FEMNodeID(iNode)) + PartNodeSource(:,iNode,iPart)
#endif /*USE_MPI*/
  END DO

#if USE_LOADBALANCE
  ! Cell is on current proc, assign load to this cell
  IF (ElementOnProc(iElem)) nDeposPerElem(iElem-offsetElem) = nDeposPerElem(iElem-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
END DO ! iPart = 1,PDM%ParticleVecLength

#if USE_MPI
! Locally completes at the origin all outstanding RMA operations initiated by the calling
! process to the target process specified by rank on the specified window. For example,
! after this routine completes, the user may reuse any buffers provided to put, get, or
! accumulate operations.
! CALL MPI_WIN_FLUSH_LOCAL(0                                                                 &
!                         ,FEMNodeSource_Shared_Win                                          &
!                         ,iError)

! Finalize all RMA operation
! CALL MPI_WIN_FLUSH( 0,FEMNodeSource_Shared_Win,iError)

!> Complete the epoch - this will block until MPI is complete
!> > ACTIVE SYNCHRONIZATION
! CALL MPI_WIN_FENCE(    0                                                                 &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
! ! All done with the window - tell MPI there are no more epochs
! CALL MPI_WIN_FENCE(    MPI_MODE_NOSUCCEED                                                &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
!> > PASSIVE SYNCHRONIZATION
CALL MPI_WIN_UNLOCK_ALL(FEMNodeSource_Shared_Win,iError)

CALL MPI_WIN_LOCK_ALL(0,FEMNodeSource_Shared_Win,iError)
CALL BARRIER_AND_SYNC(FEMNodeSource_Shared_Win,MPI_COMM_SHARED)
! Communicate results between CN roots
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_IALLREDUCE(MPI_IN_PLACE,FEMNodeSource_Shared,PP_nVar*nUniqueFEMNodes,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,MPI_DEPO_REQUEST,iError)
END IF ! myComputeNodeRank.EQ.0
#endif /*USE_MPI*/

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DEPOSITION,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositionMethod_CVLM


!===================================================================================================================================
! Deposition method 'step'
! > Deposition on all cell associated with BGM cell
!===================================================================================================================================
SUBROUTINE DepositionMethod_Step()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars             ,ONLY: wGPVol
USE MOD_Mesh_Vars                ,ONLY: nElems,offsetElem
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Mesh_Vars       ,ONLY: GEO,FIBGM_nElems,FIBGM_offsetElem,FIBGM_Element
USE MOD_Particle_Mesh_Vars       ,ONLY: Elem_xGP_Shared
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBPauseTime
USE MOD_LoadBalance_Vars         ,ONLY: nDeposPerElem
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Globals         ,ONLY: ElementOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iPart
INTEGER             :: ElemID,CNElemID
INTEGER             :: i,j,k
INTEGER             :: iBGM,jBGM,kBGM
INTEGER             :: imin,jmin,kmin
INTEGER             :: imax,jmax,kmax
REAL                :: r
REAL                :: Vol
REAL                :: Source(PP_nVar)
! Timers
#if USE_LOADBALANCE
REAL                :: tLBStart
#endif /*USE_LOADBALANCE*/
!==================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! Loop all particles and deposit their force contribution
DO iPart = 1,PDM%ParticleVecLength
  ! Particle not inside
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Particle has no RHS
  IF (Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER        .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TCONVERGENCE  .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_HPCONVERGENCE .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER) CYCLE

  CALL CalcPartSourceTerm(iPart,Source)
  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING,TRACING)
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(:,iPart),PEM%Element(iPart))
    END SELECT

  ! --- get background mesh cell of point
  imin = FLOOR((  PartState(1,iPart)-r-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
  imin = MAX(GEO%FIBGMimin,imin)
  imax = CEILING((PartState(1,iPart)+r-GEO%xminglob)/GEO%FIBGMdeltas(1))
  imax = MIN(GEO%FIBGMimax,imax)
  jmin = FLOOR((  PartState(2,iPart)-r-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
  jmin = MAX(GEO%FIBGMjmin,jmin)
  jmax = CEILING((PartState(2,iPart)+r-GEO%yminglob)/GEO%FIBGMdeltas(2))
  jmax = MIN(GEO%FIBGMjmax,jmax)
  kmin = FLOOR((  PartState(3,iPart)-r-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
  kmin = MAX(GEO%FIBGMkmin,kmin)
  kmax = CEILING((PartState(3,iPart)+r-GEO%zminglob)/GEO%FIBGMdeltas(3))
  kmax = MIN(GEO%FIBGMkmax,kmax)

  ! Radius of particle
  r   = PartState(PART_DIAM,iPart)*0.5
  Vol = 0.
  PartSource_tmp = 0.

  DO kBGM = kmin,kmax; DO jBGM = jmin,jmax; DO iBGM = imin,imax
    !--- check all cells associated with this background mesh cell
    DO iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)
      ElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iElem)
      CNElemID = GetCNElemID(ElemID)
      ElemID   = ElemID - offsetElem

      ! Cell not on the local processor
      IF (ElemID.LT.1 .OR. ElemID.GT.nElems) CYCLE

      DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
        IF (NORM2(Elem_xGP_Shared(1:3,i,j,k,CNElemID)-PartState(PART_POS1:PART_POS3,iPart)) .LE. r) THEN
          ! TODO: sJ_shared
          Vol = Vol + wGPVol(i,j,k)/sJ(i,j,k,ElemID,0)
          PartSource_tmp(:,i,j,k) = Source(:)
        END IF
      END DO; END DO; END DO ! i,j,k

      PartSource_tmp             = PartSource_tmp / Vol
      PartSource(:,:,:,:,ElemID) = PartSource(:,:,:,:,ElemID) + PartSource_tmp
    END DO ! iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)

#if USE_LOADBALANCE
    ! Cell is on current proc, assign load to this cell
    IF (ElementOnProc(ElemID)) nDeposPerElem(ElemID-offsetElem) = nDeposPerElem(ElemID-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
  END DO; END DO; END DO ! iBGM,jBGM,kBGM

END DO ! iPart = 1,PDM%ParticleVecLength

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DEPOSITION,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositionMethod_Step


!===================================================================================================================================
! Deposition method 'dirac'
! > Deposition as Dirac pulse on nearest Gauss point
!===================================================================================================================================
SUBROUTINE DepositionMethod_Dirac()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars             ,ONLY: wGPVol
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Mesh_Vars                ,ONLY: Elem_xGP
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBPauseTime
USE MOD_LoadBalance_Vars         ,ONLY: nDeposPerElem
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Globals         ,ONLY: ElementOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iPart
INTEGER             :: i,j,k
INTEGER             :: ijk(3)
REAL                :: dist_tmp,dist_loc
REAL                :: Source(PP_nVar)
! Timers
#if USE_LOADBALANCE
REAL                :: tLBStart
#endif /*USE_LOADBALANCE*/
!==================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! Loop all particles and deposit their force contribution
DO iPart = 1,PDM%ParticleVecLength
  ! Particle not inside
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Particle has no RHS
  IF (Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER        .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TCONVERGENCE  .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_HPCONVERGENCE .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER) CYCLE

  CALL CalcPartSourceTerm(iPart,Source)
  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING,TRACING)
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(:,iPart),PEM%Element(iPart))
    END SELECT

  iElem = PEM%Element(iPart) - offsetElem

  ! Search for nearest Gauss point, beginning in the lower left corner
  dist_tmp = VECNORM(Elem_xGP(:,0,0,0,iElem)-PartState(1:3,iPart))
  ijk(:)   = 0
  DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
    dist_loc = VECNORM(Elem_xGP(:,i,j,k,iElem) - PartState(1:3,iPart))
    IF (dist_loc .LT. dist_tmp) THEN; ijk(:) = (/i,j,k/); dist_tmp = dist_loc; END IF
  END DO; END DO; END DO ! i,j,k

  PartSource(:,ijk(1),ijk(2),ijk(3),iElem) = PartSource(:,ijk(1),ijk(2),ijk(3),iElem) + &
                                              Source*sJ(ijk(1),ijk(2),ijk(3),iElem,0)/wGPVol(ijk(1),ijk(2),ijk(3))

#if USE_LOADBALANCE
  ! Cell is on current proc, assign load to this cell
  IF (ElementOnProc(iElem)) nDeposPerElem(iElem-offsetElem) = nDeposPerElem(iElem-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
END DO ! iPart = 1,PDM%ParticleVecLength

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DEPOSITION,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositionMethod_Dirac


!===================================================================================================================================
! Deposition method 'shapefunction_Gauss'
! > Deposition with Gaussian shape function
!===================================================================================================================================
SUBROUTINE DepositionMethod_SF_Gauss()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars       ,ONLY: VertexInfo_Shared,VertexVol_Shared
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
USE MOD_Particle_Vars            ,ONLY: PartNodeSource
#if USE_MPI
USE MOD_MPI_Shared               ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars          ,ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars          ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SHARED
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBPauseTime
USE MOD_LoadBalance_Vars         ,ONLY: nDeposPerElem
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Globals         ,ONLY: ElementOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER      :: sigma = 1.
REAL                :: kernel
REAL                :: Vol

INTEGER             :: iElem,iPart,iNode
INTEGER             :: FEMNodeID(   1:8)
REAL                :: PartDistDepo(1:8)
REAL                :: alpha(3)
REAL                :: Source(PP_nVar)
! Timers
#if USE_LOADBALANCE
REAL                :: tLBStart
#endif /*USE_LOADBALANCE*/
!==================================================================================================================================
! Nullify
#if USE_MPI
IF (myComputeNodeRank.EQ.0) FEMNodeSource_Shared(:,:) = 0.
CALL BARRIER_AND_SYNC(FEMNodeSource_Shared_Win,MPI_COMM_SHARED)
CALL MPI_WIN_UNLOCK_ALL(FEMNodeSource_Shared_Win,iError)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

!> Start an RMA exposure epoch
! MPI_WIN_POST must complete first as per https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node281.htm
! > "The call to MPI_WIN_START can block until the matching call to MPI_WIN_POST occurs at all target processes."
! No local operations prior to this epoch, so give an assertion
!> > ACTIVE SYNCHRONIZATION
! CALL MPI_WIN_FENCE(    MPI_MODE_NOPRECEDE                                                &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
!> > PASSIVE SYNCHRONIZATION
CALL MPI_WIN_LOCK_ALL(0,FEMNodeSource_Shared_Win,iError)
#else
FEMNodeSource_Shared(:,:) = 0.
#endif /*USE_MPI*/

! Loop all particles and deposit their force contribution
DO iPart = 1,PDM%ParticleVecLength
  ! Particle not inside
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Particle has no RHS
  IF (Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER        .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TCONVERGENCE  .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_HPCONVERGENCE .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER) CYCLE

  CALL CalcPartSourceTerm(iPart,Source)
  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING,TRACING)
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(:,iPart),PEM%Element(iPart))
    END SELECT

  ! Distance to the 8 corner nodes
  alpha           = 0.5*(PartPosRef(1:3,iPart)+1.0)
  PartDistDepo(1) = VECNORM((/0.,0.,0./)-alpha) ! (0,0,0)
  PartDistDepo(2) = VECNORM((/1.,0.,0./)-alpha) ! (1,0,0)
  PartDistDepo(3) = VECNORM((/1.,1.,0./)-alpha) ! (1,1,0)
  PartDistDepo(4) = VECNORM((/0.,1.,0./)-alpha) ! (0,1,0)
  PartDistDepo(5) = VECNORM((/0.,0.,1./)-alpha) ! (0,0,1)
  PartDistDepo(6) = VECNORM((/1.,0.,1./)-alpha) ! (1,0,1)
  PartDistDepo(7) = VECNORM((/1.,1.,1./)-alpha) ! (1,1,1)
  PartDistDepo(8) = VECNORM((/0.,1.,1./)-alpha) ! (0,1,1)

  Vol   = 0.
  iElem = PEM%Element(iPart)

  DO iNode = 1,8
    FEMNodeID(iNode) = VertexInfo_Shared(1,(iElem-1)*8 + iNode)
    ! FEMNodeSource(MOMV,FEMNodeID(iNode)) = FEMNodeSource(MOMV,FEMNodeID(iNode)) + Source(MOMV)*PartDistDepo(iNode)
    ! FEMNodeSource(ENER,FEMNodeID(iNode)) = FEMNodeSource(ENER,FEMNodeID(iNode)) + Source(ENER)*PartDistDepo(iNode)
    kernel = EXP(-0.5*PartDistDepo(iNode)**2/sigma**2)
    Vol    = Vol + kernel*VertexVol_Shared(FEMNodeID(iNode))
    PartNodeSource(MOMV,iNode,iPart) = PartNodeSource(MOMV,iNode,iPart) + Source(MOMV)*PartDistDepo(iNode)
    PartNodeSource(ENER,iNode,iPart) = PartNodeSource(ENER,iNode,iPart) + Source(ENER)*PartDistDepo(iNode)
  END DO ! iNode = 1,8

  ! Normalize deposition with volume
  PartNodeSource(:,:,iPart) = PartNodeSource(:,:,iPart) / Vol

  ! Loop over all nodes and add normalized deposition to shared array
  DO iNode = 1,8
    ! FEMNodeID(iNode) = VertexInfo_Shared(1,(iElem-1)*8 + iNode)
#if USE_MPI
    CALL MPI_ACCUMULATE(PartNodeSource(:,iNode,iPart)                                                                        ,&
                        PP_nVar, MPI_DOUBLE_PRECISION, 0,       INT(PP_nVar*(FEMNodeID(iNode)-1)*SIZE_REAL,MPI_ADDRESS_KIND) ,&
                        PP_nVar, MPI_DOUBLE_PRECISION, MPI_SUM, FEMNodeSource_Shared_Win, iError)

#else
    FEMNodeSource_Shared(:,FEMNodeID(iNode)) = FEMNodeSource_Shared(:,FEMNodeID(iNode)) + PartNodeSource(:,iNode,iPart)
#endif /*USE_MPI*/
  END DO

#if USE_LOADBALANCE
  ! Cell is on current proc, assign load to this cell
  IF (ElementOnProc(iElem)) nDeposPerElem(iElem-offsetElem) = nDeposPerElem(iElem-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
END DO ! iPart = 1,PDM%ParticleVecLength

#if USE_MPI
! Locally completes at the origin all outstanding RMA operations initiated by the calling
! process to the target process specified by rank on the specified window. For example,
! after this routine completes, the user may reuse any buffers provided to put, get, or
! accumulate operations.
! CALL MPI_WIN_FLUSH_LOCAL(0                                                                 &
!                         ,FEMNodeSource_Shared_Win                                          &
!                         ,iError)

! Finalize all RMA operation
! CALL MPI_WIN_FLUSH( 0,FEMNodeSource_Shared_Win,iError)

!> Complete the epoch - this will block until MPI is complete
!> > ACTIVE SYNCHRONIZATION
! CALL MPI_WIN_FENCE(    0                                                                 &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
! ! All done with the window - tell MPI there are no more epochs
! CALL MPI_WIN_FENCE(    MPI_MODE_NOSUCCEED                                                &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
!> > PASSIVE SYNCHRONIZATION
CALL MPI_WIN_UNLOCK_ALL(FEMNodeSource_Shared_Win,iError)

CALL MPI_WIN_LOCK_ALL(0,FEMNodeSource_Shared_Win,iError)
CALL BARRIER_AND_SYNC(FEMNodeSource_Shared_Win,MPI_COMM_SHARED)
! Communicate results between CN roots
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_IALLREDUCE(MPI_IN_PLACE,FEMNodeSource_Shared,PP_nVar*nUniqueFEMNodes,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,MPI_DEPO_REQUEST,iError)
END IF ! myComputeNodeRank.EQ.0
#endif /*USE_MPI*/

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DEPOSITION,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositionMethod_SF_Gauss


!===================================================================================================================================
! Deposition method 'shapefunction_Gauss'
! > Deposition with Gaussian shape function, respecting PP_N
!===================================================================================================================================
SUBROUTINE DepositionMethod_SF_Poly()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars       ,ONLY: VertexInfo_Shared,VertexVol_Shared
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
USE MOD_Particle_Vars            ,ONLY: PartNodeSource
#if USE_MPI
USE MOD_MPI_Shared               ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars          ,ONLY: myComputeNodeRank
USE MOD_MPI_Shared_Vars          ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SHARED
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBPauseTime
USE MOD_LoadBalance_Vars         ,ONLY: nDeposPerElem
USE MOD_Mesh_Vars                ,ONLY: offsetElem
USE MOD_Particle_Globals         ,ONLY: ElementOnProc
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: kernel,sigma
REAL                :: Vol

INTEGER             :: iElem,iPart,iNode
INTEGER             :: FEMNodeID(   1:8)
REAL                :: PartDistDepo(1:8)
REAL                :: alpha(3)
REAL                :: Source(PP_nVar)
! Timers
#if USE_LOADBALANCE
REAL                :: tLBStart
#endif /*USE_LOADBALANCE*/
!==================================================================================================================================
! Nullify
#if USE_MPI
IF (myComputeNodeRank.EQ.0) FEMNodeSource_Shared(:,:) = 0.
CALL BARRIER_AND_SYNC(FEMNodeSource_Shared_Win,MPI_COMM_SHARED)
CALL MPI_WIN_UNLOCK_ALL(FEMNodeSource_Shared_Win,iError)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

!> Start an RMA exposure epoch
! MPI_WIN_POST must complete first as per https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node281.htm
! > "The call to MPI_WIN_START can block until the matching call to MPI_WIN_POST occurs at all target processes."
! No local operations prior to this epoch, so give an assertion
!> > ACTIVE SYNCHRONIZATION
! CALL MPI_WIN_FENCE(    MPI_MODE_NOPRECEDE                                                &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
!> > PASSIVE SYNCHRONIZATION
CALL MPI_WIN_LOCK_ALL(0,FEMNodeSource_Shared_Win,iError)
#else
FEMNodeSource_Shared(:,:) = 0.
#endif /*USE_MPI*/

sigma = PP_N+1

! Loop all particles and deposit their force contribution
DO iPart = 1,PDM%ParticleVecLength
  ! Particle not inside
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE

  ! Particle has no RHS
  IF (Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER        .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TCONVERGENCE  .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_HPCONVERGENCE .OR. &
      Species(PartSpecies(iPart))%RHSMethod .EQ. RHS_TRACER) CYCLE

  CALL CalcPartSourceTerm(iPart,Source)
  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING,TRACING)
      CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(:,iPart),PEM%Element(iPart))
    END SELECT

  ! Distance to the 8 corner nodes
  alpha           = 0.5*(PartPosRef(1:3,iPart)+1.0)
  PartDistDepo(1) = VECNORM((/0.,0.,0./)-alpha) ! (0,0,0)
  PartDistDepo(2) = VECNORM((/1.,0.,0./)-alpha) ! (1,0,0)
  PartDistDepo(3) = VECNORM((/1.,1.,0./)-alpha) ! (1,1,0)
  PartDistDepo(4) = VECNORM((/0.,1.,0./)-alpha) ! (0,1,0)
  PartDistDepo(5) = VECNORM((/0.,0.,1./)-alpha) ! (0,0,1)
  PartDistDepo(6) = VECNORM((/1.,0.,1./)-alpha) ! (1,0,1)
  PartDistDepo(7) = VECNORM((/1.,1.,1./)-alpha) ! (1,1,1)
  PartDistDepo(8) = VECNORM((/0.,1.,1./)-alpha) ! (0,1,1)

  Vol   = 0.
  iElem = PEM%Element(iPart)

  DO iNode = 1,8
    FEMNodeID(iNode) = VertexInfo_Shared(1,(iElem-1)*8 + iNode)
    ! FEMNodeSource(MOMV,FEMNodeID(iNode)) = FEMNodeSource(MOMV,FEMNodeID(iNode)) + Source(MOMV)*PartDistDepo(iNode)
    ! FEMNodeSource(ENER,FEMNodeID(iNode)) = FEMNodeSource(ENER,FEMNodeID(iNode)) + Source(ENER)*PartDistDepo(iNode)
    kernel = EXP(-0.5*PartDistDepo(iNode)**2/sigma**2)
    Vol    = Vol + kernel*VertexVol_Shared(FEMNodeID(iNode))
    Vol = Vol + VertexVol_Shared(FEMNodeID(iNode))
    PartNodeSource(MOMV,iNode,iPart) = Source(MOMV)*PartDistDepo(iNode)
    PartNodeSource(ENER,iNode,iPart) = Source(ENER)*PartDistDepo(iNode)
  END DO ! iNode = 1,8

  ! Normalize deposition with volume
  PartNodeSource(:,:,iPart) = PartNodeSource(:,:,iPart) / Vol

  ! Loop over all nodes and add normalized deposition to shared array
  DO iNode = 1,8
    ! FEMNodeID(iNode) = VertexInfo_Shared(1,(iElem-1)*8 + iNode)
#if USE_MPI
    CALL MPI_ACCUMULATE(PartNodeSource(:,iNode,iPart)                                                                        ,&
                        PP_nVar, MPI_DOUBLE_PRECISION, 0,       INT(PP_nVar*(FEMNodeID(iNode)-1)*SIZE_REAL,MPI_ADDRESS_KIND) ,&
                        PP_nVar, MPI_DOUBLE_PRECISION, MPI_SUM, FEMNodeSource_Shared_Win, iError)

#else
    FEMNodeSource_Shared(:,FEMNodeID(iNode)) = FEMNodeSource_Shared(:,FEMNodeID(iNode)) + PartNodeSource(:,iNode,iPart)
#endif /*USE_MPI*/
  END DO

#if USE_LOADBALANCE
  ! Cell is on current proc, assign load to this cell
  IF (ElementOnProc(iElem)) nDeposPerElem(iElem-offsetElem) = nDeposPerElem(iElem-offsetElem) + 1
#endif /*USE_LOADBALANCE*/
END DO ! iPart = 1,PDM%ParticleVecLength

#if USE_MPI
! Locally completes at the origin all outstanding RMA operations initiated by the calling
! process to the target process specified by rank on the specified window. For example,
! after this routine completes, the user may reuse any buffers provided to put, get, or
! accumulate operations.
! CALL MPI_WIN_FLUSH_LOCAL(0                                                                 &
!                         ,FEMNodeSource_Shared_Win                                          &
!                         ,iError)

! Finalize all RMA operation
! CALL MPI_WIN_FLUSH( 0,FEMNodeSource_Shared_Win,iError)

!> Complete the epoch - this will block until MPI is complete
!> > ACTIVE SYNCHRONIZATION
! CALL MPI_WIN_FENCE(    0                                                                 &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
! ! All done with the window - tell MPI there are no more epochs
! CALL MPI_WIN_FENCE(    MPI_MODE_NOSUCCEED                                                &
!                   ,    FEMNodeSource_Shared_Win                                          &
!                   ,    iError)
!> > PASSIVE SYNCHRONIZATION
CALL MPI_WIN_UNLOCK_ALL(FEMNodeSource_Shared_Win,iError)

!> Return to shared memory paradigm
CALL MPI_WIN_LOCK_ALL(0,FEMNodeSource_Shared_Win,iError)
CALL BARRIER_AND_SYNC(FEMNodeSource_Shared_Win,MPI_COMM_SHARED)
! Communicate results between CN roots
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_IALLREDUCE(MPI_IN_PLACE,FEMNodeSource_Shared,PP_nVar*nUniqueFEMNodes,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,MPI_DEPO_REQUEST,iError)
END IF ! myComputeNodeRank.EQ.0
#endif /*USE_MPI*/

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DEPOSITION,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositionMethod_SF_Poly


!==================================================================================================================================
!> Compute source terms for particles
!==================================================================================================================================
SUBROUTINE CalcPartSourceTerm(PartID,Source)
! MODULES
USE MOD_Particle_Vars       ,ONLY: Species,PartSpecies,PartState,Pt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID
REAL,INTENT(INOUT)  :: Source(PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

Source       = 0

! Calculate particle force
Source(MOMV) = - Pt(1:3,PartID)*Species(PartSpecies(PartID))%MassIC
! Calculate the work
Source(ENER) = DOT_PRODUCT(Source(MOMV),PartState(PART_VELV,PartID))

#if USE_PARTTEMP
! Calculate particle heat exchange (if activated)
Source(ENER) = Source(ENER) - Pt(4,PartID) * Species(PartSpecies(PartID))%MassIC * Species(PartSpecies(PartID))%SpecificHeatIC
#endif

END SUBROUTINE CalcPartSourceTerm

#endif /*PARTICLES_COUPLING >= 2*/

END MODULE MOD_Particle_Deposition_Method
