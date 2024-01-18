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

INTERFACE InitDepositionMethod
  MODULE PROCEDURE InitDepositionMethod
END INTERFACE

PUBLIC :: InitDepositionMethod
!==================================================================================================================================

CONTAINS

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
SELECT CASE(TRIM(DepositionType))
  CASE('cell_vol')
    DepositionMethod => DepositionMethod_CV
  CASE('cell_vol_linear')
    DepositionMethod => DepositionMethod_CVLM
  CASE('step')
    DepositionMethod => DepositionMethod_Step
  CASE('dirac')
    DepositionMethod => DepositionMethod_Dirac
  CASE('shapefunc_gauss')
    DepositionMethod => DepositionMethod_SF_Gauss
  CASE('shapefunc_poly')
    DepositionMethod => DepositionMethod_SF_Poly
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'Unknown DepositionMethod! '//TRIM(DepositionType))
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
!==================================================================================================================================

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

END DO ! iPart = 1,PDM%ParticleVecLength

END SUBROUTINE DepositionMethod_CV


!===================================================================================================================================
! Deposition method 'cell_vol_linear'
! > Linear deposition with deposition radius equal to cell size
!===================================================================================================================================
SUBROUTINE DepositionMethod_CVLM()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                ,ONLY: nElems,offsetElem
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Mesh_Vars       ,ONLY: VertexInfo_Shared
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iPart,iNode
INTEGER             :: i,j,k
INTEGER             :: FEMNodeID(   1:8)
REAL                :: PartDistDepo(1:8)
REAL                :: alpha1,alpha2,alpha3
REAL                :: Source(PP_nVar)
!==================================================================================================================================

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

  DO iNode = 1,8
    FEMNodeID(iNode) = VertexInfo_Shared(1,(PEM%Element(iPart)-1)*8 + iNode)
    FEMNodeSource(MOMV,FEMNodeID(iNode)) = FEMNodeSource(MOMV,FEMNodeID(iNode)) + Source(MOMV)*PartDistDepo(iNode)
    FEMNodeSource(ENER,FEMNodeID(iNode)) = FEMNodeSource(ENER,FEMNodeID(iNode)) + Source(ENER)*PartDistDepo(iNode)
  END DO ! iNode = 1,8

! #if USE_LOADBALANCE
!   CALL LBElemSplitTime(PEM%LocalElemID(iPart),tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to iElem
! #endif /*USE_LOADBALANCE*/
END DO ! iPart = 1,PDM%ParticleVecLength

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,FEMNodeSource,nUniqueFEMNodes,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

! Interpolate node source values to volume polynomial
DO iElem = 1,nElems

  ! Fill the NodeID array
  DO iNode = 1,8
    FEMNodeID(iNode) = VertexInfo_Shared(1,(offsetElem+iElem-1)*8 + iNode)
  END DO

  DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
    alpha1 = CellVolWeightFac(i)
    alpha2 = CellVolWeightFac(j)
    alpha3 = CellVolWeightFac(k)
    PartSource(1:PP_nVar,i,j,k,iElem) = FEMNodeSource(1:PP_nVar,FEMNodeID(1)) * (1-alpha1)*(1-alpha2)*(1-alpha3) + &
                                        FEMNodeSource(1:PP_nVar,FEMNodeID(2)) * (  alpha1)*(1-alpha2)*(1-alpha3) + &
                                        FEMNodeSource(1:PP_nVar,FEMNodeID(4)) * (1-alpha1)*  (alpha2)*(1-alpha3) + &
                                        FEMNodeSource(1:PP_nVar,FEMNodeID(3)) * (  alpha1)*  (alpha2)*(1-alpha3) + &
                                        FEMNodeSource(1:PP_nVar,FEMNodeID(5)) * (1-alpha1)*(1-alpha2)*  (alpha3) + &
                                        FEMNodeSource(1:PP_nVar,FEMNodeID(6)) * (  alpha1)*(1-alpha2)*  (alpha3) + &
                                        FEMNodeSource(1:PP_nVar,FEMNodeID(8)) * (1-alpha1)*  (alpha2)*  (alpha3) + &
                                        FEMNodeSource(1:PP_nVar,FEMNodeID(7)) * (  alpha1)*  (alpha2)*  (alpha3)

  END DO; END DO; END DO
END DO ! iElem

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
!==================================================================================================================================

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
  r = PartState(PART_DIAM,iPart)*0.5

  DO kBGM = kmin,kmax; DO jBGM = jmin,jmax; DO iBGM = imin,imax
    !--- check all cells associated with this background mesh cell
    DO iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)
      ElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iElem)
      CNElemID = GetCNElemID(ElemID)
      ElemID   = ElemID - offsetElem

      ! Cell not on the local processor
      IF (ElemID.LT.1 .OR. ElemID.GT.nElems) CYCLE

      PartSource_tmp = 0.

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
  END DO; END DO; END DO ! iBGM,jBGM,kBGM

END DO ! iPart = 1,PDM%ParticleVecLength

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
REAL                :: Vol
REAL                :: Source(PP_nVar)
!==================================================================================================================================

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
  ijk(:)   =           0
  DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
    dist_loc = VECNORM(Elem_xGP(:,i,j,k,iElem) - PartState(1:3,iPart))
    IF (dist_loc .LT. dist_tmp) THEN; ijk(:) = (/i,j,k/); dist_tmp = dist_loc; END IF
  END DO; END DO; END DO ! i,j,k

  Vol = Vol + wGPVol(ijk(1),ijk(2),ijk(3))/sJ(ijk(1),ijk(2),ijk(3),iElem,0)
  PartSource_tmp(:,i,j,k) = Source(:)

  PartSource_tmp            = PartSource_tmp / Vol
  PartSource(:,:,:,:,iElem) = PartSource(:,:,:,:,iElem) + PartSource_tmp

END DO ! iPart = 1,PDM%ParticleVecLength

END SUBROUTINE DepositionMethod_Dirac


!===================================================================================================================================
! Deposition method 'shapefunction_Gauss'
! > Deposition with Gaussian shape function
!===================================================================================================================================
SUBROUTINE DepositionMethod_SF_Gauss()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars             ,ONLY: wGPVol
USE MOD_Mesh_Vars                ,ONLY: nElems,offsetElem
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars       ,ONLY: GEO,FIBGM_nElems,FIBGM_offsetElem,FIBGM_Element
USE MOD_Particle_Mesh_Vars       ,ONLY: Elem_xGP_Shared
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER      :: sigma = 1.
INTEGER             :: iElem,iPart
INTEGER             :: ElemID,CNElemID
INTEGER             :: i,j,k
INTEGER             :: iBGM,jBGM,kBGM
INTEGER             :: imin,jmin,kmin
INTEGER             :: imax,jmax,kmax
REAL                :: dist_tmp
REAL                :: kernel
REAL                :: r
REAL                :: Vol
REAL                :: Source(PP_nVar)
!==================================================================================================================================

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
  r = PartState(PART_DIAM,iPart)*0.5

  DO kBGM = kmin,kmax; DO jBGM = jmin,jmax; DO iBGM = imin,imax
    !--- check all cells associated with this background mesh cell
    DO iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)
      ElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iElem)
      CNElemID = GetCNElemID(ElemID)
      ElemID   = ElemID - offsetElem

      ! Cell not on the local processor
      IF (ElemID.LT.1 .OR. ElemID.GT.nElems) CYCLE

      PartSource_tmp = 0.

      DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
        ! Check if point is within the particle radius
        dist_tmp = VECNORM(Elem_xGP_Shared(:,i,j,k,CNElemID)-PartState(1:3,iPart))

        IF (dist_tmp .LE. r) THEN
          ! TODO: sJ_shared
          kernel = EXP(-0.5*dist_tmp**2/sigma**2)
          Vol    = Vol + kernel*wGPVol(i,j,k)/sJ(i,j,k,ElemID,0)
          !kernel = 1./(2*PP_PI*sigma**2) * EXP(-0.5*dist_tmp**2/sigma**2)
          PartSource_tmp(:,i,j,k) = Source(:)
        END IF
      END DO; END DO; END DO ! i,j,k

      PartSource_tmp             = PartSource_tmp / Vol
      PartSource(:,:,:,:,ElemID) = PartSource(:,:,:,:,ElemID) + PartSource_tmp
    END DO ! iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)
  END DO; END DO; END DO ! iBGM,jBGM,kBGM

END DO ! iPart = 1,PDM%ParticleVecLength

END SUBROUTINE DepositionMethod_SF_Gauss


!===================================================================================================================================
! Deposition method 'shapefunction_Gauss'
! > Deposition with Gaussian shape function, respecting PP_N
!===================================================================================================================================
SUBROUTINE DepositionMethod_SF_Poly()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars             ,ONLY: wGPVol
USE MOD_Mesh_Vars                ,ONLY: nElems,offsetElem
USE MOD_Mesh_Vars                ,ONLY: sJ
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
USE MOD_Particle_Deposition_Vars
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Mesh_Vars       ,ONLY: GEO,FIBGM_nElems,FIBGM_offsetElem,FIBGM_Element
USE MOD_Particle_Mesh_Vars       ,ONLY: Elem_xGP_Shared
USE MOD_Particle_Mesh_Tools      ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Vars            ,ONLY: PDM,PEM,PartPosRef,PartState,Species,PartSpecies
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
REAL                :: dist_tmp
REAL                :: kernel
REAL                :: sigma
REAL                :: r
REAL                :: Vol
REAL                :: Source(PP_nVar)
!==================================================================================================================================

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
  r = PartState(PART_DIAM,iPart)*0.5

  DO kBGM = kmin,kmax; DO jBGM = jmin,jmax; DO iBGM = imin,imax
    !--- check all cells associated with this background mesh cell
    DO iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)
      ElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iElem)
      CNElemID = GetCNElemID(ElemID)
      ElemID   = ElemID - offsetElem

      ! Cell not on the local processor
      IF (ElemID.LT.1 .OR. ElemID.GT.nElems) CYCLE

      PartSource_tmp = 0.

      DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
        ! Check if point is within the particle radius
        dist_tmp = VECNORM(Elem_xGP_Shared(:,i,j,k,CNElemID)-PartState(1:3,iPart))

        IF (dist_tmp .LE. r) THEN
          ! TODO: sJ_shared
          kernel = EXP(-0.5*dist_tmp**2/sigma**2)
          Vol    = Vol + kernel*wGPVol(i,j,k)/sJ(i,j,k,ElemID,0)
          !kernel = 1./(2*PP_PI*sigma**2) * EXP(-0.5*dist_tmp**2/sigma**2)
          PartSource_tmp(:,i,j,k) = Source(:)
        END IF
      END DO; END DO; END DO ! i,j,k

      PartSource_tmp             = PartSource_tmp / Vol
      PartSource(:,:,:,:,ElemID) = PartSource(:,:,:,:,ElemID) + PartSource_tmp
    END DO ! iElem = 1, FIBGM_nElems(iBGM,jBGM,kBGM)
  END DO; END DO; END DO ! iBGM,jBGM,kBGM

END DO ! iPart = 1,PDM%ParticleVecLength

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

END MODULE MOD_Particle_Deposition_Method
