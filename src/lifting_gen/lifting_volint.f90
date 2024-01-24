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
#if PARABOLIC
#include "flexi.h"
#include "eos.h"

MODULE MOD_Lifting_VolInt_gen
IMPLICIT NONE
PRIVATE

INTERFACE Lifting_VolInt_gen
   MODULE PROCEDURE Lifting_VolInt_Nonconservative
   MODULE PROCEDURE Lifting_VolInt_Conservative
END INTERFACE

PUBLIC::Lifting_VolInt_gen
!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> \brief Computes the volume integral of the BR1 or BR2 scheme in the non conservative way for all directions in strong form.
!>
!> Requires lifting in strong form
!> In the non conservative form of the volume integral in BR1/2 we first differentiate the flux (which is the solution in BR1) and
!> then apply the metric terms. This is the fastest implementation of the volume integral and only available in strong form.
!==================================================================================================================================
PPURE SUBROUTINE Lifting_VolInt_Nonconservative(nVarIn,nVarOut,UPrim,gradUx,gradUy,gradUz)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: D_T
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde
#if (PP_dim==3)
USE MOD_Mesh_Vars    ,ONLY: Metrics_hTilde
#endif
USE MOD_Mesh_Vars    ,ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
USE MOD_FV_Vars      ,ONLY: FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ
#if (PP_dim==3)
USE MOD_FV_Vars      ,ONLY: FV_Metrics_hTilde_sJ
#endif
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                              :: nVarIn
INTEGER,INTENT(IN)                              :: nVarOut
REAL,INTENT(IN)                                 :: UPrim( nVarIn ,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution
REAL,INTENT(OUT)                                :: gradUx(nVarOut,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< gradients in x-direction
REAL,INTENT(OUT)                                :: gradUy(nVarOut,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< gradients in y-direction
REAL,INTENT(OUT)                                :: gradUz(nVarOut,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< gradients in z-direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(nVarOut)                         :: gradUxi,gradUeta
#if (PP_dim==3)
REAL,DIMENSION(nVarOut)                         :: gradUzeta
#endif
INTEGER                                         :: iElem,i,j,k,l
#if FV_ENABLED
REAL,DIMENSION(nVarOut,0:PP_N,0:PP_N,0:PP_NZ)   :: gradUxi_central
REAL,DIMENSION(nVarOut,0:PP_N,0:PP_N,0:PP_NZ)   :: gradUeta_central
#if PP_dim == 3
REAL,DIMENSION(nVarOut,0:PP_N,0:PP_N,0:PP_NZ)   :: gradUzeta_central
#endif
#endif
!==================================================================================================================================
! volume integral
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.0) THEN ! DG element
#endif
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      gradUxi     = D_T(0,i)*UPrim(:,0,j,k,iElem)
      gradUeta    = D_T(0,j)*UPrim(:,i,0,k,iElem)
#if (PP_dim==3)
      gradUzeta   = D_T(0,k)*UPrim(:,i,j,0,iElem)
#endif
      DO l=1,PP_N
        gradUxi   = gradUxi   + D_T(l,i)*UPrim(:,l,j,k,iElem)
        gradUeta  = gradUeta  + D_T(l,j)*UPrim(:,i,l,k,iElem)
#if (PP_dim==3)
        gradUzeta = gradUzeta + D_T(l,k)*UPrim(:,i,j,l,iElem)
#endif
      END DO
#if (PP_dim==3)
    gradUx(:,i,j,k,iElem) = Metrics_fTilde(1,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(1,i,j,k,iElem,0)*gradUeta  &
                          + Metrics_hTilde(1,i,j,k,iElem,0)*gradUzeta
    gradUy(:,i,j,k,iElem) = Metrics_fTilde(2,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(2,i,j,k,iElem,0)*gradUeta  &
                          + Metrics_hTilde(2,i,j,k,iElem,0)*gradUzeta
    gradUz(:,i,j,k,iElem) = Metrics_fTilde(3,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(3,i,j,k,iElem,0)*gradUeta  &
                          + Metrics_hTilde(3,i,j,k,iElem,0)*gradUzeta
#else
    gradUx(:,i,j,k,iElem) = Metrics_fTilde(1,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(1,i,j,k,iElem,0)*gradUeta
    gradUy(:,i,j,k,iElem) = Metrics_fTilde(2,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(2,i,j,k,iElem,0)*gradUeta
#endif /*PP_dim==3*/
   END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
  ELSE
    CALL FV_CalcGradients(nVarIn,nVarOut,iElem,UPrim(:,:,:,:,iElem),gradUxi_central,gradUeta_central&
#if (PP_dim==3)
                                                                   ,gradUzeta_central&
#endif /*PP_dim==3*/
                                                                   )
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
#if (PP_dim==3)
      gradUx(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ(1,i,j,k,iElem)*gradUxi_central  (:,i,j,k) &
                            + FV_Metrics_gTilde_sJ(1,i,j,k,iElem)*gradUeta_central (:,i,j,k) &
                            + FV_Metrics_hTilde_sJ(1,i,j,k,iElem)*gradUzeta_central(:,i,j,k)
      gradUy(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ(2,i,j,k,iElem)*gradUxi_central  (:,i,j,k) &
                            + FV_Metrics_gTilde_sJ(2,i,j,k,iElem)*gradUeta_central (:,i,j,k) &
                            + FV_Metrics_hTilde_sJ(2,i,j,k,iElem)*gradUzeta_central(:,i,j,k)
      gradUz(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ(3,i,j,k,iElem)*gradUxi_central  (:,i,j,k) &
                            + FV_Metrics_gTilde_sJ(3,i,j,k,iElem)*gradUeta_central (:,i,j,k) &
                            + FV_Metrics_hTilde_sJ(3,i,j,k,iElem)*gradUzeta_central(:,i,j,k)
#else
      gradUx(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ(1,i,j,k,iElem)*gradUxi_central  (:,i,j,k) &
                            + FV_Metrics_gTilde_sJ(1,i,j,k,iElem)*gradUeta_central (:,i,j,k)
      gradUy(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ(2,i,j,k,iElem)*gradUxi_central  (:,i,j,k) &
                            + FV_Metrics_gTilde_sJ(2,i,j,k,iElem)*gradUeta_central (:,i,j,k)
#endif
   END DO; END DO; END DO! i,j,k=0,PP_N
  END IF
#endif /*FV_ENABLED*/
END DO ! iElem=1,nElems

END SUBROUTINE Lifting_VolInt_Nonconservative

!==================================================================================================================================
!> \brief Computes the volume integral of the BR1 scheme in the conservative way for one direction at a time (x,y,z) in weak or
!> strong form (BR1) or only in strong form (BR2)
!>
!> In the conservative form of the volume integral in BR1 we calculate the transformed solution using the Lifting_Metrics routine
!> and then integrate over the derivative of this transformed solution. A weak or strong form volume integral can be performed.
!> - Weak form: We integrate over the transformed flux times the derivative of the trial functions, which is implemented
!>   using the D_hat_T matrix. This is the transpose of the D_hat matrix which is defined as
!>   \f$ \hat{D}_{ij}=-\frac{\omega_i}{\omega_j} D_{ij} \f$ where \f$ D_{ij} \f$ is the standard polynomial derivative matrix.
!> - Strong form: We integrate over the derivative of the transformed flux times the trial function. The derivative of the
!>   flux is calculated using the transpose of the polynomial derivative matrix D_T.
!>
!> For the implementation this means we only have to decide between using the D_hat_T or D_T matrix in the volume integral at the
!> the beginning of the routine to choose between weak or strong form.
!==================================================================================================================================
SUBROUTINE Lifting_VolInt_Conservative(nVarIn,nVarOut,dir,UPrim,gradU)
! MODULES
USE MOD_PreProc
USE MOD_Lifting_Vars       ,ONLY: doWeakLifting
USE MOD_DG_Vars            ,ONLY: D_hat_T,D_T
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Mesh_Vars          ,ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars            ,ONLY: FV_Elems
USE MOD_FV_Vars            ,ONLY: FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ
USE MOD_FV_Vars            ,ONLY: gradUxi_central, gradUeta_central
#if (PP_dim==3)
USE MOD_FV_Vars            ,ONLY: FV_Metrics_hTilde_sJ
USE MOD_FV_Vars            ,ONLY: gradUzeta_central
#endif
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nVarIn
INTEGER,INTENT(IN)            :: nVarOut
INTEGER,INTENT(IN)            :: dir       !< direction (x,y,z)
REAL,INTENT(IN)               :: UPrim(nVarIn ,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution
REAL,INTENT(OUT)              :: gradU(nVarOut,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution gradient in direction dir
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                          :: DMat(0:PP_N,0:PP_N)
REAL,DIMENSION(nVarOut,0:PP_N,0:PP_N,0:PP_NZ) :: UE_f,UE_g,UE_h ! transformed gradient flux (i.e. transformed solution)
INTEGER                                       :: iElem,i,j,k,l
!==================================================================================================================================
IF(doWeakLifting)THEN
  DMat=D_hat_T
ELSE
  DMat=D_T
END IF

! volume integral
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.0) THEN ! DG element
#endif
  CALL Lifting_Metrics(nVarIn,nVarOut,dir,UPrim(:,:,:,:,iElem),&
                       Metrics_fTilde(:,:,:,:,iElem,0),&
                       Metrics_gTilde(:,:,:,:,iElem,0),&
                       Metrics_hTilde(:,:,:,:,iElem,0),&
                       UE_f,UE_g,UE_h)

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    gradU(:,i,j,k,iElem)   =                      DMat(0,i)*UE_f(:,0,j,k)+&
#if (PP_dim==3)
                                                  DMat(0,k)*UE_h(:,i,j,0)+&
#endif
                                                  DMat(0,j)*UE_g(:,i,0,k)
    DO l=1,PP_N
      gradU(:,i,j,k,iElem) = gradU(:,i,j,k,iElem)+DMat(l,i)*UE_f(:,l,j,k)+&
#if (PP_dim==3)
                                                  DMat(l,k)*UE_h(:,i,j,l)+&
#endif
                                                  DMat(l,j)*UE_g(:,i,l,k)
    END DO ! l
  END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
  ELSE
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      gradU(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ(dir,i,j,k,iElem)*gradUxi_central  (:,i,j,k,iElem) &
#if (PP_dim==3)
                           + FV_Metrics_hTilde_sJ(dir,i,j,k,iElem)*gradUzeta_central(:,i,j,k,iElem) &
#endif
                           + FV_Metrics_gTilde_sJ(dir,i,j,k,iElem)*gradUeta_central (:,i,j,k,iElem)
    END DO; END DO; END DO! i,j,k=0,PP_N
  END IF
#endif /*FV_ENABLED*/

END DO ! iElem=1,nElems
END SUBROUTINE Lifting_VolInt_Conservative

!==================================================================================================================================
!> \brief Compute the tranformed gradient fluxes
!>
!> Transform the gradient terms by multiplying them with the metrics terms
!> For the direction \f$ d \f$ the transformed gradient flux is \f$ \sum_{n=1}^3 Ja^d_n U \f$.
!==================================================================================================================================
SUBROUTINE Lifting_Metrics(nVarIn,nVarOut,dir,UPrim,Mf,Mg,Mh,UPrim_f,UPrim_g,UPrim_h)
! MODULES
USE MOD_DG_Vars,ONLY:nDOFElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarIn
INTEGER,INTENT(IN) :: nVarOut
INTEGER,INTENT(IN) :: dir                          !< direction (x,y,z)
REAL,INTENT(IN)    :: Mf(3,nDOFElem)               !< metrics in xi
REAL,INTENT(IN)    :: Mg(3,nDOFElem)               !< metrics in eta
REAL,INTENT(IN)    :: Mh(3,nDOFElem)               !< metrics in zeta
REAL,INTENT(IN)    :: UPrim  (nVarIn ,nDOFElem)    !< solution ("flux")
REAL,INTENT(OUT)   :: UPrim_f(nVarOut,nDOFElem)    !< gradient flux xi
REAL,INTENT(OUT)   :: UPrim_g(nVarOut,nDOFElem)    !< gradient flux eta
REAL,INTENT(OUT)   :: UPrim_h(nVarOut,nDOFElem)    !< gradient flux zeta
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i
!==================================================================================================================================
DO i=1,nDOFElem
  UPrim_f(:,i) = Mf(dir,i)*UPrim(:,i)
  UPrim_g(:,i) = Mg(dir,i)*UPrim(:,i)
#if (PP_dim==3)
  UPrim_h(:,i) = Mh(dir,i)*UPrim(:,i)
#endif
END DO ! i
END SUBROUTINE Lifting_Metrics

#if FV_RECONSTRUCT
!==================================================================================================================================
!> Calculate slopes in the inside and limit them using TVD limiters.
!> It is important to use physical distances (and not reference distances) to calculate the slope, since otherwise the
!> slope limiter can not be applied. (Scenario: inner-cell-stretching)
!> Additionally build central limited slopes for the computation of gradients used for the viscous fluxes.
!==================================================================================================================================
PPURE SUBROUTINE FV_CalcGradients(nVarIn,nVarOut,iElem,UPrim,gradUxi_central,gradUeta_central&
#if PP_dim == 3
  ,gradUzeta_central&
#endif
  )
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars        ,ONLY: FV_sdx_XI,FV_sdx_ETA
#if PP_dim == 3
USE MOD_FV_Vars        ,ONLY: FV_sdx_ZETA
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarIn
INTEGER,INTENT(IN) :: nVarOut
INTEGER,INTENT(IN) :: iElem
REAL,INTENT(IN)    :: UPrim            (nVarIn ,0:PP_N,0:PP_N,0:PP_NZ) !< primitive volume solution
REAL,INTENT(OUT)   :: gradUxi_central  (nVarOut,0:PP_N,0:PP_N,0:PP_NZ) !< physical slope in   xi-direction (mean value)
REAL,INTENT(OUT)   :: gradUeta_central (nVarOut,0:PP_N,0:PP_N,0:PP_NZ) !< physical slope in  eta-direction (mean value)
#if PP_dim == 3
REAL,INTENT(OUT)   :: gradUzeta_central(nVarOut,0:PP_N,0:PP_N,0:PP_NZ) !< physical slope in zeta-direction (mean value)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(nVarIn,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUxi_tmp
REAL,DIMENSION(nVarIn,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUeta_tmp
#if PP_dim == 3
REAL,DIMENSION(nVarIn,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUzeta_tmp
#endif
INTEGER                                        :: l,iVar
INTEGER                                        :: i,j,k
!==================================================================================================================================
! Strategy:
! 1. compute slopes between subcells in all 3 directions and store in temporary arrays
!    Attention concerning the storage order: last index corresponds to the respective direction (xi/eta/zeta)!
! 2. limit the slopes from the temporary array and store in gradUxi/eta/zeta

! 1. gradients of inner subcells
gradUxi_tmp   = 0.
gradUeta_tmp  = 0.
#if PP_dim == 3
gradUzeta_tmp = 0.
#endif
DO l=1,PP_N
  DO iVar=1,nVarIn
    gradUxi_tmp  (iVar,:,:,l) = (UPrim(iVar,l,:,:) - UPrim(iVar,l-1,:,:)) * FV_sdx_XI  (:,:,l,iElem)
    gradUeta_tmp (iVar,:,:,l) = (UPrim(iVar,:,l,:) - UPrim(iVar,:,l-1,:)) * FV_sdx_ETA (:,:,l,iElem)
#if PP_dim == 3
    gradUzeta_tmp(iVar,:,:,l) = (UPrim(iVar,:,:,l) - UPrim(iVar,:,:,l-1)) * FV_sdx_ZETA(:,:,l,iElem)
#endif
  END DO
END DO

! 2. limit with central limiter for viscous fluxes
DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  gradUxi_central  (:,i,j,k) = 0.5*(gradUxi_tmp  (:,j,k,i)+gradUxi_tmp  (:,j,k,i+1))
  gradUeta_central (:,i,j,k) = 0.5*(gradUeta_tmp (:,i,k,j)+gradUeta_tmp (:,i,k,j+1))
#if PP_dim == 3
  gradUzeta_central(:,i,j,k) = 0.5*(gradUzeta_tmp(:,i,j,k)+gradUzeta_tmp(:,i,j,k+1))
#endif
END DO; END DO; END DO! i,j,k=0,PP_N

END SUBROUTINE FV_CalcGradients
#endif

END MODULE MOD_Lifting_VolInt_gen
#endif /*PARABOLIC*/
