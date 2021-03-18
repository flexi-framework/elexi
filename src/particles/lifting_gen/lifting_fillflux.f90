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

MODULE MOD_Lifting_fillflux_gen
IMPLICIT NONE
PRIVATE


INTERFACE Lifting_FillFlux_gen
   MODULE PROCEDURE Lifting_FillFlux
END INTERFACE

INTERFACE Lifting_FillFlux_BC_gen
   MODULE PROCEDURE Lifting_FillFlux_BC
END INTERFACE

INTERFACE Lifting_FillFlux_NormVec_gen
   MODULE PROCEDURE Lifting_FillFlux_NormVec
END INTERFACE

PUBLIC::Lifting_FillFlux_gen, Lifting_FillFlux_BC_gen, Lifting_FillFlux_NormVec_gen

CONTAINS
!==================================================================================================================================
!> \brief Computes the untransformed (BR1) and transformated (BR2) surface fluxes
!==================================================================================================================================
PPURE SUBROUTINE Lifting_FillFlux(nVarIn,nVarOut,dir,UPrimface_master,UPrimface_slave,Flux,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars,    ONLY: doWeakLifting
USE MOD_Mesh_Vars,       ONLY: SurfElem
USE MOD_Mesh_Vars,       ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide
USE MOD_Mesh_Vars,       ONLY: firstMPISide_MINE,lastMPISide_MINE
#if FV_ENABLED
USE MOD_FV_Vars         ,ONLY: FV_Elems_Sum,FV_sVdm
USE MOD_ChangeBasisByDim,ONLY: ChangeBasisSurf
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarIn
INTEGER,INTENT(IN) :: nVarOut
INTEGER,INTENT(IN) :: dir                                             !< direction (x,y,z)
LOGICAL,INTENT(IN) :: doMPISides                                      !< = .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(INOUT) :: UPrimface_master(nVarIn ,0:PP_N,0:PP_NZ,nSides) !< Solution on master sides
REAL,INTENT(INOUT) :: UPrimface_slave( nVarOut,0:PP_N,0:PP_NZ,nSides) !< Solution on slave sides
REAL,INTENT(INOUT) :: Flux(            nVarOut,0:PP_N,0:PP_NZ,nSides) !< Untransformed lifting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID,sig
#if FV_ENABLED
REAL               :: UPrim_glob(1:nVarIn,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide
   lastSideID = lastInnerSide
END IF

! for strong form subtract solution from the inside
sig = MERGE(1, -1, doWeakLifting)

DO SideID=firstSideID,lastSideID
#if FV_ENABLED
  SELECT CASE(FV_Elems_Sum(SideID))
  CASE(0) ! both DG
#endif
    Flux(:,:,:,SideID) = sig*UPrimface_master(:,:,:,SideID) + UPrimface_slave(:,:,:,SideID)
#if FV_ENABLED
  CASE(1) ! master=FV, slave=DG
    CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_sVdm,UPrimface_master(:,:,:,SideID),UPrim_glob)
    Flux(:,:,:,SideID) = sig*UPrim_glob(:,:,:) + UPrimface_slave(:,:,:,SideID)
  CASE(2) ! master=DG, slave=FV
    CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_sVdm,UPrimface_slave(:,:,:,SideID),UPrim_glob)
    Flux(:,:,:,SideID) = sig*UPrimface_master(:,:,:,SideID) + UPrim_glob(:,:,:)
  CASE(3) ! both FV
    CYCLE
  END SELECT
#endif
  IF(dir.EQ.0)THEN
    ! BR1 uses arithmetic mean value of states for the Riemann flux
    ! Transformation with NormVec is applied in separate routine
    DO q=0,PP_NZ
      DO p=0,PP_N
        Flux(:,p,q,SideID)=0.5*Flux(:,p,q,SideID)*SurfElem(p,q,0,SideID)
      END DO
    END DO
  ELSE
    DO q=0,PP_NZ
      DO p=0,PP_N
        Flux(:,p,q,SideID)=0.5*Flux(:,p,q,SideID)*NormVec(dir,p,q,0,SideID)*SurfElem(p,q,0,SideID)
      END DO
    END DO
  END IF
END DO ! SideID

END SUBROUTINE Lifting_FillFlux

!==================================================================================================================================
!> \brief Computes the BR1 surface fluxes for boundary conditions, surfelem contribution is considered as well.
!>
!> Fills the Untransformed boundary surface fluxes for the BR1 scheme. The numerical flux in the
!> BR1 lifting is simply taken as the arithmetic mean of the solution.
!> This routine calls the equation system dependant routine Lifting_GetBoundaryFlux which will fill the flux depending on the
!> boundary condition that has to be applied. The Lifting_GetBoundaryFlux routine will also differentiate between weak and
!> strong form and already multiply the flux by the surface element.
!>
!> For the BR1, the multiplication by the normal vector is done in a separate routine.
!>
!> The flux is filled for the master side, the contribution for the slave side (which is different because the inner solution
!> is equal to \f$ U^+ \f$) is taken into account in the SurfInt routine.
!==================================================================================================================================
SUBROUTINE Lifting_FillFlux_BC(nVarIn,nVarOut,UPrim_master,Flux)!,FluxX,FluxY,FluxZ)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: nBCSides,nSides,SurfElem
!USE MOD_Mesh_Vars,       ONLY: NormVec
USE MOD_Lifting_Vars,    ONLY: doWeakLifting
#if FV_ENABLED
USE MOD_FV_Vars,         ONLY: FV_Elems_master
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nVarIn
INTEGER,INTENT(IN)          :: nVarOut
REAL,INTENT(IN)             :: UPrim_master(nVarIn,0:PP_N,0:PP_NZ,1:nSides) !< Primitive solution from the inside
REAL,OPTIONAL,INTENT(OUT)   :: Flux (nVarOut,0:PP_N,0:PP_NZ,1:nSides)       !< gradient flux
!REAL,OPTIONAL,INTENT(OUT)   :: FluxX(nVarOut,0:PP_N,0:PP_NZ,1:nSides)       !< gradient flux in x-dir
!REAL,OPTIONAL,INTENT(OUT)   :: FluxY(nVarOut,0:PP_N,0:PP_NZ,1:nSides)       !< gradient flux in y-dir
!REAL,OPTIONAL,INTENT(OUT)   :: FluxZ(nVarOut,0:PP_N,0:PP_NZ,1:nSides)       !< gradient flux in z-dir
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: SideID,p,q
!REAL                        :: Flux_loc(nVarOut,0:PP_N,0:PP_NZ)             !< local gradient flux
!==================================================================================================================================

Flux = 0.
DO SideID=1,nBCSides
#if FV_ENABLED
  IF (FV_Elems_master(SideID).GT.0) CYCLE
#endif
  ! flux = 0.5 * (UPrim_master(:,p,q,SideID) + UPrim_master(:,p,q,SideID))
  IF(PRESENT(Flux))THEN
    DO q=0,PP_NZ; DO p=0,PP_N
      Flux(:,p,q,SideID)=UPrim_master(:,p,q,SideID)
    END DO; END DO
  END IF
  !in case lifting is done in strong form
  IF(.NOT.doWeakLifting) Flux(:,:,:,SideID)=Flux(:,:,:,SideID)-UPrim_master(:,:,:,SideID)

  DO q=0,PP_NZ; DO p=0,PP_N
    Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*SurfElem(p,q,0,SideID)
  END DO; END DO
!  IF(PRESENT(FluxX).AND.PRESENT(FluxY).AND.PRESENT(FluxZ))THEN
!    DO q=0,PP_NZ; DO p=0,PP_N
!      FluxX(:,p,q,SideID)=Flux_loc(:,p,q)*NormVec(1,p,q,0,SideID)
!      FluxY(:,p,q,SideID)=Flux_loc(:,p,q)*NormVec(2,p,q,0,SideID)
!#if PP_dim == 3
!      FluxZ(:,p,q,SideID)=Flux_loc(:,p,q)*NormVec(3,p,q,0,SideID)
!#else
!      FluxZ(:,p,q,SideID)=0.
!#endif
!    END DO; END DO
!  END IF
END DO

END SUBROUTINE Lifting_FillFlux_BC

!==================================================================================================================================
!> \brief Multiplies the untransformed flux (already weighted by the surface element) by the normal vector to complete the
!>        transformation into reference space.
!>
!> The untransformed fluxes is multiplied for each gradient direction by the corresponding normal vector,
!> which is the second step of the transformation. We end up with the lifting fluxes for X/Y/Z direction.
!==================================================================================================================================
PPURE SUBROUTINE Lifting_FillFlux_NormVec(nVarOut,Flux,FluxX,FluxY,FluxZ,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,nSides,MortarType
USE MOD_Mesh_Vars,       ONLY: firstMPISide_YOUR,lastMPISide_MINE,lastMPISide_YOUR
#if FV_ENABLED
USE MOD_FV_Vars,         ONLY: FV_Elems_master,FV_Elems_Sum
USE MOD_Mesh_Vars,       ONLY: nBCSides
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarOut
LOGICAL,INTENT(IN) :: doMPISides                             !< = .TRUE. MPI_YOUR
                                                             !< =.FALSE. BC+Inner+MPI_Mine
REAL,INTENT(IN)    :: Flux( nVarOut,0:PP_N,0:PP_NZ,1:nSides) !< Untransformed Lifting boundary flux
REAL,INTENT(INOUT) :: FluxX(nVarOut,0:PP_N,0:PP_NZ,1:nSides) !< Lifting boundary flux in x direction
REAL,INTENT(INOUT) :: FluxY(nVarOut,0:PP_N,0:PP_NZ,1:nSides) !< Lifting boundary flux in y direction
REAL,INTENT(INOUT) :: FluxZ(nVarOut,0:PP_N,0:PP_NZ,1:nSides) !< Lifting boundary flux in z direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: firstSideID,lastSideID,SideID,p,q
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! apply normal vector to YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID =  lastMPISide_YOUR
ELSE
  ! apply normal vector to all MINE fluxes
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
#if FV_ENABLED
  IF (FV_Elems_master(SideID).GT.0.AND.SideID.LE.nBCSides) CYCLE ! FV BC already filled
  IF (FV_Elems_Sum(SideID).EQ.3)                           CYCLE ! FV/FV already filled
#endif
  IF(MortarType(1,SideID).GT.0)                            CYCLE ! no big mortars
  DO q=0,PP_NZ; DO p=0,PP_N
    FluxX(:,p,q,SideID)=Flux(:,p,q,SideID)*NormVec(1,p,q,0,SideID)
    FluxY(:,p,q,SideID)=Flux(:,p,q,SideID)*NormVec(2,p,q,0,SideID)
#if (PP_dim==3)
    FluxZ(:,p,q,SideID)=Flux(:,p,q,SideID)*NormVec(3,p,q,0,SideID)
#else
    FluxZ(:,p,q,SideID)=0.
#endif
  END DO; END DO
END DO

END SUBROUTINE Lifting_FillFlux_NormVec


END MODULE MOD_Lifting_fillflux_gen
