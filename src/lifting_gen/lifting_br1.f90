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

MODULE MOD_Lifting_BR1_gen
IMPLICIT NONE
PRIVATE


INTERFACE Lifting_BR1_gen
   MODULE PROCEDURE Lifting_BR1
   MODULE PROCEDURE Lifting_BR1_Euler
END INTERFACE

PUBLIC:: Lifting_BR1_gen

CONTAINS
!==================================================================================================================================
!> \brief Contains the BR1 lifting procedures (initialization and lifting operator) for computing the lifted solution gradients
!> according to Bassi & Rebay 1997. The lifted gradients are required for the viscous fluxes.
!>
!> Local gradients of the DG polynomial are not well suited to compute the viscous fluxes due to their discontinuous nature.
!> Instead, modified gradients \f$ Q \approx \nabla_x U \f$ are introduced and this equation is discretized using the DG method.
!> The equation for the gradients can be discretized in the weak or the strong form, which implies integration by parts
!> once (weak) or twice (strong), see "On the Quadrature and Weak Form Choices in Collocation Type Discontinuous Galerkin Spectral
!> Element Methods" (Gassner & Kopriva 2010) for details. If the strong form is chosen, the volume integral can be computed in a
!> conservative and a non conservative way, see "Implementing Spectral Methods for Partial Differential Equations" (Kopriva 2009)
!> for details. In the non conservative form we derive the solution and multiply by the metric terms while the conservative
!> formulation derives the solution multiplied by the metric terms.
!> The numerical flux for the BR1 procedure is simply choosen as the arithmetic mean of the solution.
!>
!> This routine will be called after the current solution U has been prolonged to the element faces.
!> To compute the lifted gradients of the solution the following steps are taken:
!> - Compute and communicate the surface fluxes on the mpi interface, do this first to use latency hiding. The fluxes will be
!>   temporarily stored in the gradUxyz_master arrays. Surface fluxes are different if strong or weak form is used.
!> - Compute the volume integral. The gradients will also get nullified in this routine.
!>   There are different versions of the VolInt routine depending on the usage of the conservative (weak or strong)
!>   or non conservative (strong only) form.
!> - The surface fluxes for all remaining sides (boundaries and inner) will be computed.
!> - The surface integral is performed, first for all inner and boundary sides. Then the communication of the fluxes is finished
!>   and the surface integral of the remaining sides is performed.
!>   In the surface integral there is a distinction between weak and strong formulation. The fluxes for the strong form are
!>   different on the slave or master sides since we substract the inner solution, this is accounted for in the SurfInt routine.
!> - The gradients are transformed back to physical space to be used by the DG routines.
!> - The computed volume gradients are prolonged to the surfaces at the end of the routine.
!==================================================================================================================================
SUBROUTINE Lifting_BR1(nVarIn,nVarOut,UPrim,UPrim_master,UPrim_slave,gradUx,gradUy,gradUz)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars,             ONLY: doWeakLifting,doConservativeLifting
USE MOD_DG_Vars,                  ONLY: L_hatMinus,L_hatPlus
USE MOD_ApplyJacobianLifting_gen, ONLY: ApplyJacobianLifting_gen
USE MOD_FillMortarLifting_gen,    ONLY: Flux_MortarLifting_gen
USE MOD_Lifting_fillflux_gen,     ONLY: Lifting_FillFlux_gen, Lifting_FillFlux_BC_gen, Lifting_FillFlux_NormVec_gen
USE MOD_Lifting_VolInt_gen,       ONLY: Lifting_VolInt_gen
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI,                      ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_Mesh_Vars,                ONLY: nSides,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarIn
INTEGER,INTENT(IN) :: nVarOut
REAL,INTENT(IN)    :: UPrim(       nVarIn ,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution vector for which lifted gradients
                                                                         !> will be computed
REAL,INTENT(INOUT) :: UPrim_master(nVarIn ,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the master sides
REAL,INTENT(INOUT) :: UPrim_slave( nVarIn ,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the slave sides
REAL,INTENT(OUT)   :: gradUx(      nVarOut,0:PP_N,0:PP_N,0:PP_NZ,nElems)
REAL,INTENT(OUT)   :: gradUy(      nVarOut,0:PP_N,0:PP_N,0:PP_NZ,nElems)
REAL,INTENT(OUT)   :: gradUz(      nVarOut,0:PP_N,0:PP_N,0:PP_NZ,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: gradUz_slave( nVarOut,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the master sides
REAL               :: gradUx_master(nVarOut,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the master sides
REAL               :: gradUy_master(nVarOut,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the master sides
REAL               :: gradUz_master(nVarOut,0:PP_N,0:PP_NZ,1:nSides)      !< solution on the master sides
INTEGER            :: DataSizeSideGrad_loc
!==================================================================================================================================
! fill the global surface flux list
! #### use gradUz_slave for storing the fluxes, NormVec is applied later ####
! ## DO NOT USE SAME STORAGE for transformed/untransformed fluxes, since NormVec can be applied before communication is finished ##
DataSizeSideGrad_loc  =nVarOut*(PP_N+1)*(PP_NZ+1)
gradUz_slave = 0.
#if USE_MPI
! Receive YOUR
CALL StartReceiveMPIData(gradUz_slave,DataSizeSideGrad_loc,1,nSides,MPIRequest_Flux(:,RECV),SendID=1)
! Compute lifting MPI fluxes
CALL Lifting_FillFlux_gen(nVarIn,nVarOut,0, UPrim_master,UPrim_slave,gradUz_slave,doMPISides=.TRUE.)
! Start Send MINE
CALL StartSendMPIData(   gradUz_slave,DataSizeSideGrad_loc,1,nSides,MPIRequest_Flux(:,SEND),SendID=1)
#endif /*USE_MPI*/


! compute volume integral contribution and add to ut
IF(doWeakLifting.OR.doConservativeLifting)THEN
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,1,UPrim,GradUx)
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,2,UPrim,GradUy)
#if (PP_dim==3)
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,3,UPrim,GradUz)
#endif
ELSE
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,UPrim,GradUx,GradUy,GradUz)
END IF

! fill the all surface fluxes on this proc
CALL Lifting_FillFlux_BC_gen(nVarIn,nVarOut,UPrim_master,Flux=gradUz_slave)
CALL Lifting_FillFlux_gen(nVarIn,nVarOut,0,UPrim_master,UPrim_slave,   gradUz_slave,doMPISides=.FALSE.)
! at this point BC, inner and MPI MINE are filled
CALL Lifting_FillFlux_NormVec_gen(nVarOut,gradUz_slave,gradUx_master,gradUy_master,gradUz_master,doMPISides=.FALSE.)

! Attention: we only have one Flux (gradUx/y/z_master) for the Lifting
!            => input it to Flux_Mortar for both fluxes (master/slave)
CALL Flux_MortarLifting_gen(gradUx_master,gradUx_master,doMPISides=.FALSE.,weak=doWeakLifting)
CALL Flux_MortarLifting_gen(gradUy_master,gradUy_master,doMPISides=.FALSE.,weak=doWeakLifting)
#if (PP_dim==3)
CALL Flux_MortarLifting_gen(gradUz_master,gradUz_master,doMPISides=.FALSE.,weak=doWeakLifting)
#endif


! compute surface integral contribution and add to ut
CALL Lifting_SurfInt_BR1(nVarOut,PP_N,gradUx_master,gradUx,.FALSE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
CALL Lifting_SurfInt_BR1(nVarOut,PP_N,gradUy_master,gradUy,.FALSE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#if (PP_dim==3)
CALL Lifting_SurfInt_BR1(nVarOut,PP_N,gradUz_master,gradUz,.FALSE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#endif
#if USE_MPI
! Complete send / receive
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux)
CALL Lifting_FillFlux_NormVec_gen(nVarOut,gradUz_slave,gradUx_master,gradUy_master,gradUz_master,doMPISides=.TRUE.)
! Attention: we only have one Flux (gradUx/y/z_master) for the Lifting
!            => input it to Flux_Mortar for both fluxes (master/slave)
CALL Flux_MortarLifting_gen(gradUx_master,gradUx_master,doMPISides=.TRUE.,weak=doWeakLifting)
CALL Flux_MortarLifting_gen(gradUy_master,gradUy_master,doMPISides=.TRUE.,weak=doWeakLifting)
#if (PP_dim==3)
CALL Flux_MortarLifting_gen(gradUz_master,gradUz_master,doMPISides=.TRUE.,weak=doWeakLifting)
#endif
! compute surface integral contribution and add to ut
CALL Lifting_SurfInt_BR1(nVarOut,PP_N,gradUx_master,gradUx,.TRUE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
CALL Lifting_SurfInt_BR1(nVarOut,PP_N,gradUy_master,gradUy,.TRUE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#if (PP_dim==3)
CALL Lifting_SurfInt_BR1(nVarOut,PP_N,gradUz_master,gradUz,.TRUE.,L_hatMinus,L_hatPlus,weak=doWeakLifting)
#endif
#endif /*USE_MPI*/

! Account for the jacobian
! The Lifting already has the right sign
! For FV elements no Jacobian is needed here, since already applied in Lifting_Volint (hidden in FV_Metrics_f/g/hTilde_sJ)
CALL ApplyJacobianLifting_gen(gradUx,toPhysical=.TRUE.,FVE=0)
CALL ApplyJacobianLifting_gen(gradUy,toPhysical=.TRUE.,FVE=0)
#if (PP_dim==3)
CALL ApplyJacobianLifting_gen(gradUz,toPhysical=.TRUE.,FVE=0)
#endif

END SUBROUTINE Lifting_BR1

SUBROUTINE Lifting_BR1_Euler(nVarIn,nVarOut,UPrim,gradUx,gradUy,gradUz)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Lifting_Vars,             ONLY: doWeakLifting,doConservativeLifting
USE MOD_ApplyJacobian,            ONLY: ApplyJacobian
USE MOD_Lifting_VolInt_gen,       ONLY: Lifting_VolInt_gen
USE MOD_Mesh_Vars,                ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarIn
INTEGER,INTENT(IN) :: nVarOut
REAL,INTENT(IN)    :: UPrim( nVarIn ,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution vector for which lifted gradients
                                                                   !> will be computed
REAL,INTENT(OUT)   :: gradUx(nVarOut,0:PP_N,0:PP_N,0:PP_NZ,nElems)
REAL,INTENT(OUT)   :: gradUy(nVarOut,0:PP_N,0:PP_N,0:PP_NZ,nElems)
REAL,INTENT(OUT)   :: gradUz(nVarOut,0:PP_N,0:PP_N,0:PP_NZ,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! compute volume integral contribution and add to ut
IF(doWeakLifting.OR.doConservativeLifting)THEN
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,1,UPrim,GradUx)
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,2,UPrim,GradUy)
#if (PP_dim==3)
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,3,UPrim,GradUz)
#endif
ELSE
  CALL Lifting_VolInt_gen(nVarIn,nVarOut,UPrim,GradUx,GradUy,GradUz)
END IF

! Account for the jacobian
! The Lifting already has the right sign
! For FV elements no Jacobian is needed here, since already applied in Lifting_Volint (hidden in FV_Metrics_f/g/hTilde_sJ)
CALL ApplyJacobian(1,gradUx,toPhysical=.TRUE.,FVE=0)
CALL ApplyJacobian(1,gradUy,toPhysical=.TRUE.,FVE=0)
#if (PP_dim==3)
CALL ApplyJacobian(1,gradUz,toPhysical=.TRUE.,FVE=0)
#endif

END SUBROUTINE Lifting_BR1_Euler

!==================================================================================================================================
!> \brief Surface integral in the BR1 scheme optimized for performance, for weak or strong formulation.
!>
!> Performs the surface integral for the BR1 routine. Uses the DoSurfInt routine from the DG operator to perform actual
!> integration.
!> If we use the strong formulation, the inner solution is substracted from the numerical flux. This is done in the FillFlux
!> routines for the master side. The flux on the master side is \f$ \frac{1}{2} (U^+ + U^-) - U^- = \frac{1}{2} (U^+ - U^-) \f$
!> since \f$ \frac{1}{2} (U^+ + U^-) \f$ is the numerical flux in the BR1 scheme. On the slave side, the flux becomes
!> \f$ \frac{1}{2} (U^+ + U^-) - U^+ = \frac{1}{2} (-U^+ + U^-) \f$ which is simply the master side flux multiplied by
!> \f$ -1 \f$. This means we don't have to flip the sign on the flux for the slave side in strong form as we normally do to
!> get the flux on the slave side.
!==================================================================================================================================
PPURE SUBROUTINE Lifting_SurfInt_BR1(nVar,Nloc,Flux,gradU,doMPISides,L_HatMinus,L_HatPlus,weak)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SurfIntLifting_gen, ONLY: DoSurfIntLifting_gen
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides,nElems
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_MINE
USE MOD_Mesh_Vars,          ONLY: S2V2
USE MOD_Mesh_Vars,          ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems_master,FV_Elems_slave
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVar                                             !< number of variables
INTEGER,INTENT(IN) :: Nloc                                             !< polynomial degree
LOGICAL,INTENT(IN) :: doMPISides                                       !< = .TRUE. only MPISides_YOUR+MPIMortar are filled
                                                                       !< =.FALSE. BCSides+(Mortar-)InnerSides+MPISides_MINE
REAL,INTENT(IN)    :: Flux(1:nVar,0:Nloc,0:ZDIM(Nloc),nSides)          !< flux to be filled
REAL,INTENT(IN)    :: L_HatPlus(0:Nloc)                                !< lagrange polynomials at xi=+1 and pre-divided by
                                                                       !< integration weight
REAL,INTENT(IN)    :: L_HatMinus(0:Nloc)                               !< lagrange polynomials at xi=-1 and pre-divided by
                                                                       !< integration weight
REAL,INTENT(INOUT) :: gradU(1:nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)!< time derivative of solution
LOGICAL,INTENT(IN) :: weak                                             !< switch for weak or strong formulation
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,nbElemID,locSideID,nblocSideID,SideID,p,q,flip
INTEGER            :: firstSideID,lastSideID
REAL               :: FluxTmp(1:nVar,0:Nloc,0:ZDIM(Nloc))
!==================================================================================================================================
IF(doMPISides)THEN
  ! MPI YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  ! inner sides and MPI mine
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ElemID      = SideToElem(S2E_ELEM_ID,   SideID)
  nbElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)

  ! master sides
  IF(ElemID.GT.0)THEN
    IF (FV_Elems_master(SideID).EQ.0) THEN ! DG element
      locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
      flip      = 0
      ! orient flux to fit flip and locSide to element local system
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        ! note: for master sides, the mapping S2V2 should be a unit matrix
        FluxTmp(:,S2V2(1,p,q,flip,locSideID),S2V2(2,p,q,flip,locSideID)) = Flux(:,p,q,SideID)
      END DO; END DO ! p,q
#if (PP_NodeType==1 || (PP_NodeType==2 && defined(EXACT_MM)))
      CALL DoSurfIntLifting_gen(Nloc,FluxTmp,L_HatMinus,   L_HatPlus,      locSideID,gradU(:,:,:,:,ElemID))
#elif (PP_NodeType==2 && !defined(EXACT_MM))
      CALL DoSurfIntLifting_gen(Nloc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),locSideID,gradU(:,:,:,:,ElemID))
#endif
    END IF
  END IF

  ! slave sides
  IF(nbElemID.GT.0)THEN
    IF (FV_Elems_slave(SideID).EQ.0) THEN ! DG element
      nblocSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
      flip        = SideToElem(S2E_FLIP,SideID)
      ! orient flux to fit flip and locSide to element local system
      IF(weak)THEN
        DO q=0,ZDIM(Nloc); DO p=0,Nloc
          ! p,q are in the master RHS system, they need to be transformed to the slave volume system using S2V2 mapping
          FluxTmp(:,S2V2(1,p,q,flip,nblocSideID),S2V2(2,p,q,flip,nblocSideID))=-Flux(:,p,q,SideID)
        END DO; END DO ! p,q
      ELSE
        ! In strong form, don't flip the sign since the slave flux is the negative of the master flux
        DO q=0,ZDIM(Nloc); DO p=0,Nloc
          ! p,q are in the master RHS system, they need to be transformed to the slave volume system using S2V2 mapping
          FluxTmp(:,S2V2(1,p,q,flip,nblocSideID),S2V2(2,p,q,flip,nblocSideID))= Flux(:,p,q,SideID)
        END DO; END DO ! p,q
      END IF
#if (PP_NodeType==1 || (PP_NodeType==2 && defined(EXACT_MM)))
      CALL DoSurfIntLifting_gen(Nloc,FluxTmp,L_HatMinus,   L_HatPlus,      nblocSideID,gradU(:,:,:,:,nbElemID))
#elif (PP_NodeType==2 && !defined(EXACT_MM))
      CALL DoSurfIntLifting_gen(Nloc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),nblocSideID,gradU(:,:,:,:,nbElemID))
#endif
    END IF
  END IF
END DO ! SideID=1,nSides
END SUBROUTINE Lifting_SurfInt_BR1

END MODULE MOD_Lifting_BR1_gen
