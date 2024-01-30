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
!> Contains routines for transformation from reference to physical space and vice versa
!===================================================================================================================================
MODULE MOD_Eval_xyz
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE GetPositionInRefElem
  MODULE PROCEDURE GetPositionInRefElem
END INTERFACE

INTERFACE TensorProductInterpolation
  MODULE PROCEDURE TensorProductInterpolation
END INTERFACE

INTERFACE EvaluateFieldAtPhysPos
  MODULE PROCEDURE EvaluateFieldAtPhysPos
END INTERFACE

INTERFACE EvaluateFieldAtRefPos
  MODULE PROCEDURE EvaluateFieldAtRefPos
END INTERFACE

#if USE_EXTEND_RHS || USE_FAXEN_CORR
INTERFACE EvaluateFieldAndGradAtPhysPos
  MODULE PROCEDURE EvaluateFieldAndGradAtPhysPos
END INTERFACE

INTERFACE EvaluateFieldAndGradAtRefPos
  MODULE PROCEDURE EvaluateFieldAndGradAtRefPos
END INTERFACE
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */

#if FV_ENABLED
INTERFACE EvaluateField_FV
  MODULE PROCEDURE EvaluateField_FV
END INTERFACE
#endif /* FV_ENABLED */

PUBLIC :: GetPositionInRefElem
PUBLIC :: TensorProductInterpolation
PUBLIC :: EvaluateFieldAtPhysPos
PUBLIC :: EvaluateFieldAtRefPos
#if USE_EXTEND_RHS || USE_FAXEN_CORR
PUBLIC :: EvaluateFieldAndGradAtPhysPos
PUBLIC :: EvaluateFieldAndGradAtRefPos
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */
#if FV_ENABLED
PUBLIC :: EvaluateField_FV
#endif /* FV_ENABLED */
!===================================================================================================================================

CONTAINS

!SUBROUTINE GetPositionInRefElem(x_in,xi,ElemID,DoReUseMap,ForceMode)
SUBROUTINE GetPositionInRefElem(x_in,xi,ElemID,ForceMode)
!===================================================================================================================================
!> Get Position within reference element (x_in -> xi=[-1,1])
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                   ONLY: LagrangeInterpolationPolys
USE MOD_Mesh_Vars,               ONLY: NGeo
USE MOD_Particle_Mesh_Tools,     ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars,      ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars,      ONLY: wBaryCL_NGeo1,XiCL_NGeo1
USE MOD_Particle_Mesh_Vars,      ONLY: BoundsOfElem_Shared,ElemCurved
#if USE_MPI
USE MOD_Particle_Mesh_Vars,      ONLY: XCL_NGeo_Shared,dXCL_NGeo_Shared
#else
USE MOD_Particle_Mesh_Vars,      ONLY: dXCL_NGeo,XCL_NGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: ElemID                                 !< Global element index
REAL,INTENT(IN)             :: x_in(3)                                !< position in physical space
! LOGICAL,INTENT(IN),OPTIONAL :: DoReUseMap                             !< flag if start values for Newton elem mapping already exists
LOGICAL,INTENT(IN),OPTIONAL :: ForceMode                              !< flag for mode change in RefElemNewton
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: xi(1:3)                                !< position in reference element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: CNElemID,iMode
REAL                       :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                       :: dXCL_NGeo1(1:3,1:3,0:1,0:1,0:1)
!===================================================================================================================================

! Check whether the particle position intersects with the element bounding box
IF (     x_in(1).LT.BoundsOfElem_Shared(1,1,ElemID) .OR. x_in(1).GT.BoundsOfElem_Shared(2,1,ElemID)  &
    .OR. x_in(2).LT.BoundsOfElem_Shared(1,2,ElemID) .OR. x_in(2).GT.BoundsOfElem_Shared(2,2,ElemID)  &
    .OR. x_in(3).LT.BoundsOfElem_Shared(1,3,ElemID) .OR. x_in(3).GT.BoundsOfElem_Shared(2,3,ElemID)) THEN
  xi = HUGE(1.)
  RETURN
END IF

#if USE_MPI
ASSOCIATE( XCL_NGeo  => XCL_NGeo_Shared     &
         ,dXCL_NGeo  => dXCL_NGeo_Shared)
#endif

IF (PRESENT(ForceMode)) THEN; iMode = 1
ELSE                        ; iMode = 2
END IF

!IF (.NOT.PRESENT(DoReUseMap)) THEN
  CALL GetRefNewtonStartValue(X_in,Xi,ElemID)
!END IF

CNElemID = GetCNElemID(ElemID)

IF (ElemCurved(CNElemID)) THEN
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID),NGeo,ElemID,Mode=iMode)
ELSE
  ! fill dummy XCL_NGeo1
  IF(NGeo.EQ.1)THEN
    CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID),NGeo,ElemID,Mode=iMode)
  ELSE
    XCL_NGeo1 (1:3,    0,0,0) = XCL_NGeo (1:3,     0  , 0  , 0  ,ElemID)
    XCL_NGeo1 (1:3,    1,0,0) = XCL_NGeo (1:3,    NGeo, 0  , 0  ,ElemID)
    XCL_NGeo1 (1:3,    0,1,0) = XCL_NGeo (1:3,     0  ,NGeo, 0  ,ElemID)
    XCL_NGeo1 (1:3,    1,1,0) = XCL_NGeo (1:3,    NGeo,NGeo, 0  ,ElemID)
    XCL_NGeo1 (1:3,    0,0,1) = XCL_NGeo (1:3,     0  , 0  ,NGeo,ElemID)
    XCL_NGeo1 (1:3,    1,0,1) = XCL_NGeo (1:3,    NGeo, 0  ,NGeo,ElemID)
    XCL_NGeo1 (1:3,    0,1,1) = XCL_NGeo (1:3,     0  ,NGeo,NGeo,ElemID)
    XCL_NGeo1 (1:3,    1,1,1) = XCL_NGeo (1:3,    NGeo,NGeo,NGeo,ElemID)
    ! fill dummy dXCL_NGeo1
    dXCL_NGeo1(1:3,1:3,0,0,0) = dXCL_NGeo(1:3,1:3, 0  , 0  , 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,1,0,0) = dXCL_NGeo(1:3,1:3,NGeo, 0  , 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,0,1,0) = dXCL_NGeo(1:3,1:3, 0  ,NGeo, 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,1,1,0) = dXCL_NGeo(1:3,1:3,NGeo,NGeo, 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,0,0,1) = dXCL_NGeo(1:3,1:3, 0  , 0  ,NGeo,ElemID)
    dXCL_NGeo1(1:3,1:3,1,0,1) = dXCL_NGeo(1:3,1:3,NGeo, 0  ,NGeo,ElemID)
    dXCL_NGeo1(1:3,1:3,0,1,1) = dXCL_NGeo(1:3,1:3, 0  ,NGeo,NGeo,ElemID)
    dXCL_NGeo1(1:3,1:3,1,1,1) = dXCL_NGeo(1:3,1:3,NGeo,NGeo,NGeo,ElemID)
    CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo1,XiCL_NGeo1,XCL_NGeo1,dXCL_NGeo1,1,ElemID,Mode=iMode)
  END IF
END IF

#if USE_MPI
END ASSOCIATE
#endif

END SUBROUTINE GetPositionInRefElem


PPURE SUBROUTINE TensorProductInterpolation(Xi_in,NVar,N_in,xGP_in,wBary_In,U_In,U_Out)
!===================================================================================================================================
!> Interpolates a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation points to the position Xi
!===================================================================================================================================
! MODULES
USE MOD_Basis,                   ONLY: LagrangeInterpolationPolys
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)            :: Xi_in(3)                              !< position in reference element
INTEGER,INTENT(IN)         :: NVar                                  !< 5 (rho,u_x,u_y,u_z,e)
INTEGER,INTENT(IN)         :: N_In                                  !< usually PP_N
REAL,INTENT(IN)            :: xGP_In(0:N_in)
REAL,INTENT(IN)            :: wBary_In(0:N_in)
REAL,INTENT(IN)            :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)     !< State in element
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)           :: U_Out(1:NVar)                         !< Interpolated state at reference position xi_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i,j,k
REAL,DIMENSION(3,0:N_in)  :: L_xi
REAL                      :: L_eta_zeta
!===================================================================================================================================
CALL LagrangeInterpolationPolys(xi_in(1),N_in,xGP_in,wBary_In,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi_in(2),N_in,xGP_in,wBary_In,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi_in(3),N_in,xGP_in,wBary_In,L_xi(3,:))

U_out(:)=0
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta=L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      U_out = U_out + U_IN(:,i,j,k)*L_xi(1,i)*L_eta_zeta
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In

END SUBROUTINE TensorProductInterpolation


SUBROUTINE EvaluateFieldAtPhysPos(x_in,NVar,N_in,U_In,NVar_out,U_Out,ElemID,PartID)
!===================================================================================================================================
!> 1) Get position within reference element (x_in -> xi=[-1,1]) by inverting the mapping
!> 2) interpolate DG solution to position (U_In -> U_Out(x_in))
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
USE MOD_Mesh_Vars,             ONLY: NGeo
USE MOD_Particle_Mesh_Tools,   ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars,    ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars,    ONLY: ElemCurved,wBaryCL_NGeo1,XiCL_NGeo1
USE MOD_Particle_Vars,         ONLY: PDM,PartState,LastPartPos
#if USE_MPI
USE MOD_Particle_Mesh_Vars,    ONLY: XCL_NGeo_Shared,dXCL_NGeo_Shared
#else
USE MOD_Particle_Mesh_Vars,    ONLY: XCL_NGeo,dXCL_NGeo
#endif /*USE_MPI*/
USE MOD_Eos,                   ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: x_in(3)                                       !< position in physical space
INTEGER,INTENT(IN)        :: NVar                                          !< 5 (rho,u_x,u_y,u_z,e)
INTEGER,INTENT(IN)        :: N_In                                          !< usually PP_N
INTEGER,INTENT(IN)        :: ElemID                                        !< Element index
REAL,INTENT(IN)           :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)             !< State in Element
INTEGER,INTENT(IN)        :: PartID                                        !< particle ID
INTEGER,INTENT(IN)        :: NVar_out                                      !< 6 (rho,u_x,u_y,u_z,p,T)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: U_Out(1:NVar_out)                             !< Interpolated state at physical position x_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: CNElemID,i,j,k
REAL                      :: xi(3)
REAL                      :: L_xi(3,0:PP_N), L_eta_zeta
REAL                      :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                      :: dXCL_NGeo1(1:3,1:3,0:1,0:1,0:1)
REAL                      :: Utmp(1:nVar)
!===================================================================================================================================

! Sanity check the values
IF (ANY(IEEE_IS_NAN(X_in))) THEN
  WRITE(UNIT_stdOut,'(A,I0,A,I0)')        ' NaN detected for PartID ', PartID, ' on proc ',myRank
  WRITE(UNIT_stdOut,'(A18,3(1X,E27.16))') ' LastPosition   ', LastPartPos(1:3,PartID)
  WRITE(UNIT_stdOut,'(A18,3(1X,E27.16))') ' Velocity       ', PartState  (4:6,PartID)
  WRITE(UNIT_stdOut,'(A)')                ' Removing particle ...'
  PDM%ParticleInside(PartID) = .FALSE.
  RETURN
END IF

#if USE_MPI
ASSOCIATE( XCL_NGeo  =>  XCL_NGeo_Shared    &
         ,dXCL_NGeo  => dXCL_NGeo_Shared)
#endif /*USE_MPI*/

CALL GetRefNewtonStartValue(X_in,Xi,ElemID)

CNElemID = GetCNElemID(ElemID)

! If the element is curved, all Gauss points are required
IF (ElemCurved(CNElemID)) THEN
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID) &
                    ,NGeo,ElemID,Mode=3,PartID=PartID)
! If the element is not curved, only the corner nodes are required
ELSE
  ! fill dummy XCL_NGeo1
  XCL_NGeo1 (1:3,    0,0,0) = XCL_NGeo (1:3,     0  , 0  , 0  ,ElemID)
  XCL_NGeo1 (1:3,    1,0,0) = XCL_NGeo (1:3,    NGeo, 0  , 0  ,ElemID)
  XCL_NGeo1 (1:3,    0,1,0) = XCL_NGeo (1:3,     0  ,NGeo, 0  ,ElemID)
  XCL_NGeo1 (1:3,    1,1,0) = XCL_NGeo (1:3,    NGeo,NGeo, 0  ,ElemID)
  XCL_NGeo1 (1:3,    0,0,1) = XCL_NGeo (1:3,     0  , 0  ,NGeo,ElemID)
  XCL_NGeo1 (1:3,    1,0,1) = XCL_NGeo (1:3,    NGeo, 0  ,NGeo,ElemID)
  XCL_NGeo1 (1:3,    0,1,1) = XCL_NGeo (1:3,     0  ,NGeo,NGeo,ElemID)
  XCL_NGeo1 (1:3,    1,1,1) = XCL_NGeo (1:3,    NGeo,NGeo,NGeo,ElemID)
  ! fill dummy dXCL_NGeo1
  dXCL_NGeo1(1:3,1:3,0,0,0) = dXCL_NGeo(1:3,1:3, 0  , 0  , 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,1,0,0) = dXCL_NGeo(1:3,1:3,NGeo, 0  , 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,0,1,0) = dXCL_NGeo(1:3,1:3, 0  ,NGeo, 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,1,1,0) = dXCL_NGeo(1:3,1:3,NGeo,NGeo, 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,0,0,1) = dXCL_NGeo(1:3,1:3, 0  , 0  ,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,1,0,1) = dXCL_NGeo(1:3,1:3,NGeo, 0  ,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,0,1,1) = dXCL_NGeo(1:3,1:3, 0  ,NGeo,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,1,1,1) = dXCL_NGeo(1:3,1:3,NGeo,NGeo,NGeo,ElemID)
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo1,XiCL_NGeo1,XCL_NGeo1,dXCL_NGeo1,1,ElemID,Mode=1,PartID=PartID)
END IF

! 2.1) get "Vandermonde" vectors
CALL LagrangeInterpolationPolys(xi(1),N_in,xGP,wBary,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi(2),N_in,xGP,wBary,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi(3),N_in,xGP,wBary,L_xi(3,:))

! "more efficient" - Quote Thomas B.
Utmp(:)=0.
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta=L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      Utmp = Utmp + U_IN(:,i,j,k)*L_xi(1,i)*L_Eta_Zeta
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In

#if USE_MPI
END ASSOCIATE
#endif

! Convert to primitive variables
IF (NVar_out.EQ.PP_nVarPrim) THEN
  CALL ConsToPrim(U_out,Utmp)
ELSE
  U_out = Utmp
END IF

END SUBROUTINE EvaluateFieldAtPhysPos


#if FV_ENABLED
SUBROUTINE EvaluateField_FV(x_in,NVar,N_in,U_In,NVar_out,U_Out,ElemID)
!===================================================================================================================================
!> 1) Get position within reference element (x_in -> xi=[-1,1]) by inverting the mapping
!> 2) interpolate FV solution to position (U_In -> U_Out(x_in))
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
#if FV_RECONSTRUCT
USE MOD_FV_Vars               ,ONLY: gradUxi,gradUeta
#if PP_dim == 3
USE MOD_FV_Vars               ,ONLY: gradUzeta
#endif /* PP_dim == 3 */
USE MOD_FV_Vars               ,ONLY: FV_X,FV_Path_XI,FV_Path_ETA,FV_Path_ZETA
USE MOD_Interpolation_Vars    ,ONLY: xGP,wGP,wBary
USE MOD_FV_Metrics            ,ONLY: Integrate_Path1D
#endif /*FV_RECONSTRUCT*/
USE MOD_Eos                   ,ONLY: ConsToPrim
USE MOD_Mesh_Vars             ,ONLY: offsetElem
USE MOD_Particle_Tracking_Vars,ONLY: TrackingMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: x_in(3)                                       !< position in physical/reference space
INTEGER,INTENT(IN)        :: NVar                                          !< 5 (rho,u_x,u_y,u_z,e)
INTEGER,INTENT(IN)        :: N_In                                          !< usually PP_N
INTEGER,INTENT(IN)        :: ElemID                                        !< Element index
REAL,INTENT(IN)           :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)             !< State in Element
INTEGER,INTENT(IN)        :: NVar_out                                      !< 6 (rho,u_x,u_y,u_z,p,T)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: U_Out(1:NVar_out)                             !< Interpolated state at physical position x_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                          :: U_Prim(1:nVar_out,0:N_In,0:N_In,0:N_In)
REAL                                          :: xi(3)
INTEGER                                       :: IJK(3)
#if FV_RECONSTRUCT
INTEGER                                       :: i,j,k
REAL                                          :: distance(3)
REAL,DIMENSION(nVar_out,0:N_In,0:N_In,0:N_In) :: gradUxi_tmp
REAL,DIMENSION(nVar_out,0:N_In,0:N_In,0:N_In) :: gradUeta_tmp
#if PP_dim == 3
REAL,DIMENSION(nVar_out,0:N_In,0:N_In,0:N_In) :: gradUzeta_tmp
#endif
#endif /*FV_RECONSTRUCT*/
!===================================================================================================================================

! Calculate position in reference space when passed as physical position
IF (TrackingMethod.NE.REFMAPPING) THEN
  CALL GetPositionInRefElem(x_in,xi,ElemID+offsetElem)
ELSE
  xi = x_in
END IF
IJK(1:3) = MIN(NINT((xi(1:3)+1.)/2*(PP_N+1)),PP_N)

#if FV_RECONSTRUCT
! Transform to primitve variables
IF (nVar_out.EQ.PP_nVarPrim) THEN
  DO k=0,ZDIM(N_In); DO j=0,N_In; DO i=0,N_In
    CALL ConsToPrim(U_Prim(:,i,j,k),U_In(:,i,j,k))
  END DO; END DO; END DO! i,j,k=0,Nloc
END IF

! Get distance to fv subcell mid point in physical space
CALL Integrate_Path1D(N_In,xGP,wGP,wBary,FV_X(IJK(1)),xi(1),FV_Path_XI  (:,:,IJK(2),IJK(3),ElemID),distance(1))
CALL Integrate_Path1D(N_In,xGP,wGP,wBary,FV_X(IJK(2)),xi(2),FV_Path_ETA (:,:,IJK(1),IJK(3),ElemID),distance(2))
#if PP_dim == 3
CALL Integrate_Path1D(N_In,xGP,wGP,wBary,FV_X(IJK(3)),xi(3),FV_Path_ZETA(:,:,IJK(1),IJK(2),ElemID),distance(3))
#endif /* PP_dim == 3 */

! Get reconstructed solution at particle position
IF (nVar_out.EQ.PP_nVarPrim) THEN
  U_out = U_Prim(:,IJK(1),IJK(2),IJK(3)) + gradUxi  (:,IJK(2),IJK(3),IJK(1),ElemID) * distance(1) &
#if PP_dim == 3
                                         + gradUzeta(:,IJK(1),IJK(2),IJK(3),ElemID) * distance(3) &
#endif /* PP_dim == 3 */
                                         + gradUeta (:,IJK(1),IJK(3),IJK(2),ElemID) * distance(2)
ELSE
  CALL FV_CalcGradients(nVar,nVar_Out,ElemID,U_In,gradUxi_tmp,gradUeta_tmp,gradUzeta_tmp)
  U_out = U_In(:,IJK(1),IJK(2),IJK(3)) + gradUxi_tmp  (:,IJK(2),IJK(3),IJK(1)) * distance(1) &
#if PP_dim == 3
                                       + gradUzeta_tmp(:,IJK(1),IJK(2),IJK(3)) * distance(3) &
#endif /* PP_dim == 3 */
                                       + gradUeta_tmp (:,IJK(1),IJK(3),IJK(2)) * distance(2)
END IF
#else
! Transform to primitve variables
IF (nVar_out.EQ.PP_nVarPrim) THEN
  CALL ConsToPrim(U_out,U_in(:,IJK(1),IJK(2),IJK(3)))
ELSE
  U_Out = U_in(:,IJK(1),IJK(2),IJK(3))
END IF
#endif /*FV_RECONSTRUCT*/

END SUBROUTINE EvaluateField_FV
#endif /*FV_ENABLED*/


PPURE SUBROUTINE EvaluateFieldAtRefPos(Xi_in,NVar,N_in,U_In,NVar_out,U_Out)
!===================================================================================================================================
!> 1) interpolate DG solution to position (U_In -> U_Out(xi_in))
!===================================================================================================================================
! MODULES
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
USE MOD_Eos,                   ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: Xi_in(3)                                      !< position in reference element
INTEGER,INTENT(IN)        :: NVar                                          !< 5 (rho,u_x,u_y,u_z,e)
INTEGER,INTENT(IN)        :: N_In                                          !< usually PP_N
REAL,INTENT(IN)           :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)             !< State in Element
INTEGER,INTENT(IN)        :: NVar_out                                      !< 6 (rho,u_x,u_y,u_z,p,T)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: U_Out(1:NVar_out)                             !< Interpolated state at reference position xi_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i,j,k
REAL                      :: L_xi(3,0:N_in), L_eta_zeta
REAL                      :: Utmp(1:nVar)
!===================================================================================================================================

! 2.1) get "Vandermonde" vectors
CALL LagrangeInterpolationPolys(xi_in(1),N_in,xGP,wBary,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi_in(2),N_in,xGP,wBary,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi_in(3),N_in,xGP,wBary,L_xi(3,:))

! "more efficient" - Quote Thomas B.
Utmp(:)=0
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta=L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      Utmp = Utmp + U_IN(:,i,j,k)*L_xi(1,i)*L_Eta_Zeta
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In

! Convert to primitive variables
IF (NVar_out.EQ.PP_nVarPrim) THEN
  CALL ConsToPrim(U_out,Utmp)
ELSE
  U_out = Utmp
END IF

END SUBROUTINE EvaluateFieldAtRefPos

#if USE_EXTEND_RHS || USE_FAXEN_CORR
SUBROUTINE EvaluateFieldAndGradAtPhysPos(x_in,NVar,N_in,U_In,NVar_out,U_Out,ElemID,PartID,gradUx,gradUy,gradUz,U_RHS,UGrad_Out)
!===================================================================================================================================
!> 1) Get position within reference element (x_in -> xi=[-1,1]) by inverting the mapping
!> 2) interpolate DG solution to position (U_In -> U_Out(x_in))
!> 3) interpolate DG gradient to position (U_In -> U_GradOut(x_in))
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
USE MOD_Mesh_Vars,             ONLY: NGeo
USE MOD_Particle_Mesh_Tools,   ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars,    ONLY: wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars,    ONLY: ElemCurved,wBaryCL_NGeo1,XiCL_NGeo1
#if USE_MPI
USE MOD_Particle_Mesh_Vars,    ONLY: XCL_NGeo_Shared,dXCL_NGeo_Shared
#else
USE MOD_Particle_Mesh_Vars,    ONLY: XCL_NGeo,dXCL_NGeo
#endif /*USE_MPI*/
USE MOD_Eos,                   ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: x_in(3)                                           !< position in physical space
INTEGER,INTENT(IN)        :: NVar                                              !< 5 (rho,u_x,u_y,u_z,e)
INTEGER,INTENT(IN)        :: N_In                                              !< usually PP_N
INTEGER,INTENT(IN)        :: ElemID                                            !< Element index
REAL,INTENT(IN)           :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)                 !< State in Element
INTEGER,INTENT(IN)        :: PartID                                            !< particle ID
INTEGER,INTENT(IN)        :: NVar_out                                          !< 6 (rho,u_x,u_y,u_z,p,T)
REAL,INTENT(IN)           :: gradUx(    RHS_LIFT,0:N_in,0:N_in,0:ZDIM(N_in))   !< Gradient in x direction
REAL,INTENT(IN)           :: gradUy(    RHS_LIFT,0:N_in,0:N_in,0:ZDIM(N_in))   !< Gradient in y direction
REAL,INTENT(IN)           :: gradUz(    RHS_LIFT,0:N_in,0:N_in,0:ZDIM(N_in))   !< Gradient in z direction
REAL,INTENT(IN)           :: U_RHS(  1:RHS_NVARS,0:N_in,0:N_in,0:ZDIM(N_in))   !< du/dt, \nabla^2 u
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: U_Out(1:NVar_out)                                 !< Interpolated state at physical position x_in
REAL,INTENT(OUT)          :: UGrad_Out(RHS_GRAD,3)                             !< Interpolated gradient at physical position x_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: CNElemID,i,j,k
REAL                      :: xi(3)
REAL                      :: L_xi(3,0:PP_N), L_eta_zeta
REAL                      :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                      :: dXCL_NGeo1(1:3,1:3,0:1,0:1,0:1)
REAL                      :: Utmp(1:nVar)
!===================================================================================================================================

#if USE_MPI
ASSOCIATE( XCL_NGeo  =>  XCL_NGeo_Shared    &
         ,dXCL_NGeo  => dXCL_NGeo_Shared)
#endif /*USE_MPI*/

CALL GetRefNewtonStartValue(X_in,Xi,ElemID)

CNElemID = GetCNElemID(ElemID)

! If the element is curved, all Gauss points are required
IF (ElemCurved(CNElemID)) THEN
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID) &
                    ,NGeo,ElemID,Mode=3,PartID=PartID)
! If the element is not curved, only the corner nodes are required
ELSE
  ! fill dummy XCL_NGeo1
  XCL_NGeo1 (1:3,    0,0,0) = XCL_NGeo (1:3,     0  , 0  , 0  ,ElemID)
  XCL_NGeo1 (1:3,    1,0,0) = XCL_NGeo (1:3,    NGeo, 0  , 0  ,ElemID)
  XCL_NGeo1 (1:3,    0,1,0) = XCL_NGeo (1:3,     0  ,NGeo, 0  ,ElemID)
  XCL_NGeo1 (1:3,    1,1,0) = XCL_NGeo (1:3,    NGeo,NGeo, 0  ,ElemID)
  XCL_NGeo1 (1:3,    0,0,1) = XCL_NGeo (1:3,     0  , 0  ,NGeo,ElemID)
  XCL_NGeo1 (1:3,    1,0,1) = XCL_NGeo (1:3,    NGeo, 0  ,NGeo,ElemID)
  XCL_NGeo1 (1:3,    0,1,1) = XCL_NGeo (1:3,     0  ,NGeo,NGeo,ElemID)
  XCL_NGeo1 (1:3,    1,1,1) = XCL_NGeo (1:3,    NGeo,NGeo,NGeo,ElemID)
  ! fill dummy dXCL_NGeo1
  dXCL_NGeo1(1:3,1:3,0,0,0) = dXCL_NGeo(1:3,1:3, 0  , 0  , 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,1,0,0) = dXCL_NGeo(1:3,1:3,NGeo, 0  , 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,0,1,0) = dXCL_NGeo(1:3,1:3, 0  ,NGeo, 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,1,1,0) = dXCL_NGeo(1:3,1:3,NGeo,NGeo, 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,0,0,1) = dXCL_NGeo(1:3,1:3, 0  , 0  ,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,1,0,1) = dXCL_NGeo(1:3,1:3,NGeo, 0  ,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,0,1,1) = dXCL_NGeo(1:3,1:3, 0  ,NGeo,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,1,1,1) = dXCL_NGeo(1:3,1:3,NGeo,NGeo,NGeo,ElemID)
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo1,XiCL_NGeo1,XCL_NGeo1,dXCL_NGeo1,1,ElemID,Mode=1,PartID=PartID)
END IF

! 2.1) get "Vandermonde" vectors
CALL LagrangeInterpolationPolys(xi(1),N_in,xGP,wBary,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi(2),N_in,xGP,wBary,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi(3),N_in,xGP,wBary,L_xi(3,:))

! "more efficient" - Quote Thomas B.
Utmp(:)        = 0.
UGrad_Out(:,:) = 0.
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta = L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      Utmp                      = Utmp                      + U_IN(:,i,j,k)  * L_xi(1,i)*L_Eta_Zeta
#if USE_EXTEND_RHS
      UGrad_out(RHS_GRADVELV  ,1) = UGrad_out(RHS_GRADVELV  ,1) + gradUx(:             ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_GRADVELV  ,2) = UGrad_out(RHS_GRADVELV  ,2) + gradUy(:             ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_GRADVELV  ,3) = UGrad_out(RHS_GRADVELV  ,3) + gradUz(:             ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_dVELdt    ,1) = UGrad_out(RHS_dVELdt    ,1) + U_RHS(RHS_dVEL1dt    ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_dVELdt    ,2) = UGrad_out(RHS_dVELdt    ,2) + U_RHS(RHS_dVEL2dt    ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_dVELdt    ,3) = UGrad_out(RHS_dVELdt    ,3) + U_RHS(RHS_dVEL3dt    ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
#endif
#if USE_FAXEN_CORR
      UGrad_out(RHS_LAPLACEVEL,1) = UGrad_out(RHS_LAPLACEVEL,1) + U_RHS(RHS_LAPLACEVEL1,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_LAPLACEVEL,2) = UGrad_out(RHS_LAPLACEVEL,2) + U_RHS(RHS_LAPLACEVEL2,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_LAPLACEVEL,3) = UGrad_out(RHS_LAPLACEVEL,3) + U_RHS(RHS_LAPLACEVEL3,i,j,k) * L_xi(1,i)*L_Eta_Zeta
#endif
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In

#if USE_MPI
END ASSOCIATE
#endif

! Convert to primitve variables
IF (NVar_out.EQ.PP_nVarPrim) THEN
  CALL ConsToPrim(U_out,Utmp)
ELSE
  U_out = Utmp
END IF

END SUBROUTINE EvaluateFieldAndGradAtPhysPos


PPURE SUBROUTINE EvaluateFieldAndGradAtRefPos(Xi_in,NVar,N_in,U_In,NVar_out,U_Out,gradUx,gradUy,gradUz,U_RHS,UGrad_Out)
!===================================================================================================================================
!> 1) interpolate DG solution to position (U_In -> U_Out(xi_in))
!> 2) interpolate DG gradient to position (U_In -> U_GradOut(xi_in))
!===================================================================================================================================
! MODULES
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
USE MOD_Eos,                   ONLY: ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: Xi_in(3)                                           !< position in reference element
INTEGER,INTENT(IN)        :: NVar                                               !< 5 (rho,u_x,u_y,u_z,e)
INTEGER,INTENT(IN)        :: N_In                                               !< usually PP_N
REAL,INTENT(IN)           :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)                  !< State in Element
INTEGER,INTENT(IN)        :: NVar_out                                           !< 6 (rho,u_x,u_y,u_z,p,T)
REAL,INTENT(IN)           :: gradUx(    RHS_LIFT,0:N_in,0:N_in,0:ZDIM(N_in))    !< Gradient in x direction
REAL,INTENT(IN)           :: gradUy(    RHS_LIFT,0:N_in,0:N_in,0:ZDIM(N_in))    !< Gradient in y direction
REAL,INTENT(IN)           :: gradUz(    RHS_LIFT,0:N_in,0:N_in,0:ZDIM(N_in))    !< Gradient in z direction
REAL,INTENT(IN)           :: U_RHS(  1:RHS_NVARS,0:N_in,0:N_in,0:ZDIM(N_in))    !< du/dt, \nabla^2 u
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: U_Out(1:NVar_out)                                  !< Interpolated state at reference position xi_in
REAL,INTENT(OUT)          :: UGrad_Out(RHS_GRAD,3)                              !< Interpolated gradient at reference position x_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i,j,k
REAL                      :: L_xi(3,0:N_in), L_eta_zeta
REAL                      :: Utmp(1:nVar)
!===================================================================================================================================

! 2.1) get "Vandermonde" vectors
CALL LagrangeInterpolationPolys(xi_in(1),N_in,xGP,wBary,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi_in(2),N_in,xGP,wBary,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi_in(3),N_in,xGP,wBary,L_xi(3,:))

! "more efficient" - Quote Thomas B.
Utmp(:)       = 0.
UGrad_Out(:,:)= 0.
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta=L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      Utmp                        = Utmp                        + U_In(:,i,j,k)                * L_xi(1,i)*L_Eta_Zeta
#if USE_EXTEND_RHS
      UGrad_out(RHS_GRADVELV  ,1) = UGrad_out(RHS_GRADVELV  ,1) + gradUx(:             ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_GRADVELV  ,2) = UGrad_out(RHS_GRADVELV  ,2) + gradUy(:             ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_GRADVELV  ,3) = UGrad_out(RHS_GRADVELV  ,3) + gradUz(:             ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_dVELdt    ,1) = UGrad_out(RHS_dVELdt    ,1) + U_RHS(RHS_dVEL1dt    ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_dVELdt    ,2) = UGrad_out(RHS_dVELdt    ,2) + U_RHS(RHS_dVEL2dt    ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_dVELdt    ,3) = UGrad_out(RHS_dVELdt    ,3) + U_RHS(RHS_dVEL3dt    ,i,j,k) * L_xi(1,i)*L_Eta_Zeta
#endif
#if USE_FAXEN_CORR
      UGrad_out(RHS_LAPLACEVEL,1) = UGrad_out(RHS_LAPLACEVEL,1) + U_RHS(RHS_LAPLACEVEL1,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_LAPLACEVEL,2) = UGrad_out(RHS_LAPLACEVEL,2) + U_RHS(RHS_LAPLACEVEL2,i,j,k) * L_xi(1,i)*L_Eta_Zeta
      UGrad_out(RHS_LAPLACEVEL,3) = UGrad_out(RHS_LAPLACEVEL,3) + U_RHS(RHS_LAPLACEVEL3,i,j,k) * L_xi(1,i)*L_Eta_Zeta
#endif
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In

! Convert to primitve variables
IF (NVar_out.EQ.PP_nVarPrim) THEN
  CALL ConsToPrim(U_out,Utmp)
ELSE
  U_out = Utmp
END IF

END SUBROUTINE EvaluateFieldAndGradAtRefPos
#endif

SUBROUTINE RefElemNewton(Xi,X_In,wBaryCL_N_In,XiCL_N_In,XCL_N_In,dXCL_N_In,N_In,ElemID,Mode,PartID)
!=================================================================================================================================
!> Newton for finding the position inside the reference element [-1,1] for an arbitrary physical point
!=================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars,            ONLY: EpsMach
USE MOD_Particle_Globals
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Particle_Mesh_Vars,      ONLY:RefMappingEps
USE MOD_Particle_Vars,           ONLY:PartState,LastPartPos,PDM
USE MOD_Particle_Tracking_Vars,  ONLY:CountNbOfLostParts,NbrOfLostParticles
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: N_In
INTEGER,INTENT(IN)               :: ElemID                   !> global ID
INTEGER,INTENT(IN)               :: Mode
INTEGER,INTENT(IN),OPTIONAL      :: PartID
REAL,INTENT(IN)                  :: X_in(3)                  !> position in physical space
REAL,INTENT(IN)                  :: XiCL_N_in(0:N_In)        !> position of CL points in reference space
REAL,INTENT(IN)                  ::  XCL_N_in(3,0:N_In,0:N_in,0:N_In)   !> position of CL points in physical space
REAL,INTENT(IN)                  :: dXCL_N_in(3,3,0:N_In,0:N_in,0:N_In) !> derivation of CL points
REAL,INTENT(IN)                  :: wBaryCL_N_in(0:N_In)     !> derivation of CL points
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)               :: Xi(3)                    !> position in reference element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: Lag(1:3,0:N_In), F(1:3),Xi_Old(1:3)
INTEGER                          :: NewTonIter,i,j,k
REAL                             :: deltaXi(1:3),deltaXi2
REAL                             :: Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
REAL                             :: buff,buff2,Norm_F,Norm_F_old,lambda
INTEGER                          :: iArmijo
!===================================================================================================================================

! initial guess
CALL LagrangeInterpolationPolys(Xi(1),N_In,XiCL_N_in,wBaryCL_N_in,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),N_In,XiCL_N_in,wBaryCL_N_in,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),N_In,XiCL_N_in,wBaryCL_N_in,Lag(3,:))
! F(xi) = x(xi) - x_in
F=-x_in ! xRp
DO k=0,N_In
  DO j=0,N_In
    buff = Lag(2,j)*Lag(3,k)
    DO i=0,N_In
      F = F+XCL_N_in(:,i,j,k)*Lag(1,i)*buff !Lag(2,j)*Lag(3,k)
    END DO !l=0,N_In
  END DO !i=0,N_In
END DO !j=0,N_In

IF(ALL(ABS(F).LT.epsMach)) THEN
  deltaXi2 = 0.
ELSE
  deltaXi2 = 1. !HUGE(1.0)
END IF

Norm_F     = DOT_PRODUCT(F,F)
Norm_F_old = Norm_F
NewtonIter = 0
!abortCrit=ElemRadiusN_in(ElemID)*ElemRadiusN_in(ElemID)*RefMappingEps
DO WHILE(deltaXi2.GT.RefMappingEps .AND. NewtonIter.LT.100)
  NewtonIter = NewtonIter+1

  ! caution, dXCL_NGeo is transposed of required matrix
  Jac = 0.
  DO k=0,N_In
    DO j=0,N_In
      buff = Lag(2,j)*Lag(3,k)
      DO i=0,N_In
        buff2 = Lag(1,i)*buff
        Jac(1,1:3) = Jac(1,1:3)+dXCL_N_in(1:3,1,i,j,k)*buff2
        Jac(2,1:3) = Jac(2,1:3)+dXCL_N_in(1:3,2,i,j,k)*buff2
        Jac(3,1:3) = Jac(3,1:3)+dXCL_N_in(1:3,3,i,j,k)*buff2
      END DO !i=0,N_In
    END DO !j=0,N_In
  END DO !k=0,N_In

  ! Compute inverse of Jacobian
  sdetJac = getDet(Jac)
  IF (sdetJac.GT.0.) THEN
   sdetJac = 1./sdetJac
  ELSE !shit
   ! Newton has not converged !?!?
   IF (Mode.EQ.1) THEN
     IPWRITE(UNIT_stdOut,'(I0,A)')        ' Particle not inside of element!'
     IPWRITE(UNIT_stdOut,'(I0,A,F12.6)')  ' sdetJac     ', sdetJac
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')     ' Newton-Iter ', NewtonIter
     IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)') ' xi          ', xi(1:3)
     IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)') ' PartPos     ', X_in
     IPWRITE(UNIT_stdOut,'(I0,A,I0)')     ' ElemID      ', ElemID
     CALL Abort(__STAMP__, 'Newton in FindXiForPartPos singular. iter,sdetJac',NewtonIter,sDetJac)
   ELSE
     Xi(1) = HUGE(1.0)
     Xi(2) = Xi(1)
     Xi(3) = Xi(1)
     RETURN
   END IF
  END IF
  sJac = getInv(Jac,sdetJac)

  ! Iterate Xi using Newton step
  ! Use FAIL
  !Xi = Xi - MATMUL(sJac,F)

  ! Armijo step size control
  deltaXi  = MATMUL(sJac,F)
  deltaXi2 = DOT_PRODUCT(deltaXi,deltaXi)
  Xi_Old   = Xi

  Norm_F_old = Norm_F
  Norm_F     = Norm_F*2.
  lambda     = 1.
  iArmijo    = 1

  DO WHILE(Norm_F.GT.Norm_F_old*(1.-0.0001*lambda) .AND. iArmijo.LE.8)

    Xi = Xi_Old - lambda*deltaXI!MATMUL(sJac,F)

    ! Compute function value
    CALL LagrangeInterpolationPolys(Xi(1),N_In,XiCL_N_in,wBaryCL_N_in,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),N_In,XiCL_N_in,wBaryCL_N_in,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),N_In,XiCL_N_in,wBaryCL_N_in,Lag(3,:))
    ! F(xi) = x(xi) - x_in
    F = -x_in ! xRp
    DO k=0,N_In
      DO j=0,N_In
        buff = Lag(2,j)*Lag(3,k)
        DO i=0,N_In
          buff2 = Lag(1,i)*buff
          F     = F+XCL_N_in(:,i,j,k)*buff2
        END DO !l=0,N_In
      END DO !i=0,N_In
    END DO !j=0,N_In
    lambda  = 0.2*lambda
    iArmijo = iArmijo+1
    Norm_F  = DOT_PRODUCT(F,F)
  END DO ! Armijo iteration

  ! check xi value for plausibility
!  IF (ANY(ABS(Xi).GT.1.5)) THEN
  IF (ANY(ABS(Xi).GT.1.5) .AND. NewtonIter.GT.2) THEN
    SELECT CASE(Mode)
      CASE(1)
        IPWRITE(UNIT_stdOut,'(I0,A)')        'ERROR: Particle not inside of element, strict mode enabled!'
        IPWRITE(UNIT_stdOut,'(I0,A,I0)')     ' Newton-Iter:   ', NewtonIter
        IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)') ' xi:            ', xi(1:3)
        IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)') ' PartPos (phys):', X_in
        IF (PRESENT(PartID)) THEN
          IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' LastPos (phys):', LastPartPos(:,PartID)
          IPWRITE(UNIT_stdOut,'(I0,A,3F12.6,A,F12.6)') ' PartVel:       ', PartState(4:6,PartID),'abs:',VECNORM(PartState(4:6,PartID)**2)
          IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' PartID:        ', PartID
        END IF
        CALL Abort(__STAMP__,'Particle not inside of Element, ElemID,',ElemID)

      CASE(3)
        IPWRITE(UNIT_stdOut,'(I0,A)')        'ERROR: Particle not inside of element, particle will be deleted!'
        IF (PRESENT(PartID)) THEN
          IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' PartID:        ', PartID
        END IF
        IPWRITE(UNIT_stdOut,'(I0,A,I0)')             ' ElemID:        ', ElemID
        IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' xi:            ', xi(1:3)
        IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' PartPos (phys):', X_in
        IF (PRESENT(PartID)) THEN
          IPWRITE(UNIT_stdOut,'(I0,A,3F12.6)')         ' LastPos (phys):', LastPartPos(:,PartID)
          IPWRITE(UNIT_stdOut,'(I0,A,3F12.6,A,F12.6)') ' PartVel (+abs):', PartState(PART_VELV,PartID),&
                                                                   'abs:', VECNORM(PartState(PART_VELV,PartID))
          PDM%ParticleInside(PartID) = .FALSE.
          IF(CountNbOfLostParts) NbrOfLostParticles = NbrOfLostParticles+1
        END IF
        EXIT

      CASE DEFAULT
        EXIT
    END SELECT
  END IF

END DO !newton


END SUBROUTINE RefElemNewton


PURE FUNCTION getDet(Mat)
!=================================================================================================================================
!> compute determinant of 3x3 matrix
!=================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3)
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getDet
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
getDet=   ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * Mat(3,3) &
        + ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * Mat(3,1) &
        + ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * Mat(3,2)
END FUNCTION getDet


PURE FUNCTION getInv(Mat,sdet)
!=================================================================================================================================
!> compute inverse of 3x3 matrix, needs sDet=1/det(Mat)
!=================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3),sDet
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getInv(3,3)
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
getInv(1,1) = ( Mat(2,2) * Mat(3,3) - Mat(2,3) * Mat(3,2) ) * sdet
getInv(1,2) = ( Mat(1,3) * Mat(3,2) - Mat(1,2) * Mat(3,3) ) * sdet
getInv(1,3) = ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * sdet
getInv(2,1) = ( Mat(2,3) * Mat(3,1) - Mat(2,1) * Mat(3,3) ) * sdet
getInv(2,2) = ( Mat(1,1) * Mat(3,3) - Mat(1,3) * Mat(3,1) ) * sdet
getInv(2,3) = ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * sdet
getInv(3,1) = ( Mat(2,1) * Mat(3,2) - Mat(2,2) * Mat(3,1) ) * sdet
getInv(3,2) = ( Mat(1,2) * Mat(3,1) - Mat(1,1) * Mat(3,2) ) * sdet
getInv(3,3) = ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * sdet
END FUNCTION getInv


PPURE SUBROUTINE GetRefNewtonStartValue(X_in,Xi,ElemID)
!===================================================================================================================================
!> Returns the initial value/ guess for the Newton's algorithm
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Interpolation_Vars,      ONLY: xGP
USE MOD_Mesh_Vars,               ONLY: NGeo!,Elem_xGP,offsetElem
USE MOD_Particle_Mesh_Vars,      ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars,      ONLY: RefMappingGuess,RefMappingEps
USE MOD_Particle_Mesh_Vars,      ONLY: XiEtaZetaBasis,slenXiEtaZetaBasis
USE MOD_Particle_Mesh_Vars,      ONLY: XiCL_NGeo
USE MOD_Particle_Mesh_Vars,      ONLY: XCL_NGeo_Shared,Elem_xGP_Shared
USE MOD_Particle_Mesh_Tools,     ONLY: GetCNElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: ElemID
REAL,INTENT(IN)                :: X_in(1:3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)             :: Xi(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: Ptild(1:3),XiLinear(1:6)
REAL                          :: Winner_Dist,Dist
REAL                          :: epsOne
INTEGER                       :: CNElemID,iDir,i,j,k
REAL                          :: dX,dY,dZ
INTEGER                       :: RefMappingGuessLoc
!===================================================================================================================================

epsOne             = 1.0+RefMappingEps
RefMappingGuessLoc = RefMappingGuess

SELECT CASE(RefMappingGuessLoc)

  ! Linear interpolation
  CASE(1)
    CNElemID = GetCNElemID(ElemID)
    Ptild    = X_in - ElemBaryNGeo(:,CNElemID)
    ! plus coord system (1-3) and minus coord system (4-6)
    DO iDir = 1,6
      XiLinear(iDir) = DOT_PRODUCT(Ptild,XiEtaZetaBasis(:,iDir,CNElemID))*slenXiEtaZetaBasis(iDir,CNElemID)
    END DO
    ! compute guess as average value
    DO iDir = 1,3
      Xi(iDir) = 0.5*(XiLinear(iDir)-XiLinear(iDir+3))
    END DO
    ! limit xi to [-1,1]
    IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi = MAX(MIN(1.0d0,Xi),-1.0d0)

  ! Compute distance on Gauss Points
  CASE(2)
    Winner_Dist = SQRT(DOT_PRODUCT((x_in(:)-Elem_xGP_Shared(:,0,0,0,ElemID)),(x_in(:)-Elem_xGP_Shared(:,0,0,0,ElemID))))
    Xi(:)       = (/xGP(0),xGP(0),xGP(0)/) ! start value

    DO i = 0,PP_N; DO j = 0,PP_N; DO k = 0,PP_N
      dX = ABS(X_in(1) - Elem_xGP_Shared(1,i,j,k,ElemID))
      IF(dX.GT.Winner_Dist) CYCLE
      dY = ABS(X_in(2) - Elem_xGP_Shared(2,i,j,k,ElemID))
      IF(dY.GT.Winner_Dist) CYCLE
      dZ = ABS(X_in(3) - Elem_xGP_Shared(3,i,j,k,ElemID))
      IF(dZ.GT.Winner_Dist) CYCLE
      Dist = SQRT(dX*dX+dY*dY+dZ*dZ)
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist = Dist
        Xi(:) = (/xGP(i),xGP(j),xGP(k)/) ! start value
      END IF
    END DO; END DO; END DO

  ! Compute distance on XCL Points
  CASE(3)
    Winner_Dist = SQRT(DOT_PRODUCT((x_in(:)-XCL_NGeo_Shared(:,0,0,0,ElemID)),(x_in(:)-XCL_NGeo_Shared(:,0,0,0,ElemID))))
    Xi(:)       = (/XiCL_NGeo(0),XiCL_NGeo(0),XiCL_NGeo(0)/) ! start value

    DO i = 0,NGeo; DO j = 0,NGeo; DO k = 0,NGeo
      dX = ABS(X_in(1) - XCL_NGeo_Shared(1,i,j,k,ElemID))
      IF (dX.GT.Winner_Dist) CYCLE
      dY = ABS(X_in(2) - XCL_NGeo_Shared(2,i,j,k,ElemID))
      IF (dY.GT.Winner_Dist) CYCLE
      dZ = ABS(X_in(3) - XCL_NGeo_Shared(3,i,j,k,ElemID))
      IF (dZ.GT.Winner_Dist) CYCLE
      Dist = SQRT(dX*dX+dY*dY+dZ*dZ)
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist = Dist
        Xi(:) = (/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
      END IF
    END DO; END DO; END DO

  CASE(4)
    ! trivial guess
    xi = 0.
END SELECT

END SUBROUTINE GetRefNewtonStartValue


#if FV_RECONSTRUCT
!==================================================================================================================================
!> Calculate slopes in the inside and limit them using TVD limiters.
!> It is important to use physical distances (and not reference distances) to calculate the slope, since otherwise the
!> slope limiter can not be applied. (Scenario: inner-cell-stretching)
!> Additionally build central limited slopes for the computation of gradients used for the viscous fluxes.
!==================================================================================================================================
SUBROUTINE FV_CalcGradients(nVarIn,nVarOut,iElem,UPrim,gradUxi,gradUeta&
#if PP_dim == 3
  ,gradUzeta&
#endif
  )
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars        ,ONLY: FV_sdx_XI,FV_sdx_ETA
#if PP_dim == 3
USE MOD_FV_Vars        ,ONLY: FV_sdx_ZETA
#endif
USE MOD_FV_Limiter     ,ONLY: FV_Limiter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVarIn
INTEGER,INTENT(IN) :: nVarOut
INTEGER,INTENT(IN) :: iElem
REAL,INTENT(IN)    :: UPrim    (nVarIn ,0:PP_N,0:PP_N,0:PP_NZ) !< primitive volume solution
REAL,INTENT(OUT)   :: gradUxi  (nVarOut,0:PP_N,0:PP_N,0:PP_NZ) !< physical slope in   xi-direction (mean value)
REAL,INTENT(OUT)   :: gradUeta (nVarOut,0:PP_N,0:PP_N,0:PP_NZ) !< physical slope in  eta-direction (mean value)
#if PP_dim == 3
REAL,INTENT(OUT)   :: gradUzeta(nVarOut,0:PP_N,0:PP_N,0:PP_NZ) !< physical slope in zeta-direction (mean value)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(nVarIn,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUxi_tmp
REAL,DIMENSION(nVarIn,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUeta_tmp
#if PP_dim == 3
REAL,DIMENSION(nVarIn,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUzeta_tmp
#endif
INTEGER                                        :: l,iVar,p,q
!==================================================================================================================================
! FIXME: right know we assume that the gradients in the outer subcells are zero
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

! 3. limit
DO l=0,PP_N
  DO q=0,PP_NZ; DO p=0,PP_N
    CALL FV_Limiter(nVarOut, gradUxi_tmp  (:,p,q,l),gradUxi_tmp  (:,p,q,l+1),gradUxi  (:,p,q,l))
    CALL FV_Limiter(nVarOut, gradUeta_tmp (:,p,q,l),gradUeta_tmp (:,p,q,l+1),gradUeta (:,p,q,l))
#if PP_dim == 3
    CALL FV_Limiter(nVarOut, gradUzeta_tmp(:,p,q,l),gradUzeta_tmp(:,p,q,l+1),gradUzeta(:,p,q,l))
#endif
  END DO; END DO ! q, p
END DO ! l

END SUBROUTINE FV_CalcGradients
#endif

END MODULE MOD_Eval_xyz
