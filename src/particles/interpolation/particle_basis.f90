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

!==================================================================================================================================
!> Routines to provide and evaluate basis function coefficients, or provide fast interpolation coefficients
!==================================================================================================================================
MODULE MOD_Particle_Basis
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE BuildBezierVdm
   MODULE PROCEDURE BuildBezierVdm
END INTERFACE

INTERFACE BuildBezierDMat
   MODULE PROCEDURE BuildBezierDMat
END INTERFACE

INTERFACE GetInverse
   MODULE PROCEDURE GetInverse
END INTERFACE GetInverse

INTERFACE DeCasteljauInterpolation
   MODULE PROCEDURE DeCasteljauInterpolation
END INTERFACE DeCasteljauInterpolation

INTERFACE BernSteinPolynomial
   MODULE PROCEDURE BernSteinPolynomial
END INTERFACE BernSteinPolynomial

INTERFACE ComputeBernSteinCoeff
   MODULE PROCEDURE ComputeBernSteinCoeff
END INTERFACE

PUBLIC :: BuildBezierVdm
PUBLIC :: BuildBezierDMat
PUBLIC :: DeCasteljauInterpolation
PUBLIC :: BernSteinPolynomial
PUBLIC :: GetInverse
PUBLIC :: ComputeBernSteinCoeff
!==================================================================================================================================

CONTAINS

SUBROUTINE BuildBezierVdm(N_In,xi_In,Vdm_Bezier,sVdm_Bezier)
!===================================================================================================================================
! build a 1D Vandermonde matrix using the Bezier basis functions of degree N_In
! todo: replace numerical recipes function gaussj() for calculation the inverse of V
! by a BLAS routine for better matrix conditioning
!===================================================================================================================================
! MODULES
!USE nr,                        ONLY : gaussj
USE MOD_Globals,                ONLY: abort
USE MOD_PreProc
USE MOD_Particle_Surfaces_Vars, ONLY: arrayNchooseK,FacNchooseK,BezierElevation,ElevationMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_In
REAL,INTENT(IN)    :: xi_In(0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Vdm_Bezier(0:N_In,0:N_In),sVdm_Bezier(0:N_In,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j, errorflag,IPIV(1:N_in+1),jStart,jEnd
REAL               :: dummy,eps
REAL               :: dummy_vec(0:N_In)
!===================================================================================================================================
! set NPartCurved to N_In (NGeo)
!IF(NPartCurved.NE.N_In)THEN
!  print*,"NPartCurved is not equal NGeo: Setting NPartCurved=NGeo=",N_In
!  !NPartCurved = N_In !CHANGETAG creates problems with supersampled intersection calculation
!END IF
! store the coefficients
ALLOCATE(arrayNchooseK(0:N_In,0:N_In))
ALLOCATE(FacNchooseK(0:N_In,0:N_In))
FacNchooseK(:,:) = 0.

ALLOCATE(ElevationMatrix(0:N_In+BezierElevation,0:N_In))
ElevationMatrix(:,:) = 0.


!Vandermonde on xi_In
DO i=0,N_In
  DO j=0,N_In
    ! 1.) evaluate the berstein polynomial at xi_in -> build vandermonde
    CALL BernsteinPolynomial(N_In,j,xi_in(i),Vdm_Bezier(i,j))
    ! 2.) build array with binomial coeffs for bezier clipping
    IF(i.GE.j)THEN!only calculate LU (for n >= k, else 0)
      arrayNchooseK(i,j)=REAL(CHOOSE(i,j))
      FacNchooseK(i,j) = (1.0/(2.0**REAL(i)))*ArrayNchooseK(i,j)
    ELSE
      arrayNchooseK(i,j)=0.
    END IF
  END DO !i
END DO !j

! 3.) build array with binomial coeffs (fractions) for elevation
IF(N_In+BezierElevation.GE.171) CALL Abort(&
__STAMP__&
,'Bezier elevation to polynomial degrees greater/equal 171 is forbiddon! exit.',171,REAL(N_In+BezierElevation))
ElevationMatrix(0,0) = 1.
ElevationMatrix(N_In+BezierElevation,N_In) = 1.
!print*,"BezierElevation",N_In,'+',BezierElevation
!print*,0," eps= ",0.
DO i=1,N_In+BezierElevation-1 ! from 0+1 to p_new-1 -> remove the edge points
  jStart = MAX(0,i-BezierElevation)
  jEnd   = MIN(N_In,i)
  DO j=jStart,jEnd
    !ElevationMatrix(i,j)=REAL(CHOOSE(N_In,j))*REAL(CHOOSE(BezierElevation,i-j)) / REAL(CHOOSE(N_In+BezierElevation,i))
    ElevationMatrix(i,j)=CHOOSE_large(N_In,j)*CHOOSE_large(BezierElevation,i-j) / CHOOSE_large(N_In+BezierElevation,i)
  END DO
  eps=ABS(SUM(ElevationMatrix(i,:))-1.0)
  !print*,i," eps= ",eps
  !print*,"------- next ------"
  IF(eps>1e-12) CALL Abort(&
__STAMP__&
,'The line of the elevation matrix does not sum to unity! 1-1=',0,eps)
END DO
!print*,N_In+BezierElevation," eps= ",0.
!STOP
! debug
 !print*,"BezierElevation",N_In,'+',BezierElevation
 !DO i=0,N_In+BezierElevation
  !write(*,"(I8,10F10.4)") i,ElevationMatrix(i,:)
  !write(*,"(I3)", ADVANCE = "NO")    i
  !write(*,"(12F10.4)", ADVANCE = "NO") ElevationMatrix(i,:)
  !write(*,"(F13.4)")                  SUM(ElevationMatrix(i,:))
 !END DO
 !STOP

!print*,arrayNchooseK !CHANGETAG
!Inverse of the Vandermonde
dummy_vec=0.
!print*,dummy_vec
!print*,SIZE(sVdm_Bezier),SIZE(dummy_vec)
!print*,"CALL gaussj(sVdm_Bezier,dummy_vec)"
!CALL gaussj(sVdm_Bezier,dummy_vec)
!CALL gaussj(Matrix,Vector)

!DO i=0,N_In
  !print*,Vdm_Bezier(i,:)
!END DO
! Invert A: Caution!!! From now on A=A^(-1)
sVdm_Bezier=Vdm_Bezier
CALL DGETRF(N_In+1,N_In+1,sVdm_Bezier,N_In+1,IPIV,errorflag)
IF (errorflag .NE. 0) CALL Abort(&
__STAMP__ &
,'LU factorisation of matrix crashed',999,999.)
CALL DGETRI(N_In+1,sVdm_Bezier,N_In+1,IPIV,dummy_vec,N_In+1,errorflag)
IF (errorflag .NE. 0) CALL Abort(&
__STAMP__ &
,'Solver crashed',999,999.)
!print*,"Matrix inverted"
!DO i=0,N_In
  !print*,sVdm_Bezier(i,:)
!END DO
!Matrix=MATMUL(sVdm_Bezier,Vdm_Bezier)
!print*,"A^-1*A"
!DO i=0,N_In
  !print*,Matrix(i,:)
!END DO
!print*,"(ABS(MATMUL(sVdm_Bezier,Vdm_Bezier))"
!print*,ABS(MATMUL(sVdm_Bezier,Vdm_Bezier))
!print*,"SUM: (ABS(MATMUL(sVdm_Bezier,Vdm_Bezier))"
!print*,SUM(ABS(MATMUL(sVdm_Bezier,Vdm_Bezier)))
!print*,"(N_In+1)"
!print*,(N_In+1)
!check (Vdm_Bezier)^(-1)*Vdm_Bezier := I
dummy=SUM(ABS(MATMUL(sVdm_Bezier,Vdm_Bezier)))-REAL(N_In+1)
!print*,dummy,PP_RealTolerance
!read*
IF(ABS(dummy).GT.1.E-13) CALL abort(&
__STAMP__&
,'problems in Bezier Vandermonde: check (Vdm_Bezier)^(-1)*Vdm_Bezier := I has a value of',999,dummy)
END SUBROUTINE BuildBezierVdm


SUBROUTINE BuildBezierDMat(N_In,xi_In,DMat)
!===================================================================================================================================
! build a 1D D matrix using the Bezier basis functions of degree N_In
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Particle_Surfaces_Vars, ONLY: FacNchooseK
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_In
REAL,INTENT(IN)    :: xi_In(0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: DMat(0:N_In,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j
REAL               :: XiPlus,XiMinus
!REAL               :: rtmp1,rtmp2
!===================================================================================================================================

DMat=0.
! for j=0
! (n over 0) (t+1)^0 (1-t)^n
DO i=0,N_in
  DMat(i,0) = -facNchooseK(N_in,0)*REAL(N_in)*(1-Xi_In(i))**(N_in-1)
END DO ! i=0,N_in
! for j=N_in
DO i=0,N_in
  DMat(i,N_in) = facNchooseK(N_in,N_in)*REAL(N_in)*(Xi_In(i)+1)**(N_in-1)
END DO ! i=0,N_in

! all inner values
DO j=1,N_in-1
  DO i=0,N_in
    XiPlus  =Xi_In(i)+1.0
    XiMinus =1.0-Xi_In(i)
    DMat(i,j) =facNchooseK(N_In,j)*(REAL(j)  *(XiPlus**(j-1))*XiMinus**(N_in-j)   &
                                   -REAL(N_in-j)*(XiPlus**j    )*XiMinus**(N_in-j-1) )
  END DO ! i=0,N_in
END DO ! j=0,N_in

! via diff, compare wiki, farin, etc
! caution, a 1/2 is missing in facNchooseK .oO
! DMat=0.
! DO j=0,N_in
!   DO i=0,N_in
!     XiPlus  =Xi_In(i)+1.0
!     XiMinus =1.0-Xi_In(i)
!     IF((j-1).EQ.-1)THEN
!       rTmp1=0.
!     ELSE
!       rTmp1=0.5*facNchooseK(N_in-1,j-1)*XiPlus**(j-1)*XiMinus**(N_in-j)
!     END IF
!     IF((j).EQ.N_In)THEN
!       rTmp2=0.
!     ELSE
!       rTmp2=0.5*facNchooseK(N_in-1,j)*XiPlus**(j)*XiMinus**(N_in-1-j)
!     END IF
!     DMat(i,j) = N_in*(rtmp1-rtmp2)
!   END DO ! i=0,N_in
! END DO ! j=0,N_in

END SUBROUTINE BuildBezierDmat


SUBROUTINE DeCasteljauInterpolation(N_In,xi_In,SideID,xPoint)
!===================================================================================================================================
! Computes a point in a Bezier-Surface by using the DeCasteljau alogrithm
!===================================================================================================================================
! MODULES
!USE nr,                        ONLY : gaussj
USE MOD_Globals,                ONLY: abort
USE MOD_PreProc
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_In
INTEGER,INTENT(IN)                 :: SideID
REAL,INTENT(IN)                    :: xi_In(1:2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: xPoint(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3,0:N_In,0:N_In)    :: ReducedBezierControlPoints
REAL                               :: MinusXi,Xi,MinusEta,Eta
INTEGER                            :: l,p,q,iDeCasteljau
!===================================================================================================================================


Xi=0.5*(Xi_In(1)+1.)
Eta=0.5*(Xi_in(2)+1.)
MinusXi =1.0-Xi
MinusEta=1.0-Eta

ReducedBezierControlPoints=BezierControlPoints3D(:,:,:,SideID)
l=N_In-1
DO iDeCasteljau=1,N_In
  DO q=0,l
    DO p=0,l
      ReducedBezierControlPoints(:,p,q)=MinusXi*ReducedBezierControlPoints(:,p,q  )  *MinusEta & ! A
                                       +MinusXi*ReducedBezierControlPoints(:,p,q+1)  *Eta      & ! B
                                       +     Xi*ReducedBezierControlPoints(:,p+1,q)  *MinusEta & ! C
                                       +     Xi*ReducedBezierControlPoints(:,p+1,q+1)*Eta        ! D

    END DO
  END DO
  l=l-1
END DO
xPoint=ReducedBezierControlPoints(:,0,0)

END SUBROUTINE DeCasteljauInterpolation


SUBROUTINE BernsteinPolynomial(N_in,j,x,B,NChooseK)
!===================================================================================================================================
!
!===================================================================================================================================
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!input parameters
INTEGER,INTENT(IN)      :: N_in,j ! polynomial degree, (N+1) points and j-th Bernstein polynomial is selected
REAL,INTENT(IN),OPTIONAl:: NChooseK(0:N_in,0:N_in)
REAL,INTENT(IN)         :: x      ! coordinate value in the interval [-1,1]
!-----------------------------------------------------------------------------------------------------------------------------------
!output parameters
REAL,INTENT(OUT)        :: B      ! B_N(xi)
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
!INTEGER             :: iLegendre
!REAL                :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
!REAL                :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!===================================================================================================================================
IF(PRESENT(NChooseK))THEN
  B = (1./(2**N_in))*NCHOOSEK(N_in,j)*((x+1.)**j)*((1.-x)**(N_in-j))
ELSE
  B = (1./(2**N_in))*REAL(CHOOSE(N_in,j))*((x+1.)**j)*((1.-x)**(N_in-j))
END IF

END SUBROUTINE BernsteinPolynomial


SUBROUTINE ComputeBernSteinCoeff(N_In,NChooseK)
!===================================================================================================================================
! required for deposition
! build a 1D Vandermonde matrix using the Bezier basis functions of degree N_In
! todo: replace numerical recipes function gaussj() for calculation the inverse of V
! by a BLAS routine for better matrix conditioning
!===================================================================================================================================
! MODULES
!USE nr,                        ONLY : gaussj
USE MOD_PreProc
USE MOD_Globals,                ONLY:abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_In
!REAL,INTENT(IN)    :: xi_In(0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: NchooseK(0:N_In,0:N_In)
!REAL,INTENT(OUT)   :: Vdm_Bezier(0:N_In,0:N_In),sVdm_Bezier(0:N_In,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j!,errorflag
!REAL               :: Vector(3)
!REAL               :: Matrix(0:N_In,0:N_In)
!===================================================================================================================================
! store the coefficients
!ALLOCATE(NchooseK(0:N_In,0:N_In))
NchooseK(:,:) = 0.

!Vandermonde on xi_In
DO i=0,N_In
  DO j=0,N_In
!    CALL BernsteinPolynomial(N_In,j,xi_in(i),Vdm_BernSteinN_GaussN(i,j))
    ! array with binomial coeffs for bezier clipping
    IF(i.GE.j)THEN!only calculate LU (for n >= k, else 0)
      NchooseK(i,j)=REAL(CHOOSE(i,j))
    ELSE
      NchooseK(i,j)=0.
    END IF
  END DO !i
END DO !j

END SUBROUTINE ComputeBernSteinCoeff


FUNCTION CHOOSE(N_in,k)
!===================================================================================================================================
! The binomial coefficient ( n  k ) is often read as "n choose k".
!===================================================================================================================================
USE MOD_PreProc
!USE MOD_Mesh_Vars,                ONLY:NGeo
!USE MOD_Particle_Surfaces_Vars,   ONLY:BezierElevation
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!input parameters
INTEGER,INTENT(IN) :: N_in,k
!-----------------------------------------------------------------------------------------------------------------------------------
!output parameters
INTEGER(KIND=8)            :: CHOOSE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
!===================================================================================================================================
IF((k.EQ.0).OR.(N_in.EQ.k))THEN
  CHOOSE = 1
ELSE
  !IF(NGeo+BezierElevation.GE.20)THEN
  IF(N_in.GE.21)THEN
    ! genauer
    !print*,"N_in.GE.21"
    !CHOOSE = INT(CHOOSELARGE(N_in,k),8)
    !PRINT*,CHOOSE
    ! weniger genau
    !CHOOSE = INT(FACTORIALLARGE(REAL(N_in)) / (FACTORIALLARGE(REAL(k)) * FACTORIALLARGE(REAL(N_in-k))),8)
    CHOOSE = INT(FACTORIAL_REAL(N_in) / (FACTORIAL_REAL(k) * FACTORIAL_REAL(N_in-k)))
  ELSE
    !print*,"else"
    CHOOSE = FACTORIAL(N_in) / (FACTORIAL(k) * FACTORIAL(N_in-k))
  END IF
END IF
!IF(CHOOSE.LT.0) CALL abort(__STAMP__&
  !'CHOOSE is negative. This is not allowed! ',999,REAL(CHOOSE))
END FUNCTION CHOOSE

FUNCTION FACTORIAL(N_in)
!===================================================================================================================================
! "In mathematics, the factorial of a non-negative integer n, denoted by n!, is the product of all positive integers less than or
! equal to n."
! CAUTION: THIS FUNCTION CAN ONLY HANDLE FACTORIALS UP TO N_in=20 !
!===================================================================================================================================
USE MOD_Globals
USE MOD_PreProc
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!input parameters
INTEGER,INTENT(IN) :: N_in
!-----------------------------------------------------------------------------------------------------------------------------------
!output parameters
INTEGER(KIND=8)    :: FACTORIAL
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER(KIND=8)    :: I
!===================================================================================================================================
!print*,"stop"
!stop
IF(N_in.LT.0) CALL abort(&
__STAMP__&
,'FACTORIAL of a negative integer number not allowed! ',999,REAL(N_in))
IF(N_in.EQ.0)THEN
  FACTORIAL = 1 !! debug, should be one!!!!
ELSE
  FACTORIAL = PRODUCT((/(I, I = 1, N_in)/))
END IF
IF(FACTORIAL.LT.0) CALL abort(&
__STAMP__&
,'FACTORIAL is negative. This is not allowed! ',999,REAL(FACTORIAL))
END FUNCTION FACTORIAL


FUNCTION FACTORIAL_REAL(N_in)
!===================================================================================================================================
! "In mathematics, the factorial of a non-negative integer n, denoted by n!, is the product of all positive integers less than or
! equal to n."
! CAUTION: THIS FUNCTION CAN ONLY HANDLE FACTORIALS UP TO N_in=20 !
!===================================================================================================================================
USE MOD_Globals
USE MOD_PreProc
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!input parameters
INTEGER,INTENT(IN) :: N_in
!-----------------------------------------------------------------------------------------------------------------------------------
!output parameters
REAL(KIND=8)    :: FACTORIAL_REAL
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER         :: I
!===================================================================================================================================
!print*,"stop"
!stop
IF(N_in.LT.0) CALL abort(&
__STAMP__&
,'FACTORIAL of a negative integer number not allowed! ',999,REAL(N_in))
FACTORIAL_REAL=1.
DO I=2,N_in
  FACTORIAL_REAL=FACTORIAL_REAL*REAL(I,8)
END DO
IF(FACTORIAL_REAL.LT.0) CALL abort(&
__STAMP__&
,'FACTORIAL is negative. This is not allowed! ',999,FACTORIAL_REAL)
END FUNCTION FACTORIAL_REAL

FUNCTION GetInverse(dim1,A) RESULT(Ainv)
!============================================================================================================================
! invert a matrix (dependant in LAPACK Routines)
!============================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
INTEGER, INTENT(IN) :: dim1   !size of matrix a
REAL,INTENT(IN)     :: A(dim1,dim1)
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
REAL                :: Ainv(dim1,dim1)
!----------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER            :: IPIV(dim1),INFO,lwork,i,j
REAL               :: WORK(dim1*dim1)
!============================================================================================================================
! Store A in Ainv to prevent it from being overwritten by LAPACK
  DO j=1,dim1
    DO i=1,dim1
      Ainv(i,j) = A(i,j)
    END DO
  END DO

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  CALL DGETRF(dim1, dim1, Ainv, dim1, IPIV, INFO)

  IF (INFO /= 0) THEN
    CALL abort(&
__STAMP__&
,' Matrix is numerically singular!')
  END IF

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  lwork=dim1*dim1
  CALL DGETRI(dim1, Ainv, dim1, IPIV, WORK, lwork , INFO)

  IF (INFO /= 0) THEN
    CALL abort(&
__STAMP__&
,' Matrix inversion failed!')
  END IF
END FUNCTION GetInverse

FUNCTION CHOOSE_large(N_in,k)
!===================================================================================================================================
! The binomial coefficient ( n  k ) is often read as "n choose k".
!===================================================================================================================================
USE MOD_PreProc
!USE MOD_Mesh_Vars,                ONLY:NGeo
!USE MOD_Particle_Surfaces_Vars,   ONLY:BezierElevation
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!input parameters
INTEGER,INTENT(IN) :: N_in,k
!-----------------------------------------------------------------------------------------------------------------------------------
!output parameters
REAL(KIND=8)       :: CHOOSE_large
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
!===================================================================================================================================
IF((k.EQ.0).OR.(N_in.EQ.k))THEN
  CHOOSE_large = 1
ELSE
  !IF(NGeo+BezierElevation.GE.20)THEN
  IF(N_in.GE.21)THEN
    ! genauer
    !print*,"N_in.GE.21"
    !CHOOSE_large = INT(CHOOSE_largeLARGE(N_in,k),8)
    !PRINT*,CHOOSE_large
    ! weniger genau
    !CHOOSE_large = INT(FACTORIALLARGE(REAL(N_in)) / (FACTORIALLARGE(REAL(k)) * FACTORIALLARGE(REAL(N_in-k))),8)
    CHOOSE_large = FACTORIAL_REAL(N_in) / (FACTORIAL_REAL(k) * FACTORIAL_REAL(N_in-k))
  ELSE
    !print*,"else"
    CHOOSE_large = FACTORIAL(N_in) / (FACTORIAL(k) * FACTORIAL(N_in-k))
  END IF
END IF
!IF(CHOOSE_large.LT.0) CALL abort(__STAMP__&
  !'CHOOSE_large is negative. This is not allowed! ',999,REAL(CHOOSE_large))
END FUNCTION CHOOSE_large

END MODULE MOD_Particle_Basis
