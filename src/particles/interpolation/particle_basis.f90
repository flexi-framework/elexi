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

SUBROUTINE ComputeBernSteinCoeff(N_In,NChooseK)
!===================================================================================================================================
! required for deposition
! build a 1D Vandermonde matrix using the Bezier basis functions of degree N_In
! todo: replace numerical recipes function gaussj() for calculation the inverse of V
! by a BLAS routine for better matrix conditioning
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,                ONLY:abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_In
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: NchooseK(0:N_In,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j
!===================================================================================================================================
! store the coefficients
NchooseK(:,:) = 0.

!Vandermonde on xi_In
DO i = 0,N_In
  DO j = 0,N_In
    !only calculate LU (for n >= k, else 0)
    IF (i.GE.j) THEN
      NchooseK(i,j) = REAL(CHOOSE(i,j))
    ELSE
      NchooseK(i,j) = 0.
    END IF
  END DO !i
END DO !j

END SUBROUTINE ComputeBernSteinCoeff


SUBROUTINE BuildBezierVdm(N_In,xi_In,Vdm_Bezier,sVdm_Bezier)
!===================================================================================================================================
! build a 1D Vandermonde matrix using the Bezier basis functions of degree N_In
! todo: replace numerical recipes function gaussj() for calculation the inverse of V
! by a BLAS routine for better matrix conditioning
!===================================================================================================================================
! MODULES
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

! store the coefficients
ALLOCATE(arrayNchooseK(0:N_In,0:N_In)   &
        ,FacNchooseK  (0:N_In,0:N_In))
FacNchooseK(:,:) = 0.

ALLOCATE(ElevationMatrix(0:N_In+BezierElevation,0:N_In))
ElevationMatrix(:,:) = 0.

! Vandermonde on xi_In
DO i = 0,N_In
  DO j = 0,N_In
    ! 1.) evaluate the Bernstein polynomial at xi_in -> build Vandermonde
    CALL BernsteinPolynomial(N_In,j,xi_in(i),Vdm_Bezier(i,j))

    ! 2.) build array with binomial coeffs for bezier clipping
    !>> only calculate LU (for n >= k, else 0)
    IF (i.GE.j) THEN
      arrayNchooseK(i,j) = REAL(CHOOSE(i,j))
      FacNchooseK(  i,j) = (1.0/(2.0**REAL(i)))*ArrayNchooseK(i,j)
    ELSE
      arrayNchooseK(i,j) = 0.
    END IF
  END DO !i
END DO !j

! 3.) build array with binomial coeffs (fractions) for elevation
IF (N_In+BezierElevation.GE.171) &
  CALL ABORT(__STAMP__,'Bezier elevation to polynomial degrees greater/equal 171 is forbiddon! exit.',171,REAL(N_In+BezierElevation))

ElevationMatrix(0                   ,0)    = 1.
ElevationMatrix(N_In+BezierElevation,N_In) = 1.

! from 0+1 to p_new-1 -> remove the edge points
DO i = 1,N_In+BezierElevation-1
  jStart = MAX(0,i-BezierElevation)
  jEnd   = MIN(N_In,i)

  DO j = jStart,jEnd
    ElevationMatrix(i,j)=CHOOSE_large(N_In,j)*CHOOSE_large(BezierElevation,i-j) / CHOOSE_large(N_In+BezierElevation,i)
  END DO

  eps=ABS(SUM(ElevationMatrix(i,:))-1.0)

  IF (eps>1e-12) &
    CALL ABORT(__STAMP__,'The line of the elevation matrix does not sum to unity! 1-1=',0,eps)
END DO

dummy_vec=0.

! Invert A: Caution!!! From now on A=A^(-1)
sVdm_Bezier=Vdm_Bezier

CALL DGETRF(N_In+1,N_In+1,sVdm_Bezier,N_In+1,IPIV,errorflag)

IF (errorflag .NE. 0) &
  CALL ABORT(__STAMP__,'LU factorisation of matrix crashed',999,999.)

CALL DGETRI(N_In+1,sVdm_Bezier,N_In+1,IPIV,dummy_vec,N_In+1,errorflag)

IF (errorflag .NE. 0) &
  CALL ABORT(__STAMP__,'Solver crashed',999,999.)

dummy=SUM(ABS(MATMUL(sVdm_Bezier,Vdm_Bezier)))-REAL(N_In+1)

IF(ABS(dummy).GT.1.E-13) &
  CALL ABORT(__STAMP__,'problems in Bezier Vandermonde: check (Vdm_Bezier)^(-1)*Vdm_Bezier := I has a value of',999,dummy)

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
DO i = 0,N_in
  DMat(i,0)    = -facNchooseK(N_in,0  )*REAL(N_in)*(1-Xi_In(i))**(N_in-1)
END DO ! i=0,N_in

! for j=N_in
DO i = 0,N_in
  DMat(i,N_in) = facNchooseK(N_in,N_in)*REAL(N_in)*(Xi_In(i)+1)**(N_in-1)
END DO ! i=0,N_in

! all inner values
DO j = 1,N_in-1
  DO i = 0,N_in
    XiPlus  = Xi_In(i)+1.0
    XiMinus = 1.0-Xi_In(i)
    DMat(i,j) = facNchooseK(N_In   ,j)*(REAL(j     )*(XiPlus**(j-1))*XiMinus**(N_in-j  )   &
                                       -REAL(N_in-j)*(XiPlus** j   )*XiMinus**(N_in-j-1))
  END DO ! i=0,N_in
END DO ! j=0,N_in

END SUBROUTINE BuildBezierDmat


SUBROUTINE DeCasteljauInterpolation(N_In,xi_In,SideID,xPoint)
!===================================================================================================================================
! Computes a point in a Bezier-Surface by using the DeCasteljau alogrithm
!===================================================================================================================================
! MODULES
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

Xi   = 0.5*(Xi_In(1)+1.)
Eta  = 0.5*(Xi_in(2)+1.)
MinusXi  = 1.0-Xi
MinusEta = 1.0-Eta

ReducedBezierControlPoints=BezierControlPoints3D(:,:,:,SideID)
l = N_In-1
DO iDeCasteljau = 1,N_In
  DO q = 0,l
    DO p = 0,l
      ReducedBezierControlPoints(:,p,q) = MinusXi*ReducedBezierControlPoints(:,p,q    )*MinusEta & ! A
                                        + MinusXi*ReducedBezierControlPoints(:,p,q+1  )*Eta      & ! B
                                        +      Xi*ReducedBezierControlPoints(:,p+1,q  )*MinusEta & ! C
                                        +      Xi*ReducedBezierControlPoints(:,p+1,q+1)*Eta        ! D

    END DO
  END DO
  l = l-1
END DO
xPoint = ReducedBezierControlPoints(:,0,0)

END SUBROUTINE DeCasteljauInterpolation


SUBROUTINE BernsteinPolynomial(N_in,j,x,B,NChooseK)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: N_in,j ! polynomial degree, (N+1) points and j-th Bernstein polynomial is selected
REAL,INTENT(IN),OPTIONAl:: NChooseK(0:N_in,0:N_in)
REAL,INTENT(IN)         :: x      ! coordinate value in the interval [-1,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)        :: B      ! B_N(xi)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(PRESENT(NChooseK))THEN
  B = (1./(2**N_in))*NCHOOSEK(N_in,j)*((x+1.)**j)*((1.-x)**(N_in-j))
ELSE
  B = (1./(2**N_in))*REAL(CHOOSE(N_in,j))*((x+1.)**j)*((1.-x)**(N_in-j))
END IF

END SUBROUTINE BernsteinPolynomial


FUNCTION GetInverse(dim1,A) RESULT(Ainv)
!============================================================================================================================
! invert a matrix (dependant in LAPACK Routines)
!============================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: dim1   !size of matrix a
REAL,INTENT(IN)     :: A(dim1,dim1)
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: Ainv(dim1,dim1)
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: IPIV(dim1),INFO,lwork,i,j
REAL               :: WORK(dim1*dim1)
!============================================================================================================================

! Store A in Ainv to prevent it from being overwritten by LAPACK
DO j = 1,dim1
  DO i = 1,dim1
    Ainv(i,j) = A(i,j)
  END DO
END DO

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
CALL DGETRF(dim1, dim1, Ainv, dim1, IPIV, INFO)

IF (INFO /= 0) &
  CALL abort(__STAMP__,' Matrix is numerically singular!')

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
lwork = dim1*dim1
CALL DGETRI(dim1, Ainv, dim1, IPIV, WORK, lwork , INFO)

IF (INFO /= 0) &
  CALL abort(__STAMP__,' Matrix inversion failed!')

END FUNCTION GetInverse



FUNCTION CHOOSE(N_in,k)
!===================================================================================================================================
! The binomial coefficient ( n  k ) is often read as "n choose k".
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_in,k
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8)            :: CHOOSE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF((k.EQ.0).OR.(N_in.EQ.k))THEN
  CHOOSE = 1
ELSE
  IF(N_in.GE.21)THEN
    ! Better for large N_in
    CHOOSE = INT(FACTORIAL_REAL(N_in) / (FACTORIAL_REAL(k) * FACTORIAL_REAL(N_in-k)))
  ELSE
    CHOOSE = FACTORIAL(N_in) / (FACTORIAL(k) * FACTORIAL(N_in-k))
  END IF
END IF

END FUNCTION CHOOSE


FUNCTION CHOOSE_large(N_in,k)
!===================================================================================================================================
! The binomial coefficient ( n  k ) is often read as "n choose k".
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_in,k
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(KIND=8)       :: CHOOSE_large
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF ((k.EQ.0).OR.(N_in.EQ.k)) THEN
  CHOOSE_large = 1
ELSE
  IF (N_in.GE.21) THEN
! Better for large N_in
    CHOOSE_large = FACTORIAL_REAL(N_in) / (FACTORIAL_REAL(k) * FACTORIAL_REAL(N_in-k))
  ELSE
    CHOOSE_large = FACTORIAL(N_in) / (FACTORIAL(k) * FACTORIAL(N_in-k))
  END IF
END IF

END FUNCTION CHOOSE_large


FUNCTION FACTORIAL(N_in)
!===================================================================================================================================
! "In mathematics, the factorial of a non-negative integer n, denoted by n!, is the product of all positive integers less than or
! equal to n."
! CAUTION: THIS FUNCTION CAN ONLY HANDLE FACTORIALS UP TO N_in=20 !
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8)    :: FACTORIAL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)    :: I
!===================================================================================================================================
IF (N_in.LT.0) &
  CALL abort(__STAMP__,'FACTORIAL of a negative integer number not allowed! ',999,REAL(N_in))

IF (N_in.EQ.0) THEN
   ! debug, should be one!
  FACTORIAL = 1
ELSE
  FACTORIAL = PRODUCT((/(I, I = 1, N_in)/))
END IF

IF (FACTORIAL.LT.0) &
  CALL abort(__STAMP__,'FACTORIAL is negative. This is not allowed! ',999,REAL(FACTORIAL))

END FUNCTION FACTORIAL


FUNCTION FACTORIAL_REAL(N_in)
!===================================================================================================================================
! "In mathematics, the factorial of a non-negative integer n, denoted by n!, is the product of all positive integers less than or
! equal to n."
! CAUTION: THIS FUNCTION CAN ONLY HANDLE FACTORIALS UP TO N_in=20 !
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(KIND=8)    :: FACTORIAL_REAL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: I
!===================================================================================================================================

IF (N_in.LT.0) &
  CALL abort(__STAMP__,'FACTORIAL of a negative integer number not allowed! ',999,REAL(N_in))

FACTORIAL_REAL=1.

DO I = 2,N_in
  FACTORIAL_REAL=FACTORIAL_REAL*REAL(I,8)
END DO

IF(FACTORIAL_REAL.LT.0) &
  CALL abort(__STAMP__,'FACTORIAL is negative. This is not allowed! ',999,FACTORIAL_REAL)

END FUNCTION FACTORIAL_REAL

END MODULE MOD_Particle_Basis
