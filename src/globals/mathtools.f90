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
!> Routines providing general math functions
!==================================================================================================================================
MODULE MOD_Mathtools
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE INVERSE
   MODULE PROCEDURE INVERSE
END INTERFACE

INTERFACE INVERSE_LU
  MODULE PROCEDURE INVERSE_LU
END INTERFACE

INTERFACE CROSS
  MODULE PROCEDURE CROSS
END INTERFACE CROSS

PUBLIC::INVERSE
PUBLIC::INVERSE_LU
PUBLIC::CROSS
PUBLIC::GlobalVectorDotProduct
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes matrix inverse using LAPACK
!> Input matrix should be a square matrix
!==================================================================================================================================
FUNCTION INVERSE(A) RESULT(AINV)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: A(:,:)                      !< input matrix
REAL             :: AINV(SIZE(A,1),SIZE(A,2))   !< result: inverse of A
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI
! LOCAL VARIABLES
REAL    :: sdet
REAL    :: work(SIZE(A,1))  ! work array for LAPACK
INTEGER :: ipiv(SIZE(A,1))  ! pivot indices
INTEGER :: n,info
!==================================================================================================================================
! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

IF (n.EQ.2) THEN
  sdet=1./(A(1,1) * A(2,2) - A(1,2)*A(2,1))
  AINV(1,1) = (  A(2,2) ) * sdet
  AINV(1,2) = (- A(1,2) ) * sdet
  AINV(2,1) = (- A(2,1) ) * sdet
  AINV(2,2) = (  A(1,1) ) * sdet
ELSE IF (n.EQ.3) THEN
  sdet = 1./ (  ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) * A(3,3) &
              + ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) * A(3,1) &
              + ( A(1,3) * A(2,1) - A(1,1) * A(2,3) ) * A(3,2))
  AINV(1,1) = ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) * sdet
  AINV(1,2) = ( A(1,3) * A(3,2) - A(1,2) * A(3,3) ) * sdet
  AINV(1,3) = ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) * sdet
  AINV(2,1) = ( A(2,3) * A(3,1) - A(2,1) * A(3,3) ) * sdet
  AINV(2,2) = ( A(1,1) * A(3,3) - A(1,3) * A(3,1) ) * sdet
  AINV(2,3) = ( A(1,3) * A(2,1) - A(1,1) * A(2,3) ) * sdet
  AINV(3,1) = ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) * sdet
  AINV(3,2) = ( A(1,2) * A(3,1) - A(1,1) * A(3,2) ) * sdet
  AINV(3,3) = ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) * sdet
ELSE
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  CALL DGETRF(n, n, Ainv, n, ipiv, info)

  IF(info.NE.0)THEN
    STOP 'Matrix is numerically singular!'
  END IF

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

  IF(info.NE.0)THEN
    STOP 'Matrix inversion failed!'
  END IF
END IF
END FUNCTION INVERSE


!==================================================================================================================================
!> Computes matrix inverse using Doolittle LU factorization for Ax=b
!> Input matrix should be a square matrix
!==================================================================================================================================
FUNCTION INVERSE_LU(A) RESULT(AINV)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: A(:,:)                      !< Input matrix
REAL             :: AINV(SIZE(A,1),SIZE(A,2))   !< Result: inverse of A
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, DIMENSION(SIZE(A,1),SIZE(A,2)) :: U,L,A2  !< Upper and Lower part of A and copy
REAL, DIMENSION(SIZE(A,1))           :: b,d,x   !< RHS and aux vectors
INTEGER                              :: i,j,k   !< Loop indices
INTEGER                              :: n       !< Size of first dimension of A, assume that A is a square matirx
REAL                                 :: c       !< Auxiliary coefficient
!==================================================================================================================================
n = size(A,1)

IF (n.LE.3) THEN
  ! Call other routine that calculates the exact inverse
  AINV = INVERSE(A)
ELSE
  ! 0.) Store A in A2 to prevent it from being overwritten
  ! and init lower and upper matrices L and U and RHS b
  A2 = A
  L  = 0.
  U  = 0.
  b  = 0.

  ! 1.) Forward elimination
  DO k=1, n-1
    DO i=k+1,n
      c      = A2(i,k)/A2(k,k)
      L(i,k) = c
      DO j=k+1,n
        A2(i,j) = A2(i,j)-c*A2(k,j)
      END DO
    END DO
  END DO

  ! 2.) Set lower L and upper U matrices
  ! L matrix: A2 matrix of the elimination coefficient and diagonal elements are unity
  DO i=1,n
    L(i,i) = 1.
  END DO

  ! U matrix: Upper part of A2
  DO j=1,n
    DO i=1,j
      U(i,j) = A2(i,j)
    END DO
  END DO

  ! 3.) Columns of the inverse matrix AINV
  DO k=1,n
    b(k)=1.
    d(1) = b(1)

    ! 4.) Forward substitution: Solve L*d = b
    DO i=2,n
      d(i)=b(i)
      DO j=1,i-1
        d(i) = d(i) - L(i,j)*d(j)
      END DO
    END DO

    ! 5.) Backward substitution: Solve U*x = d
    x(n)=d(n)/U(n,n)
    DO i = n-1,1,-1
      x(i) = d(i)
      DO j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
      END DO
      x(i) = x(i)/u(i,i)
    END DO

    ! 6.) Copy x into AINV
    DO i=1,n
      AINV(i,k) = x(i)
    END DO
    b(k)=0.
  END DO
END IF
END FUNCTION INVERSE_LU


!FUNCTION INV(A) RESULT(AINV)
!!===================================================================================================================================
!! Computes matrix inverse using lapack
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)  :: A(:,:)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL             :: AINV(SIZE(A,1),SIZE(A,2))
!!-----------------------------------------------------------------------------------------------------------------------------------
!! External procedures defined in LAPACK
!EXTERNAL DGETRF
!EXTERNAL DGETRI
!! LOCAL VARIABLES
!REAL    :: work(SIZE(A,1))  ! work array for lapack
!INTEGER :: ipiv(SIZE(A,1))  ! pivot indices
!INTEGER :: n,info
!!===================================================================================================================================
!! Store A in Ainv to prevent it from being overwritten by LAPACK
!Ainv = A
!n = size(A,1)

!! DGETRF computes an LU factorization of a general M-by-N matrix A
!! using partial pivoting with row interchanges.
!CALL DGETRF(n, n, Ainv, n, ipiv, info)

!IF(info.NE.0)THEN
!    CALL abort(&
!__STAMP__&
!,' Matrix is numerically singular!')
!END IF

!! DGETRI computes the inverse of a matrix using the LU factorization
!! computed by DGETRF.
!CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

!IF(info.NE.0)THEN
!    CALL abort(&
!__STAMP__&
!,' Matrix inversion failed!')
!END IF
!END FUNCTION INV

!SUBROUTINE INV33(M,MInv,detM)
!!===================================================================================================================================
!! Computes the inverse of a 3x3 matrix
!!===================================================================================================================================
!! MODULES
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)     :: M(3,3)  ! ?
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)    :: MInv(3,3),detM  ! ?
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!!===================================================================================================================================
!detM =   M(1,1)*M(2,2)*M(3,3)  &
!       - M(1,1)*M(2,3)*M(3,2)  &
!       - M(1,2)*M(2,1)*M(3,3)  &
!       + M(1,2)*M(2,3)*M(3,1)  &
!       + M(1,3)*M(2,1)*M(3,2)  &
!       - M(1,3)*M(2,2)*M(3,1)

!IF(ABS(detM).LE.1.E-12*SUM(ABS(M)))THEN
!   MInv=0.
!   detM=0.
!   RETURN
!END IF

!MInv(1,1) =  (M(2,2)*M(3,3)-M(2,3)*M(3,2))
!MInv(2,1) = -(M(2,1)*M(3,3)-M(2,3)*M(3,1))
!MInv(3,1) =  (M(2,1)*M(3,2)-M(2,2)*M(3,1))
!MInv(1,2) = -(M(1,2)*M(3,3)-M(1,3)*M(3,2))
!MInv(2,2) =  (M(1,1)*M(3,3)-M(1,3)*M(3,1))
!MInv(3,2) = -(M(1,1)*M(3,2)-M(1,2)*M(3,1))
!MInv(1,3) =  (M(1,2)*M(2,3)-M(1,3)*M(2,2))
!MInv(2,3) = -(M(1,1)*M(2,3)-M(1,3)*M(2,1))
!MInv(3,3) =  (M(1,1)*M(2,2)-M(1,2)*M(2,1))
!MInv=MInv/detM

!END SUBROUTINE INV33


! !===================================================================================================================================
! !> Bruce Dawson quote:
! !> "There is no silver bullet. You have to choose wisely."
! !>    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
! !>      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
! !>      to your calculation. Maybe."
! !>    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
! !>      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
! !>      An absolute epsilon could be used if you knew exactly what number you were comparing against."
! !>    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
! !>      Good luck and God speed."
! !>
! !>      NOTE: The functions below are implemented as preprocessor macros, which are by definition
! !>            inlined and are thus beneficial in terms of performance and accuracy as they do not
! !>            depend on the data type.
! !>
! !===================================================================================================================================
! PURE FUNCTION ALMOSTEQUALRELATIVE(x,y,tol)
! ! MODULES
! IMPLICIT NONE
! !----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT/OUTPUT VARIABLES
! REAL,INTENT(IN) :: x                !< (IN)  first scalar to be compared
! REAL,INTENT(IN) :: y                !< (IN)  second scalar to be compared
! REAL,INTENT(IN) :: tol              !< (IN) relative epsilon value as input
! LOGICAL         :: ALMOSTEQUALRELATIVE
! !==================================================================================================================================
! ALMOSTEQUALRELATIVE=(ABS(x-y).LE.MAX(ABS(x),ABS(y))*tol)
! END FUNCTION ALMOSTEQUALRELATIVE
! PURE FUNCTION ALMOSTEQUALABSOLUTE(x,y,tol)
! ! MODULES
! IMPLICIT NONE
! !----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT/OUTPUT VARIABLES
! REAL,INTENT(IN) :: x                !< (IN)  first scalar to be compared
! REAL,INTENT(IN) :: y                !< (IN)  second scalar to be compared
! REAL,INTENT(IN) :: tol              !< (IN) relative epsilon value as input
! LOGICAL         :: ALMOSTEQUALABSOLUTE
! !==================================================================================================================================
! ALMOSTEQUALRELATIVE=(ABS(x-y).LE.tol)
!
! END FUNCTION ALMOSTEQUALABSOLUTE


!==================================================================================================================================
!> computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
!==================================================================================================================================
PPURE FUNCTION CROSS(v1,v2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    !< input vector 1
REAL,INTENT(IN) :: v2(3)    !< input vector 2
REAL            :: CROSS(3) !< cross product of vectors
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CROSS=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS

!===================================================================================================================================
!> Computes the global dot Product for vectors a and b: resu=a.b. Each processor computes the local dot product, then an allreduce
!> is called.
!===================================================================================================================================
SUBROUTINE GlobalVectorDotProduct(A,B,nDOFProc,Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: A(nDOFProc)
REAL,INTENT(IN)    :: B(nDOFProc)
INTEGER,INTENT(IN) :: nDOFProc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
!===================================================================================================================================
Resu=0.
DO i=1,nDOFProc
  Resu=Resu + A(i)*B(i)
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif

END SUBROUTINE GlobalVectorDotProduct

END MODULE MOD_Mathtools
