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
!> Provides particle parameters, used globally (please use EXTREMLY carefully!)
!============================================================================================================   ===================
MODULE MOD_Particle_Globals
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Memory
#if USE_MPI
USE MPI
#endif /*MPI*/
IMPLICIT NONE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                                  :: PI         = ACOS(-1.0D0)
REAL                                  :: epsMach    = epsilon(0.)
REAL                                  :: TwoEpsMach = 2.d0 * epsilon(0.)

! Keep nElems and PP_nElems separate for now
INTEGER                               :: PP_nElems
!=================================================================================================================================

INTERFACE CROSSNORM
  MODULE PROCEDURE CROSSNORM
END INTERFACE CROSSNORM

INTERFACE ALMOSTZERO
  MODULE PROCEDURE ALMOSTZERO
END INTERFACE ALMOSTZERO

INTERFACE ALMOSTEQUAL
  MODULE PROCEDURE ALMOSTEQUAL
END INTERFACE ALMOSTEQUAL

INTERFACE LESSEQUALTOLERANCE
  MODULE PROCEDURE LESSEQUALTOLERANCE
END INTERFACE

INTERFACE DOTPRODUCT
  MODULE PROCEDURE DOTPRODUCT
END INTERFACE DOTPRODUCT

INTERFACE OrthoNormVec
  MODULE PROCEDURE OrthoNormVec
END INTERFACE OrthoNormVec

INTERFACE UnitVector
  MODULE PROCEDURE UnitVector
END INTERFACE UnitVector

INTERFACE FindLinIndependentVectors
  MODULE PROCEDURE FindLinIndependentVectors
END INTERFACE FindLinIndependentVectors

INTERFACE Find2DNormIndependentVectors
  MODULE PROCEDURE Find2DNormIndependentVectors
END INTERFACE Find2DNormIndependentVectors

INTERFACE GETFREEUNIT
  MODULE PROCEDURE GETFREEUNIT
END INTERFACE GETFREEUNIT

INTERFACE RandNormal
  MODULE PROCEDURE RandNormal
END INTERFACE

INTERFACE StringBeginsWith
  MODULE PROCEDURE StringBeginsWith
END INTERFACE

PUBLIC :: PI
PUBLIC :: CROSSNORM
PUBLIC :: VECNORM
PUBLIC :: ALMOSTZERO
PUBLIC :: ALMOSTEQUAL
PUBLIC :: LESSEQUALTOLERANCE
PUBLIC :: DOTPRODUCT
PUBLIC :: UnitVector
PUBLIC :: OrthoNormVec
PUBLIC :: FindLinIndependentVectors
PUBLIC :: Find2DNormIndependentVectors
PUBLIC :: RandNormal
PUBLIC :: StringBeginsWith
!===================================================================================================================================

CONTAINS


PURE FUNCTION CROSSNORM(v1,v2)
!===================================================================================================================================
! computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
! and normalizes the vector
!===================================================================================================================================
! MODULES

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    !
REAL,INTENT(IN) :: v2(3)    !
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: CROSSNORM(3) !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: length
!===================================================================================================================================
CROSSNORM = (/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
length    = SQRT(CROSSNORM(1)*CROSSNORM(1)+CROSSNORM(2)*CROSSNORM(2)+CROSSNORM(3)*CROSSNORM(3))
CROSSNORM = CROSSNORM/length
END FUNCTION CROSSNORM


PURE FUNCTION VECNORM(v1)
!===================================================================================================================================
! Computes the Euclidean norm (length) of a vector
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    ! Vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: VECNORM  ! Euclidean norm (length) of the vector v1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
VECNORM = SQRT(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
END FUNCTION VECNORM


PURE FUNCTION ALMOSTZERO(num)
!===================================================================================================================================
! Performe an almost zero check. But ...
! Bruce Dawson quote:
! "There is no silver bullet. You have to choose wisely."
!    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
!      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
!      to your calculation. Maybe."
!    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
!      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
!      An absolute epsilon could be used if you knew exactly what number you were comparing against."
!    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
!      Good luck and God speed."
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: Num ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL         :: AlmostZero
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

AlmostZero=.FALSE.
IF(ABS(Num).LE.EpsMach) AlmostZero=.TRUE.

END FUNCTION ALMOSTZERO


PURE FUNCTION ALMOSTEQUAL(num1,num2)
!===================================================================================================================================
! Bruce Dawson quote:
! "There is no silver bullet. You have to choose wisely."
!    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
!      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
!      to your calculation. Maybe."
!    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
!      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
!      An absolute epsilon could be used if you knew exactly what number you were comparing against."
!    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
!      Good luck and God speed."
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: Num1,Num2      ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL         :: ALMOSTEQUAL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(Num1-Num2).LE.MAX(ABS(Num1),ABS(Num2))*TwoEpsMach*1.01)THEN
 ALMOSTEQUAL=.TRUE.
ELSE
 ALMOSTEQUAL=.FALSE.
END IF
END FUNCTION ALMOSTEQUAL


!===================================================================================================================================
!> Check if a <= b or a is almost equal to b via ALMOSTEQUALRELATIVE
!> Catch tolerance issues when b is only an epsilon smaller than a but the inquiry should be that they are equal
!===================================================================================================================================
PPURE LOGICAL FUNCTION LESSEQUALTOLERANCE(a,b,tol)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: a,b !< Two real numbers for comparison
REAL,INTENT(IN) :: tol !< fix for tolerance issues
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF((a.LE.b).OR.(ALMOSTEQUALRELATIVE(a,b,tol)))THEN
  LESSEQUALTOLERANCE = .TRUE.
ELSE
  LESSEQUALTOLERANCE = .FALSE.
END IF
END FUNCTION LESSEQUALTOLERANCE


PURE FUNCTION DOTPRODUCT(v1)
!===================================================================================================================================
! Computes the dot product of a vector with itself
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)       ! Input 3D vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: DOTPRODUCT  ! Dot product of v1 with itself
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DOTPRODUCT=v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)
END FUNCTION DOTPRODUCT


PURE SUBROUTINE OrthoNormVec(v1,v2,v3)
!===================================================================================================================================
!> computes orthonormal basis from a given vector v1 (v1 must be normalized)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: v1(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: v2(3), v3(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!IF(ABS(v1(3)).LT.100*EpsMach)THEN
!  v2=(/-v1(2)-v1(3) , v1(1) , v1(1)       /)
!ELSE
!  v2=(/ v1(3)       , v1(3) ,-v1(1)-v1(2) /)
!END IF
IF (ABS(v1(1)).LT.100*EpsMach) THEN
  v2 = (/ 0.    , -v1(3), v1(2)/)
ELSEIF (ABS(v1(2)).LT.100*EpsMach) THEN
  v2 = (/ -v1(3), 0.    , v1(1)/)
ELSEIF (ABS(v1(3)).LT.100*EpsMach) THEN
  v2 = (/ -v1(2), v1(1) , 0.   /)
ELSE
  v2 = (/ v1(3) , v1(3) , -v1(1)-v1(2) /)
END IF
v2    = UNITVECTOR(v2)
v3(:) = CROSSNORM(v1,v2)
END SUBROUTINE OrthoNormVec


PURE FUNCTION UNITVECTOR(v1)
!===================================================================================================================================
! compute  a unit vector from a given vector
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    !
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: UNITVECTOR(3)
REAL            :: invL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
invL=SQRT(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
invL=1./invL
UNITVECTOR=v1*invL
END FUNCTION UNITVECTOR


SUBROUTINE FindLinIndependentVectors(NormalVector, Vector1, Vector2)
!===================================================================================================================================
!> Finds two linear vectors of a normal vector around a base point
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: NormalVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT) :: Vector1(3), Vector2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Find the second vector which is in the normal plane
IF (NormalVector(1).NE.0) THEN
  Vector1(1) = (0 - NormalVector(2) - NormalVector(3)) / NormalVector(1)
  Vector1(2) = 1
  Vector1(3) = 1
ELSE IF (NormalVector(2).NE.0) THEN
  Vector1(1) = 1
  Vector1(2) = (0 - NormalVector(1) - NormalVector(3)) / NormalVector(2)
  Vector1(3) = 1
ELSE IF (NormalVector(3).NE.0) THEN
  Vector1(1) = 1
  Vector1(2) = 1
  Vector1(3) = (0 - NormalVector(1) - NormalVector(2)) / NormalVector(3)
ELSE
  CALL abort(__STAMP__,'The normal direction vector can not be (0,0,0)')
END IF

! Find the third vecord vector with the cross product
Vector2(1) = NormalVector(2)*Vector1(3) - NormalVector(3)*Vector1(2)
Vector2(2) = NormalVector(3)*Vector1(1) - NormalVector(1)*Vector1(3)
Vector2(3) = NormalVector(1)*Vector1(2) - NormalVector(2)*Vector1(1)
END SUBROUTINE FindLinIndependentVectors


SUBROUTINE Find2DNormIndependentVectors(NormalVector, Vector1, Vector2)
!===================================================================================================================================
!> Finds two linear vectors of a normal vector around a base point
!> > Inspired by Efremov, A.; Havran, V. & Seidel, H.-P., Robust and numerically stable Bézier clipping method for ray tracing NURBS
!> > surfaces, Proceedings of the 21st spring conference on Computer graphics - SCCG '05, ACM Press, 2005
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: NormalVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT) :: Vector1(3),Vector2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Find the first vector which is in the normal plane
IF (NormalVector(1).NE.0) THEN
  Vector1 = (/-NormalVector(2), NormalVector(1) , 0.             /)

ELSE IF (NormalVector(2).NE.0) THEN
  Vector1 = (/-NormalVector(2), NormalVector(1) , 0.             /)

ELSE IF (NormalVector(3).NE.0) THEN
  Vector1 = (/-NormalVector(3),               0., NormalVector(1)/)

ELSE
  CALL abort(__STAMP__,'The normal direction vector can not be (0,0,0)')
END IF

Vector1 = UNITVECTOR(Vector1)
!  The second vector has to be perpendicular to the first vector and the NormalVector
Vector2 = CROSSNORM(NormalVector,Vector1)

END SUBROUTINE Find2DNormIndependentVectors


FUNCTION GETFREEUNIT()
!===================================================================================================================================
! Get unused file unit number
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetFreeUnit ! File unit number
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: connected
!===================================================================================================================================
GetFreeUnit=55
INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
IF(connected)THEN
  DO
    GetFreeUnit=GetFreeUnit+1
    INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
    IF(.NOT.connected)EXIT
  END DO
END IF
END FUNCTION GETFREEUNIT


FUNCTION RandNormal(mean_opt,deviation_opt)
!===================================================================================================================================
! Computes probability density of the Gaussian normal distribution with standard deviation unity with the Box—Muller—Wiener (BMW)
! algorithm. Might be inefficient due to the high number of trigonometric and square function calls
! > Toral and Chakrabarti, "Generation of Gaussian distributed random numbers by using a numerical inversion method", 1993
! > Code based on https://people.sc.fsu.edu/~jburkardt/f_src/normal/normal.f90 (r4_normal_ab) and
! >               https://sukhbinder.wordpress.com/fortran-random-number-generation/
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL :: mean_opt
REAL,INTENT(IN),OPTIONAL :: deviation_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLE
REAL            :: RandNormal
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: mean
REAL            :: deviation
REAL            :: random(2),r,theta
REAL            :: Pi
!===================================================================================================================================
IF(PRESENT(mean_opt)) THEN
    mean = mean_opt
ELSE
    mean = 0.
END IF
IF(PRESENT(deviation_opt)) THEN
    deviation = deviation_opt
ELSE
    deviation = 1.
END IF

Pi = ACOS(-1.)
CALL RANDOM_NUMBER(random)

! Box—Muller—Wiener (BMW) algorithm
r          = SQRT(-2. * LOG(random(1)))
theta      = 2. * PI * random(2)
RandNormal = mean + deviation*r*COS(THETA)

END FUNCTION RandNormal


PURE LOGICAL FUNCTION StringBeginsWith(MainString,SubString)
!===================================================================================================================================
! Check if the string MainString starts with the string SubString
! Note that if one of the strings is of length zero, the result will be false and if both are zero the result will be true
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: MainString !< String in which the substring is looked for
CHARACTER(LEN=*),INTENT(IN) :: SubString  !< String which might be in MainString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: MainStringLength,SubStringLength
!===================================================================================================================================
MainStringLength = LEN(TRIM(ADJUSTL(MainString)))
SubStringLength  = LEN(TRIM(ADJUSTL(SubString)))
IF(SubStringLength.GT.0.AND.MainStringLength.GT.0)THEN
  StringBeginsWith = TRIM(MainString(1:MIN(SubStringLength,LEN(TRIM(ADJUSTL(MainString)))))).EQ.TRIM(ADJUSTL(SubString))
ELSEIF(SubStringLength.EQ.0.AND.MainStringLength.EQ.0)THEN
  StringBeginsWith = .TRUE.
ELSE
  StringBeginsWith = .FALSE.
END IF ! SubStringLength.GT.0.AND.MainStringLength.GT.0
END FUNCTION StringBeginsWith


END MODULE MOD_Particle_Globals
