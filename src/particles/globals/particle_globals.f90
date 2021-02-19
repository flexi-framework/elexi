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
#if USE_MPI
USE MPI
#endif /*MPI*/
IMPLICIT NONE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                                  :: PI
REAL                                  :: epsMach    = epsilon(0.)
REAL                                  :: TwoEpsMach = 2.d0 * epsilon(0.)

! Keep nElems and PP_nElems separate for now
INTEGER                               :: PP_nElems

!- LOADBALANCE ---------------------------------------------------------------------------------------------------------------------
REAL                         :: WallTime                              !> Wall time needed by a simulation (is not reset by
                                                                      !> performing a load balance step, only by user restart)
REAL                         :: InitializationWallTime                !> Wall time needed to initialize a simulation (or
                                                                      !> re-initialize a simulation by performing a load balance
                                                                      !>  step)
REAL                         :: SimulationEfficiency                  !> relates the simulated time to the used CPUh (SIMULATION
                                                                      !> TIME PER CALCULATION in [s]/[CPUh])
REAL                         :: PID                                   !> Performance index: (CalcTimeEnd-CalcTimeStart)*nProcessors/
                                                                      !> (nGlobalElems*(PP_N+1)**3*iter_loc)

!=================================================================================================================================

INTERFACE CROSSNORM
  MODULE PROCEDURE CROSSNORM
END INTERFACE CROSSNORM

INTERFACE AlmostZero
  MODULE PROCEDURE AlmostZero
END INTERFACE AlmostZero

INTERFACE AlmostEqual
  MODULE PROCEDURE AlmostEqual
END INTERFACE AlmostEqual

INTERFACE DOTPRODUCT
  MODULE PROCEDURE DOTPRODUCT
END INTERFACE DOTPRODUCT

INTERFACE UnitVector
  MODULE PROCEDURE UnitVector
END INTERFACE UnitVector

INTERFACE OrthoNormVec
  MODULE PROCEDURE OrthoNormVec
END INTERFACE OrthoNormVec

INTERFACE GETFREEUNIT
  MODULE PROCEDURE GETFREEUNIT
END INTERFACE GETFREEUNIT

INTERFACE RandNormal
  MODULE PROCEDURE RandNormal
END INTERFACE

INTERFACE StringBeginsWith
  MODULE PROCEDURE StringBeginsWith
END INTERFACE

INTERFACE
  SUBROUTINE processmemusage(memUsed,memAvail,memTotal) BIND(C, name='processmemusage')
    USE ISO_C_BINDING,   ONLY : c_double
    real(c_double) :: memUsed
    real(c_double) :: memAvail
    real(c_double) :: memTotal
  END SUBROUTINE processmemusage
END INTERFACE

PUBLIC :: PI
PUBLIC :: CROSSNORM
PUBLIC :: VECNORM
PUBLIC :: AlmostZero
PUBLIC :: AlmostEqual
PUBLIC :: DOTPRODUCT
PUBLIC :: UnitVector
PUBLIC :: RandNormal
PUBLIC :: StringBeginsWith
PUBLIC :: processmemusage
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


PURE FUNCTION AlmostZero(Num)
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

END FUNCTION AlmostZero


PURE FUNCTION AlmostEqual(Num1,Num2)
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
END FUNCTION AlmostEqual


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
REAL,INTENT(IN) :: v1(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: v2(3), v3(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(v1(3)).LT.100*EpsMach)THEN
  v2=(/-v1(2)-v1(3) , v1(1) , v1(1)       /)
ELSE
  v2=(/ v1(3)       , v1(3) ,-v1(1)-v1(2) /)
END IF
v2    = UNITVECTOR(v2)
v3(:) = CROSSNORM(v1,v2)
END SUBROUTINE OrthoNormVec


PURE FUNCTION UNITVECTOR(v1)
!===================================================================================================================================
! compute  a unit vector from a given vector
!===================================================================================================================================
! MODULES
USE Mod_Globals
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
