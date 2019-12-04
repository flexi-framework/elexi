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

INTERFACE InitParticleGlobals
  MODULE PROCEDURE InitParticleGlobals
END INTERFACE InitParticleGlobals

INTERFACE CROSS
  MODULE PROCEDURE CROSS
END INTERFACE CROSS

INTERFACE CROSSNORM
  MODULE PROCEDURE CROSSNORM
END INTERFACE CROSSNORM

INTERFACE AlmostZero
  MODULE PROCEDURE AlmostZero
END INTERFACE AlmostZero

INTERFACE AlmostEqual
  MODULE PROCEDURE AlmostEqual
END INTERFACE AlmostEqual

INTERFACE UnitVector
  MODULE PROCEDURE UnitVector
END INTERFACE UnitVector

INTERFACE GETFREEUNIT
  MODULE PROCEDURE GETFREEUNIT
END INTERFACE GETFREEUNIT

INTERFACE RandNormal
  MODULE PROCEDURE RandNormal
END INTERFACE

PUBLIC :: PI
PUBLIC :: CROSS
PUBLIC :: CROSSNORM
PUBLIC :: AlmostZero
PUBLIC :: AlmostEqual
PUBLIC :: UnitVector
PUBLIC :: InitParticleGlobals
PUBLIC :: RandNormal
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleGlobals()
!===================================================================================================================================
! Pre-compute required constants
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GLOBALS ...'

PI=ACOS(-1.0D0)

END SUBROUTINE InitParticleGlobals

PURE FUNCTION CROSS(v1,v2)
!===================================================================================================================================
! computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
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
REAL            :: CROSS(3) !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CROSS=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS

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
CROSSNORM=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
length=SQRT(CROSSNORM(1)*CROSSNORM(1)+CROSSNORM(2)*CROSSNORM(2)+CROSSNORM(3)*CROSSNORM(3))
CROSSNORM=CROSSNORM/length
END FUNCTION CROSSNORM

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
FUNCTION AlmostZero(Num)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL            :: Num ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL         :: AlmostZero
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

AlmostZero=.FALSE.
IF(ABS(Num).LE.EpsMach) AlmostZero=.TRUE.

END FUNCTION AlmostZero

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
FUNCTION AlmostEqual(Num1,Num2)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL            :: Num1,Num2      ! Number
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

FUNCTION UNITVECTOR(v1)
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

END MODULE MOD_Particle_Globals
