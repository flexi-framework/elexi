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

!===================================================================================================================================
! Contains tools for particles
!===================================================================================================================================
MODULE MOD_Part_Tools
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE UpdateNextFreePosition
  MODULE PROCEDURE UpdateNextFreePosition
END INTERFACE

INTERFACE DiceUnitVector
  MODULE PROCEDURE DiceUnitVector
END INTERFACE

PUBLIC :: UpdateNextFreePosition
PUBLIC :: DiceUnitVector
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition()
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars            ,ONLY: offsetElem
USE MOD_Particle_Vars        ,ONLY: PDM,PEM,useLinkedList
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: counter1,i,n
INTEGER            :: ElemID
!===================================================================================================================================

IF(PDM%maxParticleNumber.EQ.0) RETURN

counter1 = 1

IF (useLinkedList) THEN
  PEM%pNumber(:) = 0
END IF

n = PDM%ParticleVecLength !PDM%maxParticleNumber
PDM%ParticleVecLength    = 0
PDM%insideParticleNumber = 0

IF (useLinkedList) THEN
 DO i=1,n
   IF (.NOT.PDM%ParticleInside(i)) THEN
     PDM%nextFreePosition(counter1) = i
     counter1 = counter1 + 1
   ELSE
     ElemID = PEM%Element(i)
      IF (PEM%pNumber(ElemID).EQ.0) THEN
        PEM%pStart(ElemID) = i                          ! Start of Linked List for Particles in Elem
     ELSE
        PEM%pNext(PEM%pEnd(ElemID)) = i                 ! Next Particle of same Elem (Linked List)
     END IF
     PEM%pEnd(ElemID)      = i
     PEM%pNumber(ElemID)   = PEM%pNumber(ElemID) + 1    ! Number of Particles in Element
     PDM%ParticleVecLength = i
   END IF
 END DO

! no linked list
ELSE
 DO i=1,n
   IF (.NOT.PDM%ParticleInside(i)) THEN
     PDM%nextFreePosition(counter1) = i
     counter1 = counter1 + 1
   ELSE
     PDM%ParticleVecLength = i
   END IF
 END DO
ENDIF

PDM%insideParticleNumber    = PDM%ParticleVecLength - counter1+1
PDM%CurrentNextFreePosition = 0

DO i = n+1,PDM%maxParticleNumber
 PDM%nextFreePosition(counter1) = i
 counter1 = counter1 + 1
END DO

PDM%nextFreePosition(counter1:PDM%MaxParticleNumber) = 0 ! exists if MaxParticleNumber is reached!!!
IF (counter1.GT.PDM%MaxParticleNumber) PDM%nextFreePosition(PDM%MaxParticleNumber)=0

RETURN
END SUBROUTINE UpdateNextFreePosition

FUNCTION DiceUnitVector()
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Particle_Globals ,ONLY: Pi
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: DiceUnitVector(3)
REAL                     :: iRan, bVec, aVec
!===================================================================================================================================
CALL RANDOM_NUMBER(iRan)
bVec              = 1. - 2.*iRan
aVec              = SQRT(1. - bVec**2.)
DiceUnitVector(3) = bVec
CALL RANDOM_NUMBER(iRan)
bVec              = Pi *2. * iRan
DiceUnitVector(1) = aVec * COS(bVec)
DiceUnitVector(2) = aVec * SIN(bVec)

END FUNCTION DiceUnitVector

END MODULE MOD_Part_Tools
