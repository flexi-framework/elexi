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

!===================================================================================================================================
! Contains utils required by xy- modules
!===================================================================================================================================
MODULE MOD_Particle_Utils
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

! INTERFACE BubbleSortID
!   MODULE PROCEDURE BubbleSortID
! END INTERFACE BubbleSortID

INTERFACE InsertionSort
  MODULE PROCEDURE InsertionSortInt
  MODULE PROCEDURE InsertionSortReal
END INTERFACE InsertionSort

! PUBLIC :: BubbleSortID
PUBLIC :: InsertionSort
!==================================================================================================================================

CONTAINS

! SUBROUTINE BubbleSortID(a,id,len)
! !===================================================================================================================================
! ! bubble sort, taken from rosetta-wiki and modified for own use
! !===================================================================================================================================
! ! MODULES
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! INTEGER,INTENT(IN)                :: len
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! REAL,INTENT(INOUT)                :: a(len)
! INTEGER,INTENT(INOUT),OPTIONAL    :: id(len)
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! REAL                              :: temp
! INTEGER                           :: iloop,jloop, temp2
! LOGICAL                           :: swapped = .TRUE.
! !===================================================================================================================================

! IF(PRESENT(id))THEN
!   DO jloop=len-1,1,-1
!     swapped = .FALSE.
!     DO iloop=1,jloop
!       IF (a(iloop).GT.a(iloop+1))THEN
!         ! switch entries
!         temp=a(iloop)
!         a(iloop) = a(iloop+1)
!         a(iloop+1) = temp
!         ! switch ids
!         temp2=id(iloop)
!         id(iloop) = id(iloop+1)
!         id(iloop+1) = temp2
!         swapped = .TRUE.
!       END IF
!     END DO ! iloop
!     IF (.NOT. swapped) EXIT
!   END DO ! jloop
! ELSE
!   DO jloop=len-1,1,-1
!     swapped = .FALSE.
!     DO iloop=1,jloop
!       IF (a(iloop).GT.a(iloop+1))THEN
!         ! switch entries
!         temp=a(iloop)
!         a(iloop) = a(iloop+1)
!         a(iloop+1) = temp
!         swapped = .TRUE.
!       END IF
!     END DO ! iloop
!     IF (.NOT. swapped) EXIT
!   END DO ! jloop
! END IF
! END SUBROUTINE BubbleSortID


SUBROUTINE InsertionSortInt(a,id,len)
!===================================================================================================================================
! Insertion sort
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: len
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)             :: a(1:len)
INTEGER,INTENT(INOUT),OPTIONAL    :: id(1:len)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: tmpR
INTEGER                           :: i, j, tmpI
!===================================================================================================================================

IF(PRESENT(ID))THEN
  DO i=2,len
    j=i-1
    tmpR=a(i)
    tmpI=ID(i)
    DO WHILE (j.GE.1) !(j.GE.1 .AND. a(j).GT.tmpR)
      IF (a(j).LE.tmpR) EXIT
      a (j+1) = a(j)
      ID(j+1) = ID(j)
      j=j-1
    END DO
    a (j+1) =tmpR
    ID(j+1) =tmpI
  END DO ! i
ELSE
  DO i=2,len
    j=i-1
    tmpR=a(i)
    DO WHILE (j.GE.1) !(j.GE.1 .AND. a(j).GT.tmpR)
      IF (a(j).LE.tmpR) EXIT
      a (j+1) = a(j)
      j=j-1
    END DO
    a (j+1) =tmpR
  END DO ! i
END IF

END SUBROUTINE InsertionSortInt


SUBROUTINE InsertionSortReal(a,id,len)
!===================================================================================================================================
! Insertion sort
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: len
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                :: a(1:len)
INTEGER,INTENT(INOUT),OPTIONAL    :: id(1:len)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: tmpR
INTEGER                           :: i, j, tmpI
!===================================================================================================================================

IF(PRESENT(ID))THEN
  DO i=2,len
    j=i-1
    tmpR=a(i)
    tmpI=ID(i)
    DO WHILE (j.GE.1) !(j.GE.1 .AND. a(j).GT.tmpR)
      IF (a(j).LE.tmpR) EXIT
      a (j+1) = a(j)
      ID(j+1) = ID(j)
      j=j-1
    END DO
    a (j+1) =tmpR
    ID(j+1) =tmpI
  END DO ! i
ELSE
  DO i=2,len
    j=i-1
    tmpR=a(i)
    DO WHILE (j.GE.1) !(j.GE.1 .AND. a(j).GT.tmpR)
      IF (a(j).LE.tmpR) EXIT
      a (j+1) = a(j)
      j=j-1
    END DO
    a (j+1) =tmpR
  END DO ! i
END IF

END SUBROUTINE InsertionSortReal

END MODULE MOD_Particle_Utils
