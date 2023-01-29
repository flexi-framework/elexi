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

! INTERFACE StoreLostParticleProperties
!   MODULE PROCEDURE StoreLostParticleProperties
! END INTERFACE

PUBLIC :: UpdateNextFreePosition
PUBLIC :: DiceUnitVector
! PUBLIC :: StoreLostParticleProperties
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition()
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars        ,ONLY: PDM,PEM,PartSpecies,useLinkedList
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: counter,i,n
INTEGER            :: ElemID
!===================================================================================================================================

IF(PDM%maxParticleNumber.EQ.0) RETURN

IF (useLinkedList) THEN
  PEM%pNumber(:) = 0
END IF

n = PDM%ParticleVecLength !PDM%maxParticleNumber
PDM%ParticleVecLength = 0
counter = 0

IF (useLinkedList) THEN
 DO i = 1,n
   IF (.NOT.PDM%ParticleInside(i)) THEN
     counter = counter + 1
     PDM%nextFreePosition(counter) = i
   ELSE
      ! Sanity check corrupted particle list (some or all entries of a particle become zero, including the species ID)
      IF(PartSpecies(i).LE.0) CALL Abort(__STAMP__,'Species ID is zero for PartID=',IntInfo=i)

      ElemID = PEM%Element(i)
        ! Start of linked list for particles in elem
        IF (PEM%pNumber(ElemID).EQ.0) THEN
          PEM%pStart(ElemID)          = i
        ! Next particle of same elem (linked list)
        ELSE
          PEM%pNext(PEM%pEnd(ElemID)) = i
     END IF
     PEM%pEnd(   ElemID)   = i
     PEM%pNumber(ElemID)   = PEM%pNumber(ElemID) + 1
     PDM%ParticleVecLength = i
   END IF
 END DO

! no linked list
ELSE
 DO i=1,n
   IF (.NOT.PDM%ParticleInside(i)) THEN
     counter = counter + 1
     PDM%nextFreePosition(counter) = i
   ELSE
     PDM%ParticleVecLength = i
   END IF
 END DO
ENDIF

PDM%CurrentNextFreePosition = 0

! Positions after ParticleVecLength in freePosition
DO i = n+1,PDM%maxParticleNumber
  counter = counter + 1
  PDM%nextFreePosition(counter) = i
END DO

! Set nextFreePosition for occupied slots to zero
PDM%nextFreePosition(counter+1:PDM%maxParticleNumber) = 0
! If maxParticleNumber are inside, counter is greater than maxParticleNumber
IF (counter+1.GT.PDM%MaxParticleNumber) PDM%nextFreePosition(PDM%MaxParticleNumber) = 0

END SUBROUTINE UpdateNextFreePosition


!SUBROUTINE StoreLostParticleProperties(iPart,ElemID,UsePartState_opt,PartMissingType_opt)
!!----------------------------------------------------------------------------------------------------------------------------------!
!! Store information of a lost particle (during restart and during the simulation)
!!----------------------------------------------------------------------------------------------------------------------------------!
!! MODULES                                                                                                                          !
!USE MOD_Globals                ,ONLY: Abort,MyRank
!USE MOD_Particle_Vars          ,ONLY: PartSpecies,Species,PartState,LastPartPos
!USE MOD_Particle_Tracking_Vars ,ONLY: PartStateLost,PartLostDataSize,PartStateLostVecLength
!USE MOD_TimeDisc_Vars          ,ONLY: t
!!----------------------------------------------------------------------------------------------------------------------------------!
!! insert modules here
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT / OUTPUT VARIABLES
!INTEGER,INTENT(IN)          :: iPart
!INTEGER,INTENT(IN)          :: ElemID ! Global element index
!LOGICAL,INTENT(IN),OPTIONAL :: UsePartState_opt
!INTEGER,INTENT(IN),OPTIONAL :: PartMissingType_opt ! 0: lost, 1: missing & found once, >1: missing & multiply found
!INTEGER                     :: dims(2)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!! Temporary arrays
!REAL, ALLOCATABLE    :: PartStateLost_tmp(:,:)   ! (1:11,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,MPF,time,ElemID,iPart
!!                                                !                 2nd index: 1 to number of lost particles
!INTEGER              :: ALLOCSTAT
!LOGICAL              :: UsePartState_loc
!!===================================================================================================================================
!UsePartState_loc = .FALSE.
!IF (PRESENT(UsePartState_opt)) UsePartState_loc = UsePartState_opt
!MPF = Species(PartSpecies(iPart))%MacroParticleFactor

!dims = SHAPE(PartStateLost)

!ASSOCIATE( iMax => PartStateLostVecLength )
!  ! Increase maximum number of boundary-impact particles
!  iMax = iMax + 1

!  ! Check if array maximum is reached.
!  ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
!  IF(iMax.GT.dims(2))THEN

!    ! --- PartStateLost ---
!    ALLOCATE(PartStateLost_tmp(1:PartLostDataSize,1:dims(2)), STAT=ALLOCSTAT)
!    IF (ALLOCSTAT.NE.0) CALL Abort(&
!          __STAMP__&
!          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateLost_tmp temporary array!')
!    ! Save old data
!    PartStateLost_tmp(1:PartLostDataSize,1:dims(2)) = PartStateLost(1:PartLostDataSize,1:dims(2))

!    ! Re-allocate PartStateLost to twice the size
!    DEALLOCATE(PartStateLost)
!    ALLOCATE(PartStateLost(1:PartLostDataSize,1:2*dims(2)), STAT=ALLOCSTAT)
!    IF (ALLOCSTAT.NE.0) CALL Abort(&
!          __STAMP__&
!          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateLost array!')
!    PartStateLost(1:PartLostDataSize,        1:  dims(2)) = PartStateLost_tmp(1:PartLostDataSize,1:dims(2))
!    PartStateLost(1:PartLostDataSize,dims(2)+1:2*dims(2)) = 0.

!  END IF

!  ! 1-3: Particle position (last valid position)
!  IF(UsePartState_loc)THEN
!    PartStateLost(1:3,iMax) = PartState(1:3,iPart)
!  ELSE
!    PartStateLost(1:3,iMax) = LastPartPos(1:3,iPart)
!  END IF ! UsePartState_loc
!  ! 4-6: Particle velocity
!  PartStateLost(4:6  ,iMax) = PartState(4:6,iPart)
!  ! 7: SpeciesID
!  PartStateLost(7    ,iMax) = REAL(PartSpecies(iPart))
!  ! 8: Macro particle factor
!  PartStateLost(8    ,iMax) = MPF
!  ! 9: time of loss
!  PartStateLost(9    ,iMax) = time
!  ! 10: Global element ID
!  PartStateLost(10   ,iMax) = REAL(ElemID)
!  ! 11: Particle ID
!  PartStateLost(11   ,iMax) = REAL(iPart)
!  ! 12-14: Particle position (position of loss)
!  PartStateLost(12:14,iMax) = PartState(1:3,iPart)
!  ! 15: myrank
!  PartStateLost(15,iMax) = myrank
!  ! 16: missing type, i.e., 0: lost, 1: missing & found once, >1: missing & multiply found
!  IF(PRESENT(PartMissingType_opt))THEN ! when particles go missing during restart (maybe they are found on other procs or lost)
!    PartStateLost(16,iMax) = PartMissingType_opt
!  ELSE ! simply lost during the simulation
!    PartStateLost(16,iMax) = 0
!  END IF ! PRESENT(PartMissingType_opt)
!END ASSOCIATE

!END SUBROUTINE StoreLostParticleProperties


FUNCTION DiceUnitVector()
!===================================================================================================================================
!> Calculates random unit vector
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
REAL                     :: rRan, cos_scatAngle, sin_scatAngle, rotAngle
!===================================================================================================================================
CALL RANDOM_NUMBER(rRan)

cos_scatAngle     = 2.*rRan-1.
sin_scatAngle     = SQRT(1. - cos_scatAngle ** 2.)
DiceUnitVector(1) = cos_scatAngle

CALL RANDOM_NUMBER(rRan)
rotAngle          = 2. * Pi * rRan

DiceUnitVector(2) = sin_scatAngle * COS(rotAngle)
DiceUnitVector(3) = sin_scatAngle * SIN(rotAngle)

END FUNCTION DiceUnitVector

END MODULE MOD_Part_Tools
