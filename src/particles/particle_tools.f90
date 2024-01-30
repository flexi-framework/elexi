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

!===================================================================================================================================
! Contains tools for particles
!===================================================================================================================================
MODULE MOD_Particle_Tools
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

INTERFACE GetOffsetAndGlobalNumberOfParts
  MODULE PROCEDURE GetOffsetAndGlobalNumberOfParts
END INTERFACE

PUBLIC :: UpdateNextFreePosition
PUBLIC :: DiceUnitVector
! PUBLIC :: StoreLostParticleProperties
PUBLIC :: GetOffsetAndGlobalNumberOfParts
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


FUNCTION DiceUnitVector()
!===================================================================================================================================
!> Calculates random unit vector
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals_Vars,      ONLY: PI
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


!===================================================================================================================================
!> Calculate the particle offset and global number of particles across all processors
!> In this routine the number are calculated using integer KIND=8, but are returned with integer KIND=ICC in order to test if using
!> integer KIND=8 is required for total number of particles, particle boundary state, lost particles or clones
!===================================================================================================================================
SUBROUTINE GetOffsetAndGlobalNumberOfParts(CallingRoutine,offsetnPart,globnPart,locnPart,GetMinMaxNbrOfParticles)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: CallingRoutine
INTEGER(KIND=IK),INTENT(IN)  :: locnPart
LOGICAL,INTENT(IN)           :: GetMinMaxNbrOfParticles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=IK),INTENT(OUT) :: offsetnPart
INTEGER(KIND=IK),INTENT(OUT) :: globnPart(6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER(KIND=8)              :: globnPart8                         ! always integer KIND=8
INTEGER(KIND=8)              :: locnPart8,locnPart8Recv            ! always integer KIND=8
INTEGER(KIND=IK)             :: SimNumSpecMin,SimNumSpecMax
#endif
!===================================================================================================================================
#if USE_MPI
locnPart8     = INT(locnPart,8)
locnPart8Recv = 0_IK
CALL MPI_EXSCAN(locnPart8,locnPart8Recv,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_FLEXI,iError)
offsetnPart   = INT(locnPart8Recv,KIND=IK)

! Last proc calculates the global number and broadcasts it
IF(myRank.EQ.nProcessors-1) locnPart8=locnPart8Recv+locnPart8
CALL MPI_BCAST(locnPart8,1,MPI_INTEGER8,nProcessors-1,MPI_COMM_FLEXI,iError)

! Global numbers
globnPart8   = locnPart8
LOGWRITE(*,*) TRIM(CallingRoutine)//'offsetnPart,locnPart,globnPart8',offsetnPart,locnPart,globnPart8

! Sanity check: Add up all particles with integer KIND=8 and compare
IF (MPIRoot) THEN
  ! Check if offsetnPart is kind=8 is the number of particles is larger than integer KIND=4
  IF (globnPart8.GT.INT(HUGE(offsetnPart),8)) THEN
    WRITE(UNIT_stdOut,'(A,I0)') '\n\n\nTotal number of particles  : ',globnPart8
    WRITE(UNIT_stdOut,'(A,I0)')       'Maximum number of particles: ',HUGE(offsetnPart)
    CALL Abort(__STAMP__,TRIM(CallingRoutine)//' has encountered more than integer KIND=4 particles!')
  END IF
END IF ! MPIRoot

! Get min/max number of particles
SimNumSpecMin = 0
SimNumSpecMax = 0
IF(GetMinMaxNbrOfParticles)THEN
  IF (MPIRoot) THEN
    CALL MPI_REDUCE(locnPart,SimNumSpecMin,1,MPI_INTEGER_INT_KIND,MPI_MIN,0,MPI_COMM_FLEXI,IERROR)
    CALL MPI_REDUCE(locnPart,SimNumSpecMax,1,MPI_INTEGER_INT_KIND,MPI_MAX,0,MPI_COMM_FLEXI,IERROR)
  ELSE
    CALL MPI_REDUCE(locnPart,0            ,1,MPI_INTEGER_INT_KIND,MPI_MIN,0,MPI_COMM_FLEXI,IERROR)
    CALL MPI_REDUCE(locnPart,0            ,1,MPI_INTEGER_INT_KIND,MPI_MAX,0,MPI_COMM_FLEXI,IERROR)
  END IF
END IF ! GetMinMaxNbrOfParticles

! Cast to Kind=IK before returning the number
globnPart(1) = INT(SimNumSpecMin , KIND = IK)
globnPart(2) = INT(SimNumSpecMax , KIND = IK)
globnPart(3) = INT(globnPart8    , KIND = IK)
#else
offsetnPart=0_IK
globnPart(1:3)=INT(locnPart,KIND=IK)

! Suppress compiler warning
NO_OP(CallingRoutine)
#endif

! Get extrema over the complete simulation only during WriteParticleToHDF5
IF(GetMinMaxNbrOfParticles)THEN
  globnPart(4) = MIN(globnPart(1),globnPart(4))
  globnPart(5) = MAX(globnPart(2),globnPart(5))
  globnPart(6) = MAX(globnPart(3),globnPart(6))
END IF ! GetMinMaxNbrOfParticles

END SUBROUTINE GetOffsetAndGlobalNumberOfParts

END MODULE MOD_Particle_Tools
