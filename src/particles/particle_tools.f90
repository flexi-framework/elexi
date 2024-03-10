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
! Contains tools for particle related operations. This routine uses MOD_Particle_Boundary_Tools, but not vice versa!
!===================================================================================================================================
MODULE MOD_Particle_Tools
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE UpdateNextFreePosition
  MODULE PROCEDURE UpdateNextFreePosition
END INTERFACE

INTERFACE GetNextFreePosition
  MODULE PROCEDURE GetNextFreePosition
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

INTERFACE IncreaseMaxParticleNumber
  MODULE PROCEDURE IncreaseMaxParticleNumber
END INTERFACE

INTERFACE ReduceMaxParticleNumber
  MODULE PROCEDURE ReduceMaxParticleNumber
END INTERFACE

PUBLIC :: UpdateNextFreePosition
PUBLIC :: GetNextFreePosition
PUBLIC :: DiceUnitVector
! PUBLIC :: StoreLostParticleProperties
PUBLIC :: GetOffsetAndGlobalNumberOfParts
PUBLIC :: IncreaseMaxParticleNumber
PUBLIC :: ReduceMaxParticleNumber
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


FUNCTION GetNextFreePosition(Offset)
!===================================================================================================================================
!> Returns the next free position in the particle vector, if no space is available it increases the maximum particle number
!> ATTENTION: If optional argument is used, the PDM%CurrentNextFreePosition will not be updated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars        ,ONLY: PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,OPTIONAL,INTENT(IN) :: Offset
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                     :: GetNextFreePosition
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: i
!===================================================================================================================================
IF (PRESENT(Offset)) THEN
  IF (PDM%CurrentNextFreePosition+Offset.GT.PDM%MaxParticleNumber) THEN
    CALL IncreaseMaxParticleNumber()

    IF (PDM%CurrentNextFreePosition.GT.PDM%MaxParticleNumber) THEN
      ! This only happens if PDM%CurrentNextFreePosition+Offset is way off (which shouldn't happen)
      IPWRITE(UNIT_stdOut,'(A)') "WARNING: PDM%CurrentNextFreePosition+Offset is way off in particle_tools.f90 GetNextFreePosition(Offset), 1"
      CALL IncreaseMaxParticleNumber(CEILING((PDM%CurrentNextFreePosition+Offset)*(1+PDM%MaxPartNumIncrease)-PDM%MaxParticleNumber))
    END IF
  END IF

  GetNextFreePosition = PDM%nextFreePosition(PDM%CurrentNextFreePosition+Offset)
  ! If next free position is equal 0, determine how much more particles are needed to get a position within the particle vector
  IF (GetNextFreePosition.EQ.0) THEN
    CALL IncreaseMaxParticleNumber()

    GetNextFreePosition = PDM%nextFreePosition(PDM%CurrentNextFreePosition+Offset)
    IF (GetNextFreePosition.EQ.0) THEN
      ! This only happens if PDM%CurrentNextFreePosition+Offset is way off (which shouldn't happen)
      IPWRITE(UNIT_stdOut,'(A)') "WARNING: PDM%CurrentNextFreePosition+Offset is way off in particle_tools.f90 GetNextFreePosition(Offset), 2"
      IF (PDM%nextFreePosition(1).EQ.0) THEN
        i = 0
      ELSE
        i = PDM%CurrentNextFreePosition+Offset
        DO WHILE(PDM%nextFreePosition(i).EQ.0.AND.i.GT.0)
          i = i - 1
        END DO
      END IF
      ! Increase the maxpartnum + margin
      CALL IncreaseMaxParticleNumber(CEILING((PDM%CurrentNextFreePosition+Offset-i)*(1+PDM%MaxPartNumIncrease)+PDM%maxParticleNumber*PDM%MaxPartNumIncrease))
      GetNextFreePosition = PDM%nextFreePosition(PDM%CurrentNextFreePosition+Offset)
    END IF
  END IF
ELSE
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
  IF(PDM%CurrentNextFreePosition.GT.PDM%MaxParticleNumber) THEN
    CALL IncreaseMaxParticleNumber()

    IF (PDM%CurrentNextFreePosition.GT.PDM%MaxParticleNumber) THEN
      ! This only happens if PDM%CurrentNextFreePosition is way off (which shouldn't happen)
      IPWRITE(UNIT_stdOut,'(A)') "WARNING: PDM%CurrentNextFreePosition is way off in particle_tools.f90 GetNextFreePosition(), 1"
      CALL IncreaseMaxParticleNumber(CEILING((PDM%CurrentNextFreePosition)*(1+PDM%MaxPartNumIncrease)-PDM%MaxParticleNumber))
    END IF
  END IF

  GetNextFreePosition = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
  ! If next free position is equal 0, determine how much more particles are needed to get a position within the particle vector
  IF(GetNextFreePosition.EQ.0) THEN
    CALL IncreaseMaxParticleNumber()
    GetNextFreePosition = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
    IF (GetNextFreePosition.EQ.0) THEN
      ! This only happens if PDM%CurrentNextFreePosition is way off (which shouldn't happen)
      IPWRITE(UNIT_stdOut,'(A)') "WARNING: PDM%CurrentNextFreePosition is way off in particle_tools.f90 GetNextFreePosition(), 2"
      IF (PDM%nextFreePosition(1).EQ.0) THEN
        i = 0
      ELSE
        i = PDM%CurrentNextFreePosition
        DO WHILE(PDM%nextFreePosition(i).EQ.0.AND.i.GT.0)
          i = i - 1
        END DO
      END IF
      ! Increase the maxpartnum + margin
      CALL IncreaseMaxParticleNumber(CEILING((PDM%CurrentNextFreePosition-i)*(1+PDM%MaxPartNumIncrease)+PDM%maxParticleNumber*PDM%MaxPartNumIncrease))
      GetNextFreePosition = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
    END IF
  END IF

  IF(PDM%ParticleInside(GetNextFreePosition)) CALL Abort(__STAMP__,'This Particle is already in use',IntInfo=GetNextFreePosition)
  IF(GetNextFreePosition.GT.PDM%ParticleVecLength) PDM%ParticleVecLength = GetNextFreePosition
END IF
IF(GetNextFreePosition.EQ.0) CALL Abort(__STAMP__,'This should not happen, PDM%MaxParticleNumber reached',IntInfo=PDM%MaxParticleNumber)

END FUNCTION GetNextFreePosition


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


SUBROUTINE IncreaseMaxParticleNumber(Amount)
!===================================================================================================================================
! Increases MaxParticleNumber and increases size of all depended arrays
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Array_Operations            ,ONLY: ChangeSizeArray
USE MOD_Particle_Interpolation_Vars ,ONLY: FieldAtParticle
USE MOD_Particle_Vars
#if USE_MPI
USE MOD_Particle_MPI_Vars           ,ONLY: PartTargetProc
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: Amount
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: NewSize,i,ii
!===================================================================================================================================
IF (PRESENT(Amount)) THEN
  IF (Amount.EQ.0) RETURN
  NewSize = PDM%MaxParticleNumber+Amount
  IF (NewSize.GT.PDM%maxAllowedParticleNumber) &
    CALL Abort(__STAMP__,'More Particles needed than allowed in PDM%maxAllowedParticleNumber',IntInfo=NewSize)
ELSE
  NewSize = MAX(CEILING(PDM%MaxParticleNumber*(1+PDM%MaxPartNumIncrease)),PDM%MaxParticleNumber+1)
  IF (PDM%MaxParticleNumber.GE.PDM%maxAllowedParticleNumber) &
    CALL Abort(__STAMP__,'More Particles needed than allowed in PDM%maxAllowedParticleNumber',IntInfo=NewSize)
  NewSize = MIN(NewSize,PDM%maxAllowedParticleNumber)
END IF

IF(ALLOCATED(PEM%Element))        CALL ChangeSizeArray(PEM%Element       ,PDM%maxParticleNumber,NewSize)
IF(ALLOCATED(PEM%lastElement))    CALL ChangeSizeArray(PEM%lastElement   ,PDM%maxParticleNumber,NewSize)
! IF(ALLOCATED(PEM%pNext))          CALL ChangeSizeArray(PEM%pNext         ,PDM%maxParticleNumber,NewSize)

IF(ALLOCATED(PDM%ParticleInside)) CALL ChangeSizeArray(PDM%ParticleInside,PDM%maxParticleNumber,NewSize,.FALSE.)
IF(ALLOCATED(PDM%IsNewPart))      CALL ChangeSizeArray(PDM%IsNewPart     ,PDM%maxParticleNumber,NewSize,.FALSE.)

IF(ALLOCATED(PartState))          CALL ChangeSizeArray(PartState         ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(TurbPartState))      CALL ChangeSizeArray(TurbPartState     ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(LastPartPos))        CALL ChangeSizeArray(LastPartPos       ,PDM%maxParticleNumber,NewSize)
IF(ALLOCATED(PartPosRef))         CALL ChangeSizeArray(PartPosRef        ,PDM%maxParticleNumber,NewSize,-888.)
IF(ALLOCATED(PartReflCount))      CALL ChangeSizeArray(PartReflCount     ,PDM%maxParticleNumber,NewSize,0)
IF(ALLOCATED(PartSpecies))        CALL ChangeSizeArray(PartSpecies       ,PDM%maxParticleNumber,NewSize,0)
IF(ALLOCATED(PartIndex))          CALL ChangeSizeArray(PartIndex         ,PDM%maxParticleNumber,NewSize,0)

IF(ALLOCATED(Pt))                 CALL ChangeSizeArray(Pt                ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(Pt_temp))            CALL ChangeSizeArray(Pt_temp           ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(FieldAtParticle))    CALL ChangeSizeArray(FieldAtParticle   ,PDM%maxParticleNumber,NewSize)

#if USE_BASSETFORCE
IF(ALLOCATED(durdt))              CALL ChangeSizeArray(durdt             ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(bIter))              CALL ChangeSizeArray(bIter             ,PDM%maxParticleNumber,NewSize,0)
IF(ALLOCATED(Fbdt))               CALL ChangeSizeArray(Fbdt              ,PDM%maxParticleNumber,NewSize,0.)
#endif /* USE_BASSETFORCE */

#if USE_EXTEND_RHS && ANALYZE_RHS
IF(ALLOCATED(Pt_ext))             CALL ChangeSizeArray(Pt_ext            ,PDM%maxParticleNumber,NewSize,0.)
#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

#if USE_MPI
IF(ALLOCATED(PartTargetProc))     CALL ChangeSizeArray(PartTargetProc    ,PDM%maxParticleNumber,NewSize)
#endif

IF(ALLOCATED(PDM%nextFreePosition)) THEN
  CALL ChangeSizeArray(PDM%nextFreePosition,PDM%maxParticleNumber,NewSize,0)

  ! Search for first entry where new poition is available
  i = 1
  DO WHILE(PDM%nextFreePosition(i).NE.0)
    i = i+1
  END DO
  i = i-1
  ! Fill the free spots with the new entrys
  DO ii = 1,NewSize-PDM%MaxParticleNumber
    PDM%nextFreePosition(i+ii) = ii+PDM%MaxParticleNumber
  END DO
END IF

PDM%MaxParticleNumber = NewSize

END SUBROUTINE IncreaseMaxParticleNumber


SUBROUTINE ReduceMaxParticleNumber()
!===================================================================================================================================
! Reduces MaxParticleNumber and increases size of all depended arrays
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Array_Operations            ,ONLY: ChangeSizeArray
USE MOD_Particle_Interpolation_Vars ,ONLY: FieldAtParticle
USE MOD_Particle_Vars
#if USE_MPI
USE MOD_Particle_MPI_Vars           ,ONLY: PartTargetProc
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES

!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: NewSize,i,ii,nPart
!===================================================================================================================================

nPart = 0
DO i = 1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) nPart = nPart + 1
END DO

! Reduce Arrays only for at least PDM%maxParticleNumber*PDM%MaxPartNumIncrease free spots
IF (nPart.GE.PDM%maxParticleNumber/(1.+PDM%MaxPartNumIncrease)**2) RETURN

! Maintain nPart*PDM%MaxPartNumIncrease free spots
Newsize = MAX(CEILING(nPart*(1.+PDM%MaxPartNumIncrease)),1)
IF (Newsize.EQ.PDM%maxParticleNumber) RETURN

IF (.NOT.PDM%RearrangePartIDs) THEN
  ! Search for highest occupied particle index and set Newsize to this Value
  i = PDM%maxParticleNumber
  DO WHILE(.NOT.PDM%ParticleInside(i).OR.i.EQ.NewSize)
    i = i-1
  END DO
  NewSize = i
ELSE
  ! Rearrange particles with IDs>NewSize to lower IDs
  DO i = NewSize+1,PDM%maxParticleNumber
    IF (PDM%ParticleInside(i)) THEN
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
      ii = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
      IF(ii.EQ.0.OR.ii.GT.NewSize) THEN
        CALL UpdateNextFreePosition()
        PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
        ii = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
        IF (ii.EQ.0.OR.ii.GT.NewSize) CALL Abort(__STAMP__,'Error in particle re-arrange!')
      END IF
      IF (PDM%ParticleVecLength.LT.ii) PDM%ParticleVecLength = ii
      CALL ChangePartID(i,ii)
    END IF
  END DO
END IF

IF(ALLOCATED(PEM%Element))        CALL ChangeSizeArray(PEM%Element       ,PDM%maxParticleNumber,NewSize)
IF(ALLOCATED(PEM%lastElement))    CALL ChangeSizeArray(PEM%lastElement   ,PDM%maxParticleNumber,NewSize)
! IF(ALLOCATED(PEM%pNext))          CALL ChangeSizeArray(PEM%pNext         ,PDM%maxParticleNumber,NewSize)

IF(ALLOCATED(PDM%ParticleInside)) CALL ChangeSizeArray(PDM%ParticleInside,PDM%maxParticleNumber,NewSize,.FALSE.)
IF(ALLOCATED(PDM%IsNewPart))      CALL ChangeSizeArray(PDM%IsNewPart     ,PDM%maxParticleNumber,NewSize,.FALSE.)

IF(ALLOCATED(PartState))          CALL ChangeSizeArray(PartState         ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(TurbPartState))      CALL ChangeSizeArray(TurbPartState     ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(LastPartPos))        CALL ChangeSizeArray(LastPartPos       ,PDM%maxParticleNumber,NewSize)
IF(ALLOCATED(PartPosRef))         CALL ChangeSizeArray(PartPosRef        ,PDM%maxParticleNumber,NewSize,-888.)
IF(ALLOCATED(PartReflCount))      CALL ChangeSizeArray(PartReflCount     ,PDM%maxParticleNumber,NewSize,0)
IF(ALLOCATED(PartSpecies))        CALL ChangeSizeArray(PartSpecies       ,PDM%maxParticleNumber,NewSize,0)
IF(ALLOCATED(PartIndex))          CALL ChangeSizeArray(PartIndex         ,PDM%maxParticleNumber,NewSize,0)

IF(ALLOCATED(Pt))                 CALL ChangeSizeArray(Pt                ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(Pt_temp))            CALL ChangeSizeArray(Pt_temp           ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(FieldAtParticle))    CALL ChangeSizeArray(FieldAtParticle   ,PDM%maxParticleNumber,NewSize)

#if USE_BASSETFORCE
IF(ALLOCATED(durdt))              CALL ChangeSizeArray(durdt             ,PDM%maxParticleNumber,NewSize,0.)
IF(ALLOCATED(bIter))              CALL ChangeSizeArray(bIter             ,PDM%maxParticleNumber,NewSize,0)
IF(ALLOCATED(Fbdt))               CALL ChangeSizeArray(Fbdt              ,PDM%maxParticleNumber,NewSize,0.)
#endif /* USE_BASSETFORCE */

#if USE_EXTEND_RHS && ANALYZE_RHS
IF(ALLOCATED(Pt_ext))             CALL ChangeSizeArray(Pt_ext            ,PDM%maxParticleNumber,NewSize,0.)
#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

#if USE_MPI
IF(ALLOCATED(PartTargetProc))     CALL ChangeSizeArray(PartTargetProc    ,PDM%maxParticleNumber,NewSize)
#endif

IF(ALLOCATED(PDM%nextFreePosition)) THEN
  CALL ChangeSizeArray(PDM%nextFreePosition,PDM%maxParticleNumber,NewSize,0)

  ! Set all NextFreePositions to zero which points to a partID>NewSize
  DO i = 1,NewSize
    IF(PDM%nextFreePosition(i).GT.NewSize) PDM%nextFreePosition(i) = 0
  END DO
END IF

IF (PDM%ParticleVecLength.GT.NewSize) PDM%ParticleVecLength = NewSize
PDM%MaxParticleNumber = NewSize

CALL UpdateNextFreePosition()

END SUBROUTINE ReduceMaxParticleNumber


SUBROUTINE ChangePartID(OldID,NewID)
!===================================================================================================================================
! Change PartID from OldID to NewID
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars                   ,ONLY: nElems,offsetElem
USE MOD_Particle_Interpolation_Vars ,ONLY: FieldAtParticle
USE MOD_Particle_Vars
#if USE_MPI
USE MOD_Particle_MPI_Vars           ,ONLY: PartTargetProc
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: OldID
INTEGER,INTENT(IN)        :: NewID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i,TempPartID,LocElemID
!===================================================================================================================================

IF(ALLOCATED(PEM%Element))     PEM%Element(NewID)     = PEM%Element(OldID)
IF(ALLOCATED(PEM%lastElement)) PEM%lastElement(NewID) = PEM%lastElement(OldID)
IF(ALLOCATED(PEM%pNext)) THEN
  PEM%pNext(NewID) = PEM%pNext(OldID)
  ! Update pNext onto this particle
  LocElemID  = PEM%Element(OldID) - offsetElem
  ! Ignore mapping for particle to be sent
  IF (LocElemID.GE.1 .AND. LocElemID.LE.nElems) THEN
    TempPartID = PEM%pStart(LocElemID)
    IF (TempPartID.EQ.OldID) THEN
      PEM%pStart(LocElemID) = NewID
    ELSE
      DO i=1,PEM%pNumber(LocElemID)
        IF(PEM%pNext(TempPartID).EQ.OldID) THEN
          PEM%pNext(TempPartID) = NewID
          EXIT
        END IF
        TempPartID = PEM%pNext(TempPartID)
      END DO
    END IF
  END IF
END IF

IF(ALLOCATED(PDM%ParticleInside)) THEN
  PDM%ParticleInside(NewID) = PDM%ParticleInside(OldID)
  PDM%ParticleInside(OldID) = .FALSE.
END IF
IF(ALLOCATED(PDM%IsNewPart)) THEN
  PDM%IsNewPart(NewID)      = PDM%IsNewPart(OldID)
  PDM%IsNewPart(OldID)      = .FALSE.
END IF

IF(ALLOCATED(PartState)) THEN
  PartState(:,NewID)        = PartState(:,OldID)
  PartState(:,OldID)        = 0.
END IF
IF(ALLOCATED(TurbPartState)) THEN
  TurbPartState(:,NewID)    = TurbPartState(:,OldID)
  TurbPartState(:,OldID)    = 0.
END IF
IF(ALLOCATED(LastPartPos)) THEN
  LastPartPos(:,NewID)      = LastPartPos(:,OldID)
END IF
IF(ALLOCATED(PartPosRef)) THEN
  PartPosRef(:,NewID)       = PartPosRef(:,OldID)
  PartPosRef(:,OldID)       = -888.
END IF
IF(ALLOCATED(PartReflCount)) THEN
  PartReflCount(NewID)      = PartReflCount(OldID)
  PartReflCount(OldID)      = 0
END IF
IF(ALLOCATED(PartSpecies)) THEN
  PartSpecies(NewID)        = PartSpecies(OldID)
  PartSpecies(OldID)        = 0
END IF
IF(ALLOCATED(PartIndex)) THEN
  PartIndex(NewID)          = PartIndex(OldID)
  PartIndex(OldID)          = 0
END IF

IF(ALLOCATED(Pt)) THEN
  Pt(:,NewID)               = Pt(:,OldID)
  Pt(:,OldID)               = 0.
END IF
IF(ALLOCATED(Pt_temp)) THEN
  Pt_temp(:,NewID)          = Pt_temp(:,OldID)
  Pt_temp(:,OldID)          = 0.
END IF
IF(ALLOCATED(FieldAtParticle)) THEN
  FieldAtParticle(:,NewID)  = FieldAtParticle(:,OldID)
END IF

#if USE_BASSETFORCE
IF(ALLOCATED(durdt)) THEN
  durdt(:,NewID)            = durdt(:,OldID)
  durdt(:,OldID)            = 0.
END IF
IF(ALLOCATED(bIter)) THEN
  bIter(NewID)              = bIter(OldID)
  bIter(OldID)              = 0
END IF
IF(ALLOCATED(Fbdt)) THEN
  Fbdt(:,NewID)             = Fbdt(:,OldID)
  Fbdt(:,OldID)             = 0.
END IF
#endif /* USE_BASSETFORCE */

#if USE_EXTEND_RHS && ANALYZE_RHS
IF(ALLOCATED(Pt_ext)) THEN
  Pt_ext(:,NewID)           = Pt_ext(:,OldID)
  Pt_ext(:,OldID)           = 0.
END IF
#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

#if USE_MPI
IF(ALLOCATED(PartTargetProc)) PartTargetProc(NewID) = PartTargetProc(OldID)
#endif

END SUBROUTINE ChangePartID


END MODULE MOD_Particle_Tools
