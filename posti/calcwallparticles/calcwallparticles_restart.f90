!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz 
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
!> \brief Routines that handle restart capabilities.
!> 
!> With this feature a simulation can be resumed from a state file that has been created during a previous
!> simulation (restart file). The restart file is passed to FLEXI as a second command line argument.
!> The restart can also be performed from a file with a different polynomial degree or node type than the current simulation.
!==================================================================================================================================
MODULE MOD_CalcWallParticles_Restart
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CalcWallParticlesRestart
  MODULE PROCEDURE CalcWallParticlesRestart
END INTERFACE

PUBLIC :: CalcWallParticlesRestart
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE CalcWallParticlesRestart(boundary_opt,reflCount_opt)
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Restart_Vars   
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_Restart_Vars,            ONLY:RestartTime,RestartFile
USE MOD_Particle_Vars,           ONLY:PartState, PartSpecies, PEM, PDM, PartReflCount,LastPartPos
USE MOD_Part_Tools,              ONLY:UpdateNextFreePosition
USE MOD_Eval_XYZ,                ONLY:TensorProductInterpolation
USE MOD_Erosionpoints_Vars,      ONLY:EP_Impacts,EPDataSize
USE MOD_Particle_Mesh,           ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Tracking,       ONLY:ParticleTracing,ParticleRefTracking
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
USE MOD_Posti_CalcWallParticles_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL      :: boundary_opt
INTEGER,INTENT(IN),OPTIONAL      :: reflCount_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: COUNTER, EP_restart, EP_considered
INTEGER                  :: i, EP_glob, ErosionDim
REAL,ALLOCATABLE         :: PartData(:,:)
LOGICAL                  :: ErosionDataExists
!===================================================================================================================================

EP_Impacts = 0
COUNTER    = 0

! Open the restart file and search for erosionData
CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL DatasetExists(File_ID,'ErosionData',ErosionDataExists)

CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)

IF(ErosionDataExists) THEN
    CALL GetDataSize(File_ID,'ErosionData',ErosionDim,HSize)
!    CHECKSAFEINT(HSize(1),4)
    EP_glob    = HSize(1)
!    SWRITE(UNIT_stdOut,'(A3,A30,A3,I33)')' | ','Number of impacts',' | ',EP_glob
    
    ! We lost the impact <-> proc association, so fill the entire array
    ALLOCATE(PartData(1:EP_glob,1:EPDataSize))
    CALL ReadArray(ArrayName='ErosionData', rank=2,&
                 nVal=      (/EP_glob  ,EPDataSize/),&
                 offset_in  = 0,&
                 offset_dim = 1,&
                 RealArray  = PartData(1:EP_glob,1:EPDataSize))
    ! Pretend all impacts happened on MPI_ROOT, so we can write out
    EP_Impacts = EP_glob
    WRITE(UNIT_stdOut,'(A3,A30,A3,I33)')' | ','Found impacts',' | ',EP_glob
    
    IF (PDM%maxParticleNumber.LE.EP_Impacts) THEN
        CALL Abort(&
          __STAMP__,&
          'Insufficient array size. Please raise Part-maxParticleNumber!')
    END IF
    
    DO i = 1,EP_Impacts
        PartState(i,1:3) = PartData(i,1:3)                    ! Position
        PartState(i,4:6) = PartData(i,4:6)                    ! Velocity
        PartSpecies(i)   = PartData(i,7)                      ! Species
        PartReflCount(i) = PartData(i,10)                     ! ReflectionCounter
        ! We have to move the particle minimally away from the wall for the code to find them. Use opposite velocity vector and user
        ! time step to do so (only half so we safely collide in the next step)
        PartState(i,1:3) = PartState(i,1:3) - dt_remap/2. * PartState(i,4:6)
        ! Assume the particles are still inside
        IF (PRESENT(boundary_opt)) THEN
            IF (PartData(i,8).EQ.boundary_opt) THEN
                PDM%ParticleInside(i) = .TRUE.
            ELSE
                PDM%ParticleInside(i) = .FALSE.
            END IF
        ELSE
                PDM%ParticleInside(i) = .TRUE.
        END IF
        ! Check if particles have the right reflection counter
        IF (PRESENT(reflCount_opt)) THEN
            IF (PartData(i,10).EQ.reflCount_opt) THEN
                PDM%ParticleInside(i) = .TRUE.
            ELSE
                PDM%ParticleInside(i) = .FALSE.
            END IF
        ELSE
                PDM%ParticleInside(i) = .TRUE.
        END IF
    END DO
    
    ! Get rid of the restart arrays and get the particle arrays ready
    SDEALLOCATE(PartData)
    PDM%ParticleVecLength = EP_Impacts
    EP_restart            = EP_Impacts
    EP_Impacts            = 0
    EP_considered         = 0
    CALL UpdateNextFreePosition()
        
    ! Now try to locate the element
    DO i = 1,PDM%ParticleVecLength
        IF (PDM%ParticleInside(i)) THEN
            IF (DoRefMapping) THEN
                CALL SingleParticleToExactElement(i,doHALO=.FALSE.,initFix=.FALSE.,doRelocate=.FALSE.)
            ELSE
                CALL SingleParticleToExactElementNoMap(i,doHALO=.FALSE.,doRelocate=.FALSE.)
            END IF
            ! Check if we lost any particles
            IF (.NOT.PDM%ParticleInside(i)) THEN
                COUNTER = COUNTER + 1
            ELSE
                PEM%LastElement(i) = PEM%Element(i)
                EP_considered      = EP_considered + 1
            END IF
        END IF
    END DO
    
    IF (COUNTER.EQ.0) THEN
        WRITE(UNIT_stdOut,'(A3,A30,A3,I33,A)')' | ','Found matching elements for',' | ', EP_considered+COUNTER,' impacts.'
    ELSE
        WRITE(UNIT_stdOut,'(A,I5,A,I5,A)')' Lost ', COUNTER,' out of ', EP_considered+COUNTER, ' expected impact locations.'
    END IF
    WRITE(UNIT_StdOut,'(132("-"))')
    
    CALL UpdateNextFreePosition()
    
    ! All particles positioned correctly, so save their position as LastPartPos
    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
    
    ! Now collide the particles again with the wall
    DO i = 1,PDM%ParticleVecLength
        PartState(i,1) = PartState(i,1) + PartState(i,4)*10*dt_remap
        PartState(i,2) = PartState(i,2) + PartState(i,5)*10*dt_remap
        PartState(i,3) = PartState(i,3) + PartState(i,6)*10*dt_remap
    END DO
    
    IF (DoRefMapping) THEN
        CALL ParticleRefTracking()
    ELSE
        CALL ParticleTracing()
    END IF
    
    IF (EP_considered.NE.EP_Impacts) THEN
        WRITE(UNIT_stdOut,'(A,I5,A,I5,A)')' Warning: ', EP_considered-EP_Impacts,' out of ', EP_considered, &
                                          ' impacts could not be recreated.'
    END IF
ELSE
    WRITE(UNIT_stdOut,'(A3,A30,A3)')' | ',' No PartData in restart file.',' | '
    WRITE(UNIT_StdOut,'(132("-"))')
END IF ! PartDataExists

CALL CloseDataFile()


END SUBROUTINE CalcWallParticlesRestart

END MODULE MOD_CalcWallParticles_Restart
