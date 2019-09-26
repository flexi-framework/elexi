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

!===================================================================================================================================
! Module for DSMC Sampling and Output
!===================================================================================================================================
MODULE MOD_Particle_Erosion_Analyze
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE CalcSurfaceValues
  MODULE PROCEDURE CalcSurfaceValues
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: CalcSurfaceValues
!===================================================================================================================================

CONTAINS


SUBROUTINE CalcSurfaceValues(restart_opt,remap_opt)
!===================================================================================================================================
!> Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Restart_Vars               ,ONLY: DoRestart,RestartTime
USE MOD_Analyze_Vars               ,ONLY: Analyze_dt
USE MOD_Mesh_Vars                  ,ONLY: MeshFile
USE MOD_Timedisc_Vars              ,ONLY: t
USE MOD_Restart_Vars               ,ONLY: RestartTime
USE MOD_DSMC_Vars                  ,ONLY: MacroSurfaceVal ,MacroSurfaceSpecVal
USE MOD_Particle_Analyze_Vars      ,ONLY: TimeSample
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfMesh,nSurfSample,SampWall
USE MOD_Particle_Boundary_Sampling ,ONLY: WriteSurfSampleToHDF5
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Vars              ,ONLY: nSpecies,WriteMacroSurfaceValues,MacroValSampTime
USE MOD_CalcWallParticles_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL      :: restart_opt   !routine was called during tstep (i.e. before iter=iter+1, t=t+dt...)
CHARACTER(LEN=*),INTENT(IN),OPTIONAL::remap_opt     !routine was called from posti. Change output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iSpec,iSurfSide,p,q,nShift,nShiftRHS
REAL                               :: ActualTime
!===================================================================================================================================

! Set ActualTime to current time
ActualTime = t

IF (WriteMacroSurfaceValues) THEN
    TimeSample       = t - MacroValSampTime
    MacroValSampTime = t
END IF

! Update values if we are called from a restart
IF (PRESENT(restart_opt)) THEN
    IF (restart_opt) THEN
        TimeSample = Analyze_dt
        t          = MERGE(RestartTime,0.,DoRestart)
        ActualTime = t
    END IF
! Avoid division by zero if no simulation time has passed
ELSE
    IF(ALMOSTZERO(TimeSample)) RETURN
END IF

! Do not try to record impacts if there are no walls on the current proc
IF(.NOT.SurfMesh%SurfOnProc) RETURN

! Allocate N+1 Species to have space for average
IF (nSpecies.EQ.1) THEN
    ALLOCATE(MacroSurfaceVal(nErosionVars-1,1:nSurfSample,1:nSurfSample,SurfMesh%nSides))
ELSE
    ALLOCATE(MacroSurfaceVal((nErosionVars-1)*(nSpecies+1),1:nSurfSample,1:nSurfSample,SurfMesh%nSides))
END IF
ALLOCATE(MacroSurfaceSpecVal(1,1:nSurfSample,1:nSurfSample,SurfMesh%nSides,nSpecies))

MacroSurfaceVal    = 0.
MacroSurfaceSpecVal= 0.

!> Erosion tracking
iSpec = 1
!---- Only one species. Only total values necessary
!===================================================================================================================================
DO iSurfSide=1,SurfMesh%nSides
    DO q=1,nSurfSample
        DO p=1,nSurfSample
        !---- 1. - .. / Impact Counter
        MacroSurfaceVal(1,p,q,iSurfSide) = SampWall(iSurfSide)%State(1,p,q)
        MacroSurfaceSpecVal(1,p,q,iSurfSide,iSpec) = SampWall(iSurfSide)%State(1,p,q) / TimeSample
        !---- 2. - .. / Impact Counter per AREA
        MacroSurfaceVal(2,p,q,iSurfSide) = SampWall(iSurfSide)%State(1,p,q) / SurfMesh%SurfaceArea(p,q,iSurfSide)
        !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
        MacroSurfaceVal(3,p,q,iSurfSide) = SampWall(iSurfSide)%State(2,p,q)
        MacroSurfaceVal(4,p,q,iSurfSide) = SampWall(iSurfSide)%State(3,p,q)
        MacroSurfaceVal(5,p,q,iSurfSide) = SampWall(iSurfSide)%State(4,p,q)
        MacroSurfaceVal(6,p,q,iSurfSide) = SampWall(iSurfSide)%State(6,p,q)
        !---- 7. - 10 / Impact angle (mean, min, max, variance)
        MacroSurfaceVal(7,p,q,iSurfSide) = SampWall(iSurfSide)%State(7,p,q)
        MacroSurfaceVal(8,p,q,iSurfSide) = SampWall(iSurfSide)%State(8,p,q)
        MacroSurfaceVal(9,p,q,iSurfSide) = SampWall(iSurfSide)%State(9,p,q)
        MacroSurfaceVal(10,p,q,iSurfSide)= SampWall(iSurfSide)%State(11,p,q)
        !---- 11 - 13 / Sampling Current Forces at walls
        MacroSurfaceVal(11,p,q,iSurfSide) = SampWall(iSurfSide)%State(12,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(12,p,q,iSurfSide) = SampWall(iSurfSide)%State(13,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(13,p,q,iSurfSide) = SampWall(iSurfSide)%State(14,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        SampWall(iSurfSide)%State(12,p,q) = 0.
!        SampWall(iSurfSide)%State(13,p,q) = 0.
!        SampWall(iSurfSide)%State(14,p,q) = 0.
        !---- 14 - 16 / Sampling Average Forces at walls
        MacroSurfaceVal(14,p,q,iSurfSide) = SampWall(iSurfSide)%State(15,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(15,p,q,iSurfSide) = SampWall(iSurfSide)%State(16,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(16,p,q,iSurfSide) = SampWall(iSurfSide)%State(17,p,q) / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
    END DO
  END DO
END DO

!---- Multiple species. All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
IF (nSpecies.GT.1) THEN
    DO iSurfSide=1,SurfMesh%nSides
    DO q=1,nSurfSample
    DO p=1,nSurfSample
    DO iSpec=1,nSpecies
        nShift    = iSpec * (nErosionVars-1)
        nShiftRHS = iSpec * nErosionVars
        !---- 1. - .. / Impact Counter
        MacroSurfaceVal(1+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(1+nShiftRHS,p,q)
        MacroSurfaceSpecVal(1,p,q,iSurfSide,iSpec)= SampWall(iSurfSide)%State(1+nShiftRHS,p,q) / TimeSample
        !---- 2. - .. / Impact Counter per AREA
        MacroSurfaceVal(2+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(1+nShiftRHS,p,q) / SurfMesh%SurfaceArea(p,q,iSurfSide)
        !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
        MacroSurfaceVal(3+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(2+nShiftRHS,p,q)
        MacroSurfaceVal(4+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(3+nShiftRHS,p,q)
        MacroSurfaceVal(5+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(4+nShiftRHS,p,q)
        MacroSurfaceVal(6+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(6+nShiftRHS,p,q)
        !---- 7. - 10 / Impact angle (mean, min, max, variance)
        MacroSurfaceVal(7+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(7+nShiftRHS,p,q)
        MacroSurfaceVal(8+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(8+nShiftRHS,p,q)
        MacroSurfaceVal(9+nShift,p,q,iSurfSide) = SampWall(iSurfSide)%State(9+nShiftRHS,p,q)
        MacroSurfaceVal(10+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(11+nShiftRHS,p,q)
        !---- 11 - 13 / Sampling Current Forces at walls
        MacroSurfaceVal(11+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(12+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(12+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(13+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        MacroSurfaceVal(13+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(14+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSample)
        !>> Set current forces to zero for new sampling run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        SampWall(iSurfSide)%State(12+nShiftRHS,p,q) = 0.
!        SampWall(iSurfSide)%State(13+nShiftRHS,p,q) = 0.
!        SampWall(iSurfSide)%State(14+nShiftRHS,p,q) = 0.
        !---- 14 - 16 / Sampling Average Forces at walls
        MacroSurfaceVal(14+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(15+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(15+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(16+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        MacroSurfaceVal(16+nShift,p,q,iSurfSide)= SampWall(iSurfSide)%State(17+nShiftRHS,p,q)                                      &
                                                  / (SurfMesh%SurfaceArea(p,q,iSurfSide) * t)
        END DO
      END DO
    END DO
  END DO
END IF

CALL WriteSurfSampleToHDF5(TRIM(MeshFile),ActualTime,remap_opt)

! Only deallocate if we don't need the values for wall calculations
IF (.NOT.doCalcWallParticles) THEN
    DEALLOCATE(MacroSurfaceVal,MacroSurfaceSpecVal)
ELSE
    DEALLOCATE(MacroSurfaceSpecVal)
END IF

END SUBROUTINE CalcSurfaceValues

END MODULE MOD_Particle_Erosion_Analyze
