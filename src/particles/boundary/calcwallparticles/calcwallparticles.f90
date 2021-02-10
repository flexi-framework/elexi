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
#include "particle.h"

!===================================================================================================================================
! Module for calculation of particle forces on surfaces for convergence analysis
!===================================================================================================================================
MODULE MOD_CalcWallParticles
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CalcWallParticles
  MODULE PROCEDURE CalcWallParticles
END INTERFACE

PUBLIC :: CalcWallParticles
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Control routine for CalcWallParticles
!==================================================================================================================================
SUBROUTINE CalcWallParticles(AlphaMean,AlphaVar,EkinMean,EkinVar,partForce,maxForce)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_AnalyzeEquation_Vars     ,ONLY: isWall
USE MOD_Mesh_Vars                ,ONLY: nBCs
USE MOD_Particle_Analyze_Vars    ,ONLY: TimeSample
USE MOD_Particle_Boundary_Vars   ,ONLY: SurfOnNode,nComputeNodeSurfSides
USE MOD_Particle_Boundary_Vars   ,ONLY: SurfSideArea,nSurfSample
USE MOD_Particle_Boundary_Vars   ,ONLY: SurfSide2GlobalSide
USE MOD_Particle_Boundary_Vars   ,ONLY: MacroSurfaceVal
USE MOD_Particle_Globals
USE MOD_Particle_Mesh_Vars       ,ONLY: SideInfo_Shared
#if USE_MPI
USE MOD_Particle_Boundary_Vars   ,ONLY: SampWallState_Shared
USE MOD_Particle_MPI_Shared_Vars ,ONLY: MPI_COMM_LEADERS_SURF,myComputeNodeRank
#else
USE MOD_Particle_Boundary_Vars   ,ONLY: SampWallState
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)               :: AlphaMean(nBCs)              !< average impact angle per wall BC
REAL,INTENT(OUT)               :: AlphaVar (nBCs)              !< average impact angle variance per wall BC
REAL,INTENT(OUT)               :: EkinMean (nBCs)              !< average kinetic energy per wall BC
REAL,INTENT(OUT)               :: EkinVar  (nBCs)              !< average kinetic energy variance per wall BC
REAL,INTENT(OUT)               :: partForce(nBCs)              !< integrated impact force per wall BC
REAL,INTENT(OUT)               :: maxForce (nBCs)              !< max impact force per wall BC
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: BCID,iSurfSide,SideID,p,q
REAL,ALLOCATABLE               :: SideArea(:)
REAL                           :: tmpForce,tmpForceX,tmpForceY,tmpForceZ
!==================================================================================================================================

! Do not try to record impacts if there are no walls on the current compute-node
IF(.NOT.SurfOnNode) RETURN

#if USE_MPI
! Only compute-node roots reduce the information
IF (myComputeNodeRank.NE.0) RETURN

ASSOCIATE(SampWallState => SampWallState_Shared)
#endif /*USE_MPI*/

ALLOCATE(SideArea(nBCs))
SideArea  = 0.
tmpForceX = 0.
tmpForceY = 0.
tmpForceZ = 0.
tmpForce  = 0.

AlphaMean = 0.
AlphaVar  = 0.
EkinMean  = 0.
EkinVar   = 0.
partForce = 0.
maxForce  = 0.

DO iSurfSide = 1,nComputeNodeSurfSides
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
  BCID   = SideInfo_Shared(    SIDE_BCID  ,SideID)

  IF(.NOT.isWall(BCID)) CYCLE

  DO q = 1,nSurfSample
    DO p = 1,nSurfSample
      ! Only consider walls with actual impacts
      IF (SampWallState(1,p,q,iSurfSide).EQ.0) CYCLE
      ! Average kinetic energy
      EkinMean(BCID)  = EkinMean(BCID)  + MacroSurfaceVal(3,p,q,iSurfSide) * SurfSideArea(p,q,iSurfSide)
      EkinVar (BCID)  = EkinVar (BCID)  + MacroSurfaceVal(6,p,q,iSurfSide) * SurfSideArea(p,q,iSurfSide)
      ! Average impact angle
      AlphaMean(BCID) = AlphaMean(BCID) + MacroSurfaceVal(7,p,q,iSurfSide) * SurfSideArea(p,q,iSurfSide)
      AlphaVar (BCID) = AlphaVar (BCID) + MacroSurfaceVal(10,p,q,iSurfSide)* SurfSideArea(p,q,iSurfSide)
      ! Integrated impact force
      IF (.NOT.ALMOSTZERO(TimeSample)) THEN
        tmpForceX   = MacroSurfaceVal(11,p,q,iSurfSide)
        tmpForceY   = MacroSurfaceVal(12,p,q,iSurfSide)
        tmpForceZ   = MacroSurfaceVal(13,p,q,iSurfSide)
      END IF
      tmpForce    = SQRT(tmpForceX**2 + tmpForceY**2 + tmpForceZ**2)

      partForce(BCID) = partForce(BCID) + tmpForce
      ! Max impact force
      maxForce(BCID)  = MAX(maxForce(BCID),tmpForce)

      SideArea(BCID)  = SideArea(BCID) + SurfSideArea(p,q,iSurfSide)
    END DO
  END DO
END DO

#if USE_MPI
END ASSOCIATE

CALL MPI_REDUCE(MPI_IN_PLACE,EkinMean ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_LEADERS_SURF,iError)
CALL MPI_REDUCE(MPI_IN_PLACE,EkinVar  ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_LEADERS_SURF,iError)
CALL MPI_REDUCE(MPI_IN_PLACE,AlphaMean,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_LEADERS_SURF,iError)
CALL MPI_REDUCE(MPI_IN_PLACE,AlphaVar ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_LEADERS_SURF,iError)
CALL MPI_REDUCE(MPI_IN_PLACE,partForce,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_LEADERS_SURF,iError)
CALL MPI_REDUCE(MPI_IN_PLACE,maxForce ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_LEADERS_SURF,iError)
CALL MPI_REDUCE(MPI_IN_PLACE,SideArea ,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_LEADERS_SURF,iError)
#endif /*USE_MPI*/

DO BCID = 1,nBCs
  IF(.NOT.isWall(BCID)) CYCLE

  IF(ALMOSTZERO(SideArea(BCID))) CYCLE

  EkinMean(BCID)  = EkinMean(BCID)  / SideArea(BCID)
  EkinVar (BCID)  = EkinVar (BCID)  / SideArea(BCID)
  AlphaMean(BCID) = AlphaMean(BCID) / SideArea(BCID)
  AlphaVar (BCID) = AlphaVar (BCID) / SideArea(BCID)
END DO

! Deallocate Macrosurfaces we couldn't deallocate earlier
SDEALLOCATE(MacroSurfaceVal)

END SUBROUTINE CalcWallParticles

END MODULE MOD_CalcWallParticles
