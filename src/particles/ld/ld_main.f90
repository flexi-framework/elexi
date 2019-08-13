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

!===================================================================================================================================
! module including low diffusion model
!===================================================================================================================================
MODULE MOD_LD
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE LD_main
  MODULE PROCEDURE LD_main
END INTERFACE
INTERFACE LD_reposition
  MODULE PROCEDURE LD_reposition
END INTERFACE
INTERFACE LD_PerfectReflection
  MODULE PROCEDURE LD_PerfectReflection
END INTERFACE
INTERFACE LD_SetParticlePosition
  MODULE PROCEDURE LD_SetParticlePosition
END INTERFACE

PUBLIC :: LD_main, LD_reposition, LD_PerfectReflection, LD_SetParticlePosition
!===================================================================================================================================

CONTAINS


SUBROUTINE LD_main()
!===================================================================================================================================
! Main LD routine
!===================================================================================================================================
! MODULES
USE MOD_LD_Vars
USE MOD_Mesh_Vars,             ONLY : nElems, nSides
USE MOD_Particle_Vars,         ONLY : PDM, PEM
USE MOD_LD_mean_cell,          ONLY : CalcMacCellLDValues
USE MOD_LD_lag_velo,           ONLY : CalcSurfLagVelo
USE MOD_LD_reassign_part_prop, ONLY : LD_reassign_prop
USE MOD_LD_part_treat,         ONLY : LDPartTreament
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER           :: iElem
!--------------------------------------------------------------------------------------------------!
  LD_RHS(1:PDM%ParticleVecLength,1) = 0.0
  LD_RHS(1:PDM%ParticleVecLength,2) = 0.0
  LD_RHS(1:PDM%ParticleVecLength,3) = 0.0
  IsDoneLagVelo(1:nSides) = .FALSE.
  CALL CalcMacCellLDValues
  CALL CalcSurfLagVelo
  IF(LD_RepositionFak.NE. 0) THEN
    CALL LD_reposition()
  END IF
  DO iElem = 1, nElems
#if (PP_TimeDiscMethod==1001)
  IF((BulkValues(iElem)%CellType.EQ.3).OR.(BulkValues(iElem)%CellType.EQ.4)) THEN  ! --- LD Cell ?
#endif
    IF (PEM%pNumber(iElem).GT. 1) THEN
      CALL LD_reassign_prop(iElem)
      CALL LDPartTreament(iElem)
    END IF
#if (PP_TimeDiscMethod==1001)
  END IF  ! --- END LD Cell?
#endif
  END DO
  
END SUBROUTINE LD_main
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!                          _   _     _   _   _   _   _                                             !
!                         / \ / \   / \ / \ / \ / \ / \                                            !
!                        ( L | D ) ( T | O | O | L | S )                                           !
!                         \_/ \_/   \_/ \_/ \_/ \_/ \_/                                            !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
SUBROUTINE LD_reposition
!===================================================================================================================================
! reposition calculation for LD
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars,         ONLY : PartState, PEM
  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_Eval_xyz,              ONLY : EvaluateFieldAtRefPos
  USE MOD_Mesh_Vars,             ONLY : NGeo!
  USE MOD_Particle_Mesh_Vars,       ONLY:XCL_NGeo,XiCL_Ngeo,wBaryCL_NGeo
  USE MOD_LD_Vars,               ONLY : LD_RepositionFak
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration  
INTEGER               :: iElem
INTEGER               :: iPart, iPartIndx,nPart
REAL                  :: RandVec(3), iRan
!--------------------------------------------------------------------------------------------------!

  DO iElem = 1, nElems
#if (PP_TimeDiscMethod==1001)
  IF((BulkValues(iElem)%CellType.EQ.3).OR.(BulkValues(iElem)%CellType.EQ.4)) THEN  ! --- LD Cell ?
#endif
    nPart     = PEM%pNumber(iElem)
    iPartIndx = PEM%pStart(iElem)
    DO iPart = 1, nPart
      CALL RANDOM_NUMBER(iRan)
      IF (iRan.LT. LD_RepositionFak) THEN     
        CALL RANDOM_NUMBER(RandVec)
        RandVec = RandVec * 2.0 - 1.0
        CALL EvaluateFieldAtRefPos(RandVec,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),PartState(iPartIndx,1:3))
        iPartIndx = PEM%pNext(iPartIndx)
      END IF
    END DO
#if (PP_TimeDiscMethod==1001)
  END IF  ! --- END LD Cell?
#endif
  END DO

END SUBROUTINE LD_reposition

!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_PerfectReflection(nx,ny,nz,xNod,yNod,zNod,PoldStarX,PoldStarY,PoldStarZ,i)
!===================================================================================================================================
! reflection at wall for LD case
!===================================================================================================================================
! MODULES
  USE MOD_LD_Vars
  USE MOD_Particle_Vars,         ONLY : lastPartPos
  USE MOD_TimeDisc_Vars,         ONLY : dt
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! Local variable declaration                                                                       !
   REAL                             :: PnewX, PnewY, PnewZ                                  !
   REAL                             :: bx,by,bz, ax,ay,az, dist                                    !
   REAL                             :: PnewStarX, PnewStarY, PnewStarZ, Velo                       !
   REAL                             :: VelX, VelY, VelZ, NewVelocity                               !
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  INTEGER, INTENT(IN)               :: i
  REAL, INTENT(IN)                  :: nx,ny,nz,xNod,yNod,zNod,PoldStarX,PoldStarY,PoldStarZ
!--------------------------------------------------------------------------------------------------!

   PnewX = lastPartPos(i,1) + PartStateBulkValues(i,1) * dt
   PnewY = lastPartPos(i,2) + PartStateBulkValues(i,2) * dt
   PnewZ = lastPartPos(i,3) + PartStateBulkValues(i,3) * dt

   bx = PnewX - xNod
   by = PnewY - yNod
   bz = PnewZ - zNod

   ax = bx - nx * (bx * nx + by * ny + bz * nz)
   ay = by - ny * (bx * nx + by * ny + bz * nz)
   az = bz - nz * (bx * nx + by * ny + bz * nz)

   dist = SQRT(((ay * bz - az * by) * (ay * bz - az * by) +   &
        (az * bx - ax * bz) * (az * bx - ax * bz) +   &
        (ax * by - ay * bx) * (ax * by - ay * bx))/   &
        (ax * ax + ay * ay + az * az))

!   If vector from old point to new point goes through the node, a will be zero
!   dist is then simply length of vector b instead of |axb|/|a|
   IF (dist.NE.dist) dist = SQRT(bx*bx+by*by+bz*bz)

   PnewStarX = PnewX - 2 * dist * nx
   PnewStarY = PnewY - 2 * dist * ny
   PnewStarZ = PnewZ - 2 * dist * nz

   !---- Calculate new velocity vector

   Velo = SQRT(PartStateBulkValues(i,1) * PartStateBulkValues(i,1) + &
               PartStateBulkValues(i,2) * PartStateBulkValues(i,2) + &
               PartStateBulkValues(i,3) * PartStateBulkValues(i,3))

   VelX = PnewStarX - PoldStarX
   VelY = PnewStarY - PoldStarY
   VelZ = PnewStarZ - PoldStarZ

   NewVelocity = SQRT(VelX * VelX + VelY * VelY + VelZ * VelZ)

   VelX = VelX/NewVelocity * Velo
   VelY = VelY/NewVelocity * Velo
   VelZ = VelZ/NewVelocity * Velo

   !---- Assign new values to "old" variables to continue loop

   PartStateBulkValues(i,1)   = VelX 
   PartStateBulkValues(i,2)   = VelY
   PartStateBulkValues(i,3)   = VelZ

END SUBROUTINE LD_PerfectReflection

!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_SetParticlePosition(chunkSize,particle_positions_Temp,iSpec,iInit)
!===================================================================================================================================
! modified particle emmission for LD case
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars,         ONLY:Species,PDM
  USE MOD_Mesh_Vars,             ONLY:nElems
  USE MOD_Particle_Mesh_Vars,    ONLY:XiCL_NGeo,wBaryCL_NGeo
  USE MOD_Particle_Mesh_Vars,    ONLY:GEO
  USE MOD_Eval_xyz,              ONLY:EvaluateFieldAtRefPos
  USE MOD_Mesh_Vars,             ONLY:NGeo
  USE MOD_Particle_Mesh_Vars,    ONLY:XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! Local variable declaration                                                                       !
!--------------------------------------------------------------------------------------------------!
INTEGER               :: iElem, ichunkSize
INTEGER               :: iPart, nPart
REAL                  :: RandVec(3), iRan, RandomPos(3)
REAL                  :: PartDens, FractNbr
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)               :: iSpec
INTEGER, INTENT(IN)               :: iInit
!--------------------------------------------------------------------------------------------------!
! INOUTPUT VARIABLES
INTEGER, INTENT(INOUT)           :: chunkSize
!--------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,ALLOCATABLE, INTENT(OUT)    :: particle_positions_Temp(:)
!--------------------------------------------------------------------------------------------------!

  ALLOCATE(particle_positions_Temp(3*PDM%maxParticleNumber))
  particle_positions_Temp=0.
  PartDens = Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor   ! numerical Partdensity is needed
  ichunkSize = 1
  DO iElem = 1, nElems
    FractNbr = PartDens * GEO%Volume(iElem) - AINT(PartDens * GEO%Volume(iElem))
    CALL RANDOM_NUMBER(iRan)
    IF (iRan .GT. FractNbr) THEN
      nPart = INT(AINT(PartDens * GEO%Volume(iElem)))
    ELSE
      nPart = INT(AINT(PartDens * GEO%Volume(iElem))) + 1
    END IF
    DO iPart = 1, nPart   
      CALL RANDOM_NUMBER(RandVec)
      RandVec = RandVec * 2.0 - 1.0
      CALL EvaluateFieldAtRefPos(RandVec,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(:,:,:,:,iElem),RandomPos)
      !CALL TensorProductInterpolation(RandomPos,xi,iElem)
      !IF(ANY(ABS(Xi).GT.1.0)) THEN
      !  print*,'upppsss'
      !  print*,'randvec',randvec
      !  print*,'xi',xi
      !  read*
      !END IF
      particle_positions_Temp(ichunkSize*3-2) = RandomPos(1)
      particle_positions_Temp(ichunkSize*3-1) = RandomPos(2)
      particle_positions_Temp(ichunkSize*3)   = RandomPos(3)
      ichunkSize = ichunkSize + 1
    END DO
  END DO
  chunkSize = ichunkSize - 1

END SUBROUTINE LD_SetParticlePosition

END MODULE MOD_LD
