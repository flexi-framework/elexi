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

#define SIGMOID(x)         REAL(1.1/(1+EXP(-x))-0.1,4)
#define SILU(x,beta)       REAL(x/(1+EXP(-beta*x)),4)
#define LOGNORM(x,maxi)    REAL(LOG(ABS(x)+1.)/LOG(maxi+1.),4)
#define LOGNORMINV(x,maxi) REAL((maxi+1.)**x-1.,4)

!===================================================================================================================================
!> Determines how particles interact with a given boundary condition
!===================================================================================================================================
MODULE MOD_Particle_Fracture
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ParticleInsertingSingle
  MODULE PROCEDURE ParticleInsertingSingle
END INTERFACE

INTERFACE DiffuseReflectionFracture
  MODULE PROCEDURE DiffuseReflectionFracture
END INTERFACE

PUBLIC :: ParticleInsertingSingle
PUBLIC :: DiffuseReflectionFracture
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleInsertingSingle(FractNbr,iInit,NbrOfParticle,Particle_state,ParticleIndexNbr)
!===================================================================================================================================
! Insert a single particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_Particle_Vars          ,ONLY: doPartIndex,PartIndex
USE MOD_Particle_Vars          ,ONLY: Species,PDM,PartState,sumOfMatchedParticlesSpecies,PartSpecies,PartReflCount
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: FractNbr
INTEGER,INTENT(IN)    :: iInit
INTEGER,INTENT(INOUT) :: NbrOfParticle
REAL,INTENT(IN)       :: Particle_state(1:PP_nVarPart)
INTEGER,INTENT(INOUT) :: ParticleIndexNbr
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER               :: InitGroup
#endif
!===================================================================================================================================

Species(FractNbr)%Init(iInit)%sumOfRequestedParticles = 1
Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles = 0
Species(FractNbr)%Init(iInit)%sumOfMatchedParticles   = 0

#if USE_MPI
! emission group communicator for the current iInit
InitGroup = Species(FractNbr)%Init(iInit)%InitCOMM
IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
  NbrofParticle=0
  RETURN
END IF
#endif /*USE_MPI*/

! Find a free position in the PDM array
ParticleIndexNbr = PDM%nextFreePosition(Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles + 1 + PDM%CurrentNextFreePosition)
IF (ParticleIndexNbr.NE.0) THEN
  PartState(PART_POSV,ParticleIndexNbr) = Particle_state(PART_POSV)
  PDM%ParticleInside( ParticleIndexNbr) = .TRUE.
  CALL LocateParticleInElement(ParticleIndexNbr,doHALO=.TRUE.)
  IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
    Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles = Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles + 1
    PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
  END IF
  PartState(PART_VELV,ParticleIndexNbr) = Particle_state(PART_VELV)
  PartState(PART_DIAM,ParticleIndexNbr) = Particle_state(PART_DIAM)
  PartSpecies(ParticleIndexNbr)         = FractNbr
  PartReflCount(ParticleIndexNbr)       = 0
ELSE
  CALL Abort(__STAMP__,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
END IF

#if USE_MPI
! Start communicating matched particles. This routine is finished in particle_emission.f90
CALL MPI_IALLREDUCE( Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles &
                   , Species(FractNbr)%Init(iInit)%sumOfMatchedParticles   &
                   , 1                                                     &
                   , MPI_INTEGER                                           &
                   , MPI_SUM                                               &
!                   , 0                                                     &
                   , PartMPI%InitGroup(InitGroup)%COMM                     &
                   , PartMPI%InitGroup(InitGroup)%Request                  &
                   , IERROR)

IF (doPartIndex) THEN
  Species(FractNbr)%Init(iInit)%nPartsPerProc=0
  CALL MPI_IEXSCAN( Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles &
                  , Species(FractNbr)%Init(iInit)%nPartsPerProc           &
                  , 1                                                     &
                  , MPI_INTEGER                                           &
                  , MPI_SUM                                               &
                  , PartMPI%InitGroup(InitGroup)%COMM                     &
                  , PartMPI%InitGroup(InitGroup)%RequestIndex             &
                  , IERROR)
END IF
#else
! in the seriell case, particles are only emitted on the current proc
Species(FractNbr)%Init(iInit)%sumOfMatchedParticles = Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles
! Assign PartIndex
IF (doPartIndex) Species(FractNbr)%Init(iInit)%nPartsPerProc = 0
#endif /*USE_MPI*/

! Return the *local* NbrOfParticle so that the following Routines only fill in
! the values for the local particles
NbrOfParticle = Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles

IF(NbrofParticle.EQ.0) RETURN

! update number of particles on proc and find next free position in particle array
PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
CALL UpdateNextFreePosition()

#if USE_MPI
InitGroup = Species(FractNbr)%Init(iInit)%InitCOMM
IF (PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL .AND. Species(FractNbr)%Init(iInit)%sumOfRequestedParticles.GT.0) THEN
  CALL MPI_WAIT(PartMPI%InitGroup(InitGroup)%Request, MPI_STATUS_IGNORE, iError)
  IF(doPartIndex) CALL MPI_WAIT(PartMPI%InitGroup(InitGroup)%RequestIndex, MPI_STATUS_IGNORE, iError)
END IF
#endif

IF (doPartIndex) THEN
  PartIndex(ParticleIndexNbr)  = sumOfMatchedParticlesSpecies + Species(FractNbr)%Init(iInit)%nPartsPerProc + 1
  sumOfMatchedParticlesSpecies = sumOfMatchedParticlesSpecies + Species(FractNbr)%Init(iInit)%sumOfMatchedParticles
END IF

END SUBROUTINE ParticleInsertingSingle


SUBROUTINE DiffuseReflectionFracture(PartTrajectory,lengthPartTrajectory,alpha,PartID,SideID,n_loc,tang1,tang2,&
    eps_n,eps_t1,eps_t2,dp_old,WallVelo)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
#if USE_PARTROT
USE MOD_Mathtools                  ,ONLY: CROSS
#endif
USE MOD_Mesh_Vars                  ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars     ,ONLY: doParticleReflectionTrack
USE MOD_Particle_Boundary_Vars     ,ONLY: LowVeloRemove
USE MOD_Particle_Boundary_Vars     ,ONLY: doParticleImpactTrack
USE MOD_Particle_Boundary_Tracking ,ONLY: StoreBoundaryParticleProperties
USE MOD_Particle_Vars              ,ONLY: PartState,LastPartPos,Species,PartSpecies,PartReflCount
USE MOD_Particle_Vars              ,ONLY: PDM,PEM
USE MOD_Part_Operations            ,ONLY: CreateParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(INOUT)                :: n_loc(1:3),tang1(1:3),tang2(1:3)
REAL,INTENT(IN)                   :: eps_n,eps_t1,eps_t2,dp_old,WallVelo(3)
INTEGER,INTENT(IN)                :: PartID,SideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: v_old(1:3)
REAL                              :: PartFaceAngle,PartFaceAngle_old
REAL                              :: PartTrajectoryTang1(3),PartTrajectoryTang2(3),PartTrajectoryNorm(3)
REAL                              :: v_norm(3),v_tang1(3),v_tang2(3)
REAL                              :: intersecRemain,v_magnitude
#if USE_PARTROT
REAL                              :: rot_old(1:3)
#endif
REAL                              :: PartState_tmp(1:PP_nVarPart),PartTrajectory_tmp(1:3)
INTEGER                           :: NbrofParticle,NewPartID
!===================================================================================================================================

! Calculate wall normal and tangential velocities, impact angle
v_norm   = DOT_PRODUCT(PartState(PART_VELV,PartID),n_loc)*n_loc
v_tang1  = DOT_PRODUCT(PartState(PART_VELV,PartID),tang1)*tang1
v_tang2  = DOT_PRODUCT(PartState(PART_VELV,PartID),tang2)*tang2
PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

! Make sure we have the old values safe
IF (doParticleImpactTrack) PartFaceAngle_old = PartFaceAngle
v_old = PartState(PART_VELV,PartID)
#if USE_PARTROT
rot_old = PartState(PART_AMOMV,PartID)
#endif

!--> flip trajectory and move the remainder of the particle push
! Calculate wall normal and tangential velocity components after impact, rescale to uniform length
PartTrajectoryNorm (1:3) = eps_n  * (DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc)
PartTrajectoryTang1(1:3) = eps_t1 * (DOT_PRODUCT(PartTrajectory(1:3),tang1)*tang1)
PartTrajectoryTang2(1:3) = eps_t2 * (DOT_PRODUCT(PartTrajectory(1:3),tang2)*tang2)
PartTrajectory_tmp(1:3)  = PartTrajectoryTang1(1:3) + PartTrajectoryTang2(1:3) - PartTrajectoryNorm(1:3)

! Rescale the remainder to the new length
intersecRemain = (lengthPartTrajectory - alpha)
intersecRemain = SQRT(eps_n*eps_n + eps_t1*eps_t1 + eps_t2*eps_t2)/SQRT(3.) * (lengthPartTrajectory - alpha)
! Compute the remainder of the new particle trajectory
PartTrajectory_tmp = intersecRemain * PartTrajectory_tmp/SQRT(SUM(PartTrajectory_tmp**2.))

! Compute moved particle || rest of movement. PartTrajectory_tmp has already been updated
PartState_tmp(PART_POSV) = LastPartPos(1:3,PartID) + PartTrajectory_tmp(1:3)

! Compute new particle velocity, modified with coefficents of restitution
PartState_tmp(PART_VELV) = eps_t1 * v_tang1 + eps_t2 * v_tang2 - eps_n * v_norm + WallVelo

! Calculate new particle diameter
PartState_tmp(PART_DIAM) = VOL_SPHERE_INV((VOL_SPHERE(dp_old)-VOL_SPHERE(PartState(PART_DIAM,PartID))))

#if USE_PARTROT
! rotation: I_2*w_2-I_1*w_1 = - r x J , J=m_2*v_2-m_1*v_1, r=d/2*n_loc
! for a constant particle volume: I_2=I_1, m_1=m_2
PartState_tmp(PART_AMOMV) = (dp_old/PartState_tmp(PART_DIAM))**5*rot_old - &
                            0.5*dp_old*Species(PartSpecies(PartID))%DensityIC*PI/6*&
                            CROSS(n_loc,(PartState_tmp(PART_DIAM)**3*PartState_tmp(PART_VELV)-v_old*dp_old**3))
#endif

! Insert new particle
NbrofParticle = 1
!CALL ParticleInsertingSingle(PartSpecies(PartID),0,NbrofParticle,PartState_tmp,NewPartID)
CALL CreateParticle(PartSpecies(PartID),PartState_tmp(:),PEM%Element(PartID),PartID,&
                    LastPartPos(1:3,PartID),PEM%Element(PartID),NewPartID)

! No new particle is inserted
IF (NbrofParticle .EQ. 0) THEN
  PartState(PART_DIAM,PartID) = dp_old
  RETURN
END IF

IF (doParticleImpactTrack) THEN
  PartFaceAngle = ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,n_loc)))

  CALL StoreBoundaryParticleProperties(BCSideID          = SideInfo_Shared(SIDE_BCID,SideID)             &
                                      ,PartID            = NewPartID                                     &
                                      ,PartFaceAngle     = PartFaceAngle                                 &
                                      ,v_old             = v_old                                         &
                                      ,PartFaceAngle_old = PartFaceAngle_old                             &
                                      ,PartReflCount     = PartReflCount(NewPartID)                      &
                                      ,alpha             = alpha                                         &
                                      ,dp_old            = dp_old                                        &
#if USE_PARTROT
                                      ,rot_old           = rot_old                                       &
#endif
                                      )
END IF

! increase reflection counter
IF (doParticleReflectionTrack) PartReflCount(NewPartID) = PartReflCount(NewPartID) + 1

! Correction for Runge-Kutta. Reflect or delete force history / acceleration
! Coefficients of Restitution cause non-differentiable jump in velocity. Always erase force history and reconstruction in timedisc during push
PDM%IsNewPart(NewPartID) = .TRUE.

! Remove sliding low velocity particles
IF (LowVeloRemove) THEN
  v_magnitude = SQRT(DOT_PRODUCT(PartState(PART_VELV,NewPartID),PartState(PART_VELV,NewPartID)))

  IF ((Species(PartSpecies(NewPartID))%LowVeloThreshold.NE.0).AND.&
    (v_magnitude.LT.Species(PartSpecies(NewPartID))%LowVeloThreshold)) THEN
    Species(PartSpecies(NewPartID))%LowVeloCounter = Species(PartSpecies(NewPartID))%LowVeloCounter + 1
    PDM%ParticleInside(NewPartID) = .FALSE.
    IPWRITE(UNIT_stdOut,*) ' Low velocity particle removed after impact. Velocity after reflection:', v_magnitude
  END IF
END IF

END SUBROUTINE DiffuseReflectionFracture



END MODULE MOD_Particle_Fracture
