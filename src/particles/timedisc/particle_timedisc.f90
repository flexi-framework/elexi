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
!> Module for the particle temporal discretization
!==================================================================================================================================
MODULE MOD_Particle_TimeDisc
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Particle_TimeStepByEuler
  MODULE PROCEDURE Particle_TimeStepByEuler
END INTERFACE

INTERFACE Particle_TimeStepByLSERK
  MODULE PROCEDURE Particle_TimeStepByLSERK
END INTERFACE Particle_TimeStepByLSERK

INTERFACE Particle_TimeStepByLSERK_RHS
  MODULE PROCEDURE Particle_TimeStepByLSERK_RHS
END INTERFACE Particle_TimeStepByLSERK_RHS

INTERFACE Particle_TimeStepByLSERK_RK
  MODULE PROCEDURE Particle_TimeStepByLSERK_RK
END INTERFACE Particle_TimeStepByLSERK_RK

INTERFACE Particle_TimeStepByLSERK_RK_RHS
  MODULE PROCEDURE Particle_TimeStepByLSERK_RK_RHS
END INTERFACE Particle_TimeStepByLSERK_RK_RHS

PUBLIC::Particle_TimeStepByEuler
PUBLIC::Particle_TimeStepByLSERK
PUBLIC::Particle_TimeStepByLSERK_RHS
PUBLIC::Particle_TimeStepByLSERK_RK
PUBLIC::Particle_TimeStepByLSERK_RK_RHS

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Euler particle time integration:
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByEuler(dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_TimeDisc_Vars,           ONLY: t
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Part_RHS,                ONLY: CalcPartRHS
USE MOD_Part_Tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Interpolation,  ONLY: InterpolateFieldToParticle
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping,TriaTracking
USE MOD_Particle_Vars,           ONLY: Species, PartSpecies, PartState, Pt, LastPartPos, DelayTime, PEM, PDM
#if USE_RW
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
#endif
#if USE_LOADBALANCE
USE MOD_Particle_Localization,   ONLY: CountPartsPerElem
USE MOD_LoadBalance_Tools,       ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: dt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_MPI
#if USE_LOADBALANCE
! Needed for scaling parts and load balance
CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.)
#endif

IF (t.GE.DelayTime) THEN
  PartMPIExchange%nMPIParticles=0
END IF
#endif /*MPI*/

! set last particle position and element
LastPartPos(  1,1:PDM%ParticleVecLength) = PartState(1,1:PDM%ParticleVecLength)
LastPartPos(  2,1:PDM%ParticleVecLength) = PartState(2,1:PDM%ParticleVecLength)
LastPartPos(  3,1:PDM%ParticleVecLength) = PartState(3,1:PDM%ParticleVecLength)
PEM%lastElement(1:PDM%ParticleVecLength) = PEM%Element(1:PDM%ParticleVecLength)

! forces on particles
IF (t.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  CALL InterpolateFieldToParticle()
#if USE_RW
  CALL ParticleRandomWalk(t)
#endif
  CALL CalcPartRHS()
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF

! particle push using Euler
IF (t.GE.DelayTime) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      !-- Tracer Particles
      IF (TRIM(Species(PartSpecies(iPart))%RHSMethod).EQ.'Tracer') THEN
        PartState(1,iPart) = PartState(1,iPart) + PartState(4,iPart)*dt
        PartState(2,iPart) = PartState(2,iPart) + PartState(5,iPart)*dt
        PartState(3,iPart) = PartState(3,iPart) + PartState(6,iPart)*dt
        PartState(4,iPart) = Pt       (1,iPart)
        PartState(5,iPart) = Pt       (2,iPart)
        PartState(6,iPart) = Pt       (3,iPart)
      !-- Normal particles
      ELSE
        !-- Particle Push
        !-> Sanity Check Particle Pusher / WARNING: Might Cause Slowdowns
        !IF (ANY(ISNAN(Pt(:,iPart)))) THEN
        !    IPWRITE(UNIT_stdOut,*) 'Found invalid particle push, ignoring. PartID:', iPart
        !    Pt(iPart,:) = 0
        !ENDIF

        PartState(1,iPart) = PartState(1,iPart) + PartState(4,iPart)*dt
        PartState(2,iPart) = PartState(2,iPart) + PartState(5,iPart)*dt
        PartState(3,iPart) = PartState(3,iPart) + PartState(6,iPart)*dt
        PartState(4,iPart) = PartState(4,iPart) + Pt       (1,iPart)*dt
        PartState(5,iPart) = PartState(5,iPart) + Pt       (2,iPart)*dt
        PartState(6,iPart) = PartState(6,iPart) + Pt       (3,iPart)*dt

        !-- Sanity Check Particle / WARNING: Might Cause Slowdowns
        !IF (ANY(ISNAN(PartState(iPart,:)))) THEN
        !    PDM%ParticleInside(iPart) = .FALSE.
        !    IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        !ENDIF

      ENDIF !< Tracer
    ENDIF !< ParticleInside
  END DO
END IF
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

! Locate and communicate particles (only if particle have changed)
IF (t.GE.DelayTime) THEN
#if USE_MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif
  ! track new particle position
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    IF (TriaTracking) THEN
      CALL ParticleTriaTracking()
    ELSE
      CALL ParticleTracing()
    END IF
  END IF
  ! emitt particles inserted in current time step
  CALL ParticleInserting()
#if USE_MPI
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! receive particles, locate and finish communication
  CALL MPIParticleRecv()
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif
  ! find next free position in particle array
  CALL UpdateNextFreePosition()
END IF

END SUBROUTINE Particle_TimeStepByEuler

!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> Calculate the right hand side before updating the field solution. Can be used to hide sending of number of particles.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK_RHS(t)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles,MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
USE MOD_Part_RHS,                ONLY: CalcPartRHS
USE MOD_Particle_Interpolation,  ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY: PartState,DelayTime,LastPartPos,PDM,PEM
#if USE_RW
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
#endif
#if USE_LOADBALANCE
USE MOD_Particle_Localization,   ONLY: CountPartsPerElem
USE MOD_LoadBalance_Tools,       ONLY: LBStartTime,LBPauseTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
#if USE_LOADBALANCE
! Needed for scaling parts and load balance
CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.)
#endif

IF (t.GE.DelayTime) THEN
  PartMPIExchange%nMPIParticles=0
END IF
#endif /*USE_MPI*/

! set last particle position and element
LastPartPos(  1,1:PDM%ParticleVecLength) = PartState(1,1:PDM%ParticleVecLength)
LastPartPos(  2,1:PDM%ParticleVecLength) = PartState(2,1:PDM%ParticleVecLength)
LastPartPos(  3,1:PDM%ParticleVecLength) = PartState(3,1:PDM%ParticleVecLength)
PEM%lastElement(1:PDM%ParticleVecLength) = PEM%Element(1:PDM%ParticleVecLength)

IF (t.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! forces on particle
  ! can be used to hide sending of number of particles
  !--> Interpolate fluid field to particle position
  CALL InterpolateFieldToParticle()
#if USE_RW
  !--> Calculate the random walk push
  CALL ParticleRandomWalk(t)
#endif
  !--> Calculate the particle right hand side and push
  CALL CalcPartRHS()
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF

END SUBROUTINE Particle_TimeStepByLSERK_RHS


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK(t,b_dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_TimeDisc_Vars,           ONLY: nRKStages
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_PruettDamping,           ONLY: TempFilterTimeDeriv
USE MOD_HDF5_Output,             ONLY: WriteState
#if FV_ENABLED
USE MOD_FV,                      ONLY: FV_Switch
USE MOD_FV_Vars,                 ONLY: FV_toDGinRK
#endif
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Part_Tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping,TriaTracking
USE MOD_Particle_Vars,           ONLY: PartState, Pt, Pt_temp, DelayTime, PDM
#if USE_RW
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Tools,       ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: b_dt(1:nRKStages)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: part_err
INTEGER                       :: iPart
!REAL                          :: v_magnitude
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

IF (t.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  part_err = .FALSE.

  ! particle push for first RK stage
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      ! Pt is always known at this position, change isNewPart to false
      PDM%IsNewPart(iPart)=.FALSE.

      !-- Particle Push
      ! Sanity Check Particle Pusher / WARNING: Might Cause Slowdowns
      !IF (ANY(ISNAN(Pt(iPart,:)))) THEN
      !    IPWRITE(UNIT_stdOut,*) 'Found invalid particle push, ignoring. PartID:', iPart
      !    Pt(iPart,:) = 0
      !ENDIF

      Pt_temp  (1,iPart) = PartState(4,iPart)
      Pt_temp  (2,iPart) = PartState(5,iPart)
      Pt_temp  (3,iPart) = PartState(6,iPart)
      Pt_temp  (4,iPart) = Pt       (1,iPart)
      Pt_temp  (5,iPart) = Pt       (2,iPart)
      Pt_temp  (6,iPart) = Pt       (3,iPart)
      PartState(1,iPart) = PartState(1,iPart) + PartState(4,iPart)*b_dt(1)
      PartState(2,iPart) = PartState(2,iPart) + PartState(5,iPart)*b_dt(1)
      PartState(3,iPart) = PartState(3,iPart) + PartState(6,iPart)*b_dt(1)
      PartState(4,iPart) = PartState(4,iPart) + Pt       (1,iPart)*b_dt(1)
      PartState(5,iPart) = PartState(5,iPart) + Pt       (2,iPart)*b_dt(1)
      PartState(6,iPart) = PartState(6,iPart) + Pt       (3,iPart)*b_dt(1)

      ! Sanity Check Particle / WARNING: Might Cause Slowdowns
      !IF (ANY(ISNAN(PartState(:,iPart)))) THEN
      !    PDM%ParticleInside(iPart) = .FALSE.
      !    IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
      !ENDIF

      ! Try to find particles with too high velocity
      !v_magnitude   = SQRT(DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart)))

      !IF ((Species(PartSpecies(iPart))%HighVeloThreshold.NE.0).AND.(v_magnitude.GT.Species(PartSpecies(iPart))%HighVeloThreshold))THEN
      !  part_err = .TRUE.
      !  IPWRITE(UNIT_stdOut,*) ' High velocity particle detected. Writing error state and removing particle ...'
      !  IPWRITE(UNIT_stdout,*) ' LastPos:',  PartState(1:3,iPart)
      !  IPWRITE(UNIT_stdout,*) ' Velocity:', PartState(4:6,iPart)
      !  PDM%ParticleInside(iPart) = .FALSE.
      !END IF

    END IF
  END DO

  !IF (part_err) THEN
  !    CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,FutureTime=tWriteData,isErrorFile=.TRUE.)
  !END IF
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

#if USE_MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif
  ! track new particle position
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    IF (TriaTracking) THEN
      CALL ParticleTriaTracking()
    ELSE
      CALL ParticleTracing()
    END IF
  END IF
  ! emitt particles inserted in current time step
  CALL ParticleInserting()
#if USE_MPI
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! find next free position in particle array
  CALL MPIParticleRecv()
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif
END IF

END SUBROUTINE Particle_TimeStepByLSERK

!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> Calculate the right hand side before updating the field solution. Can be used to hide sending of number of particles.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK_RK_RHS(t)
! MODULES
USE MOD_Globals
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
USE MOD_Particle_Interpolation,  ONLY: InterpolateFieldToParticle
USE MOD_Part_RHS,                ONLY: CalcPartRHS
USE MOD_Particle_Vars,           ONLY: PartState,DelayTime,LastPartPos,PDM,PEM
#if USE_RW
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
#endif
#if USE_LOADBALANCE
USE MOD_Particle_Localization,   ONLY: CountPartsPerElem
USE MOD_LoadBalance_Tools,       ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_MPI
#if USE_LOADBALANCE
! Needed for scaling parts and load balance
CALL CountPartsPerElem(ResetNumberOfParticles=.FALSE.)
#endif

IF (t.GE.DelayTime) THEN
  PartMPIExchange%nMPIParticles=0
END IF
#endif /*USE_MPI*/

! set last particle position and element
LastPartPos(  1,1:PDM%ParticleVecLength) = PartState(1,1:PDM%ParticleVecLength)
LastPartPos(  2,1:PDM%ParticleVecLength) = PartState(2,1:PDM%ParticleVecLength)
LastPartPos(  3,1:PDM%ParticleVecLength) = PartState(3,1:PDM%ParticleVecLength)
PEM%lastElement(1:PDM%ParticleVecLength) = PEM%Element(1:PDM%ParticleVecLength)

IF (t.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! forces on particle
  ! can be used to hide sending of number of particles
  !--> Interpolate fluid field to particle position
  CALL InterpolateFieldToParticle()   ! forces on particles
#if USE_RW
  !--> Calculate the random walk push
  CALL ParticleRandomWalk(t)
#endif
  !--> Calculate the particle right hand side and push
  CALL CalcPartRHS()
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF

END SUBROUTINE Particle_TimeStepByLSERK_RK_RHS

!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version. Inner RK iteration
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK_RK(t,iStage,b_dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_TimeDisc_Vars,           ONLY: RKA,nRKStages
USE MOD_HDF5_Output,             ONLY: WriteState
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Part_RHS,                ONLY: CalcPartRHS
USE MOD_Part_Tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Interpolation,  ONLY: InterpolateFieldToParticle
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_Vars,  ONLY: DoRefMapping,TriaTracking
USE MOD_Particle_Vars,           ONLY: PartState, Pt, Pt_temp, DelayTime, PDM
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
#if USE_RW
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Tools,       ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
INTEGER,INTENT(IN)            :: iStage
REAL,INTENT(IN)               :: b_dt(1:nRKStages)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: part_err
INTEGER                       :: iPart, iStage_loc
REAL                          :: RandVal!,v_magnitude
REAL                          :: Pa_rebuilt_coeff(1:nRKStages),Pa_rebuilt(1:3,1:nRKStages),Pv_rebuilt(1:3,1:nRKStages),v_rebuilt(1:3,0:nRKStages-1)
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! Rebuild Pt_tmp-coefficients assuming F=const. (value at wall) in previous stages
DO iStage_loc=1,nRKStages
  IF (iStage_loc.EQ.1) THEN
    Pa_rebuilt_coeff(iStage_loc) = 1.
  ELSE
    Pa_rebuilt_coeff(iStage_loc) = 1. - RKA(iStage_loc)*Pa_rebuilt_coeff(iStage_loc-1)
  END IF
END DO

IF (t.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  part_err = .FALSE.

  ! particle push for nth RK stage
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
        ! "normal" particles are pushed with whole timestep
      IF (.NOT.PDM%IsNewPart(iPart)) THEN
        Pt_temp(1,iPart) = PartState(4,iPart) - RKA(iStage) * Pt_temp(1,iPart)
        Pt_temp(2,iPart) = PartState(5,iPart) - RKA(iStage) * Pt_temp(2,iPart)
        Pt_temp(3,iPart) = PartState(6,iPart) - RKA(iStage) * Pt_temp(3,iPart)
        Pt_temp(4,iPart) = Pt       (1,iPart) - RKA(iStage) * Pt_temp(4,iPart)
        Pt_temp(5,iPart) = Pt       (2,iPart) - RKA(iStage) * Pt_temp(5,iPart)
        Pt_temp(6,iPart) = Pt       (3,iPart) - RKA(iStage) * Pt_temp(6,iPart)

        ! Sanity Check Particle Pusher / WARNING: Might Cause Slowdowns
        !IF (ANY(ISNAN(Pt(:,iPart)))) THEN
        !  IPWRITE(UNIT_stdOut,*) 'Found invalid particle push, ignoring. PartID:', iPart
        !  Pt(:,iPart) = 0
        !ENDIF

        PartState(1,iPart) = PartState(1,iPart) + Pt_temp(1,iPart)*b_dt(iStage)
        PartState(2,iPart) = PartState(2,iPart) + Pt_temp(2,iPart)*b_dt(iStage)
        PartState(3,iPart) = PartState(3,iPart) + Pt_temp(3,iPart)*b_dt(iStage)
        PartState(4,iPart) = PartState(4,iPart) + Pt_temp(4,iPart)*b_dt(iStage)
        PartState(5,iPart) = PartState(5,iPart) + Pt_temp(5,iPart)*b_dt(iStage)
        PartState(6,iPart) = PartState(6,iPart) + Pt_temp(6,iPart)*b_dt(iStage)

        ! Sanity Check Particle / WARNING: Might Cause Slowdowns
        !IF (ANY(ISNAN(PartState(:,iPart)))) THEN
        !  PDM%ParticleInside(iPart) = .FALSE.
        !  IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        !ENDIF

      !IsNewPart: no Pt_temp history available. Either because of emissionType = 1 or because of reflection with almost zero wallVelo
      ELSE
        RandVal         = 1.
        Pa_rebuilt(:,:) = 0.
        DO iStage_loc=1,iStage
          Pa_rebuilt(1:3,iStage_loc)=Pa_rebuilt_coeff(iStage_loc)*Pt(1:3,iPart)
        END DO
        v_rebuilt(:,:)=0.
        DO iStage_loc=iStage-1,0,-1
          IF (iStage_loc.EQ.iStage-1) THEN
            v_rebuilt(1:3,iStage_loc) = PartState(4:6,iPart) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
          ELSE
            v_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc+1) - b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
          END IF
        END DO
        Pv_rebuilt(:,:)=0.
        DO iStage_loc=1,iStage
          IF (iStage_loc.EQ.1) THEN
            Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,0)
          ELSE
            Pv_rebuilt(1:3,iStage_loc) = v_rebuilt(1:3,iStage_loc-1) - RKA(iStage_loc)*Pv_rebuilt(1:3,iStage_loc-1)
          END IF
        END DO

        ! Pt_temp is rebuilt, do particle push
        Pt_temp  (1:3,iPart) = Pv_rebuilt(1:3,iStage)
        Pt_temp  (4:6,iPart) = Pa_rebuilt(1:3,iStage)
        PartState(1  ,iPart) = PartState(1,iPart) + Pt_temp(1,iPart)*b_dt(iStage)*RandVal
        PartState(2  ,iPart) = PartState(2,iPart) + Pt_temp(2,iPart)*b_dt(iStage)*RandVal
        PartState(3  ,iPart) = PartState(3,iPart) + Pt_temp(3,iPart)*b_dt(iStage)*RandVal
        PartState(4  ,iPart) = PartState(4,iPart) + Pt_temp(4,iPart)*b_dt(iStage)*RandVal
        PartState(5  ,iPart) = PartState(5,iPart) + Pt_temp(5,iPart)*b_dt(iStage)*RandVal
        PartState(6  ,iPart) = PartState(6,iPart) + Pt_temp(6,iPart)*b_dt(iStage)*RandVal

        PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...

        ! Sanity Check Particle / WARNING: Might Cause Slowdowns
        !IF (ANY(ISNAN(PartState(:,iPart)))) THEN
        !  PDM%ParticleInside(iPart) = .FALSE.
        !  IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        !ENDIF

      END IF !IsNewPart

      ! Try to find particles with too high velocity
      !v_magnitude   = SQRT(DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart)))

      !IF ((Species(PartSpecies(iPart))%HighVeloThreshold.NE.0).AND.(v_magnitude.GT.Species(PartSpecies(iPart))%HighVeloThreshold))THEN
      !  part_err = .TRUE.
      !  IPWRITE(UNIT_stdOut,*) ' High velocity particle detected. Writing error state and removing particle ...'
      !  IPWRITE(UNIT_stdout,*) ' LastPos:',  PartState(1:3,iPart)
      !  IPWRITE(UNIT_stdout,*) ' Velocity:', PartState(4:6,iPart)
      !  PDM%ParticleInside(iPart) = .FALSE.
      !END IF
    END IF
  END DO

  !IF (part_err) THEN
  !    CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
  !                        FutureTime=tWriteData,isErrorFile=.TRUE.)
  !END IF

#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

  ! particle tracking
#if USE_MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    IF (TriaTracking) THEN
      CALL ParticleTriaTracking()
    ELSE
      CALL ParticleTracing()
    END IF
  END IF
  ! emitt particles inserted in current time step
  CALL ParticleInserting()
#if USE_MPI
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! find next free position in particle array
  CALL MPIParticleRecv()
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif
END IF

! This >MIGHT< be only needed at the last RK stage ...
!IF (iStage.EQ.nRKStages) THEN
CALL UpdateNextFreePosition()
!END IF

END SUBROUTINE Particle_TimeStepByLSERK_RK

END MODULE MOD_Particle_TimeDisc
