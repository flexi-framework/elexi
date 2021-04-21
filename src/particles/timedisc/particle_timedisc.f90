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
#include "eos.h"

!==================================================================================================================================
!> Module for the particle temporal discretization
!==================================================================================================================================
MODULE MOD_Particle_TimeDisc
! MODULES
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
SUBROUTINE Particle_TimeStepByEuler(t,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                 ONLY: U
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Part_RHS,                ONLY: CalcPartRHS
USE MOD_Part_Tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze,        ONLY: TrackingParticlePath
USE MOD_Particle_Analyze_Vars,   ONLY: doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Interpolation,  ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars,  ONLY: FieldAtParticle
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_vars,  ONLY: TrackingMethod
USE MOD_Particle_Vars,           ONLY: Species, PartSpecies, PartState, Pt, LastPartPos, DelayTime, PEM, PDM
USE MOD_Particle_SGS,            ONLY: ParticleSGS
USE MOD_Particle_SGS_Vars,       ONLY: SGSinUse
USE MOD_Particle_Surface_Flux,   ONLY: ParticleSurfaceflux
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
#if USE_RW
USE MOD_DG_Vars,                 ONLY: UTurb
USE MOD_Equation_Vars,           ONLY: nVarTurb
USE MOD_Particle_Interpolation_Vars,  ONLY: TurbFieldAtParticle
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
USE MOD_Restart_Vars,            ONLY: RestartTurb
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBPauseTime,LBSplitTime
USE MOD_Particle_Localization,   ONLY: CountPartsPerElem
#endif
#if USE_EXTEND_RHS
USE MOD_Particle_Interpolation_Vars,ONLY: GradAtParticle
USE MOD_Lifting_Vars,            ONLY: gradUx,gradUy,gradUz
USE MOD_Part_RHS,                ONLY: tauRHS
USE MOD_Mesh_Vars,               ONLY: nElems
#endif /* USE_EXTEND_RHS */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: dt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
#if USE_EXTEND_RHS
REAL                          :: divtau(1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL                          :: gradp( 1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
#endif /* USE_EXTEND_RHS */
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
LastPartPos(1:3,1:PDM%ParticleVecLength) = PartState(1:3,1:PDM%ParticleVecLength)
PEM%lastElement(1:PDM%ParticleVecLength) = PEM%Element(  1:PDM%ParticleVecLength)

! forces on particles
IF (t.GE.DelayTime) THEN
  CALL ParticleSurfaceflux()
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if  USE_EXTEND_RHS
  ! Calculate tau
  CALL tauRHS(U,divtau,gradp)
#endif
  CALL InterpolateFieldToParticle(PP_nVar,U     ,PP_nVarPrim,FieldAtParticle&
#if  USE_EXTEND_RHS
    ,gradUx(RHS_LIFTVARS,:,:,:,:),gradUy(RHS_LIFTVARS,:,:,:,:),gradUz(RHS_LIFTVARS,:,:,:,:),divtau,gradp,GradAtParticle&
#endif
    )
#if USE_RW
  IF (RestartTurb) CALL InterpolateFieldToParticle(nVarTurb,UTurb,nVarTurb,TurbFieldAtParticle)
  CALL ParticleRandomWalk(t)
#endif
  CALL CalcPartRHS(&
#if USE_BASSETFORCE
    dt)
#else
    )
#endif /* USE_BASSETFORCE */
  IF (SGSinUse) CALL ParticleSGS(1,dt)
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF

! particle push using Euler
IF (t.GE.DelayTime) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      !-- Tracer Particles
      IF (Species(PartSpecies(iPart))%RHSMethod.EQ.RHS_TRACER) THEN
        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart)*dt
        PartState(4:6,iPart) = Pt       (1:3,iPart)
      !-- Normal particles
      ELSE
        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart)*dt
        PartState(4:6,iPart) = PartState(4:6,iPart) + Pt       (1:3,iPart)*dt
      ENDIF !< Tracer
    ENDIF !< ParticleInside
  END DO
END IF

! No BC interaction expected, so path can be calculated here. Periodic BCs are ignored purposefully
IF (doParticleDispersionTrack.OR.doParticlePathTrack) CALL TrackingParticlePath

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
  SELECT CASE (TrackingMethod)
    CASE(REFMAPPING)
      CALL ParticleRefTracking()
    CASE(TRACING)
      CALL ParticleTracing()
    CASE(TRIATRACKING)
      CALL ParticleTriaTracking()
  END SELECT
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
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
SUBROUTINE Particle_TimeStepByLSERK_RHS(t,iStage,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Part_RHS,                ONLY: CalcPartRHS
USE MOD_Particle_Interpolation,  ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars,  ONLY: FieldAtParticle
USE MOD_Particle_Vars,           ONLY: PartState,DelayTime,LastPartPos,PDM,PEM
USE MOD_Particle_SGS,            ONLY: ParticleSGS
USE MOD_Particle_SGS_Vars,       ONLY: SGSinUse
USE MOD_Particle_Surface_Flux,   ONLY: ParticleSurfaceflux
! USE MOD_Timedisc_Vars,           ONLY: nRKStages
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles,MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
#if USE_RW
USE MOD_DG_Vars,                 ONLY: UTurb
USE MOD_Equation_Vars,           ONLY: nVarTurb
USE MOD_Particle_Interpolation_Vars,  ONLY: TurbFieldAtParticle
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
USE MOD_Restart_Vars,            ONLY: RestartTurb
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBPauseTime
USE MOD_Particle_Localization,   ONLY: CountPartsPerElem
#endif
#if USE_EXTEND_RHS
USE MOD_Particle_Interpolation_Vars,ONLY: GradAtParticle
USE MOD_Lifting_Vars,            ONLY: gradUx,gradUy,gradUz
USE MOD_Part_RHS,                ONLY: tauRHS
USE MOD_Mesh_Vars,               ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
INTEGER,INTENT(IN)            :: iStage
REAL,INTENT(IN)               :: dt
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
#if USE_EXTEND_RHS
REAL                          :: divtau(1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL                          :: gradp( 1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
#endif /* USE_EXTEND_RHS */
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
LastPartPos(1:3,1:PDM%ParticleVecLength) = PartState(1:3,1:PDM%ParticleVecLength)
PEM%lastElement(1:PDM%ParticleVecLength) = PEM%Element(  1:PDM%ParticleVecLength)

IF (t.GE.DelayTime) THEN
  CALL ParticleSurfaceflux()
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_EXTEND_RHS
  ! Calculate tau
  CALL tauRHS(U,divtau,gradp)
#endif /* USE_EXTEND_RHS */
  CALL InterpolateFieldToParticle(PP_nVar,U,PP_nVarPrim,FieldAtParticle&
#if USE_EXTEND_RHS
    ,gradUx(RHS_LIFTVARS,:,:,:,:),gradUy(RHS_LIFTVARS,:,:,:,:),gradUz(RHS_LIFTVARS,:,:,:,:),divtau,gradp,GradAtParticle&
#endif
  )
#if USE_RW
  IF (RestartTurb) CALL InterpolateFieldToParticle(nVarTurb,UTurb,nVarTurb,TurbFieldAtParticle)
  !--> Calculate the random walk push
  CALL ParticleRandomWalk(t)
#endif
  !--> Calculate the particle right hand side and push
  CALL CalcPartRHS(&
#if USE_BASSETFORCE
  dt,iStage)
#else
  )
#endif /* USE_BASSETFORCE */
  IF (SGSinUse) CALL ParticleSGS(iStage,dt)
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
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_PruettDamping,           ONLY: TempFilterTimeDeriv
USE MOD_TimeDisc_Vars,           ONLY: nRKStages
#if FV_ENABLED
USE MOD_FV,                      ONLY: FV_Switch
#endif
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Part_Tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze,        ONLY: TrackingParticlePath
USE MOD_Particle_Analyze_Vars,   ONLY: doParticleDispersionTrack,doParticlePathTrack,RecordPart
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_vars,  ONLY: TrackingMethod
USE MOD_Particle_Vars,           ONLY: PartState,Pt,Pt_temp,DelayTime,PDM,PartSpecies,Species
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
USE MOD_Particle_Analyze_Tools,  ONLY: ParticleRecord
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: b_dt(1:nRKStages)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
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

  ! particle push for first RK stage
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      ! Pt is always known at this position, change isNewPart to false
      PDM%IsNewPart(iPart) = .FALSE.

      IF (Species(PartSpecies(iPart))%RHSMethod.EQ.RHS_TRACER)THEN
        Pt_temp  (1:3,iPart) = PartState(4:6,iPart)
        PartState(4:6,iPart) = Pt(       1:3,iPart)

        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart)*b_dt(1)
      ELSE
        Pt_temp  (1:3,iPart) = PartState(4:6,iPart)
        Pt_temp  (4:6,iPart) = Pt       (1:3,iPart)

        PartState(1:3,iPart) = PartState(1:3,iPart) + PartState(4:6,iPart)*b_dt(1)
        PartState(4:6,iPart) = PartState(4:6,iPart) + Pt       (1:3,iPart)*b_dt(1)
      END IF
    END IF
  END DO

  ! No BC interaction expected, so path can be calculated here. Periodic BCs are ignored purposefully
  IF (doParticleDispersionTrack.OR.doParticlePathTrack) CALL TrackingParticlePath

  IF (RecordPart.GT.0) CALL ParticleRecord(t)

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
  SELECT CASE (TrackingMethod)
    CASE(REFMAPPING)
      CALL ParticleRefTracking()
    CASE(TRACING)
      CALL ParticleTracing()
    CASE(TRIATRACKING)
      CALL ParticleTriaTracking()
  END SELECT
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
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

CALL UpdateNextFreePosition()

END SUBROUTINE Particle_TimeStepByLSERK


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> Calculate the right hand side before updating the field solution. Can be used to hide sending of number of particles.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK_RK_RHS(t,iStage,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Particle_Interpolation,  ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars,  ONLY: FieldAtParticle
USE MOD_Part_RHS,                ONLY: CalcPartRHS
USE MOD_Particle_Vars,           ONLY: PartState,DelayTime,LastPartPos,PDM,PEM
USE MOD_Particle_SGS,            ONLY: ParticleSGS
USE MOD_Particle_SGS_Vars,       ONLY: SGSinUse
USE MOD_Particle_Surface_Flux,   ONLY: ParticleSurfaceflux
! USE MOD_Timedisc_Vars,           ONLY: nRKStages
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
#if USE_RW
USE MOD_DG_Vars,                 ONLY: UTurb
USE MOD_Equation_Vars,           ONLY: nVarTurb
USE MOD_Particle_Interpolation_Vars,  ONLY: TurbFieldAtParticle
USE MOD_Particle_RandomWalk,     ONLY: ParticleRandomWalk
USE MOD_Restart_Vars,            ONLY: RestartTurb
#endif
#if USE_LOADBALANCE
USE MOD_Particle_Localization,   ONLY: CountPartsPerElem
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
#if USE_EXTEND_RHS
USE MOD_Particle_Interpolation_Vars,ONLY: GradAtParticle
USE MOD_Lifting_Vars,            ONLY: gradUx,gradUy,gradUz
USE MOD_Part_RHS,                ONLY: tauRHS
USE MOD_Mesh_Vars,               ONLY: nElems
#endif /* USE_EXTEND_RHS */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                   :: t
INTEGER,INTENT(IN)                :: iStage
REAL,INTENT(IN)                   :: dt
#if USE_LOADBALANCE
REAL                              :: tLBStart
#endif /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_EXTEND_RHS
REAL                              :: divtau(1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL                              :: gradp( 1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
#endif /* USE_EXTEND_RHS */
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
LastPartPos(1:3,1:PDM%ParticleVecLength) = PartState(1:3,1:PDM%ParticleVecLength)
PEM%lastElement(1:PDM%ParticleVecLength) = PEM%Element(  1:PDM%ParticleVecLength)

IF (t.GE.DelayTime) THEN
  CALL ParticleSurfaceflux()
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_EXTEND_RHS
  ! Calculate tau
  CALL tauRHS(U,divtau,gradp)
#endif /* USE_EXTEND_RHS */
  CALL InterpolateFieldToParticle(PP_nVar,U,PP_nVarPrim,FieldAtParticle&
#if USE_EXTEND_RHS
    ,gradUx(RHS_LIFTVARS,:,:,:,:),gradUy(RHS_LIFTVARS,:,:,:,:),gradUz(RHS_LIFTVARS,:,:,:,:),divtau,gradp,GradAtParticle&
#endif /* USE_EXTEND_RHS */
    )
#if USE_RW
  IF (RestartTurb) CALL InterpolateFieldToParticle(nVarTurb,UTurb,nVarTurb,TurbFieldAtParticle)
  CALL ParticleRandomWalk(t)
#endif
  !--> Calculate the particle right hand side and push
  CALL CalcPartRHS(&
#if USE_BASSETFORCE
  dt,iStage)
#else
  )
#endif /* USE_BASSETFORCE */
  IF (SGSinUse) CALL ParticleSGS(iStage,dt)
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
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Part_Tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze,        ONLY: TrackingParticlePath
USE MOD_Particle_Analyze_Vars,   ONLY: doParticleDispersionTrack,doParticlePathTrack,RecordPart
USE MOD_Particle_TimeDisc_Vars,  ONLY: Pa_rebuilt,Pa_rebuilt_coeff,Pv_rebuilt,v_rebuilt
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_Vars,  ONLY: TrackingMethod
USE MOD_Particle_Vars,           ONLY: PartState,Pt,Pt_temp,DelayTime,PDM,PartSpecies,Species
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
USE MOD_Particle_Analyze_Tools,  ONLY: ParticleRecord
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
INTEGER,INTENT(IN)            :: iStage
REAL,INTENT(IN)               :: b_dt(1:nRKStages)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,iStage_loc
REAL,PARAMETER                :: RandVal = 1.                         ! Random time increment for new parts to accomplish temporal
                                                                      ! dispersion within on time step. Currently disabled
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

IF (t.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

  ! particle push for nth RK stage
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      ! "normal" particles are pushed with whole timestep
      IF (.NOT.PDM%IsNewPart(iPart)) THEN
        IF (Species(PartSpecies(iPart))%RHSMethod.EQ.RHS_TRACER) THEN
          Pt_temp(1:3,iPart) = PartState(4:6,iPart) - RKA(iStage) * Pt_temp(1:3,iPart)

          PartState(1:3,iPart) = PartState(1:3,iPart) + Pt_temp(1:3,iPart)*b_dt(iStage)
          PartState(4:6,iPart) = Pt(       1:3,iPart)
        ELSE
          Pt_temp(1:3,iPart) = PartState(4:6,iPart) - RKA(iStage) * Pt_temp(1:3,iPart)
          Pt_temp(4:6,iPart) = Pt       (1:3,iPart) - RKA(iStage) * Pt_temp(4:6,iPart)

          PartState(1:6,iPart) = PartState(1:6,iPart) + Pt_temp(1:6,iPart)*b_dt(iStage)
        END IF

      !IsNewPart: no Pt_temp history available. Either because of emissionType = 1 or because of reflection with almost zero wallVelo
      ELSE
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
        PartState(1:6,iPart) = PartState( 1:6,iPart) + Pt_temp(1:6,iPart)*b_dt(iStage)*RandVal

        PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...
      END IF !IsNewPart
    END IF
  END DO

  ! No BC interaction expected, so path can be calculated here. Periodic BCs are ignored purposefully
  IF (doParticleDispersionTrack.OR.doParticlePathTrack) CALL TrackingParticlePath()

  IF (RecordPart.GT.0) CALL ParticleRecord(t)

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
  SELECT CASE (TrackingMethod)
    CASE(REFMAPPING)
      CALL ParticleRefTracking()
    CASE(TRACING)
      CALL ParticleTracing()
    CASE(TRIATRACKING)
      CALL ParticleTriaTracking()
  END SELECT
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
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

CALL UpdateNextFreePosition()

END SUBROUTINE Particle_TimeStepByLSERK_RK

END MODULE MOD_Particle_TimeDisc
