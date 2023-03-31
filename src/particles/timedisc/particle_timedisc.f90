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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE Particle_InitTimeDisc
  MODULE PROCEDURE Particle_InitTimeDisc
END INTERFACE

INTERFACE ParticleTimeRHS
  MODULE PROCEDURE ParticleTimeRHS
END INTERFACE

! > Dummy interface for time step function pointer
ABSTRACT INTERFACE
  SUBROUTINE ParticleTimeStepPointer(t,dt)
    REAL,INTENT(IN)          :: t
    REAL,INTENT(IN)          :: dt
  END SUBROUTINE
END INTERFACE

! > Dummy interface for time step function pointer
ABSTRACT INTERFACE
  SUBROUTINE ParticleTimeStepRKPointer(t,dt,CurrentStage)
    REAL,INTENT(IN)    :: t
    REAL,INTENT(IN)    :: dt
    INTEGER,INTENT(IN) :: CurrentStage
  END SUBROUTINE
END INTERFACE

INTERFACE Particle_FinalizeTimeDisk
  MODULE PROCEDURE Particle_FinalizeTimeDisk
END INTERFACE

PROCEDURE(ParticleTimeStepPointer),  POINTER :: ParticleTimeStep     !< Point to the particle time step routine to be used
PROCEDURE(ParticleTimeStepRKPointer),POINTER :: ParticleTimeStepRK   !< Point to the particle RK time step routine to be used

PUBLIC :: Particle_InitTimeDisc
PUBLIC :: ParticleTimeRHS
PUBLIC :: ParticleTimeStep
PUBLIC :: ParticleTimeStepRK
PUBLIC :: Particle_FinalizeTimeDisk

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Determine particle time stepping and set required pointers
!===================================================================================================================================
SUBROUTINE Particle_InitTimeDisc()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETSTR,GETREAL
USE MOD_TimeDisc_Vars             ,ONLY: TimeStep,TimeDiscType
USE MOD_Particle_TimeDisc_Vars    ,ONLY: ParticleTimeDiscMethod,UseManualTimeStep,ManualTimeStep
USE MOD_Particle_Vars             ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SELECT CASE(TimeDiscType)
  CASE('LSERKW2','LSERKK3')
    ! Do nothing
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'Particle tracking only supported with DG TimeDiscType=LSERKW')
END SELECT

ParticleTimeDiscMethod = GETSTR('ParticleTimeDiscMethod','Runge-Kutta')
! Check if we are running a steady state tracking
LBWRITE(UNIT_stdOut,'(66("-"))')
!--- Read Manual Time Step
useManualTimeStep = GETLOGICAL('Part-SteadyState'   )
ManualTimeStep    = GETREAL   ('Part-ManualTimeStep')
IF (useManualTimeStep .OR. ManualTimeStep.GT.0.) THEN
  useManualTimeStep = .TRUE.
  TimeStep => TimeStepSteadyState
END IF

! Select the time disc method
SELECT CASE (TRIM(ParticleTimeDiscMethod))
  CASE('Runge-Kutta')
    ParticleTimeStep   => Particle_TimeStepByLSERK
    ParticleTimeStepRK => Particle_TimeStepByLSERK_RK
  CASE('Euler')
    ParticleTimeStep   => Particle_TimeStepByEuler
    ParticleTimeStepRK => Particle_TimeStepByEuler_RK
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,&
                    'Unknown method of particle time discretization: '//TRIM(ParticleTimeDiscMethod))
END SELECT

! Switch to dummy routines if running particle code without any species
IF (nSpecies.LE.0) THEN
  ParticleTimeStep   => Particle_TimeStepDummy
  ParticleTimeStepRK => Particle_TimeStepDummy_RK
END IF

END SUBROUTINE Particle_InitTimeDisc


!===================================================================================================================================
!> Calculate the right hand side before updating the field solution. Can be used to hide sending of number of particles.
!===================================================================================================================================
SUBROUTINE ParticleTimeRHS(t,dt,iStage)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,                     ONLY: U
USE MOD_Part_RHS,                    ONLY: CalcPartRHS
USE MOD_Particle_Interpolation,      ONLY: InterpolateFieldToParticle
USE MOD_Particle_Interpolation_Vars, ONLY: FieldAtParticle
USE MOD_Particle_Vars,               ONLY: PartState,LastPartPos,PDM,PEM
USE MOD_Particle_SGS,                ONLY: ParticleSGS
USE MOD_Particle_SGS_Vars,           ONLY: SGSinUse
USE MOD_Particle_Surface_Flux,       ONLY: ParticleSurfaceflux
#if USE_MPI
USE MOD_Particle_MPI_Vars,           ONLY: PartMPIExchange
#endif /*MPI*/
#if USE_RW
USE MOD_DG_Vars,                     ONLY: UTurb
USE MOD_Equation_Vars,               ONLY: nVarTurb
USE MOD_Particle_Interpolation_Vars, ONLY: TurbFieldAtParticle
USE MOD_Particle_RandomWalk,         ONLY: ParticleRandomWalk
USE MOD_Restart_Vars,                ONLY: RestartTurb
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,          ONLY: LBStartTime,LBPauseTime
USE MOD_Particle_Localization,       ONLY: CountPartsPerElem
#endif
#if USE_EXTEND_RHS || USE_FAXEN_CORR
USE MOD_DG_Vars,                     ONLY: Ut,UPrim
USE MOD_Particle_Interpolation_Vars, ONLY: GradAtParticle
USE MOD_Lifting_Vars,                ONLY: gradUx,gradUy,gradUz
USE MOD_Part_RHS,                    ONLY: extRHS
USE MOD_Mesh_Vars,                   ONLY: nElems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: dt
INTEGER,INTENT(IN),OPTIONAL   :: iStage
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
#if USE_EXTEND_RHS || USE_FAXEN_CORR
REAL                          :: U_RHS(1:RHS_NVARS,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR  */
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
! #if USE_LOADBALANCE
! ! Needed for scaling parts and load balance
! CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.)
! #endif

PartMPIExchange%nMPIParticles=0
#endif /*USE_MPI*/

! set last particle position and element
LastPartPos(1:3,1:PDM%ParticleVecLength) = PartState(PART_POSV,1:PDM%ParticleVecLength)
PEM%lastElement(1:PDM%ParticleVecLength) = PEM%Element(        1:PDM%ParticleVecLength)

CALL ParticleSurfaceflux()
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_EXTEND_RHS || USE_FAXEN_CORR
! Calculate tau
CALL extRHS(UPrim,Ut,U_RHS)
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */
CALL InterpolateFieldToParticle(PP_nVar,U,PP_nVarPrim,FieldAtParticle&
#if USE_EXTEND_RHS || USE_FAXEN_CORR
  ,gradUx(RHS_LIFTVARS,:,:,:,:),gradUy(RHS_LIFTVARS,:,:,:,:),gradUz(RHS_LIFTVARS,:,:,:,:),U_RHS,GradAtParticle&
#endif /* USE_EXTEND_RHS || USE_FAXEN_CORR */
)
#if USE_RW
IF (RestartTurb) CALL InterpolateFieldToParticle(nVarTurb,UTurb,nVarTurb,TurbFieldAtParticle)
!--> Calculate the random walk push
CALL ParticleRandomWalk(t)
#endif
!--> Calculate the particle right hand side and push
CALL CalcPartRHS(&
#if USE_BASSETFORCE || ANALYZE_RHS
t,dt,iStage)
#else
)
#endif /* USE_BASSETFORCE || ANALYZE_RHS */
IF (SGSinUse) CALL ParticleSGS(dt,iStage)
#if USE_LOADBALANCE
CALL LBPauseTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/

! Suppress compiler warning
NO_OP(t)

END SUBROUTINE ParticleTimeRHS


!===================================================================================================================================
!> Euler particle time integration:
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepDummy(t,dt)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: dt
!===================================================================================================================================
NO_OP(t)
NO_OP(dt)
END SUBROUTINE Particle_TimeStepDummy


!===================================================================================================================================
!> Euler particle time integration:
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepDummy_RK(t,dt,iStage)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: dt
INTEGER,INTENT(IN)            :: iStage
!===================================================================================================================================
! Suppress compiler warning
NO_OP(t)
NO_OP(dt)
NO_OP(iStage)
END SUBROUTINE Particle_TimeStepDummy_RK


!===================================================================================================================================
!> Euler particle time integration:
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByEuler(t,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Particle_Analyze_Tools,  ONLY: TrackingParticlePath
USE MOD_Particle_Analyze_Vars,   ONLY: doParticleDispersionTrack,doParticlePathTrack
USE MOD_Particle_Tracking,       ONLY: PerformTracking
USE MOD_Particle_Vars,           ONLY: Species,PartSpecies,PartState,Pt,PDM
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles,MPIParticleSend,MPIParticleRecv,SendNbOfParticles
#endif /*MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBSplitTime
#endif
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
!===================================================================================================================================

CALL ParticleTimeRHS(t,dt)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! particle push using Euler
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    !-- Tracer Particles
    IF (Species(PartSpecies(iPart))%RHSMethod.EQ.RHS_TRACER) THEN
      PartState(PART_POSV,iPart) = PartState(PART_POSV,iPart) + PartState(PART_VELV,iPart)*dt
      PartState(PART_VELV,iPart) = Pt       (1:3      ,iPart)
    !-- Normal particles
    ELSE
      PartState(PART_POSV,iPart)                  = PartState(PART_POSV                 ,iPart) + PartState(PART_VELV,iPart)*dt
      PartState(PART_VEL1:3+PP_nVarPartRHS,iPart) = PartState(PART_VEL1:3+PP_nVarPartRHS,iPart) + Pt       (:        ,iPart)*dt
    ENDIF !< Tracer
  ENDIF !< ParticleInside
END DO

! No BC interaction expected, so path can be calculated here. Periodic BCs are ignored purposefully
IF (doParticleDispersionTrack.OR.doParticlePathTrack) CALL TrackingParticlePath

#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

! Locate and communicate particles (only if particle have changed)
#if USE_MPI
! open receive buffer for number of particles
CALL IRecvNbofParticles()
#if USE_LOADBALANCE
CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif
! track new particle position
CALL PerformTracking()
#if USE_LOADBALANCE
CALL LBSplitTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
! emitt particles inserted in current time step
CALL ParticleInserting()

! Suppress compiler warning
NO_OP(t)

END SUBROUTINE Particle_TimeStepByEuler


!===================================================================================================================================
!> Euler particle time integration:
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByEuler_RK(t,dt,iStage)
! MODULES
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: dt
INTEGER,INTENT(IN)            :: iStage
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Suppress compiler warning
NO_OP(t)
NO_OP(dt)
NO_OP(iStage)
! Necessary to avoid further if statements in dg operator
#if USE_MPI
! open receive buffer for number of particles
CALL IRecvNbofParticles()
#endif
END SUBROUTINE Particle_TimeStepByEuler_RK


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK(t,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Particle_Analyze_Tools,  ONLY: ParticleRecord,ParticleRecordPath,TrackingParticlePath
USE MOD_Particle_Analyze_Vars,   ONLY: doParticleDispersionTrack,doParticlePathTrack,RecordPart
USE MOD_Particle_Tracking,       ONLY: PerformTracking
USE MOD_Particle_Vars,           ONLY: PartState,Pt,Pt_temp,PDM,PartSpecies,Species
USE MOD_TimeDisc_Vars,           ONLY: b_dt
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles
#endif /*MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: dt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart
!REAL                          :: v_magnitude
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

CALL ParticleTimeRHS(t,dt,1)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! particle push for first RK stage
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! Pt is always known at this position, change isNewPart to false
    PDM%IsNewPart(iPart) = .FALSE.

    IF (Species(PartSpecies(iPart))%RHSMethod.EQ.RHS_TRACER)THEN
      Pt_temp  (1:3      ,iPart) = PartState(PART_VELV,iPart)
      PartState(PART_VELV,iPart) = Pt(       1:3      ,iPart)

      PartState(PART_POSV,iPart) = PartState(PART_POSV,iPart) + PartState(PART_VELV,iPart)*b_dt(1)
    ELSE
      Pt_temp  (1:3,iPart)             = PartState(PART_VELV       ,iPart)
      Pt_temp  (4:3+PP_nVarPartRHS,iPart) = Pt       (1:PP_nVarPartRHS,iPart)

      PartState(PART_POSV,iPart)                  = PartState(PART_POSV                 ,iPart) + PartState(PART_VELV,iPart)*b_dt(1)
      PartState(PART_VEL1:3+PP_nVarPartRHS,iPart) = PartState(PART_VEL1:3+PP_nVarPartRHS,iPart) + Pt       (:        ,iPart)*b_dt(1)
    END IF
  END IF
END DO

! No BC interaction expected, so path can be calculated here. Periodic BCs are ignored purposefully
IF (doParticleDispersionTrack.OR.doParticlePathTrack) CALL TrackingParticlePath

IF (RecordPart.GT.0) CALL ParticleRecordPath()

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
CALL PerformTracking()
#if USE_LOADBALANCE
CALL LBSplitTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
! emitt particles inserted in current time step
CALL ParticleInserting()

IF (RecordPart.GT.0) CALL ParticleRecord(t)

! Suppress compiler warning
NO_OP(dt)

END SUBROUTINE Particle_TimeStepByLSERK


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version. Inner RK iteration
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK_RK(t,dt,iStage)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_Part_Emission,           ONLY: ParticleInserting
USE MOD_Particle_Analyze_Tools,  ONLY: ParticleRecord,ParticleRecordPath,TrackingParticlePath
USE MOD_Particle_Analyze_Vars,   ONLY: doParticleDispersionTrack,doParticlePathTrack,RecordPart
USE MOD_Particle_TimeDisc_Vars,  ONLY: Pa_rebuilt,Pa_rebuilt_coeff,Pv_rebuilt,v_rebuilt
USE MOD_Particle_Tracking,       ONLY: PerformTracking
USE MOD_Particle_Vars,           ONLY: PartState,Pt,Pt_temp,PDM,PartSpecies,Species
USE MOD_TimeDisc_Vars,           ONLY: RKA,b_dt
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles
#endif /*MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,      ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
REAL,INTENT(IN)               :: dt
INTEGER,INTENT(IN)            :: iStage
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,iStage_loc
REAL,PARAMETER                :: RandVal = 1.                         ! Random time increment for new parts to accomplish temporal
                                                                      ! dispersion within on time step. Currently disabled
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

CALL ParticleTimeRHS(t,dt,iStage)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! particle push for nth RK stage
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! "normal" particles are pushed with whole timestep
    IF (.NOT.PDM%IsNewPart(iPart)) THEN
      IF (Species(PartSpecies(iPart))%RHSMethod.EQ.RHS_TRACER) THEN
        Pt_temp(1:3,iPart) = PartState(PART_VELV,iPart) - RKA(iStage) * Pt_temp(1:3,iPart)

        PartState(PART_POSV,iPart) = PartState(PART_POSV,iPart) + Pt_temp(1:3,iPart)*b_dt(iStage)
        PartState(PART_VELV,iPart) = Pt(       1:3      ,iPart)
      ELSE
        Pt_temp(1:3,iPart)                  = PartState(PART_VELV       ,iPart) - RKA(iStage) * Pt_temp(1:3               ,iPart)
        Pt_temp(4:3+PP_nVarPartRHS,iPart)   = Pt       (1:PP_nVarPartRHS,iPart) - RKA(iStage) * Pt_temp(4:3+PP_nVarPartRHS,iPart)

        PartState(1:3+PP_nVarPartRHS,iPart) = PartState(1:3+PP_nVarPartRHS,iPart) + Pt_temp(:,iPart)*b_dt(iStage)
      END IF

    !IsNewPart: no Pt_temp history available. Either because of emissionType = 1 or because of reflection with almost zero wallVelo
    ELSE
      Pa_rebuilt(:,:) = 0.
      DO iStage_loc=1,iStage
        Pa_rebuilt(1:PP_nVarPartRHS,iStage_loc) = Pa_rebuilt_coeff(iStage_loc)*Pt(1:PP_nVarPartRHS,iPart)
      END DO
      v_rebuilt(:,:)=0.
      DO iStage_loc=iStage-1,0,-1
        IF (iStage_loc.EQ.iStage-1) THEN
          v_rebuilt(1:3,iStage_loc) = PartState(PART_VELV,iPart) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
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
      Pt_temp  (1:3,iPart)                = Pv_rebuilt(1:3,iStage)
      Pt_temp  (4:3+PP_nVarPartRHS,iPart) = Pa_rebuilt(1:PP_nVarPartRHS,iStage)
      PartState(1:3+PP_nVarPartRHS,iPart) = PartState( 1:3+PP_nVarPartRHS,iPart) + Pt_temp(:,iPart)*b_dt(iStage)*RandVal

      PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...
    END IF !IsNewPart
  END IF
END DO

! No BC interaction expected, so path can be calculated here. Periodic BCs are ignored purposefully
IF (doParticleDispersionTrack.OR.doParticlePathTrack) CALL TrackingParticlePath()

IF (RecordPart.GT.0) CALL ParticleRecordPath()

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
CALL PerformTracking()
#if USE_LOADBALANCE
CALL LBSplitTime(LB_TRACK,tLBStart)
#endif /*USE_LOADBALANCE*/
! emitt particles inserted in current time step
CALL ParticleInserting()

IF (RecordPart.GT.0) CALL ParticleRecord(t)

END SUBROUTINE Particle_TimeStepByLSERK_RK


#if USE_PARTICLES
!===================================================================================================================================
!> Low-Storage Runge-Kutta integration with frozen fluid: 2 register version
!> This procedure takes the current time t, the time step dt and the intial solution.
!> Only particle push is integrated in time
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE TimeStepSteadyState(t)
! MODULES
USE MOD_Globals               ,ONLY: CollectiveStop
USE MOD_PreProc
USE MOD_TimeDisc_Vars         ,ONLY: dt,b_dt,RKb,RKc,nRKStages,CurrentStage
USE MOD_Part_Tools            ,ONLY: UpdateNextFreePosition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers    ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
USE MOD_LoadBalance_Vars      ,ONLY: PerformLBSample
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_Particle_MPI          ,ONLY: IRecvNbOfParticles,MPIParticleSend,MPIParticleRecv,SendNbOfParticles
#endif /* USE_MPI */
!#if USE_MPI_SHARED
!USE MOD_Particle_MPI_Shared,ONLY:UpdateDGShared
!#endif /* MPI_SHARED */
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: t                                     !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tStage
INTEGER                         :: iStage
#if USE_LOADBALANCE
REAL                            :: tLBStart,tDGStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! Premultiply with dt
b_dt = RKb*dt

#if USE_LOADBALANCE
! Add a minimal compute load to DG elements
IF (PerformLBSample) THEN
  CALL LBStartTime(tDGStart)
END IF
#endif /*USE_LOADBALANCE*/

! First evaluation of DG operator already done in timedisc
CurrentStage = 1
tStage       = t

CALL ParticleTimeStep(tStage,dt)

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
#endif /*USE_MPI*/

! find next free position in particle array
CALL UpdateNextFreePosition()

! Following steps
DO iStage = 2,nRKStages
  CurrentStage = iStage
  tStage       = t+dt*RKc(iStage)

  CALL ParticleTimeStepRK(tStage,dt,iStage)

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
#endif /*USE_MPI*/

  ! find next free position in particle array
  CALL UpdateNextFreePosition()
END DO

CurrentStage = 1

#if USE_LOADBALANCE
! Add a minimal compute load to DG elements
IF (PerformLBSample) THEN
  CALL LBStartTime(tLBStart)
  ! 20% is empirical value which seems as sufficient DG load for load balancing
  tLBStart = tLBStart - 0.2*(tLBStart - tDGStart)
  CALL LBSplitTime(LB_DG,tLBStart)
END IF
#endif /*USE_LOADBALANCE*/

END SUBROUTINE TimeStepSteadyState
#endif

!===================================================================================================================================
!> Finalize particle time stepping and free variables
!===================================================================================================================================
SUBROUTINE Particle_FinalizeTimeDisk
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE Particle_FinalizeTimeDisk

END MODULE MOD_Particle_TimeDisc
