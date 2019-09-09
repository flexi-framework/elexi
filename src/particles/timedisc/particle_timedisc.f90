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
!> Module for the GTS Temporal discretization
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

INTERFACE Particle_TimeDisc
  MODULE PROCEDURE Particle_TimeDisc
END INTERFACE

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

!INTERFACE Particle_FinalizeTimeDisc
!  MODULE PROCEDURE Particle_FinalizeTimeDisc
!END INTERFACE

PUBLIC::Particle_InitTimeDisc
PUBLIC::Particle_Timedisc
PUBLIC::Particle_TimeStepByEuler
PUBLIC::Particle_TimeStepByLSERK
PUBLIC::Particle_TimeStepByLSERK_RHS
PUBLIC::Particle_TimeStepByLSERK_RK
PUBLIC::Particle_TimeStepByLSERK_RK_RHS

!===================================================================================================================================

CONTAINS

SUBROUTINE Particle_InitTimeDisc()
!===================================================================================================================================
! Get information for end time and max time steps from ini file
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars
USE MOD_ReadInTools         ,ONLY:GETREAL,GETINT,GETSTR
USE MOD_StringTools         ,ONLY:LowCase,StripSpaces
USE MOD_Particle_Vars
#if USE_MPI
USE MOD_PICDepo_Vars        ,ONLY: DepositionType
USE MOD_Particle_MPI        ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

! compute ratio of dt for particle surface flux emission (example a). Additionally, the ratio of RK_inflow(iStage)/RK_c(iStage) 
! gives the ! maximum distance, the Surface Flux particles can during the initial implicit step (example b).
! example a)
!   particle number for stage 4: dt*RK_inflow(4) 
! or: generate all particles for dt, but only particles with random number R < RK_c(iStage) participate in current stage.
! This results in dt*RK_inflow(iStage) particles in the current stage, hence, we can decide if we generate all particles or
! only the particles per stage. Currently, all particles are generated prior to the RK stages.
! example b)
! ESDIRKO4 from kennedy and carpenter without an acting force. Assume again stage 4. The initial particles during stage 2 are 
! moved in stage 3 without tracking (because it could be dropped out of the domain and the negative increment of the time level. 
! The new particles in this  stage could travel a distance up to RK_c(4)=SUM(ESDIRKA(4,:)). Now, the new particles are pushed 
! further into the domain than the particles of the second stage has moved. This would create a non-uniform particle distribution.
! this is prevented by reducing their maximum emission/initial distance by RK_inflow(4)/RK_c(4).
! Note: A small overlap is possible, but this is required. See the charts in the docu folder.

! init
#if USE_MPI
IF ((TRIM(DepositionType).EQ."shape_function")             &
.OR.(TRIM(DepositionType).EQ."shape_function_1d")          &
.OR.(TRIM(DepositionType).EQ."shape_function_spherical")   &
.OR.(TRIM(DepositionType).EQ."shape_function_simple")      &
.OR.(TRIM(DepositionType).EQ."shape_function_cylindrical"))THEN
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
END IF
#endif /*MPI PARTICLES*/

END SUBROUTINE Particle_InitTimeDisc

!===================================================================================================================================
! GTS Temporal discretization 
!===================================================================================================================================
SUBROUTINE Particle_TimeDisc(iter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Mesh_Vars
USE MOD_TimeDisc_Vars       ,ONLY: dt,maxIter
USE MOD_Analyze             ,ONLY: Analyze
USE MOD_TestCase            ,ONLY: AnalyzeTestCase,CalcForcing
USE MOD_CalcTimeStep        ,ONLY: CalcTimeStep
USE MOD_Output              ,ONLY: Visualize,PrintStatusLine
USE MOD_HDF5_Output         ,ONLY: WriteState,WriteBaseFlow
USE MOD_DG                  ,ONLY: DGTimeDerivative_weakForm
USE MOD_Overintegration     ,ONLY: Overintegration
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
USE MOD_RecordPoints        ,ONLY: RecordPoints,WriteRP
#if FV_ENABLED
USE MOD_FV
#endif
use MOD_IO_HDF5
USE MOD_Particle_Mesh       ,ONLY: CountPartsPerElem
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_Particle_MPI        ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
USE MOD_Particle_Output     ,ONLY: Visualize_Particles
USE MOD_Particle_Tracking_vars
USE MOD_Particle_Vars
USE MOD_ReadInTools
USE MOD_Particle_HDF5_output,ONLY: WriteParticleToHDF5

IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)   :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: dt_Min
INTEGER                      :: iPart
LOGICAL                      :: NoPartInside
INTEGER                      :: nLostPartsTot
REAL                         :: vMax,vMaxx,vMaxy,vMaxz
!===================================================================================================================================

IF (iter.LE.maxIter) THEN
    dt_max_particles = dt ! initial evolution of field with maxwellts
ELSE
  NoPartInside=.TRUE.
  DO 
    vMaxx = 0.
    vMaxy = 0.
    vMaxz = 0.
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        vMaxx = MAX( vMaxx , ABS(PartState(iPart, 4) + dt*Pt(iPart,1)) )
        vMaxy = MAX( vMaxy , ABS(PartState(iPart, 5) + dt*Pt(iPart,2)) )
        vMaxz = MAX( vMaxz , ABS(PartState(iPart, 6) + dt*Pt(iPart,3)) )
        NoPartInside=.FALSE. 
      END IF
    END DO
vMax = MAX(vMaxx,vMaxy,vMaxz,1.0) 

#if USE_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,vMax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,NoPartInside,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,iError)
#endif /*MPI*/
    IF (NoPartInside) THEN
      dt_max_particles = dt
      EXIT
    ELSE
      dt_max_particles =  dt
    END IF
    dt = (dt_max_particles+dt)/2
    IF((dt.GE.dt_max_particles*0.95).AND.(dt.LE.dt_max_particles*1.05)) EXIT
  END DO
END IF

dt_Min = dt_max_particles

IF(CountNbOfLostParts)THEN
#if USE_MPI
    IF(MPIRoot) THEN
      CALL MPI_REDUCE(nLostParts,nLostPartsTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    ELSE ! no Root
      CALL MPI_REDUCE(nLostParts,nLostPartsTot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    END IF
#else
    nLostPartsTot=nLostParts
#endif /*MPI*/
END IF

IF(CountNbOfLostParts)THEN
    WRITE(UNIT_stdOut,'(A,I12)')' NbOfLostParticle : ',nLostPartsTot
END IF

END SUBROUTINE Particle_TimeDisc

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
USE MOD_Particle_MPI_Vars,       ONLY: DoExternalParts
USE MOD_Particle_Mesh,           ONLY: CountPartsPerElem
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
USE MOD_Particle_MPI_Vars,       ONLY: ExtPartState,ExtPartSpecies,ExtPartToFIBGM
#endif /*MPI*/
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_PICInterpolation
USE MOD_PICDepo
USE MOD_Part_tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping,TriaTracking
USE MOD_Particle_Vars,           ONLY: Species, PartSpecies, PartState, Pt, Pt_temp, LastPartPos, DelayTime, PEM, PDM
#if EQNSYSNR == 4
USE MOD_Particle_RandomWalk,     ONLY: Particle_RandomWalk
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: dt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart
!===================================================================================================================================
#if USE_MPI
CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB
#endif /*MPI*/

IF (t.GE.DelayTime) THEN
  ! communicate shape function particles
#if USE_MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! send number of particles
    CALL SendNbOfParticles()
  END IF
#endif /*MPI*/
CALL Deposition(doInnerParts=.TRUE.) ! because of emmision and UpdateParticlePosition
#if USE_MPI
  IF(DoExternalParts)THEN
    ! finish communication
    CALL MPIParticleRecv()
  END IF
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
END IF
!#endif /*PARTICLES*/

! set last data already here, since surfaceflux moved before interpolation
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)   ! forces on particles
  CALL CalcPartRHS()
#if EQNSYSNR == 4
  CALL Particle_RandomWalk()
#endif
END IF

IF (t.GE.DelayTime) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      !-- Tracer Particles
      IF (TRIM(Species(PartSpecies(iPart))%RHSMethod).EQ.'Tracer') THEN
        PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4)*dt
        PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5)*dt
        PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6)*dt
        PartState(iPart,4) = Pt(iPart,1)
        PartState(iPart,5) = Pt(iPart,2)
        PartState(iPart,6) = Pt(iPart,3)
      ELSE ! Normal particles
      !-- Particle Push
        ! Sanity Check Particle Pusher / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(Pt(iPart,:)))) THEN
            IPWRITE(UNIT_stdOut,*) 'Found invalid particle push, ignoring. PartID:', iPart
            Pt(iPart,:) = 0
        ENDIF
          
        Pt_temp(iPart,1) = PartState(iPart,4)
        Pt_temp(iPart,2) = PartState(iPart,5)
        Pt_temp(iPart,3) = PartState(iPart,6)
        Pt_temp(iPart,4) = Pt(iPart,1)
        Pt_temp(iPart,5) = Pt(iPart,2)
        Pt_temp(iPart,6) = Pt(iPart,3)
        PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4)*dt
        PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5)*dt
        PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6)*dt
        PartState(iPart,4) = PartState(iPart,4) + Pt(iPart,1)*dt
        PartState(iPart,5) = PartState(iPart,5) + Pt(iPart,2)*dt
        PartState(iPart,6) = PartState(iPart,6) + Pt(iPart,3)*dt
        
        ! Sanity Check Particle / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(PartState(iPart,:)))) THEN
            PDM%ParticleInside(iPart) = .FALSE.
            IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        ENDIF
        
      ENDIF !< Tracer
    ENDIF
  END DO
END IF

! Communicate particles
IF (t.GE.DelayTime) THEN ! removed .OR.(iter.EQ.0) because particles have not moved
#if USE_MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
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
  CALL ParticleInserting()
#if USE_MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()   ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()   ! finish communication
#endif
  CALL UpdateNextFreePosition()
END IF
  
END SUBROUTINE Particle_TimeStepByEuler

!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> Calculate the right hand side before updating the field solution. Can be used to hide sending of number of particles.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK_RHS(t,iStage,b_dt)
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars,           ONLY: nRKStages
#if USE_MPI
USE MOD_Particle_MPI_Vars,       ONLY: DoExternalParts
USE MOD_Particle_Mesh,           ONLY: CountPartsPerElem
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles,MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
USE MOD_Particle_MPI_Vars,       ONLY: ExtPartState,ExtPartSpecies,ExtPartToFIBGM
#endif /*MPI*/
USE MOD_PICInterpolation
USE MOD_PICDepo
USE MOD_Part_RHS,                ONLY: CalcPartRHS
#if EQNSYSNR == 4
USE MOD_Particle_RandomWalk,     ONLY: Particle_RandomWalk
#endif
USE MOD_Particle_Vars,           ONLY: PartState,DelayTime,LastPartPos,PDM,PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
INTEGER,INTENT(INOUT)         :: iStage
REAL,INTENT(IN)               :: b_dt(1:nRKStages)
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
CALL CountPartsPerElem(ResetNumberOfParticles=.TRUE.) !for scaling of tParts of LB
#endif /*MPI*/

IF (t.GE.DelayTime) THEN
  ! communicate shape function particles
#if USE_MPI
  PartMPIExchange%nMPIParticles=0
  IF(DoExternalParts)THEN
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! send number of particles
    CALL SendNbOfParticles()
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
  END IF
#endif /*MPI*/

  ! because of emmission and UpdateParticlePosition
  CALL Deposition(doInnerParts=.TRUE.) 

#if USE_MPI
  IF(DoExternalParts)THEN
    ! finish communication
    CALL MPIParticleRecv()
  END IF
  ! here: finish deposition with delta kernel
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*USE_MPI*/

  CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
END IF
  
  ! set last data already here, since surfaceflux moved before interpolation
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)

IF (t.GE.DelayTime) THEN
  ! forces on particle
  ! can be used to hide sending of number of particles
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
  CALL CalcPartRHS()
#if EQNSYSNR == 4
  CALL Particle_RandomWalk()
#endif
END IF

END SUBROUTINE Particle_TimeStepByLSERK_RHS

    
!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK(t,iStage,b_dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_TimeDisc_Vars,           ONLY: RKa,nRKStages
USE MOD_Mesh_Vars,               ONLY: MeshFile
USE MOD_DG,                      ONLY: DGTimeDerivative_weakForm
USE MOD_PruettDamping,           ONLY: TempFilterTimeDeriv
USE MOD_Analyze_Vars,            ONLY: tWriteData
USE MOD_HDF5_Output,             ONLY: WriteState
#if FV_ENABLED
USE MOD_FV,                      ONLY: FV_Switch
USE MOD_FV_Vars,                 ONLY: FV_toDGinRK
#endif
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*MPI*/
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_Part_tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping,TriaTracking
USE MOD_Particle_Vars,           ONLY: PartState, Pt, Pt_temp, LastPartPos, DelayTime, PEM, PDM, Species,PartSpecies
#if EQNSYSNR == 4
USE MOD_Particle_RandomWalk,     ONLY: Particle_RandomWalk
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
INTEGER,INTENT(INOUT)         :: iStage
REAL,INTENT(IN)               :: b_dt(1:nRKStages)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: part_err
INTEGER                       :: iPart, iStage_loc
REAL                          :: RandVal,v_magnitude
REAL                          :: Pa_rebuilt_coeff(1:nRKStages),Pa_rebuilt(1:3,1:nRKStages),Pv_rebuilt(1:3,1:nRKStages),v_rebuilt(1:3,0:nRKStages-1)
!===================================================================================================================================

DO iStage_loc=1,nRKStages
  ! Rebuild Pt_tmp-coefficients assuming F=const. (value at wall) in previous stages
  IF (iStage_loc.EQ.1) THEN
    Pa_rebuilt_coeff(iStage_loc) = 1.
  ELSE
    Pa_rebuilt_coeff(iStage_loc) = 1. - RKA(iStage_loc)*Pa_rebuilt_coeff(iStage_loc-1)
  END IF
END DO

IF (t.GE.DelayTime) THEN
  part_err = .FALSE.
  
  ! particle step
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      !-- Pt is not known only for new Surfaceflux-Parts -> change IsNewPart back to F for other Parts
        IF (.NOT.PDM%dtFracPush(iPart)) PDM%IsNewPart(iPart)=.FALSE.
      !-- Particle Push
        ! Sanity Check Particle Pusher / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(Pt(iPart,:)))) THEN
            IPWRITE(UNIT_stdOut,*) 'Found invalid particle push, ignoring. PartID:', iPart
            Pt(iPart,:) = 0
        ENDIF
          
      IF (.NOT.PDM%IsNewPart(iPart)) THEN
        Pt_temp(iPart,1) = PartState(iPart,4)
        Pt_temp(iPart,2) = PartState(iPart,5)
        Pt_temp(iPart,3) = PartState(iPart,6)
        Pt_temp(iPart,4) = Pt(iPart,1)
        Pt_temp(iPart,5) = Pt(iPart,2)
        Pt_temp(iPart,6) = Pt(iPart,3)
        PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4)*b_dt(1)
        PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5)*b_dt(1)
        PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6)*b_dt(1)
        PartState(iPart,4) = PartState(iPart,4) + Pt(iPart,1)*b_dt(1)
        PartState(iPart,5) = PartState(iPart,5) + Pt(iPart,2)*b_dt(1)
        PartState(iPart,6) = PartState(iPart,6) + Pt(iPart,3)*b_dt(1)

        ! Sanity Check Particle / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(PartState(iPart,:)))) THEN
            PDM%ParticleInside(iPart) = .FALSE.
            IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        ENDIF
        
      ELSE !IsNewPart: no Pt_temp history available!
        CALL RANDOM_NUMBER(RandVal)
        Pa_rebuilt(:,:)=0.
        DO iStage_loc=1,iStage
          Pa_rebuilt(1:3,iStage_loc)=Pa_rebuilt_coeff(iStage_loc)*Pt(iPart,1:3)
        END DO
        v_rebuilt(:,:)=0.
        DO iStage_loc=iStage-1,0,-1
          IF (iStage_loc.EQ.iStage-1) THEN
            v_rebuilt(1:3,iStage_loc) = PartState(iPart,4:6) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
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
        Pt_temp(iPart,1:3) = Pv_rebuilt(1:3,iStage)
        Pt_temp(iPart,4:6) = Pa_rebuilt(1:3,iStage)
        PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)*RandVal
        PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)*RandVal
        PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)*RandVal
        PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)*RandVal
        PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)*RandVal
        PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)*RandVal
        PDM%dtFracPush(iPart) = .FALSE.
        PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...
        
!        ! Sanity Check Particle / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(PartState(iPart,:)))) THEN
            PDM%ParticleInside(iPart) = .FALSE.
            IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        ENDIF
        
      END IF !IsNewPart
    
        ! Try to find particles with too high velocity
        v_magnitude   = SQRT(DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6)))
        
    IF ((Species(PartSpecies(iPart))%HighVeloThreshold.NE.0).AND.(v_magnitude.GT.Species(PartSpecies(iPart))%HighVeloThreshold))THEN
        part_err = .TRUE.
        IPWRITE(UNIT_stdOut,*) ' High velocity particle detected. Writing error state and removing particle ...'
        IPWRITE(UNIT_stdout,*) ' LastPos:',  PartState(iPart,1:3)
        IPWRITE(UNIT_stdout,*) ' Velocity:', PartState(iPart,4:6)
        PDM%ParticleInside(iPart) = .FALSE.
        END IF
      END IF
    END DO
    
    IF (part_err) THEN
        CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
                            FutureTime=tWriteData,isErrorFile=.TRUE.)
    END IF
    
#if USE_MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
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
  CALL ParticleInserting()
#if USE_MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()   ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()   ! finish communication
#endif
!  CALL ParticleCollectCharges()
END IF
  
END SUBROUTINE Particle_TimeStepByLSERK

!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> Calculate the right hand side before updating the field solution. Can be used to hide sending of number of particles.
!===================================================================================================================================
SUBROUTINE Particle_TimeStepByLSERK_RK_RHS(t,iStage,b_dt)
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars,           ONLY: RKA,nRKStages
#if USE_MPI
!USE MOD_Particle_MPI_Vars,       ONLY: DoExternalParts
USE MOD_Particle_Mesh,           ONLY: CountPartsPerElem
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
USE MOD_PICInterpolation
USE MOD_PICDepo
USE MOD_Part_RHS,                ONLY: CalcPartRHS
#if EQNSYSNR == 4
USE MOD_Particle_RandomWalk,     ONLY: Particle_RandomWalk
#endif
USE MOD_Particle_Vars,           ONLY: PartState,DelayTime,LastPartPos,PDM,PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t
INTEGER,INTENT(INOUT)         :: iStage
REAL,INTENT(IN)               :: b_dt(1:nRKStages)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_MPI
CALL CountPartsPerElem(ResetNumberOfParticles=.FALSE.) !for scaling of tParts of LB
#endif
  
! deposition
IF (t.GE.DelayTime) THEN
  CALL Deposition(doInnerParts=.TRUE.) ! because of emission and UpdateParticlePosition
#if USE_MPI
  ! here: finish deposition with delta kernel
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*USE_MPI*/
  CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
END IF

LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)

IF (t.GE.DelayTime) THEN
  ! forces on particle
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)   ! forces on particles
  CALL CalcPartRHS()
#if EQNSYSNR == 4
  CALL Particle_RandomWalk()
#endif
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
USE MOD_Mesh_Vars,               ONLY: MeshFile
USE MOD_TimeDisc_Vars,           ONLY: RKA,nRKStages
USE MOD_PruettDamping,           ONLY: TempFilterTimeDeriv
USE MOD_Analyze_Vars,            ONLY: tWriteData
USE MOD_HDF5_Output,             ONLY: WriteState
#if FV_ENABLED
USE MOD_FV,                      ONLY: FV_Switch
USE MOD_FV_Vars,                 ONLY: FV_toDGinRK
#endif
USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping,TriaTracking
USE MOD_PICDepo,                 ONLY: Deposition
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY: PartState, Pt, Pt_temp, LastPartPos, DelayTime, PEM, PDM, Species,PartSpecies
USE MOD_part_RHS,                ONLY: CalcPartRHS
USE MOD_Particle_Tracking,       ONLY: ParticleTracing,ParticleRefTracking,ParticleTriaTracking
USE MOD_part_emission,           ONLY: ParticleInserting
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
#if USE_MPI
USE MOD_Particle_Mesh,           ONLY: CountPartsPerElem
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif /*MPI*/
#if EQNSYSNR == 4
USE MOD_Particle_RandomWalk,     ONLY: Particle_RandomWalk
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
REAL                          :: RandVal,v_magnitude
REAL                          :: Pa_rebuilt_coeff(1:nRKStages),Pa_rebuilt(1:3,1:nRKStages),Pv_rebuilt(1:3,1:nRKStages),v_rebuilt(1:3,0:nRKStages-1)
!===================================================================================================================================

DO iStage_loc=1,nRKStages
  ! Rebuild Pt_tmp-coefficients assuming F=const. (value at wall) in previous stages
  IF (iStage_loc.EQ.1) THEN
    Pa_rebuilt_coeff(iStage_loc) = 1.
  ELSE
    Pa_rebuilt_coeff(iStage_loc) = 1. - RKA(iStage_loc)*Pa_rebuilt_coeff(iStage_loc-1)
  END IF
END DO

IF (t.GE.DelayTime) THEN
  part_err = .FALSE.

  ! particle step
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      IF (.NOT.PDM%IsNewPart(iPart)) THEN
        Pt_temp(iPart,1) = PartState(iPart,4) - RKA(iStage) * Pt_temp(iPart,1)
        Pt_temp(iPart,2) = PartState(iPart,5) - RKA(iStage) * Pt_temp(iPart,2)
        Pt_temp(iPart,3) = PartState(iPart,6) - RKA(iStage) * Pt_temp(iPart,3)
        Pt_temp(iPart,4) = Pt(iPart,1) - RKA(iStage) * Pt_temp(iPart,4)
        Pt_temp(iPart,5) = Pt(iPart,2) - RKA(iStage) * Pt_temp(iPart,5)
        Pt_temp(iPart,6) = Pt(iPart,3) - RKA(iStage) * Pt_temp(iPart,6)

        !< Use current particle pusher since we already have it  
!          Pt_temp(iPart,1) = PartState(iPart,4)
!          Pt_temp(iPart,2) = PartState(iPart,5)
!          Pt_temp(iPart,3) = PartState(iPart,6)
!          Pt_temp(iPart,4) = Pt(iPart,1)
!          Pt_temp(iPart,5) = Pt(iPart,2)
!          Pt_temp(iPart,6) = Pt(iPart,3)

        ! Sanity Check Particle Pusher / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(Pt(iPart,:)))) THEN
          IPWRITE(UNIT_stdOut,*) 'Found invalid particle push, ignoring. PartID:', iPart
          Pt(iPart,:) = 0
        ENDIF

        PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)
        PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)
        PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)
        PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)
        PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)
        PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)

        ! Sanity Check Particle / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(PartState(iPart,:)))) THEN
          PDM%ParticleInside(iPart) = .FALSE.
          IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        ENDIF

      ELSE !IsNewPart: no Pt_temp history available!
        RandVal=1. !"normal" particles (i.e. not from SurfFlux) are pushed with whole timestep!
        Pa_rebuilt(:,:)=0.
        DO iStage_loc=1,iStage
          Pa_rebuilt(1:3,iStage_loc)=Pa_rebuilt_coeff(iStage_loc)*Pt(iPart,1:3)
        END DO
        v_rebuilt(:,:)=0.
        DO iStage_loc=iStage-1,0,-1
          IF (iStage_loc.EQ.iStage-1) THEN
            v_rebuilt(1:3,iStage_loc) = PartState(iPart,4:6) + (RandVal-1.)*b_dt(iStage_loc+1)*Pa_rebuilt(1:3,iStage_loc+1)
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
        Pt_temp(iPart,1:3) = Pv_rebuilt(1:3,iStage)
        Pt_temp(iPart,4:6) = Pa_rebuilt(1:3,iStage)
        PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)*RandVal
        PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)*RandVal
        PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)*RandVal
        PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)*RandVal
        PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)*RandVal
        PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)*RandVal
        PDM%IsNewPart(iPart) = .FALSE. !change to false: Pt_temp is now rebuilt...

        ! Sanity Check Particle / WARNING: Might Cause Slowdowns
        IF (ANY(ISNAN(PartState(iPart,:)))) THEN
          PDM%ParticleInside(iPart) = .FALSE.
          IPWRITE(UNIT_stdOut,*) 'Found invalid particle, removing. PartID:', iPart
        ENDIF

      END IF !IsNewPart

      ! Try to find particles with too high velocity
      v_magnitude   = SQRT(DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6)))

  IF ((Species(PartSpecies(iPart))%HighVeloThreshold.NE.0).AND.(v_magnitude.GT.Species(PartSpecies(iPart))%HighVeloThreshold))THEN
      part_err = .TRUE.
      IPWRITE(UNIT_stdOut,*) ' High velocity particle detected. Writing error state and removing particle ...'
      IPWRITE(UNIT_stdout,*) ' LastPos:',  PartState(iPart,1:3)
      IPWRITE(UNIT_stdout,*) ' Velocity:', PartState(iPart,4:6)
      PDM%ParticleInside(iPart) = .FALSE.
      END IF
    END IF
  END DO

  IF (part_err) THEN
      CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
                          FutureTime=tWriteData,isErrorFile=.TRUE.)
  END IF

  ! particle tracking
#if USE_MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
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
  CALL ParticleInserting()

#if USE_MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()   ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()   ! finish communication
#endif
END IF
  
! <<<<<
! AB HIEER vielleicht nur 1x am Ende der RK stage
IF (iStage.EQ.nRKStages) THEN
!#if USE_MPI
!  PartMPIExchange%nMPIParticles=0 ! and set number of received particles to zero for deposition
!#endif
  IF ((t.GE.DelayTime).OR.(t.EQ.0)) THEN
    CALL UpdateNextFreePosition()
  END IF

  CALL UpdateNextFreePosition()
END IF

END SUBROUTINE Particle_TimeStepByLSERK_RK

END MODULE MOD_Particle_TimeDisc
