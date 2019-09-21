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
! Contains the Particles' variables (general for all modules: PIC, DSMC, FP)
!===================================================================================================================================
MODULE MOD_Particle_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL, PARAMETER       :: BoltzmannConst=1.380648813E-23                      ! Boltzmann constant [J/K] SI-Unit! in m^2/(s^2*K)
REAL                  :: ManualTimeStep                                      ! Manual TimeStep
LOGICAL               :: useManualTimeStep                                   ! Logical Flag for manual timestep. For consistency
                                                                             ! with IAG programming style
LOGICAL               :: AllowLoosing                                        ! Flag if a lost particle should abort the program
LOGICAL               :: KeepWallParticles                                   ! Flag for tracking of adsorbed Particles
LOGICAL               :: printRandomSeeds                                    ! print random seeds or not
REAL                  :: dt_max_particles                                    ! Maximum timestep for particles (for static fields!)
!INTEGER               :: WeirdElems                                          ! Number of Weird Elements (=Elements which are folded
                                                                              ! into themselves)

REAL    , ALLOCATABLE :: PartState(:,:)                                      ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
#if USE_RW
REAL    , ALLOCATABLE :: PartTurbState(:,:)                                  ! (1:NParts,1:4) with 2nd index vx',vy',vz',t_remaining
#endif
REAL    , ALLOCATABLE :: PartPosRef(:,:)                                     ! (1:3,1:NParts) particles pos mapped to -1|1 space
INTEGER , ALLOCATABLE :: PartPosGauss(:,:)                                   ! (1:NParts,1:3) Gauss point localization of particles
REAL    , ALLOCATABLE :: Pt(:,:)                                             ! Derivative of PartState (vx,xy,vz) only
INTEGER , ALLOCATABLE :: PartReflCount(:)                                    ! Counter of number of reflections

REAL    , ALLOCATABLE :: Pt_temp(:,:)                                        ! LSERK4 additional derivative of PartState                                                          ! (1:NParts,1:6) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE :: LastPartPos(:,:)                                    ! (1:NParts,1:3) with 2nd index: x,y,z
INTEGER , ALLOCATABLE :: PartSpecies(:)                                      ! (1:NParts)
INTEGER               :: PartRHSMethod
REAL                  :: PartGravity(3)
INTEGER               :: nrSeeds                                             ! Number of Seeds for Random Number Generator
INTEGER , ALLOCATABLE :: seeds(:)                        !        =>NULL()   ! Seeds for Random Number Generator

TYPE tExcludeRegion
  CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  REAL                                   :: RadiusIC                         ! Radius for IC circle
  REAL                                   :: Radius2IC                        ! Radius2 for IC cylinder (ring)
  REAL                                   :: NormalIC(3)                      ! Normal / Orientation of cylinder (altern. to BV1/2)
  REAL                                   :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  REAL                                   :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  REAL                                   :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  REAL                                   :: CuboidHeightIC                   ! third measure of cuboid
                                                                             ! (set 0 for flat rectangle),
                                                                             ! negative value = opposite direction
  REAL                                   :: CylinderHeightIC                 ! third measure of cylinder
                                                                             ! (set 0 for flat circle),
                                                                             ! negative value = opposite direction
  REAL                                   :: ExcludeBV_lenghts(2)                    ! lenghts of BV1/2 (to be calculated)
END TYPE

TYPE tInit                                                                   ! Particle Data for each init emission for each species
  !Specific Emission/Init values
  LOGICAL                                :: UseForInit                       ! Use Init/Emission for init.?
  LOGICAL                                :: UseForEmission                   ! Use Init/Emission for emission?
  CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  INTEGER(8)                             :: initialParticleNumber            ! Number of Particles at time 0.0
  REAL                                   :: RadiusIC                         ! Radius for IC circle
  REAL                                   :: Radius2IC                        ! Radius2 for IC cylinder (ring)
  REAL                                   :: RadiusICGyro                     ! Radius for Gyrotron gyro radius
  INTEGER                                :: Rotation                         ! direction of rotation, similar to TE-mode
  INTEGER                                :: VelocitySpreadMethod             ! method to compute the velocity spread
  REAL                                   :: InflowRiseTime                   ! time to ramp the number of inflow particles
                                                                             ! linearly from zero to unity
  REAL                                   :: VelocitySpread                   ! velocity spread in percent
  REAL                                   :: NormalIC(3)                      ! Normal / Orientation of circle
  REAL                                   :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  REAL                                   :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  REAL                                   :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  REAL                                   :: CuboidHeightIC                   ! third measure of cuboid
                                                                             ! (set 0 for flat rectangle),
                                                                             ! negative value = opposite direction
  REAL                                   :: CylinderHeightIC                 ! third measure of cylinder
                                                                             ! (set 0 for flat rectangle),
                                                                             ! negative value = opposite direction
  LOGICAL                                :: CalcHeightFromDt                 ! Calc. cuboid/cylinder height from v and dt?
  REAL                                   :: VeloIC                           ! velocity for inital Data
  REAL                                   :: VeloTurbIC                       ! turbulent velocity fluctuation for inital Data
  REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                                   :: Amplitude                        ! Amplitude for sin-deviation initiation.
  REAL                                   :: WaveNumber                       ! WaveNumber for sin-deviation initiation.
  INTEGER                                :: maxParticleNumberX               ! Maximum Number of all Particles in x direction
  INTEGER                                :: maxParticleNumberY               ! Maximum Number of all Particles in y direction
  INTEGER                                :: maxParticleNumberZ               ! Maximum Number of all Particles in z direction
  REAL                                   :: MJxRatio                         ! x direction portion of velocity for Maxwell-Juettner
  REAL                                   :: MJyRatio                         ! y direction portion of velocity for Maxwell-Juettner
  REAL                                   :: MJzRatio                         ! z direction portion of velocity for Maxwell-Juettner
  REAL                                   :: WeibelVeloPar                    ! Parallel velocity component for Weibel
  REAL                                   :: WeibelVeloPer                    ! Perpendicular velocity component for Weibel
  REAL                                   :: OneDTwoStreamVelo                ! Stream Velocity for the Two Stream Instability
  REAL                                   :: OneDTwoStreamTransRatio          ! Ratio between perpendicular and parallel velocity
  REAL                                   :: Alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: ConstantPressure                 ! Pressure for an Area with a Constant Pressure
  REAL                                   :: ConstPressureRelaxFac            ! RelaxFac. for ConstPressureSamp
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3) for (vpi_)cub./
                                                                             ! cyl. as alternative to Part.Emis. in Type1
  INTEGER                                :: ParticleEmissionType             ! Emission Type 1 = emission rate in 1/s,
                                                                             !               2 = emission rate 1/iteration
                                                                             !               3 = user def. emission rate
                                                                             !               4 = const. cell pressure
                                                                             !               5 = cell pres. w. complete part removal
                                                                             !               6 = outflow BC (characteristics method)
  REAL                                   :: ParticleEmission                 ! Emission in [1/s] or [1/Iteration]
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=4)                        :: InsertedParticleMisMatch=0       ! error in number of inserted particles of last step
  REAL                                   :: Nsigma                           ! sigma multiple of maxwell for virtual insert length
  LOGICAL                                :: VirtPreInsert                    ! virtual Pre-Inserting region (adapted SetPos/Velo)?
  CHARACTER(40)                          :: vpiDomainType                    ! specifying Keyword for virtual Pre-Inserting region
                                                                             ! implemented: - perpendicular_extrusion (default)
                                                                             !              - freestream
                                                                             !              - orifice
                                                                             !              - ...more following...
  LOGICAL                                :: vpiBVBuffer(4)                   ! incl. buffer region in -BV1/+BV1/-BV2/+BV2 direction?
  !TYPE(tConstPressure)                   :: ConstPress!(:)           =>NULL() !
  INTEGER                                :: NumberOfExcludeRegions           ! Number of different regions to be excluded
  TYPE(tExcludeRegion), ALLOCATABLE      :: ExcludeRegion(:)
#if USE_MPI
  INTEGER                                :: InitComm                          ! number of init-communicator
#endif /*MPI*/
END TYPE tInit

TYPE tSurfFluxSubSideData
  REAL                                   :: projFak                          ! VeloVecIC projected to inwards normal
  REAL                                   :: a_nIn                            ! speed ratio projected to inwards normal
  REAL                                   :: Velo_t1                          ! Velo comp. of first orth. vector
  REAL                                   :: Velo_t2                          ! Velo comp. of second orth. vector
  REAL                                   :: nVFR                             ! normal volume flow rate through subside
  REAL                                   :: Dmax                             ! maximum Jacobian determinant of subside for opt. ARM
  REAL,ALLOCATABLE                       :: BezierControlPoints2D(:,:,:)     ! BCP of SubSide projected to VeloVecIC
                                                                             ! (1:2,0:NGeo,0:NGeo)
END TYPE tSurfFluxSubSideData

TYPE typeSurfaceflux
  INTEGER                                :: BC                               ! PartBound to be emitted from
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  REAL                                   :: VeloIC                           ! velocity for inital Data
  REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3)
  LOGICAL                                :: VeloIsNormal                     ! VeloIC is in Surf-Normal instead of VeloVecIC
  LOGICAL                                :: ReduceNoise                      ! reduce stat. noise by global calc. of PartIns
  LOGICAL                                :: AcceptReject                     ! perform ARM for skewness of RefMap-positioning
  INTEGER                                :: ARM_DmaxSampleN                  ! number of sample intervals in xi/eta for Dmax-calc.
  REAL                                   :: VFR_total                        ! Total Volumetric flow rate through surface
  REAL                     , ALLOCATABLE :: VFR_total_allProcs(:)            ! -''-, all values for root in ReduceNoise-case
  REAL                                   :: VFR_total_allProcsTotal          !     -''-, total
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=8)                        :: tmpInsertedParticle              ! tmp Number of all already inserted Particles
  INTEGER(KIND=8)                        :: tmpInsertedParticleSurplus       ! tmp accumulated "negative" number of inserted Particles
  TYPE(tSurfFluxSubSideData), ALLOCATABLE :: SurfFluxSubSideData(:,:,:)      ! SF-specific Data of Sides (1:N,1:N,1:SideNumber)
  INTEGER, ALLOCATABLE                   :: SurfFluxSideRejectType(:)        ! Type if parts in side can be rejected (1:SideNumber)
  LOGICAL                                :: SimpleRadialVeloFit              !fit of veloR/veloTot=-r*(A*exp(B*r)+C)
  REAL                                   :: preFac                           !A
  REAL                                   :: powerFac                         !B
  REAL                                   :: shiftFac                         !C
  INTEGER                                :: dir(3)                           ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                                   :: origin(2)                        ! origin in orth. coordinates of polar system
  REAL                                   :: rmax                             ! max radius of to-be inserted particles
  REAL                                   :: rmin                             ! min radius of to-be inserted particles
  REAL                                   :: PressureFraction
END TYPE

TYPE tSpecies                                                                ! Particle Data for each Species
  !General Species Values
  TYPE(tInit), ALLOCATABLE               :: Init(:)  !     =>NULL()          ! Particle Data for each Initialisation
  CHARACTER(40)                          :: RHSMethod                        ! specifying Keyword for RHS method
  REAL                                   :: ChargeIC                         ! Particle Charge (without MPF)
  REAL                                   :: MassIC                           ! Particle Mass (without MPF)
  REAL                                   :: DensityIC                        ! Particle Density (without MPF)
  REAL                                   :: MacroParticleFactor              ! Number of Microparticle per Macroparticle
  INTEGER                                :: NumberOfInits                    ! Number of different initial particle placements
  INTEGER                                :: StartnumberOfInits               ! 0 if old emit defined (array is copied into 0. entry)
  TYPE(typeSurfaceflux),ALLOCATABLE      :: Surfaceflux(:)                   ! Particle Data for each SurfaceFlux emission
  INTEGER                                :: nSurfacefluxBCs                  ! Number of SF emissions
  REAL                                   :: LowVeloThreshold                 ! Threshold value for removal of low velocity particles
  REAL                                   :: HighVeloThreshold                ! Threshold value for removal of high velocity particle
  INTEGER                                :: LowVeloCounter                   ! Counter how many low velocity particles were removed
  !Random Walk method
  CHARACTER(40)                          :: RWModel                          ! specifying Keyword for RW model
  ! Bons particle rebound model
  REAL                                   :: YoungIC                          ! Young's modulus
  REAL                                   :: PoissonIC                        ! Poisson's ration for transverse strain under ax. comp
  REAL                                   :: YieldCoeff                       ! Yield strength coefficient
END TYPE

REAL, ALLOCATABLE                        :: Adaptive_MacroVal(:,:,:)
REAL,ALLOCATABLE                         :: MacroRestartData_tmp(:,:,:,:)

INTEGER                                  :: nSpecies                         ! number of species
INTEGER                                  :: nMacroRestartFiles               ! number of macroscopic restart files used for particles
TYPE(tSpecies), ALLOCATABLE              :: Species(:)      !      => NULL() ! Species Data Vector

TYPE tParticleElementMapping
  INTEGER                , ALLOCATABLE   :: Element(:)      !      =>NULL()  ! Element number allocated to each Particle
  INTEGER                , ALLOCATABLE   :: lastElement(:)  !      =>NULL()  ! Element number allocated
!----------------------------------------------------------------------------!----------------------------------
                                                                             ! Following vectors are assigned in
                                                                             ! SUBROUTINE UpdateNextFreePosition
                                                                             ! IF (PIC%withDSMC .OR. PIC%withFP)
  INTEGER                , ALLOCATABLE    :: pStart(:)         !   =>NULL()  ! Start of Linked List for Particles in Element
                                                               !             ! pStart(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: pNumber(:)        !   =>NULL()  ! Number of Particles in Element
                                                               !             ! pStart(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: wNumber(:)        !   =>NULL()  ! Number of Wall-Particles in Element
                                                                             ! pStart(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: pEnd(:)           !   =>NULL()  ! End of Linked List for Particles in Element
                                                               !             ! pEnd(1:PIC%nElem)
  INTEGER                , ALLOCATABLE    :: pNext(:)          !   =>NULL()  ! Next Particle in same Element (Linked List)
                                                                             ! pStart(1:PIC%maxParticleNumber)
END TYPE

TYPE(tParticleElementMapping)            :: PEM

TYPE tParticleDataManagement
  INTEGER                                :: CurrentNextFreePosition         ! Index of nextfree index in nextFreePosition-Array
  INTEGER                                :: maxParticleNumber               ! Maximum Number of all Particles
  INTEGER                                :: ParticleVecLength               ! Vector Length for Particle Push Calculation
  INTEGER                                :: insideParticleNumber            ! Number of all recent Particles inside
  INTEGER , ALLOCATABLE                  :: PartInit(:)                     ! (1:NParts), initial emission condition number
                                                                            ! the calculation area
  INTEGER ,ALLOCATABLE                   :: nextFreePosition(:)  ! =>NULL() ! next_free_Position(1:max_Particle_Number)
                                                                            ! List of free Positon
  LOGICAL ,ALLOCATABLE                   :: ParticleInside(:)    ! =>NULL() ! Particle_inside(1:Particle_Number)
  LOGICAL , ALLOCATABLE                  :: ParticleAtWall(:)               ! Particle_adsorbed_on_to_wall(1:Particle_number)
  INTEGER , ALLOCATABLE                  :: PartAdsorbSideIndx(:,:)         ! Surface index on which Particle i adsorbed
                                                                            ! (1:3,1:PDM%maxParticleNumber)
                                                                            ! 1: surface index ElemToSide(i,localsideID,ElementID)
                                                                            ! 2: p
                                                                            ! 3: q
  LOGICAL ,ALLOCATABLE                   :: dtFracPush(:)                   ! Push random fraction only
  LOGICAL ,ALLOCATABLE                   :: IsNewPart(:)                    ! Reconstruct RK-scheme in next stage
END TYPE

TYPE (tParticleDataManagement)           :: PDM

REAL                                     :: DelayTime

LOGICAL                                  :: ParticlesInitIsDone=.FALSE.

LOGICAL                                  :: WRITEMacroValues = .FALSE.
LOGICAL                                  :: WriteMacroVolumeValues =.FALSE.   ! Output of macroscopic values in volume
LOGICAL                                  :: WriteMacroSurfaceValues=.FALSE.   ! Output of macroscopic values on surface
INTEGER                                  :: MacroValSamplIterNum              ! Number of iterations for sampling
                                                                              ! macroscopic values
REAL                                     :: MacroValSampTime                  ! Sampling time for WriteMacroVal. (e.g., for td201)
INTEGER                                  :: NumRanVec                         ! Number of predefined random vectors
REAL, ALLOCATABLE                        :: RandomVec(:,:)                    ! Random Vectos (NumRanVec, direction)
LOGICAL                                  :: OutputVpiWarnings                 ! Flag for warnings for rejected v if VPI+PartDensity
LOGICAL                                  :: DoSurfaceFlux                     ! Flag for emitting by SurfaceFluxBCs
LOGICAL                                  :: DoPoissonRounding                 ! Perform Poisson sampling instead of random rounding
LOGICAL                                  :: DoTimeDepInflow                   ! Insertion and SurfaceFlux w simple random rounding
LOGICAL                                  :: DoZigguratSampling                ! Sample normal randoms with Ziggurat method
LOGICAL                                  :: FindNeighbourElems=.FALSE.
LOGICAL                                  :: RepWarn = .FALSE.                 ! Warning for Reynolds limit of particle model

INTEGER(8)                               :: nTotalPart
INTEGER(8)                               :: nTotalHalfPart

INTEGER :: nCollectChargesBCs
INTEGER :: nDataBC_CollectCharges

!===================================================================================================================================
END MODULE MOD_Particle_Vars
