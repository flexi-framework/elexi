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

!===================================================================================================================================
! Contains general particles variables
!===================================================================================================================================
MODULE MOD_Particle_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: useLinkedList                               ! Flag to trigger the start of the linked list for output tools

REAL    , ALLOCATABLE         :: PartState(:,:)                              ! (1:6,1:NParts) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE         :: TurbPartState(:,:)                          ! (1:4,1:NParts) with 2nd index vx',vy',vz',t_remaining
REAL    , ALLOCATABLE         :: PartPosRef(:,:)                             ! (1:3,1:NParts) particles pos mapped to -1|1 space
INTEGER , ALLOCATABLE         :: PartPosGauss(:,:)                           ! (1:3,1:NParts) Gauss point localization of particles
REAL    , ALLOCATABLE         :: Pt(:,:)                                     ! Derivative of PartState (vx,xy,vz) only
INTEGER , ALLOCATABLE         :: PartReflCount(:)                            ! Counter of number of reflections

REAL    , ALLOCATABLE         :: Pt_temp(:,:)                                ! LSERK4 additional derivative of PartState
                                                                             ! (1:6,1:NParts) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE         :: TurbPt_temp(:,:)                            ! LSERK4 additional derivative of TurbPartState
                                                                             ! (1:3,1:NParts) with 2nd index: vx',vy',vz'
REAL    , ALLOCATABLE         :: LastPartPos(:,:)                            ! (1:3,1:NParts) with 2nd index: x,y,z
INTEGER , ALLOCATABLE         :: PartSpecies(:)                              ! (1:NParts)
INTEGER(KIND=8), ALLOCATABLE  :: PartIndex(:)                                ! (1:NParts)
INTEGER                       :: sumOfMatchedParticlesSpecies                ! previous sumOfMatchedParticles for all species
INTEGER                       :: PartRHSMethod
REAL                          :: PartGravity(3)
INTEGER                       :: nrSeeds                                     ! Number of Seeds for Random Number Generator
INTEGER , ALLOCATABLE         :: seeds(:)                                    ! =>NULL()   ! Seeds for Random Number Generator

LOGICAL                       :: doPartIndex                                 ! Flag to give particles an unique (or not) index

#if USE_BASSETFORCE
REAL    , ALLOCATABLE         :: durdt(:,:)                                  ! Old dur/dt for Basset force
INTEGER                       :: N_Basset                                    ! Number of old dur/dt terms used for the Basset force
INTEGER                       :: bIter
INTEGER                       :: nBassetVars
#endif

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
  REAL                                   :: ExcludeBV_lenghts(2)             ! lenghts of BV1/2 (to be calculated)
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
  REAL                                   :: InflowRiseTime                   ! time to ramp the number of inflow particles
                                                                             ! linearly from zero to unity
  REAL                                   :: NormalIC(3)                      ! Normal / Orientation of circle
  REAL                                   :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  REAL                                   :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  REAL                                   :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  REAL                                   :: BaseVariance                     ! variance for Gaussian distribtion
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
  REAL                                   :: Alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3)
  INTEGER                                :: ParticleEmissionType             ! Emission Type 1 = emission rate in 1/s,
  INTEGER                                :: CountIndex
                                                                             !               2 = emission rate 1/iteration
  REAL                                   :: ParticleEmission                 ! Emission in [1/s] or [1/Iteration]
  REAL                                   :: ParticleEmissionTime             ! Scale emission time for emission type 1
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=4)                        :: InsertedParticleMisMatch=0       ! error in number of inserted particles of last step
  INTEGER                                :: NumberOfExcludeRegions           ! Number of different regions to be excluded
  TYPE(tExcludeRegion), ALLOCATABLE      :: ExcludeRegion(:)
#if USE_MPI
  INTEGER                                :: InitComm                         ! number of init-communicator
#endif /*USE_MPI*/
  INTEGER                                :: sumOfMatchedParticles            ! Sum of matched particles on all procs
  INTEGER                                :: sumOfRequestedParticles          ! Sum of requested particles on all procs
  INTEGER                                :: mySumOfMatchedParticles          ! Sum of matched particles on current proc
  INTEGER                                :: nPartsPerProc                    ! Particle displacement on current proc
  REAL, ALLOCATABLE                      :: PartField(:,:,:)
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
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3)
  LOGICAL                                :: VeloIsNormal                     ! VeloIC is in Surf-Normal instead of VeloVecIC
  LOGICAL                                :: ReduceNoise                      ! reduce stat. noise by global calc. of PartIns
  LOGICAL                                :: AcceptReject                     ! perform ARM for skewness of RefMap-positioning
  INTEGER                                :: ARM_DmaxSampleN                  ! number of sample intervals in xi/eta for Dmax-calc.
  REAL                                   :: VFR_total                        ! Total Volumetric flow rate through surface
  REAL                                   :: VFR_total_allProcsTotal          !     -''-, total
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=8)                        :: tmpInsertedParticle              ! tmp Number of all already inserted Particles
  INTEGER(KIND=8)                        :: tmpInsertedParticleSurplus       ! tmp accumulated "negative" number of inserted Particles
  TYPE(tSurfFluxSubSideData),ALLOCATABLE :: SurfFluxSubSideData(:,:,:)       ! SF-specific Data of Sides (1:N,1:N,1:SideNumber)
  LOGICAL                                :: CircularInflow                   ! Circular region, which can be used to define small
                                                                             ! geometry features on large boundaries
  INTEGER                                :: dir(3)                           ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                                   :: origin(2)                        ! origin in orth. coordinates of polar system
  REAL                                   :: rmax                             ! max radius of to-be inserted particles
  REAL                                   :: rmin                             ! min radius of to-be inserted particles
  INTEGER, ALLOCATABLE                   :: SurfFluxSideRejectType(:)        ! Type if parts in side can be rejected (1:SideNumber)
  REAL, ALLOCATABLE                      :: ConstMassflowWeight(:,:,:)       ! Adaptive, Type 4: Weighting factor for SF-sides to
                                                                             ! insert the right amount of particles
  REAL                                   :: SampledMassflow                  ! Actual mass flow rate through a surface flux boundary
  REAL, ALLOCATABLE                      :: nVFRSub(:,:)                     ! normal volume flow rate through subsubside
END TYPE

TYPE tSpecies                                                                ! Particle Data for each Species
  !General Species Values
  TYPE(tInit), ALLOCATABLE               :: Init(:)  !     =>NULL()          ! Particle Data for each Initialisation
  INTEGER                                :: RHSMethod                        ! specifying Keyword for RHS method
  REAL                                   :: MassIC                           ! Particle mass (without MPF)
  REAL                                   :: DiameterIC                       ! Particle diameter (without MPF)
  REAL                                   :: DensityIC                        ! Particle density (without MPF)
  INTEGER                                :: NumberOfInits                    ! Number of different initial particle placements
  TYPE(typeSurfaceflux),ALLOCATABLE      :: Surfaceflux(:)                   ! Particle Data for each SurfaceFlux emission
  INTEGER                                :: nSurfacefluxBCs                  ! Number of SF emissions
  INTEGER                                :: StartnumberOfInits               ! 0 if old emit defined (array is copied into 0. entry)
  REAL                                   :: LowVeloThreshold                 ! Threshold value for removal of low velocity particles
  REAL                                   :: HighVeloThreshold                ! Threshold value for removal of high velocity particle
  INTEGER                                :: LowVeloCounter                   ! Counter how many low velocity particles were removed
  REAL                                   :: SphericityIC                     ! Particle sphericity
  ! Bons particle rebound model
  REAL                                   :: YoungIC                          ! Young's modulus
  REAL                                   :: PoissonIC                        ! Poisson's ration for transverse strain under ax. comp
  REAL                                   :: YieldCoeff                       ! Yield strength coefficient
  INTEGER                                :: CountIndex                       ! Count number of particles
#if USE_EXTEND_RHS
  LOGICAL                                :: CalcLiftForce                    ! Calculate the lift (Saffman) force
  LOGICAL                                :: CalcUndisturbedFlow              ! Calculate the undisturbed flow force
  LOGICAL                                :: CalcVirtualMass                  ! Calculate the virtual mass force
  LOGICAL                                :: CalcBassetForce                  ! Calculate the Basset force
#endif /* USE_EXTEND_RHS */
END TYPE

INTEGER                                  :: nSpecies                         ! number of species
TYPE(tSpecies), ALLOCATABLE              :: Species(:)        !   => NULL()  ! Species Data Vector

TYPE tParticleElementMapping
  INTEGER                , ALLOCATABLE   :: Element(:)        !   =>NULL()   ! Element number allocated to each Particle
  INTEGER                , ALLOCATABLE   :: lastElement(:)    !   =>NULL()   ! Element number allocated
!----------------------------------------------------------------------------!----------------------------------
                                                                             ! Following vectors are assigned in
                                                                             ! SUBROUTINE UpdateNextFreePosition
  INTEGER                , ALLOCATABLE   :: pStart(:)         !   =>NULL()   ! Start of Linked List for Particles in Element
                                                                             ! pStart(1:PEM%nElem)
  INTEGER                , ALLOCATABLE   :: pNumber(:)        !   =>NULL()   ! Number of Particles in Element
                                                                             ! pStart(1:PEM%nElem)
  INTEGER                , ALLOCATABLE   :: wNumber(:)        !   =>NULL()   ! Number of Wall-Particles in Element
                                                                             ! pStart(1:PEM%nElem)
  INTEGER                , ALLOCATABLE   :: pEnd(:)           !   =>NULL()   ! End of Linked List for Particles in Element
                                                                             ! pEnd(1:PEM%nElem)
  INTEGER                , ALLOCATABLE   :: pNext(:)          !   =>NULL()   ! Next Particle in same Element (Linked List)
                                                                             ! pStart(1:PEM%maxParticleNumber)
END TYPE

TYPE(tParticleElementMapping)            :: PEM

TYPE tParticleDataManagement
  INTEGER                                :: CurrentNextFreePosition          ! Index of nextfree index in nextFreePosition-Array
  INTEGER                                :: maxParticleNumber                ! Maximum Number of all Particles
  INTEGER                                :: ParticleVecLength                ! Vector Length for Particle Push Calculation
  INTEGER                                :: insideParticleNumber             ! Number of all recent Particles inside
  INTEGER , ALLOCATABLE                  :: PartInit(:)                      ! (1:NParts), initial emission condition number
                                                                             ! the calculation area
  INTEGER ,ALLOCATABLE                   :: nextFreePosition(:)  ! =>NULL()  ! next_free_Position(1:max_Particle_Number)
                                                                             ! List of free Positon
  LOGICAL ,ALLOCATABLE                   :: ParticleInside(:)    ! =>NULL()  ! Particle_inside(1:Particle_Number)
  LOGICAL ,ALLOCATABLE                   :: IsNewPart(:)                     ! Reconstruct RK-scheme in next stage
END TYPE

TYPE (tParticleDataManagement)           :: PDM

REAL                                     :: DelayTime

LOGICAL                                  :: ParticlesInitIsDone=.FALSE.

LOGICAL                                  :: DoSurfaceFlux                    ! Flag for emitting by SurfaceFluxBCs
LOGICAL                                  :: UseCircularInflow                !
INTEGER, ALLOCATABLE                     :: CountCircInflowType(:,:,:)
LOGICAL                                  :: DoPoissonRounding                ! Perform Poisson sampling instead of random rounding
LOGICAL                                  :: DoTimeDepInflow                  ! Insertion and SurfaceFlux w simple random rounding
LOGICAL                                  :: RepWarn = .FALSE.                ! Warning for Reynolds limit of particle model

#if USE_EXTEND_RHS && ANALYZE_RHS
REAL                                     :: tWriteRHS
REAL                                     :: dtWriteRHS
CHARACTER(LEN=255)                       :: Filename_RHS
#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

!===================================================================================================================================
END MODULE MOD_Particle_Vars
