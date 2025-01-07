!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: offsetComputeNodeParts                      !
INTEGER                       :: nComputeNodeParts                           !
INTEGER                       :: nComputeNodeTotalParts                      !
INTEGER                       :: nGlobalParts                                !
INTEGER,ALLOCATABLE           :: offsetPartMPI(:)                            !< gives offsetposotion of particles of all procs,
                                                                             !< allocated from 0:nProcessors set

LOGICAL                       :: useLinkedList                               ! Flag to trigger the start of the linked list for output tools

REAL    , ALLOCATABLE         :: PartState(:,:)                              ! (1:6,1:NParts) with 2nd index: x,y,z,vx,vy,vz
REAL    , ALLOCATABLE         :: TurbPartState(:,:)                          ! (1:4,1:NParts) with 2nd index vx',vy',vz',t_remaining
REAL    , ALLOCATABLE         :: PartPosRef(:,:)                             ! (1:3,1:NParts) particles pos mapped to -1|1 space
! INTEGER , ALLOCATABLE         :: PartPosGauss(:,:)                           ! (1:3,1:NParts) Gauss point localization of particles
INTEGER , ALLOCATABLE         :: PartReflCount(:)                            ! Counter of number of reflections
REAL    , ALLOCATABLE         :: LastPartPos(:,:)                            ! (1:3,1:NParts) with 2nd index: x,y,z
INTEGER , ALLOCATABLE         :: PartSpecies(:)                              ! (1:NParts)
INTEGER(KIND=DP), ALLOCATABLE :: PartIndex(:)                                ! (1:NParts)
REAL    , ALLOCATABLE         :: Pt(:,:)                                     ! Derivative of PartState (vx,xy,vz) only
REAL    , ALLOCATABLE         :: Pt_temp(:,:)                                ! LSERK4 additional derivative of PartState
                                                                             ! (1:6,1:NParts) with 2nd index: x,y,z,vx,vy,vz


! Particle information ordered along SFC
REAL    , ALLOCATABLE         :: PartData(:,:)                               ! PartState ordered along SFC, particle number per
                                                                             ! element given in PartInt
INTEGER , ALLOCATABLE         :: PartInt(:,:)                                ! Particle number per element
REAL    , ALLOCATABLE         :: TurbPartData(:,:)                           ! Turbulent particle data
REAL    , ALLOCATABLE         :: TurbPt_temp(:,:)                            ! LSERK4 additional derivative of TurbPartState
                                                                             ! (1:3,1:NParts) with 2nd index: vx',vy',vz'
INTEGER                       :: PartDataSize     = 0                        ! number of entries in each line of PartData
INTEGER                       :: TurbPartDataSize = 0                        ! number of turbulent properties with current setup

INTEGER                       :: sumOfMatchedParticlesSpecies                ! previous sumOfMatchedParticles for all species
REAL                          :: PartGravity(3)
! INTEGER                       :: nrSeeds                                     ! Number of Seeds for Random Number Generator
INTEGER , ALLOCATABLE         :: seeds(:)                                    ! =>NULL()   ! Seeds for Random Number Generator

LOGICAL                       :: doPartIndex                                 ! Flag to give particles an unique (or not) index
LOGICAL                       :: doCalcPartSource                            ! Flag to enable two-way coupling (NN)
LOGICAL                       :: doCalcPartCollision                         ! Flag to enable four-way coupling (NN)
LOGICAL                       :: doWritePartDiam                             ! Flag to enable writeout of particle diameter
LOGICAL                       :: doRandomPartDiam                            ! Flag to enable random particle diameter
#if USE_SPHERICITY
LOGICAL                       :: doRandomSphericity                          ! Flag to enable random sphericity
#endif

#if PARTICLES_COUPLING >= 2
REAL    , ALLOCATABLE         :: PartNodeSource(:,:,:)                       ! Source term for DG solution (1:PP_nVar,1:8,1:NParts)
#endif /*PARTICLES_COUPLING >= 2*/

#if USE_BASSETFORCE
REAL    , ALLOCATABLE         :: durdt(:,:)                                  ! Old dur/dt for Basset force
REAL    , ALLOCATABLE         :: FbCoeff(:)                                  ! Coeff. for Basset force
REAL    , ALLOCATABLE         :: Fbi(:,:,:)
REAL    , ALLOCATABLE         :: Fbdt(:,:)
INTEGER                       :: N_Basset                                    ! Number of old dur/dt terms used for the Basset force
INTEGER , ALLOCATABLE         :: bIter(:)
INTEGER                       :: nBassetVars
INTEGER                       :: FbCoeffm
REAL    , ALLOCATABLE         :: FbCoeffa(:),FbCoefft(:)
#endif

CHARACTER(255)                :: FilenameRecordPart

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
END TYPE tExcludeRegion

TYPE tInit                                                                   ! Particle Data for each init emission for each species
  !Specific Emission/Init values
  LOGICAL                                :: UseForInit                       ! Use Init/Emission for init.?
  LOGICAL                                :: UseForEmission                   ! Use Init/Emission for emission?
  CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
#if USE_PARTTEMP
  CHARACTER(30)                          :: tempDistribution                 ! specifying keyword for temperature distribution
#endif /*USE_PARTTEMP*/
  REAL                                   :: initialParticleNumber            ! Number of Particles at time 0.0
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
#if USE_PARTTEMP
  REAL                                   :: TemperatureIC                    ! Particle initial temperature
#endif /*USE_PARTTEMP*/
  REAL                                   :: Alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3)
  INTEGER                                :: ParticleEmissionType             ! Emission Type 1 = emission rate in 1/s,
                                                                             !               2 = emission rate 1/iteration
  REAL                                   :: ParticleEmission                 ! Emission in [1/s] or [1/Iteration]
  REAL                                   :: ParticleEmissionTime             ! Scale emission time for emission type 1
  INTEGER(KIND=DP)                       :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=DP)                       :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=SP)                       :: InsertedParticleMisMatch=0       ! error in number of inserted particles of last step
  INTEGER                                :: NumberOfExcludeRegions           ! Number of different regions to be excluded
  TYPE(tExcludeRegion), ALLOCATABLE      :: ExcludeRegion(:)
#if USE_MPI
  INTEGER                                :: InitComm                         ! number of init-communicator
#endif /*USE_MPI*/
  INTEGER                                :: sumOfMatchedParticles            ! Sum of matched particles on all procs
  INTEGER                                :: sumOfRequestedParticles          ! Sum of requested particles on all procs
  INTEGER                                :: mySumOfMatchedParticles          ! Sum of matched particles on current proc
  INTEGER                                :: nPartsPerProc                    ! Particle displacement on current proc
  REAL, ALLOCATABLE                      :: PartField(:,:)
  INTEGER                                :: nPartField
END TYPE tInit

TYPE tSurfFluxSubSideData
  REAL                                   :: projFak                          ! VeloVecIC projected to inwards normal
  ! REAL                                   :: a_nIn                            ! speed ratio projected to inwards normal
  REAL                                   :: Velo_t1                          ! Velo comp. of first orth. vector
  REAL                                   :: Velo_t2                          ! Velo comp. of second orth. vector
  REAL                                   :: nVFR                             ! normal volume flow rate through subside
  REAL                                   :: Dmax                             ! maximum Jacobian determinant of subside for opt. ARM
  REAL,ALLOCATABLE                       :: BezierControlPoints2D(:,:,:)     ! BCP of SubSide projected to VeloVecIC
END TYPE tSurfFluxSubSideData

TYPE tSurfaceFlux
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
  INTEGER(KIND=DP)                       :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=DP)                       :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  INTEGER(KIND=DP)                       :: tmpInsertedParticle              ! tmp Number of all already inserted Particles
  INTEGER(KIND=DP)                       :: tmpInsertedParticleSurplus       ! tmp accumulated "negative" number of inserted Particles
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
  REAL                                   :: DGMeanPrimState(1:PP_nVarPrim)   ! mean DG flux through boundary
END TYPE tSurfaceFlux

ABSTRACT INTERFACE
  FUNCTION DragFactorInt(Rep,SphericityIC,Mp) RESULT(f)
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    REAL,INTENT(IN) :: Rep,SphericityIC,Mp
    REAL            :: f
  END FUNCTION DragFactorInt
END INTERFACE

TYPE type_F
  PROCEDURE(DragFactorInt),POINTER,NOPASS:: op
END TYPE type_F

TYPE tSpecies                                                                ! Particle Data for each Species
  ! General Species Values
  TYPE(tInit), ALLOCATABLE               :: Init(:)  !     =>NULL()          ! Particle Data for each Initialisation
  INTEGER                                :: RHSMethod                        ! specifying Keyword for RHS method
  TYPE(type_F)                           :: DragFactor_pointer               ! Pointer defining the drag factor
  REAL                                   :: MassIC                           ! Particle mass (without MPF)
  REAL                                   :: DiameterIC                       ! Particle diameter (without MPF)
  REAL                                   :: DensityIC                        ! Particle density (without MPF)
  REAL                                   :: StokesIC                         ! Particle Stokes number (without MPF)
#if USE_PARTTEMP
  REAL                                   :: SpecificHeatIC                   ! Particle specific heat
#endif /*USE_PARTTEMP*/
  INTEGER                                :: NumberOfInits                    ! Number of different initial particle placements
  ! SurfaceFlux
  TYPE(tSurfaceFlux),ALLOCATABLE      :: SurfaceFlux(:)                   ! Particle Data for each SurfaceFlux emission
  INTEGER                                :: nSurfacefluxBCs                  ! Number of SF emissions
  INTEGER                                :: StartnumberOfInits               ! 0 if old emit defined (array is copied into 0. entry)
  REAL                                   :: LowVeloThreshold                 ! Threshold value for removal of low velocity particles
  INTEGER                                :: LowVeloCounter                   ! Counter how many low velocity particles were removed
  REAL                                   :: SphericityIC                     ! Particle sphericity
  ! Random particle diameter
  REAL                                   :: PartDiamVarianceIC               ! Variance of random particle diameter
  REAL                                   :: ScalePartDiam                    ! Scale particle diameter if doRandomPart==T
#if USE_SPHERICITY
  ! Random particle sphericity
  REAL                                   :: PartSpheVarianceIC               ! Variance of random sphericity
#endif
  ! Bons particle rebound model
  REAL                                   :: YoungIC                          ! Young's modulus
  REAL                                   :: PoissonIC                        ! Poisson's ration for transverse strain under ax. comp
  REAL                                   :: Whitaker_a                       ! Constant which scales the velocity to obtain the yield coeff
  REAL                                   :: YieldCoeff                       ! Yield strength coefficient
  LOGICAL                                :: doRoughWallModelling             ! Walls are modelled as rough walls
#if USE_EXTEND_RHS
  LOGICAL                                :: CalcSaffmanForce                 ! Calculate the lift (Saffman) force
  LOGICAL                                :: CalcMagnusForce                  ! Calculate the Magnus force
  LOGICAL                                :: CalcUndisturbedFlow              ! Calculate the undisturbed flow force
  LOGICAL                                :: CalcVirtualMass                  ! Calculate the virtual mass force
  LOGICAL                                :: CalcBassetForce                  ! Calculate the Basset force
#endif /* USE_EXTEND_RHS */
END TYPE tSpecies

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
END TYPE tParticleElementMapping

TYPE(tParticleElementMapping)            :: PEM

TYPE tParticleDataManagement
  INTEGER                                :: CurrentNextFreePosition          ! Index of nextfree index in nextFreePosition-Array
  INTEGER                                :: maxParticleNumber                ! Maximum Number of all Particles
  INTEGER                                :: maxAllowedParticleNumber         ! Maximum allowed number of PDM%maxParticleNumber
  LOGICAL                                :: RearrangePartIDs                 ! Rearrange PartIDs during shrinking maxPartNum
  INTEGER                                :: ParticleVecLength=0              ! Vector Length for Particle Push Calculation
  REAL                                   :: MaxPartNumIncrease               ! How much shall the PDM%MaxParticleNumber be incresed if it is full
  INTEGER , ALLOCATABLE                  :: PartInit(:)                      ! (1:NParts), initial emission condition number
                                                                             ! the calculation area
  INTEGER ,ALLOCATABLE                   :: nextFreePosition(:)  ! =>NULL()  ! next_free_Position(1:max_Particle_Number)
                                                                             ! List of free Positon
  LOGICAL ,ALLOCATABLE                   :: ParticleInside(:)    ! =>NULL()  ! Particle_inside(1:Particle_Number)
  LOGICAL ,ALLOCATABLE                   :: IsNewPart(:)                     ! Reconstruct RK-scheme in next stage
END TYPE tParticleDataManagement

TYPE (tParticleDataManagement)           :: PDM

LOGICAL                                  :: ParticlesInitIsDone=.FALSE.

LOGICAL                                  :: DoSurfaceFlux                    ! Flag for emitting by SurfaceFluxBCs
LOGICAL                                  :: UseCircularInflow                !
INTEGER, ALLOCATABLE                     :: CountCircInflowType(:,:,:)
LOGICAL                                  :: DoPoissonRounding                ! Perform Poisson sampling instead of random rounding
LOGICAL                                  :: DoTimeDepInflow                  ! Insertion and SurfaceFlux w simple random rounding
LOGICAL                                  :: RepWarn = .FALSE.                ! Warning for Reynolds limit of particle model

REAL,ALLOCATABLE                         :: gradUx2(:,:,:,:,:,:),gradUy2(:,:,:,:,:,:),gradUz2(:,:,:,:,:,:)
#if USE_FAXEN_CORR
REAL,ALLOCATABLE                         :: gradUx_master_loc(:,:,:,:), gradUx_slave_loc(:,:,:,:)
REAL,ALLOCATABLE                         :: gradUy_master_loc(:,:,:,:), gradUy_slave_loc(:,:,:,:)
REAL,ALLOCATABLE                         :: gradUz_master_loc(:,:,:,:), gradUz_slave_loc(:,:,:,:)
!REAL,ALLOCATABLE                         :: U_local(:,:,:,:,:)
!REAL,ALLOCATABLE                         :: gradp_local(:,:,:,:,:,:)
#endif /* USE_FAXEN_CORR */

#if USE_EXTEND_RHS && ANALYZE_RHS
REAL                                     :: tWriteRHS
REAL                                     :: dtWriteRHS
CHARACTER(LEN=255)                       :: Filename_RHS
REAL,ALLOCATABLE                         :: Pt_ext(:,:)
#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

!===================================================================================================================================
END MODULE MOD_Particle_Vars
