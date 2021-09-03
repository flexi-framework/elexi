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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Vars
! MODULES
#if USE_MPI
USE MPI
#endif /*USE_MPI*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

LOGICAL                                 :: SurfOnNode
INTEGER                                 :: SurfSampSize                  !> Impact number + Energy + Angle + Force
REAL,ALLOCPOINT,DIMENSION(:,:,:)        :: SurfSideArea                  !> Area of supersampled surface side
! ====================================================================
! Mesh info
INTEGER                                 :: nSurfTotalSides

INTEGER                                 :: nComputeNodeSurfSides         !> Number of surface sampling sides on compute node
INTEGER                                 :: nComputeNodeSurfTotalSides    !> Number of surface sampling sides on compute node (including halo region)
INTEGER                                 :: offsetComputeNodeSurfSide     !> elem offset of compute-node root

! ====================================================================
! Impact statistics
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: SampWallState                 !> 1     Impact Counter
                                                                         !> 2-5   E_kin (mean, min, max, M2, variance)
                                                                         !> 6-9   Impact Angle (mean, min, max, deltaE)
                                                                         !> 10-12 Current Forces in x, y, z direction
                                                                         !> 13-15 Average Forces in x, y, z direction
                                                                         !> 10    Impact angle
                                                                         !> 11    Wall-Collision counter
REAL,ALLOCPOINT,DIMENSION(:,:,:,:)      :: SampWallState_Shared

! ====================================================================
! Impact tracking
LOGICAL                                 :: ImpactTrackInitIsDone  = .FALSE. !< mark if impact tracking init routine is finished
LOGICAL                                 :: doParticleImpactTrack  = .FALSE. !< mark if impacts tracking should be performed during computation
LOGICAL                                 :: ImpactSideOnProc       = .FALSE. !< marks if current proc has impact tracking side
INTEGER                                 :: PartStateBoundaryVecLength    !< Impacts on current proc
REAL,ALLOCATABLE,DIMENSION(:,:)         :: PartStateBoundary             !< solution evaluated at EPs (nvar,nEP,nSamples)
INTEGER                                 :: ImpactDataSize                !< number of variables stored per impact
INTEGER                                 :: ImpactnGlob                   !< Global number of occured impacts
INTEGER                                 :: ImpactnLoc                    !< Local  number of occured impacts
INTEGER                                 :: ImpactOffset                  !< Offset number of occured impacts

!----------------------------------------------------------------------------------------------------------------------------------
! MPI Communicator for EPs
!----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
INTEGER                                 :: myImpactRank                  !< rank within impact tracking communicator
INTEGER                                 :: MPI_COMM_IMPACT=MPI_COMM_NULL !< MPI impact tracking communicator
INTEGER                                 :: nImpactProcs                  !< number of procs with impact tracking
#endif /* USE_MPI */

! ====================================================================
! MPI3 shared variables
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: GlobalSide2SurfSide           ! Mapping Global Side ID to Surf Side ID
                                                                         !> 1 - Surf SideID
                                                                         !> 2 - Surf Side proc global rank
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: SurfSide2GlobalSide           ! Inverse mapping
                                                                         !> 1 - Surf SideID
                                                                         !> 2 - Surf Side proc global rank
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: GlobalSide2SurfSide_Shared
INTEGER,ALLOCPOINT,DIMENSION(:,:)       :: SurfSide2GlobalSide_Shared

#if USE_MPI
REAL,POINTER,DIMENSION(:,:,:)           :: SurfSideArea_Shared           !> Area of supersampled surface side
INTEGER                                 :: SurfSideArea_Shared_Win

!INTEGER,ALLOCATABLE,DIMENSION(:,:)      :: GlobalSide2SurfHaloSide       ! Mapping Global Side ID to Surf Halo Side ID (exists only on leader procs)
!                                                                         !> 1st dim: leader rank
!                                                                         !> 2nd dim: Surf SideID
!INTEGER,ALLOCATABLE,DIMENSION(:,:)      :: SurfHaloSide2GlobalSide       ! Inverse mapping  (exists only on leader procs)
!                                                                         !> 1st dim: leader rank
!                                                                         !> 2nd dim: Surf SideID

INTEGER                                 :: GlobalSide2SurfSide_Shared_Win
INTEGER                                 :: SurfSide2GlobalSide_Shared_Win

TYPE tSurfaceMapping
  INTEGER,ALLOCATABLE                   :: RecvSurfGlobalID(:)
  INTEGER,ALLOCATABLE                   :: SendSurfGlobalID(:)
  INTEGER                               :: nSendSurfSides
  INTEGER                               :: nRecvSurfSides
END TYPE
TYPE (tSurfaceMapping),ALLOCATABLE      :: SurfMapping(:)

INTEGER                                 :: SampWallState_Shared_Win
#else /*USE_MPI*/
INTEGER                                 :: mySurfRank
#endif /* USE_MPI */

LOGICAL                                 :: doParticleReflectionTrack = .TRUE.      ! Flag if reflections should be counted
LOGICAL                                 :: doParticleImpactSample                  ! Flag if impact data should be tracked
LOGICAL                                 :: WriteMacroSurfaceValues   = .FALSE.     ! Output of macroscopic values on surface
REAL                                    :: MacroValSampTime                        ! Sampling time for WriteMacroVal
LOGICAL                                 :: ImpactRestart                           ! Flag if we are restarting surface sampling

INTEGER                                 :: nImpactVars                             ! Number of Vars = nImpactVars * (nSpecies + 1)

REAL,ALLOCATABLE                        :: MacroSurfaceVal(:,:,:,:)                ! variables,p,q,sides
REAL,ALLOCATABLE                        :: MacroSurfaceSpecVal(:,:,:,:,:)          ! Macrovalues for Species specific surface output


!-----------------------------------------------------------------------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: NSurfSample                   ! polynomial degree of particle BC sampling
CHARACTER(LEN=255),ALLOCATABLE          :: SurfSampleBCs(:)              ! name of additional surface sampling BCs
REAL,ALLOCATABLE                        :: XiEQ_SurfSample(:)            ! position of XiEQ_SurfSample
REAL                                    :: dXiEQ_SurfSample              ! deltaXi in [-1,1]
INTEGER                                 :: OffSetSurfSide                ! offset of local surf side
INTEGER                                 :: nSurfBC                       ! number of surface side BCs
CHARACTER(LEN=255),ALLOCATABLE          :: SurfBCName(:)                 ! names of belonging surface BC

! Boundary
TYPE tPartBoundary
  INTEGER                                :: InternalBC              = 0      !
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryBC              = 10     ! = 10 (s.u.) Boundary Condition Integer Definition
!  INTEGER                                :: AnalyzeBC               = 100    ! = 100 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)    ! Link part 1 for mapping Boltzplatz BCs to Particle BC
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundType(:)    ! Link part 2 for mapping Boltzplatz BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)    ! Link part 3 for mapping Boltzplatz BCs to Particle BC
  REAL    , ALLOCATABLE                  :: WallTemp(:)
  REAL    , ALLOCATABLE                  :: WallVelo(:,:)
!  LOGICAL , ALLOCATABLE                  :: AmbientCondition(:)
!  LOGICAL , ALLOCATABLE                  :: AmbientConditionFix(:)
!  REAL    , ALLOCATABLE                  :: AmbientTemp(:)
!  REAL    , ALLOCATABLE                  :: AmbientVelo(:,:)
!  REAL    , ALLOCATABLE                  :: AmbientDens(:)
!  REAL    , ALLOCATABLE                  :: AmbientDynamicVisc(:) ! dynamic viscousity
  CHARACTER(LEN=255),ALLOCATABLE         :: WallModel(:)
  CHARACTER(LEN=255),ALLOCATABLE         :: WallCoeffModel(:)
  ! Bons particle rebound model
  REAL    , ALLOCATABLE                  :: Young(:)                ! Young's modulus
  REAL    , ALLOCATABLE                  :: Poisson(:)              ! Poisson's ration for transverse strain under axial compression
  ! Fong coefficient of restitution
  REAL    , ALLOCATABLE                  :: CoR(:)                  ! Coefficient of restitution for normal velocity component
  ! Rough wall modelling
  LOGICAL , ALLOCATABLE                  :: doRoughWallModelling(:) ! Flag if walls are modelled as rough walls
  REAL    , ALLOCATABLE                  :: RoughMeanIC(:)          ! Mean of the Gaussian distribution
  REAL    , ALLOCATABLE                  :: RoughVarianceIC(:)      ! Standard deviation of the Gaussian distribution
END TYPE

INTEGER                                  :: nPartBound              ! number of particle boundaries
TYPE(tPartBoundary)                      :: PartBound               ! Boundary Data for Particles

INTEGER                                  :: nAuxBCs                 ! number of aux. BCs that are checked during tracing
LOGICAL                                  :: UseAuxBCs               ! number of aux. BCs that are checked during tracing
CHARACTER(LEN=200), ALLOCATABLE          :: AuxBCType(:)            ! type of BC (plane, ...)
INTEGER           , ALLOCATABLE          :: AuxBCMap(:)             ! index of AuxBC in respective Type

TYPE tAuxBC_plane
  REAL                                   :: r_vec(3)
  REAL                                   :: n_vec(3)
  REAL                                   :: radius
END TYPE tAuxBC_plane
TYPE(tAuxBC_plane), ALLOCATABLE          :: AuxBC_plane(:)

TYPE tAuxBC_cylinder
  REAL                                   :: r_vec(3)
  REAL                                   :: axis(3)
  REAL                                   :: radius
  REAL                                   :: lmin
  REAL                                   :: lmax
  LOGICAL                                :: inwards
END TYPE tAuxBC_cylinder
TYPE(tAuxBC_cylinder), ALLOCATABLE       :: AuxBC_cylinder(:)

TYPE tAuxBC_cone
  REAL                                   :: r_vec(3)
  REAL                                   :: axis(3)
  REAL                                   :: halfangle
  REAL                                   :: lmin
  REAL                                   :: lmax
  REAL                                   :: geomatrix(3,3)
  REAL                                   :: rotmatrix(3,3)
  LOGICAL                                :: inwards
END TYPE tAuxBC_cone
TYPE(tAuxBC_cone), ALLOCATABLE       :: AuxBC_cone(:)

TYPE tAuxBC_parabol
  REAL                                   :: r_vec(3)
  REAL                                   :: axis(3)
  REAL                                   :: zfac
  REAL                                   :: lmin
  REAL                                   :: lmax
  REAL                                   :: geomatrix4(4,4)
  REAL                                   :: rotmatrix(3,3)
  LOGICAL                                :: inwards
END TYPE tAuxBC_parabol
TYPE(tAuxBC_parabol), ALLOCATABLE       :: AuxBC_parabol(:)

TYPE tPartAuxBC
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER , ALLOCATABLE                  :: TargetBoundCond(:)
  REAL    , ALLOCATABLE                  :: MomentumACC(:)
  REAL    , ALLOCATABLE                  :: WallTemp(:)
  REAL    , ALLOCATABLE                  :: WallVelo(:,:)
END TYPE
TYPE(tPartAuxBC)                         :: PartAuxBC                         ! auxBC Data for Particles

LOGICAL                                  :: LowVeloRemove                 !Flag if low velocity particles should be removed
!===================================================================================================================================

END MODULE MOD_Particle_Boundary_Vars
