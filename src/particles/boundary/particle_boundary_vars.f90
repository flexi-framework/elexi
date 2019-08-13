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

!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: NSurfSample                   ! polynomial degree of particle BC sampling   
REAL,ALLOCATABLE                        :: XiEQ_SurfSample(:)            ! position of XiEQ_SurfSample
REAL                                    :: dXiEQ_SurfSample              ! deltaXi in [-1,1]
INTEGER                                 :: OffSetSurfSide                ! offset of local surf side
INTEGER                                 :: nSurfBC                       ! number of surface side BCs
CHARACTER(LEN=255),ALLOCATABLE          :: SurfBCName(:)                 ! names of belonging surface BC
#if USE_MPI
INTEGER,ALLOCATABLE                     :: OffSetSurfSideMPI(:)          ! integer offset for particle boundary sampling            
#endif /*MPI*/

#if USE_MPI
TYPE tSurfaceSendList
  INTEGER                               :: NativeProcID
  INTEGER,ALLOCATABLE                   :: SendList(:)                   ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: RecvList(:)                   ! list containing surfsideid of sides to recv from proc
  
  INTEGER,ALLOCATABLE                   :: SurfDistSendList(:)           ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: SurfDistRecvList(:)           ! list containing surfsideid of sides to recv from proc
  INTEGER,ALLOCATABLE                   :: CoverageSendList(:)           ! list containing surfsideid of sides to send to proc
  INTEGER,ALLOCATABLE                   :: CoverageRecvList(:)           ! list containing surfsideid of sides to recv from proc
  
END TYPE
#endif /*MPI*/

TYPE tSurfaceCOMM
  LOGICAL                               :: MPIRoot                       ! if root of mpi communicator
  INTEGER                               :: MyRank                        ! local rank in new group
  INTEGER                               :: nProcs                        ! number of processes
  LOGICAL                               :: MPIOutputRoot                 ! if root of mpi communicator
  INTEGER                               :: MyOutputRank                  ! local rank in new group
  INTEGER                               :: nOutputProcs                  ! number of output processes
#if USE_MPI
  INTEGER                               :: COMM                          ! communicator
  INTEGER                               :: nMPINeighbors                 ! number of processes to communicate with
  TYPE(tSurfaceSendList),ALLOCATABLE    :: MPINeighbor(:)                ! list containing all mpi neighbors
  INTEGER                               :: OutputCOMM                    ! communicator for output
#endif /*MPI*/
END TYPE
TYPE (tSurfaceCOMM)                     :: SurfCOMM

TYPE tSurfaceMesh
  INTEGER                               :: SampSize                      ! integer of sampsize
  LOGICAL                               :: SurfOnProc                    ! flag if reflective boundary condition is on proc
  INTEGER                               :: nSides                        ! Number of Sides on Surface (reflective)
  INTEGER                               :: nTotalSides                   ! Number of Sides on Surface incl. HALO sides
  INTEGER                               :: nGlobalSides                  ! Global number of Sides on Surfaces (reflective)
  INTEGER,ALLOCATABLE                   :: SideIDToSurfID(:)             ! Mapping form the SideID to shorter side list
  REAL, ALLOCATABLE                     :: SurfaceArea(:,:,:)            ! Area of Surface 
  INTEGER,ALLOCATABLE                   :: SurfSideToGlobSideMap(:)      ! map of surfside ID to global Side ID
END TYPE

TYPE (tSurfaceMesh)                     :: SurfMesh

TYPE tSampWall             ! DSMC sample for Wall                                             
  ! easier to communicate
  ! Data structure is repeated for every species + average
  REAL,ALLOCATABLE                      :: State(:,:,:)                ! 1     Impact Counter
                                                                       ! 2-5   E_kin (mean, min, max, M2, variance)
                                                                       ! 6-9   Impact Angle (mean, min, max, deltaE)
                                                                       ! 10-12 Current Forces in x, y, z direction
                                                                       ! 13-15 Average Forces in x, y, z direction
                                                                       ! 10    Impact angle
                                                                       ! 11 Wall-Collision counter
  !REAL, ALLOCATABLE                    :: Force(:,:,:)                ! x, y, z direction
  !REAL, ALLOCATABLE                    :: Counter(:,:,:)              ! Wall-Collision counter
END TYPE
TYPE(tSampWall), ALLOCATABLE            :: SampWall(:)             ! Wall sample array (number of BC-Sides)

TYPE tSurfColl
  INTEGER                               :: NbrOfSpecies           ! Nbr. of Species to be counted for wall collisions (def. 0: all)
  LOGICAL,ALLOCATABLE                   :: SpeciesFlags(:)        ! Species counted for wall collisions (def.: all species=T)
  LOGICAL                               :: OnlySwaps              ! count only wall collisions being SpeciesSwaps (def. F)
  LOGICAL                               :: Only0Swaps             ! count only wall collisions being delete-SpeciesSwaps (def. F)
  LOGICAL                               :: Output                 ! Print sums of all counted wall collisions (def. F)
  LOGICAL                               :: AnalyzeSurfCollis      ! Output of collided/swaped particles 
                                                                  ! during Sampling period? (def. F)
END TYPE
TYPE (tSurfColl)                        :: CalcSurfCollis
  
TYPE tAnalyzeSurfCollis 
  INTEGER                               :: maxPartNumber          ! max. number of collided/swaped particles during Sampling
  REAL, ALLOCATABLE                     :: Data(:,:)              ! Output of collided/swaped particles during Sampling period
                                                                  ! (Species,Particles,Data(x,y,z,u,v,w)
  INTEGER, ALLOCATABLE                  :: Spec(:)                ! Species of Particle in Data-array
  INTEGER, ALLOCATABLE                  :: BCid(:)                ! ID of PartBC from crossing of Particle in Data-array
  INTEGER, ALLOCATABLE                  :: Number(:)              ! collided/swaped particles per Species during Sampling period
  !REAL, ALLOCATABLE                     :: Rate(:)                ! collided/swaped particles/s per Species during Sampling period
  INTEGER                               :: NumberOfBCs            ! Nbr of BC to be analyzed (def.: 1)
  INTEGER, ALLOCATABLE                  :: BCs(:)                 ! BCs to be analyzed (def.: 0 = all)

END TYPE tAnalyzeSurfCollis
TYPE(tAnalyzeSurfCollis)                :: AnalyzeSurfCollis

TYPE tPartBoundary
  INTEGER                                :: OpenBC                  = 1      ! = 1 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: ReflectiveBC            = 2      ! = 2 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: PeriodicBC              = 3      ! = 3 (s.u.) Boundary Condition Integer Definition
  INTEGER                                :: SymmetryBC              = 10     ! = 10 (s.u.) Boundary Condition Integer Definition
!  INTEGER                                :: AnalyzeBC               = 100    ! = 100 (s.u.) Boundary Condition Integer Definition
  CHARACTER(LEN=200)   , ALLOCATABLE     :: SourceBoundName(:)    ! Link part 1 for mapping Boltzplatz BCs to Particle BC
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)    ! Link part 2 for mapping Boltzplatz BCs to Particle BC
!  INTEGER              , ALLOCATABLE     :: Map(:)               ! Map from Boltzplatz BCindex to Particle BC
  INTEGER              , ALLOCATABLE     :: MapToPartBC(:)        ! Map from Boltzplatz BCindex to Particle BC (NOT TO TYPE!)
  REAL    , ALLOCATABLE                  :: MomentumACC(:)      
  REAL    , ALLOCATABLE                  :: WallTemp(:)     
  REAL    , ALLOCATABLE                  :: WallVelo(:,:) 
  LOGICAL , ALLOCATABLE                  :: SolidState(:)         ! flag defining if reflective BC is solid or liquid
  INTEGER , ALLOCATABLE                  :: SolidSpec(:)
  REAL    , ALLOCATABLE                  :: SolidPartDens(:)
  INTEGER , ALLOCATABLE                  :: LiquidSpec(:)         ! Species of Liquid Boundary
  LOGICAL , ALLOCATABLE                  :: AmbientCondition(:)
  LOGICAL , ALLOCATABLE                  :: AmbientConditionFix(:)
  REAL    , ALLOCATABLE                  :: AmbientTemp(:)
  REAL    , ALLOCATABLE                  :: AmbientMeanPartMass(:)
  REAL    , ALLOCATABLE                  :: AmbientBeta(:)
  REAL    , ALLOCATABLE                  :: AmbientVelo(:,:)
  REAL    , ALLOCATABLE                  :: AmbientDens(:)
  REAL    , ALLOCATABLE                  :: AmbientDynamicVisc(:) ! dynamic viscousity
  LOGICAL , ALLOCATABLE                  :: Adaptive(:)
  INTEGER , ALLOCATABLE                  :: AdaptiveType(:)
  INTEGER , ALLOCATABLE                  :: AdaptiveMacroRestartFileID(:)
  REAL    , ALLOCATABLE                  :: AdaptivePressure(:)
  REAL    , ALLOCATABLE                  :: AdaptiveTemp(:)
!  LOGICAL , ALLOCATABLE                  :: UseForQCrit(:)       ! Use Boundary for Q-Criterion ?
  LOGICAL , ALLOCATABLE                  :: Resample(:)           ! Resample Equilibirum Distribution with reflection
  CHARACTER(LEN=255),ALLOCATABLE         :: WallModel(:)
  CHARACTER(LEN=255),ALLOCATABLE         :: WallCoeffModel(:)
  ! Bons particle rebound model
  REAL    , ALLOCATABLE                  :: Young(:)              ! Young's modulus
  REAL    , ALLOCATABLE                  :: Poisson(:)            ! Poisson's ration for transverse strain under axial compression
END TYPE

INTEGER                                  :: nPartBound            ! number of particle boundaries
INTEGER                                  :: nAdaptiveBC
TYPE(tPartBoundary)                      :: PartBound             ! Boundary Data for Particles

INTEGER                                  :: nAuxBCs               ! number of aux. BCs that are checked during tracing
LOGICAL                                  :: UseAuxBCs             ! number of aux. BCs that are checked during tracing
CHARACTER(LEN=200), ALLOCATABLE          :: AuxBCType(:)          ! type of BC (plane, ...)
INTEGER           , ALLOCATABLE          :: AuxBCMap(:)           ! index of AuxBC in respective Type

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
  !REAL                                   :: geomatrix2(3,3)
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
  INTEGER              , ALLOCATABLE     :: TargetBoundCond(:)
  REAL    , ALLOCATABLE                  :: MomentumACC(:)      
  REAL    , ALLOCATABLE                  :: WallTemp(:)     
  REAL    , ALLOCATABLE                  :: WallVelo(:,:) 
  INTEGER , ALLOCATABLE                  :: NbrOfSpeciesSwaps(:)          !Number of Species to be changed at wall
  REAL    , ALLOCATABLE                  :: ProbOfSpeciesSwaps(:)         !Probability of SpeciesSwaps at wall
  INTEGER , ALLOCATABLE                  :: SpeciesSwaps(:,:,:)           !Species to be changed at wall (in, out), out=0: delete
  LOGICAL , ALLOCATABLE                  :: Resample(:)                      !Resample Equilibirum Distribution with reflection
END TYPE
TYPE(tPartAuxBC)                         :: PartAuxBC                         ! auxBC Data for Particles

LOGICAL                                  :: LowVeloRemove                 !Flag if low velocity particles should be removed
!===================================================================================================================================

END MODULE MOD_Particle_Boundary_Vars
