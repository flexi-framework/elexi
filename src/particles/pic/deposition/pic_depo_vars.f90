!!=================================================================================================================================
!! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz 
!! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
!! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!!
!! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!!
!! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!!=================================================================================================================================
!
!!===================================================================================================================================
!! Contains the variables for the particle deposition
!!===================================================================================================================================
!MODULE MOD_PICDepo_Vars 
!! MODULES
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!PUBLIC
!SAVE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! GLOBAL VARIABLES 
!!-----------------------------------------------------------------------------------------------------------------------------------
!LOGICAL                               :: DoDeposition       ! flag to switch deposition on/off
!REAL,ALLOCATABLE                      :: PartSource(:,:,:,:,:)  ! PartSource(1:4,PP_N,PP_N,PP_N,nElems) current and charge density
!REAL,ALLOCATABLE                      :: GaussBorder(:)     ! 1D coords of gauss points in -1|1 space
!INTEGER,ALLOCATABLE                   :: GaussBGMIndex(:,:,:,:,:) ! Background mesh index of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
!REAL,ALLOCATABLE                      :: GaussBGMFactor(:,:,:,:,:) ! BGM factor of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
!REAL,ALLOCATABLE                      :: GPWeight(:,:,:,:,:,:,:) ! Weights for splines deposition (check pic_depo for details)
!CHARACTER(LEN=256)                    :: DepositionType     ! Type of Deposition-Method
!INTEGER,ALLOCATABLE                   :: PartToFIBGM(:,:)   ! Mapping form Particle to FIBGM
!REAL,ALLOCATABLE                      :: ElemRadius2_SF(:)  ! elem radius plus radius_sf
!REAL, ALLOCATABLE                     :: BGMSource(:,:,:,:)
!REAL                                  :: r_sf               ! cutoff radius of shape function
!REAL                                  :: r2_sf              ! cutoff radius of shape function * cutoff radius of shape function
!REAL                                  :: r2_sf_inv          ! 1/cutoff radius of shape function * cutoff radius of shape function
!REAL                                  :: w_sf               ! shapefuntion weight
!REAL                                  :: r_sf0              ! minimal shape function radius
!REAL                                  :: r_sf_scale         ! scaling of shape function radius
!REAL                                  :: BetaFac            ! betafactor of shape-function || integral =1
!INTEGER                               :: sf1d_dir           ! direction of 1D shape function 
!INTEGER                               :: NDepo              ! polynomial degree of delta distri
!REAL,ALLOCATABLE                      :: tempcharge(:)      ! temp-charge for epo. kernal
!REAL,ALLOCATABLE                      :: NDepoChooseK(:,:)               ! array n over n
!REAL,ALLOCATABLE                      :: wBaryNDepo(:)      ! barycentric weights for deposition
!REAL,ALLOCATABLE                      :: swGPNDepo(:)       ! integration weights for deposition
!REAL,ALLOCATABLE                      :: XiNDepo(:)         ! gauss position of barycenters
!REAL,ALLOCATABLE                      :: Vdm_NDepo_GaussN(:,:) ! VdM between different polynomial degrees
!REAL,ALLOCATABLE                      :: DDMassinv(:,:,:,:) ! inverse mass-matrix for deposition
!LOGICAL                               :: DeltaDistriChangeBasis   ! Change polynomial degree
!LOGICAL                               :: DoSFEqui           ! use equidistant points for SF
!INTEGER                               :: SfRadiusInt        ! radius integer for cylindrical and spherical shape function
!REAL,ALLOCATABLE                      :: ElemDepo_xGP(:,:,:,:,:)  ! element xGPs for deposition 
!REAL,ALLOCATABLE                      :: Vdm_EquiN_GaussN(:,:)  ! Vdm from equidistant points to Gauss Points
!INTEGER                               :: alpha_sf           ! shapefuntion exponent 
!REAL                                  :: BGMdeltas(3)       ! Backgroundmesh size in x,y,z
!REAL                                  :: FactorBGM(3)       ! Divider for BGM (to allow real numbers)
!REAL                                  :: BGMVolume          ! Volume of a BGM Cell
!INTEGER                               :: BGMminX            ! Local minimum BGM Index in x
!INTEGER                               :: BGMminY            ! Local minimum BGM Index in y
!INTEGER                               :: BGMminZ            ! Local minimum BGM Index in z
!INTEGER                               :: BGMmaxX            ! Local maximum BGM Index in x
!INTEGER                               :: BGMmaxY            ! Local maximum BGM Index in y
!INTEGER                               :: BGMmaxZ            ! Local maximum BGM Index in z
!LOGICAL                               :: Periodic_Depo      ! Flag for periodic treatment for deposition
!INTEGER                               :: DeltaType          ! Flag
!INTEGER                               :: NKnots
!REAL,ALLOCATABLE                      :: Knots(:)
!REAL,ALLOCATABLE                      :: CellVolWeightFac(:)
!INTEGER                               :: VolIntOrder
!REAL,ALLOCATABLE                      :: VolInt_X(:)
!REAL,ALLOCATABLE                      :: VolInt_W(:)
!REAL,ALLOCATABLE                      :: CellVolWeight_Volumes(:,:,:,:)
!INTEGER                               :: NbrOfSFdepoFixes                  !Number of fixes for shape func depo at planar BCs
!REAL    , ALLOCATABLE                 :: SFdepoFixesGeo(:,:,:)             !1:nFixes;1:2(base,normal);1:3(x,y,z) normal outwards!!!
!REAL    , ALLOCATABLE                 :: SFdepoFixesBounds(:,:,:)          !1:nFixes;1:2(min,max);1:3(x,y,z)
!REAL    , ALLOCATABLE                 :: SFdepoFixesChargeMult(:)          !multiplier for mirrored charges (wall: -1.0, sym: 1.0)
!LOGICAL , ALLOCATABLE                 :: SFdepoFixesPartOfLink(:)          !this fix is part of a link
!REAL                                  :: SFdepoFixesEps                    !epsilon for defined planes
!INTEGER                               :: NbrOfSFdepoFixLinks               !Number of linked SFdepoFixes
!INTEGER , ALLOCATABLE                 :: SFdepoFixLinks(:,:)               !1:nLinks;1:3 (2 fixes are linked with each other!)
!                                                                           !             (:,3 is fraction of 180 deg)
!INTEGER                               :: NbrOfSFdepoLayers                 !Number of const. source layer for sf-depo at planar BCs
!LOGICAL                               :: PrintSFDepoWarnings               ! flag to print the warnings
!REAL    , ALLOCATABLE                 :: SFdepoLayersGeo(:,:,:)            !1:nFixes;1:2(base,normal);1:3(x,y,z) normal outwards!!!
!REAL    , ALLOCATABLE                 :: SFdepoLayersBounds(:,:,:)         !1:nFixes;1:2(min,max);1:3(x,y,z)
!LOGICAL , ALLOCATABLE                 :: SFdepoLayersUseFixBounds(:)       !use alls planes of SFdepoFixes as additional bounds?
!CHARACTER(LEN=256),ALLOCATABLE        :: SFdepoLayersSpace(:)              !name of space (cuboid or cylinder)
!REAL    , ALLOCATABLE                 :: SFdepoLayersBaseVector(:,:,:)     !1:nFixes;1:2;1:3(x,y,z)
!INTEGER , ALLOCATABLE                 :: SFdepoLayersSpec(:)               !species of particles for respective layer
!REAL    , ALLOCATABLE                 :: SFdepoLayersPartNum(:)            !number of particles in volume
!REAL    , ALLOCATABLE                 :: SFdepoLayersRadius(:)             !radius for cylinder-space
!LOGICAL                               :: SFResampleAnalyzeSurfCollis
!TYPE tLastAnalyzeSurfCollis
!  INTEGER                             :: PartNumberSamp                    !number of parts from last sampling
!  INTEGER                             :: PartNumberReduced                 !max. allowed number of parts to be saved
!  LOGICAL                             :: ReducePartNumber                  !reduce PartNumberSamp to PartNumberReduced
!  INTEGER                             :: PartNumberDepo                    !number of parts to be inserted in depo
!  REAL, ALLOCATABLE                   :: WallState(:,:)                    !Pos at wall and velocities from last sampling
!  INTEGER, ALLOCATABLE                :: Species(:)                        !Spec of parts
!  REAL                                :: pushTimeStep                      !timestep for (fractional) euler push from wall into ghost domain
!  INTEGER                             :: PartNumThreshold                  !Threshold for checking inserted parts per depo (otherwise abort)
!  REAL                                :: NormVecOfWall(3)                  !normVec for pushTimeStep
!  REAL                                :: Bounds(1:2,1:3)                   !bounds after push for preventing parts outside of...
!                                                                           !...extruded domain 1:2(min,max);1:3(x,y,z)
!  LOGICAL                             :: UseFixBounds                      !use alls planes of SFdepoFixes as additional bounds?
!  LOGICAL                             :: Restart                           !read-in old DSMCSurfCollis-file for restart
!  CHARACTER(LEN=256)                  :: DSMCSurfCollisRestartFile
!  INTEGER                             :: NumberOfBCs            ! Nbr of BC to be analyzed (def.: 1)
!  INTEGER, ALLOCATABLE                :: BCs(:)                 ! BCs to be analyzed (def.: 0 = all)
!END TYPE
!TYPE(tLastAnalyzeSurfCollis)          :: LastAnalyzeSurfCollis
!!REAL,ALLOCATABLE                      :: Vdm_BernSteinN_GaussN(:,:)
!!REAL,ALLOCATABLE                      :: sVdm_BernSteinN_GaussN(:,:)
!!===================================================================================================================================
!END MODULE MOD_PICDepo_Vars
