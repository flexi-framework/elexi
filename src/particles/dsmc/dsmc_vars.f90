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
! Contains the DSMC variables
!===================================================================================================================================
MODULE MOD_DSMC_Vars
! MODULES
#if USE_MPI
USE MOD_Particle_MPI_Vars, ONLY: tPartMPIConnect
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: useDSMC

TYPE tDSMC 
  INTEGER                       :: ElectronSpecies          ! Species of the electron
  REAL                          :: EpsElecBin               ! percentage parameter of electronic energy level merging
  REAL                          :: GammaQuant               ! GammaQuant for zero point energy in Evib (perhaps also Erot), 
                                                            ! should be 0.5 or 0
  INTEGER(KIND=8), ALLOCATABLE  :: NumColl(:)               ! Number of Collision for each case + entire Collision number
  REAL                          :: TimeFracSamp=0.          ! %-of simulation time for sampling
  INTEGER                       :: SampNum                  ! number of Samplingsteps
  INTEGER                       :: NumOutput                ! number of Outputs
  REAL                          :: DeltaTimeOutput          ! Time intervall for Output
  LOGICAL                       :: ReservoirSimu            ! Flag for reservoir simulation
  LOGICAL                       :: ReservoirSimuRate        ! Does not performe the collision.
                                                            ! Switch to enable to create reaction rates curves
  LOGICAL                       :: ReservoirSurfaceRate     ! Switch enabling surface rate output without changing surface coverages                                                          
  LOGICAL                       :: ReservoirRateStatistic   ! if false, calculate the reaction coefficient rate by the probability
                                                            ! Default Value is false
  INTEGER                       :: VibEnergyModel           ! Model for vibration Energy: 
                                                            !       0: SHO (default value!)
                                                            !       1: TSHO 
  LOGICAL                       :: DoTEVRRelaxation         ! Flag for T-V-E-R or more simple T-V-R T-E-R relaxation
  INTEGER                       :: PartNumOctreeNode        ! Max Number of Particles per Octree Node
  INTEGER                       :: PartNumOctreeNodeMin     ! Min Number of Particles per Octree Node
  LOGICAL                       :: UseOctree                ! Flag for Octree
  LOGICAL                       :: CalcSurfaceVal           ! Flag for calculation of surfacevalues like heatflux or force at walls
  LOGICAL                       :: CalcSurfaceTime          ! Flag for sampling in time-domain or iterations
  REAL                          :: CalcSurfaceSumTime       ! Flag for sampling in time-domain or iterations
  REAL                          :: CollProbMean             ! Summation of collision probability
  REAL                          :: CollProbMax              ! Maximal collision probability per cell
  INTEGER                       :: CollProbMeanCount        ! counter of possible collision pairs
  INTEGER                       :: CollSepCount             ! counter of actual collision pairs
  REAL                          :: CollSepDist              ! Summation of mean collision separation distance
  LOGICAL                       :: CalcQualityFactors       ! Enables/disables the calculation and output of flow-field variables
  REAL, ALLOCATABLE             :: QualityFactors(:,:)      ! Quality factors for DSMC
                                                            !     1: Maximal collision prob
                                                            !     2: Time-averaged mean collision prob
                                                            !     3: Mean collision separation distance over mean free path
  REAL, ALLOCATABLE             :: QualityFacSamp(:,:)      ! Sampling of quality factors
                                                            !     1: Time-averaged mean collision prob
                                                            !     2: Mean collision separation distance over mean free path
  LOGICAL                       :: ElectronicModel          ! Flag for Electronic State of atoms and molecules
  CHARACTER(LEN=64)             :: ElectronicModelDatabase  ! Name of Electronic State Database | h5 file
  INTEGER                       :: NumPolyatomMolecs        ! Number of polyatomic molecules
  LOGICAL                       :: OutputMeshInit           ! Write Outputmesh (for const. pressure BC) at Init.
  LOGICAL                       :: OutputMeshSamp           ! Write Outputmesh (for const. pressure BC) 
                                                            ! with sampling values at t_analyze
  INTEGER                       :: WallModel                ! Model for wall interaction
                                                            ! 0 perfect/diffusive reflection
                                                            ! 1 adsorption (Kisluik) / desorption (Polanyi Wigner)
                                                            ! 2 adsorption/desorption + chemical interaction (UBI-QEP)
  REAL                          :: RotRelaxProb             ! Model for calculation of rotational relaxation probability, ini_1
                                                            !    0-1: constant probability  (0: no relaxation)
                                                            !    2: Boyd's model
                                                            !    3: Nonequilibrium Direction Dependent model (Zhang,Schwarzentruber)
  REAL                          :: VibRelaxProb             ! Model for calculation of vibrational relaxation probability, ini_1
                                                            !    0-1: constant probability (0: no relaxation)
                                                            !    2: Boyd's model, with correction from Abe
  REAL                          :: ElecRelaxProb            ! electronic relaxation probability
  LOGICAL                       :: PolySingleMode           ! Separate relaxation of each vibrational mode of a polyatomic in a
                                                            ! loop over all vibrational modes (every mode has its own corrected
                                                            ! relaxation probability, comparison with the same random number
                                                            ! while the previous probability is added to the next)
  REAL, ALLOCATABLE             :: InstantTransTemp(:)      ! Instantaneous translational temprerature for each cell (nSpieces+1)
  LOGICAL                       :: BackwardReacRate         ! Enables the automatic calculation of the backward reaction rate
                                                            ! coefficient with the equilibrium constant by partition functions
  REAL                          :: PartitionMaxTemp         ! Temperature limit for pre-stored partition function (DEF: 20 000K)
  REAL                          :: PartitionInterval        ! Temperature interval for pre-stored partition function (DEF: 10K)
END TYPE tDSMC

TYPE(tDSMC)                     :: DSMC

REAL,ALLOCATABLE                :: MacroSurfaceVal(:,:,:,:)      ! variables,p,q,sides
REAL,ALLOCATABLE                :: MacroSurfaceSpecVal(:,:,:,:,:)! Macrovalues for Species specific surface output
                                                                   ! (4,p,q,nSurfSides,nSpecies)
                                                                   ! 1: Surface Collision Counter
                                                                   ! 2: Accomodation
                                                                   ! 3: Coverage
                                                                   ! 4: Recombination Coefficient

INTEGER(KIND=8)                 :: iter_macsurfvalout
!===================================================================================================================================
END MODULE MOD_DSMC_Vars
