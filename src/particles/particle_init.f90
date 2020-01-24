!=================================================================================================================================
! Copyright (c) 2010-2020  Prof. Claus-Dieter Munz
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
! Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_ParticleInit
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersParticles
  MODULE PROCEDURE DefineParametersParticles
END INTERFACE

INTERFACE InitParticles
  MODULE PROCEDURE InitParticles
END INTERFACE

INTERFACE FinalizeParticles
  MODULE PROCEDURE FinalizeParticles
END INTERFACE

PUBLIC :: DefineParametersParticles
PUBLIC :: InitParticles
PUBLIC :: FinalizeParticles
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particles
!==================================================================================================================================
SUBROUTINE DefineParametersParticles()
! MODULES
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE Mod_Particle_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection('Particle')

CALL prms%CreateIntOption('BezierClipLineVectorMethod',     '' ,                                            '2')
CALL prms%CreateIntOption('NbrOfRegions',                   'Number of regions to be mapped to Elements',   '0')
CALL prms%CreateIntOption('Part-nAuxBCs',                   'Number of auxillary BCs that are checked during tracing',  '0')
!CALL prms%CreateIntOption('Part-nBounds',                   'Number of particle boundaries.',               '1')
CALL prms%CreateIntOption('Part-nPeriodicVectors',          'Number of the periodic vectors j=1,...,n.'                          //&
                                                            ' Value has to be the same as defined in preprog.ini',      '0')
CALL prms%CreateIntOption('Part-nSpecies',                  'Number of species in part',                    '1')
CALL prms%CreateIntOption('Part-maxParticleNumber',         'Max number of particles in part',              '1')
CALL prms%CreateIntOption('Part-NumberOfRandomSeeds',       'Number of random seeds for particle random number generator','0')

CALL prms%CreateLogicalOption('Part-SteadyState',           'Only update particle position while keeping fluid state frozen',      &
                                                            '.FALSE.')
CALL prms%CreateLogicalOption('Part-TrackPosition',         'Track particle position',                      '.FALSE.')
CALL prms%CreateLogicalOption('Part-WriteMacroSurfaceValues','Set [T] to activate iteration dependant sampling and h5 output'    //&
                                                            ' surfaces.',                                                          &
                                                                                                            '.FALSE.')
CALL prms%CreateLogicalOption('Part-KeepWallParticles',     'Flag if particles stuck to walls through fouling should be removed',  &
                                                                                                            '.FALSE.')
CALL prms%CreateLogicalOption('printRandomSeeds',           'Flag defining if random seeds are written.',   '.FALSE.')
CALL prms%CreateLogicalOption('Part-AllowLoosing',          'Flag if a lost particle should abort the programm','.FALSE.')
CALL prms%CreateLogicalOption('Part-LowVeloRemove',         'Flag if low velocity particles should be removed', '.FALSE.')

CALL prms%CreateRealOption('Part-DelayTime',                "During delay time the particles won't be moved so "                 //&
                                                            "the fluid field can be evolved",               '0.')
CALL prms%CreateRealOption('Part-SafetyFactor',             'Factor to scale the halo region with MPI',     '1.')
CALL prms%CreateRealOption('Part-SteadyTimeStep',           'Manual time step routine for frozen fluid state', '0.')
CALL prms%CreateRealOption('Particles-HaloEpsVelo',         'Maximum velocity to be considered for halo region', '0.')
CALL prms%CreateRealOption('Particles-ManualTimeStep',      'Manual time step routine for frozen fluid state', '0.')

CALL prms%CreateRealArrayOption('Part-Gravity',             'Gravitational acceleration as vector',         '0. , 0. , 0.')

#if USE_RW
CALL prms%SetSection("Particle Random Walk")
!===================================================================================================================================
! >>> Values in this section only apply for turbulence models providing turbulent kinetic energy and a turbulent length/time scale
!===================================================================================================================================
CALL prms%CreateStringOption(   'Part-RWModel' &
                                , 'Random walk model used for steady-state calculations.\n'//&
                                ' - Gosman \n'//&
                                ' - Dehbi \n'//&
                                ' - Langevin \n'&
                                , 'none')
CALL prms%CreateStringOption(   'Part-RWTime'  &
                                , 'Time stepping used for random walk model.\n'//&
                                ' - RK \n'//&
                                ' - RW'&
                                , 'RW')
#endif /* USE_RW */

!===================================================================================================================================
! >>> Options for particle SGS model
!===================================================================================================================================
CALL prms%CreateStringOption(   'Part-SGSModel' &
                                , 'SGS model used for reconstruction of SGS influence on particle\n'//&
                                ' - Breuer \n'//&
                                ' - none'&
                                , 'none')
CALL prms%CreateIntOption(      'Part-SGSNFilter' &
                                , 'Number of cut-off modes in the high-pass SGS filter','2')

!===================================================================================================================================
! > Species
! >>> Values in this section appear multiple times
!===================================================================================================================================
CALL prms%SetSection("Particle Species")
! species inits
CALL prms%CreateIntOption(      'Part-Species[$]-nInits'  &
                                , 'Number of different initial particle placements for Species [$]', '0', numberedmulti=.TRUE.)

CALL prms%CreateLogicalOption(  'Part-Species[$]-CalcHeightFromDt'  &
                                , 'Calculated cuboid/cylinder height from v and dt?'&
                                , '.FALSE.', numberedmulti=.TRUE.)

CALL prms%CreateStringOption(   'Part-Species[$]-RHSMethod' &
                                , 'Particle model used for forces calculation.\n'//&
                                ' - Wang \n'//&
                                ' - Jacobs \n'//&
                                ' - Vinkovic \n'&
                                , 'none' , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Species[$]-SpaceIC'  &
                                , 'Specifying Keyword for particle space condition of species [$] in case of one init.\n'//&
                                ' - point \n'//&
                                ' - line_with_equidistant_distribution \n'//&
                                ' - line \n'//&
                                ' - disc \n'//&
                                ' - circle_equidistant \n'//&
                                ' - cuboid \n'//&
                                ' - cylinder \n'&
                              , 'cuboid', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Species[$]-velocityDistribution'  &
                                , 'Used velocity distribution.\n'//&
                                  '   constant: all particles have the same velocity defined in VeloVecIC\n'//&
                                  '   fluid:    particles have local fluid velocity\n'&
                                  , 'constant', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-CuboidHeightIC'  &
                                , 'Height of cuboid if SpaceIC=cuboid', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-CylinderHeightIC'  &
                                , 'Third measure of cylinder  (set 0 for flat rectangle),'//&
                                  ' negative value = opposite direction', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-InflowRiseTime' &
                                , 'Time to ramp the number of inflow particles linearly from zero to unity'&
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-ParticleEmission' &
                                , 'Emission rate in part/s or part/iteration.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MacroParticleFactor' &
                                , 'Number of Microparticle per Macroparticle for species [$]', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MassIC'  &
                                , 'Particle Mass of species [$] [kg]', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-DensityIC'  &
                                , 'Particle density of species [$] [kg/m^3]', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-RadiusIC'  &
                                , 'Radius for IC circle', '1.', numberedmulti=.TRUE.)
!CALL prms%CreateRealOption(     'Part-Species[$]-VeloIC'  &
!                                , 'Absolute value of initial velocity. (ensemble velocity) ', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-VeloTurbIC'  &
                                , 'Turbulent fluctuation of initial velocity. (ensemble velocity) ', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-LowVeloThreshold' &
                                , 'Threshold velocity of particles after reflection. Slower particles are deleted [$] [m/s]' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-HighVeloThreshold' &
                                , 'Threshold velocity of particles in the entire field. Faster particles are deleted [$] [m/s]' &
                                , '0.', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Part-Species[$]-initialParticleNumber'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Initial particle number', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-NumberOfExcludeRegions'  &
                                , 'Number of different regions to be excluded', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-ParticleEmissionType'  &
                                , 'Define Emission Type for particles (volume emission)\n'//&
                                  '1 = emission rate in part/s,\n'//&
                                  '2 = emission rate part/iteration\n'&
                                  ,'2', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption('Part-Species[$]-BasePointIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Base point for IC cuboid and IC sphere', '0. , 0. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-BaseVector1IC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'First base vector for IC cuboid', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-BaseVector2IC'  &
                                , 'Second base vector for IC cuboid', '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-BaseVariance'  &
                                , 'Variance for Gaussian distribtution','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-NormalIC'  &
                                , 'Normal orientation of circle.', '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-VeloVecIC '  &
                                , 'Velocity vector for given species', '0. , 0. , 0.', numberedmulti=.TRUE.)

CALL prms%SetSection("Particle Species Ninits")
CALL prms%CreateLogicalOption(  'Part-Species[$]-UseForEmission' &
                                , 'Flag to use Init/Emission for emission', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-UseForInit' &
                                , 'Flag to use Init/Emission for init', '.TRUE.', numberedmulti=.TRUE.)

! Bons particle rebound model
CALL prms%SetSection("Particle Rebound Model")
CALL prms%CreateRealOption(     'Part-Species[$]-YoungIC'       &
                                , "Young's modulus defining stiffness of particle material",'0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PoissonIC'     &
                                , "Poisson ratio defining relation of transverse to axial strain",'0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-YieldCoeff'    &
                                , "Yield strength defining elastic deformation",'0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-Young'  &
                                , "Young's modulus defining stiffness of wall material", '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-Poisson'  &
                                , "Poisson ratio defining relation of transverse to axial strain", '0.', numberedmulti=.TRUE.)


! Fong coefficient of restitution
CALL prms%CreateRealOption(     'Part-Boundary[$]-CoR'  &
                                , "Coefficent of restitution for normal velocity component", numberedmulti=.TRUE.)

!===================================================================================================================================
! > Boundaries
! >>> Values in this section appear multiple times
!===================================================================================================================================
CALL prms%SetSection("Particle Boundaries")

CALL prms%CreateStringOption(   'Part-Boundary[$]-Type'  &
                                , 'Used boundary condition for boundary.\n'//&
                                  '- open\n'//&
                                  '- reflective\n'//&
                                  '- periodic\n'//&
                                 'If condition=open, the following parameters are'//&
                                  ' used: (Part-Boundary-=PB) PB-Ambient ,PB-AmbientTemp,PB-AmbientMeanPartMass,'//&
                                  'PB-AmbientVelo,PB-AmbientDens,PB-AmbientDynamicVisc,PB-AmbientThermalCond,PB-Voltage\n'//&
                                 'If condition=reflective: PB-MomentumACC,PB-WallTemp,PB-TransACC,PB-VibACC,PB-RotACC,'//&
                                  'PB-WallVelo,Voltage,SpeciesSwaps.If condition=periodic:Part-nPeriodicVectors,'//&
                                  'Part-PeriodicVector', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Boundary[$]-Name'  &
                                , 'Source Name of Boundary. Has to be same name as defined in preproc tool', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption('Part-PeriodicVector[$]'      , 'TODO-DEFINE-PARAMETER\nVector for periodic boundaries.'//&
                                                                   'Has to be the same as defined in preproc.ini in their'//&
                                                                   ' respective order. ', '1.,0.,0.', numberedmulti=.TRUE.)

! Wall model =======================================================================================================================
CALL prms%CreateStringOption(   'Part-Boundary[$]-WallModel' &
                                , 'Wall model to be used. Options:.\n'//&
                                  'perfRef  - perfect reflection\n'//&
                                  'coeffRes - Coefficient of restitution',numberedmulti=.TRUE.)
CALL prms%CreateStringOption('Part-Boundary[$]-WallCoeffModel' &
                                , 'Coefficients to be used. Options:.\n'//&
                                  'Tabakoff1981 \n'//&
                                  'Bons2017 \n'    //&
                                  'Fong2019',numberedmulti=.TRUE.)

! Ambient condition=================================================================================================================
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-AmbientCondition'  &
                                , 'Use ambient condition (condition "behind" boundary).', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-AmbientConditionFix'  &
                                , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption('Part-Boundary[$]-AmbientVelo'  &
                                , 'Ambient velocity', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-WallVelo'  &
                                , 'Velocity (global x,y,z in [m/s]) of reflective particle boundary.' &
                                , numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Boundary[$]-AmbientDens', 'Ambient density', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-AmbientDynamicVisc'  &
                                , 'Ambient dynamic viscosity', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Particles-RandomSeed[$]'     , 'Seed [$] for Random Number Generator', '1', numberedmulti=.TRUE.)

END SUBROUTINE DefineParametersParticles

!===================================================================================================================================
! Glue Subroutine for particle initialization
!===================================================================================================================================
SUBROUTINE InitParticles(ManualTimeStep_opt)
! MODULES
USE MOD_Globals
USE Mod_Particle_Globals
USE MOD_ReadInTools
USE MOD_IO_HDF5,                    ONLY: AddToElemData
USE MOD_Part_Emission,              ONLY: InitializeParticleEmission
USE MOD_Particle_Boundary_Sampling, ONLY: InitParticleBoundarySampling
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Restart,           ONLY: ParticleRestart
USE MOD_Particle_SGS,               ONLY: ParticleSGS
USE MOD_Particle_Vars,              ONLY: ParticlesInitIsDone,WriteMacroSurfaceValues
#if USE_MPI
USE MOD_Particle_MPI,               ONLY: InitParticleCommSize
#endif
#if USE_RW
USE MOD_Particle_RandomWalk,        ONLY: ParticleInitRandomWalk
#endif
USE MOD_Particle_SGS,               ONLY: ParticleInitSGS
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL       :: ManualTimeStep_opt                                             !> ManualTimeStep coming from Posti
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ParticlesInitIsDone)THEN
   SWRITE(*,*) "InitParticles already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLES ...'

!IF(.NOT.ALLOCATED(nPartsPerElem))THEN
!  ALLOCATE(nPartsPerElem(1:nElems))
!  nPartsPerElem=0
!  CALL AddToElemData(ElementOut,'nPartsPerElem',LongIntArray=nPartsPerElem(:))
!END IF

CALL InitParticleGlobals()
CALL InitializeVariables(ManualTimeStep_opt)
! InitRandomWalk must be called after InitializeVariables to know the size of TurbPartState
#if USE_RW
CALL ParticleInitRandomWalk()
#endif
! InitSGS must be called after InitRandomWalk and will abort FLEXI if an RW model is defined
CALL ParticleInitSGS()

! Restart particles here, otherwise we can not know if we need to have an initial emission
CALL ParticleRestart()
! Initialize emission. If no particles are present, assume restart from pure fluid and perform initial inserting
CALL InitializeParticleEmission()

! Initialize surface sampling
IF (WriteMacroSurfaceValues) THEN
  CALL InitParticleBoundarySampling()
END IF

#if USE_MPI
! has to be called AFTER InitializeVariables because we need to read the parameter file to know the CommSize
CALL InitParticleCommSize()
#endif

ParticlesInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticles

!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
SUBROUTINE InitializeVariables(ManualTimeStep_opt)
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_ReadInTools
!USE MOD_Equation_Vars         ,ONLY:IniRefState, RefStatePrim
USE MOD_Mesh_Vars,             ONLY:BoundaryName,BoundaryType, nBCs
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars,ONLY:PartBound,nPartBound,PartAuxBC,LowVeloRemove
USE MOD_Particle_Boundary_Vars,ONLY:nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol,UseAuxBCs
USE MOD_Particle_Interpolation,ONLY:InitParticleInterpolation
USE MOD_Particle_Mesh,         ONLY:InitFIBGM,MarkAuxBCElems
USE MOD_Particle_Mesh_Vars,    ONLY:GEO
USE MOD_Particle_MPI_Vars,     ONLY:SafetyFactor,halo_eps_velo
#if USE_MPI
USE MOD_Particle_MPI,          ONLY:InitEmissionComm
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL       :: ManualTimeStep_opt                                             !> ManualTimeStep coming from Posti
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iInit, iPartBound, iSeed
INTEGER               :: SeedSize, iBC, iExclude
INTEGER               :: iAuxBC, nAuxBCplanes, nAuxBCcylinders, nAuxBCcones, nAuxBCparabols
INTEGER               :: ALLOCSTAT
CHARACTER(32)         :: tmpStr, tmpStr2, tmpStr3
CHARACTER(200)        :: tmpString
CHARACTER(200), ALLOCATABLE        :: tmpStringBC(:)
LOGICAL               :: TrueRandom
INTEGER,ALLOCATABLE   :: iSeeds(:)
REAL                  :: n_vec(3), cos2, rmax
REAL, DIMENSION(3,1)  :: n,n1,n2
REAL, DIMENSION(3,3)  :: rot1, rot2
REAL                  :: alpha1, alpha2
!===================================================================================================================================
! Read print flags
printRandomSeeds      = GETLOGICAL(  'printRandomSeeds','.FALSE.')

! Read basic particle parameter
PDM%maxParticleNumber = GETINT(      'Part-maxParticleNumber','1')
PartGravity           = GETREALARRAY('Part-Gravity'          ,3  ,'0. , 0. , 0.')

! Allocate array to hold particle properties
ALLOCATE(PartState(       1:6,1:PDM%maxParticleNumber),    &
         PartReflCount(       1:PDM%maxParticleNumber),    &
         LastPartPos(     1:3,1:PDM%maxParticleNumber),    &
         PartPosRef(      1:3,1:PDM%MaxParticleNumber),    &
! Allocate array for Runge-Kutta time stepping
         Pt(              1:3,1:PDM%maxParticleNumber),    &
         Pt_temp(         1:6,1:PDM%maxParticleNumber),    &
         PartSpecies(         1:PDM%maxParticleNumber),    &
         PDM%ParticleInside(  1:PDM%maxParticleNumber),    &
#if USE_SM
         PDM%ParticleInsideSM(1:PDM%maxParticleNumber),    &
#endif
! Allocate array for particle position in reference coordinates
         PDM%nextFreePosition(1:PDM%maxParticleNumber),    &
         PDM%IsNewPart(       1:PDM%maxParticleNumber),    &
! Allocate particle-to-element-mapping (PEM) arrays
         PEM%Element(         1:PDM%maxParticleNumber),    &
         PEM%lastElement(     1:PDM%maxParticleNumber),    &
#if USE_SM
         PEM%hasCrossedSM(    1:PDM%maxParticleNumber),    &
#endif
         STAT=ALLOCSTAT)

IF (ALLOCSTAT.NE.0) &
  CALL abort(__STAMP__,'ERROR in particle_init.f90: Cannot allocate particle arrays!')

PDM%ParticleInside(1:PDM%maxParticleNumber)  = .FALSE.
#if USE_SM
PDM%ParticleInsideSM(1:PDM%maxParticleNumber)= .FALSE.
#endif
PDM%IsNewPart(     1:PDM%maxParticleNumber)  = .FALSE.
LastPartPos(   1:3,1:PDM%maxParticleNumber)  = 0.
PartState                                    = 0.
PartReflCount                                = 0.
Pt                                           = 0.
PartSpecies                                  = 0
PDM%nextFreePosition(1:PDM%maxParticleNumber)= 0
Pt_temp                                      = 0
PartPosRef                                   =-888.
#if USE_SM
PEM%hasCrossedSM                             = .FALSE.
#endif

! Initialize tracking of particle trapped by fouling
KeepWallParticles       = GETLOGICAL('Part-KeepWallParticles','.FALSE.')

! Output of macroscopic values
WriteMacroSurfaceValues = GETLOGICAL('Part-WriteMacroSurfaceValues','.FALSE.')

! Number of species
nSpecies                = GETINT(    'Part-nSpecies','1')
! Abort if running particle code without any species
IF (nSpecies.LE.0) &
  CALL abort(__STAMP__,'ERROR: nSpecies .LE. 0:', nSpecies)
! Allocate species array
ALLOCATE(Species(1:nSpecies))

! Loop over all species and get requested data
DO iSpec = 1, nSpecies
  WRITE(UNIT=tmpStr,FMT='(I2)') iSpec

  ! Get number of requested inits
  CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(tmpStr))//'-nInits',"")
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(ADJUSTL(tmpStr))//'-nInits','0')
  ALLOCATE(Species(iSpec)%Init(0:Species(iSpec)%NumberOfInits))

  ! Loop over all inits and get requested data
  DO iInit = 0, Species(iSpec)%NumberOfInits
    ! set help characters to read strings from parameter file
    IF(iInit.EQ.0)THEN
      tmpStr2=TRIM(ADJUSTL(tmpStr))
    ELSE ! iInit >0
      WRITE(UNIT=tmpStr2,FMT='(I0)') iInit
      tmpStr2=TRIM(ADJUSTL(tmpStr))//'-Init'//TRIM(ADJUSTL(tmpStr2))
    END IF ! iInit

    ! get species values // only once
    !--> General Species Values
    Species(iSpec)%RHSMethod             = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-RHSMethod','none'))
    Species(iSpec)%MassIC                = GETREAL(    'Part-Species'//TRIM(tmpStr2)         //'-MassIC','0.')
    Species(iSpec)%DensityIC             = GETREAL(    'Part-Species'//TRIM(tmpStr2)         //'-DensityIC','0.')
    Species(iSpec)%LowVeloThreshold      = GETREAL(    'Part-Species'//TRIM(tmpStr2)         //'-LowVeloThreshold','0.')
    Species(iSpec)%HighVeloThreshold     = GETREAL(    'Part-Species'//TRIM(tmpStr2)         //'-HighVeloThreshold','0.')

    !--> Bons particle rebound model
    Species(iSpec)%YoungIC               = GETREAL(    'Part-Species'//TRIM(tmpStr2)         //'-YoungIC','0.')
    Species(iSpec)%PoissonIC             = GETREAL(    'Part-Species'//TRIM(tmpStr2)         //'-PoissonIC','0.')
    Species(iSpec)%YieldCoeff            = GETREAL(    'Part-Species'//TRIM(tmpStr2)         //'-YieldCoeff','0.')

    ! Emission and init data
    Species(iSpec)%Init(iInit)%UseForInit            = GETLOGICAL('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-UseForInit','.TRUE.')
    Species(iSpec)%Init(iInit)%UseForEmission        = GETLOGICAL('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-UseForEmission','.TRUE.')
    Species(iSpec)%Init(iInit)%SpaceIC               = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-SpaceIC','cuboid'))
    Species(iSpec)%Init(iInit)%velocityDistribution  = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-velocityDistribution','constant'))
    Species(iSpec)%Init(iInit)%InflowRiseTime        = GETREAL('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-InflowRiseTime','0.')
    Species(iSpec)%Init(iInit)%initialParticleNumber = GETINT('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-initialParticleNumber','0')
    Species(iSpec)%Init(iInit)%RadiusIC              = GETREAL('Part-Species'//TRIM(tmpStr2)//'-RadiusIC','1.')
    Species(iSpec)%Init(iInit)%NormalIC              = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-NormalIC',3,'0. , 0. , 1.')
    Species(iSpec)%Init(iInit)%BasePointIC           = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVariance          = GETREAL('Part-Species'//TRIM(tmpStr2)//'-BaseVariance','1.')
    Species(iSpec)%Init(iInit)%BaseVector1IC         = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector1IC',3,'1. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVector2IC         = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector2IC',3,'0. , 1. , 0.')
    Species(iSpec)%Init(iInit)%CuboidHeightIC        = GETREAL('Part-Species'//TRIM(tmpStr2)//'-CuboidHeightIC','1.')
    Species(iSpec)%Init(iInit)%CylinderHeightIC      = GETREAL('Part-Species'//TRIM(tmpStr2)//'-CylinderHeightIC','1.')
    Species(iSpec)%Init(iInit)%CalcHeightFromDt      = GETLOGICAL('Part-Species'//TRIM(tmpStr2)//'-CalcHeightFromDt','.FALSE.')
!   Species(iSpec)%Init(iInit)%VeloIC                = GETREAL('Part-Species'//TRIM(tmpStr2)//'-VeloIC','0.')
    Species(iSpec)%Init(iInit)%VeloTurbIC            = GETREAL('Part-Species'//TRIM(tmpStr2)//'-VeloTurbIC','0.')
    Species(iSpec)%Init(iInit)%VeloVecIC             = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-VeloVecIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-ParticleEmissionType','2')
    Species(iSpec)%Init(iInit)%ParticleEmission      = GETREAL('Part-Species'//TRIM(ADJUSTL(tmpStr2))//'-ParticleEmission','0.')
    Species(iSpec)%Init(iInit)%NumberOfExcludeRegions= GETINT('Part-Species'//TRIM(tmpStr2)//'-NumberOfExcludeRegions','0')

    ! Nullify additional init data here
    Species(iSpec)%Init(iInit)%InsertedParticle      = 0
    Species(iSpec)%Init(iInit)%InsertedParticleSurplus = 0

    ! Get absolute value of particle velocity vector and normalize the VeloVecIC vector
    Species(iSpec)%Init(iInit)%VeloIC                = SQRT(Species(iSpec)%Init(iInit)%VeloVecIC(1)**2.                      &
                                                           +Species(iSpec)%Init(iInit)%VeloVecIC(2)**2.                      &
                                                           +Species(iSpec)%Init(iInit)%VeloVecIC(3)**2.)

    ! Only normalize if the vector does not have zero length. If it has, our job is done
    IF (Species(iSpec)%Init(iInit)%VeloIC.NE.0) THEN
      Species(iSpec)%Init(iInit)%VeloVecIC             = Species(iSpec)%Init(iInit)%VeloVecIC/Species(iSpec)%Init(iInit)%VeloIC
    END IF

    !----------- various checks/calculations after read-in of Species(i)%Init(iInit)%-data ----------------------------------!

    !--- Check if Initial ParticleInserting is really used
    IF   (((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.1).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2)) &
      .AND.(Species(iSpec)%Init(iInit)%UseForInit) ) THEN
      IF (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0) THEN
        Species(iSpec)%Init(iInit)%UseForInit=.FALSE.
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled as no ParticleNumber"
        SWRITE(*,*) " detected for Species, Init ", iSpec, iInit
      END IF
    END IF

    !--- cuboid-/cylinder-height calculation from v and dt
    IF (.NOT.Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
        ! flag is initialized with -1, compatibility issue
        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CuboidHeightIC,-1.)) THEN
          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.
          SWRITE(*,*) "WARNING: Cuboid height will be calculated from v and dt!"
        END IF
      ELSE IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') THEN
        ! flag is initialized with -1, compatibility issue
        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CylinderHeightIC,-1.)) THEN
          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.
          SWRITE(*,*) "WARNING: Cylinder height will be calculated from v and dt!"
        END IF
      END IF
    END IF
    IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cuboid') &
          .AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cylinder')) &
        CALL abort(__STAMP__,' Calculating height from v and dt is only supported for cuboid or cylinder!')
      IF (Species(iSpec)%Init(iInit)%UseForInit) &
        CALL abort(__STAMP__,' Calculating height from v and dt is not supported for initial ParticleInserting!')
    END IF

    !--- read ExcludeRegions and normalize/calculate corresponding vectors
    IF (Species(iSpec)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      ALLOCATE(Species(iSpec)%Init(iInit)%ExcludeRegion(1:Species(iSpec)%Init(iInit)%NumberOfExcludeRegions))
      IF  ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') &
       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
        DO iExclude=1,Species(iSpec)%Init(iInit)%NumberOfExcludeRegions
          WRITE(UNIT=tmpStr3,FMT='(I0)') iExclude
          tmpStr3=TRIM(tmpStr2)//'-ExcludeRegion'//TRIM(tmpStr3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC              &
            = TRIM(GETSTR('Part-Species'//TRIM(tmpStr3)//'-SpaceIC','cuboid'))
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC             &
            = GETREAL(    'Part-Species'//TRIM(tmpStr3)//'-RadiusIC','1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC            &
            = GETREAL(    'Part-Species'//TRIM(tmpStr3)//'-Radius2IC','0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC             &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-NormalIC',3,'0. , 0. , 1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC          &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-BasePointIC',3,'0. , 0. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC        &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-BaseVector1IC',3,'1. , 0. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC        &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-BaseVector2IC',3,'0. , 1. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC       &
            = GETREAL(     'Part-Species'//TRIM(tmpStr3)//'-CuboidHeightIC','1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC     &
            = GETREAL(     'Part-Species'//TRIM(tmpStr3)//'-CylinderHeightIC','1.')

          !--check and normalize data
          IF ((TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).EQ.'cuboid') .OR.             &
               ((((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1),1.)   &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2),0.))  &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3),0.))  &
            .OR. ((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1),0.)   &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2),1.))  &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3),0.))) &
            .AND.    (((ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1),0.))       &
              .AND.    (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2),0.)))      &
              .AND.    (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3),1.))))) THEN
            !-- cuboid; or BV are non-default and NormalIC is default: calc. NormalIC for ExcludeRegions from BV1/2
            !   (for def. BV and non-def. NormalIC; or all def. or non-def.: Use User-defined NormalIC when ExclRegion is cylinder)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)          &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)          &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)          &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)
          ELSE IF ( (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cuboid')    .AND. &
                    (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cylinder')) THEN
            CALL abort(__STAMP__,'Error in ParticleInit, ExcludeRegions must be cuboid or cylinder!')
          END IF

          IF (Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2 +           &
              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2 +           &
              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2 .GT. 0.) THEN
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC                     &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC                 &
              / SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2      &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2           &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1)         &
              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2)**2      &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)**2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2)         &
              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)**2      &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)**2)
          ELSE
            CALL abort(__STAMP__,'Error in ParticleInit, NormalIC Vector must not be zero!')
          END IF
        END DO !iExclude

      ! invalid combination of SpaceIC and exclude region
      ELSE
        CALL abort(__STAMP__,'Error in ParticleInit, ExcludeRegions are currently only implemented for the SpaceIC cuboid or cylinder!')
      END IF
    END IF

    !--- determine StartnumberOfInits (start loop index, e.g., for emission loops)
    !--> two options:
    !->> old option: format Part-Species(i)-***
    !->> new option: format Part-Species(i)-Init(i)-***
    IF(iInit.EQ.0)THEN
      ! only new style paramaters defined
      IF   (((Species(iSpec)%Init(iInit)%initialParticleNumber .EQ. 0 )  &
        .AND.(Species(iSpec)%Init(iInit)%ParticleEmission      .EQ. 0.)) &
        .AND.(Species(iSpec)%NumberOfInits                     .GT. 0 )) THEN
        Species(iSpec)%StartnumberOfInits = 1
      ELSE
        ! old style parameters has been defined for inits/emissions
        Species(iSpec)%StartnumberOfInits = 0
      END IF
      SWRITE(*,*) "StartnumberOfInits of Species ", iSpec, " = ", Species(iSpec)%StartnumberOfInits
    END IF ! iInit .EQ.0

  END DO ! iInit
END DO ! iSpec

! Read in boundary parameters
! Leave out this check in FLEXI even though we should probably do it
!dummy_int  = CountOption('Part-nBounds')       ! check if Part-nBounds is present in .ini file

! get number of particle boundaries
nPartBound = CountOption('Part-BoundaryName')

! Allocate arrays for particle boundaries
ALLOCATE(PartBound%SourceBoundName    (1:nBCs))
ALLOCATE(PartBound%SourceBoundType    (1:nBCs))
ALLOCATE(PartBound%TargetBoundCond    (1:nBCs))
ALLOCATE(PartBound%MomentumACC        (1:nBCs))
ALLOCATE(PartBound%WallTemp           (1:nBCs))
ALLOCATE(PartBound%WallVelo       (1:3,1:nBCs))
ALLOCATE(PartBound%WallModel          (1:nBCs))
ALLOCATE(PartBound%WallCoeffModel     (1:nBCs))
ALLOCATE(PartBound%AmbientCondition   (1:nBCs))
ALLOCATE(PartBound%AmbientConditionFix(1:nBCs))
ALLOCATE(PartBound%AmbientTemp        (1:nBCs))
ALLOCATE(PartBound%AmbientVelo    (1:3,1:nBCs))
ALLOCATE(PartBound%AmbientDens        (1:nBCs))
ALLOCATE(PartBound%AmbientDynamicVisc (1:nBCs))

! Bons particle rebound model
ALLOCATE(PartBound%Young              (1:nBCs))
ALLOCATE(PartBound%Poisson            (1:nBCs))

! Fong coefficent of restitution
ALLOCATE(PartBound%CoR                (1:nBCs))

ALLOCATE(tmpStringBC                  (1:nBCs))

! Loop over all particle boundaries and get information
!DO iPartBound=1,nPartBound
!  tmpStringBC(iPartBound) = TRIM(GETSTR('Part-BoundaryName'))
!END DO

GEO%nPeriodicVectors=0
DO iBC=1,nBCs
  IF (BoundaryType(iBC,1).EQ.0) THEN
    PartBound%TargetBoundCond(iBC) = -1
    SWRITE(*,*)"... PartBound",iBC,"is internal bound, no mapping needed"
    CYCLE
  END IF

  SELECT CASE (BoundaryType(iBC, BC_TYPE))
  CASE(1)
    tmpString='periodic'
    GEO%nPeriodicVectors=GEO%nPeriodicVectors+1
!  CASE(2,12,22)
!    tmpString='open'
  CASE(3,4)
    tmpString='reflective'
  CASE(9)
    tmpString='symmetry'
  CASE DEFAULT
    tmpString='open'
  END SELECT
!  DO iPartBound=1,nPartBound
!    IF (TRIM(BoundaryName(iBC)).EQ.tmpStringBC(iPartBound)) THEN
!      PartBound%SourceBoundType(iBC) = TRIM(GETSTR('Part-BoundaryType'))
!      PartBound%SourceBoundName(iBC) = tmpStringBC(iPartBound)
!      ! check if PartBC is not periodic, but DGBC
!      IF ((BoundaryType(iBC, BC_TYPE).EQ.1) .AND. (TRIM(PartBound%SourceBoundType(iBC)).NE.'periodic')) &
!        GEO%nPeriodicVectors=GEO%nPeriodicVectors-1
!    END IF
!  END DO
  WRITE(UNIT=tmpStr,FMT='(I0)') iBC
  PartBound%SourceBoundType(iBC) = TRIM(GETSTR('Part-Boundary'//TRIM(tmpStr)//'-Type', tmpString))
  PartBound%SourceBoundName(iBC) = TRIM(GETSTR('Part-Boundary'//TRIM(tmpStr)//'-Name', BoundaryName(iBC))) !tmpStringBC(iPartBound)
  ! Select boundary condition for particles
  SELECT CASE (PartBound%SourceBoundType(iBC))
  ! Inflow / outflow
  CASE('open')
    PartBound%TargetBoundCond(iBC)      = PartBound%OpenBC
    PartBound%AmbientCondition(iBC)     = GETLOGICAL(  'Part-Boundary'//TRIM(tmpStr)//'-AmbientCondition','.FALSE.')
    IF(PartBound%AmbientCondition(iBC)) THEN
      PartBound%AmbientConditionFix(iBC)= GETLOGICAL(  'Part-Boundary'//TRIM(tmpStr)//'-AmbientConditionFix','.TRUE.')
      PartBound%AmbientVelo(1:3,iBC)    = GETREALARRAY('Part-Boundary'//TRIM(tmpStr)//'-AmbientVelo',3,'0., 0., 0.')
      PartBound%AmbientDens(iBC)        = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-AmbientDens','0')
      PartBound%AmbientDynamicVisc(iBC) = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-AmbientDynamicVisc','1.72326582572253E-5') ! N2:T=288K
    END IF

  ! Reflective (wall)
  CASE('reflective')
    PartBound%TargetBoundCond(iBC)      = PartBound%ReflectiveBC
    PartBound%WallVelo(1:3,iBC)         = GETREALARRAY('Part-Boundary'//TRIM(tmpStr)//'-WallVelo',3,'0. , 0. , 0.')
    PartBound%WallModel(iBC)            = GETSTR(      'Part-Boundary'//TRIM(tmpStr)//'-WallModel','perfRef')

    ! Non-perfect reflection
    IF (PartBound%WallModel(iBC).EQ.'coeffRes') THEN
        PartBound%WallCoeffModel(iBC)   = GETSTR(      'Part-Boundary'//TRIM(tmpStr)//'-WallCoeffModel','Tabakoff1981')

        ! Bons particle rebound model
        IF (PartBound%WallCoeffModel(iBC).EQ.'Bons2017') THEN
            PartBound%Young(iBC)        = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-Young')
            PartBound%Poisson(iBC)      = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-Poisson')
       ELSEIF (PartBound%WallCoeffModel(iBC).EQ.'Fong2019') THEN
            PartBound%CoR(iBC)          = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-CoR','1.')
        END IF
    END IF

  ! Periodic
  CASE('periodic')
    PartBound%TargetBoundCond(iBC)      = PartBound%PeriodicBC

  ! Invalid boundary option
  CASE DEFAULT
    SWRITE(*,*) ' Boundary does not exists: ', TRIM(PartBound%SourceBoundType(iBC))
    CALL abort(__STAMP__,'Particle Boundary Condition does not exist')
  END SELECT
END DO

GEO%nPeriodicVectors=GEO%nPeriodicVectors/2


!--- Read Manual Time Step
useManualTimeStep = .FALSE.
!> Check if we're running Posti
IF (.NOT.PRESENT(ManualTimeStep_opt)) ManualTimeStep    = GETREAL('Particles-ManualTimeStep', '0.0')
IF (ManualTimeStep.GT.0.0)            useManualTimeStep = .True.

!--- initialize randomization (= random if one or more seeds are 0 or random is wanted)
nrSeeds = GETINT('Part-NumberOfRandomSeeds','0')
IF (nrSeeds.GT.0) THEN
   ALLOCATE(seeds(1:nrSeeds))
   DO iSeed = 1, nrSeeds
      WRITE(UNIT=tmpStr,FMT='(I0)') iSeed
      seeds(iSeed) = GETINT('Particles-RandomSeed'//TRIM(tmpStr),'0')
   END DO
END IF

CALL RANDOM_SEED(Size = SeedSize)                ! Check for number of needed Seeds
TrueRandom = .FALSE.                             ! FALSE for defined random seed

! Check if the correct number of random seeds was given
IF (nrSeeds.GT.0) THEN
   IF (nrSeeds.NE.SeedSize) THEN
      IF(printRandomSeeds)THEN
        IPWRITE(UNIT_StdOut,*) 'Error: Number of seeds for RNG must be ',SeedSize
        IPWRITE(UNIT_StdOut,*) 'Random RNG seeds are used'
      END IF
      TrueRandom = .TRUE.
   END IF
   DO iSeed = 1, nrSeeds
      IF (Seeds(iSeed).EQ.0) THEN
         IF(printRandomSeeds)THEN
           IPWRITE(UNIT_StdOut,*) 'Error: ',SeedSize,' seeds for RNG must be defined'
           IPWRITE(UNIT_StdOut,*) 'Random RNG seeds are used'
         END IF
         TrueRandom = .TRUE.
      END IF
   END DO
ELSE
   TrueRandom = .TRUE.
END IF

! Get true random number if requested
IF (TrueRandom) THEN
! Changed in FLEXI to get true random numbers FOR EACH PROC
   SDEALLOCATE(seeds)
   ALLOCATE(seeds(SeedSize))
   CALL RANDOM_SEED(GET = seeds(1:SeedSize))
#if USE_MPI
   Seeds(1:SeedSize)    = Seeds(1:SeedSize)+PartMPI%MyRank
#endif
   CALL RANDOM_SEED(PUT = Seeds(1:SeedSize))
ELSE
#if USE_MPI
   Seeds(1:SeedSize)    = Seeds(1:SeedSize)+PartMPI%MyRank
#endif
   CALL RANDOM_SEED(PUT = Seeds(1:SeedSize))
END IF

! Get final random seeds and print if requested
ALLOCATE(iseeds(SeedSize))
iseeds(:)=0
CALL RANDOM_SEED(GET = iseeds(1:SeedSize))
IF(printRandomSeeds)THEN
  IPWRITE(UNIT_StdOut,*) 'Random seeds in Particle_Interpolation_init:'
  DO iSeed = 1,SeedSize
     IPWRITE(UNIT_StdOut,*) iseeds(iSeed)
  END DO
END IF
DEALLOCATE(iseeds)

! Time delay before initial particle insering
DelayTime         = GETREAL(   'Part-DelayTime'    ,'0.')

! Flag if a lost particle should abort the program
AllowLoosing      = GETLOGICAL('Part-AllowLoosing' ,'.FALSE.')

! Flag if low velocity particles should be removed
LowVeloRemove     = GETLOGICAL('Part-LowVeloRemove','.FALSE.')

! AuxBCs
nAuxBCs=GETINT('Part-nAuxBCs','0')

!--> Only read anything if there are auxiliary BCs
IF (nAuxBCs.GT.0) THEN
  UseAuxBCs=.TRUE.
  ALLOCATE (AuxBCType(1:nAuxBCs) &
            ,AuxBCMap(1:nAuxBCs) )
  AuxBCMap=0

  !- Read in BC parameters
  ALLOCATE(PartAuxBC%TargetBoundCond  (1:nAuxBCs))
  ALLOCATE(PartAuxBC%MomentumACC      (1:nAuxBCs))
  ALLOCATE(PartAuxBC%WallTemp         (1:nAuxBCs))
  ALLOCATE(PartAuxBC%WallVelo     (1:3,1:nAuxBCs))

  ! Get auxiliary boundary types
  DO iPartBound=1,nAuxBCs
    WRITE(UNIT=tmpStr,FMT='(I0)') iPartBound
    tmpString = TRIM(GETSTR('Part-AuxBC'//TRIM(tmpStr)//'-Condition','open'))

    SELECT CASE (TRIM(tmpString))

    ! Inflow / outflow
    CASE('open')
      PartAuxBC%TargetBoundCond(iPartBound) = PartAuxBC%OpenBC

    ! Reflective (wall)
    CASE('reflective')
      PartAuxBC%TargetBoundCond(iPartBound) = PartAuxBC%ReflectiveBC
      PartAuxBC%MomentumACC(iPartBound)     = GETREAL(     'Part-AuxBC'//TRIM(tmpStr)//'-MomentumACC','0')
      PartAuxBC%WallTemp(iPartBound)        = GETREAL(     'Part-AuxBC'//TRIM(tmpStr)//'-WallTemp','0')
      PartAuxBC%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-WallVelo',3,'0. , 0. , 0.')

    ! Unknown auxiliary boundary type
    CASE DEFAULT
      SWRITE(*,*) ' AuxBC Condition does not exists: ', TRIM(tmpString)
      CALL abort(__STAMP__,'AuxBC Condition does not exist')

    END SELECT
  END DO

  !- read and count types
  nAuxBCplanes    = 0
  nAuxBCcylinders = 0
  nAuxBCcones     = 0
  nAuxBCparabols  = 0

  DO iAuxBC=1,nAuxBCs
    WRITE(UNIT=tmpStr,FMT='(I0)') iAuxBC
    AuxBCType(iAuxBC) = TRIM(GETSTR('Part-AuxBC'//TRIM(tmpStr)//'-Type','plane'))

    SELECT CASE (TRIM(AuxBCType(iAuxBC)))

    CASE ('plane')
      nAuxBCplanes     = nAuxBCplanes + 1
      AuxBCMap(iAuxBC) = nAuxBCplanes

    CASE ('cylinder')
      nAuxBCcylinders  = nAuxBCcylinders + 1
      AuxBCMap(iAuxBC) = nAuxBCcylinders

    CASE ('cone')
      nAuxBCcones      = nAuxBCcones + 1
      AuxBCMap(iAuxBC) = nAuxBCcones

    CASE ('parabol')
      nAuxBCparabols   = nAuxBCparabols + 1
      AuxBCMap(iAuxBC) = nAuxBCparabols

    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(__STAMP__,'AuxBC does not exist')
    END SELECT
  END DO

  !- allocate type-specifics
  IF (nAuxBCplanes.GT.0) THEN
    ALLOCATE (AuxBC_plane(   1:nAuxBCplanes))
  END IF

  IF (nAuxBCcylinders.GT.0) THEN
    ALLOCATE (AuxBC_cylinder(1:nAuxBCcylinders))
  END IF

  IF (nAuxBCcones.GT.0) THEN
    ALLOCATE (AuxBC_cone(    1:nAuxBCcones))
  END IF

  IF (nAuxBCparabols.GT.0) THEN
    ALLOCATE (AuxBC_parabol( 1:nAuxBCparabols))
  END IF

  !- read type-specifics
  DO iAuxBC=1,nAuxBCs
    WRITE(UNIT=tmpStr,FMT='(I0)') iAuxBC
    SELECT CASE (TRIM(AuxBCType(iAuxBC)))

    CASE ('plane')
      AuxBC_plane(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-r_vec',3,'0. , 0. , 0.')
      WRITE(UNIT=tmpStr2,FMT='(G0)') HUGE(AuxBC_plane(AuxBCMap(iAuxBC))%radius)
      AuxBC_plane(AuxBCMap(iAuxBC))%radius= GETREAL(     'Part-AuxBC'//TRIM(tmpStr)//'-radius',TRIM(tmpStr2))
      n_vec                               = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-n_vec',3,'1. , 0. , 0.')
      ! Check if normal vector is zero
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(__STAMP__,'Part-AuxBC-n_vec is zero for AuxBC',iAuxBC)
      ! If not, scale vector
      ELSE
        AuxBC_plane(AuxBCMap(iAuxBC))%n_vec = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF

    CASE ('cylinder')
      AuxBC_cylinder(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                                  = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-axis',3,'1. , 0. , 0.')
      ! Check if normal vector is zero
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(__STAMP__,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ! If not, scale vector
      ELSE
        AuxBC_cylinder(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF

      AuxBC_cylinder(AuxBCMap(iAuxBC))%radius  = GETREAL(   'Part-AuxBC'//TRIM(tmpStr)//'-radius','1.')
      WRITE(UNIT=tmpStr2,FMT='(G0)') -HUGE(AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin    = GETREAL(   'Part-AuxBC'//TRIM(tmpStr)//'-lmin',TRIM(tmpStr2))
      WRITE(UNIT=tmpStr2,FMT='(G0)') HUGE(AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cylinder(AuxBCMap(iAuxBC))%lmax    = GETREAL(   'Part-AuxBC'//TRIM(tmpStr)//'-lmax',TRIM(tmpStr2))
      AuxBC_cylinder(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(tmpStr)//'-inwards','.TRUE.')

    CASE ('cone')
      AuxBC_cone(AuxBCMap(iAuxBC))%r_vec     = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                                  = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-axis',3,'1. , 0. , 0.')
      ! Check if normal vector is zero
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(__STAMP__,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ! If not, scale vector
      ELSE
        AuxBC_cone(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF

      AuxBC_cone(AuxBCMap(iAuxBC))%lmin  = GETREAL('Part-AuxBC'//TRIM(tmpStr)//'-lmin','0.')
      IF (AuxBC_cone(AuxBCMap(iAuxBC))%lmin.LT.0.) &
        CALL abort(__STAMP__,'Part-AuxBC-lminis .lt. zero for AuxBC',iAuxBC)

      WRITE(UNIT=tmpStr2,FMT='(G0)') HUGE(AuxBC_cone(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cone(AuxBCMap(iAuxBC))%lmax  = GETREAL('Part-AuxBC'//TRIM(tmpStr)//'-lmax',TRIM(tmpStr2))
      rmax                               = GETREAL('Part-AuxBC'//TRIM(tmpStr)//'-rmax','0.')

      ! either define rmax at lmax or the halfangle
      IF (rmax.EQ.0.) THEN
        AuxBC_cone(AuxBCMap(iAuxBC))%halfangle  = GETREAL('Part-AuxBC'//TRIM(tmpStr)//'-halfangle','45.')*PI/180.
      ELSE
        AuxBC_cone(AuxBCMap(iAuxBC))%halfangle  = ATAN(rmax/AuxBC_cone(AuxBCMap(iAuxBC))%lmax)
      END IF

      IF (AuxBC_cone(AuxBCMap(iAuxBC))%halfangle.LE.0.) &
        CALL abort(__STAMP__,'Part-AuxBC-halfangle is .le. zero for AuxBC',iAuxBC)

      AuxBC_cone(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(tmpStr)//'-inwards','.TRUE.')
      cos2 = COS(AuxBC_cone(AuxBCMap(iAuxBC))%halfangle)**2
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,1) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(1)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/cos2,0.,0./)
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,2) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(2)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/0.,cos2,0./)
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,3) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(3)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/0.,0.,cos2/)

    CASE ('parabol')
      AuxBC_parabol(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                                 = GETREALARRAY('Part-AuxBC'//TRIM(tmpStr)//'-axis',3,'1. , 0. , 0.')
      ! Check if normal vector is zero
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(__STAMP__,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ! If not, scale vector
      ELSE
        AuxBC_parabol(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF

      AuxBC_parabol(AuxBCMap(iAuxBC))%lmin  = GETREAL(     'Part-AuxBC'//TRIM(tmpStr)//'-lmin','0.')
      IF (AuxBC_parabol(AuxBCMap(iAuxBC))%lmin.LT.0.) &
        CALL abort(__STAMP__,'Part-AuxBC-lmin is .lt. zero for AuxBC',iAuxBC)

      WRITE(UNIT=tmpStr2,FMT='(G0)') HUGE(AuxBC_parabol(AuxBCMap(iAuxBC))%lmin)
      AuxBC_parabol(AuxBCMap(iAuxBC))%lmax  = GETREAL(     'Part-AuxBC'//TRIM(tmpStr)//'-lmax',TRIM(tmpStr2))
      AuxBC_parabol(AuxBCMap(iAuxBC))%zfac  = GETREAL(     'Part-AuxBC'//TRIM(tmpStr)//'-zfac','1.')
      AuxBC_parabol(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(tmpStr)//'-inwards','.TRUE.')

      n(:,1)=AuxBC_parabol(AuxBCMap(iAuxBC))%axis

      ! Check if normal vector is colliniar with y?
      IF (.NOT.ALMOSTZERO(SQRT(n(1,1)**2+n(3,1)**2))) THEN
        alpha1 = ATAN2(n(1,1),n(3,1))
        CALL roty(rot1,alpha1)
        n1     = MATMUL(rot1,n)
      ELSE
        alpha1 = 0.
        CALL ident(rot1)
        n1     = n
      END IF

      ! Check if normal vector is colliniar with x?
      IF (.NOT.ALMOSTZERO(SQRT(n1(2,1)**2+n1(3,1)**2))) THEN
        alpha2 = -ATAN2(n1(2,1),n1(3,1))
        CALL rotx(rot2,alpha2)
        n2     = MATMUL(rot2,n1)
      ELSE
        CALL abort(__STAMP__,'Vector is collinear with x-axis. this should not be possible... AuxBC:',iAuxBC)
      END IF

      AuxBC_parabol(AuxBCMap(iAuxBC))%rotmatrix(:,:)  = MATMUL(rot2,rot1)
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(:,:) = 0.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(1,1) = 1.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(2,2) = 1.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(3,3) = 0.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(3,4) = -0.5*AuxBC_parabol(AuxBCMap(iAuxBC))%zfac
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(4,3) = -0.5*AuxBC_parabol(AuxBCMap(iAuxBC))%zfac

    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(__STAMP__,'AuxBC does not exist for AuxBC',iAuxBC)

    END SELECT
  END DO

  ! Mark elements with auxiliary BCs
  CALL MarkAuxBCElems()
ELSE
  ! Flag if AuxBCs are used
  UseAuxBCs=.FALSE.
END IF

! Initialize interpolation and particle-in-cell for field -> particle coupling
!--> Could not be called earlier because a halo region has to be build depending on the given BCs
CALL InitParticleInterpolation()

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT FIBGM...'
SafetyFactor  = GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo = GETREAL('Particles-HaloEpsVelo','0')

! Initialize Fast-Init BackGround Mesh (FIBGM)
CALL InitFIBGM()

! Initialize MPI communicator for emmission procs
#if USE_MPI
CALL InitEmissionComm()
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/

SDEALLOCATE(tmpStringBC)

SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitializeVariables


!===================================================================================================================================
! finalize particle variables
!===================================================================================================================================
SUBROUTINE FinalizeParticles()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Interpolation_Vars
#if USE_RW
USE MOD_Particle_RandomWalk,    ONLY: ParticleFinalizeRandomWalk
#endif
USE MOD_Particle_SGS,           ONLY: ParticleFinalizeSGS
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Pt_temp)
SDEALLOCATE(PartPosRef)
SDEALLOCATE(PartState)
SDEALLOCATE(LastPartPos)
SDEALLOCATE(PartSpecies)
SDEALLOCATE(Pt)
SDEALLOCATE(PDM%ParticleInside)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%IsNewPart)
SDEALLOCATE(Species)
SDEALLOCATE(PartBound%SourceBoundName)
SDEALLOCATE(PartBound%SourceBoundType)
SDEALLOCATE(PartBound%TargetBoundCond)
SDEALLOCATE(PartBound%MomentumACC)
SDEALLOCATE(PartBound%WallTemp)
SDEALLOCATE(PartBound%WallVelo)
SDEALLOCATE(PartBound%AmbientCondition)
SDEALLOCATE(PartBound%AmbientConditionFix)
SDEALLOCATE(PartBound%AmbientTemp)
SDEALLOCATE(PartBound%AmbientVelo)
SDEALLOCATE(PartBound%AmbientDens)
SDEALLOCATE(PartBound%AmbientDynamicVisc)
SDEALLOCATE(PartBound%WallModel)
SDEALLOCATE(PartBound%WallCoeffModel)
SDEALLOCATE(PartBound%Young)
SDEALLOCATE(PartBound%Poisson)
SDEALLOCATE(PartBound%CoR)
SDEALLOCATE(PEM%Element)
SDEALLOCATE(PEM%lastElement)
SDEALLOCATE(PEM%pStart)
SDEALLOCATE(PEM%pNumber)
SDEALLOCATE(PEM%pEnd)
SDEALLOCATE(PEM%pNext)
SDEALLOCATE(FieldAtParticle)
! Sliding Mesh
#if USE_SM
SDEALLOCATE(PEM%hasCrossedSM)
SDEALLOCATE(PDM%ParticleInsideSM)
! SGS/RW turbulence model
#endif
#if USE_RW
SDEALLOCATE(TurbFieldAtParticle)
CALL ParticleFinalizeRandomWalk()
#endif
CALL ParticleFinalizeSGS()

END SUBROUTINE FinalizeParticles


!===================================================================================================================================
SUBROUTINE rotx(mat,a)
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
REAL, INTENT(IN)                  :: a
!===================================================================================================================================
mat(:,1)=(/1.0 , 0.     , 0.  /)
mat(:,2)=(/0.0 , COS(a) ,-SIN(a)/)
mat(:,3)=(/0.0 , SIN(a) , COS(a)/)
END SUBROUTINE


!===================================================================================================================================
SUBROUTINE roty(mat,a)
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
REAL, INTENT(IN)                  :: a
!===================================================================================================================================
mat(:,1)=(/ COS(a) , 0., SIN(a)/)
mat(:,2)=(/ 0.     , 1., 0.  /)
mat(:,3)=(/-SIN(a) , 0., COS(a)/)
END SUBROUTINE


!===================================================================================================================================
SUBROUTINE ident(mat)
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
INTEGER                           :: j
!===================================================================================================================================

mat                      = 0.
FORALL(j = 1:3) mat(j,j) = 1.
END SUBROUTINE

END MODULE MOD_ParticleInit
