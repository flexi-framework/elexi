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
! Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_ParticleInit
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

! Private Part --------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersParticles
  MODULE PROCEDURE DefineParametersParticles
END INTERFACE

INTERFACE InitParticles
  MODULE PROCEDURE InitParticles
END INTERFACE

INTERFACE FinalizeParticles
  MODULE PROCEDURE FinalizeParticles
END INTERFACE

PUBLIC::DefineParametersParticles
PUBLIC::InitParticles,FinalizeParticles
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
CALL prms%CreateIntOption('Part-IterationForMacroVal',      'Set number of iterations used for sampling if Part-WriteMacroValues'//&
                                                            ' is set true.',                                '1')
CALL prms%CreateIntOption('Part-nAuxBCs',                   'Number of auxillary BCs that are checked during tracing',  '0')
CALL prms%CreateIntOption('Part-nBounds',                   'Number of particle boundaries.',               '1')
CALL prms%CreateIntOption('Part-nPeriodicVectors',          'Number of the periodic vectors j=1,...,n.'                          //&
                                                            ' Value has to be the same as defined in preprog.ini',      '0')
CALL prms%CreateIntOption('Part-nSpecies',                  'Number of species in part',                    '1')
CALL prms%CreateIntOption('Part-maxParticleNumber',         'Max number of particles in part',              '1')
CALL prms%CreateIntOption('Part-NumberOfRandomSeeds',       'Number of random seeds for particle random number generator','0')
CALL prms%CreateIntOption('Particles-NumberOfRandomVectors','Number of random vectors',                     '0')

CALL prms%CreateLogicalOption('OutputSurfaceFluxLinked',    'Flag to print the SurfaceFlux-linked Info' ,   '.FALSE.')
CALL prms%CreateLogicalOption('Part-SteadyState',           'Only update particle position while keeping fluid state frozen',      &
                                                            '.FALSE.')
CALL prms%CreateLogicalOption('Part-TrackPosition',         'Track particle position',                      '.FALSE.')
CALL prms%CreateLogicalOption('Part-WriteMacroValues',      'Set [T] to activate ITERATION DEPENDANT h5 output of '              //&
                                                            'macroscopic values sampled every [Part-IterationForMacroVal] '      //&
                                                            'iterations from particles. Sampling starts from simulation'         //&
                                                            ' start. Can not be enabled together with Part-TimeFracForSampling',   &
                                                                                                            '.FALSE.')
CALL prms%CreateLogicalOption('Part-WriteMacroVolumeValues','Similar to Part-WriteMacroValues. Set [T] to activate iteration'    //&
                                                            ' dependant sampling and h5 output for each element.',                 &
                                                                                                            '.FALSE.')
CALL prms%CreateLogicalOption('Part-WriteMacroSurfaceValues','Similar to Part-WriteMacroValues. Set [T] to activate iteration'   //&
                                                            ' dependant sampling and h5 output on surfaces.',                      &
                                                                                                            '.FALSE.')
CALL prms%CreateLogicalOption('Part-WriteFieldsToVTK',      'Not in Code anymore, but read-in has to be deleted'                 //&
                                                            ' in particle_init.f90',                        '.FALSE.')
!CALL prms%CreateLogicalOption('PIC-SFResampleAnalyzeSurfCollis', '',                                        '.FALSE.')
CALL prms%CreateLogicalOption('printRandomSeeds',           'Flag defining if random seeds are written.',   '.FALSE.')
!CALL prms%CreateLogicalOption('PrintSFDepoWarnings',        'Print the shapefunction warnings',             '.FALSE.')
CALL prms%CreateLogicalOption('Particles-OutputVpiWarnings','Flag for warnings for rejected'                                     //&
                                                            ' v if VPI+PartDensity',                        '.FALSE.')
CALL prms%CreateLogicalOption('Part-AllowLoosing',          'Flag if a lost particle should abort the programm','.FALSE.')
CALL prms%CreateLogicalOption('Part-LowVeloRemove',         'Flag if low velocity particles should be removed', '.FALSE.')

CALL prms%CreateRealOption('Part-DelayTime',                "During delay time the particles won't be moved so "                 //&
                                                            "the fluid field can be evolved",               '0.')
CALL prms%CreateRealOption('Part-SafetyFactor',             'Factor to scale the halo region with MPI',     '1.')
CALL prms%CreateRealOption('Part-SteadyTimeStep',           'Manual time step routine for frozen fluid state', '0.')
CALL prms%CreateRealOption('Particles-HaloEpsVelo',         'Maximum velocity to be considered for halo region', '0.')
CALL prms%CreateRealOption('Particles-ManualTimeStep',      'Manual time step routine for frozen fluid state', '0.')

CALL prms%CreateRealArrayOption('Part-Gravity',             'Gravitational acceleration as vector',         '0. , 0. , 0.')

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
                                  '   constant: all particles have the same defined velocity.(VeloIC, VeloVec)\n'//&
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
CALL prms%CreateRealOption(     'Part-Species[$]-PartDensity' &
                                , 'Define particle density for species [$]. PartDensity (real particles per m^3).\n'//&
                                  'Used for DSMC with (vpi_)cuboid/cylinder and cell_local initial inserting.\n'// &
                                  'Also for (vpi_)cub./cyl. / cell_local as alternative to Part.Emis. in Type1', '0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-ParticleEmission' &
                                , 'Emission rate in part/s or part/iteration.', '0.', numberedmulti=.TRUE.)
    
CALL prms%CreateRealOption(     'Part-Species[$]-ChargeIC'  &
                                , 'Particle Charge (without MPF) of species[$] dim' &
                                , '0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-MacroParticleFactor' &
                                , 'Number of Microparticle per Macroparticle for species [$]', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MassIC'  &
                                , 'Particle Mass (without MPF) of species [$] [kg]', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-DensityIC'  &
                                , 'Particle density (without MPF) of species [$] [kg/m^3]', '0.', numberedmulti=.TRUE.)    
CALL prms%CreateRealOption(     'Part-Species[$]-RadiusIC'  &
                                , 'Radius for IC circle', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-VeloIC'  &
                                , 'Absolute value of initial velocity. (ensemble velocity) ', '0.', numberedmulti=.TRUE.)
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
CALL prms%CreateIntOption(      'Part-Species[$]-nSurfacefluxBCs'&
                                , 'Number of SF emissions', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-ParticleEmissionType'  &
                                , 'Define Emission Type for particles (volume emission)\n'//&
                                  '1 = emission rate in part/s,\n'//&
                                  '2 = emission rate part/iteration\n'//&
                                  '3 = user def. emission rate\n', '2', numberedmulti=.TRUE.)
    
CALL prms%CreateRealArrayOption('Part-Species[$]-BasePointIC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Base point for IC cuboid and IC sphere', '0. , 0. , 0.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-BaseVector1IC'  &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'First base vector for IC cuboid', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-BaseVector2IC'  &
                                , 'Second base vector for IC cuboid', '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-NormalIC'  &
                                , 'Normal orientation of circle.', '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-VeloVecIC '  &
                                , 'Normalized velocity vector for given VeloIC', '0. , 0. , 0.', numberedmulti=.TRUE.)
                                
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
                                
!===================================================================================================================================
! > Random Walk (Subgroup of Part-Species)  
! >>> Values in this section only apply for turbulence models providing turbulent kinetic energy and a turbulent length/time scale
!===================================================================================================================================
#if EQNSYSNR == 4
CALL prms%CreateStringOption(   'Part-Species[$]-RWModel' &
                                , 'Random walk model used for RANS calculations.\n'//&
                                ' - Gosman \n'//&
                                ' - Dehbi \n'//&
                                ' - Langevin \n'&
                                , 'none' , numberedmulti=.TRUE.)
#endif
                                
!===================================================================================================================================
! > Boundaries       
! >>> Values in this section appear multiple times
!===================================================================================================================================
CALL prms%SetSection("Particle Boundaries")
                                
CALL prms%CreateStringOption(   'Part-Boundary[$]-Condition'  &
                                , 'Used boundary condition for boundary[$].\n'//&
                                  '- open\n'//&
                                  '- reflective\n'//&
                                  '- periodic\n'//&
                                 'If condition=open, the following parameters are'//&
                                  ' used: (Part-Boundary[$]-=PB) PB-Ambient ,PB-AmbientTemp,PB-AmbientMeanPartMass,'//&
                                  'PB-AmbientVelo,PB-AmbientDens,PB-AmbientDynamicVisc,PB-AmbientThermalCond,PB-Voltage\n'//&
                                 'If condition=reflective: PB-MomentumACC,PB-WallTemp,PB-TransACC,PB-VibACC,PB-RotACC,'//&
                                  'PB-WallVelo,Voltage,SpeciesSwaps.If condition=periodic:Part-nPeriodicVectors,'//&
                                  'Part-PeriodicVector[$]', 'open', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Boundary[$]-SourceName'  &
                                , 'No Default. Source Name of Boundary[i]. Has to be selected for all'//&
                                  'nBounds. Has to be same name as defined in preproc tool', numberedmulti=.TRUE.)
                                  
CALL prms%CreateRealArrayOption('Part-PeriodicVector[$]'      , 'TODO-DEFINE-PARAMETER\nVector for periodic boundaries.'//&
                                                                   'Has to be the same as defined in preproc.ini in their'//&
                                                                   ' respective order. ', '1. , 0. , 0.', numberedmulti=.TRUE.)
                                  
! Wall model =======================================================================================================================
CALL prms%CreateStringOption(   'Part-Boundary[$]-WallModel' &
                                , 'Wall model to be used. Options:.\n'//&
                                  'perfRef  - perfect reflection\n'//&
                                  'coeffRes - Coefficient of restitution','perfRef', numberedmulti=.TRUE.)
CALL prms%CreateStringOption('Part-Boundary[$]-WallCoeffModel' &
                                , 'Coefficients to be used. Options:.\n'//&
                                  'Tabakoff1981','Tabakoff1981', numberedmulti=.TRUE.)
                                  
! Ambient condition=================================================================================================================
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-AmbientCondition'  &
                                , 'Use ambient condition (condition "behind" boundary).', '.FALSE.'&
                                , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Boundary[$]-AmbientConditionFix'  &
                                , 'TODO-DEFINE-PARAMETER', '.TRUE.', numberedmulti=.TRUE.)
    
CALL prms%CreateRealArrayOption('Part-Boundary[$]-AmbientVelo'  &
                                , 'Ambient velocity', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Boundary[$]-WallVelo'  &
                                , 'Velocity (global x,y,z in [m/s]) of reflective particle boundary [$].' &
                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
    
CALL prms%CreateRealOption(     'Part-Boundary[$]-AmbientDens'  &
                                , 'Ambient density', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Boundary[$]-AmbientDynamicVisc'  &
                                , 'Ambient dynamic viscosity', '1.72326582572253E-5', numberedmulti=.TRUE.)
    
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
USE MOD_Particle_Vars,              ONLY: ParticlesInitIsDone,WriteMacroSurfaceValues
USE MOD_Part_Emission,              ONLY: InitializeParticleEmission, InitializeParticleSurfaceflux
USE MOD_DSMC_Vars,                  ONLY: DSMC
USE MOD_Particle_Boundary_Sampling, ONLY: InitParticleBoundarySampling
USE MOD_Particle_Erosion_Vars
#if USE_MPI
USE MOD_Particle_MPI,               ONLY: InitParticleCommSize
#endif
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

Call InitParticleGlobals

CALL InitializeVariables(ManualTimeStep_opt)
CALL InitializeParticleEmission()
CALL InitializeParticleSurfaceflux()

! Initialize surface sampling
IF (WriteMacroSurfaceValues.OR.DSMC%CalcSurfaceVal) THEN
  CALL InitParticleBoundarySampling()
END IF

#if USE_MPI
! has to be called AFTER InitializeVariables and InitDSMC 
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
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars,ONLY:PartBound,nPartBound,PartAuxBC,nAdaptiveBC,LowVeloRemove
USE MOD_Particle_Boundary_Vars,ONLY:nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol,UseAuxBCs
USE MOD_Particle_Mesh_Vars    ,ONLY:NbrOfRegions,RegionBounds
USE MOD_Mesh_Vars,             ONLY:BoundaryName,BoundaryType, nBCs
USE MOD_Particle_Surfaces_Vars,ONLY:BCdata_auxSF
USE MOD_DSMC_Vars,             ONLY:DSMC
USE MOD_PICInterpolation,      ONLY:InitializeInterpolation
USE MOD_PICInit,               ONLY:InitPIC
USE MOD_Particle_Mesh,         ONLY:InitFIBGM,MapRegionToElem,MarkAuxBCElems
USE MOD_Particle_Tracking_Vars,ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,     ONLY:SafetyFactor,halo_eps_velo!,PartMPI
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
INTEGER               :: iSpec, iInit, iPartBound, iSeed!, iCC
INTEGER               :: SeedSize, iPBC, iBC, iSwaps, iRegions, iExclude
INTEGER               :: iAuxBC, nAuxBCplanes, nAuxBCcylinders, nAuxBCcones, nAuxBCparabols
INTEGER               :: ALLOCSTAT
INTEGER               :: dummy_int
CHARACTER(32)         :: hilf, hilf2, hilf3
CHARACTER(200)        :: tmpString
LOGICAL               :: PartDens_OnlyInit, TrueRandom
INTEGER,ALLOCATABLE   :: iSeeds(:)
REAL                  :: lineVector(3), v_drift_line, A_ins, n_vec(3), cos2, rmax
INTEGER               :: MaxNbrOfSpeciesSwaps
REAL, DIMENSION(3,1)  :: n,n1,n2
REAL, DIMENSION(3,3)  :: rot1, rot2
REAL                  :: alpha1, alpha2
LOGICAL,ALLOCATABLE   :: MacroRestartFileUsed(:)
INTEGER               :: FileID
!===================================================================================================================================
! Read print flags
printRandomSeeds = GETLOGICAL('printRandomSeeds','.FALSE.')

! Read basic particle parameter
PDM%maxParticleNumber = GETINT('Part-maxParticleNumber','1')
PartGravity           = GETREALARRAY('Part-Gravity',3,'0. , 0. , 0.')
ALLOCATE(Pt_temp(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT)  
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF
Pt_temp=0.

!IF(DoRefMapping)THEN
! We got rid of pic_depo, so we need to allocate it here
  ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,' Cannot allocate partposref!')
  PartPosRef=-888.
!END IF

! predefine random vectors
NumRanVec = GETINT('Particles-NumberOfRandomVectors','100000')

ALLOCATE(PartState(1:PDM%maxParticleNumber,1:6)       , &
         PartReflCount(1:PDM%maxParticleNumber)      , &
         LastPartPos(1:PDM%maxParticleNumber,1:3)     , &
         Pt(1:PDM%maxParticleNumber,1:3)              , &
         PartSpecies(1:PDM%maxParticleNumber)         , &
         PDM%ParticleInside(1:PDM%maxParticleNumber)  , &
         PDM%nextFreePosition(1:PDM%maxParticleNumber), &
         PDM%dtFracPush(1:PDM%maxParticleNumber)      , &
         PDM%IsNewPart(1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR in particle_init.f90: Cannot allocate Particle arrays!')
END IF

! output of macroscopic values
WriteMacroValues        = GETLOGICAL('Part-WriteMacroValues',       '.FALSE.')
WriteMacroVolumeValues  = GETLOGICAL('Part-WriteMacroVolumeValues', '.FALSE.')
WriteMacroSurfaceValues = GETLOGICAL('Part-WriteMacroSurfaceValues','.FALSE.')
IF(WriteMacroValues)THEN
  WriteMacroVolumeValues  = .TRUE.
  WriteMacroSurfaceValues = .TRUE.
ELSE IF((WriteMacroVolumeValues.AND.WriteMacroSurfaceValues).AND.(.NOT.WriteMacroValues))THEN
  WriteMacroValues        = .TRUE.
END IF
MacroValSamplIterNum    = GETINT(    'Part-IterationForMacroVal',   '1')
DSMC%TimeFracSamp       = GETREAL(   'Part-TimeFracForSampling',    '0.0')
DSMC%CalcSurfaceVal     = GETLOGICAL('Particles-CalcSurfaceVal',    '.FALSE.')

IF(WriteMacroVolumeValues.OR.WriteMacroSurfaceValues)THEN
  IF(DSMC%TimeFracSamp.GT.0.0) CALL abort(&
__STAMP__&
    ,'ERROR: Init Macrosampling: WriteMacroValues and Time fraction sampling enabled at the same time')
  IF(WriteMacroSurfaceValues.AND.(.NOT.DSMC%CalcSurfaceVal)) DSMC%CalcSurfaceVal = .TRUE.
END IF

PDM%ParticleInside(1:PDM%maxParticleNumber)  = .FALSE.
PDM%dtFracPush(1:PDM%maxParticleNumber)      = .FALSE.
PDM%IsNewPart(1:PDM%maxParticleNumber)       = .FALSE.
LastPartPos(1:PDM%maxParticleNumber,1:3)     = 0.
PartState                                    = 0.
PartReflCount                                = 0.
Pt                                           = 0.
PartSpecies                                  = 0
PDM%nextFreePosition(1:PDM%maxParticleNumber)= 0

! Initialize of adsorbed particle tracking 
KeepWallParticles = .FALSE.

nSpecies = GETINT('Part-nSpecies','1')

! Read particle species data
!nSpecies = CNTSTR('Part-Species-SpaceIC')

IF (nSpecies.LE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nSpecies .LE. 0:', nSpecies)
END IF

IF (nSpecies.LE.0) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nSpecies .LE. 0:', nSpecies)
END IF

ALLOCATE(Species(1:nSpecies))

DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  
  ! Define GETINT
  CALL prms%CreateIntOption('Part-Species'//TRIM(ADJUSTL(hilf))//'-nInits',"")
  
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(ADJUSTL(hilf))//'-nInits','0')
  ALLOCATE(Species(iSpec)%Init(0:Species(iSpec)%NumberOfInits)) 
  DO iInit = 0, Species(iSpec)%NumberOfInits
    ! set help characters
    IF(iInit.EQ.0)THEN
      hilf2=TRIM(ADJUSTL(hilf))
    ELSE ! iInit >0
      WRITE(UNIT=hilf2,FMT='(I0)') iInit
      hilf2=TRIM(ADJUSTL(hilf))//'-Init'//TRIM(ADJUSTL(hilf2))
    END IF ! iInit
    ! get species values // only once

    ! General Species Values
    Species(iSpec)%RHSMethod             = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(hilf2))//'-RHSMethod',''))
    Species(iSpec)%ChargeIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-ChargeIC','0.')
    Species(iSpec)%MassIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-MassIC','0.')
    Species(iSpec)%DensityIC             = GETREAL('Part-Species'//TRIM(hilf2)//'-DensityIC','0.')
    Species(iSpec)%MacroParticleFactor   = GETREAL('Part-Species'//TRIM(hilf2)//'-MacroParticleFactor','1.')
    Species(iSpec)%LowVeloThreshold      = GETREAL('Part-Species'//TRIM(hilf2)//'-LowVeloThreshold','0.')
    Species(iSpec)%HighVeloThreshold     = GETREAL('Part-Species'//TRIM(hilf2)//'-HighVeloThreshold','0.')
    
    ! Bons particle rebound model
    Species(iSpec)%YoungIC               = GETREAL('Part-Species'//TRIM(hilf2)//'-YoungIC')
    Species(iSpec)%PoissonIC             = GETREAL('Part-Species'//TRIM(hilf2)//'-PoissonIC')
    Species(iSpec)%YieldCoeff            = GETREAL('Part-Species'//TRIM(hilf2)//'-YieldCoeff')
    
    ! Random Walk model
#if EQNSYSNR == 4
    Species(iSpec)%RWModel               = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(hilf2))//'-RWModel',''))
#endif

    ! get emission and init data
    Species(iSpec)%Init(iInit)%UseForInit            = GETLOGICAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForInit','.TRUE.')
    Species(iSpec)%Init(iInit)%UseForEmission        = GETLOGICAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-UseForEmission','.TRUE.')
    Species(iSpec)%Init(iInit)%SpaceIC               = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(hilf2))//'-SpaceIC','cuboid'))
    Species(iSpec)%Init(iInit)%velocityDistribution  = TRIM(GETSTR('Part-Species'//TRIM(ADJUSTL(hilf2))//'-velocityDistribution','constant'))
    Species(iSpec)%Init(iInit)%InflowRiseTime        = GETREAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-InflowRiseTime','0.')
    Species(iSpec)%Init(iInit)%initialParticleNumber = GETINT('Part-Species'//TRIM(ADJUSTL(hilf2))//'-initialParticleNumber','0')
    Species(iSpec)%Init(iInit)%RadiusIC              = GETREAL('Part-Species'//TRIM(hilf2)//'-RadiusIC','1.')
    Species(iSpec)%Init(iInit)%NormalIC              = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-NormalIC',3,'0. , 0. , 1.')
    Species(iSpec)%Init(iInit)%BasePointIC           = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BasePointIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVector1IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector1IC',3,'1. , 0. , 0.')
    Species(iSpec)%Init(iInit)%BaseVector2IC         = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-BaseVector2IC',3,'0. , 1. , 0.')
    Species(iSpec)%Init(iInit)%CuboidHeightIC        = GETREAL('Part-Species'//TRIM(hilf2)//'-CuboidHeightIC','1.')
    Species(iSpec)%Init(iInit)%CylinderHeightIC      = GETREAL('Part-Species'//TRIM(hilf2)//'-CylinderHeightIC','1.')
    Species(iSpec)%Init(iInit)%CalcHeightFromDt      = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-CalcHeightFromDt','.FALSE.')
    Species(iSpec)%Init(iInit)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC','0.')
    Species(iSpec)%Init(iInit)%VeloTurbIC            = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloTurbIC','0.')
    Species(iSpec)%Init(iInit)%VeloVecIC             = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3,'0. , 0. , 0.')
    Species(iSpec)%Init(iInit)%PartDensity           = GETREAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-PartDensity','0.')
    Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmissionType','2')
    Species(iSpec)%Init(iInit)%ParticleEmission      = GETREAL('Part-Species'//TRIM(ADJUSTL(hilf2))//'-ParticleEmission','0.')
    Species(iSpec)%Init(iInit)%NumberOfExcludeRegions= GETINT('Part-Species'//TRIM(hilf2)//'-NumberOfExcludeRegions','0')
    Species(iSpec)%Init(iInit)%InsertedParticle      = 0
    Species(iSpec)%Init(iInit)%InsertedParticleSurplus = 0
    
    !----------- various checks/calculations after read-in of Species(i)%Init(iInit)%-data ----------------------------------!
    !--- Check if Initial ParticleInserting is really used
    IF ( ((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.1).OR.(Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2)) &
      .AND.(Species(iSpec)%Init(iInit)%UseForInit) ) THEN
      IF ( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0) &
      .AND. (ABS(Species(iSpec)%Init(iInit)%PartDensity).LE.0.) ) THEN
        Species(iSpec)%Init(iInit)%UseForInit=.FALSE.
        SWRITE(*,*) "WARNING: Initial ParticleInserting disabled as neither ParticleNumber"
        SWRITE(*,*) "nor PartDensity detected for Species, Init ", iSpec, iInit
      END IF
    END IF
    !--- cuboid-/cylinder-height calculation from v and dt
    IF (.NOT.Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CuboidHeightIC,-1.)) THEN ! flag is initialized with -1, compatibility issue 
          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.                 
          SWRITE(*,*) "WARNING: Cuboid height will be calculated from v and dt!"
        END IF
      ELSE IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder') THEN
        IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CylinderHeightIC,-1.)) THEN !flag is initialized with -1, compatibility issue 
          Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.                   
          SWRITE(*,*) "WARNING: Cylinder height will be calculated from v and dt!"
        END IF
      END IF
    END IF
    IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) .AND. (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2) ) &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is only supported for EmiType1 or EmiType2(=default)!')
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cuboid') &
          .AND.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).NE.'cylinder')) &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is only supported for cuboid or cylinder!')
      IF (Species(iSpec)%Init(iInit)%UseForInit) &
        CALL abort(&
__STAMP__&
          ,' Calculating height from v and dt is not supported for initial ParticleInserting!')
    END IF
    !--- virtual pre-insertion (vpi) checks and calculations
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
      .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
      IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) .AND. (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.2) ) &
        CALL abort(&
__STAMP__&
        ,' Wrong emission-type for virtual Pre-Inserting region!')
      IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).NE.'maxwell_lpn') &
        CALL abort(&
__STAMP__&
        ,' Only maxwell_lpn is implemened as velocity-distribution for virtual Pre-Inserting region!')
      IF (Species(iSpec)%Init(iInit)%UseForInit) &
        CALL abort(&
__STAMP__&
          ,' virtual Pre-Inserting is not supported for initial ParticleInserting. Use additional Init!')
      !-- Virtual Pre-Inserting is used correctly !
      Species(iSpec)%Init(iInit)%VirtPreInsert = .TRUE.
      SWRITE(*,*) "Virtual Pre-Inserting is used for Species, Init ", iSpec, iInit
      IF (Species(iSpec)%Init(iInit)%PartDensity .EQ. 0.) THEN
        SWRITE(*,*) "WARNING: If VPI-BC is open, a backflow might not be compensated"
        SWRITE(*,*) "         (use PartDensity instead of ParticleEmission)!"
      END IF
      Species(iSpec)%Init(iInit)%vpiDomainType = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-vpiDomainType','perpendicular_extrusion'))
      SELECT CASE ( TRIM(Species(iSpec)%Init(iInit)%vpiDomainType) )
      CASE ( 'freestream' )
        IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .NE. 'cuboid_vpi' ) THEN
          CALL abort(&
__STAMP__&
            ,' Only cuboid_vpi is supported for a freestream vpiDomainType! (Use default vpiDomainType for cylinder.)')
        ELSE
          Species(iSpec)%Init(iInit)%vpiBVBuffer(1) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV1BufferNeg','.TRUE.')
          Species(iSpec)%Init(iInit)%vpiBVBuffer(2) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV1BufferPos','.TRUE.')
          Species(iSpec)%Init(iInit)%vpiBVBuffer(3) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV2BufferNeg','.TRUE.')
          Species(iSpec)%Init(iInit)%vpiBVBuffer(4) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-vpiBV2BufferPos','.TRUE.')
        END IF
      CASE ( 'orifice' )
        Species(iSpec)%Init(iInit)%vpiBVBuffer = .TRUE.
        IF ( ABS(Species(iSpec)%Init(iInit)%Radius2IC) .GT. 0. ) THEN
          CALL abort(&
__STAMP__&
            ,' Annular orifice is not implemented yet!')
        END IF
      CASE ( 'perpendicular_extrusion' )
        Species(iSpec)%Init(iInit)%vpiBVBuffer = .TRUE. !dummy
      CASE DEFAULT
        CALL abort(&
__STAMP__&
,'vpiDomainType is not implemented!')
      END SELECT
      !--
    ELSE
      Species(iSpec)%Init(iInit)%VirtPreInsert = .FALSE.
    END IF
    !--- integer check for ParticleEmissionType 2
    IF((Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2).AND. &
         ((Species(iSpec)%Init(iInit)%ParticleEmission-INT(Species(iSpec)%Init(iInit)%ParticleEmission)).NE.0)) THEN
       CALL abort(&
__STAMP__&
       ,' If ParticleEmissionType = 2 (parts per iteration), ParticleEmission has to be an integer number')
    END IF
    !--- flag for cell-based constant-pressure-EmiTypes
    !--- normalize VeloVecIC and NormalIC (and BaseVector 1 & 2 IC for cylinder) for Inits
    IF (.NOT. ALL(Species(iSpec)%Init(iInit)%VeloVecIC(:).eq.0.)) THEN
      Species(iSpec)%Init(iInit)%VeloVecIC = Species(iSpec)%Init(iInit)%VeloVecIC            / &
        SQRT(Species(iSpec)%Init(iInit)%VeloVecIC(1)*Species(iSpec)%Init(iInit)%VeloVecIC(1) + &
        Species(iSpec)%Init(iInit)%VeloVecIC(2)*Species(iSpec)%Init(iInit)%VeloVecIC(2)      + &
        Species(iSpec)%Init(iInit)%VeloVecIC(3)*Species(iSpec)%Init(iInit)%VeloVecIC(3))
    END IF
    Species(iSpec)%Init(iInit)%NormalIC = Species(iSpec)%Init(iInit)%NormalIC /                 &
      SQRT(Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1) + &
      Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2) + &
      Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
    IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')&
        .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi')) THEN
        Species(iSpec)%Init(iInit)%BaseVector1IC =&
                  Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector1IC /     &
        SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)*Species(iSpec)%Init(iInit)%BaseVector1IC(1) + &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2)*Species(iSpec)%Init(iInit)%BaseVector1IC(2) + &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3)*Species(iSpec)%Init(iInit)%BaseVector1IC(3))
        Species(iSpec)%Init(iInit)%BaseVector2IC =&
                   Species(iSpec)%Init(iInit)%RadiusIC * Species(iSpec)%Init(iInit)%BaseVector2IC /    &
        SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)*Species(iSpec)%Init(iInit)%BaseVector2IC(1) + &
        Species(iSpec)%Init(iInit)%BaseVector2IC(2)*Species(iSpec)%Init(iInit)%BaseVector2IC(2)      + &
        Species(iSpec)%Init(iInit)%BaseVector2IC(3)*Species(iSpec)%Init(iInit)%BaseVector2IC(3))
    END IF
    !--- read stuff for ExcludeRegions and normalize/calculate corresponding vectors
    IF (Species(iSpec)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      ALLOCATE(Species(iSpec)%Init(iInit)%ExcludeRegion(1:Species(iSpec)%Init(iInit)%NumberOfExcludeRegions)) 
      IF (((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') &
       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) &
      .OR.((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid_vpi') &
       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder_vpi'))) THEN
        DO iExclude=1,Species(iSpec)%Init(iInit)%NumberOfExcludeRegions
          WRITE(UNIT=hilf3,FMT='(I0)') iExclude
          hilf3=TRIM(hilf2)//'-ExcludeRegion'//TRIM(hilf3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC             &
            = TRIM(GETSTR('Part-Species'//TRIM(hilf3)//'-SpaceIC','cuboid'))
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC             &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-RadiusIC','1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC            &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-Radius2IC','0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC             &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-NormalIC',3,'0. , 0. , 1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC          &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BasePointIC',3,'0. , 0. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC        &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector1IC',3,'1. , 0. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC        &
            = GETREALARRAY('Part-Species'//TRIM(hilf3)//'-BaseVector2IC',3,'0. , 1. , 0.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC       &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-CuboidHeightIC','1.')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC     &
            = GETREAL('Part-Species'//TRIM(hilf3)//'-CylinderHeightIC','1.')
          !--normalize and stuff
          IF ((TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).EQ.'cuboid') .OR. &
               ((((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1),1.) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2),0.)) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3),0.)) &
            .OR. ((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1),0.) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2),1.)) &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3),0.))) &
            .AND. (((ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1),0.)) &
              .AND. (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2),0.))) &
              .AND. (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3),1.))))) THEN
            !-- cuboid; or BV are non-default and NormalIC is default: calc. NormalIC for ExcludeRegions from BV1/2
            !   (for def. BV and non-def. NormalIC; or all def. or non-def.: Use User-defined NormalIC when ExclRegion is cylinder)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3) &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)
          ELSE IF ( (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cuboid') .AND. &
                    (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cylinder') )THEN
            CALL abort(&
__STAMP__&
,'Error in ParticleInit, ExcludeRegions must be cuboid or cylinder!')
          END IF
          IF (Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2 + &
              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2 + &
              Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2 .GT. 0.) THEN
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC &
              / SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)**2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1) &
              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)**2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2) &
              = SQRT(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)**2 &
              + Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)**2)
          ELSE
            CALL abort(&
__STAMP__&
,'Error in ParticleInit, NormalIC Vector must not be zero!')
          END IF
        END DO !iExclude
      ELSE
        CALL abort(&
__STAMP__&
,'Error in ParticleInit, ExcludeRegions are currently only implemented for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
      END IF
    END IF
    !--- stuff for calculating ParticleEmission/InitialParticleNumber from PartDensity when this value is not used for LD-stuff
    !                                                                                  (additional not-LD checks might be necassary)
    PartDens_OnlyInit=.FALSE.
    IF ((Species(iSpec)%Init(iInit)%PartDensity.GT.0.)) THEN
      IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.1) THEN
        IF ( (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.2) .AND. (Species(iSpec)%Init(iInit)%UseForInit) ) THEN
          PartDens_OnlyInit=.TRUE.
        ELSE
          CALL abort(&
__STAMP__&
            , 'PartDensity is only supported for EmiType1 or initial ParticleInserting with EmiType1/2!')
        END IF
      END IF
      IF ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid').OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
        IF  ((((TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'constant') &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell') ) &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'maxwell_lpn') ) &
          .OR.(TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).EQ.'emmert') ) THEN
          IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
            CALL abort(&
__STAMP__&
            ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
          END IF
          !---calculation of Base-Area and corresponding component of VeloVecIC
          lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
            Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
          lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
            Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
          lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
            Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
          A_ins = lineVector(1)*lineVector(1) + lineVector(2)*lineVector(2) + lineVector(3)*lineVector(3)
          IF (A_ins .GT. 0.) THEN
            A_ins = SQRT(A_ins)
            lineVector = lineVector / A_ins
            IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
              v_drift_line = Species(iSpec)%Init(iInit)%VeloIC * &
                ( Species(iSpec)%Init(iInit)%VeloVecIC(1)*lineVector(1) + Species(iSpec)%Init(iInit)%VeloVecIC(2)*lineVector(2) &
                + Species(iSpec)%Init(iInit)%VeloVecIC(3)*lineVector(3) ) !lineVector component of drift-velocity
            ELSE
              v_drift_line = 0.
              IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
                PartDens_OnlyInit=.TRUE.
              ELSE
                CALL abort(&
__STAMP__&
                  ,'PartDensity is only supported for CalcHeightFromDt, vpi, or initial ParticleInserting!')
              END IF
            END IF
            IF ( TRIM(Species(iSpec)%Init(iInit)%SpaceIC) .EQ. 'cylinder' ) THEN
              A_ins = Pi * (Species(iSpec)%Init(iInit)%RadiusIC**2-Species(iSpec)%Init(iInit)%Radius2IC**2)
            END IF
            !---calculation of particle flow (macroparticles/s) through boundary
            IF (.NOT.PartDens_OnlyInit) THEN
              Species(iSpec)%Init(iInit)%ParticleEmission &
                = Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor * v_drift_line * A_ins
            END IF
            !---calculation of initial (macro)particle number
            IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
              IF (Species(iSpec)%Init(iInit)%initialParticleNumber .GT. 0) THEN
                CALL abort(&
__STAMP__&
                  ,'Either initialParticleNumber or PartDensity can be defined for selected parameters, not both!')
              END IF
              IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') THEN
                Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
                  * Species(iSpec)%Init(iInit)%CuboidHeightIC * A_ins)
              ELSE !cylinder
                Species(iSpec)%Init(iInit)%initialParticleNumber &
                  = INT(Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor &
                  * Species(iSpec)%Init(iInit)%CylinderHeightIC * A_ins)
              END IF
            END IF
          ELSE
            CALL abort(&
__STAMP__&
              ,'BaseVectors are parallel or zero!')
          END IF
        ELSE
          CALL abort(&
__STAMP__&
          ,'Only const. or maxwell(_lpn) is supported as velocityDistr. for PartDensity without LD!')
        END IF
      ELSE IF (Species(iSpec)%Init(iInit)%VirtPreInsert) THEN
        IF (Species(iSpec)%Init(iInit)%ParticleEmission .GT. 0.) THEN
               CALL abort(&
__STAMP__&
          ,'Either ParticleEmission or PartDensity can be defined for selected emission parameters, not both!')
        ELSE
          SWRITE(*,*) "PartDensity is used for VPI of Species, Init ", iSpec, iInit !Value is calculated inside SetParticlePostion!
        END IF
      ELSE
        CALL abort(&
__STAMP__&
        ,'PartDensity is only supported for the SpaceIC cuboid(_vpi) or cylinder(_vpi)!')
      END IF
    END IF
    !--- determine StartnumberOfInits (start loop index, e.g., for emission loops)
    IF(iInit.EQ.0)THEN
      !!!for new case: check if to be included!!!
      IF((( (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0)&
        .AND.(Species(iSpec)%Init(iInit)%ParticleEmission.EQ.0.) )  &
        .AND.(Species(iSpec)%Init(iInit)%PartDensity.EQ.0.) )       &
        .AND.(Species(iSpec)%Init(iInit)%ConstantPressure.EQ.0.)    &
        .AND.(Species(iSpec)%NumberOfInits.GT.0))       THEN 
        Species(iSpec)%StartnumberOfInits = 1 ! only new style paramaters defined (Part-Species(i)-Init(iInit)-***)
      ELSE
        Species(iSpec)%StartnumberOfInits = 0 ! old style parameters has been defined for inits/emissions (Part-Species(i)-***)
      END IF
      SWRITE(*,*) "StartnumberOfInits of Species ", iSpec, " = ", Species(iSpec)%StartnumberOfInits
    END IF ! iInit .EQ.0

  END DO ! iInit
END DO ! iSpec 

! Read in boundary parameters
! Leave out this check in FLEXI even though we should probably do it
dummy_int  = CountOption('Part-nBounds')       ! check if Part-nBounds is present in .ini file

nPartBound = GETINT('Part-nBounds','1.')  ! get number of particle boundaries
IF ((nPartBound.LE.0).OR.(dummy_int.LT.0)) THEN
  CALL abort(&
__STAMP__&
  ,'ERROR: nPartBound .LE. 0:', nPartBound)
END IF
ALLOCATE(PartBound%SourceBoundName(1:nPartBound))
ALLOCATE(PartBound%TargetBoundCond(1:nPartBound))
ALLOCATE(PartBound%MomentumACC(1:nPartBound))
ALLOCATE(PartBound%WallTemp(1:nPartBound))
ALLOCATE(PartBound%WallVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%WallModel(1:nPartBound))
ALLOCATE(PartBound%WallCoeffModel(1:nPartBound))
ALLOCATE(PartBound%AmbientCondition(1:nPartBound))
ALLOCATE(PartBound%AmbientConditionFix(1:nPartBound))
ALLOCATE(PartBound%AmbientTemp(1:nPartBound))
ALLOCATE(PartBound%AmbientVelo(1:3,1:nPartBound))
ALLOCATE(PartBound%AmbientDens(1:nPartBound))
ALLOCATE(PartBound%SolidState(1:nPartBound))
PartBound%SolidState(1:nPartBound)=.FALSE.

ALLOCATE(PartBound%Adaptive(1:nPartBound))
ALLOCATE(PartBound%AdaptiveType(1:nPartBound))
ALLOCATE(PartBound%AdaptiveTemp(1:nPartBound))
ALLOCATE(PartBound%AdaptivePressure(1:nPartBound))
nAdaptiveBC = 0
PartBound%Adaptive(:) = .FALSE.
PartBound%AdaptiveType(:) = -1
PartBound%AdaptiveTemp(:) = -1.
PartBound%AdaptivePressure(:) = -1.

! Bons particle rebound model
ALLOCATE(PartBound%Young(1:nPartBound))
ALLOCATE(PartBound%Poisson(1:nPartBound))

!--
DO iPartBound=1,nPartBound
  WRITE(UNIT=hilf,FMT='(I0)') iPartBound
  tmpString = TRIM(GETSTR('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-Condition','open'))
  SELECT CASE (TRIM(tmpString))
  CASE('open')
     PartBound%TargetBoundCond(iPartBound) = PartBound%OpenBC          ! definitions see typesdef_pic
     PartBound%AmbientCondition(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientCondition','.FALSE.')
     IF(PartBound%AmbientCondition(iPartBound)) THEN
       PartBound%AmbientConditionFix(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientConditionFix','.TRUE.')
       PartBound%AmbientVelo(1:3,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientVelo',3,'0. , 0. , 0.')
       PartBound%AmbientDens(iPartBound)     = GETREAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientDens','0')
       PartBound%AmbientDynamicVisc(iPartBound)=&
           GETREAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-AmbientDynamicVisc','1.72326582572253E-5') ! N2:T=288K
     END IF
  CASE('reflective')
     PartBound%TargetBoundCond(iPartBound) = PartBound%ReflectiveBC
     PartBound%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-WallVelo',3,'0. , 0. , 0.')
     PartBound%WallModel(iPartBound)       = GETSTR('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-WallModel','perfRef')
     IF (PartBound%WallModel(iPartBound).EQ.'coeffRes') THEN
         PartBound%WallCoeffModel(iPartBound) = GETSTR('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-WallCoeffModel','Tabakoff1981')
         ! Bons particle rebound model
         IF (PartBound%WallCoeffModel(iPartBound).EQ.'Bons2017') THEN
             PartBound%Young(iPartBound)   = GETREAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-Young')
             PartBound%Poisson(iPartBound) = GETREAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-Poisson')
         END IF
     END IF
  CASE('periodic')
     PartBound%TargetBoundCond(iPartBound) = PartBound%PeriodicBC
  CASE DEFAULT
     SWRITE(*,*) ' Boundary does not exists: ', TRIM(tmpString)
     CALL abort(&
__STAMP__&
         ,'Particle Boundary Condition does not exist')
  END SELECT
  PartBound%SourceBoundName(iPartBound) = TRIM(GETSTR('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-SourceName'))
!  PartBound%UseForQCrit(iPartBound) = GETLOGICAL('Part-Boundary'//TRIM(ADJUSTL(hilf))//'-UseForQCrit','.TRUE.')
!  SWRITE(*,*)"PartBound",iPartBound,"is used for the Q-Criterion"
END DO

IF (nMacroRestartFiles.GT.0) THEN
  IF (ALL(.NOT.MacroRestartFileUsed(:))) CALL abort(&
__STAMP__&
,'None of defined Macro-Restart-Files used for any init!')
  DO FileID = 1,nMacroRestartFiles
    IF (.NOT.MacroRestartFileUsed(FileID)) THEN
      SWRITE(*,*) "WARNING: MacroRestartFile: ",FileID," not used for any Init"
    END IF
  END DO
END IF

! Set mapping from field boundary to particle boundary index
ALLOCATE(PartBound%MapToPartBC(1:nBCs))
PartBound%MapToPartBC(:)=-10
DO iPBC=1,nPartBound
  DO iBC = 1, nBCs
    IF (BoundaryType(iBC,1).EQ.0) THEN
      PartBound%MapToPartBC(iBC) = -1 !there are no internal BCs in the mesh, they are just in the name list!
      SWRITE(*,*)"... PartBound",iPBC,"is internal bound, no mapping needed"
    END IF
    IF (TRIM(BoundaryName(iBC)).EQ.TRIM(PartBound%SourceBoundName(iPBC))) THEN
      PartBound%MapToPartBC(iBC) = iPBC !PartBound%TargetBoundCond(iPBC)
      SWRITE(*,*)"... Mapped PartBound",iPBC,"on FieldBound",BoundaryType(iBC,1),",i.e.:",TRIM(BoundaryName(iBC))
    END IF
  END DO
END DO

! Errorhandler for PartBound-Types that could not be mapped to the FieldBound-Types.
DO iBC = 1,nBCs
  IF (PartBound%MapToPartBC(iBC).EQ.-10) THEN
    CALL abort(&
__STAMP__&
    ,' PartBound%MapToPartBC for Boundary is not set. iBC: :',iBC)
  END IF
END DO

ALLOCATE(PEM%Element(1:PDM%maxParticleNumber), PEM%lastElement(1:PDM%maxParticleNumber), STAT=ALLOCSTAT) 
IF (ALLOCSTAT.NE.0) THEN
 CALL abort(&
__STAMP__&
  ,' Cannot allocate PEM arrays!')
END IF

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
      WRITE(UNIT=hilf,FMT='(I0)') iSeed
      seeds(iSeed) = GETINT('Particles-RandomSeed'//TRIM(hilf),'0')
   END DO
END IF

CALL RANDOM_SEED(Size = SeedSize)                ! Check for number of needed Seeds
TrueRandom = .FALSE.                             ! FALSE for defined random seed

! to be stored in HDF5-state file!!!
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

IF (TrueRandom) THEN
! CALL RANDOM_SEED()
! Changed in FLEXI to get true random numbers FOR EACH PROC
   SDEALLOCATE(seeds)
   ALLOCATE(seeds(SeedSize))
   CALL RANDOM_SEED(GET = seeds(1:SeedSize))
#if USE_MPI
   Seeds(1:SeedSize) = Seeds(1:SeedSize)+PartMPI%MyRank
#endif
   CALL RANDOM_SEED(PUT = Seeds(1:SeedSize))
ELSE
#if USE_MPI
   Seeds(1:SeedSize) = Seeds(1:SeedSize)+PartMPI%MyRank
#endif
   CALL RANDOM_SEED(PUT = Seeds(1:SeedSize))
END IF

ALLOCATE(iseeds(SeedSize))
iseeds(:)=0
CALL RANDOM_SEED(GET = iseeds(1:SeedSize))
! to be stored in HDF5-state file!!!
IF(printRandomSeeds)THEN
  IPWRITE(UNIT_StdOut,*) 'Random seeds in PIC_init:'
  DO iSeed = 1,SeedSize
     IPWRITE(UNIT_StdOut,*) iseeds(iSeed)
  END DO
END IF
DEALLOCATE(iseeds)

DelayTime = GETREAL('Part-DelayTime','0.')

!-- Read Flag if warnings to be displayed for rejected velocities when virtual Pre-Inserting region (vpi) is used with PartDensity
OutputVpiWarnings = GETLOGICAL('Particles-OutputVpiWarnings','.FALSE.')

!-- Read Flag if a lost particle should abort the program
AllowLoosing      = GETLOGICAL('Part-AllowLoosing','.FALSE.')

!-- Read Flag if low velocity particles should be removed
LowVeloRemove     = GETLOGICAL('Part-LowVeloRemove','.FALSE.')

!! init interpolation
CALL InitializeInterpolation() ! not any more required ! has to be called earliear
CALL InitPIC()
! always, because you have to construct a halo_eps region around each bc element

#if USE_MPI
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT FIBGM...' 
SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','0')

!-- AuxBCs
nAuxBCs=GETINT('Part-nAuxBCs','0')
IF (nAuxBCs.GT.0) THEN
  UseAuxBCs=.TRUE.
  ALLOCATE (AuxBCType(1:nAuxBCs) &
            ,AuxBCMap(1:nAuxBCs) )
  AuxBCMap=0
  !- Read in BC parameters
  ALLOCATE(PartAuxBC%TargetBoundCond(1:nAuxBCs))
  ALLOCATE(PartAuxBC%MomentumACC(1:nAuxBCs))
  ALLOCATE(PartAuxBC%WallTemp(1:nAuxBCs))
  ALLOCATE(PartAuxBC%Resample(1:nAuxBCs))
  ALLOCATE(PartAuxBC%WallVelo(1:3,1:nAuxBCs))
  ALLOCATE(PartAuxBC%NbrOfSpeciesSwaps(1:nAuxBCs))
  !--determine MaxNbrOfSpeciesSwaps for correct allocation
  MaxNbrOfSpeciesSwaps=0
  DO iPartBound=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iPartBound
    PartAuxBC%NbrOfSpeciesSwaps(iPartBound)= GETINT('Part-AuxBC'//TRIM(hilf)//'-NbrOfSpeciesSwaps','0')
    MaxNbrOfSpeciesSwaps=max(PartAuxBC%NbrOfSpeciesSwaps(iPartBound),MaxNbrOfSpeciesSwaps)
  END DO
  IF (MaxNbrOfSpeciesSwaps.gt.0) THEN
    ALLOCATE(PartAuxBC%ProbOfSpeciesSwaps(1:nAuxBCs))
    ALLOCATE(PartAuxBC%SpeciesSwaps(1:2,1:MaxNbrOfSpeciesSwaps,1:nAuxBCs))
  END IF
  !--
  DO iPartBound=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iPartBound
    tmpString = TRIM(GETSTR('Part-AuxBC'//TRIM(hilf)//'-Condition','open'))
    SELECT CASE (TRIM(tmpString))
    CASE('open')
      PartAuxBC%TargetBoundCond(iPartBound) = PartAuxBC%OpenBC          ! definitions see typesdef_pic
    CASE('reflective')
      PartAuxBC%TargetBoundCond(iPartBound) = PartAuxBC%ReflectiveBC
      PartAuxBC%MomentumACC(iPartBound)     = GETREAL('Part-AuxBC'//TRIM(hilf)//'-MomentumACC','0')
      PartAuxBC%WallTemp(iPartBound)        = GETREAL('Part-AuxBC'//TRIM(hilf)//'-WallTemp','0')
      PartAuxBC%Resample(iPartBound)        = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-Resample','.FALSE.')
      PartAuxBC%WallVelo(1:3,iPartBound)    = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-WallVelo',3,'0. , 0. , 0.')
      IF (PartAuxBC%NbrOfSpeciesSwaps(iPartBound).gt.0) THEN
        !read Species to be changed at wall (in, out), out=0: delete
        PartAuxBC%ProbOfSpeciesSwaps(iPartBound)= GETREAL('Part-AuxBC'//TRIM(hilf)//'-ProbOfSpeciesSwaps','1.')
        DO iSwaps=1,PartAuxBC%NbrOfSpeciesSwaps(iPartBound)
          WRITE(UNIT=hilf2,FMT='(I0)') iSwaps
          PartAuxBC%SpeciesSwaps(1:2,iSwaps,iPartBound) = &
            GETINTARRAY('Part-AuxBC'//TRIM(hilf)//'-SpeciesSwaps'//TRIM(hilf2),2,'0. , 0.')
        END DO
      END IF
    CASE DEFAULT
      SWRITE(*,*) ' AuxBC Condition does not exists: ', TRIM(tmpString)
      CALL abort(&
        __STAMP__&
        ,'AuxBC Condition does not exist')
    END SELECT
  END DO
  !- read and count types
  nAuxBCplanes = 0
  nAuxBCcylinders = 0
  nAuxBCcones = 0
  nAuxBCparabols = 0
  DO iAuxBC=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iAuxBC
    AuxBCType(iAuxBC) = TRIM(GETSTR('Part-AuxBC'//TRIM(hilf)//'-Type','plane'))
    SELECT CASE (TRIM(AuxBCType(iAuxBC)))
    CASE ('plane')
      nAuxBCplanes = nAuxBCplanes + 1
      AuxBCMap(iAuxBC) = nAuxBCplanes
    CASE ('cylinder')
      nAuxBCcylinders = nAuxBCcylinders + 1
      AuxBCMap(iAuxBC) = nAuxBCcylinders
    CASE ('cone')
      nAuxBCcones = nAuxBCcones + 1
      AuxBCMap(iAuxBC) = nAuxBCcones
    CASE ('parabol')
      nAuxBCparabols = nAuxBCparabols + 1
      AuxBCMap(iAuxBC) = nAuxBCparabols
    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(&
        __STAMP__&
        ,'AuxBC does not exist')
    END SELECT
  END DO
  !- allocate type-specifics
  IF (nAuxBCplanes.GT.0) THEN
    ALLOCATE (AuxBC_plane(1:nAuxBCplanes))
  END IF
  IF (nAuxBCcylinders.GT.0) THEN
    ALLOCATE (AuxBC_cylinder(1:nAuxBCcylinders))
  END IF
  IF (nAuxBCcones.GT.0) THEN
    ALLOCATE (AuxBC_cone(1:nAuxBCcones))
  END IF
  IF (nAuxBCparabols.GT.0) THEN
    ALLOCATE (AuxBC_parabol(1:nAuxBCparabols))
  END IF
  !- read type-specifics
  DO iAuxBC=1,nAuxBCs
    WRITE(UNIT=hilf,FMT='(I0)') iAuxBC
    SELECT CASE (TRIM(AuxBCType(iAuxBC)))
    CASE ('plane')
      AuxBC_plane(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_plane(AuxBCMap(iAuxBC))%radius)
      AuxBC_plane(AuxBCMap(iAuxBC))%radius= GETREAL('Part-AuxBC'//TRIM(hilf)//'-radius',TRIM(hilf2))
      n_vec                               = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-n_vec',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-n_vec is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_plane(AuxBCMap(iAuxBC))%n_vec = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
    CASE ('cylinder')
      AuxBC_cylinder(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                                  = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-axis',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_cylinder(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
      AuxBC_cylinder(AuxBCMap(iAuxBC))%radius  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-radius','1.')
      WRITE(UNIT=hilf2,FMT='(G0)') -HUGE(AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmin',TRIM(hilf2))
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cylinder(AuxBCMap(iAuxBC))%lmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmax',TRIM(hilf2))
      AuxBC_cylinder(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-inwards','.TRUE.')
    CASE ('cone')
      AuxBC_cone(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                              = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-axis',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_cone(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
      AuxBC_cone(AuxBCMap(iAuxBC))%lmin  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmin','0.')
      IF (AuxBC_cone(AuxBCMap(iAuxBC))%lmin.LT.0.) CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-lminis .lt. zero for AuxBC',iAuxBC)
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_cone(AuxBCMap(iAuxBC))%lmin)
      AuxBC_cone(AuxBCMap(iAuxBC))%lmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmax',TRIM(hilf2))
      rmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-rmax','0.')
      ! either define rmax at lmax or the halfangle
      IF (rmax.EQ.0.) THEN
        AuxBC_cone(AuxBCMap(iAuxBC))%halfangle  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-halfangle','45.')*PI/180.
      ELSE
        AuxBC_cone(AuxBCMap(iAuxBC))%halfangle  = ATAN(rmax/AuxBC_cone(AuxBCMap(iAuxBC))%lmax)
      END IF
      IF (AuxBC_cone(AuxBCMap(iAuxBC))%halfangle.LE.0.) CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-halfangle is .le. zero for AuxBC',iAuxBC)
      AuxBC_cone(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-inwards','.TRUE.')
      cos2 = COS(AuxBC_cone(AuxBCMap(iAuxBC))%halfangle)**2
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,1) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(1)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/cos2,0.,0./)
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,2) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(2)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/0.,cos2,0./)
      AuxBC_cone(AuxBCMap(iAuxBC))%geomatrix(:,3) &
        = AuxBC_cone(AuxBCMap(iAuxBC))%axis(3)*AuxBC_cone(AuxBCMap(iAuxBC))%axis - (/0.,0.,cos2/)
    CASE ('parabol')
      AuxBC_parabol(AuxBCMap(iAuxBC))%r_vec = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-r_vec',3,'0. , 0. , 0.')
      n_vec                              = GETREALARRAY('Part-AuxBC'//TRIM(hilf)//'-axis',3,'1. , 0. , 0.')
      IF (DOT_PRODUCT(n_vec,n_vec).EQ.0.) THEN
        CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-axis is zero for AuxBC',iAuxBC)
      ELSE !scale vector
        AuxBC_parabol(AuxBCMap(iAuxBC))%axis = n_vec/SQRT(DOT_PRODUCT(n_vec,n_vec))
      END IF
      AuxBC_parabol(AuxBCMap(iAuxBC))%lmin  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmin','0.')
      IF (AuxBC_parabol(AuxBCMap(iAuxBC))%lmin.LT.0.) CALL abort(&
          __STAMP__&
          ,'Part-AuxBC-lmin is .lt. zero for AuxBC',iAuxBC)
      WRITE(UNIT=hilf2,FMT='(G0)') HUGE(AuxBC_parabol(AuxBCMap(iAuxBC))%lmin)
      AuxBC_parabol(AuxBCMap(iAuxBC))%lmax  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-lmax',TRIM(hilf2))
      AuxBC_parabol(AuxBCMap(iAuxBC))%zfac  = GETREAL('Part-AuxBC'//TRIM(hilf)//'-zfac','1.')
      AuxBC_parabol(AuxBCMap(iAuxBC))%inwards = GETLOGICAL('Part-AuxBC'//TRIM(hilf)//'-inwards','.TRUE.')

      n(:,1)=AuxBC_parabol(AuxBCMap(iAuxBC))%axis
      IF (.NOT.ALMOSTZERO(SQRT(n(1,1)**2+n(3,1)**2))) THEN !collinear with y?
        alpha1=ATAN2(n(1,1),n(3,1))
        CALL roty(rot1,alpha1)
        n1=MATMUL(rot1,n)
      ELSE
        alpha1=0.
        CALL ident(rot1)
        n1=n
      END IF
      !print*,'alpha1=',alpha1/PI*180.,'n1=',n1
      IF (.NOT.ALMOSTZERO(SQRT(n1(2,1)**2+n1(3,1)**2))) THEN !collinear with x?
        alpha2=-ATAN2(n1(2,1),n1(3,1))
        CALL rotx(rot2,alpha2)
        n2=MATMUL(rot2,n1)
      ELSE
        CALL abort(&
          __STAMP__&
          ,'vector is collinear with x-axis. this should not be possible... AuxBC:',iAuxBC)
      END IF
      !print*,'alpha2=',alpha2/PI*180.,'n2=',n2
      AuxBC_parabol(AuxBCMap(iAuxBC))%rotmatrix(:,:)=MATMUL(rot2,rot1)
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(:,:)=0.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(1,1)=1.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(2,2)=1.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(3,3)=0.
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(3,4)=-0.5*AuxBC_parabol(AuxBCMap(iAuxBC))%zfac
      AuxBC_parabol(AuxBCMap(iAuxBC))%geomatrix4(4,3)=-0.5*AuxBC_parabol(AuxBCMap(iAuxBC))%zfac
    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(&
        __STAMP__&
        ,'AuxBC does not exist for AuxBC',iAuxBC)
    END SELECT
  END DO
  CALL MarkAuxBCElems()
ELSE
  UseAuxBCs=.FALSE.
END IF

!-- Finalizing InitializeVariables
CALL InitFIBGM()
!!CALL InitSFIBGM()
#if USE_MPI
CALL InitEmissionComm()
#endif /*MPI*/
#if USE_MPI
CALL MPI_BARRIER(PartMPI%COMM,IERROR)
#endif /*MPI*/

SWRITE(UNIT_StdOut,'(132("-"))')

!-- Read parameters for particle-data on region mapping

!-- Read parameters for region mapping
NbrOfRegions = GETINT('NbrOfRegions','0')
IF (NbrOfRegions .GT. 0) THEN
  ALLOCATE(RegionBounds(1:6,1:NbrOfRegions))
  DO iRegions=1,NbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I0)') iRegions
    RegionBounds(1:6,iRegions) = GETREALARRAY('RegionBounds'//TRIM(hilf2),6,'0. , 0. , 0. , 0. , 0. , 0.')
  END DO
END IF

IF (NbrOfRegions .GT. 0) THEN
  CALL MapRegionToElem()
END IF

!-- Floating Potential
ALLOCATE(BCdata_auxSF(1:nPartBound))
DO iPartBound=1,nPartBound
  BCdata_auxSF(iPartBound)%SideNumber=-1 !init value when not used
  BCdata_auxSF(iPartBound)%GlobalArea=0.
  BCdata_auxSF(iPartBound)%LocalArea=0.
END DO
nDataBC_CollectCharges=0

END SUBROUTINE InitializeVariables

!----------------------------------------------------------------------------------------------------------------------------------!
! finalize particle variables
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE FinalizeParticles() 
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!SDEALLOCATE(Pt_temp)
SDEALLOCATE(PartPosRef)
SDEALLOCATE(RandomVec)
SDEALLOCATE(PartState)
SDEALLOCATE(LastPartPos)
SDEALLOCATE(PartSpecies)
SDEALLOCATE(Pt)
SDEALLOCATE(PDM%ParticleInside)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%dtFracPush)
SDEALLOCATE(PDM%IsNewPart)
SDEALLOCATE(Species)
SDEALLOCATE(PartBound%SourceBoundName)
SDEALLOCATE(PartBound%TargetBoundCond)
SDEALLOCATE(PartBound%MomentumACC)
SDEALLOCATE(PartBound%WallTemp)
SDEALLOCATE(PartBound%WallVelo)
SDEALLOCATE(PartBound%AmbientCondition)
SDEALLOCATE(PartBound%AmbientConditionFix)
SDEALLOCATE(PartBound%AmbientTemp)
SDEALLOCATE(PartBound%AmbientVelo)
SDEALLOCATE(PartBound%AmbientDens)
SDEALLOCATE(PEM%Element)
SDEALLOCATE(PEM%lastElement)
SDEALLOCATE(PEM%pStart)
SDEALLOCATE(PEM%pNumber)
SDEALLOCATE(PEM%pEnd)
SDEALLOCATE(PEM%pNext)
END SUBROUTINE FinalizeParticles


SUBROUTINE rotx(mat,a)
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
REAL, INTENT(IN) :: a
mat(:,1)=(/1.0 , 0.     , 0.  /)
mat(:,2)=(/0.0 , COS(a) ,-SIN(a)/)
mat(:,3)=(/0.0 , SIN(a) , COS(a)/)
END SUBROUTINE


SUBROUTINE roty(mat,a)
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
REAL, INTENT(IN) :: a
mat(:,1)=(/ COS(a) , 0., SIN(a)/)
mat(:,2)=(/ 0.     , 1., 0.  /)
mat(:,3)=(/-SIN(a) , 0., COS(a)/)
END SUBROUTINE


SUBROUTINE ident(mat)
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
REAL, INTENT(OUT), DIMENSION(3,3) :: mat
INTEGER :: j
mat = 0.
FORALL(j = 1:3) mat(j,j) = 1.
END SUBROUTINE

END MODULE MOD_ParticleInit
