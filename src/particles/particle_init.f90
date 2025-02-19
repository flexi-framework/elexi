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
#include "flexi.h"
#include "eos.h"
#include "particle.h"

!==================================================================================================================================
! Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_Particle_Init
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE PortabilityGetPID
  FUNCTION GetPID_C() BIND (C, name='getpid')
    !GETPID() is an intrinstic compiler function in gnu. This routine ensures the portability with other compilers.
    USE ISO_C_BINDING,         ONLY: PID_T => C_INT
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    INTEGER(KIND=PID_T)        :: GetPID_C
  END FUNCTION GetPID_C
END INTERFACE

PUBLIC:: DefineParametersParticles
PUBLIC:: InitParticleGlobals
PUBLIC:: InitParticles
PUBLIC:: FinalizeParticles
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particles
!==================================================================================================================================
SUBROUTINE DefineParametersParticles()
! MODULES
USE MOD_ReadInTools
USE MOD_Particle_Analyze,           ONLY:DefineParametersParticleAnalyze,InitParticleAnalyze
USE MOD_Particle_Boundary_Sampling, ONLY:DefineParametersParticleBoundarySampling
USE MOD_Particle_Boundary_Tracking, ONLY:DefineParametersParticleBoundaryTracking
USE MOD_Particle_Globals
USE MOD_Particle_Interpolation,     ONLY:DefineParametersParticleInterpolation
USE MOD_Particle_Mesh,              ONLY:DefineParametersParticleMesh
USE MOD_Particle_Surface_Flux,      ONLY:DefineParametersParticleSurfaceFlux
USE MOD_Particle_Vars
#if PARTICLES_COUPLING >= 2
USE MOD_Particle_Deposition_Method, ONLY:DefineParametersDepositionMethod
#endif /*PARTICLES_COUPLING >= 2*/
#if PARTICLES_COUPLING == 4
USE MOD_Particle_Collision,         ONLY:DefineParametersCollision
#endif /*PARTICLES_COUPLING == 4*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection('Tracking')

CALL prms%CreateIntFromStringOption('TrackingMethod', "Define Method that is used for tracking of particles:\n"                  //&
                                                      "refmapping   (1): reference mapping of particle position with (bi-)linear\n"//&
                                                      "                  and bezier (curved) description of sides.\n"            //&
                                                      "tracing      (2): tracing of particle path with (bi-)linear and bezier\n" //&
                                                      "                  (curved) description of sides.\n"                       //&
                                                      "triatracking (3): tracing of particle path with triangle-aproximation\n"  //&
                                                      "                  of (bi-)linear sides.\n",                                 &
                                                      "triatracking")
CALL addStrListEntry(               'TrackingMethod', 'refmapping'      ,REFMAPPING)
CALL addStrListEntry(               'TrackingMethod', 'tracing'         ,TRACING)
CALL addStrListEntry(               'TrackingMethod', 'triatracking'    ,TRIATRACKING)
CALL addStrListEntry(               'TrackingMethod', 'default'         ,REFMAPPING)

CALL prms%CreateLogicalOption(      'TriaSurfaceFlux'         , 'Using Triangle-aproximation [T] or (bi-)linear and bezier '     //&
                                                                '(curved) description [F] of sides for surfaceflux.'               &
                                                              , '.TRUE.')

CALL prms%CreateLogicalOption(      'CountNbOfLostParts'       , 'Count number of lost particles during tracking that can not '  //&
                                                                 'be found with fallbacks.'                                        &
                                                               , '.FALSE.')
CALL prms%CreateLogicalOption(      'DisplayLostParticles'     , 'Display position, velocity, species and host element of '      //&
                                                                 'particles lost during particle tracking (TrackingMethod = '    //&
                                                                 'triatracking, tracing)'                                          &
                                                               ,'.FALSE.')

CALL prms%CreateIntOption(          'BezierClipLineVectorMethod','TODO-DEFINE-PARAMETER'                                           &
                                                               , '2')
CALL prms%CreateIntOption(          'NbrOfRegions'             , 'Number of regions to be mapped to Elements'                      &
                                                               , '0')


CALL prms%CreateRealOption(         'Part-MaxParticleNumber'   , 'Maximum allowed number of Particles for the whole simulation. '//&
                                                                 'Note that this property is a global value and will be divided '//&
                                                                 'by the number of processors before it is used for '            //&
                                                                 'processor-local array allocation sizes during initialization. '//&
                                                                 '0 corresponds to no limit!')
CALL prms%CreateRealOption(         'Part-MaxPartNumIncrease'  , 'How much shall the PDM%MaxParticleNumber be increased if it '  //&
                                                                 'is full'                                                         &
                                                               , '0.1')
CALL prms%CreateLogicalOption(      'Part-RearrangePartIDs'    , 'Rearrange PartIDs in the process of reducing maxPartNum to '   //&
                                                                 'allow lower memory usage'                                        &
                                                               , '.TRUE.')
CALL prms%CreateIntOption(          'Part-NumberOfRandomSeeds' , 'Number of random seeds for particle random number generator'     &
                                                               , '0')
CALL prms%CreateIntOption(          'Part-RandomSeed[$]'       , 'Seed [$] for Random Number Generator'                            &
                                                               , '1'        , numberedmulti=.TRUE.)

! Timedisc
CALL prms%CreateLogicalOption(      'Part-SteadyState'         , 'Only update particle position while keeping fluid state frozen'  &
                                                               , '.FALSE.')
CALL prms%CreateRealOption(         'Part-ManualTimestep'      , 'Manual time step routine for frozen fluid state'                 &
                                                               , '0.')
CALL prms%CreateLogicalOption(      'Part-LowVeloRemove'       , 'Flag if low velocity particles should be removed'                &
                                                               , '.FALSE.')
CALL prms%CreateLogicalOption(      'Part-WritePartDiam'       , 'Flag to enable writeout of particle diameter'                    &
                                                               , '.FALSE.')
CALL prms%CreateLogicalOption(      'Part-RandomPartDiam'      , 'Flag to enable random particle diameter with a certain variance' &
                                                               , '.FALSE.')
CALL prms%CreateLogicalOption(      'Part-RandomSphericity'    , 'Flag to enable random particle sphericity with a certain var'    &
                                                               , '.FALSE.')

CALL prms%CreateRealArrayOption(    'Part-Gravity'             , 'Gravitational acceleration as vector'                            &
                                                               , '0. , 0. , 0.')
CALL prms%CreateRealOption(         'Part-charactLength'       , 'Characteristic length scale for Stokes number'                   &
                                                                , '1.')
CALL prms%CreateIntOption(          'Part-IniRefState'         , 'IniRefState for the calculation of the Stokes number'            &
                                                                , '1')
CALL prms%CreateIntOption(          'Part-NBasset    '         , 'Number of previous steps used in Basset force'                   &
                                                                , '20')

#if PARABOLIC
#if USE_RW
CALL prms%SetSection("Particle Random Walk")
!===================================================================================================================================
! >>> Values in this section only apply for turbulence models providing turbulent kinetic energy and a turbulent length/time scale
!===================================================================================================================================
CALL prms%CreateStringOption(       'Part-RWModel'  , 'Random walk model used for steady-state calculations.\n'                  //&
                                                      ' - Gosman   : Gosman (1983)\n'                                            //&
                                                      ' - Mofakham : Modified Mofakham (2020)\n'                                 //&
                                                      ' - none     : LES/DNS mode. Assume fully resolved turbulence'               &
                                                    , 'none')
CALL prms%CreateStringOption(       'Part-RWTime'   , 'Time stepping used for random walk model.\n'                              //&
                                                      ' - RK       : Update random walk every RK stage\n'                        //&
                                                      ' - RW       : Update random walk after min(eddy time scale, transit time scale)'&
                                                    , 'RW')
#endif /* USE_RW */

!===================================================================================================================================
! >>> Options for particle SGS model
!===================================================================================================================================
CALL prms%CreateStringOption(       'Part-SGSModel' , 'SGS model used for reconstruction of SGS influence on particle\n'         //&
                                                      ' - Amiri           : Amiri (2006)\n'                                      //&
                                                      ' - Breuer          : Breuer (2017) model, first option\n'                 //&
                                                      ' - Breuer-Analytic : Breuer (2017) model, second option\n'                //&
                                                      ' - Fukagata        : Fukagata (2004)\n'                                   //&
                                                      ' - Jin             : Jin (2010)\n'                                        //&
                                                      ' - Minier          : Minier and Peirano (2004)\n'                         //&
                                                      ' - Sommerfeld      : Sommerfeld (2001)\n'                                 //&
                                                      ' - none            : DNS mode. Assume fully resolved turbulence'            &
                                                    , 'none')
CALL prms%CreateIntOption(          'Part-SGSNFilter','Number of cut-off modes in the high-pass SGS filter'                        &
                                                    , '2')
#endif

!===================================================================================================================================
! > Species
! >>> Values in this section appear multiple times
!===================================================================================================================================
CALL prms%SetSection("Particle Species")
! species inits and properties
CALL prms%CreateIntOption(          'Part-nSpecies'             , 'Number of species in part'                                      &
                                                                , '0')
CALL prms%CreateIntOption(          'Part-Species[$]-nInits'    , 'Number of different initial particle placements for Species [$]'&
                                                                , '0'        , numberedmulti=.TRUE.)
CALL prms%CreateIntFromStringOption('Part-Species[$]-RHSMethod' , 'Particle model used for calculation of the drag force.\n'     //&
                                                                  ' - none          : no coupling between field and particles\n' //&
                                                                  ' - tconvergence  : special case for time-convergence testing\n'//&
                                                                  ' - hpconvergence : special case for hp-convergence testing\n' //&
                                                                  ' - inertia       : particles are convected as initerial particles\n'//&
                                                                  ' - tracer        : particles act as ideal tracers\n'              &
                                                                , 'none'     , numberedmulti=.TRUE.)
CALL addStrListEntry(               'Part-Species[$]-RHSMethod' , 'none',            RHS_NONE)
CALL addStrListEntry(               'Part-Species[$]-RHSMethod' , 'tracer',          RHS_TRACER)
CALL addStrListEntry(               'Part-Species[$]-RHSMethod' , 'tconvergence',    RHS_TCONVERGENCE)
CALL addStrListEntry(               'Part-Species[$]-RHSMethod' , 'hpconvergence',   RHS_HPCONVERGENCE)
CALL addStrListEntry(               'Part-Species[$]-RHSMethod' , 'inertia',         RHS_INERTIA)
CALL addStrListEntry(               'Part-Species[$]-RHSMethod' , 'inertiaeuler',    RHS_INERTIA_EULER)
CALL prms%CreateIntFromStringOption('Part-Species[$]-DragFactor', 'Particle model used for calculation of the drag factor.\n'    //&
                                                                  ' - stokes      : Stokes (1851)\n'                             //&
                                                                  ' - schiller    : Schiller and Naumann (1933)\n'               //&
                                                                  ' - putman      : Putnam et al. (1961)\n'                      //&
                                                                  ' - haider      : Haider and Levenspiel (1989)\n'              //&
                                                                  ' - hoelzer     : Hoelzer et al. (2008)\n'                     //&
                                                                  ' - loth        : Loth (2008)\n'                               //&
                                                                  ' - ganser      : Ganser (2008)\n'                               &
                                                                , 'none'     , numberedmulti=.TRUE.)
CALL addStrListEntry(               'Part-Species[$]-DragFactor', 'stokes',          DF_PART_STOKES)
CALL addStrListEntry(               'Part-Species[$]-DragFactor', 'schiller',        DF_PART_SCHILLER)
CALL addStrListEntry(               'Part-Species[$]-DragFactor', 'putnam',          DF_PART_PUTNAM)
CALL addStrListEntry(               'Part-Species[$]-DragFactor', 'haider',          DF_PART_HAIDER)
CALL addStrListEntry(               'Part-Species[$]-DragFactor', 'hoelzer',         DF_PART_HOELZER)
CALL addStrListEntry(               'Part-Species[$]-DragFactor', 'loth',            DF_PART_LOTH)
CALL addStrListEntry(               'Part-Species[$]-DragFactor', 'ganser',          DF_PART_GANSER)
CALL prms%CreateRealOption(         'Part-Species[$]-MassIC'    , 'Particle mass of species [$] [kg]'                              &
                                                                , '0.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-DiameterIC', 'Particle diameter of species [$] [m]'                           &
                                                                , '0.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-DensityIC' , 'Particle density of species [$] [kg/m^3]'                       &
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-StokesIC'  , 'Particle Stokes number of species [$]'                          &
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(       'Part-Species[$]-velocityDistribution', 'Used velocity distribution.\n'                      //&
                                                                  ' - constant          : Emission with constant velocity\n'     //&
                                                                  ' - constant_turbulent: Emission with constant velocity plus'  //&
                                                                                        ' Gaussian random fluctuations\n'        //&
                                                                  ' - radial_constant   : Emission with constant velocity scaled'//&
                                                                                        ' with radial position\n'                //&
                                                                  ' - load_from_file    : Emission with velocity from file\n'    //&
                                                                  ' - fluid             : Emission with local fluid velocity'      &
                                                                , 'constant', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-VeloIC'    , 'Absolute value of initial velocity. (ensemble velocity) '       &
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-VeloVecIC ', 'Velocity vector for given species'                              &
                                                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-VeloTurbIC', 'Turbulent fluctuation of initial velocity. (ensemble velocity) '&
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-LowVeloThreshold', 'Threshold velocity of particles after reflection.'      //&
                                                                  ' Slower particles are deleted [$] [m/s]'                        &
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-SphericityIC', 'Particle sphericity of species [$] [m]'                       &
                                                                , '1.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-PartDiamVarianceIC', 'Particle diameter variance of species [$] [-]'          &
                                                                , '0.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-ScalePartDiam', 'Scale particle diameter of species [$] [m]'                  &
                                                                , '1.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-PartSpheVarianceIC', 'Particle sphericity variance of species [$] [-]'        &
                                                                , '0.'       , numberedmulti=.TRUE.)
#if USE_PARTTEMP
CALL prms%CreateRealOption(         'Part-Species[$]-SpecificHeatIC'  , 'Particle specific heat [$] [J/kg/K]'                      &
                                                                      , '0.'       , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(       'Part-Species[$]-tempDistribution', 'Particle temperature distribtution. \n'                 //&
                                                                  ' - constant          : Emission with constant velocity\n'     //&
                                                                  ' - fluid             : Emission with local fluid velocity'      &
                                                                , 'constant', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-TemperatureIC'   , 'Absolute value of initial temperature.  '                 &
                                                                , '0.'      , numberedmulti=.TRUE.)
#endif /*USE_PARTTEMP*/
#if USE_EXTEND_RHS
CALL prms%CreateLogicalOption(      'Part-Species[$]-CalcSaffmanForce', 'Flag to calculate the Saffman lift force'                 &
                                                                , '.FALSE.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-CalcMagnusForce', 'Flag to calculate the Magnus force'                        &
                                                                , '.FALSE.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-CalcVirtualMass', 'Flag to calculate the virtual mass force'                  &
                                                                , '.FALSE.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-CalcUndisturbedFlow', 'Flag to calculate the undisturbed flow force'          &
                                                                , '.FALSE.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-CalcBassetForce', 'Flag to calculate the (famous) Basset force'               &
                                                                , '.FALSE.' , numberedmulti=.TRUE.)
#endif /* USE_EXTEND_RHS */


! emission time
CALL prms%SetSection('Particle Species Emission')
CALL prms%CreateLogicalOption(      'Part-Species[$]-UseForEmission', 'Flag to use Init/Emission for emission'                     &
                                                                , '.TRUE.'  , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-UseForInit', 'Flag to use Init/Emission for init'                             &
                                                                , '.TRUE.'  , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-initialParticleNumber', 'Initial particle number'                             &
                                                                , '0'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-InflowRiseTime', 'Time to ramp the number of inflow particles linearly from'//&
                                                                  ' zero to unity'                                                 &
                                                                , '0.'       , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(          'Part-Species[$]-ParticleEmissionType', 'Define Emission Type for particles (volume'         //&
                                                                  ' emission)\n'                                                 //&
                                                                  '1 = emission rate in part/s,\n'                               //&
                                                                  '2 = emission rate in part/iteration\n'                             &
                                                                , '2'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-ParticleEmission', 'Emission rate in part/s or part/iteration.'               &
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-ParticleEmissionTime', 'Scale emission time for EmissionType==1.'             &
                                                                , '1.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-PartDensity', 'PartDensity (real particles per m^3) or (vpi_)cub./cyl.'     //&
                                                                   'as alternative to Part.Emis. in Type1 '                        &
                                                                 , '0.'     , numberedmulti=.TRUE.)

! emission region
CALL prms%CreateStringOption(       'Part-Species[$]-SpaceIC'   , 'Specifying Keyword for particle space condition of species '  //&
                                                                  '[$] in case of one init.\n'                                   //&
                                                                  ' - point\n'                                                   //&
                                                                  ' - line_with_equidistant_distribution\n'                      //&
                                                                  ' - line\n'                                                    //&
                                                                  ' - Gaussian\n'                                                //&
                                                                  ' - plane\n'                                                   //&
                                                                  ' - disc\n'                                                    //&
                                                                  ' - cross\n'                                                   //&
                                                                  ' - circle\n'                                                  //&
                                                                  ' - circle_equidistant\n'                                      //&
                                                                  ' - cuboid\n'                                                  //&
                                                                  ' - cylinder\n'                                                //&
                                                                  ' - sphere\n'                                                  //&
                                                                  ' - load_from_file\n'                                            &
                                                                , 'cuboid'  , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-BasePointIC','Base point for IC'                                              &
                                                                 , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-BaseVector1IC','First base vector for IC'                                     &
                                                                 , '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-BaseVector2IC', 'Second base vector for IC'                                   &
                                                                 , '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-BaseVariance','Variance for Gaussian distribtution'                           &
                                                                 ,'1.'           , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-NormalIC'   , 'Normal orientation of IC.'                                     &
                                                                 , '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-CalcHeightFromDt', 'Calculated cuboid/cylinder height from v and dt'          &
                                                                , '.FALSE.'  , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-CuboidHeightIC'  , 'Height of cuboid if SpaceIC=cuboid'                       &
                                                                , '1.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-CylinderHeightIC', 'Third measure of cylinder  (set 0 for flat rectangle),' //&
                                                                  ' negative value = opposite direction'                           &
                                                                , '1.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-RadiusIC'  , 'Radius for IC'                                                  &
                                                                , '1.'       , numberedmulti=.TRUE.)

! Exclude regions
CALL prms%CreateIntOption(          'Part-Species[$]-NumberOfExcludeRegions', 'Number of different regions to be excluded'         &
                                                                , '0'       , numberedmulti=.TRUE.)

CALL prms%SetSection("Particle Species nInits")
! if nInit > 0 some variables have to be defined twice
CALL prms%CreateStringOption(       'Part-Species[$]-Init[$]-velocityDistribution', 'Used velocity distribution.\n'              //&
                                                                  ' - constant          : Emission with constant velocity\n'     //&
                                                                  ' - constant_turbulent: Emission with constant velocity plus'  //&
                                                                                        ' Gaussian random fluctuations\n'        //&
                                                                  ' - radial_constant   : Emission with constant velocity scaled'//&
                                                                                        ' with radial position\n'                //&
                                                                  ' - load_from_file    : Emission with velocity from file\n'    //&
                                                                  ' - fluid             : Emission with local fluid velocity'      &
                                                                , 'constant', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-VeloIC'    , 'Absolute value of initial velocity. (ensemble velocity)'&
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-Init[$]-VeloVecIC ', 'Velocity vector for given species'                      &
                                                                , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-VeloTurbIC', 'Turbulent fluctuation of initial velocity. (ensemble velocity)'&
                                                                , '0.'      , numberedmulti=.TRUE.)
#if USE_PARTTEMP
CALL prms%CreateStringOption(       'Part-Species[$]-Init[$]-tempDistribution', 'Particle temperature distribtution. \n'         //&
                                                                  ' - constant          : Emission with constant velocity\n'     //&
                                                                  ' - fluid             : Emission with local fluid velocity'      &
                                                                , 'constant', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-TemperatureIC'   , 'Absolute value of initial temperature.  '         &
                                                                , '0.'      , numberedmulti=.TRUE.)
#endif /*USE_PARTTEMP*/

! emission time
CALL prms%SetSection('Particle Species nInits Emission')
CALL prms%CreateLogicalOption(      'Part-Species[$]-Init[$]-UseForEmission', 'Flag to use Init/Emission for emission'             &
                                                                , '.TRUE.'  , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-Init[$]-UseForInit', 'Flag to use Init/Emission for init'                     &
                                                                , '.TRUE.'  , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-initialParticleNumber', 'Initial particle number'                     &
                                                                , '0'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-InflowRiseTime', 'Time to ramp the number of inflow particles linearly from' //&
                                                                  ' zero to unity'                                                 &
                                                                , '0.'       , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(          'Part-Species[$]-Init[$]-ParticleEmissionType', 'Define Emission Type for particles (volume' //&
                                                                  ' emission)\n'                                                 //&
                                                                  '1 = emission rate in part/s,\n'                               //&
                                                                  '2 = emission rate in part/iteration\n'                             &
                                                                , '1'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-ParticleEmission', 'Emission rate in part/s or part/iteration.'       &
                                                                , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-ParticleEmissionTime', 'Scale emission time for EmissionType==1.'     &
                                                                , '1.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-PartDensity', 'PartDensity (real particles per m^3) or (vpi_)cub./cyl.'//&
                                                                   ' as alternative to Part.Emis. in Type1 '                        &
                                                                 , '0.'     , numberedmulti=.TRUE.)

! emission region
CALL prms%CreateStringOption(       'Part-Species[$]-Init[$]-SpaceIC'   , 'Specifying Keyword for particle space condition of species'//&
                                                                  ' [$] in case of one init.\n'                                   //&
                                                                  ' - point\n'                                                   //&
                                                                  ' - line_with_equidistant_distribution\n'                      //&
                                                                  ' - line\n'                                                    //&
                                                                  ' - Gaussian\n'                                                //&
                                                                  ' - plane\n'                                                   //&
                                                                  ' - disc\n'                                                    //&
                                                                  ' - cross\n'                                                   //&
                                                                  ' - circle\n'                                                  //&
                                                                  ' - circle_equidistant\n'                                      //&
                                                                  ' - cuboid\n'                                                  //&
                                                                  ' - cylinder\n'                                                //&
                                                                  ' - sphere\n'                                                  //&
                                                                  ' - load_from_file\n'                                            &
                                                                , 'cuboid'  , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-Init[$]-BasePointIC','Base point for IC'                                      &
                                                                 , '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-Init[$]-BaseVector1IC','First base vector for IC'                             &
                                                                 , '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-Init[$]-BaseVector2IC', 'Second base vector for IC'                           &
                                                                 , '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-BaseVariance','Variance for Gaussian distribtution'                   &
                                                                 ,'1.'           , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Species[$]-Init[$]-NormalIC'   , 'Normal orientation of IC.'                             &
                                                                 , '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-Init[$]-CalcHeightFromDt', 'Calculated cuboid/cylinder height from v and dt'  &
                                                                , '.FALSE.'  , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-CuboidHeightIC'  , 'Height of cuboid if SpaceIC=cuboid'               &
                                                                , '1.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-CylinderHeightIC', 'Third measure of cylinder  (set 0 for flat rectangle),'//&
                                                                  ' negative value = opposite direction'                           &
                                                                , '1.'       , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Init[$]-RadiusIC'  , 'Radius for IC'                                          &
                                                                , '1.'       , numberedmulti=.TRUE.)

! Exclude regions
CALL prms%CreateIntOption(          'Part-Species[$]-Init[$]-NumberOfExcludeRegions', 'Number of different regions to be excluded' &
                                                                , '0'       , numberedmulti=.TRUE.)

! Surface Flux

CALL prms%CreateLogicalOption(      'OutputSurfaceFluxLinked'         , 'Flag to print the SurfaceFlux-linked Info'                &
                                                                      , '.FALSE.')

!===================================================================================================================================
! > Boundaries
!===================================================================================================================================
CALL prms%SetSection("Particle Boundaries")
CALL prms%CreateStringOption(       'Part-Boundary[$]-Type'     , 'Used boundary condition for boundary.\n'                      //&
                                                                  '- open\n'                                                     //&
                                                                  '- reflective\n'                                               //&
                                                                  '- periodic\n'                                                   &
                                                                            , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(       'Part-Boundary[$]-Name'     , 'Source name of boundary. Has to be same name as defined in'   //&
                                                                ' preproc tool'                                                    &
                                                                            , numberedmulti=.TRUE.)

! Wall model =======================================================================================================================
CALL prms%SetSection("Particle Rebound Model")
CALL prms%CreateStringOption(       'Part-Boundary[$]-WallModel', 'Wall model to be used. Options:.\n'                           //&
                                                                  ' - perfRef      : perfect reflection\n'                       //&
                                                                  ' - coeffRes     : Coefficient of restitution'                   &
                                                                  ,'perfRef', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(       'Part-Boundary[$]-WallCoeffModel', 'Coefficients to be used. Options:.\n'                    //&
                                                                  ' - Grant1975    : Grant and Tabakoff (1975)\n'                //&
                                                                  ' - Tabakoff1981 : Tabaoff and Wakeman (1981)\n'               //&
                                                                  ' - Bons2017     : Bons et al. (2017)\n'                       //&
                                                                  ' - Whitaker2018 : Whitaker and Bons (2018)\n'                 //&
                                                                  ' - Fong2019     : Fong et al. (2019)\n'                       //&
                                                                  ' - RebANN       : Rebound  artificial neural network\n'       //&
                                                                  ' - FracANN      : Fracture artificial neural network'           &
                                                                            , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Boundary[$]-RoughWall'  ,'Rough wall modelling is used.'                                 &
                                                                  ,'.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Boundary[$]-RoughMeanIC', 'Mean of Gaussian dist..'                                      &
                                                                  ,'0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Boundary[$]-RoughVarianceIC', 'Mean of Gaussian dist..'                                  &
                                                                  ,'0.04', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Boundary[$]-Young'      , "Young's modulus defining stiffness of wall material"          &
                                                                  , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Boundary[$]-Poisson'    , "Poisson ratio defining relation of transverse to axial strain"&
                                                                  , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Boundary[$]-FricCoeff'  , "Friction coeff"                                               &
                                                                  , '0.3'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Boundary[$]-CoR'        , "Coefficent of restitution for normal velocity component"      &
                                                                  , '1.'      , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(       'Part-Boundary[$]-ANNModel'   , 'Specifying file which contains the weights of an MLP.'        &
                                                                  , 'model/weights.dat', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-YoungIC'     , "Young's modulus of particle defining stiffness of particle material"       &
                                                                  , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-Whitaker_a'  , "Scale factor which defines Young's modulus of particle"       &
                                                                  , '1.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-PoissonIC'   , "Poisson ratio of particle defining relation of transverse to axial strain" &
                                                                  , '0.'      , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(         'Part-Species[$]-YieldCoeff'  , "Yield strength defining elastic deformation"                  &
                                                                  , '0.'       , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(      'Part-Species[$]-RoughWall'   , "Enables rough wall modelling"                                 &
                                                                  ,'.TRUE.'   , numberedmulti=.TRUE.)

! Ambient condition ================================================================================================================
! CALL prms%CreateLogicalOption(      'Part-Boundary[$]-AmbientCondition', 'Use ambient condition (condition "behind" boundary).'    &
!                                                                             , numberedmulti=.TRUE.)
! CALL prms%CreateLogicalOption(      'Part-Boundary[$]-AmbientConditionFix', 'TODO-DEFINE-PARAMETER'                                &
!                                                                             , numberedmulti=.TRUE.)
! CALL prms%CreateRealArrayOption(    'Part-Boundary[$]-AmbientVelo', 'Ambient velocity'                                             &
!                                                                   ,'0. , 0. , 0.', , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(    'Part-Boundary[$]-WallVelo'   , 'Velocity (global x,y,z in [m/s]) of reflective particle'    //&
                                                                    ' boundary.'                                                   &
                                                                  ,'0. , 0. , 0.', numberedmulti=.TRUE.)

! CALL prms%CreateRealOption(         'Part-Boundary[$]-AmbientDens', 'Ambient density'                                              &
!                                                                             , numberedmulti=.TRUE.)
! CALL prms%CreateRealOption(         'Part-Boundary[$]-AmbientDynamicVisc' , 'Ambient dynamic viscosity'                            &
!                                                                             , numberedmulti=.TRUE.)
!#if USE_EXTEND_RHS && ANALYZE_RHS
!CALL prms%CreateRealOption(         'Part-tWriteRHS'              , 'Output time for RHS', '0.'                                    )
!#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

! Call every other DefineParametersParticle routine
#if PARTICLES_COUPLING >= 2
CALL DefineParametersDepositionMethod()
#endif /*PARTICLES_COUPLING >= 2*/
#if PARTICLES_COUPLING == 4
CALL DefineParametersCollision()
#endif /*PARTICLES_COUPLING == 4*/
CALL DefineParametersParticleAnalyze()
CALL DefineParametersParticleBoundarySampling()
CALL DefineParametersParticleBoundaryTracking()
CALL DefineParametersParticleInterpolation()
CALL DefineParametersParticleMesh()
CALL DefineParametersParticleSurfaceFlux()

END SUBROUTINE DefineParametersParticles


!===================================================================================================================================
! Global particle parameters needed for other particle inits
!===================================================================================================================================
SUBROUTINE InitParticleGlobals()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools,                ONLY: GETREAL,GETLOGICAL,GETINTFROMSTR,CountOption
USE MOD_Particle_Interpolation_Vars,ONLY: DoInterpolation
USE MOD_Particle_Tracking_Vars,     ONLY: TrackingMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Read basic particle parameter
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GLOBALS...'

! Find tracking method immediately, a lot of the later variables depend on it
TrackingMethod  = GETINTFROMSTR('TrackingMethod')
SELECT CASE(TrackingMethod)
  CASE(TRIATRACKING,TRACING,REFMAPPING)
    ! Valid tracking method, do nothing
  CASE DEFAULT
    SWRITE(UNIT_stdOut,'(A)')' TrackingMethod not implemented! Select refmapping (1), tracing (2) or triatracking (3).'
    CALL CollectiveStop(__STAMP__,'TrackingMethod not implemented! TrackingMethod=',IntInfo=TrackingMethod)
END SELECT

DoInterpolation       = GETLOGICAL('Part-DoInterpolation')
! IF (.NOT.DoInterpolation) &
!   CALL CollectiveStop(__STAMP__,'Simulation without particle interpolation currently not supported!')

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE GLOBALS DONE!'

END SUBROUTINE InitParticleGlobals


!===================================================================================================================================
! Glue Subroutine for particle initialization
!===================================================================================================================================
!SUBROUTINE InitParticles(ManualTimeStep_opt)
SUBROUTINE InitParticles(doLoadBalance_opt)
! MODULES
USE MOD_Globals
USE Mod_Particle_Globals
USE MOD_ReadInTools
USE MOD_IO_HDF5,                    ONLY: AddToElemData
USE MOD_Particle_Tools,             ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze,           ONLY: InitParticleAnalyze
USE MOD_Particle_Boundary_Sampling, ONLY: RestartParticleBoundarySampling
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Emission,          ONLY: InitializeParticleEmission
USE MOD_Particle_Mesh_Build,        ONLY: InitElemVolumes
USE MOD_Particle_Mesh_Vars,         ONLY: LocalVolume,MeshVolume
USE MOD_Particle_Restart,           ONLY: ParticleRestart
USE MOD_Particle_Surfaces,          ONLY: InitParticleSurfaces
USE MOD_Particle_Surface_Flux,      ONLY: InitializeParticleSurfaceFlux
USE MOD_Particle_TimeDisc,          ONLY: Particle_InitTimeDisc
USE MOD_Particle_Tracking_Vars,     ONLY: TrackingMethod
USE MOD_Particle_Vars,              ONLY: ParticlesInitIsDone,nSpecies,PDM
USE MOD_ReadInTools,                ONLY: PrintOption
#if USE_LOADBALANCE
USE MOD_LoadBalance,                ONLY:InitLoadBalanceTracking
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_Particle_MPI,               ONLY: InitParticleCommSize
#endif /*USE_MPI*/
#if PARABOLIC
#if USE_RW
USE MOD_Particle_RandomWalk,        ONLY: ParticleInitRandomWalk
#endif /*USE_RW*/
USE MOD_Particle_SGS,               ONLY: ParticleInitSGS
#endif /*PARABOLIC*/
#if PARTICLES_COUPLING >= 2
USE MOD_Particle_Vars,              ONLY: doCalcPartSource
USE MOD_Particle_Deposition        ,ONLY: InitializeDeposition
#endif /*PARTICLES_COUPLING >= 2*/
#if PARTICLES_COUPLING == 4
USE MOD_Particle_Collision         ,ONLY: InitializeCollision
USE MOD_Particle_Vars,              ONLY: doCalcPartCollision
#endif /*PARTICLES_COUPLING == 4*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN),OPTIONAL       :: ManualTimeStep_opt                                             !> ManualTimeStep coming from Posti
LOGICAL,INTENT(IN),OPTIONAL      :: doLoadBalance_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: maxParticleNumberGlobal,maxParticleNumberUniform
INTEGER           :: minParticleNumberLocal,maxParticleNumberLocal,sumParticleNumberLocal
CHARACTER(32)     :: hilf
!===================================================================================================================================

IF(ParticlesInitIsDone)THEN
   SWRITE(UNIT_stdOut,'(A)') "InitParticles already called."
   RETURN
END IF

SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLES...'

CALL InitElemVolumes()

! Read global number of particles as datatype real to allow number inputs in the format of, e.g., 5e6
SELECT CASE(CountOption('Part-MaxParticleNumber'))
  CASE(1)
    WRITE(UNIT=hilf,FMT='(I0)') HUGE(PDM%maxAllowedParticleNumber)
    maxParticleNumberGlobal  = GETREAL('Part-MaxParticleNumber')
    PDM%MaxPartNumIncrease   = GETREAL('Part-MaxPartNumIncrease','0.1')
    PDM%RearrangePartIDs     = GETLOGICAL('Part-RearrangePartIDs','.TRUE.')
    ! Divide by number of processors, but use at least 2 (at least one next free position in the case of MPI)
    maxParticleNumberUniform = MAX(maxParticleNumberGlobal/nProcessors,2.)
    ! Increase the max particle number according to the local volume to account for element size differences
    IF (maxParticleNumberUniform * MAX(1.,LocalVolume/(MeshVolume/nProcessors)).GT. HUGE(INT(1,KIND=4)) .OR. &
        maxParticleNumberUniform * MAX(1.,LocalVolume/(MeshVolume/nProcessors)).LT.-HUGE(INT(1,KIND=4)))     &
      CALL CollectiveStop(__STAMP__,'maxParticleNumber too big for current number of processors. Decrease maxParticleNumber or increase nProcs!')
    PDM%maxAllowedParticleNumber = INT(maxParticleNumberUniform * MAX(1.,LocalVolume/(MeshVolume/nProcessors)))

    minParticleNumberLocal = PDM%maxAllowedParticleNumber
    maxParticleNumberLocal = PDM%maxAllowedParticleNumber
    sumParticleNumberLocal = PDM%maxAllowedParticleNumber
#if USE_MPI
    IF (MPIRoot) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,minParticleNumberLocal,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_FLEXI,iERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,maxParticleNumberLocal,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_FLEXI,iERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,sumParticleNumberLocal,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
    ELSE
      CALL MPI_REDUCE(minParticleNumberLocal,0           ,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_FLEXI,iERROR)
      CALL MPI_REDUCE(maxParticleNumberLocal,0           ,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_FLEXI,iERROR)
      CALL MPI_REDUCE(sumParticleNumberLocal,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
    END IF
#endif /*USE_MPI*/
    CALL PrintOption('Particle NUMBER (Min/Max/Glob)','CALC',IntArrayOpt=(/minParticleNumberLocal &
                                                                          ,maxParticleNumberLocal &
                                                                          ,sumParticleNumberLocal/))
    PDM%maxParticleNumber    = 1
    nGlobalNbrOfParticles    = 0
    nGlobalNbrOfParticles(4) = HUGE(nGlobalNbrOfParticles(4))
  CASE(0)
    PDM%maxParticleNumber    = 0
    nGlobalNbrOfParticles    = 0
    nGlobalNbrOfParticles(4) = HUGE(nGlobalNbrOfParticles(4))
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'Part-maxParticleNumber must only be specified once!')
END SELECT

IF(TrackingMethod.NE.TRIATRACKING) THEN
  CALL InitParticleSurfaces()
END IF

!CALL InitializeVariables(ManualTimeStep_opt)
CALL InitializeVariables()
! Set particle timedisc pointer
CALL Particle_InitTimeDisc()

#if USE_LOADBALANCE
CALL InitLoadBalanceTracking()
#endif /*USE_LOADBALANCE*/

! Return if running particle code without any species
IF (nSpecies.LE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLES DONE!'
  SWRITE(UNIT_stdOut,'(132("-"))')
  RETURN
END IF

! InitRandomWalk must be called after InitializeVariables to know the size of TurbPartState
#if PARABOLIC
#if USE_RW
CALL ParticleInitRandomWalk()
#endif /*USE_RW*/
! InitSGS must be called after InitRandomWalk and will abort FLEXI if an RW model is defined
CALL ParticleInitSGS()
#endif /*PARABOLIC*/

! Requires information about initialized variables
CALL InitParticleAnalyze()

! Restart particles here, otherwise we can not know if we need to have an initial emission
IF (PRESENT(doLoadBalance_opt)) THEN
  PDM%ParticleVecLength = 0
  CALL UpdateNextFreePosition()

  CALL ParticleRestart(doFlushFiles=.FALSE.)
ELSE
  CALL ParticleRestart()
END IF
! Initialize emission. If no particles are present, assume restart from pure fluid and perform initial inserting
CALL InitializeParticleEmission()

#if USE_MPI
! has to be called AFTER InitializeVariables because we need to read the parameter file to know the CommSize
CALL InitParticleCommSize()
#endif /*USE_MPI*/

#if PARTICLES_COUPLING >= 2
! Fluid-particle coupling
doCalcPartSource        = GETLOGICAL('Part-CalcSource'    )
IF (doCalcPartSource)    CALL InitializeDeposition()
#endif /*PARTICLES_COUPLING >= 2*/

#if PARTICLES_COUPLING == 4
! Particle-particle coupling
doCalcPartCollision     = GETLOGICAL('Part-CalcCollision'    )
IF (doCalcPartCollision) CALL InitializeCollision()
#endif /*PARTICLES_COUPLING == 4*/

! Rebuild previous sampling on surfMesh
CALL RestartParticleBoundarySampling()

! Initialize surface flux
CALL InitializeParticleSurfaceFlux()

ParticlesInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLES DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitParticles


!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
!SUBROUTINE InitializeVariables(ManualTimeStep_opt)
SUBROUTINE InitializeVariables()
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Analyze_Vars      ,ONLY: RPP_MaxBufferSize, RPP_Plane, RecordPart, RPP_nVarNames
USE MOD_Particle_Analyze_Vars      ,ONLY: RPP_Records,RPP_Records_Glob
USE MOD_Particle_Boundary_Sampling ,ONLY: InitParticleBoundarySampling
USE MOD_Particle_Boundary_Tracking ,ONLY: InitParticleBoundaryTracking
USE MOD_Particle_Boundary_Vars     ,ONLY: LowVeloRemove,doParticleReflectionTrack
USE MOD_Particle_Interpolation     ,ONLY: InitParticleInterpolation
USE MOD_Particle_Interpolation_Vars,ONLY: DoInterpolation
USE MOD_Particle_Mesh              ,ONLY: InitParticleMesh
USE MOD_ReadInTools
#if USE_MPI
USE MOD_Particle_Analyze_Vars      ,ONLY: RPP_MPI_Request
USE MOD_Particle_MPI_Emission      ,ONLY: InitEmissionComm
USE MOD_Particle_MPI_Halo          ,ONLY: IdentifyPartExchangeProcs
#endif /*USE_MPI*/
!#if USE_EXTEND_RHS && ANALYZE_RHS
!USE MOD_Output_Vars                ,ONLY: ProjectName
!USE MOD_Output                     ,ONLY: InitOutputToFile
!USE MOD_Restart_Vars               ,ONLY: DoRestart,RestartTime
!USE MOD_Analyze_Vars               ,ONLY: analyze_dt
!#endif /* USE_EXTEND_RHS && ANALYZE_RHS */
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN),OPTIONAL       :: ManualTimeStep_opt                                             !> ManualTimeStep coming from Posti
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: RPP_maxMemory, jP
CHARACTER(30)         :: tmpStr
!===================================================================================================================================
doPartIndex             = GETLOGICAL('Part-PartIndex'     )
IF(doPartIndex) sumOfMatchedParticlesSpecies = 0

doWritePartDiam         = GETLOGICAL('Part-WritePartDiam' )
doRandomPartDiam        = GETLOGICAL('Part-RandomPartDiam')
#if USE_SPHERICITY
doRandomSphericity      = GETLOGICAL('Part-RandomSphericity')
#endif

CALL AllocateParticleArrays()
CALL InitializeVariablesRandomNumbers()

! gravitational acceleration
PartGravity             = GETREALARRAY('Part-Gravity'      ,3  )
FilenameRecordPart      = GETSTR('Part-FilenameRecordPart'     )

! Number of species
nSpecies                = GETINT(     'Part-nSpecies')
! Return if running particle code without any species
IF (nSpecies.LE.0) THEN
  DoInterpolation = .FALSE.
  RETURN
END IF

! Check if maxParticleNumber is set
IF (PDM%maxParticleNumber.LE.0) &
  CALL CollectiveStop(__STAMP__,'nSpecies > 0 but maxParticleNumber <= 0!')

! Allocate species array
ALLOCATE(Species(1:nSpecies))
CALL InitializeVariablesSpeciesInits()

CALL InitializeVariablesPartBoundary()
! Flag if low velocity particles should be removed
LowVeloRemove       = GETLOGICAL('Part-LowVeloRemove')

! Initialize record plane of particles
RecordPart          = GETINT('Part-nRPs')
IF (RecordPart.GT.0) THEN
  RPP_nVarNames = 9
  IF(doPartIndex)               RPP_nVarNames = RPP_nVarNames + 1
  IF(doParticleReflectionTrack) RPP_nVarNames = RPP_nVarNames + 1
#if USE_SPHERICITY
  RPP_nVarNames = RPP_nVarNames + 1
#endif
#if ANALYZE_RHS
  RPP_nVarNames = RPP_nVarNames + 18
#endif
  CALL SYSTEM('mkdir -p recordpoints')
  ! Get size of buffer array
  RPP_maxMemory     = GETINT('Part-RPMemory')       ! Max buffer (100MB)
  RPP_MaxBufferSize = RPP_MaxMemory*131072/6        != size in bytes/(real*RPP_maxMemory)
  ! Readin record planes
  ALLOCATE(RPP_Plane(       RecordPart))
  ALLOCATE(RPP_Records(     RecordPart))
  ALLOCATE(RPP_Records_Glob(RecordPart))
  RPP_Records_Glob = 0
#if USE_MPI
  RPP_MPI_Request = MPI_REQUEST_NULL
#endif /*USE_MPI*/
  DO jP = 1,RecordPart
    WRITE(UNIT=tmpStr,FMT='(I2)') jP
    ALLOCATE(RPP_Plane(jP)%RPP_Data(RPP_nVarNames,RPP_MaxBufferSize))
    RPP_Plane(jP)%RPP_Data = 0.
    RPP_Plane(jP)%pos  = GETREALARRAY('Part-RPBasePoint'//TRIM(ADJUSTL(tmpStr)),3)
    RPP_Plane(jP)%dir  = GETREALARRAY('Part-RPNormVec'  //TRIM(ADJUSTL(tmpStr)),3)
    RPP_Plane(jp)%dist = DOT_PRODUCT(RPP_Plane(jP)%pos,RPP_Plane(jp)%dir)
    RPP_Plane(jP)%RPP_Records = 0.
  END DO
END IF

! CALL InitializeVariablesTimeStep(ManualTimeStep_opt)
CALL InitializeVariablesTimeStep()

! Build BGM and initialize particle mesh
CALL InitParticleMesh()
#if USE_MPI
!-- Build MPI communication
CALL IdentifyPartExchangeProcs()
#endif

! Initialize surface sampling
CALL InitParticleBoundarySampling()

! Initialize impact recording
CALL InitParticleBoundaryTracking()

! Initialize interpolation and particle-in-cell for field -> particle coupling
!--> Could not be called earlier because a halo region has to be build depending on the given BCs
CALL InitParticleInterpolation()

! Initialize MPI communicator for emmission procs
#if USE_MPI
CALL InitEmissionComm()
CALL MPI_BARRIER(MPI_COMM_FLEXI,IERROR)
#endif /*MPI*/

!#if USE_EXTEND_RHS && ANALYZE_RHS
!WRITE(tmpStr,FMT='(I0)') myRank
!FileName_RHS = TRIM(ProjectName)//'_'//TRIM(ADJUSTL(tmpStr))//'_RHS'
!WRITE(tmpStr,FMT='(E8.2)') analyze_dt
!dtWriteRHS    = GETREAL('Part-tWriteRHS',TRIM(ADJUSTL(tmpStr)))
!IF(dtWriteRHS.GT.0.0)THEN
!  tWriteRHS     = MERGE(dtWriteRHS+RestartTime,dtWriteRHS,doRestart)
!
!  CALL InitOutputToFile(FileName_RHS,'RHS',23,&
!    [CHARACTER(6)::"Spec","Fx","Fy","Fz","Fdmx","Fdmy","Fdmz","Flmx","Flmy","Flmz","Fmmx","Fmmy","Fmmz",&
!                   "Fumx","Fumy","Fumz","Fvmx","Fvmy","Fvmz","Fbmx","Fbmy","Fbmz","nIndex"],WriteRootOnly=.FALSE.)
!END IF
!#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitializeVariables


SUBROUTINE AllocateParticleArrays()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Memory                 ,ONLY: VerifyMemUsage
USE MOD_ReadInTools
USE MOD_Particle_Vars
#if USE_BASSETFORCE
USE MOD_Equation_Vars          ,ONLY: s43
USE MOD_Restart_Vars           ,ONLY: RestartTime
#endif
#if USE_FAXEN_CORR
USE MOD_Mesh_Vars              ,ONLY: nElems,nSides
#endif /* USE_FAXEN_CORR */
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)               :: ArraySize
#if USE_BASSETFORCE
INTEGER                       :: k
REAL                          :: s32
#endif /* USE_BASSETFORCE */
!===================================================================================================================================
! Allocate array to hold particle properties
! CALL Allocate_Safe(PartState    ,(/PP_nVarPart   ,PDM%maxParticleNumber/))
! CALL Allocate_Safe(PartReflCount,(/               PDM%maxParticleNumber/))
! CALL Allocate_Safe(LastPartPos  ,(/3             ,PDM%maxParticleNumber/))
! CALL Allocate_Safe(PartPosRef   ,(/3             ,PDM%MaxParticleNumber/))
! CALL Allocate_Safe(PartSpecies  ,(/               PDM%maxParticleNumber/))
! ! Allocate array for Runge-Kutta time stepping
! CALL Allocate_Safe(Pt           ,(/PP_nVarPartRHS,PDM%maxParticleNumber/))
! CALL Allocate_Safe(Pt_temp      ,(/PP_nVarPart   ,PDM%maxParticleNumber/))
! ! Allocate array for particle position in reference coordinates
! CALL Allocate_Safe(PDM%ParticleInside  ,(/        PDM%maxParticleNumber/))
! CALL Allocate_Safe(PDM%nextFreePosition,(/        PDM%maxParticleNumber/))
! CALL Allocate_Safe(PDM%IsNewPart       ,(/        PDM%maxParticleNumber/))
! ! Allocate particle-to-element-mapping (PEM) arrays
! CALL Allocate_Safe(PEM%Element         ,(/        PDM%maxParticleNumber/))
! CALL Allocate_Safe(PEM%lastElement     ,(/        PDM%maxParticleNumber/))

ArraySize = INT(((2.*REAL(PP_nVarPart) + REAL(PP_nVarPartRHS) + 6.)*SIZE_REAL + &
                  5.*SIZE_INT          + 2.*SIZE_LOG)              *PDM%maxParticleNumber,KIND=8)
IF (.NOT.VerifyMemUsage(ArraySize)) &
  CALL Abort(__STAMP__,'ERROR in particle_init.f90: Cannot allocate particle arrays. Array size too large!')

ALLOCATE(PartState(     1:PP_nVarPart,1:PDM%maxParticleNumber), &
         PartReflCount(               1:PDM%maxParticleNumber), &
         LastPartPos(             1:3,1:PDM%maxParticleNumber), &
         PartPosRef(              1:3,1:PDM%MaxParticleNumber), &
         PartSpecies(                 1:PDM%maxParticleNumber), &
! Allocate array for Runge-Kutta time stepping
         Pt(       1:PP_nVarPartRHS  ,1:PDM%maxParticleNumber), &
         Pt_temp(  1:PP_nVarPartRHS+3,1:PDM%maxParticleNumber), &
! Allocate array for particle position in reference coordinates
         PDM%ParticleInside(          1:PDM%maxParticleNumber), &
         PDM%nextFreePosition(        1:PDM%maxParticleNumber), &
         PDM%IsNewPart(               1:PDM%maxParticleNumber), &
! Allocate particle-to-element-mapping (PEM) arrays
         PEM%Element(                 1:PDM%maxParticleNumber), &
         PEM%lastElement(             1:PDM%maxParticleNumber)  &
#if PARTICLES_COUPLING >= 2
        ,PartNodeSource(1:PP_nVar,1:8,1:PDM%maxParticleNumber)  &
#endif /*PARTICLES_COUPLING >= 2*/
         )

PDM%ParticleInside(1:PDM%maxParticleNumber)  = .FALSE.
PDM%IsNewPart(     1:PDM%maxParticleNumber)  = .FALSE.
LastPartPos(   1:3,1:PDM%maxParticleNumber)  = 0.
PartState                                    = 0.
PartReflCount                                = 0.
Pt                                           = 0.
PartSpecies                                  = 0
PDM%nextFreePosition(1:PDM%maxParticleNumber)= 0
PEM%Element(         1:PDM%maxParticleNumber)= 0
PEM%lastElement(     1:PDM%maxParticleNumber)= 0
Pt_temp                                      = 0
#if PARTICLES_COUPLING >= 2
PartNodeSource                               = 0
#endif /*PARTICLES_COUPLING >= 2*/
PartPosRef                                   =-888.

IF(doPartIndex) THEN
  IF (.NOT.VerifyMemUsage(INT(PDM%maxParticleNumber*SIZE_INT,KIND=8))) &
    CALL Abort(__STAMP__,'ERROR in particle_init.f90: Cannot allocate particle arrays. Array size too large!')

  ALLOCATE(PartIndex(1:PDM%maxParticleNumber))
  PartIndex                                  = 0
END IF

! Extended RHS
#if USE_FAXEN_CORR
ArraySize = INT((9*3*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*nElems)*SIZE_REAL*PDM%maxParticleNumber,KIND=8)
IF (.NOT.VerifyMemUsage(ArraySize)) &
  CALL Abort(__STAMP__,'ERROR in particle_init.f90: Cannot allocate particle arrays. Array size too large!')

ALLOCATE(gradUx2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems),  &
         gradUy2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems),  &
         gradUz2(1:3,1:3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems),  &
!         U_local   (PRIM,0:PP_N,0:PP_N,0:PP_NZ,1:nElems),  &
!         gradp_local(1,3,0:PP_N,0:PP_N,0:PP_NZ,1:nElems),  &
         gradUx_master_loc( 1:3,0:PP_N,0:PP_NZ,1:nSides),  &
         gradUx_slave_loc(  1:3,0:PP_N,0:PP_NZ,1:nSides),  &
         gradUy_master_loc( 1:3,0:PP_N,0:PP_NZ,1:nSides),  &
         gradUy_slave_loc(  1:3,0:PP_N,0:PP_NZ,1:nSides),  &
         gradUz_master_loc( 1:3,0:PP_N,0:PP_NZ,1:nSides),  &
         gradUz_slave_loc(  1:3,0:PP_N,0:PP_NZ,1:nSides))
#endif /* USE_FAXEN_CORR */

! Basset force
#if USE_BASSETFORCE
N_Basset = GETINT('Part-NBasset')
s32 = 3./2.
nBassetVars = INT((N_Basset+1) * 3)
ALLOCATE(durdt(1:nBassetVars,1:PDM%maxParticleNumber))
durdt = 0.
ALLOCATE(bIter(1:PDM%maxParticleNumber))
bIter = 0
ALLOCATE(Fbdt(1:N_Basset+2,1:PDM%maxParticleNumber))
Fbdt = 0.
Fbdt(1,:) = RestartTime
! Window kernel
ALLOCATE(FbCoeff(1:2*N_Basset-1))
DO k=1,N_Basset-1
  FbCoeff(k) = ((k+s43)/((k+1)*SQRT(REAL(k+1))+(k+s32)*SQRT(REAL(k)))+(k-s43)/((k-1)*SQRT(REAL(k-1))+(k-s32)*SQRT(REAL(k))))
  FbCoeff(N_Basset+k-1) = (k-s43)/((k-1)*SQRT(REAL(k-1))+(k-s32)*SQRT(REAL(k)))
END DO
FbCoeff(N_Basset) = (N_Basset-s43)/((N_Basset-1)*SQRT(REAL(N_Basset-1))+(N_Basset-s32)*SQRT(REAL(N_Basset)))
! Exponential kernel
!FbCoeffm = 10
!ALLOCATE(Fbi(1:3,1:FbCoeffm,1:PDM%maxParticleNumber))
!Fbi = 0.
!ALLOCATE(FbCoeffa(FbCoeffm),FbCoefft(FbCoeffm))
!FbCoeffa = (/0.23477481312586,0.28549576238194,0.28479416718255,0.26149775537574,0.32056200511938,0.35354490689146,&
!  0.39635904496921,0.42253908596514,0.48317384225265,0.63661146557001/)
!FbCoefft = (/0.1,0.3,1.,3.,10.,40.,190.,1000.,6500.,50000./)
#endif /* USE_BASSETFORCE */

#if USE_EXTEND_RHS && ANALYZE_RHS
ALLOCATE(Pt_ext(18,PDM%maxParticleNumber))
Pt_ext = 0.
#endif /* USE_EXTEND_RHS && ANALYZE_RHS */

END SUBROUTINE AllocateParticleArrays


SUBROUTINE InitializeVariablesRandomNumbers()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSeed, nRandomSeeds,SeedSize
CHARACTER(32)         :: hilf
!===================================================================================================================================
!--- initialize randomization
! Read print flags
nRandomSeeds = GETINT('Part-NumberOfRandomSeeds')
! Specifies compiler specific minimum number of seeds
CALL RANDOM_SEED(Size = SeedSize)

ALLOCATE(Seeds(SeedSize))
! to ensure a solid run when an unfitting number of seeds is provided in ini
Seeds(:) = 1
SELECT CASE(nRandomSeeds)
  CASE(-1)
    ! Ensures different random numbers through irreproducible random seeds (via System_clock)
    CALL InitRandomSeed(nRandomSeeds,SeedSize,Seeds)

  CASE(0)
    ! IF (Restart) THEN
    !   CALL !numbers from state file
    ! ELSE IF (.NOT.Restart) THEN
    CALL InitRandomSeed(nRandomSeeds,SeedSize,Seeds)
    ! END IF

  CASE(1 :)
    ! read in numbers from ini
    IF (nRandomSeeds.NE.SeedSize) THEN
      SWRITE (UNIT_stdOut,'(A,I0,A,I0,A)') ' | Expected ',SeedSize,' seeds. Provided ', nRandomSeeds ,&
                                           '. Computer uses default value for all unset values.'
    END IF

    DO iSeed = 1,MIN(SeedSize,nRandomSeeds)
      WRITE(UNIT=hilf,FMT='(I0)') iSeed
      Seeds(iSeed) = GETINT('Part-RandomSeed'//TRIM(hilf))
    END DO

    IF (ALL(Seeds(:).EQ.0)) CALL CollectiveStop(__STAMP__,'Not all seeds can be set to zero ')
    CALL InitRandomSeed(nRandomSeeds,SeedSize,Seeds)

  CASE DEFAULT
    SWRITE (UNIT_stdOut,'(A)') ' Error: nRandomSeeds not defined.'        //&
                               ' Choose nRandomSeeds'                     //&
                               ' =-1    pseudo random'                    //&
                               ' = 0    hard-coded deterministic numbers' //&
                               ' > 0    numbers from ini. Expected ',SeedSize,'seeds.'

END SELECT

END SUBROUTINE InitializeVariablesRandomNumbers


SUBROUTINE InitializeVariablesSpeciesInits()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_ISO_VARYING_STRING
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: PI
USE MOD_Equation_Vars         ,ONLY: RefStatePrim
USE MOD_Particle_Globals
USE MOD_Particle_RHS          ,ONLY: InitRHS
USE MOD_Particle_Vars
USE MOD_ReadInTools
USE MOD_Utils                 ,ONLY: ALMOSTEQUAL
USE MOD_Viscosity
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec,iInit,iExclude
CHARACTER(32)         :: tmpStr,tmpStr2,tmpStr3
CHARACTER(255)        :: Filename_loc             ! specifying keyword for velocity distribution
INTEGER               :: PartField_shape(2)
INTEGER               :: drag_factor,RefStatePart
REAL                  :: tmp,charactLength,mu,prim(1:PP_nVarPrim)
!===================================================================================================================================
! Loop over all species and get requested data
RefStatePart            = GETINT('Part-IniRefState'   )
! characteristic length scale
charactLength           = GETREAL('Part-charactLength')
DO iSpec = 1, nSpecies
  LBWRITE(UNIT_stdOut,'(66("-"))')
  WRITE(UNIT=tmpStr,FMT='(I2)') iSpec

  ! Get number of requested inits
  Species(iSpec)%NumberOfInits         = GETINT('Part-Species'//TRIM(ADJUSTL(tmpStr))//'-nInits')
  ALLOCATE(Species(iSpec)%Init(0:Species(iSpec)%NumberOfInits))

  ! get species values // only once
  !--> General Species Values
  LBWRITE(UNIT_StdOut,'(A,I0,A,I0)') ' | Reading general  particle properties for Species',iSpec
  Species(iSpec)%RHSMethod             = GETINTFROMSTR('Part-Species'//TRIM(ADJUSTL(tmpStr))//'-RHSMethod'      )
  SELECT CASE (Species(iSpec)%RHSMethod)
    CASE (RHS_INERTIA)
#if PARABOLIC
      drag_factor                      = GETINTFROMSTR('Part-Species'//TRIM(ADJUSTL(tmpStr))//'-DragFactor'     )
#else
      drag_factor                      = DF_PART_STOKES
#endif
      CALL InitRHS(drag_factor, Species(iSpec)%DragFactor_pointer)
    CASE (RHS_INERTIA_EULER)
      drag_factor                      = DF_PART_STOKES
      CALL InitRHS(drag_factor, Species(iSpec)%DragFactor_pointer)
    CASE DEFAULT
      ! do nothing
  END SELECT
  Species(iSpec)%DensityIC             = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-DensityIC'      )
  Species(iSpec)%StokesIC              = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-StokesIC'       )
#if USE_PARTTEMP
  Species(iSpec)%SpecificHeatIC        = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-SpecificHeatIC' )
  IF (Species(iSpec)%SpecificHeatIC .LE. 0.) CALL CollectiveStop(__STAMP__, 'Invalid particle specific heat given, Species=',IntInfo=iSpec)
#endif /*USE_PARTTEMP*/

  IF (Species(iSpec)%RHSMethod .EQ. RHS_TCONVERGENCE) THEN
    IF (Species(iSpec)%StokesIC .EQ. 0) CALL COLLECTIVESTOP(__STAMP__,'Stokes number is zero!')
  ELSEIF (Species(iSpec)%StokesIC .GT. 0.) THEN
    IF (Species(iSpec)%RHSMethod .EQ. RHS_INERTIA_EULER) &
      CALL COLLECTIVESTOP(__STAMP__,'Stokes number greater than zero are not allowed with RHS_INERTIA_EULER!')
#if PARABOLIC
    ! dyn. viscosity
    prim = RefStatePrim(:,RefStatePart)
    mu   = VISCOSITY_PRIM(prim)
    tmp  = 18 * mu / (Species(iSpec)%DensityIC) * charactLength
    !tmp = 18 * mu * charactLength / (Species(iSpec)%DensityIC * NORM2(RefStatePrim(VELV,RefStatePart)))
    Species(iSpec)%DiameterIC = SQRT(Species(iSpec)%StokesIC * tmp)
    LBWRITE(UNIT_stdOut,'(A,I0,A,E16.5)') ' | Diameter of species (spherical) ', iSpec, ' = ', Species(iSpec)%DiameterIC
    Species(iSpec)%MassIC = MASS_SPHERE(Species(iSpec)%DensityIC, Species(iSpec)%DiameterIC)
    LBWRITE(UNIT_stdOut,'(A,I0,A,E16.5)') ' | Mass of species (spherical)     ', iSpec, ' = ', Species(iSpec)%MassIC
#endif
  ELSE
    Species(iSpec)%MassIC                = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-MassIC'         )
    Species(iSpec)%DiameterIC            = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-DiameterIC'     )
    Species(iSpec)%ScalePartDiam         = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-ScalePartDiam'  )
    IF (Species(iSpec)%DensityIC .EQ. 0.) THEN
      IF (Species(iSpec)%DiameterIC .EQ. 0.) CALL COLLECTIVESTOP(__STAMP__,'Particle density and diameter both zero!')
      IF (Species(iSpec)%MassIC     .EQ. 0.) CALL COLLECTIVESTOP(__STAMP__,'Particle density and mass     both zero!')
      Species(iSpec)%DensityIC = DENS_SPHERE(Species(iSpec)%MassIC, Species(iSpec)%DiameterIC)
    ELSEIF (Species(iSpec)%MassIC .EQ. 0.) THEN
      IF (Species(iSpec)%DiameterIC .EQ. 0.) CALL COLLECTIVESTOP(__STAMP__,'Particle mass and diameter both zero!')
      IF (Species(iSpec)%DensityIC  .EQ. 0.) CALL COLLECTIVESTOP(__STAMP__,'Particle mass and density  both zero!')
      Species(iSpec)%MassIC = MASS_SPHERE(Species(iSpec)%DensityIC, Species(iSpec)%DiameterIC)
      LBWRITE(UNIT_stdOut,'(A,I0,A,E16.5)') ' | Mass of species (spherical)     ', iSpec, ' = ', Species(iSpec)%MassIC
    ELSEIF (Species(iSpec)%DiameterIC .EQ. 0.) THEN
      IF (Species(iSpec)%DensityIC .EQ. 0.) CALL COLLECTIVESTOP(__STAMP__,'Particle density and diameter both zero!')
      Species(iSpec)%DiameterIC = DIAM_SPHERE(Species(iSpec)%DensityIC, Species(iSpec)%MassIC)
      LBWRITE(UNIT_stdOut,'(A,I0,A,E16.5)') ' | Diameter of species (spherical) ', iSpec, ' = ', Species(iSpec)%DiameterIC
    END IF
  END IF

  Species(iSpec)%LowVeloThreshold      = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-LowVeloThreshold' )
  Species(iSpec)%SphericityIC          = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-SphericityIC'     )
  ! Warn when outside valid range of Haider model
  IF ((drag_factor .EQ. DF_PART_HAIDER) .AND. (Species(iSpec)%SphericityIC .LT. 0.670)) THEN
    LBWRITE(UNIT_stdOut,*) 'WARNING: SphericityIC ', Species(iSpec)%SphericityIC,'< 0.670, drag coefficient may not be accurate.'
  END IF
  IF ((drag_factor .EQ. DF_PART_HOELZER) .AND. (Species(iSpec)%SphericityIC .EQ. 1.0)) THEN
    LBWRITE(UNIT_stdOut,*) 'INFO: SphericityIC ', Species(iSpec)%SphericityIC,'== 1, drag coefficient of Schiller and Naumann is used.'
    drag_factor = DF_PART_SCHILLER
    CALL InitRHS(drag_factor, Species(iSpec)%DragFactor_pointer)
  END IF
  IF ((drag_factor .NE. DF_PART_HAIDER .AND. drag_factor .NE. DF_PART_HOELZER) .AND. (Species(iSpec)%SphericityIC .LT. 1.0)) THEN
    LBWRITE(UNIT_stdOut,*) 'WARNING: SphericityIC is ignored for ', iSpec,', as the wrong drag force was chosen.'
  END IF
  Species(iSpec)%PartDiamVarianceIC    = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-PartDiamVarianceIC' )
#if USE_SPHERICITY
  Species(iSpec)%PartSpheVarianceIC    = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-PartSpheVarianceIC' )
#endif
#if USE_EXTEND_RHS
  Species(iSpec)%CalcSaffmanForce      = GETLOGICAL(   'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-CalcSaffmanForce'   )
  Species(iSpec)%CalcMagnusForce       = GETLOGICAL(   'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-CalcMagnusForce'    )
  Species(iSpec)%CalcVirtualMass       = GETLOGICAL(   'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-CalcVirtualMass'    )
  Species(iSpec)%CalcUndisturbedFlow   = GETLOGICAL(   'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-CalcUndisturbedFlow')
  Species(iSpec)%CalcBassetForce       = GETLOGICAL(   'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-CalcBassetForce    ')
#endif

  !--> Bons particle rebound model
  LBWRITE(UNIT_stdOut,'(A,I0,A,I0)') ' | Reading rebound  particle properties for Species',iSpec
  Species(iSpec)%YoungIC               = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-YoungIC'            )
  Species(iSpec)%PoissonIC             = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-PoissonIC'          )
  Species(iSpec)%YieldCoeff            = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-YieldCoeff'         )
  Species(iSpec)%Whitaker_a            = GETREAL(      'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-Whitaker_a'         )

  !--> Rough wall modelling
  Species(iSpec)%doRoughWallModelling  = GETLOGICAL(   'Part-Species'//TRIM(ADJUSTL(tmpStr))//'-RoughWall'          )

  !--> Check if particles have valid mass/density
  IF (Species(iSpec)%MassIC .LE. 0.    .AND..NOT.(Species(iSpec)%RHSMethod.EQ.RHS_NONE          &
                                       .OR.       Species(iSpec)%RHSMethod.EQ.RHS_TRACER        &
                                       .OR.       Species(iSpec)%RHSMethod.EQ.RHS_HPCONVERGENCE &
                                       .OR.       Species(iSpec)%RHSMethod.EQ.RHS_TCONVERGENCE)) &
    CALL CollectiveStop(__STAMP__, 'Invalid particle mass given, Species=',IntInfo=iSpec)

  IF (Species(iSpec)%DensityIC .LE. 0. .AND..NOT.(Species(iSpec)%RHSMethod.EQ.RHS_NONE          &
                                       .OR.       Species(iSpec)%RHSMethod.EQ.RHS_TRACER        &
                                       .OR.       Species(iSpec)%RHSMethod.EQ.RHS_HPCONVERGENCE &
                                       .OR.       Species(iSpec)%RHSMethod.EQ.RHS_TCONVERGENCE)) &
    CALL CollectiveStop(__STAMP__, 'Invalid particle density given, Species=',IntInfo=iSpec)

  ! Loop over all inits and get requested data
  DO iInit = 0, Species(iSpec)%NumberOfInits
    ! set help characters to read strings from parameter file
    IF(iInit.EQ.0)THEN
      tmpStr2=TRIM(ADJUSTL(tmpStr))
    ELSE ! iInit >0
      WRITE(UNIT=tmpStr2,FMT='(I0)') iInit
      tmpStr2=TRIM(ADJUSTL(tmpStr))//'-Init'//TRIM(ADJUSTL(tmpStr2))
    END IF ! iInit

    ! Emission and init data
    LBWRITE(UNIT_stdOut,'(A,I0,A,I0)') ' | Reading emission particle properties for Species',iSpec,'-Init',iInit
    Species(iSpec)%Init(iInit)%UseForInit            = GETLOGICAL(  'Part-Species'//TRIM(tmpStr2)//'-UseForInit'            )
    Species(iSpec)%Init(iInit)%UseForEmission        = GETLOGICAL(  'Part-Species'//TRIM(tmpStr2)//'-UseForEmission'        )
    Species(iSpec)%Init(iInit)%initialParticleNumber = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-initialParticleNumber' )
    Species(iSpec)%Init(iInit)%ParticleEmissionType  = GETINT(      'Part-Species'//TRIM(tmpStr2)//'-ParticleEmissionType'  )
    Species(iSpec)%Init(iInit)%ParticleEmission      = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-ParticleEmission'      )
    Species(iSpec)%Init(iInit)%ParticleEmissionTime  = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-ParticleEmissionTime'  )
    Species(iSpec)%Init(iInit)%PartDensity           = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-PartDensity'           )
    Species(iSpec)%Init(iInit)%SpaceIC               = TRIM(ADJUSTL(GETSTR( 'Part-Species'//TRIM(tmpStr2)//'-SpaceIC'      )))
    Species(iSpec)%Init(iInit)%VeloVecIC             = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-VeloVecIC'           ,3)
    Species(iSpec)%Init(iInit)%VeloIC                = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-VeloIC'                )
    Species(iSpec)%Init(iInit)%VeloTurbIC            = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-VeloTurbIC'            )
    Species(iSpec)%Init(iInit)%velocityDistribution  = TRIM(ADJUSTL(GETSTR( 'Part-Species'//TRIM(tmpStr2)//'-velocityDistribution')))
#if USE_PARTTEMP
    Species(iSpec)%Init(iInit)%TemperatureIC         = GETREAL(      'Part-Species'//TRIM(tmpStr2)//'-TemperatureIC' )
    IF (Species(iSpec)%Init(iInit)%TemperatureIC .LT. 0.) CALL CollectiveStop(__STAMP__, 'Invalid particle initial temperature given, Species=',IntInfo=iSpec)
    Species(iSpec)%Init(iInit)%tempDistribution  = TRIM(ADJUSTL(GETSTR( 'Part-Species'//TRIM(tmpStr2)//'-tempDistribution')))
#endif /*USE_PARTTEMP*/
    IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution) .EQ. "load_from_file") THEN
      ! load from binary file
      IF (MPIRoot) THEN
        WRITE(FileName_loc,"(A200,I2,A4)") FilenameRecordPart, iSpec-1, ".dat"
        Filename_loc = TRIM(REPLACE(Filename_loc," ","",Every=.TRUE.))
        OPEN(33, FILE=TRIM(FileName_loc),FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
        READ(33) PartField_shape
        ALLOCATE(Species(iSpec)%Init(iInit)%PartField(PartField_shape(1), PartField_shape(2)))
        READ(33) Species(iSpec)%Init(iInit)%PartField(:,:)
        CLOSE(33)
      END IF
#if USE_MPI
      CALL MPI_BCAST(PartField_shape,2,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
      Species(iSpec)%Init(iInit)%nPartField = PartField_shape(1)
      IF(.NOT.ALLOCATED(Species(iSpec)%Init(iInit)%PartField))&
        ALLOCATE(Species(iSpec)%Init(iInit)%PartField(PartField_shape(1), PartField_shape(2)))
      CALL MPI_BCAST(Species(iSpec)%Init(iInit)%PartField,PartField_shape(1)*PartField_shape(2),MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
#endif
    END IF
    Species(iSpec)%Init(iInit)%InflowRiseTime        = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-InflowRiseTime'        )
    Species(iSpec)%Init(iInit)%NumberOfExcludeRegions= GETINT(      'Part-Species'//TRIM(tmpStr2)//'-NumberOfExcludeRegions')

    !> Read only emission properties required for SpaceIC
    !>>> Parameters must be set to false to allow conformity checks afterwards
    Species(iSpec)%Init(iInit)%CalcHeightFromDt = .FALSE.

    SELECT CASE(Species(iSpec)%Init(iInit)%SpaceIC)
      CASE('point')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
      CASE('line','line_with_equidistant_distribution')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%BaseVector1IC     = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector1IC' ,3)
      CASE('plane')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%BaseVector1IC     = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector1IC' ,3)
        Species(iSpec)%Init(iInit)%BaseVector2IC     = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector2IC' ,3)
      CASE('disc')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%RadiusIC          = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-RadiusIC'      )
        Species(iSpec)%Init(iInit)%NormalIC          = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-NormalIC'      ,3)
        IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution) .EQ. "load_from_file") THEN
          Species(iSpec)%Init(iInit)%RadiusIC  = &
            MIN(MAXVAL(SQRT(Species(iSpec)%Init(iInit)%PartField(:,1)**2+Species(iSpec)%Init(iInit)%PartField(:,2)**2)),Species(iSpec)%Init(iInit)%RadiusIC)
          LBWRITE(UNIT_stdOut,'(A,I0,A,E16.5)') ' | RadiusIC of species ', iSpec, ' = ', Species(iSpec)%Init(iInit)%RadiusIC
        END IF
      CASE('circle', 'circle_equidistant')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%RadiusIC          = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-RadiusIC'        )
        Species(iSpec)%Init(iInit)%NormalIC          = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-NormalIC'      ,3)
      CASE('cuboid')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%BaseVector1IC     = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector1IC' ,3)
        Species(iSpec)%Init(iInit)%BaseVector2IC     = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector2IC' ,3)
        Species(iSpec)%Init(iInit)%CuboidHeightIC    = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-CuboidHeightIC'  )
        Species(iSpec)%Init(iInit)%CalcHeightFromDt  = GETLOGICAL(  'Part-Species'//TRIM(tmpStr2)//'-CalcHeightFromDt')
      CASE('cylinder')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%RadiusIC          = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-RadiusIC'        )
        Species(iSpec)%Init(iInit)%BaseVector1IC     = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector1IC' ,3)
        Species(iSpec)%Init(iInit)%BaseVector2IC     = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BaseVector2IC' ,3)
        Species(iSpec)%Init(iInit)%CylinderHeightIC  = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-CylinderHeightIC')
      CASE('sphere')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%RadiusIC          = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-RadiusIC'        )
      CASE('Gaussian')
        Species(iSpec)%Init(iInit)%BasePointIC       = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-BasePointIC'   ,3)
        Species(iSpec)%Init(iInit)%BaseVariance      = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-BaseVariance'    )
        Species(iSpec)%Init(iInit)%RadiusIC          = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-RadiusIC'        )
        Species(iSpec)%Init(iInit)%NormalIC          = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-NormalIC'      ,3)
!      CASE('sin_deviation')
        ! Currently not implemented
      CASE('cell_local')
        ! No readin
      CASE DEFAULT
        CALL CollectiveStop(__STAMP__,'Unknown particle emission type')
    END SELECT

    ! Nullify additional init data here
    Species(iSpec)%Init(iInit)%InsertedParticle        = 0
    Species(iSpec)%Init(iInit)%InsertedParticleSurplus = 0

    ! Get absolute value of particle velocity vector and normalize the VeloVecIC vector
    IF (Species(iSpec)%Init(iInit)%VeloIC.EQ.0.) THEN
      Species(iSpec)%Init(iInit)%VeloIC                = VECNORM(Species(iSpec)%Init(iInit)%VeloVecIC(1:3))
    END IF

    ! Only normalize if the vector does not have zero length. If it has, our job is done
    IF (VECNORM(Species(iSpec)%Init(iInit)%VeloVecIC(1:3)).NE.0) THEN
      Species(iSpec)%Init(iInit)%VeloVecIC             = Species(iSpec)%Init(iInit)%VeloVecIC/VECNORM(Species(iSpec)%Init(iInit)%VeloVecIC(1:3))
    END IF

    !----------- various checks/calculations after read-in of Species(i)%Init(iInit)%-data ----------------------------------!

    !--- Check if Initial ParticleInserting is really used
    IF (Species(iSpec)%Init(iInit)%UseForInit                   .AND. &
       (Species(iSpec)%Init(iInit)%PartDensity          .EQ.0)  .AND. &
       (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0)) THEN
      Species(iSpec)%Init(iInit)%UseForInit = .FALSE.
      LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A)',ADVANCE='NO') ' | WARNING: Species',iSpec,'-Init',iInit,' - no ParticleNumber or PartDensity detected,'
      LBWRITE(UNIT_stdOut,'(A)')                        ' disabling initial particle inserting!'
    END IF

    !--- Check if ParticleEmission is really used
    IF (Species(iSpec)%Init(iInit)%UseForEmission         .AND. &
       (Species(iSpec)%Init(iInit)%ParticleEmission.EQ.0)) THEN
      Species(iSpec)%Init(iInit)%UseForEmission = .FALSE.
      LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A)',ADVANCE='NO') ' | WARNING: Species',iSpec,'-Init',iInit,' - no emission rate                 detected,'
      LBWRITE(UNIT_stdOut,'(A)')                        ' disabling particle emission!'
    END IF

    !--- Check if cell_local is used with correct emission type
    IF (Species(iSpec)%Init(iInit)%UseForInit                   .AND. &
       (Species(iSpec)%Init(iInit)%PartDensity          .EQ.0)  .AND. &
       (Species(iSpec)%Init(iInit)%initialParticleNumber.EQ.0)) THEN
      ! Species(iSpec)%Init(iInit)%UseForInit = .FALSE.
      LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A)',ADVANCE='NO') ' | WARNING: Species',iSpec,'-Init',iInit,' - no ParticleNumber or PartDensity detected,'
      LBWRITE(UNIT_stdOut,'(A)')                        ' disabling initial particle inserting!'
    END IF

    !--- cuboid-/cylinder-height calculation from v and dt
    SELECT CASE(Species(iSpec)%Init(iInit)%SpaceIC)
      CASE('cuboid')
        IF (.NOT.Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
          IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CuboidHeightIC,-1.)) THEN
            Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.
            LBWRITE(UNIT_stdOut,'(A)') " | WARNING: Cuboid height will be calculated from v and dt!"
          END IF
        END IF

      CASE('cylinder')
        IF (.NOT.Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
          IF (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%CylinderHeightIC,-1.)) THEN
            Species(iSpec)%Init(iInit)%CalcHeightFromDt=.TRUE.
            LBWRITE(UNIT_stdOut,'(A)') " | WARNING: Cylinder height will be calculated from v and dt!"
          END IF
        END IF

      CASE('cell_local')
        IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.NE.0) THEN
          Species(iSpec)%Init(iInit)%UseForInit = .FALSE.
          LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A)',ADVANCE='NO') ' | WARNING: Species',iSpec,'-Init',iInit,' - no cell_local & EmissionType !=0 detected,'
          LBWRITE(UNIT_stdOut,'(A)')                        ' disabling initial particle inserting!'
        END IF

      CASE DEFAULT
        IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
          CALL CollectiveStop(__STAMP__,' Calculating height from v and dt is only supported for cuboid or cylinder!')
        END IF
    END SELECT

    IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN
      IF (Species(iSpec)%Init(iInit)%UseForInit) THEN
        CALL CollectiveStop(__STAMP__,' Calculating height from v and dt is not supported for initial ParticleInserting!')
      END IF
    END IF

    !--- read ExcludeRegions and normalize/calculate corresponding vectors
    IF (Species(iSpec)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      ALLOCATE(Species(iSpec)%Init(iInit)%ExcludeRegion(1:Species(iSpec)%Init(iInit)%NumberOfExcludeRegions))

      IF  ((TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cuboid') &
       .OR.(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN

        ! Read in information for exclude regions
        DO iExclude=1,Species(iSpec)%Init(iInit)%NumberOfExcludeRegions
          WRITE(UNIT=tmpStr3,FMT='(I0)') iExclude
          tmpStr3=TRIM(tmpStr2)//'-ExcludeRegion'//TRIM(tmpStr3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC              &
            = TRIM(GETSTR('Part-Species'//TRIM(tmpStr3)//'-SpaceIC','cuboid'))
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC             &
            = GETREAL(    'Part-Species'//TRIM(tmpStr3)//'-RadiusIC')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC            &
            = GETREAL(    'Part-Species'//TRIM(tmpStr3)//'-Radius2IC')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC             &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-NormalIC',3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC          &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-BasePointIC',3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC        &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-BaseVector1IC',3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC        &
            = GETREALARRAY('Part-Species'//TRIM(tmpStr3)//'-BaseVector2IC',3)
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC       &
            = GETREAL(     'Part-Species'//TRIM(tmpStr3)//'-CuboidHeightIC')
          Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC     &
            = GETREAL(     'Part-Species'//TRIM(tmpStr3)//'-CylinderHeightIC')

          !--check and normalize data
          IF ((TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).EQ.'cuboid') .OR.             &
               ((((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1),1.)   &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2),0.))  &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3),0.))  &
            .OR. ((.NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1),0.)   &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2),1.))  &
              .OR. .NOT.ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3),0.))) &
            .AND.    (((ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC     (1),0.))  &
              .AND.    (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC     (2),0.))) &
              .AND.    (ALMOSTEQUAL(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC     (3),1.))))) THEN
            !-- cuboid; or BV are non-default and NormalIC is default: calc. NormalIC for ExcludeRegions from BV1/2
            !   (for def. BV and non-def. NormalIC; or all def. or non-def.: Use User-defined NormalIC when ExclRegion is cylinder)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1)                  &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2)         &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)         &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)         &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2)                  &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)         &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)         &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1)         &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
            Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)                  &
              = Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1)         &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2)         &
              - Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2)         &
              * Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1)
          ELSE IF ( (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cuboid')    .AND. &
                    (TRIM(Species(iSpec)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC).NE.'cylinder')) THEN
            CALL CollectiveStop(__STAMP__,'Error in ParticleInit, ExcludeRegions must be cuboid or cylinder!')
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
            CALL CollectiveStop(__STAMP__,'Error in ParticleInit, NormalIC Vector must not be zero!')
          END IF
        END DO !iExclude

      ! invalid combination of SpaceIC and exclude region
      ELSE
        CALL CollectiveStop(__STAMP__,'Error in ParticleInit, ExcludeRegions are currently only implemented for the SpaceIC cuboid or cylinder!')
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
      LBWRITE(UNIT_stdOut,'(A,I0,A,I0)') ' | StartNumberOfInits of Species ', iSpec, ' = ', Species(iSpec)%StartnumberOfInits
    END IF ! iInit .EQ.0

  END DO ! iInit
END DO ! iSpec

END SUBROUTINE InitializeVariablesSpeciesInits


SUBROUTINE InitializeVariablesPartBoundary()
!===================================================================================================================================
! Read in boundary parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: BoundaryName,BoundaryType,nBCs
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_Particle_Mesh_Vars     ,ONLY: MeshHasPeriodic
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF
USE MOD_Particle_Vars
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iBC
CHARACTER(32)         :: tmpStr
CHARACTER(200)        :: tmpString
!===================================================================================================================================
! get number of particle boundaries
nPartBound = CountOption('Part-BoundaryName')

! if no part boundaries are given, assume the same number as for the DG part
IF (nPartBound.EQ.0) THEN
  nPartBound = nBCs
END IF

LBWRITE(UNIT_stdOut,'(132("."))')
LBWRITE(UNIT_stdOut,'(A)') ' | Reading particle boundary properties'

! Allocate arrays for particle boundaries
ALLOCATE(PartBound%SourceBoundName     (1:nBCs))
ALLOCATE(PartBound%SourceBoundType     (1:nBCs))
ALLOCATE(PartBound%TargetBoundCond     (1:nBCs))
ALLOCATE(PartBound%WallTemp            (1:nBCs))
ALLOCATE(PartBound%WallVelo        (1:3,1:nBCs))
ALLOCATE(PartBound%WallModel           (1:nBCs))
ALLOCATE(PartBound%WallCoeffModel      (1:nBCs))
ALLOCATE(PartBound%doRoughWallModelling(1:nBCs))
ALLOCATE(PartBound%RoughMeanIC         (1:nBCs))
ALLOCATE(PartBound%RoughVarianceIC     (1:nBCs))
!ALLOCATE(PartBound%AmbientCondition   (1:nBCs))
!ALLOCATE(PartBound%AmbientConditionFix(1:nBCs))
!ALLOCATE(PartBound%AmbientTemp        (1:nBCs))
!ALLOCATE(PartBound%AmbientVelo    (1:3,1:nBCs))
!ALLOCATE(PartBound%AmbientDens        (1:nBCs))
!ALLOCATE(PartBound%AmbientDynamicVisc (1:nBCs))
PartBound%SourceBoundName = ''
PartBound%SourceBoundType = ''
PartBound%TargetBoundCond = -1
PartBound%WallTemp        = 0.
PartBound%WallVelo        = 0.
PartBound%WallModel       = ''
PartBound%WallCoeffModel  = ''
PartBound%doRoughWallModelling= .FALSE.
PartBound%RoughMeanIC     = 0.
PartBound%RoughVarianceIC = 0.

! Bons particle rebound model
ALLOCATE(PartBound%Young               (1:nBCs))
ALLOCATE(PartBound%Poisson             (1:nBCs))
ALLOCATE(PartBound%FricCoeff           (1:nBCs))

! Fong coefficent of restitution
ALLOCATE(PartBound%CoR                 (1:nBCs))

! Periodic BCs
MeshHasPeriodic = .FALSE.

! Surface Flux
ALLOCATE(BCdata_auxSF                  (1:nBCs))
DO iBC=1,nBCs
  ! init value when not used
  BCdata_auxSF(iBC)%SideNumber = -1
  BCdata_auxSF(iBC)%GlobalArea = 0.
  BCdata_auxSF(iBC)%LocalArea  = 0.
END DO

! Loop over all boundaries and get information
DO iBC = 1,nBCs
  IF (BoundaryType(iBC,1).EQ.0) THEN
    PartBound%TargetBoundCond(iBC) = PartBound%InternalBC
    LBWRITE(UNIT_stdOut,'(A,I0,A)') " ... PartBound ",iBC," is internal bound, no mapping needed"
    CYCLE
  END IF

  ! Set DG boundary as default input
  SELECT CASE (BoundaryType(iBC, BC_TYPE))
   CASE(1)
     tmpString = 'periodic'
  ! CASE(2,12,22)
  !   tmpString='open'
   CASE(3,4)
     tmpString = 'reflective'
   CASE(9)
     tmpString = 'symmetry'
   CASE DEFAULT
     tmpString = 'open'
  END SELECT

  ! Check if boundary is overwritten for particles
  WRITE(UNIT=tmpStr,FMT='(I0)') iBC
  PartBound%SourceBoundType(iBC) = TRIM(GETSTR('Part-Boundary'//TRIM(tmpStr)//'-Type', tmpString))
  PartBound%SourceBoundName(iBC) = TRIM(GETSTR('Part-Boundary'//TRIM(tmpStr)//'-Name', BoundaryName(iBC)))
  ! Select boundary condition for particles
  SELECT CASE (PartBound%SourceBoundType(iBC))
    ! Inflow / outflow
    CASE('open')
      PartBound%TargetBoundCond(iBC)      = PartBound%OpenBC
!      PartBound%AmbientCondition(iBC)     = GETLOGICAL(  'Part-Boundary'//TRIM(tmpStr)//'-AmbientCondition'   )
!      IF(PartBound%AmbientCondition(iBC)) THEN
!        PartBound%AmbientConditionFix(iBC)= GETLOGICAL(  'Part-Boundary'//TRIM(tmpStr)//'-AmbientConditionFix')
!        PartBound%AmbientVelo(1:3,iBC)    = GETREALARRAY('Part-Boundary'//TRIM(tmpStr)//'-AmbientVelo'      ,3)
!        PartBound%AmbientDens(iBC)        = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-AmbientDens'        )
!        PartBound%AmbientDynamicVisc(iBC) = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-AmbientDynamicVisc' ) ! N2:T=288K
!      END IF

    ! Reflective (wall)
    CASE('reflective')
      PartBound%TargetBoundCond(iBC)      = PartBound%ReflectiveBC
      PartBound%WallVelo(1:3,iBC)         = GETREALARRAY('Part-Boundary'//TRIM(tmpStr)//'-WallVelo'         ,3)
      PartBound%WallModel(iBC)            = GETSTR(      'Part-Boundary'//TRIM(tmpStr)//'-WallModel'          )

      ! Rough wall modelling
      PartBound%doRoughWallModelling(iBC) = GETLOGICAL(  'Part-Boundary'//TRIM(tmpStr)//'-RoughWall'          )
      IF (PartBound%doRoughWallModelling(iBC)) THEN
        PartBound%RoughMeanIC(iBC)        = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-RoughMeanIC'        )
        ! Variance in radii
        PartBound%RoughVarianceIC(iBC)    = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-RoughVarianceIC'    )
        IF (PartBound%RoughVarianceIC(iBC).GT.1.0) &
          CALL CollectiveSTOP(__STAMP__,'Rough variance is high, disable this message only if you are sure that you are correct!')
      END IF
      ! Non-perfect reflection
      IF (PartBound%WallModel(iBC).EQ.'coeffRes') THEN
          PartBound%WallCoeffModel(iBC)   = GETSTR(      'Part-Boundary'//TRIM(tmpStr)//'-WallCoeffModel'     )

          SELECT CASE(PartBound%WallCoeffModel(iBC))
            CASE ('Bons2017','Whitaker2018')
              PartBound%Young(iBC)        = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-Young')
              PartBound%Poisson(iBC)      = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-Poisson')
              PartBound%FricCoeff(iBC)    = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-FricCoeff'          )

            ! Fong coefficient of reflection
            CASE('Fong2019')
              PartBound%CoR(iBC)          = GETREAL(     'Part-Boundary'//TRIM(tmpStr)//'-CoR'                )

            CASE('Tabakoff1981','Grant1975')
              ! nothing to readin

            ! Rebound ANN valid for v_{air} \in [150,350] [m/s]
            CASE('RebANN')
              ! Readin weights
              CALL LoadWeightsANN(GETSTR(      'Part-Boundary'//TRIM(tmpStr)//'-ANNModel'))

            CASE('FracANN')
              ! Readin weights
              CALL LoadWeightsANN(GETSTR(      'Part-Boundary'//TRIM(tmpStr)//'-ANNModel'))
              doWritePartDiam = .TRUE.
              doRandomPartDiam = .TRUE.

            CASE DEFAULT
              CALL CollectiveSTOP(__STAMP__,'Unknown wall model given!')

          END SELECT
      END IF

    ! Periodic
    CASE('periodic')
      PartBound%TargetBoundCond(iBC)      = PartBound%PeriodicBC
      MeshHasPeriodic                     = .TRUE.

    CASE('symmetry')
      PartBound%TargetBoundCond(iBC)      = PartBound%SymmetryBC

    ! Invalid boundary option
    CASE DEFAULT
      SWRITE(UNIT_stdOut,'(A,A)') ' Boundary does not exists: ', TRIM(PartBound%SourceBoundType(iBC))
      CALL CollectiveStop(__STAMP__,'Particle Boundary Condition does not exist')
  END SELECT
END DO

END SUBROUTINE InitializeVariablesPartBoundary


!SUBROUTINE InitializeVariablesTimeStep(ManualTimeStep_opt)
SUBROUTINE InitializeVariablesTimeStep()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Particle_TimeDisc_Vars,  ONLY: Pa_rebuilt,Pa_rebuilt_coeff,Pv_rebuilt,v_rebuilt
USE MOD_TimeDisc_Vars,           ONLY: RKA,nRKStages
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN),OPTIONAL   :: ManualTimeStep_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iStage_loc
!===================================================================================================================================
! Rebuild Pt_tmp-coefficients assuming F=const. (value at wall) in previous stages
ALLOCATE(Pa_rebuilt_coeff(1:nRKStages) &
        ,Pa_rebuilt  (1:PP_nVarPartRHS,1:nRKStages) &
        ,Pv_rebuilt  (1:3,1:nRKStages) &
        ,v_rebuilt   (1:3,0:nRKStages-1))

DO iStage_loc=1,nRKStages
  IF (iStage_loc.EQ.1) THEN
    Pa_rebuilt_coeff(iStage_loc) = 1.
  ELSE
    Pa_rebuilt_coeff(iStage_loc) = 1. - RKA(iStage_loc)*Pa_rebuilt_coeff(iStage_loc-1)
  END IF
END DO

END SUBROUTINE InitializeVariablesTimeStep


!===================================================================================================================================
! finalize particle variables
!===================================================================================================================================
SUBROUTINE FinalizeParticles()
! MODULES
USE MOD_Globals
USE MOD_Particle_Analyze,           ONLY: FinalizeParticleAnalyze
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Boundary_Sampling, ONLY: FinalizeParticleBoundarySampling
USE MOD_Particle_Boundary_Tracking, ONLY: FinalizeParticleBoundaryTracking
USE MOD_Particle_Interpolation,     ONLY: FinalizeParticleInterpolation
USE MOD_Particle_Mesh,              ONLY: FinalizeParticleMesh
USE MOD_Particle_Surfaces,          ONLY: FinalizeParticleSurfaces
USE MOD_Particle_TimeDisc,          ONLY: Particle_FinalizeTimeDisk
USE MOD_Particle_TimeDisc_Vars,     ONLY: Pa_rebuilt,Pa_rebuilt_coeff,Pv_rebuilt,v_rebuilt
USE MOD_Particle_Vars
#if USE_MPI
USE MOD_Particle_MPI_Emission,      ONLY: FinalizeEmissionComm
USE MOD_Particle_MPI_Halo,          ONLY: FinalizePartExchangeProcs
#endif /*USE_MPI*/
#if PARABOLIC
#if USE_RW
USE MOD_Particle_RandomWalk,        ONLY: ParticleFinalizeRandomWalk
#endif /*USE_RW*/
USE MOD_Particle_SGS,               ONLY: ParticleFinalizeSGS
#endif /*PARABOLIC*/
#if PARTICLES_COUPLING >= 2
USE MOD_Particle_Deposition,        ONLY: FinalizeDeposition
#endif /*PARTICLES_COUPLING >= 2*/
#if PARTICLES_COUPLING == 4
USE MOD_Particle_Collision,         ONLY: FinalizeCollision
#endif /*PARTICLES_COUPLING == 4*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_MPI
! Must be finalized before Species is deallocated
CALL FinalizeEmissionComm
#endif

! particle properties
SDEALLOCATE(Species)
SDEALLOCATE(PartState)
SDEALLOCATE(PartReflCount)
SDEALLOCATE(LastPartPos)
SDEALLOCATE(PartPosRef)
SDEALLOCATE(PartSpecies)
SDEALLOCATE(PartIndex)

! Runge-Kutta time stepping
SDEALLOCATE(Pt)
SDEALLOCATE(Pt_temp)
SDEALLOCATE(Pa_rebuilt)
SDEALLOCATE(Pa_rebuilt_coeff)
SDEALLOCATE(Pv_rebuilt)
SDEALLOCATE(v_rebuilt)

! particle position in reference coordinates
SDEALLOCATE(PDM%ParticleInside)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%nextFreePosition)
SDEALLOCATE(PDM%IsNewPart)

! particle boundary information
SDEALLOCATE(PartBound%SourceBoundName)
SDEALLOCATE(PartBound%SourceBoundType)
SDEALLOCATE(PartBound%TargetBoundCond)
SDEALLOCATE(PartBound%WallTemp)
SDEALLOCATE(PartBound%WallVelo)
!SDEALLOCATE(PartBound%AmbientCondition)
!SDEALLOCATE(PartBound%AmbientConditionFix)
!SDEALLOCATE(PartBound%AmbientTemp)
!SDEALLOCATE(PartBound%AmbientVelo)
!SDEALLOCATE(PartBound%AmbientDens)
!SDEALLOCATE(PartBound%AmbientDynamicVisc)
SDEALLOCATE(PartBound%WallModel)
SDEALLOCATE(PartBound%WallCoeffModel)
SDEALLOCATE(PartBound%doRoughWallModelling)
SDEALLOCATE(PartBound%RoughMeanIC)
SDEALLOCATE(PartBound%RoughVarianceIC)
SDEALLOCATE(PartBound%Young)
SDEALLOCATE(PartBound%Poisson)
SDEALLOCATE(PartBound%FricCoeff)
SDEALLOCATE(PartBound%CoR)
SDEALLOCATE(PartBoundANN%nN)
SDEALLOCATE(PartBoundANN%w)
SDEALLOCATE(PartBoundANN%b)
SDEALLOCATE(PartBoundANN%beta)
SDEALLOCATE(PartBoundANN%output)
SDEALLOCATE(PartBoundANN%max_in)
SDEALLOCATE(PartBoundANN%max_out)

! particle-to-element-mapping (PEM) arrays
SDEALLOCATE(PEM%Element)
SDEALLOCATE(PEM%lastElement)
SDEALLOCATE(PEM%pStart)
SDEALLOCATE(PEM%pNumber)
SDEALLOCATE(PEM%pEnd)
SDEALLOCATE(PEM%pNext)

! Particle deposition
#if PARTICLES_COUPLING >= 2
SDEALLOCATE(PartNodeSource)
#endif /*PARTICLES_COUPLING >= 2*/

! Extended RHS
#if USE_FAXEN_CORR
SDEALLOCATE(gradUx2)
SDEALLOCATE(gradUy2)
SDEALLOCATE(gradUz2)
SDEALLOCATE(gradUx_master_loc)
SDEALLOCATE(gradUx_slave_loc)
SDEALLOCATE(gradUy_master_loc)
SDEALLOCATE(gradUy_slave_loc)
SDEALLOCATE(gradUz_master_loc)
SDEALLOCATE(gradUz_slave_loc)
!SDEALLOCATE(U_local)
!SDEALLOCATE(gradp_local)
#endif /* USE_FAXEN_CORR */

! Basset force
#if USE_BASSETFORCE
SDEALLOCATE(durdt)
SDEALLOCATE(bIter)
SDEALLOCATE(Fbdt)
SDEALLOCATE(FbCoeff)
!SDEALLOCATE(FbCoeffa)
!SDEALLOCATE(FbCoefft)
!SDEALLOCATE(Fbi)
#endif /* USE_BASSETFORCE */

#if ANALYZE_RHS
SDEALLOCATE(Pt_ext)
#endif /* ANALYZE_RHS */

! interpolation
CALL FinalizeParticleInterpolation

! random walk
#if USE_RW
CALL ParticleFinalizeRandomWalk()
#endif /*USE_RW*/

! time disk
CALL Particle_FinalizeTimeDisk()

!
SDEALLOCATE(Seeds)

! subgrid-scale model
#if PARABOLIC
CALL ParticleFinalizeSGS()
#endif

! particle impact tracking
CALL FinalizeParticleBoundaryTracking()

! particle surface sampling
CALL FinalizeParticleBoundarySampling()

#if USE_MPI
! particle MPI halo exchange
CALL FinalizePartExchangeProcs()
#endif

CALL FinalizeParticleAnalyze()
CALL FinalizeParticleSurfaces()
CALL FinalizeParticleMesh()

#if PARTICLES_COUPLING >= 2
IF (doCalcPartSource) CALL FinalizeDeposition()
#endif /*PARTICLES_COUPLING >= 2*/

#if PARTICLES_COUPLING == 4
IF (doCalcPartCollision) CALL FinalizeCollision()
#endif /*PARTICLES_COUPLING == 4*/

ParticlesInitIsDone = .FALSE.

END SUBROUTINE FinalizeParticles


! !===================================================================================================================================
! SUBROUTINE rotx(mat,a)
! ! MODULES                                                                                                                          !
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! ! INPUT VARIABLES
! !----------------------------------------------------------------------------------------------------------------------------------!
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! REAL, INTENT(OUT), DIMENSION(3,3) :: mat
! REAL, INTENT(IN)                  :: a
! !===================================================================================================================================
! mat(:,1)=(/1.0 , 0.     , 0.  /)
! mat(:,2)=(/0.0 , COS(a) ,-SIN(a)/)
! mat(:,3)=(/0.0 , SIN(a) , COS(a)/)
! END SUBROUTINE
!
!
! !===================================================================================================================================
! SUBROUTINE roty(mat,a)
! ! MODULES                                                                                                                          !
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! ! INPUT VARIABLES
! !----------------------------------------------------------------------------------------------------------------------------------!
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! REAL, INTENT(OUT), DIMENSION(3,3) :: mat
! REAL, INTENT(IN)                  :: a
! !===================================================================================================================================
! mat(:,1)=(/ COS(a) , 0., SIN(a)/)
! mat(:,2)=(/ 0.     , 1., 0.  /)
! mat(:,3)=(/-SIN(a) , 0., COS(a)/)
! END SUBROUTINE
!
!
! !===================================================================================================================================
! SUBROUTINE ident(mat)
! ! MODULES                                                                                                                          !
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! ! INPUT VARIABLES
! !----------------------------------------------------------------------------------------------------------------------------------!
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! REAL, INTENT(OUT), DIMENSION(3,3) :: mat
! INTEGER                           :: j
! !===================================================================================================================================
! mat                      = 0.
! FORALL(j = 1:3) mat(j,j) = 1.
! END SUBROUTINE


SUBROUTINE InitRandomSeed(nRandomSeeds,SeedSize,Seeds)
!===================================================================================================================================
!> Initialize pseudo random numbers: Create Random_seed array
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Globals,               ONLY:myRank
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nRandomSeeds
INTEGER,INTENT(IN)             :: SeedSize
INTEGER,INTENT(INOUT)          :: Seeds(SeedSize)
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                        :: iSeed,DateTime(8),ProcessID,iStat,OpenFileID,GoodSeeds
INTEGER(KIND=8)                :: Clock,AuxilaryClock
LOGICAL                        :: uRandomExists
!==================================================================================================================================

uRandomExists=.FALSE.
IF (nRandomSeeds.NE.-1) THEN
  Clock     = 1536679165842_8
  ProcessID = 3671
ELSE
! First try if the OS provides a random number generator
  OPEN(NEWUNIT=OpenFileID, FILE="/dev/urandom", ACCESS="stream", &
       FORM="unformatted", ACTION="read", STATUS="old", IOSTAT=iStat)
  IF (iStat.EQ.0) THEN
    READ(OpenFileID) Seeds
    CLOSE(OpenFileID)
    uRandomExists=.TRUE.
  ELSE
    ! Fallback to XOR:ing the current time and pid. The PID is
    ! useful in case one launches multiple instances of the same
    ! program in parallel.
    CALL SYSTEM_CLOCK(COUNT=Clock)
    IF (Clock .EQ. 0) THEN
      CALL DATE_AND_TIME(values=DateTime)
      Clock =(DateTime(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
      + DateTime(2) * 31_8 * 24 * 60 * 60 * 1000 &
      + DateTime(3) * 24_8 * 60 * 60 * 1000 &
      + DateTime(5) * 60 * 60 * 1000 &
      + DateTime(6) * 60 * 1000 &
      + DateTime(7) * 1000 &
      + DateTime(8)
    END IF
    ProcessID = GetPID_C()
  END IF
END IF

IF(.NOT. uRandomExists) THEN
  Clock = IEOR(Clock, INT(ProcessID, KIND(Clock)))
  AuxilaryClock=Clock
  DO iSeed = 1, SeedSize
#if USE_MPI
    IF (nRandomSeeds.EQ.0) THEN
      AuxilaryClock=AuxilaryClock+myRank
    ELSE IF(nRandomSeeds.GT.0) THEN
      AuxilaryClock=AuxilaryClock+(myRank+1)*INT(Seeds(iSeed),8)*37
    END IF
#else
    IF (nRandomSeeds.GT.0) THEN
      AuxilaryClock=AuxilaryClock+INT(Seeds(iSeed),8)*37
    END IF
#endif
    IF (AuxilaryClock .EQ. 0) THEN
      AuxilaryClock = 104729
    ELSE
      AuxilaryClock = MOD(AuxilaryClock, 4294967296_8)
    END IF
    AuxilaryClock = MOD(AuxilaryClock * 279470273_8, 4294967291_8)
    GoodSeeds = INT(MOD(AuxilaryClock, INT(HUGE(0),KIND=8)), KIND(0))
    Seeds(iSeed) = GoodSeeds
  END DO
END IF
CALL RANDOM_SEED(PUT=Seeds)

END SUBROUTINE InitRandomSeed


SUBROUTINE LoadWeightsANN(Filename_loc)
!===================================================================================================================================
! Read in parameters of ANN
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: PartBoundANN
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(255),INTENT(IN)     :: Filename_loc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: tmp(2),iLayer
!===================================================================================================================================

IF (MPIRoot) THEN
  !Filename_loc = TRIM(REPLACE(Filename_loc," ","",Every=.TRUE.))
  OPEN(34, FILE=TRIM(FileName_loc),FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
  ! Read number of layers
  READ(34) PartBoundANN%nLayer
  ! Read maximum number of neurons
  READ(34) tmp(1)
  PartBoundANN%hgs = MERGE(tmp(1), PartBoundANN%hgs, tmp(1).LE.PartBoundANN%hgs)
  ! Allocate arrays and nullify
  ALLOCATE(PartBoundANN%nN(PartBoundANN%nLayer+2))
  ALLOCATE(PartBoundANN%w(1:PartBoundANN%hgs,1:PartBoundANN%hgs,1:PartBoundANN%nLayer+1))
  ALLOCATE(PartBoundANN%b(1:PartBoundANN%hgs,1:PartBoundANN%nLayer+1))
  ALLOCATE(PartBoundANN%beta(1:PartBoundANN%hgs,1:PartBoundANN%nLayer))
  PartBoundANN%w = 0.
  PartBoundANN%b = 0.
  PartBoundANN%beta = 0.
  ! Read weights
  DO iLayer=1,PartBoundANN%nLayer+1
    READ(34) tmp
    PartBoundANN%nN(iLayer) = tmp(1)
    READ(34) PartBoundANN%w(1:tmp(1),1:tmp(2),iLayer)
  END DO
  PartBoundANN%nN(PartBoundANN%nLayer+2) = tmp(2)
  DO iLayer=1,PartBoundANN%nLayer+1
    READ(34) tmp(1)
    READ(34) PartBoundANN%b(1:tmp(1),iLayer)
  END DO
  DO iLayer=1,PartBoundANN%nLayer
    READ(34) tmp(1)
    READ(34) PartBoundANN%beta(1:tmp(1),iLayer)
  END DO
  ALLOCATE(PartBoundANN%max_in(1:PartBoundANN%nN(1)))
  READ(34) PartBoundANN%max_in
  ALLOCATE(PartBoundANN%max_out(1:PartBoundANN%nN(PartBoundANN%nLayer+2)))
  READ(34) PartBoundANN%max_out
  CLOSE(34)
END IF
#if USE_MPI
CALL MPI_BCAST(PartBoundANN%nLayer,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
IF (.NOT.ALLOCATED(PartBoundANN%nN)) ALLOCATE(PartBoundANN%nN(1:PartBoundANN%nLayer+2))
CALL MPI_BCAST(PartBoundANN%nN,PartBoundANN%nLayer+2,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
tmp(1) = MAXVAL(PartBoundANN%nN)
PartBoundANN%hgs = MERGE(tmp(1), PartBoundANN%hgs, tmp(1).LE.PartBoundANN%hgs)
ALLOCATE(PartBoundANN%output(1:PartBoundANN%hgs))
PartBoundANN%output = 0.
! Allocate arrays and nullify
IF (.NOT.ALLOCATED(PartBoundANN%w)) THEN
  ALLOCATE(PartBoundANN%w(1:PartBoundANN%hgs,1:PartBoundANN%hgs,1:PartBoundANN%nLayer+1))
  ALLOCATE(PartBoundANN%b(1:PartBoundANN%hgs,1:PartBoundANN%nLayer+1))
  ALLOCATE(PartBoundANN%beta(1:PartBoundANN%hgs,1:PartBoundANN%nLayer))
  ALLOCATE(PartBoundANN%max_in(1:PartBoundANN%nN(1)))
  ALLOCATE(PartBoundANN%max_out(1:PartBoundANN%nN(PartBoundANN%nLayer+2)))
END IF
CALL MPI_BCAST(PartBoundANN%w,PartBoundANN%hgs*PartBoundANN%hgs*(PartBoundANN%nLayer+1),MPI_FLOAT,0,MPI_COMM_FLEXI,iError)
CALL MPI_BCAST(PartBoundANN%b,PartBoundANN%hgs*(PartBoundANN%nLayer+1),MPI_FLOAT,0,MPI_COMM_FLEXI,iError)
CALL MPI_BCAST(PartBoundANN%beta,PartBoundANN%hgs*PartBoundANN%nLayer,MPI_FLOAT,0,MPI_COMM_FLEXI,iError)
CALL MPI_BCAST(PartBoundANN%max_in,PartBoundANN%nN(1),MPI_FLOAT,0,MPI_COMM_FLEXI,iError)
CALL MPI_BCAST(PartBoundANN%max_out,PartBoundANN%nN(PartBoundANN%nLayer+2),MPI_FLOAT,0,MPI_COMM_FLEXI,iError)
#endif

END SUBROUTINE LoadWeightsANN

END MODULE MOD_Particle_Init
