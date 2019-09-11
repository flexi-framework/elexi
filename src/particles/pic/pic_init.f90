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
! Includes PIC Init
!===================================================================================================================================
MODULE MOD_PICInit
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitPIC
  MODULE PROCEDURE InitPIC
END INTERFACE

INTERFACE DefineParametersPIC
  MODULE PROCEDURE DefineParametersPIC
END INTERFACE


PUBLIC::InitPIC
PUBLIC::DefineParametersPIC

!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for PIC
!==================================================================================================================================
SUBROUTINE DefineParametersPIC()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("PIC")

CALL prms%CreateStringOption(   'PIC-Interpolation-Type'      , "Type of Interpolation-Method to calculate"                      //&
                                                                " the EM field's value for the particle", 'particle_position')
CALL prms%CreateLogicalOption(  'PIC-InterpolationElemLoop'   , 'Interpolate with outer iElem-loop (not'                         //&
                                                                'for many Elems per proc!)', '.TRUE.')
CALL prms%CreateRealArrayOption('PIC-externalField'           , 'External field is added to the'                                 //&
                                                                'maxwell-solver-field', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'PIC-scaleexternalField'      , '', '1.0')
CALL prms%CreateLogicalOption(  'PIC-DoInterpolation'         , "Compute the self field's influence on the Particle", '.TRUE.')
CALL prms%CreateLogicalOption(  'PIC-BG-Field'                , 'BGField data (1:x,0:NBG,0:NBG,0:NBG,'                           //&
                                                                '1:PP_nElems) \n'                                                //&
                                                                'If PIC-BG-Field=T\n'                                            //&
                                                                'Define:\n'                                                      //&
                                                                'PIC-BGFilename\n'                                               //&
                                                                'PIC-BGFieldScaling\n'                                           //&
                                                                'PIC-NBG', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-BGFileName'              , 'File name for background field ([character].h5)', 'none')
CALL prms%CreateIntOption(      'PIC-NBG'                     , 'Polynomial degree that shall be used '                          //&
                                                                'for background field during simulation', '1')
CALL prms%CreateRealOption(     'PIC-BGFieldScaling'          , 'Space scaling of background field', '1.')
CALL prms%CreateStringOption(   'PIC-curvedexternalField'     , 'File to curved external field data.','none')
CALL prms%CreateStringOption(   'PIC-variableexternalField'   , 'File containing the external field CSV table', 'none')

CALL prms%CreateRealArrayOption('PIC-NormVecOfWall'  , 'Normal vector for pushTimeStep', '1. , 0. , 0.')
CALL prms%CreateIntOption(      'PIC-DeltaType'      , 'Flag ', '1')
CALL prms%CreateIntOption(      'PIC-DeltaType-N'    , 'Polynomial degree of delta distribution', '1')
CALL prms%CreateRealArrayOption('PIC-BGMdeltas'      , 'Dimensions of PIC background mesh', '0. , 0. , 0.')
CALL prms%CreateRealArrayOption('PIC-FactorBGM'      , 'Denominator of PIC-BGMdeltas', '1. , 1. , 1.')

CALL prms%SetSection("PIC Deposition")

!CALL prms%CreateLogicalOption(  'PIC-DoDeposition'         , 'Switch deposition on/off', '.TRUE.')
!CALL prms%CreateStringOption(   'PIC-Deposition-Type'      , '(HALOWIKI:)\n'                                                     //&
!                                                                    'If Deposition-Type=shape_function\n'                        //&
!                                                             'Define:\n'                                                         //&
!                                                             'PIC-shapefunction-radius\n'                                        //&
!                                                             'PIC-shapefunction-alpha.\n'                                        //&
!                                                             'If Deposition-Type =(cartmesh_volumeweighting/ cartmesh_splines)\n'//&
!                                                             'Define:\n'                                                         //&
!                                                             'PIC-BGMdeltas\n'                                                   //&
!                                                             'PIC-FactorBGM', 'nearest-blurrycenter')
!CALL prms%CreateStringOption(   'PIC-TimeAverageFile'      , 'TODO-DEFINE-PARAMETER', 'none')
!
!CALL prms%CreateRealOption(     'PIC-epanechnikov-radius'  , 'TODO-DEFINE-PARAMETER', '1.')
!CALL prms%CreateRealOption(     'PIC-shapefunction-radius' , 'Radius of shape function', '1.')
!CALL prms%CreateIntOption(      'PIC-shapefunction-alpha'  , 'Exponent of shape function', '2')
!CALL prms%CreateLogicalOption(  'PIC-shapefunction-equi'   , 'Use equidistant points for shapefunction', '.FALSE.')
!CALL prms%CreateIntOption(      'PIC-shapefunction1d-direction'  , 'Direction of 1D shape function', '1')
!CALL prms%CreateRealOption(     'PIC-shapefunction-radius0', 'Minimal shape function radius', '1.')
!CALL prms%CreateRealOption(     'PIC-shapefunction-scale'  , 'Scaling factor of shape function radius', '0.')
!CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoFixes'     , 'Number of fixes for shape func depo at planar BCs', '0')
!CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoLayers'  ,    'Number of const. source layer for sf-depo at planar BCs', '0')

END SUBROUTINE DefineParametersPIC

SUBROUTINE InitPIC()
!===================================================================================================================================
! PIC Init
!===================================================================================================================================
! MODULES
USE MOD_Globals
!USE MOD_PICInterpolation_Vars,  ONLY: externalField
USE MOD_PIC_Vars ,              ONLY: PICInitIsDone!, PIC
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(PICInitIsDone)THEN
   SWRITE(*,*) "InitPIC already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PIC ...'

! So far, nothing to do here...
!IF (externalField(6).NE.0) PIC%GyroVecDirSIGN = -externalField(6)/(ABS(externalField(6)))

PICInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PIC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPIC

END MODULE MOD_PICInit
