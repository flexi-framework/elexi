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
!================================================================================================================================= ! Here, preprocessor variables for particle subsystem are defined
!=================================================================================================================================

! Calculate GCC version
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

! Safe allocation
! #define SALLOCATE(A,S) CALL Allocate_Safe(A,S) */

! Boundaries for Particles
#define PLANAR_RECT    0
#define PLANAR_NONRECT 1
#define BILINEAR       2
#define PLANAR_CURVED  3
#define CURVED         4

! Entry position in ElemBCSides
#define ELEM_NBR_BCSIDES   1
#define ELEM_FIRST_BCSIDE  2

! Entry position in FIBGMToProc
#define FIBGM_FIRSTPROCIND 1
#define FIBGM_NPROCS       2
#define FIBGM_NLOCALPROCS  3

! Entry position in SideBCMetrics
#define BCSIDE_SIDEID      1
#define BCSIDE_ELEMID      2
#define BCSIDE_DISTANCE    3
#define BCSIDE_RADIUS      4
#define BCSIDE_ORIGINX     5
#define BCSIDE_ORIGINY     6
#define BCSIDE_ORIGINZ     7

! entries for PartExchange
#define EXCHANGE_PROC_SIZE 2
#define EXCHANGE_PROC_TYPE 1
#define EXCHANGE_PROC_RANK 2

#define NATIVE_ELEM_ID  1
#define NATIVE_PROC_ID  2
#define LOCAL_PROC_ID   3
!#define NATIVE_SIDE_ID  1
#define LOCAL_SEND_ID   4

! surface sampling entries
#define SURF_SIDEID       1
#define SURF_RANK         2
#define SURF_LEADER       3

! Tracking method
#define REFMAPPING    1
#define TRACING       2
#define TRIATRACKING  3

! RHS method
#define RHS_NONE           1
#define RHS_TRACER         2
#define RHS_TCONVERGENCE   3
#define RHS_HPCONVERGENCE  4
#define RHS_INERTIA        5
#define RHS_INERTIA_EULER  6
#define RHS_SGS1           7
#define RHS_SGS1_TRACER    8
#define RHS_SGS2           9

! Drag factor model
#define DF_PART_STOKES    1
#define DF_PART_SCHILLER  2
#define DF_PART_PUTNAM    3
#define DF_PART_HAIDER    4
#define DF_PART_HOELZER   5
#define DF_PART_LOTH      6
#define DF_PART_GANSER    7

! Deposition method
#define DEPO_CV           1
#define DEPO_CVLM         2
#define DEPO_STEP         3
#define DEPO_DIRAC        4
#define DEPO_SF_GAUSS     5
#define DEPO_SF_POLY      6

#if USE_EXTEND_RHS || USE_FAXEN_CORR
! Velocity and pressure for extended RHS
#define RHS_LIFTVARS    (/LIFT_VEL1,LIFT_VEL2,LIFT_VEL3/)
#define RHS_LIFT        3
#if USE_EXTEND_RHS
#define RHS_GRADVEL1    1
#define RHS_GRADVEL2    2
#define RHS_GRADVEL3    3
#define RHS_GRADVELV    RHS_GRADVEL1:RHS_GRADVEL3
#define RHS_dVELdt      4
#define RHS_dVEL1dt     1
#define RHS_dVEL2dt     2
#define RHS_dVEL3dt     3
#define RHS_dVELVdt     RHS_dVEL1dt:RHS_dVEL3dt
#if USE_FAXEN_CORR
#define RHS_LAPLACEVEL  5
#define RHS_GRAD        5
#define RHS_GRADVARS    (/1,2,3,4,5/)
#define RHS_LAPLACEVEL1 4
#define RHS_LAPLACEVEL2 5
#define RHS_LAPLACEVEL3 6
#define RHS_NVARS       6
#else
#define RHS_GRAD        4
#define RHS_GRADVARS    (/1,2,3,4/)
#define RHS_NVARS       3
#endif
#elif USE_FAXEN_CORR
#define RHS_LAPLACEVEL  1
#define RHS_GRAD        1
#define RHS_GRADVARS    (/1/)
#define RHS_LAPLACEVEL1 1
#define RHS_LAPLACEVEL2 2
#define RHS_LAPLACEVEL3 3
#define RHS_NVARS       3
#endif
#endif

#define PART_POS1       1
#define PART_POS2       2
#define PART_POS3       3
#define PART_POSV       PART_POS1:PART_POS3
#define PART_VEL1       4
#define PART_VEL2       5
#define PART_VEL3       6
#define PART_VELV       PART_VEL1:PART_VEL3
! Only defined for PartData_Shared
#define PART_OLD1       4
#define PART_OLD2       5
#define PART_OLD3       6
#define PART_OLDV       PART_OLD1:PART_OLD3
#if USE_PARTTEMP
#define PART_TEMP       7
#endif
#if USE_PARTROT
#if USE_PARTTEMP
#define PART_AMOM1      8
#else
#define PART_AMOM1      7
#endif
#define PART_AMOM2      PART_AMOM1+1
#define PART_AMOM3      PART_AMOM2+1
#define PART_AMOMV      PART_AMOM1:PART_AMOM3
#endif
#define PART_DIAM       PP_nVarPart
#if USE_SPHERICITY
#define PART_SPHE       PP_nVarPart-1
#endif

#define ENERGY_ROTATION(rho,dp,vor)   0.5*pi/60*rho*dp**5*DOT_PRODUCT(vor,vor)
#define ENERGY_KINETIC(rho,dp,v)      0.5*rho*pi/6.*dp**3*DOT_PRODUCT(v,v)
#define MOM_OF_INERTIA(rho,dp)        pi/60.*rho*dp**5
#define MASS_SPHERE(rho,dp)           rho*pi/6.*dp**3
#define DENS_SPHERE(mass,dp)          mass/(pi/6.*dp**3)
#define DIAM_SPHERE(rho,mass)         (mass/rho*6./pi)**(1./3.)
#define VOL_SPHERE(dp)                pi/6.*dp**3
#define VOL_SPHERE_INV(Vol)           (Vol*6./pi)**(1./3.)

! formats
! print to std out like  "    1.41421356237310E+000   -1.41421356237310E+000   -1.41421356237310E+000"
! (looks good and prevents the first digit of being a zero)
#define WRITEFORMAT '(ES25.14E3)'
! print to csv file like "0.1414213562373095E+001,-.1414213562373095E+001,-.1414213562373095E+001"
! (does not look that good but it saves disk space)
#define CSVFORMAT '(A1,E23.16E3)'
