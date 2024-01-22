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
! Here, preprocessor variables for different equation systems and abbreviations for specific expressions are defined
!===================================================================================================================================

! Abbrevations
#ifdef SUN
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SUN COMPILER'
#  define IEEE_ISNAN
#elif SX
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SX COMPILER'
#elif PGI
#  define NO_ISNAN
#endif
#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#define SIZEOF_F(x) (STORAGE_SIZE(x)/8)
#define NO_OP(x)    ASSOCIATE( x => x ); END ASSOCIATE

! Size of data types
#define SIZE_LOG  KIND(.TRUE.)
#define SIZE_INT  KIND(INT(1))
#define SIZE_INT4 4
#define SIZE_INT8 8
#define SIZE_REAL KIND(REAL(1))
#define SIZE_CHAR KIND('a')
#ifdef INTKIND8
#define MPI_INTEGER_INT_KIND MPI_INTEGER8
#else
#define MPI_INTEGER_INT_KIND MPI_INTEGER
#endif

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(INT( 1,KIND=k)).OR.x<-HUGE(INT( 1,KIND=k))) CALL Abort(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(REAL(1,KIND=k)).OR.x<-HUGE(REAL(1,KIND=k))) CALL Abort(__STAMP__,'Real conversion failed: out of range!')
#elif CRAY
#define CHECKSAFEINT(x,k)
#define CHECKSAFEREAL(x,k)
#else
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  ## k).OR.x<-HUGE(1_  ## k)) CALL Abort(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ ## k).OR.x<-HUGE(1._ ## k)) CALL Abort(__STAMP__,'Real conversion failed: out of range!')
#endif

! Time Step Minimum: dt_Min
#define DT_NVAR       3
#define DT_MIN        1
#define DT_ANALYZE    2
#define DT_END        3

! Test for equality: read description in mathtools.f90 for further infos
#define ALMOSTEQUALABSOLUTE(x,y,tol)  (ABS((x)-(y)).LE.(tol))
#define ALMOSTEQUALRELATIVE(x,y,tol)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(tol))
#define ALMOSTEQUALABSORREL(x,y,tol)  (ALMOSTEQUALABSOLUTE(x,y,tol) .OR.  ALMOSTEQUALRELATIVE(x,y,tol))
#define ALMOSTEQUALABSANDREL(x,y,tol) (ALMOSTEQUALABSOLUTE(x,y,tol) .AND. ALMOSTEQUALRELATIVE(x,y,tol))

! Define MPI specific write shortcuts
#if USE_MPI
#  define SWRITE IF(MPIRoot) WRITE
#if USE_LOADBALANCE
#  define LBWRITE IF(MPIRoot.AND..NOT.PerformLoadBalance) WRITE
#else /*USE_LOADBALANCE*/
#  define LBWRITE IF(MPIRoot) WRITE
#endif /*USE_LOADBALANCE*/
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#  define GETTIME(a) a=MPI_WTIME()
#else /*USE_MPI*/
#  define SWRITE WRITE
#  define LBWRITE WRITE
#  define IPWRITE WRITE
#  define GETTIME(a) CALL CPU_TIME(a)
#endif /*USE_MPI*/
#define ERRWRITE(a,b)  CALL CreateErrFile(); IF(ErrorFiles) WRITE(UNIT_errOut,b)
#define LOGWRITE(a,b)  IF(Logging) WRITE(UNIT_logOut,b)
#define SDEALLOCATE(A) IF(ALLOCATED(A))  DEALLOCATE(A)
#define ADEALLOCATE(A) IF(ASSOCIATED(A)) DEALLOCATE(A)
#define PDEALLOCATE(A) IF(PRESENT(A)) THEN; IF(ALLOCATED(A)) DEALLOCATE(A); ENDIF

! Shared Memory
#if USE_MPI
#define ALLOCPOINT POINTER
#define MDEALLOCATE(A) IF(ASSOCIATED(A)) NULLIFY(A)
#else
#define ALLOCPOINT ALLOCATABLE
#define MDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)
#endif

! Debug memory
#if DEBUG_MEMORY
#define Allocate_Shared(a,b,c)   Allocate_Shared_DEBUG(a,b,c,'b')
#endif
#if USE_MPI
#define LWRITE IF(myComputeNodeRank.EQ.0) WRITE
#else
#define LWRITE WRITE
#endif

! Define OpenMP specific shortcuts
#if USE_OPENMP
#  define OMP_FLEXITIME() OMP_GET_WTIME()
#else
#  define OMP_FLEXITIME() FLEXITIME()
#endif

! Loop variables
#define PP_IJK     i,j,k
#define PP_ij      i,j

! Predefined "PARAMETER-like" variables
#define XI_MINUS   5
#define XI_PLUS    3
#define ETA_MINUS  2
#define ETA_PLUS   4
#define ZETA_MINUS 1
#define ZETA_PLUS  6

! Entry position in SideToElem
#define S2E_ELEM_ID        1
#define S2E_NB_ELEM_ID     2
#define S2E_LOC_SIDE_ID    3
#define S2E_NB_LOC_SIDE_ID 4
#define S2E_FLIP           5

! Entry position in ElemToSide
#define E2S_SIDE_ID 1
#define E2S_FLIP    2

! Entry position in BC
#define MI_SIDEID 1
#define MI_FLIP   2

! Entry position in BC
#define BC_TYPE  1
#define BC_STATE 2
#define BC_ALPHA 3

! Entry position in BC
#define SEND 1
#define RECV 2

! Shared memory mesh arrays
! number of entry in each line of ElemInfo
#define ELEMINFOSIZE_H5   6
!#if USE_MPI
#define ELEMINFOSIZE      9
!#else
!#define ELEMINFOSIZE      6
!#endif /* USE_MPI*/
! ElemInfo in H5 file
#define ELEM_TYPE         1
#define ELEM_ZONE         2
#define ELEM_FIRSTSIDEIND 3
#define ELEM_LASTSIDEIND  4
#define ELEM_FIRSTNODEIND 5
#define ELEM_LASTNODEIND  6
! ElemInfo for shared array
#define ELEM_RANK         7
#define ELEM_HALOFLAG     8
#define ELEM_HASMORTAR    9

! number of entries in each line of SideInfo
#define SIDEINFOSIZE_H5   5
#define SIDEINFOSIZE      8
#define SIDE_TYPE         1
#define SIDE_ID           2
#define SIDE_NBELEMID     3
#define SIDE_FLIP         4
#define SIDE_BCID         5
#define SIDE_ELEMID       6
#define SIDE_LOCALID      7
#define SIDE_NBSIDEID     8
#define SIDE_NBELEMTYPE   9

! number of entries in each line of FEMInfo
#define FEMELEMINFOSIZE   4
#define VERTEXINFOSIZE    3

#if USE_LOADBALANCE
! Load Balance (LB) position in array for measuring the time that is spent on specific operations
#define LB_DG            1
#define LB_DGANALYZE     2
#define LB_DGCOMM        3
#define LB_FV            4
#define LB_FVCOMM        5
#define LB_EMISSION      6
#define LB_TRACK         7
#define LB_INTERPOLATION 8
#define LB_PUSH          9
#define LB_PARTCOMM      10
#define LB_SURFFLUX      11

#define LB_NTIMES        11

! Macros for LB timings
#define MeasureStartTime()               CALL LBStartTime(tLBStart)
#define MeasurePauseTime_DG()            CALL LBPauseTime(LB_DG           ,tLBStart)
#define MeasurePauseTime_DGANALYZE()     CALL LBPauseTime(LB_DGANALYZE    ,tLBStart)
#define MeasurePauseTime_DGCOMM()        CALL LBPauseTime(LB_DGCOMM       ,tLBStart)
#define MeasurePauseTime_FV()            CALL LBPauseTime(LB_FV           ,tLBStart)
#define MeasurePauseTime_FVCOMM()        CALL LBPauseTime(LB_FVCOMM       ,tLBStart)
#define MeasurePauseTime_EMISSION()      CALL LBPauseTime(LB_EMISSION     ,tLBStart)
#define MeasurePauseTime_TRACK()         CALL LBPauseTime(LB_TRACK        ,tLBStart)
#define MeasurePauseTime_INTERPOLATION() CALL LBPauseTime(LB_INTERPOLATION,tLBStart)
#define MeasurePauseTime_PUSH()          CALL LBPauseTime(LB_PUSH         ,tLBStart)
#define MeasurePauseTime_PARTCOMM()      CALL LBPauseTime(LB_PARTCOMM     ,tLBStart)
#define MeasurePauseTime_SURFFLUX()      CALL LBPauseTime(LB_SURFFLUX     ,tLBStart)
#define MeasureSplitTime_DG()            CALL LBSplitTime(LB_DG           ,tLBStart)
#define MeasureSplitTime_DGANALYZE()     CALL LBSplitTime(LB_DGANALYZE    ,tLBStart)
#define MeasureSplitTime_DGCOMM()        CALL LBSplitTime(LB_DGCOMM       ,tLBStart)
#define MeasureSplitTime_FV()            CALL LBSplitTime(LB_FV           ,tLBStart)
#define MeasureSplitTime_FVCOMM()        CALL LBSplitTime(LB_FVCOMM       ,tLBStart)
#define MeasureSplitTime_EMISSION()      CALL LBSplitTime(LB_EMISSION     ,tLBStart)
#define MeasureSplitTime_TRACK()         CALL LBSplitTime(LB_TRACK        ,tLBStart)
#define MeasureSplitTime_INTERPOLATION() CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#define MeasureSplitTime_PUSH()          CALL LBSplitTime(LB_PUSH         ,tLBStart)
#define MeasureSplitTime_PARTCOMM()      CALL LBSplitTime(LB_PARTCOMM     ,tLBStart)
#define MeasureSplitTime_SURFFLUX()      CALL LBSplitTime(LB_SURFFLUX     ,tLBStart)
#define MeasureElemPauseTime()           CALL LBElemPauseTime(ElemID      ,tLBStart)
#define MeasureElemSplitTime()           CALL LBElemSplitTime(ElemID      ,tLBStart)
#else
! Macros for LB timings
! #define MeasureStartTime()               CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_DG()            CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_DGANALYZE()     CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_DGCOMM()        CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_FV()            CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_FVCOMM()        CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_EMISSION()      CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_TRACK()         CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_INTERPOLATION() CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_PUSH()          CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_PARTCOMM()      CALL NOP_INLINE(0) ! do nothing
! #define MeasurePauseTime_SURFFLUX()      CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_DG()            CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_DGANALYZE()     CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_DGCOMM()        CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_FV()            CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_FVCOMM()        CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_EMISSION()      CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_TRACK()         CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_INTERPOLATION() CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_PUSH()          CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_PARTCOMM()      CALL NOP_INLINE(0) ! do nothing
! #define MeasureSplitTime_SURFFLUX()      CALL NOP_INLINE(0) ! do nothing
! #define MeasureElemPauseTime()           CALL NOP_INLINE(0) ! do nothing
! #define MeasureElemSplitTime()           CALL NOP_INLINE(0) ! do nothing
#define MeasureStartTime()               ! do nothing
#define MeasurePauseTime_DG()            ! do nothing
#define MeasurePauseTime_DGANALYZE()     ! do nothing
#define MeasurePauseTime_DGCOMM()        ! do nothing
#define MeasurePauseTime_FV()            ! do nothing
#define MeasurePauseTime_FVCOMM()        ! do nothing
#define MeasurePauseTime_EMISSION()      ! do nothing
#define MeasurePauseTime_TRACK()         ! do nothing
#define MeasurePauseTime_INTERPOLATION() ! do nothing
#define MeasurePauseTime_PUSH()          ! do nothing
#define MeasurePauseTime_PARTCOMM()      ! do nothing
#define MeasurePauseTime_SURFFLUX()      ! do nothing
#define MeasureSplitTime_DG()            ! do nothing
#define MeasureSplitTime_DGANALYZE()     ! do nothing
#define MeasureSplitTime_DGCOMM()        ! do nothing
#define MeasureSplitTime_FV()            ! do nothing
#define MeasureSplitTime_FVCOMM()        ! do nothing
#define MeasureSplitTime_EMISSION()      ! do nothing
#define MeasureSplitTime_TRACK()         ! do nothing
#define MeasureSplitTime_INTERPOLATION() ! do nothing
#define MeasureSplitTime_PUSH()          ! do nothing
#define MeasureSplitTime_PARTCOMM()      ! do nothing
#define MeasureSplitTime_SURFFLUX()      ! do nothing
#define MeasureElemPauseTime()           ! do nothing
#define MeasureElemSplitTime()           ! do nothing
#endif /*USE_LOADBALANCE*/

!#define DEBUGMESH

#if FV_ENABLED
#define FV_SIZE 1
#else
#define FV_SIZE 0
#endif

#if !(FV_ENABLED)
#define FV_Elems(x) 0
#define FV_Elems_master(x) 0
#define FV_Elems_slave(x) 0
#define FV_Elems_Sum(x) 0
#endif

! Compute viscous contributions in volume integral
! NOT if FV-Blending or if non-parabolic
#if (FV_ENABLED==2) || !PARABOLIC
#define VOLINT_VISC 0
#else
#define VOLINT_VISC 1
#endif

#define KILL(x) SWRITE(*,*) __FILE__,__LINE__,x; stop

! overintegration
#define OVERINTEGRATIONTYPE_NONE       0
#define OVERINTEGRATIONTYPE_CUTOFF     1
#define OVERINTEGRATIONTYPE_CONSCUTOFF 2

! filter
#define FILTERTYPE_NONE   0
#define FILTERTYPE_CUTOFF 1
#define FILTERTYPE_MODAL  2
#define FILTERTYPE_LAF    3

! PURE debug switch
#if DEBUG
#define PPURE
#else
#define PPURE PURE
#endif

!2d functionality
#if (PP_dim==2)
#define ZDIM(a) 0
#define PP_NZ   0
#define DIMV    1:2
#else
#define ZDIM(a) a
#define PP_NZ   PP_N
#define DIMV    1:3
#endif
