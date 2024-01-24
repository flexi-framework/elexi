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

!===================================================================================================================================
!> Module contains the timers for load balancing
!===================================================================================================================================
MODULE MOD_LoadBalance_Timers
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_LOADBALANCE
INTERFACE LBStartTime
  MODULE PROCEDURE LBStartTime
END INTERFACE

INTERFACE LBSplitTime
  MODULE PROCEDURE LBSplitTime
END INTERFACE

INTERFACE LBPauseTime
  MODULE PROCEDURE LBPauseTime
END INTERFACE

INTERFACE LBElemSplitTime
  MODULE PROCEDURE LBElemSplitTime
END INTERFACE

INTERFACE LBElemPauseTime
  MODULE PROCEDURE LBElemPauseTime
END INTERFACE

PUBLIC::LBStartTime
PUBLIC::LBSplitTime
PUBLIC::LBPauseTime
PUBLIC::LBElemSplitTime
PUBLIC::LBElemPauseTime
!===================================================================================================================================

CONTAINS

SUBROUTINE LBStartTime(tLBStart)
!===================================================================================================================================
!> calculates and sets start time for Loadbalance.
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_LoadBalance_Vars ,ONLY: PerformLBSample
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)       :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBStart = LOCALTIME() ! LB Time Start
END SUBROUTINE LBStartTime


SUBROUTINE LBSplitTime(LB_index,tLBStart)
!===================================================================================================================================
!> Splits the time and resets LB_start. Adds time to tcurrent(LB_index) for current proc
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_LoadBalance_Vars ,ONLY: PerformLBSample,tCurrent
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: LB_index
REAL,INTENT(INOUT)  :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(LB_index) = tCurrent(LB_index) + tLBEnd-tLBStart
tLBStart = tLBEnd !LOCALTIME() ! LB Time Start
END SUBROUTINE LBSplitTime


SUBROUTINE LBPauseTime(LB_index,tLBStart)
!===================================================================================================================================
!> calculates end time and adds time to tcurrent(LB_index) for current proc
!> does not reset tLBstart
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_LoadBalance_Vars ,ONLY: PerformLBSample,tCurrent
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: LB_index
REAL,INTENT(IN)     :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(LB_index) = tCurrent(LB_index) + tLBEnd - tLBStart
END SUBROUTINE LBPauseTime


SUBROUTINE LBElemSplitTime(ElemID,tLBStart)
!===================================================================================================================================
!> Measure particle-related times for specific elements. Splits the time and resets LB_start.
!> Adds time to Elemtime(ElemID) and ElemTimePart.
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,PerformLBSample
#if USE_PARTICLES
USE MOD_LoadBalance_Vars ,ONLY: ElemTimePart
#endif /*USE_PARTICLES*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: ElemID   ! Element ID
REAL,INTENT(INOUT)  :: tLBStart ! Time when timer was started
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd,DeltaTime
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd           = LOCALTIME() ! LB Time End
DeltaTime        = tLBEnd-tLBStart
ElemTime(ELemID) = ElemTime(ElemID) + DeltaTime
#if USE_PARTICLES
ElemTimePart     = ElemTimePart     + DeltaTime
#endif /*USE_PARTICLES*/
tLBStart         = tLBEnd !LOCALTIME() ! LB Time Start
END SUBROUTINE LBElemSplitTime


SUBROUTINE LBElemPauseTime(ElemID,tLBStart)
!===================================================================================================================================
!> Measure particle-related times for specific elements.
!> Calculates end time and adds time to Elemtime(ElemID) and ElemTimePart.
!> Does not reset tLBstart
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,PerformLBSample
#if USE_PARTICLES
USE MOD_LoadBalance_Vars ,ONLY: ElemTimePart
#endif /*USE_PARTICLES*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: ElemID
REAL,INTENT(IN)     :: tLBStart
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                :: tLBEnd,DeltaTime
!===================================================================================================================================
IF(.NOT. PerformLBSample) RETURN
tLBEnd = LOCALTIME() ! LB Time End
DeltaTime        = tLBEnd-tLBStart
ElemTime(ELemID) = ElemTime(ElemID) + DeltaTime
#if USE_PARTICLES
ElemTimePart     = ElemTimePart     + DeltaTime
#endif /*USE_PARTICLES*/
END SUBROUTINE LBElemPauseTime


FUNCTION LOCALTIME()
!===================================================================================================================================
! Calculates current local time (own function because of a laterMPI implementation)
!===================================================================================================================================
! MODULES
USE MPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: LocalTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
LocalTime = MPI_WTIME()
#else
CALL CPU_TIME(LocalTime)
#endif
END FUNCTION LOCALTIME

#endif /* LOADBALANCE */
END MODULE MOD_LoadBalance_Timers
