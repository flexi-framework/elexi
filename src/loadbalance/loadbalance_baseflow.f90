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
!> Module contains the routines for load balancing
!> This file is needed only in FLEXI to avoid a circular depencendy
!===================================================================================================================================
MODULE MOD_LoadBalance_BaseFlow
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: BaseFlowRestart
PUBLIC:: SpongeRestart
!===================================================================================================================================

CONTAINS

SUBROUTINE BaseFlowRestart()
!===================================================================================================================================
! routine perfoming the field restart
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_BaseFlow_Vars          ,ONLY: doBaseFlow,BaseFlow,BaseFlowFiltered
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_Mesh_Vars              ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                :: BaseFlowTmp(:,:,:,:,:)
REAL,ALLOCATABLE                :: BaseFlowFilteredTmp(:,:,:,:,:)
! Timers
REAL                            :: StartT,EndT
!===================================================================================================================================

IF (.NOT. doBaseFlow) RETURN

StartT = MPI_WTIME()
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' REDISTRIBUTING BASEFLOW DURING LOADBALANCE...'

ALLOCATE(BaseFlowTmp(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ASSOCIATE (&
        counts_send  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemSend     ) ,&
        disp_send    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemSend) ,&
        counts_recv  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemRecv     ) ,&
        disp_recv    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemRecv))
  ! Communicate PartInt over MPI
  CALL MPI_ALLTOALLV(BaseFlow,counts_send,disp_send,MPI_DOUBLE_PRECISION,BaseFlowTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
END ASSOCIATE
CALL MOVE_ALLOC(BaseFlowTmp,BaseFlow)

ALLOCATE(BaseFlowFilteredTmp(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ASSOCIATE (&
        counts_send  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemSend     ) ,&
        disp_send    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemSend) ,&
        counts_recv  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemRecv     ) ,&
        disp_recv    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemRecv))
  ! Communicate PartInt over MPI
  CALL MPI_ALLTOALLV(BaseFlowFiltered,counts_send,disp_send,MPI_DOUBLE_PRECISION,BaseFlowFilteredTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
END ASSOCIATE
CALL MOVE_ALLOC(BaseFlowFilteredTmp,BaseFlowFiltered)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.TRUE.)

END SUBROUTINE BaseFlowRestart


!==================================================================================================================================
!> In case of load balancing, all dimensions match. Only shift the solution along the SFC!
!==================================================================================================================================
SUBROUTINE SpongeRestart()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Sponge_Vars,        ONLY: SpRefState
USE MOD_LoadBalance_Vars,   ONLY: PerformLoadBalance
USE MOD_LoadBalance_Vars,   ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: StartT,EndT
REAL,ALLOCATABLE              :: SpRefStateTmp(:,:,:,:,:)
! ==================================================================================================================================

StartT = MPI_WTIME()
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' REDISTRIBUTING SPONGE REFERENCE STATE DURING LOADBALANCE...'

ALLOCATE(SpRefStateTmp(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ASSOCIATE (&
        counts_send  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemSend     ) ,&
        disp_send    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemSend) ,&
        counts_recv  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemRecv     ) ,&
        disp_recv    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemRecv))
  ! Communicate PartInt over MPI
  CALL MPI_ALLTOALLV(SpRefState,counts_send,disp_send,MPI_DOUBLE_PRECISION,SpRefStateTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
END ASSOCIATE
CALL MOVE_ALLOC(SpRefStateTmp,SpRefState)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.TRUE.)

END SUBROUTINE SpongeRestart

END MODULE MOD_LoadBalance_BaseFlow
