!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
#include "particle.h"

!===================================================================================================================================
!> Module contains the routines for load balancing
!> This file is needed only in FLEXI to avoid a circular depencendy
!===================================================================================================================================
MODULE MOD_LoadBalance_Restart
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_LOADBALANCE
INTERFACE FieldRestart
  MODULE PROCEDURE FieldRestart
END INTERFACE

PUBLIC :: FieldRestart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

CONTAINS

#if USE_LOADBALANCE
SUBROUTINE FieldRestart()
!===================================================================================================================================
! routine perfoming the field restart
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars                ,ONLY: U
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
USE MOD_LoadBalance_Vars       ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_Mesh_Vars              ,ONLY: nElems
#if FV_ENABLED
USE MOD_FV                     ,ONLY: FV_ProlongFVElemsToFace
USE MOD_FV_Vars                ,ONLY: FV_Elems
USE MOD_Indicator_Vars         ,ONLY: IndValue
#endif /*FV_ENABLED*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                :: UTmp(:,:,:,:,:)
#if FV_ENABLED
REAL,ALLOCATABLE                :: IndValueTmp(:)
INTEGER,ALLOCATABLE             :: FV_ElemsTmp(:)
#endif /*FV_ENABLED*/
!===================================================================================================================================

ALLOCATE(UTmp(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ASSOCIATE (&
        counts_send  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemSend     ) ,&
        disp_send    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemSend) ,&
        counts_recv  => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPInElemRecv     ) ,&
        disp_recv    => (PP_nVar*(PP_N+1)*(PP_N+1)*(PP_NZ+1)*MPIoffsetElemRecv))
  ! Communicate PartInt over MPI
  CALL MPI_ALLTOALLV(U,counts_send,disp_send,MPI_DOUBLE_PRECISION,UTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
END ASSOCIATE
CALL MOVE_ALLOC(UTmp,U)

#if FV_ENABLED
ALLOCATE(FV_ElemsTmp(nElems))
ASSOCIATE (&
        counts_send  => (MPInElemSend     ) ,&
        disp_send    => (MPIoffsetElemSend) ,&
        counts_recv  => (MPInElemRecv     ) ,&
        disp_recv    => (MPIoffsetElemRecv))
  ! Communicate PartInt over MPI
  CALL MPI_ALLTOALLV(FV_Elems,counts_send,disp_send,MPI_INTEGER,FV_ElemsTmp,counts_recv,disp_recv,MPI_INTEGER,MPI_COMM_FLEXI,iError)
END ASSOCIATE
CALL MOVE_ALLOC(FV_ElemsTmp,FV_Elems)

ALLOCATE(IndValueTmp(nElems))
ASSOCIATE (&
        counts_send  => (MPInElemSend     ) ,&
        disp_send    => (MPIoffsetElemSend) ,&
        counts_recv  => (MPInElemRecv     ) ,&
        disp_recv    => (MPIoffsetElemRecv))
  ! Communicate PartInt over MPI
  CALL MPI_ALLTOALLV(IndValue,counts_send,disp_send,MPI_DOUBLE_PRECISION,IndValueTmp,counts_recv,disp_recv,MPI_DOUBLE_PRECISION,MPI_COMM_FLEXI,iError)
END ASSOCIATE
CALL MOVE_ALLOC(IndValueTmp,IndValue)

! Update the FV face information
CALL FV_ProlongFVElemsToFace()
#endif /*FV_ENABLED*/

END SUBROUTINE FieldRestart
#endif /*USE_LOADBALANCE*/

END MODULE MOD_LoadBalance_Restart
