!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "flexi.h"
#if USE_PARTICLES
#include "particle.h"
#endif /*USE_PARTICLES*/

MODULE MOD_LoadBalance_Metrics
!===================================================================================================================================
!> \brief This module contains routines for computing the geometries volume and surface metric terms.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: MoveCoords
PUBLIC:: MoveMetrics
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine rearranges the coordinates along the space-filling curve
!==================================================================================================================================
SUBROUTINE MoveCoords()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP,nElems
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                   :: Elem_xGP_LB(:,:,:,:,:)
! Custom data type
INTEGER                            :: MPI_LENGTH(1),MPI_TYPE(1),MPI_STRUCT
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
! Timer
REAL                               :: StartT,EndT,WallTime
!==================================================================================================================================

IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  ! SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' Shift volume coordinates during loadbalance...'
  GETTIME(StartT)

  ALLOCATE(Elem_xGP_LB   (3,0:PP_N,0:PP_N,0:PP_N,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate Elem_xGP over MPI
    MPI_LENGTH       = 3*(PP_N+1)*(PP_N+1)*(PP_NZ+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(Elem_xGP,counts_send,disp_send,MPI_STRUCT,Elem_xGP_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(Elem_xGP)
  CALL MOVE_ALLOC(Elem_xGP_LB,Elem_xGP)

  GETTIME(EndT)
  WallTime = EndT-StartT
  ! CALL DisplayMessageAndTime(WallTime,'DONE',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
END IF

END SUBROUTINE MoveCoords


SUBROUTINE MoveMetrics()
!===================================================================================================================================
!> This routine rearranges the metrics along the space-filling curve
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
! USE MOD_LoadBalance_Vars   ,ONLY: MPInSideSend,MPInSideRecv,MPIoffsetSideSend,MPIoffsetSideRecv
USE MOD_Mesh_Vars          ,ONLY: nElems,NGeoRef!,Elem_xGP,nSides
USE MOD_Mesh_Vars          ,ONLY: JaCL_N,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,XCL_N,dXCL_N
! USE MOD_Mesh_Vars          ,ONLY: Face_xGP,NormVec,TangVec1,TangVec2,SurfElem,Ja_Face
USE MOD_Mesh_Vars          ,ONLY: sJ!,detJac_Ref
USE MOD_Mesh_Vars          ,ONLY: detJac_N,detJac_Ref
#if USE_PARTICLES
USE MOD_Mesh_Vars          ,ONLY: NGeo
USE MOD_Particle_Mesh_Vars ,ONLY: XCL_NGeo,dXCL_NGeo
#endif /*USE_PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Elements
REAL,ALLOCATABLE                   :: XCL_N_LB(         :,:,:,:,:)
REAL,ALLOCATABLE                   :: dXCL_N_LB(      :,:,:,:,:,:)
REAL,ALLOCATABLE                   :: JaCL_N_LB(      :,:,:,:,:,:)
REAL,ALLOCATABLE                   :: Metrics_fTilde_LB(:,:,:,:,:,:)
REAL,ALLOCATABLE                   :: Metrics_gTilde_LB(:,:,:,:,:,:)
REAL,ALLOCATABLE                   :: Metrics_hTilde_LB(:,:,:,:,:,:)
REAL,ALLOCATABLE                   :: sJ_LB            (  :,:,:,:,:)
REAL,ALLOCATABLE                   :: DetJac_N_LB(      :,:,:,:,:)
REAL,ALLOCATABLE                   :: DetJac_Ref_LB(    :,:,:,:,:)
#if USE_PARTICLES
REAL,ALLOCATABLE                   :: XCL_NGeo_LB(      :,:,:,:,:)
REAL,ALLOCATABLE                   :: dXCL_NGeo_LB(   :,:,:,:,:,:)
#endif /*USE_PARTICLES*/
! REAL,ALLOCATABLE                   :: DetJac_Ref_LB    (:,:,:,:,:)
! Sides
! REAL,ALLOCATABLE                   :: Face_xGP_LB    (  :,:,:,:)
! REAL,ALLOCATABLE                   :: NormVec_LB     (  :,:,:,:)
! REAL,ALLOCATABLE                   :: TangVec1_LB    (  :,:,:,:)
! REAL,ALLOCATABLE                   :: TangVec2_LB    (  :,:,:,:)
! REAL,ALLOCATABLE                   :: SurfElem_LB    (    :,:,:)
! REAL,ALLOCATABLE                   ::      Ja_Face_LB(:,:,:,:,:)
! Custom data type
INTEGER                            :: MPI_LENGTH(1),MPI_TYPE(1),MPI_STRUCT
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
! Timer
REAL                               :: StartT,EndT,WallTime
#if FV_ENABLED
INTEGER                            :: iFV
#else
INTEGER,PARAMETER                  :: iFV = 0
#endif
! !===================================================================================================================================
!
IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' | Shift mesh metrics during loadbalance...'
  GETTIME(StartT)

  ! volume data
  ALLOCATE(       XCL_N_LB(  3,0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate XCL_N over MPI
    MPI_LENGTH       = 3*(PP_N+1)*(PP_N+1)*(PP_NZ+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(XCL_N,counts_send,disp_send,MPI_STRUCT,XCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(XCL_N)
  CALL MOVE_ALLOC(XCL_N_LB,XCL_N)

  ALLOCATE(      dXCL_N_LB(3,3,0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate dXCL_N over MPI
    MPI_LENGTH       = 3*3*(PP_N+1)*(PP_N+1)*(PP_NZ+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(dXCL_N,counts_send,disp_send,MPI_STRUCT,dXCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(dXCL_N)
  CALL MOVE_ALLOC(dXCL_N_LB,dXCL_N)

#if USE_PARTICLES
  ALLOCATE(XCL_NGeo_LB(1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate JaCL_N over MPI
    MPI_LENGTH       = 3*(NGeo+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(XCL_NGeo,counts_send,disp_send,MPI_STRUCT,XCL_NGeo_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_WORLD,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(XCL_NGeo)
  CALL MOVE_ALLOC(XCL_NGeo_LB,XCL_NGeo)

  ALLOCATE(dXCL_NGeo_LB(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate JaCL_N over MPI
    MPI_LENGTH       = 3*3*(NGeo+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(dXCL_NGeo,counts_send,disp_send,MPI_STRUCT,dXCL_NGeo_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_WORLD,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(dXCL_NGeo)
  CALL MOVE_ALLOC(dXCL_NGeo_LB,dXCL_NGeo)
#endif /*USE_PARTICLES*/

  ALLOCATE(      JaCL_N_LB(3,3,0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate JaCL_N over MPI
    MPI_LENGTH       = 3*3*(PP_N+1)*(PP_N+1)*(PP_NZ+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(JaCL_N,counts_send,disp_send,MPI_STRUCT,JaCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(JaCL_N)
  CALL MOVE_ALLOC(JaCL_N_LB,JaCL_N)

  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate metrics over MPI
    MPI_LENGTH       = 3*(PP_N+1)*(PP_N+1)*(PP_NZ+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    ALLOCATE(Metrics_fTilde_LB(3,0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems, 0:FV_SIZE))
#if FV_ENABLED
    DO iFV = 0,FV_SIZE
#endif
    CALL MPI_ALLTOALLV(Metrics_fTilde(:,:,:,:,:,iFV),counts_send,disp_send,MPI_STRUCT,Metrics_fTilde_LB(:,:,:,:,:,iFV),counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
#if FV_ENABLED
    END DO
#endif
    DEALLOCATE(Metrics_fTilde)
    CALL MOVE_ALLOC(Metrics_fTilde_LB,Metrics_fTilde)

    ALLOCATE(Metrics_gTilde_LB(3,0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems, 0:FV_SIZE))
#if FV_ENABLED
    DO iFV = 0,FV_SIZE
#endif
    CALL MPI_ALLTOALLV(Metrics_gTilde(:,:,:,:,:,iFV),counts_send,disp_send,MPI_STRUCT,Metrics_gTilde_LB(:,:,:,:,:,iFV),counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
#if FV_ENABLED
    END DO
#endif
    DEALLOCATE(Metrics_gTilde)
    CALL MOVE_ALLOC(Metrics_gTilde_LB,Metrics_gTilde)

    ALLOCATE(Metrics_hTilde_LB(3,0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems, 0:FV_SIZE))
#if FV_ENABLED
    DO iFV = 0,FV_SIZE
#endif
    CALL MPI_ALLTOALLV(Metrics_hTilde(:,:,:,:,:,iFV),counts_send,disp_send,MPI_STRUCT,Metrics_hTilde_LB(:,:,:,:,:,iFV),counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
#if FV_ENABLED
    END DO
#endif
    DEALLOCATE(Metrics_hTilde)
    CALL MOVE_ALLOC(Metrics_hTilde_LB,Metrics_hTilde)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE

  ALLOCATE(sJ_LB            (    0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems, 0:FV_SIZE))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate scaled Jacobians over MPI
    MPI_LENGTH       = (PP_N+1)*(PP_N+1)*(PP_NZ+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

#if FV_ENABLED
    DO iFV = 0,FV_SIZE
#endif
    CALL MPI_ALLTOALLV(sJ(:,:,:,:,0),counts_send,disp_send,MPI_STRUCT,sJ_LB(:,:,:,:,0),counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
#if FV_ENABLED
    END DO
#endif
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(sJ)
  CALL MOVE_ALLOC(sJ_LB,sJ)

  ALLOCATE(DetJac_N_LB            (1,  0:PP_N   ,0:PP_N   ,0:PP_NZ  ,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate scaled Jacobians over MPI
    MPI_LENGTH       = (PP_N+1)*(PP_N+1)*(PP_NZ+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(DetJac_N,counts_send,disp_send,MPI_STRUCT,DetJac_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(DetJac_N)
  CALL MOVE_ALLOC(DetJac_N_LB,DetJac_N)

  ALLOCATE(DetJac_Ref_LB            (1,  0:NGeoRef   ,0:NGeoRef   ,0:ZDIM(NGeoRef)  ,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate scaled Jacobians over MPI
    MPI_LENGTH       = (NGeoRef+1)*(NGeoRef+1)*(ZDIM(NGeoRef)+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(DetJac_Ref,counts_send,disp_send,MPI_STRUCT,DetJac_Ref_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(DetJac_Ref)
  CALL MOVE_ALLOC(DetJac_Ref_LB,DetJac_Ref)

  ! ! surface data
  ! ASSOCIATE (&
  !         counts_send  => INT(MPInSideSend     ) ,&
  !         disp_send    => INT(MPIoffsetSideSend) ,&
  !         counts_recv  => INT(MPInSideRecv     ) ,&
  !         disp_recv    => INT(MPIoffsetSideRecv))
  !   ! Communicate scaled Jacobians over MPI
  !   IPWRITE(*,*) 'MPInSideSend,MPIoffsetSideSend:', MPInSideSend,MPIoffsetSideSend
  !   IPWRITE(*,*) 'MPInSideRecv,MPIoffsetSideRecv:', MPInSideSend,MPIoffsetSideSend
  !   MPI_LENGTH       = 3*(PP_N+1)**2
  !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
  !   MPI_TYPE         = MPI_DOUBLE_PRECISION
  !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
  !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
  !
  !   ! Communicate side metrics over MPI
  !   ALLOCATE(Face_xGP_LB      (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(Face_xGP,counts_send,disp_send,MPI_STRUCT,Face_xGP_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
  !   DEALLOCATE(Face_xGP)
  !   CALL MOVE_ALLOC(Face_xGP_LB,Face_xGP)
  !
  !   ALLOCATE(NormVec_LB       (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(NormVec,counts_send,disp_send,MPI_STRUCT,NormVec_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
  !   DEALLOCATE(NormVec)
  !   CALL MOVE_ALLOC(NormVec_LB,NormVec)
  !
  !   ALLOCATE(TangVec1_LB      (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(TangVec1,counts_send,disp_send,MPI_STRUCT,TangVec1_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
  !   DEALLOCATE(TangVec1)
  !   CALL MOVE_ALLOC(TangVec1_LB,TangVec1)
  !
  !   ALLOCATE(TangVec2_LB      (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(TangVec2,counts_send,disp_send,MPI_STRUCT,TangVec2_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
  !   DEALLOCATE(TangVec2)
  !   CALL MOVE_ALLOC(TangVec2_LB,TangVec2)
  ! END ASSOCIATE
  !
  ! ASSOCIATE (&
  !         counts_send  => INT(MPInSideSend     ) ,&
  !         disp_send    => INT(MPIoffsetSideSend) ,&
  !         counts_recv  => INT(MPInSideRecv     ) ,&
  !         disp_recv    => INT(MPIoffsetSideRecv))
  !   ! Communicate scaled Jacobians over MPI
  !   MPI_LENGTH       = (PP_N+1)**2
  !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
  !   MPI_TYPE         = MPI_DOUBLE_PRECISION
  !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
  !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
  !   ALLOCATE(SurfElem_LB      (  0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(SurfElem,counts_send,disp_send,MPI_STRUCT,SurfElem_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
  !   DEALLOCATE(SurfElem)
  !   CALL MOVE_ALLOC(SurfElem_LB,SurfElem)
  ! END ASSOCIATE
  !
  ! ASSOCIATE (&
  !         counts_send  => INT(MPInSideSend     ) ,&
  !         disp_send    => INT(MPIoffsetSideSend) ,&
  !         counts_recv  => INT(MPInSideRecv     ) ,&
  !         disp_recv    => INT(MPIoffsetSideRecv))
  !   ! Communicate scaled Jacobians over MPI
  !   MPI_LENGTH       = 3*3*(PP_N+1)**2
  !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
  !   MPI_TYPE         = MPI_DOUBLE_PRECISION
  !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
  !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
  !   ALLOCATE(     Ja_Face_LB(3,3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(Ja_Face,counts_send,disp_send,MPI_STRUCT,Ja_Face_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_FLEXI,iError)
  !   DEALLOCATE(Ja_Face)
  !   CALL MOVE_ALLOC(Ja_Face_LB,Ja_Face)
  ! END ASSOCIATE

  GETTIME(EndT)
  WallTime = EndT-StartT
  CALL DisplayMessageAndTime(WallTime,'DONE!',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
END IF

END SUBROUTINE MoveMetrics

END MODULE MOD_LoadBalance_Metrics
