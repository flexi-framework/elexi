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

!==================================================================================================================================
!> Contains control routines to read high-order meshes, provide mesh data to the solver, build the metrics, partition the domain.
!==================================================================================================================================
MODULE MOD_Mesh_Check
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CheckMesh
  MODULE PROCEDURE CheckMesh
END INTERFACE

PUBLIC::CheckMesh
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Sanity checks the mesh connectivity by comparing the unique node IDs
!==================================================================================================================================
SUBROUTINE CheckMesh()
! MODULES
USE MOD_Globals
USE MOD_PreProc
! USE MOD_HDF5_Input          ,ONLY: File_ID,OpenDataFile,CloseDataFile,ReadAttribute
! USE MOD_Mesh_Vars           ,ONLY: UseCurveds,MeshFile
USE MOD_Mesh_Vars           ,ONLY: firstInnerSide,lastInnerSide,lastMPISide_MINE
USE MOD_Mesh_Vars           ,ONLY: nSides,nBCSides,nMPISides!,nMPIPeriodics
USE MOD_Mesh_Vars           ,ONLY: SideToElem,S2V2
USE MOD_Mesh_Vars           ,ONLY: BoundaryType
USE MOD_Mesh_Vars           ,ONLY: AnalyzeSide,meshHasMortars
! USE MOD_Mesh_Vars           ,ONLY: Elem_xGP
USE MOD_Mesh_Vars           ,ONLY: NodeCoords,NGeo
USE MOD_ReadInTools         ,ONLY: GETLOGICAL
#if USE_MPI
USE MOD_Mesh_Vars           ,ONLY: nGlobalSides,nGlobalBCSides,nGlobalPeriodicSides
USE MOD_Mesh_Vars           ,ONLY: firstMPISide_MINE,firstMPISide_YOUR
USE MOD_MPI_Vars            ,ONLY: MPIRequest_U,nNbProcs
USE MOD_MPI                 ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#else
USE MOD_Mesh_Vars           ,ONLY: nBCSides,nPeriodicSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars    ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER               :: SIDE_BCID=5
! REAL,PARAMETER                  :: RelTol=1.E-9
LOGICAL                         :: doCheckMesh!,doCheckRel
INTEGER                         :: NGeo_HDF5
INTEGER                         :: p,q
INTEGER                         :: ElemID,nbElemID,BCID
INTEGER                         :: SideID,firstSideID,lastSideID,locSide,nbLocSide,flip
INTEGER                         :: nCheckSides,nCheckedSides,nFailedAbsSides,nFailedTolSides,nSkippedPeriodicSides
! REAL                            :: Side_xGP(       3,0:PP_N,0:PP_NZ)
! REAL                            :: Side_xGP_master(3,0:PP_N,0:PP_NZ,1:nSides)
! REAL                            :: Side_xGP_slave( 3,0:PP_N,0:PP_NZ,1:nSides)
REAL                            :: Side_NodeCoords(       3,0:NGeo,0:ZDIM(NGeo))
REAL                            :: Side_NodeCoords_master(3,0:NGeo,0:ZDIM(NGeo),1:nSides)
REAL                            :: Side_NodeCoords_slave( 3,0:NGeo,0:ZDIM(NGeo),1:nSides)
#if USE_MPI
INTEGER                         :: DataSizeSide
#endif /*USE_MPI*/
! Timers
REAL                            :: StartT,EndT
!==================================================================================================================================

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

LBWRITE(UNIT_stdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)')' CHECK SIDE CONNECTIVITY...'

doCheckMesh = GETLOGICAL('meshCheckConnectivity')
IF (.NOT.doCheckMesh .OR. meshHasMortars) THEN
  IF (meshHasMortars) THEN
    LBWRITE(UNIT_stdOut,'(A)')' Mortar elements encountered. Force-disabling check side connectivity!'
  END IF
  LBWRITE(UNIT_stdOut,'(A)')' CHECK SIDE CONNECTIVITY DONE'
  LBWRITE(UNIT_stdOut,'(132("-"))')
  RETURN
END IF
GETTIME(StartT)

! doCheckRel = .FALSE.
! IF (.NOT.UseCurveds) THEN
!   CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
!   CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo_HDF5)
!   CALL CloseDataFile()

!   IF (NGeo_HDF5.GT.1) THEN
!     CALL PrintWarning('NGeo_HDF5>1 but UseCurveds=F. Mesh connectivity is underresolved and might not be accurate!')
!     doCheckRel = .TRUE.
!   END IF
! END IF

! Nullify arrays
! Side_xGP        = 0
! Side_xGP_master = 0
! Side_xGP_slave  = 0
Side_NodeCoords        = 0
Side_NodeCoords_master = 0
Side_NodeCoords_slave  = 0

#if USE_MPI
! DataSizeSide = 3*(PP_N+1)*(PP_NZ+1)
DataSizeSide = 3*(NGeo+1)*(ZDIM(NGeo)+1)

! Fill the arrays for all slave MPI sides
! CALL StartReceiveMPIData(Side_xGP_slave,DataSizeSide,1,nSides,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE / Side_xGP_slave: slave -> master
CALL StartReceiveMPIData(Side_NodeCoords_slave,DataSizeSide,1,nSides,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE / Side_xGP_slave: slave -> master

firstSideID = firstMPISide_YOUR
 lastSideID = nSides

DO SideID = firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID    ,SideID)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID ,SideID)

  ! master sides
  IF (ElemID.GT.0) THEN
    CALL Abort(__STAMP__,'Encountered master side in MPISide_YOUR!')
  END IF ! ElemID.GT.0

  ! slave side (ElemID,locSide and flip =-1 if not existing)
  IF (nbElemID.GT.0) THEN
    nbLocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    flip      = SideToElem(S2E_FLIP          ,SideID)

    SELECT CASE(nbLocSide)
      CASE(XI_MINUS)
        ! Side_xGP = Elem_xGP(:,0   ,:   ,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,0   ,:   ,:   ,nbElemID)
      CASE(ETA_MINUS)
        ! Side_xGP = Elem_xGP(:,:   ,0   ,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,0   ,:   ,nbElemID)
      CASE(ZETA_MINUS)
        ! Side_xGP = Elem_xGP(:,:   ,:   ,0   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,:   ,0   ,nbElemID)
      CASE(XI_PLUS)
        ! Side_xGP = Elem_xGP(:,PP_N,:   ,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,NGeo,:   ,:   ,nbElemID)
      CASE(ETA_PLUS)
        ! Side_xGP = Elem_xGP(:,:   ,PP_N,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,NGeo,:   ,nbElemID)
      CASE(ZETA_PLUS)
        ! Side_xGP = Elem_xGP(:,:   ,:   ,PP_N,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,:   ,NGeo,nbElemID)
    END SELECT

    ! DO q=0,ZDIM(PP_N); DO p=0,PP_N
      ! Side_xGP_slave( :,p,q,SideID) = Side_xGP(:,S2V2(1,p,q,flip,nblocSide),S2V2(2,p,q,flip,nblocSide))
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      Side_NodeCoords_slave( :,p,q,SideID) = Side_NodeCoords(:,S2V2(1,p,q,flip,nblocSide),S2V2(2,p,q,flip,nblocSide))
    END DO; END DO
  END IF
END DO

! TODO: Mortars
! CALL U_MortarCons(U_master,Side_xGP_slave,doMPISides=.TRUE.)
! CALL StartSendMPIData(   Side_xGP_slave,DataSizeSide,1,nSides,MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR / Side_xGP_slave: slave -> master
CALL StartSendMPIData(   Side_NodeCoords_slave,DataSizeSide,1,nSides,MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR / Side_xGP_slave: slave -> master
#endif /*USE_MPI*/

! Fill the array for all remaining sides
firstSideID = 1
 lastSideID =  lastMPISide_MINE

DO SideID = firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID    ,SideID)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID ,SideID)

  ! master sides
  IF (ElemID.GT.0) THEN
    locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
    flip    = 0
    ! ElemID_master(SideID) = ElemID

    SELECT CASE(locSide)
      CASE(XI_MINUS)
        ! Side_xGP = Elem_xGP(:,0   ,:   ,:   ,ElemID)
        Side_NodeCoords = NodeCoords(:,0   ,:   ,:   ,ElemID)
      CASE(ETA_MINUS)
        ! Side_xGP = Elem_xGP(:,:   ,0   ,:   ,ElemID)
        Side_NodeCoords = NodeCoords(:,:   ,0   ,:   ,ElemID)
      CASE(ZETA_MINUS)
        ! Side_xGP = Elem_xGP(:,:   ,:   ,0   ,ElemID)
        Side_NodeCoords = NodeCoords(:,:   ,:   ,0   ,ElemID)
      CASE(XI_PLUS)
        ! Side_xGP = Elem_xGP(:,PP_N,:   ,:   ,ElemID)
        Side_NodeCoords = NodeCoords(:,NGeo,:   ,:   ,ElemID)
      CASE(ETA_PLUS)
        ! Side_xGP = Elem_xGP(:,:   ,PP_N,:   ,ElemID)
        Side_NodeCoords = NodeCoords(:,:   ,NGeo,:   ,ElemID)
      CASE(ZETA_PLUS)
        ! Side_xGP = Elem_xGP(:,:   ,:   ,PP_N,ElemID)
        Side_NodeCoords = NodeCoords(:,:   ,:   ,NGeo,ElemID)
    END SELECT

    ! DO q=0,ZDIM(PP_N); DO p=0,PP_N
    !   Side_xGP_master(:,p,q,SideID) = Side_xGP(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      Side_NodeCoords_master(:,p,q,SideID) = Side_NodeCoords(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
    END DO; END DO
  END IF ! ElemID.GT.0

  ! slave side (ElemID,locSide and flip =-1 if not existing)
  IF (nbElemID.GT.0) THEN
    nbLocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    flip      = SideToElem(S2E_FLIP          ,SideID)
    ! ElemID_slave(SideID) = ElemID

    SELECT CASE(nbLocSide)
      CASE(XI_MINUS)
        ! Side_xGP = Elem_xGP(:,0   ,:   ,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,0   ,:   ,:   ,nbElemID)
      CASE(ETA_MINUS)
        ! Side_xGP = Elem_xGP(:,:   ,0   ,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,0   ,:   ,nbElemID)
      CASE(ZETA_MINUS)
        ! Side_xGP = Elem_xGP(:,:   ,:   ,0   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,:   ,0   ,nbElemID)
      CASE(XI_PLUS)
        ! Side_xGP = Elem_xGP(:,PP_N,:   ,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,NGeo,:   ,:   ,nbElemID)
      CASE(ETA_PLUS)
        ! Side_xGP = Elem_xGP(:,:   ,PP_N,:   ,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,NGeo,:   ,nbElemID)
      CASE(ZETA_PLUS)
        ! Side_xGP = Elem_xGP(:,:   ,:   ,PP_N,nbElemID)
        Side_NodeCoords = NodeCoords(:,:   ,:   ,NGeo,nbElemID)
    END SELECT

    ! DO q=0,ZDIM(PP_N); DO p=0,PP_N
    !   Side_xGP_slave( :,p,q,SideID) = Side_xGP(:,S2V2(1,p,q,flip,nbLocSide),S2V2(2,p,q,flip,nbLocSide))
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      Side_NodeCoords_slave( :,p,q,SideID) = Side_NodeCoords(:,S2V2(1,p,q,flip,nbLocSide),S2V2(2,p,q,flip,nbLocSide))
    END DO; END DO
  END IF
END DO

! TODO: Mortars
! CALL U_MortarCons(U_master,U_slave,doMPISides=.FALSE.)

#if USE_MPI
! Complete send / receive of side data
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! Side_xGP_slave: slave -> master
#endif /*USE_MPI*/

nCheckSides           = nSides-nMPISides/2-nBCSides!-nPeriodicSides+nMPIPeriodics/2
nCheckedSides         = 0
nFailedAbsSides       = 0
nFailedTolSides       = 0
nSkippedPeriodicSides = 0

#if USE_MPI
DO SideID = firstMPISide_MINE,lastMPISide_MINE
  nCheckedSides = nCheckedSides + 1

  ! Periodic sides are included in analyzeSide
  BCID = AnalyzeSide(SideID)
  IF (BCID.GT.0) THEN
    ! Ignore periodic sides as they cannot match
    IF (BoundaryType(BCID,BC_TYPE).EQ.1) THEN
      nSkippedPeriodicSides = nSkippedPeriodicSides + 1
      CYCLE
    END IF
  END IF

  ! IF (ANY(Side_xGP_master(:,:,:,SideID).NE.Side_xGP_slave(:,:,:,SideID))) THEN
  IF (ANY(Side_NodeCoords_master(:,:,:,SideID).NE.Side_NodeCoords_slave(:,:,:,SideID))) THEN
    nFailedAbsSides = nFailedAbsSides + 1

    ! IF (.NOT.doCheckRel) CYCLE
    ! IF (.NOT.ALL(ALMOSTEQUALABSANDREL(Side_xGP_master(:,:,:,SideID),Side_xGP_slave(:,:,:,SideID),RelTol))) THEN
      ! nFailedTolSides = nFailedTolSides + 1
    ! END IF
  END IF
END DO ! SideID = 1,nSides
#endif /*USE_MPI*/

DO SideID = firstInnerSide,lastInnerSide
  nCheckedSides = nCheckedSides + 1

  ! Periodic sides are included in analyzeSide
  BCID = AnalyzeSide(SideID)
  IF (BCID.GT.0) THEN
    ! Ignore periodic sides as they cannot match
    IF (BoundaryType(BCID,BC_TYPE).EQ.1) THEN
      nSkippedPeriodicSides = nSkippedPeriodicSides + 1
      CYCLE
    END IF
  END IF

  ! IF (ANY(Side_xGP_master(:,:,:,SideID).NE.Side_xGP_slave(:,:,:,SideID))) THEN
  IF (ANY(Side_NodeCoords_master(:,:,:,SideID).NE.Side_NodeCoords_slave(:,:,:,SideID))) THEN
    nFailedAbsSides = nFailedAbsSides + 1

    ! IF (.NOT.doCheckRel) CYCLE
    ! IF (.NOT.ALL(ALMOSTEQUALABSANDREL(Side_xGP_master(:,:,:,SideID),Side_xGP_slave(:,:,:,SideID),RelTol))) THEN
    !   nFailedTolSides = nFailedTolSides + 1
    ! END IF
  END IF
END DO ! SideID = 1,nSides

#if USE_MPI
IF (MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nCheckSides          ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nCheckedSides        ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nFailedAbsSides      ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  ! CALL MPI_REDUCE(MPI_IN_PLACE,nFailedTolSides      ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nSkippedPeriodicSides,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
ELSE
  CALL MPI_REDUCE(nCheckSides          ,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  CALL MPI_REDUCE(nCheckedSides        ,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  CALL MPI_REDUCE(nFailedAbsSides      ,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  ! CALL MPI_REDUCE(nFailedTolSides      ,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
  CALL MPI_REDUCE(nSkippedPeriodicSides,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iERROR)
END IF
#endif /*USE_MPI*/

#if !USE_MPI
ASSOCIATE(nGlobalSides         => nSides         ,&
          nGlobalPeriodicSides => nPeriodicSides ,&
          nGlobalBCSides       => nBCSides)
#endif
nCheckSides     = nGlobalSides-nGlobalBCSides!-nPeriodicSides+nMPIPeriodics/2

IF (MPIRoot .AND. nGlobalPeriodicSides.NE.nSkippedPeriodicSides) &
  CALL Abort(__STAMP__,'Not all periodic sides could be skipped!')

LBWRITE(UNIT_stdOut,'(A,A34,I0,A,I0)')' |','nSides (requested/checked) | ',nCheckSides,'/',nCheckedSides
LBWRITE(UNIT_stdOut,'(A,A34,I0,A,I0)')' |','nSides (periodic /skipped) | ',nSkippedPeriodicSides,'/',nGlobalPeriodicSides
LBWRITE(UNIT_stdOut,'(A,A34,I0)'     )' |','nSides (failed,abs.)       | ',nFailedAbsSides
! IF (doCheckRel) THEN
! LBWRITE(UNIT_stdOut,'(A,A34,I0)'     )' |','nSides (failed,tol.)       | ',nFailedTolSides
! END IF

#if !USE_MPI
END ASSOCIATE
#endif

! Abort if we encounter incorrectly connected sides
IF (MPIRoot) THEN
  IF (nFailedAbsSides.GT.0) THEN
    ! IF (.NOT.doCheckRel) CALL Abort(__STAMP__,'Aborting due to incorrectly connected sides')
                         CALL Abort(__STAMP__,'Aborting due to incorrectly connected sides')
  END IF

  ! IF (nFailedTolSides.GT.0) THEN
                         ! CALL Abort(__STAMP__,'Aborting due to incorrectly connected sides')
  ! END IF
END IF

GETTIME(EndT)
SWRITE(UNIT_stdOut,'(A,F0.3,A)')' CHECK SIDE CONNECTIVITY DONE! [',EndT-StartT,'s]'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE CheckMesh
END MODULE MOD_Mesh_Check
