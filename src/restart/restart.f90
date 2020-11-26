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
!> \brief Routines that handle restart capabilities.
!>
!> With this feature a simulation can be resumed from a state file that has been created during a previous
!> simulation (restart file). The restart file is passed to FLEXI as a second command line argument.
!> The restart can also be performed from a file with a different polynomial degree or node type than the current simulation.
!==================================================================================================================================
MODULE MOD_Restart
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitRestart
  MODULE PROCEDURE InitRestart
END INTERFACE

INTERFACE Restart
  MODULE PROCEDURE Restart
END INTERFACE

INTERFACE FinalizeRestart
  MODULE PROCEDURE FinalizeRestart
END INTERFACE

PUBLIC :: DefineParametersRestart
PUBLIC :: InitRestart,FinalizeRestart
PUBLIC :: Restart
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE DefineParametersRestart()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Restart")
CALL prms%CreateLogicalOption('ResetTime',      "Override solution time to t=0 on restart.", '.FALSE.')
CALL prms%CreateLogicalOption('RestartMean',    "Flag to denote restart from time-averaged file.", '.FALSE.')
CALL prms%CreateLogicalOption('RestartTurb',    "Flag to denote restart from time-averaged filecontaining turbulent quantities.", '.FALSE.')
#if EQNSYSNR == 3
CALL prms%CreateRealOption(   'RestartMuTilda', "Constant mu_tilda to be applied throughout the domain.", '0.')
#endif
#if FV_ENABLED
CALL prms%CreateIntOption(    'NFVRestartSuper',"Polynomial degree for equidistant supersampling of FV subcells when restarting&
                                                 &on a different polynomial degree. Default 2*MAX(N,NRestart).")
#endif
END SUBROUTINE DefineParametersRestart

!==================================================================================================================================
!> \brief Initialize all necessary information to perform the restart.
!>
!> The routine checks if two arguments have been passed to FLEXI on the command line. If so, the second one is supposed
!> to be the restart state. If only one argument has been passed, no restart will be performed.
!> - In the restart case, it is checked if the restart file exists at all. If so, the properties of the restart file will be
!>   read (polynomial degree and node type are needed) and stored to use later.The flag DoRestart indicating the restart is set
!>   to be used by other routines. Also the simulation time of the restart is read.
!>   A optional parameter ResetTime can be used to set the restart time to 0.
!> - If no restart is performed, the RestartTime is set to 0.
!>
!> The routine also checks if the node type and polynomial degree of the restart file is the same than in the current simulation.
!> If not, a flag InterpolateSolution is set. This will be used by the actual Restart routine.
!==================================================================================================================================
SUBROUTINE InitRestart(RestartFile_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Restart_Vars
#if FV_ENABLED
USE MOD_StringTools,        ONLY: INTTOSTR
USE MOD_ReadInTools,        ONLY: GETINT
#endif
USE MOD_HDF5_Input,         ONLY: ISVALIDHDF5FILE
USE MOD_Interpolation_Vars, ONLY: InterpolationInitIsDone,NodeType
USE MOD_HDF5_Input,         ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_ReadInTools,        ONLY: GETLOGICAL,GETREAL,GETREALARRAY
USE MOD_Mesh_Vars,          ONLY: nGlobalElems,NGeo
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: RestartFile_in !< state file to restart from
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: ResetTime,validHDF5
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.RestartInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,'InitRestart not ready to be called or already called.')
END IF

IF (PRESENT(RestartFile_in)) RestartFile = RestartFile_in

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'

! Check if we want to perform a restart
IF (LEN_TRIM(RestartFile).GT.0) THEN
  SWRITE(UNIT_StdOut,'(A,A,A)')' | Restarting from file "',TRIM(RestartFile),'":'
  ! Check if restart file is a valid state
  validHDF5 = ISVALIDHDF5FILE(RestartFile)
  IF(.NOT.validHDF5) &
      CALL CollectiveStop(__STAMP__,'ERROR - Restart file not a valid state file.')
  ! Set flag indicating a restart to other routines
  DoRestart = .TRUE.
  ! Read in parameters of restart solution
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  ! Read this here because we need to change array names in case of a mean file
  RestartTurb=GETLOGICAL('RestartTurb','.FALSE.')
  IF (.NOT.RestartTurb.OR.postiMode) RestartMean=GETLOGICAL('RestartMean','.FALSE.')
  ! Read in attributes
  IF (.NOT.RestartMean.OR.postiMode) THEN
    CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart)
  ELSE
    CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart,'Mean')
  END IF
  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  ! Option to set the calculation time to 0 even tho performing a restart
  ResetTime=GETLOGICAL('ResetTime','.FALSE.')
#if EQNSYSNR == 3
  MuTilda  =GETREAL   ('RestartMuTilda','0.')
#endif
  IF(ResetTime) RestartTime=0.
  CALL CloseDataFile()

  IF (nElems_Restart.NE.nGlobalElems) THEN
    CALL CollectiveStop(__STAMP__, "Restart File has different number of elements!")
  END IF
  IF (nVar_Restart.GT.PP_nVar) THEN
    IF (RestartMean.AND..NOT.postiMode) THEN
#if EQNSYSNR == 3
      SWRITE(UNIT_StdOut,'(A)') ' Restart from mean file. Parameter for mu_tilda required!'
#endif
    ELSEIF(RestartTurb.AND.(nVar_Restart.NE.PP_nVar+2)) THEN
      CALL CollectiveStop(__STAMP__, "Restart File has wrong number of variables for restart with turbulent quantities!")
    ELSE
      SWRITE(UNIT_StdOut,'(A)') ' Restart file has more variables than current equation system, will be truncated!'
    END IF
  END IF
  IF (nVar_Restart.LT.PP_nVar) THEN
#if EQNSYSNR == 3
    SWRITE(UNIT_StdOut,'(A)') ' Restart file has less variables than current equation system. Parameter for mu_tilda required!'
#endif
  END IF
ELSE
  ! No restart
  RestartTime = 0.
  SWRITE(UNIT_StdOut,'(A)')' | No restart wanted, doing a fresh computation!'
END IF

! Check if we need to interpolate the restart file to our current polynomial degree and node type
IF(DoRestart .AND. ((N_Restart.NE.PP_N) .OR. (TRIM(NodeType_Restart).NE.TRIM(NodeType))))THEN
  InterpolateSolution=.TRUE.
  IF(MIN(N_Restart,PP_N).LT.NGeo) &
    CALL PrintWarning('The geometry is or was underresolved and will potentially change on restart!')
#if FV_ENABLED
  NFVRestartSuper = GETINT('NFVRestartSuper',INTTOSTR(2*MAX(PP_N,N_Restart)))
#endif
ELSE
  InterpolateSolution=.FALSE.
END IF

RestartInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RESTART DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRestart


!==================================================================================================================================
!> \brief This routine performs the actual restart. It is called from the FLEXI main program just before the TimeDisc routine.
!>
!> For the restart there are 2 cases, depending on the flag InterpolateSolution set by InitRestart:
!> - No interpolation is necessary. Then simply read in the DG_Solution array from the restart file and store it in the solution
!>   array U.
!> - We need to interpolate the restart solution. If the polynomial degree of our computation is lower than in the restart file,
!>   a simple change basis is used to get the current solution U. If the polynomial degree is higher than in the restart file,
!>   special care is taken to ensure a conservative projection of the restart solution. To do this, the restart solution is
!>   transformed to reference space using the Jacobian built on a polynomial degree of 3*NGeo (so it can be represented exactly),
!>   then the change basis is applied. The resulting solution U is then transformed back to physical space.
!>
!> All state files that would be re-written by the simulation (with the same project name and a time stamp after the restart time
!> or all files with the same project name if no restart is performed) are deleted at the end of the routine.
!==================================================================================================================================
SUBROUTINE Restart(doFlushFiles)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Restart_Vars
USE MOD_DG_Vars,            ONLY: U
USE MOD_Mesh_Vars,          ONLY: offsetElem,detJac_Ref,Ngeo
USE MOD_Mesh_Vars,          ONLY: nElems,nGlobalElems
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_HDF5_Input,         ONLY: OpenDataFile,CloseDataFile,ReadArray,GetArrayAndName,GetVarNames
USE MOD_HDF5_Output,        ONLY: FlushFiles
USE MOD_Interpolation,      ONLY: GetVandermonde
USE MOD_ApplyJacobianCons,  ONLY: ApplyJacobianCons
USE MOD_Interpolation_Vars, ONLY: NodeType
#if FV_ENABLED
USE MOD_FV,                 ONLY: FV_ProlongFVElemsToFace
USE MOD_FV_Vars,            ONLY: FV_Elems
USE MOD_Indicator_Vars,     ONLY: IndValue
USE MOD_StringTools,        ONLY: STRICMP
#endif
#if PP_dim == 3
USE MOD_2D,                 ONLY: ExpandArrayTo3D
#else
USE MOD_2D,                 ONLY: to2D_rank5
#endif
USE MOD_IO_HDF5
USE MOD_HDF5_Input,         ONLY: GetDataSize
#if USE_RW
USE MOD_DG_Vars,            ONLY: Uturb
USE MOD_Equation_Vars,      ONLY: nVarTurb
#endif
USE MOD_EOS,                ONLY: PrimToCons
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: doFlushFiles !< flag to delete old state files
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE   :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE   :: U_localNVar(:,:,:,:,:)
#if PP_dim == 3
REAL,ALLOCATABLE   :: U_local2(:,:,:,:,:)
#endif
INTEGER            :: iElem,i,j,k
INTEGER            :: HSize_proc(5)
REAL,ALLOCATABLE   :: JNR(:,:,:,:)
REAL               :: Vdm_NRestart_N(0:PP_N,0:N_Restart)
REAL               :: Vdm_3Ngeo_NRestart(0:N_Restart,0:3*NGeo)
LOGICAL            :: doFlushFiles_loc
#if FV_ENABLED
INTEGER            :: nVal(15),iVar
REAL,ALLOCATABLE   :: ElemData(:,:),tmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesElemData(:)
#endif
#if GCL
REAL,ALLOCATABLE   :: Jac_local(:,:,:,:,:)
#if PP_dim == 3
REAL,ALLOCATABLE   :: Jac_local2(:,:,:,:,:)
#endif
LOGICAL            :: foundJac
INTEGER            :: HSize_procJac(5)
#endif
LOGICAL            :: VarNamesExist
CHARACTER(LEN=255),ALLOCATABLE  :: varnames_tmp(:)
!==================================================================================================================================
IF (PRESENT(doFlushFiles)) THEN
  doFlushFiles_loc = doFlushFiles
ELSE
  doFlushFiles_loc = .TRUE.
END IF

IF(DoRestart)THEN
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' PERFORMING RESTART...'

  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#if FV_ENABLED
  ! Read FV element distribution and indicator values from elem data array if possible
  CALL GetArrayAndName('ElemData','VarNamesAdd',nVal,tmp,VarNamesElemData)
  ALLOCATE(ElemData(nVal(1),nVal(2)))
  ElemData = RESHAPE(tmp,(/nVal(1),nVal(2)/))
  ! search for FV_Elems and IndValue
  FV_Elems=0
  IndValue=0.
  DO iVar=1,nVal(1)
    IF (STRICMP(VarNamesElemData(iVar),"FV_Elems")) THEN
      FV_Elems = INT(ElemData(iVar,:))
    END IF
    IF (STRICMP(VarNamesElemData(iVar),"IndValue")) THEN
      IndValue = ElemData(iVar,:)
    END IF
  END DO
  DEALLOCATE(ElemData,VarNamesElemData,tmp)
  CALL FV_ProlongFVElemsToFace()
#endif

  ! Mean files only have a dummy DG_Solution, we have to pick the "Mean" array in this case
  IF (RestartMean.AND..NOT.postiMode) THEN
    SWRITE(UNIT_stdOut,'(A)') 'Trying to restart from a mean file...'
    SWRITE(UNIT_StdOut,'(A)') 'WARNING: Make absolutely sure the mean file has ALL conservative or primitve variables available!'
    SWRITE(UNIT_StdOut,'(132("-"))')
    CALL GetDataSize(File_ID,'Mean',nDims,HSize)
  ELSE
    CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
  END IF
! Sanity check
!  IF ((HSize(1).LT.PP_nVar).OR.(HSize(2).NE.N_Restart+1).OR.(HSize(3).NE.N_Restart+1).OR.(HSize(5).NE.nGlobalElems)) THEN
!    CALL Abort(__STAMP__,"Dimensions of restart file do not match!")
!  END IF
  IF ((HSize(2).NE.N_Restart+1).OR.(HSize(3).NE.N_Restart+1).OR.(HSize(5).NE.nGlobalElems)) THEN
    CALL Abort(__STAMP__,"Dimensions of restart file do not match!")
  END IF
! Allow restart from LES to RANS solution
  IF ((HSize(1).LT.PP_nVar).OR.(HSize(2).NE.N_Restart+1).OR.(HSize(3).NE.N_Restart+1).OR.(HSize(5).NE.nGlobalElems)) THEN
#if EQNSYSNR == 3
    IF (MuTilda.NE.0) THEN
        SWRITE(UNIT_StdOut,'(A)')'Dimensions of restart file do not match! Proceeding because I assume you want a RANS restart.'
    ELSE
        CALL Abort(__STAMP__,"Dimensions of restart file do not match! If RANS restart wanted, provide RestartMuTilda.")
    END IF
#endif
  END IF
  HSize_proc = INT(HSize)
  HSize_proc(5) = nElems
  ! Allocate larger array for restart with turbulent quantitites
  ALLOCATE(U_local(nVar_Restart,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
  ! Mean files only have a dummy DG_Solution, we have to pick the "Mean" array in this case
  IF (RestartMean.AND..NOT.postiMode) THEN
    VarNamesExist=.FALSE.
    SDEALLOCATE(varnames_tmp)
    CALL GetVarNames("VarNames_Mean",varnames_tmp,VarNamesExist)
    IF((TRIM(varnames_tmp(2)).NE.'VelocityX').AND.(TRIM(varnames_tmp(2)).NE.'MomentumX'))THEN
      CALL Abort(__STAMP__,"Mean file has not all primitve and conservative variables.")
    ELSE
      CALL ReadArray('Mean',5,HSize_proc,OffsetElem,5,RealArray=U_local)
    END IF
  ELSE
    CALL ReadArray('DG_Solution',5,HSize_proc,OffsetElem,5,RealArray=U_local)
  END IF

  ! Truncate the solution if we read a restart file from a different equation system or from a time-averaged file
  IF (PP_nVar.LT.nVar_Restart) THEN
    ALLOCATE(U_localNVar(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
!    U_localNVar(1:PP_nVar,:,:,:,:) = U_local(1:PP_nVar,:,:,:,:)
    IF(TRIM(varnames_tmp(2)).NE.'MomentumX')THEN
      CALL PrimToCons(HSize_proc(2)-1,U_local(1:PP_nVarPrim,:,:,:,:),U_localNVar(1:PP_nVar,:,:,:,:))
    ELSE
      U_localNVar(1:PP_nVar,:,:,:,:) = U_local(1:PP_nVar,:,:,:,:)
    END IF
#if USE_RW
    IF (RestartTurb) THEN
      ALLOCATE(Uturb(nVarTurb,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
      Uturb(1:nVarTurb,:,:,:,:) = U_local(PP_nVar+1:PP_nVar+nVarTurb,:,:,:,:)
    END IF
#endif
    DEALLOCATE(U_local)
    ALLOCATE(U_local(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
    U_local = U_localNVar
    DEALLOCATE(U_localNVar)
  END IF

  ! Expand the solution if we read a restart file from a different equation system (RANS --> Navier-Stokes)
  IF (PP_nVar.GT.nVar_Restart) THEN
    ALLOCATE(U_localNVar(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
    U_localNVar(1:nVar_Restart,:,:,:,:) = U_local(1:nVar_Restart,:,:,:,:)
#if EQNSYSNR == 3
    ! Plugging in the constant mu_tilda
    U_localNVar(PP_nVar,:,:,:,:)   = MuTilda
#endif
    DEALLOCATE(U_local)
    ALLOCATE(U_local(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
    U_local = U_localNVar
    DEALLOCATE(U_localNVar)
  END IF

  ! Read in state
  IF(.NOT. InterpolateSolution)THEN
    ! No interpolation needed, read solution directly from file
#if PP_dim == 3
    IF (HSize_proc(4).EQ.1) THEN
      ! FLEXI compiled 3D, but data is 2D => expand third space dimension
      CALL ExpandArrayTo3D(5,(/PP_nVar,PP_N+1,PP_N+1,1,nElems/),4,PP_N+1,U_local,U)
    ELSE
      ! FLEXI compiled 3D + data 3D
      U = U_local
    END IF
#else
    IF (HSize_proc(4).EQ.1) THEN
      ! FLEXI compiled 2D + data 2D
      U = U_local
    ELSE
      ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
      CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVar,PP_N,PP_N,PP_N,nElems/),4,U_local)
      U = U_local
    END IF
#endif
#if GCL
#if PP_dim == 3
    IF (HSize_procJac(4).EQ.1) THEN
      ! FLEXI compiled 3D, but data is 2D => expand third space dimension
      IF (foundJac) CALL ExpandArrayTo3D(5,(/1,PP_N+1,PP_N+1,1,nElems/),4,PP_N+1,Jac_local,Jac)
    ELSE
      ! FLEXI compiled 3D + data 3D
      IF (foundJac) Jac = Jac_local
    END IF
#else
    IF (HSize_procJac(4).EQ.1) THEN
      ! FLEXI compiled 2D + data 2D
      IF (foundJac) Jac = Jac_local
    ELSE
      ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
      IF (foundJac) THEN
        CALL to2D_rank5((/1,0,0,0,1/),(/1,PP_N,PP_N,PP_N,nElems/),4,Jac_local)
        Jac = Jac_local
      END IF
    END IF
#endif /*PP_dim == 3*/
#endif /*GCL*/
  ELSE ! InterpolateSolution
    IF (RestartTurb) CALL CollectiveStop(__STAMP__,'Interpolation not supported for turbulent quantities. Non-linear operation!')
    ! We need to interpolate the solution to the new computational grid
    SWRITE(UNIT_stdOut,*)'Interpolating solution from restart grid with N=',N_restart,' to computational grid with N=',PP_N

    CALL GetVandermonde(N_Restart, NodeType_Restart,PP_N,      NodeType,         &
                        Vdm_NRestart_N,     modal=.TRUE.)
    CALL GetVandermonde(3*Ngeo,    NodeType,        N_Restart, NodeType_Restart, &
                        Vdm_3Ngeo_NRestart, modal=.TRUE.)

#if PP_dim == 3
    IF (HSize_proc(4).EQ.1) THEN
      ! FLEXI compiled 3D, but data is 2D => expand third space dimension
      ! use temporary array 'U_local2' to store 3D data
      ALLOCATE(U_local2(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,nElems))
      CALL ExpandArrayTo3D(5,HSize_proc,4,N_Restart,U_local,U_local2)
      ! Reallocate 'U_local' to 3D and mv data from U_local2 to U_local
      DEALLOCATE(U_local)
      ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,nElems))
      U_local = U_local2
      DEALLOCATE(U_local2)
    END IF
#else
    IF (HSize_proc(4).NE.1) THEN
      ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
      CALL to2D_rank5((/1,0,0,0,1/),(/PP_nVar,N_Restart,N_Restart,N_Restart,nElems/),4,U_local)
    END IF
#if GCL
    IF (HSize_procJac(4).NE.1) THEN
      ! FLEXI compiled 2D, but data is 3D => reduce third space dimension
      IF (foundJac) THEN
        CALL to2D_rank5((/1,0,0,0,1/),(/1,N_Restart,N_Restart,N_Restart,nElems/),4,Jac_local)
      END IF
    END IF
#endif /*GCL*/
#endif /*PP_dim == 3*/

#if GCL
    ! Transform Jacobian
    DO iElem=1,nElems
      CALL ChangeBasisVolume(1,N_Restart,PP_N,Vdm_NRestart_N,Jac_local(:,:,:,:,iElem),Jac(:,:,:,:,iElem))
    END DO
    ! Set inverse of Jacobian
    sJ(:,:,:,:,0) = 1./Jac(1,:,:,:,:)
#endif
    ! Transform solution to refspace and project solution to N
    ! For conservativity deg of detJac should be identical to EFFECTIVE polynomial deg of solution
    ! (e.g. beware when filtering the jacobian )
    IF(N_Restart.GT.PP_N)THEN
      ALLOCATE(JNR(1,0:N_Restart,0:N_Restart,0:N_Restart*(PP_dim-2)))
      DO iElem=1,nElems
        IF (FV_Elems(iElem).EQ.0) THEN ! DG element
          CALL ChangeBasisVolume(1,3*Ngeo,N_Restart,Vdm_3Ngeo_NRestart,detJac_Ref(:,:,:,:,iElem),JNR)
          DO k=0,N_Restart*(PP_dim-2); DO j=0,N_Restart; DO i=0,N_Restart
            U_local(:,i,j,k,iElem)=U_local(:,i,j,k,iElem)*JNR(1,i,j,k)
          END DO; END DO; END DO
          CALL ChangeBasisVolume(PP_nVar,N_Restart,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
#if FV_ENABLED
        ELSE ! FV element
          CALL SupersampleFVCell(U_local(:,:,:,:,iElem),U(:,:,:,:,iElem),N_Restart,PP_N,NFVRestartSuper)
#endif
        END IF
      END DO
      DEALLOCATE(JNR)
      ! Transform back
      CALL ApplyJacobianCons(U,toPhysical=.TRUE.,FVE=0)
    ELSE
      DO iElem=1,nElems
        IF (FV_Elems(iElem).EQ.0) THEN ! DG element
          CALL ChangeBasisVolume(PP_nVar,N_Restart,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
#if FV_ENABLED
        ELSE ! FV element
          CALL SupersampleFVCell(U_local(:,:,:,:,iElem),U(:,:,:,:,iElem),N_Restart,PP_N,NFVRestartSuper)
#endif
        END IF
      END DO
    END IF

    DEALLOCATE(U_local)
    SWRITE(UNIT_stdOut,*)'DONE!'
  END IF
  CALL CloseDataFile()

  IF (RestartMean.AND..NOT.postiMode) THEN
    SWRITE(UNIT_stdOut,'(A)') 'Restart from mean file successful.'
  END IF

#if !(USE_PARTICLES)
  ! Delete all files that will be rewritten --> moved to particle_restart.f90 since we need it there
  IF (doFlushFiles_loc) CALL FlushFiles(RestartTime)
#endif
ELSE
#if !(USE_PARTICLES)
  ! Delete all files since we are doing a fresh start --> moved to particle_restart.f90 since we need it there ! Delete all files since we are doing a fresh start --> moved to particle_restart.f90 since we need it there
  IF (doFlushFiles_loc) CALL FlushFiles()
#endif
END IF
END SUBROUTINE Restart


#if FV_ENABLED
!==================================================================================================================================
!> This routine will take a FV element with a certain number of subcells and convert it to a different number of subcells.
!> Is used during the restart from one polynomial degree to another.
!> The procedure is as follows: The old solution will be supersampled with an adjustable number of equidistant points.
!> The new solution is then simply calculated by taking the mean value of the supersampling points inside of the new subcell.
!> Attention: This procedure is not conservative! We do not take the Jacobian into account.
!==================================================================================================================================
SUBROUTINE SupersampleFVCell(UOld,UNew,NOld,NNew,NSuper)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: UOld(1:PP_nVar,0:NOld,0:NOld,0:ZDIM(NOld)) !< One FV element on NOld
REAL,INTENT(OUT)   :: UNew(1:PP_nVar,0:NNew,0:NNew,0:ZDIM(NNew)) !< FV Element on NNew
INTEGER,INTENT(IN) :: NOld                                       !< Old polynomial degree
INTEGER,INTENT(IN) :: NNew                                       !< New polynomial degree
INTEGER,INTENT(IN) :: NSuper                                     !< Polynomial degree for supersampling
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: deltaXiOld,deltaXi,deltaXiSuper,xiSuper
REAL                :: U_mean(PP_nVar)
INTEGER             :: iOld,jOld,kOld,iSuper,jSuper,kSuper
INTEGER             :: i,j,k
!==================================================================================================================================
! Supersample the old FV solution (with (NSuper+1)**dim superampling points in each new
! sub cell), then take the mean value
deltaXiOld = 2.0/(REAL(NOld)+1.)       ! Length (in reference space) of a FV element in the old element
deltaXi          = 2.0/(REAL(NNew)+1.)       ! Length (in reference space) of a FV element in the new element
deltaXiSuper     = deltaXi/(REAL(NSuper)+1.) ! Spacing (in reference space) between supersampling points
DO k=0,ZDIM(NNew); DO j=0,NNew; DO i=0,NNew
  U_mean = 0.
  DO kSuper=0,ZDIM(NSuper); DO jSuper=0,NSuper; DO iSuper=0,NSuper
    ! Calculate the index that the current super sampling point has in the old solution
    xiSuper = (REAL(i)*deltaXi) + (REAL(iSuper)+0.5)*deltaXiSuper
    iOld = INT(xiSuper/deltaXiOld)
    xiSuper = (REAL(j)*deltaXi) + (REAL(jSuper)+0.5)*deltaXiSuper
    jOld = INT(xiSuper/deltaXiOld)
#if PP_dim == 3
    xiSuper = (REAL(k)*deltaXi) + (REAL(kSuper)+0.5)*deltaXiSuper
    kOld = INT(xiSuper/deltaXiOld)
#else
    kOld = 0
#endif
    ! Calculate sum for mean value
    U_mean(:) = U_mean(:)+UOld(:,iOld,jOld,kOld)
  END DO; END DO; END DO! iSuper,jSuper,kSuper=0,NOld
  UNew(:,i,j,k) = U_mean(:)/REAL((NSuper+1)**PP_dim)
END DO; END DO; END DO! i,j,k=0,NNew
END SUBROUTINE SupersampleFVCell
#endif

!==================================================================================================================================
!> Finalizes variables necessary for restart subroutines
!==================================================================================================================================
SUBROUTINE FinalizeRestart()
! MODULES
USE MOD_Restart_Vars
#if USE_RW
USE MOD_DG_Vars,     ONLY:Uturb
#endif
IMPLICIT NONE
!==================================================================================================================================
#if USE_RW
SDEALLOCATE(Uturb)
#endif

RestartInitIsDone = .FALSE.
END SUBROUTINE FinalizeRestart

END MODULE MOD_Restart
