!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
!> Routines that handle restart capabilities.
!==================================================================================================================================
MODULE MOD_CalcTurb_Restart
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

INTERFACE GetDataProps
  MODULE PROCEDURE GetDataProps
END INTERFACE

PUBLIC :: InitRestart,FinalizeRestart,Restart,GetDataProps
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Initialize all necessary information to perform the restart.
!>
!> The routine also checks if the node type and polynomial degree of the restart file is the same than in the current simulation.
!> For TurbMode 1 & 2, a flag InterpolateSolution is set. This will be used by the actual Restart routine.
!> For TurbMode 3 & 4, interpolation results in artefacts and the tool will abort here.
!==================================================================================================================================
SUBROUTINE InitRestart(RestartFile_in,ArrayName_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CalcTurb_Vars,      ONLY: nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,TurbMode
USE MOD_HDF5_Input,         ONLY: ISVALIDHDF5FILE
USE MOD_HDF5_Input,         ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID
USE MOD_Interpolation_Vars, ONLY: NodeType
USE MOD_Mesh_Vars,          ONLY: nGlobalElems
USE MOD_ReadInTools,        ONLY: GETLOGICAL,GETREAL,GETREALARRAY
USE MOD_Restart_Vars,       ONLY: RestartFile,RestartTime,InterpolateSolution,RestartInitIsDone!,nElems_Restart
USE MOD_Restart_Vars,       ONLY: DoRestart
!USE MOD_Restart_Vars
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: ArrayName_in      !< array name to use in state file
CHARACTER(LEN=*),INTENT(IN)          :: RestartFile_in    !< state file to restart from
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                              :: validHDF5
CHARACTER(LEN=255)                   :: ArrayName
!==================================================================================================================================
! TKE means we have to read the 'Mean' array
IF (ArrayName_in.EQ.'TKE') THEN
    ArrayName = 'Mean'
ELSE
    ArrayName = ArrayName_in
END IF

RestartFile = RestartFile_in

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'
SWRITE(UNIT_StdOut,'(A,A,A)')' | Restarting from file "',TRIM(RestartFile),'":'

! Check if restart file is a valid state
validHDF5 = ISVALIDHDF5FILE(RestartFile)
IF(.NOT.validHDF5) CALL CollectiveStop(__STAMP__,'ERROR - Restart file not a valid state file.')

! Set flag indicating a restart to other routines
DoRestart = .TRUE.

! Read in parameters of restart solution
CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Read in attributes
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,ArrayName)

! Read in time from restart file
CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
CALL CloseDataFile()

! Check number of elements
IF (nElems_HDF5.NE.nGlobalElems) THEN
  CALL CollectiveStop(__STAMP__, "Restart File has different number of elements!")
END IF

! Check if we need to interpolate the restart file to our current polynomial degree and node type
IF ((N_HDF5.NE.PP_N) .OR. (TRIM(NodeType_HDF5).NE.TRIM(NodeType))) THEN
  ! Interpolation possible, set flag
  IF (TurbMode .LE. 2) THEN
    InterpolateSolution=.TRUE.
  ELSE
    CALL ABORT(__STAMP__,'Interpolation not supported for TurbMode=',TurbMode)
  END IF
END IF

RestartInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RESTART DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitRestart


!==================================================================================================================================
!> This routine performs the actual restart. It is called from the FLEXI main program just before the TimeDisc routine.
!>
!> For the restart there are 2 cases, depending on the flag InterpolateSolution set by InitRestart:
!> - No interpolation is necesary. Then simply read in the DG_Solution array from the restart file and store it in the solution
!>   array U.
!> - We need to interpolate the restart solution. If the polynomial degree of our computation is lower than in the restart file,
!>   a simple change basis is used to get the current solution U. If the polynomial degree is higher than in the restart file,
!>   special care is taken to ensure a conservative projection of the restart solution. To do this, the restart solution is
!>   transformed to reference space using the Jacobian built on a polynomial degree of 3*NGeo (so it can be represented exactly),
!>   then the change basis is applied. The resulting solution U is then transformed back to physical space.
!==================================================================================================================================
SUBROUTINE Restart(ArrayName_in )
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CalcTurb_Vars
USE MOD_ApplyJacobianCons,  ONLY: ApplyJacobianCons
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_DG_Vars,            ONLY: U
USE MOD_HDF5_Input,         ONLY: OpenDataFile,CloseDataFile,ReadArray,GetArrayAndName,ReadAttribute,GetDataSize
USE MOD_Interpolation,      ONLY: GetVandermonde
USE MOD_Interpolation_Vars, ONLY: NodeType
USE MOD_IO_HDF5
USE MOD_Mesh_Vars,          ONLY: offsetElem,detJac_Ref,Ngeo
USE MOD_Mesh_Vars,          ONLY: nElems,nGlobalElems
USE MOD_Restart_Vars
USE MOD_StringTools,        ONLY: STRICMP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)        :: ArrayName_in   !< array name to use in state file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                              :: FoundFluc = .FALSE.
INTEGER                              :: iElem,iFluc,i,j,k
INTEGER                              :: HSize_proc(5)
INTEGER                              :: nVal(15)
REAL                                 :: Vdm_NRestart_N(0:PP_N,0:N_HDF5)
REAL                                 :: Vdm_3Ngeo_NRestart(0:N_HDF5,0:3*NGeo)
REAL,ALLOCATABLE                     :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE                     :: U_localTKE(:,:,:,:,:)
REAL,ALLOCATABLE                     :: U_localNVar(:,:,:,:,:)
REAL,ALLOCATABLE                     :: JNR(:,:,:,:),tmp(:)
CHARACTER(LEN=255),ALLOCATABLE       :: VarNamesFluc(:)
CHARACTER(LEN=255)                   :: ArrayName
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING RESTART...'

! TKE means we have to read the 'Mean' array
IF (ArrayName_in.EQ.'TKE') THEN
    ArrayName = 'Mean'
ELSE
    ArrayName = ArrayName_in
END IF

CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataSize(File_ID,ArrayName,nDims,HSize)

! Reduced sanity check. Keep the entire array and sort it out later to allow for more diverse use cases
!> Check if nGlobalElems match
HSize_proc    = INT(HSize)
IF (HSize_proc(5).NE.nGlobalElems) CALL Abort(__STAMP__,"Dimensions of restart file do not match!")
HSize_proc(5) = nElems

ALLOCATE(U_local(nVar_HDF5,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
CALL ReadArray(ArrayName,5,HSize_proc,OffsetElem,5,RealArray=U_local)

! Now pick the right array to read
SELECT CASE(ArrayName_in)
    CASE('DG_Solution')
        ! Truncate the solution if we read a restart file from a different equation system
        IF (PP_nVar.LT.nVar_HDF5) THEN
            ALLOCATE(U_localNVar(PP_nVar,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
            U_localNVar(1:PP_nVar,:,:,:,:) = U_local(1:PP_nVar,:,:,:,:)
            DEALLOCATE(U_local)
            ALLOCATE(U_local(PP_nVar,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
            U_local = U_localNVar
            DEALLOCATE(U_localNVar)
        END IF

    CASE('Mean')
        ! Truncate the array after the mean DG solution
        IF (PP_nVar.NE.nVar_HDF5) THEN
            ALLOCATE(U_localNVar(PP_nVar,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
            U_localNVar(1:PP_nVar,:,:,:,:) = U_local(1:PP_nVar,:,:,:,:)
            DEALLOCATE(U_local)
            ALLOCATE(U_local(PP_nVar,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
            U_local = U_localNVar
            DEALLOCATE(U_localNVar)
        END IF

!    CASE('Fluc')
!        ! Find out which names we have available
!        CALL GetArrayAndName('Fluc'      ,'VarNames_Fluc'      ,nVal          ,tmp         ,VarNamesFluc)
!        ! Fill the array with the fluc equivalent of the DG solution
!        ALLOCATE(U_localNVar(PP_nVar,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
!        nValFluc = INT(HSize(1))
!        ! Density from mean array
!        U_localNVar(1,:,:,:,:) = UMean(1,:,:,:,:)
!        ! EnergyStagnationDensity from mean array
!        U_localNVar(5,:,:,:,:) = UMean(5,:,:,:,:)
!        DO iFluc=1,nValFluc
!            ! Fluc variables in conservative form
!            IF(STRICMP(TRIM(VarNamesFluc(iFluc)),'MomentumX')) THEN
!                U_localNVar(2,:,:,:,:) = SQRT(U_local(iFluc,:,:,:,:))
!                FoundFluc(  1)         = .TRUE.
!            END IF
!            IF(STRICMP(TRIM(VarNamesFluc(iFluc)),'MomentumY')) THEN
!                U_localNVar(3,:,:,:,:) = SQRT(U_local(iFluc,:,:,:,:))
!                FoundFluc(  2)         = .TRUE.
!            END IF
!            IF(STRICMP(TRIM(VarNamesFluc(iFluc)),'MomentumZ')) THEN
!                U_localNVar(4,:,:,:,:) = SQRT(U_local(iFluc,:,:,:,:))
!                FoundFluc(  3)         = .TRUE.
!            END IF
!        END DO
!        ! Check if all required fluc variables exist
!        IF (.NOT.(ALL(FoundFluc.EQV..TRUE.))) CALL CollectiveStop(__STAMP__,'Missing fluc variables in timeavg file. Aborting ...')
!        DEALLOCATE(U_local)
!        ALLOCATE(U_local(PP_nVar,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
!        U_local = U_localNVar
!        DEALLOCATE(U_localNVar)

    CASE('TKE')
        ! Fill the array with the fluc equivalent of the DG solution, replace EnergyStagnationDensity with TKE
        ! >Truncate the array after the mean DG solution
        IF (PP_nVar.LE.nVar_HDF5) THEN
            ALLOCATE(U_localNVar(PP_nVar,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
            U_localNVar(1:PP_nVar-1,:,:,:,:) = U_local(1:PP_nVar-1,:,:,:,:)
        ELSE
            CALL CollectiveStop(__STAMP__,'Missing conservative variable(s) in timeavg file. Aborting ...')
        END IF

        ! Now get TKE from the Fluc array
        ! >Read in attributes
        CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,'Fluc')
        CALL GetDataSize(File_ID,'Fluc',nDims,HSize)
        HSize_proc    = INT(HSize)
        HSize_proc(5) = nElems
        ! >Read actual data
        ALLOCATE(U_localTKE(nVar_HDF5,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
        CALL GetArrayAndName('Fluc','VarNames_Fluc',nVal,tmp,VarNamesFluc)
        CALL ReadArray      ('Fluc',5,HSize_proc,OffsetElem,5,RealArray=U_localTKE)

        FoundFluc = .FALSE.
        DO iFluc=1,nVal(1)
            ! TKE variables in conservative form. Keep the square since we are really looking for TKE
            IF(STRICMP(TRIM(VarNamesFluc(iFluc)),'TKE')) THEN
                U_localNVar(5,:,:,:,:) = U_localTKE(iFluc,:,:,:,:)
                FoundFluc = .TRUE.
            END IF
        END DO
        ! Check if all required fluc variables exist
        IF (.NOT.(FoundFluc.EQV..TRUE.)) CALL CollectiveStop(__STAMP__,'Missing TKE variables in timeavg file. Aborting ...')
        DEALLOCATE(U_local)
        ALLOCATE(U_local(1,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
        U_local = U_localNVar
        DEALLOCATE(U_localNVar)

END SELECT

! Read in state
IF(.NOT.InterpolateSolution) THEN
    ! No interpolation needed, read solution directly from file
    U = U_local
    ! We need mean values later, so keep them in their own array
    IF (ArrayName_in.EQ.'Mean') THEN
        SDEALLOCATE(UMean)
        ALLOCATE(UMean(PP_nVar,0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
        UMean = U
    ELSEIF(ArrayName_in.EQ.'TKE') THEN
        ALLOCATE(TKE  (        0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
        TKE   = U(5,:,:,:,:)
    END IF
ELSE
    ! We need to interpolate the solution to the new computational grid
    SWRITE(UNIT_stdOut,*)'  Interpolating solution from restart grid with N=',N_HDF5,' to computational grid with N=',PP_N
    CALL GetVandermonde(N_HDF5,    NodeType_HDF5,  PP_N,      NodeType,         Vdm_NRestart_N,     modal=.TRUE.)
    CALL GetVandermonde(3*Ngeo,    NodeType,       N_HDF5,    NodeType_HDF5,    Vdm_3Ngeo_NRestart, modal=.TRUE.)
    ! Transform solution to refspace and project solution to N
    ! For conservativity deg of detJac should be identical to EFFECTIVE polynomial deg of solution
    ! (e.g. beware when filtering the jacobian )
    IF(N_HDF5.GT.PP_N)THEN
        ALLOCATE(JNR(1,0:N_HDF5,0:N_HDF5,0:N_HDF5*(PP_dim-2)))
        ! We need mean values later, so keep them in their own array
        IF (ArrayName_in.EQ.'Mean') THEN
            SDEALLOCATE(UMean)
            ALLOCATE(UMean(PP_nVar,0:PP_N,0:PP_N,0:ZDIM(PP_N),nElems))
        END IF
        IF (ArrayName_in.EQ.'TKE')  ALLOCATE(TKE  (        0:PP_N,0:PP_N,0:ZDIM(PP_N),nElems))
        ! Interpolate
        DO iElem=1,nElems
            CALL ChangeBasisVolume(1,3*Ngeo,N_HDF5,Vdm_3Ngeo_NRestart,detJac_Ref(:,:,:,:,iElem),JNR)
            DO k=0,N_HDF5*(PP_dim-2); DO j=0,N_HDF5; DO i=0,N_HDF5
                U_local(:,i,j,k,iElem)=U_local(:,i,j,k,iElem)*JNR(1,i,j,k)
            END DO; END DO; END DO
            CALL ChangeBasisVolume(PP_nVar,N_HDF5,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
        END DO
        DEALLOCATE(JNR)
        ! Transform back
        CALL ApplyJacobianCons(U,toPhysical=.TRUE.,FVE=0)
        IF (ArrayName_in.EQ.'Mean') UMean = U
        IF (ArrayName_in.EQ.'TKE')  TKE   = U(5,:,:,:,:)
    ELSE
        DO iElem=1,nElems
            CALL ChangeBasisVolume(PP_nVar,N_HDF5,PP_N,Vdm_NRestart_N,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
        END DO
        ! We need mean values later, so keep them in their own array
        IF (ArrayName_in.EQ.'Mean') THEN
            SDEALLOCATE(UMean)
            ALLOCATE(UMean(PP_nVar,0:PP_N,0:PP_N,0:ZDIM(PP_N),nElems))
            UMean = U
        ELSE IF(ArrayName_in.EQ.'TKE') THEN
            ALLOCATE(TKE  (       0:PP_N,0:PP_N,0:ZDIM(PP_N),nElems))
            TKE   = U(5,:,:,:,:)
        END IF
    END IF

    DEALLOCATE(U_local)
    SWRITE(UNIT_stdOut,*)'DONE!'
END IF

CALL CloseDataFile()

END SUBROUTINE Restart

!==================================================================================================================================
!> Subroutine to determine HDF5 dataset properties in Flexi terminology
!==================================================================================================================================
SUBROUTINE GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,ArrayName)
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input,         ONLY: ReadAttribute
USE MOD_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(OUT)                     :: nVar_HDF5     !< number of variables
INTEGER,INTENT(OUT)                     :: N_HDF5        !< polynomial degree
INTEGER,INTENT(OUT)                     :: nElems_HDF5   !< number of elements
CHARACTER(LEN=*),INTENT(IN)             :: ArrayName     !< array name to use in state file
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT)   :: NodeType_HDF5 !< nodetype string
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                 :: Rank
INTEGER(HID_T)                          :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,A)')' GET SIZE OF DATA IN HDF5 FILE...'

! Read in attributes
CALL H5DOPEN_F(File_ID, ArrayName, Dset_ID, iError)

! Get the data space of the dataset.
CALL H5DGET_SPACE_F(Dset_ID, FileSpace, iError)

! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, Rank, iError)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Rank of database',' | ',Rank,' | HDF5    |'

! Get size and max size of data space
Dims   =0
DimsMax=0
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Dims(1:Rank), DimsMax(1:Rank), iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(Dset_ID, iError)

IF(PRESENT(NodeType_HDF5)) THEN
  ! Read in NodeType
  CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType_HDF5)
END IF

! Display data
! nVar = first array index
!nVar_HDF5 = INT(Dims(1),4)
CHECKSAFEINT(Dims(1),4)
nVar_HDF5 = INT(Dims(1),4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Number of variables nVar',' | ',nVar_HDF5,' | HDF5    |'

! N = index 2-4 of array, is expected to have the same value for each direction
CHECKSAFEINT(Dims(2)-1,4)
N_HDF5 = INT(Dims(2)-1,4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Polynomial degree N',' | ',N_HDF5,' | HDF5    |'
IF(PRESENT(NodeType_HDF5)) THEN
  SWRITE(UNIT_stdOut,'(A3,A30,A3,A33,A13)')' | ','          Node type',' | ',TRIM(NodeType_HDF5),' | HDF5    |'
END IF

! nElems = index 5 of array
CHECKSAFEINT(Dims(5),4)
nElems_HDF5 = INT(Dims(5),4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','GeometricnElems',' | ',nElems_HDF5,' | HDF5    |'

SWRITE(UNIT_stdOut,'(A)')' DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE GetDataProps


!==================================================================================================================================
!> Finalizes variables necessary for restart subroutines
!==================================================================================================================================
SUBROUTINE FinalizeRestart()
! MODULES
USE MOD_Restart_Vars
IMPLICIT NONE
!==================================================================================================================================
RestartInitIsDone = .FALSE.
END SUBROUTINE FinalizeRestart


END MODULE MOD_CalcTurb_Restart
