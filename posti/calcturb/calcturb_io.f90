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

!===================================================================================================================================
!> Contains routines that read the requested data and interpolate it to the required degree and node type. All data is first read
!> and interpolated to NCalc on the node type the tool is compiled with. In case a different node type is requested for output, the
!> data is interpolated a second time before performing nonlinear calculations to avoid aliasing.
!===================================================================================================================================
MODULE MOD_CalcTurb_IO
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ReadStateFile
  MODULE PROCEDURE ReadStateFile
END INTERFACE

INTERFACE WriteStateFile
  MODULE PROCEDURE WriteStateFile
END INTERFACE

PUBLIC:: ReadStateFile
PUBLIC:: WriteStateFile

CONTAINS

!===================================================================================================================================
!> Read quantities for all DG elements from state file
!> 1. Read quantities on PP_N
!> 2. Calculate gradients on PP_N
!> 3. Interpolate quantities and gradients to NCalc on same node type
!>(4. Interpolate quantities and gradients to NVisu on visu node type)
!===================================================================================================================================
SUBROUTINE ReadStateFile(ParameterFile,StateFile,ArrayName)
! Modules
USE MOD_Globals       ,ONLY: MPIROOT
USE MOD_PreProc
USE MOD_CalcTurb_Restart
USE MOD_CalcTurb_Vars
USE MOD_DG            ,ONLY: InitDG,DGTimeDerivative_weakForm,FinalizeDG
USE MOD_EOS           ,ONLY: DefineParametersEos
USE MOD_Equation      ,ONLY: DefineParametersEquation,InitEquation,FinalizeEquation
USE MOD_Exactfunc     ,ONLY: DefineParametersExactFunc
USE MOD_Filter        ,ONLY: DefineParametersFilter,InitFilter,FinalizeFilter
USE MOD_Lifting       ,ONLY: DefineParametersLifting,InitLifting,FinalizeLifting
USE MOD_Indicator     ,ONLY: DefineParametersIndicator,InitIndicator,FinalizeIndicator
USE MOD_Interpolation ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5       ,ONLY: DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Mesh          ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars     ,ONLY: MeshFile
USE MOD_Mortar        ,ONLY: InitMortar,FinalizeMortar
USE MOD_MPI           ,ONLY: DefineParametersMPI
USE MOD_Output        ,ONLY: DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Overintegration,ONLY:DefineParametersOverintegration,InitOverintegration,FinalizeOverintegration
USE MOD_ReadInTools   ,ONLY: prms,FinalizeParameters
USE MOD_Restart       ,ONLY: DefineParametersRestart
USE MOD_Restart_Vars  ,ONLY: RestartTime
#if USE_MPI
USE MOD_MPI           ,ONLY: InitMPIvars,FinalizeMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: ParameterFile
CHARACTER(LEN=*),INTENT(IN)   :: StateFile
CHARACTER(LEN=*),INTENT(IN)   :: ArrayName
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! finalize everything to make sure all arrays are empty
CALL FinalizeInterpolation()
CALL FinalizeMortar()
CALL FinalizeRestart()
#if USE_MPI
! We might only need this on the first run
    CALL FinalizeMPI()
#endif
CALL FinalizeIndicator()
CALL FinalizeEquation()
CALL FinalizeDG()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
CALL FinalizeLifting()
CALL FinalizeOutput()

! read options from parameter file
CALL FinalizeParameters()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersFilter()
CALL DefineParametersOverintegration()
CALL DefineParametersIndicator()
CALL DefineParametersEos()
CALL DefineParametersEquation()
CALL DefineParametersExactFunc()
CALL DefineParametersLifting()
CALL prms%read_options(ParameterFile)

! Initialization of I/O routines
CALL InitIOHDF5()
CALL InitInterpolation(NCalc)
CALL InitMortar()
CALL InitOutput()
! We might only need this on the first run
    CALL FinalizeMesh()
    CALL InitMesh(meshMode=2,MeshFile_IN=MeshFile)
CALL InitRestart(StateFile,ArrayName)
CALL InitFilter()
CALL InitOverintegration()
CALL InitIndicator()
#if USE_MPI
! We might only need this on the first run
    CALL InitMPIvars()
#endif
CALL InitEquation()
CALL InitDG()
CALL InitLifting()

! Finally do the Restart
CALL Restart(ArrayName)
! TKE can not be lifted
IF (ArrayName.NE.'TKE') THEN
    SWRITE(*,*) "Call DGTimeDerivative_weakForm ... "
    CALL DGTimeDerivative_weakForm(RestartTime)
END IF
SWRITE(*,*) "DONE!"

! Finalize Parameters
CALL FinalizeParameters()

END SUBROUTINE ReadStateFile


!==================================================================================================================================
!> Subroutine to write the solution U to HDF5 format
!==================================================================================================================================
SUBROUTINE WriteStateFile(MeshFileName,SolutionArray,ArrayName)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CalcTurb_Vars
USE MOD_ChangeBasisByDim  ,ONLY: ChangeBasisVolume
USE MOD_HDF5_Output       ,ONLY: GatheredWriteArray,MarkWriteSuccessfull
USE MOD_Mesh_Vars         ,ONLY: offsetElem,nGlobalElems,sJ,nElems
USE MOD_Output_Vars       ,ONLY: ProjectName,NOut,Vdm_N_NOut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName   !< file name of mesh used for simulation
REAL,TARGET,INTENT(IN)         :: SolutionArray(nVarTurb,0:PP_N,0:PP_N,0:ZDIM(PP_N),nElems)
CHARACTER(LEN=*),INTENT(IN)    :: ArrayName      !< name of array to be written
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,i,j,k
INTEGER                        :: nVal(5)
REAL                           :: StartT,EndT
REAL                           :: JN(1,0:NCalc,0:NCalc,0:NCalc),JOut(1,0:NOut,0:NOut,0:ZDIM(NOut))
REAL                           :: Utmp(nVarTurb,0:NCalc,0:NCalc,0:NCalc)
REAL,POINTER                   :: UOut(:,:,:,:,:)
CHARACTER(LEN=255)             :: FileName
!==================================================================================================================================
IF(MPIRoot)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' WRITE STATE TO HDF5 FILE...'
  GETTIME(StartT)
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_calcturb_',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton(TRIM(FileName),'State',nVarTurb,NOut,varnames_loc,MeshFileName,OutputTime,ArrayName)

! Set size of output
nVal=(/nVarTurb,NOut+1,NOut+1,ZDIM(NOut)+1,nElems/)

! build output data
IF(NOut.NE.NCalc) THEN
    ! Project JU and J to NOut, compute U on Nout
    SWRITE(UNIT_stdOut,*)'  Interpolating solution from computational grid with N=',NCalc,' to output grid with N=',NOut
    ALLOCATE(UOut(nVarTurb,0:NOut,0:NOut,0:ZDIM(NOut),nElems))
    DO iElem=1,nElems
        JN(1,:,:,:)=1./sJ(:,:,:,iElem,0)
        DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
            Utmp(:,i,j,k)=SolutionArray(:,i,j,k,iElem)*JN(1,i,j,k)
        END DO; END DO; END DO
        CALL ChangeBasisVolume(nVarTurb,PP_N,NOut,Vdm_N_NOut,Utmp,UOut(1:nVarTurb,:,:,:,iElem))
        ! Jacobian
        CALL ChangeBasisVolume(1,PP_N,NOut,Vdm_N_NOut,JN,JOut)
        DO k=0,ZDIM(NOut); DO j=0,NOut; DO i=0,NOut
            UOut(:,i,j,k,iElem)=UOut(:,i,j,k,iElem)/JOut(1,i,j,k)
        END DO; END DO; END DO
    END DO !iElem
    
ELSE
    ! write state on same polynomial degree as the solution
    UOut => SolutionArray
END IF

! Reopen file and write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName=ArrayName, rank=5,&
                        nValGlobal=(/PP_nVar,NOut+1,NOut+1,NOut+1,nGlobalElems/),&
                        nVal=nVal                                              ,&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=UOut)
    
! Deallocate UOut only if we did not point to U
IF(NOut.NE.NCalc) DEALLOCATE(UOut)

IF(MPIRoot)THEN
    CALL MarkWriteSuccessfull(FileName)
    GETTIME(EndT)
    WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' DONE!  [',EndT-StartT,'s]'
END IF

END SUBROUTINE WriteStateFile


!==================================================================================================================================
!> Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!==================================================================================================================================
SUBROUTINE GenerateFileSkeleton(FileName,TypeString,nVar,NData,StrVarNames,MeshFileName,OutputTime,&
                                DataName,Dataset)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_HDF5_Output        ,ONLY: WriteAttribute,WriteHeader
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_IO_HDF5
USE MOD_Mesh_Vars          ,ONLY: nGlobalElems
USE MOD_Output_Vars        ,ONLY: userblock_total_len
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName           !< Name of file to create
CHARACTER(LEN=*),INTENT(IN)    :: TypeString         !< Type of file to be created (state,timeaverage etc.)
INTEGER,INTENT(IN)             :: nVar               !< Number of variables
INTEGER,INTENT(IN)             :: NData              !< Polynomial degree of data
CHARACTER(LEN=*)               :: StrVarNames(nVar)  !< Variabel names
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName       !< Name of mesh file
REAL,INTENT(IN)                :: OutputTime         !< Time of output
CHARACTER(LEN=*),INTENT(IN)    :: DataName           !< Name of the dataset array
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: Dataset      !< Name of the dataset
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID,FileSpace,HDF5DataType
INTEGER(HSIZE_T)               :: Dimsf(5)
CHARACTER(LEN=255)             :: Dataset_Str,Varname_Str
LOGICAL                        :: withUserblock_loc,create_loc
!==================================================================================================================================
! Create file
create_loc=.TRUE.
withUserblock_loc=.FALSE.
Dataset_Str=TRIM(DataName)
Varname_Str='VarNames'
IF(PRESENT(Dataset))THEN
  Dataset_Str=TRIM(Dataset)
  Varname_Str='VarNames_'//TRIM(DataSet)
END IF

CALL OpenDataFile(TRIM(FileName),create=create_loc,single=.TRUE.,readOnly=.FALSE.,&
                  userblockSize=MERGE(userblock_total_len,0,withUserblock_loc))

! Preallocate the data space for the dataset.
IF(output2D) THEN
  Dimsf=(/nVar,NData+1,NData+1,1,nGlobalElems/)
ELSE
  Dimsf=(/nVar,NData+1,NData+1,NData+1,nGlobalElems/)
END IF

CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_DOUBLE
CALL H5DCREATE_F(File_ID,TRIM(Dataset_Str), HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL WriteAttribute(File_ID,TRIM(Varname_Str),nVar,StrArray=StrVarNames)

! Write default attributes only if file is created
IF(create_loc)THEN

  ! Write file header
  CALL WriteHeader(TRIM(TypeString),File_ID)

  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttribute(File_ID,'N',1,IntScalar=PP_N)
  CALL WriteAttribute(File_ID,'Dimension',1,IntScalar=PP_dim)
  CALL WriteAttribute(File_ID,'Time',1,RealScalar=OutputTime)
  CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=(/MeshFileName/))
  CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttribute(File_ID,'NComputation',1,IntScalar=PP_N)
END IF

CALL CloseDataFile()

END SUBROUTINE GenerateFileSkeleton

END MODULE MOD_CalcTurb_IO
