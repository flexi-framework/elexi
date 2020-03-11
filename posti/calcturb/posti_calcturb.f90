!===================================================================================================================================
! Copyright (c) 2010-2019  Prof. Claus-Dieter Munz
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
!===================================================================================================================================
#include "flexi.h"

!===================================================================================================================================
!> This tool will take a TimeAvg (and multiple state) file and calculate turbulence quantities from Mean and MeanSquare values
!===================================================================================================================================
!
! There are quantities that require the computation of gradients. In this case the DG operator 'DGTimeDerivative_weakForm' is
! called once to fill the gradients and the reconstruction of the FV subcell method. This requires the initialization of several
! modules of FLEXI. U is read via a call of 'Restart'. In the DGTimeDerivative_weakForm the primitive quantities U_Prim and
! gradUx/y/z as well as gradUxi/eta/zeta are filled. These are used to calculate the remaining quantities.
!
! * The calculation of derived quantities is performed on a arbitrary polynomial degree NCalc and afterwards interpolated to NVisu.
!   Default is PP_N. These require to reconstruct the solution first to the visu grid and afterwards can calculate the derived
!   quantities on the NVisu_FV grid.
!
!===================================================================================================================================
PROGRAM Posti_CalcTurb
! MODULES
USE ISO_C_BINDING
USE MOD_ISO_VARYING_STRING
USE MOD_Globals               ,ONLY: CollectiveStop,UNIT_stdOut,SetStackSizeUnlimited,StartTime,FlexiTime,iError,MPIROOT
USE MOD_Globals               ,ONLY: ParameterFile
USE MOD_CalcTurb              ,ONLY: DefineCalcTurb,InitCalcTurb,FinalizeCalcTurb,ALMOSTEQUAL,ALMOSTZERO
USE MOD_CalcTurb_IO           ,ONLY: ReadStateFile,WriteStateFile
USE MOD_CalcTurb_Vars
USE MOD_Commandline_Arguments
USE MOD_DG_Vars               ,ONLY: U
USE MOD_EOS_Vars              ,ONLY: mu0,KappaM1
USE MOD_Mesh_Vars             ,ONLY: MeshFile,nElems
USE MOD_Lifting_Vars          ,ONLY: gradUx,gradUy,GradUz
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_PreProc
USE MOD_ReadInTools           ,ONLY: GetReal
USE MOD_StringTools           ,ONLY: STRICMP,GetFileExtension,set_formatting,clear_formatting
!USE MOD_Mesh_Vars             ,ONLY: Elem_xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: a,i,j,k,p,q,iArg,iElem
INTEGER                        :: skipArgs
REAL                           :: Time                                       ! Used to measure execution time
REAL                           :: Pressure,rho
REAL                           :: Vel(1:3),GradVel(1:3,1:3)
REAL                           :: S(1:3,1:3)                                 ! Strain rate tensor S (symmetric)
REAL                           :: Sd(1:3,1:3)                                ! Deviatoric part of the strain rate tensor S
REAL                           :: divU                                       ! Divergence of velocity vector
REAL                           :: eps1                                       ! Integrand: Sd:Sd
REAL                           :: eps2                                       ! Integrand: (div u)^2
REAL                           :: eps3                                       ! Integrand: p*(div u)
REAL                           :: u_tens, s_tens, sd_tens                    ! matrix : matrix product, integrands of Gradvel, S, Sd
!CHARACTER(LEN=255)             :: ParameterFile
CHARACTER(LEN=255)             :: StateFile
CHARACTER(LEN=255)             :: varnames_gen(5)
!CHARACTER(LEN=255)             :: ArrayName
#if !USE_MPI
INTEGER                        :: MPI_COMM_WORLD = 0
#endif
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()

! Abort on currently unsupported options
#if EQNSYSNR != 2
CALL CollectiveStop(__STAMP__, 'Currently only Navier-Stokes equation system supported. Please recompile with eqn=navierstokes')
#endif
#if FV_ENABLED
CALL CollectiveStop(__STAMP__, 'FV mode not yet implemented. Please disable FV in CMake and recompile.')
#endif

! Check and pass all command line arguments
CALL ParseCommandlineArguments()
! Parameter file
IF(STRICMP(GetFileExtension(Args(1)),'ini')) THEN
    ParameterFile = Args(1)
    skipArgs      = 1           ! First argument is the parameter file. Skip it for state files
ELSE IF(STRICMP(GetFileExtension(Args(1)),'h5')) THEN
    skipArgs      = 0           ! Do not skip a argument. First argument is a .h5 file
ELSE
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: calcturb [prm-file] statefile [statefiles]')
END IF

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(A)') &
".___________. __    __  .______      .______    __    __   __       _______ .__   __.   ______  _______ "
SWRITE(UNIT_stdOut,'(A)') &
"|           ||  |  |  | |   _  \     |   _  \  |  |  |  | |  |     |   ____||  \ |  |  /      ||   ____|"
SWRITE(UNIT_stdOut,'(A)') &
"`---|  |----`|  |  |  | |  |_)  |    |  |_)  | |  |  |  | |  |     |  |__   |   \|  | |  ,----'|  |__   "
SWRITE(UNIT_stdOut,'(A)') &
"    |  |     |  |  |  | |      /     |   _  <  |  |  |  | |  |     |   __|  |  . `  | |  |     |   __|  "
SWRITE(UNIT_stdOut,'(A)') &
"    |  |     |  `--'  | |  |\  \----.|  |_)  | |  `--'  | |  `----.|  |____ |  |\   | |  `----.|  |____ "
SWRITE(UNIT_stdOut,'(A)') &
"    |__|      \______/  | _| `._____||______/   \______/  |_______||_______||__| \__|  \______||_______|"
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Measure init duration
StartTime=FLEXITIME()

! Define parameters needed
CALL DefineCalcTurb

! Call init with DG_Solution since this array should be globally available, if only as dummy
StateFile = Args(1+skipArgs)
CALL InitCalcTurb(ParameterFile,StateFile,ArrayName='DG_Solution')

varnames_gen(1) = 'Density'
varnames_gen(2) = 'MomentumX'
varnames_gen(3) = 'MomentumY'
varnames_gen(4) = 'MomentumZ'
varnames_gen(5) = 'EnergyStagnationDensity'

! Loop over all state files
DO iArg=1+skipArgs,nArgs
    CALL set_formatting("green")
    SWRITE(UNIT_stdOut,"(A,I0,A,I0,A)") "Processing File ", iArg-1, " of ", nArgs-1, "..."
    CALL clear_formatting()
    SWRITE(UNIT_stdOut,'(132("-"))')

    StateFile=Args(iArg)

    SELECT CASE(TurbMode)
!> =================================================================================================================================
!> Simple conversion mode to combine LES timeavg to restart to RANS. Sets epsilon to RestartEpsilon
!> =================================================================================================================================
        CASE (1)
            ! Conservative variables plus mu_tilda
            nVarTurb = 6

            ! Read state file
            IF (iArg.EQ.(1+skipArgs)) THEN
                ! Get epsilon desired for restart
                RestartEpsilon=GETREAL('RestartEpsilon'  ,'0.')
                ! Get mean values
                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING MEAN VALUES...'
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='Mean')
            END IF

            ! Set remaining varname
            ALLOCATE(varNames_loc(nVarTurb))
            varnames_loc(1:5) = varnames_gen(1:5)
            varnames_loc(6)   = 'mutilda'

            ! Fill restart epsilon
            ALLOCATE(USolution(nVarTurb,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
            USolution(1:5,:,:,:,:) = U(:,:,:,:,:)
            USolution(6  ,:,:,:,:) = RestartEpsilon

            ! Output the state file
            SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT...'
            CALL WriteStateFile(MeshFileName=TRIM(MeshFile),SolutionArray=USolution,ArrayName='DG_Solution')
            SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT DONE'

!> =================================================================================================================================
!> Dissipation calculation according to Canuto (1997)
!> FLEXI "should" already operate on mass-averaged fluctuations for compressible flows. CHECK THIS!
!> =================================================================================================================================
        CASE(2)
            ! Conservative variables plus k and epsilon
            nVarTurb = 7

            ! We expect one timeAvg file and several time-accurate files. Start by filling the mean array ONCE
            IF (iArg.EQ.(1+skipArgs)) THEN
                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING MEAN VALUES...'
                ! Get TKE so we can overwrite with mean values later
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='TKE')

                ! Get mean values
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='Mean')

                ! Allocate and fill array to the mean gradients
                ALLOCATE(GradUMeanx (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(GradUMeany (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(GradUMeanz (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                GradUMeanx = GradUx; GradUMeany = GradUy; GradUMeanz = GradUz;

                ! Allocate more arrays ONCE to hold epsilon
                ALLOCATE(EpsilonFluc(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonTmp (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonSum (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                EpsilonSum = 0.
                ALLOCATE(EpsilonFin   (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                epsilonFin   = 0.

                ! Allocate arrays to hold the dilatation and the vorticity
                ALLOCATE(dFluc(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(Vorticity(1:3,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))

                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING MEAN VALUES DONE'
                SWRITE(UNIT_stdOut,'(132("-"))')
            ELSE
                ! Now get the time-accurate values from the remaining state files
                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING STATE VALUES...'
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='DG_Solution')

                ! Nullify EpsilonFluc
                EpsilonFluc = 0.

                ! Get fluctuation of gradients. Linear operation, so we can operate on gradients
                DO iElem=1,nElems; DO k=0,ZDIM(NCalc); DO j=0,NCalc; DO i=0,NCalc
                    ! Work on velocity DIFFERENCES due to turbulence
                    GradVel(:,1) = GradUx(2:4,i,j,k,iElem) - gradUMeanx(2:4,i,j,k,iElem)
                    GradVel(:,2) = GradUy(2:4,i,j,k,iElem) - gradUMeany(2:4,i,j,k,iElem)
                    GradVel(:,3) = GradUz(2:4,i,j,k,iElem) - gradUMeanz(2:4,i,j,k,iElem)

                    dFluc(i,j,k,iElem) = GradVel(1,1) + GradVel(2,2) + GradVel(3,3)

                    ! compute vorticity
                    Vorticity(1,i,j,k,iElem)=GradVel(3,2) - GradVel(2,3)
                    Vorticity(2,i,j,k,iElem)=GradVel(1,3) - GradVel(3,1)
                    Vorticity(3,i,j,k,iElem)=GradVel(2,1) - GradVel(1,2)
                END DO; END DO; END DO; END DO

                ! Epsilon
                !< epsilon = epsilon_s      + epsilon_d
                !<         = nu*<vort*vort> + 4/3*nu*d^2
                !<         = nu*(<vort*vort>+ 4/3*d^2)
                !> Calculate rho*epsilon first and then divide by rho to ease save loops
                DO iElem=1,nElems; DO k=0,ZDIM(NCalc);  DO j=0,NCalc; DO i=0,NCalc;
                    ! Add epsilon_s contribution
                    DO, a=1,3
                        EpsilonFluc(i,j,k,iElem) = EpsilonFluc(i,j,k,iElem) + mu0*(Vorticity(a,i,j,k,iElem) * &
                                                                                   Vorticity(a,i,j,k,iElem))
                    END DO
                    ! Add epsilon_d contribution
                        EpsilonFluc(i,j,k,iElem) = EpsilonFluc(i,j,k,iElem) + 4./3.*mu0*(dFluc(i,j,k,iElem) * &
                                                                                         dFluc(i,j,k,iElem))
                    ! Divide by rho to get epsilon
                        EpsilonFluc(i,j,k,iElem) = EpsilonFluc(i,j,k,iElem) / UMean(1,i,j,k,iElem)
                END DO; END DO; END DO; END DO

                ! Average epsilon
                !<<< Welford's algorithm for variance
                !    (count, mean, M2) = existingAggregate
                !    count = count + 1
                !    delta = newValue - mean
                !    mean = mean + delta / count
                !    delta2 = newValue - mean
                !    M2 = M2 + delta * delta2
                EpsilonTmp = EpsilonFluc - EpsilonSum
                EpsilonSum = EpsilonSum  + EpsilonTmp / (iArg-2)            ! prm-file, timeAvg-file, state-file(s)

                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING STATE VALUES DONE'
                SWRITE(UNIT_stdOut,'(132("-"))')
            END IF !/TimeAvg or. State-File

            ! Postpone output until we are done with all state files
            IF (iArg.EQ.nArgs) THEN
                SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT...'
                ! Allocate output arrays
                ALLOCATE(USolution(nVarTurb,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))

                ! Set remaining varname
                ALLOCATE(varNames_loc(nVarTurb))
                varnames_loc(1:5) = varnames_gen(1:5)
                varnames_loc(6)   = 'TKE'
                varnames_loc(7)   = 'epsilon'

                ! Fill the output array
                USolution(1:5,:,:,:,:) = UMean(:,:,:,:,:)
                USolution(6  ,:,:,:,:) = TKE
                USolution(7  ,:,:,:,:) = EpsilonFin

                ! Call output routine
                CALL WriteStateFile(MeshFileName=TRIM(MeshFile),SolutionArray=USolution,ArrayName='DG_Solution')

                SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT DONE'
            END IF

!> =================================================================================================================================
!> Dissipation calculation according to Pope (2015)
!> - Comparable to incompressible turbulent dissipation (eps1) according to Hillewaert (2013)
!> =================================================================================================================================
        CASE(3)
            ! This mode needs at least a time-average and a state-file
            IF ((nArgs-skipArgs).LT.2) &
                CALL CollectiveStop(__STAMP__, 'At least a timeAvg and a state files are requirede for Mode=3!')

            ! Conservative variables plus k and epsilon
            nVarTurb = 7

            ! We expect one timeAvg file and several time-accurate files. Start by filling the mean array ONCE
            IF (iArg.EQ.(1+skipArgs)) THEN
                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING MEAN VALUES...'
                ! Get TKE so we can overwrite with mean values later
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='TKE')

                ! Get mean values
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='Mean')

                ! Can't have dissipation without viscosity
                IF (ALMOSTEQUAL(mu0,0.)) &
                    CALL CollectiveStop(__STAMP__, 'No viscosity given. Please set mu0 in your parameter file.')

                ! Allocate and fill array to the mean gradients
                ALLOCATE(GradUMeanx (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(GradUMeany (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(GradUMeanz (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                GradUMeanx = GradUx; GradUMeany = GradUy; GradUMeanz = GradUz;

                ! Allocate arrays ONCE to hold epsilon
                ALLOCATE(EpsilonFluc(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonTmp (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonSum (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonMean(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                EpsilonSum = 0.
                ALLOCATE(EpsilonFin   (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                EpsilonFin   = 0.

                ! Calculate dissipation of mean flow
                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
                    ! Density
                    rho            = Umean(1,i,j,k,iElem)

                    ! Get gradients
                    ! Get fluctuation of gradients. Linear operation, so we can operate on gradients
                    GradVel(1:3,1) = gradUMeanx(2:4,i,j,k,iElem)
                    GradVel(1:3,2) = gradUMeany(2:4,i,j,k,iElem)
                    GradVel(1:3,3) = gradUMeanz(2:4,i,j,k,iElem)

                    ! compute tensor of velocity gradients
                    S              = 0.5*(Gradvel+TRANSPOSE(GradVel))

                    !< epsilon = 2*nu*<s_ij s_ij>
                    eps1        = 0.
                    DO p = 1,3; DO q = 1,3
                        eps1   = eps1 + 2*mu0/rho * S(p,q)*S(p,q)
                    END DO; END DO

                    EpsilonMean(i,j,k,iElem) = eps1

                END DO; END DO; END DO; END DO

                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING MEAN VALUES DONE'
                SWRITE(UNIT_stdOut,'(132("-"))')
            ELSE
                ! Now get the time-accurate values from the remaining state files
                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING STATE VALUES...'
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='DG_Solution')

                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
                    ! Density
                    rho            = U(1,i,j,k,iElem)

                    ! Get fluctuation of gradients. Linear operation, so we can operate on gradients
                    !--> Use Eps0 = EpsMean + EpsFluc instead
                    GradVel(1:3,1) = GradUx(2:4,i,j,k,iElem) !- gradUMeanx(2:4,i,j,k,iElem)
                    GradVel(1:3,2) = GradUy(2:4,i,j,k,iElem) !- gradUMeany(2:4,i,j,k,iElem)
                    GradVel(1:3,3) = GradUz(2:4,i,j,k,iElem) !- gradUMeanz(2:4,i,j,k,iElem)

                    ! compute tensor of velocity gradients
                    S              = 0.5*(Gradvel+TRANSPOSE(GradVel))

                    !< epsilon = 2*nu*<s_ij s_ij>
                    eps1        = 0.
                    DO p = 1,3; DO q = 1,3
                        eps1   = eps1 + 2*mu0/rho * S(p,q)*S(p,q)
                    END DO; END DO

                    EpsilonFluc(i,j,k,iElem) = eps1

                END DO; END DO; END DO; END DO

                ! Average epsilon
                !<<< Welford's algorithm for variance
                !    (count, mean, M2) = existingAggregate
                !    count = count + 1
                !    delta = newValue - mean
                !    mean = mean + delta / count
                !    delta2 = newValue - mean
                !    M2 = M2 + delta * delta2
                EpsilonTmp = EpsilonFluc - EpsilonSum
                EpsilonSum = EpsilonSum  + EpsilonTmp / (iArg-2)            ! prm-file, timeAvg-file, state-file(s)

                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING STATE VALUES DONE'
                SWRITE(UNIT_stdOut,'(132("-"))')
            END IF

            ! Postpone output until we are done with all state files
            IF (iArg.EQ.nArgs) THEN
                SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT...'
                ! Allocate output arrays
                ALLOCATE(USolution(nVarTurb,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))

                !--- We might have areas with laminar flow. Set epsilon to zero in the absence of TKE
                !--> stricly only required for omega, so leave epsilon alone
                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
!                    IF(ALMOSTZERO(TKE(i,j,k,iElem))) THEN
!                        EpsilonFin(i,j,k,iElem) = 0.
!                    ELSE
                        EpsilonFin(i,j,k,iElem) = EpsilonSum(i,j,k,iElem) - EpsilonMean(i,j,k,iElem)
!                    END IF
                END DO; END DO; END DO; END DO

                ! Fill the output array
                USolution(1:5,:,:,:,:)     = UMean(:,:,:,:,:)
                USolution(6  ,:,:,:,:)     = TKE  (  :,:,:,:)

                ! Exact Solution eps1, mu0 = 0.03, rho =  2
!                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
!                    USolution(6  ,i,j,k,iElem) = (27.*PP_Pi**2.*cos(PP_Pi*(Elem_xGP(1,i,j,k,iElem) + Elem_xGP(2,i,j,k,iElem) + &
!                                                                         Elem_xGP(3,i,j,k,iElem)))**2.)/40000
!                END DO; END DO; END DO; END DO

                USolution(7  ,:,:,:,:)     = EpsilonFin

                ! Set remaining varname
                ALLOCATE(varNames_loc(nVarTurb))
                varnames_loc(1:5) = varnames_gen(1:5)
                varnames_loc(6)   = 'TKE'
                varnames_loc(7)   = 'epsilon'

                ! Call output routine
                CALL WriteStateFile(MeshFileName=TRIM(MeshFile),SolutionArray=USolution,ArrayName='DG_Solution')

                SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT DONE'
            END IF

!> =================================================================================================================================
!> Dissipation calculation according to Hillewaert (2013)
!> =================================================================================================================================
        CASE(4)
            ! Conservative variables plus k and epsilon
            nVarTurb = 7

            ! We expect one timeAvg file and several time-accurate files. Start by filling the mean array ONCE
            IF (iArg.EQ.(1+skipArgs)) THEN
                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING MEAN VALUES...'
                ! Get TKE so we can overwrite with mean values later
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='TKE')

                ! Get mean values
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='Mean')

                ! Can't have dissipation without viscosity
                IF (ALMOSTEQUAL(mu0,0.)) &
                    CALL CollectiveStop(__STAMP__, 'No viscosity given. Please set mu0 in your parameter file.')

                ! Allocate and fill array to the mean gradients
                ALLOCATE(GradUMeanx (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(GradUMeany (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(GradUMeanz (PP_nVar,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                GradUMeanx = GradUx; GradUMeany = GradUy; GradUMeanz = GradUz;

                ! Allocate arrays ONCE to hold epsilon
                ALLOCATE(EpsilonFluc(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonTmp (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonSum (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                ALLOCATE(EpsilonMean(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                EpsilonSum = 0.
                ALLOCATE(EpsilonFin   (0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))
                EpsilonFin   = 0.

                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING MEAN VALUES DONE'
                SWRITE(UNIT_stdOut,'(132("-"))')
            ELSE
                ! Now get the time-accurate values from the remaining state files
                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING STATE VALUES...'
                CALL ReadStateFile(Parameterfile,Statefile,ArrayName='DG_Solution')

                ! Follow procedure for Taylor-Green-Vortex
                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
                    ! Density
                    rho         = U(1,i,j,k,iElem)

                    ! Compute primitive values (of u,v,w) at each GP
                    Vel(1:3)    = U(2:4,i,j,k,iElem)/U(1,i,j,k,iElem)

                    ! Get fluctuation of gradients. Linear operation, so we can operate on gradients
                    GradVel(1:3,1) = GradUx(2:4,i,j,k,iElem) - gradUMeanx(2:4,i,j,k,iElem)
                    GradVel(1:3,2) = GradUy(2:4,i,j,k,iElem) - gradUMeany(2:4,i,j,k,iElem)
                    GradVel(1:3,3) = GradUz(2:4,i,j,k,iElem) - gradUMeanz(2:4,i,j,k,iElem)

                    ! Pressure (WITH MEAN VELOCITY?)
                    Pressure    = KappaM1*(U(5,i,j,k,iElem)-0.5*SUM(U(2:4,i,j,k,iElem)*Vel(1:3)))

                    ! compute divergence of velocity
                    divU        = GradVel(1,1)+GradVel(2,2)+GradVel(3,3)

                    ! compute tensor of velocity gradients
                    S           = 0.5*(Gradvel+TRANSPOSE(GradVel))

                    ! deviatoric part of strain tensor
                    Sd=S
                    DO p=1,3
                        Sd(p,p) = Sd(p,p)-1./3.*divU
                    END DO

                    ! compute integrand for epsilon2, viscous distribution to dissipation
!                    eps2        = divU*divU

                    ! Matrix : Matrix product for velocity gradient tensor, S:S and Sd:Sd
                    u_tens = 0.; s_tens = 0.; sd_tens = 0.
                    eps1   = 0.; eps2   = 0.; eps3    = 0.

                    ! compute integrand for epsilon3, pressure contribution to dissipation (compressibility effect)
                    eps3        = Pressure*divU

                    DO p = 1,3
                        DO q = 1,3
                            ! dissipation rate epsilon incomp from velocity gradient tensor (Diss Fauconnier)
                            u_tens = u_tens+GradVel(p,q)*GradVel(p,q)
                            ! dissipation rate epsilon incomp from strain rate tensor (incomp) (Sagaut)
                            s_tens = s_tens+S(p,q)*S(p,q)
                            ! dissipation rate epsilon 1 from deviatoric part of strain rate tensor Sd (compressible)
                            eps1   = eps1  +Sd(p,q)*Sd(p,q)
                            ! dissipation rate epsilon 3 from pressure times div u (compressible)
                            ! eps3
                        END DO
                    END DO

                    ! Scale kinetic energy dissipation
                    u_tens      = mu0       * u_tens
                    s_tens      = 2.*mu0/rho* s_tens
                    eps1        = 2.*mu0/rho* eps1
!                    eps2        = muv/rho*eps2
                    eps3        = -1./rho   * eps3

                    ! Epsilon
                    !> epsilon = eps1 + eps2 + eps3
                    EpsilonFluc(i,j,k,iElem) = eps1 + eps2 + eps3
                END DO; END DO; END DO; END DO

                ! Average epsilon
                !<<< Welford's algorithm for variance
                !    (count, mean, M2) = existingAggregate
                !    count = count + 1
                !    delta = newValue - mean
                !    mean = mean + delta / count
                !    delta2 = newValue - mean
                !    M2 = M2 + delta * delta2
                EpsilonTmp = EpsilonFluc - EpsilonSum
                EpsilonSum = EpsilonSum  + EpsilonTmp / (iArg-2)            ! prm-file, timeAvg-file, state-file(s)

                SWRITE(UNIT_stdOut,'(A)') ' PROCESSING STATE VALUES DONE'
                SWRITE(UNIT_stdOut,'(132("-"))')
            END IF

            ! Postpone output until we are done with all state files
            IF (iArg.EQ.nArgs) THEN
                SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT...'
                ! Allocate output arrays
                ALLOCATE(USolution(nVarTurb,0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems))

                !--- We might have areas with laminar flow. Set epsilon to zero in the absence of TKE
                !--> stricly only required for omega, so leave epsilon alone
                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
!                    IF(ALMOSTZERO(TKE(i,j,k,iElem))) THEN
!                        EpsilonFin(i,j,k,iElem) = 0.
!                    ELSE
                        EpsilonFin(i,j,k,iElem) = EpsilonSum(i,j,k,iElem)
!                    END IF
                END DO; END DO; END DO; END DO

                ! Fill the output array
                USolution(1:5,:,:,:,:)     = UMean(:,:,:,:,:)
                USolution(6  ,:,:,:,:)     = TKE  (  :,:,:,:)

                ! Exact Solution eps1
!                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
!                    USolution(6  ,i,j,k,iElem) = (9.*PP_Pi**2.*cos(PP_Pi*(Elem_xGP(1,i,j,k,iElem) + Elem_xGP(2,i,j,k,iElem) + &
!                                                                         Elem_xGP(3,i,j,k,iElem)))**2.)/20000
!                END DO; END DO; END DO; END DO
                ! Exact Solution eps3
!                DO iElem=1,nElems; DO k=0,ZDIM(nCalc);  DO j=0,nCalc; DO i=0,nCalc
!                    USolution(6  ,i,j,k,iElem) = (3.*PP_Pi*cos(PP_Pi*SUM(Elem_xGP(1:3,i,j,k,iElem)))* &
!                                                ((3.*     (sin(PP_Pi*SUM(Elem_xGP(1:3,i,j,k,iElem)))/10. + 2.)*&
!                                                          (sin(PP_Pi*SUM(Elem_xGP(1:3,i,j,k,iElem)))/20. + 1.))/5. - 8./5.))/40.
!                END DO; END DO; END DO; END DO

                USolution(7  ,:,:,:,:) = EpsilonFin

                ! Set remaining varname
                ALLOCATE(varNames_loc(nVarTurb))
                varnames_loc(1:5) = varnames_gen(1:5)
                varnames_loc(6)   = 'TKE'
                varnames_loc(7)   = 'epsilon'

                ! Call output routine
                CALL WriteStateFile(MeshFileName=TRIM(MeshFile),SolutionArray=USolution,ArrayName='DG_Solution')

                SWRITE(UNIT_stdOut,'(A)') ' CALLING OUTPUT DONE'
            END IF

    END SELECT

END DO !iArg

! Measure tool duration
Time=FLEXITIME()
! Write final output
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' CALCTURBULENCE FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

SDEALLOCATE(Args)

CALL FinalizeCalcTurb()

! we also have to finalize MPI itself here
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif


END PROGRAM Posti_CalcTurb
