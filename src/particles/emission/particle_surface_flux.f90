!=================================================================================================================================
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
!=================================================================================================================================
#include "flexi.h"
#include "particle.h"

!===================================================================================================================================
! module for particle emission
!===================================================================================================================================
MODULE MOD_Surface_Flux
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitializeParticleSurfaceflux
  MODULE PROCEDURE InitializeParticleSurfaceflux
END INTERFACE

INTERFACE ParticleSurfaceflux
  MODULE PROCEDURE ParticleSurfaceflux
END INTERFACE

PUBLIC :: InitializeParticleSurfaceflux
PUBLIC :: ParticleSurfaceflux
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
SUBROUTINE ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_Particle_Boundary_Vars  ,ONLY: nPartBound
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF,BezierSampleN,TriaSurfaceFlux
USE MOD_Particle_Vars           ,ONLY: nSpecies,Species,DoPoissonRounding,DoTimeDepInflow
USE MOD_Particle_Vars           ,ONLY: UseCircularInflow
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT) :: MaxSurfacefluxBCs, nDataBC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(42)          :: tmpStr,tmpStr2,tmpStr3
INTEGER                :: iSpec,iSF
!===================================================================================================================================

DO iSpec = 1,nSpecies
  ! Find number of SurfaceFluxBC per species
  WRITE(UNIT=tmpStr,FMT='(I0)') iSpec
  Species(iSpec)%nSurfacefluxBCs = GETINT('Part-Species'//TRIM(tmpStr)//'-nSurfacefluxBCs','0')

  IF (Species(iSpec)%nSurfacefluxBCs.EQ.0) THEN
    CYCLE
  ELSE
    ! Initialize SurfaceFlux to BC mapping
    ALLOCATE(Species(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
    Species(iSpec)%Surfaceflux(:)%BC = -1

    DO iSF = 1,Species(iSpec)%nSurfacefluxBCs
      WRITE(UNIT=tmpStr2,FMT='(I0)') iSF
      tmpStr2 = TRIM(tmpStr)//'-Surfaceflux'//TRIM(tmpStr2)
      Species(iSpec)%Surfaceflux(iSF)%BC = GETINT('Part-Species'//TRIM(tmpStr2)//'-BC','0')
    END DO
  END IF

  ! Increase global counter of SurfaceFluxBC
  MaxSurfacefluxBCs = MAX(MaxSurfacefluxBCs,Species(iSpec)%nSurfacefluxBCs)

  DO iSF = 1,Species(iSpec)%nSurfacefluxBCs
    WRITE(UNIT=tmpStr2,FMT='(I0)') iSF
    tmpStr2 = TRIM(tmpStr)//'-Surfaceflux'//TRIM(tmpStr2)
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticle        = 0
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus = 0
    Species(iSpec)%Surfaceflux(iSF)%VFR_total               = 0
    Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal = 0

    ! Get SurfaceFlux data
    IF ((Species(iSpec)%Surfaceflux(iSF)%BC.LT.1) .OR. Species(iSpec)%Surfaceflux(iSF)%BC.GT.nPartBound) THEN
      CALL abort(__STAMP__, 'SurfacefluxBCs must be between 1 and nPartBound!')

    ! SurfaceFlux data not yet set
    ELSE IF (BCdata_auxSF(Species(iSpec)%Surfaceflux(iSF)%BC)%SideNumber.EQ.-1) THEN
      BCdata_auxSF(Species(iSpec)%Surfaceflux(iSF)%BC)%SideNumber = 0
      nDataBC = nDataBC+1
    END IF

    ! Get SurfaceFlux velocity distribution
    Species(iSpec)%Surfaceflux(iSF)%velocityDistribution  = &
        TRIM(GETSTR('Part-Species'//TRIM(tmpStr2)//'-velocityDistribution','constant'))
    IF (TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'constant') &
      CALL abort(__STAMP__, 'Only constant velocity distribution implemented for SurfaceFlux!')

    Species(iSpec)%Surfaceflux(iSF)%VeloIC       = GETREAL(   'Part-Species'//TRIM(tmpStr2)//'-VeloIC'      )
    Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal = GETLOGICAL('Part-Species'//TRIM(tmpStr2)//'-VeloIsNormal')

    ! Get additional SurfaceFlux velocity information
    IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      Species(iSpec)%Surfaceflux(iSF)%CircularInflow      = .FALSE.
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%VeloVecIC           = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-VeloVecIC'          ,3)
      Species(iSpec)%Surfaceflux(iSF)%CircularInflow      = GETLOGICAL(  'Part-Species'//TRIM(tmpStr2)//'-CircularInflow'       )

      !
      IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        UseCircularInflow = .TRUE.
        Species(iSpec)%Surfaceflux(iSF)%dir(1)   = GETINT('Part-Species'//TRIM(tmpStr2)//'-axialDir')

        IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.1) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2) = 2
          Species(iSpec)%Surfaceflux(iSF)%dir(3) = 3
        ELSE IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.2) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2) = 3
          Species(iSpec)%Surfaceflux(iSF)%dir(3) = 1
        ELSE IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.3) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2) = 1
          Species(iSpec)%Surfaceflux(iSF)%dir(3) = 2
        ELSE
          CALL ABORT(__STAMP__,'ERROR in init: axialDir for SFradial must be between 1 and 3!')
        END IF

        IF ( Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(Species(iSpec)%Surfaceflux(iSF)%dir(2)).NE.0. .OR. &
             Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(Species(iSpec)%Surfaceflux(iSF)%dir(3)).NE.0. ) &
          CALL abort(__STAMP__,'ERROR in init: axialDir for SFradial do not correspond to VeloVecIC!')

        Species(iSpec)%Surfaceflux(iSF)%origin   = GETREALARRAY('Part-Species'//TRIM(tmpStr2)//'-origin',2,'0. , 0.')
        WRITE(UNIT=tmpStr3,FMT='(E16.8)') HUGE(Species(iSpec)%Surfaceflux(iSF)%rmax)
        Species(iSpec)%Surfaceflux(iSF)%rmax     = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-rmax',TRIM(tmpStr3))
        Species(iSpec)%Surfaceflux(iSF)%rmin     = GETREAL(     'Part-Species'//TRIM(tmpStr2)//'-rmin','0.')
      END IF
    END IF !.NOT.VeloIsNormal

    !--- normalize VeloVecIC
    IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      IF (.NOT. ALL(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(:).EQ.0.)) THEN
        Species(iSpec)%Surfaceflux(iSF)%VeloVecIC =                   Species(iSpec)%Surfaceflux(iSF)%VeloVecIC  &
                                                    /SQRT(DOT_PRODUCT(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC, &
                                                                      Species(iSpec)%Surfaceflux(iSF)%VeloVecIC))
      END IF
    END IF

    Species(iSpec)%Surfaceflux(iSF)%PartDensity  = GETREAL(   'Part-Species'//TRIM(tmpStr2)//'-PartDensity'  ,'0.'     )
    Species(iSpec)%Surfaceflux(iSF)%ReduceNoise  = GETLOGICAL('Part-Species'//TRIM(tmpStr2)//'-ReduceNoise'  ,'.FALSE.')

    IF (DoPoissonRounding .AND. Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
      SWRITE(UNIT_StdOut,'(A)') ' WARNING: Poisson sampling not possible for noise reduction of surfacefluxes:'
      SWRITE(UNIT_StdOut,'(A)') ' switching now to Random rounding...'
      DoPoissonRounding   = .FALSE.
    END IF
    IF (DoTimeDepInflow   .AND. Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
      SWRITE(UNIT_StdOut,'(A)') ' WARNING: Time-dependent inflow is not possible for noise reduction of surfacefluxes:'
      SWRITE(UNIT_StdOut,'(A)') ' switching now to Random rounding...'
      DoTimeDepInflow   = .FALSE.
    END IF

    ! TriaSurfaceFlux
    IF (TriaSurfaceFlux) THEN
      Species(iSpec)%Surfaceflux(iSF)%AcceptReject = .FALSE.
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%AcceptReject = GETLOGICAL('Part-Species'//TRIM(tmpStr2)//'-AcceptReject','.TRUE.')
    END IF

    IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject .AND. BezierSampleN.GT.1) THEN
      SWRITE(UNIT_StdOut,'(A)') ' WARNING: BezierSampleN > 0 may not be necessary as ARM is used for SurfaceFlux!'
    ELSE IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%AcceptReject .AND. BezierSampleN.LE.NGeo .AND. .NOT.TriaSurfaceFlux) THEN
      SWRITE(UNIT_StdOut,'(A)') ' WARNING: The choosen small BezierSampleN (def.: NGeo) might result in inhom. SurfFluxes without ARM!'
    END IF

    IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
      ! 1 for linear elements, this is an arbitrary estimation for higher N!
      WRITE( tmpStr3, '(I0.2)') NGeo*NGeo*NGeo
      Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = GETINT('Part-Species'//TRIM(tmpStr2)//'-ARM_DmaxSampleN',tmpStr3)
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = 0
    END IF

  END DO ! iSF
END DO ! iSpec

END SUBROUTINE ReadInAndPrepareSurfaceFlux


!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
SUBROUTINE BCSurfMeshSideAreasandNormals()
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize,SurfMeshSubSideData,BezierSampleN,SurfMeshSideAreas,TriaSurfaceFlux
USE MOD_Mesh_Vars              ,ONLY: nBCSides, offsetElem, SideToElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas, CalcNormAndTangTriangle
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: totalArea
INTEGER                 :: BCSideID, ElemID, iLocSide, SideID, jSample, iSample
REAL                    :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_nOut(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_t1(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_t2(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
!===================================================================================================================================

! Calculate total area of SurfaceFlux BCs
totalArea = 0.

DO BCSideID = 1,nBCSides
  ! SideToElem is coming from DG and available for local elems
  ElemID = SideToElem(1,BCSideID)

  ! mortar elements MIGHT need to be treated differently. Take the side pointing back
  IF (ElemID.LT.1) THEN
    ElemID   = SideToElem(2,BCSideID)
    iLocSide = SideToElem(4,BCSideID)
  ELSE
    iLocSide = SideToElem(3,BCSideID)
  END IF

  ! Get global SideID from local SideID
  SideID = GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)

  !
  IF (TriaSurfaceFlux) THEN
    IF (SurfFluxSideSize(1).NE.1 .OR. SurfFluxSideSize(2).NE.2) &
      CALL ABORT(__STAMP__, 'SurfFluxSideSize must be 1,2 for TriaSurfaceFlux!')

    ! Calculate surfSideArea
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      CALL CalcNormAndTangTriangle( SideID = SideID                              &
                                  , nVec   = tmp_Vec_nOut  (:,iSample,jSample)   &
                                  , tang1  = tmp_Vec_t1    (:,iSample,jSample)   &
                                  , tang2  = tmp_Vec_t2    (:,iSample,jSample)   &
                                  , area   = tmp_SubSideAreas(iSample,jSample)   &
                                  , TriNum = jSample)

      SurfMeshSideAreas(BCSideID) = SurfMeshSideAreas(BCSideID) + tmp_SubSideAreas(iSample,jSample)
    END DO; END DO

  ELSE
    IF (ANY(SurfFluxSideSize.NE.BezierSampleN)) &
      CALL ABORT(__STAMP__, 'SurfFluxSideSize must be BezierSampleN,BezierSampleN for .NOT.TriaSurfaceFlux!')

    CALL GetBezierSampledAreas( SideID                      = SideID                        &
                              , BezierSampleN               = BezierSampleN                 &
                              , SurfMeshSubSideAreas        = tmp_SubSideAreas              &
                              , SurfMeshSideArea_opt        = SurfMeshSideAreas(BCSideID)   &
                              , SurfMeshSubSideVec_nOut_opt = tmp_Vec_nOut                  &
                              , SurfMeshSubSideVec_t1_opt   = tmp_Vec_t1                    &
                              , SurfMeshSubSideVec_t2_opt   = tmp_Vec_t2)
  END IF

  ! Sum up total area
  totalArea = totalArea + SurfMeshSideAreas(BCSideID)

  ! Store subside data
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn = -tmp_Vec_nOut(  :,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1  =  tmp_Vec_t1(    :,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2  =  tmp_Vec_t2(    :,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%area    =  tmp_SubSideAreas(iSample,jSample)
  END DO; END DO
END DO

#if CODE_ANALYZE
IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
IPWRITE(*,*)" totalArea       = ",totalArea
IPWRITE(*,*)" totalArea/(pi) = ",totalArea/(ACOS(-1.))
IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
#endif /*CODE_ANALYZE*/

END SUBROUTINE BCSurfMeshSideAreasandNormals


!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux(nDataBC)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nBCSides,offsetElem,BC,SideToElem
USE MOD_Particle_Boundary_Vars ,ONLY: nPartBound
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, SurfMeshSubSideData,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces_Vars ,ONLY: SideType
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow,Species,DoSurfaceFlux,nSpecies
USE MOD_ReadInTools
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                           :: nDataBC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: TmpMapToBC(1:nDataBC),TmpSideStart(1:nDataBC),TmpSideNumber(1:nDataBC),TmpSideEnd(1:nDataBC)
INTEGER               :: TmpSideNext(1:nBCSides)
INTEGER               :: countDataBC,iBC,BCSideID,currentBC,iSF,ElemID,iCount,iLocSide,SideID,iPartBound
INTEGER               :: iSample, jSample, iSpec
REAL, ALLOCATABLE     :: areasLoc(:),areasGlob(:)
LOGICAL               :: OutputSurfaceFluxLinked
!===================================================================================================================================

!-- 2.: create Side lists for applicable BCs
!--- 2a: temporary (linked) lists
OutputSurfaceFluxLinked = GETLOGICAL('OutputSurfaceFluxLinked','.FALSE.')

TmpMapToBC    = 0
TmpSideStart  = 0
TmpSideNumber = 0
TmpSideEnd    = 0
TmpSideNext   = 0

countDataBC   = 0

! Loop over all particle boundaries and add to list
DO iBC=1,nPartBound
  ! Side is not a SurfaceFlux side
  IF (BCdata_auxSF(iBC)%SideNumber.EQ. -1) CYCLE

  countDataBC = countDataBC + 1
  TmpMapToBC(countDataBC) = iBC
END DO

! Loop over all local DG boundary sides and add SurfaceFlux sides to linked list
DO BCSideID = 1,nBCSides
  currentBC = 0

  ! compare current BCID against SF BCIDs
  DO iBC = 1,countDataBC
    IF (BC(BCSideID) .EQ. TmpMapToBC(iBC)) currentBC = iBC
  END DO

  ! current BCID is not a SF BCIDs
  IF (currentBC.EQ.0) CYCLE

  ! Start of Linked List for Sides
  IF (TmpSideNumber(currentBC).EQ.0) THEN
    TmpSideStart(currentBC) = BCSideID
  ! Next Side
  ELSE
    TmpSideNext(TmpSideEnd(currentBC)) = BCSideID
  END IF

  !-- prepare for next entry in list and increment number of sides
  TmpSideEnd(currentBC)    = BCSideID
  TmpSideNumber(currentBC) = TmpSideNumber(currentBC) + 1
END DO ! BCSideID

!--- 2b: save sequential lists in BCdata_auxSF
DO iBC = 1,countDataBC
  ! store number of sides for BCID (including non-SF BCIDs)
  BCdata_auxSF(TmpMapToBC(iBC))%SideNumber = TmpSideNumber(iBC)
  IF (TmpSideNumber(iBC).EQ.0) CYCLE

  ! allocate array to hold additional SF data
  ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%SideList(1:TmpSideNumber(iBC)))

  IF (TriaSurfaceFlux) THEN
    ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
    ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(1:TmpSideNumber(iBC)))
  END IF

  DO iSpec = 1,nSpecies
    DO iSF = 1,Species(iSpec)%nSurfacefluxBCs
      ! only surfacefluxes with iBC
      IF (TmpMapToBC(iBC).EQ.Species(iSpec)%Surfaceflux(iSF)%BC) THEN
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))

        IF (UseCircularInflow .AND. (iSF .LE. Species(iSpec)%nSurfacefluxBCs)) &
          ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(1:TmpSideNumber(iBC)))

      END IF
    END DO ! iSF
  END DO ! iSpec

  BCSideID = TmpSideStart(iBC)
  iCount   = 0

  ! follow BCSideID list seq. with iCount
  DO
    iCount = iCount+1
    BCdata_auxSF(TmpMapToBC(iBC))%SideList(iCount) = BCSideID

    IF (TriaSurfaceFlux) THEN
      ! SideToElem is coming from DG and available for local elems
      ElemID = SideToElem(1,BCSideID)

      ! mortar elements MIGHT need to be treated differently. Take the side pointing back
      IF (ElemID.LT.1) THEN
        ElemID   = SideToElem(2,BCSideID)
        iLocSide = SideToElem(4,BCSideID)
      ELSE
        iLocSide = SideToElem(3,BCSideID)
      END IF

      ! Get global SideID from local SideID
      SideID = GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)

      !check that all sides are planar if TriaSurfaceFlux is used for Tracing or RefMapping
      IF (.NOT.TriaTracking) THEN
        IF (SideType(SideID).NE.PLANAR_RECT .AND. SideType(SideID).NE.PLANAR_NONRECT) &
          CALL ABORT(__STAMP__,'Every surfaceflux-sides must be planar if TriaSurfaceFlux is used for Tracing or RefMapping!')
      END IF !.NOT.TriaTracking

      ! Store additional information for TriaSurfaceFlux
      DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
        CALL CalcNormAndTangTriangle( SideID   = SideID                                                                     &
                                    , midpoint = BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%midpoint &
                                    , ndist    = BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%ndist    &
                                    , xyzNod   = BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%xyzNod                   &
                                    , Vectors  = BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%Vectors                  &
                                    , TriNum   = jSample)
      END DO; END DO
    END IF ! TriaSurfaceFlux

    !-- BC-list specific data
    ! sum up total area
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      BCdata_auxSF(TmpMapToBC(iBC))%LocalArea = BCdata_auxSF(TmpMapToBC(iBC))%LocalArea &
        + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    END DO; END DO

    !-- next Side
    IF (BCSideID .EQ. TmpSideEnd(iBC)) THEN
      IF (TmpSideNumber(iBC).NE.iCount) THEN
        CALL ABORT(__STAMP__,'Someting is wrong with TmpSideNumber of iBC',iBC,999.)
      ELSE
        IF (OutputSurfaceFluxLinked) THEN
          IPWRITE(*,'(I4,I7,A53,I0)') iCount,' Sides have been found for Surfaceflux-linked PartBC ',TmpMapToBC(iBC)
        END IF

        DoSurfaceFlux = .TRUE.
        EXIT
      END IF
    END IF

    BCSideID = TmpSideNext(BCSideID)
  END DO ! BCSideID (iCount)
END DO ! iBC

!-- communicate areas
#if USE_MPI
ALLOCATE(areasLoc (1:nPartBound) &
        ,areasGlob(1:nPartBound))
areasLoc  = 0.
areasGlob = 0.

DO iPartBound = 1,nPartBound
  areasLoc(iPartBound) = BCdata_auxSF(iPartBound)%LocalArea
END DO

CALL MPI_ALLREDUCE(areasLoc,areasGlob,nPartBound,MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,IERROR)
#endif /*USE_MPI*/

DO iPartBound = 1,nPartBound
#if USE_MPI
  BCdata_auxSF(iPartBound)%GlobalArea = areasGlob(iPartBound)
#else
  BCdata_auxSF(iPartBound)%GlobalArea = BCdata_auxSF(iPartBound)%LocalArea
#endif /*USE_MPI*/
END DO

#if USE_MPI
DEALLOCATE(areasLoc,areasGlob)
#endif /*USE_MPI*/

END SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux


SUBROUTINE DefineCircInflowRejectType(iSpec, iSF, iSide)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars                 ,ONLY: offsetElem,SideToElem
USE MOD_Particle_Mesh_Tools       ,ONLY: GetGlobalNonUniqueSideID,GetSideBoundingBoxTria
USE MOD_Particle_Surfaces         ,ONLY: GetSideBoundingBox
USE MOD_Particle_Surfaces_Vars    ,ONLY: BCdata_auxSF
USE MOD_Particle_Vars             ,ONLY: Species,CountCircInflowType
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackingMethod
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF, iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: BoundingBox(1:3,1:8), origin(2), Vector1(3), Vector2(3), Vector3(3), xyzNod(3), VecBoundingBox(3)
REAL                  :: corner(3), corners(2,4), radiusCorner(2,4), rmax, rmin, point(2), vec(2)
LOGICAL               :: r0inside, intersecExists(2,2)
INTEGER               :: dir(3), iNode, currentBC, BCSideID, ElemID, iLocSide, SideID
!===================================================================================================================================

!-- check where the sides are located relative to rmax (based on corner nodes of bounding box)
!- RejectType=0 : complete side is inside valid bounds
!- RejectType=1 : complete side is outside of valid bounds
!- RejectType=2 : side is partly inside valid bounds
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
BCSideID  = BCdata_auxSF(currentBC)%SideList(iSide)
ElemID    = SideToElem(1,BCSideID)

! mortar elements MIGHT need to be treated differently. Take the side pointing back
IF (ElemID.LT.1) THEN
  ElemID   = SideToElem(2,BCSideID)
  iLocSide = SideToElem(4,BCSideID)
ELSE
  iLocSide = SideToElem(3,BCSideID)
END IF

! Get global SideID from local SideID
SideID = GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)

IF (TrackingMethod.EQ.TRIATRACKING) THEN
  CALL GetSideBoundingBoxTria(SideID,BoundingBox)
ELSE
  CALL GetSideBoundingBox(SideID,BoundingBox)
END IF

intersecExists = .FALSE.
r0inside       = .FALSE.
dir            = Species(iSpec)%Surfaceflux(iSF)%dir
origin         = Species(iSpec)%Surfaceflux(iSF)%origin
Vector1(:)     = 0.
Vector2(:)     = 0.
Vector3(:)     = 0.
xyzNod(1)      = MINVAL(BoundingBox(1,:))
xyzNod(2)      = MINVAL(BoundingBox(2,:))
xyzNod(3)      = MINVAL(BoundingBox(3,:))
VecBoundingBox(1) = MAXVAL(BoundingBox(1,:)) -MINVAL(BoundingBox(1,:))
VecBoundingBox(2) = MAXVAL(BoundingBox(2,:)) -MINVAL(BoundingBox(2,:))
VecBoundingBox(3) = MAXVAL(BoundingBox(3,:)) -MINVAL(BoundingBox(3,:))
Vector1(dir(2)) = VecBoundingBox(dir(2))
Vector2(dir(2)) = VecBoundingBox(dir(2))
Vector2(dir(3)) = VecBoundingBox(dir(3))
Vector3(dir(3)) = VecBoundingBox(dir(3))

!-- determine rmax (and corners)
DO iNode = 1,4
  SELECT CASE(iNode)
    CASE(1)
      corner = xyzNod
    CASE(2)
      corner = xyzNod + Vector1
    CASE(3)
      corner = xyzNod + Vector2
    CASE(4)
      corner = xyzNod + Vector3
  END SELECT

  corner(dir(2)) = corner(dir(2)) - origin(1)
  corner(dir(3)) = corner(dir(3)) - origin(2)
  corners(1:2,iNode) = (/corner(dir(2)),corner(dir(3))/) !coordinates of orth. dirs
  radiusCorner(1,iNode) = SQRT(corner(dir(2))**2+corner(dir(3))**2)
END DO !iNode

rmax=MAXVAL(radiusCorner(1,1:4))

!-- determine rmin
DO iNode = 1,4
  SELECT CASE(iNode)
    CASE(1)
      point = (/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
      vec   = (/Vector1(dir(2)),Vector1(dir(3))/)
    CASE(2)
      point = (/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
      vec   = (/Vector3(dir(2)),Vector3(dir(3))/)
    CASE(3)
      point = (/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
      vec   = (/-Vector1(dir(2)),-Vector1(dir(3))/)
    CASE(4)
      point = (/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
      vec   = (/-Vector3(dir(2)),-Vector3(dir(3))/)
  END SELECT

  vec = point + MIN(MAX(-DOT_PRODUCT(point,vec)/DOT_PRODUCT(vec,vec),0.),1.)*vec
  radiusCorner(2,iNode) = SQRT(DOT_PRODUCT(vec,vec)) !rmin
END DO !iNode

!-- determine if r0 is inside of bounding box
IF ((origin(1) .GE. MINVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(2),:))) .AND. &
    (origin(1) .LE. MAXVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(2),:))) .AND. &
    (origin(2) .GE. MINVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(3),:))) .AND. &
    (origin(2) .LE. MAXVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(3),:))) ) THEN
   r0inside = .TRUE.
END IF

IF (r0inside) THEN
  rmin = 0.
ELSE
  rmin = MINVAL(radiusCorner(2,1:4))
END IF

! define RejectType
IF ( (rmin .GT. Species(iSpec)%Surfaceflux(iSF)%rmax) .OR. (rmax .LT. Species(iSpec)%Surfaceflux(iSF)%rmin) ) THEN
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide) = 1

#if CODE_ANALYZE
  CountCircInflowType(2,iSF,iSpec)=CountCircInflowType(2,iSF,iSpec)+1
#endif

ELSE IF ( (rmax .LE. Species(iSpec)%Surfaceflux(iSF)%rmax) .AND. (rmin .GE. Species(iSpec)%Surfaceflux(iSF)%rmin) ) THEN
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide) = 0

#if CODE_ANALYZE
  CountCircInflowType(1,iSF,iSpec)=CountCircInflowType(1,iSF,iSpec)+1
#endif

ELSE
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide) = 2
#ifdef CODE_ANALYZE
  CountCircInflowType(3,iSF,iSpec)=CountCircInflowType(3,iSF,iSpec)+1
#endif
END IF !  (rmin > Surfaceflux-rmax) .OR. (rmax < Surfaceflux-rmin)

END SUBROUTINE DefineCircInflowRejectType


!===================================================================================================================================
! Init Particle Inserting via Surface Flux
!===================================================================================================================================
SUBROUTINE InitializeParticleSurfaceflux()
! Modules
USE MOD_Globals
USE MOD_HDF5_Input             ,ONLY: DatasetExists,ReadAttribute,ReadArray,GetDataSize
USE MOD_Mesh_Vars              ,ONLY: nBCSides,SideToElem,NGeo,offsetElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF,BezierSampleN,SurfMeshSubSideData,SurfMeshSideAreas
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas,GetSideBoundingBox,CalcNormAndTangTriangle
USE MOD_Particle_Vars          ,ONLY: Species,nSpecies,DoSurfaceFlux,DoPoissonRounding,DoTimeDepInflow
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow
USE MOD_ReadInTools
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartTime
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
INTEGER               :: iSpec,iSF,SideID,BCSideID,iSide,ElemID,iLocSide,iSample,jSample,currentBC
INTEGER               :: iCopy1, iCopy2, iCopy3
INTEGER               :: MaxSurfacefluxBCs
INTEGER               :: nDataBC                                                       ! number of different PartBounds used for SFs
REAL                  :: tmp_SubSideDmax(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL,ALLOCATABLE      :: tmp_BezierControlPoints2D(:,:,:,:,:)
REAL                  :: VFR_total
!===================================================================================================================================

ALLOCATE(SurfMeshSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:nBCSides), &
         SurfMeshSideAreas  (1:nBCSides))
SurfMeshSideAreas = 0.

! global calculations for sampling the faces for area and vector calculations (checks the integration with CODE_ANALYZE)
CALL BCSurfMeshSideAreasandNormals()

UseCircularInflow = .FALSE.
MaxSurfacefluxBCs = 0
DoSurfaceFlux     = .FALSE.

!-- 1.: read/prepare parameters and determine nec. BCs
CALL ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoPoissonRounding,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) ! set T if this is for all procs
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoTimeDepInflow  ,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) ! set T if this is for all procs
#endif /*USE_MPI*/

!
CALL CreateSideListAndFinalizeAreasSurfFlux(nDataBC)

#if CODE_ANALYZE
IF (UseCircularInflow) THEN
  ALLOCATE(CountCircInflowType(1:3,1:MaxSurfacefluxBCs,1:nSpecies))
  CountCircInflowType = 0
END IF
#endif /*CODE_ANALYZE*/

!-- 3.: initialize Surfaceflux-specific data
DO iSpec = 1,nSpecies
  DO iSF = 1,Species(iSpec)%nSurfacefluxBCs
    !--- 3a: SF-specific data of Sides
    ! go through sides if present in proc
    currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
      DO iSide = 1,BCdata_auxSF(currentBC)%SideNumber
        ! BCSideID is local SideID on current proc. SideToElem is coming from DG
        BCSideID = BCdata_auxSF(currentBC)%SideList(iSide)
        ElemID   = SideToElem(1,BCSideID)

        ! mortar elements MIGHT need to be treated differently. Take the side pointing back
        IF (ElemID.LT.1) THEN
          ElemID = SideToElem(2,BCSideID)
          iLocSide = SideToElem(4,BCSideID)
        ELSE
          iLocSide = SideToElem(3,BCSideID)
        END IF
        ! Get global SideID from local SideID
        SideID = GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)

        ! SurfaceFlux with Acceptance Rejection Method (ARM)
        IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
          CALL GetBezierSampledAreas( SideID                       = SideID                                             &
                                    , BezierSampleN                = BezierSampleN                                      &
                                    , BezierSurfFluxProjection_opt = .NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal  &
                                    , SurfMeshSubSideAreas         = tmp_SubSideAreas     &  !SubSide-areas proj. to inwards normals
                                    , DmaxSampleN_opt              = Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN    &
                                    , Dmax_opt                     = tmp_SubSideDmax                                    &
                                    , BezierControlPoints2D_opt    = tmp_BezierControlPoints2D)

        ELSE IF (.NOT.TriaSurfaceFlux) THEN
          CALL GetBezierSampledAreas( SideID                       = SideID                                             &
                                    , BezierSampleN                = BezierSampleN                                      &
                                    , BezierSurfFluxProjection_opt = .NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal  &
                                    , SurfMeshSubSideAreas         = tmp_SubSideAreas)       !SubSide-areas proj. to inwards normals

        ! TriaSurfaceFlux
        ELSE
          DO jSample = 1,SurfFluxSideSize(2); DO iSample = 1,SurfFluxSideSize(1)
            tmp_SubSideAreas(iSample,jSample) = SurfMeshSubSideData(iSample,jSample,BCSideID)%area
          END DO; END DO
        END IF

        ! Circular Inflow Init
        IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) CALL DefineCircInflowRejectType(iSpec, iSF, iSide)

        ! Init non-adaptive SF
        CALL InitNonAdaptiveSurfFlux(iSpec, iSF, iSide, tmp_SubSideAreas)

        ! Init Stuff for acceptance-rejection sampling on SF
        IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
          DO jSample = 1,SurfFluxSideSize(2); DO iSample = 1,SurfFluxSideSize(1)
            Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Dmax = tmp_SubSideDmax(iSample,jSample)
            IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
              ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
                                                                          ,iSide)%BezierControlPoints2D(1:2,0:NGeo,0:NGeo))
              DO iCopy1=0,NGeo; DO iCopy2=0,NGeo; DO iCopy3=1,2
                Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
                                                                   ,iSide)%BezierControlPoints2D(iCopy3,iCopy2,iCopy1) &
                  = tmp_BezierControlPoints2D(iCopy3,iCopy2,iCopy1,iSample,jSample)
              END DO; END DO; END DO
            END IF !.NOT.VeloIsNormal
          END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
        END IF
      END DO ! iSide

    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL ABORT(__STAMP__,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF

#if CODE_ANALYZE
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0 .AND. Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      IPWRITE(*,'(I4,A,2(x,I0),A,3(x,I0))') ' For Surfaceflux/Spec',iSF,iSpec,' are nType0,1,2: ' &
                                            , CountCircInflowType(1,iSF,iSpec),CountCircInflowType(2, iSF,iSpec) &
                                            , CountCircInflowType(3, iSF,iSpec)
    END IF
#endif /*CODE_ANALYZE*/
  END DO ! iSF
END DO ! iSpec

#if CODE_ANALYZE
SDEALLOCATE(CountCircInflowType)
#endif

! Setting variables required after a restart
IF (DoRestart) THEN
  DO iSpec = 1,nSpecies
    DO iSF = 1,Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iSF)%InsertedParticle = INT(Species(iSpec)%Init(iSF)%ParticleEmission * RestartTime,8)
    END DO

    DO iSF = 1, Species(iSpec)%nSurfacefluxBCs
      ! proc global total (for non-root: dummy!)
      IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal
      ! proc local total
      ELSE
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total
      END IF

      Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity*RestartTime / VFR_total,8)
    END DO
  END DO
END IF

#if USE_MPI
!set DoSurfaceFlux=T if at least 1 proc have SFs
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSurfaceFlux,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError)
#endif  /*USE_MPI*/

!-- no SFs defined
IF (.NOT.DoSurfaceFlux) THEN
  SWRITE(*,*) 'WARNING: No Sides for SurfacefluxBCs found! DoSurfaceFlux is now disabled!'
END IF

END SUBROUTINE InitializeParticleSurfaceflux


SUBROUTINE CalcPartInsSubSidesStandardCase(iSpec, iSF, PartInsSubSides)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Part_Emission_Tools     ,ONLY: IntegerDivide, SamplePoissonDistri
USE MOD_Particle_Globals        ,ONLY: ALMOSTEQUAL
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, BCdata_auxSF
USE MOD_Particle_TimeDisc_Vars  ,ONLY: RKdtFrac,RKdtFracTotal
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_TimeDisc_Vars           ,ONLY: t,dt
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: iSpec, iSF
INTEGER, INTENT(OUT), ALLOCATABLE   :: PartInsSubSides(:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)        :: inserted_Particle_iter,inserted_Particle_time,inserted_Particle_diff
INTEGER                :: currentBC, PartInsSF, IntSample
REAL                   :: VFR_total, PartIns, RandVal1
!===================================================================================================================================

currentBC = Species(iSpec)%Surfaceflux(iSF)%BC

! proc local total
VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total
PartIns   = Species(iSpec)%Surfaceflux(iSF)%PartDensity * dt*RKdtFrac * VFR_total
inserted_Particle_iter = INT(PartIns,8)

PartIns   = Species(iSpec)%Surfaceflux(iSF)%PartDensity * (t + dt*RKdtFracTotal) * VFR_total
!
!-- random-round the inserted_Particle_time for preventing periodicity
IF (inserted_Particle_iter.GE.1) THEN
  CALL RANDOM_NUMBER(RandVal1)
  inserted_Particle_time = INT(PartIns+RandVal1,8)

! needed, since InsertedParticleSurplus can increase and _iter>1 needs to be possible for preventing periodicity
ELSE IF (inserted_Particle_iter.GE.0) THEN
  ! dummy for procs without SFs (needed for mpi-comm, are cycled later)
  IF (ALMOSTEQUAL(PartIns,0.)) THEN
    inserted_Particle_time = INT(PartIns,8)
  ! poisson-distri of PartIns-INT(PartIns)
  ELSE
    CALL SamplePoissonDistri( PartIns-INT(PartIns), IntSample )
    ! INT(PartIns) + POISDISTRI( PartIns-INT(PartIns) )
    inserted_Particle_time = INT(INT(PartIns)+IntSample,8)
  END IF

! dummy for procs without SFs (needed for mpi-comm, are cycled later)
ELSE
  inserted_Particle_time = INT(PartIns,8)
END IF

!-- evaluate inserted_Particle_time and inserted_Particle_iter
inserted_Particle_diff = inserted_Particle_time - Species(iSpec)%Surfaceflux(iSF)%InsertedParticle &
                       - inserted_Particle_iter - Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus
Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0))
PartInsSF = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)
Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = Species(iSpec)%Surfaceflux(iSF)%InsertedParticle + INT(PartInsSF,8)

!-- calc global to-be-inserted number of parts and distribute to SubSides (proc local)
SDEALLOCATE(PartInsSubSides)
ALLOCATE(PartInsSubSides(SurfFluxSideSize(1),SurfFluxSideSize(2),1:BCdata_auxSF(currentBC)%SideNumber))
PartInsSubSides = 0

IF (BCdata_auxSF(currentBC)%SideNumber.LT.1) THEN
  IF (PartInsSF.NE.0) CALL ABORT(__STAMP__,'ERROR in ParticleSurfaceflux: Someting is wrong with PartInsSF of BC ',currentBC)
ELSE
  CALL IntegerDivide( PartInsSF                                                                                       &
                    , BCdata_auxSF(currentBC)%SideNumber*SurfFluxSideSize(1)*SurfFluxSideSize(2)                      &
                    , Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2) &
                                                                         ,1:BCdata_auxSF(currentBC)%SideNumber)%nVFR  &
                                                                         ,PartInsSubSides(1:SurfFluxSideSize(1)       &
                                                                         ,1:SurfFluxSideSize(2)                       &
                                                                         ,1:BCdata_auxSF(currentBC)%SideNumber) )
END IF

END SUBROUTINE CalcPartInsSubSidesStandardCase


!===================================================================================================================================
! Particle Inserting via Surface Flux and (if present) adaptiveBC (Surface Flux adapting part density, velocity or temperature)
!===================================================================================================================================
SUBROUTINE ParticleSurfaceflux()
! Modules
USE MOD_Globals
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Part_Emission_Tools     ,ONLY: IntegerDivide,SetParticleMass,SamplePoissonDistri
USE MOD_Part_Pos_and_Velo       ,ONLY: SetParticleVelocity
USE MOD_Part_Tools              ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Particle_Globals        ,ONLY: ALMOSTEQUAL
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF
USE MOD_Particle_Timedisc_Vars  ,ONLY: RKdtFrac
USE MOD_Particle_Vars           ,ONLY: Species,nSpecies,PDM,PEM
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos
USE MOD_Particle_Vars           ,ONLY: DoPoissonRounding,DoTimeDepInflow
USE MOD_Particle_Vars           ,ONLY: tSurfFluxLink
USE MOD_TimeDisc_Vars           ,ONLY: dt
#if CODE_ANALYZE
USE MOD_Part_Emission_Tools     ,ONLY: CalcVectorAdditionCoeffs
USE MOD_Particle_Tracking_Vars  ,ONLY: PartOut, MPIRankOut
#endif /*CODE_ANALYZE*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: nSurfacefluxPerElem
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime,LBElemSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: globElemID
INTEGER                     :: iSpec,iSF,iSide,currentBC,SideID
INTEGER                     :: NbrOfParticle
INTEGER                     :: BCSideID,ElemID,iLocSide,iSample,jSample,PartInsSubSide,iPart,iPartTotal
INTEGER                     :: ParticleIndexNbr
REAL                        :: Particle_pos(3),RandVal1,xyzNod(3)
REAL,ALLOCATABLE            :: particle_positions(:),particle_xis(:)
INTEGER,ALLOCATABLE         :: PartInsSubSides(:,:,:)
REAL                        :: xi(2)
INTEGER                     :: nReject,allowedRejections
LOGICAL                     :: AcceptPos
! variables used for sampling of of energies and impulse of emitted particles from surfaces
INTEGER                     :: PartsEmitted
REAL                        :: Vector1(3),Vector2(3),ndist(3),midpoint(3)
INTEGER                     :: Node1,Node2
! prevent inconsistency with non-triatracking by bilinear-routine (tol. might be increased)
REAL,PARAMETER              :: eps_nontria = 1.E-6
#if USE_LOADBALANCE
! load balance
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
#if CODE_ANALYZE
REAL                        :: tmpVec(3)
#endif /*CODE_ANALYZE*/
!===================================================================================================================================

DO iSpec = 1,nSpecies
  DO iSF = 1,Species(iSpec)%nSurfacefluxBCs
    PartsEmitted  = 0
    currentBC     = Species(iSpec)%Surfaceflux(iSF)%BC
    ! calculated within (sub)side-Loops!
    NbrOfParticle = 0
    iPartTotal    = 0

    ! Calc Particles for insertion in standard case
    IF ((.NOT.DoPoissonRounding).AND.(.NOT. DoTimeDepInflow)) CALL CalcPartInsSubSidesStandardCase(iSpec,iSF, PartInsSubSides)

    !----- 0.: go through (sub)sides if present in proc
    IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) THEN
      CYCLE
    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL ABORT(__STAMP__,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF

#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

    DO iSide=1,BCdata_auxSF(currentBC)%SideNumber

      !
      IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
      END IF

      BCSideID = BCdata_auxSF(currentBC)%SideList(iSide)
      ElemID   = SideToElem(1,BCSideID)

      ! mortar elements MIGHT need to be treated differently. Take the side pointing back
      IF (ElemID.LT.1) THEN
        ElemID = SideToElem(2,BCSideID)
        iLocSide = SideToElem(4,BCSideID)
      ELSE
        iLocSide = SideToElem(3,BCSideID)
      END IF

      !  Get global ElemID and SideID from local SideID
      globElemId = ElemID + offSetElem
      SideID     = GetGlobalNonUniqueSideID(globElemId,iLocSide)

      ! NodeCoords for TriaGeometry
      IF (TriaSurfaceFlux) xyzNod(1:3) = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(1:3)

      DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
          !-- compute parallelogram of triangle
        IF (TriaSurfaceFlux) THEN
          Node1 = jSample+1     ! normal = cross product of 1-2 and 1-3 for first triangle
          Node2 = jSample+2     !          and 1-3 and 1-4 for second triangle
          Vector1       = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
          Vector2       = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
          midpoint(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
          ndist   (1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)
        END IF

        !----- 1.: set positions
        !-- compute number of to be inserted particles
        IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
          PartInsSubSide = PartInsSubSides(iSample,jSample,iSide)

        ELSE IF (DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
              CALL CalcPartInsPoissonDistr(iSpec, iSF, iSample, jSample, iSide, PartInsSubSide)
        ELSE ! DoTimeDepInflow
          CALL RANDOM_NUMBER(RandVal1)
          PartInsSubSide = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity * dt*RKdtFrac                                           &
                         * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR+RandVal1)
        END IF ! DoPoissonRounding

        !-- proceed with calculated to be inserted particles
        IF (PartInsSubSide.LT.0) THEN
          CALL ABORT(__STAMP__,'ERROR in ParticleSurfaceflux: PartInsSubSide.LT.0!')
        ELSE IF (PartInsSubSide.LE.0) THEN
          CYCLE
        END IF

        NbrOfParticle  = PartInsSubSide + NbrOfParticle
        ALLOCATE( particle_positions(1:PartInsSubSide*3))

        IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
          ALLOCATE( particle_xis(1:PartInsSubSide*2))
        END IF ! VeloIsNormal

        !-- put particles in subside (rejections are used if constraint reduces actual inserted number)
        iPart   = 1
        nReject = 0
        allowedRejections = 0

        DO WHILE (iPart+allowedRejections .LE. PartInsSubSide)
          IF (TriaSurfaceFlux) THEN
            Particle_pos(1:3) = CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
          ! .NOT.TriaSurfaceFlux
          ELSE
            Particle_pos(1:3) = CalcPartPosBezier(iSpec,iSF,iSample,jSample,iSide,SideID)
          END IF !TriaSurfaceFlux

          AcceptPos=.TRUE.
          IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN !check rmax-rejection
            IF (.NOT.InSideCircularInflow(iSpec, iSF, iSide, Particle_pos)) AcceptPos=.FALSE.
          END IF ! CircularInflow

          !-- save position if accepted:
          IF (AcceptPos) THEN
            particle_positions(iPart*3-2) = Particle_pos(1)
            particle_positions(iPart*3-1) = Particle_pos(2)
            particle_positions(iPart*3  ) = Particle_pos(3)
            IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
              particle_xis(iPart*2-1) = xi(1)
              particle_xis(iPart*2  ) = xi(2)
            END IF !VeloIsNormal
            iPart=iPart+1
          ELSE
            nReject=nReject+1
            ! check rmax-rejection
            IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
              allowedRejections = allowedRejections+1
            END IF
          END IF
        END DO !put particles in subside: WHILE(iPart+allowedRejections .LE. PartInsSubSide)

        PartInsSubSide = PartInsSubSide - allowedRejections
        NbrOfParticle = NbrOfParticle - allowedRejections

        !-- Fill Particle Informations (PartState, Partelem, etc.)
        ParticleIndexNbr = 1
        DO iPart = 1,PartInsSubSide
          IF ((iPart.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) &
              ParticleIndexNbr = PDM%nextFreePosition(iPartTotal + 1 + PDM%CurrentNextFreePosition)
          IF (ParticleIndexNbr .NE. 0) THEN
            PartState(1:3,ParticleIndexNbr) = particle_positions(3*(iPart-1)+1:3*(iPart-1)+3)

            IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal.AND.(.NOT.TriaSurfaceFlux)) THEN
              ! use velo as dummy-storage for xi!
              PartState(4:5,ParticleIndexNbr) = particle_xis(2*(iPart-1)+1:2*(iPart-1)+2)
            END IF
            ! shift LastPartPos minimal into cell for fail-safe tracking
            LastPartPos(1:3,ParticleIndexNbr)    = PartState(1:3,ParticleIndexNbr)
            PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
            PDM%IsNewPart(     ParticleIndexNbr) = .TRUE.
            PEM%Element(ParticleIndexNbr)        = ElemID
            ! needed when ParticlePush is not executed, e.g. "delay"
            PEM%LastElement(ParticleIndexNbr)    = ElemID
            iPartTotal = iPartTotal + 1

#if CODE_ANALYZE
            CALL AnalyzePartPos(ParticleIndexNbr)
#endif /*CODE_ANALYZE*/
          ELSE
            CALL ABORT(__STAMP__,'ERROR in ParticleSurfaceflux: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
          END IF
        END DO

        DEALLOCATE(particle_positions)

        IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) DEALLOCATE(particle_xis)
!----- 2a.: set velocities if special for each subside
        CALL SetSurfacefluxVelocities(iSpec,iSF,iSample,jSample,BCSideID,SideID,NbrOfParticle,PartInsSubSide)

        PartsEmitted = PartsEmitted + PartInsSubSide

#if USE_LOADBALANCE
        ! used for calculating LoadBalance of tCurrent(LB_SURFFLUX) ==> "2b.: set remaining properties"
        nSurfacefluxPerElem(ElemID)=nSurfacefluxPerElem(ElemID)+PartInsSubSide
#endif /*USE_LOADBALANCE*/

      END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
#if USE_LOADBALANCE
      CALL LBElemSplitTime(ElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO ! iSide

    IF (NbrOfParticle.NE.iPartTotal) CALL ABORT(__STAMP__, 'Error 2 in ParticleSurfaceflux!')

    !----- 2b.: set remaining properties
    CALL SetParticleMass(iSpec,NbrOfParticle)
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
    PDM%ParticleVecLength       = PDM%ParticleVecLength       + NbrOfParticle

#if USE_LOADBALANCE
    CALL LBPauseTime(LB_SURFFLUX,tLBStart)
#endif /*USE_LOADBALANCE*/

    IF (NbrOfParticle.NE.PartsEmitted) THEN
      ! should be equal for including the following lines in tSurfaceFlux
      CALL ABORT(__STAMP__,'ERROR in ParticleSurfaceflux: NbrOfParticle.NE.PartsEmitted')
    END IF

  END DO !iSF
END DO !iSpec

END SUBROUTINE ParticleSurfaceflux


FUNCTION InSideCircularInflow(iSpec, iSF, iSide, Particle_pos)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES)
INTEGER, INTENT(IN)         :: iSpec,iSF,iSide
REAL, INTENT(IN)            :: Particle_pos(3)
LOGICAL                     :: InSideCircularInflow
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: point(2),radius,origin(2)
!===================================================================================================================================

origin = Species(iSpec)%Surfaceflux(iSF)%origin

SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide))
  !- RejectType=0 : complete side is inside valid bounds
  CASE(0)
    InSideCircularInflow = .TRUE.

  !- RejectType=1 : complete side is outside of valid bounds
  CASE(1)
    CALL ABORT(__STAMP__, 'Side outside of valid bounds was considered although nVFR=0...?!')

  !- RejectType=2 : side is partly inside valid bounds
  CASE(2)
    point(1) = Particle_pos(Species(iSpec)%Surfaceflux(iSF)%dir(2))-origin(1)
    point(2) = Particle_pos(Species(iSpec)%Surfaceflux(iSF)%dir(3))-origin(2)
    radius   = SQRT( (point(1))**2+(point(2))**2 )
    IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
      InSideCircularInflow = .TRUE.
    ELSE
      InSideCircularInflow = .FALSE.
    END IF

  CASE DEFAULT
    CALL ABORT(__STAMP__, 'Wrong SurfFluxSideRejectType!')
END SELECT ! SurfFluxSideRejectType

END FUNCTION InSideCircularInflow


FUNCTION CalcPartPosBezier(iSpec,iSF,iSample,jSample,iSide,SideID)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D,BezierSampleXi
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Vars           ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES)
INTEGER,INTENT(IN)          :: iSpec
INTEGER,INTENT(IN)          :: iSF
INTEGER,INTENT(IN)          :: iSample
INTEGER,INTENT(IN)          :: jSample
INTEGER,INTENT(IN)          :: iSide
INTEGER,INTENT(IN)          :: SideID
REAL                        :: CalcPartPosBezier(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: xiab(1:2,1:2),xi(2)
REAL                        :: D,E,F,G
REAL                        :: gradXiEta2D(1:2,1:2),gradXiEta3D(1:2,1:3)
REAL                        :: RandVal1,RandVal2(2)
INTEGER                     :: iLoop
!===================================================================================================================================

iLoop = 0

! ARM for xi considering the dA of the Subside in RefSpace
DO
  iLoop = iLoop+1
  CALL RANDOM_NUMBER(RandVal2)
  xiab(1,1:2) = (/BezierSampleXi(iSample-1),BezierSampleXi(iSample)/) !correct order?!?
  xiab(2,1:2) = (/BezierSampleXi(jSample-1),BezierSampleXi(jSample)/) !correct order?!?
  xi          = (xiab(:,2)-xiab(:,1))*RandVal2+xiab(:,1)

  IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
    IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      CALL EvaluateBezierPolynomialAndGradient( xi                                                                     &
                                              , NGeo                                                                   &
                                              , 2                                                                      &
                                              , Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample    &
                                              , iSide)%BezierControlPoints2D(1:2,0:NGeo,0:NGeo)                        &
                                              , Gradient = gradXiEta2D)
      E = DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(1,1:2))
      F = DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(2,1:2))
      G = DOT_PRODUCT(gradXiEta2D(2,1:2),gradXiEta2D(2,1:2))
    ELSE
      CALL EvaluateBezierPolynomialAndGradient( xi                                                                     &
                                              , NGeo                                                                   &
                                              , 3                                                                      &
                                              , BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID)                        &
                                              , Gradient = gradXiEta3D)
      E = DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
      F = DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
      G = DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
    END IF ! .NOT.VeloIsNormal

    D = SQRT(E*G-F*F)
    ! scaled Jacobian of xi
    D = D/Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Dmax

    ! arbitrary warning threshold
    IF (D .GT. 1.01) THEN
      IPWRITE(*,'(I4,x,A28,I0,A9,I0,A22,I0)') 'WARNING: ARM of SurfaceFlux ',iSF,' of Spec ',iSpec,' has inaccurate Dmax! ',D
    END IF

    CALL RANDOM_NUMBER(RandVal1)
    ! accept xi
    IF (RandVal1.LE.D) THEN
      EXIT
    ELSE
      ! arbitrary warning threshold
      IF (MOD(iLoop,100).EQ.0) THEN
        IPWRITE(*,'(I4,x,A28,I0,A9,I0,A18,I0)') 'WARNING: ARM of SurfaceFlux ',iSF,' of Spec ',iSpec,' has reached loop ',iLoop
        IPWRITE(*,'(I4,x,A19,2(x,E16.8))')      '         R, D/Dmax:',RandVal1,D
      END IF
    END IF

  ! no ARM -> accept xi
  ELSE
    EXIT
  END IF
END DO ! Jacobian-based ARM-loop

IF (MINVAL(XI).LT.-1.) THEN
  IPWRITE(UNIT_StdOut,'(I0,A,E16.8)') ' Xi<-1',XI
END IF

IF (MAXVAL(XI).GT. 1.) THEN
  IPWRITE(UNIT_StdOut,'(I0,A,E16.8)') ' Xi>1',XI
END IF

CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=CalcPartPosBezier)

END FUNCTION CalcPartPosBezier


FUNCTION CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)            :: xyzNod(3), Vector1(3), Vector2(3), ndist(3), midpoint(3)
REAL                        :: CalcPartPosTriaSurface(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal2(2), PartDistance
REAL, PARAMETER             :: eps_nontria=1.0E-6
!===================================================================================================================================
CALL RANDOM_NUMBER(RandVal2)

! prevent inconsistency with non-triatracking by bilinear-routine (tol. might be increased)
IF (.NOT.TrackingMethod.EQ.TRIATRACKING) THEN
  ! shift randVal off from 0 and 1
  RandVal2 = RandVal2 + eps_nontria*(1. - 2.*RandVal2)
  ! sum must not be 1, since this corresponds to third edge
  DO WHILE (ABS(RandVal2(1)+RandVal2(2)-1.0).LT.eps_nontria)
    CALL RANDOM_NUMBER(RandVal2)
    RandVal2 = RandVal2 + eps_nontria*(1. - 2.*RandVal2)
  END DO
END IF

CalcPartPosTriaSurface = xyzNod + Vector1 * RandVal2(1)
CalcPartPosTriaSurface = CalcPartPosTriaSurface + Vector2 * RandVal2(2)
! Distance from v1-v2
PartDistance = ndist(1)*(CalcPartPosTriaSurface(1)-midpoint(1)) &
             + ndist(2)*(CalcPartPosTriaSurface(2)-midpoint(2)) &
             + ndist(3)*(CalcPartPosTriaSurface(3)-midpoint(3))

! flip into right triangle if outside
IF (PartDistance.GT.0.) THEN
  CalcPartPosTriaSurface(1:3) = 2.*midpoint(1:3)-CalcPartPosTriaSurface(1:3)
END IF

END FUNCTION CalcPartPosTriaSurface


SUBROUTINE CalcPartInsPoissonDistr(iSpec, iSF, iSample, jSample, iSide, PartInsSubSide)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars           ,ONLY: dt
USE MOD_Part_Emission_Tools     ,ONLY: SamplePoissonDistri
USE MOD_Particle_TimeDisc_Vars  ,ONLY: RKdtFrac
USE MOD_Particle_Vars           ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, iSample, jSample, iSide
INTEGER, INTENT(OUT)        :: PartInsSubSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: PartIns
!===================================================================================================================================
PartIns = Species(iSpec)%Surfaceflux(iSF)%PartDensity * dt*RKdtFrac  &
        * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR

IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
  CALL ABORT(__STAMP__, 'ERROR in ParticleSurfaceflux: flux is too large for poisson sampling!')

! poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
ELSE
  CALL SamplePoissonDistri( PartIns , PartInsSubSide )
END IF

END SUBROUTINE CalcPartInsPoissonDistr


!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
SUBROUTINE SetSurfacefluxVelocities(FractNbr,iSF,iSample,jSample,BCSideID,SideID,NbrOfParticle,PartIns)
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals         ,ONLY: VECNORM
USE MOD_Particle_Surfaces        ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars   ,ONLY: SurfMeshSubSideData,TriaSurfaceFlux
USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr
INTEGER,INTENT(IN)               :: iSF
INTEGER,INTENT(IN)               :: iSample
INTEGER,INTENT(IN)               :: jSample
INTEGER,INTENT(IN)               :: BCSideID
INTEGER,INTENT(IN)               :: SideID
INTEGER,INTENT(IN)               :: NbrOfParticle
INTEGER,INTENT(IN)               :: PartIns
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,PositionNbr
REAL                             :: Vec3D(3),vec_nIn(1:3),vec_t1(1:3),vec_t2(1:3)
CHARACTER(30)                    :: velocityDistribution
REAL                             :: VeloVecIC(1:3)
!===================================================================================================================================

! Return if no particles are to be inserted
IF (PartIns.LT.1) RETURN

SELECT CASE(TRIM(Species(FractNbr)%Surfaceflux(iSF)%velocityDistribution))
  CASE('constant')
    velocityDistribution='constant'
    VeloVecIC(1:3) = Species(FractNbr)%Surfaceflux(iSF)%VeloVecIC(1:3)
    VeloVecIC(1:3) = VeloVecIC(1:3) / VECNORM(VeloVecIC(1:3))
  CASE DEFAULT
    CALL ABORT(__STAMP__,'Wrong VelocityDistribution!')
END SELECT

IF (.NOT.Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal) THEN
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  vec_t1 (1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
  vec_t2 (1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
END IF ! .NOT.VeloIsNormal

!-- set velocities
SELECT CASE(TRIM(velocityDistribution))
  CASE('constant')

    DO i = NbrOfParticle-PartIns+1,NbrOfParticle
      PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
      IF (PositionNbr .NE. 0) THEN
        !-- In case of side-normal velocities: calc n-vector at particle position, xi was saved in PartState(4:5)
        IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
          vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
          vec_t1 (1:3) = 0. ! dummy
          vec_t2 (1:3) = 0. ! dummy
        ELSE IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal) THEN
          CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
          vec_nIn(1:3) = -vec_nIn(1:3)
          vec_t1 (1:3) = 0. ! dummy
          vec_t2 (1:3) = 0. ! dummy
        ELSE
        vec_nIn(1:3) = VeloVecIC(1:3)
        END IF ! VeloIsNormal

        !-- build complete velo-vector
        Vec3D(1:3) = vec_nIn(1:3) * Species(FractNbr)%Surfaceflux(iSF)%VeloIC
        PartState(4:6,PositionNbr) = Vec3D(1:3)
      END IF ! PositionNbr .NE. 0
    END DO ! i = ...NbrOfParticle

  CASE DEFAULT
    CALL ABORT(__STAMP__,'Wrong VelocityDistribution!')
END SELECT

END SUBROUTINE SetSurfacefluxVelocities

END MODULE MOD_Surface_Flux
