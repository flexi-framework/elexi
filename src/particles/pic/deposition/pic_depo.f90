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

!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
MODULE  MOD_PICDepo                                                                                
 IMPLICIT NONE                                                                                   
 PRIVATE                                                                                         
!===================================================================================================================================
INTERFACE Deposition
  MODULE PROCEDURE Deposition
END INTERFACE

INTERFACE InitializeDeposition
  MODULE PROCEDURE InitializeDeposition
END INTERFACE

INTERFACE FinalizeDeposition
  MODULE PROCEDURE FinalizeDeposition
END INTERFACE

PUBLIC:: Deposition, InitializeDeposition, FinalizeDeposition
!===================================================================================================================================

CONTAINS                                                                                           
                                                                                                   
SUBROUTINE InitializeDeposition
!===================================================================================================================================
! Initialize the deposition variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_Mesh_Vars,              ONLY:nElems, Elem_xGP, sJ,nGlobalElems
USE MOD_Particle_Mesh_Vars
USE MOD_Interpolation_Vars,     ONLY:xGP,wBary,wGP
USE MOD_Particle_Basis,         ONLY:ComputeBernsteinCoeff
USE MOD_Basis,                  ONLY:BarycentricWeights,InitializeVandermonde
USE MOD_Basis,                  ONLY:LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_ChangeBasis,            ONLY:ChangeBasis3D
USE MOD_PreProc,                ONLY:PP_N
USE MOD_ReadInTools,            ONLY:GETREAL,GETINT,GETLOGICAL,GETSTR,GETREALARRAY,GETINTARRAY
USE MOD_PICInterpolation_Vars,  ONLY:InterpolationType
USE MOD_Eval_xyz,               ONLY:TensorProductInterpolation
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
#if USE_MPI
USE MOD_Particle_MPI_Vars,      ONLY:DoExternalParts
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE          :: wBary_tmp(:),Vdm_GaussN_EquiN(:,:), Vdm_GaussN_NDepo(:,:)
REAL,ALLOCATABLE          :: dummy(:,:,:,:),dummy2(:,:,:,:),xGP_tmp(:),wGP_tmp(:)
INTEGER                   :: ALLOCSTAT, iElem, i, j, k, m, dir, weightrun, mm, r, s, t, iSFfix, iPoint, iVec, iBC
REAL                      :: MappedGauss(1:PP_N+1), xmin, ymin, zmin, xmax, ymax, zmax, x0
REAL                      :: auxiliary(0:3),weight(1:3,0:3), nTotalDOF, VolumeShapeFunction, r_sf_average
REAL                      :: DetLocal(1,0:PP_N,0:PP_N,0:PP_N), DetJac(1,0:1,0:1,0:1)
REAL, ALLOCATABLE         :: Vdm_tmp(:,:)
REAL                      :: SFdepoLayersCross(3),BaseVector(3),BaseVector2(3),SFdepoLayersChargedens,n1(3),n2(3),nhalf
REAL                      :: NormVecCheck(3),eps,diff,BoundPoints(3,8),BVlengths(2),DistVec(3),Dist(3)
CHARACTER(32)             :: hilf, hilf2
CHARACTER(255)            :: TimeAverageFile
LOGICAL                   :: LayerOutsideOfBounds, ChangeOccured
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE DEPOSITION...'

DoDeposition = GETLOGICAL('PIC-DoDeposition','T')
IF(.NOT.DoDeposition) THEN
  ! fill deposition type with empty string
  DepositionType='NONE'
  RETURN
END IF
DepositionType = GETSTR('PIC-Deposition-Type','nearest_blurrycenter')
! check for interpolation type incompatibilities (cannot be done at interpolation_init
! because DepositionType is not known yet)
IF((TRIM(InterpolationType).EQ.'nearest_gausspoint').AND. &
   (TRIM(DepositionType).NE.'nearest_gausspoint')) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR in pic_depo.f90: Interpolation type nearest_gausspoint only allowed with same deposition type!')
END IF
!--- Allocate arrays for charge density collection and initialize
ALLOCATE(PartSource(1:4,0:PP_N,0:PP_N,0:PP_N,nElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR in pic_depo.f90: Cannot allocate PartSource!')
END IF
PartSource=0.

!--- check if chargedensity is computed from TimeAverageFile
TimeAverageFile = GETSTR('PIC-TimeAverageFile','none')
IF (TRIM(TimeAverageFile).NE.'none') THEN
  CALL ReadTimeAverage(TimeAverageFile)
  DoDeposition=.FALSE.
  DepositionType='NONE'
  RETURN
END IF

!--- init DepositionType-specific vars
SELECT CASE(TRIM(DepositionType))
    
CASE('nearest_blurrycenter')
    
CASE('nearest_blurycenter')
  DepositionType = 'nearest_blurrycenter'
  
CASE('cell_volweight') 
  ALLOCATE(CellVolWeightFac(0:PP_N),wGP_tmp(0:PP_N) , xGP_tmp(0:PP_N))
  ALLOCATE(CellVolWeight_Volumes(0:1,0:1,0:1,nElems))
  CellVolWeightFac(0:PP_N) = xGP(0:PP_N)
  CellVolWeightFac(0:PP_N) = (CellVolWeightFac(0:PP_N)+1.0)/2.0
  CALL LegendreGaussNodesAndWeights(1,xGP_tmp,wGP_tmp)
  ALLOCATE( Vdm_tmp(0:1,0:PP_N))
  CALL InitializeVandermonde(PP_N,1,wBary,xGP,xGP_tmp,Vdm_tmp)
  DO iElem=1, nElems
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          DetLocal(1,i,j,k)=1./sJ(i,j,k,iElem,0)
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
    CALL ChangeBasis3D(1,PP_N, 1,Vdm_tmp, DetLocal(:,:,:,:),DetJac(:,:,:,:))
    DO k=0,1
      DO j=0,1
        DO i=0,1
          CellVolWeight_Volumes(i,j,k,iElem) = DetJac(1,i,j,k)*wGP_tmp(i)*wGP_tmp(j)*wGP_tmp(k)
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
  END DO
  DEALLOCATE(Vdm_tmp)
  DEALLOCATE(wGP_tmp, xGP_tmp)
  
CASE('epanechnikov') 
  r_sf     = GETREAL('PIC-epanechnikov-radius','1.')
  r2_sf = r_sf * r_sf 
  ALLOCATE( tempcharge(1:nElems))
CASE('nearest_gausspoint')
  ! Allocate array for particle positions in -1|1 space (used for deposition as well as interpolation)
  ! only if NOT DoRefMapping
  IF(.NOT.DoRefMapping)THEN
    ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,' Cannot allocate partposref!')
  END IF
  ! compute the borders of the virtual volumes around the gauss points in -1|1 space
  ALLOCATE(GaussBorder(1:PP_N),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in pic_depo.f90: Cannot allocate Mapped Gauss Border Coords!')
  END IF
  DO i=0,PP_N
    ! bullshit here, use xGP
    !CALL TensorProductInterpolation(Elem_xGP(:,i,1,1,1),Temp(:),1)
    !MappedGauss(i+1) = Temp(1)
    MappedGauss(i+1) = xGP(i)
  END DO
  DO i = 1,PP_N
    GaussBorder(i) = (MappedGauss(i+1) + MappedGauss(i))/2
  END DO
  ! allocate array for saving the gauss points of particles for nearest_gausspoint interpolation
  IF(TRIM(InterpolationType).EQ.'nearest_gausspoint')THEN
    SDEALLOCATE(PartPosGauss)
    ALLOCATE(PartPosGauss(1:PDM%maxParticleNumber,1:3),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR in pic_depo.f90: Cannot allocate Part Pos Gauss!')
    END IF
  END IF
  
CASE('shape_function','shape_function_simple')
  r_sf     = GETREAL('PIC-shapefunction-radius','1.')
  alpha_sf = GETINT('PIC-shapefunction-alpha','2')
  DoSFEqui = GETLOGICAL('PIC-shapefunction-equi','F')
  BetaFac = beta(1.5, REAL(alpha_sf) + 1.)
  w_sf = 1./(2. * BetaFac * REAL(alpha_sf) + 2 * BetaFac) &
                        * (REAL(alpha_sf) + 1.)/(PI*(r_sf**3))
  r2_sf = r_sf * r_sf 
  r2_sf_inv = 1./r2_sf
  !-- fixes for shape func depo at planar BCs
  NbrOfSFdepoFixes = GETINT('PIC-NbrOfSFdepoFixes','0')
  PrintSFDepoWarnings=GETLOGICAL('PrintSFDepoWarnings','.FALSE.')
  IF (NbrOfSFdepoFixes.GT.0) THEN
    SFdepoFixesEps = GETREAL('PIC-SFdepoFixesEps','0.')
    SDEALLOCATE(SFdepoFixesGeo)
    SDEALLOCATE(SFdepoFixesChargeMult)
    SDEALLOCATE(SFdepoFixesBounds)
    ALLOCATE(SFdepoFixesGeo(1:NbrOfSFdepoFixes,1:2,1:3),STAT=ALLOCSTAT)  !1:nFixes;1:2(base,normal);1:3(x,y,z)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR in pic_depo.f90: Cannot allocate SFdepoFixesGeo!')
    END IF
    ALLOCATE(SFdepoFixesChargeMult(1:NbrOfSFdepoFixes),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR in pic_depo.f90: Cannot allocate SFdepoFixesChargeMult!')
    END IF
    ALLOCATE(SFdepoFixesPartOfLink(0:NbrOfSFdepoFixes),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoFixesPartOfLink!')
    END IF
    SFdepoFixesPartOfLink=.FALSE.
    ALLOCATE(SFdepoFixesBounds(1:NbrOfSFdepoFixes,1:2,1:3),STAT=ALLOCSTAT)  !1:nFixes;1:2(min,max);1:3(x,y,z)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR in pic_depo.f90: Cannot allocate SFdepoFixesBounds!')
    END IF
    DO iSFfix=1,NbrOfSFdepoFixes
      WRITE(UNIT=hilf,FMT='(I0)') iSFfix
      SFdepoFixesGeo(iSFfix,1,1:3) = &
        GETREALARRAY('PIC-SFdepoFixes'//TRIM(hilf)//'-Basepoint',3,'0. , 0. , 0.')
      SFdepoFixesGeo(iSFfix,2,1:3) = &
        GETREALARRAY('PIC-SFdepoFixes'//TRIM(hilf)//'-Normal',3,'1. , 0. , 0.') !directed outwards
      IF (SFdepoFixesGeo(iSFfix,2,1)**2 + SFdepoFixesGeo(iSFfix,2,2)**2 + SFdepoFixesGeo(iSFfix,2,3)**2 .GT. 0.) THEN
        SFdepoFixesGeo(iSFfix,2,1:3) = SFdepoFixesGeo(iSFfix,2,1:3) &
          / SQRT(SFdepoFixesGeo(iSFfix,2,1)**2 + SFdepoFixesGeo(iSFfix,2,2)**2 + SFdepoFixesGeo(iSFfix,2,3)**2)
      ELSE
        CALL abort(&
      __STAMP__&
      ,' SFdepoFixesXX-Normal is zero for Fix ',iSFfix)
      END IF
      SFdepoFixesChargeMult(iSFfix) = &
        GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-ChargeMult','1.')
      WRITE(UNIT=hilf2,FMT='(E16.8)') -HUGE(1.0)
      SFdepoFixesBounds(iSFfix,1,1)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-xmin',TRIM(hilf2))
      SFdepoFixesBounds(iSFfix,1,2)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-ymin',TRIM(hilf2))
      SFdepoFixesBounds(iSFfix,1,3)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-zmin',TRIM(hilf2))
      WRITE(UNIT=hilf2,FMT='(E16.8)') HUGE(1.0)
      SFdepoFixesBounds(iSFfix,2,1)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-xmax',TRIM(hilf2))
      SFdepoFixesBounds(iSFfix,2,2)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-ymax',TRIM(hilf2))
      SFdepoFixesBounds(iSFfix,2,3)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-zmax',TRIM(hilf2))
    END DO
    NbrOfSFdepoFixLinks = GETINT('PIC-NbrOfSFdepoFixLinks','0')
    IF (NbrOfSFdepoFixLinks.GT.0) THEN
      ALLOCATE(SFdepoFixLinks(1:NbrOfSFdepoFixLinks,1:3),STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) THEN
        CALL abort(&
          __STAMP__&
          ,'ERROR in pic_depo.f90: Cannot allocate SFdepoFixesChargeMult!')
      END IF
      DO iSFfix=1,NbrOfSFdepoFixLinks
        WRITE(UNIT=hilf,FMT='(I0)') iSFfix
        SFdepoFixLinks(iSFfix,1:2) = &
          GETINTARRAY('PIC-SFdepoFixLink'//TRIM(hilf),2,'1 , 2')
        IF (   SFdepoFixLinks(iSFfix,1).GT.NbrOfSFdepoFixes &
          .OR. SFdepoFixLinks(iSFfix,2).GT.NbrOfSFdepoFixes &
          .OR. SFdepoFixLinks(iSFfix,1).LE.0 &
          .OR. SFdepoFixLinks(iSFfix,2).LE.0) THEN
          CALL abort(&
            __STAMP__&
            ,' SFdepoFixes not defined for Link ',iSFfix)
        ELSE
          SFdepoFixesPartOfLink(SFdepoFixLinks(iSFfix,1))=.TRUE.
          n1=SFdepoFixesGeo(SFdepoFixLinks(iSFfix,1),2,1:3)
          SFdepoFixesPartOfLink(SFdepoFixLinks(iSFfix,2))=.TRUE.
          n2=SFdepoFixesGeo(SFdepoFixLinks(iSFfix,2),2,1:3)
          nhalf=ACOS(-(DOT_PRODUCT(n1,n2))) !n already normalized
          IF (nhalf.EQ.0.) THEN
            CALL abort(__STAMP__, &
              'ERROR in pic_depo.f90: angle between vectors of SFdepoFixLinks is zero!')
          ELSE
            nhalf=PI/nhalf
          END IF
          IF ( ABS(nhalf-1.5).LT.SFdepoFixesEps ) THEN !120 deg
            SFdepoFixLinks(iSFfix,3)=-2 !negative as flag for special case (120 deg), 2 since abs(-2) is needed for respective loop
            SWRITE(*,*) 'SFdepoFixLink ',iSFfix,' was determined to have an angle of 120 deg...'
          ELSE IF ( ABS(nhalf-NINT(nhalf)).GT.SFdepoFixesEps .OR. NINT(nhalf).LT.2 ) THEN !check if integer fraction (>2) of 180 deg
            CALL abort(__STAMP__, &
              'ERROR in pic_depo.f90: angle between vectors of SFdepoFixLink is neither 120 deg nor a fraction of 180 deg!', &
              NINT(nhalf),ABS(nhalf-NINT(nhalf)))
          ELSE !integer fraction of 180 deg
            SFdepoFixLinks(iSFfix,3) = NINT(nhalf)
            SWRITE(*,*) 'SFdepoFixLink ',iSFfix,' was determined to divide 180 deg into ',NINT(nhalf),' parts...'
          END IF
        END IF
      END DO
    END IF !NbrOfSFdepoFixLinks>0
    DO iSFfix=0,NbrOfSFdepoFixes
      SWRITE(*,*) 'SFdepoFix ',iSFfix,' is part of link: ',SFdepoFixesPartOfLink(iSFfix)
    END DO
  END IF !NbrOfSFdepoFixes>0

  !-- const. PartSource layer for shape func depo at planar BCs
  NbrOfSFdepoLayers = GETINT('PIC-NbrOfSFdepoLayers','0')
  IF (NbrOfSFdepoLayers.GT.0) THEN
    SDEALLOCATE(SFdepoLayersGeo)
    SDEALLOCATE(SFdepoLayersBounds)
    SDEALLOCATE(SFdepoLayersUseFixBounds)
    SDEALLOCATE(SFdepoLayersSpace)
    SDEALLOCATE(SFdepoLayersBaseVector)
    SDEALLOCATE(SFdepoLayersSpec)
    SDEALLOCATE(SFdepoLayersPartNum)
    SDEALLOCATE(SFdepoLayersRadius)
    ALLOCATE(SFdepoLayersGeo(1:NbrOfSFdepoLayers,1:2,1:3),STAT=ALLOCSTAT)  !1:nLayers;1:2(base,normal);1:3(x,y,z)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersGeo!')
    END IF
    ALLOCATE(SFdepoLayersBounds(1:NbrOfSFdepoLayers,1:2,1:3),STAT=ALLOCSTAT)  !1:nLayers;1:2(min,max);1:3(x,y,z)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersBounds!')
    END IF
    ALLOCATE(SFdepoLayersUseFixBounds(1:NbrOfSFdepoLayers),STAT=ALLOCSTAT)  !1:nLayers;1:2(min,max);1:3(x,y,z)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersUseFixBounds!')
    END IF
    ALLOCATE(SFdepoLayersSpace(1:NbrOfSFdepoLayers),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersSpace!')
    END IF
    ALLOCATE(SFdepoLayersBaseVector(1:NbrOfSFdepoLayers,1:2,1:3),STAT=ALLOCSTAT)  !1:nLayers;1:2(base,normal);1:3(x,y,z)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersBaseVector!')
    END IF
    ALLOCATE(SFdepoLayersSpec(1:NbrOfSFdepoLayers),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersSpec!')
    END IF
    ALLOCATE(SFdepoLayersPartNum(1:NbrOfSFdepoLayers),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersPartNum!')
    END IF
    ALLOCATE(SFdepoLayersRadius(1:NbrOfSFdepoLayers),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
        'ERROR in pic_depo.f90: Cannot allocate SFdepoLayersRadius!')
    END IF
    SFdepoLayersRadius=0. !not used for cuboid...
    DO iSFfix=1,NbrOfSFdepoLayers
#if !(defined (PP_HDG) && (PP_nVar==1))
      CALL abort(__STAMP__, &
        ' NbrOfSFdepoLayers are only implemented for electrostatic HDG!')
#endif
      WRITE(UNIT=hilf,FMT='(I0)') iSFfix
      SFdepoLayersGeo(iSFfix,1,1:3) = &
        GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-Basepoint',3,'0. , 0. , 0.')
      SFdepoLayersGeo(iSFfix,2,1:3) = &
        GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-Normal',3,'1. , 0. , 0.') !directed outwards
      IF (SFdepoLayersGeo(iSFfix,2,1)**2 + SFdepoLayersGeo(iSFfix,2,2)**2 + SFdepoLayersGeo(iSFfix,2,3)**2 .GT. 0.) THEN
        SFdepoLayersGeo(iSFfix,2,1:3) = SFdepoLayersGeo(iSFfix,2,1:3) &
          / SQRT(SFdepoLayersGeo(iSFfix,2,1)**2 + SFdepoLayersGeo(iSFfix,2,2)**2 + SFdepoLayersGeo(iSFfix,2,3)**2)
      ELSE
        CALL abort(__STAMP__&
          ,' SFdepoLayersXX-Normal is zero for Layer ',iSFfix)
      END IF
      WRITE(UNIT=hilf2,FMT='(E16.8)') -HUGE(1.0)
      SFdepoLayersBounds(iSFfix,1,1)   = MAX(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-xmin',TRIM(hilf2)),GEO%xmin-r_sf)
      SFdepoLayersBounds(iSFfix,1,2)   = MAX(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-ymin',TRIM(hilf2)),GEO%ymin-r_sf)
      SFdepoLayersBounds(iSFfix,1,3)   = MAX(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-zmin',TRIM(hilf2)),GEO%zmin-r_sf)
      WRITE(UNIT=hilf2,FMT='(E16.8)') HUGE(1.0)
      SFdepoLayersBounds(iSFfix,2,1)   = MIN(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-xmax',TRIM(hilf2)),GEO%xmax+r_sf)
      SFdepoLayersBounds(iSFfix,2,2)   = MIN(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-ymax',TRIM(hilf2)),GEO%ymax+r_sf)
      SFdepoLayersBounds(iSFfix,2,3)   = MIN(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-zmax',TRIM(hilf2)),GEO%zmax+r_sf)
      IF (NbrOfSFdepoFixes.EQ.0) THEN
        SFdepoLayersUseFixBounds(iSFfix) = .FALSE.
          ELSE
        SFdepoLayersUseFixBounds(iSFfix) = GETLOGICAL('PIC-SFdepoLayers'//TRIM(hilf)//'-UseFixBounds','.TRUE.')
      END IF
      SFdepoLayersSpace(iSFfix) = &
        GETSTR('PIC-SFdepoLayers'//TRIM(hilf)//'-Space','cuboid')
      LayerOutsideOfBounds = .FALSE.
      ChangeOccured = .FALSE.
      SELECT CASE (TRIM(SFdepoLayersSpace(iSFfix)))
      CASE('cuboid')
        SFdepoLayersBaseVector(iSFfix,1,1:3)   = GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-BaseVector1',3,'0. , 1. , 0.')
        SFdepoLayersBaseVector(iSFfix,2,1:3)   = GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-BaseVector2',3,'0. , 0. , 1.')
        !--check if BV1/2 are perpendicular to NormalVec (need to be for domain reduction with BoundPoints, otherwise should not matter)
        eps=1.0E-10
        NormVecCheck=CROSS(SFdepoLayersBaseVector(iSFfix,1,1:3),SFdepoLayersBaseVector(iSFfix,2,1:3))
        IF (NormVecCheck(1)**2 + NormVecCheck(2)**2 + NormVecCheck(3)**2 .GT. 0.) THEN
          NormVecCheck(1:3) = NormVecCheck(1:3) / SQRT(NormVecCheck(1)**2 + NormVecCheck(2)**2 + NormVecCheck(3)**2)
          diff = SQRT(DOT_PRODUCT(NormVecCheck-SFdepoLayersGeo(iSFfix,2,:),NormVecCheck-SFdepoLayersGeo(iSFfix,2,:))) !NormVecCheck - SFdepoLayersGeo(iSFfix,2,1) should be (/0,0,0/) when identical with same sign
          IF (diff.GT.eps) THEN
            diff = MIN(diff,SQRT(DOT_PRODUCT(NormVecCheck+SFdepoLayersGeo(iSFfix,2,:),NormVecCheck+SFdepoLayersGeo(iSFfix,2,:)))) !NormVecCheck + SFdepoLayersGeo(iSFfix,2,1) should be (/0,0,0/) when identical with opposite sign
            IF (diff.GT.eps) THEN !yes, the first diff could be checked again...
              CALL abort(__STAMP__&
                ,' SFdepoLayersBaseVectors are not perpendicular to Normal for Layer ',iSFfix,diff)
            END IF
          END IF
        ELSE
          CALL abort(__STAMP__&
            ,' SFdepoLayersBaseVectors are parallel for Layer ',iSFfix)
        END IF
        !--calculate area and check if BV1/2 are orthgonal (need to be for domain reduction with BoundPoints, otherwise should not matter)
        SFdepoLayersCross(1) = SFdepoLayersBaseVector(iSFfix,1,2)*SFdepoLayersBaseVector(iSFfix,2,3) &
                             - SFdepoLayersBaseVector(iSFfix,1,3)*SFdepoLayersBaseVector(iSFfix,2,2)
        SFdepoLayersCross(2) = SFdepoLayersBaseVector(iSFfix,1,3)*SFdepoLayersBaseVector(iSFfix,2,1) &
                             - SFdepoLayersBaseVector(iSFfix,1,1)*SFdepoLayersBaseVector(iSFfix,2,3)
        SFdepoLayersCross(3) = SFdepoLayersBaseVector(iSFfix,1,1)*SFdepoLayersBaseVector(iSFfix,2,2) &
                             - SFdepoLayersBaseVector(iSFfix,1,2)*SFdepoLayersBaseVector(iSFfix,2,1)
        SFdepoLayersPartNum(iSFfix)=ABS(DOT_PRODUCT(SFdepoLayersCross,SFdepoLayersGeo(iSFfix,2,1:3))) !area of parallelogram
        BVlengths(1)=SQRT(DOT_PRODUCT(SFdepoLayersBaseVector(iSFfix,1,:),SFdepoLayersBaseVector(iSFfix,1,:)))
        BVlengths(2)=SQRT(DOT_PRODUCT(SFdepoLayersBaseVector(iSFfix,2,:),SFdepoLayersBaseVector(iSFfix,2,:)))
        IF (.NOT.ALMOSTEQUAL( SFdepoLayersPartNum(iSFfix),BVlengths(1)*BVlengths(2) )) THEN !area of assumed rectangle
          CALL abort(__STAMP__&
            ,' SFdepoLayersBaseVectors are not orthogonal for Layer ',iSFfix)
        END IF
        !--domain reduction based on BoundPoints
        BoundPoints(:,1) = (/SFdepoLayersBounds(iSFfix,1,1),SFdepoLayersBounds(iSFfix,1,2),SFdepoLayersBounds(iSFfix,1,3)/)
        BoundPoints(:,2) = (/SFdepoLayersBounds(iSFfix,2,1),SFdepoLayersBounds(iSFfix,1,2),SFdepoLayersBounds(iSFfix,1,3)/)
        BoundPoints(:,3) = (/SFdepoLayersBounds(iSFfix,1,1),SFdepoLayersBounds(iSFfix,2,2),SFdepoLayersBounds(iSFfix,1,3)/)
        BoundPoints(:,4) = (/SFdepoLayersBounds(iSFfix,2,1),SFdepoLayersBounds(iSFfix,2,2),SFdepoLayersBounds(iSFfix,1,3)/)
        BoundPoints(:,5) = (/SFdepoLayersBounds(iSFfix,1,1),SFdepoLayersBounds(iSFfix,1,2),SFdepoLayersBounds(iSFfix,2,3)/)
        BoundPoints(:,6) = (/SFdepoLayersBounds(iSFfix,2,1),SFdepoLayersBounds(iSFfix,1,2),SFdepoLayersBounds(iSFfix,2,3)/)
        BoundPoints(:,7) = (/SFdepoLayersBounds(iSFfix,1,1),SFdepoLayersBounds(iSFfix,2,2),SFdepoLayersBounds(iSFfix,2,3)/)
        BoundPoints(:,8) = (/SFdepoLayersBounds(iSFfix,2,1),SFdepoLayersBounds(iSFfix,2,2),SFdepoLayersBounds(iSFfix,2,3)/)
        !-1: shift BasePoint if applicable
        Dist=HUGE(1.)
        DO iPoint=1,8
          DistVec = BoundPoints(:,iPoint) - SFdepoLayersGeo(iSFfix,1,:) !vec from Basepoint to BoundPoint
          DO iVec=1,2
            Dist(iVec) = MIN(Dist(iVec),MAX(0.,DOT_PRODUCT(DistVec,SFdepoLayersBaseVector(iSFfix,iVec,:))/BVlengths(iVec))) !DistVec projected on BaseVec, if > 0
          END DO
          iVec=3
          Dist(iVec) = MIN(Dist(iVec),MAX(0.,DOT_PRODUCT(DistVec,SFdepoLayersGeo(iSFfix,2,:)))) !smallest DistVec projected on (already normalized) NormVec, if > 0
        END DO
        DO iVec=1,3 !shift of BP is independently for the 3 vectors possible because they are orthogonal!
          IF ( Dist(iVec).LT.HUGE(1.) .AND. Dist(iVec).GT.0. ) THEN !not-initial and positive dist was found (negative would not be reduction of domain)
            ChangeOccured = .TRUE.
            IF (iVec.NE.3) THEN
              IF (Dist(iVec).GE.BVlengths(iVec)) THEN !new BVlength would be < 0
                LayerOutsideOfBounds = .TRUE.
                EXIT
              ELSE !shift BP and reduce BVlength so that end of domain is at same position
                SFdepoLayersGeo(iSFfix,1,:) &
                  = SFdepoLayersGeo(iSFfix,1,:) + Dist(iVec)/BVlengths(iVec)*SFdepoLayersBaseVector(iSFfix,iVec,:)
                SFdepoLayersBaseVector(iSFfix,iVec,:) &
                  = SFdepoLayersBaseVector(iSFfix,iVec,:) * (BVlengths(iVec)-Dist(iVec))/BVlengths(iVec)
                BVlengths(iVec)=BVlengths(iVec)-Dist(iVec)
              END IF
            ELSE !NormalVec (r_SF is fix)
              IF (Dist(iVec).GE.r_SF) THEN
                LayerOutsideOfBounds = .TRUE.
                EXIT
              END IF
            END IF
          END IF
        END DO
        !-2: shorten BaseVectors if applicable (r_SF is fix)
        IF (.NOT.LayerOutsideOfBounds) THEN
          Dist=-HUGE(1.)
          DO iPoint=1,8
            DistVec = BoundPoints(:,iPoint) - SFdepoLayersGeo(iSFfix,1,:) !vec from Basepoint to BoundPoint
            DO iVec=1,2
              Dist(iVec) = MAX(Dist(iVec),DOT_PRODUCT(DistVec,SFdepoLayersBaseVector(iSFfix,iVec,:))/BVlengths(iVec)) !largest DistVec projected on BaseVec
            END DO
            iVec=3
            Dist(iVec) = MAX(Dist(iVec),DOT_PRODUCT(DistVec,SFdepoLayersGeo(iSFfix,2,:)))
          END DO
          DO iVec=1,3
            IF ( Dist(iVec).LE.0.) THEN !completely outside
              LayerOutsideOfBounds = .TRUE.
              EXIT
            ELSE IF ( iVec.NE.3 ) THEN
              IF (Dist(iVec).LT.BVlengths(iVec)) THEN !positive dist < BVlength was found (> BVlength would not be reduction of domain)
                ChangeOccured = .TRUE.
                SFdepoLayersBaseVector(iSFfix,iVec,:) = SFdepoLayersBaseVector(iSFfix,iVec,:) * (Dist(iVec))/BVlengths(iVec)
                BVlengths(iVec)=Dist(iVec)
              END IF
            END IF
          END DO
        END IF
        SFdepoLayersPartNum(iSFfix)=BVlengths(1)*BVlengths(2)
      CASE('cylinder')
        !calc BaseVectors from Normalvector
        IF (SFdepoLayersGeo(iSFfix,2,3).NE.0) THEN
          BaseVector(1) = 1.0
          BaseVector(2) = 1.0
          BaseVector(3) = -(SFdepoLayersGeo(iSFfix,2,1)+SFdepoLayersGeo(iSFfix,2,2))/ &
            SFdepoLayersGeo(iSFfix,2,3)
        ELSE
          IF (SFdepoLayersGeo(iSFfix,2,2).NE.0) THEN
            BaseVector(1) = 1.0
            BaseVector(3) = 1.0
            BaseVector(2) = -(SFdepoLayersGeo(iSFfix,2,1)+SFdepoLayersGeo(iSFfix,2,3))/ &
              SFdepoLayersGeo(iSFfix,2,2)
          ELSE
            IF (SFdepoLayersGeo(iSFfix,2,1).NE.0) THEN
              BaseVector(2) = 1.0
              BaseVector(3) = 1.0
              BaseVector(1) = -(SFdepoLayersGeo(iSFfix,2,2)+SFdepoLayersGeo(iSFfix,2,3))/ &
                SFdepoLayersGeo(iSFfix,2,1)
            END IF
          END IF
        END IF
        BaseVector = BaseVector / SQRT(BaseVector(1) * BaseVector(1) + BaseVector(2) * &
          BaseVector(2) + BaseVector(3) * BaseVector(3))
        BaseVector2(1) = SFdepoLayersGeo(iSFfix,2,2) * BaseVector(3) - &
          SFdepoLayersGeo(iSFfix,2,3) * BaseVector(2)
        BaseVector2(2) = SFdepoLayersGeo(iSFfix,2,3) * BaseVector(1) - &
          SFdepoLayersGeo(iSFfix,2,1) * BaseVector(3)
        BaseVector2(3) = SFdepoLayersGeo(iSFfix,2,1) * BaseVector(2) - &
          SFdepoLayersGeo(iSFfix,2,2) * BaseVector(1)
        BaseVector2 = BaseVector2 / SQRT(BaseVector2(1) * BaseVector2(1) + BaseVector2(2) * &
          BaseVector2(2) + BaseVector2(3) * BaseVector2(3))
        SFdepoLayersBaseVector(iSFfix,1,1:3)=BaseVector
        SFdepoLayersBaseVector(iSFfix,2,1:3)=BaseVector2
        SFdepoLayersRadius(iSFfix) = &
          GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-SFdepoLayersRadius','1.')
        SFdepoLayersPartNum(iSFfix)=PI*SFdepoLayersRadius(iSFfix)*SFdepoLayersRadius(iSFfix) !volume scaled with height (r_sf)
      CASE DEFAULT
        CALL abort(__STAMP__, &
          ' Wrong Space for SFdepoLayer: only cuboid and cylinder implemented!')
      END SELECT
      SFdepoLayersChargedens = GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-Chargedens','1.')
      SFdepoLayersSpec(iSFFix) = GETINT('PIC-SFdepoLayers'//TRIM(hilf)//'-Spec','1')
      IF (.NOT.LayerOutsideOfBounds) THEN
        SFdepoLayersGeo(iSFfix,2,:) = SFdepoLayersGeo(iSFfix,2,:)*r_sf
        IF (SFdepoLayersPartNum(iSFfix).LE.0.) CALL abort(__STAMP__&
          ,' Volume of SFdepoLayersXX is zero for Layer ',iSFfix)
        SFdepoLayersPartNum(iSFfix) = SFdepoLayersPartNum(iSFfix)*r_sf &
          * SFdepoLayersChargedens/Species(SFdepoLayersSpec(iSFFix))%MacroParticleFactor
        SWRITE(*,'(E12.5,A,I0)') SFdepoLayersPartNum(iSFfix), &
          ' additional particles will be inserted for SFdepoLayer ',iSFfix
        IF (ChangeOccured) THEN
          IF(PrintSFDepoWarnings)THEN
            IPWRITE(*,'(I4,A,I0,A)') ' WARNING: SFdepoLayer ',iSFfix,' was changed!'
            IPWRITE(*,'(I4,A,3(x,E12.5))')                 '          xyz minBounds:',SFdepoLayersBounds(iSFfix,1,:)
            IPWRITE(*,'(I4,A,3(x,E12.5))')                 '          xyz maxBounds:',SFdepoLayersBounds(iSFfix,2,:)
            IPWRITE(*,'(I4,A,3(x,E12.5))')                 '          New Basepoint:',SFdepoLayersGeo(iSFfix,1,:)
            IPWRITE(*,'(I4,A,3(x,E12.5))')                 '          New BaseVec01:',SFdepoLayersBaseVector(iSFfix,1,:)
            IPWRITE(*,'(I4,A,3(x,E12.5))')                 '          New BaseVec02:',SFdepoLayersBaseVector(iSFfix,2,:)
          END IF
        END IF
      ELSE !LayerOutsideOfBounds
        SFdepoLayersPartNum(iSFfix)=0
        IF(PrintSFDepoWarnings)THEN
          IPWRITE(*,'(I4,A,I0,A)') ' WARNING: SFdepoLayer ',iSFfix,' was disabled!'
          IPWRITE(*,'(I4,A,3(x,E12.5))')                 '          xyz minBounds:',SFdepoLayersBounds(iSFfix,1,:)
          IPWRITE(*,'(I4,A,3(x,E12.5))')                 '          xyz maxBounds:',SFdepoLayersBounds(iSFfix,2,:)
        END IF
      END IF
    END DO
  END IF !NbrOfSFdepoLayers>0

  !-- ResampleAnalyzeSurfCollis
  SFResampleAnalyzeSurfCollis = GETLOGICAL('PIC-SFResampleAnalyzeSurfCollis','.FALSE.')
  IF (SFResampleAnalyzeSurfCollis) THEN
    LastAnalyzeSurfCollis%PartNumberSamp = 0
    LastAnalyzeSurfCollis%PartNumberDepo = 0
    LastAnalyzeSurfCollis%ReducePartNumber = GETLOGICAL('PIC-SFResampleReducePartNumber','.FALSE.')
    LastAnalyzeSurfCollis%PartNumThreshold = GETINT('PIC-PartNumThreshold','0')
    IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN
      WRITE(UNIT=hilf,FMT='(I0)') LastAnalyzeSurfCollis%PartNumThreshold
      LastAnalyzeSurfCollis%PartNumberReduced = GETINT('PIC-SFResamplePartNumberReduced',TRIM(hilf)) !def. PartNumThreshold
    END IF
    WRITE(UNIT=hilf,FMT='(E16.8)') -HUGE(1.0)
    LastAnalyzeSurfCollis%Bounds(1,1)   = MAX(GETREAL('PIC-SFResample-xmin',TRIM(hilf)),GEO%xmin-r_sf)
    LastAnalyzeSurfCollis%Bounds(1,2)   = MAX(GETREAL('PIC-SFResample-ymin',TRIM(hilf)),GEO%ymin-r_sf)
    LastAnalyzeSurfCollis%Bounds(1,3)   = MAX(GETREAL('PIC-SFResample-zmin',TRIM(hilf)),GEO%zmin-r_sf)
    WRITE(UNIT=hilf,FMT='(E16.8)') HUGE(1.0)
    LastAnalyzeSurfCollis%Bounds(2,1)   = MIN(GETREAL('PIC-SFResample-xmax',TRIM(hilf)),GEO%xmax+r_sf)
    LastAnalyzeSurfCollis%Bounds(2,2)   = MIN(GETREAL('PIC-SFResample-ymax',TRIM(hilf)),GEO%ymax+r_sf)
    LastAnalyzeSurfCollis%Bounds(2,3)   = MIN(GETREAL('PIC-SFResample-zmax',TRIM(hilf)),GEO%zmax+r_sf)
    IF (NbrOfSFdepoFixes.EQ.0) THEN
      LastAnalyzeSurfCollis%UseFixBounds = .FALSE.
    ELSE
      LastAnalyzeSurfCollis%UseFixBounds = GETLOGICAL('PIC-SFResample-UseFixBounds','.TRUE.')
    END IF
    LastAnalyzeSurfCollis%NormVecOfWall = GETREALARRAY('PIC-NormVecOfWall',3,'1. , 0. , 0.')  !directed outwards
    IF (DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%NormVecOfWall).GT.0.) THEN
      LastAnalyzeSurfCollis%NormVecOfWall = LastAnalyzeSurfCollis%NormVecOfWall &
        / SQRT( DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%NormVecOfWall) )
    END IF
    LastAnalyzeSurfCollis%Restart = GETLOGICAL('PIC-SFResampleRestart','.FALSE.')
    IF (LastAnalyzeSurfCollis%Restart) THEN
      LastAnalyzeSurfCollis%DSMCSurfCollisRestartFile = GETSTR('PIC-SFResampleRestartFile','dummy')
    END IF
    LastAnalyzeSurfCollis%NumberOfBCs = GETINT('PIC-SFResampleNumberOfBCs','1')
    ALLOCATE(LastAnalyzeSurfCollis%BCs(1:LastAnalyzeSurfCollis%NumberOfBCs))
    IF (LastAnalyzeSurfCollis%NumberOfBCs.EQ.1) THEN !already allocated
      LastAnalyzeSurfCollis%BCs = GETINT('PIC-SFResampleSurfCollisBC','0') ! 0 means all...
    ELSE     
      hilf2=''
      DO iBC=1,LastAnalyzeSurfCollis%NumberOfBCs !build default string: 0,0,0,...
        WRITE(UNIT=hilf,FMT='(I0)') 0
        hilf2=TRIM(hilf2)//TRIM(hilf)
        IF (iBC.NE.LastAnalyzeSurfCollis%NumberOfBCs) hilf2=TRIM(hilf2)//','
      END DO
      LastAnalyzeSurfCollis%BCs = GETINTARRAY('PIC-SFResampleSurfCollisBC',LastAnalyzeSurfCollis%NumberOfBCs,hilf2)
    END IF
  END IF

  VolumeShapeFunction=4./3.*PI*r_sf*r2_sf
  nTotalDOF=REAL(nGlobalElems)*REAL((PP_N+1)**3)
  IF(MPIRoot)THEN
    IF(VolumeShapeFunction.GT.GEO%MeshVolume) &
      CALL abort(&
      __STAMP__&
      ,'ShapeFunctionVolume > MeshVolume')
  END IF
  
  SWRITE(UNIT_stdOut,'(A,F12.6)') ' | Average DOFs in Shape-Function |                      ' , &
      nTotalDOF*VolumeShapeFunction/GEO%MeshVolume

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1. 
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

#if USE_MPI
  DoExternalParts=.TRUE.
#endif /*MPI*/

CASE('shape_function_1d')
  r_sf     = GETREAL('PIC-shapefunction-radius','1.')
  alpha_sf = GETINT ('PIC-shapefunction-alpha','2')
  sf1d_dir = GETINT ('PIC-shapefunction1d-direction','1')
  DoSFEqui = GETLOGICAL('PIC-shapefunction-equi','F')
  r2_sf = r_sf * r_sf 
  r2_sf_inv = 1./r2_sf
  SELECT CASE(alpha_sf)
  CASE(2)
    w_sf=16.*r_sf/15.
  CASE(3)
    w_sf=32.*r_sf/35.
  CASE(4)
    w_sf=256.*r_sf/315.
  CASE(5)
    w_sf=512.*r_sf/693.
  CASE(6)
    w_sf=2048.*r_sf/3003.
  CASE(7)
    w_sf=4096.*r_sf/6435.
  CASE(8)
    w_sf=65536.*r_sf/109395.
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,' Correct 1D weight not precomputed!')
  END SELECT

  IF(sf1d_dir.EQ.1)THEN
    w_sf=w_sf*(GEO%ymaxglob-GEO%yminglob)*(GEO%zmaxglob-GEO%zminglob)
  ELSE IF (sf1d_dir.EQ.2)THEN
    w_sf=w_sf*(GEO%xmaxglob-GEO%xminglob)*(GEO%zmaxglob-GEO%zminglob)
  ELSE IF (sf1d_dir.EQ.3)THEN
    w_sf=w_sf*(GEO%xmaxglob-GEO%xminglob)*(GEO%ymaxglob-GEO%yminglob)
  ELSE
    w_sf=2*GEO%MeshVolume
  END IF
  VolumeShapeFunction=w_sf
  w_sf=1.0/w_sf
  nTotalDOF=REAL(nGlobalElems)*REAL((PP_N+1)**3)
  IF(MPIRoot)THEN
    IF(VolumeShapeFunction.GT.GEO%MeshVolume) &
      CALL abort(&
      __STAMP__&
      ,'ShapeFunctionVolume > MeshVolume')
  END IF
  
  SWRITE(UNIT_stdOut,'(A,F12.6)') ' | Shape function volume          |                      ' , &
      VolumeShapeFunction
  SWRITE(UNIT_stdOut,'(A,F12.6)') ' | Average DOFs in Shape-Function |                      ' , &
      nTotalDOF*VolumeShapeFunction/GEO%MeshVolume

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1. 
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

#if USE_MPI
  DoExternalParts=.TRUE.
#endif /*MPI*/

CASE('shape_function_cylindrical','shape_function_spherical')
  IF(TRIM(DepositionType).EQ.'shape_function_cylindrical')THEN
    SfRadiusInt=2
  ELSE
    SfRadiusInt=3
  END IF
  r_sf       = GETREAL   ('PIC-shapefunction-radius','1.')
  r_sf0      = GETREAL   ('PIC-shapefunction-radius0','1.')
  r_sf_scale = GETREAL   ('PIC-shapefunction-scale','0.')
  alpha_sf   = GETINT    ('PIC-shapefunction-alpha','2')
  DoSFEqui   = GETLOGICAL('PIC-shapefunction-equi','F')
  BetaFac    = beta(1.5, REAL(alpha_sf) + 1.)
  w_sf = 1./(2. * BetaFac * REAL(alpha_sf) + 2 * BetaFac) &
                        * (REAL(alpha_sf) + 1.)!/(PI*(r_sf**3))
  r2_sf = r_sf * r_sf 
  r2_sf_inv = 1./r2_sf

  r_sf_average=0.5*(r_sf+r_sf0)
  VolumeShapeFunction=4./3.*PI*r_sf_average*r_sf_average
  nTotalDOF=REAL(nGlobalElems)*REAL((PP_N+1)**3)
  IF(MPIRoot)THEN
    IF(VolumeShapeFunction.GT.GEO%MeshVolume) &
      CALL abort(&
      __STAMP__&
      ,'ShapeFunctionVolume > MeshVolume')
  END IF
  
  SWRITE(UNIT_stdOut,'(A,F12.6)') ' | Average DOFs in Shape-Function |                      ' , &
      nTotalDOF*VolumeShapeFunction/GEO%MeshVolume

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1. 
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

#if USE_MPI
  DoExternalParts=.TRUE.
#endif /*MPI*/

CASE('delta_distri')
  ! Allocate array for particle positions in -1|1 space (used for deposition as well as interpolation)
  IF(.NOT.DoRefMapping)THEN
    ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,' Cannot allocate partposref!')
  END IF
  DeltaType = GETINT('PIC-DeltaType','1')
  WRITE(hilf,'(I0)') PP_N
  NDepo     = GETINT('PIC-DeltaType-N',hilf)
  SELECT CASE(DeltaType)
  CASE(1)
    SWRITE(UNIT_stdOut,'(A)') ' Lagrange-Polynomial'
  CASE(2)
    SWRITE(UNIT_stdOut,'(A)') ' Bernstein-Polynomial'
    !CALL BuildBernSteinVdm(PP_N,xGP)
    ALLOCATE(NDepochooseK(0:NDepo,0:NDepo))
    CALL ComputeBernsteinCoeff(NDepo,NDepoChooseK)
  CASE(3)
    SWRITE(UNIT_stdOut,'(A)') ' Uniform B-Spline '
    NKnots=(NDepo+1)*2-1
    ALLOCATE( Knots(0:NKnots)  )
    ! knots distance is 2 [-1,1)
    X0=-1.-REAL(NDepo)*2.0
    DO i=0,nKnots
      Knots(i) =X0+2.0*REAL(i)
    END DO ! i
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,' Wrong Delta-Type')
  END SELECT
  IF(NDepo.GT.PP_N) CALL abort(&
    __STAMP__&
    ,' NDepo must be smaller than N!')
  DeltaDistriChangeBasis=.FALSE.
  IF(NDepo.NE.PP_N) DeltaDistriChangeBasis=.TRUE.
  ALLOCATE(DDMassInv(0:NDepo,0:NDepo,0:NDepo,1:PP_nElems) &
          , XiNDepo(0:NDepo)                              &
          , swGPNDepo(0:NDepo)                            &
          , wBaryNDepo(0:NDepo)                           )
  DDMassInv=0.
#if (PP_NodeType==1)
  CALL LegendreGaussNodesAndWeights(NDepo,XiNDepo,swGPNDepo)
#elif (PP_NodeType==2)
  CALL LegGaussLobNodesAndWeights(NDepo,XiNDepo,swGPNDepo)
#endif
  CALL BarycentricWeights(NDepo,XiNDepo,wBaryNDepo)

  IF(DeltaDistriChangeBasis)THEN
    ALLOCATE( Vdm_NDepo_GaussN(0:PP_N,0:NDepo)             &
            , Vdm_GaussN_NDepo(0:NDepo,0:PP_N)             &
            , dummy(1,0:PP_N,0:PP_N,0:PP_N)                &
            , dummy2(1,0:NDepo,0:NDepo,0:NDepo)            )
    CALL InitializeVandermonde(NDepo,PP_N,wBaryNDepo,XiNDepo,xGP,Vdm_NDepo_GaussN)
    ! and inverse of mass matrix
    DO i=0,NDepo
      swGPNDepo(i)=1.0/swGPNDepo(i)
    END DO ! i=0,PP_N
    CALL InitializeVandermonde(PP_N,NDepo,wBary,xGP,xiNDepo,Vdm_GaussN_NDepo)
    DO iElem=1,PP_nElems
      dummy(1,:,:,:)=sJ(:,:,:,iElem,0)
      CALL ChangeBasis3D(1,PP_N,NDepo,Vdm_GaussN_NDepo,dummy,dummy2)
      DO k=0,NDepo
        DO j=0,NDepo
          DO i=0,NDepo
            DDMassInv(i,j,k,iElem)=dummy2(1,i,j,k)*swGPNDepo(i)*swGPNDepo(j)*swGPNDepo(k)
          END DO ! i
        END DO ! j
      END DO ! k
    END DO ! iElem=1,PP_nElems
    !DEALLOCATE(Vdm_NDepo_GaussN)
    DEALLOCATE(Vdm_GaussN_NDepo)
    DEALLOCATE(dummy,dummy2)
  ELSE
    DO i=0,NDepo
      swGPNDepo(i)=1.0/wGP(i)
    END DO ! i=0,PP_N
    DO iElem=1,PP_nElems
      DO k=0,NDepo
        DO j=0,NDepo
          DO i=0,NDepo
            DDMassInv(i,j,k,iElem)=sJ(i,j,k,iElem,0)*swGPNDepo(i)*swGPNDepo(j)*swGPNDepo(k)
          END DO ! i
        END DO ! j
      END DO ! k
    END DO ! iElem=1,PP_nElems
  END IF
  DEALLOCATE(swGPNDepo)

CASE('cartmesh_volumeweighting')
    
  ! read in background mesh size
  BGMdeltas(1:3) = GETREALARRAY('PIC-BGMdeltas',3,'0. , 0. , 0.')
  FactorBGM(1:3) = GETREALARRAY('PIC-FactorBGM',3,'1. , 1. , 1.')
  BGMdeltas(1:3) = 1./FactorBGM(1:3)*BGMdeltas(1:3)
  IF (ANY(BGMdeltas.EQ.0.0)) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR: PIC-BGMdeltas: No size for the cartesian background mesh definded.')
  END IF
  ! calc min and max coordinates for local mesh
  xmin = 1.0E200
  ymin = 1.0E200
  zmin = 1.0E200
  xmax = -1.0E200
  ymax = -1.0E200
  zmax = -1.0E200
  DO iElem = 1, nElems
    xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
    xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
    ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
    ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
    zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
    zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  END DO
  ! define minimum and maximum backgroundmesh index, compute volume
  BGMVolume = BGMdeltas(1)*BGMdeltas(2)*BGMdeltas(3)
  BGMminX = FLOOR(xmin/BGMdeltas(1)-0.0001)
  BGMminY = FLOOR(ymin/BGMdeltas(2)-0.0001)
  BGMminZ = FLOOR(zmin/BGMdeltas(3)-0.0001)
  BGMmaxX = CEILING(xmax/BGMdeltas(1)+0.0001)
  BGMmaxY = CEILING(ymax/BGMdeltas(2)+0.0001)
  BGMmaxZ = CEILING(zmax/BGMdeltas(3)+0.0001)
  ! mapping from gausspoints to BGM
  ALLOCATE(GaussBGMIndex(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in pic_depo.f90: Cannot allocate GaussBGMIndex!') 
  END IF
  ALLOCATE(GaussBGMFactor(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in pic_depo.f90: Cannot allocate GaussBGMFactor!')
  END IF
  DO iElem = 1, nElems
    DO j = 0, PP_N
      DO k = 0, PP_N
        DO m = 0, PP_N
          GaussBGMIndex(1,j,k,m,iElem) = FLOOR(Elem_xGP(1,j,k,m,iElem)/BGMdeltas(1))
          GaussBGMIndex(2,j,k,m,iElem) = FLOOR(Elem_xGP(2,j,k,m,iElem)/BGMdeltas(2))
          GaussBGMIndex(3,j,k,m,iElem) = FLOOR(Elem_xGP(3,j,k,m,iElem)/BGMdeltas(3))
          GaussBGMFactor(1,j,k,m,iElem) = (Elem_xGP(1,j,k,m,iElem)/BGMdeltas(1))-REAL(GaussBGMIndex(1,j,k,m,iElem))
          GaussBGMFactor(2,j,k,m,iElem) = (Elem_xGP(2,j,k,m,iElem)/BGMdeltas(2))-REAL(GaussBGMIndex(2,j,k,m,iElem))
          GaussBGMFactor(3,j,k,m,iElem) = (Elem_xGP(3,j,k,m,iElem)/BGMdeltas(3))-REAL(GaussBGMIndex(3,j,k,m,iElem))
        END DO
      END DO
    END DO
  END DO
#if USE_MPI
  CALL MPIBackgroundMeshInit()
#else
  IF(GEO%nPeriodicVectors.GT.0)THEN
    ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
    ALLOCATE(GEO%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
    IF (allocStat .NE. 0) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR in MPIBackgroundMeshInit: cannot allocate GEO%PeriodicBGMVectors!')
    END IF
    DO i = 1, GEO%nPeriodicVectors
      GEO%PeriodicBGMVectors(1,i) = NINT(GEO%PeriodicVectors(1,i)/BGMdeltas(1))
      IF(ABS(GEO%PeriodicVectors(1,i)/BGMdeltas(1)-REAL(GEO%PeriodicBGMVectors(1,i))).GT.1E-10)THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(2,i) = NINT(GEO%PeriodicVectors(2,i)/BGMdeltas(2))
      IF(ABS(GEO%PeriodicVectors(2,i)/BGMdeltas(2)-REAL(GEO%PeriodicBGMVectors(2,i))).GT.1E-10)THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(3,i) = NINT(GEO%PeriodicVectors(3,i)/BGMdeltas(3))
      IF(ABS(GEO%PeriodicVectors(3,i)/BGMdeltas(3)-REAL(GEO%PeriodicBGMVectors(3,i))).GT.1E-10)THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
    END DO
  END IF
#endif
  ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
CASE('cartmesh_splines')
  BGMdeltas(1:3) = GETREALARRAY('PIC-BGMdeltas',3,'0. , 0. , 0.')
  FactorBGM(1:3) = GETREALARRAY('PIC-FactorBGM',3,'1. , 1. , 1.')
  BGMdeltas(1:3) = 1./FactorBGM(1:3)*BGMdeltas(1:3)
  IF (ANY(BGMdeltas.EQ.0)) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR: PIC-BGMdeltas: No size for the cartesian background mesh definded.')
  END IF
  ! calc min and max coordinates for local mesh
  xmin = 1.0E200
  ymin = 1.0E200
  zmin = 1.0E200
  xmax = -1.0E200
  ymax = -1.0E200
  zmax = -1.0E200
  DO iElem = 1, nElems
    xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
    xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
    ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
    ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
    zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
    zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  END DO
  ! define minimum and maximum backgroundmesh index, compute volume
  BGMVolume = BGMdeltas(1)*BGMdeltas(2)*BGMdeltas(3)
  BGMminX = FLOOR(xmin/BGMdeltas(1)-0.0001)-2
  BGMminY = FLOOR(ymin/BGMdeltas(2)-0.0001)-2
  BGMminZ = FLOOR(zmin/BGMdeltas(3)-0.0001)-2
  BGMmaxX = CEILING(xmax/BGMdeltas(1)+0.0001)+2
  BGMmaxY = CEILING(ymax/BGMdeltas(2)+0.0001)+2
  BGMmaxZ = CEILING(zmax/BGMdeltas(3)+0.0001)+2
  ! mapping from gausspoints to BGM
  ALLOCATE(GaussBGMIndex(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in pic_depo.f90: Cannot allocate GaussBGMIndex!')
  END IF
  ALLOCATE(GaussBGMFactor(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in pic_depo.f90: Cannot allocate GaussBGMFactor!')
  END IF
  DO iElem = 1, nElems
    DO j = 0, PP_N
      DO k = 0, PP_N
        DO m = 0, PP_N
          GaussBGMIndex(1,j,k,m,iElem) = FLOOR(Elem_xGP(1,j,k,m,iElem)/BGMdeltas(1))
          GaussBGMIndex(2,j,k,m,iElem) = FLOOR(Elem_xGP(2,j,k,m,iElem)/BGMdeltas(2))
          GaussBGMIndex(3,j,k,m,iElem) = FLOOR(Elem_xGP(3,j,k,m,iElem)/BGMdeltas(3))
        END DO
      END DO
    END DO
  END DO
  ! pre-build weights for BGM to gausspoint interpolation
  ALLOCATE(GPWeight(1:nElems,0:PP_N,0:PP_N,0:PP_N,1:4,1:4,1:4),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in pic_depo.f90: Cannot allocate GPWeight!')
  END IF
  DO iElem = 1, nElems
    DO j = 0, PP_N
      DO k = 0, PP_N
        DO m = 0, PP_N
          DO dir = 1,3               ! x,y,z direction
            DO weightrun = 0,3
              DO mm = 0, 3
                IF (mm.EQ.weightrun) then
                  auxiliary(mm) = 1.0
                ELSE
                  auxiliary(mm) = 0.0
                END IF
              END DO
              CALL DeBoor(GaussBGMIndex(dir,j,k,m,iElem),auxiliary,Elem_xGP(dir,j,k,m,iElem),weight(dir,weightrun),dir)
            END DO
          END DO
          DO r = 1,4
            DO s = 1,4
              DO t = 1,4
                GPWeight(iElem,j,k,m,r,s,t) = weight(1,-(r-4)) * weight(2,-(s-4)) * weight(3,-(t-4))
              END DO !t
            END DO !s
          END DO !r
        END DO !m
      END DO !k
    END DO !j
  END DO !iElem
#if USE_MPI
  CALL MPIBackgroundMeshInit()
#else
  IF(GEO%nPeriodicVectors.GT.0)THEN
    ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
    ALLOCATE(GEO%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
    IF (allocStat .NE. 0) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR in MPIBackgroundMeshInit: cannot allocate GEO%PeriodicBGMVectors!')
    END IF
    DO i = 1, GEO%nPeriodicVectors
      GEO%PeriodicBGMVectors(1,i) = NINT(GEO%PeriodicVectors(1,i)/BGMdeltas(1))
      IF(ABS(GEO%PeriodicVectors(1,i)/BGMdeltas(1)-REAL(GEO%PeriodicBGMVectors(1,i))).GT.1E-10)THEN
        CALL abort(&
        __STAMP__ &
        ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(2,i) = NINT(GEO%PeriodicVectors(2,i)/BGMdeltas(2))
      IF(ABS(GEO%PeriodicVectors(2,i)/BGMdeltas(2)-REAL(GEO%PeriodicBGMVectors(2,i))).GT.1E-10)THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(3,i) = NINT(GEO%PeriodicVectors(3,i)/BGMdeltas(3))
      IF(ABS(GEO%PeriodicVectors(3,i)/BGMdeltas(3)-REAL(GEO%PeriodicBGMVectors(3,i))).GT.1E-10)THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
    END DO
  END IF
#endif
  ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Unknown DepositionType in pic_depo.f90')
END SELECT

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DEPOSITION DONE!'

END SUBROUTINE InitializeDeposition


SUBROUTINE Deposition(doInnerParts,doParticle_In)
!============================================================================================================================
! This subroutine performes the deposition of the particle charge and current density to the grid
! following list of distribution methods are implemted
! - nearest blurrycenter (barycenter of hexahedra)
! - nearest Gauss Point  (only volome of IP - higher resolution than nearest blurrycenter )
! - shape function       (only one type implemented)
! - delta distributio
! useVMPF added, therefore, this routine contains automatically the use of variable mpfs
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_PreProc
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Mesh_Vars,              ONLY:nElems, Elem_xGP, sJ
USE MOD_ChangeBasis,            ONLY:ChangeBasis3D
USE MOD_Interpolation_Vars,     ONLY:wGP
USE MOD_PICInterpolation_Vars,  ONLY:InterpolationType
USE MOD_Eval_xyz,               ONLY:TensorProductInterpolation
USE MOD_Basis,                  ONLY:LagrangeInterpolationPolys
USE MOD_Particle_Basis,         ONLY:BernSteinPolynomial
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,casematrix, NbrOfCases
#if USE_MPI
USE MOD_Particle_MPI_Vars,      ONLY:ExtPartState,ExtPartSpecies,ExtPartToFIBGM,NbrOfExtParticles
USE MOD_Particle_MPI_Vars,      ONLY:PartMPIExchange
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
#endif  /*MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT variable declaration                                                                       
LOGICAL,INTENT(IN)               :: doInnerParts
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT variable declaration                                                                       
!-----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration                                                                       
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                          :: doParticle(1:PDM%ParticleVecLength)
INTEGER                          :: firstPart,lastPart
INTEGER                          :: i,j, k, l, m, iElem, iPart, iPart2, iSFfix
LOGICAL                          :: chargedone(1:nElems)        
LOGICAL                          :: SAVE_GAUSS
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax                           
INTEGER                          :: kk, ll, mm, ppp                                              
INTEGER                          :: ElemID, iCase, ind
REAL                             :: radius2, S, S1, Fac(4)
REAL                             :: dx,dy,dz
REAL, ALLOCATABLE                :: BGMSourceCellVol(:,:,:,:,:), tempsource(:,:,:), tempgridsource(:)
REAL                             :: Vec1(1:3), Vec2(1:3), Vec3(1:3), ShiftedPart(1:3)
INTEGER                          :: a,b, ii, expo
REAL                             :: ElemSource(nElems,1:4)
REAL                             :: Charge, TSource(1:4), auxiliary(0:3),weight(1:3,0:3), locweight
REAL                             :: alpha1, alpha2, alpha3, TempPartPos(1:3), alpha
INTEGER                          :: PosInd(3),r,ss,t,u,v,w, dir, weightrun
INTEGER                          :: iLayer, layerParts
REAL,DIMENSION(3,0:NDepo)        :: L_xi
REAL                             :: DeltaIntCoeff,prefac
REAL                             :: local_r_sf, local_r2_sf, local_r2_sf_inv
REAL                             :: RandVal, RandVal2(2), layerPartPos(3), PartRadius, FractPush(3), SFfixDistance
LOGICAL                          :: DoCycle
!#if USE_MPI
! load balance
!REAL                             :: tLBStart,tLBEnd
!#endif /*MPI*/
!============================================================================================================================

IF(PRESENT(DoParticle_IN))THEN
  DoParticle=PDM%ParticleInside(1:PDM%ParticleVecLength).AND.DoParticle_In
ELSE
  DoParticle(1:PDM%ParticleVecLength)=PDM%ParticleInside(1:PDM%ParticleVecLength)
END IF

IF(doInnerParts)THEN
  IF(.NOT.DoDeposition) RETURN
  PartSource=0.0
  firstPart=1
  lastPart =PDM%ParticleVecLength
  !IF(firstPart.GT.lastPart) RETURN
ELSE
  IF(.NOT.DoDeposition) RETURN
#if USE_MPI
  firstPart=PDM%ParticleVecLength-PartMPIExchange%nMPIParticles+1
  lastPart =PDM%ParticleVecLength
#else
  firstPart=1
  lastPart =0
#endif /*MPI*/
END IF


SELECT CASE(TRIM(DepositionType))
    
CASE('nearest_blurrycenter')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  ElemSource=0.0
  DO iElem=1,PP_nElems
!#if USE_MPI
!    tLBStart = FLEXITIME() ! LB Time Start
!#endif /*MPI*/
    DO iPart=firstPart,lastPart
      !IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      IF(.NOT.DoParticle(iPart)) CYCLE
      IF(PEM%Element(iPart).EQ.iElem)THEN
         ElemSource(iElem,1:3) = ElemSource(iElem,1:3)+ &
                PartState(iPart,4:6)* Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
         ElemSource(iElem,4) = ElemSource(iElem,4) + & 
              Species(PartSpecies(iPart))%ChargeIC* Species(PartSpecies(iPart))%MacroParticleFactor
      END IF ! Element(iPart).EQ.iElem
    END DO ! iPart
    PartSource(1,:,:,:,iElem) = PartSource(1,:,:,:,iElem)+ElemSource(iElem,1) 
    PartSource(2,:,:,:,iElem) = PartSource(2,:,:,:,iElem)+ElemSource(iElem,2) 
    PartSource(3,:,:,:,iElem) = PartSource(3,:,:,:,iElem)+ElemSource(iElem,3)                                        
    PartSource(4,:,:,:,iElem) = PartSource(4,:,:,:,iElem)+ElemSource(iElem,4) 
  END DO ! iElem=1,PP_nElems
  IF(.NOT.doInnerParts)THEN
    DO iElem=1,PP_nElems
      PartSource(1:4,:,:,:,iElem) = PartSource(1:4,:,:,:,iElem) / GEO%Volume(iElem)
    END DO ! iElem=1,PP_nElems
  END IF ! .NOT. doInnerParts
  
CASE('cell_volweight') 
  ALLOCATE(BGMSourceCellVol(0:1,0:1,0:1,1:nElems,1:4))
  BGMSourceCellVol(:,:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
    !IF (PDM%ParticleInside(iPart)) THEN
    IF (DoParticle(iPart)) THEN
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
      iElem = PEM%Element(iPart)
      IF(DoRefMapping)THEN
        TempPartPos(1:3)=PartPosRef(1:3,iPart)
      ELSE
        CALL TensorProductInterpolation(PartState(iPart,1:3),TempPartPos,iElem,ForceMode=.TRUE.)
      END IF
      TSource(:) = 0.0
      TSource(1) = PartState(iPart,4)*Charge
      TSource(2) = PartState(iPart,5)*Charge
      TSource(3) = PartState(iPart,6)*Charge
      TSource(4) = Charge
      alpha1=(TempPartPos(1)+1.0)/2.0
      alpha2=(TempPartPos(2)+1.0)/2.0
      alpha3=(TempPartPos(3)+1.0)/2.0
      BGMSourceCellVol(0,0,0,iElem,1:4) = BGMSourceCellVol(0,0,0,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
      BGMSourceCellVol(0,0,1,iElem,1:4) = BGMSourceCellVol(0,0,1,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(alpha3))
      BGMSourceCellVol(0,1,0,iElem,1:4) = BGMSourceCellVol(0,1,0,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(alpha2)*(1-alpha3))
      BGMSourceCellVol(0,1,1,iElem,1:4) = BGMSourceCellVol(0,1,1,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(alpha2)*(alpha3))
      BGMSourceCellVol(1,0,0,iElem,1:4) = BGMSourceCellVol(1,0,0,iElem,1:4) + (TSource(1:4)*(alpha1)*(1-alpha2)*(1-alpha3))
      BGMSourceCellVol(1,0,1,iElem,1:4) = BGMSourceCellVol(1,0,1,iElem,1:4) + (TSource(1:4)*(alpha1)*(1-alpha2)*(alpha3))
      BGMSourceCellVol(1,1,0,iElem,1:4) = BGMSourceCellVol(1,1,0,iElem,1:4) + (TSource(1:4)*(alpha1)*(alpha2)*(1-alpha3))
      BGMSourceCellVol(1,1,1,iElem,1:4) = BGMSourceCellVol(1,1,1,iElem,1:4) + (TSource(1:4)*(alpha1)*(alpha2)*(alpha3))   
    END IF
  END DO

  DO iElem=1, nElems
    BGMSourceCellVol(0,0,0,iElem,:) = BGMSourceCellVol(0,0,0,iElem,1:4)/CellVolWeight_Volumes(0,0,0,iElem)
    BGMSourceCellVol(0,0,1,iElem,:) = BGMSourceCellVol(0,0,1,iElem,1:4)/CellVolWeight_Volumes(0,0,1,iElem)
    BGMSourceCellVol(0,1,0,iElem,:) = BGMSourceCellVol(0,1,0,iElem,1:4)/CellVolWeight_Volumes(0,1,0,iElem)
    BGMSourceCellVol(0,1,1,iElem,:) = BGMSourceCellVol(0,1,1,iElem,1:4)/CellVolWeight_Volumes(0,1,1,iElem)
    BGMSourceCellVol(1,0,0,iElem,:) = BGMSourceCellVol(1,0,0,iElem,1:4)/CellVolWeight_Volumes(1,0,0,iElem)
    BGMSourceCellVol(1,0,1,iElem,:) = BGMSourceCellVol(1,0,1,iElem,1:4)/CellVolWeight_Volumes(1,0,1,iElem)
    BGMSourceCellVol(1,1,0,iElem,:) = BGMSourceCellVol(1,1,0,iElem,1:4)/CellVolWeight_Volumes(1,1,0,iElem)
    BGMSourceCellVol(1,1,1,iElem,:) = BGMSourceCellVol(1,1,1,iElem,1:4)/CellVolWeight_Volumes(1,1,1,iElem)   
  END DO

  DO iElem = 1, nElems
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
         alpha1 = CellVolWeightFac(kk)
         alpha2 = CellVolWeightFac(ll)
         alpha3 = CellVolWeightFac(mm)
         PartSource(1:4,kk,ll,mm,iElem) =PartSource(1:4,kk,ll,mm,iElem) +&
              BGMSourceCellVol(0,0,0,iElem,1:4) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSourceCellVol(0,0,1,iElem,1:4) * (1-alpha1) * (1-alpha2) * (alpha3) + &
              BGMSourceCellVol(0,1,0,iElem,1:4) * (1-alpha1) * (alpha2) * (1-alpha3) + &
              BGMSourceCellVol(0,1,1,iElem,1:4) * (1-alpha1) * (alpha2) * (alpha3) + &
              BGMSourceCellVol(1,0,0,iElem,1:4) * (alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSourceCellVol(1,0,1,iElem,1:4) * (alpha1) * (1-alpha2) * (alpha3) + &
              BGMSourceCellVol(1,1,0,iElem,1:4) * (alpha1) * (alpha2) * (1-alpha3) + &
              BGMSourceCellVol(1,1,1,iElem,1:4) * (alpha1) * (alpha2) * (alpha3)
       END DO !mm
     END DO !ll
   END DO !kk
 END DO !iEle
 DEALLOCATE(BGMSourceCellVol)
 
CASE('epanechnikov') 
  ALLOCATE(tempsource(0:PP_N,0:PP_N,0:PP_N))
  IF(DoInnerParts)  tempcharge= 0.0
  DO iPart = firstPart, lastPart
    !IF (PDM%ParticleInside(iPart)) THEN
    IF (DoParticle(iPart)) THEN
      Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
      iElem = PEM%Element(iPart)
      alpha = 0.0 
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
            radius2 = (PartState(iPart,1) - Elem_xGP(1,kk,ll,mm,iElem)) * (PartState(iPart,1) - Elem_xGP(1,kk,ll,mm,iElem)) &
                    + (PartState(iPart,2) - Elem_xGP(2,kk,ll,mm,iElem)) * (PartState(iPart,2) - Elem_xGP(2,kk,ll,mm,iElem)) &
                    + (PartState(iPart,3) - Elem_xGP(3,kk,ll,mm,iElem)) * (PartState(iPart,3) - Elem_xGP(3,kk,ll,mm,iElem))
           IF (radius2 .LT. r2_sf) THEN
             tempsource(kk,ll,mm) = r2_sf - radius2
             alpha = alpha + tempsource(kk,ll,mm)               
           ELSE
             tempsource(kk,ll,mm) = 0.0
           END IF         
         END DO !mm
       END DO !ll
      END DO !kk
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
           PartSource(1:3,kk,ll,mm,iElem) = PartSource(1:3,kk,ll,mm,iElem)  + 1./alpha*tempsource(kk,ll,mm)*PartState(iPart,4:6) &
                                          * Species(PartSpecies(iPart))%ChargeIC &
                                          * Species(PartSpecies(iPart))%MacroParticleFactor
           PartSource(4,kk,ll,mm,iElem) = PartSource(4,kk,ll,mm,iElem)  + 1./alpha*tempsource(kk,ll,mm) & 
                                          * Species(PartSpecies(iPart))%ChargeIC &
                                          * Species(PartSpecies(iPart))%MacroParticleFactor
         END DO !mm
       END DO !ll
      END DO !kk
    END IF
  END DO

  IF(.NOT.DoInnerParts)THEN
    ALLOCATE(tempgridsource(1:nElems))
    tempgridsource= 0.0
    ! seemps to be finalize
    DO iElem=1,nElems
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
            alpha = wGP(kk)*wGP(ll)*wGP(mm)/sJ(kk,ll,mm,iElem,0)
            PartSource(1:4,kk,ll,mm,iElem) = 1./SQRT(alpha) * PartSource(1:4,kk,ll,mm,iElem)*(PP_N+1)**3/ (GEO%Volume(iElem))
            tempgridsource(iElem) = tempgridsource(iElem) + PartSource(4,kk,ll,mm,iElem)*alpha
         END DO !mm
       END DO !ll
      END DO !kk
      ! possible ABS???
      !IF (tempgridsource(iElem).GT.0.0) THEN
      IF (ABS(tempgridsource(iElem)).GT.0.0) THEN
        alpha = tempcharge(iElem)/tempgridsource(iElem)
        PartSource(1:4,:,:,:,iElem) = PartSource(1:4,:,:,:,iElem)*alpha
      END IF
    END DO
    DEALLOCATE(tempgridsource)
  END IF
  DEALLOCATE(tempsource)

CASE('shape_function','shape_function_simple')
  !-- "normal" particles
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
    DO iPart=firstPart,LastPart
      IF (DoParticle(iPart)) THEN
        CALL calcSfSource(4,Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor*w_sf &
          ,Vec1,Vec2,Vec3,PartState(iPart,1:3),iPart,PartVelo=PartState(iPart,4:6))
      END IF ! DoParticle
    END DO ! iPart
!  END IF ! usevMPF
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
#if USE_MPI
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF
#endif /*MPI*/

    !-- layer particles (only once, i.e., during call with .NOT.DoInnerParts)
    DO iLayer=1,NbrOfSFdepoLayers
      CALL RANDOM_NUMBER(RandVal)
      IF (SFdepoLayersPartNum(iLayer).GT.0.) THEN
        layerParts=INT(SFdepoLayersPartNum(iLayer)+RandVal)
      ELSE
        layerParts=0
      END IF
      DO iPart=1,layerParts
        SELECT CASE (TRIM(SFdepoLayersSpace(iLayer)))
        CASE('cuboid')
          CALL RANDOM_NUMBER(RandVal2)
          layerPartPos = SFdepoLayersGeo(iLayer,1,:) + &
            (RandVal2(1)*SFdepoLayersBaseVector(iLayer,1,:) + RandVal2(2)*SFdepoLayersBaseVector(iLayer,2,:))
        CASE('cylinder')
          PartRadius = SFdepoLayersRadius(iLayer) + 1
          DO WHILE(PartRadius.GT.SFdepoLayersRadius(iLayer))
            CALL RANDOM_NUMBER(RandVal2)
            RandVal2 = RandVal2 * 2. - 1.
            layerPartPos = SFdepoLayersGeo(iLayer,1,:) + SFdepoLayersRadius(iLayer) * &
              (RandVal2(1)*SFdepoLayersBaseVector(iLayer,1,:) + RandVal2(2)*SFdepoLayersBaseVector(iLayer,2,:))
            PartRadius = SQRT( (layerPartPos(1)-SFdepoLayersGeo(iLayer,1,1)) * &
              (layerPartPos(1)-SFdepoLayersGeo(iLayer,1,1)) + &
              (layerPartPos(2)-SFdepoLayersGeo(iLayer,1,2)) * &
              (layerPartPos(2)-SFdepoLayersGeo(iLayer,1,2)) + &
              (layerPartPos(3)-SFdepoLayersGeo(iLayer,1,3)) * &
              (layerPartPos(3)-SFdepoLayersGeo(iLayer,1,3)) )
          END DO
        CASE DEFAULT
          CALL abort(__STAMP__, &
            ' Wrong Space for SFdepoLayer: only cuboid and cylinder implemented!')
        END SELECT
        CALL RANDOM_NUMBER(RandVal)
        layerPartPos = layerPartPos + RandVal*SFdepoLayersGeo(iLayer,2,:)
        IF ( SFdepoLayersUseFixBounds(iLayer) ) THEN
          DoCycle=.FALSE.
          DO iSFfix=1,NbrOfSFdepoFixes
            SFfixDistance = SFdepoFixesGeo(iSFfix,2,1)*(layerPartPos(1)-SFdepoFixesGeo(iSFfix,1,1)) &
              + SFdepoFixesGeo(iSFfix,2,2)*(layerPartPos(2)-SFdepoFixesGeo(iSFfix,1,2)) &
              + SFdepoFixesGeo(iSFfix,2,3)*(layerPartPos(3)-SFdepoFixesGeo(iSFfix,1,3))
            IF (SFfixDistance .GT. 0.) THEN !outside of plane
              DoCycle=.TRUE.
              EXIT
            END IF
          END DO
          IF (DoCycle) CYCLE
        ELSE IF ( SFdepoLayersBounds(iLayer,1,1).GT.layerPartPos(1) .OR. layerPartPos(1).GT.SFdepoLayersBounds(iLayer,2,1) .OR. &
          SFdepoLayersBounds(iLayer,1,2).GT.layerPartPos(2) .OR. layerPartPos(2).GT.SFdepoLayersBounds(iLayer,2,2) .OR. &
          SFdepoLayersBounds(iLayer,1,3).GT.layerPartPos(3) .OR. layerPartPos(3).GT.SFdepoLayersBounds(iLayer,2,3) ) THEN
          CYCLE !outside of bounds
        END IF
        CALL calcSfSource(1 &
          ,Species(SFdepoLayersSpec(iLayer))%ChargeIC*Species(SFdepoLayersSpec(iLayer))%MacroParticleFactor*w_sf &
          ,Vec1,Vec2,Vec3,layerPartPos,iPart)
      END DO ! iPart
    END DO ! iLayer=1,NbrOfSFdepoLayers

    !--SFResampleAnalyzeSurfCollis 
    IF (SFResampleAnalyzeSurfCollis) THEN
      iPart=0
      DO iPart2=1,LastAnalyzeSurfCollis%PartNumberDepo
        !get random (equal!) position between [1,PartNumberSamp]
        CALL RANDOM_NUMBER(RandVal)
        iPart=MIN(1+INT(RandVal*REAL(LastAnalyzeSurfCollis%PartNumberSamp)),LastAnalyzeSurfCollis%PartNumberSamp)
        !perform surfaceflux-like push into sf-layer outside of mesh
        CALL RANDOM_NUMBER(RandVal)
        FractPush = RandVal*LastAnalyzeSurfCollis%pushTimeStep*LastAnalyzeSurfCollis%WallState(4:6,iPart)
        IF ( DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,FractPush).LE.r_SF  ) THEN
          layerPartPos = LastAnalyzeSurfCollis%WallState(1:3,iPart) + FractPush
          IF ( LastAnalyzeSurfCollis%UseFixBounds ) THEN
            DoCycle=.FALSE.
            DO iSFfix=1,NbrOfSFdepoFixes
              SFfixDistance = SFdepoFixesGeo(iSFfix,2,1)*(layerPartPos(1)-SFdepoFixesGeo(iSFfix,1,1)) &
                + SFdepoFixesGeo(iSFfix,2,2)*(layerPartPos(2)-SFdepoFixesGeo(iSFfix,1,2)) &
                + SFdepoFixesGeo(iSFfix,2,3)*(layerPartPos(3)-SFdepoFixesGeo(iSFfix,1,3))
              IF (SFfixDistance .GT. 0.) THEN !outside of plane
                DoCycle=.TRUE.
                EXIT
              END IF
            END DO
            IF (DoCycle) CYCLE
          ELSE IF ( LastAnalyzeSurfCollis%Bounds(1,1).GT.layerPartPos(1) .OR. &
            layerPartPos(1).GT.LastAnalyzeSurfCollis%Bounds(2,1) .OR. &
            LastAnalyzeSurfCollis%Bounds(1,2).GT.layerPartPos(2) .OR. &
            layerPartPos(2).GT.LastAnalyzeSurfCollis%Bounds(2,2) .OR. &
            LastAnalyzeSurfCollis%Bounds(1,3).GT.layerPartPos(3) .OR. &
            layerPartPos(3).GT.LastAnalyzeSurfCollis%Bounds(2,3) ) THEN
            CYCLE !outside of bounds
          END IF
        ELSE
          CYCLE !outside of r_SF
        END IF
        CALL calcSfSource(4 &
          ,Species(LastAnalyzeSurfCollis%Species(iPart))%ChargeIC &
          *Species(LastAnalyzeSurfCollis%Species(iPart))%MacroParticleFactor*w_sf &
          ,Vec1,Vec2,Vec3,layerPartPos,iPart2,PartVelo=LastAnalyzeSurfCollis%WallState(4:6,iPart))
      END DO ! iPart2
    END IF !SFResampleAnalyzeSurfCollis

    !-- external particles
#if USE_MPI     
      DO iPart=1,NbrOfextParticles  !external Particles
        CALL calcSfSource(4 &
          ,Species(ExtPartSpecies(iPart))%ChargeIC*Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf &
          ,Vec1,Vec2,Vec3,ExtPartState(iPart,1:3),iPart,PartVelo=ExtPartState(iPart,4:6))
      END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    NbrOfExtParticles=0
#endif /*MPI*/
  END IF !.NOT.DoInnerParts

  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map PartSource from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('shape_function_1d')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  DO iPart=firstPart,LastPart
    !IF (PDM%ParticleInside(iPart)) THEN
    IF (DoParticle(iPart)) THEN
      Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
      !IF(fac(4).GT.0.) print*,'charge pos'
      Fac(1:3) = PartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      chargedone(:) = .FALSE.
      DO iCase = 1, NbrOfCases
        DO ind = 1,3
          ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        IF(sf1d_dir.EQ.1)THEN
          kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
          kmax = MIN(kmax,GEO%FIBGMimax)
          kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
          kmin = MAX(kmin,GEO%FIBGMimin)
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSEIF(sf1d_dir.EQ.2)THEN
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
          lmax = MIN(lmax,GEO%FIBGMjmax)
          lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
          lmin = MAX(lmin,GEO%FIBGMjmin)
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSE
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
          mmax = MIN(mmax,GEO%FIBGMkmax)
          mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
          mmin = MAX(mmin,GEO%FIBGMkmin)
        END IF
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#if USE_MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                    radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                            * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. r2_sf) THEN
                      S = 1. - r2_sf_inv * radius2
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO ! iCase (periodicity)
    END IF ! inside
  END DO ! i
#if USE_MPI           
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF
    
    DO iPart=1,NbrOfextParticles  !external Particles
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
      Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
      chargedone(:) = .FALSE.
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        DO ind = 1,3
          ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        IF(sf1d_dir.EQ.1)THEN
          kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
          kmax = MIN(kmax,GEO%FIBGMimax)
          kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
          kmin = MAX(kmin,GEO%FIBGMimin)
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSEIF(sf1d_dir.EQ.2)THEN
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
          lmax = MIN(lmax,GEO%FIBGMjmax)
          lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
          lmin = MAX(lmin,GEO%FIBGMjmin)
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSE
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
          mmax = MIN(mmax,GEO%FIBGMkmax)
          mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
          mmin = MAX(mmin,GEO%FIBGMkmin)
        END IF
        !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain)
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#if USE_MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                      radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                              * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                      !-- calculate charge and current density at ip point using a shape function
                      !-- currently only one shapefunction available, more to follow (including structure change)
                      IF (radius2 .LT. r2_sf) THEN
                        S = 1. - r2_sf_inv * radius2
                        !IF(S.LT.0.) print*,'dist neg '
                      !radius2=GaussDistance(k,l,m)
                      !IF (radius2 .LT. 1.0) THEN
                      !  S = 1 -  radius2
                        S1 = S*S
                        DO expo = 3, alpha_sf
                          S1 = S*S1
                        END DO
                        PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                        PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                      END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO
    END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    NbrOfExtParticles=0
  END IF
#endif /*MPI*/

  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map source from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('shape_function_cylindrical','shape_function_spherical')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  DO iPart=firstPart,LastPart
    !IF (PDM%ParticleInside(iPart)) THEN
    IF (DoParticle(iPart)) THEN
      ! compute local radius
      local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(iPart,1:SfRadiusInt),PartState(iPart,1:SfRadiusInt)))
      local_r2_sf=local_r_sf*local_r_sf
      local_r2_sf_inv=1./local_r2_sf
      Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
      !IF(fac(4).GT.0.) print*,'charge pos'
      Fac(1:3) = PartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        chargedone(:) = .FALSE.
        DO ind = 1,3
          ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#if USE_MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                    dX = ABS(ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID))
                    IF(dX.GT.local_r_sf) CYCLE
                    dY = ABS(ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID))
                    IF(dY.GT.local_r_sf) CYCLE
                    dZ = ABS(ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                    IF(dZ.GT.local_r_sf) CYCLE
                    radius2 = dX*dX+dY*dY+dZ*dZ
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. local_r2_sf) THEN
                      S = 1. - local_r2_sf_inv * radius2
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO ! iCase (periodicity)
    END IF ! inside
  END DO ! i
#if USE_MPI           
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF
    
    DO iPart=1,NbrOfextParticles  !external Particles
      local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(iPart,1:SfRadiusInt),PartState(iPart,1:SfRadiusInt)))
      local_r2_sf=local_r_sf*local_r_sf
      local_r2_sf_inv=1./local_r2_sf
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC & 
              * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
      Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        chargedone(:) = .FALSE.
        DO ind = 1,3
          ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#if USE_MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                      dX = ABS(ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID))
                      IF(dX.GT.local_r_sf) CYCLE
                      dY = ABS(ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID))
                      IF(dY.GT.local_r_sf) CYCLE
                      dZ = ABS(ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                      IF(dZ.GT.local_r_sf) CYCLE
                      radius2 = dX*dX+dY*dY+dZ*dZ
                      !-- calculate charge and current density at ip point using a shape function
                      !-- currently only one shapefunction available, more to follow (including structure change)
                      IF (radius2 .LT. local_r2_sf) THEN
                        S = 1. - local_r2_sf_inv * radius2
                        !IF(S.LT.0.) print*,'dist neg '
                      !radius2=GaussDistance(k,l,m)
                      !IF (radius2 .LT. 1.0) THEN
                      !  S = 1 -  radius2
                        S1 = S*S
                        DO expo = 3, alpha_sf
                          S1 = S*S1
                        END DO
                        PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                        PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                      END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO
    END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    NbrOfExtParticles=0
  END IF
#endif /*MPI*/
  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map PartSource from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('delta_distri')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  DO iElem=1,PP_nElems
!#if USE_MPI
!    tLBStart = FLEXITIME() ! LB Time Start
!#endif /*MPI*/
    DO iPart=firstPart,LastPart
      ! IF (PDM%ParticleInside(iPart)) THEN
      IF (DoParticle(iPart)) THEN
        IF(PEM%Element(iPart).EQ.iElem)THEN
          prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
          ! Map Particle to -1|1 space (re-used in interpolation)
          IF(.NOT.DoRefMapping)THEN
            CALL TensorProductInterpolation(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem)
          END IF
          ! get value of test function at particle position
          SELECT CASE(DeltaType)
          CASE(1)
            ! xi   -direction
            CALL LagrangeInterpolationPolys(PartPosRef(1,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(1,:))
            ! eta  -direction                                 
            CALL LagrangeInterpolationPolys(PartPosRef(2,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(2,:))
            ! zeta -direction                                 
            CALL LagrangeInterpolationPolys(PartPosRef(3,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(3,:))
          CASE(2)
            DO i=0,NDepo
              ! xi   -direction
               CALL BernsteinPolynomial(NDepo,i,PartPosRef(1,iPart),L_xi(1,i),NDepoChooseK)
              ! eta  -direction
               CALL BernsteinPolynomial(NDepo,i,PartPosRef(2,iPart),L_xi(2,i),NDepoChooseK)
              ! zeta  -direction
               CALL BernsteinPolynomial(NDepo,i,PartPosRef(3,iPart),L_xi(3,i),NDepoChooseK)
            END DO ! i
          CASE(3)
            ! xi - direction
            CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(1,iPart),L_xi(1,:))
            ! eta - direction                                  
            CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(2,iPart),L_xi(2,:))
            ! zeta - direction                                 
            CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(3,iPart),L_xi(3,:))
          END SELECT
          DO k=0,NDepo
            DO j=0,NDepo
              DO i=0,NDepo
           !     print*,'i,j,k,L',i,j,k,L_xi(1,i)* L_xi(2,j)* L_xi(3,k)
                DeltaIntCoeff = L_xi(1,i)* L_xi(2,j)* L_xi(3,k)*prefac 
                PartSource(1:3,i,j,k,iElem) = PartSource(1:3,i,j,k,iElem) + DeltaIntCoeff*PartState(iPart,4:6)
                PartSource( 4 ,i,j,k,iElem) = PartSource( 4 ,i,j,k,iElem) + DeltaIntCoeff
              END DO ! i
            END DO ! j
          END DO ! k
        END IF ! Particle in Element
      END IF ! ParticleInside of domain
    END DO ! ParticleVecLength
  END DO ! iElem
  IF(.NOT.DoInnerParts)THEN
    DO iElem=1,PP_nElems
      DO k=0,NDepo
        DO j=0,NDepo
          DO i=0,NDepo
            PartSource( : ,i,j,k,iElem) = PartSource( : ,i,j,k,iElem) * DDMassInv(i,j,k,iElem)
          END DO ! i
        END DO ! j
      END DO ! k
      IF(DeltaDistriChangeBasis)THEN
        CALL ChangeBasis3D(4,NDepo,PP_N,Vdm_NDepo_GaussN,PartSource(1:4,0:NDepo,0:NDepo,0:NDepo,iElem)&
                                                        ,PartSource(1:4,0:PP_N ,0:PP_N ,0:PP_N, iElem))
      END IF
    END DO ! loop over all elems
  END IF ! DoInnerParts
  
CASE('nearest_gausspoint')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  SAVE_GAUSS = .FALSE.
  IF(TRIM(InterpolationType).EQ.'nearest_gausspoint') SAVE_GAUSS = .TRUE.
  IF(MOD(PP_N,2).EQ.0) THEN
    a = PP_N/2
    b = a
  ELSE
    a = (PP_N+1)/2
    b = a-1
  END IF
  DO iElem=1,PP_nElems
    DO iPart=firstPart,LastPart
     ! IF (PDM%ParticleInside(iPart)) THEN
      IF (DoParticle(iPart)) THEN
        IF(PEM%Element(iPart).EQ.iElem)THEN
          prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
          ! Map Particle to -1|1 space (re-used in interpolation)
!          IF(.NOT.DoRefMapping)THEN
!            CALL TensorProductInterpolation(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem,iPart)
!          END IF
          ! Find out which gausspoint is closest and add up charges and currents
          !! x-direction
          IF(.NOT.SAVE_GAUSS) THEN
            k = a
            DO ii = 0,b-1
              IF(ABS(PartPosRef(1,iPart)).GE.GaussBorder(PP_N-ii))THEN
                k = PP_N-ii
                EXIT
              END IF
            END DO
            k = NINT((PP_N+SIGN(2.0*k-PP_N,PartPosRef(1,iPart)))/2)
            !! y-direction
            l = a
            DO ii = 0,b-1
              IF(ABS(PartPosRef(2,iPart)).GE.GaussBorder(PP_N-ii))THEN
                l = PP_N-ii
                EXIT
              END IF
            END DO
            l = NINT((PP_N+SIGN(2.0*l-PP_N,PartPosRef(2,iPart)))/2)
            !! z-direction
            m = a
            DO ii = 0,b-1
              IF(ABS(PartPosRef(3,iPart)).GE.GaussBorder(PP_N-ii))THEN
                m = PP_N-ii
                EXIT
              END IF
            END DO
            m = NINT((PP_N+SIGN(2.0*m-PP_N,PartPosRef(3,iPart)))/2)
          END IF
          PartSource(1:3,k,l,m,iElem) = PartSource(1:3,k,l,m,iElem) + PartState(iPart,4:6) * prefac
          PartSource( 4 ,k,l,m,iElem) = PartSource( 4 ,k,l,m,iElem) + prefac
        END IF ! Element .EQ. iElem
      END IF ! Particle inside
    END DO ! iPart
  END DO ! iElem=1,PP_nElems
  IF(.NOT.DoInnerParts)THEN
    DO iElem=1,PP_nElems
      DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
      ! get densities by dividing by gauss volume
        PartSource(1:4,k,l,m,iElem) = PartSource(1:4,k,l,m,iElem) * sJ(k,l,m,iElem,0)/(wGP(k)*wGP(l)*wGP(m))
      END DO; END DO; END DO
    END DO ! iElem=1,PP_nElems
  END IF
  
CASE('cartmesh_volumeweighting')
  ! Step 1: Deposition of all particles onto background mesh -> densities
  ! IF(DoInnerParts) BGMSource=0.0 ! not possible due to periodic stuff --> two communications
  BGMSource(:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
    !IF (PDM%ParticleInside(iPart)) THEN
    IF (DoParticle(iPart)) THEN
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
      k = FLOOR(PartState(iPart,1)/BGMdeltas(1))
      l = FLOOR(PartState(iPart,2)/BGMdeltas(2))
      m = FLOOR(PartState(iPart,3)/BGMdeltas(3))
      alpha1 = (PartState(iPart,1) / BGMdeltas(1)) - k
      alpha2 = (PartState(iPart,2) / BGMdeltas(2)) - l
      alpha3 = (PartState(iPart,3) / BGMdeltas(3)) - m
      TSource(:) = 0.0
      TSource(1) = PartState(iPart,4)*Charge
      TSource(2) = PartState(iPart,5)*Charge
      TSource(3) = PartState(iPart,6)*Charge
      TSource(4) = Charge

      BGMSource(k,l,m,1:4)       = BGMSource(k,l,m,1:4) + (TSource * (1-alpha1)*(1-alpha2)*(1-alpha3))
      BGMSource(k,l,m+1,1:4)     = BGMSource(k,l,m+1,1:4) + (TSource * (1-alpha1)*(1-alpha2)*(alpha3))
      BGMSource(k,l+1,m,1:4)     = BGMSource(k,l+1,m,1:4) + (TSource * (1-alpha1)*(alpha2)*(1-alpha3))
      BGMSource(k,l+1,m+1,1:4)   = BGMSource(k,l+1,m+1,1:4) + (TSource * (1-alpha1)*(alpha2)*(alpha3))
      BGMSource(k+1,l,m,1:4)     = BGMSource(k+1,l,m,1:4) + (TSource * (alpha1)*(1-alpha2)*(1-alpha3))
      BGMSource(k+1,l,m+1,1:4)   = BGMSource(k+1,l,m+1,1:4) + (TSource * (alpha1)*(1-alpha2)*(alpha3))
      BGMSource(k+1,l+1,m,1:4)   = BGMSource(k+1,l+1,m,1:4) + (TSource * (alpha1)*(alpha2)*(1-alpha3))
      BGMSource(k+1,l+1,m+1,1:4) = BGMSource(k+1,l+1,m+1,1:4) + (TSource * (alpha1)*(alpha2)*(alpha3))
    END IF
  END DO
  BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#if USE_MPI
  ! should be treated in this way, unfortunately, we would neglect the periodic stuff
  !IF(.NOT.DoInnerParts)
  CALL MPISourceExchangeBGM()
#else
  IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif

  ! Step 2: Interpolation of densities onto grid
  DO iElem = 1, nElems
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
         k = GaussBGMIndex(1,kk,ll,mm,iElem)
         l = GaussBGMIndex(2,kk,ll,mm,iElem)
         m = GaussBGMIndex(3,kk,ll,mm,iElem)
         alpha1 = GaussBGMFactor(1,kk,ll,mm,iElem)
         alpha2 = GaussBGMFactor(2,kk,ll,mm,iElem)
         alpha3 = GaussBGMFactor(3,kk,ll,mm,iElem)
         DO i = 1,3
           PartSource(i,kk,ll,mm,iElem) = PartSource(i,kk,ll,mm,iElem)            + &
                BGMSource(k,l,m,i) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k,l,m+1,i) * (1-alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k,l+1,m,i) * (1-alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k,l+1,m+1,i) * (1-alpha1) * (alpha2) * (alpha3) + &
                BGMSource(k+1,l,m,i) * (alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k+1,l,m+1,i) * (alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k+1,l+1,m,i) * (alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k+1,l+1,m+1,i) * (alpha1) * (alpha2) * (alpha3)
         END DO
           PartSource(4,kk,ll,mm,iElem) = PartSource(4,kk,ll,mm,iElem)          + &
              BGMSource(k,l,m,4) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k,l,m+1,4) * (1-alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k,l+1,m,4) * (1-alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k,l+1,m+1,4) * (1-alpha1) * (alpha2) * (alpha3) + &
              BGMSource(k+1,l,m,4) * (alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k+1,l,m+1,4) * (alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k+1,l+1,m,4) * (alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k+1,l+1,m+1,4) * (alpha1) * (alpha2) * (alpha3)
       END DO !mm
     END DO !ll
   END DO !kk
 END DO !iElem
 !DEALLOCATE(BGMSource)
 
CASE('cartmesh_splines')
  ! Step 1: Deposition of all particles onto background mesh -> densities
  !ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
  ! IF(DoInnerParts) BGMSource=0. not possible due to periodic stuff
  BGMSource(:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
    !IF (PDM%ParticleInside(iPart)) THEN
    IF (DoParticle(iPart)) THEN
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
      PosInd(1) = FLOOR(PartState(iPart,1)/BGMdeltas(1))
      PosInd(2) = FLOOR(PartState(iPart,2)/BGMdeltas(2))
      PosInd(3) = FLOOR(PartState(iPart,3)/BGMdeltas(3))
      !print*,'posind(1:3),charge',posInd,charge
      DO dir = 1,3               ! x,y,z direction
        DO weightrun = 0,3
          DO mm = 0, 3
            IF (mm.EQ.weightrun) then
              auxiliary(mm) = 1.0
            ELSE
              auxiliary(mm) = 0.0
            END IF
          END DO
          CALL DeBoor(PosInd(dir),auxiliary,PartState(iPart,dir),weight(dir,weightrun),dir)
        END DO
      END DO
      DO k = PosInd(1)-1, PosInd(1)+2
        kk = abs(k - PosInd(1) - 2)
        DO l = PosInd(2)-1, PosInd(2)+2
          ll = abs(l - PosInd(2) - 2)
          DO m = PosInd(3)-1, PosInd(3)+2
            mm = abs(m - PosInd(3) - 2)
            locweight = weight(1,kk)*weight(2,ll)*weight(3,mm)*charge
            BGMSource(k,l,m,1) = BGMSource(k,l,m,1) + PartState(iPart,4)* locweight
            BGMSource(k,l,m,2) = BGMSource(k,l,m,2) + PartState(iPart,5)* locweight
            BGMSource(k,l,m,3) = BGMSource(k,l,m,3) + PartState(iPart,6)* locweight
            BGMSource(k,l,m,4) = BGMSource(k,l,m,4) + locweight
         !   print*,'BMGSOURCE4',BGMSOURCE(k,l,m,4)
          END DO
        END DO
      END DO
    END IF
  END DO
  BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#if USE_MPI
  !IF(.NOT.DoInnerParts)THEN has to be communicated each time :(
  CALL MPISourceExchangeBGM()
#else
  IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif
  ! Step 2: Interpolation of densities onto grid
  DO iElem = 1, nElems
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
          k = GaussBGMIndex(1,kk,ll,mm,iElem)
          l = GaussBGMIndex(2,kk,ll,mm,iElem)
          m = GaussBGMIndex(3,kk,ll,mm,iElem)
          DO r = k-1,k+2
            u = r-k+2
            DO ss = l-1,l+2
              v = ss-l+2
              DO t = m-1,m+2
                w = t-m+2     
                PartSource(1:4,kk,ll,mm,iElem) = PartSource(1:4,kk,ll,mm,iElem)+BGMSource(r,ss,t,1:4)*GPWeight(iElem,kk,ll,mm,u,v,w)
              END DO !t
            END DO !s
          END DO !r
        END DO !mm
      END DO !ll
    END DO !kk
  END DO !iElem
 !DEALLOCATE(BGMSource)
  
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Unknown DepositionType in pic_depo.f90')
END SELECT

RETURN
END SUBROUTINE Deposition


#if USE_MPI
SUBROUTINE MPISourceExchangeBGM()            
!=================================================================================================================================
! Exchange sources in periodic case for MPI
!==================================================================================================================================
! use MODULES                                                         
USE MOD_Particle_MPI_Vars,  ONLY: PartMPI,tMPIMEssage
USE MOD_Particle_Mesh_Vars, ONLY: GEO
USE MOD_PICDepo_Vars
USE MOD_Globals
USE MOD_Particle_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tMPIMessage)           :: send_message(0:PartMPI%nProcs-1)                              
TYPE(tMPIMessage)           :: recv_message(0:PartMPI%nProcs-1)                              
INTEGER                     :: send_request(0:PartMPI%nProcs-1)                              
INTEGER                     :: recv_request(0:PartMPI%nProcs-1)                              
INTEGER                     :: send_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
INTEGER                     :: recv_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
INTEGER                     :: iProc,k,l,m,n, ppp, Counter, MsgLength(0:PartMPI%nProcs-1), iPer
INTEGER                     :: SourceLength(0:PartMPI%nProcs-1)                              
INTEGER                     :: RecvLength(0:PartMPI%nProcs-1)                                
INTEGER                     :: allocStat, Counter2                                     
INTEGER                     :: messageCounterS, messageCounterR                                
INTEGER                     :: myRealKind, k2,l2,m2                                          
REAL                        :: myRealTestValue                                                
!-----------------------------------------------------------------------------------------------------------------------------------

myRealKind = KIND(myRealTestValue)
IF (myRealKind.EQ.4) THEN
 myRealKind = MPI_REAL
ELSE IF (myRealKind.EQ.8) THEN
 myRealKind = MPI_DOUBLE_PRECISION
ELSE
 myRealKind = MPI_REAL
END IF
     
!--- Determine which Sources actually need to be sent (<> 0) and build corresponding true/false list
!    One list per process
DO iProc = 0,PartMPI%nProcs-1
   MsgLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      MsgLength(iProc) = (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3) + 1)
   END IF
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO k = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         MsgLength(iProc) = MsgLength(iProc) + &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,1) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,1) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,2) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,2) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,3) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,3) + 1)
      END DO
   END IF
   IF((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor))THEN
      ALLOCATE(send_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
      ALLOCATE(recv_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
      END IF
   END IF
   !--- check which sources are <> 0
   Counter = 0
   SourceLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
               send_message(iProc)%content_log(Counter) = .TRUE.
               SourceLength(iProc) = SourceLength(iProc) + 1
            ELSE
               send_message(iProc)%content_log(Counter) = .FALSE.
            END IF
          END DO
        END DO
      END DO
   END IF
   !--- same for periodic
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
          DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                 PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
            DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                   PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
               Counter = Counter + 1
               IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
                  send_message(iProc)%content_log(Counter) = .TRUE.
                  SourceLength(iProc) = SourceLength(iProc) + 1
               ELSE
                  send_message(iProc)%content_log(Counter) = .FALSE.
               END IF
            END DO
          END DO
         END DO
      END DO
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%Content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF

!--- Assemble actual sources to send
 DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ALLOCATE(send_message(iProc)%content(1:SourceLength(iProc)*4), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
   END IF
   Counter = 0
   Counter2 = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
             Counter = Counter + 1
             IF (send_message(iProc)%content_log(Counter)) THEN
                Counter2 = Counter2 + 1
                DO n = 1,4
                   send_message(iProc)%content((Counter2-1)*4 +n) = BGMSource(k,l,m,n)
                END DO
             END IF
          END DO
        END DO
      END DO
   END IF
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
           DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
             DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                Counter = Counter + 1
                IF (send_message(iProc)%content_log(Counter)) THEN
                   Counter2 = Counter2 + 1
                   DO ppp = 1,4
                      send_message(iProc)%content((Counter2-1)*4 +ppp) = BGMSource(k,l,m,ppp)
                   END DO
                END IF
             END DO
           END DO
         END DO
      END DO
   END IF
END DO

!--- allocate actual PartSource receive buffer
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      Counter = 0
      DO k = 1, MsgLength(iProc)
         IF (recv_message(iProc)%Content_log(k)) THEN
            Counter = Counter + 1
         END IF
      END DO
      RecvLength(iProc) = Counter
      IF (RecvLength(iProc).GT.0) THEN
         ALLOCATE(recv_message(iProc)%content(1:Counter*4), STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
         END IF
      END IF
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content,SourceLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(RecvLength(iProc).GT.0)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%content,RecvLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF
!--- Deallocate Send Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      DEALLOCATE(send_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
      END IF
      IF (SourceLength(iProc).GT.0) THEN
         DEALLOCATE(send_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
         END IF
      END IF
   END IF
END DO

!--- add selfperiodic sources, if any (needs to be done after send message is compiled and before
!---           received sources have been added!
IF ((GEO%nPeriodicVectors.GT.0).AND.(GEO%SelfPeriodic)) THEN
   DO iPer = 1, GEO%nPeriodicVectors
      DO k = BGMminX, BGMmaxX
         k2 = k + GEO%PeriodicBGMVectors(1,iPer)
        DO l = BGMminY, BGMmaxY
          l2 = l + GEO%PeriodicBGMVectors(2,iPer)
          DO m = BGMminZ, BGMmaxZ
             m2 = m + GEO%PeriodicBGMVectors(3,iPer)
             IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
             IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
             IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
                BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
                BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
             END IF
             END IF
             END IF
          END DO
        END DO
      END DO
   END DO
END IF

!--- Add Sources and Deallocate Receive Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF (RecvLength(iProc).GT.0) THEN
      Counter = 0
      Counter2 = 0
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
         DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
         DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
         DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF(recv_message(iProc)%content_log(Counter))THEN
               Counter2 = Counter2 + 1
               DO n = 1,4
                 BGMSource(k,l,m,n) = BGMSource(k,l,m,n) + recv_message(iProc)%content((Counter2-1)*4+n)
               END DO
            END IF
         END DO
         END DO
         END DO
      END IF
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
         DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
           DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
             DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
               DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                      PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                  Counter = Counter + 1
                  IF(recv_message(iProc)%content_log(Counter))THEN
                     Counter2 = Counter2 + 1
                     DO ppp = 1,4
                        BGMSource(k,l,m,ppp) = BGMSource(k,l,m,ppp) + recv_message(iProc)%content((Counter2-1)*4+ppp)
                     END DO
                  END IF
               END DO
             END DO
           END DO
         END DO
      END IF
      IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
         DEALLOCATE(recv_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
         END IF
      END IF
   END IF
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
      DEALLOCATE(recv_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
      END IF
   END IF
END DO

END SUBROUTINE MPISourceExchangeBGM


SUBROUTINE MPIBackgroundMeshInit()  
!==================================================================================================================================
! initialize MPI background mesh
!==================================================================================================================================
! use MODULES                                                                   
USE MOD_PICDepo_Vars
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Mesh_Vars,  ONLY:GEO
USE MOD_Particle_MPI_Vars,   ONLY:PartMPI
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iProc
INTEGER                     :: k,m,n,m0,n0
INTEGER                     :: localminmax(6), maxofmin, minofmax                              
INTEGER                     :: completeminmax(6*PartMPI%nProcs)                              
INTEGER                     :: allocStat, NeighCount                                                                 
INTEGER                     :: TempBorder(1:2,1:3)            
INTEGER                     :: Periodicminmax(6), coord, PeriodicVec(1:3)                   
INTEGER                     :: TempPeriBord(1:26,1:2,1:3)                                      
LOGICAL                     :: CHECKNEIGHBOR     
!-----------------------------------------------------------------------------------------------------------------------------------

!Periodic Init stuff
IF(GEO%nPeriodicVectors.GT.0)THEN
  ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
  ALLOCATE(GEO%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
  IF (allocStat .NE. 0) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR in MPIBackgroundMeshInit: cannot allocate GEO%PeriodicBGMVectors!')
  END IF
  DO iProc = 1, GEO%nPeriodicVectors
    GEO%PeriodicBGMVectors(1,iProc) = NINT(GEO%PeriodicVectors(1,iProc)/BGMdeltas(1))
    IF(ABS(GEO%PeriodicVectors(1,iProc)/BGMdeltas(1)-REAL(GEO%PeriodicBGMVectors(1,iProc))).GT.1E-10)THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
    GEO%PeriodicBGMVectors(2,iProc) = NINT(GEO%PeriodicVectors(2,iProc)/BGMdeltas(2))
    IF(ABS(GEO%PeriodicVectors(2,iProc)/BGMdeltas(2)-REAL(GEO%PeriodicBGMVectors(2,iProc))).GT.1E-10)THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
    GEO%PeriodicBGMVectors(3,iProc) = NINT(GEO%PeriodicVectors(3,iProc)/BGMdeltas(3))
    IF(ABS(GEO%PeriodicVectors(3,iProc)/BGMdeltas(3)-REAL(GEO%PeriodicBGMVectors(3,iProc))).GT.1E-10)THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
  END DO
  ! Check whether process is periodic with itself
  GEO%SelfPeriodic = .FALSE.
  !--- virtually move myself according to periodic vectors in order to find overlapping areas
  !--- 26 possibilities,
  localminmax(1) = BGMminX
  localminmax(2) = BGMminY
  localminmax(3) = BGMminZ
  localminmax(4) = BGMmaxX
  localminmax(5) = BGMmaxY
  localminmax(6) = BGMmaxZ
  DO k = -1,1
    DO m = -1,1
      DO n = -1,1
        PeriodicVec = k*GEO%PeriodicBGMVectors(:,1)
        IF (GEO%nPeriodicVectors.GT.1) THEN
          PeriodicVec = PeriodicVec + m*GEO%PeriodicBGMVectors(:,2)
        END IF
        IF (GEO%nPeriodicVectors.GT.2) THEN
          PeriodicVec = PeriodicVec + n*GEO%PeriodicBGMVectors(:,3)
        END IF
        IF (ALL(PeriodicVec(:).EQ.0)) CYCLE
        periodicminmax(1) = localminmax(1) + PeriodicVec(1)
        periodicminmax(2) = localminmax(2) + PeriodicVec(2)
        periodicminmax(3) = localminmax(3) + PeriodicVec(3)
        periodicminmax(4) = localminmax(4) + PeriodicVec(1)
        periodicminmax(5) = localminmax(5) + PeriodicVec(2)
        periodicminmax(6) = localminmax(6) + PeriodicVec(3)
        !--- find overlap
        DO coord = 1,3  !          x y z direction
          maxofmin = MAX(periodicminmax(coord),localminmax(coord))
          minofmax = MIN(periodicminmax(3+coord),localminmax(3+coord))
          IF (maxofmin.LE.minofmax) GEO%SelfPeriodic = .TRUE.     !  overlapping
        END DO
      END DO
    END DO
  END DO
END IF

! --- send and receive min max indices to and from all processes

!--- enter local min max vector (xmin, ymin, zmin, xmax, ymax, zmax)
localminmax(1) = BGMminX
localminmax(2) = BGMminY
localminmax(3) = BGMminZ
localminmax(4) = BGMmaxX
localminmax(5) = BGMmaxY
localminmax(6) = BGMmaxZ
!--- do allgather into complete min max vector
CALL MPI_ALLGATHER(localminmax,6,MPI_INTEGER,completeminmax,6,MPI_INTEGER,PartMPI%COMM,IERROR)
! Allocate MPIConnect
SDEALLOCATE(PartMPI%DepoBGMConnect)
ALLOCATE(PartMPI%DepoBGMConnect(0:PartMPI%nProcs-1),STAT=allocStat)
IF (allocStat .NE. 0) THEN
  CALL abort(&
  __STAMP__&
  ,' Cannot allocate PartMPI%DepoBGMConnect')
END IF

!--- determine borders indices (=overlapping BGM mesh points) with each process
DO iProc=0,PartMPI%nProcs-1
  PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
  PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = 0
   IF (iProc.EQ.PartMPI%MyRank) THEN
      PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .FALSE.
   ELSE
      PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .TRUE.
      DO k = 1,3          !  x y z direction
         maxofmin = MAX(localminmax(k),completeminmax((iProc*6)+k))
         minofmax = MIN(localminmax(3+k),completeminmax((iProc*6)+3+k))
         IF (maxofmin.LE.minofmax) THEN          !  overlapping
            TempBorder(1,k) = maxofmin
            TempBorder(2,k) = minofmax
         ELSE
            PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .FALSE.
         END IF
      END DO          
   END IF
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)THEN
      SDEALLOCATE(PartMPI%DepoBGMConnect(iProc)%BGMBorder)
      ALLOCATE(PartMPI%DepoBGMConnect(iProc)%BGMBorder(1:2,1:3),STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,' Cannot allocate PartMPI%DepoMPIConnect%BGMBorder')
      END IF
      PartMPI%DepoBGMConnect(iProc)%BGMBorder(1:2,1:3) = TempBorder(1:2,1:3)
   END IF
END DO ! iProc=0,PartMPI%nProcs-1

!--- determine border indices for periodic meshes  
IF (GEO%nPeriodicVectors.GT.0) THEN   
  !--- m-/-n-loops must be executed just once (with 0) if respective PV is not present
  !    (otherwise the unnec. 0 are sent and the 0-0-0-case occurs two 2 more which are not cycled, see below for more info)!
  IF (GEO%nPeriodicVectors.GT.1) THEN
    m0=1
  ELSE
    m0=0
  END IF
  IF (GEO%nPeriodicVectors.GT.2) THEN
    n0=1
  ELSE
    n0=0
  END IF
  DO iProc = 0,PartMPI%nProcs-1
    IF (iProc.EQ.PartMPI%MyRank) THEN
      PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
      PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = 0
    ELSE
      !--- virtually move myself according to periodic vectors in order to find overlapping areas
      !--- 26 possibilities, processes need to work through them in opposite direction in order
      !--- to get matching areas.
      !--- Example for 2D:  I am process #3, I compare myself with #7
      !--- Periodic Vectors are p1 and p2.
      !--- I check p1, p2, p1+p2, p1-p2, -p1+p2, -p1-p2, -p2, -p1
      !--- #7 has to check -p1, -p2, -p1-p2, -p1+p2, p1-p2, p1+p1, p2, p1
      !--- This is done by doing 3 loops from -1 to 1 (for the higher process number)
      !--- or 1 to -1 (for the lower process number) and multiplying
      !--- these numbers to the periodic vectors
      NeighCount = 0  ! -- counter: how often is the process my periodic neighbor?
      PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
      DO k = -SIGN(1,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
        DO m = -SIGN(m0,PartMPI%MyRank-iProc),SIGN(m0,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
          DO n = -SIGN(n0,PartMPI%MyRank-iProc),SIGN(n0,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
            IF ((k.EQ.0).AND.(m.EQ.0).AND.(n.EQ.0)) CYCLE ! this is not periodic and already done above
            CHECKNEIGHBOR = .TRUE.
            PeriodicVec = k*GEO%PeriodicBGMVectors(:,1)
            IF (GEO%nPeriodicVectors.GT.1) THEN
              PeriodicVec = PeriodicVec + m*GEO%PeriodicBGMVectors(:,2)
            END IF
            IF (GEO%nPeriodicVectors.GT.2) THEN
              PeriodicVec = PeriodicVec + n*GEO%PeriodicBGMVectors(:,3)
            END IF
            periodicminmax(1) = localminmax(1) + PeriodicVec(1)
            periodicminmax(2) = localminmax(2) + PeriodicVec(2)
            periodicminmax(3) = localminmax(3) + PeriodicVec(3)
            periodicminmax(4) = localminmax(4) + PeriodicVec(1)
            periodicminmax(5) = localminmax(5) + PeriodicVec(2)
            periodicminmax(6) = localminmax(6) + PeriodicVec(3)
            !--- find overlap
            DO coord = 1,3           ! x y z direction
              maxofmin = MAX(periodicminmax(coord),completeminmax((iProc*6)+coord))
              minofmax = MIN(periodicminmax(3+coord),completeminmax((iProc*6)+3+coord))
              IF (maxofmin.LE.minofmax) THEN         !   overlapping
                TempBorder(1,coord) = maxofmin
                TempBorder(2,coord) = minofmax
              ELSE
                CHECKNEIGHBOR = .FALSE.
              END IF
            END DO
            IF(CHECKNEIGHBOR)THEN
              NeighCount = NeighCount + 1
              TempBorder(:,1) = TempBorder(:,1) - PeriodicVec(1)
              TempBorder(:,2) = TempBorder(:,2) - PeriodicVec(2)
              TempBorder(:,3) = TempBorder(:,3) - PeriodicVec(3)
              TempPeriBord(NeighCount,1:2,1:3) = TempBorder(1:2,1:3)
              PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .TRUE.
            END IF
          END DO
        END DO
      END DO
      PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = NeighCount
      ALLOCATE(PartMPI%DepoBGMConnect(iProc)%Periodic(1:PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount),STAT=allocStat)
      IF (allocStat .NE. 0) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR in MPIBackgroundMeshInit: cannot allocate PartMPI%DepoBGMConnect')
      END IF
      DO k = 1,NeighCount
        ALLOCATE(PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1:2,1:3),STAT=allocStat)
        IF (allocStat .NE. 0) THEN
          CALL abort(&
          __STAMP__&
          ,'ERROR in MPIBackgroundMeshInit: cannot allocate PartMPI%DepoBGMConnect')
        END IF
        PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1:2,1:3) = TempPeriBord(k,1:2,1:3)
      END DO
    END IF
  END DO
ELSE
  !--- initialize to FALSE for completely non-periodic cases
  DO iProc = 0,PartMPI%nProcs-1
    PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
  END DO
END IF

END SUBROUTINE MPIBackgroundMeshInit
#else
SUBROUTINE PeriodicSourceExchange()    
!============================================================================================================================
! Exchange sources in periodic case
!============================================================================================================================
! use MODULES                                                    
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars,  ONLY: GEO
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT)         :: BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                           
INTEGER                     :: i,k,l,m,k2,l2,m2
!-----------------------------------------------------------------------------------------------------------------------------------

DO i = 1,GEO%nPeriodicVectors
  DO k = BGMminX, BGMmaxX
    k2 = k + GEO%PeriodicBGMVectors(1,i)
    DO l = BGMminY, BGMmaxY
      l2 = l + GEO%PeriodicBGMVectors(2,i)
      DO m = BGMminZ, BGMmaxZ
        m2 = m + GEO%PeriodicBGMVectors(3,i)
        IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
          IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
            IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
              BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
              BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
END DO

END SUBROUTINE PeriodicSourceExchange
#endif /*MPI*/

SUBROUTINE DeBoor(PosInd, aux, coord, results, dir)       
!============================================================================================================================
! recursive function for evaluating a b-spline basis function
!============================================================================================================================
! use MODULES 
   USE MOD_PICDepo_Vars                                                          
!-----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
   INTEGER, INTENT(IN)                          :: PosInd, dir  
   REAL, INTENT(IN)                             :: coord                                              
   REAL, INTENT(INOUT)                          :: aux(0:3), results                     
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES  
   INTEGER                          :: i,k,jL,jR                                                  
   REAL                              :: hlp1,hlp2                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
    DO i = 0, 2
       DO k = 0, 2-i
          jL = PosInd - k
          jR = jL + 3 - i
          
          hlp1 = jR * BGMdeltas(dir) - coord
          hlp2 = coord - jL * BGMdeltas(dir)
         
          aux(k) = (hlp1 * aux(k+1) + hlp2 * aux(k)) / (hlp1+hlp2)
       ENDDO
    ENDDO
    results = aux(0)

    RETURN
END SUBROUTINE DeBoor


SUBROUTINE DeBoorRef(N_in,NKnots,Knots,Xi_in,UBspline)
!===================================================================================================================================
! DeBoor algorithms for uniform B-Splines
! rule of DeBoor algorithm
! N_i^0 (x) = 1 if x in [u_i, u_i+1); else 0     
! N_i^n (x) = (x-u_i) /(u_i+n-u_i) N_i^n-1(x) + (u_i+n+1 - x)/(u_i+n+1 - u_i+1) N_i+1^n-1(x)
! this algorithm evaluates the complete 1D basis fuction, because certain knots can be reused
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: N_in,Nknots
REAL,INTENT(IN)         :: Xi_in
REAL,INTENT(OUT)        :: UBspline(0:N_in)
REAL,INTENT(IN)         :: knots(0:Nknots) ! range of parameter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,n, last
REAL                    :: tmpArray(0:NKnots-1,0:NKnots-1)
REAL                    :: sDxiN
!===================================================================================================================================

 
! init, first layer (constant)
! zero order
n=0
tmpArray=0.
last=nknots-1
DO i=0,last
  IF((knots(i).LE.Xi_in).AND.(Xi_in.LT.knots(i+1)))THEN
    tmpArray(i,n)=1.0
  END IF ! select
END DO ! i

DO n=1,N_in
  last=last-1
  DO i=0,last
    ! standard
!    tmpArray(i,n) = (Xi_in-knots(i))/(knots(i+n)-knots(i))*tmpArray(i,n-1) &
!                  + (knots(i+n+1)-Xi_in)/(knots(i+n+1)-knots(i+1))*tmpArray(i+1,n-1)
    ! optimized
    sDxiN=0.5/REAL(n)
    tmpArray(i,n) = sDxiN*( (Xi_in-knots(i) )*tmparray(i,n-1)+(knots(i+n+1)-Xi_in)*tmpArray(i+1,n-1))
  END DO ! i
END DO ! n

! move back to correct range
UBSpline(0:N_in)=tmpArray(0:N_in,N_in)

END SUBROUTINE DeBoorRef


FUNCTION beta(z,w)  
!============================================================================================================================
! calculates the beta function
!============================================================================================================================
! use MODULES                                                                                               
!-----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
   REAL, INTENT(IN)                             :: w, z                                             
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   REAL                                          :: beta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                    
!----------------------------------------------------------------------------------------------------------------------------------
   beta = GAMMA(z)*GAMMA(w)/GAMMA(z+w)                                                                    
END FUNCTION beta 


#ifdef donotcompilethis
SUBROUTINE ComputeGaussDistance(N_In,scaleR,X_in,Elem_xGP,GaussDistance) 
!----------------------------------------------------------------------------------------------------------------------------------!
! compute all distance between X_in and given array 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES 
INTEGER,INTENT(IN)      :: N_In
REAL,INTENT(IN)         :: X_in(1:3)
REAL,INTENT(IN)         :: scaleR
REAL,INTENT(IN)         :: Elem_xGP(1:3,0:N_in,0:N_in,0:N_in)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)        :: GaussDistance(0:N_in,0:N_in,0:N_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k
!===================================================================================================================================

DO k=0,N_in
  DO j=0,N_in
    DO i=0,N_in
      GaussDistance(i,j,k) =((X_in(1)-Elem_xGP(1,i,j,k))*(X_in(1)-Elem_xGP(1,i,j,k)) &
                            +(X_in(2)-Elem_xGP(2,i,j,k))*(X_in(2)-Elem_xGP(2,i,j,k)) &
                            +(X_in(3)-Elem_xGP(3,i,j,k))*(X_in(3)-Elem_xGP(3,i,j,k)))*scaleR
    END DO ! i=0,N_in
  END DO !  j=0,N_in
END DO ! k=0,N_in

END SUBROUTINE ComputeGaussDistance
#endif


SUBROUTINE calcSfSource(SourceSize_in,ChargeMPF,Vec1,Vec2,Vec3,PartPos,PartIdx,PartVelo)
!============================================================================================================================
! deposit charges on DOFs via shapefunction including periodic displacements and mirroring with SFdepoFixes
!============================================================================================================================
! use MODULES                                                                                               
USE MOD_PICDepo_Vars,           ONLY:r_sf,DepositionType
USE MOD_PICDepo_Vars,           ONLY:NbrOfSFdepoFixes,SFdepoFixesGeo,SFdepoFixesBounds,SFdepoFixesChargeMult
USE MOD_PICDepo_Vars,           ONLY:SFdepoFixesPartOfLink,SFdepoFixesEps,NbrOfSFdepoFixLinks,SFdepoFixLinks
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Mesh_Vars,     ONLY:casematrix,NbrOfCases
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: SourceSize_in,PartIdx
REAL, INTENT(IN)                 :: ChargeMPF,PartPos(3),Vec1(3),Vec2(3),Vec3(3)
REAL, INTENT(IN), OPTIONAL       :: PartVelo(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                    
INTEGER                          :: SourceSize
REAL                             :: Fac(1:SourceSize_in), Fac2(1:SourceSize_in)
INTEGER                          :: iCase, ind
REAL                             :: ShiftedPart(1:3), caseShiftedPart(1:3)
INTEGER                          :: iSFfix, LinkLoopEnd(2), iSFfixLink, iTwin, iLinkRecursive, SFfixIdx, SFfixIdx2
LOGICAL                          :: DoCycle, DoNotDeposit
REAL                             :: SFfixDistance, SFfixDistance2
LOGICAL , ALLOCATABLE            :: SFdepoFixDone(:)
!----------------------------------------------------------------------------------------------------------------------------------
SourceSize=SourceSize_in
IF (SourceSize.EQ.1) THEN
  Fac2= ChargeMPF
ELSE IF (SourceSize.EQ.4) THEN
  Fac2(1:3) = PartVelo*ChargeMPF
  Fac2(4)= ChargeMPF
ELSE
  CALL abort(&
__STAMP__ &
,'SourceSize has to be either 1 or 4!',SourceSize)
END IF
IF (NbrOfSFdepoFixes.EQ.0) THEN
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      ShiftedPart(ind) = PartPos(ind) + casematrix(iCase,1)*Vec1(ind) + &
        casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    Fac = Fac2
    SELECT CASE(TRIM(DepositionType))
    CASE('shape_function')
        ! Do nothing in FLEXI
    CASE('shape_function_simple')
        ! Do nothing in FLEXI
    END SELECT
  END DO ! iCase (periodicity)
ELSE ! NbrOfSFdepoFixes.NE.0
  ALLOCATE(SFdepoFixDone(0:NbrOfSFdepoFixes))
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      caseShiftedPart(ind) = PartPos(ind) + casematrix(iCase,1)*Vec1(ind) + &
        casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    SFdepoFixDone=.FALSE.
    DO iSFfix=0,NbrOfSFdepoFixes
      IF (SFdepoFixesPartOfLink(iSFfix)) CYCLE !this SFfix will be already covered as part of a FixLink
      IF (iSFfix.EQ.0) THEN
        LinkLoopEnd(1)=NbrOfSFdepoFixLinks !non-mirrored position: consider FixLinks
      ELSE
        LinkLoopEnd(1)=0
      END IF
      DO iSFfixLink=0,LinkLoopEnd(1)
        DoCycle=.FALSE.
        !-- strategy for SFdepoFixLink:
        !-- 1.: create 2 identities by iTwin for iLinkRecursive=0 (non-mirror and mirror at SFdepoFixLinks(*,1))
        !-- 2.: mirror recursevely further for iLinkRecursive>0 (alternating at SFdepoFixLinks(*,2) and SFdepoFixLinks(*,1))
        DO iTwin=1,2
          IF (iSFfixLink.EQ.0 .AND. iTwin.EQ.2) EXIT !no SFdepoFixLink
          Fac = Fac2
          ShiftedPart=caseShiftedPart
          IF (iSFfixLink.EQ.0) THEN
            LinkLoopEnd(2)=0 !no SFdepoFixLink
          ELSE
            LinkLoopEnd(2)=ABS(SFdepoFixLinks(iSFfixLink,3))-1 !might be negative as flag for special case
          END IF
          DO iLinkRecursive=0,LinkLoopEnd(2)
            IF (iLinkRecursive.EQ.0) THEN !(first!)
              IF (iTwin.EQ.1) THEN
                SFfixIdx = iSFfix
              ELSE
                SFfixIdx =SFdepoFixLinks(iSFfixLink,1)
              END IF
            ELSE
              IF (MOD(iLinkRecursive,2).NE.0) THEN !uneven (second and later: 1, 3, ...)
                SFfixIdx =SFdepoFixLinks(iSFfixLink,2)
              ELSE !even (third and later: 2, 4, ...)
                SFfixIdx =SFdepoFixLinks(iSFfixLink,1)
              END IF
            END IF
            DoNotDeposit=.FALSE.
            IF (iLinkRecursive.EQ.0 .OR. (iLinkRecursive.EQ.1 .AND. iTwin.EQ.1)) THEN !single or one-time-mirrored identity
              IF (SFdepoFixDone(SFfixIdx)) DoNotDeposit=.TRUE. !do not deposit this charge (but save position for further mirroring)
              SFdepoFixDone(SFfixIdx) = .TRUE.
            ELSE IF (iSFfixLink.GT.0) THEN
              IF (SFdepoFixLinks(iSFfixLink,3).LT.0) EXIT !skip double mirroring for 120 deg case
            END IF
            IF (SFfixIdx.GT.0) THEN
              IF ( SFdepoFixesBounds(SFfixIdx,1,1).GT.ShiftedPart(1) .OR. ShiftedPart(1).GT.SFdepoFixesBounds(SFfixIdx,2,1) .OR. &
                SFdepoFixesBounds(SFfixIdx,1,2).GT.ShiftedPart(2) .OR. ShiftedPart(2).GT.SFdepoFixesBounds(SFfixIdx,2,2) .OR. &
                SFdepoFixesBounds(SFfixIdx,1,3).GT.ShiftedPart(3) .OR. ShiftedPart(3).GT.SFdepoFixesBounds(SFfixIdx,2,3) ) THEN
                DoCycle=.TRUE.
                EXIT !do not shift this particle (-> go to next single SFfixIdx -> iSFfixLink-loop)
              END IF
              SFfixDistance = SFdepoFixesGeo(SFfixIdx,2,1)*(ShiftedPart(1)-SFdepoFixesGeo(SFfixIdx,1,1)) &
                + SFdepoFixesGeo(SFfixIdx,2,2)*(ShiftedPart(2)-SFdepoFixesGeo(SFfixIdx,1,2)) &
                + SFdepoFixesGeo(SFfixIdx,2,3)*(ShiftedPart(3)-SFdepoFixesGeo(SFfixIdx,1,3))
              SFfixIdx2=0 !init
              IF (SFfixDistance .GT. SFdepoFixesEps) THEN
                IPWRITE(*,'(I4,A,3(x,E12.5))')' original case-pos:',caseShiftedPart
                IPWRITE(*,'(I4,A,3(x,E12.5))')' current pos:',ShiftedPart
                IPWRITE(*,'(I4,7(A,I0))') &
                  ' iCase: ',iCase,', iSFfix:',iSFfix,', iSFfixLink:',iSFfixLink,', iTwin:',iTwin,', iLinkRec:',iLinkRecursive,&
                  ', SFfixIdx:',SFfixIdx,', PartIdx:',PartIdx
                CALL abort(&
                  __STAMP__ &
                  ,'Particle is outside of SF-Fix-Plane! (For Layer-/Resample-Parts: try -UseFixBounds)',SFfixIdx,SFfixDistance)
              ELSE IF ( (SFfixDistance.LT.-r_sf) ) THEN !.OR. (SFfixDistance.GT.0.) ) THEN !SFfixDistance>0 are particle within eps
                DoNotDeposit=.TRUE. !too far inside so that mirrored SF would not reach any DOF
              ELSE IF (iSFfixLink.GT.0) THEN !check which is the other SFfixIdx of current link
                IF (SFfixIdx.EQ.SFdepoFixLinks(iSFfixLink,1)) THEN
                  SFfixIdx2=SFdepoFixLinks(iSFfixLink,2)
                ELSE IF (SFfixIdx.EQ.SFdepoFixLinks(iSFfixLink,2)) THEN
                  SFfixIdx2=SFdepoFixLinks(iSFfixLink,1)
                ELSE
                  CALL abort(&
                    __STAMP__ &
                    ,'Something is wrong with iSFfixLink',iSFfixLink)
                END IF
              END IF
              ShiftedPart(1:3) = ShiftedPart(1:3) - 2.*SFfixDistance*SFdepoFixesGeo(SFfixIdx,2,1:3)
              Fac = Fac * SFdepoFixesChargeMult(SFfixIdx)
              IF (SFfixIdx2.NE.0) THEN !check if new position would not reach a dof because of the other plane
                SFfixDistance2 = SFdepoFixesGeo(SFfixIdx2,2,1)*(ShiftedPart(1)-SFdepoFixesGeo(SFfixIdx2,1,1)) &
                  + SFdepoFixesGeo(SFfixIdx2,2,2)*(ShiftedPart(2)-SFdepoFixesGeo(SFfixIdx2,1,2)) &
                  + SFdepoFixesGeo(SFfixIdx2,2,3)*(ShiftedPart(3)-SFdepoFixesGeo(SFfixIdx2,1,3))
                IF (SFfixDistance2 .GT. r_sf) DoNotDeposit=.TRUE. !too far outside of plane
              END IF
            END IF
            IF (DoNotDeposit) CYCLE !(-> do not deposit but save position for possible further recursive mirroring)
            !------------- actual deposition:
            SELECT CASE(TRIM(DepositionType))
            CASE('shape_function')
                ! Do nothing in FLEXI
            CASE('shape_function_simple')
                ! Do nothing in FLEXI
            END SELECT
          END DO ! iLinkRecursive
          IF (DoCycle) EXIT
        END DO ! iTwin
      END DO ! iSFfixLink
    END DO ! iSFfix
  END DO ! iCase (periodicity)
END IF !NbrOfSFdepoFixes

END SUBROUTINE calcSfSource


SUBROUTINE ReadTimeAverage(FileName)
!===================================================================================================================================
! Read in ChargeDensity and save to PartSource
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_PreProc
USE MOD_IO_HDF5
USE MOD_HDF5_Input,              ONLY:ReadArray,ReadAttribute,File_ID,OpenDataFile,CloseDataFile,DatasetExists
USE MOD_Particle_Vars,           ONLY:nSpecies
USE MOD_PICDepo_Vars,            ONLY:PartSource
USE MOD_Mesh_Vars,               ONLY:OffsetElem,nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE         :: U(:,:,:,:,:)
INTEGER                  :: iSpec, iElem, kk, ll, mm
INTEGER                  :: Rank
INTEGER                  :: nVars, iVar, N_HDF5
INTEGER,ALLOCATABLE      :: PartSourceToVar(:)
INTEGER(HID_T)           :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
LOGICAL                  :: SolutionExists
CHARACTER(LEN=255),ALLOCATABLE          :: VarNames(:)
CHARACTER(LEN=10)        :: strhelp
#if USE_MPI
REAL                     :: StartT,EndT
#endif /*MPI*/
!===================================================================================================================================

  SWRITE(UNIT_stdOut,*)'Using TimeAverage as constant PartSource(4) from file:',TRIM(FileName)
#if USE_MPI
  StartT=MPI_WTIME()
#endif

  IF(MPIRoot) THEN
    IF(.NOT.FILEEXISTS(FileName))  CALL abort(__STAMP__, &
          'TimeAverage-File "'//TRIM(FileName)//'" does not exist',999,999.)
  END IF
#if USE_MPI
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
#else
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#endif
  ! get attributes
  CALL DatasetExists(File_ID,'DG_Solution',SolutionExists)
  IF(MPIRoot)THEN
    IF(.NOT.SolutionExists)  CALL abort(&
      __STAMP__&
      ,'DG_Solution in TimeAverage-File "'//TRIM(FileName)//'" does not exist!')
  END IF
  CALL H5DOPEN_F(File_ID, 'DG_Solution', Dset_ID, iError)
  ! Get the data space of the dataset.
  CALL H5DGET_SPACE_F(Dset_ID, FileSpace, iError)
  ! Get number of dimensions of data space
  CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, Rank, iError)
  ! Get size and max size of data space
  Dims   =0
  DimsMax=0
  CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Dims(1:Rank), DimsMax(1:Rank), iError)
  CALL H5SCLOSE_F(FileSpace, iError)
  CALL H5DCLOSE_F(Dset_ID, iError)
  IF(MPIRoot)THEN
    IF(INT(Dims(Rank),4).NE.nGlobalElems)  CALL abort(&
      __STAMP__&
      ,' MeshSize and Size of TimeAverage-File "'//TRIM(FileName)//'" does not match!')
  END IF
  nVars=INT(Dims(1),4)
  ALLOCATE(VarNames(nVars))
  CALL ReadAttribute(File_ID,'VarNames',nVars,StrArray=VarNames)

  ALLOCATE(PartSourceToVar(nSpecies))
  PartSourceToVar=0
  DO iSpec=1,nSpecies
    WRITE(strhelp,'(I2.2)') iSpec
    DO iVar=1,nVars
      IF (VarNames(iVar).EQ.TRIM('ChargeDensity-Spec')//TRIM(strhelp)) THEN
        PartSourceToVar(iSpec)=iVar
        EXIT
      END IF
    END DO
  END DO
  IF (.NOT.ANY(PartSourceToVar.NE.0)) CALL abort(__STAMP__, &
    'No PartSource found in TimeAverage-File "'//TRIM(FileName)//'"!!!',999,999.)
  DEALLOCATE(VarNames)

  !-- read state
  ALLOCATE(U(nVars,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  CALL ReadAttribute(File_ID,'N',1,IntScalar=N_HDF5)
  IF(N_HDF5.EQ.PP_N)THEN! No interpolation needed, read solution directly from file
    CALL ReadArray('DG_Solution',5,(/nVars,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
  ELSE
    CALL abort(__STAMP__, &
          'N_HDF5.NE.PP_N !',999,999.)
  END IF
  CALL CloseDataFile() 

  !-- save to PartSource
  PartSource(4,:,:,:,:)=0.
  DO iSpec=1,nSpecies
    IF (PartSourceToVar(iSpec).NE.0) THEN
      DO iElem=1,PP_nElems
        DO kk = 0, PP_N
          DO ll = 0, PP_N
            DO mm = 0, PP_N
              PartSource(4,mm,ll,kk,iElem)=PartSource(4,mm,ll,kk,iElem)+U(PartSourceToVar(iSpec),mm,ll,kk,iElem)
            END DO
          END DO
        END DO
      END DO
    END IF
  END DO
  DEALLOCATE(U)
  DEALLOCATE(PartSourceToVar)

#if USE_MPI
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' Readin took  [',EndT-StartT,'s].'
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE!' 
#else
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE!' 
#endif

END SUBROUTINE ReadTimeAverage


SUBROUTINE FinalizeDeposition() 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize pic deposition
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars,   ONLY:Geo
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(PartSource)
SDEALLOCATE(GaussBorder)
SDEALLOCATE(ElemDepo_xGP)
SDEALLOCATE(Vdm_EquiN_GaussN)
SDEALLOCATE(Knots)
SDEALLOCATE(GaussBGMIndex)
SDEALLOCATE(GaussBGMFactor)
SDEALLOCATE(GEO%PeriodicBGMVectors)
SDEALLOCATE(BGMSource)
SDEALLOCATE(GPWeight)
SDEALLOCATE(ElemRadius2_sf)
SDEALLOCATE(Vdm_NDepo_GaussN)
SDEALLOCATE(DDMassInv)
SDEALLOCATE(XiNDepo)
SDEALLOCATE(swGPNDepo)
SDEALLOCATE(wBaryNDepo)
SDEALLOCATE(NDepochooseK)
SDEALLOCATE(tempcharge)
SDEALLOCATE(CellVolWeightFac)
SDEALLOCATE(CellVolWeight_Volumes)
END SUBROUTINE FinalizeDeposition

END MODULE MOD_PICDepo
