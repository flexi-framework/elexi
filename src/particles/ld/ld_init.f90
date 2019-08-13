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
! Initialisation of LD variables!
!===================================================================================================================================
MODULE MOD_LD_Init
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitLD
  MODULE PROCEDURE InitLD
END INTERFACE
INTERFACE CalcDegreeOfFreedom
  MODULE PROCEDURE CalcDegreeOfFreedom
END INTERFACE

PUBLIC :: InitLD, CalcDegreeOfFreedom
!===================================================================================================================================

CONTAINS

SUBROUTINE InitLD()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_LD_Vars
USE MOD_Mesh_Vars,             ONLY : nElems, nSides, NGEO
USE MOD_Particle_Mesh_Vars,    ONLY: ElemBaryNGeo
USE MOD_Particle_Vars,         ONLY : PDM, Species, PartSpecies
USE MOD_ReadInTools
USE MOD_Particle_Tracking_Vars,ONLY: DoRefMapping
USE MOD_Particle_Mesh_Vars,    ONLY: PartElemToElemAndSide
#if USE_MPI
USE MOD_MPI_Vars
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem, trinum, iLocSide, iPart, iInit, Elem2
#if USE_MPI
INTEGER                 :: SumOfMPISides, EndOfMPINeighbor, iProc, OffsetInnerAndBCSides
#endif
!CHARACTER(32)           :: hilf
!===================================================================================================================================

  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' LD INIT ...'

  IF(.NOT.DoRefMapping) CALL abort(&
__STAMP__&
,' LD requires ref-mapping tracking.')
  IF(NGeo.NE.1) CALL abort(&
__STAMP__&
,' LD requires a linear mesh (NGeo=1)!')

  LD_SecantMeth%Guess = GETREAL('LD-InitialGuess','10')
  LD_SecantMeth%MaxIter = GETINT('LD-MaxIterNumForLagVelo','100')
  LD_SecantMeth%Accuracy = GETREAL('LD-AccuracyForLagVelo','0.001')
  LD_RepositionFak = GETREAL('LD-RepositionsFaktor','0.0')
  LD_RelaxationFak = GETREAL('LD-RelaxationsFaktor','0.0')
  LD_DSMC_RelaxationFak_BufferA = GETREAL('LD-DSMC-RelaxationsFaktorForBufferA','0.0')
  LD_CalcDelta_t=GETLOGICAL('LD-CFL-CalcDelta-t','.FALSE.')
  LD_CalcResidual=GETLOGICAL('LD_CalcResidual','.FALSE.')

  ALLOCATE(LD_Residual(nElems,6))
  LD_Residual(1:nElems,1) = 0.0
  LD_Residual(1:nElems,2) = 0.0
  LD_Residual(1:nElems,3) = 0.0
  LD_Residual(1:nElems,4) = 0.0
  LD_Residual(1:nElems,5) = 0.0
  LD_Residual(1:nElems,6) = 0.0

  ALLOCATE(TempDens(nElems))

  ALLOCATE(BulkValues(nElems))
  BulkValues(1:nElems)%CellV(1)        = 0.0
  BulkValues(1:nElems)%CellV(2)        = 0.0
  BulkValues(1:nElems)%CellV(3)        = 0.0
  BulkValues(1:nElems)%DegreeOfFreedom = 0.0 
  BulkValues(1:nElems)%Beta            = 0.0  
  BulkValues(1:nElems)%MassDens        = 0.0
  BulkValues(1:nElems)%BulkTemperature = 0.0
  BulkValues(1:nElems)%DynamicVisc     = 0.0
  BulkValues(1:nElems)%ThermalCond     = 0.0

  ALLOCATE(BulkValuesOpenBC(nElems))
  BulkValuesOpenBC(1:nElems)%CellV(1)        = 0.0
  BulkValuesOpenBC(1:nElems)%CellV(2)        = 0.0
  BulkValuesOpenBC(1:nElems)%CellV(3)        = 0.0
  BulkValuesOpenBC(1:nElems)%DegreeOfFreedom = 0.0 
  BulkValuesOpenBC(1:nElems)%Beta            = 0.0  
  BulkValuesOpenBC(1:nElems)%MassDens        = 0.0
  BulkValuesOpenBC(1:nElems)%DynamicVisc     = 0.0
  BulkValuesOpenBC(1:nElems)%ThermalCond     = 0.0

  ALLOCATE(SurfLagValues(6,nElems,2))
  SurfLagValues(1:6,1:nElems,1:2)%LagVelo    = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaM(1)  = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaM(2)  = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaM(3)  = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaE     = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagNormVec(1) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagNormVec(2) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagNormVec(3) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagTangVec(1,1) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagTangVec(1,2) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagTangVec(1,3) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagTangVec(2,1) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagTangVec(2,2) = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%LagTangVec(2,3) = 0.0

!--- calculate cellcenter distance for viscousity terms

! P.O. should be computed AFTER construction of HALO mesh DEBUG

!--- end calculate cellcenter distance for viscousity terms

  ALLOCATE(IsDoneLagVelo(nSides))
  IsDoneLagVelo(1:nSides)   = .FALSE.
  DO iElem=1, nElems
    DO iLocSide = 1, 6
      DO trinum=1, 2
        CALL CalcLagNormVec(iLocSide, iElem, trinum)
!--- calculate cellcenter distance for viscousity terms
        Elem2=PartElemToElemAndSide(1,ilocSide,iElem) 
        IF(Elem2.EQ.-1) CYCLE
        MeanSurfValues(iLocSide, iElem)%CellCentDist(1) = ElemBaryNGeo(1,iElem) - ElemBaryNGeo(1,Elem2)
        MeanSurfValues(iLocSide, iElem)%CellCentDist(2) = ElemBaryNGeo(2,iElem) - ElemBaryNGeo(2,Elem2)
        MeanSurfValues(iLocSide, iElem)%CellCentDist(3) = ElemBaryNGeo(3,iElem) - ElemBaryNGeo(3,Elem2)
!--- end calculate cellcenter distance for viscousity terms  
      END DO
      CALL SetMeanSurfValues(iLocSide, iElem)
    END DO
  END DO
  ALLOCATE(PartStateBulkValues(PDM%maxParticleNumber,5))
  ALLOCATE(LD_RHS(PDM%maxParticleNumber,3))
  LD_RHS = 0.0
  ALLOCATE(LD_DSMC_RHS(PDM%maxParticleNumber,3))
  LD_DSMC_RHS = 0.0

! Set Particle Bulk Values

  DO iPart = 1, PDM%maxParticleNumber
    IF (PDM%ParticleInside(ipart)) THEN
      iInit = PDM%PartInit(iPart)
      PartStateBulkValues(iPart,1) = Species(PartSpecies(iPart))%Init(iInit)%VeloVecIC(1) &
                                   * Species(PartSpecies(iPart))%Init(iInit)%VeloIC
      PartStateBulkValues(iPart,2) = Species(PartSpecies(iPart))%Init(iInit)%VeloVecIC(2) &
                                   * Species(PartSpecies(iPart))%Init(iInit)%VeloIC
      PartStateBulkValues(iPart,3) = Species(PartSpecies(iPart))%Init(iInit)%VeloVecIC(3) &
                                   * Species(PartSpecies(iPart))%Init(iInit)%VeloIC
      PartStateBulkValues(iPart,4) = Species(PartSpecies(iPart))%Init(iInit)%MWTemperatureIC
      PartStateBulkValues(iPart,5) = CalcDegreeOfFreedom(iPart)
    END IF
  END DO

#if USE_MPI
  SumOfMPISides = 0
  DO iProc =1, nNbProcs
    SumOfMPISides =SumOfMPISides + nMPISides_MINE_Proc(iProc) + nMPISides_YOUR_Proc(iProc)
  END DO
  OffsetInnerAndBCSides = OffsetMPISides_MINE(0) + 1
  EndOfMPINeighbor = OffsetInnerAndBCSides + SumOfMPISides - 1
  ALLOCATE(MPINeighborBulkVal(OffsetInnerAndBCSides:EndOfMPINeighbor,1:8))
#endif

  SWRITE(UNIT_stdOut,'(A)')' INIT LD DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')

  DEALLOCATE(PDM%PartInit)  ! normaly done in DSMC_ini.f90

END SUBROUTINE InitLD

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

REAL FUNCTION CalcDegreeOfFreedom(iPart)
!===================================================================================================================================
! calculation of degree of freedom per part
!===================================================================================================================================
! MODULES
  USE MOD_LD_Vars
  USE MOD_Particle_Vars,      ONLY : PartSpecies, BoltzmannConst
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
!--------------------------------------------------------------------------------------------------!
  REAL                          :: ZetaRot, ZetaVib, JToEv
  INTEGER                       :: iSpec
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN)           :: iPart
!===================================================================================================
  JToEv = 1.602176565E-19
  iSpec = PartSpecies(iPart)
  
  ZetaRot = 0.0
  ZetaVib = 0.0
  
  CalcDegreeOfFreedom = 3.0 + ZetaRot + ZetaVib

END FUNCTION CalcDegreeOfFreedom

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

REAL FUNCTION CalcTriNumArea(Vector1,Vector2)
!===================================================================================================================================
! Calculation of triangle surface area
!===================================================================================================================================
! MODULES
  !USE MOD_Particle_Vars,          ONLY : GEO
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  REAL, INTENT(IN)            :: Vector1(1:3), Vector2(1:3)
!--------------------------------------------------------------------------------------------------!

!--- 2 * Area = |cross product| of vector 1-2 and 1-3 for first triangle or
                                 ! vector 1-3 and 1-4 for second triangle
   CalcTriNumArea = 0.5*SQRT((Vector1(2)*Vector2(3)-Vector1(3)*Vector2(2))**2 &
           + (-Vector1(1)*Vector2(3)+Vector1(3)*Vector2(1))**2 &
           + (Vector1(1)*Vector2(2)-Vector1(2)*Vector2(1))**2)

  RETURN

END FUNCTION CalcTriNumArea

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE CalcLagNormVec(iLocSide, Element, trinum)
!===================================================================================================================================
! Calculation of normal vector
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_LD_Vars
  USE MOD_Mesh_Vars,              ONLY:ElemToSide,NGeo
  USE MOD_Particle_Mesh_Vars,     ONLY:XCL_NGeo,ElemBaryNGeo
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  REAL                        :: Vector1(1:3), Vector2(1:3),nVecTest
  REAL                        :: nx, ny, nz, nVal
  INTEGER                     :: flip,p,q
  REAL                        :: SideCoord(1:3,0:1,0:1)
  REAL                        :: SideCoord_tmp(1:3,0:1,0:1)
  REAL                        :: xNod1,xnod2,xnod3,ynod1,ynod2,ynod3,znod1,znod2,znod3
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  INTEGER, INTENT(IN)         :: iLocSide, Element, trinum
!--------------------------------------------------------------------------------------------------!

! normal vector for each side, pointing outside

flip=ElemToSide(E2S_FLIP,ilocSide,Element)
! master side, flip=0
! slave side,  flip=1,..,4

SELECT CASE(ilocSide)

CASE(XI_MINUS)
  DO q=0,NGeo
    DO p=0,NGeo
      SideCoord_tmp(1:3,p,q)=XCL_NGeo(1:3,0,q,p,Element)
    END DO !p
  END DO !q

CASE(XI_PLUS)
  SideCoord_tmp(1:3,:,:)=XCL_NGeo(1:3,NGeo,:,:,Element)

CASE(ETA_MINUS)
  SideCoord_tmp(1:3,:,:)=XCL_NGeo(1:3,:,0,:,Element)

CASE(ETA_PLUS)
  DO q=0,NGeo
    DO p=0,NGeo
      SideCoord_tmp(1:3,NGeo-p,q)=XCL_NGeo(1:3,p,NGeo,q,Element)
    END DO !p
  END DO !q

CASE(ZETA_MINUS)
  DO q=0,NGeo
    DO p=0,NGeo
      SideCoord_tmp(1:3,q,p)=XCL_NGeo(1:3,p,q,0,Element)
    END DO !p
  END DO !q

CASE(ZETA_PLUS)
  DO q=0,NGeo
    DO p=0,NGeo
      SideCoord_tmp(1:3,p,q)=XCL_NGeo(1:3,p,q,NGeo,Element)
    END DO !p
  END DO ! q
END SELECT

SELECT CASE(flip)
  CASE(0) ! master side
    SideCoord(:,:,:)=SideCoord_tmp
  CASE(1) ! slave side, SideID=q,jSide=p
    DO q=0,NGeo
      DO p=0,NGeo
        SideCoord(:,p,q)=SideCoord_tmp(:,q,p)
      END DO ! p
    END DO ! q
  CASE(2) ! slave side, SideID=N-p,jSide=q
    DO q=0,NGeo
      DO p=0,NGeo
        SideCoord(:,p,q)=SideCoord_tmp(:,NGeo-p,q)
      END DO ! p
    END DO ! q
  CASE(3) ! slave side, SideID=N-q,jSide=N-p
    DO q=0,NGeo
      DO p=0,NGeo
        SideCoord(:,p,q)=SideCoord_tmp(:,NGeo-q,NGeo-p)
      END DO ! p
    END DO ! q
  CASE(4) ! slave side, SideID=p,jSide=N-q
    DO q=0,NGeo
      DO p=0,NGeo
        SideCoord(:,p,q)=SideCoord_tmp(:,p,NGeo-q)
      END DO ! p
    END DO ! q
END SELECT


IF(trinum.EQ.1) THEN
!--- Node 1 ---
   xNod1 = SideCoord(1,0,0)
   yNod1 = SideCoord(2,0,0)
   zNod1 = SideCoord(3,0,0)

!--- Node 2 ---
   xNod2 = SideCoord(1,NGeo,0)
   yNod2 = SideCoord(2,NGeo,0)
   zNod2 = SideCoord(3,NGeo,0)

!--- Node 3 ---
   xNod3 = SideCoord(1,NGeo,NGeo)
   yNod3 = SideCoord(2,NGeo,NGeo)
   zNod3 = SideCoord(3,NGeo,NGeo)
ELSE
!--- Node 1 ---
   xNod1 = SideCoord(1,0,0)
   yNod1 = SideCoord(2,0,0)
   zNod1 = SideCoord(3,0,0)
!--- Node 3 ---
   xNod2 = SideCoord(1,NGeo,NGeo)
   yNod2 = SideCoord(2,NGeo,NGeo)
   zNod2 = SideCoord(3,NGeo,NGeo)
!--- Node 4 ---
   xNod3 = SideCoord(1,0,NGeo)
   yNod3 = SideCoord(2,0,NGeo)
   zNod3 = SideCoord(3,0,NGeo)
END IF

   Vector1(1) = xNod2 - xNod1
   Vector1(2) = yNod2 - yNod1
   Vector1(3) = zNod2 - zNod1

   Vector2(1) = xNod3 - xNod1
   Vector2(2) = yNod3 - yNod1
   Vector2(3) = zNod3 - zNod1

   nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2) ! n is inward normal vector
   ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
   nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

   nVal = SQRT(nx*nx + ny*ny + nz*nz)

   SurfLagValues(iLocSide, Element,trinum)%LagNormVec(1) = nx/nVal
   SurfLagValues(iLocSide, Element,trinum)%LagNormVec(2) = ny/nVal
   SurfLagValues(iLocSide, Element,trinum)%LagNormVec(3) = nz/nVal

   nVal = SQRT( Vector1(1)*Vector1(1) + Vector1(2)*Vector1(2) + Vector1(3)*Vector1(3) )
!--- first tangential Vector == Node1->Node2
   SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,1) = Vector1(1)/nVal
   SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,2) = Vector1(2)/nVal
   SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,3) = Vector1(3)/nVal
!--- second tangential Vector == |cross product| of N_Vec and Tang_2
   SurfLagValues(iLocSide, Element,trinum)%LagTangVec(2,1) = &
          SurfLagValues(iLocSide, Element,trinum)%LagNormVec(2) * SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,3) &
        - SurfLagValues(iLocSide, Element,trinum)%LagNormVec(3) * SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,2)
   SurfLagValues(iLocSide, Element,trinum)%LagTangVec(2,2) = &
        - SurfLagValues(iLocSide, Element,trinum)%LagNormVec(1) * SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,3) &
        + SurfLagValues(iLocSide, Element,trinum)%LagNormVec(3) * SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,1)
   SurfLagValues(iLocSide, Element,trinum)%LagTangVec(2,3) = &
          SurfLagValues(iLocSide, Element,trinum)%LagNormVec(1) * SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,2) &
        - SurfLagValues(iLocSide, Element,trinum)%LagNormVec(2) * SurfLagValues(iLocSide, Element,trinum)%LagTangVec(1,1)
   !IF(flip.NE.0) THEN
   !  SurfLagValues(iLocSide, Element,trinum)%LagTangVec(:,:) = -SurfLagValues(iLocSide, Element,trinum)%LagTangVec(:,:)
   !  SurfLagValues(iLocSide, Element,trinum)%LagNormVec      = -SurfLagValues(iLocSide, Element,trinum)%LagNormVec
   !END IF
   NVecTest = (SideCoord(1,0,0)-ElemBaryNGeo(1,Element) ) &
            * SurfLagValues(iLocSide, Element,trinum)%LagNormVec(1) &
            + (SideCoord(2,0,0)-ElemBaryNGeo(2,Element) ) &
            * SurfLagValues(iLocSide, Element,trinum)%LagNormVec(2) &
            + (SideCoord(3,0,0)-ElemBaryNGeo(3,Element) ) &
            * SurfLagValues(iLocSide, Element,trinum)%LagNormVec(3)
   IF (NVecTest.LT.0.0) THEN
     SurfLagValues(iLocSide, Element,trinum)%LagTangVec(:,:) = -SurfLagValues(iLocSide, Element,trinum)%LagTangVec(:,:)
     SurfLagValues(iLocSide, Element,trinum)%LagNormVec      = -SurfLagValues(iLocSide, Element,trinum)%LagNormVec
   END IF  

   SurfLagValues(iLocSide, Element, trinum)%Area = CalcTriNumArea(Vector1,Vector2)

END SUBROUTINE CalcLagNormVec
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
SUBROUTINE SetMeanSurfValues(iLocSide, Element)
!===================================================================================================================================
! Definition of surface fit
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Globals,             ONLY:CROSSNORM
  USE MOD_LD_Vars
  USE MOD_Mesh_Vars,                    ONLY:NGeo
  USE MOD_Particle_Mesh_Vars,           ONLY:XCL_NGeo,ElemBaryNGeo
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration  
  REAL                        :: Vector1(3), Vector2(3),BaseVectorS(3),nVectest
  REAL                        :: SideCoord(1:3,0:1,0:1)
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  INTEGER, INTENT(IN)         :: iLocSide, Element                                                 !
!--------------------------------------------------------------------------------------------------!

SELECT CASE(ilocSide)

CASE(XI_MINUS)
  SideCoord(1:3,:,:)=XCL_NGeo(1:3,0,:,:,Element)
CASE(XI_PLUS)
  SideCoord(1:3,:,:)=XCL_NGeo(1:3,NGeo,:,:,Element)
CASE(ETA_MINUS)
  SideCoord(1:3,:,:)=XCL_NGeo(1:3,:,0,:,Element)
CASE(ETA_PLUS)
  SideCoord(1:3,:,:)=XCL_NGeo(1:3,:,NGeo,:,Element)
CASE(ZETA_MINUS)
  SideCoord(1:3,:,:)=XCL_NGeo(1:3,:,:,0,Element)
CASE(ZETA_PLUS)
  SideCoord(1:3,:,:)=XCL_NGeo(1:3,:,:,NGeo,Element)
END SELECT

vector1(:) = SideCoord(:,NGeo,NGeo)-SideCoord(:,0,0)
vector2(:) = SideCoord(:,0,NGeo)-SideCoord(:,NGeo,0)

  BaseVectorS(1) = 0.25*SUM(SideCoord(1,:,:))
  BaseVectorS(2) = 0.25*SUM(SideCoord(2,:,:))
  BaseVectorS(3) = 0.25*SUM(SideCoord(3,:,:))

  MeanSurfValues(iLocSide, Element)%MeanNormVec(1:3) = CROSSNORM(Vector1,Vector2)


  NVecTest = DOT_PRODUCT(BaseVectorS-ElemBaryNGeo(:,Element),MeanSurfValues(ilocSide,Element)%MeanNormVec)
  IF(NVecTest.LE.0.0) MeanSurfValues(ilocSide,Element)%MeanNormVec=(-1.0)*MeanSurfValues(ilocSide,Element)%MeanNormVec

  MeanSurfValues(iLocSide, Element)%MeanBaseD = MeanSurfValues(iLocSide, Element)%MeanNormVec(1) * BaseVectorS(1) &
                                              + MeanSurfValues(iLocSide, Element)%MeanNormVec(2) * BaseVectorS(2) &
                                              + MeanSurfValues(iLocSide, Element)%MeanNormVec(3) * BaseVectorS(3)

END SUBROUTINE SetMeanSurfValues

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_Init
