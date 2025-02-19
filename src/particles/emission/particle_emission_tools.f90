!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
! helper functions for particle emission
!===================================================================================================================================
MODULE MOD_Particle_Emission_Tools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: IntegerDivide
PUBLIC:: InsideExcludeRegionCheck
PUBLIC:: SetParticleMass
PUBLIC:: SetCellLocalParticlePosition
PUBLIC:: SetParticlePositionPoint
PUBLIC:: SetParticlePositionEquidistLine
PUBLIC:: SetParticlePositionLine
PUBLIC:: SetParticlePositionPlane
PUBLIC:: SetParticlePositionDisk
PUBLIC:: SetParticlePositionCross
PUBLIC:: SetParticlePositionCircle
PUBLIC:: SetParticlePositionCuboidCylinder
PUBLIC:: SetParticlePositionSphere
! PUBLIC:: SetParticlePositionSinDeviation
PUBLIC:: SetParticlePositionGaussian
PUBLIC:: SetParticlePositionFromFile
PUBLIC:: SamplePoissonDistri
!===================================================================================================================================
CONTAINS


SUBROUTINE IntegerDivide(Ntot,length,Ai,Ni)
!===================================================================================================================================
! Divide the Integer Ntot into separate Ni inside different "areas" Ai (attention: old Ni is counted up -> needs to be initialized!)
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: Ntot, length
REAL,INTENT(IN)                  :: Ai(1:length)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: Ni(1:length)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iN, iRan, Nitemp, Nrest, Ntot0
REAL            :: Atot, Bi(0:length), RandVal1, A2i(1:length), A2tot !,Error,Nrel(1:length),Arel(1:length)
!===================================================================================================================================

IF (Ntot.EQ.0) RETURN

Atot  = 0.
Ntot0 = 0
DO iN = 1,length
  Atot  = Atot +Ai(iN)
  Ntot0 = Ntot0+Ni(iN)
END DO

!-- divide into INT-parts
Nrest = Ntot
A2tot = 0.
Bi(:) = 0.

DO iN = 1,length
  Nitemp  = INT(REAL(Ai(iN))/REAL(Atot)*Ntot)     ! INT-part
  Ni(iN)  = Ni(iN) + Nitemp
  Nrest   = Nrest  - Nitemp                       ! remaining number
  A2i(iN) = REAL(Ai(iN))/REAL(Atot)*Ntot - Nitemp ! elem weight for remaining number
  A2tot   = A2tot + A2i(iN)
  Bi(iN)  = A2tot                                 ! elem upper limit for remaining number
END DO

!-- distribute remaining number
IF (Nrest.LT.0) CALL Abort(__STAMP__,'ERROR 1 in IntegerDivide!')

IF (Nrest.GT.0) THEN
  DO iN = 1,length
    Bi(iN) = Bi(iN)/A2tot                         ! normalized upper limit
  END DO

  DO iRan = 1,Nrest
    CALL RANDOM_NUMBER(RandVal1)
    DO iN = 1,length
      IF (Bi(iN-1).LT.RandVal1 .AND. RandVal1.LE.Bi(iN) ) THEN
        Ni(iN) = Ni(iN)+1
        EXIT
      END IF
    END DO
  END DO
END IF

!-- test if remaining number was distributed
Nrest = Ntot + Ntot0
DO iN = 1,length
  Nrest = Nrest-Ni(iN)
END DO

IF (Nrest.NE.0) THEN
  IPWRITE(UNIT_stdOut,'(I0,A,I0)') 'Ntot:  ',Ntot
  IPWRITE(UNIT_stdOut,'(I0,A,I0)') 'Ntot0: ',Ntot0
  IPWRITE(UNIT_stdOut,'(I0,A,I0)') 'Nrest: ',Nrest
  CALL Abort(__STAMP__,'ERROR 2 in IntegerDivide!')
END IF

END SUBROUTINE IntegerDivide


SUBROUTINE SetParticleMass(FractNbr,NbrOfParticle)
!===================================================================================================================================
! Add particles mass and charge
!===================================================================================================================================
! MODULES
USE MOD_Globals,          ONLY: Abort
USE MOD_Particle_Vars,    ONLY: PartSpecies,Species,PartState,doRandomPartDiam
USE MOD_Particle_Tools,   ONLY: GetNextFreePosition
#if USE_SPHERICITY
USE MOD_Particle_Vars,    ONLY: doRandomSphericity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: i,PositionNbr
REAL                                     :: randnum
!===================================================================================================================================

DO i = 1,NbrOfParticle
  PositionNbr = GetNextFreePosition(i)
  PartSpecies(PositionNbr) = FractNbr
  IF (doRandomPartDiam) THEN
    CALL RANDOM_NUMBER(randnum)
    randnum = randnum*2.-1.
    PartState(PART_DIAM, PositionNbr) = (Species(FractNbr)%DiameterIC * Species(FractNbr)%ScalePartDiam + FLOOR(randnum * &
                 Species(FractNbr)%PartDiamVarianceIC*Species(FractNbr)%ScalePartDiam))/Species(FractNbr)%ScalePartDiam
  ELSE
    PartState(PART_DIAM, PositionNbr) = Species(FractNbr)%DiameterIC
  END IF
#if USE_SPHERICITY
  IF (doRandomSphericity) THEN
    CALL RANDOM_NUMBER(randnum)
    randnum = randnum*2.-1.
    PartState(PART_SPHE, PositionNbr) = (Species(FractNbr)%SphericityIC * 100 + FLOOR(randnum * &
                                         Species(FractNbr)%PartSpheVarianceIC*100))*0.01
  ELSE
    PartState(PART_SPHE, PositionNbr) = Species(FractNbr)%SphericityIC
  END IF
#endif
END DO

END SUBROUTINE SetParticleMass


SUBROUTINE SamplePoissonDistri(RealTarget,IntSample,Flag_opt)
!===================================================================================================================================
! Sample IntSample from Poisson-Distri around RealTarget (if Flag present it will be turned off at sample limit, otherwise abort)
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: RealTarget
LOGICAL,INTENT(INOUT),OPTIONAL :: Flag_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: IntSample
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL         :: Flag
INTEGER         :: Npois
REAL            :: Tpois, RandVal1
!===================================================================================================================================

IF (PRESENT(Flag_opt)) THEN; Flag = Flag_opt
ELSE                       ; Flag = .FALSE.
END IF
Npois = 0
Tpois = 1.0
CALL RANDOM_NUMBER(RandVal1)

! Continue looping until found a valid sample or ran into an error
DO
  Tpois = RandVal1*Tpois
  IF (Tpois.LT.TINY(Tpois)) THEN
    ! Turn off Poisson Sampling and "sample" by random-rounding
    IF (Flag) THEN
      IPWRITE(UNIT_stdOut,'(I0,A)') ' WARNING: target is too large for poisson sampling: switching now to Random rounding...'
      IntSample = INT(RealTarget + RandVal1)
      Flag      = .FALSE.
      EXIT
    ! Turning off not allowed: abort (RealTarget must be decreased ot PoissonSampling turned off manually)
    ELSE
      CALL Abort(__STAMP__,'ERROR in SamplePoissonDistri: RealTarget (e.g. flux) is too large for poisson sampling!')
    END IF
  END IF

  ! Invalid randVal draw, try again
  IF (Tpois.GT.EXP(-RealTarget)) THEN
    Npois = Npois+1
    CALL RANDOM_NUMBER(RandVal1)
  ELSE
    IntSample = Npois
    EXIT
  END IF
END DO

END SUBROUTINE SamplePoissonDistri


SUBROUTINE SetCellLocalParticlePosition(chunkSize,iSpec,iInit,UseExactPartNum)
!===================================================================================================================================
!> routine for inserting particles positions locally in every cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
USE MOD_Particle_Mesh_Vars     ,ONLY: BoundsOfElem_Shared,ElemVolume_Shared
USE MOD_Particle_Tracking      ,ONLY: ParticleInsideCheck
USE MOD_Particle_Vars          ,ONLY: Species,PDM,PEM,PartState
USE MOD_Particle_Tools         ,ONLY: GetNextFreePosition
USE MOD_Particle_Tools         ,ONLY: IncreaseMaxParticleNumber
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: offsetComputeNodeElem
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
INTEGER, INTENT(IN)              :: iInit
LOGICAL, INTENT(IN)              :: UseExactPartNum
INTEGER, INTENT(INOUT)           :: chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iElem, iChunkSize
INTEGER                          :: iPart, nPart
REAL                             :: iRan, RandomPos(3)
REAL                             :: PartDens
LOGICAL                          :: InsideFlag
INTEGER                          :: CellChunkSize(1:nElems)
INTEGER                          :: chunkSize_tmp,ParticleIndexNbr
INTEGER                          :: offsetElemCNProc
!===================================================================================================================================
#if USE_MPI
! ElemVolume is only built for local DG elements. Therefore, array is only filled for elements on the same compute node
offsetElemCNProc = offsetElem - offsetComputeNodeElem
#else
offsetElemCNProc = 0
#endif  /*USE_MPI*/

IF (UseExactPartNum) THEN
  IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.0) CALL IncreaseMaxParticleNumber(chunkSize)

  CellChunkSize(:) = 0
  CALL IntegerDivide(chunkSize,nElems,ElemVolume_Shared(1+offsetElemCNProc:nElems+offsetElemCNProc),CellChunkSize(:))
ELSE
  ! numerical PartDensity is needed
  PartDens      = Species(iSpec)%Init(iInit)%PartDensity
  chunkSize_tmp = INT(PartDens * LocalVolume)
  IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.0) CALL IncreaseMaxParticleNumber(chunkSize_tmp)
END IF

iChunkSize       = 1
ParticleIndexNbr = 1

DO iElem = 1,nElems
  ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem+offsetElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
    IF (UseExactPartNum) THEN
      nPart = CellChunkSize(iElem)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      nPart = INT(PartDens * ElemVolume_Shared(iElem+offsetElemCNProc) + iRan)
    END IF

    DO iPart = 1,nPart
      ParticleIndexNbr = GetNextFreePosition(iChunkSize)
      InsideFlag       = .FALSE.

      DO WHILE(.NOT.InsideFlag)
        CALL RANDOM_NUMBER(RandomPos)
        RandomPos  = Bounds(1,:) + RandomPos*(Bounds(2,:)-Bounds(1,:))
        InsideFlag = ParticleInsideCheck(RandomPos,ParticleIndexNbr,iElem+offsetElem)
      END DO

      PartState(     1:3,ParticleIndexNbr) = RandomPos(1:3)
      PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
      PDM%IsNewPart(     ParticleIndexNbr) = .TRUE.
      PEM%Element(       ParticleIndexNbr) = iElem+offsetElem
      iChunkSize = iChunkSize + 1
    END DO
  END ASSOCIATE
END DO
chunkSize = ichunkSize - 1

END SUBROUTINE SetCellLocalParticlePosition


SUBROUTINE SetParticlePositionPoint(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr,iInit,chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3)
INTEGER                 :: i
!===================================================================================================================================

Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC
DO i = 1,chunkSize
   particle_positions(i*3-2) = Particle_pos(1)
   particle_positions(i*3-1) = Particle_pos(2)
   particle_positions(i*3  ) = Particle_pos(3)
END DO

END SUBROUTINE SetParticlePositionPoint


SUBROUTINE SetParticlePositionEquidistLine(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr,iInit,chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),VectorGap(3)
INTEGER                 :: i
!===================================================================================================================================

IF (chunkSize.EQ.1) THEN
  Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + 0.5 * Species(FractNbr)%Init(iInit)%BaseVector1IC
ELSE
  VectorGap = Species(FractNbr)%Init(iInit)%BaseVector1IC/(REAL(chunkSize)-1.)
  DO i = 1,chunkSize
    Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + (i-1)*VectorGap
    particle_positions(i*3-2) = Particle_pos(1)
    particle_positions(i*3-1) = Particle_pos(2)
    particle_positions(i*3  ) = Particle_pos(3)
  END DO
END IF

END SUBROUTINE SetParticlePositionEquidistLine


SUBROUTINE SetParticlePositionLine(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr,iInit,chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVal
INTEGER                 :: i
!===================================================================================================================================

DO i = 1,chunkSize
  CALL RANDOM_NUMBER(RandVal)
  Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC*RandVal

  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)
END DO

END SUBROUTINE SetParticlePositionLine


SUBROUTINE SetParticlePositionCross(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! Modules
USE MOD_Particle_Globals       ,ONLY: VECNORM,FindLinIndependentVectors
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr,iInit,chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVec,lineVector(3),lineVector2(3),frac
INTEGER                 :: i
!===================================================================================================================================

CALL FindLinIndependentVectors(Species(FractNbr)%Init(iInit)%NormalIC(1:3),lineVector(1:3),lineVector2(1:3))
CALL GramSchmidtAlgo(Species(FractNbr)%Init(iInit)%NormalIC(1:3),lineVector(1:3),lineVector2(1:3))

! first particle is at (0,0)
frac = 1./(chunkSize*0.5-1.)

DO i = 1,chunkSize
  RandVec      = MERGE(2.*(i-1.)*frac-1.      ,2.*(i-chunkSize*0.5-1.)*frac-1.,i.LE.chunkSize*0.5)
  Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * RandVec
  Particle_pos = MERGE(Particle_pos*lineVector,Particle_pos*lineVector2       ,i.LE.chunkSize*0.5)

 particle_positions(i*3-2) = Particle_pos(1)
 particle_positions(i*3-1) = Particle_pos(2)
 particle_positions(i*3  ) = Particle_pos(3)

END DO

END SUBROUTINE SetParticlePositionCross


SUBROUTINE SetParticlePositionPlane(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)  :: chunkSize
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVal(2)
INTEGER                 :: i
!===================================================================================================================================

DO i = 1,chunkSize
  CALL RANDOM_NUMBER(RandVal)
  Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1)
  Particle_pos = Particle_pos                              + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2)

  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)
END DO

END SUBROUTINE SetParticlePositionPlane


SUBROUTINE SetParticlePositionDisk(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! Modules
USE MOD_Particle_Globals       ,ONLY: VECNORM,FindLinIndependentVectors!,RandNormal
USE MOD_Particle_Vars          ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVec(2),lineVector(3),lineVector2(3),radius
INTEGER                 :: i
!===================================================================================================================================

CALL FindLinIndependentVectors(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
CALL GramSchmidtAlgo(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))

DO i = 1,chunkSize
  radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.

  DO WHILE(radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC)
    !  RandVec(1) = MAX(MIN(RandNormal(0.0,1.0),1.0),-1.)
    !  RandVec(2) = MAX(MIN(RandNormal(0.0,1.0),1.0),-1.)
     CALL RANDOM_NUMBER(RandVec)
     RandVec      = RandVec * 2. - 1.
     Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * &
                   (RandVec(1) * lineVector + RandVec(2) *lineVector2)

     radius = VECNORM((Particle_pos(1:3)-Species(FractNbr)%Init(iInit)%BasePointIC(1:3)))
  END DO

  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)
END DO

END SUBROUTINE SetParticlePositionDisk


SUBROUTINE SetParticlePositionCircle(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_Particle_Globals       ,ONLY: FindLinIndependentVectors
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVal,lineVector(3),lineVector2(3),radius,Phi
INTEGER                 :: i
!===================================================================================================================================

CALL FindLinIndependentVectors(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
CALL GramSchmidtAlgo(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
radius = Species(FractNbr)%Init(iInit)%RadiusIC

DO i = 1,chunkSize
  IF(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle') THEN
    CALL RANDOM_NUMBER(RandVal)
    Phi = 2.*PI*RandVal
  ELSE
    Phi = 2.*PI*REAL(i)/ REAL(chunkSize)
  END IF
  Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + &
                linevector  * COS(Phi) * radius            + &
                linevector2 * SIN(Phi) * radius
  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)

END DO

END SUBROUTINE SetParticlePositionCircle


SUBROUTINE SetParticlePositionCuboidCylinder(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Timedisc_Vars ,ONLY: RKdtFrac
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Timedisc_Vars          ,ONLY: dt
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)  :: chunkSize
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVal(3),lineVector(3),radius
INTEGER                 :: i,chunkSize2
LOGICAL                 :: insideExcludeRegion
!===================================================================================================================================

lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)

! Sanity check line vectors
IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
  CALL Abort(__STAMP__,'BaseVectors are parallel!')
ELSE
  lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + lineVector(3) * lineVector(3))
END IF

chunkSize2 = 0
DO i = 1,chunkSize
  SELECT CASE (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
    CASE ('cuboid')
      CALL RANDOM_NUMBER(RandVal)
      Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1)
      Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2)
      ! Height directly calculated by timestep
      IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
      ELSE
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC       * RandVal(3)
      END IF

    CASE ('cylinder')
      radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
      DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC))!.OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
        CALL RANDOM_NUMBER(RandVal)
        Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC * (RandVal(1)*2.-1.) &
                     + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2.-1.)
        radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                       Particle_pos(2) * Particle_pos(2) + &
                       Particle_pos(3) * Particle_pos(3) )
      END DO
      Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
      ! Height directly calculated by timestep
      IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
      ELSE
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal(3)
      END IF
  END SELECT

  ! Check if emission was inside an exclude region
  IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
    CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
    ! Particle is in excluded region
    IF (insideExcludeRegion) CYCLE
  END IF

  particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
  particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
  particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)

  chunkSize2 = chunkSize2 + 1
END DO

chunkSize = chunkSize2

END SUBROUTINE SetParticlePositionCuboidCylinder


SUBROUTINE SetParticlePositionSphere(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Particle_Tools         ,ONLY: DICEUNITVECTOR
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)  :: chunkSize
REAL,INTENT(OUT)        :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVal,radius
INTEGER                 :: i,chunkSize2
LOGICAL                 :: insideExcludeRegion
!===================================================================================================================================

chunkSize2 = 0
DO i = 1,chunkSize
  CALL RANDOM_NUMBER(RandVal)
  radius       = Species(FractNbr)%Init(iInit)%RadiusIC*RandVal**(1./3.)
  Particle_pos = DICEUNITVECTOR()*radius + Species(FractNbr)%Init(iInit)%BasePointIC

  ! Check if emission was inside an exclude region
  IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
    CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
    ! Particle is in excluded region
    IF (insideExcludeRegion) CYCLE
  END IF

  particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
  particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
  particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)

  chunkSize2=chunkSize2 + 1
END DO
chunkSize = chunkSize2

END SUBROUTINE SetParticlePositionSphere


SUBROUTINE SetParticlePositionGaussian(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! position particle along line by drawing it from a Gaussian distribution
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVal(3),lineVector(3),lineVector2(3),radius
REAL                    :: argumentTheta,pdf,norm_pdf
INTEGER                 :: i
!===================================================================================================================================

! Find normal vector direction and length.
IF (Species(FractNbr)%Init(iInit)%NormalIC(1).EQ.0.) THEN
  IF (Species(FractNbr)%Init(iInit)%NormalIC(2).EQ.0.) THEN
    IF (Species(FractNbr)%Init(iInit)%NormalIC(3).EQ.0.) &
      CALL Abort(__STAMP__,'Error in SetParticlePosition, NormalIC should not be zero')
    lineVector(1:2)   = 1.0
    lineVector(3)     = 0.0
  ELSE
    IF (Species(FractNbr)%Init(iInit)%NormalIC(3).EQ.0.) THEN
      lineVector(1)   = 1.0
      lineVector(3)   = 1.0
      lineVector(2)   = 0.0
    ELSE
      lineVector(1)   = 1.0
      lineVector(2:3) = 0.0
    END IF
  END IF
ELSE
  IF (Species(FractNbr)%Init(iInit)%NormalIC(2).EQ.0.) THEN
    IF (Species(FractNbr)%Init(iInit)%NormalIC(3).EQ.0.) THEN
      lineVector(2:3) = 1.0
      lineVector(1)   = 0.0
    ELSE
      lineVector(2)   = 1.0
      lineVector(1)   = 0.0
      lineVector(3)   = 0.0
    END IF
  ELSE
    IF (Species(FractNbr)%Init(iInit)%NormalIC(3).EQ.0.) THEN
      lineVector(3)   = 1.0
      lineVector(1:2) = 0.0
    ELSE
       lineVector(2)  = 0.0
       lineVector(3)  = -Species(FractNbr)%Init(iInit)%NormalIC(1)
       lineVector(1)  = Species(FractNbr)%Init(iInit)%NormalIC(3)
    END IF
  END IF
END IF

! Normalize lineVector to unit length
lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))

! calculate lineVector cross product with normal vector initial condition
lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
                 Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
                 Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
                 Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

! Normalize lineVector2 to unit length
lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
              lineVector2(2) + lineVector2(3) * lineVector2(3))

DO i=1,chunkSize
  CALL RANDOM_NUMBER(RandVal)
  ! expand RandVal from [0,1] to [0,2PI]
  ! x: normalized Gaussian distribution
  norm_pdf = COS(2.*PI*RandVal(1))*SQRT(-2.*LOG(RandVal(2)))
  ! z = mu + std*x, mu=0.0
  pdf    = norm_pdf * Species(FractNbr)%Init(iInit)%BaseVariance
  radius = MIN(ABS(pdf*Species(FractNbr)%Init(iInit)%RadiusIC),Species(FractNbr)%Init(iInit)%RadiusIC)
  radius = MAX(0.,radius)-MIN(0.,radius)
  argumentTheta = 2.*PI*RandVal(3)

  ! position particle at random angle
  Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +  &
                 linevector  * cos(argumentTheta) * radius +  &
                 linevector2 * sin(argumentTheta) * radius

  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)
END DO

END SUBROUTINE SetParticlePositionGaussian


SUBROUTINE SetParticlePositionFromFile(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! position particle along line by drawing it from a Gaussian distribution
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_Particle_Globals       ,ONLY: FindLinIndependentVectors
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),RandVal(2),lineVector(3),lineVector2(3),radius
REAL                    :: argumentTheta
INTEGER                 :: argumentIndex
INTEGER                 :: i
!===================================================================================================================================

CALL FindLinIndependentVectors(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
CALL GramSchmidtAlgo(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))

DO i=1,chunkSize
  CALL RANDOM_NUMBER(RandVal)
  ! expand RandVal from [0,1] to [0,2PI]
  argumentTheta = RandVal(1) * 2.*PI
  ! expand RandVal from [0,1] to [0,nPartField]
  argumentIndex = 1 + FLOOR(RandVal(2)*Species(FractNbr)%Init(iInit)%nPartField)
  ! Calculate radius
  radius = SQRT((Species(FractNbr)%Init(iInit)%PartField(argumentIndex,1)-Species(FractNbr)%Init(iInit)%BasePointIC(2))**2+&
                (Species(FractNbr)%Init(iInit)%PartField(argumentIndex,2)-Species(FractNbr)%Init(iInit)%BasePointIC(3))**2)

  Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC  +  &
                 linevector   * COS(argumentTheta) * radius +  &
                 linevector2  * SIN(argumentTheta) * radius
  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)
END DO

END SUBROUTINE SetParticlePositionFromFile


SUBROUTINE InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
!===================================================================================================================================
! Subroutine for checking if calculated particle position would be inside user-defined ExcludeRegion (cuboid or cylinder)
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY : Abort
USE MOD_Particle_Vars,          ONLY : Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr, iInit
REAL,INTENT(IN)                  :: Particle_pos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: insideExcludeRegion
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: VecExclude(3),DistExclude
INTEGER                          :: iExclude
!===================================================================================================================================
insideExcludeRegion=.FALSE.

DO iExclude=1,Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions
  VecExclude = Particle_pos - Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC

  SELECT CASE (TRIM(Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC))

    CASE ('cuboid')
      !--check normal direction
      DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
                  + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
                  + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
      ! Particle is inside current ExcludeRegion based an normal dimensions
      IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC) &
      .AND.(DistExclude .GE. 0.) ) THEN
        insideExcludeRegion = .TRUE.
      ELSE
        insideExcludeRegion = .FALSE.
        CYCLE
      END IF

      !--check BV1 direction
      DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
                  + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
                  + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)
      ! Particle is inside current ExcludeRegion based an BV1 dimensions
      IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1)**2) &
      .AND.(DistExclude .GE. 0.) ) THEN
        insideExcludeRegion = .TRUE.
      ELSE
        insideExcludeRegion = .FALSE.
        CYCLE
      END IF

      !--check BV2 direction
      DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
                  + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
                  + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
      ! Particle is inside current ExcludeRegion based an all dimensions
      IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2)**2) &
      .AND.(DistExclude .GE. 0.) ) THEN
        insideExcludeRegion = .TRUE.
        RETURN
      ELSE
        insideExcludeRegion = .FALSE.
        CYCLE
      END IF

    CASE ('cylinder')
      !--check normal direction
      DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
                  + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
                  + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
      ! Particle is inside current ExcludeRegion based an normal dimensions
      IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC) &
      .AND.(DistExclude .GE. 0.) ) THEN
        insideExcludeRegion = .TRUE.
      ELSE
        insideExcludeRegion = .FALSE.
        CYCLE
      END IF

      !--check radial direction
      DistExclude = SQRT( VecExclude(1)**2 + VecExclude(2)**2 + VecExclude(3)**2 - DistExclude**2 )
      ! Particle is inside current ExcludeRegion based an all dimensions
      IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC) &
      .AND.(DistExclude .GE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC) ) THEN
        insideExcludeRegion = .TRUE.
        RETURN
      ELSE
        insideExcludeRegion = .FALSE.
        CYCLE
      END IF

    CASE DEFAULT
        CALL Abort(__STAMP__,'Wrong SpaceIC for ExcludeRegion!')

  END SELECT
END DO

END SUBROUTINE InsideExcludeRegionCheck


PURE SUBROUTINE GramSchmidtAlgo(Vector1, Vector2, Vector3)
!===================================================================================================================================
!> Contains the Gram Schmidt algorithm for an orthonormal basis
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT) :: Vector1(3), Vector2(3), Vector3(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! v1 = w1/||w1||
Vector1(:) = Vector1(:) / SQRT(Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2)

! v2 = w2 - <v1,w2>*v1
Vector2(:) = Vector2(:) - DOT_PRODUCT(Vector1, Vector2) * Vector1(:)
! v2 = v2/||v2||
Vector2(:) = Vector2(:) / SQRT(Vector2(1)**2 + Vector2(2)**2 + Vector2(3)**2)

! v3 = w3 - <v1,w3>*v1 - <v2,w3>*v2
Vector3(:) = Vector3(:) - DOT_PRODUCT(Vector1, Vector3) * Vector1(:) -&
                          DOT_PRODUCT(Vector2, Vector3) * Vector2(:)
! v3 = v3/||v3||
Vector3(:) = Vector3(:) / SQRT(Vector3(1)**2 + Vector3(2)**2 + Vector3(3)**2)

END SUBROUTINE GramSchmidtAlgo

END MODULE MOD_Particle_Emission_Tools
