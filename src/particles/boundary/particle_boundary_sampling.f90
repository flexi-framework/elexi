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
!> Samples the particle conditions upon encountering a boundary
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Sampling
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleBoundarySampling
  MODULE PROCEDURE InitParticleBoundarySampling
END INTERFACE

INTERFACE RestartParticleBoundarySampling
  MODULE PROCEDURE RestartParticleBoundarySampling
END INTERFACE

INTERFACE FinalizeParticleBoundarySampling
  MODULE PROCEDURE FinalizeParticleBoundarySampling
END INTERFACE

INTERFACE WriteSurfSampleToHDF5
  MODULE PROCEDURE WriteSurfSampleToHDF5
END INTERFACE

#if USE_MPI
INTERFACE ExchangeSurfData
  MODULE PROCEDURE ExchangeSurfData
END INTERFACE
#endif /*MPI*/

PUBLIC :: InitParticleBoundarySampling
PUBLIC :: RestartParticleBoundarySampling
PUBLIC :: FinalizeParticleBoundarySampling
PUBLIC :: WriteSurfSampleToHDF5
#if USE_MPI
PUBLIC :: ExchangeSurfData
#endif /*MPI*/
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleBoundarySampling()
!===================================================================================================================================
! init of particle boundary sampling
! default: use for sampling same polynomial degree as NGeo
! 1) mark sides for sampling
! 2) build special MPI communicator
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Basis                   ,ONLY:LegendreGaussNodesAndWeights
USE MOD_Mesh_Vars               ,ONLY:NGeo,BC,nSides,nBCSides,nBCs,BoundaryName
USE MOD_Particle_Boundary_Vars  ,ONLY:nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample,SurfMesh,SampWall,nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars  ,ONLY:SurfCOMM
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Mesh_Vars      ,ONLY:GEO
USE MOD_Particle_Surfaces       ,ONLY:EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY:DoRefMapping,TriaTracking
USE MOD_Particle_Vars           ,ONLY:nSpecies
USE MOD_ReadInTools
USE MOD_StringTools             ,ONLY:LowCase
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY:ElemInfo_Shared,SideInfo_Shared,ElemSideNodeID_Shared,NodeCoords_Shared
USE MOD_Particle_MPI_Vars       ,ONLY:PartMPI
USE MOD_Particle_MPI_Shared_Vars,ONLY:nNonUniqueGlobalSides,nComputeNodeTotalSides
#else
USE MOD_Particle_Boundary_Vars  ,ONLY:offSetSurfSide
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: p,q,iSide,SurfSideID,SideID,ElemID,LocSideID
INTEGER                                :: iSample,jSample,iBC,iSurfBC
INTEGER                                :: TriNum, Node1, Node2
INTEGER                                :: nSurfSampleBC
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,ALLOCATABLE,DIMENSION(:)          :: Xi_NGeo,wGP_NGeo
REAL                                   :: XiOut(1:2),E,F,G,D,tmp1,area,tmpI2,tmpJ2
REAL                                   :: xNod, zNod, yNod, Vector1(3), Vector2(3), nx, ny, nz
REAL                                   :: nVal, SurfaceVal
CHARACTER(20)                          :: hilf, hilfBC
CHARACTER(LEN=255),ALLOCATABLE         :: BCName(:),SurfSampleBCs(:)
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING ...'

! standard is sampling on NGeo
WRITE(UNIT=hilf,FMT='(I0)') NGeo
nSurfSample             = GETINT    ('Particles-nSurfSample',TRIM(hilf))

! get equidistant supersampling points
ALLOCATE(XiEQ_SurfSample(0:nSurfSample))
dXiEQ_SurfSample =2./REAL(nSurfSample)
DO q=0,nSurfSample
  XiEQ_SurfSample(q) = dXiEQ_SurfSample * REAL(q) - 1.
END DO

! create boundary name mapping for surfaces SurfaceBC number mapping
nSurfSampleBC  = CountOption('Particles-SurfSampleBC')
ALLOCATE(SurfSampleBCs(nSurfSampleBC))
DO iSurfBC=1,nSurfSampleBC
  SurfSampleBCs(iSurfBC) = TRIM(GETSTR('Particles-SurfSampleBC'))
END DO

! compare BCs against requested BCs
nSurfBC = 0
ALLOCATE(BCName(1:nBCs))
DO iBC=1,nBCs
  BCName=''
  ! inner side
  IF (PartBound%TargetBoundCond(iBC).EQ.-1) CYCLE
  
  ! always sample walls
  IF (PartBound%TargetBoundCond(iBC).EQ.PartBound%ReflectiveBC) THEN
    nSurfBC = nSurfBC + 1
    BCName(nSurfBC) = BoundaryName(iBC)
    CYCLE
  END IF
  
  ! check if BC is explicitly requested
  CALL lowcase(TRIM(BoundaryName(iBC)),hilfBC)
  DO iSurfBC=1,nSurfSampleBC
    IF (hilfBC.EQ.SurfSampleBCs(iSurfBC)) THEN
      nSurfBC = nSurfBC + 1
      BCName(nSurfBC) = BoundaryName(iBC)
      EXIT
    END IF
  END DO
END DO

! write back found BC names
IF (nSurfBC.GE.1) THEN
  ALLOCATE(SurfBCName(1:nSurfBC))
  DO iBC=1,nSurfBC
    SurfBCName(iBC) = BCName(iBC)
  END DO
END IF
DEALLOCATE(BCName)

! get number of BC-Sides
ALLOCATE(SurfMesh%SideIDToSurfID(1:nNonUniqueGlobalSides))
SurfMesh%SideIDToSurfID(1:nNonUniqueGlobalSides) = -1

! first sides on local proc. This comes from the DG part and is not affected by the new halo region
SurfMesh%nSides = 0
!DO iSide=1,nBCSides
!  ! inner side
!  IF(BC(iSide).EQ.0) CYCLE
!  
!  ! always sample walls
!  IF (PartBound%TargetBoundCond(BC(iSide)).EQ.PartBound%ReflectiveBC) THEN
!    SurfMesh%nSides                = SurfMesh%nSides + 1
!    SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nSides
!    CYCLE
!  END IF
!  
!  ! check if BC is explicitly requested
!  CALL lowcase(TRIM(BoundaryName(BC(iSide))),hilfBC)
!  DO iSurfBC=1,nSurfSampleBC
!    IF (hilfBC.EQ.SurfSampleBCs(iSurfBC)) THEN
!      SurfMesh%nSides                = SurfMesh%nSides + 1
!      SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nSides
!      EXIT
!    END IF
!  END DO
!END DO

! use here global side IDs as well. Simplifies things later
SurfMesh%nTotalSides = SurfMesh%nSides
DO iSide=1,nComputeNodeTotalSides
  ! inner (mortar) side
  IF(SideInfo_Shared(SIDE_BCID,iSide).EQ.0) CYCLE
  
  ! skip if element is not proc (i.e. later added to global sides)
  IF (ElemInfo_Shared(ELEM_RANK,SideInfo_Shared(SIDE_ELEMID,iSide)).NE.myRank) CYCLE
  
  ! always sample walls
  IF (PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,iSide)).EQ.PartBound%ReflectiveBC) THEN
    SurfMesh%nTotalSides           = SurfMesh%nTotalSides + 1
    SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nTotalSides
    CYCLE
  END IF
  
  ! check if BC is explicitly requested
  CALL lowcase(TRIM(BoundaryName(SideInfo_Shared(SIDE_BCID,iSide))),hilfBC)
  DO iSurfBC=1,nSurfSampleBC
    IF (hilfBC.EQ.SurfSampleBCs(iSurfBC)) THEN
      SurfMesh%nTotalSides           = SurfMesh%nTotalSides + 1
      SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nTotalSides
      EXIT
    END IF
  END DO
END DO

! then sides in halo region. These currently include all BC sides on node (incl. halo region). Massive overkill, change to ALLREDUCE
! later
SurfMesh%nTotalSides = SurfMesh%nSides
DO iSide=1,nComputeNodeTotalSides
  ! inner (mortar) side
  IF(SideInfo_Shared(SIDE_BCID,iSide).EQ.0) CYCLE
  
  ! skip if element is on proc (i.e. already added to local side)
  IF (ElemInfo_Shared(ELEM_RANK,SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.myRank) CYCLE
  
  ! always sample walls
  IF (PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,iSide)).EQ.PartBound%ReflectiveBC) THEN
    SurfMesh%nTotalSides           = SurfMesh%nTotalSides + 1
    SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nTotalSides
    CYCLE
  END IF
  
  ! check if BC is explicitly requested
  CALL lowcase(TRIM(BoundaryName(SideInfo_Shared(SIDE_BCID,iSide))),hilfBC)
  DO iSurfBC=1,nSurfSampleBC
    IF (hilfBC.EQ.SurfSampleBCs(iSurfBC)) THEN
      SurfMesh%nTotalSides           = SurfMesh%nTotalSides + 1
      SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nTotalSides
      EXIT
    END IF
  END DO
END DO

! build reverse mapping
ALLOCATE(SurfMesh%SurfSideToGlobSideMap(1:SurfMesh%nTotalSides))
SurfMesh%SurfSideToGlobSideMap(:) = -1
DO iSide = 1,nComputeNodeTotalSides
  IF (SurfMesh%SideIDToSurfID(iSide).LE.0) CYCLE
  SurfMesh%SurfSideToGlobSideMap(SurfMesh%SideIDToSurfID(iSide)) = iSide
END DO

! check if local proc as BC sides to sample
SurfMesh%SurfOnProc=.FALSE.
IF(SurfMesh%nTotalSides.GT.0) SurfMesh%SurfOnProc=.TRUE.

! communicate all local sides to get the total number of surf sampling sides
#if USE_MPI
CALL MPI_ALLREDUCE(SurfMesh%nSides,SurfMesh%nGlobalSides,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError)
#else
SurfMesh%nGlobalSides=SurfMesh%nSides
#endif
SWRITE(UNIT_stdOut,'(A,I8)') ' nGlobalSurfSides ', SurfMesh%nGlobalSides

! allocate arrays to hold data. If nSpecies > 1, one more species is added to hold the global average
nErosionVars        = 17
IF (nSpecies.EQ.1) THEN
    SurfMesh%SampSize   = nErosionVars
ELSE
    SurfMesh%SampSize   = nErosionVars*(nSpecies+1)
END IF

#if USE_MPI
! split communitator
CALL InitSurfCommunicator()
IF(SurfMesh%SurfOnProc) THEN
!  CALL GetHaloSurfMapping()
END IF
#else
SurfCOMM%MyRank        = 0
SurfCOMM%MPIRoot       = .TRUE.
SurfCOMM%nProcs        = 1
SurfCOMM%MyOutputRank  = 0
SurfCOMM%MPIOutputRoot = .TRUE.
SurfCOMM%nOutputProcs  = 1
OffSetSurfSide         = 0
#endif /*USE_MPI*/

! nothing more to do for procs without sampling surfaces on them
IF(.NOT.SurfMesh%SurfOnProc) RETURN

! allocate arrays for surface data
ALLOCATE(SampWall(1:SurfMesh%nTotalSides))
DO iSide=1,SurfMesh%nTotalSides
  ALLOCATE(SampWall(iSide)%State(1:SurfMesh%SampSize,1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%State=0.
END DO

! compute area of sub-faces
ALLOCATE(SurfMesh%SurfaceArea(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides))
SurfMesh%SurfaceArea = 0.

ALLOCATE(Xi_NGeo( 0:NGeo)  &
        ,wGP_NGeo(0:NGeo) )
CALL LegendreGaussNodesAndWeights(NGeo,Xi_NGeo,wGP_NGeo)
tmp1 = dXiEQ_SurfSample/2.0 !(b-a)/2

DO iSide = 1,nComputeNodeTotalSides
  SurfSideID = SurfMesh%SideIDToSurfID(iSide)
  
  ! ignore sides that are not sampled
  IF(SurfSideID.EQ.-1) CYCLE
  
  IF (TriaTracking) THEN
    ! TODO: This has to come from ElemInfo_Shared
    ElemID    = SideInfo_Shared(SIDE_ELEMID ,iSide)
    LocSideID = SideInfo_Shared(SIDE_LOCALID,iSide)

    SurfaceVal = 0.
    xNod = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,LocSideID,ElemID))
    yNod = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,LocSideID,ElemID))
    zNod = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,LocSideID,ElemID))

    DO TriNum = 1,2
      Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
      Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle
        Vector1(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node1,LocSideID,ElemID)) - xNod
        Vector1(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node1,LocSideID,ElemID)) - yNod
        Vector1(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node1,LocSideID,ElemID)) - zNod
        Vector2(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node2,LocSideID,ElemID)) - xNod
        Vector2(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node2,LocSideID,ElemID)) - yNod
        Vector2(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node2,LocSideID,ElemID)) - zNod
      nx = - Vector1(2) * Vector2(3) + Vector1(3) * Vector2(2) !NV (inwards)
      ny = - Vector1(3) * Vector2(1) + Vector1(1) * Vector2(3)
      nz = - Vector1(1) * Vector2(2) + Vector1(2) * Vector2(1)
      nVal = SQRT(nx*nx + ny*ny + nz*nz)
      SurfaceVal = SurfaceVal + nVal/2.
    END DO
    SurfMesh%SurfaceArea(1,1,SurfSideID) = SurfaceVal
  ELSE
    ! call here stephens algorithm to compute area
    DO jSample=1,nSurfSample
      DO iSample=1,nSurfSample
        area=0.
        tmpI2=(XiEQ_SurfSample(iSample-1)+XiEQ_SurfSample(iSample))/2. ! (a+b)/2
        tmpJ2=(XiEQ_SurfSample(jSample-1)+XiEQ_SurfSample(jSample))/2. ! (a+b)/2
        DO q=0,NGeo
          DO p=0,NGeo
            XiOut(1)=tmp1*Xi_NGeo(p)+tmpI2
            XiOut(2)=tmp1*Xi_NGeo(q)+tmpJ2
            CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
                                                    ,Gradient=gradXiEta3D)
            ! calculate first fundamental form
            E=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
            F=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
            G=DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
            D=SQRT(E*G-F*F)
            area = area+tmp1*tmp1*D*wGP_NGeo(p)*wGP_NGeo(q)
          END DO
        END DO
        SurfMesh%SurfaceArea(iSample,jSample,SurfSideID) = area
      END DO ! iSample=1,nSurfSample
    END DO ! jSample=1,nSurfSample
  END IF
END DO ! iSide=1,nComputeNodeTotalSides

! get the full area of the BC's of all faces on local proc
Area=0.
DO iSide=1,nComputeNodeTotalSides
  ! skip if element is not proc
  IF (ElemInfo_Shared(ELEM_RANK,SideInfo_Shared(SIDE_ELEMID,iSide)).NE.myRank) CYCLE
  
  SurfSideID = SurfMesh%SideIDToSurfID(iSide)
  IF(SurfSideID.EQ.-1) CYCLE
  Area = Area + SUM(SurfMesh%SurfaceArea(:,:,SurfSideID))
END DO ! iSide=1,nComputeNodeTotalSides

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Area,1,MPI_DOUBLE_PRECISION,MPI_SUM,SurfCOMM%COMM,iError)
#endif /*USE_MPI*/

! clean-up
DEALLOCATE(Xi_NGeo,wGP_NGeo)

SWRITE(UNIT_stdOut,'(A,E25.14E3)') ' Surface-Area: ', Area
SWRITE(UNIT_stdOut,'(A)') ' ... DONE.'

END SUBROUTINE InitParticleBoundarySampling

#if USE_MPI
! TODO: rebuild for new halo region
SUBROUTINE InitSurfCommunicator()
!===================================================================================================================================
! Creates two new subcommunicators.
! SurfCOMM%COMM contains all MPI-Ranks which have sampling boundary faces in their halo-region and process which have sampling
! boundary faces in their origin region. This communicator is used to communicate the wall-sampled values of halo-faces to the
! origin face
! SurfCOMM%OutputCOMM is another subset. This communicator contains only the processes with origin surfaces. It is used to perform
! collective writes of the surf-sampled values.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars   ,ONLY:SurfMesh
USE MOD_Particle_Boundary_Vars   ,ONLY:SurfCOMM
USE MOD_Particle_MPI_Vars        ,ONLY:PartMPI
USE MOD_Particle_Boundary_Vars   ,ONLY:OffSetSurfSideMPI,OffSetSurfSide
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: color,iProc
INTEGER                   :: noSurfrank,Surfrank
LOGICAL                   :: hasSurf
INTEGER,ALLOCATABLE       :: countSurfSideMPI(:)
LOGICAL                   :: OutputOnProc
!===================================================================================================================================

! define color for procs with sampling surfaces
color=MPI_UNDEFINED
IF(SurfMesh%SurfOnProc) color=1001

! create ranks for sampling communicator
! TODO: shouldn't this just be an MPI_EXSCAN? But then, root wouldn't know the total number ...
IF(PartMPI%MPIRoot) THEN
  SurfRank        = -1
  noSurfrank      = -1
  SurfCOMM%MyRank = 0
  IF(SurfMesh%SurfOnProc) THEN
    SurfRank   = 0
  ELSE
    noSurfRank = 0
  END IF
  ! root communicates with all procs to get global numbers and informs the proc
  DO iProc=1,nProcessors-1
    CALL MPI_RECV(hasSurf,1,MPI_LOGICAL,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
    IF(hasSurf) THEN
      SurfRank   = SurfRank   + 1
      CALL MPI_SEND(SurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    ELSE
      noSurfRank = noSurfRank + 1
      CALL MPI_SEND(noSurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    END IF
  END DO
ELSE
  ! all other procs send their information to root
  CALL MPI_SEND(SurfMesh%SurfOnProc,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,iError)
  CALL MPI_RECV(SurfCOMM%MyRank    ,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,MPIstatus,iError)
END IF

! create new SurfMesh communicator for SurfMesh communication
CALL MPI_COMM_SPLIT(PartMPI%COMM,color,SurfCOMM%MyRank,SurfCOMM%COMM,iError)
IF(SurfMesh%SurfOnPRoc) THEN
  CALL MPI_COMM_SIZE(SurfCOMM%COMM,SurfCOMM%nProcs,iError)
ELSE
  SurfCOMM%nProcs = 0
END IF

! root flags itself as SurfCOMM%MPIROOT
SurfCOMM%MPIRoot = .FALSE.
IF(SurfCOMM%MyRank.EQ.0 .AND. SurfMesh%SurfOnProc) THEN
  SurfCOMM%MPIRoot=.TRUE.
END IF

! now, create output communicator when proc has local sides
OutputOnProc = .FALSE.
color        = MPI_UNDEFINED
IF (SurfMesh%nSides.GT.0) THEN
  OutputOnProc = .TRUE.
  color        = 1002
END IF

! create ranks for output communicator
! TODO: shouldn't this just be an MPI_EXSCAN? But then, root wouldn't know the total number ...
IF(PartMPI%MPIRoot) THEN
  SurfRank              = -1
  noSurfRank            = -1
  SurfCOMM%MyOutputRank = 0
  IF (SurfMesh%nSides.GT.0) THEN
    SurfRank   = 0
  ELSE
    noSurfRank = 0
  END IF
  ! root communicates with all procs to get global numbers and informs the proc
  DO iProc=1,nProcessors-1
    CALL MPI_RECV(hasSurf,1,MPI_LOGICAL,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
    IF(hasSurf) THEN
      SurfRank   = SurfRank+1
      CALL MPI_SEND(SurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    ELSE
      noSurfRank = noSurfRank+1
      CALL MPI_SEND(noSurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    END IF
  END DO
ELSE
  ! all other procs send their information to root
  CALL MPI_SEND(OutputOnProc,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,iError)
  CALL MPI_RECV(SurfCOMM%MyOutputRank,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,MPIstatus,iError)
END IF

! create new SurfMesh Output-communicator
CALL MPI_COMM_SPLIT(PartMPI%COMM, color, SurfCOMM%MyOutputRank, SurfCOMM%OutputCOMM,iError)
IF(OutputOnPRoc)THEN
  CALL MPI_COMM_SIZE(SurfCOMM%OutputCOMM, SurfCOMM%nOutputProcs,iError)
ELSE
  SurfCOMM%nOutputProcs = 0
END IF

! root flags itself as SurfCOMM%MPIOutputRoot
SurfCOMM%MPIOutputRoot  = .FALSE.
IF(SurfCOMM%MyOutputRank.EQ.0 .AND. OutputOnProc) THEN
  SurfCOMM%MPIOutputRoot = .TRUE.
END IF

IF(SurfMesh%nSides.EQ.0) RETURN
! get correct offsets for output of hdf5 file (master sides, i.e. proc local sides)
ALLOCATE(offsetSurfSideMPI(0:SurfCOMM%nOutputProcs)   &
        ,countSurfSideMPI (0:SurfCOMM%nOutputProcs-1))
offsetSurfSideMPI = 0
countSurfSideMPI  = 0

! SurfCOMM%MPIOutputRoot gathers all information and sends it to the proc
CALL MPI_GATHER(SurfMesh%nSides,1,MPI_INTEGER,countSurfSideMPI,1,MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)
IF (SurfCOMM%MPIOutputRoot) THEN
  DO iProc=1,SurfCOMM%nOutputProcs-1
    offsetSurfSideMPI(iProc) = SUM(countSurfSideMPI(0:iProc-1))
  END DO
  offsetSurfSideMPI(SurfCOMM%nOutputProcs) = SUM(countSurfSideMPI(:))
END IF
CALL MPI_BCAST (offsetSurfSideMPI,size(offsetSurfSideMPI),MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)
offsetSurfSide=offsetSurfSideMPI(SurfCOMM%MyOutputRank)

END SUBROUTINE InitSurfCommunicator

! TODO: rebuild for new halo region
!SUBROUTINE InitSurfCommunication()
!!===================================================================================================================================
!
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Particle_Globals
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------!
!! INPUT/OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                           :: iProc,GlobalProcID
!INTEGER                           :: firstSide,lastSide
!INTEGER                           :: SurfToGlobal(0:SurfCOMM%nProcs-1)
!!===================================================================================================================================
!
!! get mapping from local rank to global rank...
!CALL MPI_ALLGATHER(MyRank,1, MPI_INTEGER, SurfToGlobal(0:SurfCOMM%nProcs-1), 1, MPI_INTEGER, SurfCOMM%COMM, IERROR)
!
!! find all procs with surface sampling sides to communicate with
!isMPINeighbor=.FALSE.
!IF (SurfMesh%nTotalSides.GT.SurfMesh%nSides) THEN
!  DO iProc=0,SurfCOMM%nProcs-1
!    ! ignore myself
!    IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!    
!    ! check if the proc has surface sampling sides
!    GlobalProcID = SurfToGlobal(iProc)
!    firstSide    = ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElemMPI(GlobalProcID)+1)
!    lastSide     = ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElemMPI(GlobalProcID+1))
!    
!    DO iSide = firstSide,lastSide
!      SurfSideID = SurfMesh%SideIDToSurfID(iSide)
!      
!      ! ignore sides that are not sampled
!      IF(SurfSideID.EQ.-1) CYCLE
!      
!      ! flag once a surface sampling side is found
!      IF (.NOT.isMPINeighbor(iProc)) THEN
!        isMPINeighbor(iProc) = .TRUE.
!        EXIT
!      END IF
!    END DO ! iSide = firstSide,lastSide
!  END DO ! iProc=0,SurfCOMM%nProcs-1
!END IF
!
!! Make sure SurfCOMM%MPINeighbor is consistent
!ALLOCATE(RECV_STATUS_LIST(1:MPI_STATUS_SIZE,0:SurfCOMM%nProcs-1))
!! receive and send connection information
!DO iProc=0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  CALL MPI_IRECV( RecvMPINeighbor(iProc)                    &
!                , 1                                         &
!                , MPI_LOGICAL                               &
!                , iProc                                     &
!                , 1001                                      &
!                , SurfCOMM%COMM                             &
!                , RecvRequest(iProc)                        &
!                , IERROR )
!END DO ! iProc
!DO iProc=0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  CALL MPI_ISEND( isMPINeighbor(iProc)                      &
!                , 1                                         &
!                , MPI_LOGICAL                               &
!                , iProc                                     &
!                , 1001                                      &
!                , SurfCOMM%COMM                             &
!                , SendRequest(iProc)                        &
!                , IERROR )
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
!END DO ! iProc
!
!! finish communication
!DO iProc=0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
!  CALL MPI_WAIT(RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__,' MPI Communication error', IERROR)
!END DO ! iProc
!DEALLOCATE(RECV_STATUS_LIST)
!
!! Make sure we agree on how many sides need to be exchanged
!!>> TODO: if the halo region is reduced, we might encounter the case where some procs on the other nodes are NOT within our halo
!!>> region. Things will badly break then
!ALLOCATE(SurfExchange%nSidesSend(1:SurfCOMM%nMPINeighbors) &
!        ,SurfExchange%nSidesRecv(1:SurfCOMM%nMPINeighbors) &
!        ,SurfExchange%SendRequest(SurfCOMM%nMPINeighbors)  &
!        ,SurfExchange%RecvRequest(SurfCOMM%nMPINeighbors))
!
!SurfExchange%nSidesSend=0
!SurfExchange%nSidesRecv=0
!! loop over all neighbors
!DO iProc=1,SurfCOMM%nMPINeighbors
!  !
!END DO
!  
!
!END SUBROUTINE InitSurfCommunication


!SUBROUTINE GetHaloSurfMapping()
!!===================================================================================================================================
!! build all missing stuff for surface-sampling communicator, like
!! offSetMPI
!! MPI-neighbor list
!! PartHaloSideToProc
!! only receiving process knows to which local side the sending information is going, the sending process does not know the final
!! sideid
!! if only a processes has to send his data to another, the whole structure is build, but not filled. only if the nsidesrecv,send>0
!! communication is performed
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Particle_Globals
!USE MOD_Preproc
!USE MOD_Mesh_Vars                   ,ONLY:nSides,nBCSides
!USE MOD_Particle_Boundary_Vars      ,ONLY:SurfMesh,SurfComm,nSurfSample
!USE MOD_Particle_MPI_Vars           ,ONLY:PartHaloSideToProc,PartHaloElemToProc,SurfSendBuf,SurfRecvBuf,SurfExchange
!USE MOD_Particle_Mesh_Vars          ,ONLY:nTotalSides,PartSideToElem,PartElemToSide
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!LOGICAL                           :: isMPINeighbor(0:SurfCOMM%nProcs-1)
!LOGICAL                           :: RecvMPINeighbor(0:SurfCOMM%nProcs-1)
!INTEGER                           :: nDOF,ALLOCSTAT,SideID,TargetProc,TargetElem
!INTEGER                           :: iProc, GlobalProcID,iSide,ElemID,SurfSideID,LocalProcID,iSendSide,iRecvSide,iPos
!INTEGER,ALLOCATABLE               :: recv_status_list(:,:)
!INTEGER                           :: RecvRequest(0:SurfCOMM%nProcs-1),SendRequest(0:SurfCOMM%nProcs-1)
!INTEGER                           :: SurfToGlobal(0:SurfCOMM%nProcs-1)
!INTEGER                           :: NativeElemID, NativeLocSideID
!!===================================================================================================================================
!
!nDOF=(nSurfSample)*(nSurfSample)
!
!! get mapping from local rank to global rank...
!CALL MPI_ALLGATHER(MyRank,1, MPI_INTEGER, SurfToGlobal(0:SurfCOMM%nProcs-1), 1, MPI_INTEGER, SurfCOMM%COMM, IERROR)
!
!! get list of mpi surf-comm neighbors in halo-region
!isMPINeighbor=.FALSE.
!IF(SurfMesh%nTotalSides.GT.SurfMesh%nSides)THEN
!  ! PartHaloSideToProc has the mapping from the local sideid to the corresponding
!  ! element-id, localsideid, native-proc-id and local-proc-id
!  ! caution:
!  ! native-proc-id is id in global list || PartMPI%COMM
!  ! local-proc-id is the nth neighbor in SurfCOMM%COMM
!  ! caution2:
!  !  this mapping is done only for reflective bcs, thus, each open side or non-reflective side
!  !  points to -1
!  ! get all MPI-neighbors to communicate with
!  DO iProc=0,SurfCOMM%nProcs-1
!    IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!    GlobalProcID=SurfToGlobal(iProc)
!    DO iSide=nSides+1,nTotalSides
!      SurfSideID=SurfMesh%SideIDToSurfID(iSide)
!      IF(SurfSideID.EQ.-1) CYCLE
!      ! get elemid
!      ! TODO: Move this to new halo region
!      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
!      IF(ElemID.LE.PP_nElems)THEN
!        CALL abort(&
!__STAMP__&
!,' Error in PartSideToElem. Halo-Side cannot be connected to local element', ElemID  )
!      END IF
!      IF(GlobalProcID.EQ.TargetProc)THEN
!        IF(.NOT.isMPINeighbor(iProc))THEN
!          isMPINeighbor(iProc)=.TRUE.
!        END IF
!      END IF
!    END DO ! iSide=nSides+1,nTotalSides
!  END DO ! iProc = 0, SurfCOMM%nProcs-1
!END IF
!
!! Make sure SurfCOMM%MPINeighbor is consistent
!! 1) communication of found is neighbor
!ALLOCATE(RECV_STATUS_LIST(1:MPI_STATUS_SIZE,0:SurfCOMM%nProcs-1))
!! receive and send connection information
!DO iProc=0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  CALL MPI_IRECV( RecvMPINeighbor(iProc)                    &
!                , 1                                         &
!                , MPI_LOGICAL                               &
!                , iProc                                     &
!                , 1001                                      &
!                , SurfCOMM%COMM                             &
!                , RecvRequest(iProc)                        &
!                , IERROR )
!END DO ! iProc
!DO iProc=0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  CALL MPI_ISEND( isMPINeighbor(iProc)                      &
!                , 1                                         &
!                , MPI_LOGICAL                               &
!                , iProc                                     &
!                , 1001                                      &
!                , SurfCOMM%COMM                             &
!                , SendRequest(iProc)                        &
!                , IERROR )
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!    __STAMP__&
!    ,' MPI Communication error', IERROR)
!END DO ! iProc
!
!! finish communication
!DO iProc=0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!    __STAMP__&
!    ,' MPI Communication error', IERROR)
!  CALL MPI_WAIT(RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!    __STAMP__&
!    ,' MPI Communication error', IERROR)
!END DO ! iProc
!DEALLOCATE(RECV_STATUS_LIST)
!
!! 2) finalize MPI consistency check.
!DO iProc=0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  IF (RecvMPINeighbor(iProc).NEQV.isMPINeighbor(iProc)) THEN
!    IF(.NOT.isMPINeighbor(iProc))THEN
!      isMPINeighbor(iProc)=.TRUE.
!      ! it is a debug output
!      !IPWRITE(UNIT_stdOut,*) ' Found missing mpi-neighbor'
!    END IF
!  END IF
!END DO ! iProc
!
!! build the mapping
!SurfCOMM%nMPINeighbors=0
!IF(ANY(IsMPINeighbor))THEN
!  ! PartHaloSideToProc has the mapping from the local sideid to the corresponding
!  ! element-id, localsideid, native-proc-id and local-proc-id
!  ! caution:
!  ! native-proc-id is id in global list || PartMPI%COMM
!  ! local-proc-id is the nth neighbor in SurfCOMM%COMM
!  ! caution2:
!  !  this mapping is done only for reflective bcs, thus, each open side or non-reflective side
!  !  points to -1
!  IF(SurfMesh%nTotalSides.GT.SurfMesh%nSides)THEN
!    ALLOCATE(PartHaloSideToProc(1:4,nSides+1:nTotalSides))
!    PartHaloSideToProc=-1
!    ! get all MPI-neighbors to communicate with
!    DO iProc=0,SurfCOMM%nProcs-1
!      IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!      IF(isMPINeighbor(iProc))THEN
!        SurfCOMM%nMPINeighbors=SurfCOMM%nMPINeighbors+1
!      ELSE
!        CYCLE
!      END IF
!      GlobalProcID=SurfToGlobal(iProc)
!      DO iSide=nSides+1,nTotalSides
!        SurfSideID=SurfMesh%SideIDToSurfID(iSide)
!        IF(SurfSideID.EQ.-1) CYCLE
!        ! get elemid
!        ! TODO: Move this to new halo region
!        ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
!        IF(ElemID.LE.PP_nElems)THEN
!        CALL abort(&
!__STAMP__&
!,' Error in PartSideToElem. Halo-Side cannot be connected to local element', ElemID  )
!        END IF
!        IF(GlobalProcID.EQ.TargetProc)THEN
!          ! caution:
!          ! native-proc-id is id in global list || PartMPI%COMM
!          ! local-proc-id is the nth neighbor in SurfCOMM%COMM
!          PartHaloSideToProc(NATIVE_PROC_ID,iSide)=GlobalProcID
!          PartHaloSideToProc(LOCAL_PROC_ID ,iSide)=SurfCOMM%nMPINeighbors
!        END IF
!      END DO ! iSide=nSides+1,nTotalSides
!    END DO ! iProc = 0, SurfCOMM%nProcs-1
!  ELSE
!    ! process receives only surface data from other processes, but does not send data.
!    DO iProc=0,SurfCOMM%nProcs-1
!      IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!      IF(isMPINeighbor(iProc)) SurfCOMM%nMPINeighbors=SurfCOMM%nMPINeighbors+1
!    END DO ! iProc = 0, SurfCOMM%nProcs-1
!  END IF
!END IF
!
!! build SurfMesh exchange information
!! next, allocate SurfCOMM%MPINeighbor
!ALLOCATE(SurfCOMM%MPINeighbor(1:SurfCOMM%nMPINeighbors))
!! set native proc-id of each SurfCOMM-MPI-Neighbor
!LocalProcID=0
!DO iProc = 0,SurfCOMM%nProcs-1
!  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
!  IF(isMPINeighbor(iProc))THEN
!    LocalProcID=LocalProcID+1
!    ! map from local proc id to global
!    SurfCOMM%MPINeighbor(LocalProcID)%NativeProcID=iProc !PartMPI%MPINeighbor(iProc)
!  END IF
!END DO ! iProc=1,PartMPI%nMPINeighbors
!
!! array how many data has to be communicated
!! number of Sides
!ALLOCATE(SurfExchange%nSidesSend(1:SurfCOMM%nMPINeighbors) &
!        ,SurfExchange%nSidesRecv(1:SurfCOMM%nMPINeighbors) &
!        ,SurfExchange%SendRequest(SurfCOMM%nMPINeighbors)  &
!        ,SurfExchange%RecvRequest(SurfCOMM%nMPINeighbors)  )
!
!SurfExchange%nSidesSend=0
!SurfExchange%nSidesRecv=0
!! loop over all neighbors
!DO iProc=1,SurfCOMM%nMPINeighbors
!  ! proc-id in SurfCOMM%nProcs
!  DO iSide=nSides+1,nTotalSides
!    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
!    IF(SurfSideID.EQ.-1) CYCLE
!    ! get elemid
!    ! TODO: Move this to new halo region
!    IF(iProc.EQ.PartHaloSideToProc(LOCAL_PROC_ID,iSide))THEN
!      SurfExchange%nSidesSend(iProc)=SurfExchange%nSidesSend(iProc)+1
!      PartHaloSideToProc(LOCAL_SEND_ID,iSide) =SurfExchange%nSidesSend(iProc)
!    END IF
!  END DO ! iSide=nSides+1,nTotalSides
!END DO ! iProc=1,SurfCOMM%nMPINeighbors
!
!! open receive number of send particles
!ALLOCATE(RECV_STATUS_LIST(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors))
!DO iProc=1,SurfCOMM%nMPINeighbors
!  CALL MPI_IRECV( SurfExchange%nSidesRecv(iProc)            &
!                , 1                                         &
!                , MPI_INTEGER                               &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID  &
!                , 1001                                      &
!                , SurfCOMM%COMM                             &
!                , SurfExchange%RecvRequest(iProc)           &
!                , IERROR )
!END DO ! iProc
!
!DO iProc=1,SurfCOMM%nMPINeighbors
!  ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SendList(SurfExchange%nSidesSend(iProc)))
!  SurfCOMM%MPINeighbor(iProc)%SendList=0
!  CALL MPI_ISEND( SurfExchange%nSidesSend(iProc)           &
!                , 1                                        &
!                , MPI_INTEGER                              &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
!                , 1001                                     &
!                , SurfCOMM%COMM                            &
!                , SurfExchange%SendRequest(iProc)          &
!                , IERROR )
!END DO ! iProc
!
!
!! 4) Finish Received number of particles
!DO iProc=1,SurfCOMM%nMPINeighbors
!  CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!,' MPI Communication error', IERROR)
!  CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
!  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!END DO ! iProc
!
!! allocate send and receive buffer
!ALLOCATE(SurfSendBuf(SurfCOMM%nMPINeighbors))
!ALLOCATE(SurfRecvBuf(SurfCOMM%nMPINeighbors))
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).GT.0)THEN
!    ALLOCATE(SurfSendBuf(iProc)%content(2*SurfExchange%nSidesSend(iProc)),STAT=ALLOCSTAT)
!    SurfSendBuf(iProc)%Content=-1
!  END IF
!  IF(SurfExchange%nSidesRecv(iProc).GT.0)THEN
!    ALLOCATE(SurfRecvBuf(iProc)%content(2*SurfExchange%nSidesRecv(iProc)),STAT=ALLOCSTAT)
!  END IF
!END DO ! iProc=1,PartMPI%nMPINeighbors
!
!! open receive buffer
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
!  CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
!                , 2*SurfExchange%nSidesRecv(iProc)             &
!                , MPI_DOUBLE_PRECISION                         &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
!                , 1004                                         &
!                , SurfCOMM%COMM                                &
!                , SurfExchange%RecvRequest(iProc)              &
!                , IERROR )
!END DO ! iProc
!
!! build message
!! after this message, the receiving process knows to which of his sides the sending process will send the
!  ! surface data
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
!  iSendSide=0
!  iPos=1
!  DO iSide=nSides+1,nTotalSides
!    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
!    IF(SurfSideID.EQ.-1) CYCLE
!    ! get elemid
!    IF(iProc.EQ.PartHaloSideToProc(LOCAL_PROC_ID,iSide))THEN
!      iSendSide=iSendSide+1
!      ! get elemid
!      ! TODO: Move this to new halo region
!      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
!!       IF(ElemID.LE.PP_nElems)THEN
!!         IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
!!       END IF
!!       IF(ElemID.LE.1)THEN
!!         IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
!!       END IF
!      SurfCOMM%MPINeighbor(iProc)%SendList(iSendSide)=SurfSideID
!!       IPWRITE(*,*) 'negative elem id',PartHaloElemToProc(NATIVE_ELEM_ID,ElemID),PartSideToElem(S2E_LOC_SIDE_ID,iSide)
!      SurfSendBuf(iProc)%content(iPos  )= REAL(PartHaloElemToProc(NATIVE_ELEM_ID,ElemID))
!      SurfSendBuf(iProc)%content(iPos+1)= REAL(PartSideToElem(S2E_LOC_SIDE_ID,iSide))
!      iPos=iPos+2
!    END IF
!  END DO ! iSide=nSides+1,nTotalSides
!  IF(iSendSide.NE.SurfExchange%nSidesSend(iProc)) CALL abort(&
!__STAMP__&
!          ,' Message too short!',iProc)
!  IF(ANY(SurfSendBuf(iProc)%content.LE.0))THEN
!    IPWRITE(UNIT_stdOut,*) ' nSendSides', SurfExchange%nSidesSend(iProc), ' to Proc ', iProc
!    CALL abort(&
!__STAMP__&
!          ,' Sent NATIVE_ELEM_ID or LOCSIDEID is zero!')
!  END IF
!END DO
!
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
!  CALL MPI_ISEND( SurfSendBuf(iProc)%content               &
!                , 2*SurfExchange%nSidesSend(iProc)         &
!                , MPI_DOUBLE_PRECISION                     &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
!                , 1004                                     &
!                , SurfCOMM%COMM                            &
!                , SurfExchange%SendRequest(iProc)          &
!                , IERROR )
!END DO ! iProc
!
!! 4) Finish Received number of particles
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
!    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
!    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!END DO ! iProc
!
!! fill list with received side ids
!! store the receiving data
!DO iProc=1,SurfCOMM%nMPINeighbors
!  ALLOCATE(SurfCOMM%MPINeighbor(iProc)%RecvList(SurfExchange%nSidesRecv(iProc)))
!  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
!  iPos=1
!  DO iRecvSide=1,SurfExchange%nSidesRecv(iProc)
!    NativeElemID   =INT(SurfRecvBuf(iProc)%content(iPos))
!    NativeLocSideID=INT(SurfRecvBuf(iProc)%content(iPos+1))
!    IF(NativeElemID.GT.PP_nElems)THEN
!     CALL abort(&
!__STAMP__&
!          ,' Cannot send halo-data to other procs. big error! ', ElemID, REAL(PP_nElems))
!    END IF
!    SideID=PartElemToSide(E2S_SIDE_ID,NativeLocSideID,NativeElemID)
!    IF(SideID.GT.nSides)THEN
!      IPWRITE(UNIT_stdOut,*) ' Received wrong sideid. Is not a BC side! '
!      IPWRITE(UNIT_stdOut,*) ' SideID, nBCSides, nSides ', SideID, nBCSides, nSides
!      IPWRITE(UNIT_stdOut,*) ' ElemID, locsideid        ', NativeElemID, NativeLocSideID
!      IPWRITE(UNIT_stdOut,*) ' Sending process has error in halo-region! ', SurfCOMM%MPINeighbor(iProc)%NativeProcID
!     CALL abort(&
!__STAMP__&
!          ,' Big error in halo region! NativeLocSideID ', NativeLocSideID )
!    END IF
!    SurfSideID=SurfMesh%SideIDToSurfID(SideID)
!    IF(SurfSideID.EQ.-1)THEN
!     CALL abort(&
!__STAMP__&
!          ,' Side is not even a reflective BC side! ', SideID )
!    END IF
!    SurfCOMM%MPINeighbor(iProc)%RecvList(iRecvSide)=SideID
!    iPos=iPos+2
!  END DO ! RecvSide=1,SurfExchange%nSidesRecv(iProc)-1,2
!END DO ! iProc
!
!DO iProc=1,SurfCOMM%nMPINeighbors
!  SDEALLOCATE(SurfSendBuf(iProc)%content)
!  SDEALLOCATE(SurfRecvBuf(iProc)%content)
!  IF(SurfExchange%nSidesSend(iProc).GT.0) THEN
!    ALLOCATE(SurfSendBuf(iProc)%content((SurfMesh%SampSize)*nDOF*SurfExchange%nSidesSend(iProc)))
!    SurfSendBuf(iProc)%content=0.
!  END IF
!  IF(SurfExchange%nSidesRecv(iProc).GT.0) THEN
!    ALLOCATE(SurfRecvBuf(iProc)%content((SurfMesh%SampSize)*nDOF*SurfExchange%nSidesRecv(iProc)))
!    SurfRecvBuf(iProc)%content=0.
!  END IF
!END DO ! iProc
!DEALLOCATE(recv_status_list)
!
!CALL MPI_BARRIER(SurfCOMM%Comm,iError)
!
!END SUBROUTINE GetHaloSurfMapping


SUBROUTINE ExchangeSurfData()
!===================================================================================================================================
! exchange the surface data
! only processes with sampling sides in their halo region and the original process participate on the communication
! structure is similar to particle communication
! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
! the receiving process adds the new data to his own sides
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY:nSpecies
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfMesh,SurfComm,nSurfSample,SampWall
USE MOD_Particle_MPI_Vars           ,ONLY:SurfSendBuf,SurfRecvBuf,SurfExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize,nValues,iSurfSide,SurfSideID
INTEGER                         :: iPos,p,q,iProc,iSpec,nShift
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
REAL,ALLOCATABLE                :: SampWallTmp(:)
!===================================================================================================================================

nValues = SurfMesh%SampSize*nSurfSample**2
! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesRecv(iProc)*nValues
  CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                , MessageSize                                  &
                , MPI_DOUBLE_PRECISION                         &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1009                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%RecvRequest(iProc)              &
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  SurfSendBuf(iProc)%content = 0.
  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SurfSendBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)= SampWall(SurfSideID)%State(:,p,q)
        iPos=iPos+SurfMesh%SampSize
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
!   This is the State in the Halo region. We can safely delete it
    SampWall(SurfSideID)%State(:,:,:)=0.
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  CALL MPI_ISEND( SurfSendBuf(iProc)%content               &
                , MessageSize                              &
                , MPI_DOUBLE_PRECISION                     &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
                , 1009                                     &
                , SurfCOMM%COMM                            &
                , SurfExchange%SendRequest(iProc)          &
                , IERROR )
END DO ! iProc

! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
END DO ! iProc

! add data do my list
ALLOCATE(SampWallTmp(1:SurfMesh%SampSize))

DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
!       Treat each variable individual, depending on if they can be added
        SampWallTmp(:)                    = SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)

!---- All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
!---- 1. - .. / Impact Counter
        SampWall(SurfSideID)%State(1,p,q) = SampWall(SurfSideID)%State(1,p,q) + SampWallTmp(1)
!===================================================================================================================================
!---- 2. - 6. / Kinetic energy on impact (mean, min, max, M2, variance)
        SampWall(SurfSideID)%State(2,p,q) = (SampWall(SurfSideID)%State(2,p,q)*(SampWall(SurfSideID)%State(1,p,q) - SampWallTmp(1))&
                                            + SampWallTmp(2)*SampWallTmp(1)) / SampWall(SurfSideID)%State(1,p,q)
        IF (SampWallTmp(3).LT.SampWall(SurfSideID)%State(3,p,q)) THEN
            SampWall(SurfSideID)%State(3,p,q) = SampWallTmp(3)
        END IF
        IF (SampWallTmp(4).GT.SampWall(SurfSideID)%State(4,p,q)) THEN
            SampWall(SurfSideID)%State(4,p,q) = SampWallTmp(4)
        END IF
! 5-6 are M2 and variance. Nothing to be done about those as we need data after each particle impact
!===================================================================================================================================
!---- 7. - 11 / Impact angle (mean, min, max, M2, variance)
        SampWall(SurfSideID)%State(7,p,q) = (SampWall(SurfSideID)%State(7,p,q)*(SampWall(SurfSideID)%State(1,p,q) - SampWallTmp(1))&
                                            + SampWallTmp(7)*SampWallTmp(1)) / SampWall(SurfSideID)%State(1,p,q)
        IF (SampWallTmp(8).LT.SampWall(SurfSideID)%State(8,p,q)) THEN
            SampWall(SurfSideID)%State(8,p,q) = SampWallTmp(8)
        END IF
        IF (SampWallTmp(9).GT.SampWall(SurfSideID)%State(9,p,q)) THEN
            SampWall(SurfSideID)%State(9,p,q) = SampWallTmp(9)
        END IF
! 10-11 are M2 and variance. Nothing to be done about those as we need data after each particle impact
!===================================================================================================================================
!---- 12 - 14 / Sampling Current Forces at walls - we can simply add those
        SampWall(SurfSideID)%State(12:14,p,q) = SampWall(SurfSideID)%State(12:14,p,q) + SampWallTmp(12:14)
!===================================================================================================================================
!---- 15 - 17 / Sampling Average Forces at walls  - we can simply add those
        SampWall(SurfSideID)%State(15:17,p,q) = SampWall(SurfSideID)%State(15:17,p,q) + SampWallTmp(15:17)
!===================================================================================================================================
!<<< Repeat for specific species >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF (nSpecies.GT.1) THEN
            DO iSpec=1,nSpecies
                nShift = iSpec * nErosionVars
!---- 1. - .. / Impact Counter
                SampWall(SurfSideID)%State(1+nShift,p,q) = SampWall(SurfSideID)%State(1+nShift,p,q) + SampWallTmp(1+nShift)
!---- 2. - 6. / Kinetic energy on impact (mean, min, max, M2, variance)
                SampWall(SurfSideID)%State(2+nShift,p,q) = (SampWall(SurfSideID)%State(2+nShift,p,q)                               &
                                                           *(SampWall(SurfSideID)%State(1+nShift,p,q) - SampWallTmp(1+nShift))     &
                                                           + SampWallTmp(2+nShift)*SampWallTmp(1+nShift))                          &
                                                           / SampWall(SurfSideID)%State(1+nShift,p,q)
                IF (SampWallTmp(3+nShift).LT.SampWall(SurfSideID)%State(3+nShift,p,q)) THEN
                    SampWall(SurfSideID)%State(3+nShift,p,q) = SampWallTmp(3+nShift)
                END IF
                IF (SampWallTmp(4+nShift).GT.SampWall(SurfSideID)%State(4+nShift,p,q)) THEN
                    SampWall(SurfSideID)%State(4+nShift,p,q) = SampWallTmp(4+nShift)
                END IF
! 5-6 are M2 and variance. Nothing to be done about those as we need data after each particle impact
!===================================================================================================================================
!---- 7. - 11 / Impact angle (mean, min, max, M2, variance)
                SampWall(SurfSideID)%State(7+nShift,p,q) = (SampWall(SurfSideID)%State(7+nShift,p,q)                               &
                                                           *(SampWall(SurfSideID)%State(1+nShift,p,q) - SampWallTmp(1+nShift))     &
                                                           + SampWallTmp(7+nShift)*SampWallTmp(1+nShift))                          &
                                                           / SampWall(SurfSideID)%State(1+nShift,p,q)
                IF (SampWallTmp(8+nShift).LT.SampWall(SurfSideID)%State(8+nShift,p,q)) THEN
                    SampWall(SurfSideID)%State(8+nShift,p,q) = SampWallTmp(8+nShift)
                END IF
                IF (SampWallTmp(9+nShift).GT.SampWall(SurfSideID)%State(9+nShift,p,q)) THEN
                    SampWall(SurfSideID)%State(9+nShift,p,q) = SampWallTmp(9+nShift)
                END IF
! 10-11 are M2 and variance. Nothing to be done about those as we need data after each particle impact
!===================================================================================================================================
!---- 12 - 14 / Sampling Current Forces at walls - we can simply add those
                SampWall(SurfSideID)%State(12+nShift:14+nShift,p,q) = SampWall(SurfSideID)%State(12+nShift:14+nShift,p,q)          &
                                                                      + SampWallTmp(12+nShift:14+nShift)
!===================================================================================================================================
!---- 15 - 17 / Sampling Average Forces at walls  - we can simply add those
                SampWall(SurfSideID)%State(15+nShift:17+nShift,p,q) = SampWall(SurfSideID)%State(15+nShift:17+nShift,p,q)          &
                                                                      + SampWallTmp(15+nShift:17+nShift)
!===================================================================================================================================
            END DO
        END IF
!        SampWall(SurfSideID)%State(:,p,q) = SampWall(SurfSideID)%State(:,p,q) &
!                                           + SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)
        iPos=iPos+SurfMesh%SampSize
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  SurfRecvBuf(iProc)%content = 0.
END DO ! iProc

SDEALLOCATE(SampWallTmp)

END SUBROUTINE ExchangeSurfData
#endif /*MPI*/


SUBROUTINE RestartParticleBoundarySampling(remap_opt)
!===================================================================================================================================
!> Restart the particle boundary sampling (assuming RestartData array exists)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Input
USE MOD_HDF5_Output
USE MOD_Output_Vars,                ONLY:ProjectName
USE MOD_Restart_Vars,               ONLY:RestartTime
USE MOD_Particle_Vars,              ONLY:nSpecies,WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars,     ONLY:SurfCOMM
USE MOD_Particle_Boundary_Vars,     ONLY:nSurfSample,SurfMesh,offSetSurfSide,SampWall
USE MOD_Particle_Erosion_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: remap_opt       ! routine was called from posti. Change input file name
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString
CHARACTER(LEN=255)                  :: H5_Name
INTEGER                             :: iVal,iSide,p,q!,NRestart
LOGICAL                             :: ErosionDataExists
INTEGER                             :: RestartVarCount,NRestartFile
REAL,ALLOCATABLE                    :: RestartArray(:,:,:,:)
!#if USE_MPI
!INTEGER                             :: iProc, GlobalProcID,ElemID,SurfSideID,LocalProcID
!#endif
!===================================================================================================================================

#if USE_MPI
! Ignore procs without surfMesh on them
IF(SurfMesh%nSides.EQ.0) RETURN
! Make sure the remaining procs are synchronized. DO NOT USE THIS, BLOCKING!
! CALL MPI_BARRIER(SurfCOMM%COMM,iERROR)
IF(.NOT.WriteMacroSurfaceValues) THEN
    ErosionRestart = .FALSE.
    RETURN
END IF
#endif /*MPI*/

IF(SurfCOMM%MPIOutputRoot)THEN
!  WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' RESTARTING EROSION TRACKING FROM HDF5 FILE...'
END IF

! Get number of required erosion vars
IF (nSpecies.EQ.1) THEN
    RestartVarCount = nErosionVars
ELSE
    RestartVarCount = (nErosionVars*(nSpecies+1))
END IF

! Allocate array for restart on ALL procs
IF (nSpecies.EQ.1) THEN
  ALLOCATE(RestartArray(nErosionVars,1:nSurfSample,1:nSurfSample,SurfMesh%nGlobalSides))
ELSE
  ALLOCATE(RestartArray((nErosionVars)*(nSpecies+1),1:nSurfSample,1:nSurfSample,SurfMesh%nGlobalSides))
END IF
RestartArray = 0.

! Open restart array in general erosion file
IF (PRESENT(remap_opt)) THEN
    FileString=remap_opt
ELSE
    FileName=TIMESTAMP(TRIM(ProjectName)//'_ErosionSurfState',RestartTime)
    FileString=TRIM(FileName)//'.h5'
END IF
! Apparently we can't read it only as root and then communicate it afterwards. Read the entire array on every proc and pick later.
!#if USE_MPI
!IF(SurfCOMM%MPIOutputRoot)THEN
!#endif

  IF(.NOT.FILEEXISTS(FileString)) THEN
      IF(SurfCOMM%MPIOutputRoot)THEN
        WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Erosion File does not exist for current time. Aborting erosion restart ...'
        WRITE(UNIT_StdOut,'(132("-"))')
      END IF
      RETURN
  END IF

  ! We got the file, now open it. Be careful to open it ONLY on SurfComm procs!
#if USE_MPI
! Make sure the remaining procs are synchronized. DO NOT USE THIS, BLOCKING!
!  CALL MPI_BARRIER(SurfCOMM%OutputCOMM,iERROR)
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
#else
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*MPI*/

  ! Check if erosion file has restart data
  CALL DatasetExists(File_ID,'RestartData',ErosionDataExists)

  IF (ErosionDataExists) THEN
      ! Check if size of restart array is correct
      CALL ReadAttribute(File_ID,'NRestart',1,IntScalar=NRestartFile)
      IF (NRestartFile.NE.RestartVarCount) THEN
          WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Number of variables in erosion file does not match. Aborting erosion restart ...'
          WRITE(UNIT_StdOut,'(132("-"))')
          RETURN
      END IF

      ! Read array from restart file
      WRITE(H5_Name,'(A)') 'RestartData'
      CALL ReadArray(ArrayName=H5_Name, rank=4,&
                     nVal=      (/RestartVarCount  ,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                     offset_in  = 0,&
                     offset_dim = 1,&
                     RealArray  = RestartArray)

      ! We got the array, close the file
      CALL CloseDataFile()

      ! Only loop over sides on current proc (meaningless in single core case)
      DO iSide=1,SurfMesh%nSides
        ! Not a surfMesh side
!        IF (SurfSideID.EQ.-1) CYCLE

        ! Write array back to SampWall
        IF (nSpecies.EQ.1) THEN
          DO p=1,nSurfSample
            DO q=1,nSurfSample
              DO iVal=1,nErosionVars
                SampWall(iSide)%State(iVal,p,q) = RestartArray(iVal,p,q,offsetSurfSide+iSide)
              END DO !iVal
            END DO !q
          END DO !p
        ELSE
          DO p=1,nSurfSample
            DO q=1,nSurfSample
              DO iVal=1,(nErosionVars*(nSpecies+1))
                SampWall(iSide)%State(iVal,p,q) = RestartArray(iVal,p,q,offsetSurfSide+iSide)
              END DO !iVal
            END DO !q
          END DO !p
        END IF !nSpecies
      END DO

      ErosionRestart = .TRUE.

      IF(SurfCOMM%MPIOutputRoot)THEN
        WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' EROSION TRACKING RESTART SUCCESSFUL.'
        WRITE(UNIT_StdOut,'(132("-"))')
      END IF

      DEALLOCATE(RestartArray)
  ELSE
    IF(SurfCOMM%MPIOutputRoot)THEN
      ErosionRestart = .FALSE.

      WRITE(UNIT_stdOut,*)'RestartData does not exists in erosion tracking file.'
      WRITE(UNIT_StdOut,'(132("-"))')
    END IF
  END IF

END SUBROUTINE RestartParticleBoundarySampling


SUBROUTINE WriteSurfSampleToHDF5(MeshFileName,OutputTime,remap_opt)
!===================================================================================================================================
!> write the final values of the surface sampling to a HDF5 state file
!> additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Output
USE MOD_Output_Vars,                ONLY:ProjectName
USE MOD_Particle_Vars,              ONLY:nSpecies
USE MOD_Particle_Boundary_Vars,     ONLY:SurfCOMM,nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars,     ONLY:nSurfSample,SurfMesh,offSetSurfSide,SampWall
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_HDF5_output,       ONLY:WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: remap_opt       !routine was called from posti. Change output file name
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString,Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255)                  :: SpecID
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: nVar2D, nVar2D_Spec, nVar2D_Total, nVarCount, iSpec, nShift
INTEGER                             :: iSide,iVal,p,q,RestartVarCount,iSurfSide,nShiftRHS
REAL,ALLOCATABLE                    :: RestartArray(:,:,:,:)
REAL                                :: startT,endT
!===================================================================================================================================

#if USE_MPI
CALL MPI_BARRIER(SurfCOMM%COMM,iERROR)
IF(SurfMesh%nSides.EQ.0) RETURN
#endif /*MPI*/
IF(SurfCOMM%MPIOutputRoot)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE EROSION SURFACE STATE TO HDF5 FILE...'
  GETTIME(startT)
END IF

! Allocate array for restart on each proc
IF (nSpecies.EQ.1) THEN
    ALLOCATE(RestartArray(nErosionVars,1:nSurfSample,1:nSurfSample,SurfMesh%nGlobalSides))
ELSE
    ALLOCATE(RestartArray((nErosionVars)*(nSpecies+1),1:nSurfSample,1:nSurfSample,SurfMesh%nGlobalSides))
END IF

IF (nSpecies.EQ.1) THEN
    DO iSide=1,SurfMesh%nSides
        DO p=1,nSurfSample
            DO q=1,nSurfSample
                DO iVal=1,nErosionVars
                    RestartArray(iVal,p,q,iSide) = SampWall(iSide)%State(iVal,p,q)
                END DO !iVal
            END DO !q
        END DO !p
    END DO !iSide
ELSE
    DO iSide=1,SurfMesh%nSides
        DO p=1,nSurfSample
            DO q=1,nSurfSample
                DO iVal=1,(nErosionVars*(nSpecies+1))
                    RestartArray(iVal,p,q,iSide) = SampWall(iSide)%State(iVal,p,q)
                END DO !iVal
            END DO !q
        END DO !p
    END DO !iSide
END IF

IF (PRESENT(remap_opt)) THEN
    FileName=TIMESTAMP(TRIM(ProjectName)//remap_opt,OutputTime)
ELSE
    FileName=TIMESTAMP(TRIM(ProjectName)//'_ErosionSurfState',OutputTime)
END IF
FileString=TRIM(FileName)//'.h5'

! Create dataset attribute "SurfVarNames"
IF (nSpecies.EQ.1) THEN
    nVar2D = nErosionVars - 1
ELSE
    nVar2D = (nErosionVars - 1) * (nSpecies+1)
END IF
nVar2D_Spec=1
nVar2D_Total = nVar2D + nVar2D_Spec*nSpecies

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF(SurfCOMM%MPIOutputRoot)THEN
#endif
  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'DSMCSurfState'

  ! Write file header
  CALL WriteHDF5Header(Statedummy,File_ID)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
  CALL WriteAttributeToHDF5(File_ID,'BC_Surf',nSurfBC,StrArray=SurfBCName)
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=nSurfSample)
  NodeTypeTemp='VISU'
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeTypeTemp/))

  ALLOCATE(Str2DVarNames(1:nVar2D_Total))
  nVarCount=0

  DO iSpec=1,nSpecies
    WRITE(SpecID,'(I3.3)') iSpec
    Str2DVarNames(nVarCount+1) ='Spec'//TRIM(SpecID)//'_Counter'
    nVarCount=nVarCount+nVar2D_Spec
  END DO ! iSpec=1,nSpecies

  ! fill varnames for total values
  Str2DVarNames(nVarCount+1) ='Impacts'
  Str2DVarNames(nVarCount+2) ='ImpactsPerArea'
  Str2DVarNames(nVarCount+3) ='E_kin(mean)'
  Str2DVarNames(nVarCount+4) ='E_kin(min)'
  Str2DVarNames(nVarCount+5) ='E_kin(max)'
  Str2DVarNames(nVarCount+6) ='E_kin(variance)'
  Str2DVarNames(nVarCount+7) ='ImpactAngle(mean)'
  Str2DVarNames(nVarCount+8) ='ImpactAngle(min)'
  Str2DVarNames(nVarCount+9) ='ImpactAngle(max)'
  Str2DVarNames(nVarCount+10)='ImpactAngle(variance)'
  Str2DVarNames(nVarCount+11)='CurrentForceX'
  Str2DVarNames(nVarCount+12)='CurrentForceY'
  Str2DVarNames(nVarCount+13)='CurrentForceZ'
  Str2DVarNames(nVarCount+14)='AverageForceX'
  Str2DVarNames(nVarCount+15)='AverageForceY'
  Str2DVarNames(nVarCount+16)='AverageForceZ'

  IF (nSpecies.GT.1) THEN
      DO iSpec=1,nSpecies
          WRITE(SpecID,'(I3.3)') iSpec
          nShift = iSpec * (nErosionVars - 1)

          Str2DVarNames(nVarCount+nShift+1) ='Species'//TRIM(SpecID)//'_Impacts'
          Str2DVarNames(nVarCount+nShift+2) ='Species'//TRIM(SpecID)//'_ImpactsPerArea'
          Str2DVarNames(nVarCount+nShift+3) ='Species'//TRIM(SpecID)//'_E_kin(mean)'
          Str2DVarNames(nVarCount+nShift+4) ='Species'//TRIM(SpecID)//'_E_kin(min)'
          Str2DVarNames(nVarCount+nShift+5) ='Species'//TRIM(SpecID)//'_E_kin(max)'
          Str2DVarNames(nVarCount+nShift+6) ='Species'//TRIM(SpecID)//'_E_kin(variance)'
          Str2DVarNames(nVarCount+nShift+7) ='Species'//TRIM(SpecID)//'_ImpactAngle(mean)'
          Str2DVarNames(nVarCount+nShift+8) ='Species'//TRIM(SpecID)//'_ImpactAngle(min)'
          Str2DVarNames(nVarCount+nShift+9) ='Species'//TRIM(SpecID)//'_ImpactAngle(max)'
          Str2DVarNames(nVarCount+nShift+10)='Species'//TRIM(SpecID)//'_ImpactAngle(variance)'
          Str2DVarNames(nVarCount+nShift+11)='Species'//TRIM(SpecID)//'_CurrentForceX'
          Str2DVarNames(nVarCount+nShift+12)='Species'//TRIM(SpecID)//'_CurrentForceY'
          Str2DVarNames(nVarCount+nShift+13)='Species'//TRIM(SpecID)//'_CurrentForceZ'
          Str2DVarNames(nVarCount+nShift+14)='Species'//TRIM(SpecID)//'_AverageForceX'
          Str2DVarNames(nVarCount+nShift+15)='Species'//TRIM(SpecID)//'_AverageForceY'
          Str2DVarNames(nVarCount+nShift+16)='Species'//TRIM(SpecID)//'_AverageForceZ'
      END DO
  END IF

  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVar2D_Total,StrArray=Str2DVarNames)

  CALL CloseDataFile()
  DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF

CALL MPI_BARRIER(SurfCOMM%OutputCOMM,iERROR)
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*MPI*/

nVarCount=0
! This array is to be visualized
WRITE(H5_Name,'(A)') 'SurfaceData'
DO iSpec = 1,nSpecies
    CALL WriteArrayToHDF5(DataSetName=H5_Name, rank=4,&
                    nValGlobal=(/nVar2D_Total,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                    nVal=      (/nVar2D_Spec ,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                    offset=    (/nVarCount   ,          0,          0,offsetSurfSide/),&
                    collective=.TRUE.,  RealArray=MacroSurfaceSpecVal(:,:,:,:,iSpec))
    nVarCount = nVarCount + nVar2D_Spec
END DO
CALL WriteArrayToHDF5(DataSetName=H5_Name, rank=4,&
                    nValGlobal=(/NVar2D_Total,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                    nVal=      (/nVar2D      ,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                    offset=    (/nVarCount   ,          0,          0,offsetSurfSide/),&
                    collective=.TRUE., RealArray=MacroSurfaceVal)

CALL CloseDataFile()

IF (nSpecies.EQ.1) THEN
    RestartVarCount = nErosionVars
ELSE
    RestartVarCount = (nErosionVars*(nSpecies+1))
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF(SurfCOMM%MPIOutputRoot)THEN
#endif
  CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'RestartState'

  ! Write file header
  !CALL WriteHDF5Header(Statedummy,File_ID)

  CALL WriteAttributeToHDF5(File_ID,'NRestart',1,IntegerScalar=RestartVarCount)

  ! We do not need more attributes as we can use the ones from above

  CALL CloseDataFile()
#if USE_MPI
END IF

CALL MPI_BARRIER(SurfCOMM%OutputCOMM,iERROR)
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*MPI*/

! This array is purely for restart
WRITE(H5_Name,'(A)') 'RestartData'
CALL WriteArrayToHDF5(DataSetName=H5_Name, rank=4,&
                    nValGlobal=(/RestartVarCount  ,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                    nVal=      (/RestartVarCount  ,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                    offset=    (/0                ,          0,          0,offsetSurfSide/),&
                    collective=.TRUE., RealArray=RestartArray)

CALL CloseDataFile()

!> Only reset current forces here so we have them in case of restart
iSpec = 1
!---- Only one species. Only total values necessary
!===================================================================================================================================
DO iSurfSide=1,SurfMesh%nSides
  DO q=1,nSurfSample
    DO p=1,nSurfSample
        SampWall(iSurfSide)%State(12,p,q) = 0.
        SampWall(iSurfSide)%State(13,p,q) = 0.
        SampWall(iSurfSide)%State(14,p,q) = 0.
    END DO
  END DO
END DO

!---- Multiple species. All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
IF (nSpecies.GT.1) THEN
  DO iSurfSide=1,SurfMesh%nSides
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        DO iSpec=1,nSpecies
          nShift    = iSpec * (nErosionVars-1)
          nShiftRHS = iSpec * nErosionVars
          SampWall(iSurfSide)%State(12+nShiftRHS,p,q) = 0.
          SampWall(iSurfSide)%State(13+nShiftRHS,p,q) = 0.
          SampWall(iSurfSide)%State(14+nShiftRHS,p,q) = 0.
        END DO
      END DO
    END DO
  END DO
END IF

IF(SurfCOMM%MPIOutputRoot)THEN
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' DONE  [',EndT-StartT,'s]'
END IF

END SUBROUTINE WriteSurfSampleToHDF5


SUBROUTINE FinalizeParticleBoundarySampling()
!===================================================================================================================================
! deallocate everything
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Boundary_Vars
#if USE_MPI
USE MOD_Particle_MPI_Vars           ,ONLY:SurfSendBuf,SurfRecvBuf,SurfExchange,PartHaloSideToProc
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iSurfSide
#if USE_MPI
INTEGER :: iProc
#endif
!===================================================================================================================================

SDEALLOCATE(XiEQ_SurfSample)
SDEALLOCATE(SurfMesh%SurfaceArea)
SDEALLOCATE(SurfMesh%SideIDToSurfID)
SDEALLOCATE(SurfMesh%SurfSideToGlobSideMap)
DO iSurfSide=1,SurfMesh%nSides
  SDEALLOCATE(SampWall(iSurfSide)%State)
END DO
SDEALLOCATE(SurfBCName)
SDEALLOCATE(SampWall)
#if USE_MPI
SDEALLOCATE(PartHaloSideToProc)
SDEALLOCATE(SurfExchange%nSidesSend)
SDEALLOCATE(SurfExchange%nSidesRecv)
SDEALLOCATE(SurfExchange%SendRequest)
SDEALLOCATE(SurfExchange%RecvRequest)
DO iProc=1,SurfCOMM%nMPINeighbors
  IF (ALLOCATED(SurfSendBuf))THEN
    SDEALLOCATE(SurfSendBuf(iProc)%content)
  END IF
  IF (ALLOCATED(SurfRecvBuf))THEN
    SDEALLOCATE(SurfRecvBuf(iProc)%content)
  END IF
  IF (ALLOCATED(SurfCOMM%MPINeighbor))THEN
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%SendList)
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%RecvList)
  END IF
END DO ! iProc=1,PartMPI%nMPINeighbors
SDEALLOCATE(SurfCOMM%MPINeighbor)
SDEALLOCATE(SurfSendBuf)
SDEALLOCATE(SurfRecvBuf)
SDEALLOCATE(OffSetSurfSideMPI)
#endif /*USE_MPI*/
END SUBROUTINE FinalizeParticleBoundarySampling

END MODULE MOD_Particle_Boundary_Sampling
