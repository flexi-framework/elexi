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
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
MODULE MOD_Particle_MPI_Halo
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

#if USE_MPI
INTERFACE IdentifyHaloMPINeighborhood
  MODULE PROCEDURE IdentifyHaloMPINeighborhood
END INTERFACE

INTERFACE ExchangeHaloGeometry
  MODULE PROCEDURE ExchangeHaloGeometry
END INTERFACE

INTERFACE WriteParticlePartitionInformation
  MODULE PROCEDURE WriteParticlePartitionInformation
END INTERFACE


PUBLIC :: IdentifyHaloMPINeighborhood,ExchangeHaloGeometry
PUBLIC :: WriteParticlePartitionInformation!,WriteParticleMappingPartitionInformation

!===================================================================================================================================

CONTAINS


SUBROUTINE IdentifyHaloMPINeighborhood(iProc,SideIndex,ElemIndex)
!===================================================================================================================================
! Searches for sides in the neighborhood of MPI-neighbors
! mark all Sides for communication which are in range of halo_eps of own MPI_Sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,         ONLY:GEO,SidePeriodicType,nPartSides,PartElemToSide
USE MOD_Particle_MPI_Vars,          ONLY:PartMPI
USE MOD_Particle_Surfaces_Vars,     ONLY:BezierControlPoints3D
USE MOD_Particle_MPI_Vars,          ONLY:halo_eps
USE MOD_Mesh_Vars,                  ONLY:NGeo,firstMPISide_MINE,BC,BoundaryType
USE MOD_Particle_Mesh_Vars,         ONLY:MortarSlave2MasterInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)   :: iProc  ! MPI proc with which the local proc has to exchange boundary information
INTEGER, INTENT(INOUT):: SideIndex(nPartSides)
INTEGER, INTENT(INOUT):: ElemIndex(PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iElem,ilocSide,SideID,iSide
INTEGER                     :: DataSize
! MPI exchange
TYPE tMPISideMessage
  REAL(KIND=8), ALLOCATABLE :: BezierSides3D(:,:,:,:) ! BezierSides of boundary faces
  INTEGER(KIND=4)           :: nMPISides              ! number of Sides to send
  REAL                      :: MinMaxXYZ(1:6)
END TYPE
TYPE(tMPISideMessage)       :: SendMsg
TYPE(tMPISideMessage)       :: RecvMsg
INTEGER                     :: ALLOCSTAT,PVID,iDir
REAL                        :: MinMax(2),Vec1(1:3)
LOGICAL                     :: SideisDone(1:nPartSides)
LOGICAL                     :: SideInside
!=================================================================================================================================

! 1) Exchange Sides:
!    Each proc receives all sides that lie in FIBGM cells within an  eps-distance of FIBGM cells containing
!    any of my own MPI sides.
! 2) Go through each FIBGM cell and if there are MPI-Neighbors, search the neighbor  surrounding FIBGM
!     cells for my own MPI sides.

! Step1: find Sides of myProc that are within halo_eps distance to iProc

! here only MPI Sides.... why not INNER Sides, too????????/
! PO: correction: we send only the MPI sides, because only the MPI-Faces has to be checked, not the interior faces for
!     distance calulation. if MPI-Side is in range, then all the other sides are in range, too.
!     caution: maybe there can be a special case, in which the nBCSides are in HaloRange but not the MPI-Side
SideIndex(:)=0
SendMsg%nMPISides=0

! exchange of bounding box of process-local (MY)-region
! take the halo_eps into account
SendMsg%MinMaxXYZ(1)=GEO%xmin-halo_eps
SendMsg%MinMaxXYZ(2)=GEO%ymin-halo_eps
SendMsg%MinMaxXYZ(3)=GEO%zmin-halo_eps
SendMsg%MinMaxXYZ(4)=GEO%xmax+halo_eps
SendMsg%MinMaxXYZ(5)=GEO%ymax+halo_eps
SendMsg%MinMaxXYZ(6)=GEO%zmax+halo_eps

IF (PartMPI%MyRank.LT.iProc) THEN
  CALL MPI_SEND(SendMsg%MinMaxXYZ,6,MPI_DOUBLE_PRECISION,iProc,1001,PartMPI%COMM,IERROR)
  CALL MPI_RECV(RecvMsg%MinMaxXYZ,6,MPI_DOUBLE_PRECISION,iProc,1002,PartMPI%COMM,MPISTATUS,IERROR)
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  CALL MPI_RECV(RecvMsg%MinMaxXYZ,6,MPI_DOUBLE_PRECISION,iProc,1001,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_SEND(SendMsg%MinMaxXYZ,6,MPI_DOUBLE_PRECISION,iProc,1002,PartMPI%COMM,IERROR)
END IF

SideisDone=.FALSE.
SideisDone(1:firstMPISide_MINE-1)=.TRUE.
DO iElem=1,PP_nElems
  DO ilocSide=1,6
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
    ! if not SideID exists, ignore
    IF(SideID.LT.1) CYCLE
    IF(SideID.GT.nPartSides) CYCLE ! only MY sides are checked, no HALO sides of other processes
    IF(SideID.LT.firstMPISide_MINE)THEN ! for my-inner-sides, we have to check the periodic sides
      IF(BC(SideID).LE.0) CYCLE ! no boundary side, cycle
      IF(BoundaryType(BC(SideID),BC_TYPE).NE.1) CYCLE  ! not a periodic BC, cycle
    END IF
    ! side is already checked
    IF(SideIsDone(SideID)) CYCLE
    ! mortar side
    IF(MortarSlave2MasterInfo(SideID).NE.-1) CYCLE
    IF(SideIndex(SideID).EQ.0)THEN
      PVID=SidePeriodicType(SideID)
      Vec1=0.
      IF(PVID.NE.0) Vec1= SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
      ! check if an overlap is possible
      DO iDir=1,3 ! x,y,z
        ! check for SideInSide at the beginning of the loop
        SideInside=.FALSE.
        MinMax(1)=MINVAL(BezierControlPoints3D(iDir,0:NGeo,0:NGeo,SideID))
        MinMax(2)=MAXVAL(BezierControlPoints3D(iDir,0:NGeo,0:NGeo,SideID))
        IF(MinMax(1).LE.RecvMsg%MinMaxXYZ(iDir+3).AND.MinMax(2).GE.RecvMsg%MinMaxXYZ(iDir))THEN
          SideInside=.TRUE.
        END IF
        ! periodic sides
        IF(Vec1(iDir).NE.0)THEN
          MinMax(1)=MinMax(1)+Vec1(iDir)
          MinMax(2)=MinMax(2)+Vec1(iDir)
          IF(MinMax(1).LE.RecvMsg%MinMaxXYZ(iDir+3).AND.MinMax(2).GE.RecvMsg%MinMaxXYZ(iDir))THEN
            SideInside=.TRUE.
          END IF
        END IF
        IF(.NOT.SideInside) EXIT ! Side outside of required range, can be skipped
      END DO ! directions
      ! if Side is marked, it has to be communicated
      IF(SideInside)THEN
        SendMsg%nMPISides=SendMsg%nMPISides+1
        SideIndex(SideID)=SendMsg%nMPISides
      END IF
      SideisDone(SideID)=.TRUE.
    END IF
  END DO ! ilocSide
END DO ! iElem=1,PP_nElems

!--- NOTE: IF SENDMSG%NNODES IS 0 AT THIS POINT, THEN I SHOULD BE ABLE TO RETURN HERE!!!
!          This is not done yet because I'm not sure whether there are still inconsistencies in the code...

!--- Debugging information
LOGWRITE(*,*)' nMPISides for iProc=',iProc,':',SendMsg%nMPISides

!--- Send number of MPI sides to MPI neighbor iProc and receive number of MPI
!    sides from MPI neighbor iProc (immediate neighbor or not)
IF (PartMPI%MyRank.LT.iProc) THEN
  CALL MPI_SEND(SendMsg%nMPISides,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
  CALL MPI_RECV(RecvMsg%nMPISides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,MPISTATUS,IERROR)
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  CALL MPI_RECV(RecvMsg%nMPISides,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_SEND(SendMsg%nMPISides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,IERROR)
END IF

!--- Allocate Message
IF (SendMsg%nMPISides.GT.0) THEN
  ALLOCATE(SendMsg%BezierSides3D(1:3,0:NGeo,0:nGEO,1:SendMsg%nMPISides), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMessage%BezierSides3D ',SendMsg%nMPISides)
  END IF
  SendMsg%BezierSides3D=0.
END IF
IF (RecvMsg%nMPISides.GT.0) THEN
  ALLOCATE(RecvMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,1:RecvMsg%nMPISides), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMessage%BezierSides3D ',RecvMsg%nMPISides)
  END IF
  RecvMsg%BezierSides3D=0.
END IF

!--- Send any (corner-) nodes from the MPI-sides to the MPI-neighbor iProc
!    and receive iProc's (corner-) nodes in return
!--- fill send buffers
!--- Step 2: send myproc MPI-side-nodes to iProc

DO iSide=1,nPartSides
  IF(SideIndex(iSide).NE.0)THEN
     SendMsg%BezierSides3D(1:3,0:NGeo,0:NGeo,SideIndex(iSide))=BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)
  END IF
END DO ! iSide=1,nPartSides

!--- send and receive data
DataSize=3*(NGeo+1)*(NGeo+1)
IF(PartMPI%MyRank.LT.iProc)THEN
  IF (SendMsg%nMPISides.GT.0) CALL MPI_SEND(SendMsg%BezierSides3D       &
                                           ,SendMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1103                        &
                                           ,PartMPI%COMM                &
                                           ,IERROR                      )
  IF (RecvMsg%nMPISides.GT.0) CALL MPI_RECV(RecvMsg%BezierSides3D       &
                                           ,RecvMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1104                        &
                                           ,PartMPI%COMM                &
                                           ,MPISTATUS                   &
                                           ,IERROR                      )
ELSE IF(PartMPI%MyRank.GT.iProc)THEN
  IF (RecvMsg%nMPISides.GT.0) CALL MPI_RECV(RecvMsg%BezierSides3D       &
                                           ,RecvMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1103                        &
                                           ,PartMPI%COMM                &
                                           ,MPISTATUS                   &
                                           ,IERROR                      )
  IF (SendMsg%nMPISides.GT.0) CALL MPI_SEND(SendMsg%BezierSides3D       &
                                           ,SendMsg%nMPISides*DataSize  &
                                           ,MPI_DOUBLE_PRECISION        &
                                           ,iProc                       &
                                           ,1104                        &
                                           ,PartMPI%COMM                &
                                           ,IERROR                      )
END IF



! PO
! For each side, identifz the FIBGM cell(s) in which the side resides and
! search the surrounding nPaddingCells for neighboring elements
! For idiots: iProc and tells me which sides are on his MPI bound,
!             and I check which sides are within eps range

SideIndex(:)=0
ElemIndex(:)=0
IF (RecvMsg%nMPISides.GT.0) THEN
  CALL CheckMPINeighborhoodByFIBGM(RecvMsg%BezierSides3D,RecvMsg%nMPISides,SideIndex,ElemIndex)
END IF

!--- Deallocate Messages
IF (SendMsg%nMPISides.GT.0) THEN
  DEALLOCATE(SendMsg%BezierSides3D, STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not deallocate SendMessage%BezierSides3D proc ',iProc)
  END IF
END IF
IF (RecvMsg%nMPISides.GT.0) THEN
  DEALLOCATE(RecvMsg%BezierSides3D, STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(&
    __STAMP__&
    ,'Could not deallocate RecvMessage%BezierSides3D proc ',iProc)
  END IF
END IF

END SUBROUTINE IdentifyHaloMPINeighborhood


SUBROUTINE CheckMPINeighborhoodByFIBGM(BezierSides3D,nExternalSides,SideIndex,ElemIndex)
!===================================================================================================================================
! Compute distance of MPI-Side to my-Sides, if distance below helo_eps2 distance, mark size for MPI-Exchange
! Question: Why does one does not mark the INNER Sides, too??
!           Then, the halo region would only require the elements inside of itself and not only the MPI sides...??
!           Or is it sufficient?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,                 ONLY:NGeo,OffSetElem
USE MOD_Particle_Mesh_Vars,        ONLY:GEO, FIBGMCellPadding,NbrOfCases,casematrix,nPartSides,PartElemToSide,PartSideToElem &
                                       ,SidePeriodicType,ElemToElemGlob
USE MOD_Particle_MPI_Vars,         ONLY:halo_eps
USE MOD_Particle_Surfaces_Vars,    ONLY:BezierControlPoints3D
USE MOD_Particle_Mappings,         ONLY:SideToAdjointLocSide
USE MOD_Particle_Tracking_Vars,    ONLY:CartesianPeriodic
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,    INTENT(IN)      :: BezierSides3D(1:3,0:NGeo,0:NGeo,1:nExternalSides)
INTEGER, INTENT(IN)      :: nExternalSides
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)   :: SideIndex(nPartSides)
INTEGER, INTENT(INOUT)   :: ElemIndex(PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iSide, NbOfSides,p,q,ElemID,ilocSide,SideID,r,s,iBGMElem,iCase,NbOfElems,NBElemID,iMortar
REAL                     :: NodeX(1:3)!,xNodes(1:3,0:NGeo,0:NGeo)
INTEGER                  :: iBGM,jBGM,kBGM,iPBGM,jPBGM,kPBGM,PVID
LOGICAL                  :: leave
REAL                     :: Vec1(1:3),Vec2(1:3),Vec3(1:3)
REAL                     :: distance(1:3)
INTEGER                  :: iElem,firstBezierPoint,lastBezierPoint!,flip,AdjointLocSideID(2)
!===================================================================================================================================

! For each (NGeo+1)^2 BezierControlPoint of each side, the FIBGM cell(s) in which the side
! resides is identified and the surrouding nPaddingCells for each neighboring element
! are searched
!--- for idiots: get BezierControlPoints of myProc that are within eps distance to MPI-bound
!                of iProc
SideIndex(:)=0
ElemIndex=0
NbOfSides=0
NBOfElems=0

! check if the full mesh is required
Vec1(1) = GEO%xmaxglob-GEO%xminglob
Vec1(2) = GEO%ymaxglob-GEO%yminglob
Vec1(3) = GEO%zmaxglob-GEO%zminglob
! if full mesh is required, then mark all elements for communication
IF(SQRT(DOT_PRODUCT(Vec1,Vec1)).LE.halo_eps) THEN
  NbOfElems=PP_nElems
  DO iElem=1,PP_nElems
    ElemIndex(iElem)=iElem
  END DO ! iElem=1,PP_nElems
  RETURN
END IF

DO iSide=1,nExternalSides
  DO q=0,NGeo
    DO p=0,NGeo
      NodeX(:) = BezierSides3d(:,p,q,iSide)
      iBGM = INT((NodeX(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
      jBGM = INT((NodeX(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
      kBGM = INT((NodeX(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
      DO iPBGM = iBGM-FIBGMCellPadding(1),iBGM+FIBGMCellPadding(1)
        DO jPBGM = jBGM-FIBGMCellPadding(2),jBGM+FIBGMCellPadding(2)
          DO kPBGM = kBGM-FIBGMCellPadding(3),kBGM+FIBGMCellPadding(3)
            IF((iPBGM.GT.GEO%FIBGMimax).OR.(iPBGM.LT.GEO%FIBGMimin) .OR. &
               (jPBGM.GT.GEO%FIBGMjmax).OR.(jPBGM.LT.GEO%FIBGMjmin) .OR. &
               (kPBGM.GT.GEO%FIBGMkmax).OR.(kPBGM.LT.GEO%FIBGMkmin) ) CYCLE
            DO iBGMElem = 1, GEO%FIBGM(iPBGM,jPBGM,kPBGM)%nElem
              ElemID = GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element(iBGMElem)
              IF(ElemIndex(ElemID).GT.0) CYCLE ! element is already marked
              DO ilocSide=1,6
                SideID=PartElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
                IF(SideID.LT.1) CYCLE
                IF(SideID.GT.nPartSides) CYCLE
                ! caution, not save if corect
                leave=.FALSE.
                IF(SideIndex(SideID).EQ.0)THEN
                  DO s=0,NGeo
                    DO r=0,NGeo
                      IF(SQRT(DOT_Product(BezierControlPoints3D(:,r,s,SideID)-NodeX &
                                         ,BezierControlPoints3D(:,r,s,SideID)-NodeX )).LE.halo_eps)THEN
                        NbOfSides=NbOfSides+1
                        SideIndex(SideID)=NbOfSides
                        IF(ElemIndex(ElemID).EQ.0)THEN
                          NbOfElems=NbOfElems+1
                          ElemIndex(ElemID)=NbofElems
                        END IF
                        leave=.TRUE.
                        ! mark potential Inner elements on the other side
                        DO iMortar=1,4
                          NBElemID=INT(ElemToElemGlob(iMortar,ilocSide,offSetElem+ElemID)-offSetElem,4)
                          CHECKSAFEINT(NBElemID,4)
                          IF(NBElemID.LE.0) CYCLE
                          IF(NBElemID.GT.PP_nElems) CYCLE
                          ! check if NBElem is already marked, if not, mark it!
                          IF(ElemIndex(NbElemID).EQ.0)THEN
                            NbOfElems=NbOfElems+1
                            ElemIndex(NbElemID)=NbofElems
                          END IF
                        END DO ! iMortar=1,4
                        EXIT
                      END IF
                    END DO ! r
                    IF(leave) EXIT
                  END DO ! s
                END IF ! SideIndex(SideID).EQ.0
                IF(leave) EXIT
              END DO ! ilocSide
            END DO ! iElem
          END DO ! kPBGM
        END DO ! jPBGM
      END DO ! i PBGM
    END DO ! p
  END DO ! q
END DO ! iSide

!--- if there are periodic boundaries, they need to be taken into account as well:
IF (GEO%nPeriodicVectors.GT.0) THEN
  IF(.NOT.CartesianPeriodic)THEN
    DO iElem=1,PP_nElems
      IF(ElemIndex(iElem).GT.0)CYCLE
      leave=.FALSE.
      DO ilocSide=1,6
        SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
        IF(SideID.LT.1) CYCLE
        IF(SideID.GT.nPartSides) CYCLE
        PVID=SidePeriodicType(SideID)
        IF(PVID.EQ.0) CYCLE
        ! do not (!) reduce the number of side-BezierControlPoints which are compared to the received sides
        ! since the (moved) sides are checked individually
        firstBezierPoint=0
        lastBezierPoint=NGeo
        Vec1(1:3)= SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
        ! loop over all send sides, you do not have to check the halo mesh like in the svn-trunk
        ! the distance-x,y,z check is cheaper than a modified check here and similar fast
        DO iSide=1,nExternalSides
          DO q=firstBezierPoint,lastBezierPoint
            DO p=firstBezierPoint,lastBezierPoint
              NodeX(:) = BezierSides3D(:,p,q,iSide)
              DO s=firstBezierPoint,lastBezierPoint
                DO r=firstBezierPoint,lastBezierPoint
                  distance(1)=ABS(BezierControlPoints3D(1,r,s,SideID)+Vec1(1)-NodeX(1))
                  IF(distance(1).GT.halo_eps) CYCLE
                  distance(2)=ABS(BezierControlPoints3D(2,r,s,SideID)+Vec1(2)-NodeX(2))
                  IF(distance(2).GT.halo_eps) CYCLE
                  distance(3)=ABS(BezierControlPoints3D(3,r,s,SideID)+Vec1(3)-NodeX(3))
                  IF(distance(3).GT.halo_eps) CYCLE
                  IF(SQRT(DOT_Product(Distance,Distance)).LE.halo_eps)THEN
                    IF(SideIndex(SideID).GT.0)THEN
                      NbOfSides=NbOfSides+1
                      SideIndex(SideID)=NbOfSides
                    END IF
                    NbOfElems=NbOfElems+1
                    ElemIndex(iElem)=NbofElems
                    leave=.TRUE.
                    ! mark potential Inner elements on the other side
                    DO iMortar=1,4
                      NBElemID=INT(ElemToElemGlob(iMortar,ilocSide,offSetElem+iElem)-offSetElem,4)
                      CHECKSAFEINT(NBElemID,4)
                      IF(NBElemID.LE.0) CYCLE
                      IF(NBElemID.GT.PP_nElems) CYCLE
                      ! check if NBElem is already marked, if not, mark it!
                      IF(ElemIndex(NbElemID).EQ.0)THEN
                        NbOfElems=NbOfElems+1
                        ElemIndex(NbElemID)=NbofElems
                      END IF
                    END DO ! iMortar=1,4
                    EXIT
                  END IF ! Distance <halo_eps
                END DO ! r=firstBezierPoint,lastBezierPoint
                IF(leave) EXIT
              END DO ! s=firstBezierPoint,lastBezierPoint
              IF(leave) EXIT
            END DO ! p=firstBezierPoint,lastBezierPoint
            IF(leave) EXIT
          END DO ! q=firstBezierPoint,lastBezierPoint
          IF(leave) EXIT
        END DO ! iSide=1,nExternalSides
      END DO ! ilocSide=1,6
    END DO ! ! iElem=1,PP_nElems
  ELSE
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
    !--- check sides shifted by periodic vectors, add to SideIndex if match
    DO iSide=1,nExternalSides
      DO iCase = 1, NbrOfCases
        IF ((casematrix(iCase,1).EQ.0) .AND. &  ! DON'T DO THE UNMOVED PART, HAS BEEN DONE ABOVE
            (casematrix(iCase,2).EQ.0) .AND. &
            (casematrix(iCase,3).EQ.0)) CYCLE
        DO q=0,NGeo
          DO p=0,NGeo
            NodeX(:) = BezierSides3d(:,p,q,iSide) + &
                       casematrix(iCase,1)*Vec1(1:3) + &
                       casematrix(iCase,2)*Vec2(1:3) + &
                       casematrix(iCase,3)*Vec3(1:3)
            iBGM = INT((NodeX(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
            jBGM = INT((NodeX(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
            kBGM = INT((NodeX(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
            DO iPBGM = iBGM-FIBGMCellPadding(1),iBGM+FIBGMCellPadding(1)
              DO jPBGM = jBGM-FIBGMCellPadding(2),jBGM+FIBGMCellPadding(2)
                DO kPBGM = kBGM-FIBGMCellPadding(3),kBGM+FIBGMCellPadding(3)
                  IF((iPBGM.GT.GEO%FIBGMimax).OR.(iPBGM.LT.GEO%FIBGMimin) .OR. &
                     (jPBGM.GT.GEO%FIBGMjmax).OR.(jPBGM.LT.GEO%FIBGMjmin) .OR. &
                     (kPBGM.GT.GEO%FIBGMkmax).OR.(kPBGM.LT.GEO%FIBGMkmin) ) CYCLE
                  DO iBGMElem = 1, GEO%FIBGM(iPBGM,jPBGM,kPBGM)%nElem
                    ElemID = GEO%FIBGM(iPBGM,jPBGM,kPBGM)%Element(iBGMElem)
                    IF(ElemIndex(ElemID).GT.0) CYCLE ! element is already marked
                    DO ilocSide=1,6
                      SideID=PartElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
                      IF(SideID.LT.1) CYCLE
                      IF(SideID.GT.nPartSides) CYCLE
                      ! caution, not save if corect
                      leave=.FALSE.
                      IF(SideIndex(SideID).EQ.0)THEN
                        DO s=0,NGeo
                          DO r=0,NGeo
                            IF(SQRT(DOT_Product(BezierControlPoints3D(:,r,s,SideID)-NodeX &
                                               ,BezierControlPoints3D(:,r,s,SideID)-NodeX )).LE.halo_eps)THEN
                              NbOfSides=NbOfSides+1
                              SideIndex(SideID)=NbOfSides
                              IF(ElemIndex(ElemID).EQ.0)THEN
                                NbOfElems=NbOfElems+1
                                ElemIndex(ElemID)=NbofElems
                              END IF
                              leave=.TRUE.
                              EXIT
                            END IF
                          END DO ! r
                          IF(leave) EXIT
                        END DO ! s
                      END IF ! SideIndex(SideID).EQ.0
                      IF(leave) EXIT
                    END DO ! ilocSide
                  END DO ! iElem
                END DO ! kPBGM
              END DO ! jPBGM
            END DO ! i PBGM
          END DO ! p
        END DO ! q
      END DO ! iCase
    END DO ! iSide
  END IF
END IF  ! nperiodicvectors>0

! finally, all elements connected to this side have to be marked
! this is a sanity step and could be omitted
DO iSide=1,nPartSides
  IF(SideIndex(iSide).GT.0)THEN
    ! check both elements connected to side if they are marked
    ! master
    ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
    IF((ElemID.GT.0).AND.(ElemID.LE.PP_nElems))THEN
      IF(ElemIndex(ElemID).EQ.0)THEN
        NbOfElems=NbOfElems+1
        ElemIndex(ElemID)=NbofElems
      END IF
    END IF
    ! slave
    ElemID=PartSideToElem(S2E_NB_ELEM_ID,iSide)
    IF((ElemID.GT.0).AND.(ElemID.LE.PP_nElems))THEN
      IF(ElemIndex(ElemID).EQ.0)THEN
        NbOfElems=NbOfElems+1
        ElemIndex(ElemID)=NbofElems
      END IF
    END IF
  END IF
END DO ! iSide=1,nSides

END SUBROUTINE CheckMPINeighborhoodByFIBGM


SUBROUTINE ExchangeHaloGeometry(iProc,ElemList)
!===================================================================================================================================
! exchange of halo geometry
! including:
! BezierControlPoints3D
! ElemToSide
! BC-Type and State
! GEO%NodeCoords
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartHaloElemToProc
USE MOD_Mesh_Vars,              ONLY:nElems, nBCSides, BC,nGeo
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalSides,nTotalElems,SidePeriodicType,PartBCSideList,nPartSides,ElemBaryNGeo,CurvedElem
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartElemToElemGlob,nTotalBCSides,ElemType,ElemHasAuxBCs,nNodes
USE MOD_Particle_Mesh_Vars,     ONLY:XCL_NGeo,dXCL_NGeo,nTotalNodes
USE MOD_Mesh_Vars,              ONLY:MortarType
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,SideType,SideDistance,SideNormVec
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToElemGlob,GEO
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping,TriaTracking
USE MOD_Particle_Boundary_Vars, ONLY:nAuxBCs,UseAuxBCs
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iProc       ! MPI proc with which the local proc is to exchange boundary information
INTEGER, INTENT(INOUT)          :: ElemList(PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE tMPISideMessage
  REAL,ALLOCATABLE          :: BezierControlPoints3D(:,:,:,:)
  INTEGER,ALLOCATABLE       :: ElemToNodeID(:,:)
  INTEGER,ALLOCATABLE       :: ElemSideNodeID(:,:,:)
  REAL,ALLOCATABLE          :: NodeCoords(:,:)
  LOGICAL,ALLOCATABLE       :: ConcaveElemSide(:,:)
  REAL,ALLOCATABLE          :: XCL_NGeo (:,:,:,:,:)
  REAL,ALLOCATABLE          :: dXCL_NGeo(:,:,:,:,:,:)
  REAL,ALLOCATABLE,DIMENSION(:,:)     :: ElemBaryNGeo
  INTEGER,ALLOCATABLE       :: ElemToSide(:,:,:)
  INTEGER(KIND=8),ALLOCATABLE ::ElemToElemGlob(:,:,:)
  INTEGER,ALLOCATABLE       :: SideToElem(:,:)
  INTEGER,ALLOCATABLE       :: MortarType(:,:)
  INTEGER,ALLOCATABLE       :: SideBCType(:)
  INTEGER,ALLOCATABLE       :: BC(:)
  INTEGER,ALLOCATABLE       :: NativeElemID(:)
  REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: SideSlabNormals                 ! normal vectors of bounding slab box
  REAL,ALLOCATABLE,DIMENSION(:,:)    :: SideSlabIntervals               ! intervalls beta1, beta2, beta3
  LOGICAL,ALLOCATABLE,DIMENSION(:)   :: BoundingBoxIsEmpty
  !INTEGER,ALLOCATABLE       :: PeriodicElemSide(:,:)
  INTEGER                   :: nNodes                 ! number of nodes to send
  INTEGER                   :: nSides                 ! number of sides to send
  INTEGER                   :: nElems                 ! number of elems to send
  LOGICAL,ALLOCATABLE       :: curvedElem(:)
  INTEGER,ALLOCATABLE       :: ElemType(:)
  INTEGER,ALLOCATABLE       :: SideType(:)
  REAL   ,ALLOCATABLE       :: SideDistance(:)
  REAL   ,ALLOCATABLE       :: SideNormVec(:,:)
  LOGICAL,ALLOCATABLE       :: ElemHasAuxBCs(:,:)
END TYPE
TYPE(tMPISideMessage)       :: SendMsg
TYPE(tMPISideMessage)       :: RecvMsg
INTEGER                     :: ALLOCSTAT
INTEGER                     :: newSideID,haloSideID,newElemID
LOGICAL,ALLOCATABLE         :: isElem(:),isSide(:),isNode(:),isDone(:)
INTEGER, ALLOCATABLE        :: ElemIndex(:), SideIndex(:), NodeIndex(:), HaloInc(:)
INTEGER                     :: iElem, ilocSide,SideID,iSide,iIndex,iHaloSide,flip, iNode
INTEGER                     :: nDoubleSides,tmpnSides,tmpnElems,tmpnNodes
INTEGER                     :: datasize,datasize2,datasize3
INTEGER                     :: tmpbcsides
!===================================================================================================================================

ALLOCATE(isElem(1:nElems))
IF (.NOT.ALLOCATED(isElem)) CALL abort(&
  __STAMP__&
  ,'Could not allocate isElem')
isElem(:) = .FALSE.

ALLOCATE(isSide(1:nPartSides))
IF (.NOT.ALLOCATED(isSide)) CALL abort(&
  __STAMP__&
  ,'Could not allocate isSide')
isSide(:) = .FALSE.

ALLOCATE(isNode(1:nNodes))
IF (.NOT.ALLOCATED(isNode)) CALL abort(&
  __STAMP__&
  ,'Could not allocate isNode')
isNode(:) = .FALSE.

ALLOCATE(ElemIndex(1:nElems))
IF (.NOT.ALLOCATED(ElemIndex)) CALL abort(&
  __STAMP__&
  ,'Could not allocate ElemIndex')
ElemIndex(:) = 0

ALLOCATE(SideIndex(1:nPartSides))
IF (.NOT.ALLOCATED(SideIndex)) CALL abort(&
  __STAMP__&
  ,'Could not allocate SideIndex')
SideIndex(:) = 0

ALLOCATE(NodeIndex(1:nNodes))
IF (.NOT.ALLOCATED(NodeIndex)) CALL abort(&
  __STAMP__&
  ,'Could not allocate NodeIndex')
NodeIndex(:) = 0

!--- First, count marker node indices (nNeighborhoodNodes are within eps distance of at least one MPI-neighbor's node)
!--- For each MPI neighbor, identify the number of sides and elements to be sent
SendMsg%nElems=0
SendMsg%nSides=0
SendMsg%nNodes=0

! 1) get number of elements and sides
DO iElem=1,nElems
  IF(ElemList(iElem).NE.0)THEN
    !IF(.NOT.isElem(iElem)) THEN
    SendMsg%nElems=SendMsg%nElems+1
    ElemIndex(iElem) = SendMsg%nElems
    isElem(iElem)=.TRUE.
    !END IF ! NOT isElem
  END IF
END DO ! iElem

! 2) mark all required sides and get number of send sides
DO iElem=1,nElems
  IF(isElem(iElem))THEN
    DO ilocSide=1,6
      SideID=PartElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      IF(SideID.LT.1) CYCLE
      IF(SideID.GT.nPartSides) CYCLE
      IF(DoRefMapping)THEN
        IF(.NOT.isSide(SideID))THEN
          IF((SideID.LE.nBCSides).OR.(BC(SideID).NE.0))THEN
            ! missing: what do do with BC sides??"
            SendMsg%nSides=SendMsg%nSides+1
            SideIndex(SideID) = SendMsg%nSides
            isSide(SideID)=.TRUE.
          END IF
        END IF ! not isSide
      ELSE
        IF(.NOT.isSide(SideID)) THEN
          SendMsg%nSides=SendMsg%nSides+1
          SideIndex(SideID) = SendMsg%nSides
          isSide(SideID)=.TRUE.
        END IF ! not isSide
      END IF
    END DO ! ilocSide
    IF (TriaTracking) THEN
      !--- name nodes and update node-count
      DO iNode=1,8
        IF(.NOT.isNode(GEO%ElemToNodeID(iNode,iElem))) THEN
          SendMsg%nNodes=SendMsg%nNodes+1
          NodeIndex(GEO%ElemToNodeID(iNode,iElem)) = SendMsg%nNodes
          isNode(GEO%ElemToNodeID(iNode,iElem))=.TRUE.
        END IF
      END DO
    END IF
  END IF ! Element is marked to send
END DO ! iElem

!--- Communicate number of sides (trias,quads), elems (tets,hexas) and nodes to each MPI proc
IF (PartMPI%MyRank.LT.iProc) THEN
  CALL MPI_SEND(SendMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
  CALL MPI_SEND(SendMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,IERROR)
  CALL MPI_SEND(SendMsg%nNodes,1,MPI_INTEGER,iProc,1103,PartMPI%COMM,IERROR)
  CALL MPI_RECV(RecvMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_RECV(RecvMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_RECV(RecvMsg%nNodes,1,MPI_INTEGER,iProc,1103,PartMPI%COMM,MPISTATUS,IERROR)
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  CALL MPI_RECV(RecvMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_RECV(RecvMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_RECV(RecvMsg%nNodes,1,MPI_INTEGER,iProc,1103,PartMPI%COMM,MPISTATUS,IERROR)
  CALL MPI_SEND(SendMsg%nElems,1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
  CALL MPI_SEND(SendMsg%nSides,1,MPI_INTEGER,iProc,1102,PartMPI%COMM,IERROR)
  CALL MPI_SEND(SendMsg%nNodes,1,MPI_INTEGER,iProc,1103,PartMPI%COMM,IERROR)
END IF

! allocate send buffers for nodes, sides and elements for each MPI neighbor
! ElemToSide Mapping
IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%ElemToSide(1:2,1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate SendMsg%ElemToSide',SendMsg%nElems)
  SendMsg%ElemToSide(:,:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemToSide(1:2,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems)
  RecvMsg%ElemToSide(:,:,:)=0
END IF

! ElemToElemGlob Mapping
IF (SendMsg%nElems.GT.0) THEN       ! ElemToElem(1:4,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%ElemToElemGlob(1:4,1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not allocate SendMsg%ElemToElemGlob',SendMsg%nElems)
  SendMsg%ElemToElemGlob(:,:,:)=-1
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemToElemGlob(1:4,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%ElemToElemGlob',RecvMsg%nElems)
  RecvMsg%ElemToElemGlob(:,:,:)=-1
END IF
! BezierControlPoints3D for exchange
IF (SendMsg%nSides.GT.0) THEN       ! Beziercontrolpoints3d
  ALLOCATE(SendMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%BezierControlPoints3D',SendMsg%nSides)
  SendMsg%BezierControlPoints3D=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%BezierControlPoints3D',RecvMsg%nSides)
  RecvMsg%BezierControlPoints3D=0.
END IF
! Elem types
IF (SendMsg%nElems.GT.0) THEN
  ALLOCATE(SendMsg%CurvedElem(1:SendMsg%nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%CurvedElem',SendMsg%nElems)
  SendMsg%CurvedElem=.FALSE.
  IF (.NOT.DoRefMapping) THEN
    ALLOCATE(SendMsg%ElemType(1:SendMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%ElemType',SendMsg%nElems)
    SendMsg%ElemType=-1
  END IF
  IF (UseAuxBCs) THEN
    ALLOCATE(SendMsg%ElemHasAuxBCs(1:SendMsg%nElems , 1:nAuxBCs),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%ElemHasAuxBCs',SendMsg%nElems)
      SendMsg%ElemHasAuxBCs=.FALSE.
END IF
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%CurvedElem(1:RecvMsg%nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate ResvMsg%CurvedElem',RecvMsg%nElems)
  RecvMsg%CurvedElem=.FALSE.
  IF (.NOT.DoRefMapping) THEN
    ALLOCATE(RecvMsg%ElemType(1:RecvMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemType',RecvMsg%nElems)
    RecvMsg%ElemType=-1
  END IF
  IF (UseAuxBCs) THEN
    ALLOCATE(RecvMsg%ElemHasAuxBCs(1:RecvMsg%nElems , 1:nAuxBCs),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemHasAuxBCs',RecvMsg%nElems)
    RecvMsg%ElemHasAuxBCs=.FALSE.
END IF
END IF

IF (TriaTracking) THEN
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToNodeID
    ALLOCATE(SendMsg%ElemToNodeID(1:8,1:SendMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%ElemToNodeID',SendMsg%nElems)
    SendMsg%ElemToNodeID=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ElemToNodeID(1:8,1:RecvMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemToNodeID',RecvMsg%nElems)
    RecvMsg%ElemToNodeID=0
  END IF
  IF (SendMsg%nElems.GT.0) THEN       ! ElemSideNodeID
    ALLOCATE(SendMsg%ElemSideNodeID(1:4,1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%ElemSideNodeID',SendMsg%nElems)
    SendMsg%ElemSideNodeID=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ElemSideNodeID(1:4,1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemSideNodeID',RecvMsg%nElems)
    RecvMsg%ElemSideNodeID=0
  END IF
  IF (SendMsg%nNodes.GT.0) THEN       ! NodeCoords
    ALLOCATE(SendMsg%Nodecoords(1:3,1:SendMsg%nNodes),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%Nodecoords',SendMsg%nNodes)
    SendMsg%Nodecoords=0.
  END IF
  IF (RecvMsg%nNodes.GT.0) THEN
    ALLOCATE(RecvMsg%NodeCoords(1:3,1:RecvMsg%nNodes),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%NodeCoords',RecvMsg%nNodes)
    RecvMsg%NodeCoords=0.
  END IF
  IF (SendMsg%nElems.GT.0) THEN       ! ConcaveElemSide
    ALLOCATE(SendMsg%ConcaveElemSide(1:6,1:SendMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%ConcaveElemSide',SendMsg%nElems)
    SendMsg%ConcaveElemSide=.FALSE.
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%ConcaveElemSide(1:6,1:RecvMsg%nElems),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ConcaveElemSide',RecvMsg%nElems)
    RecvMsg%ConcaveElemSide=.FALSE.
  END IF
END IF

IF(DoRefMapping)THEN
  ! XCL_NGeo for exchange
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%XCL_NGeo',SendMsg%nElems)
    SendMsg%XCL_NGeo(:,:,:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%XCL_NGeo',RecvMsg%nElems)
    RecvMsg%XCL_NGeo(:,:,:,:,:)=0
  END IF

  ! DXCL_NGeo for exchange
  IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
    ALLOCATE(SendMsg%DXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:SendMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate SendMsg%DXCL_NGeo',SendMsg%nElems)
    SendMsg%DXCL_NGeo(:,:,:,:,:,:)=0
  END IF
  IF (RecvMsg%nElems.GT.0) THEN
    ALLOCATE(RecvMsg%DXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
    IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
      ,'Could not allocate RecvMsg%ElemToSide',RecvMsg%nElems)
    RecvMsg%DXCL_NGeo(:,:,:,:,:,:)=0
  END IF
END IF

! ElemBaryNGeo
IF (SendMsg%nElems.GT.0) THEN       ! ElemToSide(1:2,1:iLocSide,1:nElems)
  ALLOCATE(SendMsg%ElemBaryNGeo(1:3,1:SendMsg%nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
   __STAMP__&
    ,'Could not allocate SendMsg%ElemBaryNGeo',SendMsg%nElems)
  SendMsg%ElemBaryNGeo(:,:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%ElemBaryNGeo(1:3,1:RecvMsg%nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%ElemBary',RecvMsg%nElems)
  RecvMsg%ElemBaryNGeo(:,:)=0
END IF

! SideToElem Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SideToElem(1:2,1:nSides)
  ALLOCATE(SendMsg%SideToElem(1:5,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideToElem',SendMsg%nSides)
  SendMsg%SideToElem(:,:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideToElem(1:5,1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideToElem',RecvMsg%nSides)
  RecvMsg%SideToElem(:,:)=0
END IF
! side types
IF (SendMsg%nSides.GT.0) THEN
  ALLOCATE(SendMsg%SideType(1:SendMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideType',SendMsg%nSides)
  SendMsg%SideType(:)=-1
  ALLOCATE(SendMsg%SideDistance(1:SendMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideDistance',SendMsg%nSides)
  SendMsg%SideDistance(:)=0.
  ALLOCATE(SendMsg%SideNormVec(1:3,1:SendMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideNormVec',SendMsg%nSides)
  SendMsg%SideNormVec(:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideType(1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideType',RecvMsg%nSides)
  RecvMsg%SideType(:)=-1
  ALLOCATE(RecvMsg%SideDistance(1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideDistance',RecvMsg%nSides)
  RecvMsg%SideDistance(:)=0.
  ALLOCATE(RecvMsg%SideNormVec(1:3,1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideNormVec',RecvMsg%nSides)
  RecvMsg%SideNormVec(:,:)=0.
END IF

! BC Mapping
! BC(1:4,1:nSides)
! 1:BC,
! 2:NBProc,
! 3:1=Mine/2=Yours
! 4:SIDE_ID-MPI_Offset(NBProc)
IF (SendMsg%nSides.GT.0) THEN
  ALLOCATE(SendMsg%BC(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%BC',SendMsg%nSides,999.)
  SendMsg%BC(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BC(1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%BC',RecvMsg%nSides,999.)
  RecvMsg%BC(:)=0
END IF
! mortartype
IF (SendMsg%nSides.GT.0) THEN
  ALLOCATE(SendMsg%MortarType(1:2,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%MortarType',SendMsg%nSides,999.)
  SendMsg%MortarType(:,:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%MortarType(1:2,1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%MortarType',RecvMsg%nSides,999.)
  RecvMsg%MortarType(:,:)=0
END IF
! SideBCType
IF (SendMsg%nSides.GT.0) THEN
  ALLOCATE(SendMsg%SideBCType(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideBCType',SendMsg%nSides,999.)
  SendMsg%SideBCType(:)=0
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideBCType(1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideBCType',RecvMsg%nSides,999.)
  RecvMsg%SideBCType(:)=0
END IF
! NativeElemID
IF (SendMsg%nElems.GT.0) THEN
  ALLOCATE(SendMsg%NativeElemID(1:SendMsg%nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%NativeElemID',SendMsg%nElems,999.)
  SendMsg%NativeElemID(:)=0
END IF
IF (RecvMsg%nElems.GT.0) THEN
  ALLOCATE(RecvMsg%NativeElemID(1:RecvMsg%nElems),STAT=ALLOCSTAT)  ! Save E2S_SIDE_ID, E2S_FLIP
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%NativeElemID',RecvMsg%nElems,999.)
  RecvMsg%NativeElemID(:)=0
END IF
! SideSlabNormals Mapping
IF (SendMsg%nSides.GT.0) THEN
  ALLOCATE(SendMsg%SideSlabNormals(1:3,1:3,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideSlabNormals',SendMsg%nSides)
  SendMsg%SideSlabNormals(:,:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideSlabNormals(1:3,1:3,1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideSlabNormals',RecvMsg%nSides)
  RecvMsg%SideSlabNormals(:,:,:)=0.
END IF
! SideSlabIntervals Mapping
IF (SendMsg%nSides.GT.0) THEN       ! SideSlabIntervals(1:2,1:nSides)
  ALLOCATE(SendMsg%SideSlabIntervals(1:6,1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%SideSlabIntervals',SendMsg%nSides)
  SendMsg%SideSlabIntervals(:,:)=0.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%SideSlabIntervals(1:6,1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%SideSlabIntervals',RecvMsg%nSides)
  RecvMsg%SideSlabIntervals(:,:)=0.
END IF
! BoundingBoxIsEmpty Mapping
IF (SendMsg%nSides.GT.0) THEN       ! BoundingBoxIsEmpty(1:2,1:nSides)
  ALLOCATE(SendMsg%BoundingBoxIsEmpty(1:SendMsg%nSides),STAT=ALLOCSTAT)  ! see boltzplatz.h
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate SendMsg%BoundingBoxIsEmpty',SendMsg%nSides)
  SendMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF
IF (RecvMsg%nSides.GT.0) THEN
  ALLOCATE(RecvMsg%BoundingBoxIsEmpty(1:RecvMsg%nSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate RecvMsg%BoundingBoxIsEmpty',RecvMsg%nSides)
  RecvMsg%BoundingBoxIsEmpty(:)=.FALSE.
END IF

! fill send buffers with node, side and element data (including connectivity!)
! ElemtoSide
DO iElem = 1,nElems
  IF (ElemIndex(iElem).NE.0) THEN
    IF(DoRefMapping)THEN
      SendMsg%XCL_NGeo(:,:,:,:,ElemIndex(iElem))=XCL_NGeo(:,:,:,:,iElem)
      SendMsg%dXCL_NGeo(:,:,:,:,:,ElemIndex(iElem))=dXCL_NGeo(:,:,:,:,:,iElem)
    END IF
    SendMsg%CurvedElem(ElemIndex(iElem))=CurvedElem(iElem)
    IF (UseAuxBCs) SendMsg%ElemHasAuxBCs(ElemIndex(iElem),:)=ElemHasAuxBCs(iElem,:)
    IF (.NOT.DoRefMapping) THEN
      SendMsg%ElemType(ElemIndex(iElem))=ElemType(iElem)
    END IF
    IF(TriaTracking)THEN
      SendMsg%ConcaveElemSide(:,ElemIndex(iElem))=GEO%ConcaveElemSide(:,iElem)
      DO iNode=1,8
        IF(NodeIndex(GEO%ElemToNodeID(iNode,iElem)).GT.0) THEN
          SendMsg%ElemToNodeID(iNode,ElemIndex(iElem))=NodeIndex(GEO%ElemToNodeID(iNode,iElem))
    END IF
      END DO
      DO iLocSide = 1,6
        DO iNode = 1,4
          IF(NodeIndex(GEO%ElemSideNodeID(iNode,iLocSide,iElem)).GT.0) THEN
            SendMsg%ElemSideNodeID(iNode,iLocSide,ElemIndex(iElem))=NodeIndex(GEO%ElemSideNodeID(iNode,iLocSide,iElem))
          END IF
        END DO
      END DO
    END IF
    SendMsg%ElemBaryNGeo(:,ElemIndex(iElem))=ElemBaryNGeo(:,iElem)
    DO iLocSide = 1,6
      SideID=PartElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      IF(SideID.LT.1) CYCLE
      IF(SideID.GT.nPartSides) CYCLE
      IF(SideIndex(SideID).GT.0)THEN
        SendMsg%ElemToSide(1,iLocSide,ElemIndex(iElem)) = SideIndex(SideID)
        SendMsg%ElemToSide(2,iLocSide,ElemIndex(iElem)) = &
                PartElemToSide(2,iLocSide,iElem)
      END IF
    END DO
  END IF
END DO

! part of elemtoelemglob which is needed
DO iElem = 1,nElems
  IF (ElemIndex(iElem).NE.0) THEN
    SendMsg%ElemToElemGlob(1:4,1:6,ElemIndex(iElem)) = &
            PartElemToElemGlob(1:4,1:6,iElem)
  END IF
END DO

! coordinates of nodes to be sent
DO iNode = 1,nNodes
  IF (NodeIndex(iNode).NE.0) THEN
    SendMsg%NodeCoords(:,NodeIndex(iNode))=GEO%Nodecoords(:,iNode)
  END IF
END DO

! SideToElem Mapping & BezierControloints3D
DO iSide = 1,nPartSides
  IF (SideIndex(iSide).GT.0) THEN
  !IF(isSide(iSide))THEN
    DO iIndex = 1,2   ! S2E_ELEM_ID, S2E_NB_ELEM_ID
      IF (PartSideToElem(iIndex,iSide).GT.0) THEN
         SendMsg%SideToElem(iIndex,SideIndex(iSide)) = &
                 ElemIndex(PartSideToElem(iIndex,iSide))
         SendMsg%SideBCType(SideIndex(iSide)) = SidePeriodicType(iSide)
      END IF
    END DO ! S2E_LOC_SIDE_ID, S2E_NB_LOC_SIDE_ID, S2E_FLIP
    IF(DoRefMapping)THEN
      SendMsg%SideType(SideIndex(iSide)) = SideType(PartBCSideList(iSide))
      SendMsg%SideDistance(SideIndex(iSide)) = SideDistance(PartBCSideList(iSide))
      SendMsg%SideNormVec(1:3,SideIndex(iSide)) = SideNormVec(1:3,PartBCSideList(iSide))
    ELSE
      SendMsg%SideType(SideIndex(iSide)) = SideType(iSide)
      SendMsg%SideDistance(SideIndex(iSide)) = SideDistance(iSide)
      SendMsg%SideNormVec(1:3,SideIndex(iSide)) = SideNormVec(1:3,iSide)
    END IF
    SendMsg%SideToElem(3:5,SideIndex(iSide)) = &
        PartSideToElem(3:5,iSide)
    SendMsg%BezierControlPoints3D(:,:,:,SideIndex(iSide)) = &
        BezierControlPoints3D(:,:,:,iSide)
    ! slabnormals
    SendMsg%SideSlabNormals(:,:,SideIndex(iSide)) = &
        SideSlabNormals(:,:,iSide)
    ! slabintervalls
    SendMsg%SideSlabIntervals(:,SideIndex(iSide)) = &
        SideSlabIntervals(:,iSide)
    ! BoundingBoxIsEmpty
    SendMsg%BoundingBoxIsEmpty(SideIndex(iSide)) = &
        BoundingBoxIsEmpty(iSide)
  END IF
END DO
!--- BC Mapping ------------------------------------------------------!
DO iSide = 1,nPartSides  ! no need to go through all side since BC(1:nBCSides)
  IF (SideIndex(iSide).NE.0) THEN
    SendMsg%BC(SideIndex(iSide)) = BC(iSide)
    SendMsg%MortarType(:,SideIndex(iSide))=MortarType(:,iSide)
  END IF
END DO

! NativeElemID
DO iElem = 1,nElems
  IF (ElemIndex(iElem).NE.0) THEN
    SendMsg%NativeElemID(ElemIndex(iElem)) = iElem
  END IF
END DO

dataSize=3*(NGeo+1)*(NGeo+1)
dataSize2=3*(NGeo+1)*(NGeo+1)*(NGeo+1)
dataSize3=9*(NGeo+1)*(NGeo+1)*(NGeo+1)

IF (PartMPI%MyRank.LT.iProc) THEN
  ! Send:
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BezierControlPoints3D,SendMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideBCType,SendMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabNormals,SendMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabIntervals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)
  IF(DoRefMapping)THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)
  END IF
  IF(TriaTracking)THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemToNodeID,SendMsg%nElems*8,MPI_INTEGER,iProc,1114,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemSideNodeID,SendMsg%nElems*6*4,MPI_INTEGER,iProc,1115,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ConcaveElemSide,SendMsg%nElems*6,MPI_LOGICAL,iProc,1116,PartMPI%COMM,IERROR)
    IF (SendMsg%nNodes.GT.0) &
        CALL MPI_SEND(SendMsg%NodeCoords,SendMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemBaryNGeo,SendMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToElemGlob,SendMsg%nElems*24,MPI_LONG       ,iProc,1119,PartMPI%COMM,IERROR)
    IF (SendMsg%nSides.GT.0) &
        CALL MPI_SEND(SendMsg%MortarType,SendMsg%nSides*2,MPI_INTEGER       ,iProc,1120,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%CurvedElem,SendMsg%nElems,MPI_LOGICAL,iProc,1121,PartMPI%COMM,IERROR)
  IF (.NOT.DoRefMapping) THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemType  ,SendMsg%nElems,MPI_INTEGER,iProc,1122,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideType    ,SendMsg%nSides  ,MPI_INTEGER,iProc,1123,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideDistance,SendMsg%nSides  ,MPI_DOUBLE_PRECISION,iProc,1124,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideNormVec ,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1125,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0 .AND. UseAuxBCs) &
      CALL MPI_SEND(SendMsg%ElemHasAuxBCs,SendMsg%nElems*nAuxBCs,MPI_LOGICAL,iProc,1126,PartMPI%COMM,IERROR)

  ! Receive:
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
  CALL MPI_RECV(RecvMsg%BezierControlPoints3D,RecvMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,MPISTATUS,IERROR)

  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) CALL MPI_RECV(RecvMsg%SideBCType,RecvMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabNormals,RecvMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabIntervals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
  IF(DoRefMapping)THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF(TriaTracking)THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemToNodeID,RecvMsg%nElems*8,MPI_INTEGER,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemSideNodeID,RecvMsg%nElems*6*4,MPI_INTEGER,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ConcaveElemSide,RecvMsg%nElems*6,MPI_LOGICAL,iProc,1116,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nNodes.GT.0) &
        CALL MPI_RECV(RecvMsg%NodeCoords,RecvMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemBaryNGeo,RecvMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemToElemGlob,RecvMsg%nElems*24,MPI_LONG,iProc,1119,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nSides.GT.0) &
        CALL MPI_RECV(RecvMsg%MortarType,RecvMsg%nSides*2,MPI_INTEGER,iProc,1120,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%CurvedElem,RecvMsg%nElems,MPI_LOGICAL,iProc,1121,PartMPI%COMM,MPISTATUS,IERROR)
  IF (.NOT.DoRefMapping) THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemType  ,RecvMsg%nElems,MPI_INTEGER,iProc,1122,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideType    ,RecvMsg%nSides  ,MPI_INTEGER,iProc,1123,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideDistance,RecvMsg%nSides  ,MPI_DOUBLE_PRECISION,iProc,1124,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideNormVec ,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1125,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0 .AND. UseAuxBCs) &
      CALL MPI_RECV(RecvMsg%ElemHasAuxBCs,RecvMsg%nElems*nAuxBCs,MPI_LOGICAL,iProc,1126,PartMPI%COMM,MPISTATUS,IERROR)
ELSE IF (PartMPI%MyRank.GT.iProc) THEN
  ! Receive:
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%ElemToSide,RecvMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%SideToElem,RecvMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
  CALL MPI_RECV(RecvMsg%BezierControlPoints3D,RecvMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
    CALL MPI_RECV(RecvMsg%BC,RecvMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
    CALL MPI_RECV(RecvMsg%NativeElemID,RecvMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) CALL MPI_RECV(RecvMsg%SideBCType,RecvMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabNormals,RecvMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideSlabIntervals,RecvMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%BoundingBoxIsEmpty,RecvMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,MPISTATUS,IERROR)
  IF(DoRefMapping)THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%XCL_NGeo,RecvMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%dXCL_NGeo,RecvMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF(TriaTracking)THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemToNodeID,RecvMsg%nElems*8,MPI_INTEGER,iProc,1114,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemSideNodeID,RecvMsg%nElems*6*4,MPI_INTEGER,iProc,1115,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ConcaveElemSide,RecvMsg%nElems*6,MPI_LOGICAL,iProc,1116,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nNodes.GT.0) &
        CALL MPI_RECV(RecvMsg%NodeCoords,RecvMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemBaryNGeo,RecvMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%ElemToElemGlob,RecvMsg%nElems*24,MPI_LONG,iProc,1119,PartMPI%COMM,MPISTATUS,IERROR)
    IF (RecvMsg%nSides.GT.0) &
        CALL MPI_RECV(RecvMsg%MortarType,RecvMsg%nSides*2,MPI_INTEGER,iProc,1120,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0) &
      CALL MPI_RECV(RecvMsg%CurvedElem,RecvMsg%nElems,MPI_LOGICAL,iProc,1121,PartMPI%COMM,MPISTATUS,IERROR)
  IF (.NOT.DoRefMapping) THEN
    IF (RecvMsg%nElems.GT.0) &
        CALL MPI_RECV(RecvMsg%ElemType  ,RecvMsg%nElems,MPI_INTEGER,iProc,1122,PartMPI%COMM,MPISTATUS,IERROR)
  END IF
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideType    ,RecvMsg%nSides  ,MPI_INTEGER,iProc,1123,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideDistance,RecvMsg%nSides  ,MPI_DOUBLE_PRECISION,iProc,1124,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nSides.GT.0) &
      CALL MPI_RECV(RecvMsg%SideNormVec ,RecvMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1125,PartMPI%COMM,MPISTATUS,IERROR)
  IF (RecvMsg%nElems.GT.0 .AND. UseAuxBCs) &
      CALL MPI_RECV(RecvMsg%ElemHasAuxBCs,RecvMsg%nElems*nAuxBCs,MPI_LOGICAL,iProc,1126,PartMPI%COMM,MPISTATUS,IERROR)

  ! Send:
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%ElemToSide,SendMsg%nElems*2*6,MPI_INTEGER       ,iProc,1104,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideToElem,SendMsg%nSides*5,MPI_INTEGER         ,iProc,1105,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BezierControlPoints3D,SendMsg%nSides*datasize,MPI_DOUBLE_PRECISION,iProc,1106,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%BC,SendMsg%nSides,MPI_INTEGER                 ,iProc,1107,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) CALL MPI_SEND(SendMsg%NativeElemID,SendMsg%nElems,MPI_INTEGER         ,iProc,1108,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) CALL MPI_SEND(SendMsg%SideBCType,SendMsg%nSides,MPI_INTEGER,iProc,1110,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabNormals,SendMsg%nSides*9,MPI_DOUBLE_PRECISION,iProc,1111,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideSlabIntervals,SendMsg%nSides*6,MPI_DOUBLE_PRECISION,iProc,1112,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%BoundingBoxIsEmpty,SendMsg%nSides,MPI_LOGICAL,iProc,1113,PartMPI%COMM,IERROR)

  IF(DoRefMapping)THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%XCL_NGeo,SendMsg%nElems*datasize2,MPI_DOUBLE_PRECISION,iProc,1114,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%dXCL_NGeo,SendMsg%nElems*datasize3,MPI_DOUBLE_PRECISION,iProc,1115,PartMPI%COMM,IERROR)
  END IF
  IF(TriaTracking)THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemToNodeID,SendMsg%nElems*8,MPI_INTEGER,iProc,1114,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemSideNodeID,SendMsg%nElems*6*4,MPI_INTEGER,iProc,1115,PartMPI%COMM,IERROR)
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ConcaveElemSide,SendMsg%nElems*6,MPI_LOGICAL,iProc,1116,PartMPI%COMM,IERROR)
    IF (SendMsg%nNodes.GT.0) &
        CALL MPI_SEND(SendMsg%NodeCoords,SendMsg%nNodes*3,MPI_DOUBLE_PRECISION,iProc,1117,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemBaryNGeo,SendMsg%nElems*3,MPI_DOUBLE_PRECISION,iProc,1118,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%ElemToElemGlob,SendMsg%nElems*24,MPI_LONG,iProc,1119,PartMPI%COMM,IERROR)
    IF (SendMsg%nSides.GT.0) &
        CALL MPI_SEND(SendMsg%MortarType,SendMsg%nSides*2,MPI_INTEGER       ,iProc,1120,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0) &
      CALL MPI_SEND(SendMsg%CurvedElem,SendMsg%nElems,MPI_LOGICAL,iProc,1121,PartMPI%COMM,IERROR)
  IF (.NOT.DoRefMapping) THEN
    IF (SendMsg%nElems.GT.0) &
        CALL MPI_SEND(SendMsg%ElemType  ,SendMsg%nElems,MPI_INTEGER,iProc,1122,PartMPI%COMM,IERROR)
  END IF
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideType    ,SendMsg%nSides  ,MPI_INTEGER,iProc,1123,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideDistance,SendMsg%nSides  ,MPI_DOUBLE_PRECISION,iProc,1124,PartMPI%COMM,IERROR)
  IF (SendMsg%nSides.GT.0) &
      CALL MPI_SEND(SendMsg%SideNormVec ,SendMsg%nSides*3,MPI_DOUBLE_PRECISION,iProc,1125,PartMPI%COMM,IERROR)
  IF (SendMsg%nElems.GT.0 .AND. UseAuxBCs) &
      CALL MPI_SEND(SendMsg%ElemHasAuxBCs,SendMsg%nElems*nAuxBCs,MPI_LOGICAL,iProc,1126,PartMPI%COMM,IERROR)
END IF

IF ((RecvMsg%nElems.EQ.0) .AND. (RecvMsg%nSides.GT.0))THEN
    ERRWRITE(*,*)'ERROR: nElems=0 when nSides=',RecvMsg%nSides,' and nSides=',RecvMsg%nSides,'!'
    CALL abort(&
     __STAMP__&
     ,'nElems=0 while nSides=',RecvMsg%nSides)
END IF

DEALLOCATE(isElem,isSide,ElemIndex,SideIndex)

IF(DoRefMapping)THEN
  IF (RecvMsg%nElems.GT.0) THEN
    ! now, the famous reconstruction of geometry
    ! add the halo region to the existing geometry
    ! therefore, the PartElemToSide,... has to be extended
    ! multiple sides ( MPI-Sides) should be ignored

    ! get number of double sides
    ! BezierControlPoints are build in master system, therefore, the node indicies are unique and 0,0 and NGeo,NGeo are on diagonal,
    ! opposide sides of the side
    nDoubleSides=0
    ALLOCATE(isSide(1:RecvMsg%nSides))
    ALLOCATE(isDone(1:RecvMsg%nSides))
    isDone=.FALSE.
    isSide=.TRUE.
    ! get increament for each halo side
    ! 1) get increment of halo side id
    ! HaloSideID is increment of new side
    ALLOCATE(HaloInc(1:RecvMsg%nSides))
    HaloInc=0
    HaloSideID=0
    DO iHaloSide=1,RecvMsg%nSides
      IF(isSide(iHaloSide))THEN
        HaloSideID=HaloSideID+1
        HaloInc(iHaloSide)=HaloSideID
      END IF
    END DO ! iHaloSide

    ! new number of sides
    tmpnSides    =nTotalSides
    tmpnElems    =nTotalElems
    tmpBCSides   =nTotalBCSides
    nTotalSides  =nTotalSides+RecvMsg%nSides-nDoubleSides
    nTotalBCSides=nTotalBCSides+RecvMsg%nSides-nDoubleSides
    nTotalElems  =nTotalElems+RecvMsg%nElems
    CALL ResizeParticleMeshData(tmpnSides,tmpnElems,nTotalSides,nTotalElems,nOldBCSides=tmpBCSides,nTotalBCSides=nTotalBCSides)

    ! loop over all elements and add them
    DO iElem=1,RecvMsg%nElems
      newElemID=tmpnElems+iElem
      PartHaloElemToProc(NATIVE_ELEM_ID,newElemId)=RecvMsg%NativeElemID(iElem)
      PartHaloElemToProc(NATIVE_PROC_ID,newElemId)=iProc
      XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
      dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,newElemID)=RecvMsg%dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,iElem)
      ElemBaryNGeo(1:3,newElemID) = RecvMsg%ElemBaryNGeo(1:3,iElem)
      CurvedElem(newElemID)       = RecvMsg%CurvedElem(iElem)
      IF (UseAuxBCs) THEN
        ElemHasAuxBCs(newElemID,:)  = RecvMsg%ElemHasAuxBCs(iElem,:)
      END IF
    END DO

    ! loop over all sides and add them
    DO iSide=1,RecvMsg%nSides
      newSideID  =tmpnSides+iSide
      ! BezierControlPoints and BoundingBox stuff
      BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,iSide)
      ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
      SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,iSide)
      SideSlabIntervals(1:6,newSideID) =RecvMsg%SideSlabIntervals(1:6 ,iSide)
      BoundingBoxIsEmpty(newSideID) =RecvMsg%BoundingBoxIsEmpty( iSide)
      ! BC type, etc.
      BC(newSideID)             =RecvMsg%BC(iSide)
      MortarType(1:2,newSideID) =RecvMsg%MortarType(1:2,iSide)
      SidePeriodicType(newSideID)=RecvMsg%SideBCType(iSide)
      PartBCSideList(newSideID)=tmpBCSides+iSide
      newSideID = tmpBCSides+iSide
      SideType(newSideID)        =RecvMsg%SideType(iSide)
      SideDistance(newSideID)    =RecvMsg%SideDistance(iSide)
      SideNormVec(1:3,newSideID) =RecvMsg%SideNormVec(1:3,iSide)
    END DO

    ! fill lists
    ! caution: PartSideToElem is only filled for BC sides
    DO iElem=1,RecvMsg%nElems
      newElemID=tmpnElems+iElem
      DO ilocSide=1,6
        SideID=RecvMsg%ElemToSide(1,iLocSide,iElem)
        IF(SideID.GT.0)THEN
          ! fill PartElemToSide
          newSideID=tmpnSides+SideID
          PartElemToSide(E2S_SIDE_ID,iLocSide,NewElemID)=newSideID
          PartElemToSide(E2S_FLIP,iLocSide,NewElemID)=0
          ! and SideToElem
          PartSideToElem(S2E_ELEM_ID,newSideID)=newElemID
          PartSideToElem(2:5,newSideID)=RecvMsg%SideToElem(2:5,SideID)
        END IF
      END DO
      ! list from ElemToElemGlob mapped to process local element
      ! new list points from local-elem-id to global
      PartElemToElemGlob(1:4,1:6,newElemID) = RecvMsg%ElemToElemGlob(1:4,1:6,iElem)
    END DO

    IF(.NOT.PartMPI%isMPINeighbor(iProc))THEN
      PartMPI%isMPINeighbor(iProc) = .true.
      PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
    END IF
    DEALLOCATE(isSide)
    DEALLOCATE(isDone)
    DEALLOCATE(HaloInc)
    DO iSide=1,nTotalSides
      IF(SUM(ABS(SideSlabNormals(:,:,iSide))).EQ.0)THEN
        CALL abort(&
__STAMP__&
          ,' SideSlabNormals is zero!,iSide',iSide)
      END IF
    END DO
  END IF ! RecvMsg%nSides>0
ELSE ! DoRefMappping=F
 IF (RecvMsg%nSides.GT.0) THEN
    ! now, the famous reconstruction of geometry
    ! add the halo region to the existing geometry
    ! therefore, the PartElemToSide,... has to be extended
    ! multiple sides ( MPI-Sides) should are duplicated to deal easier with
    ! mortar faces, the elemtoelemandside lists are reconstructed via the global element list
    ! this loop is not optimized, hence, the MPI-faces exists twice,
    ! once for the master and the slave element

    nDoubleSides=0
    ALLOCATE(isSide(1:RecvMsg%nSides))
    ALLOCATE(isDone(1:RecvMsg%nSides))
    isDone=.FALSE.
    isSide=.TRUE.
    ! get increament for each halo side
    ! 1) get increment of halo side id
    ! HaloSideID is increment of new side
    ALLOCATE(HaloInc(1:RecvMsg%nSides))
    HaloInc=0
    HaloSideID=0
    DO iHaloSide=1,RecvMsg%nSides
      IF(isSide(iHaloSide))THEN
        HaloSideID=HaloSideID+1
        HaloInc(iHaloSide)=HaloSideID
      END IF
    END DO ! iHaloSide

    tmpnSides =nTotalSides
    tmpnElems =nTotalElems
    tmpnNodes =nTotalNodes
    nTotalSides=nTotalSides+RecvMsg%nSides-nDoubleSides
    nTotalElems=nTotalElems+RecvMsg%nElems
    nTotalNodes=nTotalNodes+RecvMsg%nNodes
    CALL ResizeParticleMeshData(tmpnSides,tmpnElems,nTotalSides,nTotalElems,nOldNodes=tmpnNodes,nTotalNodes=nTotalNodes)

    !DO iElem=tmpnElems+1,nTotalElems
    DO iElem=1,RecvMsg%nElems
      ! first, new SideID=entry of RecvMsg+tmpnSides
      newElemID=tmpnElems+iElem
      DO ilocSide=1,6
        haloSideID=RecvMsg%ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
        ! first, set new sideid
        IF(isSide(haloSideID)) THEN
          IF(HaloInc(haloSideID).EQ.0) IPWRITE(UNIT_stdOut,*) ' Warning: wrong halo inc'
          newSideID=tmpnSides+haloinc(haloSideID)
          IF(newSideID.LT.tmpnSides) IPWRITE(UNIT_stdOut,*) 'Warning: wrong new sideid', newsideid
        END IF
        ! build PartSideToElem, so much as possible
        ! get correct side out of RecvMsg%SideToElem
        IF(iElem.EQ.RecvMsg%SideToElem(S2E_ELEM_ID,haloSideID))THEN
          PartSideToElem(S2E_ELEM_ID,newSideID)    =newElemID
          PartSideToElem(S2E_LOC_SIDE_ID,newSideID)=ilocSide
          PartSideToElem(S2E_FLIP       ,newSideID)=RecvMsg%SideToElem(S2E_FLIP,haloSideID)
        ELSE IF(iElem.EQ.RecvMsg%SideToElem(S2E_NB_ELEM_ID,haloSideID))THEN
          PartSideToElem(S2E_NB_ELEM_ID,newSideID) =newElemID
          PartSideToElem(S2E_NB_LOC_SIDE_ID,newSideID) =ilocSide
          PartSideToElem(S2E_FLIP      ,newSideID) =RecvMsg%SideToElem(S2E_FLIP,haloSideID)
        ELSE ! should be found, because there should be halo sides without any connection
            CALL abort(&
 __STAMP__&
          ,'Non-Critical error in domain reconstrution. IF NOT encountered, something is terrible wrong.')
        END IF
        IF(.NOT.isDone(haloSideID))THEN
          BC(newSideID)=RecvMsg%BC(haloSideID)
          SidePeriodicType(newSideID)=RecvMsg%SideBCType(haloSideID)

          ! copy Bezier to new side id
          BezierControlPoints3D(1:3,0:NGeo,0:NGeo,newSideID)=RecvMsg%BezierControlpoints3D(1:3,0:NGeo,0:NGeo,haloSideID)
          ! SlabBoundingBox has to be sent because only BezierPoints of Slave-Sides are received
          SideSlabNormals(1:3,1:3,newSideID)=RecvMsg%SideSlabNormals(1:3,1:3,haloSideID)
          SideSlabIntervals(1:6,newSideID)  =RecvMsg%SideSlabIntervals(1:6 ,haloSideID)
          BoundingBoxIsEmpty(newSideID)     =RecvMsg%BoundingBoxIsEmpty( haloSideID)
          MortarType(1:2,newSideID)         =RecvMsg%MortarType(1:2,haloSideID)
          SideType(newSideID)               =RecvMsg%SideType(haloSideID)
          SideDistance(newSideID)           =RecvMsg%SideDistance(haloSideID)
          SideNormVec(1:3,newSideID)           =RecvMsg%SideNormVec(1:3,haloSideID)
          isDone(haloSideID)=.TRUE.
        END IF
        ! build entry to PartElemToSide
        PartElemToSide(1,iLocSide,newElemId)=newSideID
        PartElemToSide(2,ilocSide,newElemId)=RecvMsg%ElemToSide(2,ilocSide,iElem)
        IF (TriaTracking) THEN
          DO iNode = 1,4
            GEO%ElemSideNodeID(iNode,ilocSide,newElemID) = tmpnNodes + RecvMsg%ElemSideNodeID(iNode,ilocSide,iElem)
          END DO
        END IF
      END DO ! ilocSide
      IF (TriaTracking) THEN
        DO iNode = 1,8
          GEO%ElemToNodeID(iNode,newElemID) = tmpnNodes + RecvMsg%ElemToNodeID(iNode,iElem)
        END DO
      END IF
      ! set native elemID
      PartHaloElemToProc(NATIVE_ELEM_ID,newElemId)=RecvMsg%NativeElemID(iElem)
      PartHaloElemToProc(NATIVE_PROC_ID,newElemId)=iProc
      ElemBaryNGeo(1:3,newElemID) = RecvMsg%ElemBaryNGeo(1:3,iElem)
      CurvedElem(newElemID)       = RecvMsg%CurvedElem(iElem)
      ElemType(newElemID)         = RecvMsg%ElemType(iElem)
      IF (UseAuxBCs) THEN
        ElemHasAuxBCs(newElemID,:)  = RecvMsg%ElemHasAuxBCs(iElem,:)
      END IF
      IF (TriaTracking) THEN
        GEO%ConcaveElemSide(1:6,newElemID)    = RecvMsg%ConcaveElemSide(1:6,iElem)
      END IF
      ! list from ElemToElemGlob mapped to process local element
      ! new list points from local-elem-id to global
      PartElemToElemGlob(1:4,1:6,newElemID) = RecvMsg%ElemToElemGlob(1:4,1:6,iElem)
    END DO ! iElem
    IF (TriaTracking) THEN
      DO iNode = 1,RecvMsg%nNodes
        ! first, new NodeID=entry of RecvMsg+tmpnNodes
        GEO%NodeCoords(:,tmpnNodes+iNode) = RecvMsg%NodeCoords(:,iNode)
      END DO
    END IF
    ! build rest: PartElemToElem, PartLocSideID
    DO iElem=PP_nElems+1,nTotalElems
      DO ilocSide=1,6
        flip   = PartElemToSide(E2S_FLIP,ilocSide,iElem)
        SideID = PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
        ! check of sideid
      END DO ! ilocSide
    END DO ! Elem
    IF(.NOT.PartMPI%isMPINeighbor(iProc))THEN
      PartMPI%isMPINeighbor(iProc) = .true.
      PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
    END IF
    DEALLOCATE(isSide)
    DEALLOCATE(isDone)
    DEALLOCATE(HaloInc)
  END IF ! RecvMsg%nSides>0
END IF

END SUBROUTINE ExchangeHaloGeometry


SUBROUTINE ResizeParticleMeshData(nOldSides,nOldElems,nTotalSides,nTotalElems,nOldBCSides,nTotalBCSides,nOldNodes,nTotalNodes)
!===================================================================================================================================
! resize the partilce mesh data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartHaloElemToProc
USE MOD_Mesh_Vars,              ONLY:BC,nGeo,nElems,MortarType
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicType,PartBCSideList,GEO,ElemType,XCL_NGeo,ElemBaryNGeo,CurvedElem,DXCL_NGEO
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem,PartElemToElemGlob,ElemHasAuxBCs
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,SideType,SideNormVec,SideDistance
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping,TriaTracking
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Boundary_Vars, ONLY:nAuxBCs,UseAuxBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: nOldSides,nOldElems,nTotalSides,nTotalElems
INTEGER,INTENT(IN),OPTIONAL        :: nOldBCSides,nTotalBCSides
INTEGER,INTENT(IN),OPTIONAL        :: nOldNodes,nTotalNodes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: ALLOCSTAT,nLower
INTEGER,ALLOCATABLE                :: DummyElemToSide(:,:,:)
INTEGER,ALLOCATABLE                :: DummyBC(:)
REAL,ALLOCATABLE                   :: DummyBezierControlPoints3D(:,:,:,:)
REAL,ALLOCATABLE                   :: DummyXCL_NGEO (:,:,:,:,:)
REAL,ALLOCATABLE                   :: DummydXCL_NGEO(:,:,:,:,:,:)
REAL,ALLOCATABLE                   :: DummyElemBaryNGeo(:,:)
LOGICAL,ALLOCATABLE                :: DummyCurvedElem(:)
LOGICAL,ALLOCATABLE                :: DummyElemHasAuxBCs(:,:)
INTEGER,ALLOCATABLE                :: DummyElemType(:)
INTEGER,ALLOCATABLE                :: DummySideType(:)
REAL,ALLOCATABLE                   :: DummySideDistance(:)
REAL,ALLOCATABLE                   :: DummySideNormVec(:,:)

INTEGER,ALLOCATABLE                :: DummyElemToNodeID(:,:)
INTEGER,ALLOCATABLE                :: DummyElemSideNodeID(:,:,:)
REAL,ALLOCATABLE                   :: DummyNodeCoords(:,:)
LOGICAL,ALLOCATABLE                :: DummyConcaveElemSide(:,:)
INTEGER,ALLOCATABLE                :: DummyHaloToProc(:,:)
INTEGER,ALLOCATABLE                :: DummyMortarType(:,:)
INTEGER,ALLOCATABLE                :: DummySideToElem(:,:)
INTEGER,ALLOCATABLE                :: DummySideBCType(:),DummyPartBCSideList(:)
INTEGER(KIND=8),ALLOCATABLE        :: DummyElemToElem(:,:,:)
REAL,ALLOCATABLE,DIMENSION(:,:,:)  :: DummySideSlabNormals                 ! normal vectors of bounding slab box
REAL,ALLOCATABLE,DIMENSION(:,:)    :: DummySideSlabIntervals               ! intervalls beta1, beta2, beta3
LOGICAL,ALLOCATABLE,DIMENSION(:)   :: DummyBoundingBoxIsEmpty
!===================================================================================================================================

ALLOCATE(DummyElemToSide(1:2,1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemToSide)) CALL abort(&
  __STAMP__&
  ,'Could not allocate DummyElemToSide')
DummyElemToSide=PartElemToSide

DEALLOCATE(PartElemToSide)
ALLOCATE(PartElemToSide(1:2,1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'Could not reallocate PartElemToSide')
PartElemToSide=-1
PartElemToSide(:,:,1:nOldElems) =DummyElemToSide(:,:,1:nOldElems)
DEALLOCATE(DummyElemToSide)

IF(DoRefMapping)THEN
  ! XCL_NGeo
  ALLOCATE(DummyXCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyXCL_NGeo)) CALL abort(&
    __STAMP__&
    ,'Could not allocate DummyXCL_NGeo')
  DummyXCL_NGeo=XCL_NGeo
  DEALLOCATE(XCL_NGeo)
  ALLOCATE(XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not reallocate XCL_NGeo')
  XCL_NGeo=0.
  XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems) =DummyXCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems)
  DEALLOCATE(DummyXCL_NGeo)

  ! dXCL_NGeo
  ALLOCATE(DummydXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems))
  IF (.NOT.ALLOCATED(DummydXCL_NGeo)) CALL abort(&
    __STAMP__&
    ,'Could not allocate DummydXCL_NGeo')
  DummydXCL_NGeo=dXCL_NGeo
  DEALLOCATE(dXCL_NGeo)
  ALLOCATE(dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not reallocate dXCL_NGeo')
  dXCL_NGeo=0.
  dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems) =DummydXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nOldElems)
  DEALLOCATE(DummydXCL_NGeo)

  ! PartBCSideList
  ALLOCATE(DummyPartBCSideList(1:nOldSides))
  IF (.NOT.ALLOCATED(DummyPartBCSideList)) CALL abort(&
    __STAMP__&
    ,'Could not allocate DummyPartBCSideList')
  DummyPartBCSideList(1:nOldSides)=PartBCSideList(1:nOldSides)
  DEALLOCATE(PartBCSideList)
  ALLOCATE(PartBCSideList(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not reallocate PartBCSideList')
  PartBCSideList=-1 !HUGE(1)
  PartBCSideList(1:nOldSides) =DummyPartBCSideList(1:nOldSides)
  DEALLOCATE(DummyPartBCSideList)
END IF

! HaloToProc
IF(.NOT.ALLOCATED(PartHaloElemToProc))THEN
  nLower=nElems+1
  ALLOCATE(PartHaloElemToProc(1:3,nLower:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not allocate PartHaloElemToProc')
  PartHaloElemToProc=-1
ELSE
  nLower=nElems+1
  ALLOCATE(DummyHaloToProc(1:3,nLower:nOldElems))
  IF (.NOT.ALLOCATED(DummyHaloToProc)) CALL abort(&
    __STAMP__&
    ,'Could not allocate DummyHaloToProc')
  DummyHaloToProc=PartHaloElemToProc
  DEALLOCATE(PartHaloElemToProc)
  ALLOCATE(PartHaloElemToProc(1:3,nLower:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not reallocate PartHaloElemToProc')

  ! copy array to new
  PartHaloElemToProc=-1
  PartHaloElemToProc(1:3,PP_nElems+1:nOldElems)    =DummyHaloToProc(1:3,PP_nElems+1:nOldElems)
  DEALLOCATE(DummyHaloToProc)
END IF

! PartSideToElem
ALLOCATE(DummySideToElem(1:5,1:nOldSides))
IF (.NOT.ALLOCATED(DummySideToElem)) CALL abort(&
  __STAMP__&
  ,'Could not allocate DummySideToElem')
DummySideToElem=PartSideToElem
DEALLOCATE(PartSideToElem)
ALLOCATE(PartSideToElem(1:5,1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'Could not reallocate PartSideToElem')
PartSideToElem=-1
PartSideToElem(1:5,1:nOldSides  )              =DummySideToElem(1:5,1:nOldSides)
DEALLOCATE(DummySideToElem)

! PartElemToElemGlob
ALLOCATE(DummyElemToElem(1:4,1:6,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemToElem)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummySideToElem')
DummyElemToElem=PartElemToElemGlob
DEALLOCATE(PartElemToElemGlob)
ALLOCATE(PartElemToElemGlob(1:4,1:6,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate PartElemToElemGlob')
PartElemToElemGlob=-1
PartElemToElemGlob(1:4,1:6,1:nOldElems)            =DummyElemToElem(1:4,1:6,1:nOldElems)
DEALLOCATE(DummyElemToElem)

IF(DoRefMapping)THEN
  ! BezierControlPoints3D
  ALLOCATE(DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBezierControlPoints3D)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummyBezierControlPoints3D')
  DummyBezierControlPoints3D=BezierControlPoints3D
  DEALLOCATE(BezierControlPoints3D)
  ALLOCATE(BezierControlPoints3d(1:3,0:NGeo,0:NGeo,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate BezierControlPoints3d')
  BezierControlPoints3d(:,:,:,1:nOldSides) =DummyBezierControlPoints3D(:,:,:,1:nOldSides)
  DEALLOCATE(DummyBezierControlPoints3D)

  ! SideSlabNormals
  ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideSlabNormals)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideSlabNormals')
  DummySideSlabNormals=SideSlabNormals
  DEALLOCATE(SideSlabNormals)
  ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideSlabNormals')
  SideSlabNormals=0
  SideSlabNormals(1:3,1:3,1:nOldSides) =DummySideSlabNormals(1:3,1:3,1:nOldSides)
  DEALLOCATE(DummySideSlabNormals)

  ! SideSlabIntervals
  ALLOCATE(DummySideSlabIntervals(1:6,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideSlabIntervals)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideSlabIntervals')
  DummySideSlabIntervals=SideSlabIntervals
  DEALLOCATE(SideSlabIntervals)
  ALLOCATE(SideSlabIntervals(1:6,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideSlabIntervals')
  SideSlabIntervals=0
  SideSlabIntervals(1:6,1:nOldSides) =DummySideSlabIntervals(1:6,1:nOldSides)
  DEALLOCATE(DummySideSlabIntervals)

  ! BoundingBoxIsEmpty
  ALLOCATE(DummyBoundingBoxIsEmpty(1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummyBoundingBoxIsEmpty')
  DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
  DEALLOCATE(BoundingBoxIsEmpty)
  ALLOCATE(BoundingBoxIsEmpty(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate BoundingBoxIsEmpty')
  BoundingBoxIsEmpty(1:nOldSides) =DummyBoundingBoxIsEmpty(1:nOldSides)
  DEALLOCATE(DummyBoundingBoxIsEmpty)

  ! side type
  ALLOCATE(DummySideType(1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummySideType)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideType')
  DummySideType=SideType
  DEALLOCATE(SideType)
  ALLOCATE(SideType(1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideType')
  SideType=-1
  SideType(1:nOldBCSides) =DummySideType(1:nOldBCSides)

  ALLOCATE(DummySideDistance(1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummySideDistance)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideDistance')
  DummySideDistance=SideDistance
  DEALLOCATE(SideDistance)
  ALLOCATE(SideDistance(1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideDistance')
  SideDistance=0.
  SideDistance(1:nOldBCSides) =DummySideDistance(1:nOldBCSides)

  ALLOCATE(DummySideNormVec(1:3,1:nOldBCSides))
  IF (.NOT.ALLOCATED(DummySideNormVec)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideNormVec')
  DummySideNormVec=SideNormVec
  DEALLOCATE(SideNormVec)
  ALLOCATE(SideNormVec(1:3,1:nTotalBCSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideNormVec')
  SideNormVec=0.
  SideNormVec(1:3,1:nOldBCSides) =DummySideNormVec(1:3,1:nOldBCSides)

ELSE ! no mapping
  ! BezierControlPoints3D
  ALLOCATE(DummyBezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBezierControlPoints3D)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummyBezierControlPoints3D')
  DummyBezierControlPoints3D=BezierControlPoints3D
  DEALLOCATE(BezierControlPoints3D)
  ALLOCATE(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate BezierControlPoints3D')
  BezierControlPoints3D(:,:,:,1:nOldSides) =DummyBezierControlPoints3D(:,:,:,1:nOldSides)
  DEALLOCATE(DummyBezierControlPoints3D)

  ! SideSlabNormals
  ALLOCATE(DummySideSlabNormals(1:3,1:3,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideSlabNormals)) CALL abort(&
    __STAMP__&
   ,'Could not allocate DummySideSlabNormals')
  DummySideSlabNormals=SideSlabNormals
  DEALLOCATE(SideSlabNormals)
  ALLOCATE(SideSlabNormals(1:3,1:3,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideSlabNormals')
  SideSlabNormals=0
  SideSlabNormals(1:3,1:3,1:nOldSides) =DummySideSlabNormals(1:3,1:3,1:nOldSides)
  DEALLOCATE(DummySideSlabNormals)

  ! SideSlabIntervals
  ALLOCATE(DummySideSlabIntervals(1:6,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideSlabIntervals)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideSlabIntervals')
  DummySideSlabIntervals=SideSlabIntervals
  DEALLOCATE(SideSlabIntervals)
  ALLOCATE(SideSlabIntervals(1:6,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideSlabIntervals')
  SideSlabIntervals=0
  SideSlabIntervals(1:6,1:nOldSides) =DummySideSlabIntervals(1:6,1:nOldSides)
  DEALLOCATE(DummySideSlabIntervals)

  ! BoundingBoxIsEmpty
  ALLOCATE(DummyBoundingBoxIsEmpty(1:nOldSides))
  IF (.NOT.ALLOCATED(DummyBoundingBoxIsEmpty)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummyBoundingBoxIsEmpty')
  DummyBoundingBoxIsEmpty=BoundingBoxIsEmpty
  DEALLOCATE(BoundingBoxIsEmpty)
  ALLOCATE(BoundingBoxIsEmpty(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate BoundingBoxIsEmpty')
  BoundingBoxIsEmpty(1:nOldSides) =DummyBoundingBoxIsEmpty(1:nOldSides)
  DEALLOCATE(DummyBoundingBoxIsEmpty)

  ! side type
  ALLOCATE(DummySideType(1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideType)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideType')
  DummySideType=SideType
  DEALLOCATE(SideType)
  ALLOCATE(SideType(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideType')
  SideType=-1
  SideType(1:nOldSides) =DummySideType(1:nOldSides)

  ALLOCATE(DummySideDistance(1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideDistance)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideDistance')
  DummySideDistance=SideDistance
  DEALLOCATE(SideDistance)
  ALLOCATE(SideDistance(1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideDistance')
  SideDistance=0.
  SideDistance(1:nOldSides) =DummySideDistance(1:nOldSides)

  ALLOCATE(DummySideNormVec(1:3,1:nOldSides))
  IF (.NOT.ALLOCATED(DummySideNormVec)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummySideNormVec')
  DummySideNormVec=SideNormVec
  DEALLOCATE(SideNormVec)
  ALLOCATE(SideNormVec(1:3,1:nTotalSides),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate SideNormVec')
  SideNormVec=0.
  SideNormVec(1:3,1:nOldSides) =DummySideNormVec(1:3,1:nOldSides)
END IF

IF (TriaTracking) THEN
  ! Resize ElemToNodeID
  ALLOCATE(DummyElemToNodeID(1:8,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyElemToNodeID)) CALL abort(&
__STAMP__&
,'Could not allocate DummyElemToNodeID')
  DummyElemToNodeID(1:8,1:nOldElems) = GEO%ElemToNodeID(1:8,1:nOldElems)
  DEALLOCATE(GEO%ElemToNodeID)
  ALLOCATE(GEO%ElemToNodeID(1:8,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'Could not resize reallocate GEO%ElemToNodeID')
  GEO%ElemToNodeID=0
  GEO%ElemToNodeID(1:8,1:nOldElems) = DummyElemToNodeID(1:8,1:nOldElems)
  DEALLOCATE(DummyElemToNodeID)

  ! Resize ElemSideNodeID
  ALLOCATE(DummyElemSideNodeID(1:4,1:6,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyElemSideNodeID)) CALL abort(&
__STAMP__&
,'Could not allocate DummyElemSideNodeID')
  DummyElemSideNodeID(1:4,1:6,1:nOldElems) = GEO%ElemSideNodeID(1:4,1:6,1:nOldElems)
  DEALLOCATE(GEO%ElemSideNodeID)
  ALLOCATE(GEO%ElemSideNodeID(1:4,1:6,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'Could not resize reallocate GEO%ConcaveElemSide')
  GEO%ElemSideNodeID=0
  GEO%ElemSideNodeID(1:4,1:6,1:nOldElems) = DummyElemSideNodeID(1:4,1:6,1:nOldElems)
  DEALLOCATE(DummyElemSideNodeID)

  ! Resize Nodecoords
  ALLOCATE(DummyNodeCoords(1:3,1:nOldNodes))
  IF (.NOT.ALLOCATED(DummyNodeCoords)) CALL abort(&
__STAMP__&
,'Could not allocate DummyNodeCoords')
  DummyNodeCoords(1:3,1:nOldNodes) = GEO%NodeCoords(1:3,1:nOldNodes)
  DEALLOCATE(GEO%NodeCoords)
  ALLOCATE(GEO%NodeCoords(1:3,1:nTotalNodes),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'Could not resize reallocate GEO%NodeCoords')
  GEO%NodeCoords=0.
  GEO%NodeCoords(1:3,1:nOldNodes) = DummyNodeCoords(1:3,1:nOldNodes)
  DEALLOCATE(DummyNodeCoords)

  ! Resize ConcaveElemSides
  ALLOCATE(DummyConcaveElemSide(1:6,1:nOldElems))
  IF (.NOT.ALLOCATED(DummyConcaveElemSide)) CALL abort(&
__STAMP__&
,'Could not allocate DummyConcaveElemSide')
  DummyConcaveElemSide(1:6,1:nOldElems) = GEO%ConcaveElemSide(1:6,1:nOldElems)
  DEALLOCATE(GEO%ConcaveElemSide)
  ALLOCATE(GEO%ConcaveElemSide(1:6,1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'Could not resize reallocate GEO%ConcaveElemSide')
  GEO%ConcaveElemSide=.FALSE.
  GEO%ConcaveElemSide(1:6,1:nOldElems) = DummyConcaveElemSide(1:6,1:nOldElems)
  DEALLOCATE(DummyConcaveElemSide)
END IF

! ElemBaryNGeo
ALLOCATE(DummyElemBaryNGeo(1:3,1:nOldElems))
IF (.NOT.ALLOCATED(DummyElemBaryNGeo)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummyElemBaryNGeo')
DummyElemBaryNGeo=ElemBaryNGeo
DEALLOCATE(ElemBaryNGeo)
ALLOCATE(ElemBaryNGeo(1:3,1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate ElemBaryNGeo')
ElemBaryNGeo=0.
ElemBaryNGeo(1:3,1:nOldElems) =DummyElemBaryNGeo(1:3,1:nOldElems)

! curved elem
ALLOCATE(DummyCurvedElem(1:nOldElems))
IF (.NOT.ALLOCATED(DummyCurvedElem)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummyCurvedElem')
DummyCurvedElem=CurvedElem
DEALLOCATE(CurvedElem)
ALLOCATE(CurvedElem(1:nTotalElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate CurvedElem')
CurvedElem=.FALSE.
CurvedElem(1:nOldElems) =DummyCurvedElem(1:nOldElems)

! ElemHasAuxBCs elem
IF (UseAuxBCs) THEN
  ALLOCATE(DummyElemHasAuxBCs(1:nOldElems,1:nAuxBCs))
  IF (.NOT.ALLOCATED(DummyElemHasAuxBCs)) CALL abort(&
    __STAMP__&
    ,'Could not allocate DummyElemHasAuxBCs')
  DummyElemHasAuxBCs=ElemHasAuxBCs
  DEALLOCATE(ElemHasAuxBCs)
  ALLOCATE(ElemHasAuxBCs(1:nTotalElems,1:nAuxBCs),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'Could not reallocate ElemHasAuxBCs')
  ElemHasAuxBCs=.FALSE.
  ElemHasAuxBCs(1:nOldElems,:) =DummyElemHasAuxBCs(1:nOldElems,:)
END IF

! Elem Type
IF (.NOT.DoRefMapping) THEN
  ALLOCATE(DummyElemType(1:nOldElems))
  IF (.NOT.ALLOCATED(DummyElemType)) CALL abort(&
      __STAMP__&
   ,'Could not allocate DummyElemType')
  DummyElemType=ElemType
  DEALLOCATE(ElemType)
  ALLOCATE(ElemType(1:nTotalElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__&
   ,'Could not reallocate ElemType')
  ElemType=-1
  ElemType(1:nOldElems) =DummyElemType(1:nOldElems)
END IF

! SideBCType
ALLOCATE(DummySideBCType(1:nOldSides))
IF (.NOT.ALLOCATED(DummySideBCType)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummySideBCType')
DummySideBCType(1:nOldSides)=SidePeriodicType(1:nOldSides)
DEALLOCATE(SidePeriodicType)
ALLOCATE(SidePeriodicType(1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate SidePeriodicType')
SidePeriodicType=-1
SidePeriodicType(1:nOldSides) =DummySideBCType(1:nOldSides)
DEALLOCATE(DummySideBCType)

! BC
ALLOCATE(DummyBC(1:nOldSides))
IF (.NOT.ALLOCATED(DummyBC)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummyBC')

! check
DummyBC(1:nOldSides)=BC(1:nOldSides)
DEALLOCATE(BC)
ALLOCATE(BC(1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate BC')
BC=0
BC(1:nOldSides) =DummyBC(1:nOldSides)
DEALLOCATE(DummyBC)
! finished copying

! MortarType
ALLOCATE(DummyMortarType(1:2,1:nOldSides))
IF (.NOT.ALLOCATED(DummyMortarType)) CALL abort(&
    __STAMP__&
 ,'Could not allocate DummyMortarType')

! check
DummyMortarType(1:2,1:nOldSides)=MortarType(1:2,1:nOldSides)
DEALLOCATE(MortarType)
ALLOCATE(MortarType(1:2,1:nTotalSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
 ,'Could not reallocate MortarType')
MortarType=-1
MortarType(1:2,1:nOldSides) =DummyMortarType(1:2,1:nOldSides)
DEALLOCATE(DummyMortarType)

END SUBROUTINE ResizeParticleMeshData


SUBROUTINE WriteParticlePartitionInformation(nPlanar,nBilinear,nCurved,nPlanarHalo,nBilinearHalo,nCurvedHalo &
                                        ,nBCElems,nLinearElems,nCurvedElems,nBCElemsHalo,nLinearElemsHalo,nCurvedElemsHalo)
!===================================================================================================================================
! write the particle partition information to file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,              ONLY:nSides,nElems
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
USE MOD_Particle_Mesh_Vars,     ONLY:nTotalSides,nTotalElems
USE MOD_LoadBalance_Vars,       ONLY:DoLoadBalance,nLoadBalanceSteps, writePartitionInfo
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              ::nPlanar,nBilinear,nCurved,nPlanarHalo,nBilinearHalo,nCurvedHalo
INTEGER,INTENT(IN)              ::nBCElems,nLinearElems,nCurvedElems,nBCElemsHalo,nLinearElemsHalo,nCurvedElemsHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(LEN=10)          :: formatstr
INTEGER,ALLOCATABLE        :: nNBProcs_glob(:), ProcInfo_glob(:,:),NBInfo_glob(:,:), NBInfo(:), tmparray(:,:),ProcInfo(:)
REAL,ALLOCATABLE           :: tmpreal(:,:)
INTEGER                    :: nVars,nNBmax,i,j,ioUnit
CHARACTER(LEN=64)          :: filename
CHARACTER(LEN=4)           :: hilf
!===================================================================================================================================

IF(.NOT.WritePartitionInfo) RETURN

nVars=12
IF(DoRefMapping) nVars=16
Allocate(ProcInfo(nVars))

!output partitioning info
ProcInfo( 1)=nElems
ProcInfo( 2)=nSides
ProcInfo( 3)=nTotalElems-nElems ! number of halo elems
ProcInfo( 4)=nTotalSides-nSides ! number of halo sides
ProcInfo( 5)=nPlanar
ProcInfo( 6)=nBilinear
ProcInfo( 7)=nCurved
ProcInfo( 8)=nPlanarHalo
ProcInfo( 9)=nBilinearHalo
ProcInfo(10)=nCurvedHalo
ProcInfo(11)=nLinearElems
ProcInfo(12)=nCurvedElems
IF(DoRefmapping)THEN
  ProcInfo(13)=nLinearElemsHalo
  ProcInfo(14)=nCurvedElemsHalo
  ProcInfo(15)=nBCElems
  ProcInfo(16)=nBCElemsHalo
END IF

IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:PartMPI%nProcs-1))
  ALLOCATE(ProcInfo_glob(nVars,0:PartMPI%nProcs-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot
CALL MPI_GATHER(PartMPI%nMPINeighbors,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,PartMPI%COMM,iError)
CALL MPI_GATHER(ProcInfo,nVars,MPI_INTEGER,ProcInfo_glob,nVars,MPI_INTEGER,0,PartMPI%COMM,iError)
IF(MPIroot)THEN
  nNBmax=MAXVAL(nNBProcs_glob) !count, total number of columns in table
  ALLOCATE(NBinfo_glob(nNBmax,0:PartMPI%nProcs-1))
  NBinfo_glob=-77777
ELSE
  ALLOCATE(NBinfo_glob(1,1)) !dummy for debug
END IF
CALL MPI_BCAST(nNBmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
ALLOCATE(NBinfo(nNbmax))
NBinfo(1:PartMPI%nMPINeighbors)=PartMPI%MPINeighbor(1:PartMPi%nMPINeighbors)
CALL MPI_GATHER(NBinfo,nNBmax,MPI_INTEGER,NBinfo_glob,nNBmax,MPI_INTEGER,0,PartMPI%COMM,iError)
DEALLOCATE(NBinfo)
IF(MPIroot)THEN
  ioUnit=GETFREEUNIT()
  IF(DoLoadBalance)THEN
    WRITE( hilf,'(I4.4)') nLoadBalanceSteps
    filename='particlepartitionInfo-'//TRIM(hilf)//'.out'
  ELSE
    filename='particlepartitionInfo.out'
  END IF
  OPEN(UNIT=ioUnit,FILE=filename,STATUS='REPLACE')
  WRITE(ioUnit,*)'Particle Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',PartMPI%nProcs
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  IF(DoRefMapping)THEN
    WRITE(ioUnit,'(18(A20))')'Rank','nElems','nSides','nHaloElems','nHaloSides','nPlanar','nBilinear','nCurved' &
                                   ,'nPlanarHalo','nBilinearHalo','nCurvedHalo','nLinearElems','nCurvedElems'   &
                                   ,'nLinearElemsHalo','nCurvedElemsHalo','nBCElems','nBCElemsHalo','nNBProcs'
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========================================================================================='
  ELSE
    WRITE(ioUnit,'(14(A20))')'Rank','nElems','nSides','nHaloElems','nHaloSides','nPlanar','nBilinear','nCurved' &
                                   ,'nPlanarHalo','nBilinearHalo','nCurvedHalo','nLinearElems','nCurvedElems','nNBProcs'
    WRITE(ioUnit,'(A90,A90,A90,A10)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========='
  END IF

  !statistics
  ALLOCATE(tmparray(nVars+1,0:3),tmpreal(nVars+1,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=0      !max
  tmparray(:,3)=HUGE(1)   !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray( 1,0)=Procinfo_glob( 1,i) ! nElems
    tmparray( 2,0)=Procinfo_glob( 2,i) ! nSides
    tmparray( 3,0)=Procinfo_glob( 3,i) ! nHaloElems
    tmparray( 4,0)=Procinfo_glob( 4,i) ! nHaloSides
    tmparray( 5,0)=Procinfo_glob( 5,i) ! nPlanar
    tmparray( 6,0)=Procinfo_glob( 6,i) ! nBilinear
    tmparray( 7,0)=Procinfo_glob( 7,i) ! nCurved
    tmparray( 8,0)=Procinfo_glob( 8,i) ! nPlanarHalo
    tmparray( 9,0)=Procinfo_glob( 9,i) ! nBilinearHalo
    tmparray(10,0)=Procinfo_glob(10,i) ! nCurvedHalo
    tmparray(11,0)=Procinfo_glob(11,i) ! nLinearElems
    tmparray(12,0)=Procinfo_glob(12,i) ! nCurvedElems
    IF(DoRefMapping)THEN
      tmparray(13,0)=Procinfo_glob(13,i) ! nLinearElemsHalo
      tmparray(14,0)=Procinfo_glob(14,i) ! nCurvedElemsHalo
      tmparray(15,0)=Procinfo_glob(15,i) ! nBCElems
      tmparray(16,0)=Procinfo_glob(16,i) ! nBCElemsHalo
    END IF
    tmparray(nVars+1,0)=nNBProcs_glob(i)   ! nNBProcs
    DO j=1,nVars+1
      !mean
      tmparray(j,1)=tmparray(j,1)+tmparray(j,0)
      !max
      tmparray(j,2)=MAX(tmparray(j,2),tmparray(j,0))
      tmparray(j,3)=MIN(tmparray(j,3),tmparray(j,0))
    END DO !j
  END DO ! i
  tmpreal(:,1)=REAL(tmparray(:,1))/REAL(PartMPI%nProcs) !mean in REAL
  tmpreal(:,2)=0.   !RMS
  DO i=0,PartMPI%nProcs-1
    !actual proc
    tmparray( 1,0)=Procinfo_glob( 1,i) ! nElems
    tmparray( 2,0)=Procinfo_glob( 2,i) ! nSides
    tmparray( 3,0)=Procinfo_glob( 3,i) ! nHaloElems
    tmparray( 4,0)=Procinfo_glob( 4,i) ! nHaloSides
    tmparray( 5,0)=Procinfo_glob( 5,i) ! nPlanar
    tmparray( 6,0)=Procinfo_glob( 6,i) ! nBilinear
    tmparray( 7,0)=Procinfo_glob( 7,i) ! nCurved
    tmparray( 8,0)=Procinfo_glob( 8,i) ! nPlanarHalo
    tmparray( 9,0)=Procinfo_glob( 9,i) ! nBilinearHalo
    tmparray(10,0)=Procinfo_glob(10,i) ! nCurvedHalo
    tmparray(11,0)=Procinfo_glob(11,i) ! nLinearElems
    tmparray(12,0)=Procinfo_glob(12,i) ! nCurvedElems
    IF(DoRefMapping)THEN
      tmparray(13,0)=Procinfo_glob(13,i) ! nLinearElemsHalo
      tmparray(14,0)=Procinfo_glob(14,i) ! nCurvedElemsHalo
      tmparray(15,0)=Procinfo_glob(15,i) ! nBCElems
      tmparray(16,0)=Procinfo_glob(16,i) ! nBCElemsHalo
    END IF
    tmparray(nVars+1,0)=nNBProcs_glob(i)   ! nNBProcs
    DO j=1,nVars+1
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2
    END DO ! j
  END DO ! i
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(PartMPI%nProcs))
  ! output
  IF(DoRefMapping)THEN
    WRITE(ioUnit,'(A15,17(8X,F12.4))')'   MEAN        ',tmpreal(:,1)
    WRITE(ioUnit,'(A90,A90,A90,A90,A20)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '--------------------'
    WRITE(ioUnit,'(A15,17(8X,F12.4))')'   RMS         ',tmpreal(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------'
    WRITE(ioUnit,'(A15,17(8X,I12))')'   MIN         ',tmparray(:,3)
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------'
    WRITE(ioUnit,'(A15,17(8X,I12))')'   MAX         ',tmparray(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A90)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========================================================================================='
    DO i=0,PartMPI%nProcs-1
      WRITE(ioUnit,'(17(5X,I10))')i,Procinfo_glob(:,i),nNBProcs_glob(i)
      WRITE(ioUnit,'(A90,A90,A90,A90)')&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------'
    END DO! i
  ELSE
    WRITE(ioUnit,'(A15,13(8X,F12.4))')'   MEAN        ',tmpreal(:,1)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A15,13(8X,F12.4))')'   RMS         ',tmpreal(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A15,13(8X,I12))')'   MIN         ',tmparray(:,3)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A15,13(8X,I12))')'   MAX         ',tmparray(:,2)
    WRITE(ioUnit,'(A90,A90,A90,A30)')&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------------------------------------------------------------------',&
        '------------------------------'
    WRITE(ioUnit,'(A90,A90,A90,A10)')&
        '==========================================================================================',&
        '==========================================================================================',&
        '==========================================================================================',&
        '=========='
    DO i=0,PartMPI%nProcs-1
      WRITE(ioUnit,'(13(5X,I10))')i,Procinfo_glob(:,i),nNBProcs_glob(i)
      WRITE(ioUnit,'(A90,A90,A90,A10)')&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '------------------------------------------------------------------------------------------',&
          '----------'
    END DO! i
  END IF
  DEALLOCATE(tmparray,tmpreal)
  CLOSE(ioUnit)
END IF !MPIroot
DEALLOCATE(NBinfo_glob,nNBProcs_glob,ProcInfo_glob)

END SUBROUTINE WriteParticlePartitionInformation


SUBROUTINE WriteParticleMappingPartitionInformation(nPlanar,nBilinear,nCurved,nTotalBCElems)
!===================================================================================================================================
! write the particle partition information to file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,            ONLY:nSides,nElems
USE MOD_Particle_MPI_Vars,    ONLY:PartMPI
USE MOD_Particle_Mesh_Vars,   ONLY:nTotalSides,nTotalElems
USE MOD_LoadBalance_Vars,     ONLY:DoLoadBalance,nLoadBalanceSteps, writePartitionInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)         :: nPlanar,nBilinear,nCurved,nTotalBCElems
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE        :: nNBProcs_glob(:), ProcInfo_glob(:,:),NBInfo_glob(:,:), NBInfo(:), tmparray(:,:)
REAL,ALLOCATABLE           :: tmpreal(:,:)
INTEGER                    :: ProcInfo(8),nNBmax,i,j,ioUnit
CHARACTER(LEN=64)          :: filename
CHARACTER(LEN=4)           :: hilf
!===================================================================================================================================

IF(.NOT.WritePartitionInfo) RETURN

!output partitioning info
ProcInfo(1)=nElems
ProcInfo(2)=nSides
ProcInfo(3)=nTotalBCElems
ProcInfo(4)=nTotalElems-nElems
ProcInfo(5)=nTotalSides-nSides
ProcInfo(6)=nPlanar
ProcInfo(7)=nBilinear
ProcInfo(8)=nCurved
IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:PartMPI%nProcs-1))
  ALLOCATE(ProcInfo_glob(8,0:PartMPI%nProcs-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot
CALL MPI_GATHER(PartMPI%nMPINeighbors,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,PartMPI%COMM,iError)
CALL MPI_GATHER(ProcInfo,8,MPI_INTEGER,ProcInfo_glob,8,MPI_INTEGER,0,PartMPI%COMM,iError)
IF(MPIroot)THEN
  nNBmax=MAXVAL(nNBProcs_glob) !count, total number of columns in table
  ALLOCATE(NBinfo_glob(nNBmax,0:PartMPI%nProcs-1))
  NBinfo_glob=-77777
ELSE
  ALLOCATE(NBinfo_glob(1,1)) !dummy for debug
END IF
CALL MPI_BCAST(nNBmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
ALLOCATE(NBinfo(nNbmax))
NBinfo(1:PartMPI%nMPINeighbors)=PartMPI%MPINeighbor(1:PartMPi%nMPINeighbors)
CALL MPI_GATHER(NBinfo,nNBmax,MPI_INTEGER,NBinfo_glob,nNBmax,MPI_INTEGER,0,PartMPI%COMM,iError)
DEALLOCATE(NBinfo)
IF(MPIroot)THEN
  ioUnit=GETFREEUNIT()
  IF(DoLoadBalance)THEN
    WRITE( hilf,'(I4.4)') nLoadBalanceSteps
    filename='particlepartitionInfo-'//TRIM(hilf)//'.out'
  ELSE
    filename='particlepartitionInfo.out'
  END IF
  OPEN(UNIT=ioUnit,FILE=filename,STATUS='REPLACE')
  WRITE(ioUnit,*)'Particle Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',PartMPI%nProcs
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  WRITE(ioUnit,'(10(A15))')'Rank','nElems','nSides','nBCElems','nHaloElems','nHaloSides','nPlanar','nBilinear','nCurved','nNBProcs'
  WRITE(ioUnit,'(A90,A60)')&
      '==========================================================================================',&
      '============================================================'
  !statistics
  ALLOCATE(tmparray(9,0:3),tmpreal(9,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=0      !max
  tmparray(:,3)=HUGE(1)   !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i) ! nElems
    tmparray(2,0)=Procinfo_glob(2,i) ! nSides
    tmparray(3,0)=Procinfo_glob(3,i) ! nBCElems
    tmparray(4,0)=Procinfo_glob(4,i) ! nHaloElems
    tmparray(5,0)=Procinfo_glob(5,i) ! nHaloSides
    tmparray(6,0)=Procinfo_glob(6,i) ! nPlanarSides
    tmparray(7,0)=Procinfo_glob(7,i) ! nBilinearSides
    tmparray(8,0)=Procinfo_glob(8,i) ! nCurvedSides
    tmparray(9,0)=nNBProcs_glob(i)   ! nNBProcs
    DO j=1,9
      !mean
      tmparray(j,1)=tmparray(j,1)+tmparray(j,0)
      !max
      tmparray(j,2)=MAX(tmparray(j,2),tmparray(j,0))
      tmparray(j,3)=MIN(tmparray(j,3),tmparray(j,0))
    END DO !j
  END DO ! i
  tmpreal(:,1)=REAL(tmparray(:,1))/REAL(PartMPI%nProcs) !mean in REAL
  tmpreal(:,2)=0.   !RMS
  DO i=0,PartMPI%nProcs-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i)
    tmparray(2,0)=Procinfo_glob(2,i)
    tmparray(3,0)=Procinfo_glob(3,i)
    tmparray(4,0)=Procinfo_glob(4,i)
    tmparray(5,0)=Procinfo_glob(5,i)
    tmparray(6,0)=Procinfo_glob(6,i)
    tmparray(7,0)=Procinfo_glob(7,i)
    tmparray(8,0)=Procinfo_glob(8,i)
    tmparray(9,0)=nNBProcs_glob(i)
    DO j=1,9
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2
    END DO ! j
  END DO ! i
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(PartMPI%nProcs))
  WRITE(ioUnit,'(A15,9(5X,F10.2))')'   MEAN        ',tmpreal(:,1)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A15,9(5X,F10.2))')'   RMS         ',tmpreal(:,2)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A15,9(5X,I10))')'   MIN         ',tmparray(:,3)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A15,9(5X,I10))')'   MAX         ',tmparray(:,2)
  WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  WRITE(ioUnit,'(A90,A60)')&
      '==========================================================================================',&
      '============================================================'
  DO i=0,PartMPI%nProcs-1
    WRITE(ioUnit,'(10(5X,I10))')i,Procinfo_glob(:,i),nNBProcs_glob(i)
    WRITE(ioUnit,'(A90,A60)')&
      '------------------------------------------------------------------------------------------',&
      '------------------------------------------------------------'
  END DO! i
  DEALLOCATE(tmparray,tmpreal)
  CLOSE(ioUnit)
END IF !MPIroot
DEALLOCATE(NBinfo_glob,nNBProcs_glob,ProcInfo_glob)

END SUBROUTINE WriteParticleMappingPartitionInformation
#endif /*MPI*/

END MODULE MOD_Particle_MPI_Halo
