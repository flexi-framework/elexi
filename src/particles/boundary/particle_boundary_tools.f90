!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
! Tools used for boundary interactions
!===================================================================================================================================
MODULE MOD_Particle_Boundary_Tools
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE MarkAuxBCElems
  MODULE PROCEDURE MarkAuxBCElems
END INTERFACE

PUBLIC :: MarkAuxBCElems
!===================================================================================================================================

CONTAINS

SUBROUTINE MarkAuxBCElems()
!===================================================================================================================================
! check if auxBCs are inside BoundingBox of Elems
! -- plane: use plane equation f=a1*x+a2*y+a3*z+a4=0 and insert corresponding intervals of box -> fmin and fmax
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars                          ,ONLY:nElems
USE MOD_Particle_Globals
USE MOD_Particle_Mesh_Vars                 ,ONLY:ElemHasAuxBCs
USE MOD_Particle_Mesh_Vars                 ,ONLY: BoundsOfElem_Shared
USE MOD_Particle_Boundary_Vars             ,ONLY:nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone!,AuxBC_parabol
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,iAuxBC,icoord,dir(3),positiontype,positiontype_tmp
REAL                     :: r_vec(3),n_vec(3),fmin,fmax,radius,BoundsBC(1:2,1:3)
REAL                     :: lmin,lmax,deltamin,deltamax,origin(2),halfangle
LOGICAL                  :: cartesian, backwards
!===================================================================================================================================

ALLOCATE(ElemHasAuxBCs(1:nElems , 1:nAuxBCs))
ElemHasAuxBCs=.FALSE.

DO iAuxBC=1,nAuxBCs
  SELECT CASE (TRIM(AuxBCType(iAuxBC)))

    CASE ('plane')
      r_vec  = AuxBC_plane(AuxBCMap(iAuxBC))%r_vec
      n_vec  = AuxBC_plane(AuxBCMap(iAuxBC))%n_vec
      radius = AuxBC_plane(AuxBCMap(iAuxBC))%radius

      ! loop over all  elements
      DO iElem = 1,nElems
        ! 1-2: Min, Max value; 1-3: x,y,z
        ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) )

        fmin = -DOT_PRODUCT(r_vec,n_vec)
        fmax = fmin
        DO icoord=1,3
          IF (n_vec(icoord).GE.0) THEN
            fmin = fmin + n_vec(icoord)*Bounds(1,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(2,icoord)
          ELSE
            fmin = fmin + n_vec(icoord)*Bounds(2,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(1,icoord)
          END IF
        END DO

        ! plane intersects the box!
        IF ((fmin.LE.0 .AND. fmax.GT.0).OR.(fmin.LT.0 .AND. fmax.GE.0)) THEN
          !radius check needs to be implemented (compute intersection polygon and minimum radii): would sort out further elements!!!
          !quick, conservative solution: calculate bounding box of disc in space and compare with bb of element
          ElemHasAuxBCs(iElem,iAuxBC) = .TRUE.

          IF (radius .LT. 0.5*HUGE(radius)) THEN !huge was default
            BoundsBC(1,1:3) = r_vec - radius * SQRT(1.-(n_vec*n_vec))
            BoundsBC(2,1:3) = r_vec + radius * SQRT(1.-(n_vec*n_vec))
            DO icoord=1,3
              IF ( BoundsBC(2,icoord).LT.Bounds(1,icoord) .OR. BoundsBC(1,icoord).GT.Bounds(2,icoord) ) THEN
                ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
                EXIT
              END IF
            END DO
          END IF
        ! plane does not intersect the box!
        ELSE IF ((fmin.LT.0 .AND. fmax.LT.0).OR.(fmin.GT.0 .AND. fmax.GT.0)) THEN
          ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
        ! e.g. if elem has zero volume...
        ELSE
          CALL abort(__STAMP__,'Error in MarkAuxBCElems for AuxBC:',iAuxBC)
        END IF

        END ASSOCIATE
      END DO

    CASE ('cylinder','cone')
      IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
        r_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%r_vec
        n_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%axis
        radius=AuxBC_cylinder(AuxBCMap(iAuxBC))%radius
        lmin=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin
        lmax=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmax
      ELSE !cone
        r_vec=AuxBC_cone(AuxBCMap(iAuxBC))%r_vec
        n_vec=AuxBC_cone(AuxBCMap(iAuxBC))%axis
        halfangle=AuxBC_cone(AuxBCMap(iAuxBC))%halfangle
        lmin=AuxBC_cone(AuxBCMap(iAuxBC))%lmin
        lmax=AuxBC_cone(AuxBCMap(iAuxBC))%lmax
      END IF

      cartesian = .TRUE.
      backwards = .FALSE.
      IF (ABS(n_vec(1)).EQ.1.) THEN
        dir(1)=1
        dir(2)=2
        dir(3)=3
        IF (n_vec(1).LT.0.) backwards = .TRUE.
      ELSE IF (ABS(n_vec(2)).EQ.1.) THEN
        dir(1)=2
        dir(2)=3
        dir(3)=1
        IF (n_vec(2).LT.0.) backwards = .TRUE.
      ELSE IF (ABS(n_vec(3)).EQ.1.) THEN
        dir(1)=3
        dir(2)=1
        dir(3)=2
        IF (n_vec(3).LT.0.) backwards = .TRUE.
      ELSE
        cartesian=.FALSE.
        SWRITE(*,*) 'WARNING in MarkAuxBCElems: all Elems are set to ElemHasAuxBCs=.TRUE. for AuxBC:',iAuxBC
        ElemHasAuxBCs(:,iAuxBC) = .TRUE. !actual intersection with box check to-be implemented!!!
      END IF

      IF (cartesian) THEN
        IF (backwards) THEN
          deltamin = -lmax
          deltamax = -lmin
        ELSE
          deltamin = lmin
          deltamax = lmax
        END IF

        origin(1) = r_vec(dir(2))
        origin(2) = r_vec(dir(3))

        ! loop over all  elements
        DO iElem=1,nElems
          ! 1-2: Min, Max value; 1-3: x,y,z
          ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) )

          ! check for lmin and lmax
          IF ( r_vec(dir(1))+deltamax.LT.Bounds(1,dir(1)) .OR. r_vec(dir(1))+deltamin.GT.Bounds(2,dir(1)) ) THEN
            ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
          ELSE !between lmin and lmax
            IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
            ELSE !cone
              !local minimum radius
              IF (backwards) THEN
                radius = MAX(-Bounds(2,dir(1))+r_vec(dir(1)),lmin)*TAN(halfangle)
              ELSE
                radius = MAX(Bounds(1,dir(1))-r_vec(dir(1)),lmin)*TAN(halfangle)
              END IF
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype_tmp)
              !local maximum radius
              IF (backwards) THEN
                radius = MIN(-Bounds(1,dir(1))+r_vec(dir(1)),lmax)*TAN(halfangle)
              ELSE
                radius = MIN(Bounds(2,dir(1))-r_vec(dir(1)),lmax)*TAN(halfangle)
              END IF
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
              !if both are type 0 or both are type 1 than the "total" type is not 2:
              IF ( .NOT.(positiontype_tmp.EQ.0 .AND. positiontype.EQ.0) &
                .AND. .NOT.(positiontype_tmp.EQ.1 .AND. positiontype.EQ.1) ) THEN
                positiontype=2
              END IF
            END IF
            IF (positiontype.EQ.2) THEN
              ElemHasAuxBCs(iElem,iAuxBC)=.TRUE.
            ELSE
              ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
            END IF
          END IF !check for lmin and lmax

          END ASSOCIATE
        END DO !iElem
      END IF !cartesian

    CASE('parabol')
      ! to be implemented!!!
      ElemHasAuxBCs(:,iAuxBC) = .TRUE.

    CASE DEFAULT
      SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
      CALL abort(__STAMP__,'AuxBC does not exist')
  END SELECT
END DO

END SUBROUTINE MarkAuxBCElems


SUBROUTINE CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
!===================================================================================================================================
! checks how a cartesian bb is located with regard to a radius with cartesian axis (dir is cartesian axis and origin in orth. dirs)
!- positiontype=0 : complete bb is inside of radius
!- positiontype=1 : complete bb is outside of radius
!- positiontype=2 : bb is partly inside of radius
! (based on "check where the sides are located relative to rmax" in particle_emission for SimpleRadialVeloFit)
!===================================================================================================================================
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)           :: Bounds(1:2,1:3), origin(2), radius
INTEGER,INTENT(IN)        :: dir(3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)       :: positiontype
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iDir1, iDir2, iDir3, iPoint
REAL                      :: BoundingBox(1:3,1:8), point(2), pointRadius
LOGICAL                   :: done, insideBound
!===================================================================================================================================
!-- convert minmax-values to bb-points
DO iDir1 = 0,1
  DO iDir2 = 0,1
      DO iDir3 = 0,1
        BoundingBox(1,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir1+1,1)
        BoundingBox(2,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir2+1,2)
        BoundingBox(3,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir3+1,3)
      END DO
  END DO
END DO

!-- check where the points are located relative to radius
done = .FALSE.
DO iDir1 = 0,1
  IF(done) EXIT

  DO iDir2  =0,1
    IF(done) EXIT

    DO iDir3 = 0,1
      !-- coords orth. to axis of point:
      iPoint=iDir1*4 + iDir2*2 + iDir3+1
      point(1) = BoundingBox(dir(2),iPoint)-origin(1)
      point(2) = BoundingBox(dir(3),iPoint)-origin(2)
      pointRadius = SQRT( (point(1))**2+(point(2))**2 )

      IF (iPoint.EQ.1) THEN
        IF (pointRadius.LE.radius) THEN
          insideBound=.TRUE.
        ELSE !outside
          insideBound=.FALSE.
        END IF !in-/outside?

      ! iPoint.GT.1: type must be 2 if state of point if different from last point
      ELSE
        IF (pointRadius.LE.radius) THEN
          ! different from last point
          IF (.NOT.insideBound) THEN
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        ! outside
        ELSE
          ! different from last point
          IF (insideBound) THEN
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        END IF !in-/outside?
      END IF !iPoint.EQ.1
    END DO !iDir3
  END DO !iDir2
END DO !iDir1

IF (.NOT.done) THEN
  IF (insideBound) THEN
    positiontype=0
  ELSE
    ! all points are outside of radius, but when radius is smaller than box, it can intersect it:
    IF ( origin(1) + radius .GE. Bounds(1,dir(2)) .AND. &
         origin(1) - radius .LE. Bounds(2,dir(2)) .AND. &
         origin(2) + radius .GE. Bounds(1,dir(3)) .AND. &
         origin(2) - radius .LE. Bounds(2,dir(3)) ) THEN !circle completely or partly inside box
      positiontype=2
    ! points are really outside
    ELSE
      positiontype=1
    END IF
  END IF
END IF

END SUBROUTINE CheckBoundsWithCartRadius

END MODULE MOD_Particle_Boundary_Tools
