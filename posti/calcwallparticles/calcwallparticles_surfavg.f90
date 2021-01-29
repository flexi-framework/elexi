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
! Module for DSMC Sampling and Output
!===================================================================================================================================
MODULE MOD_CalcWallParticles_SurfAvg
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitCalcWallParticles_SurfAvg
  MODULE PROCEDURE InitCalcWallParticles_SurfAvg
END INTERFACE

INTERFACE CalcWallParticles_SurfAvg
  MODULE PROCEDURE CalcWallParticles_SurfAvg
END INTERFACE

INTERFACE WriteAvgSampleToHDF5
  MODULE PROCEDURE WriteAvgSampleToHDF5
END INTERFACE

INTERFACE FinalizeCalcWallParticles_SurfAvg
  MODULE PROCEDURE FinalizeCalcWallParticles_SurfAvg
END INTERFACE

PUBLIC :: InitCalcWallParticles_SurfAvg, CalcWallParticles_SurfAvg, FinalizeCalcWallParticles_SurfAvg, WriteAvgSampleToHDF5
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Init spatial averaging of surface values
!==================================================================================================================================
SUBROUTINE InitCalcWallParticles_SurfAvg()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY:SurfMesh
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Mesh_Vars      ,ONLY:nTotalSides,PartSideToElem,GEO
USE MOD_Posti_CalcWallParticles_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: iSide,SurfSideID,ElemID,LocSideID
INTEGER                                :: iMinSide,SurfMinSideID,ElemMinID,LocSideMinID
INTEGER                                :: minCounter
REAL                                   :: minDir
REAL                                   :: xSide(2),xMin(2)
INTEGER,DIMENSION(3)                   :: surfAvgVecDir
!===================================================================================================================================
ALLOCATE(SurfSideIDtoSurfAvgID(nTotalSides))
!ALLOCATE(SurfAvgIDtoSurfSideID(nTotalSides))
SurfSideIDtoSurfAvgID = -1
!SurfAvgIDtoSurfSideID = -1

! Check if averaging vector is Cartesian
IF(surfAvgDir(1).NE.0) THEN
    IF((surfAvgDir(2).NE.0).OR.(surfAvgDir(3).NE.0)) THEN
        CALL abort(&
        __STAMP__&
        ,'Averaging vector not in Cartesian direction!')
    END IF
    ! Set vector for minimum localization
    surfAvgVecDir(1) = 2
    surfAvgVecDir(2) = 3
    surfAvgVecDir(3) = 1
ELSEIF(surfAvgDir(2).NE.0) THEN
    IF((surfAvgDir(3).NE.0).OR.(surfAvgDir(1).NE.0)) THEN
        CALL abort(&
        __STAMP__&
        ,'Averaging vector not in Cartesian direction!')
    END IF
    ! Set vector for minimum localization
    surfAvgVecDir(1) = 1
    surfAvgVecDir(2) = 3
    surfAvgVecDir(3) = 2
ELSEIF(surfAvgDir(3).NE.0) THEN
    IF((surfAvgDir(1).NE.0).OR.(surfAvgDir(2).NE.0)) THEN
        CALL abort(&
        __STAMP__&
        ,'Averaging vector not in Cartesian direction!')
    END IF
    ! Set vector for minimum localization
    surfAvgVecDir(1) = 1
    surfAvgVecDir(2) = 2
    surfAvgVecDir(3) = 3
ELSE
    CALL abort(&
    __STAMP__&
    ,'No valid averaging vector found!')
END IF

! Loop over all sides to find the minimum in averaging direction
minDir     = HUGE(1.)
minCounter = 0

DO iSide=1,nTotalSides
    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
    ! Not a surfMesh side
    IF(SurfSideID.EQ.-1) CYCLE

    ElemID    = PartSideToElem(S2E_ELEM_ID,iSide)
    LocSideID = PartSideToElem(S2E_LOC_SIDE_ID,iSide)

    IF (MINVAL(GEO%NodeCoords(surfAvgVecDir(3),GEO%ElemSideNodeID(:,LocSideID,ElemID))).LT.minDir) THEN
        minDir = MINVAL(GEO%NodeCoords(surfAvgVecDir(3),GEO%ElemSideNodeID(:,LocSideID,ElemID)))
    END IF
END DO

! Loop over all sides to assign to corresponding minimum side
DO iSide=1,nTotalSides
    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
    ! Not a surfMesh side
    IF(SurfSideID.EQ.-1) CYCLE

    ElemID = PartSideToElem(S2E_ELEM_ID,iSide)
    LocSideID = PartSideToElem(S2E_LOC_SIDE_ID,iSide)

    xSide(1) = MINVAL(GEO%NodeCoords(surfAvgVecDir(1),GEO%ElemSideNodeID(:,LocSideID,ElemID)))
    xSide(2) = MINVAL(GEO%NodeCoords(surfAvgVecDir(2),GEO%ElemSideNodeID(:,LocSideID,ElemID)))

    IF (.NOT.ALMOSTEQUAL(MINVAL(GEO%NodeCoords(surfAvgVecDir(3),GEO%ElemSideNodeID(:,LocSideID,ElemID))),minDir)) THEN
        ! Loop till we found the corresponding minimum side
        DO iMinSide=1,nTotalSides
            SurfMinSideID=SurfMesh%SideIDToSurfID(iMinSide)
            ! Not a surfMesh side
            IF(SurfMinSideID.EQ.-1) CYCLE
            ! Same side always matches but we do not want this side
            IF(SurfMinSideID.EQ.SurfSideID) CYCLE

            ElemMinID    = PartSideToElem(S2E_ELEM_ID,iMinSide)
            LocSideMinID = PartSideToElem(S2E_LOC_SIDE_ID,iMinSide)

            xMin(1) = MINVAL(GEO%NodeCoords(surfAvgVecDir(1),GEO%ElemSideNodeID(:,LocSideMinID,ElemMinID)))
            xMin(2) = MINVAL(GEO%NodeCoords(surfAvgVecDir(2),GEO%ElemSideNodeID(:,LocSideMinID,ElemMinID)))

            IF(ALMOSTEQUAL(xSide(1),xMin(1)).AND.ALMOSTEQUAL(xSide(2),xMin(2))) THEN
                ! Make sure we actually found the minimum side
                IF(ALMOSTEQUAL(MINVAL(GEO%NodeCoords(surfAvgVecDir(3),GEO%ElemSideNodeID(:,LocSideMinID,ElemMinID))),minDir)) THEN
                    SurfSideIDtoSurfAvgID(SurfSideID)       = SurfMinSideID
!                    SurfAvgIDtoSurfSideID(SurfMinSideID)    = SurfSideID
                    ! Increment counter
                    minCounter                              = minCounter + 1
                    EXIT
                END IF
            END IF

        END DO
    END IF
END DO

SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') '| Found corresponding minimum side for ',minCounter,' of ',SurfMesh%nSides,' sides.'

END SUBROUTINE InitCalcWallParticles_SurfAvg


!==================================================================================================================================
!> Perform spatial averaging of surface values
!==================================================================================================================================
SUBROUTINE CalcWallParticles_SurfAvg()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Preproc                 ,ONLY:PP_N
USE MOD_Mesh_Vars               ,ONLY:Face_xGP
USE MOD_Particle_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY:SurfMesh,nSurfSample
USE MOD_DSMC_Vars               ,ONLY:MacroSurfaceVal
USE MOD_Particle_Vars           ,ONLY:nSpecies
USE MOD_Particle_Erosion_Vars
USE MOD_Particle_Mesh_Vars      ,ONLY:nTotalSides
USE MOD_Posti_CalcWallParticles_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: p,q,iSpec,nShift,nShiftRHS
INTEGER                                :: ptemp,qtemp,pcount,qcount
REAL,DIMENSION(3)                      :: surfSideVec
INTEGER                                :: surfSideAvgDir            !< Flag if the averaging direction is xi or eta
INTEGER                                :: iSide,SurfSideID
INTEGER                                :: SurfMinSideID
!===================================================================================================================================

! Loop over all sides to find sides we must move
DO iSide=1,nTotalSides
    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
    ! Not a surfMesh side
    IF(SurfSideID.EQ.-1) CYCLE
    ! Already a minimum side
    IF(SurfSideIDtoSurfAvgID(SurfSideID).EQ.-1) CYCLE

    SurfMinSideID = SurfSideIDtoSurfAvgID(SurfSideID)
    ! Sanity check SurfMinSideID
    IF (SurfMinSideID.EQ.SurfSideID) THEN
    CALL Abort(&
      __STAMP__,&
      'Error finding corresponding minimum side. This should not happen.')
    END IF

!---- Only one species. Only total values necessary
!===================================================================================================================================
    DO q=1,nSurfSample
        DO p=1,nSurfSample
            !---- 1. - .. / Impact Counter
            MacroSurfaceVal(1,p,q,SurfMinSideID) = MacroSurfaceVal(1,p,q,SurfMinSideID) + MacroSurfaceVal(1,p,q,SurfSideID)
!-----------! Ignore faces with no impacts. We can only be sure after we check both faces !----------------------------------------!
            IF (MacroSurfaceVal(1,p,q,SurfMinSideID).EQ.0) CYCLE
            !---- 2. - .. / Impact Counter per AREA
            MacroSurfaceVal(2,p,q,SurfMinSideID) = MacroSurfaceVal(2,p,q,SurfMinSideID) + MacroSurfaceVal(2,p,q,SurfSideID)
            !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
            MacroSurfaceVal(3,p,q,SurfMinSideID) = MacroSurfaceVal(3,p,q,SurfMinSideID) * (MacroSurfaceVal(1,p,q,SurfMinSideID) - &
                                                   MacroSurfaceVal(1,p,q,SurfSideID))   /  MacroSurfaceVal(1,p,q,SurfMinSideID) + &
                                                   MacroSurfaceVal(3,p,q,SurfSideID)    *  MacroSurfaceVal(1,p,q,SurfSideID)      &
                                                                                        /  MacroSurfaceVal(1,p,q,SurfMinSideID)
            MacroSurfaceVal(4,p,q,SurfMinSideID) = MIN(MacroSurfaceVal(4,p,q,SurfMinSideID),MacroSurfaceVal(4,p,q,SurfSideID))
            MacroSurfaceVal(5,p,q,SurfMinSideID) = MAX(MacroSurfaceVal(4,p,q,SurfMinSideID),MacroSurfaceVal(5,p,q,SurfSideID))
            MacroSurfaceVal(6,p,q,SurfMinSideID) = 0.                         ! Can't reconstruct variance. Erase to avoid confusion
            !---- 7. - 10 / Impact angle (mean, min, max, variance)
            MacroSurfaceVal(7,p,q,SurfMinSideID) = MacroSurfaceVal(7,p,q,SurfMinSideID) * (MacroSurfaceVal(1,p,q,SurfMinSideID) - &
                                                   MacroSurfaceVal(1,p,q,SurfSideID))   /  MacroSurfaceVal(1,p,q,SurfMinSideID) + &
                                                   MacroSurfaceVal(7,p,q,SurfSideID)    *  MacroSurfaceVal(1,p,q,SurfSideID)      &
                                                                                        /  MacroSurfaceVal(1,p,q,SurfMinSideID)
            MacroSurfaceVal(8,p,q,SurfMinSideID) = MIN(MacroSurfaceVal(8,p,q,SurfMinSideID),MacroSurfaceVal(8,p,q,SurfSideID))
            MacroSurfaceVal(9,p,q,SurfMinSideID) = MAX(MacroSurfaceVal(9,p,q,SurfMinSideID),MacroSurfaceVal(9,p,q,SurfSideID))
            MacroSurfaceVal(10,p,q,SurfMinSideID)= 0.                         ! Can't reconstruct variance. Erase to avoid confusion
            !---- 11 - 13 / Sampling Current Forces at walls
            MacroSurfaceVal(11:13,p,q,SurfMinSideID)=0.                       ! Erase for now
            !---- 14 - 16 / Sampling Average Forces at walls
            MacroSurfaceVal(14:16,p,q,SurfMinSideID)=0.                       ! Erase for now
        END DO
    END DO
!---- Multiple species. All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
    IF (nSpecies.GT.1) THEN
        DO q=1,nSurfSample
            DO p=1,nSurfSample
                DO iSpec=1,nSpecies
                    nShift    = iSpec * (nErosionVars-1)
                    nShiftRHS = iSpec * nErosionVars
                    !---- 1. - .. / Impact Counter
                    MacroSurfaceVal(1+nShift,p,q,SurfMinSideID) = MacroSurfaceVal(1+nShift,p,q,SurfMinSideID) + &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfSideID)
!-------------------! Ignore faces with no impacts. We can only be sure after we check both faces !--------------------------------!
                    IF (MacroSurfaceVal(1+nShift,p,q,SurfMinSideID).EQ.0) CYCLE

                    !---- 2. - .. / Impact Counter per AREA
                    MacroSurfaceVal(2+nShift,p,q,SurfMinSideID) = MacroSurfaceVal(2+nShift,p,q,SurfMinSideID) + &
                                                                  MacroSurfaceVal(2+nShift,p,q,SurfSideID)
                    !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
                    MacroSurfaceVal(3+nShift,p,q,SurfMinSideID) = MacroSurfaceVal(3+nShift,p,q,SurfMinSideID) * &
                                                                 (MacroSurfaceVal(1+nShift,p,q,SurfMinSideID) - &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfSideID))   / &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfMinSideID) + &
                                                                  MacroSurfaceVal(3+nShift,p,q,SurfSideID)    * &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfSideID)    / &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfMinSideID)
                    MacroSurfaceVal(4+nShift,p,q,SurfMinSideID) = MIN(MacroSurfaceVal(4+nShift,p,q,SurfMinSideID), &
                                                                      MacroSurfaceVal(4+nShift,p,q,SurfSideID))
                    MacroSurfaceVal(5+nShift,p,q,SurfMinSideID) = MAX(MacroSurfaceVal(4+nShift,p,q,SurfMinSideID), &
                                                                      MacroSurfaceVal(5+nShift,p,q,SurfSideID))
                    MacroSurfaceVal(6+nShift,p,q,SurfMinSideID) = 0.          ! Can't reconstruct variance. Erase to avoid confusion
                    !---- 7. - 10 / Impact angle (mean, min, max, variance)
                    MacroSurfaceVal(7+nShift,p,q,SurfMinSideID) = MacroSurfaceVal(7+nShift,p,q,SurfMinSideID) * &
                                                                 (MacroSurfaceVal(1+nShift,p,q,SurfMinSideID) - &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfSideID))   / &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfMinSideID) + &
                                                                  MacroSurfaceVal(7+nShift,p,q,SurfSideID)    * &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfSideID)    / &
                                                                  MacroSurfaceVal(1+nShift,p,q,SurfMinSideID)
                    MacroSurfaceVal(8+nShift,p,q,SurfMinSideID) = MIN(MacroSurfaceVal(8+nShift,p,q,SurfMinSideID), &
                                                                      MacroSurfaceVal(8+nShift,p,q,SurfSideID))
                    MacroSurfaceVal(9+nShift,p,q,SurfMinSideID) = MAX(MacroSurfaceVal(9+nShift,p,q,SurfMinSideID), &
                                                                      MacroSurfaceVal(9+nShift,p,q,SurfSideID))
                    MacroSurfaceVal(10+nShift,p,q,SurfMinSideID)= 0.          ! Can't reconstruct variance. Erase to avoid confusion
                    !---- 11 - 13 / Sampling Current Forces at walls
                    MacroSurfaceVal(11+nShift:13+nShift,p,q,SurfMinSideID)=0.                       ! Erase for now
                    !---- 14 - 16 / Sampling Average Forces at walls
                    MacroSurfaceVal(14+nShift:16+nShift,p,q,SurfMinSideID)=0.                       ! Erase for now
                END DO
            END DO
        END DO
    END IF
END DO

! Write values back so we have them on the entire surface
DO iSide=1,nTotalSides
    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
    ! Not a surfMesh side
    IF(SurfSideID.EQ.-1) CYCLE
    ! Already a minimum side
    IF(SurfSideIDtoSurfAvgID(SurfSideID).EQ.-1) CYCLE

    SurfMinSideID = SurfSideIDtoSurfAvgID(SurfSideID)
    ! Sanity check SurfMinSideID
    IF (SurfMinSideID.EQ.SurfSideID) THEN
    CALL Abort(&
      __STAMP__,&
      'Error finding corresponding minimum side. This should not happen.')
    END IF
!===================================================================================================================================
    DO q=1,nSurfSample
        DO p=1,nSurfSample
            MacroSurfaceVal(:,p,q,SurfSideID) = MacroSurfaceVal(:,p,q,SurfMinSideID)
        END DO
    END DO
END DO

! Initialize counters
pcount = 0
qcount = 0

! Also average in xi direction if nSurfSample > 1
IF (nSurfSample.GT.1) THEN
    ! Loop over all surfaces
    DO iSide=1,nTotalSides
        SurfSideID=SurfMesh%SideIDToSurfID(iSide)
        ! Not a surfMesh side
        IF(SurfSideID.EQ.-1) CYCLE
        ! Already a minimum side

        DO q=1,nSurfSample
            DO p=1,nSurfSample
                !< Find averaging direction by first testing xi-direction
                surfSideVec = Face_xGP(:,PP_N,0,0,iSide) - Face_xGP(:,0,0,0,iSide)

                IF(.NOT.ALMOSTZERO(DOT_PRODUCT(surfSideVec,surfAvgDir))) THEN
                    ! xi-direction works
                    ptemp          = 1
                    qtemp          = q
                    surfSideAvgDir = 1
                    pcount         = pcount + 1
                    ! Do not add to itself
                    IF (p.EQ.1) CYCLE
                ELSE
                    surfSideVec = Face_xGP(:,0,PP_NZ,0,iSide) - Face_xGP(:,0,0,0,iSide)

                    IF(.NOT.ALMOSTZERO(DOT_PRODUCT(surfSideVec,surfAvgDir))) THEN
                        ! eta-direction works
                        ptemp          = p
                        qtemp          = 1
                        surfSideAvgDir = 2
                        qcount         = qcount + 1
                        ! Do not add to itself
                        IF (q.EQ.1) CYCLE
                    ELSE
                        ! We could not find the averaging direction in the xi,eta plane
                        CALL Abort(&
                        __STAMP__,&
                        'Invalid averaging direction.')
                    END IF
                END IF

!---- Only one species. Only total values necessary
!===================================================================================================================================
        !---- 1. - .. / Impact Counter
        MacroSurfaceVal(1,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(1,ptemp,qtemp,SurfSideID) + MacroSurfaceVal(1,p,q,SurfSideID)
!-------! Ignore faces with no impacts. We can only be sure after we check both faces !----------------------------------------!
        IF (MacroSurfaceVal(1,ptemp,qtemp,SurfSideID).EQ.0) CYCLE
        !---- 2. - .. / Impact Counter per AREA
        MacroSurfaceVal(2,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(2,ptemp,qtemp,SurfSideID) +                   &
                                                    MacroSurfaceVal(2,p,q,SurfSideID)
        !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
        MacroSurfaceVal(3,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(3,ptemp,qtemp,SurfSideID)  * &
                                                   (MacroSurfaceVal(1,ptemp,qtemp,SurfSideID)  - &
                                                    MacroSurfaceVal(1,p,q,SurfSideID))         / &
                                                    MacroSurfaceVal(1,ptemp,qtemp,SurfSideID)  + &
                                                    MacroSurfaceVal(3,p,q,SurfSideID)          * &
                                                    MacroSurfaceVal(1,p,q,SurfSideID)          / &
                                                    MacroSurfaceVal(1,ptemp,qtemp,SurfSideID)
        MacroSurfaceVal(4,ptemp,qtemp,SurfSideID) = MIN(MacroSurfaceVal(4,ptemp,qtemp,SurfSideID),MacroSurfaceVal(4,p,q,SurfSideID))
        MacroSurfaceVal(5,ptemp,qtemp,SurfSideID) = MAX(MacroSurfaceVal(4,ptemp,qtemp,SurfSideID),MacroSurfaceVal(5,p,q,SurfSideID))
        MacroSurfaceVal(6,ptemp,qtemp,SurfSideID) = 0.                        ! Can't reconstruct variance. Erase to avoid confusion
        !---- 7. - 10 / Impact angle (mean, min, max, variance)
        MacroSurfaceVal(7,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(7,ptemp,qtemp,SurfSideID)  * &
                                                   (MacroSurfaceVal(1,ptemp,qtemp,SurfSideID)  - &
                                                    MacroSurfaceVal(1,p,q,SurfSideID))         / &
                                                    MacroSurfaceVal(1,ptemp,qtemp,SurfSideID)  + &
                                                    MacroSurfaceVal(7,p,q,SurfSideID)          * &
                                                    MacroSurfaceVal(1,p,q,SurfSideID)          / &
                                                    MacroSurfaceVal(1,ptemp,qtemp,SurfSideID)
        MacroSurfaceVal(8,ptemp,qtemp,SurfSideID) = MIN(MacroSurfaceVal(8,ptemp,qtemp,SurfSideID),MacroSurfaceVal(8,p,q,SurfSideID))
        MacroSurfaceVal(9,ptemp,qtemp,SurfSideID) = MAX(MacroSurfaceVal(9,ptemp,qtemp,SurfSideID),MacroSurfaceVal(9,p,q,SurfSideID))
        MacroSurfaceVal(10,ptemp,qtemp,SurfSideID)= 0.                        ! Can't reconstruct variance. Erase to avoid confusion
        !---- 11 - 13 / Sampling Current Forces at walls
        MacroSurfaceVal(11:13,ptemp,qtemp,SurfSideID)=0.                       ! Erase for now
        !---- 14 - 16 / Sampling Average Forces at walls
        MacroSurfaceVal(14:16,ptemp,qtemp,SurfSideID)=0.                       ! Erase for now
!---- Multiple species. All Variables are saved DOUBLE. First Total, then per SPECIES
!===================================================================================================================================
        DO iSpec=1,nSpecies
            nShift    = iSpec * (nErosionVars-1)
            nShiftRHS = iSpec *  nErosionVars
            !---- 1. - .. / Impact Counter
            MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID) + &
                                                               MacroSurfaceVal(1+nShift,p,q,SurfSideID)
!-----------! Ignore faces with no impacts. We can only be sure after we check both faces !--------------------------------!
            IF (MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID).EQ.0) CYCLE

            !---- 2. - .. / Impact Counter per AREA
            MacroSurfaceVal(2+nShift,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(2+nShift,ptemp,qtemp,SurfSideID) + &
                                                               MacroSurfaceVal(2+nShift,p,q,SurfSideID)
            !---- 3. - 6. / Kinetic energy on impact (mean, min, max, variance)
            MacroSurfaceVal(3+nShift,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(3+nShift,ptemp,qtemp,SurfSideID)  * &
                                                              (MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID)  - &
                                                               MacroSurfaceVal(1+nShift,p,q,SurfSideID))         / &
                                                               MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID)  + &
                                                               MacroSurfaceVal(3+nShift,p,q,SurfSideID)          * &
                                                               MacroSurfaceVal(1+nShift,p,q,SurfSideID)          / &
                                                               MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID)
            MacroSurfaceVal(4+nShift,ptemp,qtemp,SurfSideID) = MIN(MacroSurfaceVal(4+nShift,ptemp,qtemp,SurfSideID), &
                                                                   MacroSurfaceVal(4+nShift,p,q,SurfSideID))
            MacroSurfaceVal(5+nShift,ptemp,qtemp,SurfSideID) = MAX(MacroSurfaceVal(4+nShift,ptemp,qtemp,SurfSideID), &
                                                                   MacroSurfaceVal(5+nShift,p,q,SurfSideID))
            MacroSurfaceVal(6+nShift,ptemp,qtemp,SurfSideID) = 0.             ! Can't reconstruct variance. Erase to avoid confusion
            !---- 7. - 10 / Impact angle (mean, min, max, variance)
            MacroSurfaceVal(7+nShift,ptemp,qtemp,SurfSideID) = MacroSurfaceVal(7+nShift,ptemp,qtemp,SurfSideID)  * &
                                                              (MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID)  - &
                                                               MacroSurfaceVal(1+nShift,p,q,SurfSideID))         / &
                                                               MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID)  + &
                                                               MacroSurfaceVal(7+nShift,p,q,SurfSideID)          * &
                                                               MacroSurfaceVal(1+nShift,p,q,SurfSideID)          / &
                                                               MacroSurfaceVal(1+nShift,ptemp,qtemp,SurfSideID)
            MacroSurfaceVal(8+nShift,ptemp,qtemp,SurfSideID) = MIN(MacroSurfaceVal(8+nShift,ptemp,qtemp,SurfSideID), &
                                                                   MacroSurfaceVal(8+nShift,p,q,SurfSideID))
            MacroSurfaceVal(9+nShift,ptemp,qtemp,SurfSideID) = MAX(MacroSurfaceVal(9+nShift,ptemp,qtemp,SurfSideID), &
                                                                   MacroSurfaceVal(9+nShift,p,q,SurfSideID))
            MacroSurfaceVal(10+nShift,ptemp,qtemp,SurfSideID)= 0.             ! Can't reconstruct variance. Erase to avoid confusion
            !---- 11 - 13 / Sampling Current Forces at walls
            MacroSurfaceVal(11+nShift:13+nShift,ptemp,qtemp,SurfSideID)=0.                       ! Erase for now
            !---- 14 - 16 / Sampling Average Forces at walls
            MacroSurfaceVal(14+nShift:16+nShift,ptemp,qtemp,SurfSideID)=0.                       ! Erase for now
                END DO ! iSpec
            END Do ! p
        END DO ! q

!===================================================================================================================================
! Write values back so we have them on the entire surface
        IF (surfSideAvgDir.EQ.1) THEN
            DO p=2,nSurfSample
                    MacroSurfaceVal(:,p,:,SurfSideID) = MacroSurfaceVal(:,1,:,SurfSideID)
            END DO
        ELSE ! surfSideAvgDir = 2
            DO q=2,nSurfSample
                MacroSurfaceVal(:,:,q,SurfSideID) = MacroSurfaceVal(:,:,1,SurfSideID)
            END DO
        END IF
    END DO ! iSide
END IF

SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') '| Found corresponding xi/eta direction for ', &
                                     pcount/nSurfSample**2,'/',qcount/nSurfSample**2,' sides'

END SUBROUTINE CalcWallParticles_SurfAvg

SUBROUTINE WriteAvgSampleToHDF5(MeshFileName,OutputTime)
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
USE MOD_Particle_Boundary_Vars,     ONLY:nSurfSample,SurfMesh,offSetSurfSide
USE MOD_Particle_Erosion_Vars
USE MOD_DSMC_Vars,                  ONLY:MacroSurfaceVal,MacroSurfaceSpecVal,CollisMode
USE MOD_Particle_HDF5_output,       ONLY:WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString,Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255)                  :: SpecID
CHARACTER(LEN=255)                  :: tmp255
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: nVar2D, nVar2D_Spec, nVar2D_Total, nVarCount, iSpec, nShift
!===================================================================================================================================

IF(SurfCOMM%MPIOutputRoot)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' WRITE EROSION AVERAGE STATE TO HDF5 FILE...'
END IF

FileName=TIMESTAMP(TRIM(ProjectName)//'_ErosionAvgState',OutputTime)
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
  CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode',1,IntegerScalar=CollisMode)
  tmp255=TRIM(MeshFileName)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/tmp255/))
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

! Finaly deallocate Macrosurfaces we couldn't deallocate earlier
SDEALLOCATE(MacroSurfaceVal)
SDEALLOCATE(MacroSurfaceSpecVal)

END SUBROUTINE WriteAvgSampleToHDF5


SUBROUTINE FinalizeCalcWallParticles_SurfAvg()
!===================================================================================================================================
! deallocate everything
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Posti_CalcWallParticles_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(SurfSideIDtoSurfAvgID)
!SDEALLOCATE(SurfAvgIDtoSurfSideID)

END SUBROUTINE FinalizeCalcWallParticles_SurfAvg

END MODULE MOD_CalcWallParticles_SurfAvg
