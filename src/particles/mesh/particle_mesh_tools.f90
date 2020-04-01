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

MODULE MOD_Particle_Mesh_Tools
!===================================================================================================================================
! Contains subroutines for mappings of particle mesh
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitGetGlobalElemID
  MODULE PROCEDURE InitGetGlobalElemID
END INTERFACE

INTERFACE InitGetCNElemID
  MODULE PROCEDURE InitGetCNElemID
END INTERFACE

INTERFACE GetGlobalElemID
  PROCEDURE GetGlobalElemID
END INTERFACE

INTERFACE GetCNElemID
  PROCEDURE GetCNElemID
END INTERFACE

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetGlobalElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalElemIDInterface),POINTER :: GetGlobalElemID !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetCNElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetCNElemIDInterface)    ,POINTER :: GetCNElemID     !< pointer defining the mapping: global element ID -> compute-node element ID

PUBLIC :: InitGetGlobalElemID
PUBLIC :: InitGetCNElemID
PUBLIC :: GetGlobalElemID
PUBLIC :: GetCNElemID
!===================================================================================================================================

CONTAINS


!==================================================================================================================================!
!> Initialize GetGlobalElemID function (mapping of compute-node element ID to global element ID)
!==================================================================================================================================!
SUBROUTINE InitGetGlobalElemID()
! MODULES
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetGlobalElemID => GetGlobalElemID_iElem
ELSE
  GetGlobalElemID => GetGlobalElemID_fromTotalElem
END IF
#else
GetGlobalElemID => GetGlobalElemID_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalElemID_fromTotalElem(1)
#endif
dummy=GetGlobalElemID_iElem(1)
END SUBROUTINE InitGetGlobalElemID


!==================================================================================================================================!
!> Get the compute-node element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE FUNCTION GetGlobalElemID_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetGlobalElemID_iElem
!===================================================================================================================================
GetGlobalElemID_iElem = iElem
END FUNCTION GetGlobalElemID_iElem


#if USE_MPI
!==================================================================================================================================!
!> Get the global element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE FUNCTION GetGlobalElemID_fromTotalElem(iElem)
! MODULES
USE MOD_Particle_MPI_Shared_Vars, ONLY:CNTotalElem2GlobalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetGlobalElemID_fromTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalElemID_fromTotalElem = CNTotalElem2GlobalElem(iElem)
END FUNCTION GetGlobalElemID_fromTotalElem
#endif /*USE_MPI*/


!==================================================================================================================================!
!> Initialize GetCNElemID function (mapping of global element ID to compute-node element ID)
!==================================================================================================================================!
SUBROUTINE InitGetCNElemID()
! MODULES
#if USE_MPI
USE MOD_Particle_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetCNElemID => GetCNElemID_iElem
ELSE
  GetCNElemID => GetCNElemID_fromTotalElem
END IF
#else
GetCNElemID => GetCNElemID_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetCNElemID_fromTotalElem(1)
#endif
dummy=GetCNElemID_iElem(1)
END SUBROUTINE InitGetCNElemID


!==================================================================================================================================!
!> Get the CN element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE FUNCTION GetCNElemID_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetCNElemID_iElem
!===================================================================================================================================
GetCNElemID_iElem = iElem
END FUNCTION GetCNElemID_iElem


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE FUNCTION GetCNElemID_fromTotalElem(iElem)
! MODULES
USE MOD_Particle_MPI_Shared_Vars, ONLY:GlobalElem2CNTotalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetCNElemID_fromTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetCNElemID_fromTotalElem = GlobalElem2CNTotalElem(iElem)
END FUNCTION GetCNElemID_fromTotalElem
#endif /*USE_MPI*/


END MODULE MOD_Particle_Mesh_Tools
