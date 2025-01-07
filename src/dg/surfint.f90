!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
#include "eos.h"

!==================================================================================================================================
!> Contains the different Surface integral formulations
!> Computes the Surface integral for all faces using U and updates Ut
!> Computes only inner surface integrals!
!> Surface integrals are separated for each direction
!==================================================================================================================================
MODULE MOD_SurfInt
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
#define WITHnVars 1

PUBLIC:: SurfInt
PUBLIC:: DoSurfInt

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfInt

!==================================================================================================================================
!> Contains the surface integral for conservative quantities
!==================================================================================================================================
MODULE MOD_SurfIntCons
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE SurfIntCons
  MODULE PROCEDURE SurfInt
#if ((FV_ENABLED==2) && (PP_NodeType==1))
  MODULE PROCEDURE SurfIntBlend
#endif
END INTERFACE

INTERFACE DoSurfIntCons
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC:: SurfIntCons
PUBLIC:: DoSurfIntCons

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntCons

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_SurfIntLifting
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVarLifting

INTERFACE SurfIntLifting
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfIntLifting
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC:: SurfIntLifting
PUBLIC:: DoSurfIntLifting

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntLifting

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_SurfIntLifting_gen
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = 3

INTERFACE SurfIntLifting_gen
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfIntLifting_gen
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfIntLifting_gen,DoSurfIntLifting_gen

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntLifting_gen

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_SurfIntPrim
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE SurfIntPrim
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfIntPrim
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC:: SurfIntPrim
PUBLIC:: DoSurfIntPrim

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntPrim
