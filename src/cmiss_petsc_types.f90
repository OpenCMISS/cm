!> \file
!> \author Chris Bradley
!> \brief This module contains type definitions related to the PETSc library.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module contains types related to the PETSc library.
MODULE CMISS_PETSC_TYPES
  
  USE KINDS
  
  IMPLICIT NONE
 
  PRIVATE

#include "finclude/petsc.h"
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 1 )
#include "finclude/petscis.h"
#include "finclude/petscksp.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscsnes.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#endif
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
#include "finclude/petscts.h"
#endif
  
  !Module parameters
  
  !Module types

  TYPE PETSC_IS_TYPE
    IS :: IS_
  END TYPE PETSC_IS_TYPE

  TYPE PETSC_ISLOCALTOGLOBALMAPPING_TYPE
    ISLocalToGlobalMapping :: ISLOCALTOGLOBALMAPPING
  END TYPE PETSC_ISLOCALTOGLOBALMAPPING_TYPE

  TYPE PETSC_ISCOLORING_TYPE
    ISColoring :: ISCOLORING
  END TYPE PETSC_ISCOLORING_TYPE

  TYPE PETSC_KSP_TYPE
    KSP :: KSP_
  END TYPE PETSC_KSP_TYPE

  TYPE PETSC_MAT_TYPE
    !PetscScalar :: MAT_DATA(1)
    !PetscOffset :: MAT_OFFSET
    Mat :: MAT
  END TYPE PETSC_MAT_TYPE
  
  !TYPE PETSC_MATCOLORING_TYPE
  !  MatColoring :: MATCOLORING
  !END TYPE PETSC_MATCOLORING_TYPE

  TYPE PETSC_MATFDCOLORING_TYPE
    MatFDColoring :: MATFDCOLORING
  END TYPE PETSC_MATFDCOLORING_TYPE

  TYPE PETSC_PC_TYPE
    PC :: PC_
  END TYPE PETSC_PC_TYPE
  
  TYPE PETSC_SNES_TYPE
    SNES :: SNES_
  END TYPE PETSC_SNES_TYPE

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  TYPE PetscSnesLineSearchType
    SNESLineSearch :: snesLineSearch
  END TYPE PetscSnesLineSearchType
#else
  TYPE PetscSnesLineSearchType
    !Empty type here to avoid adding mess in types.f90
  END TYPE PetscSnesLineSearchType
#endif
  
  TYPE PETSC_TS_TYPE
    TS :: TS_
  END TYPE PETSC_TS_TYPE
  
  TYPE PETSC_VEC_TYPE
    !PetscScalar :: VEC_DATA(1)
    !PetscOffset :: VEC_OFFSET
    Vec :: VEC
  END TYPE PETSC_VEC_TYPE
  
  !Module variables

  !Interfaces
 
  PUBLIC PETSC_IS_TYPE,PETSC_ISLOCALTOGLOBALMAPPING_TYPE,PETSC_ISCOLORING_TYPE,PETSC_KSP_TYPE,PETSC_MAT_TYPE, &
    & PETSC_MATFDCOLORING_TYPE,PETSC_PC_TYPE,PETSC_SNES_TYPE,PetscSnesLineSearchType,PETSC_TS_TYPE,PETSC_VEC_TYPE

END MODULE CMISS_PETSC_TYPES
