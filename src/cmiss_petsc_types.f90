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
MODULE CmissPetscTypes
  
  USE KINDS
  
  IMPLICIT NONE
 
  PRIVATE

#include "petsc/finclude/petsc.h"

  !Module parameters
  
  !Module types

  TYPE PetscISType
    IS :: is
  END TYPE PetscISType

  TYPE PetscISLocalToGloabalMappingType
    ISLocalToGlobalMapping :: isLocalToGlobalMapping
  END TYPE PetscISLocalToGloabalMappingType

  TYPE PetscISColoringType
    ISColoring :: isColoring
  END TYPE PetscISColoringType

  TYPE PetscKspType
    KSP :: ksp
  END TYPE PetscKspType

  TYPE PetscMatType
    Mat :: mat
  END TYPE PetscMatType
  
  TYPE PetscMatColoringType
    MatColoring :: matColoring
  END TYPE PetscMatColoringType
  
  TYPE PETScMatFDColoringType
    MatFDColoring :: matFDColoring
  END TYPE PETScMatFDColoringType

  TYPE PetscPCType
    PC :: pc
  END TYPE PetscPCType
  
  TYPE PetscSnesType
    SNES :: snes
  END TYPE PetscSnesType

  TYPE PetscSnesLineSearchType
    SNESLineSearch :: snesLineSearch
  END TYPE PetscSnesLineSearchType
  
  TYPE PetscTSType
    TS :: ts
  END TYPE PetscTSType
  
  TYPE PetscVecType
    Vec :: vec
  END TYPE PetscVecType
  
  !Interfaces
 
  PUBLIC PetscISType,PetscISLocalToGloabalMappingType,PetscISColoringType,PetscKspType,PetscMatType,PetscMatColoringType, &
    & PetscMatFDColoringType,PetscPCType,PetscSnesType,PetscSnesLineSearchType,PetscTSType,PetscVecType

END MODULE CmissPetscTypes
