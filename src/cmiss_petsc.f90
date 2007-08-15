!> \file
!> $Id: cmiss_petsc.f90 27 2007-07-24 16:52:51Z cpb $
!> \author Chris Bradley
!> \brief This module is a CMISS buffer module to the PETSc library.
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
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
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

!> This module is a CMISS buffer module to the PETSc library.
MODULE CMISS_PETSC
  
  USE BASE_ROUTINES
  USE KINDS
  USE ISO_VARYING_STRING
  
  IMPLICIT NONE
 
  PRIVATE

#include "include/finclude/petsc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscviewer.h"

  !Module parameters

  !Module types

  TYPE PETSC_VEC_TYPE
    Vec :: VEC
  END TYPE PETSC_VEC_TYPE

  TYPE PETSC_MAT_TYPE
    Mat :: MAT
  END TYPE PETSC_MAT_TYPE

  !Module variables

  LOGICAL, PARAMETER :: PETSC_HANDLE_ERROR=.TRUE.

  !Interfaces

  INTERFACE

    SUBROUTINE ISLocalToGlobalMappingApply(ctx,n,in,out,ierr)
      ISLocalToGlobalMapping ctx
      PetscInt n
      PetscInt in
      PetscInt out
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingApply
    
    SUBROUTINE ISLocalToGlobalMappingApplyIS(ctx,isin,isout,ierr)
      ISLocalToGlobalMapping ctx
      IS isin
      IS isout
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingApplyIS
    
    SUBROUTINE ISLocalToGlobalMappingCreate(comm,N,globalnum,ctx,ierr)
      MPI_Comm comm
      PetscInt N
      PetscInt globalnum(*)
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingCreate
    
    SUBROUTINE ISLocalToGlobalMappingDestroy(ctx,ierr)
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingDestroy
    
    SUBROUTINE MatAssemblyBegin(A,assemblytype,ierr)
      Mat A
      MatAssemblyType assemblytype
      PetscInt ierr
    END SUBROUTINE MatAssemblyBegin
    
    SUBROUTINE MatAssemblyEnd(A,assemblytype,ierr)
      Mat A
      MatAssemblyType assemblytype
      PetscInt ierr
    END SUBROUTINE MatAssemblyEnd
    
    SUBROUTINE MatCreate(comm,A,ierr)
      MPI_Comm comm
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreate

    SUBROUTINE MatCreateMPIAIJ(comm,localm,localn,globalm,globaln,diagnumbernzperrow,diagnumbernzeachrow,offdiagnumbernzperrow, &
      & offdiagnumbernzeachrow,A,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt localn
      PetscInt globalm
      PetscInt globaln
      PetscInt diagnumbernzperrow
      PetscInt diagnumbernzeachrow(*)
      PetscInt offdiagnumbernzperrow
      PetscInt offdiagnumbernzeachrow(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateMPIAIJ
        
    SUBROUTINE MatCreateMPIDense(comm,localm,localn,globalm,globaln,matrixdata,A,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt localn
      PetscInt globalm
      PetscInt globaln
      PetscScalar matrixdata(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateMPIDense
        
    SUBROUTINE MatCreateSeqAIJ(comm,m,n,numbernzperrow,numbernzeachrow,A,ierr)
      MPI_Comm comm
      PetscInt m
      PetscInt n
      PetscInt numbernzperrow
      PetscInt numbernzeachrow(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateSeqAIJ
          
    SUBROUTINE MatCreateSeqDense(comm,m,n,matrixdata,A,ierr)
      MPI_Comm comm
      PetscInt m
      PetscInt n
      PetscScalar matrixdata(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateSeqDense

    SUBROUTINE MatGetOwnershipRange(A,firstrow,lastrow,ierr)
      Mat A
      PetscInt firstrow
      PetscInt lastrow
      PetscInt ierr
    END SUBROUTINE MatGetOwnershipRange
    
    SUBROUTINE MatSetOption(A,option,ierr)
      Mat A
      MatOption option
      PetscInt ierr
    END SUBROUTINE MatSetOption

    SUBROUTINE MatSetSizes(A,localm,localn,globalM,globalN,ierr)
      Mat A
      PetscInt localm
      PetscInt localn
      PetscInt globalM
      PetscInt globalN
      PetscInt ierr
    END SUBROUTINE MatSetSizes
    
    SUBROUTINE MatSetValues(A,m,mindices,n,nindices,values,insertmode,ierr)
      Mat A
      PetscInt m
      PetscInt mindices(*)
      PetscInt n
      PetscInt nindices(*)
      PetscScalar values(*)
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE MatSetValues
    
    SUBROUTINE PetscFinalize(ierr)
      PetscInt ierr
    END SUBROUTINE PetscFinalize
    
    SUBROUTINE PetscInitialize(file,ierr)
      PetscChar(*) file
      PetscInt ierr
    END SUBROUTINE PetscInitialize

    SUBROUTINE PetscLogPrintSummary(comm,file,ierr)
      MPI_Comm comm
      PetscChar(*) file
      PetscInt ierr
    END SUBROUTINE PetscLogPrintSummary

    SUBROUTINE VecAssemblyBegin(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecAssemblyBegin

    SUBROUTINE VecAssemblyEnd(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecAssemblyEnd

    SUBROUTINE VecCreate(comm,x,ierr)
      MPI_Comm comm
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreate

    SUBROUTINE VecCreateGhost(comm,localm,globalm,nghost,ghosts,x,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt globalm
      PetscInt nghost
      PetscInt ghosts(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateGhost

    SUBROUTINE VecCreateGhostWithArray(comm,localm,globalm,nghost,ghosts,array,x,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt globalm
      PetscInt nghost
      PetscInt ghosts(*)
      PetscScalar array(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateGhostWithArray

    SUBROUTINE VecCreateMPI(comm,localm,globalm,x,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt globalm
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateMPI

    SUBROUTINE VecCreateMPIWithArray(comm,localn,globaln,array,x,ierr)
      MPI_Comm comm
      PetscInt localn
      PetscInt globaln
      PetscScalar array(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateMPIWithArray

    SUBROUTINE VecCreateSeq(comm,m,x,ierr)
      MPI_Comm comm
      PetscInt m
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateSeq

    SUBROUTINE VecCreateSeqWithArray(comm,n,array,x,ierr)
      MPI_Comm comm
      PetscInt n
      PetscScalar array(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateSeqWithArray

    SUBROUTINE VecDestroy(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecDestroy

    SUBROUTINE VecDuplicate(old,new,ierr)
      Vec old,new
      PetscInt ierr
    END SUBROUTINE VecDuplicate

    SUBROUTINE VecGetLocalSize(x,size,ierr)
      Vec x
      PetscInt size
      PetscInt ierr
    END SUBROUTINE VecGetLocalSize

    SUBROUTINE VecGetOwnershipRange(x,low,high,ierr)
      Vec x
      PetscInt low
      PetscInt high
      PetscInt ierr
    END SUBROUTINE VecGetOwnershipRange

    SUBROUTINE VecGetSize(x,size,ierr)
      Vec x
      PetscInt size
      PetscInt ierr
    END SUBROUTINE VecGetSize

    SUBROUTINE VecGhostGetLocalForm(g,l,ierr)
      Vec g
      Vec l
      PetscInt ierr
    END SUBROUTINE VecGhostGetLocalForm

    SUBROUTINE VecGhostRestoreLocalForm(g,l,ierr)
      Vec g
      Vec l
      PetscInt ierr
    END SUBROUTINE VecGhostRestoreLocalForm

   SUBROUTINE VecGhostUpdateBegin(x,insertmode,scattermode,ierr)
      Vec x
      InsertMode insertmode
      ScatterMode scattermode
      PetscInt ierr
    END SUBROUTINE VecGhostUpdateBegin

    SUBROUTINE VecGhostUpdateEnd(x,insertmode,scattermode,ierr)
      Vec x
      InsertMode insertmode
      ScatterMode scattermode
      PetscInt ierr
    END SUBROUTINE VecGhostUpdateEnd

    SUBROUTINE VecSet(x,value,ierr)
      Vec x
      PetscScalar value
      PetscInt ierr
    END SUBROUTINE VecSet

    SUBROUTINE VecSetFromOptions(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecSetFromOptions

    SUBROUTINE VecSetLocalToGlobalMapping(v,ctx,ierr)
      Vec v
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE VecSetLocalToGlobalMapping

    SUBROUTINE VecSetSizes(x,localm,globalm,ierr)
      Vec x
      PetscInt localm,globalm
      PetscInt ierr
    END SUBROUTINE VecSetSizes

    SUBROUTINE VecSetValues(x,n,indices,values,insertmode,ierr)
      Vec x
      PetscInt n
      PetscInt indices(*)
      PetscScalar values(*)
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE VecSetValues

    SUBROUTINE VecSetValuesLocal(x,n,indices,values,insertmode,ierr)
      Vec x
      PetscInt n
      PetscInt indices(*)
      PetscScalar values(*)
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE VecSetValuesLocal

    SUBROUTINE VecView(x,v,ierr)
      Vec x
      PetscViewer v
      PetscInt ierr
    END SUBROUTINE VecView

  END INTERFACE

  PUBLIC PETSC_VEC_TYPE
  
  PUBLIC ADD_VALUES,INSERT_VALUES,PETSC_COMM_WORLD,PETSC_COMM_SELF,PETSC_DECIDE,PETSC_NULL,PETSC_NULL_CHARACTER, &
    & SCATTER_FORWARD,SCATTER_REVERSE,PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_DRAW_WORLD, &
    & PETSC_VIEWER_DRAW_SELF
  
  PUBLIC PETSC_FINALIZE,PETSC_INITIALIZE,PETSC_LOGPRINTSUMMARY

  PUBLIC PETSC_ISLOCALTOGLOBALMAPPINGAPPLY,PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS,PETSC_ISLOCALTOGLOBALMAPPINGCREATE, &
    & PETSC_ISLOCALTOGLOBALMAPPINGDESTROY

  PUBLIC PETSC_VECASSEMBLYBEGIN,PETSC_VECASSEMBLYEND,PETSC_VECCREATE,PETSC_VECCREATEGHOST,PETSC_VECCREATEGHOSTWITHARRAY, &
    & PETSC_VECCREATEMPI,PETSC_VECCREATEMPIWITHARRAY,PETSC_VECCREATESEQ,PETSC_VECCREATESEQWITHARRAY,PETSC_VECDESTROY, &
    & PETSC_VECDUPLICATE,PETSC_VECGETLOCALSIZE,PETSC_VECGETOWNERSHIPRANGE,PETSC_VECGETSIZE,PETSC_VECGHOSTGETLOCALFORM, &
    & PETSC_VECGHOSTRESTORELOCALFORM,PETSC_VECGHOSTUPDATEBEGIN,PETSC_VECGHOSTUPDATEEND,PETSC_VECSET,PETSC_VECSETFROMOPTIONS, &
    & PETSC_VECSETLOCALTOGLOBALMAPPING,PETSC_VECSETSIZES,PETSC_VECSETVALUES,PETSC_VECSETVALUESLOCAL,PETSC_VECVIEW

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_FINALIZE(ERR,ERROR,*)

    !#### Subroutine: PETSC_FINALIZE
    !###  Description:
    !###    Buffer routine to the PETSc PetscFinalize routine.

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_FINALIZE",ERR,ERROR,*999)

    CALL PetscFinalize(ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscFinalize",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_FINALIZE")
    RETURN
999 CALL ERRORS("PETSC_FINALIZE",ERR,ERROR)
    CALL EXITS("PETSC_FINALIZE")
    RETURN 1
  END SUBROUTINE PETSC_FINALIZE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_INITIALIZE(FILE,ERR,ERROR,*)

    !#### Subroutine: PETSC_INITIALIZE
    !###  Description:
    !###    Buffer routine to the PETSc PetscInitialize routine.

    !Argument Variables
    CHARACTER(LEN=*), INTENT(IN) :: FILE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_INITIALIZE",ERR,ERROR,*999)

    CALL PetscInitialize(FILE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscInitialize",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_INITIALIZE")
    RETURN
999 CALL ERRORS("PETSC_INITIALIZE",ERR,ERROR)
    CALL EXITS("PETSC_INITIALIZE")
    RETURN 1
  END SUBROUTINE PETSC_INITIALIZE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLY(CTX,N,IN,OUT,ERR,ERROR,*)

    !#### Subroutine: PETSC_ISLOCALTOGLOBALMAPPINGAPPLY
    !###  Description:
    !###    Buffer routine to the PETSc ISLocalToGlobalMappingApply routine.

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: CTX
    INTEGER(INTG), INTENT(IN) :: N
    INTEGER(INTG), INTENT(IN) :: IN
    INTEGER(INTG), INTENT(OUT) :: OUT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingApply(CTX,N,IN,OUT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingApply",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLY
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS(CTX,ISIN,ISOUT,ERR,ERROR,*)

    !#### Subroutine: PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS
    !###  Description:
    !###    Buffer routine to the PETSc ISLocalToGlobalMappingApplyIS routine.

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: CTX
    INTEGER(INTG), INTENT(IN) :: ISIN
    INTEGER(INTG), INTENT(OUT) :: ISOUT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingApplyIS(CTX,ISIN,ISOUT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingApplyIS",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGCREATE(COMMUNICATOR,N,GLOBALNUM,CTX,ERR,ERROR,*)

    !#### Subroutine: PETSC_ISLOCALTOGLOBALMAPPINGCREATE
    !###  Description:
    !###    Buffer routine to the PETSc ISLocalToGlobalMappingCreate routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: N
    INTEGER(INTG), INTENT(IN) :: GLOBALNUM(*)
    INTEGER(INTG), INTENT(INOUT) :: CTX
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingCreate(COMMUNICATOR,N,GLOBALNUM,CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGCREATE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGDESTROY(CTX,ERR,ERROR,*)

    !#### Subroutine: PETSC_ISLOCALTOGLOBALMAPPINGDESTROY
    !###  Description:
    !###    Buffer routine to the PETSc ISLocalToGlobalMappingDestroy routine.

    !Argument Variables
    INTEGER(INTG), INTENT(INOUT) :: CTX
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingDestroy(CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingDestroy",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGDESTROY
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_LOGPRINTSUMMARY(COMMUNICATOR,FILE,ERR,ERROR,*)

    !#### Subroutine: PETSC_LOGPRINTSUMMARY
    !###  Description:
    !###    Buffer routine to the PETSc PetscLogPrintSummary routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    CHARACTER(LEN=*), INTENT(IN) :: FILE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_LOGPRINTSUMMARY",ERR,ERROR,*999)

    CALL PetscLogPrintSummary(COMMUNICATOR,FILE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscLogPrintSummary",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_LOGPRINTSUMMARY")
    RETURN
999 CALL ERRORS("PETSC_LOGPRINTSUMMARY",ERR,ERROR)
    CALL EXITS("PETSC_LOGPRINTSUMMARY")
    RETURN 1
  END SUBROUTINE PETSC_LOGPRINTSUMMARY
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATASSEMBLYBEGIN(A,ASSEMBLY_TYPE,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATASSEMBLYBEGIN
    !###  Description:
    !###    Buffer routine to the PETSc MatAssemblyBegin routine.

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(IN) :: ASSEMBLY_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATASSEMBLYBEGIN",ERR,ERROR,*999)

    CALL MatAssemblyBegin(A%MAT,ASSEMBLY_TYPE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatAssemblyBegin",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATASSEMBLYBEGIN")
    RETURN
999 CALL ERRORS("PETSC_MATASSEMBLYBEGIN",ERR,ERROR)
    CALL EXITS("PETSC_MATASSEMBLYBEGIN")
    RETURN 1
  END SUBROUTINE PETSC_MATASSEMBLYBEGIN
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATASSEMBLYEND(A,ASSEMBLY_TYPE,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATASSEMBLYEND
    !###  Description:
    !###    Buffer routine to the PETSc MatAssemblyEnd routine.

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(IN) :: ASSEMBLY_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATASSEMBLYEND",ERR,ERROR,*999)

    CALL MatAssemblyEnd(A%MAT,ASSEMBLY_TYPE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatAssemblyEnd",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATASSEMBLYEND")
    RETURN
999 CALL ERRORS("PETSC_MATASSEMBLYEND",ERR,ERROR)
    CALL EXITS("PETSC_MATASSEMBLYEND")
    RETURN 1
  END SUBROUTINE PETSC_MATASSEMBLYEND
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATCREATE(COMMUNICATOR,A,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATCREATE
    !###  Description:
    !###    Buffer routine to the PETSc MatCreate routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATCREATE",ERR,ERROR,*999)

    CALL MatCreate(COMMUNICATOR,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATE")
    RETURN
999 CALL ERRORS("PETSC_MATCREATE",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATE")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATE    
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATCREATEMPIAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
    & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATCREATEMPIAIJ
    !###  Description:
    !###    Buffer routine to the PETSc MatCreateMPIAIJ routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: LOCAL_M
    INTEGER(INTG), INTENT(IN) :: LOCAL_N
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N
    INTEGER(INTG), INTENT(IN) :: DIAG_NUMBER_NZ_PERROW
    INTEGER(INTG), INTENT(IN) :: DIAG_NUMBER_NZ_EACHROW(*)
    INTEGER(INTG), INTENT(IN) :: OFFDIAG_NUMBER_NZ_PERROW
    INTEGER(INTG), INTENT(IN) :: OFFDIAG_NUMBER_NZ_EACHROW(*)
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATCREATEMPIAIJ",ERR,ERROR,*999)

    CALL MatCreateMPIAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
      & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateMPIAIJ",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATEMPIAIJ")
    RETURN
999 CALL ERRORS("PETSC_MATCREATEMPIAIJ",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATEMPIAIJ")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEMPIAIJ
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATCREATEMPIDENSE(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATCREATEMPIDENSE
    !###  Description:
    !###    Buffer routine to the PETSc MatCreateMPIDense routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: LOCAL_M
    INTEGER(INTG), INTENT(IN) :: LOCAL_N
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*)
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATCREATEMPIDENSE",ERR,ERROR,*999)

    CALL MatCreateMPIDense(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateMPIDense",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATEMPIDENSE")
    RETURN
999 CALL ERRORS("PETSC_MATCREATEMPIDENSE",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATEMPIDENSE")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEMPIDENSE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATCREATESEQAIJ(COMMUNICATOR,M,N,NUMBER_NZ_PERROW,NUMBER_NZ_EACHROW,A,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATCREATESEQAIJ
    !###  Description:
    !###    Buffer routine to the PETSc MatCreateSeqAIJ routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: M
    INTEGER(INTG), INTENT(IN) :: N
    INTEGER(INTG), INTENT(IN) :: NUMBER_NZ_PERROW
    INTEGER(INTG), INTENT(IN) :: NUMBER_NZ_EACHROW(*)
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATCREATESEQAIJ",ERR,ERROR,*999)

    CALL MatCreateSeqAIJ(COMMUNICATOR,M,N,NUMBER_NZ_PERROW,NUMBER_NZ_EACHROW,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateSeqAIJ",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATESEQAIJ")
    RETURN
999 CALL ERRORS("PETSC_MATCREATESEQAIJ",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATESEQAIJ")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATESEQAIJ
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATCREATESEQDENSE(COMMUNICATOR,M,N,MATRIX_DATA,A,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATCREATESEQDENSE
    !###  Description:
    !###    Buffer routine to the PETSc MatCreateSeqDense routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: M
    INTEGER(INTG), INTENT(IN) :: N
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*)
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATCREATESEQDENSE",ERR,ERROR,*999)

    CALL MatCreateSeqDense(COMMUNICATOR,M,N,MATRIX_DATA,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateSeqDense",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATESEQDENSE")
    RETURN
999 CALL ERRORS("PETSC_MATCREATESEQDENSE",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATESEQDENSE")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATESEQDENSE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATGETOWNERSHIPRANGE(A,FIRST_ROW,LAST_ROW,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATGETOWNERSHIPRANGE
    !###  Description:
    !###    Buffer routine to the PETSc MatGetOwnershipRange routine.

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(OUT) :: FIRST_ROW
    INTEGER(INTG), INTENT(OUT) :: LAST_ROW
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATGETOWNERSHIPRANGE",ERR,ERROR,*999)

    CALL MatGetOwnershipRange(A%MAT,FIRST_ROW,LAST_ROW,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatGetOwnershipRange",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATGETOWNERSHIPRANGE")
    RETURN
999 CALL ERRORS("PETSC_MATGETOWNERSHIPRANGE",ERR,ERROR)
    CALL EXITS("PETSC_MATGETOWNERSHIPRANGE")
    RETURN 1
  END SUBROUTINE PETSC_MATGETOWNERSHIPRANGE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATSETOPTION(A,OPTION,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATSETOPTION
    !###  Description:
    !###    Buffer routine to the PETSc MatSetOption routine.

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    MatOption, INTENT(IN) :: OPTION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATSETOPTION",ERR,ERROR,*999)

    CALL MatSetOption(A%MAT,OPTION,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetOption",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETOPTION")
    RETURN
999 CALL ERRORS("PETSC_MATSETOPTION",ERR,ERROR)
    CALL EXITS("PETSC_MATSETOPTION")
    RETURN 1
  END SUBROUTINE PETSC_MATSETOPTION
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATSETSIZES(A,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATSETSIZES
    !###  Description:
    !###    Buffer routine to the PETSc MatSetSizes routine.

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(IN) :: LOCAL_M
    INTEGER(INTG), INTENT(IN) :: LOCAL_N
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATSETSIZES",ERR,ERROR,*999)

    CALL MatSetSizes(A%MAT,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetSizes",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETSIZES")
    RETURN
999 CALL ERRORS("PETSC_MATSETSIZES",ERR,ERROR)
    CALL EXITS("PETSC_MATSETSIZES")
    RETURN 1
  END SUBROUTINE PETSC_MATSETSIZES
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_MATSETVALUES(A,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,ERR,ERROR,*)

    !#### Subroutine: PETSC_MATSETVALUES
    !###  Description:
    !###    Buffer routine to the PETSc MatSetValues routine.

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A
    INTEGER(INTG), INTENT(IN) :: M
    INTEGER(INTG), INTENT(IN) :: M_INDICES(*)
    INTEGER(INTG), INTENT(IN) :: N
    INTEGER(INTG), INTENT(IN) :: N_INDICES(*)
    REAL(DP), INTENT(IN) :: VALUES(*)
    InsertMode, INTENT(IN) :: INSERT_MODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_MATSETVALUES",ERR,ERROR,*999)

    CALL MatSetValues(A%MAT,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetValues",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETVALUES")
    RETURN
999 CALL ERRORS("PETSC_MATSETVALUES",ERR,ERROR)
    CALL EXITS("PETSC_MATSETVALUES")
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUES
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECASSEMBLYBEGIN(X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECASSEMBLYBEGIN
    !###  Description:
    !###    Buffer routine to the PETSc VecAssemblyBegin routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECASSEMBLYBEGIN",ERR,ERROR,*999)

    CALL VecAssemblyBegin(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecAssemblyBegin",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECASSEMBLYBEGIN")
    RETURN
999 CALL ERRORS("PETSC_VECASSEMBLYBEGIN",ERR,ERROR)
    CALL EXITS("PETSC_VECASSEMBLYBEGIN")
    RETURN 1
  END SUBROUTINE PETSC_VECASSEMBLYBEGIN
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECASSEMBLYEND(X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECASSEMBLYEND
    !###  Description:
    !###    Buffer routine to the PETSc VecAssemblyEnd routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECASSEMBLYEND",ERR,ERROR,*999)

    CALL VecAssemblyEnd(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecAssemblyEnd",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECASSEMBLYEND")
    RETURN
999 CALL ERRORS("PETSC_VECASSEMBLYEND",ERR,ERROR)
    CALL EXITS("PETSC_VECASSEMBLYEND")
    RETURN 1
  END SUBROUTINE PETSC_VECASSEMBLYEND
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECCREATE(COMMUNICATOR,X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECCREATE
    !###  Description:
    !###    Buffer routine to the PETSc VecCreate routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECCREATE",ERR,ERROR,*999)

    CALL VecCreate(COMMUNICATOR,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATE")
    RETURN
999 CALL ERRORS("PETSC_VECCREATE",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATE")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECCREATEGHOST(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECCREATEGHOST
    !###  Description:
    !###    Buffer routine to the PETSc VecCreateGhost routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE
    INTEGER(INTG), INTENT(IN) :: NUMBER_GHOST
    INTEGER(INTG), INTENT(IN) :: GHOSTS(*)
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEGHOST",ERR,ERROR,*999)

    CALL VecCreateGhost(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateGhost",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEGHOST")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEGHOST",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEGHOST")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEGHOST
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECCREATEGHOSTWITHARRAY(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,ARRAY,X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECCREATEGHOSTWITHARRAY
    !###  Description:
    !###    Buffer routine to the PETSc VecCreateGhostWithArray routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE
    INTEGER(INTG), INTENT(IN) :: NUMBER_GHOST
    INTEGER(INTG), INTENT(IN) :: GHOSTS(*)
    REAL(DP), INTENT(OUT) :: ARRAY(*)
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEGHOSTWITHARRAY",ERR,ERROR,*999)

    CALL VecCreateGhostWithArray(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,ARRAY,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateGhostWithArray",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEGHOSTWITHARRAY")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEGHOSTWITHARRAY",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEGHOSTWITHARRAY")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEGHOSTWITHARRAY
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECCREATEMPI(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECCREATEMPI
    !###  Description:
    !###    Buffer routine to the PETSc VecCreateMPI routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEMPI",ERR,ERROR,*999)

    CALL VecCreateMPI(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateMPI",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEMPI")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEMPI",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEMPI")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEMPI
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECCREATEMPIWITHARRAY(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,ARRAY,X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECCREATEMPIWITHARRAY
    !###  Description:
    !###    Buffer routine to the PETSc VecCreateMPIWithArray routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE
    REAL(DP), INTENT(OUT) :: ARRAY(*)
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEMPIWITHARRAY",ERR,ERROR,*999)

    CALL VecCreateMPIWithArray(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,ARRAY,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateMPIWithArray",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEMPIWITHARRAY")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEMPIWITHARRAY",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEMPIWITHARRAY")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEMPIWITHARRAY
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECCREATESEQ(COMMUNICATOR,SIZE,X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECCREATESEQ
    !###  Description:
    !###    Buffer routine to the PETSc VecCreateSeq routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: SIZE
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECCREATESEQ",ERR,ERROR,*999)

    CALL VecCreateSeq(COMMUNICATOR,SIZE,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateSeq",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATESEQ")
    RETURN
999 CALL ERRORS("PETSC_VECCREATESEQ",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATESEQ")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATESEQ
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECCREATESEQWITHARRAY(COMMUNICATOR,SIZE,ARRAY,X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECCREATESEQWITHARRAY
    !###  Description:
    !###    Buffer routine to the PETSc VecCreateSeqWithArray routine.

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(IN) :: SIZE
    REAL(DP), INTENT(OUT) :: ARRAY(*)
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECCREATESEQWITHARRAY",ERR,ERROR,*999)

    CALL VecCreateSeqWithArray(COMMUNICATOR,SIZE,ARRAY,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateSeqWithArray",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATESEQWITHARRAY")
    RETURN
999 CALL ERRORS("PETSC_VECCREATESEQWITHARRAY",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATESEQWITHARRAY")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATESEQWITHARRAY
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECDESTROY(X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECDESTROY
    !###  Description:
    !###    Buffer routine to the PETSc VecDestroy routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECDESTROY",ERR,ERROR,*999)

    CALL VecDestroy(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecDestroy",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECDESTROY")
    RETURN
999 CALL ERRORS("PETSC_VECDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_VECDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_VECDESTROY
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECDUPLICATE(OLD,NEW,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECDUPLICATE
    !###  Description:
    !###    Buffer routine to the PETSc VecDuplicate routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: OLD,NEW
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECDUPLICATE",ERR,ERROR,*999)

    CALL VecDuplicate(OLD%VEC,NEW%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecDuplicate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECDUPLICATE")
    RETURN
999 CALL ERRORS("PETSC_VECDUPLICATE",ERR,ERROR)
    CALL EXITS("PETSC_VECDUPLICATE")
    RETURN 1
  END SUBROUTINE PETSC_VECDUPLICATE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECGETLOCALSIZE(X,SIZE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECGETLOCALSIZE
    !###  Description:
    !###    Buffer routine to the PETSc VecGetLocalSize routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: SIZE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECGETLOCALSIZE",ERR,ERROR,*999)

    CALL VecGetLocalSize(X%VEC,SIZE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGetLocalSize",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGETLOCALSIZE")
    RETURN
999 CALL ERRORS("PETSC_VECGETLOCALSIZE",ERR,ERROR)
    CALL EXITS("PETSC_VECGETLOCALSIZE")
    RETURN 1
  END SUBROUTINE PETSC_VECGETLOCALSIZE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECGETOWNERSHIPRANGE(X,LOW,HIGH,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECGETOWNERSHIPRANGE
    !###  Description:
    !###    Buffer routine to the PETSc VecGetOwnershipRange routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: LOW
    INTEGER(INTG), INTENT(OUT) :: HIGH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECGETOWNERSHIPRANGE",ERR,ERROR,*999)

    CALL VecGetOwnershipRange(X%VEC,LOW,HIGH,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
        ENDIF
      CALL FLAG_ERROR("PETSc error in VecGetOwnershipRange",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGETOWNERSHIPRANGE")
    RETURN
999 CALL ERRORS("PETSC_VECGETOWNERSHIPRANGE",ERR,ERROR)
    CALL EXITS("PETSC_VECGETOWNERSHIPRANGE")
    RETURN 1
  END SUBROUTINE PETSC_VECGETOWNERSHIPRANGE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECGETSIZE(X,SIZE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECGETSIZE
    !###  Description:
    !###    Buffer routine to the PETSc VecGetSize routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: SIZE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECGETSIZE",ERR,ERROR,*999)

    CALL VecGetSize(X%VEC,SIZE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGetSize",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGETSIZE")
    RETURN
999 CALL ERRORS("PETSC_VECGETSIZE",ERR,ERROR)
    CALL EXITS("PETSC_VECGETSIZE")
    RETURN 1
  END SUBROUTINE PETSC_VECGETSIZE
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECGHOSTGETLOCALFORM(G,L,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECGHOSTGETLOCALFORM
    !###  Description:
    !###    Buffer routine to the PETSc VecGhostGetLocalForm routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: G
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: L
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTGETLOCALFORM",ERR,ERROR,*999)

    CALL VecGhostGetLocalForm(G%VEC,L%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostGetLocalForm",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTGETLOCALFORM")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTGETLOCALFORM",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTGETLOCALFORM")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTGETLOCALFORM
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECGHOSTRESTORELOCALFORM(G,L,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECGHOSTRESTORELOCALFORM
    !###  Description:
    !###    Buffer routine to the PETSc VecGhostRestoreLocalForm routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: G
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: L
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTRESTORELOCALFORM",ERR,ERROR,*999)

    CALL VecGhostRestoreLocalForm(G%VEC,L%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostRestoreLocalForm",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTRESTORELOCALFORM")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTRESTORELOCALFORM",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTRESTORELOCALFORM")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTRESTORELOCALFORM
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECGHOSTUPDATEBEGIN(X,INSERT_MODE,SCATTER_MODE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECGHOSTUPDATEBEGIN
    !###  Description:
    !###    Buffer routine to the PETSc VecGhostUpdateBegin routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    InsertMode, INTENT(IN) :: INSERT_MODE
    ScatterMode, INTENT(IN) :: SCATTER_MODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTUPDATEBEGIN",ERR,ERROR,*999)

    CALL VecGhostUpdateBegin(X%VEC,INSERT_MODE,SCATTER_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostUpdateBegin",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTUPDATEBEGIN")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTUPDATEBEGIN",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTUPDATEBEGIN")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTUPDATEBEGIN
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECGHOSTUPDATEEND(X,INSERT_MODE,SCATTER_MODE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECGHOSTUPDATEEND
    !###  Description:
    !###    Buffer routine to the PETSc VecGhostUpdateEnd routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    InsertMode, INTENT(IN) :: INSERT_MODE
    ScatterMode, INTENT(IN) :: SCATTER_MODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTUPDATEEND",ERR,ERROR,*999)

    CALL VecGhostUpdateEnd(X%VEC,INSERT_MODE,SCATTER_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostUpdateEnd",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTUPDATEEND")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTUPDATEEND",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTUPDATEEND")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTUPDATEEND
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECSET(X,VALUE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECSET
    !###  Description:
    !###    Buffer routine to the PETSc VecSet routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    REAL(DP), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECSET",ERR,ERROR,*999)

    CALL VecSet(X%VEC,VALUE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSet",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSET")
    RETURN
999 CALL ERRORS("PETSC_VECSET",ERR,ERROR)
    CALL EXITS("PETSC_VECSET")
    RETURN 1
  END SUBROUTINE PETSC_VECSET
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECSETFROMOPTIONS(X,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECSETFROMOPTIONS
    !###  Description:
    !###    Buffer routine to the PETSc VecSetFromOptions routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECSETFROMOPTIONS",ERR,ERROR,*999)

    CALL VecSetFromOptions(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetFromOptions",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETFROMOPTIONS")
    RETURN
999 CALL ERRORS("PETSC_VECSETFROMOPTIONS",ERR,ERROR)
    CALL EXITS("PETSC_VECSETFROMOPTIONS")
    RETURN 1
  END SUBROUTINE PETSC_VECSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECSETLOCALTOGLOBALMAPPING(X,CTX,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECSETLOCALTOGLOBALMAPPING
    !###  Description:
    !###    Buffer routine to the PETSc VecSetLocalToGlobalMapping routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(IN) :: CTX
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECSETLOCALTOGLOBALMAPPING",ERR,ERROR,*999)

    CALL VecSetLocalToGlobalMapping(X%VEC,CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetLocalToGlobalMapping",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETLOCALTOGLOBALMAPPING")
    RETURN
999 CALL ERRORS("PETSC_VECSETLOCALTOGLOBALMAPPING",ERR,ERROR)
    CALL EXITS("PETSC_VECSETLOCALTOGLOBALMAPPING")
    RETURN 1
  END SUBROUTINE PETSC_VECSETLOCALTOGLOBALMAPPING
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECSETVALUES(X,N,INDICES,VALUES,INSERT_MODE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECSETVALUES
    !###  Description:
    !###    Buffer routine to the PETSc VecSetValues routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(IN) :: N
    INTEGER(INTG), INTENT(IN) :: INDICES(*)
    REAL(DP), INTENT(IN) :: VALUES(*)
    InsertMode, INTENT(IN) :: INSERT_MODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECSETVALUES",ERR,ERROR,*999)

    CALL VecSetValues(X%VEC,N,INDICES,VALUES,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetValues",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETVALUES")
    RETURN
999 CALL ERRORS("PETSC_VECSETVALUES",ERR,ERROR)
    CALL EXITS("PETSC_VECSETVALUES")
    RETURN 1
  END SUBROUTINE PETSC_VECSETVALUES
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECSETVALUESLOCAL(X,N,INDICES,VALUES,INSERT_MODE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECSETVALUESLOCAL
    !###  Description:
    !###    Buffer routine to the PETSc VecSetValuesLocal routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(IN) :: N
    INTEGER(INTG), INTENT(IN) :: INDICES(*)
    REAL(DP), INTENT(IN) :: VALUES(*)
    InsertMode :: INSERT_MODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECSETVALUESLOCAL",ERR,ERROR,*999)

    CALL VecSetValuesLocal(X%VEC,N,INDICES,VALUES,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetValuesLocal",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETVALUESLOCAL")
    RETURN
999 CALL ERRORS("PETSC_VECSETVALUESLOCAL",ERR,ERROR)
    CALL EXITS("PETSC_VECSETVALUESLOCAL")
    RETURN 1
  END SUBROUTINE PETSC_VECSETVALUESLOCAL
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECSETSIZES(X,LOCAL_SIZE,GLOBAL_SIZE,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECSETSIZES
    !###  Description:
    !###    Buffer routine to the PETSc VecSetSizes routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE,GLOBAL_SIZE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECSETSIZES",ERR,ERROR,*999)

    CALL VecSetSizes(X%VEC,LOCAL_SIZE,GLOBAL_SIZE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetSizes",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETSIZES")
    RETURN
999 CALL ERRORS("PETSC_VECSETSIZES",ERR,ERROR)
    CALL EXITS("PETSC_VECSETSIZES")
    RETURN 1
  END SUBROUTINE PETSC_VECSETSIZES
    
  !
  !================================================================================================================================
  !

  SUBROUTINE PETSC_VECVIEW(X,V,ERR,ERROR,*)

    !#### Subroutine: PETSC_VECVIEW
    !###  Description:
    !###    Buffer routine to the PETSc VecView routine.

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X
    INTEGER(PTR), INTENT(IN) :: V
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PETSC_VECVIEW",ERR,ERROR,*999)

    CALL VecView(X%VEC,V,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecView",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECVIEW")
    RETURN
999 CALL ERRORS("PETSC_VECVIEW",ERR,ERROR)
    CALL EXITS("PETSC_VECVIEW")
    RETURN 1
  END SUBROUTINE PETSC_VECVIEW
    
  !
  !================================================================================================================================
  !

END MODULE CMISS_PETSC
    
