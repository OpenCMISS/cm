!> \file
!> $Id$
!> \author Caton Little
!> \brief This module handles writing out FieldML files.
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

!> Temporary IO routines for fluid mechanics

MODULE FIELDML_OUTPUT_ROUTINES

  USE BASE_ROUTINES
  USE INPUT_OUTPUT
  USE KINDS
  USE FIELDML_API
  USE BASIS_ROUTINES
  USE UTIL_ARRAY
  USE CMISS

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  INTEGER(INTG), PARAMETER :: BUFFER_SIZE = 1024

  !INTEGER(INTG), PARAMETER ::
  !INTEGER(INTG), PARAMETER ::
  !INTEGER(INTG), PARAMETER ::

  TYPE(VARYING_STRING) :: errorString

  !Interfaces
  INTERFACE FieldmlOutput_WriteRawData
    MODULE PROCEDURE FieldmlOutput_WriteRawData_Int2
    MODULE PROCEDURE FieldmlOutput_WriteRawData_Real1
    MODULE PROCEDURE FieldmlOutput_WriteRawData_Real2
  END INTERFACE FieldmlOutput_WriteRawData

  INTERFACE

  END INTERFACE

  PUBLIC :: FieldmlOutput_WriteRawData

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_WriteRawData_Int2( parseHandle, parametersHandle, array, append, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: parseHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    INTEGER(C_INT), INTENT(INOUT) :: array(:,:)
    INTEGER(C_INT), INTENT(IN) :: append
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: i, indexCount, count1, count2, handle1, handle2
    INTEGER(C_INT), TARGET :: dummy(0)
    INTEGER(C_INT), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: writer
    
    indexCount = Fieldml_GetSemidenseIndexCount( parseHandle, parametersHandle, 1 )
    IF( indexCount /= 0 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    indexCount = Fieldml_GetSemidenseIndexCount( parseHandle, parametersHandle, 0 )
    IF( indexCount /= 2 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    handle1 = Fieldml_GetSemidenseIndex( parseHandle, parametersHandle, 1, 0 )
    handle2 = Fieldml_GetSemidenseIndex( parseHandle, parametersHandle, 2, 0 )
    IF( ( handle1 == FML_INVALID_HANDLE ) .OR. ( handle2 == FML_INVALID_HANDLE ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    count1 = Fieldml_GetEnsembleDomainElementCount( parseHandle, handle1 )
    count2 = Fieldml_GetEnsembleDomainElementCount( parseHandle, handle2 )

    ALLOCATE( buffer( count1 ) )

    writer = Fieldml_OpenWriter( parseHandle, parametersHandle, append )
    IF( .NOT. C_ASSOCIATED( writer ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    DO i = 1, count2
      buffer( 1:count1 ) = array( i, 1:count1 )
      err = Fieldml_WriteIntSlice( parseHandle, writer, C_LOC(dummy), C_LOC(buffer) )
    ENDDO
    
    err = Fieldml_CloseWriter( parseHandle, writer )

    DEALLOCATE( buffer )
    
  END SUBROUTINE FieldmlOutput_WriteRawData_Int2

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_WriteRawData_Real1( parseHandle, parametersHandle, array, append, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: parseHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    REAL(C_DOUBLE), INTENT(IN), ALLOCATABLE, TARGET :: array(:)
    INTEGER(C_INT), INTENT(IN) :: append
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: i, indexCount, count1, handle1
    INTEGER(C_INT), TARGET :: dummy(0)
    TYPE(C_PTR) :: writer
    
    indexCount = Fieldml_GetSemidenseIndexCount( parseHandle, parametersHandle, 1 )
    IF( indexCount /= 0 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    indexCount = Fieldml_GetSemidenseIndexCount( parseHandle, parametersHandle, 0 )
    IF( indexCount /= 1 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    handle1 = Fieldml_GetSemidenseIndex( parseHandle, parametersHandle, 1, 0 )
    IF( handle1 == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    count1 = Fieldml_GetEnsembleDomainElementCount( parseHandle, handle1 )

    writer = Fieldml_OpenWriter( parseHandle, parametersHandle, append )
    IF( .NOT. C_ASSOCIATED( writer ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    err = Fieldml_WriteDoubleSlice( parseHandle, writer, C_LOC(dummy), C_LOC(array) )
    
    err = Fieldml_CloseWriter( parseHandle, writer )
    
  END SUBROUTINE FieldmlOutput_WriteRawData_Real1

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_WriteRawData_Real2( parseHandle, parametersHandle, array, append, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: parseHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    REAL(C_DOUBLE), INTENT(INOUT), TARGET :: array(:,:)
    INTEGER(C_INT), INTENT(IN) :: append
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: i, indexCount, count1, count2, handle1, handle2
    INTEGER(C_INT), TARGET :: dummy(0)
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: writer
    
    indexCount = Fieldml_GetSemidenseIndexCount( parseHandle, parametersHandle, 1 )
    IF( indexCount /= 0 ) THEN
      WRITE(*,'("BAH1")')
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    indexCount = Fieldml_GetSemidenseIndexCount( parseHandle, parametersHandle, 0 )
    IF( indexCount /= 2 ) THEN
      WRITE(*,'("BAH1")')
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    handle1 = Fieldml_GetSemidenseIndex( parseHandle, parametersHandle, 1, 0 )
    handle2 = Fieldml_GetSemidenseIndex( parseHandle, parametersHandle, 2, 0 )
    IF( ( handle1 == FML_INVALID_HANDLE ) .OR. ( handle2 == FML_INVALID_HANDLE ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    count1 = Fieldml_GetEnsembleDomainElementCount( parseHandle, handle1 )
    count2 = Fieldml_GetEnsembleDomainElementCount( parseHandle, handle2 )

    ALLOCATE( buffer( count1 ) )

    writer = Fieldml_OpenWriter( parseHandle, parametersHandle, append )
    IF( .NOT. C_ASSOCIATED( writer ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    DO i = 1, count2
      buffer( 1:count1 ) = array( i, 1:count1 )
      err = Fieldml_WriteDoubleSlice( parseHandle, writer, C_LOC(dummy), C_LOC(buffer) )
    ENDDO
    
    err = Fieldml_CloseWriter( parseHandle, writer )

    DEALLOCATE( buffer )
    
  END SUBROUTINE FieldmlOutput_WriteRawData_Real2

  !
  !================================================================================================================================
  !

END MODULE FIELDML_OUTPUT_ROUTINES
