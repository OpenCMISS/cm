!> \file
!> $Id: interface_matrices_routines.f90 690 2009-09-30 23:27:16Z chrispbradley $
!> \author Chris Bradley
!> \brief This module contains all interface matrices routines.
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

!>This module contains all interface matrices routines.
MODULE INTERFACE_MATRICES_ROUTINES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the interface matrices for the interface equations
  SUBROUTINE INTERFACE_MATRICES_CREATE_FINISH(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<The pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,NUMBER_OF_NON_ZEROS
    INTEGER(INTG), POINTER :: ROW_INDICES(:),COLUMN_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)

    NULLIFY(ROW_DOMAIN_MAP)
    NULLIFY(COLUMN_DOMAIN_MAP)

    CALL ENTERS("INTERFACE_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Interface matrices have already been finished.",ERR,ERROR,*998)
      ELSE
        INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          IF(ASSOCIATED(ROW_DOMAIN_MAP)) THEN
            !Now create the individual interface matrices
            DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
              INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
                IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
                  !Create the distributed equations matrix
                  CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,INTERFACE_MATRICES% &
                    & MATRICES(matrix_idx)%PTR%MATRIX,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(INTERFACE_MATRIX%MATRIX,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(INTERFACE_MATRIX%MATRIX,INTERFACE_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                  !Calculate and set the matrix structure/sparsity pattern
                  IF(INTERFACE_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                    & INTERFACE_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                    CALL INTERFACE_MATRIX_STRUCTURE_CALCULATE(INTERFACE_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                      & ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(INTERFACE_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(INTERFACE_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                      & ERR,ERROR,*999)
                    IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                  ENDIF
                  CALL DISTRIBUTED_MATRIX_CREATE_FINISH(INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="Column domain map for interface matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Interface matrix for matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx                
            !Finish up
            INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Row domain map is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    CALL INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MATRICES_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of the interface matrices and rhs for the interface equations
  SUBROUTINE INTERFACE_MATRICES_CREATE_START(INTERFACE_EQUATIONS,INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<The pointer to the interface equations to create the interface equations matrices for
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<On return, a pointer to the interface matrices being created. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables

    CALL ENTERS("INTERFACE_MATRICES_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN      
      IF(INTERFACE_EQUATIONS%EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
          CALL FLAG_ERROR("Interface matrices is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(INTERFACE_MATRICES)
          !Initialise the interface matrices
          CALL INTERFACE_MATRICES_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
          !Return the pointer
          INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRICES_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_CREATE_START")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MATRICES_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the interface matrices
  SUBROUTINE INTERFACE_MATRICES_DESTROY(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer the interface matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      CALL INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("INTERFACE_MATRICES_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_DESTROY",ERR,ERROR)    
    CALL EXITS("INTERFACE_MATRICES_DESTROY")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MATRICES_DESTROY
  !
  !================================================================================================================================
  !

  !>Finalise a interface matrix and deallocate all memory
  SUBROUTINE INTERFACE_MATRIX_FINALISE(INTERFACE_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX !<A pointer to the interface matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
      IF(ASSOCIATED(EQUATIONS_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
      CALL INTERFACE_MATRICES_ELEMENT_MATRIX_FINALISE(INTERFACE_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_MATRIX)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalise the interface matrices and deallocate all memory.
  SUBROUTINE INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   
    CALL ENTERS("INTERFACE_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      DEALLOCATE(INTERFACE_MATRICES)
    ENDIF
       
    CALL EXITS("INTERFACE_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the interface matrices for the interface equations.
  SUBROUTINE INTERFACE_MATRICES_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*)
    
     !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to initialise the interface matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("INTERFACE_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MATRICES)) THEN
        CALL FLAG_ERROR("Interface matrices is already associated for this interface equations.",ERR,ERROR,*998)
      ELSE
        INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
            ALLOCATE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface equations interface matrices.",ERR,ERROR,*999)
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED=.FALSE.
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_MAPPING=>INTERFACE_MAPPING
            NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MATRICES%SOLVER_MAPPING)
            EQUATIONS%EQUATIONS_MATRICES%NUMBER_OF_ROWS=EQUATIONS_MAPPING%NUMBER_OF_ROWS
            EQUATIONS%EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS=EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
            EQUATIONS%EQUATIONS_MATRICES%NUMBER_OF_GLOBAL_ROWS=EQUATIONS_MAPPING%NUMBER_OF_GLOBAL_ROWS
           ELSE
            CALL FLAG_ERROR("Interface mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface equations interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MATRICES_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS%EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_INITIALISE



END MODULE INTERFACE_MATRICES_ROUTINES
