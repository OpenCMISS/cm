!> \file
!> $Id: equations_mapping_routines.f90 28 2007-07-27 08:35:14Z cpb $
!> \author Chris Bradley
!> \brief This module handles all equations mapping routines.
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

!>This module handles all equations mapping routines.
MODULE EQUATIONS_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES
 
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC EQUATIONS_MAPPING_CREATE_FINISH,EQUATIONS_MAPPING_CREATE_START,EQUATIONS_MAPPING_DESTROY, &
    & EQUATIONS_MAPPING_MATRICES_NUMBER_SET,EQUATIONS_MAPPING_MATRICES_VARIABLE_TYPES_SET, &
    & EQUATIONS_MAPPING_MATRIX_COEFFICIENTS_SET,EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the equations/dofs mapping.
  SUBROUTINE EQUATIONS_MAPPING_CALCULATE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,dof_idx,matrix_idx,NUMBER_OF_GLOBAL_DOFS,NUMBER_OF_LOCAL_DOFS,NUMBER_OF_ROWS,row_idx, &
      & TOTAL_NUMBER_OF_ROWS,variable_idx,variable_type
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,SOURCE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
          IF(ASSOCIATED(EQUATIONS_SET)) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
              IF(ASSOCIATED(DECOMPOSITION)) THEN
                !Allocate and initialise the variable type maps
                ALLOCATE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping variable to equations map.",ERR,ERROR,*999)
                DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                  CALL EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_INITIALISE(EQUATIONS_MAPPING% &
                    & VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type),ERR,ERROR,*999)
                  EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE_INDEX=variable_type
                  EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE_TYPE=variable_type
                  EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE=>DEPENDENT_FIELD% &
                    & VARIABLE_TYPE_MAP(variable_type)%PTR
                ENDDO !variable_type
                IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                  !Calculate the number of variable type maps and initialise
                  DO matrix_idx=1,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                    variable_type=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)
                    EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES= &
                      & EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES+1
                  ENDDO !matrix_idx
                  IF(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE/=0) THEN
                    EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE)% &
                      & NUMBER_OF_EQUATIONS_MATRICES=-1
                    EQUATIONS_MAPPING%RHS_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE)%PTR
                  ENDIF
                ELSE
                  IF(EQUATIONS_MAPPING%RESIDUAL_VARIABLE_TYPE/=0) THEN
                    EQUATIONS_MAPPING%RESIDUAL_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING% &
                      & RESIDUAL_VARIABLE_TYPE)%PTR
                    EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(EQUATIONS_MAPPING%RESIDUAL_VARIABLE_TYPE)% &
                      & NUMBER_OF_EQUATIONS_MATRICES=0
                    IF(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE/=0) THEN                      
                      EQUATIONS_MAPPING%RHS_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE)%PTR
                      EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE)% &
                        & NUMBER_OF_EQUATIONS_MATRICES=0
                    ELSE
                      CALL FLAG_ERROR("Equations mapping RHS variable type is not set.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations mapping residual variable type is not set.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES=0
                !Allocate and initialise the variable to equations matrices maps
                DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                    NUMBER_OF_LOCAL_DOFS=DEPENDENT_VARIABLE%NUMBER_OF_DOFS
                    NUMBER_OF_GLOBAL_DOFS=DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS                  
                    IF(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES==-1) THEN
                      ALLOCATE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP( &
                        & NUMBER_OF_LOCAL_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to equations matrices maps dof to rows map.", &
                        & ERR,ERROR,*999)
                      EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP=0
                      EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES=EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES+1
                    ELSE IF(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                      & NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                      ALLOCATE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%EQUATIONS_MATRIX_NUMBERS( &
                        & EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES), &
                        & STAT=ERR)
                      IF(ERR/=0) &
                        & CALL FLAG_ERROR("Could not allocate variable to equations matrices maps equations matrix numbers.", &
                        & ERR,ERROR,*999)
                      ALLOCATE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_COLUMNS_MAPS( &
                        & EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES), &
                        & STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to equations matrices maps dof to columns map.", &
                        & ERR,ERROR,*999)                
                      EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%EQUATIONS_MATRIX_NUMBERS=0
                      EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES=0
                      DO matrix_idx=1,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                        IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)==variable_type) THEN
                          EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES= &
                            & EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES+1
                          EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%EQUATIONS_MATRIX_NUMBERS( &
                            & EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES) = &
                            & matrix_idx
                          ALLOCATE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_COLUMNS_MAPS( &
                            & EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES)% &
                            & COLUMN_DOF(NUMBER_OF_LOCAL_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to columns map column dof.",ERR,ERROR,*999)
                          DO dof_idx=1,NUMBER_OF_LOCAL_DOFS
                            !1-1 mapping for now
                            EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_COLUMNS_MAPS( &
                              & EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                              & NUMBER_OF_EQUATIONS_MATRICES)%COLUMN_DOF(dof_idx)=dof_idx
                          ENDDO !dof_idx
                        ENDIF
                      ENDDO !matrix_idx
                      ALLOCATE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP( &
                        & NUMBER_OF_LOCAL_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to equations matrices maps dof to rows map.", &
                        & ERR,ERROR,*999)
                      DO dof_idx=1,NUMBER_OF_LOCAL_DOFS
                        !1-1 mappings for now.
                        EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP(dof_idx)=dof_idx
                      ENDDO !dof_idx
                      EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES=EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES+1
                    ENDIF
                  ENDIF
                ENDDO !variable_type
                !Allocate and initialise the variable types            
                ALLOCATE(EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping matrix variable types.",ERR,ERROR,*999)
                EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES=0
                DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                  IF(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                    EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES=EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES+1
                    EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES)=variable_type
                  ENDIF
                ENDDO !variable_type
                !Allocate and initialise the equations matrix to variable maps types
                ALLOCATE(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping equations matrix to variable maps.",ERR,ERROR,*999)
                !Create the individual matrix maps and column maps
                DO matrix_idx=1,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                  variable_type=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  NUMBER_OF_LOCAL_DOFS=DEPENDENT_VARIABLE%NUMBER_OF_DOFS
                  NUMBER_OF_GLOBAL_DOFS=DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS                  
                  CALL EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_INITIALISE(EQUATIONS_MAPPING% &
                    & EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx),ERR,ERROR,*999)
                  EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%MATRIX_NUMBER=matrix_idx
                  EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE_TYPE=variable_type
                  EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE=>DEPENDENT_VARIABLE
                  EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%NUMBER_OF_COLUMNS=NUMBER_OF_GLOBAL_DOFS
                  EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%MATRIX_COEFFICIENT=EQUATIONS_MAPPING% &
                    CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                  ALLOCATE(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%COLUMN_TO_DOF_MAP( &
                    & NUMBER_OF_GLOBAL_DOFS),STAT=ERR)                  
                  DO column_idx=1,NUMBER_OF_GLOBAL_DOFS
                    !1-1 mapping for now
                    EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%COLUMN_TO_DOF_MAP(column_idx)=column_idx
                  ENDDO !column_idx
                  EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING=> &
                    & DEPENDENT_VARIABLE%DOMAIN_MAPPING
                ENDDO !matrix_idx
                IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                  !Check that the number of rows are consistent across the matrices and RHS vector.
                  NUMBER_OF_ROWS=EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(1)%VARIABLE%NUMBER_OF_DOFS
                  TOTAL_NUMBER_OF_ROWS=EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(1)%VARIABLE%TOTAL_NUMBER_OF_DOFS
                  DO matrix_idx=2,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                    IF(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The number of rows in equations matrix 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the number of rows in equations matrix "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" ("// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE% &
                        & NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE%TOTAL_NUMBER_OF_DOFS/= &
                      & TOTAL_NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The total number of rows in equations matrix 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the total number of rows in equations matrix "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" ("// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE% &
                        & TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !matrix_idx
                ELSE
                  IF(EQUATIONS_MAPPING%RESIDUAL_VARIABLE_TYPE/=0) THEN
                    NUMBER_OF_ROWS=EQUATIONS_MAPPING%RESIDUAL_VARIABLE%NUMBER_OF_DOFS
                    TOTAL_NUMBER_OF_ROWS=EQUATIONS_MAPPING%RESIDUAL_VARIABLE%TOTAL_NUMBER_OF_DOFS
                  ELSE
                    CALL FLAG_ERROR("Residual variable type has not been set.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE/=0) THEN
                  IF(EQUATIONS_MAPPING%RHS_VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                    LOCAL_ERROR="Invalid equations set up. The number of rows in the equations matrix 1 ("// &
                      & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                      & ") does not match the number of rows in the RHS vector ("// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%RHS_VARIABLE%NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  IF(EQUATIONS_MAPPING%RHS_VARIABLE%TOTAL_NUMBER_OF_DOFS/=TOTAL_NUMBER_OF_ROWS) THEN
                    LOCAL_ERROR="Invalid equations set up. The total number of rows in the equations matrix 1 ("// &
                      & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                      & ") does not match the total number of rows in the RHS vector ("// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%RHS_VARIABLE%TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDIF
                !Calculate the row maps. For now assume 1-1 mapping.
                ALLOCATE(EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(NUMBER_OF_ROWS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping equations row to variables maps.",ERR,ERROR,*999)
                DO row_idx=1,NUMBER_OF_ROWS
                  ALLOCATE(EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(row_idx)%ROW_TO_DOFS_MAP(EQUATIONS_MAPPING% &
                    & NUMBER_OF_MATRIX_VARIABLES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to variables maps row to dofs map.",ERR,ERROR,*999)
                  DO variable_idx=1,EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES
                    EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(row_idx)%ROW_TO_DOFS_MAP(variable_idx)=row_idx
                  ENDDO !variable_idx
                  IF(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE/=0) EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(row_idx)% &
                    & ROW_TO_RHS_DOF=row_idx
                ENDDO !row_idx
                EQUATIONS_MAPPING%NUMBER_OF_ROWS=NUMBER_OF_ROWS
                EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS=TOTAL_NUMBER_OF_ROWS
                !Calcuate the source row maps if there is a source field
                IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
                  SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
                  IF(ASSOCIATED(SOURCE_FIELD)) THEN
                    CALL EQUATIONS_MAPPING_SOURCE_MAPPINGS_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*999)
!!TODO: source mappings!
                  ELSE
                    CALL FLAG_ERROR("Source field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                  !For now make the rows dofs domain mapping the same as the first matrix variable domain mapping
                  EQUATIONS_MAPPING%ROW_DOFS_MAPPING=>EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(1)%VARIABLE%DOMAIN_MAPPING
                ELSE
                  !For now make the rows dofs domain mapping the same as the residual variable domain mapping
                  EQUATIONS_MAPPING%ROW_DOFS_MAPPING=>EQUATIONS_MAPPING%RESIDUAL_VARIABLE%DOMAIN_MAPPING
                ENDIF
              ELSE
                CALL FLAG_ERROR("Dependent field decomposition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations equations set is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations mapping equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_CALCULATE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_CALCULATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_CALCULATE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equations mapping
  SUBROUTINE EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping has already been finished.",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          CALL EQUATIONS_MAPPING_CALCULATE(EQUATIONS_MAPPING,ERR,ERROR,*999)
          CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
          EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED=.TRUE.            
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_CREATE_FINISH")
    RETURN
999 CALL EQUATIONS_MAPPING_FINALISE(EQUATIONS_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equations mapping for a problem solution
  SUBROUTINE EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equation to create the equations mapping from.
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<On return, a pointer to the equations mapping. This must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    
    CALL ENTERS("EQUATIONS_MAPPING_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Equations has already been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          CALL FLAG_ERROR("Equations mapping is already assocaited",ERR,ERROR,*999)
        ELSE
          NULLIFY(EQUATIONS_MAPPING)
          CALL EQUATIONS_MAPPING_INITIALISE(EQUATIONS,ERR,ERROR,*999)
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises an equations mapping create values cache and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations mapping create values cache 
  SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Equations mapping create values cache is already associated",ERR,ERROR,*998)
      ELSE
        EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
          IF(ASSOCIATED(EQUATIONS_SET)) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              IF(DEPENDENT_FIELD%NUMBER_OF_VARIABLES==1) THEN
                IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                  CALL FLAG_ERROR("Dependent field only has one variable which cannot be mapped to both a global matrix "// &
                    & "and rhs vector",ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Dependent field only has one variable which cannot be ampped to both the residual "// &
                    & "and rhs vector",ERR,ERROR,*999)
                ENDIF           
              ELSE IF(DEPENDENT_FIELD%NUMBER_OF_VARIABLES==2) THEN
                IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                  EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES=1
                  EQUATIONS_MAPPING%RHS_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(2)%PTR%VARIABLE_TYPE
                ELSE
                  EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES=0
                  EQUATIONS_MAPPING%RESIDUAL_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(1)%PTR%VARIABLE_TYPE
                  EQUATIONS_MAPPING%RHS_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(2)%PTR%VARIABLE_TYPE
                ENDIF
              ELSE
                IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                  EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES=DEPENDENT_FIELD%NUMBER_OF_VARIABLES-1
                  EQUATIONS_MAPPING%RHS_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(2)%PTR%VARIABLE_TYPE
                ELSE
                  EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES=DEPENDENT_FIELD%NUMBER_OF_VARIABLES-2
                  EQUATIONS_MAPPING%RESIDUAL_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(1)%PTR%VARIABLE_TYPE
                  EQUATIONS_MAPPING%RHS_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(2)%PTR%VARIABLE_TYPE
                ENDIF
              ENDIF
              ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping create values cache",ERR,ERROR,*999)
              ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING% &
                & NUMBER_OF_EQUATIONS_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping create values cache matrix variable types", &
                & ERR,ERROR,*999)
              ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(EQUATIONS_MAPPING% &
                & NUMBER_OF_EQUATIONS_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping create values cache matrix coefficients", &
                & ERR,ERROR,*999)
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=0
              IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1)=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(1)%PTR% &
                  & VARIABLE_TYPE
                DO matrix_idx=2,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)= &
                    & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR%VARIABLE_TYPE
                ENDDO !matrix_idx
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
              ELSE
                DO matrix_idx=1,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)= &
                    & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR%VARIABLE_TYPE
                ENDDO !matrix_idx
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The equations equations set is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The equations mapping equations is not associated",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Destroy an equations mapping.
  SUBROUTINE EQUATIONS_MAPPING_DESTROY(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer the equations mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      CALL EQUATIONS_MAPPING_FINALISE(EQUATIONS_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("EQUATIONS_MAPPING_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_DESTROY",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise an equations matrix to variable maps and deallocate all memory.
  SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_FINALISE(EQUATIONS_MATRIX_TO_VARIABLE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TO_VARIABLE_MAP_TYPE) :: EQUATIONS_MATRIX_TO_VARIABLE_MAP !<The equations matrix to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_MATRIX_TO_VARIABLE_MAP%COLUMN_TO_DOF_MAP)) &
      & DEALLOCATE(EQUATIONS_MATRIX_TO_VARIABLE_MAP%COLUMN_TO_DOF_MAP)
    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an equations matrix to variable maps.
  SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_INITIALISE(EQUATIONS_MATRIX_TO_VARIABLE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TO_VARIABLE_MAP_TYPE) :: EQUATIONS_MATRIX_TO_VARIABLE_MAP !<The equations matrix to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_MATRIX_TO_VARIABLE_MAP%MATRIX_NUMBER=0
    EQUATIONS_MATRIX_TO_VARIABLE_MAP%VARIABLE_TYPE=0
    NULLIFY(EQUATIONS_MATRIX_TO_VARIABLE_MAP%VARIABLE)
    EQUATIONS_MATRIX_TO_VARIABLE_MAP%NUMBER_OF_COLUMNS=0
    EQUATIONS_MATRIX_TO_VARIABLE_MAP%MATRIX_COEFFICIENT=1.0_DP !Matrices in an equation set are added by default
    NULLIFY(EQUATIONS_MATRIX_TO_VARIABLE_MAP%COLUMN_DOFS_MAPPING)
    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise an equations row to variable map and deallocate all memory
  SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_ROW_TO_VARIABLE_MAP_FINALISE(EQUATIONS_ROW_TO_VARIABLE_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_ROW_TO_VARIABLE_MAP_TYPE) :: EQUATIONS_ROW_TO_VARIABLE_MAP !<The equations row to variable map to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_EQUATIONS_ROW_TO_VARIABLE_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_ROW_TO_VARIABLE_MAP%ROW_TO_DOFS_MAP)) &
      & DEALLOCATE(EQUATIONS_ROW_TO_VARIABLE_MAP%ROW_TO_DOFS_MAP)
    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_ROW_TO_VARIABLE_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_EQUATIONS_ROW_TO_VARIABLE_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_ROW_TO_VARIABLE_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_ROW_TO_VARIABLE_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variables and the equations matrices
  SUBROUTINE EQUATIONS_MAPPING_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,NUMBER_OF_EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set the number of matrices for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_EQUATIONS_MATRICES !<The number of equations matrices for the mapping.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_VARIABLE_TYPES(:)
    REAL(DP), ALLOCATABLE :: OLD_MATRIX_COEFFICIENTS(:)
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
            IF(ASSOCIATED(EQUATIONS_SET)) THEN            
              !Check number of matrices to create is valid
              IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                IF(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE==0) THEN
                  IF(NUMBER_OF_EQUATIONS_MATRICES<1.OR.NUMBER_OF_EQUATIONS_MATRICES>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    LOCAL_ERROR="The requested number of matrices ("// &
                      & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_EQUATIONS_MATRICES,"*",ERR,ERROR))// &
                      & ") is invalid. For problems without a equations set RHS the number must be between >= 1 and <= "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE            
                  IF(NUMBER_OF_EQUATIONS_MATRICES<1.OR.NUMBER_OF_EQUATIONS_MATRICES>=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    LOCAL_ERROR="The requested number of matrices ("// &
                      & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_EQUATIONS_MATRICES,"*",ERR,ERROR))// &
                      & ") is invalid. For problems with a equations set RHS the number must be between >= 1 and < "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                IF(NUMBER_OF_EQUATIONS_MATRICES<0.OR.NUMBER_OF_EQUATIONS_MATRICES>FIELD_NUMBER_OF_VARIABLE_TYPES-2) THEN
                  LOCAL_ERROR="The requested number of matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_EQUATIONS_MATRICES,"*",ERR,ERROR))// &
                    & ") is invalid. For nonlinear problems the number must be between >= 0 and <= "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES-2,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
              !If we need to reallocate and reset all the create_values cache arrays and change the number of matrices
              IF(NUMBER_OF_EQUATIONS_MATRICES/=EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES) THEN
                ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types",ERR,ERROR,*999)
                ALLOCATE(OLD_MATRIX_COEFFICIENTS(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix coefficients",ERR,ERROR,*999)
                OLD_MATRIX_VARIABLE_TYPES=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
                OLD_MATRIX_COEFFICIENTS=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS
                DEALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
                DEALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(NUMBER_OF_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types",ERR,ERROR,*999)
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(NUMBER_OF_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix coefficients",ERR,ERROR,*999)
                IF(NUMBER_OF_EQUATIONS_MATRICES>EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES) THEN
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES)= &
                    & OLD_MATRIX_VARIABLE_TYPES
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES+1: &
                    & NUMBER_OF_EQUATIONS_MATRICES)=OLD_MATRIX_VARIABLE_TYPES(1)
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES)= &
                    & OLD_MATRIX_COEFFICIENTS
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES+1: &
                    & NUMBER_OF_EQUATIONS_MATRICES)=OLD_MATRIX_COEFFICIENTS(1)
                ELSE
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:NUMBER_OF_EQUATIONS_MATRICES)= &
                    & OLD_MATRIX_VARIABLE_TYPES(1:NUMBER_OF_EQUATIONS_MATRICES)
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:NUMBER_OF_EQUATIONS_MATRICES)= &
                    & OLD_MATRIX_COEFFICIENTS(1:NUMBER_OF_EQUATIONS_MATRICES)
                ENDIF
                EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES=NUMBER_OF_EQUATIONS_MATRICES
                IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
                IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations equations set is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations mapping equations is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MAPPING_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)    
    IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)    
    CALL ERRORS("EQUATIONS_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_MATRICES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variable types and the equations matrices
  SUBROUTINE EQUATIONS_MAPPING_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,MATRIX_VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(IN) :: MATRIX_VARIABLE_TYPES(:) !<The matrix variable types to map to each equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MAPPING_MATRICES_VARIABLE_TYPES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(SIZE(MATRIX_VARIABLE_TYPES,1)==EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES) THEN
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  !Check input values
                  DO matrix_idx=1,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                    !Check to see if the variable numbers are defined on the dependent field
                    IF(MATRIX_VARIABLE_TYPES(matrix_idx)>=1.OR. &
                      & MATRIX_VARIABLE_TYPES(matrix_idx)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                      IF(.NOT.ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(MATRIX_VARIABLE_TYPES(matrix_idx))%PTR)) THEN
                        LOCAL_ERROR="The matrix variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(MATRIX_VARIABLE_TYPES(matrix_idx),"*",ERR,ERROR))//" for matrix "// &
                          & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                          & " is not defined on the dependent field."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The matrix variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(MATRIX_VARIABLE_TYPES(matrix_idx),"*",ERR,ERROR))//" for matrix "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                        & " is invalid. The variable types must be >= 1 and <= "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    !Check to see if the matrix variable number is not already being used for the RHS vector
                    IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)== &
                      & EQUATIONS_MAPPING%RHS_VARIABLE_TYPE) THEN
                      LOCAL_ERROR="The matrix variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(MATRIX_VARIABLE_TYPES(matrix_idx),"*",ERR,ERROR))//" for matrix  "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                        & " is the same as the current RHS global variable type."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
                    ENDIF
                  ENDDO !matrix_idx
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=MATRIX_VARIABLE_TYPES
                ELSE
                  CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations mapping equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Invalid size of matrix variable types. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_VARIABLE_TYPES,1),"*",ERR,ERROR))// &
              & ") must match the number of equations matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_MATRICES_VARIABLE_TYPES_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_MATRICES_VARIABLE_TYPES_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_MATRICES_VARIABLE_TYPES_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_MATRICES_VARIABLE_TYPES_SET

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set rhs vector.
  SUBROUTINE EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,RHS_VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: RHS_VARIABLE_TYPE !<The variable type associated with the equations set rhs vector. If the problem does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping have been finished.",ERR,ERROR,*999)
      ELSE
        IF(RHS_VARIABLE_TYPE==0) THEN
          EQUATIONS_MAPPING%RHS_VARIABLE_TYPE=0
        ELSE
          IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  !Check that the given RHS variable number is not already being used for a matrix
                  DO matrix_idx=1,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
                    IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)==RHS_VARIABLE_TYPE) THEN
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " is the same as the variable type for equations matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//"."
                    ENDIF
                  ENDDO !matrix_idx        
                  !Check the RHS variable number is defined on the dependent field
                  IF(RHS_VARIABLE_TYPE>=1.AND.RHS_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(RHS_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%RHS_VARIABLE_TYPE=RHS_VARIABLE_TYPE
                    ELSE
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " is not defined on the dependent field."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The specified RHS variable type of "//TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The number must either be zero or >= 1 and <= "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations mapping equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations mapping create values cache is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Finalise an equations mapping source mappings and deallocate all memory
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE(SOURCE_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOURCE_EQUATIONS_MATRICES_MAP_TYPE), POINTER :: SOURCE_MAPPINGS !<A pointer to the source equations matrices mappings to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOURCE_MAPPINGS)) THEN
      IF(ALLOCATED(SOURCE_MAPPINGS%EQUATIONS_ROW_TO_SOURCE_DOF_MAP)) &
        & DEALLOCATE(SOURCE_MAPPINGS%EQUATIONS_ROW_TO_SOURCE_DOF_MAP)
      IF(ALLOCATED(SOURCE_MAPPINGS%SOURCE_DOF_TO_EQUATIONS_ROW_MAP)) &
        & DEALLOCATE(SOURCE_MAPPINGS%SOURCE_DOF_TO_EQUATIONS_ROW_MAP)
      DEALLOCATE(SOURCE_MAPPINGS)
    ENDIF
    
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an equations mapping source mappings.
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPINGS_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the source mappings for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%SOURCE_MAPPINGS)) THEN
        CALL FLAG_ERROR("Equations mapping source mappings is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%SOURCE_MAPPINGS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mappings source mappings.",ERR,ERROR,*999)        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE(EQUATIONS_MAPPING%SOURCE_MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_INITIALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPINGS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPINGS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise an equations mapping equations matrix map.
  SUBROUTINE EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_COLUMN_MAP_FINALISE(VARIABLE_TO_EQUATIONS_COLUMN_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VARIABLE_TO_EQUATIONS_COLUMN_MAP_TYPE) :: VARIABLE_TO_EQUATIONS_COLUMN_MAP !<The variable dof to equations column map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_COLUMN_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(VARIABLE_TO_EQUATIONS_COLUMN_MAP%COLUMN_DOF)) &
      & DEALLOCATE(VARIABLE_TO_EQUATIONS_COLUMN_MAP%COLUMN_DOF)
    
    CALL EXITS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_COLUMN_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_COLUMN_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_COLUMN_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_COLUMN_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations matrices map and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_FINALISE(VARIABLE_TO_EQUATIONS_MATRICES_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VARIABLE_TO_EQUATIONS_MATRICES_MAP_TYPE) :: VARIABLE_TO_EQUATIONS_MATRICES_MAP !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(VARIABLE_TO_EQUATIONS_MATRICES_MAP%EQUATIONS_MATRIX_NUMBERS)) &
      & DEALLOCATE(VARIABLE_TO_EQUATIONS_MATRICES_MAP%EQUATIONS_MATRIX_NUMBERS)
    IF(ALLOCATED(VARIABLE_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS)) THEN
      DO matrix_idx=1,SIZE(VARIABLE_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS,1)
        CALL EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_COLUMN_MAP_FINALISE(VARIABLE_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS( &
          & matrix_idx),ERR,ERROR,*999)
      ENDDO !matrix_idx
      DEALLOCATE(VARIABLE_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS)
    ENDIF
    IF(ALLOCATED(VARIABLE_TO_EQUATIONS_MATRICES_MAP%DOF_TO_ROWS_MAP)) &
      & DEALLOCATE(VARIABLE_TO_EQUATIONS_MATRICES_MAP%DOF_TO_ROWS_MAP)
    
    CALL EXITS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an equations mapping equations matrix map.
  SUBROUTINE EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_INITIALISE(VARIABLE_TO_EQUATIONS_MATRICES_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VARIABLE_TO_EQUATIONS_MATRICES_MAP_TYPE) :: VARIABLE_TO_EQUATIONS_MATRICES_MAP !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_INITIALISE",ERR,ERROR,*999)

    VARIABLE_TO_EQUATIONS_MATRICES_MAP%VARIABLE_INDEX=0
    VARIABLE_TO_EQUATIONS_MATRICES_MAP%VARIABLE_TYPE=0
    NULLIFY(VARIABLE_TO_EQUATIONS_MATRICES_MAP%VARIABLE)
    VARIABLE_TO_EQUATIONS_MATRICES_MAP%NUMBER_OF_EQUATIONS_MATRICES=0
    
    CALL EXITS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_FINALISE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,row_idx,variable_type

    CALL ENTERS("EQUATIONS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ALLOCATED(EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES)) DEALLOCATE(EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS)) THEN
        DO variable_type=1,SIZE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS,1)
          CALL EQUATIONS_MAPPING_VARIABLE_TO_EQUATIONS_MATRICES_MAP_FINALISE(EQUATIONS_MAPPING% &
            & VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type),ERR,ERROR,*999)
        ENDDO !variable_type
        DEALLOCATE(EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS)        
      ENDIF
      IF(ALLOCATED(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS)) THEN
        DO matrix_idx=1,SIZE(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS,1)
          CALL EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VARIABLE_MAP_FINALISE(EQUATIONS_MAPPING% &
            & EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS)
      ENDIF
      CALL EQUATIONS_MAPPING_SOURCE_MAPPINGS_FINALISE(EQUATIONS_MAPPING%SOURCE_MAPPINGS,ERR,ERROR,*999)
      IF(ALLOCATED(EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS)) THEN
        DO row_idx=1,SIZE(EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS,1)
          CALL EQUATIONS_MAPPING_EQUATIONS_ROW_TO_VARIABLE_MAP_FINALISE(EQUATIONS_MAPPING% &
            & EQUATIONS_ROW_TO_VARIABLES_MAPS(row_idx),ERR,ERROR,*999)
        ENDDO !row_idx
        DEALLOCATE(EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS)
      ENDIF
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(EQUATIONS_MAPPING%ROW_DOFS_MAPPING,ERR,ERROR,*999)
      CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS_MAPPING)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_INITIALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to initialise the equations mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MAPPING)) THEN
        CALL FLAG_ERROR("Equations mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS%EQUATIONS_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations equations mapping.",ERR,ERROR,*999)
        EQUATIONS%EQUATIONS_MAPPING%EQUATIONS=>EQUATIONS
        EQUATIONS%EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED=.FALSE.
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%EQUATIONS_MATRICES)
        EQUATIONS%EQUATIONS_MAPPING%NUMBER_OF_ROWS=0
        EQUATIONS%EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS=0
        EQUATIONS%EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES=0
        EQUATIONS%EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES=0
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%SOURCE_MAPPINGS)
        EQUATIONS%EQUATIONS_MAPPING%RHS_VARIABLE_TYPE=0
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%RHS_VARIABLE)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%RHS_VARIABLE_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%ROW_DOFS_MAPPING)
        CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE(EQUATIONS%EQUATIONS_MAPPING,ERR,ERROR,*999)        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_FINALISE(EQUATIONS%EQUATIONS_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the equations matrices in an equation set. 
  SUBROUTINE EQUATIONS_MAPPING_MATRIX_COEFFICIENTS_SET(EQUATIONS_MAPPING,MATRIX_COEFFICIENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping.
    REAL(DP), INTENT(IN) :: MATRIX_COEFFICIENTS(:) !<The matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_MATRIX_COEFFICIENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping is finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN          
          IF(SIZE(MATRIX_COEFFICIENTS,1)==EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=MATRIX_COEFFICIENTS          
          ELSE
            LOCAL_ERROR="Invalid size of matrix coefficeints. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_COEFFICIENTS,1),"*",ERR,ERROR))// &
              & ") must match the number of equations matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MAPPING_MATRIX_COEFFICIENTS_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_MATRIX_COEFFICIENTS_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_MATRIX_COEFFICIENTS_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_MATRIX_COEFFICIENTS_SET

  !
  !================================================================================================================================
  !

END MODULE EQUATIONS_MAPPING_ROUTINES
