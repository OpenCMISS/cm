!> \file
!> $Id$
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
  USE INPUT_OUTPUT
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
    & EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFICIENTS_SET,EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET, &
    & EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET,EQUATIONS_MAPPING_RESIDUAL_COEFFICIENT_SET, &
    & EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET,EQUATIONS_MAPPING_RHS_COEFFICIENT_SET, &
    & EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET,EQUATIONS_MAPPING_SOURCE_COEFFICIENT_SET, &
    & EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET
  
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
    INTEGER(INTG) :: column_idx,dof_idx,matrix_idx,NUMBER_OF_ROWS,NUMBER_OF_GLOBAL_ROWS,row_idx, &
      & TOTAL_NUMBER_OF_ROWS,variable_idx,variable_type
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,SOURCE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,SOURCE_VARIABLE
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
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE/=0) THEN
                IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
                  SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
                  IF(.NOT.ASSOCIATED(SOURCE_FIELD)) THEN
                    CALL FLAG_ERROR("Source field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF                  
              !Calculate the number of rows in the equations set
              IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                !Linear equations set
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>=1) THEN
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & MATRIX_VARIABLE_TYPES(1))%PTR
                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                    NUMBER_OF_ROWS=DEPENDENT_VARIABLE%NUMBER_OF_DOFS
                    TOTAL_NUMBER_OF_ROWS=DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    EQUATIONS_MAPPING%ROW_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                    IF(ASSOCIATED(EQUATIONS_MAPPING%ROW_DOFS_MAPPING)) THEN
                      NUMBER_OF_GLOBAL_ROWS=EQUATIONS_MAPPING%ROW_DOFS_MAPPING%NUMBER_OF_GLOBAL
                    ELSE
                      CALL FLAG_ERROR("Dependent variable domain mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The dependent variable mapped to the first linear matrix is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The number of linear equations matrices is zero of less for a linear equations set.", &
                    & ERR,ERROR,*999)
                ENDIF
                !Check that the number of rows is consistent across the remaining linear matrices
                DO matrix_idx=2,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & MATRIX_VARIABLE_TYPES(matrix_idx))%PTR
                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                    IF(DEPENDENT_VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The number of rows in equations matrix number 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the number of rows in equations matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS/=TOTAL_NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The total number of rows in equations matrix number 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the total number of rows in equations matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The dependent variable mapped to linear matrix number "// &
                      & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx                  
                !Check that the number of rows are consistent with the RHS vector if it exists
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) THEN
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & MATRIX_VARIABLE_TYPES(1))%PTR
                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                    IF(DEPENDENT_VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The number of rows in the equations matrix number 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the number of rows in the RHS vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS/=TOTAL_NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The total number of rows in the equations matrix number 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the total number of rows in the RHS vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The dependent variable mapped to the RHS vector is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                !Check that the number of rows are consistent with the source vector if it exists
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE/=0) THEN
                  SOURCE_VARIABLE=>SOURCE_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & SOURCE_VARIABLE_TYPE)%PTR
                  IF(ASSOCIATED(SOURCE_VARIABLE)) THEN
                    IF(SOURCE_VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The number of rows in equations matrix number 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the number of rows in the source vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE%NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(SOURCE_VARIABLE%TOTAL_NUMBER_OF_DOFS/=TOTAL_NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The total number of rows in equations matrix number 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the total number of rows in the source vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE%TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The source variable mapped to the source vector is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                !Non linear equations set
                DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                  & RESIDUAL_VARIABLE_TYPE)%PTR
                IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                  NUMBER_OF_ROWS=DEPENDENT_VARIABLE%NUMBER_OF_DOFS
                  TOTAL_NUMBER_OF_ROWS=DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                  EQUATIONS_MAPPING%ROW_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  IF(ASSOCIATED(EQUATIONS_MAPPING%ROW_DOFS_MAPPING)) THEN
                    NUMBER_OF_GLOBAL_ROWS=EQUATIONS_MAPPING%ROW_DOFS_MAPPING%NUMBER_OF_GLOBAL
                  ELSE
                    CALL FLAG_ERROR("Dependent variable domain mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Residual variable is not associated.",ERR,ERROR,*999)
                ENDIF
                !Check that the number of rows is consistent across the linear matrices (if any)
                DO matrix_idx=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & MATRIX_VARIABLE_TYPES(matrix_idx))%PTR
                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                    IF(DEPENDENT_VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The number of rows in the residual vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the number of rows in equations matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS/=TOTAL_NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The total number of rows in the residual vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the total number of rows in equations matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The dependent variable mapped to linear matrix number "// &
                      & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
                !Check that the number of rows are consistent with the RHS vector if it exists
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) THEN
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & RHS_VARIABLE_TYPE)%PTR
                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                    IF(DEPENDENT_VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The number of rows in the residual vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the number of rows in the RHS vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS/=TOTAL_NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The total number of rows in the residual vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the total number of rows in the RHS vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The dependent variable mapped to the RHS vector is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                !Check that the number of rows are consistent with the source vector if it exists
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE/=0) THEN
                  SOURCE_VARIABLE=>SOURCE_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & SOURCE_VARIABLE_TYPE)%PTR
                  IF(ASSOCIATED(SOURCE_VARIABLE)) THEN
                    IF(SOURCE_VARIABLE%NUMBER_OF_DOFS/=NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The number of rows in the residual vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the number of rows in the source vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE%NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(SOURCE_VARIABLE%TOTAL_NUMBER_OF_DOFS/=TOTAL_NUMBER_OF_ROWS) THEN
                      LOCAL_ERROR="Invalid equations set up. The total number of rows in the residual vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(TOTAL_NUMBER_OF_ROWS,"*",ERR,ERROR))// &
                        & ") does not match the total number of rows in the source vector ("// &
                        & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE%TOTAL_NUMBER_OF_DOFS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The source variable mapped to the source vector is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDIF
              EQUATIONS_MAPPING%NUMBER_OF_ROWS=NUMBER_OF_ROWS
              EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS=TOTAL_NUMBER_OF_ROWS
              EQUATIONS_MAPPING%NUMBER_OF_GLOBAL_ROWS=NUMBER_OF_GLOBAL_ROWS 
              !Calculate linear mappings
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN                  
                CALL EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*999)
                LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                  LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  !Allocate and initialise the variable type maps
                  ALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping variable to equations map.",ERR,ERROR,*999)
                  DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE( &
                      & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type),ERR,ERROR,*999)
                    LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE_INDEX=variable_type
                    LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE_TYPE=variable_type
                    LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE=>DEPENDENT_FIELD% &
                      & VARIABLE_TYPE_MAP(variable_type)%PTR
                  ENDDO !variable_type
                  !Calculate the number of variable type maps and initialise
                  DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    variable_type=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)
                    LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES= &
                      & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES+1
                  ENDDO !matrix_idx
                  IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) LINEAR_MAPPING% &
                    & VAR_TO_EQUATIONS_MATRICES_MAPS(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE)% &
                    & NUMBER_OF_LINEAR_EQUATIONS_MATRICES=-1                    
                  LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES=0
                  !Allocate and initialise the variable to equations matrices maps
                  DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                    IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                      IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                        & NUMBER_OF_LINEAR_EQUATIONS_MATRICES==-1) THEN
!!TODO: check if this can be removed and just allocate those variables that are actually used
                        ALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP( &
                          & DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to equations matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP=0
                        LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES=LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES+1
                      ELSE IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                        & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
                        ALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%EQUATIONS_MATRIX_NUMBERS( &
                          & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                          & NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                        IF(ERR/=0) &
                          & CALL FLAG_ERROR("Could not allocate variable to equations matrices maps equations matrix numbers.", &
                          & ERR,ERROR,*999)
                        ALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_COLUMNS_MAPS( &
                          & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                          & NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to equations matrices maps dof to columns map.", &
                          & ERR,ERROR,*999)                
                        LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%EQUATIONS_MATRIX_NUMBERS=0
                        LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                        DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)==variable_type) THEN
                            LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES= &
                              & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES+1
                            LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%EQUATIONS_MATRIX_NUMBERS( &
                              & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES)=matrix_idx
                            ALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_COLUMNS_MAPS( &
                              & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES)%COLUMN_DOF(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to columns map column dof.", &
                              & ERR,ERROR,*999)
                            DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                              !1-1 mapping for now
                              column_idx=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                              LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_COLUMNS_MAPS( &
                                & LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                & NUMBER_OF_LINEAR_EQUATIONS_MATRICES)%COLUMN_DOF(dof_idx)=column_idx
                            ENDDO !dof_idx
                          ENDIF
                        ENDDO !matrix_idx
                        ALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP( &
                          & DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to equations matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          !1-1 mappings for now.
                          row_idx=dof_idx
                          LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%DOF_TO_ROWS_MAP(dof_idx)=row_idx
                        ENDDO !dof_idx
                        LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES=LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES+1
                      ENDIF
                    ENDIF
                  ENDDO !variable_type
                  !Allocate and initialise the variable types
                  ALLOCATE(LINEAR_MAPPING%MATRIX_VARIABLE_TYPES(LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping matrix variable types.",ERR,ERROR,*999)
                  LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES=0
                  DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    IF(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
                      LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES=LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES+1
                      LINEAR_MAPPING%MATRIX_VARIABLE_TYPES(LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES)=variable_type
                    ENDIF
                  ENDDO !variable_type
                  !Allocate and initialise the equations matrix to variable maps types
                  ALLOCATE(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES), &
                    & STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping equations matrix to variable maps.", &
                    & ERR,ERROR,*999)
                  !Create the individual matrix maps and column maps
                  DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    variable_type=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)
                    DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                    CALL EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE(LINEAR_MAPPING% &
                      & EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx),ERR,ERROR,*999)
                    LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%MATRIX_NUMBER=matrix_idx
                    LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=variable_type
                    LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE=>DEPENDENT_VARIABLE
                    LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_COLUMNS=DEPENDENT_VARIABLE% &
                      & DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                    LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT=EQUATIONS_MAPPING% &
                      CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                    ALLOCATE(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_TO_DOF_MAP( &
                      & DEPENDENT_VARIABLE%DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)                  
                    LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_TO_DOF_MAP=0
                    DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                      !1-1 mapping for now
                      column_idx=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                      LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_TO_DOF_MAP(column_idx)=dof_idx
                    ENDDO !dof_idx
                    LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING=> &
                      & DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  ENDDO !matrix_idx
                  !Allocate the row mappings
                  ALLOCATE(LINEAR_MAPPING%EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS, &
                    & LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to variable dof maps.",ERR,ERROR,*999)
                  !Set up the row mappings
                  DO variable_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                    DO row_idx=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
                      !1-1 mapping for now
                      LINEAR_MAPPING%EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(row_idx,variable_idx)=row_idx
                    ENDDO !row_idx
                  ENDDO !variable_idx
                ELSE
                  CALL FLAG_ERROR("Linear mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              !Calculate non-linear mappings
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE/=0) THEN                  
                CALL EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*999)
                NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                  CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE(NONLINEAR_MAPPING% &
                    & VAR_TO_JACOBIAN_MAP,ERR,ERROR,*999)
                  NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%VARIABLE_TYPE=EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & RESIDUAL_VARIABLE_TYPE
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & RESIDUAL_VARIABLE_TYPE)%PTR
                  NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%VARIABLE=>DEPENDENT_VARIABLE
                  !Allocate and set dof to Jacobian columns map
                  ALLOCATE(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(DEPENDENT_VARIABLE% &
                    & TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to Jacobian map dof to columns map.",ERR,ERROR,*999)
                  DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    column_idx=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                    NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP(dof_idx)=column_idx
                  ENDDO !dof_idx
                  !Allocate and set dof to Jacobian rows map
                  ALLOCATE(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_ROWS_MAP(DEPENDENT_VARIABLE% &
                    & TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to Jacobian map dof to columns map.",ERR,ERROR,*999)
                  DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    row_idx=dof_idx
                    NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_ROWS_MAP(dof_idx)=row_idx
                  ENDDO !dof_idx
                  CALL EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE(NONLINEAR_MAPPING% &
                    & JACOBIAN_TO_VAR_MAP,ERR,ERROR,*999)
                  NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE_TYPE=EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & RESIDUAL_VARIABLE_TYPE
                  NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE=>DEPENDENT_VARIABLE
                  NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS=DEPENDENT_VARIABLE%DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                  NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT=EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & RESIDUAL_COEFFICIENT
                  ALLOCATE(NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP( &
                    & DEPENDENT_VARIABLE%DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
                  NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP=0
                  DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    column_idx=DEPENDENT_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                    NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP(column_idx)=dof_idx
                  ENDDO !dof_idx
                  NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%COLUMN_DOFS_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  !Set up the residual vector mapping
                  NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE
                  NONLINEAR_MAPPING%RESIDUAL_VARIABLE=>DEPENDENT_VARIABLE
                  NONLINEAR_MAPPING%RESIDUAL_VARIABLE_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  !Set up the row mappings
                  ALLOCATE(NONLINEAR_MAPPING%RESIDUAL_DOF_TO_EQUATIONS_ROW_MAP(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate residual dof to equations row map.",ERR,ERROR,*999)
                  ALLOCATE(NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP(TOTAL_NUMBER_OF_ROWS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to residual dof map.",ERR,ERROR,*999)
                  DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    row_idx=dof_idx
                    NONLINEAR_MAPPING%RESIDUAL_DOF_TO_EQUATIONS_ROW_MAP(dof_idx)=row_idx
                  ENDDO !dof_idx
                  DO row_idx=1,TOTAL_NUMBER_OF_ROWS
                    !1-1 mapping for now
                    dof_idx=row_idx
                    NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP(row_idx)=dof_idx
                  ENDDO !row_idx
                ELSE
                  CALL FLAG_ERROR("Nonlinear mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              !Calculate RHS mappings
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) THEN                  
                CALL EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*999)
                RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                IF(ASSOCIATED(RHS_MAPPING)) THEN
                  RHS_MAPPING%RHS_VARIABLE_TYPE=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE
                  DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & RHS_VARIABLE_TYPE)%PTR
                  RHS_MAPPING%RHS_VARIABLE=>DEPENDENT_VARIABLE
                  RHS_MAPPING%RHS_VARIABLE_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                  RHS_MAPPING%RHS_COEFFICIENT=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT
                  !Allocate and set up the row mappings
                  ALLOCATE(RHS_MAPPING%RHS_DOF_TO_EQUATIONS_ROW_MAP(DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate rhs dof to equations row map.",ERR,ERROR,*999)
                  ALLOCATE(RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(TOTAL_NUMBER_OF_ROWS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to dof map.",ERR,ERROR,*999)
                  DO dof_idx=1,DEPENDENT_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    row_idx=dof_idx
                    RHS_MAPPING%RHS_DOF_TO_EQUATIONS_ROW_MAP(dof_idx)=row_idx
                  ENDDO !dof_idx
                  DO row_idx=1,TOTAL_NUMBER_OF_ROWS
                    !1-1 mapping for now
                    dof_idx=row_idx
                    RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(row_idx)=dof_idx
                  ENDDO !row_idx
                ELSE
                  CALL FLAG_ERROR("RHS mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              !Calcuate the source mappings
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE/=0) THEN                  
                CALL EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*999)
                SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
                IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                  SOURCE_MAPPING%SOURCE_VARIABLE_TYPE=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE
                  SOURCE_VARIABLE=>SOURCE_FIELD%VARIABLE_TYPE_MAP(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                    & SOURCE_VARIABLE_TYPE)%PTR
                  SOURCE_MAPPING%SOURCE_VARIABLE=>SOURCE_VARIABLE
                  SOURCE_MAPPING%SOURCE_VARIABLE_MAPPING=>SOURCE_VARIABLE%DOMAIN_MAPPING
                  SOURCE_MAPPING%SOURCE_COEFFICIENT=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_COEFFICIENT
                  !Allocate and set up the row mappings
                  ALLOCATE(SOURCE_MAPPING%SOURCE_DOF_TO_EQUATIONS_ROW_MAP(SOURCE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate source dof to equations row map.",ERR,ERROR,*999)
                  ALLOCATE(SOURCE_MAPPING%EQUATIONS_ROW_TO_SOURCE_DOF_MAP(TOTAL_NUMBER_OF_ROWS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations row to source map.",ERR,ERROR,*999)
                  DO dof_idx=1,SOURCE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    row_idx=dof_idx
                    SOURCE_MAPPING%SOURCE_DOF_TO_EQUATIONS_ROW_MAP(dof_idx)=row_idx
                  ENDDO !dof_idx
                  DO row_idx=1,TOTAL_NUMBER_OF_ROWS
                    !1-1 mapping for now
                    dof_idx=row_idx
                    SOURCE_MAPPING%EQUATIONS_ROW_TO_SOURCE_DOF_MAP(row_idx)=dof_idx
                  ENDDO !row_idx
                ELSE
                  CALL FLAG_ERROR("Source mapping is not associated.",ERR,ERROR,*999)
                ENDIF
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

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Equations mappings:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",EQUATIONS_MAPPING%NUMBER_OF_ROWS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total umber of rows = ",EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global rows = ",EQUATIONS_MAPPING%NUMBER_OF_GLOBAL_ROWS, &
        & ERR,ERROR,*999)
      LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
      IF(ASSOCIATED(LINEAR_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Linear mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear equations matrices = ",LINEAR_MAPPING% &
          & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear matrix variables = ",LINEAR_MAPPING% &
          & NUMBER_OF_LINEAR_MATRIX_VARIABLES,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES,4,4, &
          & LINEAR_MAPPING%MATRIX_VARIABLE_TYPES,'("    Matrix variable types :",4(X,I12))','(27X,4(X,I12))', &
          & ERR,ERROR,*999) 
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",ERR,ERROR,*999)
        DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ",variable_type,ERR,ERROR,*999)
          IF(ASSOCIATED(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE)) THEN
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",LINEAR_MAPPING% &
              & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%VARIABLE%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of linear equations matrices = ",LINEAR_MAPPING% &
              & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,4,4,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
              & EQUATIONS_MATRIX_NUMBERS,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',ERR,ERROR,*999) 
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      DOF to column maps :",ERR,ERROR,*999)
            DO matrix_idx=1,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix number : ",matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                & variable_type)%VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                & variable_type)%DOF_TO_COLUMNS_MAPS(matrix_idx)%COLUMN_DOF, &
                & '("        Column numbers :",5(X,I13))','(24X,5(X,I13))',ERR,ERROR,*999) 
            ENDDO !matrix_idx
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
              & VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
              & DOF_TO_ROWS_MAP,'("      DOF to row maps  :",5(X,I13))','(24X,5(X,I13))',ERR,ERROR,*999)
          ENDIF
        ENDDO !variable_type
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",ERR,ERROR,*999)
        DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrix_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",LINEAR_MAPPING% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",LINEAR_MAPPING% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",LINEAR_MAPPING% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)% &
            & NUMBER_OF_COLUMNS,5,5,LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_TO_DOF_MAP, &
            & '("        Column to DOF maps :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999) 
        ENDDO !matrix_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",ERR,ERROR,*999)
        DO variable_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Variable number : ",variable_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS,5,5, &
            & LINEAR_MAPPING%EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(:,variable_idx), &
            & '("    Row to DOF maps :",5(X,I13))','(21X,5(X,I13))',ERR,ERROR,*999) 
        ENDDO !variable_idx
      ENDIF
      NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
      IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Nonlinear mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Residual variable type = ",NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of residual DOFs = ",NONLINEAR_MAPPING%RESIDUAL_VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Residual coefficient = ",NONLINEAR_MAPPING%RESIDUAL_COEFFICIENT, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Residual row mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NONLINEAR_MAPPING%RESIDUAL_VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5, &
          & NONLINEAR_MAPPING%RESIDUAL_DOF_TO_EQUATIONS_ROW_MAP,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))', &
          & ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS,5,5, &
          & NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))', &
          & ERR,ERROR,*999) 
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Jacobian mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to Jacobian mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian variable type = ",NONLINEAR_MAPPING% &
          & VAR_TO_JACOBIAN_MAP%VARIABLE_TYPE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of Jacobain DOFs = ",NONLINEAR_MAPPING% &
          & VAR_TO_JACOBIAN_MAP%VARIABLE%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,5,5,NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP, &
          & '("      DOF to column map :",5(X,I13))','(26X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,5,5,NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP%DOF_TO_ROWS_MAP, &
          & '("      DOF to row map    :",5(X,I13))','(26X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Jacobian to variable mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian variable type = ",NONLINEAR_MAPPING% &
          & JACOBIAN_TO_VAR_MAP%VARIABLE_TYPE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",NONLINEAR_MAPPING% &
          & JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian coefficient = ",NONLINEAR_MAPPING% &
          & JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS, &
          & 5,5,NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP, &
          & '("      Column to DOF map :",5(X,I13))','(26X,5(X,I13))',ERR,ERROR,*999) 
      ENDIF
      RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
      IF(ASSOCIATED(RHS_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  RHS mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    RHS variable type = ",RHS_MAPPING%RHS_VARIABLE_TYPE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of RHS DOFs = ",RHS_MAPPING%RHS_VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    RHS coefficient = ",RHS_MAPPING%RHS_COEFFICIENT,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,RHS_MAPPING%RHS_VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5, &
          & RHS_MAPPING%RHS_DOF_TO_EQUATIONS_ROW_MAP,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS,5,5, &
          & RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))',ERR,ERROR,*999) 
       ENDIF
      SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
      IF(ASSOCIATED(SOURCE_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Source mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Source variable type = ",SOURCE_MAPPING%SOURCE_VARIABLE_TYPE, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of source DOFs = ",SOURCE_MAPPING%SOURCE_VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Source coefficient = ",SOURCE_MAPPING%SOURCE_COEFFICIENT,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOURCE_MAPPING%SOURCE_VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5, &
          & SOURCE_MAPPING%SOURCE_DOF_TO_EQUATIONS_ROW_MAP,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))', &
          & ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS,5,5, &
          & SOURCE_MAPPING%EQUATIONS_ROW_TO_SOURCE_DOF_MAP,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))', &
          & ERR,ERROR,*999) 
      ENDIF
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
              !Allocate and initialise the create values cache
              ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping create values cache",ERR,ERROR,*999)
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENT=1.0_DP
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=1.0_DP
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_COEFFICIENT=1.0_DP
              !Set the default equations mapping in the create values cache
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
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=1
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                      & VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE=DEPENDENT_FIELD% &
                      & VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                      & VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=DEPENDENT_FIELD%NUMBER_OF_VARIABLES-1
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                      & VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=DEPENDENT_FIELD%NUMBER_OF_VARIABLES-2
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE=DEPENDENT_FIELD% &
                      & VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                      & VARIABLE_TYPE_MAP(FIELD_NORMAL_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDIF
               ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING% &
                & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping create values cache matrix variable types", &
                & ERR,ERROR,*999)
              ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(EQUATIONS_MAPPING% &
                & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping create values cache matrix coefficients", &
                & ERR,ERROR,*999)
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=0
              IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR)) THEN
                  EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1)=DEPENDENT_FIELD% &
                    & VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                ELSE
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                ENDIF
                DO matrix_idx=2,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)= &
                      & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
              ELSE
                DO matrix_idx=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)= &
                      & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR%VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  ENDIF
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

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE(EQUATIONS_JACOBIAN_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_JACOBIAN_TO_VAR_MAP_TYPE) :: EQUATIONS_JACOBIAN_TO_VAR_MAP !<The equations Jacobian to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(EQUATIONS_JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP)) &
      & DEALLOCATE(EQUATIONS_JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP)
    
    CALL EXITS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map.
  SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE(EQUATIONS_JACOBIAN_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_JACOBIAN_TO_VAR_MAP_TYPE) :: EQUATIONS_JACOBIAN_TO_VAR_MAP !<The equations Jacobian to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE",ERR,ERROR,*999)
    
    EQUATIONS_JACOBIAN_TO_VAR_MAP%VARIABLE_TYPE=0
    NULLIFY(EQUATIONS_JACOBIAN_TO_VAR_MAP%VARIABLE)
    NULLIFY(EQUATIONS_JACOBIAN_TO_VAR_MAP%JACOBIAN)
    EQUATIONS_JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS=0
    EQUATIONS_JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT=0
    NULLIFY(EQUATIONS_JACOBIAN_TO_VAR_MAP%COLUMN_DOFS_MAPPING)    
    
    CALL EXITS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise an equations matrix to variable maps and deallocate all memory.
  SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE(EQUATIONS_MATRIX_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TO_VAR_MAP_TYPE) :: EQUATIONS_MATRIX_TO_VAR_MAP !<The equations matrix to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_MATRIX_TO_VAR_MAP%COLUMN_TO_DOF_MAP)) &
      & DEALLOCATE(EQUATIONS_MATRIX_TO_VAR_MAP%COLUMN_TO_DOF_MAP)
    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an equations matrix to variable maps.
  SUBROUTINE EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE(EQUATIONS_MATRIX_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TO_VAR_MAP_TYPE) :: EQUATIONS_MATRIX_TO_VAR_MAP !<The equations matrix to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_MATRIX_TO_VAR_MAP%MATRIX_NUMBER=0
    EQUATIONS_MATRIX_TO_VAR_MAP%VARIABLE_TYPE=0
    NULLIFY(EQUATIONS_MATRIX_TO_VAR_MAP%VARIABLE)
    EQUATIONS_MATRIX_TO_VAR_MAP%NUMBER_OF_COLUMNS=0
    EQUATIONS_MATRIX_TO_VAR_MAP%MATRIX_COEFFICIENT=1.0_DP !Matrices in an equation set are added by default
    NULLIFY(EQUATIONS_MATRIX_TO_VAR_MAP%COLUMN_DOFS_MAPPING)
    
    CALL EXITS("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE

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

    CALL ENTERS("EQUATIONS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(EQUATIONS_MAPPING%ROW_DOFS_MAPPING,ERR,ERROR,*999)
      CALL EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%LINEAR_MAPPING,ERR,ERROR,*999)
      CALL EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%NONLINEAR_MAPPING,ERR,ERROR,*999)
      CALL EQUATIONS_MAPPING_RHS_MAPPING_FINALISE(EQUATIONS_MAPPING%RHS_MAPPING,ERR,ERROR,*999)      
      CALL EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE(EQUATIONS_MAPPING%SOURCE_MAPPING,ERR,ERROR,*999)      
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
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%ROW_DOFS_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%LINEAR_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%SOURCE_MAPPING)
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

  !>Finalises the equations mapping linear mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE(LINEAR_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING !<A pointer to the linear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,variable_type
 
    CALL ENTERS("EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
      IF(ALLOCATED(LINEAR_MAPPING%MATRIX_VARIABLE_TYPES)) DEALLOCATE(LINEAR_MAPPING%MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS)) THEN
        DO variable_type=1,SIZE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS,1)
          CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE(LINEAR_MAPPING% &
            & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type),ERR,ERROR,*999)
        ENDDO !variable_type
        DEALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS)        
      ENDIF
      IF(ALLOCATED(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS)) THEN
        DO matrix_idx=1,SIZE(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS,1)
          CALL EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE(LINEAR_MAPPING% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS)
      ENDIF
      IF(ALLOCATED(LINEAR_MAPPING%EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS)) DEALLOCATE(LINEAR_MAPPING%EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS)
      DEALLOCATE(LINEAR_MAPPING)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping linear mapping
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the linear mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%LINEAR_MAPPING)) THEN
        CALL FLAG_ERROR("Equations mapping linear mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%LINEAR_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping linear mapping.",ERR,ERROR,*999)
        EQUATIONS_MAPPING%LINEAR_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING       
        EQUATIONS_MAPPING%LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
        EQUATIONS_MAPPING%LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%LINEAR_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the linear equations matrices in an equation set. 
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFICIENTS_SET(EQUATIONS_MAPPING,MATRIX_COEFFICIENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping.
    REAL(DP), INTENT(IN) :: MATRIX_COEFFICIENTS(:) !<The matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFICIENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping is finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN          
          IF(SIZE(MATRIX_COEFFICIENTS,1)==EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=MATRIX_COEFFICIENTS          
          ELSE
            LOCAL_ERROR="Invalid size of matrix coefficeints. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_COEFFICIENTS,1),"*",ERR,ERROR))// &
              & ") must match the number of equations matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFICIENTS_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFICIENTS_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFICIENTS_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFICIENTS_SET

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variables and the linear equations matrices
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set the number of matrices for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_LINEAR_EQUATIONS_MATRICES !<The number of linear equations matrices for the mapping.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_VARIABLE_TYPES(:)
    REAL(DP), ALLOCATABLE :: OLD_MATRIX_COEFFICIENTS(:)
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET",ERR,ERROR,*999)

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
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE==0) THEN
                  IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<1.OR. &
                    & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    LOCAL_ERROR="The requested number of matrices ("// &
                      & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",ERR,ERROR))// &
                      & ") is invalid. For linear problems without a equations set RHS the number must be between >= 1 and <= "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE            
                  IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<1.OR. &
                    & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    LOCAL_ERROR="The requested number of matrices ("// &
                      & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",ERR,ERROR))// &
                      & ") is invalid. For linear problems with a equations set RHS the number must be between >= 1 and < "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<0.OR. &
                  & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>FIELD_NUMBER_OF_VARIABLE_TYPES-2) THEN
                  LOCAL_ERROR="The requested number of matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",ERR,ERROR))// &
                    & ") is invalid. For nonlinear problems the number must be between >= 0 and <= "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES-2,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
              !If we need to reallocate and reset all the create_values cache arrays and change the number of matrices
              IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES/= &
                & EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
                ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types",ERR,ERROR,*999)
                ALLOCATE(OLD_MATRIX_COEFFICIENTS(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix coefficients",ERR,ERROR,*999)
                OLD_MATRIX_VARIABLE_TYPES=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
                OLD_MATRIX_COEFFICIENTS=EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS
                DEALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
                DEALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types",ERR,ERROR,*999)
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix coefficients",ERR,ERROR,*999)
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
                  IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES> &
                    & EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES)=OLD_MATRIX_VARIABLE_TYPES
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES+1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)=OLD_MATRIX_VARIABLE_TYPES(1)
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES)=OLD_MATRIX_COEFFICIENTS
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES+1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)=OLD_MATRIX_COEFFICIENTS(1)
                  ELSE
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)= &
                      & OLD_MATRIX_VARIABLE_TYPES(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)= &
                      & OLD_MATRIX_COEFFICIENTS(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)
                  ENDIF
                ELSE
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=0
                    IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_LINEAR) THEN
                      IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR)) THEN
                        EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1)=DEPENDENT_FIELD% &
                          & VARIABLE_TYPE_MAP(FIELD_STANDARD_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                      ELSE
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      ENDIF
                      DO matrix_idx=2,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                        IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR)) THEN
                          EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)= &
                            & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR%VARIABLE_TYPE
                        ELSE
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO !matrix_idx
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
                    ELSE
                      DO matrix_idx=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                        IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR)) THEN
                          EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)= &
                            & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR%VARIABLE_TYPE
                        ELSE
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO !matrix_idx
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES
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
    
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)    
    IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)    
    CALL ERRORS("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variable types and the linear equations matrices
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,MATRIX_VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(IN) :: MATRIX_VARIABLE_TYPES(:) !<The matrix variable types to map to each linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(SIZE(MATRIX_VARIABLE_TYPES,1)==EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  !Check input values
                  DO matrix_idx=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
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
                      & EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE) THEN
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
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping nonlinear mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE(NONLINEAR_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING !<A pointer to the nonlinear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
      CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP,ERR,ERROR,*999)
      CALL EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE(NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP,ERR,ERROR,*999)
      IF(ALLOCATED(NONLINEAR_MAPPING%RESIDUAL_DOF_TO_EQUATIONS_ROW_MAP)) &
        & DEALLOCATE(NONLINEAR_MAPPING%RESIDUAL_DOF_TO_EQUATIONS_ROW_MAP)
      IF(ALLOCATED(NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP)) &
        & DEALLOCATE(NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP)
      DEALLOCATE(NONLINEAR_MAPPING)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping nonlinear mapping
  SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the nonlinear mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%NONLINEAR_MAPPING)) THEN
        CALL FLAG_ERROR("Equations mapping nonlinear mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%NONLINEAR_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping nonlinear mapping.",ERR,ERROR,*999)
        EQUATIONS_MAPPING%NONLINEAR_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING
        CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE(EQUATIONS_MAPPING%NONLINEAR_MAPPING% &
          & VAR_TO_JACOBIAN_MAP,ERR,ERROR,*999)
        CALL EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE(EQUATIONS_MAPPING%NONLINEAR_MAPPING% &
          & JACOBIAN_TO_VAR_MAP,ERR,ERROR,*999)
        EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE=0
        NULLIFY(EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLE)
        NULLIFY(EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLE_MAPPING)
        EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%NONLINEAR_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set residual vector.
  SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_COEFFICIENT_SET(EQUATIONS_MAPPING,RESIDUAL_COEFFICIENT,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: RESIDUAL_COEFFICIENT!<The coefficient applied to the equations set residual vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_RESIDUAL_COEFFICIENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE/=0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENT=RESIDUAL_COEFFICIENT
          ELSE
            CALL FLAG_ERROR("The equations mapping residual variable type has not been set.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_RESIDUAL_COEFFICIENT_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_RESIDUAL_COEFFICIENT_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_RESIDUAL_COEFFICIENT_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_COEFFICIENT_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set residual vector.
  SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,RESIDUAL_VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: RESIDUAL_VARIABLE_TYPE !<The variable type associated with the equations set residual vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_NONLINEAR) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN                 
                !Check the residual variable type is not being used for the RHS vector
                  IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE==RESIDUAL_VARIABLE_TYPE) THEN
                    LOCAL_ERROR="The specified residual variable type of "// &
                      & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " is the same as the variable type for the RHS vector."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  !Check the residual variable number is defined on the dependent field
                  IF(RESIDUAL_VARIABLE_TYPE>=1.AND.RESIDUAL_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(RESIDUAL_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE=RESIDUAL_VARIABLE_TYPE
                    ELSE
                      LOCAL_ERROR="The specified residual variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " is not defined on the dependent field."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The specified residual variable type of "// &
                      & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " is invalid. The number must either be zero or >= 1 and <= "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The equations set is not a nonlinear equations set.",ERR,ERROR,*999)
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
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set RHS vector.
  SUBROUTINE EQUATIONS_MAPPING_RHS_COEFFICIENT_SET(EQUATIONS_MAPPING,RHS_COEFFICIENT,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: RHS_COEFFICIENT!<The coefficient applied to the equations set RHS vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_RHS_COEFFICIENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=RHS_COEFFICIENT
          ELSE
            CALL FLAG_ERROR("The equations mapping RHS variable type has not been set.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_RHS_COEFFICIENT_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_RHS_COEFFICIENT_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_RHS_COEFFICIENT_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_COEFFICIENT_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping RHS mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_FINALISE(RHS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_MAPPING_RHS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(RHS_MAPPING)) THEN
      IF(ALLOCATED(RHS_MAPPING%RHS_DOF_TO_EQUATIONS_ROW_MAP)) DEALLOCATE(RHS_MAPPING%RHS_DOF_TO_EQUATIONS_ROW_MAP)
      IF(ALLOCATED(RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP)) DEALLOCATE(RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP)
      DEALLOCATE(RHS_MAPPING)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_RHS_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_RHS_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping RHS mapping
  SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%RHS_MAPPING)) THEN
        CALL FLAG_ERROR("Equations mapping RHS mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%RHS_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping RHS mapping.",ERR,ERROR,*999)
        EQUATIONS_MAPPING%RHS_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING        
        EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE=0
        NULLIFY(EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE)
        NULLIFY(EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_MAPPING)
        EQUATIONS_MAPPING%RHS_MAPPING%RHS_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_RHS_MAPPING_FINALISE(EQUATIONS_MAPPING%RHS_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE

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
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(RHS_VARIABLE_TYPE==0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=0
          ELSE
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  !Check that the given RHS variable type is not already being used for a linear matrix
                  DO matrix_idx=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)==RHS_VARIABLE_TYPE) THEN
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " is the same as the variable type for linear equations matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !matrix_idx
                  !If we have a non-linear problem check the RHS variable type is not being used for the residual vector
                  IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_NONLINEAR) THEN
                    IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPE==RHS_VARIABLE_TYPE) THEN
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " is the same as the variable type for the residual vector."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                  !Check the RHS variable number is defined on the dependent field
                  IF(RHS_VARIABLE_TYPE>=1.AND.RHS_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(RHS_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=RHS_VARIABLE_TYPE
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
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated",ERR,ERROR,*999)
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

  !>Sets the coefficient applied to the equations set source vector.
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_COEFFICIENT_SET(EQUATIONS_MAPPING,SOURCE_COEFFICIENT,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: SOURCE_COEFFICIENT!<The coefficient applied to the equations set source vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_SOURCE_COEFFICIENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE/=0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_COEFFICIENT=SOURCE_COEFFICIENT
          ELSE
            CALL FLAG_ERROR("The equations mapping source variable type has not been set.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_COEFFICIENT_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_SOURCE_COEFFICIENT_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_COEFFICIENT_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_COEFFICIENT_SET
  
 !
  !================================================================================================================================
  !

  !>Finalises the equations mapping source mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE(SOURCE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING !<A pointer to the SOURCE mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOURCE_MAPPING)) THEN
      IF(ALLOCATED(SOURCE_MAPPING%SOURCE_DOF_TO_EQUATIONS_ROW_MAP)) DEALLOCATE(SOURCE_MAPPING%SOURCE_DOF_TO_EQUATIONS_ROW_MAP)
      IF(ALLOCATED(SOURCE_MAPPING%EQUATIONS_ROW_TO_SOURCE_DOF_MAP)) DEALLOCATE(SOURCE_MAPPING%EQUATIONS_ROW_TO_SOURCE_DOF_MAP)
      DEALLOCATE(SOURCE_MAPPING)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping source mapping
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE(EQUATIONS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the source mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%SOURCE_MAPPING)) THEN
        CALL FLAG_ERROR("Equations mapping source mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%SOURCE_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations mapping source mapping.",ERR,ERROR,*999)
        EQUATIONS_MAPPING%SOURCE_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING        
        EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_VARIABLE_TYPE=0
        NULLIFY(EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_VARIABLE)
        NULLIFY(EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_VARIABLE_MAPPING)
        EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE(EQUATIONS_MAPPING%SOURCE_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a source field variable and the equations set source vector.
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,SOURCE_VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: SOURCE_VARIABLE_TYPE !<The variable type associated with the equations set source vector. If the problem does not have a source vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Equations mapping have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(SOURCE_VARIABLE_TYPE==0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=0
          ELSE
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
                  SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
                  IF(ASSOCIATED(SOURCE_FIELD)) THEN                    
                    !Check the source variable type is defined on the source field
                    IF(SOURCE_VARIABLE_TYPE>=1.AND.SOURCE_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                      IF(ASSOCIATED(SOURCE_FIELD%VARIABLE_TYPE_MAP(SOURCE_VARIABLE_TYPE)%PTR)) THEN
                        EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=SOURCE_VARIABLE_TYPE
                      ELSE
                        LOCAL_ERROR="The specified source variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " is not defined on the source field."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The specified source variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " is invalid. The number must either be zero or >= 1 and <= "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Source field is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations equations set is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations mapping equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Finalise an equations mapping equations matrix map.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE(VAR_TO_EQUATIONS_COLUMN_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_COLUMN_MAP_TYPE) :: VAR_TO_EQUATIONS_COLUMN_MAP !<The variable dof to equations column map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(VAR_TO_EQUATIONS_COLUMN_MAP%COLUMN_DOF)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_COLUMN_MAP%COLUMN_DOF)
    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE(VAR_TO_EQUATIONS_JACOBIAN_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_JACOBIAN_MAP_TYPE) :: VAR_TO_EQUATIONS_JACOBIAN_MAP !<The variable to equations Jacobian map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(VAR_TO_EQUATIONS_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP)
    IF(ALLOCATED(VAR_TO_EQUATIONS_JACOBIAN_MAP%DOF_TO_ROWS_MAP)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_JACOBIAN_MAP%DOF_TO_ROWS_MAP)
    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE(VAR_TO_EQUATIONS_JACOBIAN_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_JACOBIAN_MAP_TYPE) :: VAR_TO_EQUATIONS_JACOBIAN_MAP !<The variable to equations Jacobian map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE",ERR,ERROR,*999)
    
    VAR_TO_EQUATIONS_JACOBIAN_MAP%VARIABLE_TYPE=0
    NULLIFY(VAR_TO_EQUATIONS_JACOBIAN_MAP%VARIABLE)
    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations matrices map and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE(VAR_TO_EQUATIONS_MATRICES_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_MATRICES_MAP_TYPE) :: VAR_TO_EQUATIONS_MATRICES_MAP !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(VAR_TO_EQUATIONS_MATRICES_MAP%EQUATIONS_MATRIX_NUMBERS)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_MATRICES_MAP%EQUATIONS_MATRIX_NUMBERS)
    IF(ALLOCATED(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS)) THEN
      DO matrix_idx=1,SIZE(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS,1)
        CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS( &
          & matrix_idx),ERR,ERROR,*999)
      ENDDO !matrix_idx
      DEALLOCATE(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS)
    ENDIF
    IF(ALLOCATED(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_ROWS_MAP)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_ROWS_MAP)
    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an equations mapping equations matrix map.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE(VAR_TO_EQUATIONS_MATRICES_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_MATRICES_MAP_TYPE) :: VAR_TO_EQUATIONS_MATRICES_MAP !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE",ERR,ERROR,*999)

    VAR_TO_EQUATIONS_MATRICES_MAP%VARIABLE_INDEX=0
    VAR_TO_EQUATIONS_MATRICES_MAP%VARIABLE_TYPE=0
    NULLIFY(VAR_TO_EQUATIONS_MATRICES_MAP%VARIABLE)
    VAR_TO_EQUATIONS_MATRICES_MAP%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE EQUATIONS_MAPPING_ROUTINES
