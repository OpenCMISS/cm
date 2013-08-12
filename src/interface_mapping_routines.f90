!> \file
!> \author Chris Bradley
!> \brief This module contains all interface mapping routines.
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

!>This module contains all interface mapping routines.
MODULE INTERFACE_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  PUBLIC INTERFACE_MAPPING_CREATE_FINISH,INTERFACE_MAPPING_CREATE_START

  PUBLIC INTERFACE_MAPPING_DESTROY

  PUBLIC INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET

  PUBLIC INTERFACE_MAPPING_MATRICES_COEFFS_SET

  PUBLIC INTERFACE_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET,INTERFACE_MAPPING_MATRICES_ROW_MESH_INDICES_SET

  PUBLIC INTERFACE_MAPPING_MATRICES_NUMBER_SET

  PUBLIC INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET

  PUBLIC INTERFACE_MAPPING_RHS_COEFF_SET

  PUBLIC INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates interface mapping
  SUBROUTINE INTERFACE_MAPPING_CALCULATE(INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to calculate the mapping for
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,dof_idx,matrix_idx,mesh_idx,variable_idx,number_of_interface_matrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,LAGRANGE_VARIABLE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
      IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
        INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
            IF(ASSOCIATED(LAGRANGE)) THEN
              INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
              IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                !Set the Lagrange variable information
                LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                NULLIFY(LAGRANGE_VARIABLE)
                CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                  & ERR,ERROR,*999)
                INTERFACE_MAPPING%LAGRANGE_VARIABLE_TYPE=CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE
                INTERFACE_MAPPING%LAGRANGE_VARIABLE=>LAGRANGE_VARIABLE
                !Set the number of columns in the interface matrices
                INTERFACE_MAPPING%NUMBER_OF_COLUMNS=LAGRANGE_VARIABLE%NUMBER_OF_DOFS
                INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS=LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                INTERFACE_MAPPING%NUMBER_OF_GLOBAL_COLUMNS=LAGRANGE_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                !Set the column dofs mapping
                INTERFACE_MAPPING%COLUMN_DOFS_MAPPING=>LAGRANGE_VARIABLE%DOMAIN_MAPPING
                ALLOCATE(INTERFACE_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Lagrange dof to column map.",ERR,ERROR,*999)
                !1-1 mapping for now
                DO dof_idx=1,LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                  column_idx=LAGRANGE_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                  INTERFACE_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(dof_idx)=column_idx
                ENDDO
                !Set the number of interface matrices
                INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES=CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES
                ALLOCATE(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface matrix rows to variable maps.",ERR,ERROR,*999)
                !Loop over the interface matrices and calculate the row mappings
                !The pointers below have been checked for association above.
                SELECT CASE(INTERFACE_CONDITION%METHOD)
                CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                  number_of_interface_matrices=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                  !Number of interface matrices whose rows/columns are related to Dependent/Lagrange variables and not Lagrange/Lagrange variables (last interface matrix is Lagrange/Lagrange (Penalty matrix)
                  number_of_interface_matrices=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES-1 
                ENDSELECT
                DO matrix_idx=1,number_of_interface_matrices
                  !Initialise and setup the interface matrix
                  CALL INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE(INTERFACE_MAPPING,matrix_idx,ERR,ERROR,*999)
                  mesh_idx=CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)
                  NULLIFY(EQUATIONS_SET)
                  NULLIFY(FIELD_VARIABLE)
                  DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    IF(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)==mesh_idx) THEN
                      EQUATIONS_SET=>INTERFACE_DEPENDENT%EQUATIONS_SETS(variable_idx)%PTR
                      FIELD_VARIABLE=>INTERFACE_DEPENDENT%FIELD_VARIABLES(variable_idx)%PTR
                      EXIT
                    ENDIF
                  ENDDO !variable_idx
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%EQUATIONS_SET=>EQUATIONS_SET
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MESH_INDEX=mesh_idx
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT=INTERFACE_MAPPING% &
                        & CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=INTERFACE_MAPPING% &
                        & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                       !Set the number of rows
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE%NUMBER_OF_DOFS
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                        & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                        & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                      !Set the row mapping
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                        & FIELD_VARIABLE%DOMAIN_MAPPING
                      ALLOCATE(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                        & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                      !1-1 mapping for now
                      DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                        INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
                      ENDDO !dof_idx
                    ELSE
                      LOCAL_ERROR="Dependent variable for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                        & " could not be found."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Equations set for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                      & " could not be found."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx

                !The pointers below have been checked for association above.
                SELECT CASE(INTERFACE_CONDITION%METHOD)
                CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                  !Sets up the Lagrange-(Penalty) interface matrix mapping and calculate the row mappings
                  matrix_idx = INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES !last of the interface matrices
                  !Initialise and setup the interface matrix
                  CALL INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE(INTERFACE_MAPPING,matrix_idx,ERR,ERROR,*999)
                  mesh_idx=CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)
                  NULLIFY(LAGRANGE_VARIABLE)
                  CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                    & ERR,ERROR,*999)
                  NULLIFY(INTERFACE_EQUATIONS)
                  NULLIFY(FIELD_VARIABLE)
                  FIELD_VARIABLE=>LAGRANGE_VARIABLE
                  INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
                  IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MESH_INDEX=mesh_idx
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT=INTERFACE_MAPPING% &
                        & CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=INTERFACE_MAPPING% &
                        & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                        !Set the number of rows
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE%NUMBER_OF_DOFS
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                        & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                        & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                      !Set the row mapping
                      INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                        & FIELD_VARIABLE%DOMAIN_MAPPING
                      ALLOCATE(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                        & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                      !1-1 mapping for now
                      DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                        INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
                      ENDDO !dof_idx
                    ELSE
                      LOCAL_ERROR="Lagrange variable for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                        & " could not be found."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Interface Equations for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                      & " could not be found."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDSELECT

                !Calculate RHS mappings
                IF(CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE/=0) THEN
                  CALL INTERFACE_MAPPING_RHS_MAPPING_INITIALISE(INTERFACE_MAPPING,ERR,ERROR,*999)
                  RHS_MAPPING=>INTERFACE_MAPPING%RHS_MAPPING
                  IF(ASSOCIATED(RHS_MAPPING)) THEN
                    RHS_MAPPING%RHS_VARIABLE_TYPE=CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE
                    LAGRANGE_VARIABLE=>LAGRANGE_FIELD%VARIABLE_TYPE_MAP(CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE)%PTR
                    RHS_MAPPING%RHS_VARIABLE=>LAGRANGE_VARIABLE
                    RHS_MAPPING%RHS_VARIABLE_MAPPING=>LAGRANGE_VARIABLE%DOMAIN_MAPPING
                    RHS_MAPPING%RHS_COEFFICIENT=CREATE_VALUES_CACHE%RHS_COEFFICIENT
                    !Allocate and set up the row mappings
                    ALLOCATE(RHS_MAPPING%RHS_DOF_TO_INTERFACE_ROW_MAP(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate rhs dof to interface row map.",ERR,ERROR,*999)
                    ALLOCATE(RHS_MAPPING%INTERFACE_ROW_TO_RHS_DOF_MAP(INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface row to dof map.",ERR,ERROR,*999)
                    DO dof_idx=1,LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                      !1-1 mapping for now
                      column_idx=dof_idx
                      RHS_MAPPING%RHS_DOF_TO_INTERFACE_ROW_MAP(dof_idx)=column_idx
                    ENDDO !dof_idx
                    DO column_idx=1,INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS
                      !1-1 mapping for now
                      dof_idx=column_idx
                      RHS_MAPPING%INTERFACE_ROW_TO_RHS_DOF_MAP(column_idx)=dof_idx
                    ENDDO !column_idx
                  ELSE
                    CALL FLAG_ERROR("RHS mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_CALCULATE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_CALCULATE")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of interface mapping
  SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH(INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to finish the creation of.
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        !Calculate the equations mapping and clean up
        CALL INTERFACE_MAPPING_CALCULATE(INTERFACE_MAPPING,ERR,ERROR,*999)
        CALL INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE(INTERFACE_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
        INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED=.TRUE.
      ENDIF
   ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the process of creating the interface mapping for interface equations.
  SUBROUTINE INTERFACE_MAPPING_CREATE_START(INTERFACE_EQUATIONS,INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to create the mapping for.
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<On exit, a pointer to the created interface mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          CALL FLAG_ERROR("Interface mapping is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(INTERFACE_MAPPING)
          CALL INTERFACE_MAPPING_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
          INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_CREATE_START")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises an interface mapping create values cache and deallocates all memory
  SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%HAS_TRANSPOSE)) DEALLOCATE(CREATE_VALUES_CACHE%HAS_TRANSPOSE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES))  &
        & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COL_FIELD_VARIABLE_INDICES)) &
        & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COL_FIELD_VARIABLE_INDICES)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface mapping create values cache 
  SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to create the create values cache for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_idx,variable_type_idx,variable_type_idx2
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ASSOCIATED(INTERFACE_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Interface mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
          IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
            !Allocate and initialise the create values cache
            ALLOCATE(INTERFACE_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface mapping create values cache.",ERR,ERROR,*999)
            INTERFACE_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES=0
            INTERFACE_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=0
            INTERFACE_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=0
            INTERFACE_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=0.0_DP
            !Set the default interface mapping in the create values cache
            !First calculate how many interface matrices we have and set the variable types
            SELECT CASE(INTERFACE_CONDITION%METHOD)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
              LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
              IF(ASSOCIATED(LAGRANGE)) THEN
                LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                  INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
                  IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                    !The pointers below have been checked for association above.
                    SELECT CASE(INTERFACE_CONDITION%METHOD)
                    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      !Default the number of interface matrices to the number of added dependent variables
                      INTERFACE_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES= &
                      INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                      !Default the number of interface matrices to the number of added dependent variables plus a single Lagrange variable
                      INTERFACE_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES= &
                      INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1
                    END SELECT
                    !Default the Lagrange variable to the first Lagrange variable
                    INTERFACE_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=0
                    DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                      IF(ASSOCIATED(LAGRANGE_FIELD%VARIABLE_TYPE_MAP(variable_type_idx)%PTR)) THEN
                        INTERFACE_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=variable_type_idx
                        EXIT
                      ENDIF
                    ENDDO !variable_type_idx
                    IF(INTERFACE_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE==0) &
                      & CALL FLAG_ERROR("Could not find a Lagrange variable type in the Lagrange field.",ERR,ERROR,*999)
                    !Default the RHS Lagrange variable to the second Lagrange variable
                    DO variable_type_idx2=variable_type_idx+1,FIELD_NUMBER_OF_VARIABLE_TYPES
                      IF(ASSOCIATED(LAGRANGE_FIELD%VARIABLE_TYPE_MAP(variable_type_idx2)%PTR)) THEN
                        INTERFACE_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=variable_type_idx2
                        EXIT
                      ENDIF
                    ENDDO !variable_type_idx2
                    IF(INTERFACE_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE==0) &
                      & CALL FLAG_ERROR("Could not find a RHS Lagrange variable type in the Lagrange field.",ERR,ERROR,*999)
                    ALLOCATE(INTERFACE_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(INTERFACE_MAPPING% &
                      & CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix coefficients.",ERR,ERROR,*999)
                    !Default the interface matrices coefficients to add.
                    INTERFACE_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=1.0_DP
                    INTERFACE_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=1.0_DP
                    ALLOCATE(INTERFACE_MAPPING%CREATE_VALUES_CACHE%HAS_TRANSPOSE(INTERFACE_MAPPING% &
                      & CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache has transpose.",ERR,ERROR,*999)
                    !Default the interface matrices to all have a transpose
                    INTERFACE_MAPPING%CREATE_VALUES_CACHE%HAS_TRANSPOSE=.TRUE.
                    ALLOCATE(INTERFACE_MAPPING%CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(INTERFACE_MAPPING% &
                      & CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix row field variable indexes.", &
                      & ERR,ERROR,*999)
                    !Default the interface matrices to be in mesh index order.
                    DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                      INTERFACE_MAPPING%CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(variable_idx)=variable_idx
                    ENDDO !variable_idx
                    !The pointers below have been checked for association above.
                    SELECT CASE(INTERFACE_CONDITION%METHOD)
                    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                      !Default the interface matrix (Penalty) to have no transpose
                      INTERFACE_MAPPING%CREATE_VALUES_CACHE%HAS_TRANSPOSE(INTERFACE_MAPPING% &
                        & CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)=.FALSE.
                      !Default the interface matrices to be in mesh index order (and set Penalty matrix (last interface matrix)to be the first Lagrange variable).
                      INTERFACE_MAPPING%CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(INTERFACE_DEPENDENT% &
                        & NUMBER_OF_DEPENDENT_VARIABLES+1)=1
                    END SELECT
                  ELSE
                    CALL FLAG_ERROR("Interface condition depdendent is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Interface condition Lagrange field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface equations method of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE(INTERFACE_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Destroys an interface mapping.
  SUBROUTINE INTERFACE_MAPPING_DESTROY(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer the interface mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      CALL INTERFACE_MAPPING_FINALISE(INTERFACE_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("INTERFACE_MAPPING_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the interface mapping and deallocates all memory.
  SUBROUTINE INTERFACE_MAPPING_FINALISE(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("INTERFACE_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ALLOCATED(INTERFACE_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)) DEALLOCATE(INTERFACE_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)
      IF(ALLOCATED(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS)) THEN
        DO matrix_idx=1,SIZE(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS,1)
          CALL INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_FINALISE(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx), &
            & ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS)
      ENDIF
      CALL INTERFACE_MAPPING_RHS_MAPPING_FINALISE(INTERFACE_MAPPING%RHS_MAPPING,ERR,ERROR,*999)
      CALL INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE(INTERFACE_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_MAPPING)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_FINALISE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mapping and deallocates all memory.
  SUBROUTINE INTERFACE_MAPPING_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to initialise the interface mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("INTERFACE_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MAPPING)) THEN
        CALL FLAG_ERROR("Interface mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(INTERFACE_EQUATIONS%INTERFACE_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface equations interface mapping.",ERR,ERROR,*999)
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED=.FALSE.
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%LAGRANGE_VARIABLE_TYPE=0
        NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MAPPING%LAGRANGE_VARIABLE)
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%NUMBER_OF_COLUMNS=0
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS=0
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%NUMBER_OF_GLOBAL_COLUMNS=0
        NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MAPPING%COLUMN_DOFS_MAPPING)
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES=0
        NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MAPPING%RHS_MAPPING)
        NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MAPPING%CREATE_VALUES_CACHE)
        CALL INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE(INTERFACE_EQUATIONS%INTERFACE_MAPPING,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_INITIALISE")
    RETURN
999 CALL INTERFACE_MAPPING_FINALISE(INTERFACE_EQUATIONS%INTERFACE_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_INITIALISE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the Lagrange variable type for the interface mapping. 
  SUBROUTINE INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET(INTERFACE_MAPPING,LAGRANGE_VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: LAGRANGE_VARIABLE_TYPE !<The Lagrange variable type to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: LAGRANGE_VARIABLE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
                IF(ASSOCIATED(LAGRANGE)) THEN
                  IF(LAGRANGE%LAGRANGE_FINISHED) THEN
                    LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                    NULLIFY(LAGRANGE_VARIABLE)
                    CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE,ERR,ERROR,*999)
                    CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=LAGRANGE_VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Interface condition Lagrange field has not been finished.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
                ENDIF
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Finalises an interface matrix to variable map and deallocates all memory.
  SUBROUTINE INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_FINALISE(INTERFACE_MATRIX_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRIX_TO_VAR_MAP_TYPE) :: INTERFACE_MATRIX_TO_VAR_MAP !<The interface matrix to var map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    CALL ENTERS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(INTERFACE_MATRIX_TO_VAR_MAP%VARIABLE_DOF_TO_ROW_MAP)) &
      & DEALLOCATE(INTERFACE_MATRIX_TO_VAR_MAP%VARIABLE_DOF_TO_ROW_MAP)
    
    CALL EXITS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface matrix to variable map.
  SUBROUTINE INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE(INTERFACE_MAPPING,matrix_idx,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to initialise the matrix to variable map for a given matrix index.
    INTEGER(INTG), INTENT(IN) :: matrix_idx !<The matrix index to intialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
  
    CALL ENTERS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ALLOCATED(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS)) THEN
        IF(matrix_idx>0.AND.matrix_idx<=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES) THEN
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_NUMBER=matrix_idx
          NULLIFY(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%INTERFACE_MATRIX)
          NULLIFY(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%EQUATIONS_SET)
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=0
          NULLIFY(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE)
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MESH_INDEX=0
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT=0.0_DP
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=.FALSE.
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS=0
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS=0
          INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS=0
          NULLIFY(INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING)          
        ELSE
          LOCAL_ERROR="The specified matrix index of "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is invalid. The index must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface mapping matrix rows to var maps is not allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the interface matrices. 
  SUBROUTINE INTERFACE_MAPPING_MATRICES_COEFFS_SET(INTERFACE_MAPPING,MATRIX_COEFFICIENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    REAL(DP), INTENT(IN) :: MATRIX_COEFFICIENTS(:) !<The interface matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_MATRICES_COEFFS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN          
           INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                !Check that the number of supplied coefficients matches the number of interface matrices
                IF(SIZE(MATRIX_COEFFICIENTS,1)==CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES) THEN
                  CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                    & MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)
                ELSE
                  LOCAL_ERROR="Invalid size of matrix coefficeints. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_COEFFICIENTS,1),"*",ERR,ERROR))// &
                    & ") must match the number of interface matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_MATRICES_COEFFS_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_MATRICES_COEFFS_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_MATRICES_COEFFS_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_COEFFS_SET

  !
  !================================================================================================================================
  !

  !>Sets the column mesh indices for the interface matrices. 
  SUBROUTINE INTERFACE_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET(INTERFACE_MAPPING,COLUMN_MESH_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: COLUMN_MESH_INDICES(:) !<The interface matrix column mesh indices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                CALL FLAG_ERROR("Can not set the column mesh indices when using the Lagrange multipliers "// &
                  "interface condition method.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET

  !
  !================================================================================================================================
  !

  !>Sets the row mesh indices for the interface matrices. 
  SUBROUTINE INTERFACE_MAPPING_MATRICES_ROW_MESH_INDICES_SET(INTERFACE_MAPPING,ROW_MESH_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: ROW_MESH_INDICES(:) !<The interface matrix mesh indices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_idx2,mesh_idx3
    LOGICAL :: FOUND
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_MATRICES_ROW_MESH_INDICES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                !Check the size of the mesh indicies array
                IF(SIZE(ROW_MESH_INDICES,1)==CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES) THEN
                  !Check that mesh indices are valid.
                  INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
                  IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                    DO mesh_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES
                      FOUND=.FALSE.
                      DO mesh_idx2=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                        IF(ROW_MESH_INDICES(mesh_idx)==INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(mesh_idx2)) THEN
                          FOUND=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO !mesh_idx2
                      IF(FOUND) THEN
                        !Check that the mesh index has not been repeated.
                        DO mesh_idx3=mesh_idx+1,CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES
                          IF(ROW_MESH_INDICES(mesh_idx)==ROW_MESH_INDICES(mesh_idx3)) THEN
                            LOCAL_ERROR="The supplied mesh index of "// &
                              & TRIM(NUMBER_TO_VSTRING(ROW_MESH_INDICES(mesh_idx),"*",ERR,ERROR))// &
                              & " at position "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                              & " has been repeated at position "//TRIM(NUMBER_TO_VSTRING(mesh_idx3,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                        ENDDO !mesh_idx3
                        !Set the mesh indices
                        CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                          & ROW_MESH_INDICES(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)
                      ELSE
                        LOCAL_ERROR="The supplied mesh index of "// &
                          & TRIM(NUMBER_TO_VSTRING(ROW_MESH_INDICES(mesh_idx),"*",ERR,ERROR))// &
                          & " at position "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                          & " has not been added as a dependent variable to the interface condition."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !mesh_idx
                  ELSE
                    CALL FLAG_ERROR("Interface condition dependent is not assocaited.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Invalid size of mesh indices. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_MESH_INDICES,1),"*",ERR,ERROR))// &
                    & ") must match the number of interface matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_MATRICES_ROW_MESH_INDICES_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_MATRICES_ROW_MESH_INDICES_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_MATRICES_ROW_MESH_INDICES_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_ROW_MESH_INDICES_SET

  !
  !================================================================================================================================
  !

  !>Sets the number of interface matrices for an interface mapping.
  SUBROUTINE INTERFACE_MAPPING_MATRICES_NUMBER_SET(INTERFACE_MAPPING,NUMBER_OF_INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to set the number of linear matrices for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_INTERFACE_MATRICES !<The number of interface matrices to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,matrix_idx2,variable_idx,number_of_dependent_variables
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(:)
    LOGICAL :: FOUND
    LOGICAL, ALLOCATABLE :: OLD_MATRIX_TRANSPOSE(:)
    REAL(DP), ALLOCATABLE :: OLD_MATRIX_COEFFICIENTS(:)
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                !Check the number of interface matrices
                IF(NUMBER_OF_INTERFACE_MATRICES>0) THEN
                  INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
                  IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                    SELECT CASE(INTERFACE_CONDITION%METHOD)
                    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      number_of_dependent_variables=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                      number_of_dependent_variables=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1
                    END SELECT
                    IF(NUMBER_OF_INTERFACE_MATRICES<=number_of_dependent_variables) THEN
                      !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
                      IF(NUMBER_OF_INTERFACE_MATRICES/=CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES) THEN
                        ALLOCATE(OLD_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix coefficients.",ERR,ERROR,*999)
                        ALLOCATE(OLD_MATRIX_TRANSPOSE(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix transpose.",ERR,ERROR,*999)
                        ALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix row field indexes.",ERR,ERROR,*999)
                        OLD_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)
                        OLD_MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                          & CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)
                        OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                          & CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE% &
                          & NUMBER_OF_INTERFACE_MATRICES)
                        IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)) &
                          & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
                        IF(ALLOCATED(CREATE_VALUES_CACHE%HAS_TRANSPOSE)) &
                          & DEALLOCATE(CREATE_VALUES_CACHE%HAS_TRANSPOSE)
                        IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES)) &
                          & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES)
                        ALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix coefficients.",ERR,ERROR,*999)
                        ALLOCATE(CREATE_VALUES_CACHE%HAS_TRANSPOSE(NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix tranpose.",ERR,ERROR,*999)
                        ALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix row field variable indexes.",ERR,ERROR,*999)
                        IF(NUMBER_OF_INTERFACE_MATRICES>CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES) THEN
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                            & OLD_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES+1: &
                            & NUMBER_OF_INTERFACE_MATRICES)=1.0_DP
                          CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                            & OLD_MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)
                          CREATE_VALUES_CACHE%HAS_TRANSPOSE(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES+1: &
                            & NUMBER_OF_INTERFACE_MATRICES)=.TRUE.
                          CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE% &
                            & NUMBER_OF_INTERFACE_MATRICES)=OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE% &
                            & NUMBER_OF_INTERFACE_MATRICES)
                          !Loop through in mesh index order and set the default matrix to variable map to be in mesh index order
                          DO matrix_idx=CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES+1,NUMBER_OF_INTERFACE_MATRICES
                            CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)=0
                            DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                              FOUND=.FALSE.
                              DO matrix_idx2=1,CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES
                                IF(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)==CREATE_VALUES_CACHE% &
                                  MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx2)) THEN
                                  FOUND=.TRUE.
                                  EXIT
                                ENDIF
                              ENDDO !matrix_idx2
                              IF(.NOT.FOUND) THEN
                                CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)=INTERFACE_DEPENDENT% &
                                  & VARIABLE_MESH_INDICES(variable_idx)
                              ENDIF
                            ENDDO !variable_idx2
                            IF(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)==0) THEN
                              LOCAL_ERROR="Could not map an interface mesh index for interface matrix "// &
                                & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//"."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !matrix_idx
                        ELSE
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:NUMBER_OF_INTERFACE_MATRICES)= &
                            & OLD_MATRIX_COEFFICIENTS(1:NUMBER_OF_INTERFACE_MATRICES)
                          CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:NUMBER_OF_INTERFACE_MATRICES)= &
                            & OLD_MATRIX_TRANSPOSE(1:NUMBER_OF_INTERFACE_MATRICES)
                          CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(1:NUMBER_OF_INTERFACE_MATRICES)= &
                            & OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:NUMBER_OF_INTERFACE_MATRICES)
                        ENDIF
                        IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
                        IF(ALLOCATED(OLD_MATRIX_TRANSPOSE)) DEALLOCATE(OLD_MATRIX_TRANSPOSE)
                        IF(ALLOCATED(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)) DEALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The specified number of interface matrices of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))// &
                        & " is invalid. The number must be <= the number of added dependent variables of "// &
                        & TRIM(NUMBER_TO_VSTRING(number_of_dependent_variables,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified number of interface matrices of "// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))// &
                    & " is invalid. The number must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
    IF(ALLOCATED(OLD_MATRIX_TRANSPOSE)) DEALLOCATE(OLD_MATRIX_TRANSPOSE)
    IF(ALLOCATED(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)) DEALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)
    CALL ERRORS("INTERFACE_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_MATRICES_NUMBER_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Sets the transpose flag for the interface matrices. 
  SUBROUTINE INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET(INTERFACE_MAPPING,MATRIX_TRANSPOSE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    LOGICAL, INTENT(IN) :: MATRIX_TRANSPOSE(:) !<MATRIX_TRANSPOSE(matrix_idx). The interface matrix transpose flag for the matrix_idx'th interface matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
           INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              SELECT CASE(INTERFACE_CONDITION%METHOD)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                !Check that the number of supplied coefficients matches the number of interface matrices
                IF(SIZE(MATRIX_TRANSPOSE,1)==CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES) THEN
                  CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)= &
                    MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES)
                ELSE
                  LOCAL_ERROR="Invalid size of matrix tranpose. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_TRANSPOSE,1),"*",ERR,ERROR))// &
                    & ") must match the number of interface matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interface condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the interface RHS vector.
  SUBROUTINE INTERFACE_MAPPING_RHS_COEFF_SET(INTERFACE_MAPPING,RHS_COEFFICIENT,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to set the RHS coefficent for
    REAL(DP), INTENT(IN) :: RHS_COEFFICIENT !<The coefficient applied to the interface RHS vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_RHS_COEFF_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(INTERFACE_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE/=0) THEN
            INTERFACE_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=RHS_COEFFICIENT
          ELSE
            CALL FLAG_ERROR("The interface mapping RHS Lagrange variable type has not been set.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_RHS_COEFF_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_RHS_COEFF_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_RHS_COEFF_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_COEFF_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the interface mapping RHS mapping and deallocates all memory
  SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_FINALISE(RHS_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("INTERFACE_MAPPING_RHS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(RHS_MAPPING)) THEN
      IF(ALLOCATED(RHS_MAPPING%RHS_DOF_TO_INTERFACE_ROW_MAP)) DEALLOCATE(RHS_MAPPING%RHS_DOF_TO_INTERFACE_ROW_MAP)
      IF(ALLOCATED(RHS_MAPPING%INTERFACE_ROW_TO_RHS_DOF_MAP)) DEALLOCATE(RHS_MAPPING%INTERFACE_ROW_TO_RHS_DOF_MAP)
      DEALLOCATE(RHS_MAPPING)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_RHS_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_RHS_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_RHS_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mapping RHS mapping
  SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_INITIALISE(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("INTERFACE_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ASSOCIATED(INTERFACE_MAPPING%RHS_MAPPING)) THEN
        CALL FLAG_ERROR("Interface mapping RHS mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(INTERFACE_MAPPING%RHS_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface mapping RHS mapping.",ERR,ERROR,*999)
        INTERFACE_MAPPING%RHS_MAPPING%INTERFACE_MAPPING=>INTERFACE_MAPPING        
        INTERFACE_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE=0
        NULLIFY(INTERFACE_MAPPING%RHS_MAPPING%RHS_VARIABLE)
        NULLIFY(INTERFACE_MAPPING%RHS_MAPPING%RHS_VARIABLE_MAPPING)
        INTERFACE_MAPPING%RHS_MAPPING%RHS_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN
999 CALL INTERFACE_MAPPING_RHS_MAPPING_FINALISE(INTERFACE_MAPPING%RHS_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a Lagrange field variable and the interface rhs vector.
  SUBROUTINE INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET(INTERFACE_MAPPING,RHS_VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to set the RHS variable type for.
    INTEGER(INTG), INTENT(IN) :: RHS_VARIABLE_TYPE !<The variable type associated with the interface rhs vector. If the interface condition does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: INTERFACE_LAGRANGE
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Interface mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>INTERFACE_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          IF(RHS_VARIABLE_TYPE==0) THEN
            CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=0
          ELSE
            INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
            IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
              INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
              IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                SELECT CASE(INTERFACE_CONDITION%METHOD)
                CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                  INTERFACE_LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
                  IF(ASSOCIATED(INTERFACE_LAGRANGE)) THEN
                    LAGRANGE_FIELD=>INTERFACE_LAGRANGE%LAGRANGE_FIELD
                    IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                      !Check the RHS variable type is not being by the interface matrices
                      IF(CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE==RHS_VARIABLE_TYPE) THEN
                        LOCAL_ERROR="The specified RHS variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " is the same as the Lagrange variable type for the interface matrices."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Check the RHS variable number is defined on the Lagrange field
                      IF(RHS_VARIABLE_TYPE>=1.AND.RHS_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                        IF(ASSOCIATED(LAGRANGE_FIELD%VARIABLE_TYPE_MAP(RHS_VARIABLE_TYPE)%PTR)) THEN
                          CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=RHS_VARIABLE_TYPE
                        ELSE
                          LOCAL_ERROR="The specified RHS variable type of "// &
                            & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " is not defined on the Lagrange field."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="The specified RHS variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " is invalid. The number must either be zero or >= 1 and <= "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Lagrange field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface Lagrange is not associated.",ERR,ERROR,*999)
                  ENDIF
                CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The interface condition method of "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

END MODULE INTERFACE_MAPPING_ROUTINES
