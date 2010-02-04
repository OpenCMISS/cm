!> \file
!> $Id$
!> \author Ting Yu
!> \brief This module handles all analytic analysis routines.
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

!>This module handles all analytic analysis routines.
MODULE ANALYTIC_ANALYSIS_ROUTINES

  USE BASIS_ROUTINES
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE MPI
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup ANALYTIC_ANALYSIS_ROUTINES_ErrorTypes ANALYTIC_ANALYSIS_ROUTINES::ErrorTypes
  !> \brief errors definition type parameters
  !> \see ANALYTIC_ANALYSIS_ROUTINES,OPENCMISS_ErrorTypes
  !>@{
  INTEGER(INTG), PARAMETER :: ABSOLUTE_ERROR_TYPE=1 !<The absolute type \see ANALYTIC_ANALYSIS_ROUTINES_ErrorTypes,ANALYTIC_ANALYSIS_ROUTINES
  INTEGER(INTG), PARAMETER :: PERCENTAGE_ERROR_TYPE=2 !<The percentage type \see ANALYTIC_ANALYSIS_ROUTINES_ErrorTypes,ANALYTIC_ANALYSIS_ROUTINES
  INTEGER(INTG), PARAMETER :: RELATIVE_ERROR_TYPE=3 !<The relative type \see ANALYTIC_ANALYSIS_ROUTINES_ErrorTypes,ANALYTIC_ANALYSIS_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC ANALYTIC_ANALYSIS_OUTPUT
  
  PUBLIC ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE,ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE, &
    & ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE

  PUBLIC ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT,ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT, &
    & ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT

  PUBLIC ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_CONSTANT,ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_CONSTANT, &
    & ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_CONSTANT
 
  PUBLIC ANALYTIC_ANALYSIS_RMS_ERROR_GET_NODE,ANALYTIC_ANALYSIS_RMS_ERROR_GET_ELEMENT
  
  PUBLIC ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET,ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET, &
    & ANALYTIC_ANALYSIS_INTEGRAL_PERCENTAGE_ERROR_GET,ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET, &
    & ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET
    
CONTAINS  

  !
  !================================================================================================================================
  !  

  !>Output the analytic error analysis for a dependent field compared to the analytic values parameter set. \see OPENCMISS::CMISSAnalyticAnalytisOutput
  SUBROUTINE ANALYTIC_ANALYSIS_OUTPUT(FIELD,FILENAME,ERR,ERROR,*)
  
    !Argument variables 
    TYPE(FIELD_TYPE), INTENT(IN), POINTER :: FIELD !<A pointer to the dependent field to calculate the analytic error analysis for
    CHARACTER(LEN=*) :: FILENAME !<If not empty, the filename to output the analytic analysis to. If empty, the analysis will be output to the standard output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,element_idx,GHOST_NUMBER(8),local_ny,MESH_COMPONENT,MPI_IERROR,node_idx, &
      & NUMBER(8),OUTPUT_ID,var_idx,variable_type
    REAL(DP) :: GHOST_RMS_ERROR_PER(8),GHOST_RMS_ERROR_ABS(8),GHOST_RMS_ERROR_REL(8),RMS_ERROR_PER(8),RMS_ERROR_ABS(8), &
      & RMS_ERROR_REL(8),VALUES(5)
    REAL(DP), POINTER :: ANALYTIC_VALUES(:),NUMERICAL_VALUES(:)
    REAL(DP), ALLOCATABLE :: INTEGRAL_ERRORS(:,:),GHOST_INTEGRAL_ERRORS(:,:)
    CHARACTER(LEN=40) :: FIRST_FORMAT
    CHARACTER(LEN=MAXSTRLEN) :: FILE_NAME
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: ELEMENTS_DECOMPOSITION
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: DECOMPOSITION_TOPOLOGY
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES_DOMAIN
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_STRING
    
    CALL ENTERS("ANALYTIC_ANALYSIS_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(FIELD%DEPENDENT_TYPE==FIELD_DEPENDENT_TYPE) THEN
          IF(LEN_TRIM(FILENAME)>=1) THEN
!!TODO \todo have more general ascii file mechanism
            IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>1) THEN
              WRITE(FILE_NAME,'(A,".opanal.",I0)') FILENAME(1:LEN_TRIM(FILENAME)),COMPUTATIONAL_ENVIRONMENT% &
                & MY_COMPUTATIONAL_NODE_NUMBER
            ELSE
              FILE_NAME=FILENAME(1:LEN_TRIM(FILENAME))//".opanal"
            ENDIF
            OUTPUT_ID=IO1_FILE_UNIT
            OPEN(UNIT=OUTPUT_ID,FILE=FILE_NAME(1:LEN_TRIM(FILE_NAME)),STATUS="REPLACE",FORM="FORMATTED",IOSTAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Error opening analysis output file.",ERR,ERROR,*999)            
          ELSE
            OUTPUT_ID=GENERAL_OUTPUT_TYPE
          ENDIF
          DECOMPOSITION=>FIELD%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            DECOMPOSITION_TOPOLOGY=>DECOMPOSITION%TOPOLOGY
            IF(ASSOCIATED(DECOMPOSITION_TOPOLOGY)) THEN
              CALL WRITE_STRING(OUTPUT_ID,"Analytic error analysis:",ERR,ERROR,*999)
              CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
              LOCAL_STRING="Field "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" : "//FIELD%LABEL
              IF(ERR/=0) GOTO 999
              CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
              !Loop over the variables
              DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
                variable_type=FIELD%VARIABLES(var_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
                  LOCAL_STRING="Variable "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))//" : "// &
                    & FIELD_VARIABLE%VARIABLE_LABEL
                  IF(ERR/=0) GOTO 999
                  CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                  CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
                  !Get the dependent and analytic parameter sets
                  CALL FIELD_PARAMETER_SET_DATA_GET(FIELD,variable_type,FIELD_VALUES_SET_TYPE,NUMERICAL_VALUES,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_DATA_GET(FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ANALYTIC_VALUES, &
                    & ERR,ERROR,*999)
                  !Loop over the components
                  DO component_idx=1,FIELD%VARIABLES(var_idx)%NUMBER_OF_COMPONENTS
                    MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER
                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                    IF(ASSOCIATED(DOMAIN)) THEN
                      DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                      IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                        LOCAL_STRING="Component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" : "// &
                          & FIELD_VARIABLE%COMPONENTS(component_idx)%COMPONENT_LABEL
                        IF(ERR/=0) GOTO 999
                        CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                        CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
                        SELECT CASE(FIELD%VARIABLES(var_idx)%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                        CASE(FIELD_CONSTANT_INTERPOLATION)
                          CALL WRITE_STRING(OUTPUT_ID,"Constant errors:",ERR,ERROR,*999)
                          LOCAL_STRING="                       Numerical      Analytic       % error  Absolute err  Relative err"
                          CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                          VALUES(1)=NUMERICAL_VALUES(local_ny)
                          VALUES(2)=ANALYTIC_VALUES(local_ny)
                          VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                          VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                          VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                          CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,"(20X,5(2X,E12.5))","(20X,5(2X,E12.5))", &
                            & ERR,ERROR,*999)
                        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                          ELEMENTS_DOMAIN=>DOMAIN_TOPOLOGY%ELEMENTS
                          IF(ASSOCIATED(ELEMENTS_DOMAIN)) THEN
                            DECOMPOSITION=>DOMAIN%DECOMPOSITION
                            IF(ASSOCIATED(DECOMPOSITION)) THEN
                              DECOMPOSITION_TOPOLOGY=>DECOMPOSITION%TOPOLOGY
                              IF(ASSOCIATED(DECOMPOSITION_TOPOLOGY)) THEN
                                ELEMENTS_DECOMPOSITION=>DECOMPOSITION_TOPOLOGY%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_DECOMPOSITION)) THEN
                                  NUMBER=0
                                  RMS_ERROR_PER=0.0_DP
                                  RMS_ERROR_ABS=0.0_DP
                                  RMS_ERROR_REL=0.0_DP
                                  GHOST_NUMBER=0
                                  GHOST_RMS_ERROR_PER=0.0_DP
                                  GHOST_RMS_ERROR_ABS=0.0_DP
                                  GHOST_RMS_ERROR_REL=0.0_DP
                                  CALL WRITE_STRING(OUTPUT_ID,"Element errors:",ERR,ERROR,*999)
                                  LOCAL_STRING= &
                                    & "  Element#             Numerical      Analytic       % error  Absolute err  Relative err"
                                  CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                  DO element_idx=1,ELEMENTS_DOMAIN%NUMBER_OF_ELEMENTS
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(element_idx)
                                    VALUES(1)=NUMERICAL_VALUES(local_ny)
                                    VALUES(2)=ANALYTIC_VALUES(local_ny)
                                    VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                                    VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                                    VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                                    !Accumlate the RMS errors
                                    NUMBER(1)=NUMBER(1)+1
                                    RMS_ERROR_PER(1)=RMS_ERROR_PER(1)+VALUES(3)*VALUES(3)
                                    RMS_ERROR_ABS(1)=RMS_ERROR_ABS(1)+VALUES(4)*VALUES(4)
                                    RMS_ERROR_REL(1)=RMS_ERROR_REL(1)+VALUES(5)*VALUES(5)
                                    WRITE(FIRST_FORMAT,"(A,I10,A)") "('",ELEMENTS_DECOMPOSITION%ELEMENTS(element_idx)%USER_NUMBER, &
                                      & "',20X,3(2X,E12.5))"
                                    CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))", &
                                      & ERR,ERROR,*999)
                                  ENDDO !element_idx
                                  DO element_idx=ELEMENTS_DOMAIN%NUMBER_OF_ELEMENTS+1,ELEMENTS_DOMAIN%TOTAL_NUMBER_OF_ELEMENTS
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(element_idx)
                                    VALUES(1)=NUMERICAL_VALUES(local_ny)
                                    VALUES(2)=ANALYTIC_VALUES(local_ny)
                                    VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                                    VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                                    VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                                    !Accumlate the RMS errors
                                    GHOST_NUMBER(1)=GHOST_NUMBER(1)+1
                                    GHOST_RMS_ERROR_PER(1)=GHOST_RMS_ERROR_PER(1)+VALUES(3)*VALUES(3)
                                    GHOST_RMS_ERROR_ABS(1)=GHOST_RMS_ERROR_ABS(1)+VALUES(4)*VALUES(4)
                                    GHOST_RMS_ERROR_REL(1)=GHOST_RMS_ERROR_REL(1)+VALUES(5)*VALUES(5)
                                    WRITE(FIRST_FORMAT,"(A,I10,A)") "('",ELEMENTS_DECOMPOSITION%ELEMENTS(element_idx)%USER_NUMBER, &
                                      & "',20X,3(2X,E12.5))"
                                    CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))", &
                                      & ERR,ERROR,*999)
                                  ENDDO !node_idx
                                  !Output RMS errors                  
                                  CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
                                  IF(NUMBER(1)>0) THEN
                                    IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>1) THEN
                                      !Local elements only
                                      CALL WRITE_STRING(OUTPUT_ID,"Local RMS errors:",ERR,ERROR,*999)
                                      LOCAL_STRING= &
                                        & "                                                     % error  Absolute err  Relative err"
                                      CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                      VALUES(1)=SQRT(RMS_ERROR_PER(deriv_idx)/NUMBER(deriv_idx))
                                      VALUES(2)=SQRT(RMS_ERROR_ABS(deriv_idx)/NUMBER(deriv_idx))
                                      VALUES(3)=SQRT(RMS_ERROR_REL(deriv_idx)/NUMBER(deriv_idx))
                                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))", &
                                        & ERR,ERROR,*999)
                                      !Local and ghost nodes
                                      CALL WRITE_STRING(OUTPUT_ID,"Local + Ghost RMS errors:",ERR,ERROR,*999)
                                      LOCAL_STRING= &
                                        & "                                                     % error  Absolute err  Relative err"
                                      CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                      VALUES(1)=SQRT((RMS_ERROR_PER(1)+GHOST_RMS_ERROR_PER(1))/(NUMBER(1)+GHOST_NUMBER(1)))
                                      VALUES(2)=SQRT((RMS_ERROR_ABS(1)+GHOST_RMS_ERROR_ABS(1))/(NUMBER(1)+GHOST_NUMBER(1)))
                                      VALUES(3)=SQRT((RMS_ERROR_REL(1)+GHOST_RMS_ERROR_REL(1))/(NUMBER(1)+GHOST_NUMBER(1)))
                                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))", &
                                        & ERR,ERROR,*999)
                                      !Global RMS values
                                      !Collect the values across the ranks
                                      CALL MPI_ALLREDUCE(MPI_IN_PLACE,NUMBER,1,MPI_INTEGER,MPI_SUM, &
                                        & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                                      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                      CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR_PER,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                        & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                                      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                      CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR_ABS,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                        & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                                      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                      CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR_REL,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                        & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                                      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                      CALL WRITE_STRING(OUTPUT_ID,"Global RMS errors:",ERR,ERROR,*999)
                                      LOCAL_STRING= &
                                        & "                                                     % error  Absolute err  Relative err"
                                      CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                      VALUES(1)=SQRT(RMS_ERROR_PER(1)/NUMBER(1))
                                      VALUES(2)=SQRT(RMS_ERROR_ABS(1)/NUMBER(1))
                                      VALUES(3)=SQRT(RMS_ERROR_REL(1)/NUMBER(1))
                                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))", &
                                        & ERR,ERROR,*999)
                                    ELSE
                                      CALL WRITE_STRING(OUTPUT_ID,"RMS errors:",ERR,ERROR,*999)
                                      LOCAL_STRING= &
                                        & "                                                     % error  Absolute err  Relative err"
                                      CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                      VALUES(1)=SQRT(RMS_ERROR_PER(deriv_idx)/NUMBER(deriv_idx))
                                      VALUES(2)=SQRT(RMS_ERROR_ABS(deriv_idx)/NUMBER(deriv_idx))
                                      VALUES(3)=SQRT(RMS_ERROR_REL(deriv_idx)/NUMBER(deriv_idx))
                                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))", &
                                        & ERR,ERROR,*999)
                                    ENDIF
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Decomposition topology elements is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Decomposition topology is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain decomposition is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Elements domain topology is not associated.",ERR,ERROR,*999)
                          ENDIF                     
                       CASE(FIELD_NODE_BASED_INTERPOLATION)
                          NODES_DOMAIN=>DOMAIN_TOPOLOGY%NODES
                          IF(ASSOCIATED(NODES_DOMAIN)) THEN
                            NUMBER=0
                            RMS_ERROR_PER=0.0_DP
                            RMS_ERROR_ABS=0.0_DP
                            RMS_ERROR_REL=0.0_DP
                            GHOST_NUMBER=0
                            GHOST_RMS_ERROR_PER=0.0_DP
                            GHOST_RMS_ERROR_ABS=0.0_DP
                            GHOST_RMS_ERROR_REL=0.0_DP
                            CALL WRITE_STRING(OUTPUT_ID,"Nodal errors:",ERR,ERROR,*999)
                            LOCAL_STRING="     Node#  Deriv#     Numerical      Analytic       % error  Absolute err  Relative err"
                            CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                            DO node_idx=1,NODES_DOMAIN%NUMBER_OF_NODES
                              DO deriv_idx=1,NODES_DOMAIN%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                                VALUES(1)=NUMERICAL_VALUES(local_ny)
                                VALUES(2)=ANALYTIC_VALUES(local_ny)
                                VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                                VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                                VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                                !Accumlate the RMS errors
                                NUMBER(deriv_idx)=NUMBER(deriv_idx)+1
                                RMS_ERROR_PER(deriv_idx)=RMS_ERROR_PER(deriv_idx)+VALUES(3)*VALUES(3)
                                RMS_ERROR_ABS(deriv_idx)=RMS_ERROR_ABS(deriv_idx)+VALUES(4)*VALUES(4)
                                RMS_ERROR_REL(deriv_idx)=RMS_ERROR_REL(deriv_idx)+VALUES(5)*VALUES(5)
                                IF(deriv_idx==1) THEN
                                  WRITE(FIRST_FORMAT,"(A,I10,A,I6,A)") "('",NODES_DOMAIN%NODES(node_idx)%USER_NUMBER,"',2X,'", &
                                    & deriv_idx,"',5(2X,E12.5))"
                                ELSE
                                  WRITE(FIRST_FORMAT,"(A,I6,A)") "(12X,'",deriv_idx,"',5(2X,E12.5))"
                                ENDIF
                                CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                            DO node_idx=NODES_DOMAIN%NUMBER_OF_NODES+1,NODES_DOMAIN%TOTAL_NUMBER_OF_NODES
                              DO deriv_idx=1,NODES_DOMAIN%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                                VALUES(1)=NUMERICAL_VALUES(local_ny)
                                VALUES(2)=ANALYTIC_VALUES(local_ny)
                                VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                                VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                                VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                                !Accumlate the RMS errors
                                GHOST_NUMBER(deriv_idx)=GHOST_NUMBER(deriv_idx)+1
                                GHOST_RMS_ERROR_PER(deriv_idx)=GHOST_RMS_ERROR_PER(deriv_idx)+VALUES(3)*VALUES(3)
                                GHOST_RMS_ERROR_ABS(deriv_idx)=GHOST_RMS_ERROR_ABS(deriv_idx)+VALUES(4)*VALUES(4)
                                GHOST_RMS_ERROR_REL(deriv_idx)=GHOST_RMS_ERROR_REL(deriv_idx)+VALUES(5)*VALUES(5)
                                IF(deriv_idx==1) THEN
                                  WRITE(FIRST_FORMAT,"(A,I10,A,I6,A)") "('",NODES_DOMAIN%NODES(node_idx)%USER_NUMBER, &
                                    & "',2X,'",deriv_idx,"',5(2X,E12.5))"
                                ELSE
                                  WRITE(FIRST_FORMAT,"(A,I6,A)") "(12X,'",deriv_idx,"',5(2X,E12.5))"
                                ENDIF
                                CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                            !Output RMS errors                  
                            CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
                            IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>1) THEN
                              IF(ANY(NUMBER>0)) THEN
                                !Local nodes only
                                CALL WRITE_STRING(OUTPUT_ID,"Local RMS errors:",ERR,ERROR,*999)
                                LOCAL_STRING= &
                                  & "            Deriv#                                   % error  Absolute err  Relative err"
                                CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                DO deriv_idx=1,8
                                  IF(NUMBER(deriv_idx)>0) THEN
                                    VALUES(1)=SQRT(RMS_ERROR_PER(deriv_idx)/NUMBER(deriv_idx))
                                    VALUES(2)=SQRT(RMS_ERROR_ABS(deriv_idx)/NUMBER(deriv_idx))
                                    VALUES(3)=SQRT(RMS_ERROR_REL(deriv_idx)/NUMBER(deriv_idx))
                                    WRITE(FIRST_FORMAT,"(A,I6,A)") "(12X,'",deriv_idx,"',28X,3(2X,E12.5))"
                                    CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,FIRST_FORMAT,"(46X,3(2X,E12.5))", &
                                      & ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !deriv_idx
                                !Local and ghost nodes
                                CALL WRITE_STRING(OUTPUT_ID,"Local + Ghost RMS errors:",ERR,ERROR,*999)
                                LOCAL_STRING= &
                                  & "            Deriv#                                   % error  Absolute err  Relative err"
                                CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                DO deriv_idx=1,8
                                  IF(NUMBER(deriv_idx)>0) THEN
                                    VALUES(1)=SQRT((RMS_ERROR_PER(deriv_idx)+GHOST_RMS_ERROR_PER(deriv_idx))/ &
                                      & (NUMBER(deriv_idx)+GHOST_NUMBER(deriv_idx)))
                                    VALUES(2)=SQRT((RMS_ERROR_ABS(deriv_idx)+GHOST_RMS_ERROR_ABS(deriv_idx))/ &
                                      & (NUMBER(deriv_idx)+GHOST_NUMBER(deriv_idx)))
                                    VALUES(3)=SQRT((RMS_ERROR_REL(deriv_idx)+GHOST_RMS_ERROR_REL(deriv_idx))/ &
                                      & (NUMBER(deriv_idx)+GHOST_NUMBER(deriv_idx)))
                                    WRITE(FIRST_FORMAT,"(A,I6,A)") "(12X,'",deriv_idx,"',28X,3(2X,E12.5))"
                                    CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,FIRST_FORMAT,"(46X,3(2X,E12.5))", &
                                      & ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !deriv_idx
                                !Global RMS values
                                !Collect the values across the ranks
                                CALL MPI_ALLREDUCE(MPI_IN_PLACE,NUMBER,8,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
                                  & MPI_IERROR)
                                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR_PER,8,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                  & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR_ABS,8,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                  & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR_REL,8,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                  & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                                CALL WRITE_STRING(OUTPUT_ID,"Global RMS errors:",ERR,ERROR,*999)
                                LOCAL_STRING= &
                                  & "            Deriv#                                   % error  Absolute err  Relative err"
                                CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                DO deriv_idx=1,8
                                  IF(NUMBER(deriv_idx)>0) THEN
                                    VALUES(1)=SQRT(RMS_ERROR_PER(deriv_idx)/NUMBER(deriv_idx))
                                    VALUES(2)=SQRT(RMS_ERROR_ABS(deriv_idx)/NUMBER(deriv_idx))
                                    VALUES(3)=SQRT(RMS_ERROR_REL(deriv_idx)/NUMBER(deriv_idx))
                                    WRITE(FIRST_FORMAT,"(A,I6,A)") "(12X,'",deriv_idx,"',28X,3(2X,E12.5))"
                                    CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,FIRST_FORMAT,"(46X,3(2X,E12.5))", &
                                      & ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !deriv_idx
                              ENDIF
                            ELSE
                              IF(ANY(NUMBER>0)) THEN
                                CALL WRITE_STRING(OUTPUT_ID,"RMS errors:",ERR,ERROR,*999)
                               LOCAL_STRING= &
                                 & "            Deriv#                                   % error  Absolute err  Relative err"
                                CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                                DO deriv_idx=1,8
                                  IF(NUMBER(deriv_idx)>0) THEN
                                    VALUES(1)=SQRT(RMS_ERROR_PER(deriv_idx)/NUMBER(deriv_idx))
                                    VALUES(2)=SQRT(RMS_ERROR_ABS(deriv_idx)/NUMBER(deriv_idx))
                                    VALUES(3)=SQRT(RMS_ERROR_REL(deriv_idx)/NUMBER(deriv_idx))
                                    WRITE(FIRST_FORMAT,"(A,I6,A)") "(12X,'",deriv_idx,"',28X,3(2X,E12.5))"
                                    CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,3,3,3,VALUES,FIRST_FORMAT,"(46X,3(2X,E12.5))", &
                                      & ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !deriv_idx
                              ENDIF
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Nodes domain topology is not associated.",ERR,ERROR,*999)
                          ENDIF                     
                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The interpolation type of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD%VARIABLES(var_idx)%COMPONENTS(component_idx)% &
                            & INTERPOLATION_TYPE,"*",ERR,ERROR))//" for component number "// &
                            & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of variable type "// &
                            & TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))//" of field number "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                    ENDIF
                    CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
                  ENDDO !component_idx
                  !Restore the dependent and analytic parameter sets
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(FIELD,variable_type,FIELD_VALUES_SET_TYPE,NUMERICAL_VALUES, &
                    & ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ANALYTIC_VALUES, &
                    & ERR,ERROR,*999)
                  !Allocated the integral errors
                  ALLOCATE(INTEGRAL_ERRORS(6,FIELD_VARIABLE%NUMBER_OF_COMPONENTS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate integral errors.",ERR,ERROR,*999)
                  ALLOCATE(GHOST_INTEGRAL_ERRORS(6,FIELD_VARIABLE%NUMBER_OF_COMPONENTS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ghost integral errors.",ERR,ERROR,*999)
                  CALL ANALYTIC_ANALYSIS_INTEGRAL_ERRORS(FIELD_VARIABLE,INTEGRAL_ERRORS,GHOST_INTEGRAL_ERRORS,ERR,ERROR,*999)
                  IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>1) THEN
                    CALL WRITE_STRING(OUTPUT_ID,"Local Integral errors:",ERR,ERROR,*999)
                    LOCAL_STRING="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      VALUES(1)=INTEGRAL_ERRORS(1,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Intg  ',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(2,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                    ENDDO !component_idx
                    LOCAL_STRING="Component#             Numerical      Analytic           NID        NID(%)"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      VALUES(1)=INTEGRAL_ERRORS(5,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                     WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Diff  ',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(6,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                    ENDDO !component_idx
                    CALL WRITE_STRING(OUTPUT_ID,"Local + Ghost Integral errors:",ERR,ERROR,*999)
                    LOCAL_STRING="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      VALUES(1)=INTEGRAL_ERRORS(1,component_idx)+GHOST_INTEGRAL_ERRORS(1,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)+GHOST_INTEGRAL_ERRORS(3,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Intg  ',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(2,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                    ENDDO !component_idx
                    LOCAL_STRING="Component#             Numerical      Analytic           NID        NID(%)"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      VALUES(1)=INTEGRAL_ERRORS(5,component_idx)+GHOST_INTEGRAL_ERRORS(5,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)+GHOST_INTEGRAL_ERRORS(3,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                      WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Diff  ',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(6,component_idx)+GHOST_INTEGRAL_ERRORS(6,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)+GHOST_INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                    ENDDO !component_idx
                    !Collect the values across the ranks
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,INTEGRAL_ERRORS,6*FIELD_VARIABLE%NUMBER_OF_COMPONENTS,MPI_DOUBLE_PRECISION, &
                      & MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL WRITE_STRING(OUTPUT_ID,"Global Integral errors:",ERR,ERROR,*999)
                    LOCAL_STRING="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(1,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Intg  ',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(2,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999) 
                    ENDDO !component_idx
                    LOCAL_STRING="Component#             Numerical      Analytic           NID        NID(%)"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      VALUES(1)=INTEGRAL_ERRORS(5,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                     WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Diff  ',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(6,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                    ENDDO !component_idx
                  ELSE
                    CALL WRITE_STRING(OUTPUT_ID,"Integral errors:",ERR,ERROR,*999)
                    LOCAL_STRING="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      VALUES(1)=INTEGRAL_ERRORS(1,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Intg  ',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(2,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(VALUES(1),VALUES(2))
                      VALUES(5)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(VALUES(1),VALUES(2))
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,5,5,5,VALUES,FIRST_FORMAT,"(20X,5(2X,E12.5))",ERR,ERROR,*999)
                    ENDDO !component_idx
                    LOCAL_STRING="Component#             Numerical      Analytic           NID        NID(%)"
                    CALL WRITE_STRING(OUTPUT_ID,LOCAL_STRING,ERR,ERROR,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      VALUES(1)=INTEGRAL_ERRORS(5,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(3,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                     WRITE(FIRST_FORMAT,"(A,I10,A)") "('",component_idx,"',2X,'Diff  ',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                      VALUES(1)=INTEGRAL_ERRORS(6,component_idx)
                      VALUES(2)=INTEGRAL_ERRORS(4,component_idx)
                      VALUES(3)=ANALYTIC_ANALYSIS_NID_ERROR(VALUES(1),VALUES(2))
                      VALUES(4)=VALUES(3)*100.0_DP
                      WRITE(FIRST_FORMAT,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
                      CALL WRITE_STRING_VECTOR(OUTPUT_ID,1,1,4,4,4,VALUES,FIRST_FORMAT,"(20X,4(2X,E12.5))",ERR,ERROR,*999)
                    ENDDO !component_idx
                    CALL WRITE_STRING(OUTPUT_ID,"",ERR,ERROR,*999)
                  ENDIF
                  IF(ALLOCATED(INTEGRAL_ERRORS)) DEALLOCATE(INTEGRAL_ERRORS)
                  IF(ALLOCATED(GHOST_INTEGRAL_ERRORS)) DEALLOCATE(GHOST_INTEGRAL_ERRORS)
                ELSE
                  CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !var_idx
            ELSE
              CALL FLAG_ERROR("Decomposition topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
          IF(LEN_TRIM(FILENAME)>=1) THEN
            CLOSE(UNIT=OUTPUT_ID,IOSTAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Error closing analysis output file.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
            & " is not a dependent field."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
          
    CALL EXITS("ANALYTIC_ANALYSIS_OUTPUT")
    RETURN
999 IF(ALLOCATED(INTEGRAL_ERRORS)) DEALLOCATE(INTEGRAL_ERRORS)
    IF(ALLOCATED(GHOST_INTEGRAL_ERRORS)) DEALLOCATE(GHOST_INTEGRAL_ERRORS)
    CALL ERRORS("ANALYTIC_ANALYSIS_OUTPUT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_OUTPUT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_OUTPUT

  !
  !================================================================================================================================
  !  

  !>Calculates the absolute error between a numeric value and an analytic value.
  PURE FUNCTION ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(NUMERIC_VALUE,ANALYTIC_VALUE)

    !Argument variables
    REAL(DP), INTENT(IN) :: NUMERIC_VALUE !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: ANALYTIC_VALUE !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: ANALYTIC_ANALYSIS_ABSOLUTE_ERROR

    ANALYTIC_ANALYSIS_ABSOLUTE_ERROR=ABS(ANALYTIC_VALUE-NUMERIC_VALUE)

  END FUNCTION ANALYTIC_ANALYSIS_ABSOLUTE_ERROR
  
  !
  !================================================================================================================================
  !  

  !>Calculates the percentage error between a numeric value and an analytic value.
  PURE FUNCTION ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(NUMERIC_VALUE,ANALYTIC_VALUE)

    !Argument variables
    REAL(DP), INTENT(IN) :: NUMERIC_VALUE !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: ANALYTIC_VALUE !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: ANALYTIC_ANALYSIS_PERCENTAGE_ERROR

    IF(ABS(ANALYTIC_VALUE)>ZERO_TOLERANCE) THEN
      ANALYTIC_ANALYSIS_PERCENTAGE_ERROR=(ANALYTIC_VALUE-NUMERIC_VALUE)/ANALYTIC_VALUE*100.0_DP
    ELSE
      ANALYTIC_ANALYSIS_PERCENTAGE_ERROR=0.0_DP
    ENDIF

  END FUNCTION ANALYTIC_ANALYSIS_PERCENTAGE_ERROR
  
  !
  !================================================================================================================================
  !  

  !>Calculates the relative error between a numeric value and an analytic value.
  PURE FUNCTION ANALYTIC_ANALYSIS_RELATIVE_ERROR(NUMERIC_VALUE,ANALYTIC_VALUE)

    !Argument variables
    REAL(DP), INTENT(IN) :: NUMERIC_VALUE !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: ANALYTIC_VALUE !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: ANALYTIC_ANALYSIS_RELATIVE_ERROR

    IF(ABS(1.0_DP+ANALYTIC_VALUE)>ZERO_TOLERANCE) THEN
      ANALYTIC_ANALYSIS_RELATIVE_ERROR=(ANALYTIC_VALUE-NUMERIC_VALUE)/(1.0_DP+ANALYTIC_VALUE)
    ELSE
      ANALYTIC_ANALYSIS_RELATIVE_ERROR=0.0_DP
    ENDIF

  END FUNCTION ANALYTIC_ANALYSIS_RELATIVE_ERROR
  
  !
  !================================================================================================================================
  !  

  !>Calculates the Normalised Integral Difference (NID) error with a numeric value and an analytic value.
  PURE FUNCTION ANALYTIC_ANALYSIS_NID_ERROR(NUMERIC_VALUE,ANALYTIC_VALUE)

    !Argument variables
    REAL(DP), INTENT(IN) :: NUMERIC_VALUE !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: ANALYTIC_VALUE !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: ANALYTIC_ANALYSIS_NID_ERROR

    IF(ABS(ANALYTIC_VALUE)>ZERO_TOLERANCE) THEN
      ANALYTIC_ANALYSIS_NID_ERROR=(ANALYTIC_VALUE-NUMERIC_VALUE)/ANALYTIC_VALUE
    ELSE
      ANALYTIC_ANALYSIS_NID_ERROR=(ANALYTIC_VALUE-NUMERIC_VALUE)/(1.0_DP+ANALYTIC_VALUE)
    ENDIF

  END FUNCTION ANALYTIC_ANALYSIS_NID_ERROR
  
  !
  !================================================================================================================================
  !  

  !>Calculates the relative error between a numeric value and an analytic value.
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ERRORS(FIELD_VARIABLE,INTEGRAL_ERRORS,GHOST_INTEGRAL_ERRORS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable to calculate the integral errors for
    REAL(DP), INTENT(OUT) :: INTEGRAL_ERRORS(:,:) !<On exit, the integral errors for the local elements
    REAL(DP), INTENT(OUT) :: GHOST_INTEGRAL_ERRORS(:,:) !<On exit, the integral errors for the ghost elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,component_idx,gauss_idx,parameter_idx,variable_type
    REAL(DP) :: ANALYTIC_INT,NUMERICAL_INT,RWG
    TYPE(BASIS_TYPE), POINTER :: BASIS,DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS1,DOMAIN_ELEMENTS2,DOMAIN_ELEMENTS3
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERP_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERP_POINT_METRICS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: ANALYTIC_INTERP_PARAMETERS,GEOMETRIC_INTERP_PARAMETERS, &
      & NUMERICAL_INTERP_PARAMETERS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(GEOMETRIC_INTERP_POINT)
    NULLIFY(GEOMETRIC_INTERP_POINT_METRICS)
    NULLIFY(ANALYTIC_INTERP_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERP_PARAMETERS)
    NULLIFY(NUMERICAL_INTERP_PARAMETERS)

    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_ERRORS",ERR,ERROR,*999)
    
    INTEGRAL_ERRORS=0.0_DP
    GHOST_INTEGRAL_ERRORS=0.0_DP
    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
      IF(SIZE(INTEGRAL_ERRORS,1)>=6.AND.SIZE(INTEGRAL_ERRORS,2)>=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
        IF(SIZE(GHOST_INTEGRAL_ERRORS,1)>=6.AND.SIZE(GHOST_INTEGRAL_ERRORS,2)>=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
          variable_type=FIELD_VARIABLE%VARIABLE_TYPE
          DEPENDENT_FIELD=>FIELD_VARIABLE%FIELD
          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
            DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
            IF(ASSOCIATED(DECOMPOSITION)) THEN
              GEOMETRIC_FIELD=>DEPENDENT_FIELD%GEOMETRIC_FIELD
              IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                IF(ASSOCIATED(GEOMETRIC_VARIABLE)) THEN
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DEPENDENT_FIELD,variable_type,NUMERICAL_INTERP_PARAMETERS, &
                    & ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DEPENDENT_FIELD,variable_type,ANALYTIC_INTERP_PARAMETERS, &
                    & ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATED_POINT_INITIALISE(GEOMETRIC_INTERP_PARAMETERS,GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(GEOMETRIC_INTERP_POINT,GEOMETRIC_INTERP_POINT_METRICS, &
                    & ERR,ERROR,*999)
                  DOMAIN_ELEMENTS1=>FIELD_VARIABLE%COMPONENTS(DECOMPOSITION%MESH_COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%ELEMENTS
                  DOMAIN_ELEMENTS2=>GEOMETRIC_VARIABLE%COMPONENTS(DECOMPOSITION%MESH_COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%ELEMENTS
                  DO element_idx=1,DOMAIN_ELEMENTS1%NUMBER_OF_ELEMENTS
                    DEPENDENT_BASIS=>DOMAIN_ELEMENTS1%ELEMENTS(element_idx)%BASIS
                    GEOMETRIC_BASIS=>DOMAIN_ELEMENTS2%ELEMENTS(element_idx)%BASIS
                    QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
                      & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
                      & NUMERICAL_INTERP_PARAMETERS,ERR,ERROR,*999)
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_ANALYTIC_VALUES_SET_TYPE,element_idx, &
                      & ANALYTIC_INTERP_PARAMETERS,ERR,ERROR,*999)
                    DO gauss_idx=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
                      CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                        & GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
                      CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERP_POINT_METRICS, &
                        & ERR,ERROR,*999)
                      RWG=GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN*QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                        DOMAIN_ELEMENTS3=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
                        BASIS=>DOMAIN_ELEMENTS3%ELEMENTS(element_idx)%BASIS
                        NUMERICAL_INT=0.0_DP
                        ANALYTIC_INT=0.0_DP
                        SELECT CASE(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE)
                        CASE(FIELD_NO_SCALING)
                          DO parameter_idx=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            NUMERICAL_INT=NUMERICAL_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & NUMERICAL_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)
                            ANALYTIC_INT=ANALYTIC_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & ANALYTIC_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)
                          ENDDO !parameter_idx
                        CASE(FIELD_UNIT_SCALING,FIELD_ARC_LENGTH_SCALING,FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
                          DO parameter_idx=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            NUMERICAL_INT=NUMERICAL_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & NUMERICAL_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)* &
                              & NUMERICAL_INTERP_PARAMETERS%SCALE_FACTORS(parameter_idx,component_idx)
                            ANALYTIC_INT=ANALYTIC_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & ANALYTIC_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)* &
                              & ANALYTIC_INTERP_PARAMETERS%SCALE_FACTORS(parameter_idx,component_idx)
                          ENDDO !parameter_idx
                        CASE DEFAULT
                          LOCAL_ERROR="The dependent field scaling type of "// &
                            & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                        INTEGRAL_ERRORS(1,component_idx)=INTEGRAL_ERRORS(1,component_idx)+NUMERICAL_INT*RWG
                        INTEGRAL_ERRORS(2,component_idx)=INTEGRAL_ERRORS(2,component_idx)+NUMERICAL_INT**2*RWG
                        INTEGRAL_ERRORS(3,component_idx)=INTEGRAL_ERRORS(3,component_idx)+ANALYTIC_INT*RWG
                        INTEGRAL_ERRORS(4,component_idx)=INTEGRAL_ERRORS(4,component_idx)+ANALYTIC_INT**2*RWG
                        INTEGRAL_ERRORS(5,component_idx)=INTEGRAL_ERRORS(5,component_idx)+(ANALYTIC_INT-NUMERICAL_INT)*RWG
                        INTEGRAL_ERRORS(6,component_idx)=INTEGRAL_ERRORS(6,component_idx)+(ANALYTIC_INT-NUMERICAL_INT)**2*RWG
                      ENDDO !component_idx
                    ENDDO !gauss_idx
                  ENDDO !element_idx
                  DO element_idx=DOMAIN_ELEMENTS1%NUMBER_OF_ELEMENTS+1,DOMAIN_ELEMENTS1%TOTAL_NUMBER_OF_ELEMENTS
                    DEPENDENT_BASIS=>DOMAIN_ELEMENTS1%ELEMENTS(element_idx)%BASIS
                    GEOMETRIC_BASIS=>DOMAIN_ELEMENTS2%ELEMENTS(element_idx)%BASIS
                    QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx,GEOMETRIC_INTERP_PARAMETERS, &
                      & ERR,ERROR,*999)
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx,NUMERICAL_INTERP_PARAMETERS, &
                      & ERR,ERROR,*999)
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_ANALYTIC_VALUES_SET_TYPE,element_idx, &
                      & ANALYTIC_INTERP_PARAMETERS,ERR,ERROR,*999)
                    DO gauss_idx=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
                      CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                        & GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
                      CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERP_POINT_METRICS, &
                        & ERR,ERROR,*999)
                      RWG=GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN*QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
                      DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                        DOMAIN_ELEMENTS3=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
                        BASIS=>DOMAIN_ELEMENTS3%ELEMENTS(element_idx)%BASIS
                        NUMERICAL_INT=0.0_DP
                        ANALYTIC_INT=0.0_DP
                        SELECT CASE(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE)
                        CASE(FIELD_NO_SCALING)
                          DO parameter_idx=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            NUMERICAL_INT=NUMERICAL_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & NUMERICAL_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)
                            ANALYTIC_INT=ANALYTIC_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & ANALYTIC_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)
                          ENDDO !parameter_idx
                        CASE(FIELD_UNIT_SCALING,FIELD_ARC_LENGTH_SCALING,FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
                          DO parameter_idx=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            NUMERICAL_INT=NUMERICAL_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & NUMERICAL_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)* &
                              & NUMERICAL_INTERP_PARAMETERS%SCALE_FACTORS(parameter_idx,component_idx)
                            ANALYTIC_INT=ANALYTIC_INT+QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)* &
                              & ANALYTIC_INTERP_PARAMETERS%PARAMETERS(parameter_idx,component_idx)* &
                              & ANALYTIC_INTERP_PARAMETERS%SCALE_FACTORS(parameter_idx,component_idx)
                          ENDDO !parameter_idx
                        CASE DEFAULT
                          LOCAL_ERROR="The dependent field scaling type of "// &
                            & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                        GHOST_INTEGRAL_ERRORS(1,component_idx)=GHOST_INTEGRAL_ERRORS(1,component_idx)+NUMERICAL_INT*RWG
                        GHOST_INTEGRAL_ERRORS(2,component_idx)=GHOST_INTEGRAL_ERRORS(2,component_idx)+NUMERICAL_INT**2*RWG
                        GHOST_INTEGRAL_ERRORS(3,component_idx)=GHOST_INTEGRAL_ERRORS(3,component_idx)+ANALYTIC_INT*RWG
                        GHOST_INTEGRAL_ERRORS(4,component_idx)=GHOST_INTEGRAL_ERRORS(4,component_idx)+ANALYTIC_INT**2*RWG
                        GHOST_INTEGRAL_ERRORS(5,component_idx)=GHOST_INTEGRAL_ERRORS(5,component_idx)+ &
                          & (ANALYTIC_INT-NUMERICAL_INT)*RWG
                        GHOST_INTEGRAL_ERRORS(6,component_idx)=GHOST_INTEGRAL_ERRORS(6,component_idx)+ &
                          & (ANALYTIC_INT-NUMERICAL_INT)**2*RWG
                      ENDDO !component_idx
                    ENDDO !gauss_idx
                  ENDDO !element_idx
                  CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATED_POINT_FINALISE(GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(ANALYTIC_INTERP_PARAMETERS,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(NUMERICAL_INTERP_PARAMETERS,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Geometric field variable is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Field geometric field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Field variable field is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Invalid size for GHOST_INTEGRAL_ERRORS. The size is ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GHOST_INTEGRAL_ERRORS,1),"*",ERR,ERROR))//","// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(GHOST_INTEGRAL_ERRORS,2),"*",ERR,ERROR))//") and it needs to be at least (6,"// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Invalid size for INTEGRAL_ERRORS. The size is ("// &
          & TRIM(NUMBER_TO_VSTRING(SIZE(INTEGRAL_ERRORS,1),"*",ERR,ERROR))//","// &
          & TRIM(NUMBER_TO_VSTRING(SIZE(INTEGRAL_ERRORS,2),"*",ERR,ERROR))//") and it needs to be at least (6,"// &
          & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//")."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ERRORS")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_ERRORS",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ERRORS")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ERRORS

  !
  !================================================================================================================================
  !

  !>Get integral absolute error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTEGRAL_ERROR, &
    & GHOST_INTEGRAL_ERROR,ERR,ERROR,*)

    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    REAL(DP), INTENT(OUT) :: INTEGRAL_ERROR(2) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: GHOST_INTEGRAL_ERROR(2) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP), ALLOCATABLE :: INTEGRAL_ERRORS(:,:) !<the integral errors for the local elements
    REAL(DP), ALLOCATABLE :: GHOST_INTEGRAL_ERRORS(:,:) !<the integral errors for the ghost elements

    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ERRORS(FIELD_VARIABLE,INTEGRAL_ERRORS,GHOST_INTEGRAL_ERRORS,ERR,ERROR,*999)
      INTEGRAL_ERROR(1)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(INTEGRAL_ERRORS(1,COMPONENT_NUMBER),INTEGRAL_ERRORS(3,COMPONENT_NUMBER))
      INTEGRAL_ERROR(2)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(INTEGRAL_ERRORS(2,COMPONENT_NUMBER),INTEGRAL_ERRORS(4,COMPONENT_NUMBER))
      GHOST_INTEGRAL_ERROR(1)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(GHOST_INTEGRAL_ERRORS(1,COMPONENT_NUMBER), &
        & GHOST_INTEGRAL_ERRORS(3,COMPONENT_NUMBER))
      GHOST_INTEGRAL_ERROR(2)=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(GHOST_INTEGRAL_ERRORS(2,COMPONENT_NUMBER), &
        & GHOST_INTEGRAL_ERRORS(4,COMPONENT_NUMBER))
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET
  
   !
  !================================================================================================================================
  !

  !>Get integral analytic value for the field TODO should we use analytical formula to calculate the integration?
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTEGRAL_ERROR, &
    & GHOST_INTEGRAL_ERROR,ERR,ERROR,*)
  
    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    REAL(DP), INTENT(OUT) :: INTEGRAL_ERROR(2) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: GHOST_INTEGRAL_ERROR(2) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP), ALLOCATABLE :: INTEGRAL_ERRORS(:,:) !<the integral errors for the local elements
    REAL(DP), ALLOCATABLE :: GHOST_INTEGRAL_ERRORS(:,:) !<the integral errors for the ghost elements

    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ERRORS(FIELD_VARIABLE,INTEGRAL_ERRORS,GHOST_INTEGRAL_ERRORS,ERR,ERROR,*999)
      INTEGRAL_ERROR(1)=INTEGRAL_ERRORS(3,COMPONENT_NUMBER)
      INTEGRAL_ERROR(2)=INTEGRAL_ERRORS(4,COMPONENT_NUMBER)
      GHOST_INTEGRAL_ERROR(1)=GHOST_INTEGRAL_ERRORS(3,COMPONENT_NUMBER)
      GHOST_INTEGRAL_ERROR(2)=GHOST_INTEGRAL_ERRORS(4,COMPONENT_NUMBER)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET
  
  !
  !================================================================================================================================
  !

  !>Get integral numerical value for the field, TODO check integral calculation
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTEGRAL_ERROR, &
    & GHOST_INTEGRAL_ERROR,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    REAL(DP), INTENT(OUT) :: INTEGRAL_ERROR(2) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: GHOST_INTEGRAL_ERROR(2) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP), ALLOCATABLE :: INTEGRAL_ERRORS(:,:) !<the integral errors for the local elements
    REAL(DP), ALLOCATABLE :: GHOST_INTEGRAL_ERRORS(:,:) !<the integral errors for the ghost elements

    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ERRORS(FIELD_VARIABLE,INTEGRAL_ERRORS,GHOST_INTEGRAL_ERRORS,ERR,ERROR,*999)
      INTEGRAL_ERROR(1)=INTEGRAL_ERRORS(1,COMPONENT_NUMBER)
      INTEGRAL_ERROR(2)=INTEGRAL_ERRORS(2,COMPONENT_NUMBER)
      GHOST_INTEGRAL_ERROR(1)=GHOST_INTEGRAL_ERRORS(1,COMPONENT_NUMBER)
      GHOST_INTEGRAL_ERROR(2)=GHOST_INTEGRAL_ERRORS(2,COMPONENT_NUMBER)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET
  
  !
  !================================================================================================================================
  !

  !>Get integral percentage error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_PERCENTAGE_ERROR_GET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTEGRAL_ERROR, &
    & GHOST_INTEGRAL_ERROR,ERR,ERROR,*)

    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    REAL(DP), INTENT(OUT) :: INTEGRAL_ERROR(2) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: GHOST_INTEGRAL_ERROR(2) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP), ALLOCATABLE :: INTEGRAL_ERRORS(:,:) !<the integral errors for the local elements
    REAL(DP), ALLOCATABLE :: GHOST_INTEGRAL_ERRORS(:,:) !<the integral errors for the ghost elements

    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENTAGE_ERROR_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ERRORS(FIELD_VARIABLE,INTEGRAL_ERRORS,GHOST_INTEGRAL_ERRORS,ERR,ERROR,*999)
      INTEGRAL_ERROR(1)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(INTEGRAL_ERRORS(1,COMPONENT_NUMBER),INTEGRAL_ERRORS(3,COMPONENT_NUMBER))
      INTEGRAL_ERROR(2)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(INTEGRAL_ERRORS(2,COMPONENT_NUMBER),INTEGRAL_ERRORS(4,COMPONENT_NUMBER))
      GHOST_INTEGRAL_ERROR(1)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(GHOST_INTEGRAL_ERRORS(1,COMPONENT_NUMBER), &
        & GHOST_INTEGRAL_ERRORS(3,COMPONENT_NUMBER))
      GHOST_INTEGRAL_ERROR(2)=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(GHOST_INTEGRAL_ERRORS(2,COMPONENT_NUMBER), &
        & GHOST_INTEGRAL_ERRORS(4,COMPONENT_NUMBER))
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENTAGE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENTAGE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENTAGE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_PERCENTAGE_ERROR_GET
  
   !
  !================================================================================================================================
  !

  !>Get integral relative error value for the field.
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTEGRAL_ERROR, &
    & GHOST_INTEGRAL_ERROR,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    REAL(DP), INTENT(OUT) :: INTEGRAL_ERROR(2) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: GHOST_INTEGRAL_ERROR(2) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP), ALLOCATABLE :: INTEGRAL_ERRORS(:,:) !<the integral errors for the local elements
    REAL(DP), ALLOCATABLE :: GHOST_INTEGRAL_ERRORS(:,:) !<the integral errors for the ghost elements

    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ERRORS(FIELD_VARIABLE,INTEGRAL_ERRORS,GHOST_INTEGRAL_ERRORS,ERR,ERROR,*999)
      INTEGRAL_ERROR(1)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(INTEGRAL_ERRORS(1,COMPONENT_NUMBER),INTEGRAL_ERRORS(3,COMPONENT_NUMBER))
      INTEGRAL_ERROR(2)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(INTEGRAL_ERRORS(2,COMPONENT_NUMBER),INTEGRAL_ERRORS(4,COMPONENT_NUMBER))
      GHOST_INTEGRAL_ERROR(1)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(GHOST_INTEGRAL_ERRORS(1,COMPONENT_NUMBER), &
        & GHOST_INTEGRAL_ERRORS(3,COMPONENT_NUMBER))
      GHOST_INTEGRAL_ERROR(2)=ANALYTIC_ANALYSIS_RELATIVE_ERROR(GHOST_INTEGRAL_ERRORS(2,COMPONENT_NUMBER), &
        & GHOST_INTEGRAL_ERRORS(4,COMPONENT_NUMBER))
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get integral absolute error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the integral absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE,h,k
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx,dev_idx,NUMBER_OF_SURROUNDING_ELEMENTS
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,SIZE(DERIVATIVE_NUMBERS)
            CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_NUMBER,FIELD_VALUES_SET_TYPE,DERIVATIVE_NUMBERS(dev_idx), &
              & node_idx,comp_idx,NUMERICAL_VALUE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_NUMBER,FIELD_ANALYTIC_VALUES_SET_TYPE,DERIVATIVE_NUMBERS(dev_idx), &
              & node_idx,comp_idx,ANALYTIC_VALUE,ERR,ERROR,*999)
            NUMBER_OF_SURROUNDING_ELEMENTS=FIELD%VARIABLES(1)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES%NODES(node_idx)% &
                & NUMBER_OF_SURROUNDING_ELEMENTS
            ! Uses Trapezoidal 2D Rule
!TODO implement integration calculation for 3D
!TODO what numerical method does cm use??
            REGULAR_MESH=>FIELD%REGION%MESHES%MESHES(1)%PTR%GENERATED_MESH%REGULAR_MESH
            h=REGULAR_MESH%MAXIMUM_EXTENT(1)/REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(1)
            k=REGULAR_MESH%MAXIMUM_EXTENT(2)/REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(2)
            SELECT CASE(POWER)
            CASE(1)
              VALUE=VALUE+(NUMERICAL_VALUE-ANALYTIC_VALUE)*NUMBER_OF_SURROUNDING_ELEMENTS*h*k/4
            CASE(2)
              VALUE=VALUE+(NUMERICAL_VALUE-ANALYTIC_VALUE)**2*NUMBER_OF_SURROUNDING_ELEMENTS*h*k/4
            CASE DEFAULT
              CALL FLAG_ERROR("Not valid power number",ERR,ERROR,*999)
            END SELECT
          ENDDO !dev_idx
        ENDDO !node_idx
      ENDDO !comp_idx
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get integral absolute error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_NIDS_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the integral absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, NUMERICAL_ERROR
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NIDS_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      CALL ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,NUMERICAL_ERROR,ERR,ERROR,*999)
!      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,NUMERICAL_VALUE,ERR, &
!        & ERROR,*999)
      SELECT CASE(POWER)
      CASE(1)
        VALUE=NUMERICAL_ERROR/NUMERICAL_VALUE
      CASE(2)
        VALUE=SQRT(NUMERICAL_ERROR/NUMERICAL_VALUE)
      CASE DEFAULT
        CALL FLAG_ERROR("Not valid power number",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NIDS_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NIDS_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NIDS_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NIDS_ERROR_GET

  !
  !================================================================================================================================
  !

  !>Get absolute error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
        & COMPONENT_NUMBER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
        & COMPONENT_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(NUMERICAL_VALUE,ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE

  !
  !================================================================================================================================
  !

  !>Get percentage error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
        & COMPONENT_NUMBER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
        & COMPONENT_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(NUMERICAL_VALUE, ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE



  !
  !================================================================================================================================
  !

  !>Get relative error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
        & COMPONENT_NUMBER,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_NODE(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
        & COMPONENT_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_RELATIVE_ERROR(NUMERICAL_VALUE, ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE

  !
  !================================================================================================================================
  !

  !>Get absolute error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER,VALUE,ERR, &
    & ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_ELEMENT(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
        & NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_ELEMENT(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,USER_ELEMENT_NUMBER, &
        & COMPONENT_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(NUMERICAL_VALUE,ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT

  !
  !================================================================================================================================
  !

  !>Get percentage error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE, &
    & USER_ELEMENT_NUMBER,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_ELEMENT(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
        & NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_ELEMENT(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,USER_ELEMENT_NUMBER, &
        & COMPONENT_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(NUMERICAL_VALUE, ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT


  !
  !================================================================================================================================
  !

  !>Get relative error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER,VALUE,ERR, &
    & ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_ELEMENT(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
        & NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_ELEMENT(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,USER_ELEMENT_NUMBER, &
        & COMPONENT_NUMBER,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_RELATIVE_ERROR(NUMERICAL_VALUE, ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT

    !
  !================================================================================================================================
  !

  !>Get absolute error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,VALUE,ERR, &
    & ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_CONSTANT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,NUMERICAL_VALUE,ERR,ERROR, &
        & *999)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,COMPONENT_NUMBER,ANALYTIC_VALUE, &
        & ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_ABSOLUTE_ERROR(NUMERICAL_VALUE,ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_CONSTANT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_CONSTANT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_CONSTANT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_CONSTANT

  !
  !================================================================================================================================
  !

  !>Get percentage error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_CONSTANT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,NUMERICAL_VALUE,ERR,ERROR, &
        & *999)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,COMPONENT_NUMBER,ANALYTIC_VALUE, &
        & ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_PERCENTAGE_ERROR(NUMERICAL_VALUE, ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_CONSTANT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_CONSTANT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_CONSTANT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_CONSTANT


  !
  !================================================================================================================================
  !

  !>Get relative error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,VALUE,ERR, &
    & ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE

    CALL ENTERS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_CONSTANT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,NUMERICAL_VALUE,ERR,ERROR, &
        & *999)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(FIELD,VARIABLE_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,COMPONENT_NUMBER,ANALYTIC_VALUE, &
        & ERR,ERROR,*999)
      VALUE=ANALYTIC_ANALYSIS_RELATIVE_ERROR(NUMERICAL_VALUE, ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_CONSTANT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_CONSTANT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_CONSTANT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_CONSTANT

  !
  !================================================================================================================================
  !

  !>Get rms error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERROR_TYPE,LOCAL_RMS,LOCAL_GHOST_RMS, &
    & GLOBAL_RMS,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component index
    INTEGER(INTG), INTENT(IN) :: ERROR_TYPE !<error type
    REAL(DP), INTENT(OUT) :: LOCAL_RMS(8) !<On return, the local rms percentage error
    REAL(DP), INTENT(OUT) :: LOCAL_GHOST_RMS(8) !<On return, the local + ghost rms percentage error
    REAL(DP), INTENT(OUT) :: GLOBAL_RMS(8) !<On return, the global rms percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: ERROR_VALUE
    INTEGER(INTG) :: GHOST_NUMBER(8),NUMBER(8),MPI_IERROR
    REAL(DP) :: RMS_ERROR(8),GHOST_RMS_ERROR(8)
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES_DOMAIN
    INTEGER(INTG) :: node_idx,deriv_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_PERCENTAGE_ERROR_GET_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      NODES_DOMAIN=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES
      IF(ASSOCIATED(NODES_DOMAIN)) THEN
        NUMBER=0
        RMS_ERROR=0.0_DP
        GHOST_NUMBER=0
        GHOST_RMS_ERROR=0.0_DP
        DO node_idx=1,NODES_DOMAIN%NUMBER_OF_NODES
          DO deriv_idx=1,NODES_DOMAIN%NODES(node_idx)%NUMBER_OF_DERIVATIVES
            SELECT CASE(ERROR_TYPE)
            CASE(ABSOLUTE_ERROR_TYPE)
              CALL ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,deriv_idx,node_idx,COMPONENT_NUMBER, &
                & ERROR_VALUE,ERR,ERROR,*999)
            CASE(PERCENTAGE_ERROR_TYPE)
              CALL ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,deriv_idx,node_idx,COMPONENT_NUMBER, &
                & ERROR_VALUE,ERR,ERROR,*999)
            CASE(RELATIVE_ERROR_TYPE)
              CALL ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,deriv_idx,node_idx,COMPONENT_NUMBER, &
                & ERROR_VALUE,ERR,ERROR,*999)
            CASE DEFAULT
              CALL FLAG_ERROR("The error type is not valid!",ERR,ERROR,*999)
            END SELECT
            !Accumlate the RMS errors
            NUMBER(deriv_idx)=NUMBER(deriv_idx)+1
            RMS_ERROR(deriv_idx)=RMS_ERROR(deriv_idx)+ERROR_VALUE*ERROR_VALUE
          ENDDO !deriv_idx
        ENDDO !node_idx
        DO node_idx=NODES_DOMAIN%NUMBER_OF_NODES+1,NODES_DOMAIN%TOTAL_NUMBER_OF_NODES
          DO deriv_idx=1,NODES_DOMAIN%NODES(node_idx)%NUMBER_OF_DERIVATIVES
            SELECT CASE(ERROR_TYPE)
            CASE(ABSOLUTE_ERROR_TYPE)
              CALL ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,deriv_idx,node_idx,COMPONENT_NUMBER, &
                & ERROR_VALUE,ERR,ERROR,*999)
            CASE(PERCENTAGE_ERROR_TYPE)
              CALL ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,deriv_idx,node_idx,COMPONENT_NUMBER, &
                & ERROR_VALUE,ERR,ERROR,*999)
            CASE(RELATIVE_ERROR_TYPE)
              CALL ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_NODE(FIELD,VARIABLE_TYPE,deriv_idx,node_idx,COMPONENT_NUMBER, &
                & ERROR_VALUE,ERR,ERROR,*999)
            CASE DEFAULT
              CALL FLAG_ERROR("The error type is not valid!",ERR,ERROR,*999)
            END SELECT
            !Accumlate the RMS errors
            GHOST_NUMBER(deriv_idx)=GHOST_NUMBER(deriv_idx)+1
            GHOST_RMS_ERROR(deriv_idx)=GHOST_RMS_ERROR(deriv_idx)+ERROR_VALUE*ERROR_VALUE
          ENDDO !deriv_idx
        ENDDO !node_idx

        IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>1) THEN
          IF(ANY(NUMBER>0)) THEN
            DO deriv_idx=1,8
              IF(NUMBER(deriv_idx)>0) THEN
                LOCAL_RMS(deriv_idx)=SQRT(RMS_ERROR(deriv_idx)/NUMBER(deriv_idx))
              ENDIF
            ENDDO !deriv_idx
            DO deriv_idx=1,8
              IF(NUMBER(deriv_idx)>0) THEN
                LOCAL_GHOST_RMS(deriv_idx)=SQRT((RMS_ERROR(deriv_idx)+GHOST_RMS_ERROR(deriv_idx))/(NUMBER(deriv_idx) &
                      & +GHOST_NUMBER(deriv_idx)))
              ENDIF
            ENDDO !deriv_idx
            !Global RMS values
            !Collect the values across the ranks
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,NUMBER,8,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR,8,MPI_DOUBLE_PRECISION,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
              & MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
            DO deriv_idx=1,8
              IF(NUMBER(deriv_idx)>0) THEN
                GLOBAL_RMS(deriv_idx)=SQRT(RMS_ERROR(deriv_idx)/NUMBER(deriv_idx))
              ENDIF
            ENDDO !deriv_idx
          ENDIF
        ELSE
          IF(ANY(NUMBER>0)) THEN
            DO deriv_idx=1,8
              IF(NUMBER(deriv_idx)>0) THEN
                LOCAL_RMS(deriv_idx)=SQRT(RMS_ERROR(deriv_idx)/NUMBER(deriv_idx))
                GLOBAL_RMS(deriv_idx)=LOCAL_RMS(deriv_idx)
              ENDIF
            ENDDO !deriv_idx
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Nodes domain topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ERROR_GET_NODE")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_ERROR_GET_NODE",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ERROR_GET_NODE")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_ERROR_GET_NODE

  !
  !================================================================================================================================
  !

  !>Get rms error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERROR_TYPE,LOCAL_RMS,LOCAL_GHOST_RMS, &
    & GLOBAL_RMS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<variable type
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component index
    INTEGER(INTG), INTENT(IN) :: ERROR_TYPE !<error type
    REAL(DP), INTENT(OUT) :: LOCAL_RMS !<On return, the local rms percentage error
    REAL(DP), INTENT(OUT) :: LOCAL_GHOST_RMS !<On return, the local + ghost rms percentage error
    REAL(DP), INTENT(OUT) :: GLOBAL_RMS !<On return, the global rms percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: ERROR_VALUE
    INTEGER(INTG) :: GHOST_NUMBER,NUMBER,MPI_IERROR
    REAL(DP) :: RMS_ERROR,GHOST_RMS_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: ELEMENTS_DECOMPOSITION
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: DECOMPOSITION_TOPOLOGY
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_DOMAIN
    INTEGER(INTG) :: element_idx

    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_ERROR_GET_ELEMENT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      DOMAIN=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR%COMPONENTS(COMPONENT_NUMBER)%DOMAIN
      ELEMENTS_DOMAIN=>DOMAIN%TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(ELEMENTS_DOMAIN)) THEN
        DECOMPOSITION=>DOMAIN%DECOMPOSITION
        IF(ASSOCIATED(DECOMPOSITION)) THEN
          DECOMPOSITION_TOPOLOGY=>DECOMPOSITION%TOPOLOGY
          IF(ASSOCIATED(DECOMPOSITION_TOPOLOGY)) THEN
            ELEMENTS_DECOMPOSITION=>DECOMPOSITION_TOPOLOGY%ELEMENTS
            IF(ASSOCIATED(ELEMENTS_DECOMPOSITION)) THEN
              NUMBER=0
              RMS_ERROR=0.0_DP
              GHOST_NUMBER=0
              GHOST_RMS_ERROR=0.0_DP
              DO element_idx=1,ELEMENTS_DOMAIN%NUMBER_OF_ELEMENTS
                SELECT CASE(ERROR_TYPE)
                CASE(ABSOLUTE_ERROR_TYPE)
                  CALL ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,element_idx,COMPONENT_NUMBER,ERROR_VALUE, &
                    & ERR,ERROR,*999)
                CASE(PERCENTAGE_ERROR_TYPE)
                  CALL ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,element_idx,COMPONENT_NUMBER, &
                    & ERROR_VALUE,ERR,ERROR,*999)
                CASE(RELATIVE_ERROR_TYPE)
                  CALL ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,element_idx,COMPONENT_NUMBER,ERROR_VALUE, &
                    & ERR,ERROR,*999)
                CASE DEFAULT
                  CALL FLAG_ERROR("The error type is not valid!",ERR,ERROR,*999)
                END SELECT
                NUMBER=NUMBER+1
                RMS_ERROR=RMS_ERROR+ERROR_VALUE*ERROR_VALUE
              ENDDO !element_idx
              DO element_idx=ELEMENTS_DOMAIN%NUMBER_OF_ELEMENTS+1,ELEMENTS_DOMAIN%TOTAL_NUMBER_OF_ELEMENTS
                SELECT CASE(ERROR_TYPE)
                CASE(ABSOLUTE_ERROR_TYPE)
                  CALL ANALYTIC_ANALYSIS_ABSOLUTE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,element_idx,COMPONENT_NUMBER,ERROR_VALUE, &
                    & ERR,ERROR,*999)
                CASE(PERCENTAGE_ERROR_TYPE)
                  CALL ANALYTIC_ANALYSIS_PERCENTAGE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,element_idx,COMPONENT_NUMBER, &
                    & ERROR_VALUE,ERR,ERROR,*999)
                CASE(RELATIVE_ERROR_TYPE)
                  CALL ANALYTIC_ANALYSIS_RELATIVE_ERROR_GET_ELEMENT(FIELD,VARIABLE_TYPE,element_idx,COMPONENT_NUMBER,ERROR_VALUE, &
                    & ERR,ERROR,*999)
                CASE DEFAULT
                  CALL FLAG_ERROR("The error type is not valid!",ERR,ERROR,*999)
                END SELECT
                GHOST_NUMBER=GHOST_NUMBER+1
                GHOST_RMS_ERROR=GHOST_RMS_ERROR+ERROR_VALUE*ERROR_VALUE
              ENDDO !element_idx
              IF(NUMBER>0) THEN
                IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>1) THEN
                  !Local elements only
                  LOCAL_RMS=SQRT(RMS_ERROR/NUMBER)
                  !Local and ghost elements
                  LOCAL_GHOST_RMS=SQRT((RMS_ERROR+GHOST_RMS_ERROR)/(NUMBER+GHOST_NUMBER))
                  !Global RMS values
                  !Collect the values across the ranks
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,NUMBER,1,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,RMS_ERROR,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
                    & MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  GLOBAL_RMS=SQRT(RMS_ERROR/NUMBER)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Decomposition topology elements is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Decomposition topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Domain decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Elements domain topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ERROR_GET_ELEMENT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_ERROR_GET_ELEMENT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ERROR_GET_ELEMENT")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_ERROR_GET_ELEMENT

  !
  !================================================================================================================================
  !  
  
 
END MODULE ANALYTIC_ANALYSIS_ROUTINES
