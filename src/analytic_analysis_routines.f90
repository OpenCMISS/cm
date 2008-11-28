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

!>This module handles all analytic analysis routines.
MODULE ANALYTIC_ANALYSIS_ROUTINES

  USE CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC ANALYTIC_ANALYSIS_EXPORT
  
  PUBLIC ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET,ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET, &
    & ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET,ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET,ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET
 
  PUBLIC ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET,ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET,ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET
  
  PUBLIC ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET,ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET, &
    & ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET,ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET, &
    & ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET
    
CONTAINS  

  !
  !================================================================================================================================
  !  

  !>Calculate the analytic analysis data. 
  SUBROUTINE ANALYTIC_ANALYSIS_CALCULATE(FIELD,FILE_ID,ERR,ERROR,*)
  
    !Argument variables 
    TYPE(FIELD_TYPE), INTENT(IN), POINTER :: FIELD  
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: var_idx,comp_idx,node_idx,dev_idx,NUM_OF_NODAL_DEV, pow_idx
    TYPE(VARYING_STRING) :: STRING_DATA
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    REAL(DP) :: VALUE_BUFFER(5)
    REAL(DP) :: RMS_PERCENT, RMS_ABSOLUTE, RMS_RELATIVE, INTEGRAL_NUM, INTEGRAL_ANA
    
    CALL ENTERS("ANALYTIC_ANALYSIS_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE(FIELD,FILE_ID,(/1/),ERR,ERROR,*999)
!TODO: for Bicubic Hermit and regular mesh only      
      IF (FIELD%VARIABLES(1)%COMPONENTS(1)%DOMAIN%TOPOLOGY%NODES%NODES(1)%NUMBER_OF_DERIVATIVES==4)THEN
        CALL ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE(FIELD,FILE_ID,(/2,3/),ERR,ERROR,*999)
        CALL ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE(FIELD,FILE_ID,(/4/),ERR,ERROR,*999) 
      ENDIF
      DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
        IF(var_idx==1) THEN
          STRING_DATA="Dependent variable integral error"
          CALL WRITE_STRING(FILE_ID,STRING_DATA,ERR,ERROR,*999)
          STRING_DATA="Node #              Numerical      Analytic      % error      Absolute error Relative error"
          CALL WRITE_STRING(FILE_ID,STRING_DATA,ERR,ERROR,*999)
          DO pow_idx=1,2 
            CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,var_idx,pow_idx,(/1/),VALUE_BUFFER(1),ERR,ERROR,*999)
            CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,var_idx,pow_idx,(/1/),VALUE_BUFFER(2),ERR,ERROR,*999)
            CALL ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET(FIELD,var_idx,pow_idx,(/1/),VALUE_BUFFER(3),ERR,ERROR,*999)
            CALL ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET(FIELD,var_idx,pow_idx,(/1/),VALUE_BUFFER(4),ERR,ERROR,*999)
            CALL ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET(FIELD,var_idx,pow_idx,(/1/),VALUE_BUFFER(5),ERR,ERROR,*999)
            SELECT CASE(pow_idx)
            CASE(1)
              CALL WRITE_STRING_VECTOR(FILE_ID,1,1,5,5,5,VALUE_BUFFER,'("  Intgl         ",5(X,D13.4))','(20X,5(X,D13.4))', &
                & ERR,ERROR,*999)
            CASE(2)
              CALL WRITE_STRING_VECTOR(FILE_ID,1,1,5,5,5,VALUE_BUFFER,'("  Int^2         ",5(X,D13.4))','(20X,5(X,D13.4))', &
                & ERR,ERROR,*999)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid power value!",ERR,ERROR,*999)     
            END SELECT
          ENDDO ! pow_idx
          
          STRING_DATA="Node #              Numerical           NID"
          CALL WRITE_STRING(FILE_ID,STRING_DATA,ERR,ERROR,*999)
          DO pow_idx=1,2 
            CALL ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET(FIELD,var_idx,pow_idx,(/1/),VALUE_BUFFER(1),ERR,ERROR,*999)
            CALL ANALYTIC_ANALYSIS_NIDS_ERROR_GET(FIELD,var_idx,pow_idx,(/1/),VALUE_BUFFER(2),ERR,ERROR,*999)
            SELECT CASE(pow_idx)
            CASE(1)
              CALL WRITE_STRING_VECTOR(FILE_ID,1,1,2,2,2,VALUE_BUFFER(1:2),'("  Diff.         ",2(X,D13.4))','(20X,2(X,D13.4))', &
                & ERR,ERROR,*999)
            CASE(2)
              CALL WRITE_STRING_VECTOR(FILE_ID,1,1,2,2,2,VALUE_BUFFER(1:2),'("  Dif^2         ",2(X,D13.4))','(20X,2(X,D13.4))', &
                & ERR,ERROR,*999)
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid power value!",ERR,ERROR,*999)     
            END SELECT
          ENDDO ! pow_idx
        ENDIF
      ENDDO !var_idx
    ELSE
       CALL FLAG_ERROR("The field is not associated!",ERR,ERROR,*999)     
    ENDIF
          
    CALL EXITS("ANALYTIC_ANALYSIS_CALCULATE")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_CALCULATE",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_CALCULATE")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_CALCULATE
  
  !
  !================================================================================================================================
  !  

  !>Calculate the analytic analysis data. 
  SUBROUTINE ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE(FIELD,FILE_ID,DERIVATIVE_NUMBERS,ERR,ERROR,*)
  
    !Argument variables 
    TYPE(FIELD_TYPE), INTENT(IN), POINTER :: FIELD  
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: var_idx,comp_idx,node_idx,dev_idx,NUM_OF_NODAL_DEV, pow_idx
    TYPE(VARYING_STRING) :: STRING_DATA
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    REAL(DP) :: VALUE_BUFFER(5)
    REAL(DP) :: RMS_PERCENT, RMS_ABSOLUTE, RMS_RELATIVE, INTEGRAL_NUM, INTEGRAL_ANA
    
    CALL ENTERS("ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(FIELD)) THEN
      DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
        IF(var_idx==1) THEN
          IF(DERIVATIVE_NUMBERS(1)==1) STRING_DATA="Dependent variable"
          IF(DERIVATIVE_NUMBERS(1)==2) STRING_DATA="Arc Length Derivative"
          IF(DERIVATIVE_NUMBERS(1)==4) STRING_DATA="Arc Length Cross Derivative"
          RMS_PERCENT=0.0_DP
          RMS_ABSOLUTE=0.0_DP
          RMS_RELATIVE=0.0_DP
          INTEGRAL_NUM=0.0_DP
          INTEGRAL_ANA=0.0_DP
          CALL WRITE_STRING(FILE_ID,STRING_DATA,ERR,ERROR,*999)
          STRING_DATA="Node #              Numerical      Analytic      % error      Absolute error Relative error"
          CALL WRITE_STRING(FILE_ID,STRING_DATA,ERR,ERROR,*999)
          DO comp_idx=1,FIELD%VARIABLES(var_idx)%NUMBER_OF_COMPONENTS
            DOMAIN_NODES=>FIELD%VARIABLES(var_idx)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
              NUM_OF_NODAL_DEV=SIZE(DERIVATIVE_NUMBERS)
              DO dev_idx=1,NUM_OF_NODAL_DEV
                CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,var_idx,VALUE_BUFFER(1),ERR,ERROR, &
                  & *999)
                CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,var_idx,VALUE_BUFFER(2),ERR,ERROR,*999)
                CALL ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,var_idx,VALUE_BUFFER(3),ERR,ERROR,*999)
                CALL ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,var_idx,VALUE_BUFFER(4),ERR,ERROR,*999)
                CALL ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,var_idx,VALUE_BUFFER(5),ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(FILE_ID,1,1,5,5,5,VALUE_BUFFER, &
                  & CHAR('("     '//NUMBER_TO_VSTRING(node_idx,"*",ERR,ERROR)//'",5(X,D13.4))'),'(20X,5(X,D13.4))', &
                  & ERR,ERROR,*999)
              ENDDO !dev_idx
            ENDDO !node_idx
          ENDDO !comp_idx
        
          CALL ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET(FIELD,var_idx,DERIVATIVE_NUMBERS,RMS_PERCENT,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(FILE_ID,"RMS error (Percent) = ",RMS_PERCENT,ERR,ERROR,*999)
          CALL ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET(FIELD,var_idx,DERIVATIVE_NUMBERS,RMS_ABSOLUTE,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(FILE_ID,"RMS error (Absolute) = ",RMS_ABSOLUTE,ERR,ERROR,*999)
          CALL ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET(FIELD,var_idx,DERIVATIVE_NUMBERS,RMS_RELATIVE,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(FILE_ID,"RMS error (Relative) = ",RMS_RELATIVE,ERR,ERROR,*999)
        ENDIF
      ENDDO !var_idx
    ELSE
       CALL FLAG_ERROR("The field is not associated!",ERR,ERROR,*999)     
    ENDIF
          
    CALL EXITS("ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODAL_RMS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Export analytic information \see{ANALYTIC_ANALYSIS_ROUTINES::ANALYTIC_ANALYSIS_EXPORT}.                 
  SUBROUTINE ANALYTIC_ANALYSIS_EXPORT(FIELD,FILE_NAME, METHOD, ERR,ERROR,*)
    !Argument variables       
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field object
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<file name
    TYPE(VARYING_STRING), INTENT(IN):: METHOD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: FILE_STATUS 
    INTEGER(INTG) :: FILE_ID

    CALL ENTERS("ANALYTIC_ANALYSIS_EXPORT", ERR,ERROR,*999)    

    IF(METHOD=="FORTRAN") THEN
       FILE_STATUS="REPLACE"
       FILE_ID=1245+FIELD%GLOBAL_NUMBER
       CALL ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
       CALL ANALYTIC_ANALYSIS_CALCULATE(FIELD,FILE_ID,ERR,ERROR,*999)
    ELSE IF(METHOD=="MPIIO") THEN
       CALL FLAG_ERROR("MPI IO has not been implemented yet!",ERR,ERROR,*999)
    ELSE 
       CALL FLAG_ERROR("Unknown method!",ERR,ERROR,*999)   
    ENDIF   
    
    CALL EXITS("ANALYTIC_ANALYSIS_EXPORT")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_EXPORT",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_EXPORT")
    RETURN 1  
  END SUBROUTINE ANALYTIC_ANALYSIS_EXPORT
  
  !
  !================================================================================================================================
  !

  !>Open a file using Fortran. TODO should we use method in FIELD_IO??     
  SUBROUTINE ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<the name of file.
    TYPE(VARYING_STRING), INTENT(IN) :: FILE_STATUS !<status for opening a file
    INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
        
    CALL ENTERS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN",ERR,ERROR,*999)       

    !CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"OPEN FILE",ERR,ERROR,*999)

    OPEN(UNIT=FILE_ID, FILE=CHAR(FILE_NAME), STATUS=CHAR(FILE_STATUS), FORM="FORMATTED", ERR=999)   
        
    
    CALL EXITS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_FORTRAN_FILE_OPEN
  
  !
  !================================================================================================================================
  !

  !>Get integral absolute error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the integral absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_ABSOLUTE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ABS(ANALYTIC_VALUE-NUMERICAL_VALUE)
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
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) !<derivative number   
    REAL(DP), INTENT(OUT) :: VALUE !<the integral analytic value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: ANALYTIC_VALUE,h,k
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx,dev_idx,NUMBER_OF_SURROUNDING_ELEMENTS
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,SIZE(DERIVATIVE_NUMBERS)
	        CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,VARIABLE_NUMBER,ANALYTIC_VALUE, &
	          & ERR,ERROR,*999)
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
              VALUE=VALUE+ANALYTIC_VALUE*NUMBER_OF_SURROUNDING_ELEMENTS*h*k/4
            CASE(2)
              VALUE=VALUE+ANALYTIC_VALUE**2*NUMBER_OF_SURROUNDING_ELEMENTS*h*k/4
            CASE DEFAULT
              CALL FLAG_ERROR("Not valid power number",ERR,ERROR,*999)
            END SELECT
	      ENDDO !dev_idx
        ENDDO !node_idx
      ENDDO !comp_idx
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
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the integral numerical value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE,h,k
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx,dev_idx,NUMBER_OF_SURROUNDING_ELEMENTS
    TYPE(GENERATED_MESH_REGULAR_TYPE), POINTER :: REGULAR_MESH
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,SIZE(DERIVATIVE_NUMBERS)
            CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,VARIABLE_NUMBER,NUMERICAL_VALUE, &
              & ERR,ERROR,*999)
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
              VALUE=VALUE+NUMERICAL_VALUE*NUMBER_OF_SURROUNDING_ELEMENTS*h*k/4
            CASE(2)
              VALUE=VALUE+NUMERICAL_VALUE**2*NUMBER_OF_SURROUNDING_ELEMENTS*h*k/4
            CASE DEFAULT
              CALL FLAG_ERROR("Not valid power number",ERR,ERROR,*999)
            END SELECT
          ENDDO !dev_idx
        ENDDO !node_idx
      ENDDO !comp_idx
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
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the integral percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ABS(ANALYTIC_VALUE)>ZERO_TOLERANCE) THEN
        VALUE=(ANALYTIC_VALUE-NUMERICAL_VALUE)/ANALYTIC_VALUE*100
      ELSE
        VALUE=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_PERCENT_ERROR_GET
  
   !
  !================================================================================================================================
  !

  !>Get integral relative error value for the field.
  SUBROUTINE ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: POWER !<power
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) !<derivative number
    REAL(DP), INTENT(OUT) :: VALUE !<the integral relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_INTEGRAL_RELATIVE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_INTEGRAL_ANALYTIC_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ABS(ANALYTIC_VALUE+1)>ZERO_TOLERANCE) THEN
        VALUE=(NUMERICAL_VALUE-ANALYTIC_VALUE)/(1+ANALYTIC_VALUE)
      ELSE
        VALUE=0.0_DP
      ENDIF
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
    REAL(DP), INTENT(OUT) :: VALUE !<the integral absolute error
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
            CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,VARIABLE_NUMBER,NUMERICAL_VALUE, &
              & ERR,ERROR,*999)
            CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,VARIABLE_NUMBER,ANALYTIC_VALUE, &
              & ERR,ERROR,*999)
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
    REAL(DP), INTENT(OUT) :: VALUE !<the integral absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, NUMERICAL_ERROR
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NIDS_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      CALL ANALYTIC_ANALYSIS_NIDS_NUMERICAL_ERROR_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,NUMERICAL_ERROR,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_INTEGRAL_NUMERICAL_VALUE_GET(FIELD,VARIABLE_NUMBER,POWER,DERIVATIVE_NUMBERS,NUMERICAL_VALUE,ERR,ERROR,*999)
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
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER,VALUE, &
    & ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, &
        & NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, &
        & ANALYTIC_VALUE,ERR,ERROR,*999)
      VALUE=ABS(NUMERICAL_VALUE-ANALYTIC_VALUE)
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get Analytic value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER,VALUE, &
    & ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the analytic value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:)
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: NUMBER_OF_DERIVATIVES
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET(FIELD,FIELD_ANALYTIC_SET_TYPE,ANALYTIC_PARAMETERS,ERR,ERROR,*999)
      IF(ASSOCIATED(ANALYTIC_PARAMETERS)) THEN      
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES
        IF(ASSOCIATED(DOMAIN_NODES)) THEN
          NUMBER_OF_DERIVATIVES=DOMAIN_NODES%NODES(node_number)%NUMBER_OF_DERIVATIVES
          VALUE=ANALYTIC_PARAMETERS((VARIABLE_NUMBER-1)*DOMAIN_NODES%NUMBER_OF_NODES*NUMBER_OF_DERIVATIVES+(NODE_NUMBER-1)* &
            & NUMBER_OF_DERIVATIVES+DERIVATIVE_NUMBER)  
        ELSE
          CALL FLAG_ERROR("Domain nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("ANALYTIC_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF 
    
    IF(ASSOCIATED(ANALYTIC_PARAMETERS)) THEN
      NULLIFY(ANALYTIC_PARAMETERS)
    ELSE
      CALL FLAG_ERROR("ANALYTIC_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET
  
  !
  !================================================================================================================================
  !

  !>Get Numerical value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, &
    & VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the numerical value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), POINTER :: NUMERICAL_PARAMETERS(:)
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: NUMBER_OF_DERIVATIVES
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_PARAMETER_SET_GET(FIELD,FIELD_VALUES_SET_TYPE,NUMERICAL_PARAMETERS,ERR,ERROR,*999)
            
      IF(ASSOCIATED(NUMERICAL_PARAMETERS)) THEN      
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES
        IF(ASSOCIATED(DOMAIN_NODES)) THEN
          NUMBER_OF_DERIVATIVES=DOMAIN_NODES%NODES(node_number)%NUMBER_OF_DERIVATIVES
          VALUE=NUMERICAL_PARAMETERS((VARIABLE_NUMBER-1)*DOMAIN_NODES%NUMBER_OF_NODES*NUMBER_OF_DERIVATIVES+(NODE_NUMBER-1) &
            & *NUMBER_OF_DERIVATIVES+DERIVATIVE_NUMBER)  
        ELSE
          CALL FLAG_ERROR("Domain nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("NUMERICAL_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF 
    
    IF(ASSOCIATED(NUMERICAL_PARAMETERS)) THEN
      NULLIFY(NUMERICAL_PARAMETERS)
    ELSE
      CALL FLAG_ERROR("NUMERICAL_PARAMETERS is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET
  
  
  !
  !================================================================================================================================
  !

  !>Get percentage error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER,VALUE, &
    & ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, & 
        & NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, & 
        & ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ABS(ANALYTIC_VALUE)>ZERO_TOLERANCE) THEN
        VALUE=(ANALYTIC_VALUE-NUMERICAL_VALUE)/ANALYTIC_VALUE*100
      ELSE
        VALUE=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET
  
  
 
  !
  !================================================================================================================================
  !

  !>Get relative error value for the node
  SUBROUTINE ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, & 
    & VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<derivative number
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<node number
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<component number
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    REAL(DP), INTENT(OUT) :: VALUE !<the relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: NUMERICAL_VALUE, ANALYTIC_VALUE
        
    CALL ENTERS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      CALL ANALYTIC_ANALYSIS_NODE_NUMERICIAL_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, &
        & NUMERICAL_VALUE,ERR,ERROR,*999)
      CALL ANALYTIC_ANALYSIS_NODE_ANALYTIC_VALUE_GET(FIELD,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER, &
        & ANALYTIC_VALUE,ERR,ERROR,*999)
      IF(ABS(ANALYTIC_VALUE+1)>ZERO_TOLERANCE) THEN
        VALUE=ABS((NUMERICAL_VALUE-ANALYTIC_VALUE)/(1+ANALYTIC_VALUE))
      ELSE
        VALUE=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get rms percentage error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET(FIELD,VARIABLE_NUMBER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:)
    REAL(DP), INTENT(OUT) :: VALUE !<the rms percentage error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: PERCENT_ERROR_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx,dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,SIZE(DERIVATIVE_NUMBERS)
            CALL ANALYTIC_ANALYSIS_NODE_PERCENT_ERROR_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,VARIABLE_NUMBER,PERCENT_ERROR_VALUE, & 
              & ERR,ERROR,*999)
            VALUE=VALUE+PERCENT_ERROR_VALUE**2/(DOMAIN_NODES%NUMBER_OF_NODES*DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES)
          ENDDO !dev_idx
        ENDDO !node_idx
        VALUE=SQRT(VALUE)
      ENDDO !comp_idx
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_PERCENT_ERROR_GET
  
  
  !
  !================================================================================================================================
  !

  !>Get rms absolute error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET(FIELD,VARIABLE_NUMBER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:) 
    REAL(DP), INTENT(OUT) :: VALUE !<the rms absolute error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: ABSOLUTE_ERROR_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx, dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,SIZE(DERIVATIVE_NUMBERS)
            CALL ANALYTIC_ANALYSIS_NODE_ABSOLUTE_ERROR_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,VARIABLE_NUMBER,ABSOLUTE_ERROR_VALUE, &
              & ERR,ERROR,*999)
            VALUE=VALUE+ABSOLUTE_ERROR_VALUE**2/(DOMAIN_NODES%NUMBER_OF_NODES*DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES)
          ENDDO !dev_idx
        ENDDO !node_idx
        VALUE=SQRT(VALUE)
      ENDDO !comp_idx
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_ABSOLUTE_ERROR_GET
  
  !
  !================================================================================================================================
  !

  !>Get rms relative error value for the field
  SUBROUTINE ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET(FIELD,VARIABLE_NUMBER,DERIVATIVE_NUMBERS,VALUE,ERR,ERROR,*)
  
    !Argument variables   
    TYPE(FIELD_TYPE), POINTER :: FIELD !<the field.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<variable number
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBERS(:)
    REAL(DP), INTENT(OUT) :: VALUE !<the rms relative error
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: RELATIVE_ERROR_VALUE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    INTEGER(INTG) :: comp_idx,node_idx, dev_idx
        
    CALL ENTERS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET",ERR,ERROR,*999)       

    IF(ASSOCIATED(FIELD)) THEN
      VALUE=0.0_DP
      DO comp_idx=1,FIELD%VARIABLES(VARIABLE_NUMBER)%NUMBER_OF_COMPONENTS
        DOMAIN_NODES=>FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
          DO dev_idx=1,SIZE(DERIVATIVE_NUMBERS)
            CALL ANALYTIC_ANALYSIS_NODE_RELATIVE_ERROR_GET(FIELD,DERIVATIVE_NUMBERS(dev_idx),node_idx,comp_idx,VARIABLE_NUMBER,RELATIVE_ERROR_VALUE, &
              & ERR,ERROR,*999)
            VALUE=VALUE+RELATIVE_ERROR_VALUE**2/(DOMAIN_NODES%NUMBER_OF_NODES*DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES)
          ENDDO !dev_idx
        ENDDO !node_idx
        VALUE=SQRT(VALUE)
      ENDDO !comp_idx
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF 
    
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET")
    RETURN
999 CALL ERRORS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET",ERR,ERROR)
    CALL EXITS("ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET")
    RETURN 1
  END SUBROUTINE ANALYTIC_ANALYSIS_RMS_RELATIVE_ERROR_GET

  !
  !================================================================================================================================
  !  
  
 
END MODULE ANALYTIC_ANALYSIS_ROUTINES
