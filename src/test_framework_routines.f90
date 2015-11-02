!> \file
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

!>This module handles test framework routines.
MODULE TEST_FRAMEWORK_ROUTINES

  USE CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  
  INTERFACE TEST_FRAMEWORK_ASSERT_EQUALS
    MODULE PROCEDURE TEST_FRAMEWORK_ASSERT_EQUALS_INTG
    MODULE PROCEDURE TEST_FRAMEWORK_ASSERT_EQUALS_DP
    MODULE PROCEDURE TEST_FRAMEWORK_ASSERT_EQUALS_DP1
    MODULE PROCEDURE TEST_FRAMEWORK_ASSERT_EQUALS_DP2
  END INTERFACE !TEST_FRAMEWORK_ASSERT_EQUALS
  
  PUBLIC TEST_FRAMEWORK_ASSERT_EQUALS, TEST_FRAMEWORK_GRADIENT_VALUE_GET

CONTAINS  

  !
  !================================================================================================================================
  !   
  
  !>Check if the actual integer value is as expected.
  SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_INTG(EXPECTED_VALUE,ACTUAL_VALUE,ERR)
  
    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: EXPECTED_VALUE !<expected value
    INTEGER(INTG), INTENT(IN) :: ACTUAL_VALUE !<actual value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("TEST_FRAMEWORK_ASSERT_EQUALS_INTG",ERR,ERROR,*999)
    
    IF(EXPECTED_VALUE/=ACTUAL_VALUE) THEN
       CALL FlagError("Not equal",ERR,ERROR,*999)
    ENDIF
          
    EXITS("TEST_FRAMEWORK_ASSERT_EQUALS_INTG")
    RETURN
999 ERRORSEXITS("TEST_FRAMEWORK_ASSERT_EQUALS_INTG",ERR,ERROR)
    RETURN 
  END SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_INTG
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the actual real(DP) value is as expected.
  SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_DP1(EXPECTED_VALUE,ACTUAL_VALUE,ERR)
  
    !Argument variables 
    REAL(DP), INTENT(IN) :: EXPECTED_VALUE !<expected value
    REAL(DP), INTENT(IN) :: ACTUAL_VALUE !<actual value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("TEST_FRAMEWORK_ASSERT_EQUALS_DP1",ERR,ERROR,*999)
    
    CALL TEST_FRAMEWORK_ASSERT_EQUALS_DP2(EXPECTED_VALUE,ACTUAL_VALUE,ZERO_TOLERANCE,ERR)
          
    EXITS("TEST_FRAMEWORK_ASSERT_EQUALS_DP1")
    RETURN
999 ERRORSEXITS("TEST_FRAMEWORK_ASSERT_EQUALS_DP1",ERR,ERROR)
    RETURN 
  END SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_DP1
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the actual real(DP) value is as expected.
  SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_DP2(EXPECTED_VALUE,ACTUAL_VALUE,TOLERANCE,ERR)
  
    !Argument variables 
    REAL(DP), INTENT(IN) :: EXPECTED_VALUE !<expected value
    REAL(DP), INTENT(IN) :: ACTUAL_VALUE !<actual value
    REAL(DP), INTENT(IN) :: TOLERANCE !<tolerance
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("TEST_FRAMEWORK_ASSERT_EQUALS_DP2",ERR,ERROR,*999)
    
    IF(ABS(EXPECTED_VALUE-ACTUAL_VALUE)>TOLERANCE) THEN
       CALL FlagError("Not Equal",ERR,ERROR,*999)
    ENDIF
          
    EXITS("TEST_FRAMEWORK_ASSERT_EQUALS_DP2")
    RETURN
999 ERRORSEXITS("TEST_FRAMEWORK_ASSERT_EQUALS_DP1",ERR,ERROR)
    RETURN 
  END SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_DP2
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the actual real(DP) values is as expected.
  SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_DP(EXPECTED_VALUE,ACTUAL_VALUE,ERR)
  
    !Argument variables 
    REAL(DP), INTENT(IN) :: EXPECTED_VALUE(:) !<expected value
    REAL(DP), INTENT(IN) :: ACTUAL_VALUE(:) !<actual value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    
    ENTERS("TEST_FRAMEWORK_ASSERT_EQUALS_DP",ERR,ERROR,*999)
    
    IF(SIZE(EXPECTED_VALUE)==SIZE(ACTUAL_VALUE)) THEN
      DO i=1,SIZE(EXPECTED_VALUE)
        CALL TEST_FRAMEWORK_ASSERT_EQUALS(EXPECTED_VALUE(i),ACTUAL_VALUE(i),ERR)
      ENDDO !i
    ELSE
       CALL FlagError("Not Equal",ERR,ERROR,*999)
    ENDIF
          
    EXITS("TEST_FRAMEWORK_ASSERT_EQUALS_DP")
    RETURN
999 ERRORSEXITS("TEST_FRAMEWORK_ASSERT_EQUALS_DP",ERR,ERROR)
    RETURN 
  END SUBROUTINE TEST_FRAMEWORK_ASSERT_EQUALS_DP 
  
  !
  !================================================================================================================================
  !   
  
  !>Get the gradient value of two array.
  SUBROUTINE TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,GRADIENT_VALUE)
  
    !Argument variables 
    REAL(DP), INTENT(IN) :: X_VALUES(:) !<expected value
    REAL(DP), INTENT(IN) :: Y_VALUES(:) !<actual value
    REAL(DP), INTENT(OUT) :: GRADIENT_VALUE !<error message
    INTEGER(INTG) :: ERR !<The error code
    TYPE(VARYING_STRING) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i, interval_size
    
    ENTERS("TEST_FRAMEWORK_GRADIENT_VALUE_GET",ERR,ERROR,*999)

    IF(SIZE(X_VALUES)==SIZE(Y_VALUES)) THEN
      GRADIENT_VALUE=0.0_DP
      interval_size=SIZE(X_VALUES)-1
      DO i=1,SIZE(X_VALUES)
        IF(i/=1) THEN
          GRADIENT_VALUE=GRADIENT_VALUE+(Y_VALUES(i)-Y_VALUES(i-1))/(X_VALUES(i)-X_VALUES(i-1))/interval_size
        ENDIF
      ENDDO !i
    ELSE
       CALL FlagError('Size of x and Size of y do not match',ERR,ERROR,*999)     
    ENDIF

    EXITS("TEST_FRAMEWORK_GRADIENT_VALUE_GET")
    RETURN
999 ERRORSEXITS("TEST_FRAMEWORK_GRADIENT_VALUE_GET",ERR,ERROR)
    RETURN 
  END SUBROUTINE TEST_FRAMEWORK_GRADIENT_VALUE_GET 
 
END MODULE TEST_FRAMEWORK_ROUTINES
