!> \file
!> \author Chris Bradley
!> \brief This module handles all bioelectric routines.
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

!> This module handles all bioelectric class routines.
MODULE BIOELECTRIC_ROUTINES

  USE BASE_ROUTINES
  USE BIODOMAIN_EQUATION_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE MONODOMAIN_EQUATIONS_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC BIOELECTRIC_CONTROL_LOOP_POST_LOOP
  
  PUBLIC BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_GET,BIOELECTRIC_PROBLEM_CLASS_TYPE_GET

  PUBLIC BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_SET

  PUBLIC BIOELECTRIC_FINITE_ELEMENT_CALCULATE

  PUBLIC BIOELECTRIC_EQUATIONS_SET_SETUP

  PUBLIC BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC BIOELECTRIC_PROBLEM_CLASS_TYPE_SET

  PUBLIC BIOELECTRIC_PROBLEM_SETUP
  
  PUBLIC BIOELECTRIC_PRE_SOLVE,BIOELECTRIC_POST_SOLVE
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop for bioelectric problems, i.e., after each time step for a time loop
  SUBROUTINE BIOELECTRIC_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BIOELECTRIC_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          SELECT CASE(PROBLEM%TYPE)
          CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
            CALL BIODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
            !do nothing
          CASE DEFAULT
            LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
              & " is not valid for a bioelectric problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BIOELECTRIC_CONTROL_LOOP_POST_LOOP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_CONTROL_LOOP_POST_LOOP")
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a bioelectric equation set class.
  SUBROUTINE BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_GET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_TYPE !<The equation type
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SUBTYPE !<The equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF (EQUATIONS_SET%CLASS==EQUATIONS_SET_BIOELECTRICS_CLASS) THEN
        EQUATIONS_TYPE=EQUATIONS_SET%TYPE
        EQUATIONS_SUBTYPE=EQUATIONS_SET%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Equations set is not a bioelectric class.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_GET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_GET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a bioelectric equation set class.
  SUBROUTINE BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_TYPE !<The equation type
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SUBTYPE !<The equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_EQUATIONS_SET_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_TYPE)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
          & ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
          & ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a bioelectric class finite element equation set.
  SUBROUTINE BIOELECTRIC_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a bioelectric equations set class.
  SUBROUTINE BIOELECTRIC_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_EQUATIONS_SET_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectric equation set class.
  SUBROUTINE BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Perform pre-solve actions for the bioelectrics problem class.
  SUBROUTINE BIOELECTRIC_PRE_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BIOELECTRIC_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          PROBLEM=>CONTROL_LOOP%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            SELECT CASE(PROBLEM%TYPE)
            CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
              CALL BIODOMAIN_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
              CALL MONODOMAIN_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_PRE_SOLVE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Performs post solve actions for a bioelectrics problem class.
  SUBROUTINE BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BIOELECTRIC_POST_SOLVE",ERR,ERROR,*999)
  
    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          PROBLEM=>CONTROL_LOOP%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            SELECT CASE(PROBLEM%TYPE)
            CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
              !Do nothing???
            CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
              CALL MONODOMAIN_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BIOELECTRIC_POST_SOLVE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_POST_SOLVE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_POST_SOLVE")
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a bioelectric problem class.
  SUBROUTINE BIOELECTRIC_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<The problem type
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<The proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("BIOELECTRIC_PROBLEM_CLASS_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%CLASS==PROBLEM_BIOELECTRICS_CLASS) THEN
        PROBLEM_EQUATION_TYPE=PROBLEM%TYPE
        PROBLEM_SUBTYPE=PROBLEM%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Problem is not bioelectric class.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_PROBLEM_CLASS_TYPE_GET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_PROBLEM_CLASS_TYPE_GET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_PROBLEM_CLASS_TYPE_GET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_PROBLEM_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a bioelectric problem class.
  SUBROUTINE BIOELECTRIC_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE !<The problem type
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_PROBLEM_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_EQUATION_TYPE)
      CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE, &
          & ERR,ERROR,*999)
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        CALL MONODOMAIN_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)        
      CASE DEFAULT
        LOCAL_ERROR="Problem equation type "//TRIM(NUMBER_TO_VSTRING(PROBLEM_EQUATION_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_PROBLEM_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_PROBLEM_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_PROBLEM_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_PROBLEM_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a bioelectric problem class.
  SUBROUTINE BIOELECTRIC_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%TYPE)
      CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        CALL MONODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

END MODULE BIOELECTRIC_ROUTINES

