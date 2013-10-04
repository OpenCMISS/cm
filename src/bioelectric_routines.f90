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

  PUBLIC Bioelectric_EquationsSetSpecificationSet

  PUBLIC BIOELECTRIC_FINITE_ELEMENT_CALCULATE

  PUBLIC BIOELECTRIC_EQUATIONS_SET_SETUP

  PUBLIC BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC Bioelectric_ProblemSpecificationSet

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
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
            CALL FlagError("Problem specification must have at least two entries for a bioelectric problem.",err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
            CALL BIODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
            !do nothing
          CASE DEFAULT
            LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
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

  !>Sets the problem specification for a bioelectric equation set class.
  SUBROUTINE Bioelectric_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("Bioelectric_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a bioelectric class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(specification(2))
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BiodomainEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BiodomainEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE DEFAULT
        localError="Equations set equation type "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for a bioelectric equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    CALL Exits("Bioelectric_EquationsSetSpecificationSet")
    RETURN
999 CALL Errors("Bioelectric_EquationsSetSpecificationSet",err,error)
    CALL Exits("Bioelectric_EquationsSetSpecificationSet")
    RETURN 1
  END SUBROUTINE Bioelectric_EquationsSetSpecificationSet

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
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a bioelectric type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
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
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a bioelectric type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
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
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a bioelectric type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
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
            IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
              CALL FlagError("Problem specification must have at least two entries for a bioelectric problem.",err,error,*999)
            END IF
            SELECT CASE(PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
              CALL BIODOMAIN_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
              CALL MONODOMAIN_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
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
            IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
              CALL FlagError("Problem specification must have at least two entries for a bioelectric problem.",err,error,*999)
            END IF
            SELECT CASE(PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE,PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
              !Do nothing???
            CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
              CALL MONODOMAIN_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
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

  !>Sets the problem specification for a bioelectric problem class.
  SUBROUTINE Bioelectric_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType

    CALL Enters("Bioelectric_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)>=2) THEN
        problemType=problemSpecification(2)
        SELECT CASE(problemType)
        CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
          CALL BiodomainEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
          CALL Monodomain_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE DEFAULT
          localError="Problem type "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for a bioelectric problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Bioelectric problem specification must have a type set.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    CALL Exits("Bioelectric_ProblemSpecificationSet")
    RETURN
999 CALL Errors("Bioelectric_ProblemSpecificationSet",err,error)
    CALL Exits("Bioelectric_ProblemSpecificationSet")
    RETURN 1
  END SUBROUTINE Bioelectric_ProblemSpecificationSet

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
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a bioelectric problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
        CALL BIODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        CALL MONODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
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

