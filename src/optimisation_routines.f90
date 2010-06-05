!> This module handles all optimisation problem class routines.
!> reference to classical_field_routines.f90
MODULE OPTIMISATION_ROUTINES

  USE BASE_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE OPTIMISATION_KALMAN_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  
  PUBLIC OPTIMISATION_PROBLEM_SETUP, OPTIMISATION_PROBLEM_CLASS_TYPE_SET,OPTIMISATION_KALMAN_PROBLEM_SETUP
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a optimisation problem class.
  !>reference to CLASSICAL_FIELD_PROBLEM_CLASS_TYPE_SET
  SUBROUTINE OPTIMISATION_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_TYPE_INT,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: PROBLEM_TYPE_INT !<The problem type (postfix 'int',avoid name conflicts with TYPE PROBLEM_TYPE)
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("OPTIMISATION_PROBLEM_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_TYPE_INT)
       CASE(PROBLEM_OPTIMISATION_KALMAN_TYPE)
        CALL OPTIMISATION_KALMAN_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem equation type "//TRIM(NUMBER_TO_VSTRING(PROBLEM_TYPE_INT,"*",ERR,ERROR))// &
          & " is not valid for a optimisation problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("OPTIMISATION_PROBLEM_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("OPTIMISATION_PROBLEM_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("OPTIMISATION_PROBLEM_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE OPTIMISATION_PROBLEM_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a optimisation problem class.
  SUBROUTINE OPTIMISATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("OPTIMISATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%TYPE)
      CASE(PROBLEM_OPTIMISATION_KALMAN_TYPE)
        CALL OPTIMISATION_KALMAN_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("OPTIMISATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("OPTIMISATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("OPTIMISATION_PROBLEM_SETUP")
    RETURN 1

  END SUBROUTINE OPTIMISATION_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

END MODULE OPTIMISATION_ROUTINES