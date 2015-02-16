MODULE OPENCMISS

  PRIVATE

  !>An example CMISS type.
  TYPE CMISSExampleType
  END TYPE CMISSExampleType

  INTERFACE CMISSExample_SomeInterface
    MODULE PROCEDURE CMISSExample_SomeInterfaceNumber
    MODULE PROCEDURE CMISSExample_SomeInterfaceObj
  END INTERFACE !CMISSExample_SomeInterface

  INTERFACE CMISSExample_CreateStart
    MODULE PROCEDURE CMISSExample_CreateStartObj
    MODULE PROCEDURE CMISSExample_CreateStartNumber
  END INTERFACE !CMISSExample_CreateStart

  INTERFACE CMISSArrayRoutine
    MODULE PROCEDURE CMISSArrayRoutine0
    MODULE PROCEDURE CMISSArrayRoutine1
  END INTERFACE !CMISSArrayRoutine

  INTERFACE CMISSStringRoutine
    MODULE PROCEDURE CMISSStringRoutineCObj
    MODULE PROCEDURE CMISSStringRoutineVSObj
    MODULE PROCEDURE CMISSStringRoutineCNumber
    MODULE PROCEDURE CMISSStringRoutineVSNumber
  END INTERFACE !CMISSStringRoutine

  !> \addtogroup OPENCMISS_ExampleEnum OPENCMISS::ExampleEnum
  !> \brief Example of an enum
  !>@{
  INTEGER(INTG), PARAMETER :: CMISS_ENUM_ONE = FIRST_VALUE !<Description of first enum value
  INTEGER(INTG), PARAMETER :: CMISS_ENUM_TWO = SECOND_VALUE !<Description of second enum value
  INTEGER(INTG), PARAMETER :: CMISS_ENUM_THREE = THIRD_VALUE !<Description of third enum value
  !>@}

  INTEGER(INTG), PARAMETER :: UNGROUPED_CONSTANT = 1 !<Description

  INTEGER(INTG), PARAMETER :: NON_PUBLIC_CONSTANT = 1

  PUBLIC CMISSExampleType, CMISSExample_SomeInterface, CMISSExample_CreateStart, &
    & CMISSExample_Initialise, CMISSStringRoutine, CMISSArrayRoutine, CMISS_ENUM_ONE, &
    & CMISS_ENUM_TWO, CMISS_ENUM_THREE, UNGROUPED_CONSTANT

CONTAINS

  !>Doxygen comment describing subroutine
  SUBROUTINE CMISSExample_SomeInterfaceObj(Example, InputString, OutputString, InputArray, &
      & OutputArray, InputArray2D, OutputArray2D, ArrayWithSize, InputReal, OutputReal, Err)
    TYPE(CMISSExampleType), INTENT(INOUT) :: Example !<Comment for Example
    CHARACTER(LEN=*), INTENT(IN) :: InputString !<Comment for InputString
    CHARACTER(LEN=*), INTENT(OUT) :: OutputString !<Comment for OutputString
    INTEGER(INTG), INTENT(IN) :: InputArray(:) !<Comment for InputArray
    INTEGER(INTG), INTENT(OUT) :: OutputArray(:) !<Comment for OutputArray
    INTEGER(INTG), INTENT(IN) :: InputArray2D(:,:) !<Comment for InputArray2D
    INTEGER(INTG), INTENT(OUT) :: OutputArray2D(:,:) !<Comment for OutputArray2D
    INTEGER(INTG), INTENT(IN) :: ArrayWithSize(2) !<Comment for ArrayWithSize
    REAL(DP), INTENT(IN) :: InputReal !<Comment for InputReal
    REAL(DP), INTENT(OUT) :: OutputReal !<Comment for OutputReal
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExample_SomeInterfaceObj

  SUBROUTINE CMISSExample_SomeInterfaceNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExample_SomeInterfaceNumber

  SUBROUTINE CMISSStringRoutineCObj(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSStringRoutineCObj

  SUBROUTINE CMISSStringRoutineVSObj(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSStringRoutineVSObj

  SUBROUTINE CMISSStringRoutineCNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSStringRoutineCNumber

  SUBROUTINE CMISSStringRoutineVSNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSStringRoutineVSNumber

  SUBROUTINE CMISSArrayRoutine0(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSArrayRoutine0

  SUBROUTINE CMISSArrayRoutine1(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSArrayRoutine1

  SUBROUTINE CMISSExample_Initialise(Example, Err)
    TYPE(CMISSExampleType), INTENT(OUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExample_Initialise

  SUBROUTINE CMISSExample_CreateStartObj(UserNumber, Example, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    TYPE(CMISSExampleType), INTENT(INOUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExample_CreateStartObj

  SUBROUTINE CMISSExample_CreateStartNumber(UserNumber, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExample_CreateStartNumber

  SUBROUTINE CMISSNonPublicRoutine(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSNonPublicRoutine

END MODULE OPENCMISS
