MODULE OPENCMISS

  PRIVATE

  !>An example CMISS type.
  TYPE CMISSExampleType
  END TYPE CMISSExampleType

  INTERFACE CMISSExampleSomeInterface
    MODULE PROCEDURE CMISSExampleSomeInterfaceNumber
    MODULE PROCEDURE CMISSExampleSomeInterfaceObj
  END INTERFACE !CMISSExampleSomeInterface

  INTERFACE CMISSExampleCreateStart
    MODULE PROCEDURE CMISSExampleCreateStartObj
    MODULE PROCEDURE CMISSExampleCreateStartNumber
  END INTERFACE !CMISSExampleCreateStart

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
  INTEGER(INTG), PARAMETER :: CMISSEnumOne = FIRST_VALUE !<Description of first enum value
  INTEGER(INTG), PARAMETER :: CMISSEnumTwo = SECOND_VALUE !<Description of second enum value
  INTEGER(INTG), PARAMETER :: CMISSEnumThree = THIRD_VALUE !<Description of third enum value
  !>@}

  INTEGER(INTG), PARAMETER :: UngroupedConstant = 1 !<Description

  INTEGER(INTG), PARAMETER :: NonPublicConstant = 1

  PUBLIC CMISSExampleType, CMISSExampleSomeInterface, CMISSExampleCreateStart, &
    & CMISSExampleTypeInitialise, CMISSStringRoutine, CMISSArrayRoutine, CMISSEnumOne, &
    & CMISSEnumTwo, CMISSEnumThree, UngroupedConstant

CONTAINS

  !>Doxygen comment describing subroutine
  SUBROUTINE CMISSExampleSomeInterfaceObj(Example, InputString, OutputString, InputArray, &
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
  END SUBROUTINE CMISSExampleSomeInterfaceObj

  SUBROUTINE CMISSExampleSomeInterfaceNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExampleSomeInterfaceNumber

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

  SUBROUTINE CMISSExampleTypeInitialise(Example, Err)
    TYPE(CMISSExampleType), INTENT(OUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExampleTypeInitialise

  SUBROUTINE CMISSExampleCreateStartObj(UserNumber, Example, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    TYPE(CMISSExampleType), INTENT(INOUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExampleCreateStartObj

  SUBROUTINE CMISSExampleCreateStartNumber(UserNumber, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSExampleCreateStartNumber

  SUBROUTINE CMISSNonPublicRoutine(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE CMISSNonPublicRoutine

END MODULE OPENCMISS
