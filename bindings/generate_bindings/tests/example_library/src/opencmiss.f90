MODULE CMFE

  PRIVATE

  !>An example cmfe_ type.
  TYPE cmfe_ExampleType
  END TYPE cmfe_ExampleType

  INTERFACE cmfe_Example_SomeInterface
    MODULE PROCEDURE cmfe_Example_SomeInterfaceNumber
    MODULE PROCEDURE cmfe_Example_SomeInterfaceObj
  END INTERFACE !cmfe_Example_SomeInterface

  INTERFACE cmfe_Example_CreateStart
    MODULE PROCEDURE cmfe_Example_CreateStartObj
    MODULE PROCEDURE cmfe_Example_CreateStartNumber
  END INTERFACE !cmfe_Example_CreateStart

  INTERFACE cmfe_ArrayRoutine
    MODULE PROCEDURE cmfe_ArrayRoutine0
    MODULE PROCEDURE cmfe_ArrayRoutine1
  END INTERFACE !cmfe_ArrayRoutine

  INTERFACE cmfe_StringRoutine
    MODULE PROCEDURE cmfe_StringRoutineCObj
    MODULE PROCEDURE cmfe_StringRoutineVSObj
    MODULE PROCEDURE cmfe_StringRoutineCNumber
    MODULE PROCEDURE cmfe_StringRoutineVSNumber
  END INTERFACE !cmfe_StringRoutine

  !> \addtogroup OPENCMISS_ExampleEnum OPENCMISS::ExampleEnum
  !> \brief Example of an enum
  !>@{
  INTEGER(INTG), PARAMETER :: CMFE_ENUM_ONE = FIRST_VALUE !<Description of first enum value
  INTEGER(INTG), PARAMETER :: CMFE_ENUM_TWO = SECOND_VALUE !<Description of second enum value
  INTEGER(INTG), PARAMETER :: CMFE_ENUM_THREE = THIRD_VALUE !<Description of third enum value
  !>@}

  INTEGER(INTG), PARAMETER :: UNGROUPED_CONSTANT = 1 !<Description

  INTEGER(INTG), PARAMETER :: NON_PUBLIC_CONSTANT = 1

  PUBLIC cmfe_ExampleType, cmfe_Example_SomeInterface, cmfe_Example_CreateStart, &
    & cmfe_Example_Initialise, cmfe_StringRoutine, cmfe_ArrayRoutine, CMFE_ENUM_ONE, &
    & CMFE_ENUM_TWO, CMFE_ENUM_THREE, UNGROUPED_CONSTANT

CONTAINS

  !>Doxygen comment describing subroutine
  SUBROUTINE cmfe_Example_SomeInterfaceObj(Example, InputString, OutputString, InputArray, &
      & OutputArray, InputArray2D, OutputArray2D, ArrayWithSize, InputReal, OutputReal, Err)
    TYPE(cmfe_ExampleType), INTENT(INOUT) :: Example !<Comment for Example
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
  END SUBROUTINE cmfe_Example_SomeInterfaceObj

  SUBROUTINE cmfe_Example_SomeInterfaceNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_Example_SomeInterfaceNumber

  SUBROUTINE cmfe_StringRoutineCObj(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_StringRoutineCObj

  SUBROUTINE cmfe_StringRoutineVSObj(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_StringRoutineVSObj

  SUBROUTINE cmfe_StringRoutineCNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_StringRoutineCNumber

  SUBROUTINE cmfe_StringRoutineVSNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_StringRoutineVSNumber

  SUBROUTINE cmfe_ArrayRoutine0(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_ArrayRoutine0

  SUBROUTINE cmfe_ArrayRoutine1(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_ArrayRoutine1

  SUBROUTINE cmfe_Example_Initialise(Example, Err)
    TYPE(cmfe_ExampleType), INTENT(OUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_Example_Initialise

  SUBROUTINE cmfe_Example_CreateStartObj(UserNumber, Example, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    TYPE(cmfe_ExampleType), INTENT(INOUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_Example_CreateStartObj

  SUBROUTINE cmfe_Example_CreateStartNumber(UserNumber, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_Example_CreateStartNumber

  SUBROUTINE cmfe_NonPublicRoutine(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE cmfe_NonPublicRoutine

END MODULE CMFE
