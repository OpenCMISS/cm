!> \file
!> $Id: opencmiss_c.f90 542 2009-06-03 17:16:22Z chrispbradley $
!> \author Chris Bradley
!> \brief The top level OpenCMISS module for C bindings.
!>
!> \mainpage OpenCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
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
!>
!> The top level OpenCMISS module for C. This module is the buffer module between the OpenCMISS library and user C code.
MODULE OPENCMISS_C

  USE ISO_C_BINDING
  USE ISO_VARYING_STRING
  USE OPENCMISS
   
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  INTEGER(C_INT), PARAMETER :: CMISSNoError = 0
  INTEGER(C_INT), PARAMETER :: CMISSPointerIsNULL = -1
  INTEGER(C_INT), PARAMETER :: CMISSPointerNotNULL = -2
  INTEGER(C_INT), PARAMETER :: CMISSCouldNotAllocatePointer = -3
  INTEGER(C_INT), PARAMETER :: CMISSErrorConvertingPointer = -4
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISSBasisTypeFinaliseC, CMISSBasisTypeInitialiseC

  PUBLIC CMISSBoundaryConditionsTypeFinaliseC, CMISSBoundaryConditionsTypeInitialiseC

  PUBLIC CMISSControlLoopTypeFinaliseC, CMISSControlLoopTypeInitialiseC

  PUBLIC CMISSCoordinateSystemTypeFinaliseC,CMISSCoordinateSystemTypeInitialiseC

  PUBLIC CMISSDecompositionTypeFinaliseC, CMISSDecompositionTypeInitialiseC

  PUBLIC CMISSEquationsTypeFinaliseC, CMISSEquationsTypeInitialiseC

  PUBLIC CMISSEquationsSetTypeFinaliseC, CMISSEquationsSetTypeInitialiseC

  PUBLIC CMISSRegionTypeFinaliseC, CMISSRegionTypeInitialiseC
  
  PUBLIC CMISSFinaliseC,CMISSInitialiseCNum,CMISSInitialiseCPtr

CONTAINS

  !
  !================================================================================================================================
  !

  !>Copys/converts a C string (array of characters) to a Fortran String (length of characters)
  SUBROUTINE CMISSC2FString(Cstring,Fstring)
    !Argument variables
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Cstring(:)
    CHARACTER(LEN=*), INTENT(OUT) :: Fstring
    !Local variables
    INTEGER(C_INT) :: i,LENGTH

    IF(LEN(Fstring)>=SIZE(Cstring,1)-1) THEN
      LENGTH=SIZE(Cstring,1)-1
    ELSE
      LENGTH=LEN(Fstring)
    ENDIF
    Fstring=""
    DO i=1,LENGTH
      IF(Cstring(i)==C_NULL_CHAR) THEN
        EXIT
      ELSE
        Fstring(i:i)=Cstring(i)
      ENDIF
    ENDDO !i
    
    RETURN
    
  END SUBROUTINE CMISSC2FSTRING
   
  !
  !================================================================================================================================
  !

  !>Copys/converts a  Fortran String (length of characters) to a C string (array of characters)
  SUBROUTINE CMISSF2CString(Fstring,Cstring)
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: Fstring
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Cstring(:)
    !Local variables
    INTEGER(C_INT) :: i,LENGTH

    IF(SIZE(Cstring,1)>LEN_TRIM(Fstring)) THEN
      LENGTH=LEN_TRIM(Fstring)
    ELSE
      LENGTH=SIZE(Cstring,1)-1
    ENDIF
    DO i=1,LENGTH     
      Cstring(i)=Fstring(i:i)
    ENDDO !i
    !Null terminate the string
    Cstring(LENGTH+1)=C_NULL_CHAR
    
    RETURN
    
  END SUBROUTINE CMISSF2CSTRING
   
  !
  !================================================================================================================================
  !

  !>Finalises a CMISSBasisType object for C.
  FUNCTION CMISSBasisTypeFinaliseC(BasisTypePtr) BIND(C, NAME= "CMISSBasisTypeFinalise")

    !Argument Variables
    TYPE(C_PTR), INTENT(INOUT) :: BasisTypePtr !<C pointer to CMISSBasisType object to finalise.
    !Function Variable
    INTEGER(C_INT) :: CMISSBasisTypeFinaliseC !<Error code.
    !Local Variables
    TYPE(CMISSBasisType), POINTER :: BasisType

    CMISSBasisTypeFinaliseC = CMISSNoError

    IF(C_ASSOCIATED(BasisTypePtr)) THEN
      CALL C_F_POINTER(BasisTypePtr,BasisType)
      IF(ASSOCIATED(BasisType)) THEN
        CALL CMISSBasisTypeFinalise(BasisType,CMISSBasisTypeFinaliseC)
        DEALLOCATE(BasisType)
        BasisTypePtr = C_NULL_PTR
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSBasisTypeFinaliseC

  !
  !============================================================================
  !

  !>Initialises a CMISSBasisType object for C.

  FUNCTION CMISSBasisTypeInitialiseC (BasisTypePtr) BIND(C, NAME = "CMISSBasisTypeInitialise")

	!Argument variables
  TYPE(C_PTR), INTENT(INOUT) :: BasisTypePtr !<C pointer to CMISSBasisType object to initialise.
	!Function variable
  INTEGER(C_INT) :: CMISSBasisTypeInitialiseC !<Error code.
	!Local Variables
  INTEGER(C_INT) :: Err
    TYPE (CMISSBasisType), POINTER :: BasisType

    IF(C_ASSOCIATED(BasisTypePtr)) THEN
      CMISSBasisTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (BasisType)
      ALLOCATE (BasisType, STAT = Err)
      IF (Err /= 0) THEN
        CMISSBasisTypeInitialiseC =CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSBasisTypeInitialise (BasisType,CMISSBasisTypeInitialiseC)
        BasisTypePtr = C_LOC (BasisType)
        CMISSBasisTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF
    RETURN

  END FUNCTION CMISSBasisTypeInitialiseC

  !
  !============================================================================
  !

  !>Finalises a CMISSBoundaryConditionsType object for C.

  FUNCTION CMISSBoundaryConditionsTypeFinaliseC(BoundaryConditionsTypePtr) BIND (C, NAME = "CMISSBoundaryConditionsTypeFinalise")

    !Argument Variables
    TYPE (C_PTR), INTENT(INOUT) :: BoundaryConditionsTypePtr !<C pointer to CMISSBoundaryConditionsType object to finalise.
    !Function Variable
    INTEGER(C_INT) :: CMISSBoundaryConditionsTypeFinaliseC !<Error Code.
    !Local Variables
    TYPE(CMISSBoundaryConditionsType), POINTER :: BoundaryConditionsType

    CMISSBoundaryConditionsTypeFinaliseC = CMISSNoError
    IF (C_ASSOCIATED(BoundaryConditionsTypePtr)) THEN
      CALL C_F_POINTER (BoundaryConditionsTypePtr, BoundaryConditionsType)
      IF(ASSOCIATED(BoundaryConditionsType)) THEN
        CALL CMISSBoundaryConditionsTypeFinalise (BoundaryConditionsType, CMISSBoundaryConditionsTypeFinaliseC)
        DEALLOCATE (BoundaryConditionsType)
        BoundaryConditionsTypePtr = C_NULL_PTR
      ELSE
        CMISSBoundaryConditionsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSBoundaryConditionsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSBoundaryConditionsTypeFinaliseC

  !
  !============================================================================
  !

  !>Initialises a CMISSBoundaryConditionsType object for C.

  FUNCTION CMISSBoundaryConditionsTypeInitialiseC (BoundaryConditionsTypePtr) BIND (C, NAME = &
  & "CMISSBoundaryConditionsTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT (INOUT) :: BoundaryConditionsTypePtr !<C pointer to the CMISSBoundaryConditionsType object to be initialised.
    !Function variables
    INTEGER(C_INT) :: CMISSBoundaryConditionsTypeInitialiseC !<Error Code.
    !Local Variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSBoundaryConditionsType), POINTER :: BoundaryConditionsType

    IF (C_ASSOCIATED(BoundaryConditionsTypePtr)) THEN
      CMISSBoundaryConditionsTypeInitialiseC = CMISSPointerNotNull
    ELSE
      NULLIFY (BoundaryConditionsType)
      ALLOCATE(BoundaryConditionsType, STAT = Err)
      IF (Err /= 0) THEN
        CMISSBoundaryConditionsTypeInitialiseC =CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsType, CMISSBoundaryConditionsTypeInitialiseC)
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSBoundaryConditionsTypeInitialiseC

  !
  !===========================================================================
  !

  !>Finalises a CMISSControlLoopType object for C.

  FUNCTION CMISSControlLoopTypeFinaliseC (ControlLoopTypePtr) BIND (C, NAME = "CMISSControlLoopTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT (INOUT) :: ControlLoopTypePtr !<C pointer to the CMISSControlLoopType object to be finalised.
    !Function variables
    INTEGER(C_INT) :: CMISSControlLoopTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSControlLoopType), POINTER :: ControlLoopType

    CMISSControlLoopTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(ControlLoopTypePtr)) THEN
      CALL C_F_POINTER(ControlLoopTypePtr, ControlLoopType)
      IF(ASSOCIATED(ControlLoopType)) THEN
        CALL CMISSControlLoopTypeFinalise(ControlLoopType,CMISSControlLoopTypeFinaliseC)
        DEALLOCATE(ControlLoopType)
        ControlLoopTypePtr = C_NULL_PTR
      ELSE
        CMISSControlLoopTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSControlLoopTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSControlLoopTypeFinaliseC

  !
  !===========================================================================
  !

  !>Initialises a CMISSControlLoopType object for C.

  FUNCTION CMISSControlLoopTypeInitialiseC (ControlLoopTypePtr) BIND(C, NAME = "CMISSControlLoopTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT (INOUT) :: ControlLoopTypePtr !<C pointer to the CMISSControlLoopType object to be intialised.
    !Function variables
    INTEGER(C_INT) :: CMISSControlLoopTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE (CMISSControlLoopType), POINTER :: ControlLoopType

    IF(C_ASSOCIATED(ControlLoopTypePtr)) THEN
      CMISSControlLoopTypeInitialiseC=CMISSPointerNotNULL
    ELSE
      NULLIFY(ControlLoopType)
      ALLOCATE(ControlLoopType,STAT=Err)
      IF(Err/=0) THEN
        CMISSControlLoopTypeInitialiseC=CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSControlLoopTypeInitialise(ControlLoopType,CMISSControlLoopTypeInitialiseC)
        ControlLoopTypePtr=C_LOC(ControlLoopType)
        CMISSControlLoopTypeInitialiseC=CMISSNoError
      ENDIF
    ENDIF


    RETURN

  END FUNCTION CMISSControlLoopTypeInitialiseC

  !
  !===========================================================================
  !

  !>Finalises a CMISSCoordinateSystemType object for C.
  FUNCTION CMISSCoordinateSystemTypeFinaliseC(CoordinateSystemTypePtr)  BIND(C,NAME="CMISSCoordinateSystemTypeFinalise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: CoordinateSystemTypePtr !<C pointer to the CMISSCoordinateSystemType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSCoordinateSystemTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSCoordinateSystemType), POINTER :: CoordinateSystemType
    
    CMISSCoordinateSystemTypeFinaliseC=CMISSNoError
    IF(C_ASSOCIATED(CoordinateSystemTypePtr)) THEN
      CALL C_F_POINTER(CoordinateSystemTypePtr,CoordinateSystemType)
      IF(ASSOCIATED(CoordinateSystemType)) THEN
        CALL CMISSCoordinateSystemTypeFinalise(CoordinateSystemType,CMISSCoordinateSystemTypeFinaliseC)
        DEALLOCATE(CoordinateSystemType)
        CoordinateSystemTypePtr=C_NULL_PTR
      ELSE
        CMISSCoordinateSystemTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSCoordinateSystemTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSCoordinateSystemTypeFinaliseC
 
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSCoordinateSystemType object for C.
  FUNCTION CMISSCoordinateSystemTypeInitialiseC(CoordinateSystemTypePtr)  BIND(C,NAME="CMISSCoordinateSystemTypeInitialise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: CoordinateSystemTypePtr !<C pointer to the CMISSCoordinateSystemType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSCoordinateSystemTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSCoordinateSystemType), POINTER :: CoordinateSystemType
    
    IF(C_ASSOCIATED(CoordinateSystemTypePtr)) THEN
      CMISSCoordinateSystemTypeInitialiseC=CMISSPointerNotNULL
    ELSE
      NULLIFY(CoordinateSystemType)
      ALLOCATE(CoordinateSystemType,STAT=Err)
      IF(Err/=0) THEN
        CMISSCoordinateSystemTypeInitialiseC=CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystemType,CMISSCoordinateSystemTypeInitialiseC)
        CoordinateSystemTypePtr=C_LOC(CoordinateSystemType)
        CMISSCoordinateSystemTypeInitialiseC=CMISSNoError
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSCoordinateSystemTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSDecompositionType object.
  FUNCTION CMISSDecompositionTypeFinaliseC(DecompositionTypePtr) BIND(C, NAME = "CMISSDecompositionTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: DecompositionTypePtr !<C pointer to the CMISSDecompositionType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeFinaliseC !Error Code.
    !Local Variables
    TYPE(CMISSDecompositionType), POINTER :: DecompositionType

    CMISSDecompositionTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionTypePtr)) THEN
      CALL C_F_POINTER (DecompositionTypePtr, DecompositionType)
      IF(ASSOCIATED(DecompositionType)) THEN
        CALL CMISSDecompositionTypeFinalise (DecompositionType, CMISSDecompositionTypeFinaliseC)
        DEALLOCATE(DecompositionType)
        DecompositionTypePtr = C_NULL_PTR
      ELSE
        CMISSDecompositionTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSDecompositionType object.
  FUNCTION CMISSDecompositionTypeInitialiseC(DecompositionTypePtr)  BIND(C, NAME = "CMISSDecompositionTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: DecompositionTypePtr !<C pointer to the CMISSDecompositionType object to be initialised.
    !Function variables
    INTEGER(C_INT) ::  CMISSDecompositionTypeInitialiseC !<Error Code.
    !Local Variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSDecompositionType), POINTER :: DecompositionType

    IF(C_ASSOCIATED(DecompositionTypePtr)) THEN
      CMISSDecompositionTypeInitialiseC = CMISSPointerNotNull
    ELSE
      NULLIFY(DecompositionType)
      ALLOCATE(DecompositionType, STAT=Err)
      IF (Err/=0) THEN
        CMISSDecompositionTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSDecompositionTypeInitialise(DecompositionType, CMISSDecompositionTypeInitialiseC)
        DecompositionTypePtr = C_LOC(DecompositionType)
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSEquationsType object for C.
  FUNCTION CMISSEquationsTypeFinaliseC (EquationsTypePtr) BIND(C, NAME = "CMISSEquationsTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: EquationsTypePtr !<C pointer to the CMISSEquationsType object to be finalised.
    !Function variables
    INTEGER(C_INT) :: CMISSEquationsTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSEquationsType), POINTER :: EquationsType

    CMISSEquationsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(EquationsTypePtr)) THEN
      CALL C_F_POINTER(EquationsTypePtr, EquationsType)
      IF (ASSOCIATED(EquationsType)) THEN
        CALL CMISSEquationsTypeFinalise (EquationsType, CMISSEquationsTypeFinaliseC)
        DEALLOCATE(EquationsType)
        EquationsTypePtr = C_NULL_PTR
      ELSE
        CMISSEquationsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSEquationsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

      RETURN

    END FUNCTION CMISSEquationsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSEquationsType object for C.
  FUNCTION CMISSEquationsTypeInitialiseC (EquationsTypePtr) BIND(C, NAME = "CMISSEquationsTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: EquationsTypePtr  !<C pointer to the CMISSEquationsType object to be intialised.
    !Function variables
    INTEGER(C_INT) :: CMISSEquationsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSEquationsType), POINTER :: EquationsType

    IF(C_ASSOCIATED(EquationsTypePTR)) THEN
      CMISSEquationsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(EquationsType)
      ALLOCATE(EquationsType, STAT= Err)
      IF(Err/=0) THEN
        CMISSEquationsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSEquationsTypeInitialise(EquationsType, CMISSEquationsTypeInitialiseC)
        EquationsTypePtr = C_LOC(EquationsType)
        CMISSEquationsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN
  END FUNCTION CMISSEquationsTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSEquationsSetType object for C.
  FUNCTION CMISSEquationsSetTypeFinaliseC (EquationsSetTypePtr) BIND(C, NAME = "CMISSEquationsSetTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: EquationsSetTypePtr  !<C pointer to the CMISSEquationsSetType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSEquationsSetTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSEquationsSetType), POINTER :: EquationsSetType

    CMISSEquationsSetTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(EquationsSetTypePtr)) THEN
      CALL C_F_POINTER(EquationsSetTypePtr, EquationsSetType)
      IF (ASSOCIATED(EquationsSetType)) THEN
        CALL CMISSEquationsSetTypeFinalise(EquationsSetType, CMISSEquationsSetTypeFinaliseC)
        DEALLOCATE(EquationsSetType)
        EquationsSetTypePtr = C_NULL_PTR
      ELSE
        CMISSEquationsSetTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSEquationsSetTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSEquationsSetTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSEquationsSetType object for C.
  FUNCTION CMISSEquationsSetTypeInitialiseC(EquationsSetTypePtr) BIND(C, NAME = "CMISSEquationsSetTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: EquationsSetTypePtr  !<C pointer to the CMISSEquationsSetType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSEquationsSetTypeInitialiseC !<Error Code.
    !Local variable
    INTEGER(C_INT) :: Err
    TYPE(CMISSEquationsSetType), POINTER :: EquationsSetType

    IF(C_ASSOCIATED(EquationsSetTypePtr)) THEN
      CMISSEquationsSetTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (EquationsSetType)
      ALLOCATE(EquationsSetType, STAT = Err)
      IF(Err/=0) THEN
        CMISSEquationsSetTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSEquationsSetTypeInitialise(EquationsSetType, CMISSEquationsSetTypeInitialiseC)
        EquationsSetTypePtr = C_LOC(EquationsSetType)
        CMISSEquationsSetTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSEquationsSetTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSFieldType object for C.
  FUNCTION CMISSFieldTypeFinaliseC(FieldTypePtr) BIND(C,NAME="CMISSFieldTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldTypePtr !<C pointer to the CMISSFieldType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: FieldType

    CMISSFieldTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(FieldTypePtr)) THEN
      CALL C_F_POINTER (FieldTypePtr, FieldType)
      IF(ASSOCIATED(FieldType)) THEN
        CALL CMISSFieldTypeFinalise(FieldType, CMISSFieldTypeFinaliseC)
        DEALLOCATE(FieldType)
        FieldTypePtr = C_NULL_PTR
     ELSE
        CMISSFieldTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSFieldType object for C.
  FUNCTION CMISSFieldTypeInitialiseC(FieldTypePtr) BIND(C, NAME = "CMISSFieldTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: FieldTypePtr !<C pointer to the CMISSFieldType object to be finalised.
    !Function Variable
    INTEGER(C_INT) :: CMISSFieldTypeInitialiseC !<Error Code
    !Local variable
    INTEGER(C_INT) :: Err
    TYPE(CMISSFieldType), POINTER :: FieldType

    IF(C_ASSOCIATED(FieldTypePtr)) THEN
      CMISSFieldTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(FieldType)
      ALLOCATE(FieldType, STAT = Err)
      IF(Err/=0) THEN
        CMISSFieldTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSFieldTypeInitialise (FieldType, CMISSFieldTypeInitialiseC)
        FieldTypePtr = C_LOC(FieldType)
        CMISSFieldTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN
  END FUNCTION CMISSFieldTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Creates a pointer to a CMISSFieldsType object for an object reference for C.
  FUNCTION CMISSFieldsTypeCreateC(RegionPtr, FieldsPtr) BIND(C, NAME="CMISSFieldsTypeCreate")

  !Argument variables
  TYPE(C_PTR), INTENT(IN) :: RegionPtr !<C pointer to the region to get the fields from.
  TYPE(C_PTR), INTENT(INOUT) :: FieldsPtr !<C pointer to the fields attached to the specified region.
  !Function variable
  INTEGER(C_INT) :: CMISSFieldsTypeCreateC !<Error Code.
  !Local variables
  INTEGER(C_INT) :: Err
  TYPE(CMISSRegionType), POINTER :: Region
  TYPE(CMISSFieldsType), POINTER :: Fields

  IF(C_ASSOCIATED(RegionPtr)) THEN
    IF (C_ASSOCIATED(FieldsPtr)) THEN
      CMISSFieldsTypeCreateC = CMISSPointerNotNULL
    ELSE
      NULLIFY(Fields)
      ALLOCATE(Fields, STAT= Err)
      IF(Err/=0) THEN
        CMISSFieldsTypeCreateC = CMISSCouldNotAllocatePointer
      ELSE
        CALL C_F_POINTER (RegionPtr, Region)
        IF(ASSOCIATED(Region)) THEN
          CALL CMISSFieldsTypeCreate(Region, Fields, CMISSFieldsTypeCreateC)
          FieldsPtr = C_LOC(Fields)
        ELSE
          CMISSFieldsTypeCreateC = CMISSErrorConvertingPointer
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  RETURN

END FUNCTION CMISSFieldsTypeCreateC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSFieldsType object for C.

  FUNCTION CMISSFieldsTypeFinaliseC (FieldsTypePtr) BIND(C,NAME="CMISSFieldsTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldsTypePtr !<C pointer to the CMISSFieldsType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldsTypeFinaliseC !<Error code.
    !Local variables
    TYPE(CMISSFieldsType), POINTER :: FieldsType

    CMISSFieldsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(FieldsTypePtr)) THEN
      CALL C_F_POINTER(FieldsTypePtr, FieldsType)
      IF(ASSOCIATED(FieldsType)) THEN
        CALL CMISSFieldsTypeFinalise(FieldsType, CMISSFieldsTypeFinaliseC)
        DEALLOCATE(FieldsType)
        FieldsTypePtr = C_NULL_PTR
      ELSE
        CMISSFieldsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSFieldsType object for C.

  FUNCTION CMISSFieldsTypeInitialiseC(FieldsTypePtr) BIND(C,NAME = "CMISSFieldsTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldsTypePtr !<C pointer to the CMISSFieldsType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSFieldsType), POINTER :: FieldsType

    IF(C_ASSOCIATED(FieldsTypePtr)) THEN
      CMISSFieldsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (FieldsType)
      ALLOCATE(FieldsType, STAT = Err)
      IF(Err/=0) THEN
        CMISSFieldsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSFieldsTypeInitialise(FieldsType, CMISSFieldsTypeInitialiseC)
        FieldsTypePtr = C_LOC(FieldsType)
        CMISSFieldsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldsTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSGeneratedMeshType object for C.

  FUNCTION CMISSGeneratedMeshTypeFinaliseC(GeneratedMeshTypePtr) BIND(C, NAME = "CMISSGeneratedMeshTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshTypePtr !<C pointer to the CMISSGeneratedMeshType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMeshType

    CMISSGeneratedMeshTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshTypePtr)) THEN
      CALL C_F_POINTER(GeneratedMeshTypePtr, GeneratedMeshType)
      IF(ASSOCIATED(GeneratedMeshType)) THEN
        CALL CMISSGeneratedMeshTypeFinalise(GeneratedMeshType, CMISSGeneratedMeshTypeFinaliseC)
        DEALLOCATE(GeneratedMeshType)
        GeneratedMeshTypePtr = C_NULL_PTR
      ELSE
        CMISSGeneratedMeshTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSGeneratedMeshType object for C.

  FUNCTION CMISSGeneratedMeshTypeInitialiseC (GeneratedMeshTypePtr) BIND(C, NAME = "CMISSGeneratedMeshTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshTypePtr !<C pointer to the CMISSGeneratedMeshType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMeshType

    IF(C_ASSOCIATED(GeneratedMeshTypePtr)) THEN
      CMISSGeneratedMeshTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (GeneratedMeshType)
      ALLOCATE(GeneratedMeshType, STAT = Err)
      IF(Err/=0) THEN
        CMISSGeneratedMeshTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSGeneratedMeshTypeInitialise(GeneratedMeshType, CMISSGeneratedMeshTypeInitialiseC)
        GeneratedMeshTypePtr = C_LOC(GeneratedMeshType)
        CMISSGeneratedMeshTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSHistoryType object for C.

  FUNCTION CMISSHistoryTypeFinaliseC(HistoryTypePtr) BIND(C, NAME = "CMISSHistoryTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: HistoryTypePtr !<C pointer to the CMISSHistoryType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSHistoryTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSHistoryType), POINTER :: HistoryType

    CMISSHistoryTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(HistoryTypePtr)) THEN
      CALL C_F_POINTER(HistoryTypePtr, HistoryType)
      IF(ASSOCIATED(HistoryType)) THEN
        CALL CMISSHistoryTypeFinalise(HistoryType, CMISSHistoryTypeFinaliseC)
        DEALLOCATE(HistoryType)
        HistoryTypePtr = C_NULL_PTR
      ELSE
        CMISSHistoryTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSHistoryTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSHistoryTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSHistoryType object for C.
  FUNCTION CMISSHistoryTypeInitialiseC(HistoryTypePtr) BIND(C,NAME="CMISSHistoryTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: HistoryTypePtr !<C pointer to the CMISSHistoryType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSHistoryTypeInitialiseC !<Error Code.
    !Local varaibles
    INTEGER(C_INT) :: Err
    TYPE(CMISSHistoryType), POINTER :: HistoryType

    IF(C_ASSOCIATED(HistoryTypePtr)) THEN
      CMISSHistoryTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(HistoryType)
      ALLOCATE(HistoryType, STAT=Err)
      IF(Err/=0) THEN
        CMISSHistoryTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSHistoryTypeInitialise(HistoryType, CMISSHistoryTypeInitialiseC)
        HistoryTypePtr = C_LOC(HistoryType)
        CMISSHistoryTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSHistoryTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSMeshType object for C.
  FUNCTION CMISSMeshTypeFinaliseC(MeshTypePtr) BIND(C,NAME="CMISSMeshTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: MeshTypePtr !<C pointer to the CMISSMeshType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSMeshType), POINTER :: MeshType

    CMISSMeshTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(MeshTypePtr))THEN
      CALL C_F_POINTER(MeshTypePtr, MeshType)
      IF(ASSOCIATED(MeshType)) THEN
        CALL CMISSMeshTypeFinalise(MeshType, CMISSMeshTypeFinaliseC)
        DEALLOCATE(MeshType)
        MeshTypePtr = C_NULL_PTR
      ELSE
        CMISSMeshTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSMeshTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSMeshType object for C.
  FUNCTION CMISSMeshTypeInitialiseC (MeshTypePtr)BIND(C,NAME= "CMISSMeshTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: MeshTypePtr
    !Function variable
    INTEGER(C_INT) :: CMISSMeshTypeInitialiseC
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSMeshType), POINTER :: MeshType

    IF(C_ASSOCIATED(MeshTypePtr)) THEN
      CMISSMeshTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(MeshType)
      ALLOCATE(MeshType, STAT=Err)
      IF(Err/=0) THEN
        CMISSMeshTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSMeshTypeInitialise(MeshType, CMISSMeshTypeInitialiseC)
        MeshTypePtr = C_LOC(MeshType)
        CMISSMeshTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSMeshTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSMeshElementsType object for C.

  FUNCTION CMISSMeshElementsTypeFinaliseC(MeshElementsTypePtr) BIND(C, NAME="CMISSMeshElementsTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: MeshElementsTypePtr !<C pointer to the CMISSMeshElementsType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSMeshElementsType), POINTER :: MeshElementsType

    CMISSMeshElementsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsTypePtr)) THEN
      CALL C_F_POINTER(MeshElementsTypePtr, MeshElementsType)
        IF(ASSOCIATED(MeshElementsType)) THEN
          CALL CMISSMeshElementsTypeFinalise(MeshElementsType, CMISSMeshElementsTypeFinaliseC)
          DEALLOCATE(MeshElementsType)
          MeshElementsTypePtr = C_NULL_PTR
      ELSE
        CMISSMeshElementsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

      RETURN

    END FUNCTION CMISSMeshElementsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSMeshElementsType object for C.
  FUNCTION CMISSMeshElementsTypeInitialiseC(MeshElementsTypePtr) BIND(C, NAME = "CMISSMeshElementsTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT (INOUT) :: MeshElementsTypePtr !<C pointer to the CMISSMeshElementsType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSMeshElementsType), POINTER :: MeshElementsType

    IF(C_ASSOCIATED(MeshElementsTypePtr)) THEN
      CMISSMeshElementsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(MeshElementsType)
      ALLOCATE(MeshElementsType, STAT = Err)
      IF(Err/=0) THEN
        CMISSMeshElementsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSMeshElementsTypeInitialise(MeshElementsType, CMISSMeshElementsTypeInitialiseC)
        MeshElementsTypePtr = C_LOC(MeshElementsType)
        CMISSMeshElementsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSNodesType object for C.
  FUNCTION CMISSNodesTypeFinaliseC(NodesTypePtr) BIND(C, NAME="CMISSNodesTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: NodesTypePtr !<C pointer to the CMISSNodesType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSNodesType), POINTER :: NodesType

    CMISSNodesTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(NodesTypePtr)) THEN
      CALL C_F_POINTER(NodesTypePtr, NodesType)
      IF(ASSOCIATED(NodesType)) THEN
        CALL CMISSNodesTypeFinalise(NodesType, CMISSNodesTypeFinaliseC)
        DEALLOCATE(NodesType)
        NodesTypePtr = C_NULL_PTR
      ELSE
        CMISSNodesTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSNodesType object for C.
  FUNCTION CMISSNodesTypeInitialiseC(NodesTypePtr) BIND(C, NAME= "CMISSNodesTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: NodesTypePtr !<C pointer to the CMISSNodesType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSNodesType), POINTER :: NodesType

    IF(C_ASSOCIATED(NodesTypePtr)) THEN
      CMISSNodesTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(NodesType)
      ALLOCATE(NodesType, STAT = Err)
      IF(Err/=0) THEN
        CMISSNodesTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSNodesTypeInitialise(NodesType, CMISSNodesTypeInitialiseC)
        NodesTypePtr = C_NULL_PTR
        CMISSNodesTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSNodesTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSProblemType object for C.
  FUNCTION CMISSProblemTypeFinaliseC(ProblemTypePtr) BIND(C, NAME="CMISSProblemTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: ProblemTypePtr !<C pointer to the CMISSProblemType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSProblemType), POINTER :: ProblemType

    CMISSProblemTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(ProblemTypePtr)) THEN
      CALL C_F_POINTER(ProblemTypePtr, ProblemType)
      IF(ASSOCIATED(ProblemType)) THEN
        CALL CMISSProblemTypeFinalise(ProblemType, CMISSProblemTypeFinaliseC)
        DEALLOCATE(ProblemType)
        ProblemTypePtr = C_NULL_PTR
      ELSE
        CMISSProblemTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSProblemType object for C.
  FUNCTION CMISSProblemTypeInitialiseC(ProblemTypePtr) BIND(C, NAME = "CMISSProblemTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: ProblemTypePtr !<C pointer to the CMISSProblemType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSProblemType), POINTER :: ProblemType

    IF(C_ASSOCIATED(ProblemTypePtr)) THEN
      CMISSProblemTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(ProblemType)
      ALLOCATE(ProblemType, STAT=Err)
      IF(Err/=0) THEN
        CMISSProblemTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSProblemTypeInitialise(ProblemType, CMISSProblemTypeInitialiseC)
        ProblemTypePtr = C_LOC(ProblemType)
        CMISSProblemTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSProblemTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSQuadratureType object for C.
  FUNCTION CMISSQuadratureTypeFinaliseC(QuadratureTypePtr) BIND(C, NAME="CMISSQuadratureTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: QuadratureTypePtr !<C pointer to the CMISSQuadratureType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSQuadratureTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSQuadratureType), POINTER :: QuadratureType

    CMISSQuadratureTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(QuadratureTypePtr)) THEN
      CALL C_F_POINTER(QuadratureTypePtr, QuadratureType)
      IF(ASSOCIATED(QuadratureType)) THEN
        CALL CMISSQuadratureTypeFinalise(QuadratureType, CMISSQuadratureTypeFinaliseC)
        DEALLOCATE(QuadratureType)
        QuadratureTypePtr = C_NULL_PTR
      ELSE
        CMISSQuadratureTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSQuadratureTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSQuadratureTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSQuadratureType object for C.
  FUNCTION CMISSQuadratureTypeInitialiseC(QuadratureTypePtr) BIND(C,NAME="CMISSQuadratureTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: QuadratureTypePtr !<C pointer to the CMISSQuadratureType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSQuadratureTypeInitialiseC !<Error Code.
    !Local variable
    INTEGER(C_INT) :: Err
    TYPE(CMISSQuadratureType), POINTER :: QuadratureType

    IF(C_ASSOCIATED(QuadratureTypePtr)) THEN
      CMISSQuadratureTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(QuadratureType)
      ALLOCATE(QuadratureType, STAT=Err)
      IF(Err/=0) THEN
        CMISSQuadratureTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSQuadratureTypeInitialise(QuadratureType, CMISSQuadratureTypeInitialiseC)
        QuadratureTypePtr = C_LOC(QuadratureType)
        CMISSQuadratureTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

  END FUNCTION CMISSQuadratureTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSRegionType object.
  FUNCTION CMISSRegionTypeFinaliseC(RegionTypePtr) BIND(C, NAME="CMISSRegionTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: RegionTypePtr !<C pointer to the CMISSRegionType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSRegionTypeFinaliseC !<Error code.
    !Local variable
    TYPE(CMISSRegionType), POINTER :: RegionType

    CMISSRegionTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(RegionTypePtr)) THEN
      CALL C_F_POINTER(RegionTypePtr, RegionType)
      IF(ASSOCIATED(RegionType)) THEN
        CALL CMISSRegionTypeFinalise(RegionType, CMISSRegionTypeFinaliseC)
        DEALLOCATE(RegionType)
        RegionTypePtr = C_NULL_PTR
      ELSE
        CMISSRegionTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSRegionTypeFinaliseC
 
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSRegionType object for C.
  FUNCTION CMISSRegionTypeInitialiseC(RegionTypePtr) BIND(C,NAME="CMISSRegionTypeInitialise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: RegionTypePtr  !<C pointer to the CMISSRegionType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSRegionTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSRegionType), POINTER :: RegionType
    
    IF(C_ASSOCIATED(RegionTypePtr)) THEN
      CMISSRegionTypeInitialiseC=CMISSPointerNotNULL
    ELSE
      NULLIFY(RegionType)
      ALLOCATE(RegionType,STAT=Err)
      IF(Err/=0) THEN
        CMISSRegionTypeInitialiseC=CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSRegionTypeInitialise(RegionType,CMISSRegionTypeInitialiseC)
        RegionTypePtr=C_LOC(RegionType)
        CMISSRegionTypeInitialiseC=CMISSNoError
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionTypeInitialiseC
  
  !
  !================================================================================================================================
  !

  !>Finalises a CMISSSolverType object for C.
  FUNCTION CMISSSolverTypeFinaliseC(SolverTypePtr) BIND(C, NAME="CMISSSolverTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverTypePtr !<C pointer to the CMISSSolverType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSSolverType), POINTER :: SolverType

    CMISSSolverTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(SolverTypePtr)) THEN
      CALL C_F_POINTER(SolverTypePtr, SolverType)
      IF(ASSOCIATED(SolverType)) THEN
        CALL CMISSSolverTypeFinalise(SolverType, CMISSSolverTypeFinaliseC)
        DEALLOCATE(SolverType)
        SolverTypePtr = C_NULL_PTR
      ELSE
        CMISSSolverTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSSolverTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSSolverTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSSolverType object for C.
  FUNCTION CMISSSolverTypeInitialiseC(SolverTypePtr) BIND(C, NAME="CMISSSolverTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverTypePtr !<C pointer to the CMISSSolverType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSSolverType), POINTER :: SolverType

    IF(C_ASSOCIATED(SolverTypePtr)) THEN
      CMISSSolverTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(SolverType)
      ALLOCATE(SolverType, STAT=Err)
      IF(Err/=0) THEN
        CMISSSolverTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSSolverTypeInitialise(SolverType, CMISSSolverTypeInitialiseC)
        SolverTypePtr = C_LOC(SolverType)
        CMISSSolverTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSSolverTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSSolverEquationsType object for C.
  FUNCTION CMISSSolverEquationsTypeFinaliseC (SolverEquationsTypePtr) BIND(C, NAME="CMISSSolverEquationsTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverEquationsTypePtr !<C pointer to the CMISSSolverEquationsType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverEquationsTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSSolverEquationsType), POINTER :: SolverEquationsType

    CMISSSolverEquationsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(SolverEquationsTypePtr)) THEN
      CALL C_F_POINTER(SolverEquationsTypePtr, SolverEquationsType)
      IF(ASSOCIATED(SolverEquationsType)) THEN
        CALL CMISSSolverEquationsTypeFinalise(SolverEquationsType, CMISSSolverEquationsTypeFinaliseC)
        DEALLOCATE(SolverEquationsType)
        SolverEquationsTypePtr = C_NULL_PTR
      ELSE
        CMISSSolverEquationsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSSolverEquationsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSSolverEquationsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSSolverEquationsType object for C.
  FUNCTION CMISSSolverEquationsTypeInitialiseC(SolverEquationsTypePtr) BIND(C, NAME="CMISSSolverEquationsTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverEquationsTypePtr !<C pointer to the CMISSSolverEquationsType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverEquationsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSSolverEquationsType), POINTER :: SolverEquationsType

    IF(C_ASSOCIATED(SolverEquationsTypePtr)) THEN
      CMISSSolverEquationsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(SolverEquationsType)
      ALLOCATE(SolverEquationsType, STAT=Err)
      IF(Err/=0) THEN
        CMISSSolverEquationsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSSolverEquationsTypeInitialise(SolverEquationsType, CMISSSolverEquationsTypeInitialiseC)
        SolverEquationsTypePtr = C_LOC(SolverEquationsType)
        CMISSSolverEquationsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSSolverEquationsTypeInitialiseC


!!==================================================================================================================================
!!
!! ANALYTIC_ANALYSIS_ROUTINES
!!
!!==================================================================================================================================

  !>Output the analytic error analysis for a field specified by a user number compared to the analytic values parameter set for C.
  FUNCTION CMISSAnalyticAnalysisOutputNumberC(RegionUserNumber, FieldUserNumber, FileName, FileNameSize) BIND(C, NAME = &
  & "CMISSAnalyticAnalysisOutputNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to
    INTEGER(C_INT),VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to calculate the analytic error analysis for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !<The file name size, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<If not empty, the filename to output the analytic analysis to, for C. If empty, the analysis will be output to the standard output.
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisOutputNumberC !<Error Code.
    !Local variables
    CHARACTER(LEN=FileNameSize-1) :: FFileName

    CALL CMISSC2FString(FileName, FFileName)
    CALL CMISSAnalyticAnalysisOutputNumber(RegionUserNumber, FieldUserNumber, FFilename, CMISSAnalyticAnalysisOutputNumberC)

    RETURN

  END FUNCTION CMISSAnalyticAnalysisOutputNumberC

  !
  !================================================================================================================================
  !

  !>Output the analytic error analysis for a field identified by an object compared to the analytic values parameter set.
  FUNCTION CMISSAnalyticAnalysisOutputPtrC(FieldPtr,FileNameSize, FileName) BIND(C, NAME = "CMISSAnalyticAnalysisOutputObj")

    !Argument variables
    TYPE(C_PTR), INTENT(OUT) :: FieldPtr !<The dependent field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !<The file name size, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<If not empty, the filename to output the analytic analysis to, for C. If empty, the analysis will be output to the standard output.
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisOutputPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = FileNameSize-1) :: FFileName

    CMISSAnalyticAnalysisOutputPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER (FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Filename, FFileName)
        CALL CMISSAnalyticAnalysisOutputObj(Field, FFileName, CMISSAnalyticAnalysisOutputPtrC)
      ELSE
        CMISSAnalyticAnalysisOutputPtrC=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSAnalyticAnalysisOutputPtrC=CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSAnalyticAnalysisOutputPtrC


  !
  !================================================================================================================================
  !

  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************




!!==================================================================================================================================
!!
!! CMISS_ROUTINES
!!
!!==================================================================================================================================

  !>Finalises CMISS for C.
  FUNCTION CMISSFinaliseC() BIND(C,NAME="CMISSFinalise")
  
    !Argument variables
    !Function variable
    INTEGER(C_INT) :: CMISSFinaliseC
    !Local variables

    CALL CMISSFinalise(CMISSFinaliseC)

    RETURN
    
  END FUNCTION CMISSFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises CMISS for C returning a user number to the world coordinate system and region.
  FUNCTION CMISSInitialiseCNum(WorldCoordinateSystemUserNumber,WorldRegionUserNumber) BIND(C,NAME="CMISSInitialiseNum")
  
    !Argument variables
    INTEGER(C_INT), INTENT(OUT) :: WorldCoordinateSystemUserNumber
    INTEGER(C_INT), INTENT(OUT) :: WorldRegionUserNumber
    !Function variable
    INTEGER(C_INT) :: CMISSInitialiseCNum
    !Local variables
  
    CALL CMISSInitialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,CMISSInitialiseCNum)

    RETURN
    
    END FUNCTION CMISSInitialiseCNum

  !
  !================================================================================================================================
  !

  !>Initialises CMISS for C returning pointers to the world coordinate system and region.
  FUNCTION CMISSInitialiseCPtr(WorldCoordinateSystemPtr,WorldRegionPtr) BIND(C,NAME="CMISSInitialise")
  
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: WorldCoordinateSystemPtr
    TYPE(C_PTR), INTENT(INOUT) :: WorldRegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSInitialiseCPtr
    !Local variables
    TYPE(CMISSCoordinateSystemType), POINTER :: WorldCoordinateSystem
    TYPE(CMISSRegionType), POINTER :: WorldRegion

    CMISSInitialiseCPtr=CMISSCoordinateSystemTypeInitialiseC(WorldCoordinateSystemPtr)
    IF(CMISSInitialiseCPtr==CMISSNoError) THEN
      CMISSInitialiseCPtr=CMISSRegionTypeInitialiseC(WorldRegionPtr)
      IF(CMISSInitialiseCPtr==CMISSNoError) THEN
        CALL C_F_POINTER(WorldCoordinateSystemPtr,WorldCoordinateSystem)
        IF(ASSOCIATED(WorldCoordinateSystem)) THEN
          CALL C_F_POINTER(WorldRegionPtr,WorldRegion)
          IF(ASSOCIATED(WorldRegion)) THEN        
            CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,CMISSInitialiseCPtr)
          ELSE
            CMISSInitialiseCPtr=CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSInitialiseCPtr=CMISSErrorConvertingPointer
        ENDIF
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSInitialiseCPtr

!!==================================================================================================================================
!!
!! FIELD_ROUTINES
!!
!!==================================================================================================================================

  !>Returns the interpolation type for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentInterpolationGetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType) BIND(C,NAME = "CMISSFieldComponentInterpolationGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the interpolation type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the interpolation type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) ::  CMISSFieldComponentInterpolationGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldComponentInterpolationGet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType, CMISSFieldComponentInterpolationGetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationGetCNum

  !!
  !!==================================================================================================================================
  !!

  !>Returns the interpolation type for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentInterpolationGetCPtr(FieldPtr,VariableType,ComponentNumber,InterpolationType) BIND(C, NAME = &
  & "CMISSFieldComponentInterpolationGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the interpolation type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the interpolation type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentInterpolationGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentInterpolationGet(Field, VariableType,ComponentNumber,InterpolationType, &
        & CMISSFieldComponentInterpolationGetCPtr)
      ELSE
        CMISSFieldComponentInterpolationGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentInterpolationGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationGetCPtr

  !!
  !!==================================================================================================================================
  !!

  !>Sets/changes the interpolation type for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentInterpolationSetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType) BIND(C, NAME = "CMISSFieldComponentInterpolationSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the interpolation type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentInterpolationSet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,InterpolationType,&
    & CMISSFieldComponentInterpolationSetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationSetCNum


  !!
  !!==================================================================================================================================
  !!

  !>Sets/changes the interpolation type for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentInterpolationSetCPtr(FieldPtr,VariableType,ComponentNumber,InterpolationType) &
   & BIND(C, NAME = "CMISSFieldComponentInterpolationSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the interpolation type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentInterpolationSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentInterpolationSet(Field, VariableType, ComponentNumber,InterpolationType, &
        & CMISSFieldComponentInterpolationSetCPtr)
      ELSE
        CMISSFieldComponentInterpolationSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentInterpolationSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationSetCPtr


  !!
  !!==================================================================================================================================
  !!

  !>Returns the character string label for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentLabelGetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,LabelSize,Label) &
  & BIND(C, NAME = "CMISSFieldComponentLabelGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable component character string label to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelGetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CALL CMISSFieldComponentLabelGet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,FLabel, &
    & CMISSFieldComponentLabelGetCNum)
    CALL CMISSF2CString(Flabel,Label)

    RETURN

  END FUNCTION CMISSFieldComponentLabelGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentLabelGetCPtr(FieldPtr,VariableType,ComponentNumber,LabelSize,Label) BIND(C, NAME = &
  & "CMISSFieldComponentLabelGet")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable component character string label to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelGetCPtr !<Error Code.
    !Local variable
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentLabelGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentLabelGet(Field,VariableType,ComponentNumber,FLabel,CMISSFieldComponentLabelGetCPtr)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldComponentLabelGetCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentLabelGetCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentLabelGetCPtr

  !
  !================================================================================================================================
  !
  !>Sets/changes the character string label for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentLabelSetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,LabelSize,Label)&
  & BIND(C,NAME= "CMISSFieldComponentLabelSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable component character string label to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelSetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CALL CMISSC2FString(Label,Flabel)
    CALL CMISSFieldComponentLabelSet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,FLabel, &
    & CMISSFieldComponentLabelSetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentLabelSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentLabelSetCPtr(FieldPtr,VariableType,ComponentNumber,LabelSize,Label) BIND(C, NAME= &
  & "CMISSFieldComponentLabelSetC")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label Size
    CHARACTER(LEN = 1, KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable component character string label to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelSetCPtr !<Error Code.
    !Local variable
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentLabelSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldComponentLabelSet(Field,VariableType,ComponentNumber,FLabel,CMISSFieldComponentLabelSetCPtr)
      ELSE
        CMISSFieldComponentLabelSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentLabelSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSFieldComponentLabelSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentMeshComponentGetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & MeshComponent) BIND(C, NAME = "CMISSFieldComponentMeshComponentGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the mesh component number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the mesh component number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the mesh component number for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the mesh component number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponent !<The mesh component number to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentMeshComponentGet(RegionUserNumber, FieldUserNumber, VariableType, ComponentNumber, MeshComponent,&
    & CMISSFieldComponentMeshComponentGetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentMeshComponentGetCPtr(FieldPtr,VariableType,ComponentNumber,MeshComponent) BIND(C, &
  & NAME = "CMISSFieldComponentMeshComponentGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the mesh component number for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the mesh component number for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the mesh component number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponent !<The mesh component number to get, for C.
    !Function Variables
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentMeshComponentGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentMeshComponentGet(Field, VariableType, ComponentNumber, MeshComponent, &
        & CMISSFieldComponentMeshComponentGetCPtr)
      ELSE
        CMISSFieldComponentMeshComponentGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentMeshComponentGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentGetCPtr

  !
  !================================================================================================================================
  !
  !>Sets/changes the mesh component number for a field variable component for a field identified by a user number.
  FUNCTION CMISSFieldComponentMeshComponentSetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & MeshComponent) BIND(C, NAME = "CMISSFieldComponentMeshComponentSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the mesh component number to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the mesh component number to, for C.
    INTEGER(C_INT), INTENT(IN) :: MeshComponent !<The mesh component number to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentMeshComponentSet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,MeshComponent, &
    & CMISSFieldComponentMeshComponentSetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentMeshComponentSetCPtr(FieldPtr,VariableType,ComponentNumber,MeshComponent) BIND(C, &
  & NAME = "CMISSFieldComponentMeshComponentSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the mesh component number to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the mesh component number to. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponent !<The mesh component number to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentMeshComponentSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentMeshComponentSet(Field,VariableType,ComponentNumber,MeshComponent, &
        & CMISSFieldComponentMeshComponentSetCPtr)
      ELSE
        CMISSFieldComponentMeshComponentSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentMeshComponentSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentSetCPtr
  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to an integer constant value for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentValuesInitialiseIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseIntgNum")

        !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldComponentValuesInitialiseIntg(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, &
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseIntgCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseIntgCNum

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to an integer constant value for a field identified by an object for C.
  FUNCTION CMISSFieldComponentValuesInitialiseIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME ="CMISSFieldComponentValuesInitialiseIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseIntg(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseIntgCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseIntgCPtr

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a single precision constant value for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentValuesInitialiseSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for , for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseSP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, &
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseSPCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseSPCNum


  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a single precision constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldComponentValuesInitialiseSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseSPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseSP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseSPCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseSPCPtr

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a double precision constant value for a field identified by a user number, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseDP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType,&
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseDPCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseDPCNum

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a double precision constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldComponentValuesInitialiseDP")

      !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseDPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseDP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseDPCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseDPCPtr

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a logical constant value for a field identified by a user number, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseL(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, ComponentNumber,&
    & Value, CMISSFieldComponentValuesInitialiseLCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseLCNum

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a logical constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldComponentValuesInitialiseL")

    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseL(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseLCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseLCPtr

  !
  !================================================================================================================================
  !

  !>Returns the data type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDataTypeGetCNum(RegionUserNumber,FieldUserNumber,VariableType,DataType) BIND(C,NAME= &
  & "CMISSFieldDataTypeGetNum")
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the data type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DataType !<The field variable data type to get, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeGetCNum !<Error Code.
    !Local variables

    CALL  CMISSFieldDataTypeGet(RegionUserNumber,FieldUserNumber,VariableType,DataType,CMISSFieldDataTypeGetCNum)

    RETURN

  END FUNCTION  CMISSFieldDataTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the data type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDataTypeGetCPtr(FieldPtr,VariableType,DataType) BIND(C, NAME = "CMISSFieldDataTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<The field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the data type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DataType !<The field variable data type to get, for C. \see OPENCMISS_FieldDataTypes
    !Function Variables
    INTEGER(C_INT) :: CMISSFieldDataTypeGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDataTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDataTypeGet(Field, VariableType, DataType, CMISSFieldDataTypeGetCPtr)
      ELSE
        CMISSFieldDataTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDataTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDataTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDataTypeSetCNum(RegionUserNumber,FieldUserNumber,VariableType,DataType) BIND(C, NAME = &
  & "CMISSFieldDataTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the data type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DataType !<The field variable data type to set, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldDataTypeSet(RegionUserNumber, FieldUserNumber,VariableType, DataType,CMISSFieldDataTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldDataTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDataTypeSetCPtr(FieldPtr,VariableType,DataType) BIND(C, NAME = "CMISSFieldDataTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to set the data type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DataType !<The field variable data type to set, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDataTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDataTypeSet(Field, VariableType, DataType,CMISSFieldDataTypeSetCPtr)
      ELSE
        CMISSFieldDataTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDataTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDataTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the DOF order type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDOFOrderTypeGetCNum(RegionUserNumber,FieldUserNumber,VariableType,DOFOrderType) BIND( &
  & C, NAME = "CMISSFieldDOFOrderTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the DOF order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the DOF order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the DOF order type for, for C.\see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DOFOrderType !<The field variable DOF Order type to get, for C.  \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDOFOrderTypeGet(RegionUserNumber, FieldUserNumber, VariableType, DOFOrderType, &
    & CMISSFieldDOFOrderTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the DOF Order type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDOFOrderTypeGetCPtr(FieldPtr,VariableType,DOFOrderType) BIND(C, NAME = "CMISSFieldDOFOrderTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer of the field to get the DOF Order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the DOF Order type for, for C \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DOFOrderType !<The field variable DOF Order type to get, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeGetCPtr
 !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDOFOrderTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDOFOrderTypeGet(Field, VariableType, DOFOrderType,CMISSFieldDOFOrderTypeGetCPtr)
      ELSE
        CMISSFieldDOFOrderTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDOFOrderTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the DOF order type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldDOFOrderTypeSetCNum(RegionUserNumber,FieldUserNumber,VariableType,DOFOrderType) BIND(C, &
  & NAME = "CMISSFieldDOFOrderTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the DOF Order type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DOFOrderType !<The field variable DOF Order type to set, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDOFOrderTypeSet(RegionUserNumber, FieldUserNumber, VariableType, DOFOrderType, &
    & CMISSFieldDOFOrderTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the DOF Order type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldDOFOrderTypeSetCPtr(FieldPtr,VariableType,DOFOrderType) BIND(C, NAME = "CMISSFieldDOFOrderTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the DOF Order type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DOFOrderType !<The field variable DOF Order type to set, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDOFOrderTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDOFOrderTypeSet(Field, VariableType, DOFOrderType, CMISSFieldDOFOrderTypeSetCPtr)
      ELSE
        CMISSFieldDOFOrderTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDOFOrderTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a field identified by a user number for C.
  FUNCTION CMISSFieldCreateFinishCNum(RegionUserNumber,FieldUserNumber) BIND(C, NAME = "CMISSFieldCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to finish the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to finish the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateFinishCNum !<Error Code.

    CALL CMISSFieldCreateFinish(RegionUserNumber,FieldUserNumber,CMISSFieldCreateFinishCNum)

    RETURN

  END FUNCTION CMISSFieldCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a field identified by an object for C.
  FUNCTION CMISSFieldCreateFinishCPtr(FieldPtr) BIND(C, NAME = "CMISSFieldCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finish the creation of.
    !Function variable
    INTEGER(C_INT) ::CMISSFieldCreateFinishCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldCreateFinish(Field,CMISSFieldCreateFinishCPtr)
      ELSE
        CMISSFieldCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the creation of a field identified by a user number for C.
  FUNCTION CMISSFieldCreateStartCNum(FieldUserNumber,RegionUserNumber) BIND(C, NAME = "CMISFieldCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber!<The user number of the region containing the field to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the field to start the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateStartCNum !<Error Code.
    !Local variable

    CALL CMISSFieldCreateStart(FieldUserNumber, RegionUserNumber, CMISSFieldCreateStartCNum)

    RETURN

  END FUNCTION CMISSFieldCreateStartCNum


  !
  !================================================================================================================================
  !

  !>Starts the creation of a field identified by an object for C.
  FUNCTION CMISSFieldCreateStartCPtr(FieldUserNumber,RegionPtr,FieldPtr) BIND(C, NAME ="CMISSFieldCreateStart")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber  !<The user number of the field to start the creation of, for C.
    TYPE(C_PTR), INTENT(IN) :: RegionPtr !<C pointer to the region to create the field on.
    TYPE(C_PTR), INTENT(IN) :: FieldPtr !<C pointer to the created field.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateStartCPtr !<Error Code.
    !Local variable
    TYPE(CMISSRegionType), POINTER :: Region
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr, Region)
      IF(ASSOCIATED(Region)) THEN
        IF(C_ASSOCIATED(FieldPtr)) THEN
          CALL C_F_POINTER(FieldPtr, Field)
          IF(ASSOCIATED(Field)) THEN
            CALL CMISSFieldCreateStart(FieldUserNumber, Region, Field, CMISSFieldCreateStartCPtr)
          ELSE
            CMISSFieldCreateStartCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldCreateStartCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldCreateStartCPtr


  !
  !================================================================================================================================
  !

  !>Returns the dependent type for a field identified by a user number for C.
  FUNCTION CMISSFieldDependentTypeGetCNum(RegionUserNumber,FieldUserNumber,DependentType) BIND(C, NAME = &
  & "CMISSFieldDependentTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to get the dependent type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the dependent type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DependentType !<The field dependent type to get, for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldDependentTypeGet(RegionUserNumber, FieldUserNumber, DependentType, CMISSFieldDependentTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldDependentTypeGetCNum


  !
  !================================================================================================================================
  !

  !>Returns the dependent type for a field identified by an object for C.
  FUNCTION CMISSFieldDependentTypeGetCPtr(FieldPtr,DependentType) BIND(C, NAME = "CMISSFieldDependentTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the dependent type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DependentType !<The field dependent type for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDependentTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF (ASSOCIATED(Field)) THEN
        CALL CMISSFieldDependentTypeGet(Field, DependentType, CMISSFieldDependentTypeGetCPtr)
      ELSE
        CMISSFieldDependentTypeGetCPtr= CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDependentTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDependentTypeGetCPtr


  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field identified by a user number for C.
  FUNCTION CMISSFieldDependentTypeSetCNum(RegionUserNumber,FieldUserNumber,DependentType) BIND(C, NAME  = &
  & "CMISSFieldDependentTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DependentType !<The field dependent type for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldDependentTypeSet(RegionUserNumber,FieldUserNumber,DependentType,CMISSFieldDependentTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldDependentTypeSetCNum


  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field identified by an object for C.
  FUNCTION CMISSFieldDependentTypeSetCPtr(FieldPtr,DependentType) BIND(C, NAME = "CMISSFieldDependentTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DependentType !<The field dependent type, for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDependentTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDependentTypeSet(Field, DependentType, CMISSFieldDependentTypeSetCPtr)
      ELSE
        CMISSFieldDependentTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDependentTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDependentTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Destroys a field identified by a user number for C.
  FUNCTION CMISSFieldDestroyCNum(RegionUserNumber,FieldUserNumber) BIND(C, NAME = "CMISSFieldDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to destroy for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDestroyCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDestroy(RegionUserNumber, FieldUserNumber, CMISSFieldDestroyCNum)

    RETURN

  END FUNCTION CMISSFieldDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys a field identified by an object for C.
  FUNCTION CMISSFieldDestroyCPtr(FieldPtr) BIND(C,NAME = "CMISSFieldDestroy")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDestroyCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF (ASSOCIATED(Field)) THEN
        CALL CMISSFieldDestroy(Field, CMISSFieldDestroyCPtr)
      ELSE
        CMISSFieldDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns the dimension for a field identified by a user number for C.
  FUNCTION CMISSFieldDimensionGetCNum(RegionUserNumber,FieldUserNumber,VariableType,DIMENSION) BIND(C, NAME = &
  & "CMISSFieldDimensionGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the dimension for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDimensionGet(RegionUserNumber, FieldUserNumber, VariableType, DIMENSION, CMISSFieldDimensionGetCNum)

    RETURN

  END FUNCTION CMISSFieldDimensionGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the dimension for a field identified by an object for C.
  FUNCTION CMISSFieldDimensionGetCPtr(FieldPtr,VariableType,Dimension) BIND(C, NAME = "CMISSFieldDimensionGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the dimension for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDimensionGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDimensionGet(Field, VariableType, Dimension, CMISSFieldDimensionGetCPtr)
      ELSE
        CMISSFieldDimensionGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDimensionGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDimensionGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension for a field identified by a user number for C.
  FUNCTION CMISSFieldDimensionSetCNum(RegionUserNumber,FieldUserNumber,VariableType,Dimension) BIND(C, NAME = &
  & "CMISSFieldDimensionSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the dimension to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDimensionSet(RegionUserNumber, FieldUserNumber, VariableType, Dimension, CMISSFieldDimensionSetCNum)

    RETURN

  END FUNCTION CMISSFieldDimensionSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension for a field identified by an object for C.
  FUNCTION CMISSFieldDimensionSetCPtr(FieldPtr,VariableType,Dimension) BIND(C, NAME= "CMISSFieldDimensionSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the dimension to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: Dimension !<The field dimension, for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDimensionSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDimensionSet(Field, VariableType, Dimension, CMISSFieldDimensionSetCPtr)
      ELSE
        CMISSFieldDimensionSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDimensionSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDimensionSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the geometric field for a field identified by a user number for C.
  FUNCTION CMISSFieldGeometricFieldGetCNum(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber) BIND(C, NAME = &
  & "CMISSFieldGeometricFieldGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the geometric field for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the geometric field for, for C.
    INTEGER(C_INT), INTENT(OUT) :: GeometricFieldUserNumber !<The field geometric field user number, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldGeometricFieldGet(RegionUserNumber, FieldUserNumber, GeometricFieldUserNumber, &
    & CMISSFieldGeometricFieldGetCNum)

    RETURN

  END FUNCTION CMISSFieldGeometricFieldGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the geometric field for a field identified by an object for C.
  FUNCTION CMISSFieldGeometricFieldGetCPtr(FieldPtr,GeometricFieldPtr) BIND(C, NAME = "CMISSFieldGeometricFieldGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the geometric field to, for C.
    TYPE(C_PTR), INTENT(OUT) :: GeometricFieldPtr !<C pointer to the geometric field for the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSFieldType), POINTER :: GeometricField

    CMISSFieldGeometricFieldGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(GeometricFieldPtr)) THEN
          CALL C_F_POINTER(GeometricFieldPtr, GeometricField)
          IF(ASSOCIATED(GeometricField)) THEN
            CALL CMISSFieldGeometricFieldGet(Field, GeometricField, CMISSFieldGeometricFieldGetCPtr)
          ELSE
            CMISSFieldGeometricFieldGetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldGeometricFieldGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldGeometricFieldGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldGeometricFieldGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldGeometricFieldGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field identified by a user number for C.
  FUNCTION CMISSFieldGeometricFieldSetCNum(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber) BIND(C, &
  & NAME = "CMISSFieldGeometricFieldSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region corresponding to the field to set the geometric field to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !< The user number for the field to set the geometric field to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeometricFieldUserNumber !<The field geometric field user number to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldGeometricFieldSetNum(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber,&
    & CMISSFieldGeometricFieldSetCNum)

    RETURN

  END FUNCTION CMISSFieldGeometricFieldSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field identified by an object for C.
  FUNCTION CMISSFieldGeometricFieldSetCPtr(FieldPtr,GeometricFieldPtr) BIND(C, NAME = "CMISSFieldGeometricFieldSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the geometric field to, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeometricFieldPtr !<C pointer to the geometric field for the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSFieldType), POINTER :: GeometricField

    CMISSFieldGeometricFieldSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(GeometricFieldPtr)) THEN
          CALL C_F_POINTER(GeometricFieldPtr, GeometricField)
          IF(ASSOCIATED(GeometricField)) THEN
            CALL CMISSFieldGeometricFieldSet(Field, GeometricField, CMISSFieldGeometricFieldSetCPtr)
          ELSE
            CMISSFieldGeometricFieldSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldGeometricFieldSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldGeometricFieldSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
     CMISSFieldGeometricFieldSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldGeometricFieldSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field identified by a user number for C.
  FUNCTION CMISSFieldLabelGetCNum(RegionUserNumber,FieldUserNumber,LabelSize,Label) BIND(C, NAME = "CMISSFieldLabelGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label Size
    CHARACTER(LEN = 1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelGetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize - 1) :: FLabel

    CALL CMISSFieldLabelGet(RegionUserNumber,FieldUserNumber,FLabel,CMISSFieldLabelGetCNum)
    CALL CMISSF2CString(FLabel, Label)

    RETURN

  END FUNCTION CMISSFieldLabelGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field identified by an object for C.
  FUNCTION CMISSFieldLabelGetCPtr(FieldPtr,LabelSize,Label) BIND(C, NAME = "CMISSFieldLabelGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label size
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = LabelSize -1) :: FLabel

    CMISSFieldLabelGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldLabelGet(Field, FLabel, CMISSFieldLabelGetCPtr)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldLabelGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldLabelGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldLabelGetCPtr

  !!
  !================================================================================================================================
  !!

  !>Sets/changes the character string label for a field identified by a user number for C.
  FUNCTION CMISSFieldLabelSetCNum(RegionUserNumber,FieldUserNumber,LabelSize, Label) BIND(C, NAME = "CMISSFieldLabelSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size.
    CHARACTER(LEN = 1, KIND= C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelSetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSC2FString (Label, FLabel)
    CALL CMISSFieldLabelSet(RegionUserNumber, FieldUserNumber, FLabel,CMISSFieldLabelSetCNum)

    RETURN

  END FUNCTION CMISSFieldLabelSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field identified by an object for C.
  FUNCTION CMISSFieldLabelSetCPtr(FieldPtr,LabelSize, Label) BIND(C, NAME = "CMISSFieldLabelSet")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CMISSFieldLabelSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldLabelSet(Field,FLabel,CMISSFieldLabelSetCPtr)
      ELSE
        CMISSFieldLabelSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldLabelSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldLabelSetCPtr


  !
  !================================================================================================================================
  !

  !>Returns the mesh decomposition for a field identified by a user number for C.

  FUNCTION CMISSFieldMeshDecompositionGetCNum(RegionUserNumber,FieldUserNumber,DecompositionUserNumber) &
    & BIND(C, NAME = "CMISSFieldMeshDecompositionGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionUserNumber !<The field mesh decomposition user number for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldMeshDecompositionGet(RegionUserNumber,FieldUserNumber,DecompositionUserNumber,&
    & CMISSFieldMeshDecompositionGetCNum)

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the mesh decomposition for a field identified by an object for C.
  FUNCTION CMISSFieldMeshDecompositionGetCPtr(FieldPtr,MeshDecompositionPtr) BIND(C, NAME = "CMISSFieldMeshDecompositionGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the mesh decomposition for, for C.
    TYPE(C_PTR), INTENT(OUT) :: MeshDecompositionPtr !<C pointer to the field mesh decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSDecompositionType), POINTER :: MeshDecomposition

    CMISSFieldMeshDecompositionGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(MeshDecompositionPtr)) THEN
          CALL C_F_POINTER(MeshDecompositionPtr, MeshDecomposition)
          IF(ASSOCIATED(MeshDecomposition)) THEN
            CALL CMISSFieldMeshDecompositionGet(Field, MeshDecomposition, CMISSFieldMeshDecompositionGetCPtr)
          ELSE
            CMISSFieldMeshDecompositionGetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldMeshDecompositionGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldMeshDecompositionGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldMeshDecompositionGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field identified by a user number for C.
  FUNCTION CMISSFieldMeshDecompositionSetCNum(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber) &
  & BIND(C, NAME = "CMISSFieldMeshDecompositionSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number for the mesh to set the mesh decomposititon to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The field mesh decomposition user number for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldMeshDecompositionSet(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber,&
    & CMISSFieldMeshDecompositionSetCNum)

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field identified by an object for C.
  FUNCTION CMISSFieldMeshDecompositionSetCPtr(FieldPtr,MeshDecompositionPtr) BIND(C,NAME = "CMISSFieldMeshDecompositionSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the mesh decomposition to, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshDecompositionPtr !<C pointer to the mesh decomposition for the field to set.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSDecompositionType), POINTER :: MeshDecomposition

    CMISSFieldMeshDecompositionSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(MeshDecompositionPtr)) THEN
          CALL C_F_POINTER(MeshDecompositionPtr, MeshDecomposition)
          IF (ASSOCIATED(MeshDecomposition)) THEN
            CALL CMISSFieldMeshDecompositionSet(Field, MeshDecomposition, CMISSFieldMeshDecompositionSetCPtr)
          ELSE
            CMISSFieldMeshDecompositionSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldMeshDecompositionSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldMeshDecompositionSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldMeshDecompositionSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the number of componenets for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfComponentsGetCNum(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents) &
  & BIND(C, NAME = "CMISSFieldNumberOfComponentsGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsGetCNum
    !Local variables

    CALL CMISSFieldnumberOfComponentsGet(RegionUserNumber, FieldUserNumber, VariableType, NumberOfComponents, &
    & CMISSFieldNumberOfComponentsGetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the number of components for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfComponentsGetCPtr(FieldPtr,VariableType,NumberOfComponents) BIND(C, NAME = &
  & "CMISSFieldNumberOfComponentsGet")

    !Argument variables
    TYPE(C_PTR),VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the number of components for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents !<The number of components in the field variable for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfComponentsGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfComponentsGet(Field, VariableType, NumberOfComponents, CMISSFieldNumberOfComponentsGetCPtr)
      ELSE
        CMISSFieldNumberOfComponentsGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfComponentsGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsGetCPtr


  !
  !================================================================================================================================
  !

  !>Sets/changes the number of componenets for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfComponentsSetCNum(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents) &
  & BIND(C, NAME = "CMISSFieldNumberOfComponentsSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the number of components to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the field variable for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldNumberOfComponentsSet(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents, &
    & CMISSFieldNumberofComponentsSetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of components for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfComponentsSetCPtr(FieldPtr,VariableType,NumberOfComponents) BIND(C, NAME = &
  & "CMISSFieldNumberOfComponentsSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the number of components for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the field variables for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfComponentsSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfComponentsSet(Field, VariableType, NumberOfComponents, CMISSFieldNumberOfComponentsSetCPtr)
      ELSE
        CMISSFieldNumberOfComponentsSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfComponentsSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the number of variables for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfVariablesGetCNum(RegionUserNumber,FieldUserNumber,NumberOfVariables) BIND(C, NAME = &
  & "CMISSFieldNumberOfVariablesGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the number of variables for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the number of variables for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfVariables !<The number of variables in the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldNumberOfVariablesGet(RegionUserNumber, FieldUserNumber, NumberOfVariables, &
    & CMISSFieldNumberOfVariablesGetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the number of variables for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfVariablesGetCPtr(FieldPtr,NumberOfVariables) BIND(C, NAME = "CMISSFieldNumberOfVariablesGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the number of variables for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfVariables !<The number of variables in the field for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfVariablesGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfVariablesGet(Field, NumberOfVariables, CMISSFieldNumberOfVariablesGetCPtr)
      ELSE
        CMISSFieldNumberOfVariablesGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfVariablesGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfVariablesSetCNum(RegionUserNumber,FieldUserNumber,NumberOfVariables) BIND(C, &
  & NAME = "CMISSFieldNumberOfVariablesSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the number of variables to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the number of variables to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfVariables !<The number of variables set to the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldNumberofVariablesSet(RegionUserNumber, FieldUserNumber, NumberOfVariables, &
    & CMISSFieldNumberOfVariablesSetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfVariablesSetCPtr(FieldPtr,NumberOfVariables) BIND(C, NAME = "CMISSFieldNumberOfVariablesSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfVariables
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesSetCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfVariablesSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfVariablesSet(Field, NumberOfVariables, CMISSFieldNumberOfVariablesSetCPtr)
      ELSE
        CMISSFieldNumberOfVariablesSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfVariablesSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesSetCPtr


  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME= "CMISSFieldParameterSetAddConstantIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value, CMISSFieldParameterSetAddConstantIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantIntgCNum

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddConstantIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME= "CMISSFieldParameterSetAddConstantIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantIntgCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantIntgCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantIntg(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantIntgCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantIntgCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.

  FUNCTION CMISSFieldParameterSetAddConstantSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantSP (RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
    & Value,CMISSFieldParameterSetAddConstantSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantSPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddConstantSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldParameterSetAddConstantSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantSP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantSPCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantSPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
    & Value,CMISSFieldParameterSetAddConstantDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantDPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddConstantDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetAddConstantDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantDP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantDPCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantDPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetAddConstantL(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, ComponentNumber, &
    & Value, CMISSFieldParameterSetAddConstantLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantLCNum

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddConstantLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetAddConstantL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantL(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantLCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN


  END FUNCTION CMISSFieldParameterSetAddConstantLCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value, CMISSFieldParameterSetAddElementIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementIntgCNum

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an element in the given parameter set for field variable component for a field identified by an object, for C.
  FUNCTION CMISSFieldParameterSetAddElementIntgCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementIntgCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementIntg(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementIntgCPtr)
      ELSE
        CMISSFieldParameterSetAddElementIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementIntgCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an element in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddElementSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementSPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an element in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddElementSPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementSPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementSP(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementSPCPtr)
      ELSE
        CMISSFieldParameterSetAddElementSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementSPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementDPCNum


    !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an element in the given parameter set for field variable component for a field identified by an object, for C.
  FUNCTION CMISSFieldParameterSetAddElementDPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementDPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementDP(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementDPCPtr)
      ELSE
        CMISSFieldParameterSetAddElementDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementDPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementLCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementLCNum

    !
  !================================================================================================================================
  !

  !>Adds the given logical value to an element in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddElementLCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementLCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementL(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
         & CMISSFieldParameterSetAddElementLCPtr)
      ELSE
        CMISSFieldParameterSetAddElementLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementLCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an node in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddNodeIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber,&
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeIntgCNum

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeIntgCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeIntgCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeIntg(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeIntgCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeIntgCPtr



  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value  !<The single precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeSPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeSPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value  !<The single precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeSPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeSP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeSPCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeSPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value  !<The double precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeDPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an node in the given parameter set for field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddNodeDPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value  !<The double precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeDPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeDP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeDPCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeDPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeLCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetAddNodeLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeLCNum

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeLCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber,ComponentNumber, &
    & Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeLCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeL(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, ComponentNumber, &
        & Value, CMISSFieldParameterSetAddNodeLCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeLCPtr

  !
  !================================================================================================================================
  !

  !>Creates a new parameter set of type set type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetCreateCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, &
  & NAME = "CMISSFieldParameterSetCreateNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to create the parameter set on for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to create, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetCreateCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetCreate(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & CMISSFieldParameterSetCreateCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetCreateCNum

  !
  !================================================================================================================================
  !

  !>Creates a new parameter set of type set type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetCreateCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = "CMISSFieldParameterSetCreate")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to create the parameter set on for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to create, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetCreateCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetCreateCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetCreate(Field,VariableType,FieldSetType,CMISSFieldParameterSetCreateCPtr)
      ELSE
        CMISSFieldParameterSetCreateCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetCreateCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetCreateCPtr

  !
  !================================================================================================================================
  !

  !>Destroys the specified parameter set type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetDestroyCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, NAME= &
  & "CMISSFieldParameterSetDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to destroy the parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to destroy, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDestroyCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetDestroy(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & CMISSFieldParameterSetDestroyCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys the specified parameter set type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetDestroyCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = "CMISSFieldParameterSetDestroy")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to destroy the parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to destroy, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDestroyCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetDestroy(Field,VariableType,FieldSetType,CMISSFieldParameterSetDestroyCPtr)
      ELSE
        CMISSFieldParameterSetDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local integer data array for a field identified by an user number, for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersPtr, &
  & ParametersSize) BIND(C, NAME = "CMISSFieldParameterSetDataGetIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get for C. \see OPENCMISS_FieldParameterSetTypes
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetIntgCNum
    !Local variables
    INTEGER(C_INT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetIntgCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetIntgCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataGetIntg(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType, Parameters, &
        & CMISSFieldParameterSetDataGetIntgCNum)
      ParametersSize = Size(Parameters,1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local integer data array for a field identified by an object for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetIntgCPtr(FieldPtr,VariableType,FieldSetType,ParametersPtr,ParametersSize) BIND(C, NAME = &
  & "CMISSFieldParameterSetDataGetIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetIntgCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    INTEGER(C_INT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetIntgCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr, Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataGetIntg(Field, VariableType, FieldSetType, Parameters, &
            & CMISSFieldParameterSetDataGetIntgCPtr)
          ParametersSize = Size(Parameters,1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
          ELSE
          CMISSFieldParameterSetDataGetIntgCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataGetIntgCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local single precision data array for a field identified by an user number for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
  &ParametersPtr,ParametersSize) BIND(C, NAME = "CMISSFieldParameterSetDataGetSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetSPCNum !<Error code.
    !Local variables
    INTEGER(C_INT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetSPCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetSPCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataGetSP(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
      & CMISSFieldParameterSetDataGetSPCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetSPCNum

  !
  !================================================================================================================================
  !
  !>Returns a pointer to the specified field parameter set local single precision data array for a field identified by an object. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.

  FUNCTION CMISSFieldParameterSetDataGetSPCPtr(FieldPtr,VariableType,FieldSetType,ParametersPtr,ParametersSize) BIND(C, NAME = &
  & "CMISSFieldParameterSetDataGetSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetSPCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    REAL(C_DOUBLE), POINTER :: Parameters (:)

    CMISSFieldParameterSetDataGetSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetSPCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataGetSP(Field, VariableType, FieldSetType, Parameters, CMISSFieldParameterSetDataGetSPCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataGetSPCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataGetSPCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetSPCPtr



!missing code


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,Value, &
      & CMISSFieldParameterSetGetConstantIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantIntgCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantIntg(Field,VariableType,FieldSetType,ComponentNumber,Value, &
          & CMISSFieldParameterSetGetConstantIntgCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetGetConstantSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantSPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantSPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantSP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetConstantSPCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantSPCPtr


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetGetConstantDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantDPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantDPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantDP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetConstantDPCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantDPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetGetConstantLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantLCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantL(Field, VariableType, FieldSetType, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetConstantLCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantLCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementIntgCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementIntg(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetElementIntgCPtr)
      ELSE
        CMISSFieldParameterSetGetElementIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementSPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified element of a field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetGetElementSPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementSP(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetElementSPCPtr)
      ELSE
        CMISSFieldParameterSetGetElementSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementSPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementDPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementDPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementDP(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetElementDPCPtr)
      ELSE
        CMISSFieldParameterSetGetElementDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementDPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementLCNum


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementLCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementL(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value,  &
          & CMISSFieldParameterSetGetElementLCPtr)
      ELSE
        CMISSFieldParameterSetGetElementLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementLCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeIntg(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeIntgCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeIntgCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeIntg(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeIntgCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_FLOAT),INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeSP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      &UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeSPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeSPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeSPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeSP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeSPCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeSPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_DOUBLE),INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeDP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeDPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeDPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeDPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeDP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeDPCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeDPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    LOGICAL(C_BOOL),INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeL(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeLCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeLCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeL(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeLCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeLCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantIntg(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateConstantIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantIntgCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetUpdateConstantIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateConstantIntg")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantIntg(Field, VariableType,FieldSetType,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateConstantIntgCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantIntgCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantSP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
      & Value,CMISSFieldParameterSetUpdateConstantSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantSPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
    & NAME = "CMISSFieldParameterSetUpdateConstantSP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantSP(Field, VariableType,FieldSetType,ComponentNumber,Value, &
          & CMISSFieldParameterSetUpdateConstantSPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantSPCPtr
  
  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantDP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateConstantDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantDPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
    & NAME = "CMISSFieldParameterSetUpdateConstantDP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantDP(Field, VariableType,FieldSetType,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateConstantDPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantDPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantL(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetUpdateConstantLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantLCNum


  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateConstantL")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantL(Field, VariableType,FieldSetType,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateConstantLCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantLCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    INTEGER(C_INT), INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementIntg(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
      & UserElementNumber,ComponentNumber,Value,CMISSFieldParameterSetUpdateElementIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementIntgCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementIntgCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementIntg")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    INTEGER(C_INT), INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementIntg(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateElementIntgCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementIntgCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementSP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementSPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementSPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementSP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementSP(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value, &
          & CMISSFieldParameterSetUpdateElementSPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementSPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementDP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementDPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementDPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementDP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementDP(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateElementDPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementDPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementL(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementLCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementLCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementL")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementL(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value, &
          & CMISSFieldParameterSetUpdateElementLCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementLCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the parameter set update for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateFinishCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to finish the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to finish the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to finish the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to finish the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateFinishCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateFinish(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
      & CMISSFieldParameterSetUpdateFinishCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the parameter set update for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateFinishCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to finish the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to finish the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateFinishCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateFinish(Field,VariableType,FieldSetType,CMISSFieldParameterSetUpdateFinishCPtr)
      ELSE
        CMISSFieldParameterSetUpdateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeIntgCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeIntgCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeIntg(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber,Value,CMISSFieldParameterSetUpdateNodeIntgCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeIntgCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeSPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the node and derivative of the field variable component for a field identified by an object, for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeSPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeSP(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeSPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeSPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeDPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeDPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeDP(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeDPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeDPCPtr

  !
  !================================================================================================================================
  !,

  !>Updates the given parameter set with the given logical value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeLCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber,&
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeLCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeLCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeL(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeLCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeLCPtr

  !
  !================================================================================================================================
  !

  !>Starts the parameter set update for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateStartCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, &
    & NAME = "CMISSFieldParameterSetUpdateStartNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber!<The user number of the region containing the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to start the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to start the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateStartCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateStart(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,&
      & CMISSFieldParameterSetUpdateStartCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the parameter set update for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateStartCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateStart")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to start the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to start the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateStart(Field,VariableType,FieldSetType,CMISSFieldParameterSetUpdateStartCPtr)
      ELSE
        CMISSFieldParameterSetUpdateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateStartCPtr

  !
  !================================================================================================================================
  !

  !>Returns the scaling type for a field identified by a user number for C.
  FUNCTION CMISSFieldScalingTypeGetCNum(RegionUserNumber,FieldUserNumber,ScalingType) BIND(C, NAME = &
    & "CMISSFieldScalingTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the scaling type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the scaling type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ScalingType !<The field scaling type to get, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldScalingTypeGet(RegionUserNumber,FieldUserNumber,ScalingType,CMISSFieldScalingTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldScalingTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the scaling type for a field identified by an object for C.
  FUNCTION CMISSFieldScalingTypeGetCPtr(FieldPtr,ScalingType) BIND(C, NAME = "CMISSFieldScalingTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the scaling type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ScalingType !<The field scaling type to get, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeGetCPtr
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldScalingTypeGetCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldScalingTypeGet(Field, ScalingType, CMISSFieldScalingTypeGetCPtr)
      ELSE
        CMISSFieldScalingTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldScalingTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldScalingTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field identified by a user number for C.
  FUNCTION CMISSFieldScalingTypeSetCNum(RegionUserNumber,FieldUserNumber,ScalingType) BIND(C, NAME = &
    & "CMISSFieldScalingTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ScalingType !<The field scaling type to set, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldScalingTypeSet(RegionUserNumber,FieldUserNumber,ScalingType,CMISSFieldScalingTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldScalingTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field identified by an object for C.
  FUNCTION CMISSFieldScalingTypeSetCPtr(FieldPtr,ScalingType) BIND(C, NAME = "CMISSFieldScalingTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ScalingType !<The field scaling type to set, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldScalingTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldScalingTypeSet(Field, ScalingType, CMISSFieldScalingTypeSetCPtr)
      ELSE
        CMISSFieldScalingTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldScalingTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldScalingTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the field type for a field identified by a user number for C.
  FUNCTION CMISSFieldTypeGetCNum(RegionUserNumber,FieldUserNumber,FieldType) BIND(C, NAME = "CMISSFieldTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the field type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the field type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: FieldType !<The field type to get, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldTypeGet(RegionUserNumber,FieldUserNumber,FieldType,CMISSFieldTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the type for a field identified by an object for C.
  FUNCTION CMISSFieldTypeGetCPtr(FieldPtr,FieldType) BIND(C, NAME = "CMISSFieldTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: FieldType !<The field type to get, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldTypeGet(Field, FieldType, CMISSFieldTypeGetCPtr)
      ELSE
        CMISSFieldTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the field type for a field identified by a user number for C.
  FUNCTION CMISSFieldTypeSetCNum(RegionUserNumber,FieldUserNumber,FieldType) BIND(C, NAME = "CMISSFieldTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to set, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldTypeSet(RegionUserNumber,FieldUserNumber,FieldType,CMISSFieldTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a field identified by an object for C.
  FUNCTION CMISSFieldTypeSetCPtr(FieldPtr,FieldType) BIND(C, NAME = "CMISSFieldTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to set, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldTypeSet(Field, FieldType, CMISSFieldTypeSetCPtr)
      ELSE
        CMISSFieldTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable for a field identified by a user number.
  FUNCTION CMISSFieldVariableLabelGetCCNum(RegionUserNumber,FieldUserNumber,VariableType,LabelSize,Label) BIND(C, &
    & NAME = "CMISSFieldVariableLabelGetCNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelGetCCNum !<Error Code.
    !Local variable
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSFieldVariableLabelGetCNum(RegionUserNumber,FieldUserNumber,VariableType,FLabel,CMISSFieldVariableLabelGetCCNum)
    CALL CMISSF2CString(FLabel,Label)

    RETURN

  END FUNCTION CMISSFieldVariableLabelGetCCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable for a field identified by an object.
  FUNCTION CMISSFieldVariableLabelGetCPtrC(FieldPtr,VariableType,LabelSize,Label) BIND(C, NAME = "CMISSFieldVariableLabelGetCObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelGetCPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = LabelSize - 1) :: FLabel


    CMISSFieldVariableLabelGetCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldVariableLabelGetCObj(Field, VariableType,FLabel, CMISSFieldVariableLabelGetCPtrC)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldVariableLabelGetCPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableLabelGetCPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableLabelGetCPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldVariableLabelSetCNumberC(RegionUserNumber,FieldUserNumber,VariableType,LabelSize,Label) BIND(C, &
    & NAME ="CMISSFieldVariableLabelSetCNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to set the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN = 1, KIND = C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelSetCNumberC !<Error Code.
    !Local variable
    CHARACTER(LEN = LabelSize -1) :: FLabel

    CALL CMISSC2FLabel(Label, FLabel)
    CALL CMISSFieldVariableLabelSetCNumber(RegionUserNumber,FieldUserNumber,VariableType,FLabel,CMISSFieldVariableLabelSetCNumberC)

    RETURN

  END FUNCTION CMISSFieldVariableLabelSetCNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldVariableLabelSetCPtrC(FieldPtr,VariableType,LabelSize,Label) BIND(C, NAME = "CMISSFieldVariableLabelSetCObj")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelSetCPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = LabelSize - 1) :: FLabel


    CMISSFieldVariableLabelSetCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldVariableLabelSetCObj(Field, VariableType,FLabel, CMISSFieldVariableLabelSetCPtrC)
      ELSE
        CMISSFieldVariableLabelSetCPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableLabelSetCPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableLabelSetCPtrC

  !
  !================================================================================================================================
  !

  !>Returns the field variable types for a field identified by a user number for C.
  FUNCTION CMISSFieldVariableTypesGetNumberC(RegionUserNumber,FieldUserNumber,VariableTypes) BIND(C, NAME = &
    & "CMISSFieldVariableTypesGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the field variable types for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the field variable types for, for C.
    INTEGER(C_INT), INTENT(OUT) :: VariableTypes !<VariableTypes(variable_idx). The field variable types for the variable_idx'th field variable, for C. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesGetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldVariableTypesGetNumber(RegionUserNumber,FieldUserNumber,VariableTypes,CMISSFieldVariableTypesGetNumberC)

    RETURN

  END FUNCTION CMISSFieldVariableTypesGetNumberC

    !
  !================================================================================================================================
  !

  !>Returns the variable types for a field identified by an object for C.
  FUNCTION CMISSFieldVariableTypesGetPtrC(FieldPtr,VariableTypes) BIND(C, NAME = "CMISSFieldVariableTypesGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field variable types for, for C.
    INTEGER(C_INT), INTENT(OUT) :: VariableTypes !<VariableTypes(variable_idx). The field variable types for the variable_idx'th field variable, for C. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldVariableTypesGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldVariableTypesGetObj(Field, VariableTypes, CMISSFieldVariableTypesGetPtrC)
      ELSE
        CMISSFieldVariableTypesGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableTypesGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableTypesGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the field variable types for a field identified by a user number for C.
  FUNCTION CMISSFieldVariableTypesSetNumberC(RegionUserNumber,FieldUserNumber,VariableTypes) BIND(C, NAME = &
    & "CMISSFieldVariableTypesSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the field variable types to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the field variable types to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableTypes !<VariableTypes(variable_idx). The field variable types for the variable_idx'th field variable, for C. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldVariableTypesSetNumber(RegionUserNumber,FieldUserNumber,VariableTypes,CMISSFieldVariableTypesSetNumberC)

    RETURN

  END FUNCTION CMISSFieldVariableTypesSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the variable types for a field identified by an object for C.
  FUNCTION CMISSFieldVariableTypesSetPtrC(FieldPtr,VariableTypes) BIND(C, NAME = "CMISSFieldVariableTypesSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the field variable types to, for C.
    INTEGER(C_INT), INTENT(OUT) :: VariableTypes !<VariableTypes(variable_idx). The field variable types for the variable_idx'th field variable, for C. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldVariableTypesSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldVariableTypesSetObj(Field, VariableTypes, CMISSFieldVariableTypesSetPtrC)
      ELSE
        CMISSFieldVariableTypesSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableTypesSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableTypesSetPtrC

!!==================================================================================================================================
!!
!! FIELD_IO_ROUTINES
!!
!!==================================================================================================================================

  !>Export element information for fields set identified by an object for C. \todo number method
  FUNCTION CMISSFieldIOElementsExportCCPtrC(FieldsPtr,FileNameSize,FileName,MethodSize,Method) BIND(C, NAME = &
    & "CMISSFieldIOElementsExportCCObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldsPtr !<The fields to export the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !< Size of the file name to export the elements to for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<The file name to export the elements to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MethodSize !<Size of the export method name for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN):: Method(MethodSize) !<The export method to use for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldIOElementsExportCCPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldsType), POINTER :: Fields
    CHARACTER(LEN = FileNameSize -1 ) :: FFileName
    CHARACTER(LEN = MethodSize - 1) :: FMethod

    CMISSFieldIOElementsExportCCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldsPtr)) THEN
      CALL C_F_POINTER(FieldsPtr, Fields)
      IF(ASSOCIATED(Fields)) THEN
        CALL CMISSC2FString(FileName, FFileName)
        CALL CMISSC2FString(Method, FMethod)
        CALL CMISSFieldIOElementsExportCCObj(Fields,FileName,Method, CMISSFieldIOElementsExportCCPtrC)
      ELSE
        CMISSFieldIOElementsExportCCPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldIOElementsExportCCPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldIOElementsExportCCPtrC

  !
  !================================================================================================================================
  !

  !>Export nodal information for fields set identified by an object for C. \todo number method
  FUNCTION CMISSFieldIONodesExportCCPtrC(FieldsPtr,FileNameSize,FileName,MethodSize,Method) BIND(C, NAME = &
    &  "CMISSFieldIONodesExportCCObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldsPtr !<The fields to export the nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !< Size of the file name to export the nodes to for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<The file name to export the nodes to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MethodSize !<Size of the export method name for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN):: Method(MethodSize) !<The export method to use for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldIONodesExportCCPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldsType), POINTER :: Fields
    CHARACTER(LEN = FileNameSize -1 ) :: FFileName
    CHARACTER(LEN = MethodSize - 1) :: FMethod

    CMISSFieldIONodesExportCCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldsPtr)) THEN
      CALL C_F_POINTER(FieldsPtr, Fields)
      IF(ASSOCIATED(Fields)) THEN
        CALL CMISSC2FString(FileName, FFileName)
        CALL CMISSC2FString(Method, FMethod)
        CALL CMISSFieldIONodesExportCCObj(Fields,FileName,Method, CMISSFieldIONodesExportCCPtrC)
      ELSE
        CMISSFieldIONodesExportCCPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldIONodesExportCCPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldIONodesExportCCPtrC

!!==================================================================================================================================
!!
!! GENERATED_MESH_ROUTINES
!!
!!==================================================================================================================================

  !>Returns the basis for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshBasisGetNumberC(GeneratedMeshUserNumber,BasisUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshBasisGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to get the basis for, for C.
    INTEGER(C_INT), INTENT(OUT) :: BasisUserNumber !<The user number of the basis to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisGetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshBasisGetNumber(GeneratedMeshUserNumber,BasisUserNumber,CMISSGeneratedMeshBasisGetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisGetNumberC

    !
  !================================================================================================================================
  !

  !>Returns the basis for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshBasisGetPtrC(GeneratedMeshPtr,BasisPtr) BIND(C, NAME = "CMISSGeneratedMeshBasisGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr!<C pointer to the generated mesh to get the basis for.
    TYPE(C_PTR), INTENT(INOUT) :: BasisPtr !<C pointer to the basis to get.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSBasisType), POINTER :: Basis

    CMISSGeneratedMeshBasisGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSGeneratedMeshBasisGetObj(GeneratedMesh, Basis, CMISSGeneratedMeshBasisGetPtrC)
          ELSE
            CMISSGeneratedMeshBasisGetPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshBasisGetPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshBasisGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshBasisGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshBasisSetNumberC(GeneratedMeshUserNumber,BasisUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshBasisSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the basis to, for C.
    INTEGER(C_INT), INTENT(IN) :: BasisUserNumber !<The user number of the basis to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisSetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshBasisSetNumber(GeneratedMeshUserNumber,BasisUserNumber,CMISSGeneratedMeshBasisSetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshBasisSetPtrC(GeneratedMeshPtr,BasisPtr) BIND(C, NAME = "CMISSGeneratedMeshBasisSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr!<C pointer to the generated mesh to set the basis to.
    TYPE(C_PTR), INTENT(INOUT) :: BasisPtr !<C pointer to the basis to set.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSBasisType), POINTER :: Basis

    CMISSGeneratedMeshBasisSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSGeneratedMeshBasisSetObj(GeneratedMesh, Basis, CMISSGeneratedMeshBasisSetPtrC)
          ELSE
            CMISSGeneratedMeshBasisSetPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshBasisSetPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshBasisSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshBasisSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisSetPtrC

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshCreateFinishNumberC(GeneratedMeshUserNumber,MeshUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateFinishNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the basis to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to generate, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateFinishNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshCreateFinishNumber(GeneratedMeshUserNumber,MeshUserNumber,CMISSGeneratedMeshCreateFinishNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateFinishNumberC

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a generated mesh identified by an object.
  FUNCTION CMISSGeneratedMeshCreateFinishPtrC(GeneratedMeshPtr,MeshUserNumber,MeshPtr) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateFinishObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr!<C pointer to the generated mesh to finish the creation of.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to generate, for C.
    TYPE(C_PTR), INTENT(INOUT) :: MeshPtr !<C pointer to the generated mesh.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateFinishPtrC !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSGeneratedMeshCreateFinishPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(MeshPtr)) THEN
          CALL C_F_POINTER(MeshPtr, Mesh)
          IF(ASSOCIATED(Mesh)) THEN
            CALL CMISSGeneratedMeshCreateFinishObj(GeneratedMesh, MeshUserNumber, Mesh, CMISSGeneratedMeshCreateFinishPtrC)
          ELSE
            CMISSGeneratedMeshCreateFinishPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshCreateFinishPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshCreateFinishPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshCreateFinishPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateFinishPtrC

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshCreateStartNumberC(GeneratedMeshUserNumber,RegionUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateStartNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to create, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region to create the generated mesh in, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateStartNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshCreateStartNumber(GeneratedMeshUserNumber,RegionUserNumber,CMISSGeneratedMeshCreateStartNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateStartNumberC

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshCreateStartPtrC(GeneratedMeshUserNumber,RegionPtr,GeneratedMeshPtr) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateStartObj")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to create, for C.
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshPtr !<C pointer to the generated mesh to finish the creation of.
    TYPE(C_PTR), INTENT(INOUT) :: RegionPtr !<C pointer to the region to created generated mesh in.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateStartPtrC !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSregionType), POINTER :: Region

    CMISSGeneratedMeshCreateStartPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(RegionPtr)) THEN
          CALL C_F_POINTER(RegionPtr, Region)
          IF(ASSOCIATED(Region)) THEN
            CALL CMISSGeneratedMeshCreateStartObj(GeneratedMeshUserNumber, Region, GeneratedMesh, CMISSGeneratedMeshCreateStartPtrC)
          ELSE
            CMISSGeneratedMeshCreateStartPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshCreateStartPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshCreateStartPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshCreateStartPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateStartPtrC

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshDestroyNumberC(GeneratedMeshUserNumber) BIND(C, NAME = "CMISSGeneratedMeshDestroyNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to destroy, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshDestroyNumberC !<Error Code.
    !Local variableC

    CALL CMISSGeneratedMeshDestroyNumber(GeneratedMeshUserNumber,CMISSGeneratedMeshDestroyNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshDestroyNumberC

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshDestroyPtrC(GeneratedMeshPtr) BIND(C, NAME = "CMISSGeneratedMeshDestroyObj")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshPtr !<C pointer to the generated mesh to destroy.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshDestroyPtrC !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshDestroyPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshDestroyObj(GeneratedMesh, CMISSGeneratedMeshDestroyPtrC)
      ELSE
        CMISSGeneratedMeshDestroyPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshDestroyPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshDestroyPtrC


  !
  !================================================================================================================================
  !

  !>Returns the type of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshTypeGetNumberC(GeneratedMeshUserNumber,GeneratedMeshType) BIND(C, NAME = &
    & "CMISSGeneratedMeshTypeGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to get the type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: GeneratedMeshType !<On return, the type of the generated mesh to get, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeGetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshTypeGetNumber(GeneratedMeshUserNumber, GeneratedMeshType, CMISSGeneratedMeshTypeGetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the type of a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshTypeGetPtrC(GeneratedMeshPtr,GeneratedMeshType) BIND(C, NAME = "CMISSGeneratedMeshTypeGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to get the type for.
    INTEGER(C_INT), INTENT(OUT) :: GeneratedMeshType !<On return, the type of the generated mesh to get, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeGetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshTypeGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshTypeGetObj(GeneratedMesh, GeneratedMeshType, CMISSGeneratedMeshTypeGetPtrC)
      ELSE
        CMISSGeneratedMeshTypeGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshTypeGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshTypeSetNumberC(GeneratedMeshUserNumber,GeneratedMeshType) BIND(C, NAME = &
    & "CMISSGeneratedMeshTypeSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshType !<On return, the type of the generated mesh to set, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeSetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshTypeSetNumber(GeneratedMeshUserNumber, GeneratedMeshType, CMISSGeneratedMeshTypeSetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh identified by an object.
  FUNCTION CMISSGeneratedMeshTypeSetPtrC(GeneratedMeshPtr,GeneratedMeshType) BIND(C, NAME = "CMISSGeneratedMeshTypeSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to set the type to.
    INTEGER(C_INT), INTENT(OUT) :: GeneratedMeshType !<On return, the type of the generated mesh to set, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshTypeSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshTypeSetObj(GeneratedMesh, GeneratedMeshType, CMISSGeneratedMeshTypeSetPtrC)
      ELSE
        CMISSGeneratedMeshTypeSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshTypeSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeSetPtrC

  !
  !================================================================================================================================
  !

  !>Calculates and sets the geometric field parameters for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshGeometricParametersCalculateNumberC(RegionUserNumber,FieldUserNumber,GeneratedMeshUserNumber) &
    &  BIND(C, NAME = "CMISSGeneratedMeshGeometricParametersCalculateNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to calculate the geometric parameters for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to calculate the geometric parameters for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to calculate the geometric parameters for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshGeometricParametersCalculateNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshGeometricParametersCalculateNumber(RegionUserNumber,FieldUserNumber,GeneratedMeshUserNumber, &
      & CMISSGeneratedMeshGeometricParametersCalculateNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshGeometricParametersCalculateNumberC

  !
  !================================================================================================================================
  !

  !>Calculates and sets the geometric field parameters for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshGeometricParametersCalculatePtrC(FieldPtr,GeneratedMeshPtr)  BIND(C, NAME = &
    & "CMISSGeneratedMeshGeometricParametersCalculateObj")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldPtr !<The field to calculate the geometric parameters for, for C.
    TYPE(C_PTR), INTENT(IN) :: GeneratedMeshPtr !<The generated mesh to calculate the geometric parameters for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshGeometricParametersCalculatePtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshGeometricParametersCalculatePtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
          CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
          IF(ASSOCIATED(GeneratedMesh)) THEN
            CALL CMISSGeneratedMeshGeometricParametersCalculateObj(Field,GeneratedMesh, &
              & CMISSGeneratedMeshGeometricParametersCalculatePtrC)
          ELSE
            CMISSGeneratedMeshGeometricParametersCalculatePtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshGeometricParametersCalculatePtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshGeometricParametersCalculatePtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshGeometricParametersCalculatePtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshGeometricParametersCalculatePtrC

!!==================================================================================================================================
!!
!! MESH_ROUTINES
!!
!!==================================================================================================================================

  !>Finishes the creation of a domain decomposition for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionCreateFinishNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber) BIND(C, NAME = &
    & "CMISSDecompositionCreateFinishNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to finish the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to finish the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to finish for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateFinishNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionCreateFinishNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
      & CMISSDecompositionCreateFinishNumberC)

    RETURN

  END FUNCTION CMISSDecompositionCreateFinishNumberC

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionCreateFinishPtrC(DecompositionPtr) BIND(C, NAME = "CMISSDecompositionCreateFinishObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to finish creating.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateFinishPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionCreateFinishPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionCreateFinishObj(Decomposition, CMISSDecompositionCreateFinishPtrC)
      ELSE
        CMISSDecompositionCreateFinishPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionCreateFinishPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionCreateFinishPtrC

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionCreateStartNumberC(DecompositionUserNumber,RegionUserNumber,MeshUserNumber)  BIND(C, NAME = &
    & "CMISSDecompositionCreateStartNumberC")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to create for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to create the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to create the decomposition for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateStartNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionCreateStartNumber(DecompositionUserNumber,RegionUserNumber,MeshUserNumber, &
      & CMISSDecompositionCreateStartNumberC)

    RETURN

  END FUNCTION CMISSDecompositionCreateStartNumberC

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a decomposition identified by an object.
  FUNCTION CMISSDecompositionCreateStartPtrC(DecompositionUserNumber,MeshPtr,DecompositionPtr)  BIND(C, NAME = &
    & "CMISSDecompositionCreateStartObj")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to finish for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to create the decomposition for.
    TYPE(C_PTR), INTENT(OUT) :: DecompositionPtr !<C pointer to the decomposition to create.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateStartPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSDecompositionCreateStartPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        IF(C_ASSOCIATED(DecompositionPtr)) THEN
          CALL C_F_POINTER(DecompositionPtr, Decomposition)
          IF(ASSOCIATED(Decomposition)) THEN
            CALL CMISSDecompositionCreateStartObj(DecompositionUserNumber, Mesh, Decomposition, CMISSDecompositionCreateStartPtrC)
          ELSE
            CMISSDecompositionCreateStartPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSDecompositionCreateStartPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSDecompositionCreateStartPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionCreateStartPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionCreateStartPtrC

  !
  !================================================================================================================================
  !

  !>Destroys a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionDestroyNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber)  BIND(C, NAME = &
    & "CMISSDecompositionDestroyNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to destroy the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to destroy the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionDestroyNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionDestroyNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, CMISSDecompositionDestroyNumberC)

    RETURN

  END FUNCTION CMISSDecompositionDestroyNumberC

  !
  !================================================================================================================================
  !

  !>Destroys a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionDestroyPtrC(DecompositionPtr) BIND(C, NAME = "CMISSDecompositionDestroyObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to destroy.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionDestroyPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionDestroyPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionDestroyObj(Decomposition, CMISSDecompositionDestroyPtrC)
      ELSE
        CMISSDecompositionDestroyPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionDestroyPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionDestroyPtrC

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionElementDomainCalculateNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber) BIND(C, NAME =&
    & "CMISSDecompositionElementDomainCalculateNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to calculate the element domains for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to calculate the element domains, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to calculate the element domains for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainCalculateNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionElementDomainCalculateNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
      & CMISSDecompositionElementDomainCalculateNumberC)

    RETURN

  END FUNCTION CMISSDecompositionElementDomainCalculateNumberC

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionElementDomainCalculatePtrC(DecompositionPtr) BIND(C, NAME = &
    & "CMISSDecompositionElementDomainCalculateObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to calculate the element domains.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainCalculatePtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionElementDomainCalculatePtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionElementDomainCalculateObj(Decomposition, CMISSDecompositionElementDomainCalculatePtrC)
      ELSE
        CMISSDecompositionElementDomainCalculatePtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionElementDomainCalculatePtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionElementDomainCalculatePtrC

  !
  !================================================================================================================================
  !

  !>Returns the domain for a given element in a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionElementDomainGetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & ElementUserNumber,Domain) BIND(C, NAME = "CMISSDecompositionElementDomainGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the element domain for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the element domain for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the element domain for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to get the domain for, for C.
    INTEGER(C_INT), INTENT(OUT) :: Domain !<The computational domain of the element, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainGetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionElementDomainGetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,ElementUserNumber,Domain,&
      & CMISSDecompositionElementDomainGetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionElementDomainGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the domain for a given element in a decomposition identified by an object.
  FUNCTION CMISSDecompositionElementDomainGetPtrC(DecompositionPtr,ElementUserNumber,Domain)  BIND(C, NAME = &
    & "CMISSDecompositionElementDomainGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the domain for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to get the domain for, for C.
    INTEGER(C_INT), INTENT(OUT) :: Domain !<The computational domain of the element, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionElementDomainGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionElementDomainGetObj(Decomposition,ElementUserNumber,Domain,CMISSDecompositionElementDomainGetPtrC)
      ELSE
        CMISSDecompositionElementDomainGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionElementDomainGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionElementDomainGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the domain for a given element in a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionElementDomainSetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & ElementUserNumber,Domain) BIND(C, NAME = "CMISSDecompositionElementDomainSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the element domain to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the element domain to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the element domain to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to set the domain to, for C.
    INTEGER(C_INT), INTENT(IN) :: Domain !<The computational domain of the element to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainSetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionElementDomainSetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,ElementUserNumber,Domain,&
      & CMISSDecompositionElementDomainSetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionElementDomainSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the domain for a given element in a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionElementDomainSetPtrC(DecompositionPtr,ElementUserNumber,Domain)  BIND(C, NAME = &
    & "CMISSDecompositionElementDomainSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the domain to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to set the domain to, for C.
    INTEGER(C_INT), INTENT(IN) :: Domain !<The computational domain of the element to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionElementDomainSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionElementDomainSetObj(Decomposition,ElementUserNumber,Domain,CMISSDecompositionElementDomainSetPtrC)
      ELSE
        CMISSDecompositionElementDomainSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionElementDomainSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionElementDomainSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number used for the decomposition of a mesh for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionMeshComponentGetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & MeshComponentNumber) BIND(C, NAME = "CMISSDecompositionMeshComponentGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the mesh component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the mesh component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the mesh component for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponentNumber !<The mesh component number for the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentGetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionMeshComponentGetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,MeshComponentNumber,&
      & CMISSDecompositionMeshComponentGetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number used for the decomposition of a mesh for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionMeshComponentGetPtrC(DecompositionPtr,MeshComponentNumber) BIND(C, NAME = &
    & "CMISSDecompositionMeshComponentGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the mesh component for.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponentNumber !<The mesh component number for the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionMeshComponentGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionMeshComponentGetObj(Decomposition,MeshComponentNumber,CMISSDecompositionMeshComponentGetPtrC)
      ELSE
        CMISSDecompositionMeshComponentGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionMeshComponentGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number used for the decomposition of a mesh for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionMeshComponentSetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & MeshComponentNumber) BIND(C, NAME = "CMISSDecompositionMeshComponentSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the mesh component to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the mesh component to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the mesh component to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number for the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentSetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionMeshComponentSetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,MeshComponentNumber,&
      & CMISSDecompositionMeshComponentSetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number used for the decomposition of a mesh for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionMeshComponentSetPtrC(DecompositionPtr,MeshComponentNumber) BIND(C, NAME = &
    & "CMISSDecompositionMeshComponentSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the mesh component to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number for the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionMeshComponentSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionMeshComponentSetObj(Decomposition,MeshComponentNumber,CMISSDecompositionMeshComponentSetPtrC)
      ELSE
        CMISSDecompositionMeshComponentSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionMeshComponentSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the number of domains for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionNumberOfDomainsGetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & NumberOfDomains) BIND(C, NAME = "CMISSDecompositionNumberOfDomainsGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the number of domains for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the number of domains for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the number of domains for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfDomains !<The number of domains in the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsGetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionNumberOfDomainsGetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NumberOfDomains,&
      & CMISSDecompositionNumberOfDomainsGetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the number of domains for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionNumberOfDomainsGetPtrC(DecompositionPtr,NumberOfDomains) BIND(C, NAME = &
    & "CMISSDecompositionNumberOfDomainsGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the number of domains for.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfDomains !<The number of domains in the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionNumberOfDomainsGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionNumberOfDomainsGetObj(Decomposition,NumberOfDomains,CMISSDecompositionNumberOfDomainsGetPtrC)
      ELSE
        CMISSDecompositionNumberOfDomainsGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionNumberOfDomainsGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsGetPtrC

    !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionNumberOfDomainsSetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & NumberOfDomains) BIND(C, NAME = "CMISSDecompositionNumberOfDomainsSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the number of domains to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the number of domains to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the number of domains to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDomains !<The number of domains in the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsSetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionNumberOfDomainsSetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NumberOfDomains,&
      & CMISSDecompositionNumberOfDomainsSetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionNumberOfDomainsSetPtrC(DecompositionPtr,NumberOfDomains) BIND(C, NAME = &
    & "CMISSDecompositionNumberOfDomainsSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the number of domains to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDomains !<The number of domains in the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionNumberOfDomainsSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionNumberOfDomainsSetObj(Decomposition,NumberOfDomains,CMISSDecompositionNumberOfDomainsSetPtrC)
      ELSE
        CMISSDecompositionNumberOfDomainsSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionNumberOfDomainsSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the type of a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionTypeGetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType) BIND(C, NAME&
    & = "CMISSDecompositionTypeGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the decomposition type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the decomposition type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the decomposition type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionType !<The type of the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeGetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionTypeGetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType,&
      & CMISSDecompositionTypeGetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionTypeGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the type of a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionTypeGetPtrC(DecompositionPtr,DecompositionType) BIND(C, NAME = "CMISSDecompositionTypeGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the decomposition type for.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionType !<The type of the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionTypeGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionTypeGetObj(Decomposition,DecompositionType,CMISSDecompositionTypeGetPtrC)
      ELSE
        CMISSDecompositionTypeGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionTypeGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionTypeSetNumberC(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType) BIND(C, NAME&
    & = "CMISSDecompositionTypeSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the decomposition type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the decomposition type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the decomposition type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionType !<The type of the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeSetNumberC !<Error Code.
    !Local variable

    CALL CMISSDecompositionTypeSetNumber(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType,&
      &CMISSDecompositionTypeSetNumberC)

    RETURN

  END FUNCTION CMISSDecompositionTypeSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionTypeSetPtrC(DecompositionPtr,DecompositionType) BIND(C, NAME = "CMISSDecompositionTypeSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the decomposition type to.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionType !<The type of the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionTypeSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionTypeSetObj(Decomposition,DecompositionType,CMISSDecompositionTypeSetPtrC)
      ELSE
        CMISSDecompositionTypeSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionTypeSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeSetPtrC

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a mesh for a mesh identified by a user number for C.
  FUNCTION CMISSMeshCreateFinishNumberC(RegionUserNumber,MeshUserNumber) BIND(C, NAME = "CMISSMeshCreateFinishNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to finish the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to finish the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateFinishNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshCreateFinishNumber(RegionUserNumber,MeshUserNumber,CMISSMeshCreateFinishNumberC)

    RETURN

  END FUNCTION CMISSMeshCreateFinishNumberC

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a mesh for a mesh identified by an object.
  FUNCTION CMISSMeshCreateFinishPtrC(MeshPtr) BIND(C, NAME = "CMISSMeshCreateFinishObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to finish creating.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateFinishPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshCreateFinishPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshCreateFinishObj(Mesh,CMISSMeshCreateFinishPtrC)
      ELSE
        CMISSMeshCreateFinishPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshCreateFinishPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshCreateFinishPtrC

  !
  !================================================================================================================================
  !

  !>Starts the creation of a mesh for a mesh identified by a user number for C.
  FUNCTION CMISSMeshCreateStartNumberC(MeshUserNumber,RegionUserNumber,NumberOfDimensions) BIND(C, NAME = &
    &"CMISSMeshCreateStartNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDimensions !<The number of dimensions for the mesh, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateStartNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshCreateStartNumber(MeshUserNumber,RegionUserNumber,NumberOfDimensions,CMISSMeshCreateStartNumberC)

    RETURN

  END FUNCTION CMISSMeshCreateStartNumberC

  !
  !================================================================================================================================
  !

  !>Starts the creation of a mesh for a mesh identified by an object for C.
  FUNCTION CMISSMeshCreateStartPtrC(MeshUserNumber,RegionPtr,NumberOfDimensions,MeshPtr) BIND(C, NAME = "CMISSMeshCreateStartObj")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to start the creation of, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr !<C pointer to the region containing the mesh to start the creation of.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDimensions !<The number of dimensions for the mesh, for C.
    TYPE(C_PTR), INTENT(OUT) :: MeshPtr !<C pointer to the created mesh.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateStartPtrC !<Error Code.
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshCreateStartPtrC = CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr, Region)
      IF(ASSOCIATED(Region)) THEN
        IF(C_ASSOCIATED(MeshPtr)) THEN
          CALL C_F_POINTER(MeshPtr, Mesh)
          IF(ASSOCIATED(Mesh)) THEN
            CALL CMISSMeshCreateStartObj(MeshUserNumber,Region,NumberOfDimensions,Mesh,CMISSMeshCreateStartPtrC)
          ELSE
            CMISSMeshCreateStartPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshCreateStartPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshCreateStartPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshCreateStartPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshCreateStartPtrC

  !
  !================================================================================================================================
  !

  !>Destroys a mesh identified by a user number for C.
  FUNCTION CMISSMeshDestroyNumberC(RegionUserNumber,MeshUserNumber) BIND(C, NAME = "CMISSMeshDestroyNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to destroy, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to destroy, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshDestroyNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshDestroyNumber(RegionUserNumber,MeshUserNumber,CMISSMeshDestroyNumberC)

    RETURN

  END FUNCTION CMISSMeshDestroyNumberC

  !
  !================================================================================================================================
  !

  !>Destroys a mesh identified by an object for C.
  FUNCTION CMISSMeshDestroyPtrC(MeshPtr) BIND(C, NAME = "CMISSMeshDestroyObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to destroy.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshDestroyPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshDestroyPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshDestroyObj(Mesh,CMISSMeshDestroyPtrC)
      ELSE
        CMISSMeshDestroyPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshDestroyPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshDestroyPtrC

  !
  !================================================================================================================================
  !

  !>Returns the number of components in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfComponentsGetNumberC(RegionUserNumber,MeshUserNumber,NumberOfComponents) BIND(C, NAME = &
    &"CMISSMeshNumberOfComponentsGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the number of components for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents !<The number of components in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsGetNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfComponentsGetNumber(RegionUserNumber,MeshUserNumber,NumberOfComponents,&
      &CMISSMeshNumberOfComponentsGetNumberC)

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the number of components in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfComponentsGetPtrC(MeshPtr,NumberOfComponents) BIND(C, NAME = "CMISSMeshNumberOfComponentsGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to get the number of components for.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents !<The number of components in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfComponentsGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfComponentsGetObj(Mesh,NumberOfComponents,CMISSMeshNumberOfComponentsGetPtrC)
      ELSE
        CMISSMeshNumberOfComponentsGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfComponentsGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of components in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfComponentsSetNumberC(RegionUserNumber,MeshUserNumber,NumberOfComponents) BIND(C, NAME =&
    & "CMISSMeshNumberOfComponentsSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsSetNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfComponentsSetNumber(RegionUserNumber,MeshUserNumber,NumberOfComponents,&
      &CMISSMeshNumberOfComponentsSetNumberC)

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of components in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfComponentsSetPtrC(MeshPtr,NumberOfComponents) BIND(C, NAME = "CMISSMeshNumberOfComponentsSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to set the number of components to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfComponentsSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfComponentsSetObj(Mesh,NumberOfComponents,CMISSMeshNumberOfComponentsSetPtrC)
      ELSE
        CMISSMeshNumberOfComponentsSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfComponentsSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the number of elements in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfElementsGetNumberC(RegionUserNumber,MeshUserNumber,NumberOfElements) BIND(C, NAME = &
    &"CMISSMeshNumberOfElementsGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the number of elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the number of elements for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements !<The number of elements in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsGetNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfElementsGetNumber(RegionUserNumber,MeshUserNumber,NumberOfElements,CMISSMeshNumberOfElementsGetNumberC)

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the number of elements in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfElementsGetPtrC(MeshPtr,NumberOfElements) BIND(C, NAME = "CMISSMeshNumberOfElementsGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to get the number of elements for.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements !<The number of elements in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfElementsGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfElementsGetObj(Mesh,NumberOfElements,CMISSMeshNumberOfElementsGetPtrC)
      ELSE
        CMISSMeshNumberOfElementsGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfElementsGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfElementsSetNumberC(RegionUserNumber,MeshUserNumber,NumberOfElements) BIND(C, NAME = &
    &"CMISSMeshNumberOfElementsSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the number of elements to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the number of elements to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfElements !<The number of elements in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsSetNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfElementsSetNumber(RegionUserNumber,MeshUserNumber,NumberOfElements,CMISSMeshNumberOfElementsSetNumberC)

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfElementsSetPtrC(MeshPtr,NumberOfElements) BIND(C, NAME = "CMISSMeshNumberOfElementsSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to set the number of elements to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfElements !<The number of elements in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfElementsSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfElementsSetObj(Mesh,NumberOfElements,CMISSMeshNumberOfElementsSetPtrC)
      ELSE
        CMISSMeshNumberOfElementsSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfElementsSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsSetPtrC

  !
  !================================================================================================================================
  !

  !>Finishes creating elements for a mesh component of a mesh identified by a user number for C.
  FUNCTION CMISSMeshElementsCreateFinishNumberC(RegionUserNumber,MeshUserNumber,MeshComponentNumber) BIND(C, NAME = &
    &"CMISSMeshElementsCreateFinishNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to finish creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to finish creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to finish creating the elements for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateFinishNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshElementsCreateFinishNumber(RegionUserNumber,MeshUserNumber,MeshComponentNumber,&
      &CMISSMeshElementsCreateFinishNumberC)

    RETURN

  END FUNCTION CMISSMeshElementsCreateFinishNumberC

  !
  !================================================================================================================================
  !

  !>Finishes creating elements for a mesh component of a mesh identified by an object for C.
  FUNCTION CMISSMeshElementsCreateFinishPtrC(MeshElementsPtr) BIND(C, NAME = "CMISSMeshElementsCreateFinishObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<C pointer the mesh elements to finish creating.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateFinishPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements

    CMISSMeshElementsCreateFinishPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        CALL CMISSMeshElementsCreateFinishObj(MeshElements,CMISSMeshElementsCreateFinishPtrC)
      ELSE
        CMISSMeshElementsCreateFinishPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsCreateFinishPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsCreateFinishPtrC

  !
  !================================================================================================================================
  !

  !>Starts creating elements for a mesh component of a mesh identified by a user number for C.
  FUNCTION CMISSMeshElementsCreateStartNumberC(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber) BIND(C, NAME = &
    &"CMISSMeshElementsCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to start creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to start creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to start creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: BasisUserNumber !<The user number of the default basis to use for the elements, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateStartNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshElementsCreateStartNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber,&
      &CMISSMeshElementsCreateStartNumberC)

    RETURN

  END FUNCTION CMISSMeshElementsCreateStartNumberC

  !
  !================================================================================================================================
  !

  !>Starts creating elements for a mesh component of a mesh identified by an object for C.
  FUNCTION CMISSMeshElementsCreateStartPtrC(MeshPtr,MeshComponentNumber,BasisPtr,MeshElementsPtr) BIND(C, NAME = &
    &"CMISSMeshElementsCreateStartObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to start the creation of elements for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to start creating the elements for, for C.
    TYPE(C_PTR), INTENT(IN) :: BasisPtr !<C pointer to the default basis to use for the elements.
    TYPE(C_PTR), INTENT(OUT) :: MeshElementsPtr !<C pointer to the created mesh elements.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateStartPtrC !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh
    TYPE(CMISSBasisType), POINTER :: Basis
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements

    CMISSMeshElementsCreateStartPtrC = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            IF(C_ASSOCIATED(MeshElementsPtr)) THEN
              CALL C_F_POINTER(MeshElementsPtr, MeshElements)
              IF(ASSOCIATED(MeshElements)) THEN
                CALL CMISSMeshElementsCreateStartObj(Mesh, MeshComponentNumber, Basis, MeshElements, &
                  & CMISSMeshElementsCreateStartPtrC)
              ELSE
                CMISSMeshElementsCreateStartPtrC = CMISSErrorConvertingPointer
              ENDIF
            ELSE
              CMISSMeshElementsCreateStartPtrC = CMISSPointerIsNULL
            ENDIF
          ELSE
            CMISSMeshElementsCreateStartPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshElementsCreateStartPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshElementsCreateStartPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsCreateStartPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsCreateStartPtrC

  !
  !================================================================================================================================
  !

  !>Returns the basis for an element in a mesh identified by an user number for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsBasisGetNumberC(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber, &
      & BasisUserNumber) BIND(C, NAME = "CMISSMeshElementsBasisGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the basis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the basis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to get the basis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to get the basis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: BasisUserNumber !<The user number of the basis for the element, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsBasisGetNumberC !<Error Code.
    !Local variable

    CALL CMISSMeshElementsBasisGetNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber,BasisUserNumber,&
      &CMISSMeshElementsBasisGetNumberC)

    RETURN

  END FUNCTION CMISSMeshElementsBasisGetNumberC


  !
  !================================================================================================================================
  !

  !>Sets/changes the basis for an element in a mesh identified by an object. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsBasisSetCPtr(MeshElementsPtr,GlobalElementNumber,BasisPtr) BIND(C, NAME = " CMISSMeshElementsBasisSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<C pointer to the mesh elements to set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to set the basis for, for C.
    TYPE(C_PTR), INTENT(IN) :: BasisPtr !<C pointer to the basis for the element to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsBasisSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements
    TYPE(CMISSBasisType), POINTER :: Basis

    CMISSMeshElementsBasisSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSMeshElementsBasisSet(MeshElements, GlobalElementNumber, Basis, &
              & CMISSMeshElementsBasisSetCPtr)
          ELSE
            CMISSMeshElementsBasisSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshElementsBasisSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshElementsBasisSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsBasisSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsBasisSetCPtr

    !
  !================================================================================================================================
  !

  !>Returns the element nodes for an element in a mesh identified by an user number. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsNodesGetCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber, &
    & ElementUserNodes) BIND(C, NAME = "CMISSMeshElementsNodesGet")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the element nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the element nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to get the element nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to get the element nodes for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ElementUserNodes !<The user node number number of the i'th element node.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsNodesGetCNum !<Error Code.
    !Local variable

    CALL CMISSMeshElementsBasisGetNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber,ElementUserNodes, &
      &CMISSMeshElementsNodesGetCNum)

    RETURN

  END FUNCTION CMISSMeshElementsNodesGetCNum







!!==================================================================================================================================
!!
!! REGION_ROUTINES
!!
!!==================================================================================================================================

  !>Finishes the process of creating a region identified by a pointer for C.
  FUNCTION CMISSRegionCreateFinishPtrC(RegionPtr) BIND(C,NAME="CMISSRegionCreateFinish")

    !Argument variablesCMISSAnalyticAnalysisOutputNumber
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateFinishPtrC !<Error Code.
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region

    CMISSRegionCreateFinishPtrC=CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr,Region)
      IF(ASSOCIATED(Region)) THEN        
        CALL CMISSRegionCreateFinish(Region,CMISSRegionCreateFinishPtrC)
      ELSE
        CMISSRegionCreateFinishPtrC=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionCreateFinishPtrC=CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionCreateFinishPtrC

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a region identified by an user number for C.
  FUNCTION CMISSRegionCreateFinishCNum(RegionUserNumber) BIND(C,NAME="CMISSRegionCreateFinishNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateFinishCNum
    !Local variables

    CALL CMISSRegionCreateFinish(RegionUserNumber,CMISSRegionCreateFinishCNum)

    RETURN
    
  END FUNCTION CMISSRegionCreateFinishCNum

  !
  !================================================================================================================================
  !
  
  !>Starts the process of creating a region identified by a pointer for C.
  FUNCTION CMISSRegionCreateStartPtrC(RegionUserNumber,ParentRegionPtr,RegionPtr) BIND(C,NAME="CMISSRegionCreateStart")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    TYPE(C_PTR), VALUE, INTENT(IN) :: ParentRegionPtr
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateStartPtrC
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region,ParentRegion

    CMISSRegionCreateStartPtrC=CMISSNoError
    IF(C_ASSOCIATED(ParentRegionPtr)) THEN
      CALL C_F_POINTER(ParentRegionPtr,ParentRegion)
      IF(ASSOCIATED(ParentRegion)) THEN        
        IF(C_ASSOCIATED(RegionPtr)) THEN
          CALL C_F_POINTER(RegionPtr,Region)
          IF(ASSOCIATED(Region)) THEN        
            CALL CMISSRegionCreateStart(RegionUserNumber,ParentRegion,Region,CMISSRegionCreateStartPtrC)
          ELSE
            CMISSRegionCreateStartPtrC=CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSRegionCreateStartPtrC=CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSRegionCreateStartPtrC=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionCreateStartPtrC=CMISSPointerIsNULL
    ENDIF
    
    RETURN
    
  END FUNCTION CMISSRegionCreateStartPtrC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a region identified by an user number for C.
  FUNCTION CMISSRegionCreateStartCNum(RegionUserNumber,ParentRegionUserNumber) BIND(C,NAME="CMISSRegionCreateStartNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: ParentRegionUserNumber
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateStartCNum
    !Local variables

    CALL CMISSRegionCreateStart(RegionUserNumber,ParentRegionUserNumber,CMISSRegionCreateStartCNum)

    RETURN
    
  END FUNCTION CMISSRegionCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a region identified by a pointer for C.
  FUNCTION CMISSRegionLabelGetCPtr(RegionPtr,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelGet")
  
    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelGetCPtr
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSRegionType), POINTER :: Region    

    CMISSRegionLabelGetCPtr=CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr,Region)
      IF(ASSOCIATED(Region)) THEN        
        CALL CMISSRegionLabelGet(Region,FLabel,CMISSRegionLabelGetCPtr)
        CALL CMISSF2CString(Flabel,Label)
      ELSE
        CMISSRegionLabelGetCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionLabelGetCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionLabelGetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a region identified by an user number for C.
  FUNCTION CMISSRegionLabelGetCNum(RegionUserNumber,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelGetNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelGetCNum
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
 
    CALL CMISSRegionLabelGet(RegionUserNumber,FLabel,CMISSRegionLabelGetCNum)
    CALL CMISSF2CString(Flabel,Label)
 
    RETURN
    
  END FUNCTION CMISSRegionLabelGetCNum

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the label for a region identified by a pointer for C.
  FUNCTION CMISSRegionLabelSetCPtr(RegionPtr,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelSet")
  
    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelSetCPtr
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSRegionType), POINTER :: Region    

    CMISSRegionLabelSetCPtr=CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr,Region)
      IF(ASSOCIATED(Region)) THEN        
        CALL CMISSC2FString(Label,Flabel)
        CALL CMISSRegionLabelSet(Region,FLabel,CMISSRegionLabelSetCPtr)
      ELSE
        CMISSRegionLabelSetCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionLabelSetCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionLabelSetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the label for a region identified by an user number for C.
  FUNCTION CMISSRegionLabelSetCNum(RegionUserNumber,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelSetNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelSetCNum
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
 
    CALL CMISSC2FString(Label,Flabel)
    CALL CMISSRegionLabelSet(RegionUserNumber,FLabel,CMISSRegionLabelSetCNum)
    
    RETURN
    
  END FUNCTION CMISSRegionLabelSetCNum

  !
  !================================================================================================================================
  !

END MODULE OPENCMISS_C
